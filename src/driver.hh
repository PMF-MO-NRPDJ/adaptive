#pragma once

#include <cmath>
#include <iostream>
#include <string>

#include <dune/common/timer.hh>

#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

// dune-istl included by pdelab
#include <dune/pdelab/adaptivity/adaptivity.hh>

#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/variablefactories.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include "operator.hh"
#include "estimator.hh"
#include "bctype.hh"


template <typename Grid>
void driver(Grid &grid, int subsampling,  int steps,
            double alpha, double tol, std::string output)
{
  using GV = typename Grid::LeafGridView;
  GV gv = grid.leafGridView();
  const int dim = GV::dimension;
  const int degree = 1;  // stupanj polinom u prostoru KE rješenja
  using DF = typename GV::Grid::ctype; // tip za koordinate
  using RF = double;                   // tip za računanje

  typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF, RF, degree> FEM;
  FEM fem(gv);
  typedef Dune::PDELab::ConformingDirichletConstraints CON;
  typedef Dune::PDELab::ISTL::VectorBackend<> VBE;
  typedef Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE> GFS;
  GFS gfs(gv, fem);


  DirichletBdry bctype; // identifikacija Dirichletove granice
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cc;
  Dune::PDELab::constraints(bctype, gfs, cc);

  // Vektor stupnjeva slobode rješenja
  using Z = Dune::PDELab::Backend::Vector<GFS, RF>;
  Z z(gfs);
  BCExtension<GV> bcext(gv);
  // Izračunaj z iz rubnog uvjeta
  Dune::PDELab::interpolate(bcext, gfs, z);

  // Petlja unutar koje adaptiramo mrežu
  for(int i = 0; i < steps; ++i)
  {
    std::string iter = std::to_string(i);

    std::cout << "===== Iteracija no. : " << iter
              << "  Najviše profinjenje mreže: " << grid.maxLevel() << std::endl;
    std::cout << "      Br vezanih stupnjeva slobode = "  << cc.size() << " od " << gfs.globalSize()
              << std::endl;

    // Lokalni operator
    using LOP = LocalOperator<DirichletBdry,FEM>;
    LOP lop(bctype);

    using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
    MBE mbe( static_cast<int>(std::pow(1 + 2 * degree, dim)) );
    using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, CC, CC>;
    GO go(gfs, cc, gfs, cc, lop, mbe);

    using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO>;
    LS ls(100, 0);

    using  SLP =Dune::PDELab::StationaryLinearProblemSolver<GO,LS,Z>;
    SLP slp(go, ls, z, 1e-10 /* =redukcija */, 1e-99 /* =min_defect */, 0 /* =verbosity */);
    slp.apply();

    // Konstrukcija procijenitelja greške.
    using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, RF, dim>;
    auto eltype = Dune::GeometryType(Dune::GeometryTypes::simplex(dim));
    // Konstruktor za P0 elemente uzima tip elementa a ne gid view.
    P0FEM p0fem(eltype);
    using NCON = Dune::PDELab::NoConstraints;
    using P0GFS = Dune::PDELab::GridFunctionSpace<GV, P0FEM, NCON, VBE>;
    P0GFS p0gfs(gv, p0fem);
    using ESTLOP = Estimator<DirichletBdry, FEM>;
    ESTLOP estlop(bctype);
    using NCC = typename P0GFS::template ConstraintsContainer<RF>::Type;
    using ESTGO = Dune::PDELab::GridOperator<GFS,    // prostor rješenja
                                             P0GFS,  // prostor test funkcija
                                             ESTLOP, // lokalni operator
                                             MBE, RF, RF, RF,  // uobičajeno
                                             NCC, NCC          // bez ograničenja
                                             >;
    ESTGO estgo(gfs, p0gfs, estlop, mbe);

    // Procjenitelj računamo tako da izračunamo rezidual mrežnog operatora ESTGO
    using Z0 = Dune::PDELab::Backend::Vector<P0GFS, RF>;
    Z0 z0(p0gfs, 0.0);
    estgo.residual(z, z0);
    auto estimated_error = sqrt(z0.one_norm());
    std::cout << "      L2 norma greške je procijenjena na " << estimated_error << std::endl;

    // vtk output
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals{subsampling});
    typedef Dune::PDELab::DiscreteGridFunction<GFS, Z> ZDGF;
    ZDGF zdgf(gfs, z);
    typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
    vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(zdgf, "fem_sol")));
    typedef Dune::PDELab::DiscreteGridFunction<P0GFS, Z0> Z0DGF;
    Z0DGF z0dgf(p0gfs, z0);
    typedef Dune::PDELab::VTKGridFunctionAdapter<Z0DGF> VTKF0;
    vtkwriter.addCellData(std::shared_ptr<VTKF0>(new VTKF0(z0dgf, "error2")));
    vtkwriter.write(output + iter, Dune::VTK::ascii);

    // Da li je greška dovoljno mala ?
    if (estimated_error <= tol)
      break;  // prekini profinjavanje
    if (i == steps - 1)
      break; // u zadnjem koraku preskoči adaptaciju mreže

    // Greška ne zadovoljava toleranciju. Označi elemente za profinjenje.
    //demangle(z0);
    RF eta_alpha, eta_beta;
    RF beta = 0.1;
    Dune::PDELab::error_fraction(z0, alpha, beta, eta_alpha, eta_beta, 1 /* verbose */);
    if (alpha >= 1.0)
      eta_alpha = 0.0;  // uniformno profinjenje

    // označi za profinjenje
    //eta_beta = 0; // ako ne želimo okrupnjavanje
    int min_level = 2;
    int max_level = 100;
    int verbosity = 1;
    Dune::PDELab::mark_grid(grid, z0, eta_alpha, eta_beta, min_level, max_level, verbosity);
    // profini mrežu i interpoliraj vektor rješenja
    Dune::PDELab::adapt_grid(grid, gfs, z, 2 * (degree + 1));
    // ponovo izračunaj Dirichletova ograničenja
    Dune::PDELab::constraints(bctype, gfs, cc);
    // korektne rubne uvjete upiši u novi vektor
    Z znew(gfs);
    Dune::PDELab::interpolate(bcext, gfs, znew);
    // kopiraj Dirichletove rubne uvjete u interpolirani vektor rješenja
    Dune::PDELab::copy_constrained_dofs(cc, znew, z);
  }
}

