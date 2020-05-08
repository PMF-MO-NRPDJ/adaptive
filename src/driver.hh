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
//#include <dune/pdelab/localoperator/defaultimp.hh>
//#include <dune/pdelab/localoperator/flags.hh>
//#include <dune/pdelab/localoperator/pattern.hh>
//#include <dune/pdelab/localoperator/variablefactories.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include "operator.hh"
#include "estimator.hh"
#include "bctype.hh"
#include "norm.hh"

template <typename Grid>
void driver(Grid &grid, int subsampling,  int steps,
            double alpha, double beta, double tol, std::string output)
{
  using GV = typename Grid::LeafGridView;
  GV gv = grid.leafGridView();  // Drži pokazivač na grid.
  const int dim = GV::dimension;

//1. Konstrukcija Grid Operatora
  const int k = 1;  // stupanj polinom u prostoru KE rješenja
  using DF = typename GV::Grid::ctype; // tip za koordinate
  using RF = double;                   // tip za računanje

  using FEM = Dune::PDELab::PkLocalFiniteElementMap<GV, DF, RF, k>;
  FEM fem(gv);
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using GFS = Dune::PDELab::GridFunctionSpace<GV, FEM, CON, VBE>;
  GFS gfs(gv, fem);

  DirichletBdry bctype; // identifikacija Dirichletove granice
  using CC = typename GFS::template ConstraintsContainer<RF>::Type;
  CC cc;
  Dune::PDELab::constraints(bctype, gfs, cc);

  // Vektor stupnjeva slobode rješenja
  using U = Dune::PDELab::Backend::Vector<GFS, RF>;
  U  u(gfs, 0.0), utmp(gfs);
  int int_order = 2 * (k + 1); // Red integracijske formule za inverz matrice mase.
  BCExtension<GV> bcext(gv);
  // Izračunaj z iz rubnog uvjeta
  Dune::PDELab::interpolate(bcext, gfs, utmp);       // sad praktički imamo egzaktno rješenje u ztmp.
  Dune::PDELab::copy_constrained_dofs(cc, utmp, u);  // copy ztmp --> z
  // Lokalni operator
  using LOP = LocalOperator<DirichletBdry,FEM>;
  LOP lop(bctype);

  using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
  MBE mbe( static_cast<int>(std::pow(1 + 2 * k, dim)) );
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, CC, CC>;
  GO go(gfs, cc, gfs, cc, lop, mbe);

  using LS = Dune::PDELab::ISTLBackend_SEQ_CG_AMG_SSOR<GO>;
  const int max_no_iter = 1000;
  const int verbosity = 0;
  LS ls(max_no_iter, verbosity);
  using  SLP =Dune::PDELab::StationaryLinearProblemSolver<GO,LS,U>;

// 2. Konstrukcija procjenitelja greške
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, RF, dim>;
  auto eltype = Dune::GeometryType(Dune::GeometryTypes::simplex(dim));
  // Konstruktor za P0 elemente uzima tip elementa a ne gid view.
  P0FEM p0fem(eltype);
  using NCON   = Dune::PDELab::NoConstraints;
  using P0GFS  = Dune::PDELab::GridFunctionSpace<GV, P0FEM, NCON, VBE>;
  using ESTLOP = Estimator<DirichletBdry, FEM>;
  using NCC = typename P0GFS::template ConstraintsContainer<RF>::Type;
  // Procjenitelj računamo tako da izračunamo rezidual mrežnog operatora ESTGO
  using ESTGO = Dune::PDELab::GridOperator<GFS,    // prostor rješenja
                                           P0GFS,  // prostor test funkcija
                                           ESTLOP, // lokalni operator
                                           MBE, RF, RF, RF,  // uobičajeno
                                           NCC, NCC          // bez ograničenja
                                           >;
  using Z = Dune::PDELab::Backend::Vector<P0GFS, RF>;


  // VTK ispis
  using UDGF  = Dune::PDELab::DiscreteGridFunction<GFS, U>;
  using VTKU  = Dune::PDELab::VTKGridFunctionAdapter<UDGF>;

  using VTKE  = Dune::PDELab::VTKGridFunctionAdapter<ExactGF<GV>>;

  using ZDGF  = Dune::PDELab::DiscreteGridFunction<P0GFS, Z>;
  using VTKZ  = Dune::PDELab::VTKGridFunctionAdapter<ZDGF>;

  ExactGF<GV> exact(gv);

  // Petlja unutar koje adaptiramo mrežu
  for(int i = 0; i < steps; ++i){
    std::string iter = std::to_string(i);

    std::cout << "===== Iteracija no. : " << iter
              << "  Najviše profinjenje mreže: " << grid.maxLevel() << std::endl;
    std::cout << "      Br vezanih stupnjeva slobode = "  << cc.size() << " od " << gfs.globalSize()
              << std::endl;

    // Rješavanje zadaće.
    // Solver se mora "obnavljati" nakon promjene mreže
    double redukcija = 1e-10;
    SLP slp(go, ls, u, redukcija);
    slp.apply();

    // Izračunaj stvarnu L2 grešku, poslije rješavanja sustava.
    U  err(gfs);
    Dune::PDELab::interpolate(exact, gfs, err);
    err -= u;
    UDGF errdgf(gfs, err);
    std::cout << "      Stvarna L2 greška = " << L2norm(gv,errdgf) << std::endl;

    // Konstrukcija procijenitelja greške. U svakom koraku konstruiramo novi procjenitelj.
    P0GFS  p0gfs(gv, p0fem);   // pogfs ne profinjujemo već iznova konstruiramo
    ESTLOP estlop(bctype);
    ESTGO  estgo(gfs, p0gfs, estlop, mbe);
    Z     z(p0gfs, 0.0);
    estgo.residual(u, z); // z = residual(u)
    auto estimated_error = sqrt(z.one_norm());
    std::cout << "      L2 norma greške je procijenjena na " << estimated_error << std::endl;

    // VTK ispis
    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, Dune::RefinementIntervals{subsampling});
    UDGF udgf(gfs, u);
    vtkwriter.addVertexData(std::shared_ptr<VTKU>(new VTKU(udgf, "solution")));
    ZDGF zdgf(p0gfs, z);
    vtkwriter.addCellData(std::shared_ptr<VTKZ>(new VTKZ(zdgf, "indicator")));
    vtkwriter.addVertexData(std::shared_ptr<VTKE>(new VTKE(exact, "exact")));
    vtkwriter.write(output + iter, Dune::VTK::ascii);

    // Da li je greška dovoljno mala ?
    if (estimated_error <= tol)
      break;  // prekini profinjavanje
    if (i == steps - 1)
      break; // u zadnjem koraku preskoči adaptaciju mreže

    // Greška ne zadovoljava toleranciju. Označi elemente za profinjenje.
    RF eta_refine, eta_coarsen;
    Dune::PDELab::error_fraction(z, alpha, beta, eta_refine, eta_coarsen, 1 /* verbose */);

    // označi za profinjenje
    //eta_coarsen = 0; // ako ne želimo okrupnjavanje
    int min_level =2;
    int max_level = 100;
    int verbosity = 1;
    Dune::PDELab::mark_grid(grid, z, eta_refine, eta_coarsen, min_level, max_level, verbosity);
    // profini mrežu i interpoliraj vektor rješenja
    Dune::PDELab::adapt_grid(grid, gfs, u, int_order); // Profini mrežu, gfs i vektor u
    // ponovo izračunaj Dirichletova ograničenja
    Dune::PDELab::constraints(bctype, gfs, cc);
    // korektne rubne uvjete upiši u novi vektor
    U  unew(gfs);
    Dune::PDELab::interpolate(bcext, gfs, unew);
    // kopiraj Dirichletove rubne uvjete u interpolirani vektor rješenja
    Dune::PDELab::copy_constrained_dofs(cc, unew, u);  // unew --> u
  }
}

