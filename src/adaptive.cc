// Adaptivna metoda konačnih elemenata. Primjer je uzet iz
// dune-pdelab-tutorials/tutorial05.
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cmath>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/grid/uggrid.hh>

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include "driver.hh"

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);

  // Otvaranje datoteke parametara
  Dune::ParameterTree ptree;
  Dune::ParameterTreeParser ptreeparser;
  ptreeparser.readINITree("adaptive.ini", ptree);  // pročitaj parametre iz datoteke
  ptreeparser.readOptions(argc, argv, ptree);      // dodaj parametre komandne linije

  const int dim = 2;
  std::string filename = "ldomain.msh";

  // Čitanje podataka iz ulazne datoteke.
  // subsampling za ispis rješenja i procjenitelja
  int subsampling    = ptree.get<int>("output.subsampling", 1);
  // ime izlazne datoteke
  std::string output = ptree.get<std::string>("output.filename", "output");
  // Tražena tolerancija za procjenitelj
  double  tol        = ptree.get<double>("fem.tol", 1E-3);
  // alpha i beta - parametri algoritma za profinjenje mreže
  double alpha    = ptree.get<double>("fem.alpha", 0.5);
  double beta     = ptree.get<double>("fem.beta", 0.1);
  // Broj profinjenja mreže
  int steps          = ptree.get<int>("fem.steps", 5);

  // alpha = 1 uniformno profinjenje, alpha = 0 nema profinjenja
  if(alpha > 1.0 or alpha <0)
      throw std::runtime_error("alpha out of bounds!");
  // beta = 1 uniformno okrupnjavanje, beta = 0 nema okrupnjavanja
  if(beta > 1.0 or beta <0)
      throw std::runtime_error("beta out of bounds!");
  if(alpha+beta > 1.0)
      throw std::runtime_error("alpha + beta > 1!");

  // Primjer korištenja UG Grida.
  //  using Grid = Dune::UGGrid<dim>;
  //  Grid * gridp = Dune::GmshReader<Grid>::read(filename);
  //  driver<Grid>(*gridp, subsampling,  steps, fraction, tol, output);

  // Primjer korištenja ALUGrid-a.
  using Grid = Dune::ALUGrid<dim, dim, Dune::simplex, Dune::conforming>;
  std::unique_ptr<Grid> gridp{Dune::GmshReader<Grid>::read(filename)};

  driver(*gridp, subsampling, steps, alpha, beta, tol, output);

  return 0;
}
