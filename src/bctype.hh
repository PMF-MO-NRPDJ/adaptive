#pragma once

#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/common/function.hh>


#include "coefficients.hh"

// Klasa mora proširivati klasu PDELab::DirichletConstraintsParameters
// i u njoj prerađuje metodu isDirichlet() odlučuje je li neka točka
// na Dirichletovoj granici ili nije.
class DirichletBdry : public Dune::PDELab::DirichletConstraintsParameters
{
public:
  //  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
  //  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
  //  povratna vrijednost: true ako je točka na Dirichletovoj granici
  //                       false ako nije.
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::coorddimension-1> & coord
                   ) const
  {
//      auto x = intersection.geometry().global(coord);
//      if(std::abs(x[1]) < 1E-8 and x[0] >= 0.0) // y=0, x>=0
//          return false;
      return true;
  }

};

/*  Klasa koja određuje vrijednost Dirichletovog rubnog uvjeta i
    njegovo proširenje na čitavu domenu.
    Template parametri:
       GV = GridView
       RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)

       Treći parametar u GridFunctionTraits je dimenzija slike funkcije (1 jer su
       naše funkcije skalarne). Ta se dimenzija ponavlja u zadnjem parametru
       GridFunctionTraits klase.
    */
template<typename GV>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >,
                                      BCExtension<GV> > {
 // Klasa čuva referencu na GridView objekt.
  GV const & gv;
public:
  using Traits = Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >;

  // Konstruktor samo uzima referencu na  GridView objekt.
  BCExtension (const GV& gv_) : gv(gv_) {}

  // Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na
  // Dirichletovoj granici, nda funkcija daje proširenje Dirichletovog rubnog
  // uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno.
  // e     = element
  // xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
  // y      = izračunata Dirichletova vrijednost
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
      auto x = e.geometry().global(xlocal);     
      y =  exact(x);
  }

  // Vrati referencu na GridView
  inline const GV& getGridView () {return gv;}
};




