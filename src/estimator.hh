#pragma once
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>


//#include<dune/pdelab/common/referenceelements.hh> 2.6.0
#include <dune/geometry/referenceelements.hh>
//#include<dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "coefficients.hh"


/** Rezidualni procjenitelj za eliptičku zadaću:
 *
 * \f{align*}{
 *   -\Delta u(x) + a u(x) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 * Za izračunati procjenitelj treba pozvati residual() na mrežnom
 * operatoru. Lokalno se računa \f$\eta_K^2\f$ na svakom elementu \f$K\f$.
 *
 * Pretpostavke i ograničenja:
 * - Uzima se da je  LFSU jednak P_k/Q_k te da je
 *   LFSV jednak P_0.
 * - Derivacije drugog reda se zanemaruju.
 *
 */
template<typename BCType, typename FEM>
class Estimator : public Dune::PDELab::LocalOperatorDefaultFlags
{
  using LocalBasis = typename FEM::Traits::FiniteElementType::Traits::LocalBasisType;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  BCType& bctype; // parameter functions

  // Metoda za računanje dijametra ćelije
  template<class GEO>
  double diameter (const GEO& geo) const
  {
    double hmax = -1.0E00;
    for (int i=0; i<geo.corners(); i++)
      {
        auto xi = geo.corner(i);
        for (int j=i+1; j<geo.corners(); j++)
          {
            auto xj = geo.corner(j);
            xj -= xi;
            hmax = std::max(hmax,xj.two_norm());
          }
      }
    return hmax;
  }

public:
  // pattern assembly flags
  enum { doPatternVolume = false };
  enum { doPatternSkeleton = false };

  // residual assembly flags
  enum { doAlphaVolume  = true };
  enum { doAlphaSkeleton  = true };
  enum { doAlphaBoundary  = true };

  Estimator (BCType& bctype_) : bctype(bctype_)
  {}

  // volumni dio indikatora
  // h_K^2 int_K (f - au + \Delta u)^2 dx. Zanemarujemo druge derivacije pa računamo
  // h_K^2 int_K (f - au)^2 dx .
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim  = EG::Geometry::mydimension;
    const auto geo = eg.geometry();
    const int order = 2*lfsu.finiteElement().localBasis().order();  // red integracije je 2k
    const auto & rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),order);

    double sum = 0.0;
    for (const auto& ip : rule)
      {
        auto& phihat = cache.evaluateFunction(ip.position(),lfsu.finiteElement().localBasis());

        // rješenje u
        double u=0.0;
        for (size_t i=0; i<lfsu.size(); i++) u += x(lfsu,i)*phihat[i];

        // slobodni član
        auto a = reactCoeff();

        // desna strana
        auto f = RHS(eg.geometry().global(ip.position()));

        // Prostorni rezidual. Za P1 elemente laplas rješenja je jednak nuli.
        double factor =ip.weight() * geo.integrationElement(ip.position());
        sum += (f-a*u)*(f-a*u)*factor;
      }

    auto h_T = diameter(eg.geometry());
    r.accumulate(lfsv, 0, h_T*h_T * sum);
  }

  // Dio indikatora po unutarnjim stranicama elemenata.
  // Svaka se stranica obilazi samo jednom.
  // Računamo (h/2) int_sigma [grad u . n]^2 ds
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_skeleton(const IG& ig,
      const LFSU & lfsu_i, const X & x_i, const LFSV & lfsv_i,  // inside
      const LFSU & lfsu_o, const X & x_o, const LFSV & lfsv_o,  // outside
      R& r_i, R& r_o) const
  {
    // geometrije stranice u lokalnim koordinatama elemenata
    auto insidegeo  = ig.geometryInInside();
    auto outsidegeo = ig.geometryInOutside();

    // elementi inside i outside
    auto cell_inside  = ig.inside();
    auto cell_outside = ig.outside();

    // geometrije elemenata u globalnim koordinatama
    auto geo_i = cell_inside.geometry();
    auto geo_o = cell_outside.geometry();

    const int dim = IG::Entity::dimension;  // dimenzija prostora

    auto globalgeo = ig.geometry();
    const int order = 2*lfsu_i.finiteElement().localBasis().order();
    //auto rule = Dune::PDELab::quadratureRule(globalgeo,order); -- jednostavniji način
    const auto & rule = Dune::QuadratureRules<double,dim-1>::rule(globalgeo.type(),order);

    double sum = 0.0;
    for (const auto& ip : rule){
        auto xi = ip.position();
        // pozicija kvadraturne točke u lokalnim koordinatama elementa.
        // Potrebna je za račun baznih funkcija.
        auto iplocal_i = insidegeo.global(xi);
        auto iplocal_o = outsidegeo.global(xi);

        const auto n_F = ig.unitOuterNormal(xi);

        // grad u . n na inside elementu
        auto& gradphihat_i = cache.evaluateJacobian(
                                           iplocal_i, lfsu_i.finiteElement().localBasis()
                                                   );
        const auto S_i = geo_i.jacobianInverseTransposed(iplocal_i);
        double gradun_i = 0.0; // gradu . n
        for (size_t k=0; k<lfsu_i.size(); k++)
          {
            Dune::FieldVector<double,dim> v;
            S_i.mv(gradphihat_i[k][0],v); // v = S_i*gradphihat_i[k][0] = gradphi_i[k]
            gradun_i += x_i(lfsu_i,k)*(v*n_F);
          }

        // grad u . n na outside elementu
        auto& gradphihat_o = cache.evaluateJacobian(
                                            iplocal_o,lfsu_o.finiteElement().localBasis()
                                                   );
        const auto S_o = geo_o.jacobianInverseTransposed(iplocal_o);
        double gradun_o = 0.0; // grad u . n
        for (size_t k=0; k<lfsu_o.size(); k++)
          {
            Dune::FieldVector<double,dim> v;
            S_o.mv(gradphihat_o[k][0],v); // v = S_o*gradphihat_i[k][0] = gradphi_o[k]
            gradun_o += x_o(lfsu_o,k)*(v*n_F);
          }

        // integracija
        double factor = ip.weight() * globalgeo.integrationElement(xi);
        double jump = gradun_i-gradun_o;
        sum += jump*jump*factor;
      }

    // akumulacija indikatora
    auto h_T = diameter(globalgeo);
    r_i.accumulate(lfsv_i, 0, 0.5*h_T * sum);
    r_o.accumulate(lfsv_o, 0, 0.5*h_T * sum);
  }

  // Doprinos rezidualu od ruba domene. Dirichletovu granicu preskačemo.
  // h_k int_sigma (j+ grad u . n)^2 ds
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary(const IG& ig,
                       const LFSU& lfsu_i, const X& x_i, const LFSV& lfsv_i, // postoji samo inside
                       R& r_i) const
  { 
    // geometrija stranice u lokalnim koordinatama elementa
    auto insidegeo = ig.geometryInInside();

    // "inside" element (element koji sadrži stranicu)
    auto cell_inside = ig.inside();

    //  geometrije elementa koji sadrži stranicu u globalnim koordinatama
    auto geo_i = cell_inside.geometry();

    const int dim = IG::Entity::dimension; // Prostorna dimenzija

    // geometrija stranice u globalnim koordinatama
    auto globalgeo = ig.geometry();
    // Kvadraturna formula -- jednostavniji način
    const int order = 2*lfsu_i.finiteElement().localBasis().order();
    const auto & rule = Dune::PDELab::quadratureRule(globalgeo,order);

    double sum = 0.0;
    for(const auto& ip : rule){
        auto xi = ip.position();
        // preskoči Dirichletovu granicu
        if( bctype.isDirichlet(ig,xi) ) continue;

        // pozicija kvadraturne točke u elementu
        auto iplocal_i = insidegeo.global(xi);

        auto n_F = ig.unitOuterNormal(xi);

        // grad u . n
        auto& gradphihat_i = cache.evaluateJacobian(
                                    iplocal_i,lfsu_i.finiteElement().localBasis()
                                                   );
        const auto S_i = geo_i.jacobianInverseTransposed(iplocal_i);
        double gradun_i = 0.0;  // grad u . n
        for (size_t k=0; k<lfsu_i.size(); k++)
          {
            Dune::FieldVector<double,dim> v;
            S_i.mv(gradphihat_i[k][0],v);  // v=S_i*gradphihat_i[k][0]
            gradun_i += x_i(lfsu_i,k)*(v*n_F);
          }

        // Neumannov rubni uvjet
        auto j = neumannBC(globalgeo.global(xi));

        // integracija
        double factor = ip.weight() * globalgeo.integrationElement(xi);
        double jump = j + gradun_i;
        sum += jump*jump * factor;
      }

    // akumulacija indicatora
    auto h_T = diameter(globalgeo);
    r_i.accumulate(lfsv_i, 0, h_T * sum);
  }
};
