#include<vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>

#include "coefficients.hh"

/** Lokalni operator za eliptičku zadaću
 *
 * \f{align*}{
 *    -\Delta u(x) + a u(x) &=& f(x) x\in\Omega,  \\
 *                     u(x) &=& g(x) x\in\partial\Omega_D \\
 *  -\nabla u(x) \cdot n(x) &=& j(x) x\in\partial\Omega_N \\
 * \f}
 *
 */
template<typename DirichletBdry, typename FEM>
class LocalOperator :
  public Dune::PDELab::NumericalJacobianVolume<LocalOperator<DirichletBdry,FEM> >,
  public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<DirichletBdry,FEM> >,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
  DirichletBdry const & dirBdry;

public:
  enum { doPatternVolume = true };
  enum { doLambdaVolume = true };
  enum { doLambdaBoundary = true };
  enum { doAlphaVolume = true };

  LocalOperator (DirichletBdry const & dirBdry_) : dirBdry(dirBdry_) {}

  // desna strana
  template<typename EG, typename LFSV, typename R>
  void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
  {
    // kvadraturna formula
    auto geo = eg.geometry();
    const int order = 2*lfsv.finiteElement().localBasis().order();
    const auto & rule = Dune::PDELab::quadratureRule(geo,order);
    // umjesto
    // const int dim  = EG::Geometry::mydimension;
    // const auto & rule = Dune::QuadratureRules<double,dim>::rule(geo.type(),order);

    for(const auto& ip : rule){
        const auto xi = ip.position();
        const auto& phihat = cache.evaluateFunction(xi,lfsv.finiteElement().localBasis());
        const auto factor = ip.weight() * geo.integrationElement(xi);

        const auto f = RHS(geo.global(xi));
        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv, i, -f*phihat[i]*factor);
      }
  }

  // Neumannova granica
  template<typename IG, typename LFSV, typename R>
  void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
  {
    auto localgeo = ig.geometryInInside();
    // integracijska formula
    auto globalgeo = ig.geometry();
    const int order = 2*lfsv.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(globalgeo,order);

    for(const auto& ip : rule){
        auto xi = ip.position();
        // Preskoči ako je točka na Dirichletovoj granici
        if (dirBdry.isDirichlet(ig, xi))
            continue;
        // Kvadraturna točka u lokalnim koordinatama elementa.
        auto local = localgeo.global(xi);
        auto& phihat = cache.evaluateFunction(local, lfsv.finiteElement().localBasis());

        auto factor = ip.weight()*globalgeo.integrationElement(xi);
        auto j = neumannBC(globalgeo.global(xi));

        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv, i, j*phihat[i]*factor);
      }
  }

  //! bilinearni dio reziduala
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    const int dim = EG::Entity::dimension;
    using Gradient = Dune::FieldVector<double,dim>;

    // kvadraturna formula
    auto geo = eg.geometry();
    const int order = 2*lfsu.finiteElement().localBasis().order();
    auto rule = Dune::PDELab::quadratureRule(geo,order);

    for(const auto& ip : rule){
        const auto xi = ip.position();
        const auto& phihat = cache.evaluateFunction(xi, lfsu.finiteElement().localBasis());

        // rješenje
        double  u=0.0;
        for (size_t i=0; i<lfsu.size(); i++)
          u += x(lfsu,i)*phihat[i];

        const auto& gradphihat = cache.evaluateJacobian(xi, lfsu.finiteElement().localBasis());
        const auto jac = geo.jacobianInverseTransposed(xi);
        std::vector<Gradient> gradphi(lfsu.size());

        for (size_t i=0; i<lfsu.size(); i++)
          jac.mv(gradphihat[i][0],gradphi[i]);

        // grad u
        Dune::FieldVector<double,dim> gradu(0.0);
        for (size_t i=0; i<lfsu.size(); i++)
          gradu.axpy(x(lfsu,i),gradphi[i]);

        const auto factor = ip.weight()*geo.integrationElement(xi);
        const auto a = reactCoeff();

        for (size_t i=0; i<lfsv.size(); i++)
          r.accumulate(lfsv,i,( gradu*gradphi[i] + a*u*phihat[i] ) * factor);
      }
  }
};
