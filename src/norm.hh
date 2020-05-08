#pragma once
#include <dune/geometry/quadraturerules.hh>


template<class GV, class DGF>
double L2norm(GV const & gv, DGF const & dgf, std::size_t intorder = 4)
{
    double integral = 0.0;
    for(auto const & elem : elements(gv))
    {
        const auto & rule = Dune::QuadratureRules<double,GV::dimension>::rule( elem.geometry().type() , intorder );
        for(auto const & qpoint : rule)
        {
            auto const xi = qpoint.position();
            auto const dx = qpoint.weight() * elem.geometry().integrationElement(xi);

            //Evaluiram DiscreteGridFunction dgf u lokalnoj koordinati
            typename DGF::Traits::RangeType f;
            dgf.evaluate(elem, xi, f);

            integral += f.two_norm2() * dx;
        }
    }
    return std::sqrt(integral);
}


