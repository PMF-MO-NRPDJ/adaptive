#pragma once

#include <cmath>

// Nelinearnost
double reactCoeff()
{
    double eta = 0.1;
    return eta;
}

// Desna strana
template<typename V>
double RHS(const V & x)
{
  return 0.0;
}


// Dirichletov rubni uvjet
template<typename V>
double exact(const V & x)
{
     double theta = std::atan2(x[1],x[0]);
     if(theta < 0.0) theta += 2*M_PI;
     auto r = x.two_norm();
     return pow(r,2.0/3.0)*std::sin(theta*2.0/3.0);
}



// Neumannov rubni uvjet
template<typename V>
double neumannBC(const V& x)
{
  return 0.0;
}


