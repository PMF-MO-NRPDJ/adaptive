#pragma once
#include <cmath>

double reactCoeff()
{
    return 100;
}



// Harmonijska funkcija r^beta sin(beta*theta), beta = 2/3.
template<typename V>
double exact(const V & x)
{
     double theta = std::atan2(x[1],x[0]);  // u (-pi, pi)
     if(theta < 0.0) theta += 2*M_PI;
     auto r = x.two_norm2();
     return pow(r,1.0/3.0)*std::sin(theta*2.0/3.0);
}

// Neumannov rubni uvjet
template<typename V>
double neumannBC(const V& x)
{
  return 0.0;
}

// Desna strana
template<typename V>
double RHS(const V & x)
{
  return reactCoeff()*exact(x);
}
