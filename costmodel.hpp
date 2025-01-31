//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// The code in this file builds the cost model with which we choose the 'optimal' torsion (header)
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <functional>

using model_function = std::function<NTL::RR(int const&, int const&, int const&)>;

model_function cost_model(NTL::ZZ const &p);
NTL::RR karatsuba(int const &k, NTL::RR const &logp, NTL::RR const &loglogp);
NTL::RR adj_cost(NTL::RR const &cost, int const &ell);
