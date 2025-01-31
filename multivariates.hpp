///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing the necessary multivariate polynomial arithmetic (header file) 
///////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cassert>

#include "Fp2k.hpp"

////////////////////////
// Bivariates
////////////////////////

struct ZZ_pEXY{
    NTL::Vec<FpEX_elem> coeffs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an
    unsigned dX; // Degree of X in the bivariate polynomial
    unsigned dY; // Degree of Y in the bivariate polynomial

    ZZ_pEXY(NTL::Vec<FpEX_elem> coeffs);
};


std::vector<FpE_elem> _powers(FpE_elem a, int k);
FpEX_elem EvaluateBivariate(ZZ_pEXY F, FpE_elem a);

////////////////////////////////////////////////////
// Multivariate polynomials with three variables
////////////////////////////////////////////////////

struct ZZ_pEXYZ{
    std::vector<ZZ_pEXY> coeffs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an
    unsigned dZ; // Degree of Z in the trivariate polynomial

    ZZ_pEXYZ(std::vector<ZZ_pEXY> coeffs);
};

FpEX_elem EvaluateTrivariate(ZZ_pEXYZ F, FpE_elem y, FpE_elem z);
