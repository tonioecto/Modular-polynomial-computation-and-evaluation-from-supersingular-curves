///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing the necessary multivariate polynomial arithmetic for Weber variants
////     of ModEvalBigLevel and ModEvalBigChar: see getweber.cpp
////    (Note that such arithmetic is not available in NTL)
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "multivariates.hpp"

////////////////////////////
// Bivariate Polynomials
////////////////////////////

ZZ_pEXY::ZZ_pEXY(NTL::Vec<FpEX_elem> cs){
    this->coeffs = cs;

    unsigned l = cs.length();
    this->dY = l-1;

    unsigned dX = 0;
    for (size_t i = 0; i < l; i++){
        auto D = deg(cs[i]);
        if(D > dX){
            dX = D;
        }
    }
    this->dX = dX;
}

std::vector<FpE_elem> _powers(FpE_elem a, unsigned k){

    std::vector<FpE_elem> powers(k+1);

    FpE_elem A = FpE_elem(1);
    powers[0] = A;

    for(unsigned i = 1; i <= k; i++){
        A = A*a;
        powers[i] = A;
    }

    return powers;
}

FpEX_elem EvaluateBivariate(ZZ_pEXY F, FpE_elem a){

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Evaluate bivariate in variables X,Y at Y = a to get a polynomial in X.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto powers = _powers(a, F.dY+1);

    FpEX_elem f = FpEX_elem(0);

    for(unsigned i = 0; i <= F.dY; i++){
        FpEX_elem g = powers[i]*F.coeffs[i];
        f = f + g;
    }

    return f;
}

////////////////////////////////////////////////////
// Multivariate Polynomials with three variables
////////////////////////////////////////////////////

ZZ_pEXYZ::ZZ_pEXYZ(std::vector<ZZ_pEXY> cs){
    this->coeffs = cs;

    unsigned l = cs.size();
    this->dZ = l-1;
}

FpEX_elem EvaluateTrivariate(ZZ_pEXYZ F, FpE_elem y, FpE_elem z){
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Evaluate multivariate in variables X,Y,Z at Y = y and Z = z to get a polynomial in X.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    auto powers = _powers(z, F.dZ+1);

    FpEX_elem f = FpEX_elem(0);

    for(unsigned i = 0; i <= F.dZ; i++){
        ZZ_pEXY C = F.coeffs[i];
        // TODO: Could precompute the powers of y that we need for extra efficiency here
        auto c = EvaluateBivariate(C, y);
        f = f + powers[i]*c;
    }

    return f;
}
