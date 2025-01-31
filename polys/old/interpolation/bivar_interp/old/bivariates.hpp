#include <NTL/ZZ_pX.h>
#include <iostream>
#include <vector>

struct bivar_poly{
    NTL::Vec<NTL::ZZ_pX> cs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an

};

int degree(const bivar_poly F);
NTL::ZZ_pX reciprocal_kronecker_sub(const bivar_poly F, const int Ly, const int Lx, const int N);
void bivariate_mult(bivar_poly &FG, const bivar_poly F, const bivar_poly G, const int Ly, const int Lx);
