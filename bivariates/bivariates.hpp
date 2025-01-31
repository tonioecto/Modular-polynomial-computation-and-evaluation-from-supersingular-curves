#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cassert>

struct ZZ_pEXY{
    NTL::Vec<NTL::ZZ_pEX> coeffs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an
    unsigned dX; // Degree of X in the bivariate polynomial
    unsigned dY; // Degree of Y in the bivariate polynomial

    ZZ_pEXY(NTL::Vec<NTL::ZZ_pEX> coeffs);
};

struct monomial_ZZ_pE{
    unsigned dX; // Degree of X in the monomial
    unsigned dY; // Degree of Y in the monomial
    NTL::ZZ_pE c; // Coefficient of monomial
};

ZZ_pEXY operator+(ZZ_pEXY const f, ZZ_pEXY const &other);
ZZ_pEXY &operator+=(ZZ_pEXY f, ZZ_pEXY const &other);
ZZ_pEXY operator-(ZZ_pEXY const f, ZZ_pEXY const &other);
ZZ_pEXY &operator-=(ZZ_pEXY f, ZZ_pEXY const &other);
NTL::ZZ_pEX _reciprocal_kronecker_substitution(ZZ_pEXY const F, unsigned Lx, unsigned Ly, unsigned N);
ZZ_pEXY operator*(ZZ_pEXY const F, ZZ_pEXY const &other);

bool _divides(ZZ_pEXY const F, ZZ_pEXY const &other);
bool _divides(monomial_ZZ_pE const M, monomial_ZZ_pE const &other);
monomial_ZZ_pE LeadingTerm(ZZ_pEXY F);
monomial_ZZ_pE operator*(monomial_ZZ_pE const M, monomial_ZZ_pE const &other);
monomial_ZZ_pE operator/(monomial_ZZ_pE const M, monomial_ZZ_pE const &other);
ZZ_pEXY operator*(monomial_ZZ_pE const M, ZZ_pEXY const &other);
ZZ_pEXY operator+(monomial_ZZ_pE const M, ZZ_pEXY const &other);
ZZ_pEXY operator-(ZZ_pEXY const &other, monomial_ZZ_pE const M);
bool IsZero(ZZ_pEXY const F);
ZZ_pEXY operator%(ZZ_pEXY const F, ZZ_pEXY const other);