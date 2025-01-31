#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include "bivariates.hpp"



ZZ_pEXY::ZZ_pEXY(NTL::Vec<NTL::ZZ_pEX> cs){
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


ZZ_pEXY operator+(ZZ_pEXY const f, ZZ_pEXY const &other) { 
    
    
    auto dY = std::max(f.dY, other.dY);

    NTL::Vec<NTL::ZZ_pEX> cs_f;
    cs_f.SetLength(dY+1);

    for(unsigned i = 0; i <= f.dY; i++){
        cs_f[i] = f.coeffs[i];
    }
    for(unsigned i = f.dY+1; i <= dY; i++){
        cs_f[i] = NTL::ZZ_pEX(0);
    }

    NTL::Vec<NTL::ZZ_pEX> cs_other;
    cs_other.SetLength(dY+1);

    for(unsigned i = 0; i <= other.dY; i++){
        cs_other[i] = other.coeffs[i];
    }
    for(unsigned i = other.dY+1; i <= dY; i++){
        cs_other[i] = NTL::ZZ_pEX(0);
    }
    

    
    for(size_t i = 0; i <= dY; ++i) {
        cs_f[i] += cs_other[i];
    } 
    return ZZ_pEXY(cs_f); 
}

ZZ_pEXY &operator+=(ZZ_pEXY f, ZZ_pEXY const &other) { 
    return f = f + other; 
}

ZZ_pEXY operator-(ZZ_pEXY const f, ZZ_pEXY const &other) { 

    auto dY = std::max(f.dY, other.dY);

    NTL::Vec<NTL::ZZ_pEX> cs_f;
    cs_f.SetLength(dY+1);

    for(unsigned i = 0; i <= f.dY; i++){
        cs_f[i] = f.coeffs[i];
    }
    for(unsigned i = f.dY+1; i <= dY; i++){
        cs_f[i] = NTL::ZZ_pEX(0);
    }

    NTL::Vec<NTL::ZZ_pEX> cs_other;
    cs_other.SetLength(dY+1);

    for(unsigned i = 0; i <= other.dY; i++){
        cs_other[i] = other.coeffs[i];
    }
    for(unsigned i = other.dY+1; i <= dY; i++){
        cs_other[i] = NTL::ZZ_pEX(0);
    }
    
    for(size_t i = 0; i <= dY; ++i) {
        cs_f[i] -= cs_other[i];
    }

    return ZZ_pEXY(cs_f); 
}

ZZ_pEXY &operator-=(ZZ_pEXY f, ZZ_pEXY const &other) { 
    return f = f - other; 
}

ZZ_pEXY normalise(ZZ_pEXY F){
    size_t l = F.coeffs.length()-1;

    while(IsZero(F.coeffs[l])){
        l -= 1;
    }
    NTL::Vec<NTL::ZZ_pEX> cs;
    cs.SetLength(l+1);
    for(unsigned i = 0; i <= l; i++){
        cs[i] = F.coeffs[i];
    }
    return ZZ_pEXY(cs);
}

NTL::ZZ_pEX _reciprocal_kronecker_substitution(ZZ_pEXY const F, unsigned Lx, unsigned Ly, unsigned N){

    // This can be changed to make multiplication more efficient 
    // https://arxiv.org/abs/0712.4046
    

    NTL::ZZ_pEX f;
    f.SetLength((Ly-1)*N+Lx);

    ZZ_pEXY F1 = F;
    while(F1.coeffs.length() < Ly){
        (F1.coeffs).append(NTL::ZZ_pEX(0));
    }
    for(unsigned i = 0; i < Ly; i++){
        NTL::ZZ_pEX Fi = F1.coeffs[i];
        for(unsigned j = 0; j < Lx; j++){
            NTL::SetCoeff(f, i*N+j, coeff(Fi, j));
        };
    }
    return f;
}

ZZ_pEXY operator*(ZZ_pEXY const F, ZZ_pEXY const &other) { 

    // Write here assumptions on degree
    auto dX = std::max(F.dX, other.dX);
    auto dY = std::max(F.dY, other.dY);
    unsigned Lx = dX + 1;
    unsigned Ly = dY + 1;
    unsigned N = 2*Lx-1;
    

    NTL::ZZ_pEX f = _reciprocal_kronecker_substitution(F, Lx, Ly, N);    
    NTL::ZZ_pEX g = _reciprocal_kronecker_substitution(other, Lx, Ly, N);

    NTL::ZZ_pEX fg; 
    NTL::mul(fg, f, g);

    unsigned degFG = dY*dY;
    NTL::Vec<NTL::ZZ_pEX> coeffs;
    coeffs.SetLength(degFG+1);
    for(unsigned i = 0; i < degFG+1; i++){
        NTL::ZZ_pEX c;

        for(unsigned j = 0; j < N; j++){
            NTL::SetCoeff(c, j, coeff(fg, j+i*N));
            
        }
        c.normalize();
        coeffs[i] = c;
    }
    ZZ_pEXY G = ZZ_pEXY(coeffs);
    return normalise(G);
}

bool _divides(ZZ_pEXY const F, ZZ_pEXY const &other) {
    return ((F.dY >= other.dY) && (deg(F.coeffs[F.dY]) >= deg(other.coeffs[other.dY])));
}

bool _divides(monomial_ZZ_pE const M, monomial_ZZ_pE const &other) {
    return ((M.dY >= other.dY) && (M.dX >= other.dX));
}

monomial_ZZ_pE LeadingTerm(ZZ_pEXY F){
    monomial_ZZ_pE M;
    M.dY = F.dY;
    NTL::ZZ_pEX f = F.coeffs[F.dY];
    M.dX = deg(f);
    M.c = NTL::LeadCoeff(f);
    return M;
}

monomial_ZZ_pE operator*(monomial_ZZ_pE const M, monomial_ZZ_pE const &other) {
    monomial_ZZ_pE m;
    m.dY = M.dY + other.dY;
    m.dX = M.dX + other.dX;
    m.c = M.c * other.c;
    return m;
}

monomial_ZZ_pE operator/(monomial_ZZ_pE const M, monomial_ZZ_pE const &other) {
    assert(_divides(M, other));
    monomial_ZZ_pE m;
    m.dY = M.dY - other.dY;
    m.dX = M.dX - other.dX;
    m.c = M.c/other.c;
    return m;
}

ZZ_pEXY operator*(monomial_ZZ_pE const M, ZZ_pEXY const &other){
    NTL::ZZ_pEX f;
    NTL::SetCoeff(f, M.dX, M.c);
    
    NTL::Vec<NTL::ZZ_pEX> coeffs;
    coeffs.SetLength(M.dY+1);
    for(unsigned i = 0; i < M.dY; i++){
        coeffs[i] = NTL::ZZ_pEX(0);
    }
    coeffs[M.dY] = f;
    ZZ_pEXY F = ZZ_pEXY(coeffs);
    return F*other;
}

ZZ_pEXY operator+(monomial_ZZ_pE const M, ZZ_pEXY const &other){
    NTL::ZZ_pEX f;
    NTL::SetCoeff(f, M.dX, M.c);
    NTL::Vec<NTL::ZZ_pEX> coeffs;
    coeffs.SetLength(M.dY+1);
    coeffs[M.dY] = f;
    ZZ_pEXY F = ZZ_pEXY(coeffs);
    return other+F;
}

ZZ_pEXY operator+(ZZ_pEXY const &other, monomial_ZZ_pE const M){
    return M + other;
}

ZZ_pEXY operator-(ZZ_pEXY const &other, monomial_ZZ_pE const M){
    NTL::ZZ_pEX f;
    NTL::SetCoeff(f, M.dX, M.c);
    NTL::Vec<NTL::ZZ_pEX> coeffs;
    coeffs.SetLength(M.dY+1);
    coeffs[M.dY] = f;
    ZZ_pEXY F = ZZ_pEXY(coeffs);
    return other-F;
}




bool IsZero(ZZ_pEXY const F){
    for(unsigned i = 0; i < F.coeffs.length(); i++){
        if(!IsZero(F.coeffs[i])){
            return false;
        };
    }
    return true;
}

ZZ_pEXY operator%(ZZ_pEXY const F, ZZ_pEXY const other) { 
    // Source: Cox, Little, O'Shea 

    assert(_divides(F,other));

    ZZ_pEXY p = F;

    NTL::Vec<NTL::ZZ_pEX> cs;
    cs.SetLength(1);
    cs[0] = NTL::ZZ_pEX(0);
    ZZ_pEXY r = ZZ_pEXY(cs);
    ZZ_pEXY q = ZZ_pEXY(cs);
    
    bool is_zero = IsZero(p);

    while(!is_zero){
        p = normalise(p);
        monomial_ZZ_pE LTp = LeadingTerm(p);
        monomial_ZZ_pE LTo = LeadingTerm(other);
        if(_divides(LTp, LTo)){
            monomial_ZZ_pE m = LTp/LTo;
            q = q + m;
            p = p - m*other;
        }else{
            r = r + LTp;
            p = p - LTp;
        }
        is_zero = IsZero(p);
    }
    assert(IsZero(r));

    return q;
}

