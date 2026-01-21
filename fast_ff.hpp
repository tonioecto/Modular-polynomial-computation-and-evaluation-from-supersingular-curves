#pragma once

#include "fast_quaternions.hpp"
#include "hashmap.hpp"
#include "multivariates.hpp"

typedef std::pair<Fp,Fp> ffp2; 


ffp2 from_Fp2(const Fp2 &a);
void to_Fp2 (Fp2 &t, const ffp2 &a);
Fp2 Fp2_cast(const ffp2 &a);
bool is_zero(const ffp2 &a);
bool is_one(const ffp2 &a);
inline bool is_Fp_fast(const ffp2 &a) {
    return NTL::IsZero(a.second);
};
void fast_mul(ffp2 &c, const ffp2 &a, const ffp2 &b);
void fast_mul_with_precomp(ffp2 &c, const ffp2 &a, const ffp2 &b, const std::vector<NTL::mulmod_precon_t> &prep);
void fast_add(ffp2 &c, const ffp2 &a, const ffp2 &b);
void negate(ffp2 &a, const ffp2 &b);
void fast_sqr(ffp2 &s, const ffp2 &a);
void fast_inv(ffp2 &y, const ffp2 &x);
ffp2 fast_Frob(const ffp2 &x);

void fast_pow(ffp2 &y, const ffp2 &x, long e);

inline std::ostream& operator<<(std::ostream& o, ffp2 const &P)
        {
            o << "[" << P.first << " " << P.second << "] ";
            return o;
        }

void get_powers(std::vector<ffp2> &powers, const ffp2 &a, int k);

class ffp2X {
    public:
        long deg;
        std::vector<ffp2> coeffs;

        ffp2X() {this->coeffs = {}; this->deg = - 1;};
        ffp2X(int d) : deg(d), coeffs(d+1) {};
        ffp2X(const std::vector<ffp2> &co, int degree);
        ffp2X(const Fp2X &poly);
        // coeff * X^degree
        ffp2X(int degree, long coeff);

        friend std::ostream& operator<<(std::ostream& o, ffp2X const &P)
        {        
            for (int i = 0; i <= P.deg; i++) {
                o << "[" << P.coeffs[i].first << " " << P.coeffs[i].second << "] ";
            }
            return o;
        }

        

        bool is_zero_poly() const;
        void resize(int k);
        void normalize();
        ffp2 lead_coeff();
        void mul(const ffp2 &x);
        bool operator==(const ffp2X &other);

        // in place reminder of division by b
        void fast_PlainRem(const ffp2X &b);

        void Evaluate(ffp2 &x, const ffp2 &a) const;
};

void sub(ffp2X &c, const ffp2X &a, const ffp2X &b);

void fast_PlainGCD_2(ffp2X& x, const ffp2X& a, const ffp2X& b);

class ffp2XY{
    public: 
    std::vector<ffp2X> coeffs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an
    int degX; // Degree of X in the bivariate polynomial
    int degY; // Degree of Y in the bivariate polynomial

    ffp2XY() {this->coeffs = {}; this->degX = - 1; this-> degY = - 1;};
    ffp2XY(const std::vector<ffp2X> &coeffs);
    ffp2XY(const ZZ_pEXY &F);

    void Evaluate(ffp2X &x, const ffp2 &a) const;
    void Evaluate(ffp2X &x, const std::vector<ffp2> &powers) const;
    void Evaluate(ffp2 &t, const ffp2 &x, const ffp2 &y) const;
};

class ffp2XYZ{
    public: 
    std::vector<ffp2XY> coeffs; // coefficients are univariate polynomials [a0,a1,a2,...,an] where a_i \in Fp[x] and f(x,y) = a0y^n + a1*y^(n-1) + ... + an
    int degZ; // Degree of Z in the trivariate polynomial

    ffp2XYZ() {};
    ffp2XYZ(const ZZ_pEXYZ &F);
};

// fast operations over Fp2
void fast_mul(Fp2 &c, const Fp2 &a, const Fp2 &b);
void faster_mul(Fp2 &c, const Fp2 &a, const Fp2 &b, FpX &t);
bool fast_sqrt(Fp2 &s, const Fp2 &a);
void fast_sqr(Fp2 &s, const Fp2 &a);
void fast_pow(Fp2 &y, const Fp2 &x, long e);
void fast_inv(Fp2 &y, const Fp2 &x);
void fast_add(Fp2 &c, const Fp2 &a, const Fp2 &b, FpX &temp);
void fast_set(Fp2 &c, const FpX &t);

// fast GCD for Fp2 polynomials
void fast_PlainGCD(FpEX_elem& x, const FpEX_elem& a, const FpEX_elem& b);

// for comparison
void fast_mul_pair(ffp2 &c, const ffp2 &a, const ffp2 &b);
void fast_add_pair(ffp2 &c, const ffp2 &a, const ffp2&b);


// FF_p^{2k}
class ffp2k {
    
    public:
    
        ffp2X def_poly; // the polynomial representing the ffpE
        ffp2X const *mod_poly; // the polynomial of degree k defining the finite field
        
        ffp2k(const ffp2X *_modpoly) : def_poly(), mod_poly(_modpoly) {};
        ffp2k(const ffp2X &_defpoly, const ffp2X *_modpoly) {this->def_poly = _defpoly; this->mod_poly = _modpoly;};

        void normalize();

};

ffp2k random_ffp2k(const ffp2X *mod_poly);
void mul(ffp2k &c, const ffp2k &b, ffp2k &a);


// faster multiplication for Fp2k elements when we know the modulus has a specific shape (lots of zero coefficients)
void SpecialMul(FpE_elem& x, const FpE_elem& a, const FpE_elem& b);


