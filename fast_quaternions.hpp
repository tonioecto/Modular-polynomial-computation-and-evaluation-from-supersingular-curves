
#pragma once
#include <optional>
#include <cassert>
#include <functional>
#include "quaternions.hpp"
#include <time.h>

typedef long FastInteger;

typedef long double FastFloat;

struct SignedBarrettReducer {
    FastInteger b;
    FastInteger bsqr;
    FastInteger m;
    FastInteger msqr;
    FastInteger m2;
    FastInteger m2sqr;


    SignedBarrettReducer(FastInteger b) : b(b) {
        assert(b > 0);
        bsqr = b * b;
        m = (((FastInteger) 2)<< 60) / b + 1;
        msqr = (((FastInteger) 2)<< 60) / (bsqr) + 1;
        m2 = (((FastInteger) 2)<< 60) / (b << 1) + 1;
        m2sqr = (((FastInteger) 2)<< 60) / (bsqr << 1) + 1;
    };

    // mod b 
    FastInteger mod(FastInteger a) const;
    // mod bsqr = b²
    FastInteger modsqr(FastInteger a) const;
    // mod 2*b
    FastInteger mod2(FastInteger a) const;
    // mod 2*b²
    FastInteger mod2sqr(FastInteger a) const;

    // multiplication of two elements mod bsqr
    FastInteger mulmodsqr(FastInteger a, FastInteger c) const;

};


FastInteger convert(const Integer &x);

FastInteger GCD(const FastInteger &a1, const FastInteger &a2); 
FastInteger InvMod(const FastInteger &a, const SignedBarrettReducer &red);
FastInteger InvModSqr(const FastInteger &a, const SignedBarrettReducer &red);
FastInteger InvMod2(const FastInteger &a, const SignedBarrettReducer &red);
FastInteger InvMod2Sqr(const FastInteger &a, const SignedBarrettReducer &red);
FastFloat IntToFloat(FastInteger a1, FastInteger a2);
FastInteger Rounding(FastFloat f);

FastInteger BytesToFastInt(unsigned char* buffer, int BUFFER_SIZE);
void FastIntToBytes(FastInteger value, unsigned char* buffer, int BUFFER_SIZE);

bool IsOdd(const FastInteger a);



class FastQuatLat;
class FastQuat;

class FastMat4;
class FastMat3;

struct FastQuatAlg
{
    FastInteger p, q;
    FastQuatAlg(const quatalg &A);
};


class FastMat3 {
    private:
        std::array<std::array<FastInteger,3>,3> mat;

    public: 
        std::array<FastInteger,3>& operator[](size_t i) {
            return mat[i];
        }

        friend std::ostream& operator<<(std::ostream& o, FastMat3 &M) {
            return o << "[ " << M[0][0] << ", " << M[0][1] << ", " << M[0][2] << ", \n" << "  " <<  M[1][0] << ", " << M[1][1] << ", " << M[1][2] << ", \n" << "  " <<  M[2][0] << ", " << M[2][1] << ", " << M[2][2] << ", \n" << " ] \n";
        }

        bool operator==(NTL::mat_ZZ &other);
        bool operator==(FastMat3 &other);

};

class FastMat4 {
    private:
        std::array<std::array<FastInteger,4>,4> mat;

    public: 
        std::array<FastInteger,4>& operator[](size_t i) {
            return mat[i];
        }

        std::array<FastInteger,4> const &operator[](size_t i) const {
            return mat[i];
        }

        friend std::ostream& operator<<(std::ostream& o, FastMat4 &M) {
            return o << "[ " << M[0][0] << ", " << M[0][1] << ", " << M[0][2] << ", " << M[0][3] << ", \n" << "  " <<  M[1][0] << ", " << M[1][1] << ", " << M[1][2] << ", " << M[1][3] << ", \n" << "  " <<  M[2][0] << ", " << M[2][1] << ", " << M[2][2] << ", " << M[2][3] << ", \n" << "  " <<  M[3][0] << ", " << M[3][1] << ", " << M[3][2] << ", " << M[3][3] << " ] \n";
        }

};

class FastQuat {

    private:
        std::array<FastInteger, 5> coeffs;

    public:
        FastQuatAlg const &alg;

        FastQuat(FastQuat const &alpha) : alg{alpha.alg} { *this = alpha; };
        FastQuat(FastQuatAlg const &alg_) : alg{alg_} { coeffs[0] = 0; };
        FastQuat(quat const &alpha, FastQuatAlg const &alg);

        FastQuat(std::array<FastInteger, 5> const &coeffs_, FastQuatAlg const &alg_)
        : coeffs{coeffs_}, alg{alg_}
        {
            normalize();
        }

        void normalize();
        void normalize2();
        
        FastQuat& operator=(const FastQuat& other) {
            if (this != &other) {  // Avoid self-assignment
                // Copy the value
                coeffs = other.coeffs;
            }
            return *this;
        }

        FastInteger &operator[](size_t index) {
            assert(index < 5);
            return coeffs[index];
        }

        FastInteger const &operator[](size_t index) const {
            assert(index < 5);
            return coeffs[index];
        }

        friend std::ostream& operator<<(std::ostream& o, FastQuat const &alpha) {
            return o
            << alpha[0] << "/" << alpha[4] << " + "
            << alpha[1] << "/" << alpha[4] << "*i + "
            << alpha[2] << "/" << alpha[4] << "*j + "
            << alpha[3] << "/" << alpha[4] << "*k";
    
        }

        FastQuat operator*(FastQuat const &other) const;

        FastQuat operator*(FastInteger n) const { return {{n*coeffs[0], n*coeffs[1], n*coeffs[2], n*coeffs[3], coeffs[4]}, alg}; }

        FastQuat operator+(FastQuat const &other) const;

        FastInteger integral_norm() const;

        void mul_left(FastQuat const &other);

        void scalar_mul(const FastInteger N);

        FastInteger scalar_remove();

        bool is_zero() const;
};

class FastQuatLat{
    public: 
        FastQuatAlg const &alg;
        FastMat4 basis; // beware when good_type is true, basis[0][1], basis[0][2] and basis[0][3] (which are supposed to be 0) will contain some precomputation data
        FastInteger denom;
        bool good_type; // tells us if it is an ideal with a basis in HNF form of a very specific shape,
        //  in that case, some of the zero coefficients may contain some precomputation used to compute faster the right order and its gram matrix
    private :
        mutable std::pair<FastInteger, FastInteger> the_norm {};
        // mutable FastQuat gen;


    public: 
        FastQuatLat(FastMat4 const &gens, std::pair<FastInteger,FastInteger> const &norm, FastInteger const &denom_, FastQuatAlg const &alg_, bool normalize, bool good_type_);

        FastQuatLat(FastQuatLat const &L) : alg{L.alg} { *this = L; }

        FastQuatLat(quatlat const L, FastQuatAlg const &_alg);

        void normalize(bool safe = false);

        FastQuatLat& operator=(const FastQuatLat& other) {
            if (this != &other) {  // Avoid self-assignment
                // Copy the values
                if (&this->alg != &other.alg) {
                    throw std::runtime_error("LOUD ERROR WHEN ASSIGNING QUATLAT");
                }
                basis = other.basis;
                denom = other.denom;
                the_norm = other.the_norm;
                good_type = other.good_type;
            }
            return *this;
        }

        friend std::ostream& operator<<(std::ostream& o, FastQuatLat &L) {
            return o << L.basis << "/" << L.denom << "  (norm " << L.norm().first << "/" << L.norm().second << ")";
        }
        bool operator==(FastQuatLat &other);

        std::pair<FastInteger, FastInteger> norm() const;

        // red contains the reducer for this.norm * other.norm
        // red1 contains the reduced for other.norm
        // inv is the inverse of this.norm mod 2 * other.norm
        // invsqr2 is the inverse of this.norm^2 mod 2 other.norm^2 
        void _fast_intersect(FastQuatLat &other, const SignedBarrettReducer &red, const SignedBarrettReducer &red1, const FastInteger &inv, const FastInteger &inv2sqr);
        bool fast_contains(FastQuat const &alpha, FastInteger const &modnorm) const;

        FastQuatLat new_right_order();
        std::pair<FastQuatLat, FastMat3> fast_right_order_and_gram(const SignedBarrettReducer &red);

        quatlat FastQuatLat_to_quatlat(const quatalg &Bp);

        FastQuatLat copy() const;


};

FastQuatLat create_from_generator_O0(const FastQuat &gen, const FastInteger &norm);
std::vector<FastQuatLat> left_ideals_of_prime_norm_O0(FastInteger const &ell, const quatalg &Bp, const FastQuatAlg &fast_Bp);

bool commutatorfind(const std::pair<FastQuat,FastQuat> &beta1, const std::pair<FastQuat,FastQuat> &beta2, const std::pair<FastInteger,FastInteger> &n, const FastQuat &prod_quat, const FastInteger &prod_quat_norm, bool is_FP, FastQuat &quat_sol, const SignedBarrettReducer &redp);