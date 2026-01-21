#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"
#include "fast_ff.hpp"

/******   Implements Montgomery curves   ********/
// this is not generic at all, and is implemented without any fail safe for maximum performance, to use with CAUTION

class montXZ;

class mont {
    private: 
        FpE_elem _A;
        FpE_elem _C;
        FpE_elem _A24;
        FpE_elem _C24;
        FpE_elem _alpha;
        FpE_elem _beta;
        FpE_elem _sq; 

    public: 
        mont() {};

        mont(FpE_elem const &A) : _A{A}, _C{FpE_elem(1)} {this->_A24 = (A + 2); this->_C24 = FpE_elem(4);};
        mont(FpE_elem const &alpha, FpE_elem const &beta, FpE_elem const &x2);

        FpE_elem const &A() const {return this->_A;};
        FpE_elem const &C() const {return this->_C;};
        FpE_elem const &A24() const {return this->_A24;};
        FpE_elem const &C24() const {return this->_C24;};
        FpE_elem const &sq() const {return this->_sq;};
        FpE_elem const &beta() const {return this->_beta;};

        FpE_elem j_inv() const;
        montXZ to_montXZ(const ecp& P);
        
        void set_A(const FpE_elem &A) {this->_A = A;};
        void set_C(const FpE_elem &C) {this->_C = C;};
        void set_A24(const FpE_elem &A24) {this->_A24 = A24;};
        void set_C24(const FpE_elem &C24) {this->_C24 = C24;};
        void set_sq(const FpE_elem &sq) {this->_sq = sq;};
        void set_beta(const FpE_elem &beta) {this->_beta = beta;};

        void coerce_curve(Fp2k const &Fext);

        ec to_E();
        
};

mont mont_from_j(const FpE_elem &j);

mont mont_from_E(const ec &E);

class montXZ {
    private: 
        FpE_elem _x,_z;
        mont *curve;
    public:
        montXZ() {};
        montXZ(mont *E) {
            this->curve = E;
        };
        montXZ(mont *E, const FpE_elem &x, const FpE_elem &z) : _x{x}, _z{z} {this->curve = E;};
        
        FpE_elem x() const {return this->_x;};
        FpE_elem z() const {return this->_z;};
        void set_x(const FpE_elem &x) {this->_x = x;};
        void set_z(const FpE_elem &z) {this->_z = z;};

        montXZ(mont *E, const montXZ &P) : _x{P.x()}, _z{P.z()} {this->curve = E;};

        void normalize();

        void xDBL(const FpE_elem &A24, const FpE_elem &C24);
        void rec_XDBL(const int k, const FpE_elem &A24, const FpE_elem &C24); 
        bool is_zero();

};

// fast 2^k isogeny from kernel
std::vector<mont> fast_2k_isog(const mont &mE, const montXZ &P, const int k);