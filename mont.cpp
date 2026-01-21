#include "mont.hpp"
#include "Fp2k.hpp"
#include <NTL/lzz_pEXFactoring.h>
#include <vector>
#include "fast_ff.hpp"

mont mont_from_j(const FpE_elem &j) {
    
    // we find the roots of some polynomial
    FpEX_elem pj;
    FpE_elem jj = j / 256;
    pj.SetLength(4);
    NTL::SetCoeff(pj, 3);
    NTL::SetCoeff(pj, 2, - 9);
    NTL::SetCoeff(pj, 1, 27 - jj);
    NTL::SetCoeff(pj, 0, - 27 + 4 * jj);


    auto v = FindRoots(pj);
    for (int i = 0; i < v.length(); i++) {
        auto r = v[i];
        auto A = sqrt(r);
        if (A) {
            return mont(*A);
        }
    }

    return mont(FpE_elem(0));
}

mont mont_from_E(const ec &E) {
    FpEX_elem pj;
    pj.SetLength(4);
    NTL::SetCoeff(pj, 3);
    NTL::SetCoeff(pj, 2, 0);
    NTL::SetCoeff(pj, 1, E.a());
    NTL::SetCoeff(pj, 0, E.b());

    auto v = FindRoots(pj);
    for (int i = 0 ; i < v.length(); i++) {
        // auto beta = sqrt(3 * v[i] * v[i] + E.a());
        // if (beta) {
        //     assert(sqrt(*beta));
        //     assert(i < 2);
        //     return mont(v[i], *beta, v[i + 1]);
        // }
        Fp2 beta;
         if(fast_sqrt(beta, 3 * v[i] * v[i] + E.a())) {
            if (beta != 0) {
                if (i < 2) {
                    if (v[i + 1] == 0 && i + 2 < v.length()) {
                        return mont(v[i], beta, v[i + 2]);    
                    }
                    else if (v[i + 1] == 0 && i >=1) {
                        return mont(v[i], beta, v[i - 1]); 
                    }
                    else {
                        assert(v[i + 1] != 0);
                        return mont(v[i], beta, v[i + 1]);    
                    }
                }
                else {
                    assert(i == 2);
                    if (v[i - 1] == 0) {
                        assert(v[i - 2] != 0);
                        return mont(v[i], beta, v[i - 2]);    
                    }
                    else {
                        assert(v[i - 1] != 0);
                        return mont(v[i], beta, v[i - 1]);    
                    }

                }

                
            }
            
         }

    }
    assert(0);
    return (mont(FpE_elem(0)));
}

mont::mont(FpE_elem const &alpha, FpE_elem const &beta, FpE_elem const &x2) : _C{beta}, _alpha{alpha}, _beta{beta} {

    this->_A  = 3 * alpha;
    assert(beta != 0);
    this->_A24 = (3 * alpha + 2 * beta); this->_C24 = 4 * beta;  

    assert(x2 != 0);

    // this computes sqrt(A^2 - 4) = 2 x2 + A / C
    this->_sq = (2 * x2 + alpha) / beta;
    assert(this->_sq * this->_sq == (NTL::power(this->_A / this->_C, 2) - 4) ); 
}

montXZ mont::to_montXZ(const ecp &P) {
    return montXZ(this, P.get_x() - lift(this->_alpha, P.field()) * P.get_z(), lift(this->_beta, P.field()) * P.get_z() );
}

FpE_elem mont::j_inv() const {
    auto a2 = NTL::power(this->A() / this->C(), 2);
    return 256 * NTL::power((a2 - 3), 3) / (a2 - 4);
}

void montXZ::xDBL(const FpE_elem &A24, const FpE_elem &C24) {

    // assumes the correct field is pushed
    auto v1 = this->_x + this->_z;
    // v1 = v1 * v1;
    SpecialMul(v1, v1, v1);
    FpE_elem v2,v3;
    v2 = this->_x - this->_z;
    // v2 = v2 * v2;
    SpecialMul(v2, v2, v2);
    this->_x = v1; // x coordinate
    v1 = v1 - v2;
    // v2 = C24 * v2;
    SpecialMul(v2, C24, v2);
    // this->_x = this->_x * v2;
    SpecialMul(this->_x, this->_x, v2);
    // v3 = (A24 * v1);
    SpecialMul(v3, A24, v1);
    v3 = v3 + v2; 
    // this->_z = v1 * v3;
    SpecialMul(this->_z, v1, v3);
}

void montXZ::rec_XDBL(const int k, const FpE_elem &A24, const FpE_elem &C24) {
    for (int i = 0; i < k; i++) {
        this->xDBL(A24, C24);
    }
}

bool montXZ::is_zero() {
    return IsZero(this->_z);
}

void montXZ::normalize() {
    this->_x = this->_x / this->_z;
    this->_z = 1;
}

void two_isog_eval_and_dom(montXZ &P, mont &domain, const montXZ &kernel) {
    assert(kernel.z() != 0 );
    auto T1 = kernel.x();

    if ( T1 == 0 ) {
        auto A = ( (domain.A24() + domain.A24() + domain.A24() + domain.A24()) / domain.C24() - 2 );
        auto sq = &(domain.sq());
        
        // newA = - 2A / sq
        // newA24 =  (-A + sq) / (2 * sq)
        domain.set_A24( - A + *sq );
        domain.set_C24( (*sq) + (*sq) );

        // newx = ((x - 1)^2 + x (A + 2)) / (x * sq)
        
        // auto T0 = P.x() * P.z();
        FpE_elem T0, T2;
        SpecialMul(T0, P.x(), P.z());

        T1 = P.x() - P.z();
        SpecialMul(T2, T0, *sq); 
        P.set_z(T2);
        SpecialMul(T1, T1, T1);
        SpecialMul(T2, T0, A+2);
        P.set_x(T1 + T2);

    }
    else {
#ifndef NDEBUG 
    auto TT1 = kernel.x() / kernel.z();
    assert(  ( TT1 * TT1 + ( (domain.A24() + domain.A24() + domain.A24() + domain.A24()) / domain.C24() - 2 ) * TT1  + 1) == 0);
#endif 
        auto T0 = kernel.z() + T1;
        T1 = kernel.z() - T1;
        FpE_elem temp;
        SpecialMul(temp, T0, T1);
        domain.set_A24(temp);
        SpecialMul(temp, kernel.z(), kernel.z());
        domain.set_C24(temp);


        auto T2 = P.x() + P.z();
        auto T3 = P.x() - P.z();
        SpecialMul(T0, T3, T0);
        SpecialMul(T1, T2, T1);
        T2 = T0 - T1;
        T3 = T0 + T1;
        SpecialMul(T0, P.x(), T2);P.set_x(T0);
        SpecialMul(T1, P.z(), T3);P.set_z(T1); 
    }

}


void mont::coerce_curve(Fp2k const &Fext) {
    this->set_A(coerce(this->A(), Fext));
    this->set_C(coerce(this->C(), Fext));
    this->set_A24(coerce(this->A24(), Fext));
    this->set_C24(coerce(this->C24(), Fext));
}

std::vector<mont> fast_2k_isog(const mont &mE, const montXZ &P, const int k) {
    std::vector<mont> codoms(4);
    
    auto tP1 = montXZ(P);
    
    mont domain = mE;

    // assumes all lifts have been performed
    for (int i = 0; i < k ; i++ ) {
        auto tP2 = montXZ(&domain, tP1);
        tP2.rec_XDBL(k - i - 1, domain.A24(), domain.C24());
#ifndef NDEBUG 
        auto testP = tP2;
        assert(!testP.is_zero());
        testP.xDBL(domain.A24(), domain.C24());
        assert(testP.is_zero());
#endif 
        two_isog_eval_and_dom(tP1, domain, tP2);
        codoms[i] = domain; 
        // if (i == 0) {
        //     codoms[i].set_A(4 * codoms[i].A24() - 2 * codoms[i].C24());
        //     codoms[i].set_C(codoms[i].C24());
        //     codoms[i].set_beta( mE.beta() * tP2.x() / tP2.z() );
        // }
    }
    codoms[k-1].set_A(codoms[k-1].A24() + codoms[k-1].A24() + codoms[k-1].A24() + codoms[k-1].A24() - codoms[k-1].C24() - codoms[k-1].C24());
    codoms[k-1].set_C(codoms[k-1].C24());


    return codoms;

}