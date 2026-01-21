//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// The code in this file develops the class ecp, which represents an elliptic curve point,
/////// and related functions for elliptic curve arithmetic
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cassert>

#include <optional>

#include <NTL/ZZ_pE.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"
#include "fast_ff.hpp"


ecp::ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, FpE_elem const &xx, FpE_elem const &yy)
    : E{curve}, Fext{Fext}, x{xx}, y{yy}, z{FpE_elem(1)}
{
#ifndef NDEBUG
    FpE_push Push(Fext.F);
    ec E = *curve.get();
    assert(y*y == x*x*x + lift(E.a(), Fext)*x + lift(E.b(), Fext));
#endif
}

ecp::ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, FpE_elem const &xx, FpE_elem const &yy, FpE_elem const &zz)
    : E{curve}, Fext{Fext}, x{xx}, y{yy}, z{zz}
{
#ifndef NDEBUG
    FpE_push Push(Fext.F);
    assert(NTL::deg(FpE_elem::modulus()) == 2*Fext.k);
    ec E = *curve.get();
    assert(y*y*z == x*x*x + lift(E.a(), Fext)*x*z*z + lift(E.b(), Fext)*z*z*z);
#endif
}

void fast_dbl(ecp &Q, const ecp &P) {

    if (!P) {Q = P; return;}

    // assumes the correct field has already been pushed
    FpE_elem XX,ZZ, s, ss, B;
    // XX = P.get_x() * P.get_x();
    SpecialMul(XX, P.get_x(), P.get_x());
    // ZZ = P.get_z() * P.get_z();
    SpecialMul(ZZ, P.get_z(), P.get_z());
    // ZZ = lift(P.curve().a(), P.field()) * ZZ;
    SpecialMul(ZZ, ZZ, lift(P.curve().a(), P.field()));
    ZZ = ZZ + XX + XX + XX;
    // s = P.get_y() * P.get_z();
    SpecialMul(s, P.get_y(), P.get_z());
    s = s + s;
    // ss = s * s;
    SpecialMul(ss, s, s);
    // ss = ss * s;
    SpecialMul(ss, ss, s);
    Q.set_z(ss);

    // ss = P.get_y() * s;
    SpecialMul(ss, P.get_y(), s);
    B = P.get_x() + ss; 
    // B = B * B;
    SpecialMul(B, B, B);
    // ss = ss * ss;
    SpecialMul(ss, ss, ss);
    B = B - ss - XX;
    // XX = ZZ * ZZ - B - B;
    SpecialMul(XX, ZZ, ZZ);
    XX = XX - B - B;

    // Q.set_x(XX * s);
    // Q.set_y(ZZ * (B - XX) - ss - ss);
    SpecialMul(ZZ, ZZ, (B - XX));
    ZZ = ZZ - ss - ss;
    Q.set_y(ZZ);
    SpecialMul(XX, XX, s);
    Q.set_x(XX);

}

void fast_add(ecp &R, const ecp &P1, const ecp &P2, const FpE_elem &a, const FpE_elem &b3) {
    
    if (!P1) {R = P2; return;}
    if (!P2) {R = P1; return ;}

    FpE_elem t0,t1,t2,t3,t4,t5,X3,Y3,Z3;

    // first testing if P = +-Q
    
    // t0 = P1.get_x() * P2.get_x();
    // t1 = P1.get_y() * P2.get_y();
    // t2 = P1.get_z() * P2.get_z();
    SpecialMul(t0, P1.get_x(), P2.get_x());
    SpecialMul(t1, P1.get_y(), P2.get_y());
    SpecialMul(t2, P1.get_z(), P2.get_z());
    t3 = P1.get_x() + P1.get_y();
    t4 = P2.get_x() + P2.get_y();
    // t3 = t3 * t4;
    SpecialMul(t3, t3, t4);
    t4 = t0 + t1;
    t3 = t3 - t4;
    t4 = P1.get_x() + P1.get_z();
    t5 = P2.get_x() + P2.get_z();
    // t4 = t4 * t5;
    SpecialMul(t4, t4, t5);
    t5 = t0 + t2;
    t4 = t4 - t5;
    t5 = P1.get_y() + P1.get_z();
    X3 = P2.get_y() + P2.get_z();
    // t5 = t5 * X3;
    SpecialMul(t5, t5, X3);
    X3 = t1 + t2;
    t5 = t5 - X3;
    SpecialMul(Z3, a, t4);
    SpecialMul(X3, b3, t2);
    Z3 = X3 + Z3;
    X3 = t1 - Z3;
    Z3 = t1 + Z3;
    // Y3 = X3 * Z3;
    SpecialMul(Y3, X3, Z3);
    t1 = t0 + t0;
    t1 = t1 + t0;
    SpecialMul(t2, a, t2);
    SpecialMul(t4, b3, t4);
    t1 = t1 + t2;
    t2 = t0 - t2;
    SpecialMul(t2, a, t2);
    t4 = t4 + t2;
    SpecialMul(t0, t1, t4);
    Y3 = Y3 + t0;
    SpecialMul(t0, t5, t4);
    SpecialMul(X3, t3, X3);
    X3 = X3 - t0;
    SpecialMul(t0, t3, t1);
    SpecialMul(Z3, t5, Z3);
    Z3 = Z3 + t0;

    R.set_x(X3);
    R.set_y(Y3);
    R.set_z(Z3);
}

ecp operator*(NTL::ZZ k, ecp P) {

    FpE_push push(P.field().F);

    FpE_elem a = lift(P.curve().a(), P.field());
    FpE_elem b3 = lift(P.curve().b(), P.field());
    b3 = b3 + b3 + b3;

    if (k==0) {
        return ecp(P.curve_ptr(),P.field());
    }

    if (k < 0) {
        return (-k)*(-P);
    }

    ecp result(P.E, P.field());
    ecp Q = P;

    while (k > 0) {
        if (k % 2 == 1) { 
            fast_add(result, Q, result, a, b3);
        }
        fast_dbl(Q, Q);

        k >>= 1;
    }

    return result;
}


void ecp::normalize() const
{
    FpE_push push(field().F);
    if (NTL::IsZero(z)) {
        NTL::set(y);
    }
    else {
        auto invz = NTL::inv(z);
        x *= invz;
        y *= invz;
        NTL::set(z);
    }
}

std::pair<FpE_elem, FpE_elem> ecp::_lambda(ecp const &other) const
{
    //Notice this is (should?) only be called from +, where ZZ_pEPush is already done, so we dont have to do it again
    auto const &P = *this, &Q = other;
    FpE_elem const &x1 = P.x, &y1 = P.y, &z1 = P.z;
    FpE_elem const &x2 = Q.x, &y2 = Q.y, &z2 = Q.z;

    if (x1*z2 == x2*z1) {
        FpE_elem x1sq = x1 * x1;
        FpE_elem z1sq = z1 * z1;
        return {x1sq+x1sq+x1sq + lift(this->curve().a(), P.field())*z1sq, (y1+y1) * z1};
    }

    return {y2*z1 - y1*z2, x2*z1 - x1*z2};
}


ecp ecp::operator+(ecp other) const
{
#ifndef NDEBUG
    if (this->curve() != other.curve())
        throw std::logic_error("points not on the same curve");
#endif 

    if (!other) return *this;
    if (!*this) return other;

#ifndef NDEBUG
    if (&this->field() != &other.field())
        throw std::logic_error("points not over the same field");
#endif

    FpE_push push(this->field().F);

    FpE_elem const &x1 = this->x, &y1 = this->y, &z1 = this->z;
    FpE_elem const &x2 = other.x, &y2 = other.y, &z2 = other.z;

    if (x1*z2 == x2*z1 && y1*z2 != y2*z1) {
        assert(y2*z1 == -y1*z2);
        return {this->E, this->Fext};
    }

    auto lam = this->_lambda(other);

    FpE_elem const &u = lam.first, &v = lam.second;

    FpE_elem uu = u*u, vv = v*v;
    FpE_elem zz = z1 * z2;

    FpE_elem z3_ = vv*zz;
    FpE_elem x3_ = uu*zz - vv*(x1*z2 + x2*z1);

    FpE_elem vz1 = v*z1;
    FpE_elem y3 = u*(x1*z3_ - x3_*z1) - y1*v*z3_;
    FpE_elem x3 = x3_*vz1;
    FpE_elem z3 = z3_*vz1;

    return ecp(this->E, this->Fext, x3,y3,z3);
}


