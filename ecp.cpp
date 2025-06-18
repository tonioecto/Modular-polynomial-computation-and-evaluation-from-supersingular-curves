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

ecp operator*(NTL::ZZ k, ecp P) {

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
            result += Q;
        }
        Q += Q;
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
    if (this->curve() != other.curve())
        throw std::logic_error("points not on the same curve");

    if (!other) return *this;
    if (!*this) return other;


    if (&this->field() != &other.field())
        throw std::logic_error("points not over the same field");

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


