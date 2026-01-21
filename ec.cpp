///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////// The code in this file develops the class ec, which represents an elliptic curve,
/////// and related functions
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <optional>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"

ec::ec(FpE_elem const &a, FpE_elem const &b) : _a{a}, _b{b}
{
    if ((4*this->_a*this->_a*this->_a + 27*this->_b*this->_b) == 0)
        throw std::logic_error("curve is singular");
}

ecp ec::operator()(Fp2k const &Fext, FpE_elem new_x, FpE_elem new_y) const {
    return ecp(std::make_shared<const ec>(*this), Fext, new_x, new_y);
}
ecp ec::operator()(Fp2k const &Fext) const {
    return ecp(std::make_shared<const ec>(*this), Fext);
}

std::optional<ecp> ec::lift_x(Fp2k const &Fext, FpE_elem const &x) const
{
    FpE_push Push(Fext.F);

    FpE_elem rhs = (x*x + lift(this->_a, Fext))*x + lift(this->_b, Fext);
    auto y = fast_sqrt(Fext, rhs);
    if (y)
        return ecp(std::make_shared<const ec>(*this), Fext, x, *y);
    return {};
}

ecp ec::random_point(Fp2k const &Fext) const
{
    FpE_push Push(Fext.F);
    while (true) {
        auto pt = this->lift_x(Fext, random_FpE_elem());
        if (pt) {
            return *pt;
        }
    }
}

FpE_elem ec::j_invariant() const
{
    auto &a = _a, &b = _b;
    auto a3 = FpE_elem(4)*a*a*a;
    auto b2 = FpE_elem(27)*b*b;
    return FpE_elem(1728) * a3/(a3 + b2);
}

ec ec::from_j(FpE_elem const &j)
{
    if (j == 0) { return ec(FpE_elem(0), FpE_elem(1));}
    FpE_elem const f(1728);
    if (j == f) { return ec(FpE_elem(1), FpE_elem(0));}
    FpE_elem two(2), three(3);
    auto j2 = j*j;
    auto j3 = j2*j;
    FpE_elem a = three * (f*j - j2);
    FpE_elem b = two * (-j3 + two*f*j2 - f*f*j);
    auto E = ec(a, b);
    assert(E.j_invariant() == j);
    return E;
}

ecp ec::random_point_of_order(Fp2k const &Fext, NTL::ZZ cof, int ell, int e) const {
    ecp P = this->random_point(Fext);
    P = cof*P;
    while (!(NTL::power(NTL::ZZ(ell), e-1)*P)) {
        P = this->random_point(Fext);
        P = cof*P;
    }

    assert (NTL::power(NTL::ZZ(ell), e-1)*P);
    assert (!(NTL::power(NTL::ZZ(ell), e)*P));
    return P;
}

std::pair<ecp, ecp> const ec::torsionBasis(Fp2k const &Fext, int ell, int e) const
{
    // What if we are on the twist?
    NTL::ZZ cof = NTL::power(NTL::ZZ(Fp_elem::modulus()), long(Fext.k)) - NTL::power_long(long(-1), long(Fext.k % 2));

    NTL::ZZ ellcof = NTL::power(NTL::ZZ(ell), e-1);
    NTL::ZZ le = ellcof * ell;
    assert (cof % le == 0);
    cof /= le;

    ecp P = this->random_point_of_order(Fext, cof, ell, e);
    auto ellP = ellcof*P;

    while (true) {
        ecp Q = this->random_point_of_order(Fext, cof, ell, e);
        if (DLP(ellcof*Q, ellP, ell, 1) == NTL::ZZ(-1)) {
            return {P, Q};
        }
    }
}

std::optional<ecp> ec::det_lift_x(Fp2k const &Fext, FpE_elem const &x) const
{
    
    FpE_elem rhs = (x*x + lift(this->_a, Fext))*x + lift(this->_b, Fext);
    auto y = fast_sqrt(Fext, rhs);
    auto p = Fp::modulus(); 
    // FpE_elem p2 = FpE_elem(p);
    if (y) {
        if (!NTL::IsZero(*y) && (2 * NTL::conv<long>(rep(*y)[0]) >= p)) {
            *y = - *y;
        }
        return ecp(std::make_shared<const ec>(*this), Fext, x, *y);
    }
    return {};
}


ecp ec::point_of_order_det(Fp2k const &Fext, NTL::ZZ cof, int ell, int e, long &try_x) const {
    
    bool found = false;
    ecp P = ecp(std::make_shared<const ec>(*this), Fext);
    FpE_elem t = FpE_elem(1) + Fext.Fp2_gen;
    int count = 0;
    ecp P_rec = ecp(std::make_shared<const ec>(*this), Fext);
    while (!found && count < 1000) {
        try_x++;
        count++;
        t._zz_pE__rep[2 * Fext.k - 1] += try_x;
        t._zz_pE__rep[ Fext.k ] += try_x;
        auto Pt = det_lift_x(Fext, t);
        if (Pt) {
            P_rec = NTL::ZZ(1) * (*Pt);
            P = cof * (*Pt);
            found = (NTL::power(NTL::ZZ(ell), e-1) * P).get_z() != 0;
        }
    }
    if (!found) {
        std::cout << "point_order_det looped 1000 times without finding an answer, this is not normal !!\n";
        std::cout << "extension degree = " << Fext.k << " modulus = " << Fp_elem::modulus() << " polymod =  " << FpE_elem::modulus() << "\n";
        std::cout <<  "curve = " << *this << " order = " << ell << " " << e << "\n";
        std::cout << "cof = " << cof << "\n";
        P_rec.normalize();
        P.normalize();
        std::cout << "Pt = " << P_rec << "\n";
        std::cout << "P = " << P << "\n";
        // auto Pt = det_lift_x(Fext, t);
        // if (Pt) {
        //     std::cout << "Pt = " << Pt << "\n";
        //     P = cof * (*Pt);
        //     std::cout << "P =  " << P << "\n";
        // }
        // else {
        //     std::cout << "no Pt \n";
        // }
        // std::cout << try_x;
        assert(0); 
    }

    assert (NTL::power(NTL::ZZ(ell), e-1)*P);
    assert (!(NTL::power(NTL::ZZ(ell), e)*P));

    return P;
}


std::pair<ecp, ecp> const ec::torsionBasisDet(Fp2k const &Fext, int ell, int e) const
{
    // What if we are on the twist?
    NTL::ZZ cof = NTL::power(NTL::ZZ(Fp_elem::modulus()), long(Fext.k)) - NTL::power_long(long(-1), long(Fext.k % 2));

    NTL::ZZ ellcof = NTL::power(NTL::ZZ(ell), e-1);
    NTL::ZZ le = ellcof * ell;
    assert (cof % le == 0);
    cof /= le;

    long try_x = 0;

    FpE_push Push(Fext.F);

    ecp P = point_of_order_det(Fext, cof, ell, e, try_x);
    auto ellP = ellcof * P;

    int count = 0;
    while (true && count < 50) {
        count++;
        ecp Q = this->point_of_order_det(Fext, cof, ell, e, try_x);
        if (DLP(ellcof * Q, ellP, ell, 1) == NTL::ZZ(-1)) {
            return {P, Q};
        }
    }
    // we shouldn't go outside the loop unless there is a problem somewhere else
    // the bound of 50 should be more than enough for the intended sizes of prime p
    assert(0);
    return {P,P};
}

std::vector<ecp> const ec::allTorsionPoints(Fp2k const &Fext, int ell, int e) const
{
    auto Basis = this->torsionBasis(Fext, ell, e);
    assert((NTL::power_ZZ(ell,e)*(Basis.first)).is_identity());
    assert((NTL::power_ZZ(ell,e)*(Basis.second)).is_identity());

    std::vector<ecp> out = {};

    int ell_e = std::pow(ell, e);

    for (int i = 0; i < ell_e; i++) {
        for (int j = 0; j < ell_e; j++) {
            if (not((i % ell == 0) and (j % ell == 0))){
            auto R = i*Basis.first + j*Basis.second;
            out.push_back(R);}
        }
    }
    return out;
}

FpE_elem ec::disc() const
{
    auto &a = _a, &b = _b;
    auto a3 = FpE_elem(4)*a*a*a;
    auto b2 = FpE_elem(27)*b*b;
    return FpE_elem(-16) * (a3 + b2);
}
