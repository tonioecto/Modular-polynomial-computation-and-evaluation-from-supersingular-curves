///////////////////////////////////////////////////////////////////////////////////////////////
/////// The code for the elliptic curve class (header)
///////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <optional>

#include <NTL/ZZ_pE.h>
#include "Fp2k.hpp"
#include <vector>

class ecp;

class ec
{
    private:
        FpE_elem _a, _b;
    public:
        ec() {};
        ec(FpE_elem const &a, FpE_elem const &b);

        FpE_elem const &a() const { return this->_a; }
        FpE_elem const &b() const { return this->_b; }
        // Kohel: b2 = 0, b4 = 2*a, b6 = 4*b

    std::optional<ecp> lift_x(Fp2k const &Fext, FpE_elem const &x) const;

    ecp random_point(Fp2k const &Fext) const;
    ecp random_point_of_order(Fp2k const &Fext, NTL::ZZ cof, int ell, int k) const;

    bool operator==(ec const &other) const { return this == &other || (this->_a == other._a && this->_b == other._b); }
    bool operator!=(ec const &other) const { return !(*this == other); }

    ecp operator()(Fp2k const &Fext, FpE_elem new_x, FpE_elem new_y) const;
    ecp operator()(Fp2k const &Fext) const;



    FpE_elem j_invariant() const;
    FpE_elem disc() const;

    ec static from_j(FpE_elem const &j);

    std::pair<ecp, ecp> const torsionBasis(Fp2k const &Fext, int ell, int k) const;
    std::vector<ecp> const allTorsionPoints(Fp2k const &Fext, int ell, int k) const;
    friend std::ostream& operator<<(std::ostream& o, ec const &E) { return o << "{y^2 = x^3 + (" << NTL::rep(E._a) << ")*x + (" << NTL::rep(E._b) << ")}"; }

};
