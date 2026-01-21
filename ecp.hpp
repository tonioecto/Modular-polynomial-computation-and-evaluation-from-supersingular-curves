//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// Elliptic curve point class (header)
//////////////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <cassert>
#include <optional>
#include <memory> // Include the memory header for shared_ptr

#include <NTL/ZZ_pE.h>
#include <NTL/mat_ZZ_p.h>
#include "Fp2k.hpp"
#include "utils.hpp"
#include "smallint.hpp"

class ecp {
    private:
        std::shared_ptr<const ec> E;
        std::reference_wrapper<const Fp2k> Fext;
        mutable FpE_elem x, y, z;
    public:
        // Given only a curve and field extension, it creates point at infinity
        // ecp() {};
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext) : E{curve}, Fext{Fext}, x{0}, y{1}, z{0} {};
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, FpE_elem const &x, FpE_elem const &y, FpE_elem const &z);
        ecp(std::shared_ptr<const ec> curve, Fp2k const &Fext, FpE_elem const &x, FpE_elem const &y);
        

        FpE_elem const &get_x() const {return this->x;}
        FpE_elem const &get_y() const {return this->y;}
        FpE_elem const &get_z() const {return this->z;}
        void set_x(const FpE_elem &_x) {this->x = _x;}
        void set_y(const FpE_elem &_y) {this->y = _y;}
        void set_z(const FpE_elem &_z) {this->z = _z;}


        FpE_elem const aff_x() const {if (this->is_identity()) {return FpE_elem(0);} FpE_push push(this->field().F); return this->x/this->z;}
        FpE_elem const aff_y() const {if (this->is_identity()) {return FpE_elem(1);} FpE_push push(this->field().F); return this->y/this->z;}

        bool is_identity() const {return NTL::IsZero(this->z);}

    friend std::ostream& operator<<(std::ostream& o, ecp const &P)
    {
        return o << "(" << NTL::rep(P.x)
                 << " : " << NTL::rep(P.y)
                 << " : " << NTL::rep(P.z)
                 << ")";
    }

    ec const &curve() const { return *E; }
    std::shared_ptr<const ec> const &curve_ptr() const { return E; }
    Fp2k const &field() const { return this->Fext.get(); }

    void normalize() const;
    ecp const &normalized() { this->normalize(); return *this;}

    std::pair<FpE_elem,FpE_elem> _lambda(ecp const &other) const;
    ecp operator+(ecp other) const;
    friend ecp operator*(NTL::ZZ k, ecp P);
    friend ecp operator*(int k, ecp P) { return NTL::ZZ(k)*P; }
    friend ecp operator*(long k, ecp P) { return NTL::ZZ(k)*P; }

    ecp const frob() const
    {
        auto const &fld = field();
        FpE_push push(fld.F);
        return {this->E, Fext, fld.frob(x), fld.frob(y), fld.frob(z)};
    }


    bool operator==(ecp const &other) const
    {
        return (this->curve() == other.curve()
        && this->aff_x() == other.aff_x()
        && this->aff_y() == other.aff_y()
        && this->is_identity() == other.is_identity());
    }

    ecp &operator+=(ecp const &other) { return *this = *this + other; }
    ecp operator-() const { ecp Q = *this; Q.y = -Q.y; return Q; }
    ecp operator-(ecp const &other) const { return *this + (-other); }
    ecp &operator-=(ecp const &other) { return *this = *this - other; }
    operator bool() const { return !(this->is_identity());}
};



namespace std {
    template<>
    struct hash<FpE_elem>
    {
        size_t operator()(FpE_elem const &a) const {
            size_t h = 0xb00b;
            if (a == 0) {
                return h;
            }
            auto apoly = NTL::rep(a);
            size_t d = NTL::deg(apoly);
            for (size_t i = 0; i <= d; i++) {
                size_t c;
                conv(c, NTL::coeff(apoly, i));
                h ^= c;
            }
            return h;
        }
    };
}

namespace std {
    template<>
    struct hash<ecp>
    {
        size_t operator()(ecp const &P) const {
            size_t h = 0xec77;
            h += std::hash<FpE_elem>()(P.curve().a());
            h += std::hash<FpE_elem>()(P.curve().b());

            h ^= std::hash<FpE_elem>()(P.aff_x());
            h *= 1111111;
            h ^= std::hash<FpE_elem>()(P.aff_y());
            return h;
        }
    };
}

