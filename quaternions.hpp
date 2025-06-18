///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing the necessary basic quaternion arithmetic
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <cassert>
#include <array>
#include <string>
#include <sstream>
#include <functional>
#include <NTL/ZZ.h>
#include <NTL/HNF.h>
#include <list>
#include <ctime>

typedef NTL::ZZ Integer;

class quatlat;
class quat;

struct quatalg
{
    NTL::ZZ p, q;

    // returns matrices of i and j modulo ell
    std::pair<NTL::mat_ZZ_p,NTL::mat_ZZ_p> splitting(NTL::ZZ const &ell) const
    {
        if (ell == 2)
            throw;  // not implemented
        if (q % ell == 0)
            throw;  // not implemented
        assert(NTL::ZZ_p::modulus() == ell);

        NTL::ZZ_p pmod, qmod, qinv, b, a;
        NTL::conv(pmod, p);
        NTL::conv(qmod, q);
        NTL::inv(qinv, qmod);
        NTL::set(b);

        std::pair<NTL::mat_ZZ_p,NTL::mat_ZZ_p> ret;
        ret.first.SetDims(2, 2);
        ret.second.SetDims(2, 2);
        ret.first[0][1] = -qmod;
        NTL::set(ret.first[1][0]);

        for (size_t i = 0; i < 999; ++i, ++b) {
            if (NTL::IsZero(b)) {
                //Might happen for very small ell
                std::cout << "bruteforcing mod ell splitting data" << std::endl;
                for (size_t a_int = 0; a_int < ell; a_int++) {
                for (size_t b_int = 0; b_int < ell; b_int++) {
                for (size_t c_int = 0; c_int < ell; c_int++) {
                for (size_t d_int = 0; d_int < ell; d_int++) {
                    ret.second[0][0] = a_int;
                    ret.second[0][1] = b_int;
                    ret.second[1][0] = c_int;
                    ret.second[1][1] = d_int;
                    NTL::mat_ZZ_p K = ret.first*ret.second;
                    if (NTL::IsIdent(ret.second*ret.second*NTL::inv(-pmod), 2) && K == -ret.second*ret.first) {
                        assert(NTL::IsIdent(ret.first*ret.first*NTL::inv(-qmod), 2));
                        assert(NTL::IsIdent(ret.second*ret.second*NTL::inv(-pmod), 2));
                        assert(ret.first*ret.second == -ret.second*ret.first);
                        return ret;
                    }
                }
                }
                }
                }
                std::cout << "Really didn't find anything" << std::endl;
                throw;
            }
            auto a2 = -pmod - b*b*qinv;
            if (NTL::Jacobi(NTL::rep(a2), ell) == 1) {
                NTL::conv(a, NTL::SqrRootMod(NTL::rep(a2), ell));
                assert(a*a == a2);
                break;
            }
        }
        assert(!NTL::IsZero(a));

        ret.second[0][0] = a;
        ret.second[0][1] = b;
        ret.second[1][0] = (-pmod-a*a)/b;
        ret.second[1][1] = -a;
        assert(NTL::IsIdent(ret.first*ret.first*NTL::inv(-qmod), 2));
        assert(NTL::IsIdent(ret.second*ret.second*NTL::inv(-pmod), 2));
        assert(ret.first*ret.second == -ret.second*ret.first);
        return ret;
    }

    quatlat maximal_order(bool const &surface=false) const; // return a standard choice of maximal order
    quatlat maximal_order_with_quat_for_check(quat *alpha, bool const &surface=false) const;
};

class quat
{
private:
    std::array<NTL::ZZ,5> coeffs; //the four coefficients, and the denominator

public:
    quatalg const &alg;

    quat(quat const &alpha) : alg{alpha.alg} { *this = alpha; };
    quat(quatalg const &alg_) : alg{alg_} { NTL::set(coeffs[4]); };

    quat(std::array<NTL::ZZ,5> const &coeffs_, quatalg const &alg_)
        : coeffs{coeffs_}, alg{alg_}
    {
        normalize();
    }

    void normalize()
    {
        assert(!NTL::IsZero(coeffs[4]));
        if (coeffs[4] < 0)
            for (auto &c: coeffs)
                c = -c;
        NTL::ZZ g = coeffs[4];
        for (unsigned i = 0; i < 4 && !NTL::IsOne(g); ++i)
            g = NTL::GCD(g, coeffs[i]);
        assert(g > 0);
        if (!NTL::IsOne(g))
            for (auto &c: coeffs)
                c /= g;
    }

    quat const conjugate() const
    {
        return {{coeffs[0], -coeffs[1], -coeffs[2], -coeffs[3], coeffs[4]}, alg};
    }

    // in-place
    void invert()
    {
        auto conj = conjugate();
        auto nrm = norm();
        coeffs[0] = conj[0]*nrm.second;
        coeffs[1] = conj[1]*nrm.second;
        coeffs[2] = conj[2]*nrm.second;
        coeffs[3] = conj[3]*nrm.second;
        coeffs[4] = conj[4]*nrm.first;
        normalize();
 }

    NTL::ZZ &operator[](size_t index)
    {
        if (index > 4)
            throw std::out_of_range("Index out of bounds");
        return coeffs[index];
    }

    NTL::ZZ const &operator[](size_t index) const
    {
        if (index > 4)
            throw std::out_of_range("Index out of bounds");
        return coeffs[index];
    }

    quat operator+(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        auto const &[a,b,c,d,e] = coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;
        return {{a*E+A*e, b*E+B*e, c*E+C*e, d*E+D*e, e*E}, alg};
    };

    quat &operator+=(quat const &other)
    {
        auto const sum = *this + other;
        std::copy(sum.coeffs.begin(), sum.coeffs.end(), this->coeffs.begin());
        return *this;
    };

    quat operator*(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        auto const &[a,b,c,d,e] = this->coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;
        auto const r = a*A - alg.q*b*B - alg.p* (c*C + alg.q*d*D);
        auto const s = a*B + A*b + alg.p* ( c*D - C*d);
        auto const t = a*C + A*c - alg.q*b*D + alg.q*B*d;
        auto const u = a*D + A*d + b*C - B*c;
        return {{r,s,t,u, e*E}, alg};
    };

    quat operator*(NTL::ZZ const &n) const { return {{n*coeffs[0], n*coeffs[1], n*coeffs[2], n*coeffs[3], coeffs[4]}, alg}; }
    quat operator*(long n) const { return {{n*coeffs[0], n*coeffs[1], n*coeffs[2], n*coeffs[3], coeffs[4]}, alg}; }

    quat& operator=(const quat& other) {
        if (this != &other) {  // Avoid self-assignment
            // Copy the values
            if (&this->alg != &other.alg) {
                throw std::runtime_error("LOUD ERROR WHEN ASSIGNING QUATLAT");
            }
            coeffs = other.coeffs;
        }
        return *this;
    }

    std::pair<NTL::ZZ,NTL::ZZ> norm() const  // returns numerator and denominator
    {
        auto const &[a,b,c,d,e] = this->coeffs;
        return { a*a + alg.q * b*b + alg.p * ( c*c + alg.q * d * d), e*e };
    }

    std::pair<NTL::ZZ,NTL::ZZ> trace() const  // returns numerator and denominator
    {
        auto num = 2*coeffs[0], denom = coeffs[4];
        auto g = NTL::GCD(num, denom);
        if (!NTL::IsOne(g)) {
            num /= g;
            denom /= g;
        }
        return {num, denom};
    }

    std::pair<NTL::ZZ,NTL::ZZ> trace_pairing(quat const &other) const
    {
        if (&other.alg != &alg) // object identity
            throw;
        return (this->conjugate() * other).trace();
    }

    friend std::ostream& operator<<(std::ostream& o, quat const &alpha)
    {
            return o
            << alpha[0] << "/" << alpha[4] << " + "
            << alpha[1] << "/" << alpha[4] << "*i + "
            << alpha[2] << "/" << alpha[4] << "*j + "
            << alpha[3] << "/" << alpha[4] << "*k";
    }

    bool is_one() const
    {
        auto const &[a,b,c,d,e] = coeffs;
        assert(e != 0);
        return a == e && b == 0 && c == 0 && d == 0;
    }

    bool is_zero() const
    {
        auto const &[a,b,c,d,e] = coeffs;
        assert(e != 0);
        return a == 0 && b == 0 && c == 0 && d == 0;
    }
};


#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

class quatlat
{
public:
    quatalg const &alg;
    NTL::mat_ZZ basis;
    NTL::ZZ denom;

private:
    mutable std::pair<NTL::ZZ,NTL::ZZ> the_norm {};
    mutable quat gen;


public:
    quatlat(NTL::mat_ZZ const &gens, NTL::ZZ const &denom_, quatalg const &alg_);

    quatlat(NTL::mat_ZZ const &gens, std::pair<NTL::ZZ,NTL::ZZ> const &norm, NTL::ZZ const &denom_, quatalg const &alg_);

    quatlat(NTL::mat_ZZ const &gens, NTL::ZZ const &denom_, quatalg const &alg_, bool normalize);

    quatlat(NTL::mat_ZZ const &gens, std::pair<NTL::ZZ,NTL::ZZ> const &norm, NTL::ZZ const &denom_, quatalg const &alg_, bool normalize);

    quatlat(quatlat const &L) : alg{L.alg}, gen{L.gen} { *this = L; }

    void normalize(bool safe = false);

    void reset_norm() const {
        the_norm.first = 0;
        the_norm.second = 0;
    }

    void _conjugate();
    quatlat conjugate() const;
    quatlat HNF_conjugate() const;

    quatlat inverse() const;

    quatlat operator+(quatlat const &other) const;

    quatlat operator*(quatlat const &other) const;

    quatlat operator*(quat const &other) const;

    quatlat operator*(NTL::ZZ const &other) const;

    friend quatlat operator*(quat const &left, quatlat const &right)
    { return (right.conjugate()*left.conjugate()).conjugate(); }  //TODO write faster version

    friend quatlat operator*(NTL::ZZ const &left, quatlat const &right)
    { return right * left; }

    quatlat& operator=(const quatlat& other) {
        if (this != &other) {  // Avoid self-assignment
            // Copy the values
            if (&this->alg != &other.alg) {
                throw std::runtime_error("LOUD ERROR WHEN ASSIGNING QUATLAT");
            }
            basis = other.basis;
            denom = other.denom;
            the_norm = other.the_norm;
        }
        return *this;
    }

    NTL::mat_ZZ normalized_invariant() const;

    quatlat conjugated_by_j() const;

    bool operator==(quatlat const &other) const;

    std::pair<NTL::ZZ,NTL::ZZ> norm() const;  // returns numerator and denominator


    std::pair<NTL::ZZ,NTL::ZZ> norm_no_comput() const  // returns numerator and denominator
    {
        return the_norm;
    }

    // in-place
    void _intersect(quatlat const &other, bool safe = false);
    void _fast_intersect(quatlat const &other);

    // set of all x such that xL <= L or Lx <= L
    quatlat _compute_order(bool right_order) const;

    quatlat new_right_order() const;
    std::pair<quatlat,NTL::mat_ZZ>  fast_right_order_and_gram();

    quatlat left_order() const { return _compute_order(false); }
    quatlat right_order() const { return new_right_order(); }
    // NB! ^definitely buggy some times
    //{ return _compute_order(true); }

    bool is_order() const;

    // Using https://math.dartmouth.edu/~jvoight/articles/73446.pdf Lemma 7.2:
    void right_ideals_of_norm(NTL::ZZ const &ell, std::function<void(quatlat const &)> const &fun);


    std::list<quatlat> left_ideals_of_prime_norm(NTL::ZZ const &ell, const quat &first_gen, const quat &iter);

    void left_ideals_of_norm(NTL::ZZ const &ell, std::function<void(quatlat const &)> const &fun)
    {
        auto const fun2 = [&](quatlat const &lat) {
            fun(lat.conjugate());
        };
        right_ideals_of_norm(ell, fun2);
    }

    friend std::ostream& operator<<(std::ostream& o, quatlat const &L)
    {
        return o << L.basis << "/" << L.denom << "  (norm " << L.norm().first << "/" << L.norm().second << ")";
    }

    std::string sage() const
    {
        std::ostringstream s;
        s << "i.parent().ideal([";
        for (unsigned i = 0; i < 4; ++i) {
            s << "[";
            for (unsigned j = 0; j < 4; ++j)
                s << basis[i][j] << "/" << denom << ",";
            s << "],";
        }
        s << "])";
        return s.str();
    }

    NTL::ZZ get_denom() const { return denom; }

    NTL::mat_ZZ get_basis() const { return basis; }

    NTL::mat_ZZ LLL_basis() const;

    void enumerate_shortish(unsigned long bnd, std::function<bool(quat const &)> const &fun) const;
    quat compute_gen_cyclic() const;

    NTL::mat_ZZ HNF_basis() const;

    quatlat copy() const;

    void set_generator(quat const &gen);

    quat get_generator() const;

    bool contains(quat const &alpha) const;
    bool fast_contains(quat const &alpha, Integer const &modnorm) const;

    void reduce_norm_cyclic(const quatlat &left_order, const NTL::ZZ &norm);

    quatlat special_add(const quatlat &other) const;

    friend class quatlatenum;
};

quatlat create_from_generator(const quat &gen, const NTL::ZZ &norm, const quatlat &left_order);
std::list<quatlat> left_ideals_of_prime_norm_O0(NTL::ZZ const &ell, const quat &first_gen, const quat &iter);
quatlat create_from_generator_O0(const quat &gen, const NTL::ZZ &norm);

quatlat connecting_ideal(const quatlat &left_order, const quatlat &right_order);


namespace std {
    template<>
    struct hash<quatlat> {
        size_t operator()(const quatlat& L) const {
            std::hash<long> hashlong;

            auto basis_mat = L.normalized_invariant();

            size_t h = 0xffffff;
            for (unsigned i = 0; i < 4; ++i) {
                for (unsigned j = 0; j < 4; ++j) {
                    long elem_long;
                    NTL::conv(elem_long, basis_mat[i][j]);
                    h ^= hashlong(elem_long);
                }
            }

            return h;
        }
    };
}
