///////////////////////////////////////////////////////////////////////////////////////////////
////   Code to generate the shortest vectors in a quaternion lattice
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <cassert>

#include <functional>

#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include "quaternions.hpp"

#include <gmp.h>
#include <fplll/fplll.h>

void ntl2gmp(mpz_t out, NTL::ZZ const &num);
void gmp2ntl(NTL::ZZ &out, mpz_t const num);

using mat_t = fplll::ZZ_mat<mpz_t>;
using gso_t = fplll::MatGSOGram<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>>;

inline gso_t quatlat2fplll(mat_t &G, mat_t &U, quatlat const &lat)
{
    //std::cerr << "p = " << lat.alg.p << std::endl;
    //std::cerr << "q = " << lat.alg.q << std::endl;
    //std::cerr << lat.sage() << std::endl;

    std::array<quat, 4> elts {{
        {{lat.basis[0][0], lat.basis[0][1], lat.basis[0][2], lat.basis[0][3], NTL::ZZ(1)}, lat.alg},
        {{lat.basis[1][0], lat.basis[1][1], lat.basis[1][2], lat.basis[1][3], NTL::ZZ(1)}, lat.alg},
        {{lat.basis[2][0], lat.basis[2][1], lat.basis[2][2], lat.basis[2][3], NTL::ZZ(1)}, lat.alg},
        {{lat.basis[3][0], lat.basis[3][1], lat.basis[3][2], lat.basis[3][3], NTL::ZZ(1)}, lat.alg},
    }};

    G.resize(4, 4);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = i; j < 4; ++j) {
            // auto val = elts[i] * elts[j].conjugate();
            // auto pair = val[0];
            auto pair = elts[i][0] * elts[j][0] + lat.alg.q * elts[i][1] * elts[j][1] + lat.alg.p * (elts[i][2] * elts[j][2] + lat.alg.q * elts[i][3] * elts[j][3]);
            ntl2gmp(G[i][j].get_data(), pair+pair);
            if (i != j) {
                ntl2gmp(G[j][i].get_data(), pair+pair);
            }
    }

    // the Gram matrix of an integral ideal is always integral
    {
        mpz_t d;
        mpz_init(d);
        ntl2gmp(d, lat.denom);
        mpz_mul(d, d, d);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j) {
                auto &c = G[i][j].get_data();
                assert(mpz_divisible_p(c, d));
                mpz_divexact(c, c, d);
            }
        mpz_clear(d);
    }

//    std::cerr << "G:\n" << G << std::endl;

    U.gen_identity(4);
    mat_t Uinv;
    gso_t gso(G, U, Uinv, fplll::GSO_INT_GRAM);

    fplll::LLLReduction<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> lllobj(gso, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, 0);
    lllobj.lll();

    //    std::cerr << "G':\n" << G << std::endl;
    //    std::cerr << "U:\n" << UU << std::endl;
    //    std::cerr << "U^-1:\n" << U << std::endl;

    return gso;
}

class quatlatenum
{
    public:
    quatlat const &lat;

    mat_t gram;
    mat_t U;
    gso_t gso;

    std::vector<quat> basis;


    quatlatenum(quatlat const &_lat)
        : lat {_lat}, gso {quatlat2fplll(gram, U, lat)}
    {
        for (unsigned i = 0; i < 4; ++i) {
            NTL::Vec<NTL::ZZ> row;
            row.SetLength(4);
            for (unsigned j = 0; j < 4; ++j) {
                NTL::ZZ c;
                gmp2ntl(c, U[i][j].get_data());
                row += c * lat.basis[j];
            }
            basis.emplace_back(std::array<NTL::ZZ,5>{row[0], row[1], row[2], row[3], lat.denom}, lat.alg);
        }
    }

    void enumerate(double bnd, std::function<bool(quat const &)> const &fun)
    {
        auto const callback = [&](size_t n, fplll::enumf *coords, void *ctx) -> bool
        {
            (void) ctx;
            assert(n == 4);
            // std::cerr << n << ": ";
            (void) n;
            // std::cerr << coords[0] << " " << coords[1] << " " << coords[2] << " " << coords[3] << " ~> ";
            quat elt {lat.alg};
            for (unsigned i = 0; i < 4; ++i)
                elt += basis[i] * (long) coords[i];
                // std::cerr << elt << " | ";
                // std::cerr << elt.norm().first << std::endl;
            return fun(elt);
        };
        fplll::CallbackEvaluator<fplll::FP_NR<double>> evaluator(callback, nullptr, 1, fplll::EVALSTRATEGY_FIRST_N_SOLUTIONS);
        fplll::Enumeration<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> enumobj(gso, evaluator);

        fplll::FP_NR<double> maxdist;
        gso.get_r(maxdist, 3, 3);
        maxdist *= bnd;
        enumobj.enumerate(0, 4, maxdist, 0);
    }
};


inline bool isPrincipal(quatlat const &J) {
    // clock_t t = clock();
    quatlatenum Enumerator(J);
    // std::cout << "isPrinc time = " << (double) (clock() -t )/CLOCKS_PER_SEC <<  " \n";
    auto &c = Enumerator.gram[0][0].get_data();
    NTL::ZZ norm;
    gmp2ntl(norm, c);
    bool isPrincipal = (norm * J.norm().second == 2 * J.norm().first);
    return isPrincipal;
}
