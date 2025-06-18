///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements the KLPT algorithm (and variants)
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

#include "ecp.hpp"
#include "choosetorsion.hpp"
#include "utils.hpp"
#include <gmp.h>
#include <NTL/mat_ZZ.h>

#include <random>
#include "quaternions.hpp"
#include "quatlatenum.hpp"

void ntl2gmp(mpz_t out, NTL::ZZ const &num);
void gmp2ntl(NTL::ZZ &out, mpz_t const num);

#include <fplll/fplll.h>

using mat_t = fplll::ZZ_mat<mpz_t>;
using gso_new_t = fplll::MatGSO<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>>;

// enumerates small solutions (x,y) to A*x + B*y = C (mod N)
void small_solutions(NTL::ZZ const &N, NTL::ZZ const &C, NTL::ZZ const &A, NTL::ZZ const &B,
                     std::function<bool(NTL::ZZ const &x, NTL::ZZ const &y)> const &fun, double maxdist)
{
    //    std::cerr << "N = " << N << std::endl;
    //    std::cerr << "C = " << C << std::endl;
    //    std::cerr << "A = " << A << std::endl;
    //    std::cerr << "B = " << B << std::endl;

    mat_t M, U, Uinv;

    // TODO: compute a basis of this 2-dimensional lattice directly instead of LLL'ing it
    M.resize(3, 2);
    ntl2gmp(M[0][0].get_data(), B);
    ntl2gmp(M[0][1].get_data(), -A);
    ntl2gmp(M[1][0].get_data(), N);
    M[2][1] = M[1][0];
    // std::cerr << M << std::endl;
    fplll::lll_reduction(M);
    assert(M[0][0].sgn() == 0 && M[0][1].sgn() == 0);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
            M[i][j] = M[i+1][j];
    M.resize(2, 2);
    // std::cerr << M << std::endl;

    gso_new_t gso(M, U, Uinv, 0);
    gso.update_gso();

    NTL::mat_ZZ L; L.SetDims(2, 2);
    for (unsigned i = 0; i < 2; ++i)
        for (unsigned j = 0; j < 2; ++j)
            gmp2ntl(L[i][j], M[i][j].get_data());
    // std::cerr << L << std::endl;

    auto CdivB = C * NTL::InvMod(B, N) % N;
    std::vector<fplll::FP_NR<double>> target {0., 0.};
    {
        auto a = NTL::to_RR(L[0][0]);
        auto b = NTL::to_RR(L[0][1]);
        auto c = NTL::to_RR(L[1][0]);
        auto d = NTL::to_RR(L[1][1]);
        NTL::RR u, v;
        if (a == 0) {
            u = -NTL::to_RR(CdivB) / b;
            assert(v == 0); //XXX
        }
        else {
            v = NTL::to_RR(CdivB) / (c/a*b - d);
            u = -c/a * v;
        }
        target[0] = NTL::to_double(u);
        target[1] = NTL::to_double(v);
    }

    auto const callback = [&](size_t n, fplll::enumf *coords, void *ctx) -> bool
    {
        (void) ctx;
        assert(n == 2); (void) n;
        // std::cerr << coords[0] << " " << coords[1] << " ~> ";
        NTL::ZZ x = coords[0]*L[0][0] + coords[1]*L[1][0];
        NTL::ZZ y = coords[0]*L[0][1] + coords[1]*L[1][1];
        y += CdivB;  // CVP
        // std::cerr << x << " " << y << "\n";
        assert((A*x + B*y - C) % N == 0);
        return fun(x, y);
    };
    fplll::CallbackEvaluator<fplll::FP_NR<double>> evaluator(callback, nullptr, 1, fplll::EVALSTRATEGY_BEST_N_SOLUTIONS);  //TODO or EVALSTRATEGY_FIRST_N_SOLUTIONS?
    fplll::Enumeration<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> enumobj(gso, evaluator);

    fplll::FP_NR<double> md = maxdist;
    enumobj.enumerate(0, 2, md, 0, target);
}

std::optional<std::pair<NTL::ZZ, NTL::ZZ>> Cornacchia(const NTL::ZZ d,const NTL::ZZ &M) {
    ////////////////////////////////////////////////
    // Finds x, y such that x**2 + d*y**2 = M
    //      also outputs indicator of success
    //      only run on M prime
    ////////////////////////////////////////////////

    NTL::ZZ D = M - d;

    if (NTL::ProbPrime(M) && (NTL::Jacobi(D, M) == 1)) {
        NTL::ZZ r = NTL::SqrRootMod(D, M);
        NTL::ZZ rs[2];
        rs[0] = r;
        rs[1] = M-r;

        NTL::ZZ Bound = SqrRoot(M);

        for (auto r : rs) {
            NTL::ZZ n = M;
            while (r >= Bound) {
                NTL::ZZ temp = r;
                r = n % r;
                n = temp;
            }
            NTL::ZZ s = NTL::SqrRoot((M-r*r)/d);

            if (r*r + d*s*s == M) {
                return std::make_pair(r, s);
            }
        }
    }

    return {};
}

quat RepresentInteger(quatalg const &B, NTL::ZZ const &M) {// Don't need the full version
    //////////////////////////////////////////////////////////////////////
    // Finds a quaternion of the form x + iy + jz + kw = M, in B
    //////////////////////////////////////////////////////////////////////

    auto p = B.p;
    auto q = B.q;
    assert (M > 50*q*p);

    std::random_device rd;
    std::mt19937 gen(rd());
    long firstbound = 10000;
    NTL::ZZ better_bound = NTL::SqrRoot(M/(2*p*q));
    if (better_bound < firstbound) {
        NTL::conv(firstbound, better_bound);
    }
    //std::cout << "REPINT bound: " << firstbound << std::endl;
    std::uniform_int_distribution<> distribution(0, firstbound); // This shouldn't really be hardcoded

    NTL::ZZ z, w, tar;
    for (size_t i = 1 ; i < 10000 ; i++) {
        NTL::conv(z, distribution(gen));
        NTL::conv(w, distribution(gen));

        tar = M - p*z*z - p*q*w*w;

        if (tar < 0) {
            std::cout << "REPINT: TOO SMALL NUM" << std::endl; // This really shouldn't happen
            continue;
        }

        auto result = Cornacchia(q, tar);
        if (result) {

            std::array<NTL::ZZ,5> coeffs;
            coeffs[0] = result->first;
            coeffs[1] = result->second;
            coeffs[2] = z;
            coeffs[3] = w;
            coeffs[4] = NTL::ZZ(1);
            return quat(coeffs, B);
        }
    }
    //In case nothing was found
    std::cout << "REPINT: failed to find a solution!" << std::endl;
    std::array<NTL::ZZ,5> coeffs;
    coeffs[0] = NTL::ZZ(0);
    coeffs[1] = NTL::ZZ(0);
    coeffs[2] = NTL::ZZ(0);
    coeffs[3] = NTL::ZZ(0);
    coeffs[4] = NTL::ZZ(1);
    return quat(coeffs, B);
}

std::pair<NTL::ZZ, NTL::ZZ> IdealModConstraint(quatlat const &I, quat const &gamma) {


    quat i({{NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, gamma.alg});
    quat j({{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, gamma.alg});

    //NTL::mat_ZZ M = NTL::transpose(I.get_basis());
    NTL::mat_ZZ M = I.get_basis();

    auto b1 = gamma*j;
    auto b2 = gamma*i*j;

    NTL::ZZ d1 = b1[4];
    NTL::ZZ d2 = b2[4];
    NTL::ZZ d3 = I.get_denom();

    NTL::mat_ZZ system;
    system.SetDims(6, 4);

    for (size_t i = 0; i < 4; i++) {
        system[0][i] = b1[i]*d2*d3;
        system[1][i] = b2[i]*d1*d3;
        system[2][i] = M[0][i]*d1*d2;
        system[3][i] = M[1][i]*d1*d2;
        system[4][i] = M[2][i]*d1*d2;
        system[5][i] = M[3][i]*d1*d2;
    }

    assert(I.norm().second == 1);

    NTL::ZZ N = I.norm().first;

    assert(NTL::ProbPrime(N));

    NTL::ZZ_pPush push(N);
    NTL::mat_ZZ_p system_reduced;
    NTL::conv(system_reduced, system);

    NTL::mat_ZZ_p X;
    NTL::kernel(X, system_reduced);

    std::pair<NTL::ZZ, NTL::ZZ> CD(NTL::rep(X[0][0]), NTL::rep(X[0][1]));

    return CD;
}

std::optional<quat> StrongApproximation(NTL::ZZ const &lam, NTL::ZZ const &C, NTL::ZZ const &D, NTL::ZZ const &N_I, NTL::ZZ const &N_mu, quatalg const &B) {
    NTL::ZZ p = B.p;
    NTL::ZZ q = B.q;

    NTL::ZZ lamC = (lam*C) % N_I;
    NTL::ZZ lamD = (lam*D) % N_I;

    NTL::ZZ lhs = N_mu - (p*lamC*lamC + p*q*lamD*lamD);

    assert (lhs % N_I == 0);
    lhs /= N_I;

    // Petit--Smith
    //  https://crypto.iacr.org/2018/affevents/mathcrypt/medias/08-50_3.pdf
    std::optional<quat> sol;

    small_solutions(N_I, lhs * NTL::InvMod(2*B.p*lam % N_I, N_I) % N_I, C, B.q*D % N_I, [&](NTL::ZZ const &z, NTL::ZZ const &w) {
        NTL::ZZ plamCz2 = p*lamC*z*2;
        NTL::ZZ M = lhs - (plamCz2 + p*(N_I*z*z + q*(lamD*w*2 + N_I*w*w)));
        //std::cout << "M = " << M << std::endl;
        if (M > 0) {
            assert (M % N_I == 0);
            auto result = Cornacchia(q, M/(N_I));
            if (result) {
                std::array<NTL::ZZ,5> coeffs {N_I*result->first, N_I*result->second, lamC + N_I*z, lamD + N_I*w, NTL::ZZ(1)};
                sol.emplace(coeffs, B);
                return true;
            }
        }
        return false;
    }, 99*NTL::to_double(N_I));

    return sol;
}

quatlat equivPrime(quatlat const &J, NTL::ZZ const &T) {
    quatlat O0 = J.left_order();
    quatlat I = O0;

    bool found = false;

    J.enumerate_shortish(100, //TODO <- good value here?
        [&](quat const &el) -> bool {
        //std::cout << "el norm: " << el.norm().first << "/" << el.norm().second << std::endl;
        //std::cout << "el should be an element of J..." << std::endl;
        //std::cout << "J norm: " << J.norm().first << "/" << J.norm().second << std::endl;
        NTL::ZZ N = el.norm().first*J.norm().second/(el.norm().second*J.norm().first);
        assert (N > 1); //Principal ideal
        // std::cout << "outside N = " << N << std::endl;
        if ((NTL::GCD(N, T) == 1) && (NTL::ProbPrime(N)) && (NTL::GCD(N*N, el.norm().first/el.norm().second) == N)) { //Test that el is also a generator. TODO: Find new generator in this case
            // std::cout << "el used: " << std::endl;
            //std::cout << el << std::endl;
            auto el_conj = el.conjugate();
            I = create_from_generator(el_conj, N, O0);
            found = true;
            return true;
        }
        return false;
    });
    I.reset_norm();
    assert (found);
    //std::cout << "Size of equiv ideal: " << NTL::NumBits(I.norm().first) << ", p: " << NTL::NumBits(O0.alg.p) << std::endl;
    return I;
}

quatlat equivPrime_conj(quat *conj, quatlat const &J, NTL::ZZ const &T) {
    quatlat O0 = J.left_order();
    quatlat I = O0;

    bool found = false;


    J.enumerate_shortish(100, //TODO <- good value here?
        [&](quat const &el) -> bool {
        //std::cout << "el norm: " << el.norm().first << "/" << el.norm().second << std::endl;
        //std::cout << "el should be an element of J..." << std::endl;
        //std::cout << "J norm: " << J.norm().first << "/" << J.norm().second << std::endl;
        NTL::ZZ N = el.norm().first*J.norm().second/(el.norm().second*J.norm().first);
        assert (N > 1); //Principal ideal
        // std::cout << "outside N = " << N << std::endl;
        if ((NTL::GCD(N, T) == 1) && (NTL::ProbPrime(N)) && (NTL::GCD(N*N, el.norm().first/el.norm().second) == N)) { //Test that el is also a generator. TODO: Find new generator in this case
            // std::cout << "el used: " << std::endl;
            //std::cout << el << std::endl;
            auto el_conj = el.conjugate();
            I = create_from_generator(el_conj, N, O0);
            (*conj)[0] = el_conj[0];
            (*conj)[1] = el_conj[1];
            (*conj)[2] = el_conj[2];
            (*conj)[3] = el_conj[3];
            (*conj)[4] = J.norm().first;
            found = true;
            return true;
        }
        return false;
    });
    I.reset_norm();
    assert (found);
    //std::cout << "Size of equiv ideal: " << NTL::NumBits(I.norm().first) << ", p: " << NTL::NumBits(O0.alg.p) << std::endl;
    return I;
}

NTL::ZZ _choose_gamma_norm(factor_list const &fac_list, NTL::ZZ const &coprime_in, NTL::ZZ const &N_I, NTL::ZZ const &p, NTL::ZZ const &q, NTL::ZZ const &Smallest_2mod3, NTL::ZZ const &Second_smallest_2_mod_3) {
    NTL::ZZ min_size = 1000*q*p;

    NTL::ZZ N_gamma(1);
    NTL::ZZ N_gammaN_I;
    NTL::ZZ coprime = coprime_in;
    if (q == 3) {
        coprime *= Second_smallest_2_mod_3;
    }
    for (const auto& tup : fac_list) {
        //std::cout << "----" << std::endl;
        //std::cout << q << std::endl;
        //std::cout << min_size << std::endl;
        //std::cout << N_gamma << std::endl;
        NTL::ZZ ell(std::get<0>(tup));
        if ((2*q)%ell == 0) { //Add so that N_gamma | T
            continue;
        }
        if (NTL::GCD(coprime, ell) != 1) {
            continue;
        }

        N_gamma *= ell;
        N_gammaN_I = N_gamma*N_I;

        if (N_gammaN_I > min_size) {
            if (q != 3) {
                break;
            } else {
                if (N_gammaN_I % 3 == 1) {
                    break;
                } else if ((N_gammaN_I/Smallest_2mod3 > min_size) && (N_gamma % Smallest_2mod3 == 0)){
                    N_gamma /= Smallest_2mod3;
                    N_gammaN_I = N_gamma*N_I;
                    break;
                } else if ((N_gamma % Smallest_2mod3 != 0)) {
                    N_gamma *= Smallest_2mod3;
                    N_gammaN_I = N_gamma*N_I;
                    break;
                }
            }
        }
    }
    assert (N_gammaN_I > min_size);
    // Shaving off some excess
    //std::cout << "N_gamma" << N_gamma << std::endl;
    //std::cout << "Jacobi" << NTL::Jacobi(N_gamma*N_I, p) << std::endl;
    //std::cout << "min_size" << min_size << std::endl;
    NTL::ZZ d;
    NTL::div(d, N_gammaN_I, min_size);
    if (d > 10000) {
        d = NTL::ZZ(10000);
    }
    // Need to be careful with jacobi symbol mod q
    while (d > 1) {
        if (N_gamma % d != 0) {
            d--;
            continue;
        }
        //if (true) {
            // Okay to divide freely
        if (q != 3) {
            while ((N_gamma % d == 0) && ((N_gamma*N_I)/d > min_size)) {
                //std::cout << "loop 1" << std::endl;
                N_gamma /= d;
            }
        } else {
            // Only divide squares
            while ((N_gamma % (d*d) == 0) && ((N_gamma*N_I)/(d*d) > min_size)) {
                //std::cout << "loop 2" << std::endl;
                N_gamma /= (d*d);
            }
        }
        d--;
    }
    if (q == 3) {
        assert ((N_gamma*N_I) % 3 == 1);
    }

    return N_gamma;
}

std::pair<quat, NTL::ZZ> _makePrimitive(quat const &alpha, quatlat const &O, NTL::ZZ const &N_I) {
    //Make alpha primitive, and update its norm
    O.reset_norm();
    auto O_basis = O.get_basis();
    auto O_denom = O.get_denom();

    assert (alpha[4] == 1);

    NTL::vec_ZZ alpha_vec;
    alpha_vec.SetLength(4);
    for (size_t i = 0; i < 4 ; i ++) {
        alpha_vec[i] = alpha[i];
    }
    //std::cout << "O0: " << std::endl;
    //std::cout << O.sage() << std::endl;
    //std::cout << "MAKE PRIMITIVE: alpha = " << alpha << std::endl;
    NTL::vec_ZZ coeffs;
    NTL::ZZ d;
    NTL::solve(d, coeffs, O_basis, alpha_vec, true); //Use deterministic strategy that doesnt fail...
    // This solves d*alpha in the 2*q*O_basis
    // coeffs % d should be 0
    NTL::ZZ div = NTL::GCD(NTL::GCD(NTL::GCD(coeffs[0], coeffs[1]), coeffs[2]), coeffs[3]);
    //std::cout << "MAKE PRIMITIVE: d = " << d << std::endl;
    //std::cout << "MAKE PRIMITIVE: div = " << div << std::endl;
    //std::cout << "MAKE PRIMITIVE: O_denom = " << O_denom << std::endl;
    //div /= O_denom;
    div *= O_denom;
    assert (div % d == 0);
    div /= d;

    quat divide({NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), div}, O.alg);

    //std::pair<quat, NTL::ZZ> output(alpha*divide, T/(2*div)); //why tf times two?
    quat alpha_prim_conj = (alpha*divide).conjugate();
    std::pair<quat, NTL::ZZ> output(alpha_prim_conj, alpha_prim_conj.norm().first/(alpha_prim_conj.norm().second*N_I));
    return output;
}

std::pair<quat, NTL::ZZ> KLPT(quatlat const &J, factor_list const &fac_list, NTL::ZZ const &coprime) {
    //Look for equivalent ideal to J of norm dividing T
    //std::cout << "Starting KLPT..." << std::endl;
    //std::cout << "norm of ideal: " << J.norm().first/J.norm().second << std::endl;
    //std::cout << J.sage() << std::endl;
    quatlat O0 = J.left_order();
    J.reset_norm();
    NTL::ZZ N_J = J.norm().first/J.norm().second;
    //std::cout << "after J.left_order..." << std::endl;
    NTL::ZZ p = J.alg.p;
    NTL::ZZ q = J.alg.q;
    quat i({{NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, J.alg});
    quat j({{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, J.alg});

    NTL::ZZ T_full(1);
    NTL::ZZ Smallest_2_mod_3(1000000);
    NTL::ZZ Second_smallest_2_mod_3(1000000);

    for (const auto& tup : fac_list) {
        NTL::ZZ ell(std::get<0>(tup));
        if (NTL::GCD(coprime, ell) == 1) {
            T_full *= ell; // ell occurs once for each exponent....
            if (q == 3 && ell < Second_smallest_2_mod_3 && (ell % 3 == 2)) {
                if (ell < Smallest_2_mod_3) {
                    Second_smallest_2_mod_3 = Smallest_2_mod_3;
                    Smallest_2_mod_3 = ell;
                } else if (ell != Smallest_2_mod_3) {
                    Second_smallest_2_mod_3 = ell;
                }
            }
        }
    }
    if (q == 3) {
        assert (Smallest_2_mod_3 < 1000000);
        assert (Second_smallest_2_mod_3 < 1000000);
        assert (Second_smallest_2_mod_3 != Smallest_2_mod_3);
    }
    T_full *= q;

    for (size_t tries_outer = 0; tries_outer < 10000; tries_outer++) {
        //std::cout << "before equivPrime..." << std::endl;
        quatlat I = equivPrime(J, T_full); // Has to be coprime to T, q, so thats why T_full is included
        assert (I.norm().first % I.norm().second == 0);
        NTL::ZZ N_I = I.norm().first/I.norm().second;
        T_full *= N_I;
        //std::cout << "N_I: " << I.norm().first << "/" << I.norm().second << std::endl;
        //std::cout << I.sage() << std::endl;
        NTL::ZZ N_gamma = _choose_gamma_norm(fac_list, coprime, N_I, p, q, Smallest_2_mod_3, Second_smallest_2_mod_3); //TODO: Might be worth giving a list of acceptible N_gammas instead, more efficient maybe?
        NTL::ZZ norm_bound = q*N_gamma*p*N_I*N_I*N_I; //Constant 1000 is a bit arbitrary
        NTL::ZZ T(1);
        for (const auto& tup : fac_list) {
            NTL::ZZ ell(std::get<0>(tup));
            if (NTL::GCD(coprime, ell) == 1) {
                T *= ell; // ell occurs once for each exponent....
            }
            if (q != 3) {
                if (T > (2*q)*norm_bound) {
                    break;
                }
            } else {
                if (T > norm_bound && (N_I*T) % 3 == 1) {
                    break;
                } else if ((T > norm_bound/Second_smallest_2_mod_3) && (T % Second_smallest_2_mod_3 != 0)) {
                    break;
                }
                else if ((T > Second_smallest_2_mod_3*norm_bound) && (T % Second_smallest_2_mod_3 == 0)) {
                    break;
                }
            }
        }
        if (T % 2 == 0) {
            T = T/2;
        }
        if (T % q == 0) {
            T = T/q;
        }

        if ((q==3) && ((N_I*T) % 3 != 1)) {
            if (T % Second_smallest_2_mod_3 == 0) {
                T = T/Second_smallest_2_mod_3;
            } else {
                T = T*Second_smallest_2_mod_3;
            }
        }
        if (T < norm_bound) {
            std::cout << "T was too small!" << std::endl;
            std::cout << "KLPT with params: p: " << p << ", N_I: " << N_I << "; N_gamma: " << N_gamma << std::endl;
            continue;
        }
        assert (T % N_gamma == 0);
        NTL::ZZ N_mu = T/N_gamma;
        // Search in full order
        // N_mu *= SuborderIndex;
        bool past_jacobi = false;

        // q == 3 is annoying...
        if (q == 3) {
            if (N_I*N_gamma % 3 == 2) {
                std::cout << "STUPID SKIP" << std::endl;
                continue;
            }
            if (N_mu % 3 == 2) {
                std::cout << "STUPID SKIP" << std::endl;
                continue;
            }
        }

        for (size_t tries = 0; tries < 10; tries++) {
            quat gamma = RepresentInteger(J.alg, N_I*N_gamma);
            if (gamma.norm().first == 0) {
                break;
            }
            auto CD = IdealModConstraint(I, gamma);

            if ((CD.first == 0) || (CD.second == 0)) {
                continue;
            }

            // quat mu_0 = j*CD.first + j*i*CD.second;
            quat mu_0 = j*CD.first + i*j*CD.second;

            assert (mu_0.norm().second == 1);


            NTL::ZZ lam_sq = N_mu*NTL::InvMod(mu_0.norm().first % N_I, N_I) % N_I;

            if (!(NTL::Jacobi((lam_sq), N_I) == 1)) { // Can get stuck here if not enough solutions to Repr.Integer
                // std::cout << "Jacobi symbol was problem" << std::endl;
                continue;
            }

            past_jacobi = true;
            NTL::ZZ lam = NTL::SqrRootMod(lam_sq, N_I);
            auto SA_output = StrongApproximation(lam, CD.first, CD.second, N_I, N_mu, J.alg);
            if (SA_output) {
                quat beta = gamma * *SA_output;
                //std::cout << "Before make primitive?" << std::endl;
                std::pair<quat, NTL::ZZ> output = _makePrimitive(beta, O0, N_I);
                if (T % output.second == 0) {
                    assert (I.contains(output.first.conjugate()));
                    return output;
                }
            }
        }
        (void) past_jacobi;
    }

    std::ostringstream msg;
    msg << "KLPT failed for p = " << p << ", q = " << q;
    throw std::runtime_error(msg.str());
}


std::pair<quat, NTL::ZZ> KLPT_conj(quat *conj, quatlat const &J, factor_list const &fac_list, NTL::ZZ const &coprime) {
    //Look for equivalent ideal to J of norm dividing T
    //std::cout << "Starting KLPT..." << std::endl;
    //std::cout << "norm of ideal: " << J.norm().first/J.norm().second << std::endl;
    //std::cout << J.sage() << std::endl;
    quatlat O0 = J.left_order();
    J.reset_norm();
    NTL::ZZ N_J = J.norm().first/J.norm().second;
    //std::cout << "after J.left_order..." << std::endl;
    NTL::ZZ p = J.alg.p;
    NTL::ZZ q = J.alg.q;
    quat i({{NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, J.alg});
    quat j({{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, J.alg});

    NTL::ZZ T_full(1);

    NTL::ZZ Smallest_2_mod_3(1000000);
    NTL::ZZ Second_smallest_2_mod_3(1000000);

    for (const auto& tup : fac_list) {
        NTL::ZZ ell(std::get<0>(tup));
        if (NTL::GCD(coprime, ell) == 1) {
            T_full *= ell; // ell occurs once for each exponent....
            if (q == 3 && ell < Second_smallest_2_mod_3 && (ell % 3 == 2)) {
                if (ell < Smallest_2_mod_3) {
                    Second_smallest_2_mod_3 = Smallest_2_mod_3;
                    Smallest_2_mod_3 = ell;
                } else if (ell != Smallest_2_mod_3) {
                    Second_smallest_2_mod_3 = ell;
                }
            }
        }
    }
    if (q == 3) {
        assert (Smallest_2_mod_3 < 1000000);
        assert (Second_smallest_2_mod_3 < 1000000);
        assert (Second_smallest_2_mod_3 != Smallest_2_mod_3);
    }
    T_full *= q;
    T_full *= q;


    for (size_t tries_outer = 0; tries_outer < 10000; tries_outer++) {
        //std::cout << "before equivPrime..." << std::endl;
        quatlat I = equivPrime_conj(conj, J, T_full); // Has to be coprime to T, q, so thats why T_full is included
        assert (I.norm().first % I.norm().second == 0);
        NTL::ZZ N_I = I.norm().first/I.norm().second;
        T_full *= N_I;
        //std::cout << "N_I: " << I.norm().first << "/" << I.norm().second << std::endl;
        //std::cout << I.sage() << std::endl;
        NTL::ZZ N_gamma = _choose_gamma_norm(fac_list, coprime, N_I, p, q, Smallest_2_mod_3, Second_smallest_2_mod_3); //TODO: Might be worth giving a list of acceptible N_gammas instead, more efficient maybe?
        NTL::ZZ norm_bound = q*N_gamma*p*N_I*N_I*N_I; //Constant 1000 is a bit arbitrary
        NTL::ZZ T(1);
        for (const auto& tup : fac_list) {
            NTL::ZZ ell(std::get<0>(tup));
            if (NTL::GCD(coprime, ell) == 1) {
                T *= ell; // ell occurs once for each exponent....
            }
            if (q != 3) {
                if (T > (2*q)*norm_bound) {
                    break;
                }
            } else {
                if (T > norm_bound && (N_I*T) % 3 == 1) {
                    break;
                } else if ((T > norm_bound/Second_smallest_2_mod_3) && (T % Second_smallest_2_mod_3 != 0)) {
                    break;
                }
                else if ((T > Second_smallest_2_mod_3*norm_bound) && (T % Second_smallest_2_mod_3 == 0)) {
                    break;
                }
            }
        }
        //Leave extra room for 2 and q torsion for division
        if (T % 2 == 0) {
            T = T/2;
        }
        if (T % q == 0) {
            T = T/q;
        }
        if ((q==3) && ((N_I*T) % 3 != 1)) {
            if (T % Second_smallest_2_mod_3 == 0) {
                T = T/Second_smallest_2_mod_3;
            } else {
                T = T*Second_smallest_2_mod_3;
            }
        }
        if (T < norm_bound) {
            std::cout << "T was too small!" << std::endl;
            std::cout << "KLPT with params: p: " << p << ", N_I: " << N_I << "; N_gamma: " << N_gamma << std::endl;
            continue;
        }
        assert (T % N_gamma == 0);
        NTL::ZZ N_mu = T/N_gamma;
        //q == 3 is annoying...
        if (q == 3) {
            if (N_I*N_gamma % 3 == 2) {
                std::cout << "STUPID SKIP" << std::endl;
                continue;
            }
            if (N_mu % 3 == 2) {
                std::cout << "STUPID SKIP" << std::endl;
                continue;
            }
        }

        bool past_jacobi = false;

        for (size_t tries = 0; tries < 10; tries++) {
            quat gamma = RepresentInteger(J.alg, N_I*N_gamma);
            if (gamma.norm().first == 0) {
                break;
            }
            auto CD = IdealModConstraint(I, gamma);

            if ((CD.first == 0) || (CD.second == 0)) {
                continue;
            }

            //quat mu_0 = j*CD.first + j*i*CD.second;
            quat mu_0 = j*CD.first + i*j*CD.second;

            assert (mu_0.norm().second == 1);

            NTL::ZZ lam_sq = N_mu*NTL::InvMod(mu_0.norm().first % N_I, N_I) % N_I;

            if (!(NTL::Jacobi((lam_sq), N_I) == 1)) { //can get stuck here if not enough solutions to Repr.Integer
                // std::cout << "Jacobi symbol was problem" << std::endl;
                continue;
            }

            past_jacobi = true;
            NTL::ZZ lam = NTL::SqrRootMod(lam_sq, N_I);

            auto SA_output = StrongApproximation(lam, CD.first, CD.second, N_I, N_mu, J.alg);
            if (SA_output) {
                quat beta = gamma * *SA_output;
                std::pair<quat, NTL::ZZ> output = _makePrimitive(beta, O0, N_I);
                if (T % output.second == 0) {
                    assert (I.contains(output.first.conjugate()));
                    return output;
                }
            }
        }
        //std::cout << "KLPT with params failed: p: " << p << ", N_I: " << N_I << "; N_gamma: " << N_gamma << std::endl;
        //std::cout << "Made it past Jacobi? " << past_jacobi << std::endl;
        //std::cout << "Trying again with new equiv prime ideal..." << past_jacobi << std::endl;
        (void) past_jacobi;
    }

    std::ostringstream msg;
    msg << "KLPT failed for p = " << p << ", q = " << q;
    throw std::runtime_error(msg.str());
}
