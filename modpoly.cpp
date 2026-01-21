///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
//////////////// The code in this file implements the main algorithms
//////////////// in the accompanying paper, namely:
////////////////                    - SpecialSupersingularEvaluation (Algorithm 1)
////////////////                    - ModularEvaluationBigCharacteristic (Algorithm 2)
////////////////                    - SupersingularEvaluation (Algorithm 4)
////////////////                    - ModEvaluationBigLevel (Algorithm 5)
//////////////// and their Weber variants (as described in Section 3.5).
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <NTL/ZZ.h>
#include <NTL/vec_lzz_pE.h>
#include "modpoly.hpp"
#include "interpolation.hpp"
#include <unordered_set>
#include "fast_quaternions.hpp"
// #include "fast_ff.hpp"



#ifdef NDEBUG
// size_t const num_threads = 1 + std::thread::hardware_concurrency();
size_t const num_threads = 1;
#else
size_t const num_threads = 1; // For debugging
#endif

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///////// Get (small) primes to perform (Special)SupersingularEvaluation on
/////////     for classical and Weber variants
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

float log2(NTL::ZZ const x){
    int x_int = NTL::conv<float>(x);
    return logf(x_int)/logf(2);
}

std::vector<NTL::ZZ> _avail_qs(NTL::ZZ const &p, NTL::ZZ const &ell) {
    std::vector<NTL::ZZ> qs;
    if (p % 4 == 3) {
        qs.push_back(NTL::ZZ(1));
    } else if (!(ell == 3) && p % 3 == 2) {
        qs.push_back(NTL::ZZ(3));
    } else if (!(ell == 7) && (p % 7 == 3 || p % 7 == 5 || p % 7 == 6)) {
        qs.push_back(NTL::ZZ(7));
    } else if (!(ell == 11) && (p % 11 == 2 || p % 11 == 6 || p % 11 == 7 || p % 11 == 8 || p % 11 == 10)) {
        qs.push_back(NTL::ZZ(11));
    } else if (!(ell == 19) && (p % 19 == 2 || p % 19 == 3 || p % 19 == 8 || p % 19 == 10 || p % 19 == 12 || p % 19 == 13 || p % 19 == 14 || p % 19 == 15 || p % 19 == 18)) {
        qs.push_back(NTL::ZZ(19));
    }
    return qs;
}

long select_q_long(long p, long ell) {
    if (p % 4 == 3) {
        return 1;
    } else if (!(ell == 3) && p % 3 == 2) {
        return 3;
    } else if (!(ell == 7) && (p % 7 == 3 || p % 7 == 5 || p % 7 == 6)) {
        return 7;
    } else if (!(ell == 11) && (p % 11 == 2 || p % 11 == 6 || p % 11 == 7 || p % 11 == 8 || p % 11 == 10)) {
        return 11;
    } else if (!(ell == 19) && (p % 19 == 2 || p % 19 == 3 || p % 19 == 8 || p % 19 == 10 || p % 19 == 12 || p % 19 == 13 || p % 19 == 14 || p % 19 == 15 || p % 19 == 18)) {
        return 19;
    }
    return 0;
}


bool check_prime(long p, long q, long j) {
    Fp_integer p_mod;
    NTL::conv(p_mod, p);
    Fp_push push(p_mod);
    Fp_elem a, b;
    if (q == 1) {
        a = 1;
        b = 0;
    } else if (q == 3) {
        a = 0;
        b = 1;
    } else if (q == 7) {
        a = -35;
        b = 98;
    } else if (q == 11) {
        a = -1056;
        b = 13552;
    } else if (q == 19) {
        a = -152;
        b = 722;
    } else {
        throw;
    }

    auto a3 = Fp_elem(4)*a*a*a;
    auto b2 = Fp_elem(27)*b*b;
    auto disc = a3 + b2;
    //check E0 has good reduction
    if (disc == 0) {
        return false;
    }
    auto j_0 = Fp_elem(1728) * a3/(a3 + b2);
    Fp_elem j_1(j);

    // This is actually an easy case in theory, but requires a separate function
    // It's very rare anyway
    if (j_0 == j_1) {
        return false;
    }
    return true;
}

NTL::ZZ GetBoundClassicalBigLevel(NTL::ZZ const j, NTL::ZZ const l){
    // TODO: this can be done better but without float theres a big rounding error
    NTL::ZZ B;

    double lf = NTL::conv<double>(l);
    double log2l = lf*NTL::NumBits(l);
    double log2j = NTL::NumBits(j);
    double log2lp2 = log2(lf+2);

    long b = round(6*log2l + 18*lf + lf * log2j + log2lp2);

    std::cout << "The target bound is " << b << " bits" << std::endl;
    NTL::power(B, NTL::ZZ(2), b);

    return B;
}

NTL::ZZ GetBoundWeberBigLevel(NTL::ZZ const w, NTL::ZZ const l){
    // TODO: this can be done better but without float theres a big rounding error
    NTL::ZZ B;

    double lf = NTL::conv<double>(l);
    double log2l = lf*NTL::NumBits(l);
    double log2w = NTL::NumBits(w);
    double log2lp2 = log2(lf+2);

    long b = round(log2l/12 + lf/5 + lf * log2w + log2lp2);

    std::cout << "The target bound is " << b << " bits" << std::endl;
    NTL::power(B, NTL::ZZ(2), b);

    return B;

}


void GetPrimes(std::vector<NTL::ZZ> &Pl, NTL::ZZ const j_int, NTL::ZZ const B, long l){

    NTL::ZZ PP(1);

    long int j_long;
    NTL::conv(j_long, j_int);

    std::mutex mtx;

    //store all primes up to given bound, required for sieving
    long int bound_sieving_primes = 33554432; //2^25, so prime sieve can reach primes up to 2^50

    //length of sieving intervals
    long int len_sieve_interval = 65536; //2^16

    //read sieving primes from file
    std::vector<long int> sieve_primes;
    if (!init_sieve("primes_upto_2pow25.csv", sieve_primes))
            std::cerr << "Failed to load primes.\n";

    int num_primes = sieve_primes.size();

    auto const fun = [&](size_t tidx) {

        int i = tidx;

        while (PP <= B && i < num_primes) {

            long int p = sieve_primes[i];

            i += num_threads;

            long q = select_q_long(p, l);

            if (p < 200) { // seems buggy with really small p
                continue;
            }

            if (q == 0) {
                continue;
            }
            Fp_integer p_mod;
            NTL::conv(p_mod, p);
            if (sutherland(j_long, p_mod)) {
                if (check_prime(p, q, j_long)) {
                    std::lock_guard lock(mtx);
                    Pl.push_back(NTL::ZZ(p));
                    PP *= NTL::ZZ(p);
                    std::cerr << "Found a good prime: " << p << ", q = " << q << ", now have " << NTL::NumBits(PP) << " bits" << std::endl;
                }
            }
        }
    };


    auto const fun2 = [&] (size_t tidx) {

        long int L = bound_sieving_primes + tidx * len_sieve_interval;
        while (PP <= B) {
            auto primes = sieve_interval(L, len_sieve_interval, sieve_primes);

            for (auto &p: primes){

                long q = select_q_long(p, l);
                if (q == 0) {
                    continue;
                }

                Fp_integer p_mod;
                NTL::conv(p_mod, p);
                if (PP <= B && sutherland(j_long, p_mod)) {
                    if (check_prime(p, q, j_long)) {
                        std::lock_guard lock(mtx);
                        Pl.push_back(NTL::ZZ(p));
                        PP *= NTL::ZZ(p);
                        std::cerr << "Found a good prime: " << p << ", q = " << q << ", now have " << NTL::NumBits(PP) << " bits" << std::endl;
                    }
                }
            }

            L += num_threads * len_sieve_interval;
        }
    };

    //go through precomputed primes
    std::vector<std::thread> ts;
    for (size_t i = 0; i < num_threads; ++i)
        ts.emplace_back(fun, i);
    for (auto &t: ts)
        t.join();

    //sieve for more primes
    if (PP <= B){
        std::vector<std::thread> ts2;
        for (size_t i = 0; i < num_threads; ++i)
            ts2.emplace_back(fun2, i);
        for (auto &t: ts2)
            t.join();
    }

    std::sort(Pl.begin(), Pl.end());
}

void GetPrimesWeber(std::vector<NTL::ZZ> &Pl, NTL::ZZ const w_int, NTL::ZZ const B, long l){

    NTL::ZZ PP(1);
    NTL::ZZ j_int;

    std::mutex mtx;

    //store all primes up to given bound, required for sieving
    long int bound_sieving_primes = 33554432; //2^25, so prime sieve can reach primes up to 2^50

    //length of sieving intervals
    long int len_sieve_interval = 65536; //2^16

    //read sieving primes from file
    std::vector<long int> sieve_primes;
    if (!init_sieve("primes_upto_2pow25.csv", sieve_primes))
            std::cerr << "Failed to load primes.\n";

    int num_primes = sieve_primes.size();

    auto const fun = [&](size_t tidx) {

        int i = tidx;

        while (PP <= B && i < num_primes) {

            long int p = sieve_primes[i];

            i += num_threads;

            long q = select_q_long(p, l);

            if (p < 200) { // seems buggy with really small p
                continue;
            }

            if (q == 0) {
                continue;
            }

            NTL::ZZ p_ZZ(p);
            long int j_long;
            {
                NTL::ZZ_pPush push(p_ZZ);
                NTL::ZZ_p w_mod;
                NTL::conv(w_mod, w_int);
                if (w_mod == 0) {
                    continue;
                }
                NTL::ZZ_p j_mod = w_to_j_prime(w_mod);
                NTL::conv(j_long, j_mod);
            }

            Fp_integer p_mod;
            NTL::conv(p_mod, p);
            if (sutherland(j_long, p_mod)) {
                if (check_prime(p, q, j_long)) {
                    std::lock_guard lock(mtx);
                    Pl.push_back(p_ZZ);
                    PP *= p_ZZ;
                    std::cerr << "Found a good prime: " << p << ", q = " << q << ", now have " << NTL::NumBits(PP) << " bits" << std::endl;
                }
            }
        }
    };


    auto const fun2 = [&] (size_t tidx) {

        long int L = bound_sieving_primes + tidx * len_sieve_interval;
        while (PP <= B) {
            auto primes = sieve_interval(L, len_sieve_interval, sieve_primes);

            for (auto &p: primes){

                long q = select_q_long(p, l);
                if (q == 0) {
                    continue;
                }
                //Remove once bugs are fixed
                if (q == 7 || q == 19) {
                    continue;
                }

                NTL::ZZ p_ZZ(p);
                long int j_long;
                {
                    NTL::ZZ_pPush push(p_ZZ);
                    NTL::ZZ_p w_mod;
                    NTL::conv(w_mod, w_int);
                    NTL::ZZ_p j_mod = w_to_j_prime(w_mod);
                    NTL::conv(j_long, j_mod);
                }

                Fp_integer p_mod;
                NTL::conv(p_mod, p);
                if (PP <= B && sutherland(j_long, p_mod)) {
                    if (check_prime(p, q, j_long)) {
                        std::lock_guard lock(mtx);
                        Pl.push_back(p_ZZ);
                        PP *= p_ZZ;
                        std::cerr << "Found a good prime: " << p << ", q = " << q << ", now have " << NTL::NumBits(PP) << " bits" << std::endl;
                    }
                }
            }

            L += num_threads * len_sieve_interval;
        }
    };

    //go through precomputed primes
    std::vector<std::thread> ts;
    for (size_t i = 0; i < num_threads; ++i)
        ts.emplace_back(fun, i);
    for (auto &t: ts)
        t.join();

    //sieve for more primes
    if (PP <= B){
        std::vector<std::thread> ts2;
        for (size_t i = 0; i < num_threads; ++i)
            ts2.emplace_back(fun2, i);
        for (auto &t: ts2)
            t.join();
    }

    std::sort(Pl.begin(), Pl.end());
}

NTL::ZZ GetBoundClassicalBigChar(NTL::ZZ const p, NTL::ZZ const l){
    // TODO: this can be done better but without float theres a big rounding error
    NTL::ZZ B;

    double lf = NTL::conv<double>(l);
    double log2l = lf*NTL::NumBits(l);
    double log2p = NTL::NumBits(p);
    double log2lp2 = log2(lf+2);

    int b = round(6*log2l + 18*lf + log2p + log2lp2);

    std::cout << "The target bound is " << b << " bits" << std::endl;
    NTL::power(B, NTL::ZZ(2), b);

    return B;
}

NTL::ZZ GetBoundWeberBigChar( NTL::ZZ const p, NTL::ZZ const l){
    // TODO: this can be done better but without float theres a big rounding error
    NTL::ZZ B;

    double lf = NTL::conv<double>(l);
    double log2l = lf*NTL::NumBits(l);
    double log2p = NTL::NumBits(p);
    double log2lp2 = log2(lf+2);
    std::cerr << lf << " " << log2l << " " << log2p << " " << log2lp2 << std::endl;

    long b = round(log2l/12 + lf/5 + log2p + log2lp2);

    std::cout << "The target bound is " << b << " bits" << std::endl;
    NTL::power(B, NTL::ZZ(2), b);

    return B;

}


void GetPrimesBigCharWeber(std::vector<Integer> &vec, const Integer B, const Integer ell, const Integer p) {

    Integer prime = NTL::NextPrime(ell/6 + 600);
    Integer prod = Integer(1);

    while (prod < B ) {

        // TODO we could use all primes
        if (prime % 4 ==3 && prime != ell && prime != p){
            prod *= prime;
            vec.push_back(prime);
        }
        prime = NTL::NextPrime(prime+1);

    }
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////// Utility functions
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


bool isPrincipal_Compute(quat *gamma, const quatlat &J) {

    quatlatenum Enumerator(J);


    auto &c = Enumerator.gram[0][0].get_data();
    Integer norm;
    gmp2ntl(norm, c);
    bool isPrincipal = (norm * J.norm().second == 2 * J.norm().first);


    // principal, we compute the smallest element which is the generator
    if (isPrincipal) {
        NTL::Vec<NTL::ZZ> row;
        row.SetLength(4);
        for (unsigned j = 0; j < 4; ++j) {
            NTL::ZZ c;
            gmp2ntl(c, Enumerator.U[0][j].get_data());
            row += c * J.basis[j];
        }
        (*gamma)[0] = row[0];
        (*gamma)[1] = row[1];
        (*gamma)[2] = row[2];
        (*gamma)[3] = row[3];
        (*gamma)[4] = J.denom;
        assert(J.contains(*gamma));
        assert(gamma->norm().first * J.norm().second == gamma->norm().second * J.norm().first);
    }

    return isPrincipal;
}

FpE_elem conjugate(FpE_elem const &a) {
    FpX_elem conj_a = NTL::rep(a);
    conj_a[1] = -conj_a[1];
    FpE_elem conj_a_done;
    NTL::conv(conj_a_done, conj_a);
    return conj_a_done;
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////// Implementing the BIG LEVEL variant
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


std::vector<FpE_elem> SupersingularEvaluation(Fp_integer p, FpE_elem const &j, long l, bool useKLPT) {

    ///////////////////////////////////////////////////////////////////////////////////
    // Implements SupersingularEvaluation (Algorithm 4) in the accompanying paper
    ///////////////////////////////////////////////////////////////////////////////////

    std::cout << "\n\nRunning SupersingularEvaluationWeber with p = " << p << ", j = " << j << " and l = " << l << std::endl;

    auto start_full = std::chrono::steady_clock::now();

    FpE_elem j_E1 = j;
    std::cout << "Using curve j(E1) = " << j_E1 << std::endl;
    ec E1 = ec::from_j(j_E1);

    NTL::ZZ p_ZZ(p);
    NTL::ZZ l_ZZ(l);
    std::vector<FpE_elem> j_invariants;
    std::map<NTL::ZZ,std::pair<ecp,ecp>> bases;
    auto qs = _avail_qs(p_ZZ, l_ZZ);
    assert(qs.size() >= 1);
    auto q = qs.at(0);

    NTL::ZZ q_ZZ(q);
    std::cout << "Using quaternion algebra (-" << q << ", -" << p << ")" << std::endl;
    quatalg Bp {p_ZZ, q_ZZ};

    auto start_computing_fieldexts = std::chrono::steady_clock::now();
    // torsion/factor basis stuff
    NTL::ZZ T;

    assert (useKLPT); (void) useKLPT; // change later if we want to use the other one
    NTL::ZZ KLPT_const;
    NTL::conv(KLPT_const, "100000"); //1000 in for N_mu, 100 extra because it overshoots
    T = q_ZZ*p_ZZ*p_ZZ*p_ZZ*p_ZZ*KLPT_const;

    factor_list tors_list;
    tors_list = choose_torsion(p_ZZ, T, NTL::ZZ(6)*5);

    NTL::ZZ S(1);
    unsigned k_bound = 20; //Some minimum
    for (const auto& tup : tors_list) {
        if (std::get<2>(tup) > (ssize_t) k_bound) {
            k_bound = std::get<2>(tup);
        }
    }

    std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO: eventually we should generate these field extensions on the fly and cache
    //      them, similar to how we now do it for torsion bases.
    //      At the moment this fails because references to the Fp2k object are stored
    //      in all kinds of other objects and those references are invalidated when
    //      the std::map is modified. possible solution: use std::shared_ptr for Fp2k
    //      references, just like we do for ec references.

    std::map<unsigned,Fp2k> Fexts;
    {
        /*NTL::ZZ_p::init(p);
        NTL::ZZ_pX f;
        SetCoeff(f, 2);
        f[0] = NTL::ZZ_p(3); //is this really needed?
        NTL::ZZ_pE::init(f);*/
        for (unsigned k = 1; k <= k_bound; ++k) {
            std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts.emplace(std::make_pair(k, Fp2k {k}));
        }
        std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }

    std::chrono::duration<long, std::milli> duration_compute_fieldexts = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_computing_fieldexts);
    std::cout << ">>> Computing field extensions took: " << duration_compute_fieldexts.count() << " milliseconds" << std::endl;

    //std::cout << "Starting compute Endring..." << std::endl;
    auto start_comp_endring = std::chrono::steady_clock::now();
    auto start = starting_curve(Bp, Fexts, false); //This doesnt care about being on the surface
    ec E0 = start.first; // This is sometimes y^2 = x^3 - x, messes things up...
    auto j_E0 = E0.j_invariant();
    assert (j_E0 != j_E1);
    quatlat O0 = start.second;


    assert(NTL::NumBytes(l) <= (ssize_t) sizeof(long));
    auto O1 = compute_endring(E1, Bp, Fexts, l);
    auto I_conn0 = connecting_ideal(O0, O1);
    //std::cout << "O1 found: " << std::endl;
    //std::cout << O1.sage() << std::endl;

    std::chrono::duration<long, std::milli> duration_compute_starting_endring = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_comp_endring);
    std::cout << ">>> Computing endring took: " << duration_compute_starting_endring.count() << " milliseconds" << std::endl;

    std::chrono::duration<long, std::milli> duration_quaternionstuff = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());
    std::chrono::duration<long, std::milli> duration_isogenies = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());

    //std::cout << "Entering left ideals..." << std::endl;
    long lp1 = l + 1;
    long count = 0;

    std::unordered_map<quatlat, FpE_elem> ideal_to_invariant;
    std::cerr << "Processing ideal " << 0 << "/" << lp1 << std::flush;
    O1.left_ideals_of_norm(l_ZZ, [&](quatlat const &I) {
        //std::cout << "\n~~~~ New left ideal (" << count++ << " / " << (time(NULL)-t0) << "s) ~~~~~" << std::endl;
        std::cerr << "\r\x1b[KProcessing ideal " << count++ << " out of " << lp1 << std::flush;
        Fp_push push((Fp_integer(p)));
        FpE_elem j_iso;

        auto hit = ideal_to_invariant.find(I.conjugated_by_j());
        if (hit != ideal_to_invariant.end()) {
            j_iso = conjugate(hit->second);
        } else if (isPrincipal(I)) {
            j_iso = j_E1;
            ideal_to_invariant.emplace(std::make_pair(I, j_iso));
        } else {
            auto I_conn = I_conn0*I;
            I_conn.reset_norm();

            if (isPrincipal(I_conn)) {
                //std::cout << "???????????" << std::endl;
                //std::cout << "The ideal I_conn was principal!" << std::endl;
                j_iso = j_E0;
                ideal_to_invariant.emplace(std::make_pair(I, j_iso));
            } else {
                //std::cout << "============" << std::endl;
                // Find a smooth normed connecting ideal
                std::optional<FpE_elem> j;
                // Have to be careful, if O_R(I_conn) == O_0 this currently fails.
                auto start_klpt = std::chrono::steady_clock::now();
                auto gens = KLPT(I_conn, tors_list, NTL::ZZ(6));
                //std::cout << "KLPT done:" << std::endl;
                //std::cout << "beta: " << gens.first << std::endl;
                //std::cout << "Norm of ideal: " << gens.second << std::endl;
                duration_quaternionstuff = std::chrono::duration_cast<std::chrono::milliseconds>(duration_quaternionstuff + std::chrono::steady_clock::now() - start_klpt);

                auto start_isogeny = std::chrono::steady_clock::now(); (void) start_isogeny;
                auto phi_I = idealToIsogeny(gens.first, gens.second, E0, Fexts, bases);
                ec E_I = phi_I.get_codomain();
                j = E_I.j_invariant();
                j_iso = (*j);
                ideal_to_invariant.emplace(std::make_pair(I, j_iso));
                duration_isogenies = std::chrono::duration_cast<std::chrono::milliseconds>(duration_isogenies + std::chrono::steady_clock::now() - start_isogeny);
            }
        }
        std::cout << j_iso << std::endl;
        j_invariants.push_back(j_iso);
    });
    std::cerr << "\r\x1b[KProcessing ideal " << lp1 << " out of " << lp1 << " >>>> DONE!" << std::flush;


    std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~ Summary ~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_full);
    std::cout << ">>> We have used in total " << duration.count() << " seconds" << std::endl;
    std::cout << ">>> On computing field extensions " << duration_compute_fieldexts.count() << " milliseconds" << std::endl;
    std::cout << ">>> On computing endomorphism ring " << duration_compute_starting_endring.count() << " milliseconds" << std::endl;
    std::cout << ">>> On quaternion computation KLPT: " << duration_quaternionstuff.count() << " milliseconds" << std::endl;
    std::cout << ">>> On computing ideal2isogeny: " << duration_isogenies.count() << " milliseconds" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    return j_invariants;
}

FpX_big_elem ModEvalBigLevel(NTL::ZZ p, NTL::ZZ_pE const j, long l)
{
    ///////////////////////////////////////////////////////////////////////////////////
    // Implements ModularEvaluationBigLevel (Algorithm 5) in the accompanying paper
    ///////////////////////////////////////////////////////////////////////////////////

    // Initialise mod poly F
    NTL::ZZ_pX F;

    // // Get the bound for prime search
    NTL::ZZ B;
    NTL::ZZ l_ZZ(l);



    // Number of coeffs
    int l_int = NTL::conv<int>(l_ZZ);
    int Ncoeffs;

    Ncoeffs = l_int + 2;

    NTL::ZZ j_ZZ = NTL::conv<NTL::ZZ>(NTL::ConstTerm(NTL::conv<NTL::ZZ_pX>(j)));
    FpE_elem j_modq;

    B = GetBoundClassicalBigLevel(j_ZZ, l_ZZ);

    // Computing the set of primes
    std::vector<NTL::ZZ> Pl;
    std::cout << "Finding primes..." << std::endl;
    GetPrimes(Pl, j_ZZ, B, l);
    std::cout << "Done!" << std::endl;

    int Nprimes = Pl.size();

    // Initialise crt structure
    crt_info crt;
    std::cout << "Initialising CRT coeffs..." << std::endl;
    crt_init(crt, Pl, Nprimes, Ncoeffs, p);
    std::cout << "Done!" << std::endl;

    // Compute the F mod q and update crt coeffs
    std::cout << "We are working with " << Nprimes << " primes." << std::endl;

    // Note: we process the primes in reverse order since large primes will probably take longer
    std::atomic<ssize_t> next_idx = Pl.size();
    std::mutex mtx;

    auto const fun = [&]() {

        while (true) {

            ssize_t qidx = --next_idx;
            if (qidx < 0)
                break;

            NTL::ZZ const &q_ZZ = Pl[qidx];
            Fp_integer q;
            NTL::conv(q, q_ZZ);

            std::cerr << "Starting with a new prime: " << q << std::endl;

            // In this function we set this to be in ZZX to not work with two moduli in a function
            std::vector<NTL::ZZ> Fq_coeffs(crt.k);

            // Not confusing at all that p is named q, q is named qq... Makes me qqq
            Fp_elem::init(q);
            FpX_elem f;
            SetCoeff(f, 2);
            auto qqs = _avail_qs(q_ZZ, l_ZZ);
            Fp_elem qq;
            NTL::conv(qq, qqs.at(0));
            f[0] = qq;
            FpE_elem::init(f);
            NTL::conv(j_modq, j_ZZ);
            // Compute Fq
            std::cerr << "Entering Supersingular evaluation..." << std::endl;
            bool KLPT = true;

            auto j_invariants = SupersingularEvaluation(q, j_modq, l, KLPT); //Either j or weber invariants
            std::cerr << "Done!" << std::endl;

            // std::cout << "\n\n~~~~~~ recovered poly: ~~~~~~" << std::endl;
            // std::cout << Phi_l << std::endl;
            // std::cout << "~~~~~~~~~~~~~~~~~~~\n\n" << std::endl;

            std::cerr << "Entering interpolation..." << std::endl;

            FpEX_elem Fq;
            FastInterpolateFromRoots(Fq, j_invariants);
            //for (auto j : j_invariants) {
            //    std::cout << j << "\n";
            //}

            std::cout << "~~~~~~~~~~~~~~~~~~~ recovered poly: ~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << Fq << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" << std::endl;


            std::cerr << "Done!" << std::endl;

            // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
            for(int i = 0; i <= deg(Fq); i++){
                assert(is_Fp(NTL::coeff(Fq,i)));
                Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::ConstTerm(NTL::rep(NTL::coeff(Fq, i))));
            }

            //Update CRT sums
            {
                std::lock_guard lock(mtx);
                std::cerr << "Updating CRT coeffs for q=" << q << "... " << std::flush;
                crt_update(crt, qidx, Fq_coeffs, crt.k);
                std::cerr << "Done!\n\n\n" << std::endl;
            }
        }
    };

    {
        std::vector<std::thread> ts;
        for (size_t i = 0; i < num_threads; ++i)
            ts.emplace_back(fun);
        for (auto &t: ts)
            t.join();
    }

    std::cout << "Finalising CRT coeffs..." << std::endl;
    crt_finalise(crt);
    for(int i = 0; i < crt.k; i++){
        NTL::SetCoeff(F, i, NTL::conv<NTL::ZZ_p>(crt.Cdata[i]));
    }
    std::cout << "Done!" << std::endl;

    return F;
}



std::vector<FpE_elem> SupersingularEvaluationWeber(Fp_integer p, FpE_elem const w, long l) {
    //////////////////////////////////////////////////////////////////////////////
    /// Weber variant of SupersingularEvaluation
    ///     Implementation of the ideas in Section 3.5 for the BIG LEVEL variant
    //////////////////////////////////////////////////////////////////////////////

    // std::cout << "\n\nRunning SupersingularEvaluationWeber with p = " << p << ", w = " << w << " and l = " << l << std::endl;

    auto start_full = std::chrono::steady_clock::now();

    FpE_elem j_E1 = w_to_j(w);
    // std::cout << "Using curve j(E1) = " << j_E1 << std::endl;
    ec E1 = ec::from_j(j_E1);

    NTL::ZZ p_ZZ(p);
    NTL::ZZ l_ZZ(l);
    auto qs = _avail_qs(p_ZZ, l_ZZ);
    assert(qs.size() >= 1);
    auto q = qs.at(0);

    NTL::ZZ q_ZZ(q);
    // std::cout << "Using quaternion algebra (-" << q << ", -" << p << ")" << std::endl;
    quatalg Bp {p_ZZ, q_ZZ};

    auto start_computing_fieldexts = std::chrono::steady_clock::now();
    // torsion/factor basis stuff
    NTL::ZZ T;

    NTL::ZZ KLPT_const;
    NTL::conv(KLPT_const, "100000"); //1000 in for N_mu, 100 as KLPT bound, 1000 for log factors
    T = q_ZZ*p_ZZ*p_ZZ*p_ZZ*p_ZZ*KLPT_const;

    factor_list tors_list;
    tors_list = choose_torsion(p_ZZ, T, NTL::ZZ(6)*5);
    // Add 32 and 3 manually
    NTL::ZZ num_3or9;
    NTL::ZZ cofac_3;
    long val3;
    if (q == 3) {
        num_3or9 = 9;
        cofac_3 = 3;
        val3 = 2;
    } else {
        num_3or9 = 3;
        cofac_3 = 1;
        val3 = 1;
    }

    NTL::ZZ ord_32(p + 1);
    int k_2 = 1;
    while (ord_32 % 32 != 0){
        k_2 += 1;
        ord_32 = NTL::power(p_ZZ, long(k_2)) - NTL::power_long(long(-1), long(k_2 % 2));
    }
    ell_tuple tup_32;
    std::get<0>(tup_32) = 2;
    std::get<1>(tup_32) = 5;
    std::get<2>(tup_32) = k_2;
    std::get<3>(tup_32) = NTL::RR(0);
    tors_list.push_back(tup_32);

    ell_tuple tup_16;
    std::get<0>(tup_16) = 2;
    std::get<1>(tup_16) = 4;
    std::get<2>(tup_16) = k_2;
    std::get<3>(tup_16) = NTL::RR(0);
    tors_list.push_back(tup_16);

    ell_tuple tup_2;
    std::get<0>(tup_2) = 2;
    std::get<1>(tup_2) = 1;
    std::get<2>(tup_2) = k_2;
    std::get<3>(tup_2) = NTL::RR(0);
    tors_list.push_back(tup_2);

    NTL::ZZ ord_3(p + 1);
    int k_3 = 1;
    while (ord_3 % num_3or9 != 0){
        k_3 += 1;
        ord_3 = NTL::power(p_ZZ, long(k_3)) - NTL::power_long(long(-1), long(k_3 % 2));
    }
    ell_tuple tup_3;
    std::get<0>(tup_3) = 3;
    std::get<1>(tup_3) = val3;
    std::get<2>(tup_3) = k_3;
    std::get<3>(tup_3) = NTL::RR(0);
    tors_list.push_back(tup_3);

    if (val3 > 1) {
        ell_tuple tup_3_really;
        std::get<0>(tup_3_really) = 3;
        std::get<1>(tup_3_really) = 1;
        std::get<2>(tup_3_really) = k_3;
        std::get<3>(tup_3_really) = NTL::RR(0);
        tors_list.push_back(tup_3_really);
    }


    NTL::ZZ S(1);
    unsigned k_bound = 20; //Some minimum
    for (const auto& tup : tors_list) {
        if (std::get<2>(tup) > (ssize_t) k_bound) {
            k_bound = std::get<2>(tup);
        }
    }

    std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
    std::map<unsigned,Fp2k> Fexts;
    {
        for (unsigned k = 1; k <= k_bound; ++k) {
            std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts.emplace(std::make_pair(k, Fp2k {k}));
        }
        std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }

    std::chrono::duration<long, std::milli> duration_compute_fieldexts = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_computing_fieldexts);
    std::cout << ">>> Computing field extensions took: " << duration_compute_fieldexts.count() << " milliseconds" << std::endl;


    weber_enum_poly_precomp precomp = SetWeberPrecomp();

    //std::cout << "Starting compute Endring..." << std::endl;
    auto start_comp_endring = std::chrono::steady_clock::now();
    //auto [E0,O0] = starting_curve(Bp, is_surface(E1)); This syntax does not work for me?
    auto start = starting_curve(Bp, Fexts, false); //This doesnt care about being on the surface
    ec E0 = start.first; // This is sometimes y^2 = x^3 - x, messes things up...
    auto j_E0 = E0.j_invariant();
    assert (j_E0 != j_E1);
    quatlat O0 = start.second;
    //std::cout << "O0_orig" << std::endl;
    //std::cout << O0.sage() << std::endl;


    assert(NTL::NumBytes(l) <= (ssize_t) sizeof(long));
    auto O1 = compute_endring(E1, Bp, Fexts, l);
    //std::cout << "O1 found: " << std::endl;
    //std::cout << O1.sage() << std::endl;

    std::chrono::duration<long, std::milli> duration_compute_starting_endring = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_comp_endring);
    std::cout << ">>> Computing endring took: " << duration_compute_starting_endring.count() << " milliseconds" << std::endl;

    // List of weber invariants

    auto start_comp_weber = std::chrono::steady_clock::now();
    std::vector<FpE_elem> w_invariants;

    FpE_elem wd;
    wd = w;
    std::vector<std::vector<ecp>> levelstruc;
    std::vector<std::vector<ecp>> levelstruc_above;

    // Computing a basis for weber enumeration
    unsigned k = torsionToFieldDegree(num_3or9);
    auto jt = Fexts.find(k);
    assert(jt != Fexts.end());
    assert(jt->first == jt->second.k);
    auto bas3 = E0.torsionBasis(jt->second, int(3), val3);

    k = torsionToFieldDegree(Integer(32));
    jt = Fexts.find(k);
    assert(jt != Fexts.end());
    assert(jt->first == jt->second.k);
    auto bas2 = E0.torsionBasis(jt->second, int(2), 5);

    weber_bas web_above = {bas2.first, bas2.second, bas3.first, bas3.second};
    weber_bas web = {2*bas2.first, 2*bas2.second, cofac_3*bas3.first, cofac_3*bas3.second};


    size_t num_elems;
    NTL::conv(num_elems, l+1);
    w_invariants.reserve(num_elems);

    std::map<NTL::ZZ,std::pair<ecp,ecp>> bases;

    size_t count = 0;
    // size_t t0 = time(NULL);

    auto I_conn0 = connecting_ideal(O0, O1);
    quat conj(O0.alg);
    //std::cout << "Connecting ideal" << std::endl;
    //std::cout << I_conn0.sage() << std::endl;
    //std::cout << "O0" << std::endl;
    //std::cout << O0.sage() << std::endl;
    auto gens = KLPT_conj(&conj, I_conn0, tors_list, NTL::ZZ(6));
    //std::cout << "KLPT done!" << std::endl;
    //std::cout << gens.first << std::endl;
    //std::cout << gens.second << std::endl;
    //needed for weber later
    auto gen_01 = gens.first;
    auto trans = conj*gens.first;
    auto trans_inv = trans;
    trans_inv.invert();

    auto I_01 = create_from_generator(gens.first, gens.second, O0);
    // I_01.reset_norm();
    auto phi_0 = idealToIsogeny(gens.first, gens.second, E0, Fexts, bases);
    //std::cout << "Id2Iso finished" << std::endl;
    assert (phi_0.get_codomain().j_invariant() == j_E1);
    // auto webdata = EnumerateAllWeberFast(web, Fexts, &precomp);
    auto webdata = EnumerateAllWeberCoeff();
    //std::cout << "Getting compatible level structure on E0..." << std::endl;
    bool found_weber = 0;

    for (auto const &webdata_i : webdata){
        auto coeff = webdata_i;
        auto c2 = coeff[0][0];
        auto c31 = coeff[1][0];
        auto c32 = coeff[1][1];
        ecp P3 = c31.first * web_above.P3 + c31.second * web_above.Q3;
        ecp Q3 = c32.first * web_above.P3 + c32.second * web_above.Q3;
        auto c161 = coeff[2][0];
        auto c162 = coeff[2][1];
        ecp P2 = (c2.first*8)*web_above.P16 + (c2.second*8)*web_above.Q16;
        assert((4*P2).is_identity());
        assert(!(2*P2).is_identity());
        assert((num_3or9*P3).is_identity());
        assert(!(cofac_3*P3).is_identity());
        std::vector<ecp> tors3or9 = {P3, Q3, P3+Q3, P3+2*Q3};
        std::vector<ecp> tors3 = {cofac_3*tors3or9[0], cofac_3*tors3or9[1], cofac_3*tors3or9[2], cofac_3*tors3or9[3]};
        std::vector<ecp> tors32 = {c161.first*web_above.P16+c161.second*web_above.Q16,c162.first*web_above.P16+c162.second*web_above.Q16};
        std::vector<ecp> tors16 = {2*tors32[0], 2*tors32[1]};
        levelstruc_above = {{P2}, tors3or9, tors32};
        levelstruc = {{2*P2}, tors3, tors16};
        Fp2 wt;
        auto check = GetWeberOfImage_chain(&wt, phi_0, levelstruc, Fexts, true);
        if((check) && (wt == wd)){
            found_weber = true;
            break;
        }
    }
    assert (found_weber); (void) found_weber;

    std::chrono::duration<long, std::milli> duration_compute_weber = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_comp_weber);
    std::cout << ">>> Computing first weber structure took: " << duration_compute_weber.count() << " milliseconds" << std::endl;

    //std::cout << "Done with weber setup stuff!" << std::endl;

    std::chrono::duration<long, std::milli> duration_quaternionstuff = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());
    std::chrono::duration<long, std::milli> duration_isogenies = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());

    //std::cout << "Entering left ideals..." << std::endl;
    long lp1 = l + 1;
    std::unordered_map<quatlat, FpE_elem> ideal_to_invariant;
    std::cerr << "Processing ideal " << 0 << "/" << lp1 << std::flush;
    O1.left_ideals_of_norm(l_ZZ, [&](quatlat const &I) {
        //std::cout << "\n~~~~ New left ideal (" << count++ << " / " << (time(NULL)-t0) << "s) ~~~~~" << std::endl;
        std::cerr << "\r\x1b[KProcessing ideal " << count++ << " out of " << lp1 << std::flush;
        Fp_push push((Fp_integer(p)));
        FpE_elem w_iso;

        auto hit = ideal_to_invariant.find(I.conjugated_by_j());
        if (hit != ideal_to_invariant.end()) {
            w_iso = conjugate(hit->second);
        } else if (isPrincipal(I)) {
            //std::cout << "!!!!!!!!!!!!!" << std::endl;
            //std::cout << "The ideal I was principal!" << std::endl;
            auto start_weber = std::chrono::steady_clock::now();
            auto I_trans = trans_inv*I*trans;
            auto rho_ideal = I_01.norm().first * I_trans;

            quat rho(O0.alg);
            bool check_principal = isPrincipal_Compute(&rho, rho_ideal);
            //std::cout << "Done!" << std::endl;
            assert (check_principal); (void) check_principal;

            quat rho_int = rho*(2*q);
            //std::cout << rho_int << std::endl;

            std::vector<std::vector<ecp>> new_levelstruc;
            std::vector<ecp> new_tors_2 = {evalEndo(rho_int, levelstruc_above[0][0], NTL::ZZ(4))};
            std::vector<ecp> new_tors_3 = {evalEndo(rho_int, levelstruc_above[1][0], num_3or9), evalEndo(rho_int, levelstruc_above[1][1], num_3or9), evalEndo(rho_int, levelstruc_above[1][2], num_3or9), evalEndo(rho_int, levelstruc_above[1][3], num_3or9)};
            std::vector<ecp> new_tors_16 = {evalEndo(rho_int, levelstruc_above[2][0], NTL::ZZ(32)), evalEndo(rho_int, levelstruc_above[2][1], NTL::ZZ(32))};

            new_levelstruc.push_back(new_tors_2);
            new_levelstruc.push_back(new_tors_3);
            new_levelstruc.push_back(new_tors_16);

            //std::cout << "get weber...." << std::endl;
            bool check = GetWeberOfImage_chain(&w_iso, phi_0, new_levelstruc, Fexts, true);
            assert (check); (void) check;
            assert (j_E1 == NTL::power((NTL::power(w_iso,24)-16),3)/NTL::power(w_iso,24));

            ideal_to_invariant.emplace(std::make_pair(I, w_iso));
            duration_compute_weber = std::chrono::duration_cast<std::chrono::milliseconds>(duration_compute_weber + std::chrono::steady_clock::now() - start_weber);
        } else {
            //auto O2 = I.right_order();
            //assert(O2.is_order());

            // For each of those, connecting idea
            //auto I_conn = connecting_ideal(O0, O2);
            auto I_conn = I_conn0*I;
            I_conn.reset_norm();

            if (isPrincipal(I_conn)) {
                //std::cout << "The ideal I_conn was principal!" << std::endl;

                auto start_weber = std::chrono::steady_clock::now();
                auto I_trans = trans_inv*I*trans;
                auto rho_ideal = I_01 * I_trans; // This is principal

                quat rho(O0.alg);
                bool check_principal = isPrincipal_Compute(&rho, rho_ideal);
                //std::cout << "Done!" << std::endl;
                assert (check_principal); (void) check_principal;

                quat rho_int = rho*(2*q);
                //std::cout << rho_int << std::endl;

                std::vector<std::vector<ecp>> new_levelstruc;
                std::vector<ecp> new_tors_2 = {evalEndo(rho_int, levelstruc_above[0][0], NTL::ZZ(4))};
                std::vector<ecp> new_tors_3 = {evalEndo(rho_int, levelstruc_above[1][0], num_3or9), evalEndo(rho_int, levelstruc_above[1][1], num_3or9), evalEndo(rho_int, levelstruc_above[1][2], num_3or9), evalEndo(rho_int, levelstruc_above[1][3], num_3or9)};
                std::vector<ecp> new_tors_16 = {evalEndo(rho_int, levelstruc_above[2][0], NTL::ZZ(32)), evalEndo(rho_int, levelstruc_above[2][1], NTL::ZZ(32))};

                new_levelstruc.push_back(new_tors_2);
                new_levelstruc.push_back(new_tors_3);
                new_levelstruc.push_back(new_tors_16);

                //std::cout << "running first weber" << std::endl;
                bool check = GetWeberOfLevelStruct(&w_iso, E0, new_levelstruc, Fexts);
                if (!(check)) {
                    //std::cout << "first weber failed" << std::endl;
                    //std::cout << "trying enum all" << std::endl;
                    check = GetWeberOfLevelStruct_j0(&w_iso, new_levelstruc, Fexts);
                }
                assert (check);
                assert (j_E0 == NTL::power((NTL::power(w_iso,24)-16),3)/NTL::power(w_iso,24));

                ideal_to_invariant.emplace(std::make_pair(I, w_iso));
                duration_compute_weber = std::chrono::duration_cast<std::chrono::milliseconds>(duration_compute_weber + std::chrono::steady_clock::now() - start_weber);
                //std::cout << "w = " << w_iso << "\n";
            } else {
                //std::cout << "============" << std::endl;
                // Find a smooth normed connecting ideal
                std::optional<FpE_elem> j;
                // Have to be careful, if O_R(I_conn) == O_0 this currently fails.
                auto start_klpt = std::chrono::steady_clock::now();
                auto gens = KLPT(I_conn, tors_list, NTL::ZZ(6));
                //std::cout << "KLPT done:" << std::endl;
                //std::cout << "beta: " << gens.first << std::endl;
                //std::cout << "Norm of ideal: " << gens.second << std::endl;
                duration_quaternionstuff = std::chrono::duration_cast<std::chrono::milliseconds>(duration_quaternionstuff + std::chrono::steady_clock::now() - start_klpt);

                auto start_isogeny = std::chrono::steady_clock::now(); (void) start_isogeny;
                auto phi_I = idealToIsogeny(gens.first, gens.second, E0, Fexts, bases);
                ec E_I = phi_I.get_codomain();
                j = E_I.j_invariant();

                duration_isogenies = std::chrono::duration_cast<std::chrono::milliseconds>(duration_isogenies + std::chrono::steady_clock::now() - start_isogeny);
                //std::cout << "j: " << *j << std::endl;

                //std::cout << "I_conn0= " << I_conn0.sage() << std::endl;
                //std::cout << "I_01= " << I_01.sage() << std::endl;
                //std::cout << "I= " << I.sage() << std::endl;
                //std::cout << "conj= " << conj << std::endl;
                //std::cout << "gen_01= " << gen_01 << std::endl;

                //Computing endomorphism E0 -> E1 -> E2 -> E0
                auto start_weber = std::chrono::steady_clock::now();

                auto I_02 = create_from_generator(gens.first, gens.second, O0);
                I_02.reset_norm();

                I.reset_norm();
                auto I_trans = trans_inv*I*trans;
                I_trans.reset_norm();
                auto temp_id = I_01 * I_trans;
                temp_id.reset_norm();
                auto rho_ideal = I_02.conjugate() * temp_id;
                rho_ideal.reset_norm();
                quat rho(O0.alg);
                bool check_principal = isPrincipal_Compute(&rho, rho_ideal);
                //std::cout << "Done!" << std::endl;
                assert (check_principal); (void) check_principal;

                quat rho_int = rho*(2*q);
                //std::cout << rho_int << std::endl;

                //New level structure (push through)
                std::vector<std::vector<ecp>> new_levelstruc;
                std::vector<ecp> new_tors_2 = {evalEndo(rho_int, levelstruc_above[0][0], NTL::ZZ(4))};
                std::vector<ecp> new_tors_3 = {evalEndo(rho_int, levelstruc_above[1][0], num_3or9), evalEndo(rho_int, levelstruc_above[1][1], num_3or9), evalEndo(rho_int, levelstruc_above[1][2], num_3or9), evalEndo(rho_int, levelstruc_above[1][3], num_3or9)};
                std::vector<ecp> new_tors_16 = {evalEndo(rho_int, levelstruc_above[2][0], NTL::ZZ(32)), evalEndo(rho_int, levelstruc_above[2][1], NTL::ZZ(32))};

                new_levelstruc.push_back(new_tors_2);
                new_levelstruc.push_back(new_tors_3);
                new_levelstruc.push_back(new_tors_16);

                //std::cout << "Entering getweberimageofchain" << std::endl;
                bool check = GetWeberOfImage_chain(&w_iso, phi_I, new_levelstruc, Fexts, true);
                assert (check); (void) check;
                //std::cout << "okay?" << std::endl;
                assert((*j) == NTL::power((NTL::power(w_iso,24)-16),3)/NTL::power(w_iso,24));

                ideal_to_invariant.emplace(std::make_pair(I, w_iso));
                duration_compute_weber = std::chrono::duration_cast<std::chrono::milliseconds>(duration_compute_weber + std::chrono::steady_clock::now() - start_weber);
                //std::cout << ">>> We have used in total " << duration.count() << " seconds" << std::endl;
                //std::cout << ">>> On quaternion computation (KLPT/bruteforce): " << duration_quaternionstuff.count() << " milliseconds" << std::endl;
                //std::cout << ">>> On computing isogenies: " << duration_isogenies.count() << " milliseconds" << std::endl;
            }
        }
        w_invariants.push_back(w_iso);
    });
    std::cerr << "\r\x1b[KProcessing ideal " << lp1 << " out of " << lp1 << " >>>> DONE!" << std::flush;


    std::cout << "\n\n~~~~~~~~~~~~~~~~~~~~~~~ Summary ~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_full);
    std::cout << ">>> We have used in total " << duration.count() << " seconds" << std::endl;
    std::cout << ">>> On computing field extensions " << duration_compute_fieldexts.count() << " milliseconds" << std::endl;
    std::cout << ">>> On computing endomorphism ring " << duration_compute_starting_endring.count() << " milliseconds" << std::endl;
    std::cout << ">>> On quaternion computation KLPT: " << duration_quaternionstuff.count() << " milliseconds" << std::endl;
    std::cout << ">>> On computing ideal2isogeny: " << duration_isogenies.count() << " milliseconds" << std::endl;
    std::cout << ">>> On computing weber stuff: " << duration_compute_weber.count() << " milliseconds" << std::endl;
    std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;

    return w_invariants;
}

FpX_big_elem ModEvalBigLevelWeber(NTL::ZZ p, NTL::ZZ_pE const w, long l) {
    //////////////////////////////////////////////////////////////////////////////
    /// Weber variant of ModEvalBigLevel
    ///     Implementation of the ideas in Section 3.5 for the BIG LEVEL variant
    //////////////////////////////////////////////////////////////////////////////

    // Initialise mod poly F
    NTL::ZZ_pX F;

    // // Get the bound for prime search
    NTL::ZZ B;
    NTL::ZZ l_ZZ(l);

    if (!(NTL::ProbPrime(l_ZZ))) {
        throw std::invalid_argument("l must be prime!");
    }

    // Number of coeffs
    int l_int = NTL::conv<int>(l_ZZ);
    int Ncoeffs;

    Ncoeffs = l_int + 2;

    NTL::ZZ w_ZZ = NTL::conv<NTL::ZZ>(NTL::ConstTerm(NTL::rep(w))); // Don't know if this is the best way to convert between ZZ_pE to ZZ when it lies in ZZ_p
    FpE_elem w_modq;

    B = GetBoundWeberBigLevel(w_ZZ, l_ZZ);

    // Computing the set of primes
    std::vector<NTL::ZZ> Pl;
    std::cout << "Finding primes..." << std::endl;
    GetPrimesWeber(Pl, w_ZZ, B, l);
    std::cout << "Done!" << std::endl;

    int Nprimes = Pl.size();

    // Initialise crt structure
    crt_info crt;
    std::cout << "Initialising CRT coeffs..." << std::endl;
    crt_init(crt, Pl, Nprimes, Ncoeffs, p);
    std::cout << "Done!" << std::endl;

    // Compute the F mod q and update crt coeffs
    std::cout << "We are working with " << Nprimes << " primes." << std::endl;

    // Note: we process the primes in reverse order since large primes will probably take longer
    std::atomic<ssize_t> next_idx = Pl.size();
    std::mutex mtx;

    auto const fun = [&]() {

        while (true) {

            ssize_t qidx = --next_idx;
            if (qidx < 0)
                break;

            NTL::ZZ const &q_ZZ = Pl[qidx];
            Fp_integer q;
            NTL::conv(q, q_ZZ);

            std::cerr << "Starting with a new prime: " << q << std::endl;

            // In this function we set this to be in ZZX to not work with two moduli in a function
            std::vector<NTL::ZZ> Fq_coeffs(crt.k);

            // Not confusing at all that p is named q, q is named qq... Makes me qqq
            Fp_elem::init(q);
            FpX_elem f;
            SetCoeff(f, 2);
            auto qqs = _avail_qs(q_ZZ, l_ZZ);
            Fp_elem qq;
            NTL::conv(qq, qqs.at(0));
            f[0] = qq;
            FpE_elem::init(f);
            NTL::conv(w_modq, w_ZZ);
            // Compute Fq
            //std::cerr << "Entering Supersingular evaluation..." << std::endl;
            auto invariants = SupersingularEvaluationWeber(q, w_modq, l);
            //std::cerr << "Done!" << std::endl;

            // std::cout << "\n\n~~~~~~ recovered poly: ~~~~~~" << std::endl;
            // std::cout << Phi_l << std::endl;
            // std::cout << "~~~~~~~~~~~~~~~~~~~\n\n" << std::endl;

            std::cerr << "Entering interpolation..." << std::endl;

            FpEX_elem Fq;
            FastInterpolateFromRoots(Fq, invariants);
            //for (auto x : invariants) {
            //    std::cout << x << "\n";
            //}

            std::cout << "~~~~~~~~~~~~~~~~~~~ recovered poly: ~~~~~~~~~~~~~~~~~~~" << std::endl;
            std::cout << Fq << std::endl;
            std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n" << std::endl;

            std::cerr << "Done!" << std::endl;

            // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
            for(int i = 0; i <= deg(Fq); i++){
                assert(is_Fp(NTL::coeff(Fq,i)));
                Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::ConstTerm(NTL::rep(NTL::coeff(Fq, i))));
            }

            //Update CRT sums
            {
                std::lock_guard lock(mtx);
                std::cerr << "Updating CRT coeffs for q=" << q << "... " << std::flush;
                crt_update(crt, qidx, Fq_coeffs, crt.k);
                std::cerr << "Done!\n\n\n" << std::endl;
            }
        }
    };

    {
        std::vector<std::thread> ts;
        for (size_t i = 0; i < num_threads; ++i)
            ts.emplace_back(fun);
        for (auto &t: ts)
            t.join();
    }

    std::cout << "Finalising CRT coeffs..." << std::endl;
    crt_finalise(crt);
    for(int i = 0; i < crt.k; i++){
        NTL::SetCoeff(F, i, NTL::conv<NTL::ZZ_p>(crt.Cdata[i]));
    }
    std::cout << "Done!" << std::endl;

    return F;
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
///////// Implementing the BIG CHARACTERISTIC variant
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////


///*
std::vector<Fp2> SSEvalJinv(const quatlat &id, const std::unordered_map<Key, std::pair<FpE_elem,quatlat>, KeyHash, KeyEqual> &order_jinv_map, const quatlat &O0, const std::list<quatlat> &ell_id_list, const Integer &ell)
{
    (void) ell; //FIXME parameter ell is unused

    std::vector <Fp2> ell_isog_j_list = {};
    auto norm = id.norm().first/id.norm().second;
    assert(NTL::GCD(ell,norm)==1);
    for (auto ellI : ell_id_list) {
        quatlat I = ellI.copy();
        I._intersect(id);
        quatlat O = I.right_order();
        O.reset_norm();
        assert(O.is_order());
        quat ii = {{Integer(0), Integer(1), Integer(0), Integer(0), Integer(1)}, O0.alg};
        Key K = order_invariant_computation(O,&ii);
        auto [j_ell,id_ell] = order_jinv_map.find(K)->second;
        // now we need to decide if we take this one or the frobenius conjugate
        quatlat J = I.conjugate() * id_ell;
        bool a= isPrincipal(J);
        if (a) {
            ell_isog_j_list.push_back(j_ell);
        }
        else {
            quat jj = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, O0.alg};
            assert(O0.contains(jj));
            J._intersect(O0*jj);
            assert(isPrincipal(J));
            Fp_integer some_mod;
            NTL::conv(some_mod, O0.alg.p);
            Fp_push push((Fp_integer(some_mod)));
            {
                assert(rep(j_ell)[0].modulus() == O0.alg.p);
                ell_isog_j_list.push_back(
                    Frob(j_ell)
                    );
            }

        }
    }
    return ell_isog_j_list;

}


void FastSSEvalWeber(std::array<std::vector<ffp2>, 3> &ell_isog_j_list, FastQuatLat &id, const std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<FastInteger, std::pair<std::pair<FastQuat,FastQuat>,std::pair<FastInteger,FastQuat>>>, weber_full_data>>, KeyHash, KeyEqual> &order_jinv_map, const FastInteger &p, const std::vector<std::pair<VerySmallMat,VerySmallMat>> &mat0, const std::vector<FastQuatLat> &ell_id_list, const FastInteger &ell, const weber_enum &enumerator)
{

    Fp_integer some_mod = Fp_integer(p);
    Fp_push push(some_mod);

    assert(GCD(ell, id.norm().first/id.norm().second) == 1);

    assert(ell_id_list.size() == (size_t) ell + 1);

    FastInteger count = 0;

    // init 
    std::pair<VerySmallMat,VerySmallMat> new_web = {{{{0,0},{0,0}}}, {{{0,0},{0,0}}}};
    Key K;
    FastInteger id_norm = id.norm().first/id.norm().second;

    // precomputing some data for fast modular computations throughout the loop 
    SignedBarrettReducer red(ell * id_norm);
    SignedBarrettReducer redid(id_norm);
    SignedBarrettReducer redp (id.alg.p);
    FastInteger inv = InvMod2(ell, redid);
    // FastInteger invsqr = InvModSqr(ell * ell, redid);
    FastInteger inv2sqr = InvMod2Sqr(ell * ell, redid);
    FastQuat fast_commut = {{0, 0, 0, 0, 1}, id.alg};

    for (auto & ellid : ell_id_list) {
        

        FastQuatLat FastI = ellid.copy();
        FastI._fast_intersect(id, red, redid, inv, inv2sqr);
        
        // TODO all coefficients of the gram matrix are divisible by p except the first one (gram[0][0]), and its actually quite close to a multiple of p (it is of the form (1 + p (c^2 + d^2))/n^2 ) where n is the norm of the ideal, and c and d are values mod n^2
        // it would probably be worth considering dividing everything by p and using floats. But it needs to be checked... it may also help the precision requirement....
        auto [FastO, FastGram] = FastI.fast_right_order_and_gram(red);

        fast_commut[0] = 0; fast_commut[1] = 0; fast_commut[2] = 0; fast_commut[3] = 0; fast_commut[4] = 1;
        std::pair<FastQuat,FastQuat> fast_gamma_pair = {fast_commut, fast_commut};

        K = fast_order_invariant_computation_from_gram(FastO, FastGram, fast_gamma_pair);


        // finding the precomputed data corresponding to K 
        auto j_ell_it = order_jinv_map.find(K);

        // somehow the previous computation failed
        // and we go back to the slow but more stable method
        // this should only happen in VERY VERY rare weird cases that were hard to debug
        // and so for now we deal with them that way, it has a negligible impact on concrete performances anyway
        if (j_ell_it == order_jinv_map.end()) {

            // std::cout << "using slow method !\n";

            quatalg Bp = {Integer(p), Integer(1)};
            quatlat testI = ellid.copy().FastQuatLat_to_quatlat(Bp);
            testI.basis[0][2] = Integer(0);
            testI.basis[0][3] = Integer(0);
            testI.basis[1][2] = Integer(0);
            testI.basis[1][3] = Integer(0);
            quatlat testid = id.FastQuatLat_to_quatlat(Bp);
            testid.basis[0][2] = Integer(0);
            testid.basis[0][3] = Integer(0);
            testid.basis[1][2] = Integer(0);
            testid.basis[1][3] = Integer(0);
            testI._fast_intersect(testid);
            
            auto [O, Gram] = testI.fast_right_order_and_gram();
            quat gamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, Bp};
            std::pair<quat,quat> gamma_pair = {gamma, gamma};

            // std::cout << "computing the key....\n";
            K = order_invariant_computation_from_gram(O, Gram, &gamma_pair);
            j_ell_it = order_jinv_map.find(K);
            // std::cout << "done ! \n";

            FastQuatAlg Fast_Bp = FastQuatAlg(Bp); 
            fast_gamma_pair.first = FastQuat(gamma_pair.first, Fast_Bp);
            fast_gamma_pair.second = FastQuat(gamma_pair.second, Fast_Bp);


        } 


        assert(j_ell_it != order_jinv_map.end());
        auto j_ell = j_ell_it->second.first;

        auto normJ = j_ell_it->second.second.first.first;
        auto normI = FastI.norm().first/FastI.norm().second;
        
        auto is_equiv = commutatorfind(fast_gamma_pair, j_ell_it->second.second.first.second.first, {normI, normJ}, j_ell_it->second.second.first.second.second.second, j_ell_it->second.second.first.second.second.first, is_Fp(j_ell), fast_commut, redp);


        new_web = WeberBasApplicationRemoteEndo(mat0, fast_commut);
        

        // now we compute the three weber elements and write them in the result
        for (unsigned char i = 0; i < 3; i++) {
            unsigned char ind_ell = WeberGetFromEnum(enumerator[i].second, new_web);
            ffp2 w_ell = j_ell_it->second.second.second.inv_list[ind_ell];
            if (is_equiv) {
                ell_isog_j_list[i][count] = w_ell;
            }
            else {
                assert(rep(j_ell)[0].modulus() == p);
                ell_isog_j_list[i][count] = fast_Frob(w_ell);
            }
        }

        

        
        count++;   
    }

}

// old algorithm based on NTL integers, it also contains the new fixed sized implementation for debugging purpose
// the api of this function probably needs to be updated to match the most recent developpement of FastSSEvalWeber 
std::vector<std::pair<FpE_elem,std::tuple<bool, Key, std::pair<VerySmallMat, VerySmallMat>>>> SSEvalWeber(const quatlat &id, const std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<Integer, std::pair<std::pair<quat,quat>,std::pair<Integer,quat>>>, weber_full_data>>, KeyHash, KeyEqual> &order_jinv_map, const quatlat &O0, const weber_bas &bas0, const std::vector<std::pair<VerySmallMat,VerySmallMat>> &mat0, const std::list<quatlat> &ell_id_list, const Integer &ell)
{


    (void) bas0; (void) ell; //FIXME parameters bas0 and ell are unused

    std::vector <std::pair<FpE_elem,std::tuple<bool, Key, std::pair<VerySmallMat, VerySmallMat>>>> ell_isog_j_list = {};
    auto norm = id.norm().first/id.norm().second;
    assert(NTL::GCD(ell,norm)==1);

    int total_count = 0;

    clock_t tot_new = 0;
    clock_t tot_old = 0;

    FastInteger id_norm = convert(id.norm().first/id.norm().second);
    SignedBarrettReducer red(convert(ell) * id_norm);
    SignedBarrettReducer redid(id_norm);
    FastInteger inv = InvMod2(convert(ell), redid);
    FastInteger inv2sqr = InvMod2Sqr(convert(ell), redid);
    SignedBarrettReducer redp (convert(id.alg.p));

    for (auto ellI : ell_id_list) {
        total_count++;
        quatlat I = ellI.copy();

    
        FastQuatAlg fast_a = FastQuatAlg(I.alg);
        FastQuatLat FastI = FastQuatLat(I, fast_a);
        FastQuatLat Fastid = FastQuatLat(id, fast_a);

        
        I._fast_intersect(id);
        
        quat gamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        std::pair<quat,quat> gamma_pair = {gamma, gamma};

        auto [O,gram] = I.fast_right_order_and_gram();
        assert(gram[0][0] != 0);

        Key K = order_invariant_computation_from_gram(O, gram, &gamma_pair);

        clock_t tt = tic();
        auto j_ell_it = order_jinv_map.find(K)->second;
        tot_old += (tic() - tt);
        assert(order_jinv_map.find(K) != order_jinv_map.end());
        auto j_ell = j_ell_it.first;

        auto normJ = j_ell_it.second.first.first;
        auto normI = I.norm().first/I.norm().second;
        

        // finding the isomorphism between the two
        
        quat commut = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        auto is_equiv = commutatorfind(gamma_pair, j_ell_it.second.first.second.first, {normI, normJ}, j_ell_it.second.first.second.second.second, j_ell_it.second.first.second.second.first, is_Fp(j_ell), commut);
        


        std::pair<VerySmallMat,VerySmallMat> new_web;

        
        FastI._fast_intersect(Fastid, red, redid, inv, inv2sqr);
        auto [FastO, FastGram] = FastI.fast_right_order_and_gram(red);
        FastMat3 Coords;
        Coords[0][0] = 1;
        Coords[1][1] = 1;
        Coords[2][2] = 1;
        Coords[0][1] = 0;
        Coords[0][2] = 0;
        Coords[1][0] = 0;
        Coords[1][2] = 0;
        Coords[2][1] = 0;
        Coords[2][0] = 0;

#ifndef NDEBUG

        if (!(FastGram == gram)) {
            std::cout << "Gram and FastGram are not equal !\n";
            std::cout << O.basis << "\n";
            std::cout << FastO.basis << "\n";
            
            std::cout << gram << "\n";
            std::cout << FastGram << "\n";   
        }
        assert(FastGram == gram);
        auto test_gram = gram;
        GreedyReduction3(test_gram);
#endif

        FastGreedyReduction3(FastGram, Coords);
        FastQuat fast_gamma = {{0, 0, 0, 0, 1}, fast_a};
        std::pair<FastQuat,FastQuat> fast_gamma_pair = {fast_gamma, fast_gamma};
        fast_gamma_pair.first[4] = FastO.denom;
        fast_gamma_pair.second[4] = FastO.denom;

        for (unsigned i = 1; i <= 3; i++) {
            for (unsigned j = 0; j < 3; ++j) {
                fast_gamma_pair.first[i]  += 2 * Coords[0][j] * FastO.basis[j+1][i];
                fast_gamma_pair.second[i] += 2 * Coords[1][j] * FastO.basis[j+1][i]; 
    
            }
            if (FastGram[0][1] < 0) {
                fast_gamma_pair.second[i] = - fast_gamma_pair.second[i]; 
            }
        }
        fast_gamma_pair.first.normalize();
        fast_gamma_pair.second.normalize();

        std::pair<FastQuat, FastQuat> fast_pair_2 = {FastQuat(j_ell_it.second.first.second.first.first, fast_a), FastQuat(j_ell_it.second.first.second.first.second, fast_a)};
        std::pair<FastInteger, FastInteger> norm_pair = {convert(normI), convert(normJ)};
        FastQuat fast_prod = FastQuat(j_ell_it.second.first.second.second.second, fast_a);
        FastInteger fast_prod_norm = convert(j_ell_it.second.first.second.second.first);

        tt = tic();
        FastQuat fast_commut = {{0, 0, 0, 0, 1}, fast_a};
        auto fast_is_equiv = commutatorfind(fast_gamma_pair, fast_pair_2, norm_pair, fast_prod, fast_prod_norm, is_Fp(j_ell), fast_commut, redp);

        tot_new += (tic() - tt);
        // std::cout << "time up to key new " << (tic()-tt) << "\n\n";
        
#ifndef NDEBUG 

        
        Key K_test;
        NTL::BytesFromZZ(K_test.IntList[0], test_gram[0][0], LenNumberBytes);
        NTL::BytesFromZZ(K_test.IntList[1], test_gram[1][1], LenNumberBytes);
        NTL::BytesFromZZ(K_test.IntList[2], test_gram[2][2], LenNumberBytes);
        assert(is_key_equal(K_test, K));


        Key K_fast;
        FastIntToBytes(FastGram[0][0],K_fast.IntList[0],LenNumberBytes);
        FastIntToBytes(FastGram[1][1],K_fast.IntList[1],LenNumberBytes);
        FastIntToBytes(FastGram[2][2],K_fast.IntList[2],LenNumberBytes);
        assert(is_key_equal(K, K_fast));

        quat quat_test = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        quat_test[0] = Integer(fast_gamma_pair.first[0]);
        quat_test[1] = Integer(fast_gamma_pair.first[1]);
        quat_test[2] = Integer(fast_gamma_pair.first[2]);
        quat_test[3] = Integer(fast_gamma_pair.first[3]);
        quat_test[4] = Integer(fast_gamma_pair.first[4]);
        assert((quat_test + gamma_pair.first).is_zero() || (quat_test + gamma_pair.first * Integer(-1)).is_zero());
        quat_test[0] = Integer(fast_gamma_pair.second[0]);
        quat_test[1] = Integer(fast_gamma_pair.second[1]);
        quat_test[2] = Integer(fast_gamma_pair.second[2]);
        quat_test[3] = Integer(fast_gamma_pair.second[3]);
        quat_test[4] = Integer(fast_gamma_pair.second[4]);
        assert(j_ell == 0 || (quat_test + gamma_pair.second).is_zero() || (quat_test + gamma_pair.second * Integer(-1)).is_zero());
        quat_test[0] = Integer(fast_commut[0]);
        quat_test[1] = Integer(fast_commut[1]);
        quat_test[2] = Integer(fast_commut[2]);
        quat_test[3] = Integer(fast_commut[3]);
        quat_test[4] = Integer(fast_commut[4]);
        if (!(j_ell == 1728 || j_ell == 0 || (quat_test + commut).is_zero() || (quat_test + commut * Integer(-1)).is_zero())) {
            std::cout << "good = " << commut << "\nbad = " << fast_commut << "\n";
        }
        assert(j_ell == 1728 || j_ell == 0 || (quat_test + commut).is_zero() || (quat_test + commut * Integer(-1)).is_zero());
        assert(is_equiv == fast_is_equiv);
        
#endif

        (void) fast_is_equiv;

        new_web = WeberBasApplicationRemoteEndo(mat0, fast_commut);

        if (is_equiv) {
            ell_isog_j_list.push_back({j_ell, {false, K, new_web}});
        }
        else {
            Fp_integer some_mod;
            NTL::conv(some_mod, O0.alg.p);
            Fp_push push(some_mod);
            {
                assert(rep(j_ell)[0].modulus() == O0.alg.p);
                ell_isog_j_list.push_back(
                    {Frob(j_ell), {true, K, new_web}}
                    );
            }
        }   
    }

    std::cout << "mean old time " << (double) tot_old / total_count << "\n";
    // std::cout << "mean time " << (double) tot_new / total_count << "\n";
    

    return ell_isog_j_list;

}



FpX SpecialSupersingularEvaluation(const Integer &p, const Integer &ell, const std::vector<Fp_elem> eval_points)
{
    //////////////////////////////////////////////////////////////////////////////////////
    // Implements SpecialSupersingularEvaluation (Algorithm 1) in the accompanying paper
    //////////////////////////////////////////////////////////////////////////////////////

    // init
    Fp_integer p_mod;
    NTL::conv(p_mod, p);

    Fp::init(p_mod);
    FpX f;
    SetCoeff(f, 2);
    if (p%4 == 3) {
        f[0] = Fp(1);
    }
    else {
        auto qqs = _avail_qs(p, ell);
        Fp qq;
        NTL::conv(qq, qqs.at(0));
        f[0] = Fp(qq);
    }

    Fp2::init(f);

    auto qs = _avail_qs(p, ell);
    assert(qs.size() >= 1);
    auto q = qs.at(0);
    quatalg Bp {p, q};
    std::cout << "q = " << q << "\n";
    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;

    std::unordered_map<Key, std::pair<FpE_elem,quatlat>, KeyHash, KeyEqual> order_jinv_map;
    std::vector <std::pair<quatlat, std::pair<quatlat, Key>>> id_list;

    unsigned k_bound = 20; //Some minimum
    std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
    std::map<unsigned,Fp2k> Fexts;
    {
        for (unsigned k = 1; k <= k_bound; ++k) {
            std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts.emplace(std::make_pair(k, Fp2k {k}));
        }
        std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }


    order_to_jinv_full_list(order_jinv_map, id_list, p, Bp, Fexts, ell);

    // first we enumerate through the set of ell-ideals
    std::vector <Fp2> input_list = {};
    std::vector <Fp2> output_list = {};


    // finding the first generator and the iterator
    bool found = false;
    quat gamma = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, Bp};
    while(!found) {
        NTL::ZZ cofac;
        if (p != 3067) {
            cofac = 3067;
        }
        else {
            cofac = 3319;
        }
        NTL::ZZ target = cofac * ell;
        while (target < 1000 * p) {
            target = cofac * target;
        }
        quat gamma_tmp = RepresentInteger(Bp, target);
        quat testing = {{gamma_tmp[0], gamma_tmp[1], gamma_tmp[2], gamma_tmp[3], NTL::ZZ(2)*gamma_tmp[4]}, Bp};
        found = !O0.contains(testing);
        gamma[0] = gamma_tmp[0];
        gamma[1] = gamma_tmp[1];
        gamma[2] = gamma_tmp[2];
        gamma[3] = gamma_tmp[3];
        gamma[4] = gamma_tmp[4];
    }
    assert(O0.contains(gamma));
    quatlat I = O0 * gamma + O0 * ell;
    assert(std::get<0>(I.norm())/std::get<1>(I.norm()) == ell);


    std::list<int> ell_list = {NTL::conv<int>(ell)};
    // generating the iterating quaternion
    quat beta = find_quaternion_iterator(ell_list, I , O0, Bp);
    // quat betabar = beta.conjugate();


    std::list<quatlat> ell_id_list = O0.left_ideals_of_prime_norm(ell, gamma, beta);

    clock_t tot_time = 0;

    for (auto [id, order_pair] : id_list) {
        if (input_list.size() < ell + 2) {
            FpX sumX;
            Fp2 sum;
            NTL::SetCoeff( sumX, 0, Fp(0) );
            NTL::SetCoeff( sumX, 1, Fp(0) );
            sum = NTL::conv<Fp2>(sumX);
            auto [new_j, new_id] = (order_jinv_map.find(order_pair.second))->second;
            std::vector <Fp2> ell_isog_j_list = {};

            ell_isog_j_list = SSEvalJinv(id, order_jinv_map, O0, ell_id_list, ell);

            clock_t t = clock();
            Fp2X phi_ell_new_j;
            FastInterpolateFromRoots(phi_ell_new_j, ell_isog_j_list);
            // std::cout << " sum= " << sum << " after fast interpolate \n";
            assert(NTL::deg(phi_ell_new_j) == ell + 1);
            for (int i = 0; i <= ell + 1; i++) {
                // std::cout << "i = " << i << "\n";
                // std::cout  << eval_points[i] << "\n";
                sum += NTL::coeff(phi_ell_new_j, i) * eval_points[i];
            }
            // std::cout << "sum com \n";
            input_list.push_back(new_j);
            output_list.push_back(sum);
            // now if the j_inv belong to Fp2 we can also eval for the conjugate
            if (!is_Fp(new_j) && input_list.size() < ell + 2 ) {
                input_list.push_back(Frob(new_j));
                output_list.push_back(Frob(sum));
            }
            tot_time += (clock() -t);


        }
        else {
            goto endloop;
        }

    }
    endloop:

    std::cout << " \n \npoly time = " << (double) (tot_time)/CLOCKS_PER_SEC << "\n \n";

    long long_ell = NTL::conv<long>(ell);

    assert(input_list.size() == (size_t) long_ell + 2 );
    // out of the loop now we can interpolate the final result
    auto res_poly = FastInterpolate(input_list, output_list);
    assert(NTL::deg(res_poly) ==  long_ell + 1 );
    auto result = FpX(long_ell + 1, 1);

    for (long i = 0; i < long_ell + 1; i++) {
        assert(is_Fp(coeff(res_poly, i)));
        NTL::SetCoeff(result, i, NTL::coeff(rep(coeff(res_poly,i)),0));
    }

    return result;

}



FpX SpecialSupersingularEvaluationWeber(const Integer &p, const Integer &ell, const std::vector<Fp_elem> eval_points)
{
    //////////////////////////////////////////////////////////////////////////////
    /// Weber variant of SpecialSupersingularEvaluation
    ///     Implementation of the ideas in Section 3.5 for the BIG CHAR variant
    //////////////////////////////////////////////////////////////////////////////


    // init of the required finite fields
    Fp_integer p_mod;
    NTL::conv(p_mod, p);
    Fp::init(p_mod);
    FpX f;
    SetCoeff(f, 2);
    if (p%4 == 3) {
        f[0] = Fp(1);
    }
    else {
        auto qqs = _avail_qs(p, ell);
        Fp qq;
        NTL::conv(qq, qqs.at(0));
        f[0] = Fp(qq);
    }
    Fp2::init(f);
    auto qs = _avail_qs(p, ell);
    assert(qs.size() >= 1);
    auto q = qs.at(0);

    // quaternion algebra
    quatalg Bp {p, q};
    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;

    // generating of the finite fields
    unsigned k_bound = 20; //Some minimum
    // std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
    std::map<unsigned,Fp2k> Fexts;
    {
        for (unsigned k = 1; k <= k_bound; ++k) {
            // std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts.emplace(std::make_pair(k, Fp2k {k}));
        }
        // std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }


    FastInteger fast_ell = convert(ell);

    clock_t total_time = clock();

    // convert to fast format 
    FastQuatAlg fast_Bp = FastQuatAlg(O0.alg);
    
    // computation of the list of ell ideals
    std::vector<quatlat> ell_id_list = left_ideals_of_prime_norm_O0(ell, Bp);

    // TODO this could be sped up a bit by directly generating the FastQuatLat. There is a function to do so in the fast_quaternions module already, but when I tried to use it, that seemed to create an error
    std::vector<FastQuatLat> fast_ell_id_list = {};
    for (auto & id : ell_id_list) {
        fast_ell_id_list.push_back( FastQuatLat(id, fast_Bp) );
    }
    // some useful precomputation for later on
    SignedBarrettReducer redl(fast_ell);
    for (auto & id : fast_ell_id_list) {
        if (id.good_type) {
            assert(id.basis[0][2] == 0 && id.basis[0][3] == 0 && id.basis[1][2] == 0 && id.basis[1][3] == 0);
            id.basis[0][3] = redl.mod((- id.basis[2][1]) * InvMod(id.basis[2][0], redl));
            id.basis[0][2] = redl.modsqr(InvModSqr( redl.modsqr(id.basis[3][0] * id.basis[3][0] + id.basis[3][1] * id.basis[3][1] - fast_Bp.p), redl) << 1);
            id.basis[1][3] = redl.modsqr(id.basis[0][2] * id.basis[3][1]);
            id.basis[1][2] = redl.modsqr( - id.basis[0][2] * id.basis[3][0]);
        }        
    }

    std::vector <Fp2> input_list = {};
    std::vector <Fp2> output_list = {};


    clock_t tot_interpolation_time = 0;
    clock_t tot_quaternion_time = 0;
    std::unordered_set<Jinv, JinvHash, JinvEqual> input_weber_list = {};
    clock_t t_loop = clock();

    long ellmod24 = NTL::conv<long>(ell) %24;

    std::vector <std::pair<FastQuatLat, Key>> fast_id_list_test;
    std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<FastInteger, std::pair<std::pair<FastQuat,FastQuat>,std::pair<FastInteger,FastQuat>>>, weber_full_data>>, KeyHash, KeyEqual> fast_order_jinv_map_test;
    auto [bas0,mat0] = fast_order_to_weber_inv_full_list(fast_order_jinv_map_test, fast_id_list_test, p, Bp, fast_Bp, Fexts, ell);

    std::cout << "total data_enum (include data conversion to fast format) = " << (double)(clock() - total_time) / (CLOCKS_PER_SEC) << "\n";

    // initializer
    std::array<std::vector<ffp2>, 3> ell_isog_web_list = {std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)}), std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)}), std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)})};

    std::pair<FpX, FpX> phi_ell_new_w;
    phi_ell_new_w.first.SetMaxLength(fast_ell + 1);
    phi_ell_new_w.second.SetMaxLength(fast_ell + 1);

 
    // Ready to start the main loop
    for (auto & [id, id_key] : fast_id_list_test) {
        if (input_list.size() < ell + 2) {
            clock_t tt = clock();
            // TODO this seems to be the only place where we use the enumerator, thus we could store those directly in fast_id_list to reduce the size of order_jinv_map. To be checked, there may be an obstacle
            auto weber_list = (fast_order_jinv_map_test.find(id_key))->second.second.second.inv_list;
            auto weber_enumerator = (fast_order_jinv_map_test.find(id_key))->second.second.second.enumerator;

            // auto weber_enumerator_test = (fast_order_jinv_map_test.find(id_key))->second.second.second.enumerator;

            // we compute the list of ell-isogenous weber invariant (we need only 3 among the 72, the others can be deduced from the three we computed)
            // FastSSEvalWeber(ell_isog_web_list_test, id, fast_order_jinv_map_test, convert(O0.alg.p), mat0_test, fast_ell_id_list, fast_ell, weber_enumerator_test);
            FastSSEvalWeber(ell_isog_web_list, id, fast_order_jinv_map_test, convert(O0.alg.p), mat0, fast_ell_id_list, fast_ell, weber_enumerator);

            // for (int i = 0; i < 3; i++) {
            //     int num_error = 0;
            //     for (int j = 0; j < ell + 1; j++) {
            //         if (ell_isog_web_list[i][j] != ell_isog_web_list_test[i][j]) {
            //             // std::cout << "error at index " << i << " " << j << " " << ell_isog_web_list[i][j] << " " << ell_isog_web_list_test[i][j] <<  "\n";
            //             num_error++;
            //         }
                     
            //     }
            //     std::cout << "num error at i = " << i << " is " << num_error << "\n";
            // }   
            
            tot_quaternion_time += (clock() - tt);
            Fp2 w_ell;

            for (int i = 0; i < 3; i++) {
                
                Fp2 web_inv = Fp2_cast(weber_enumerator[i].first);
                // std::cout << "big3 winv " << web_inv << "\n";
                Jinv ww= JToJinv(web_inv);
                Jinv wwp = JToJinv(Frob(web_inv));
                auto a = input_weber_list.insert(ww);
                auto ap = input_weber_list.insert(wwp);
                bool proceed = a.second && (is_Fp(web_inv) || ap.second) && input_list.size() < ell + 2;

                if (proceed) {
                    clock_t t = clock();
                    // FastInterpolateFromRoots(phi_ell_new_w, ell_isog_web_list[i]);
                    FastInterpolateFromRootsKaratsubaPlusTrick(phi_ell_new_w, ell_isog_web_list[i]);
                    assert(NTL::deg(phi_ell_new_w.first) == ell + 1);

                    // now we compute the sums
                    ffp2 sum;
                    std::array<ffp2, 24> sub_sums;
                    for (int i =0; i < 24; i++) {
                        sub_sums[i] = {Fp(0), Fp(0)};
                    }

                    // now we compute the sums
                    ffp2 coeff;
                    for (int k = 0; k <= ell + 1; k++) {

                        if (k <= deg(phi_ell_new_w.second)) {
                            coeff = {phi_ell_new_w.first[k], phi_ell_new_w.second[k]};
                        }
                        else {
                            coeff = {phi_ell_new_w.first[k], Fp(0)};
                        }
                        
                        // fast_mul(coeff, coeff, {eval_points[k], Fp(0)} );
                        mul(coeff.first, coeff.first, eval_points[k]);
                        mul(coeff.second, coeff.second, eval_points[k]);

                        fast_add(sub_sums[k % 24], sub_sums[k % 24], coeff);

                    }          

                    for (int ii = 0; ii <24 && input_list.size() < ell + 2; ii++) {
                        
                        Fp2 loc_web = Fp2_cast(weber_list[24*i + ii]);
                        Jinv ww_loc= JToJinv(loc_web);
                        Jinv wwp_loc = JToJinv(Frob(loc_web));
                        auto b = input_weber_list.insert(ww_loc);
                        auto bp = input_weber_list.insert(wwp_loc);
                        bool proceed_bis = ((ii == 0) || (b.second && (is_Fp(loc_web) || bp.second))) && input_list.size() < ell + 2;
                        
                        if (proceed_bis) {
                            ffp2 inverse_web = from_Fp2(web_inv);
                            fast_inv(inverse_web, inverse_web);
                            fast_mul(inverse_web, inverse_web, from_Fp2(loc_web));

                            std::vector<ffp2> pow(25); get_powers(pow, inverse_web, 24);
                            sum = {Fp(0), Fp(0)};
                            for (int i = 0; i < 24; i++ ) {
                                fast_mul(coeff, pow[(ellmod24 * ( 25 - i) + 1) % 24], sub_sums[i]);
                                fast_add(sum, sum, coeff);
                                
                            }
                            
                            input_list.push_back(loc_web);
                            // output_list.push_back(sum_test);
                            output_list.push_back(Fp2_cast(sum));
                            // now if the w_inv belong to Fp2 we can also eval for the conjugate
                            if (!is_Fp(loc_web) && input_list.size() < ell + 2 ) {
                                input_list.push_back(Frob(loc_web));
                                output_list.push_back(Frob(Fp2_cast(sum)));
                            }

                        }
                        // else {
                            // std::cout << "not proceed bis !!!!!!!!!!!!!!!! \n";
                            // std::cout << "j = " << fast_order_jinv_map_test.find(id_key)->second.first << " w = " << loc_web << "\n";
                        // }
                    }
                    tot_interpolation_time += (clock() -t);
                }
                else if (input_list.size() >= ell + 2) {
                    goto endloop;
                }
            }


        }
        else {
            goto endloop;
        }

    }
    endloop:

    std::cout << "poly time = " << (double) (tot_interpolation_time)/CLOCKS_PER_SEC << "\n";
    std::cout << "quat time = " << (double)(tot_quaternion_time) / (CLOCKS_PER_SEC) << "\n";    
    std::cout << "total loop time : " << (double) (clock() - t_loop)/(CLOCKS_PER_SEC) << "\n";
    long long_ell = NTL::conv<long>(ell);

    assert(input_list.size() == (size_t) long_ell + 2 );

    // out of the loop now we can interpolate the final result
    auto res_poly = FastInterpolate(input_list, output_list);
    assert(NTL::deg(res_poly) ==  long_ell + 1 );
    auto result = FpX(long_ell + 1, 1);
    for (long i = 0; i < long_ell + 1; i++) {

#ifndef NDEBUG
        if (!is_Fp(coeff(res_poly, i))){
            std::cout << "not in Fp i = " << i << " \n";
            // for (size_t j = 0; j < input_list.size() ; j++) {
            //     std::cout << input_list[j] << " " << output_list[j] << "\n";
            // }
        }
#endif
        assert(is_Fp(coeff(res_poly, i)));
        NTL::SetCoeff(result, i, NTL::coeff(rep(coeff(res_poly,i)),0));
    }

    std::cout << "total computation time = " << (double)(clock() - total_time) / (CLOCKS_PER_SEC) << "\n";

    return result;

}



//
std::pair<long,std::list<std::pair<quatlat, Fp2>>> Reordering_list(const std::vector <std::pair<quatlat, std::pair<quatlat, Key>>> &id_list, const std::unordered_map<Key, std::pair<FpE_elem,quatlat>, KeyHash, KeyEqual> &order_jinv_map ){
    std::list<std::pair<quatlat, Fp2>> reorder_list = {};
    long count_fp2 = 0;
    for (auto const &[id,order_pair] : id_list) {
        auto [new_j, new_id] = (order_jinv_map.find(order_pair.second))->second;
        if (is_Fp(new_j)) {
            // then we append at the end
            reorder_list.push_back({ id, new_j });
        }
        else {
            count_fp2 ++ ;
            reorder_list.push_front({id, new_j});
        }
    }
    return {count_fp2, reorder_list};
}

Integer const_CRT( const std::vector<Integer> &values, const std::vector<Integer> &modulus, const long len) {
    Integer res = values[0];
    Integer mod = modulus[0];
    for (int i = 1; i < len; i++) {
        NTL::CRT(res, mod, values[i], modulus[i]);
    }
    return res;
}

Fp2 const_CRT_polynomials( const std::vector<Fp2> &values, const std::vector<Integer> &modulus, const long len) {
    Integer res0 = NTL::conv<Integer>(coeff(rep(values[0]),0));
    Integer res1 = NTL::conv<Integer>(coeff(rep(values[0]),1));
    Integer mod = modulus[0];
    for (int i = 1; i < len; i++) {
        Integer new_res0 = NTL::conv<Integer>(coeff(rep(values[i]),0));
        Integer new_res1 = NTL::conv<Integer>(coeff(rep(values[i]),1));
        Integer mod_temp = mod;
        NTL::CRT(res0, mod_temp, new_res0, modulus[i]);
        NTL::CRT(res1, mod, new_res1, modulus[i]);
    }
    FpX res_polX;
    NTL::SetCoeff( res_polX, 0, NTL::conv<Fp>(res0));
    NTL::SetCoeff( res_polX, 1, NTL::conv<Fp>(res1));
    return NTL::conv<Fp2>(res_polX);
}

// compute the CRT from FFpi^2 = ZZpi[X] / X^2 + 1 (given as the type ffp2 which is a pair of NTL::_zz_p element) to ZZm[X] / X^2 + 1 where m is the product of the pi (where there are len pi). The product m should be smaller than the bound NTL_SP_BOUND
// values should contain the values of the elements in FFpi^2 
// modulus should contain the pi 
// prod_mod should contain the list {p1, p1 p2 , p1 p2 p3, ....}
// inv should contain the list {p1^{-1} mod p2, (p1p2)^{-1} mod p3, ....}
// assumes that NTL::_zz_p is init with the value of m
ffp2 CRT_ffp2( const std::vector<ffp2> &values, const std::vector<FastInteger> &modulus, const std::vector<FastInteger> &prod_mod, const std::vector<FastInteger> &inv, const long len) {
    FastInteger res0  = values[0].first._zz_p__rep;
    FastInteger res1  = values[0].second._zz_p__rep;
    for (int i = 1; i < len; i++) {
        FastInteger new_res0  = values[i].first._zz_p__rep;
        FastInteger new_res1  = values[i].second._zz_p__rep;    
        res0 = (res0 + (((new_res0 - res0) * inv[i - 1]) % modulus[i] ) * prod_mod[i - 1]) % prod_mod[i]; 
        res1 = (res1 + (((new_res1 - res1) * inv[i - 1]) % modulus[i] ) * prod_mod[i - 1]) % prod_mod[i];
    }

#ifndef NDEBUG 
        for (int i = 0; i < len; i++) {
            // std::cout << "i = " << i << " " << modulus[i] << " " << (res0 + modulus[i]) % modulus[i] << " " << (values[i][array_index][j].first._zz_p__rep + modulus[i]) % modulus[i] << " " << (res1 + modulus[i]) % modulus[i] << " " << (values[i][array_index][j].second._zz_p__rep + modulus[i]) % modulus[i] << "\n";
            assert((res0 - values[i].first._zz_p__rep  ) % modulus[i] == 0);
            assert((res1 - values[i].second._zz_p__rep + modulus[i] ) % modulus[i] == 0);
        }
#endif

    return {NTL::conv<Fp>(res0), NTL::conv<Fp>(res1)};
}

// values contains a vector of len arrays of 3 vectors of length_vec
// the result will be given in the vector crt_values of length length_vec corresponding to the CRT of the vectors at indices array_index in the array
void CRT_ffp2_array_vector(std::vector<ffp2> &crt_values, const std::vector<std::array<std::vector<ffp2>,3>> &values, int array_index, const std::vector<FastInteger> &modulus, const std::vector<FastInteger> &prod_mod, const std::vector<FastInteger> &inv, const long len, const long len_vec) {
    FastInteger res0, res1, new_res0, new_res1;
    for (int j = 0; j < len_vec; j++ ) {
        res0  = values[0][array_index][j].first._zz_p__rep;
        res1  = values[0][array_index][j].second._zz_p__rep;
        for (int i = 1; i < len; i++) {
            new_res0  = values[i][array_index][j].first._zz_p__rep;
            new_res1  = values[i][array_index][j].second._zz_p__rep;    
            res0 = (res0 + (((new_res0 - res0) * inv[i - 1]) % modulus[i] ) * prod_mod[i - 1]) % prod_mod[i]; 
            res1 = (res1 + (((new_res1 - res1) * inv[i - 1]) % modulus[i] ) * prod_mod[i - 1]) % prod_mod[i];
        }

#ifndef NDEBUG 
        for (int i = 0; i < len; i++) {
            // std::cout << "i = " << i << " " << modulus[i] << " " << (res0 + modulus[i]) % modulus[i] << " " << (values[i][array_index][j].first._zz_p__rep + modulus[i]) % modulus[i] << " " << (res1 + modulus[i]) % modulus[i] << " " << (values[i][array_index][j].second._zz_p__rep + modulus[i]) % modulus[i] << "\n";
            assert((res0 - values[i][array_index][j].first._zz_p__rep  ) % modulus[i] == 0);
            assert((res1 - values[i][array_index][j].second._zz_p__rep + modulus[i] ) % modulus[i] == 0);
        }
#endif
        crt_values[j] =  {NTL::conv<Fp>(res0), NTL::conv<Fp>(res1)};
    }    
}

FpX SpecialSupersingularEvaluationCRT(const Integer &p1, const Integer &p2, const Integer ell, const std::vector<Fp_elem> eval_points) {
    //////////////////////////////////////////////////////////////////////////////
    /// Same as SpecialSupersingularEvaluation
    /// but we compute the result mod m = p1 p2
    //////////////////////////////////////////////////////////////////////////////

    Integer m = p1*p2;
    assert(m < NTL_SP_BOUND);
    std::cout << p1 << " " << p2 << "\n";

    // init
    Fp_integer p1_mod;
    NTL::conv(p1_mod, p1);
    Fp::init(p1_mod);
    FpX f1;
    SetCoeff(f1, 2);
    if (p1%4 == 3) {
        f1[0] = Fp(1);
    }
    else {
        f1[0] = Fp(3);
    }
    Fp2::init(f1);
    auto qs1 = _avail_qs(p1, ell);
    assert(qs1.size() >= 1);
    auto q1 = qs1.at(0);
    quatalg Bp1 {p1, q1};
    std::cout << "q1 = " << q1 << "\n";
    auto start1 = starting_curve(Bp1, false);
    quatlat O01 = start1.second;
    Fp_integer p2_mod;
    NTL::conv(p2_mod, p2);
    Fp::init(p2_mod);
    FpX f2;
    SetCoeff(f2, 2);
    if (p2%4 == 3) {
        f2[0] = Fp(1);
    }
    else {
        f2[0] = Fp(3);
    }
    auto qs2 = _avail_qs(p2, ell);
    assert(qs2.size() >= 1);
    auto q2 = qs2.at(0);
    quatalg Bp2 {p2, q2};
    std::cout << "q2 = " << q2 << "\n";
    Fp2::init(f2);
    auto start2 = starting_curve(Bp2, false);
    quatlat O02 = start2.second;

    Fp_integer m_mod;
    NTL::conv(m_mod, m);
    Fp::init(m_mod);
    FpX fm;
    SetCoeff(fm,2);
    std::vector<Integer> vals = {NTL::conv<Integer>(f1[0]), NTL::conv<Integer>(f2[0])};
    std::vector<Integer> modulus = {p1, p2};
    fm[0] = NTL::conv<Fp>(const_CRT(vals, modulus ,2));

    std::unordered_map<Key, std::pair<FpE_elem,quatlat>, KeyHash, KeyEqual> order_jinv_map1,order_jinv_map2;
    std::vector <std::pair<quatlat, std::pair<quatlat, Key>>> id_list1, id_list2;

    Fp::init(p1_mod);
    Fp2::init(f1);

    unsigned k_bound = 20; //Some minimum
    std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
    std::map<unsigned,Fp2k> Fexts1;
    {
        for (unsigned k = 1; k <= k_bound; ++k) {
            std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts1.emplace(std::make_pair(k, Fp2k {k}));
        }
        std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }

    order_to_jinv_full_list(order_jinv_map1, id_list1, p1, Bp1, Fexts1, ell);


    Fp::init(p2_mod);
    Fp2::init(f2);
    std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
    //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
    std::map<unsigned,Fp2k> Fexts2;
    {
        for (unsigned k = 1; k <= k_bound; ++k) {
            std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
            Fexts2.emplace(std::make_pair(k, Fp2k {k}));
        }
        std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
    }

    order_to_jinv_full_list(order_jinv_map2, id_list2, p2, Bp2, Fexts2, ell);

    // first we enumerate through the set of ell-ideals
    Fp::init(m_mod);
    Fp2::init(fm);
    std::vector <Fp2> input_list = {};
    std::vector <Fp2> output_list = {};

    // finding the first generator and the iterator
    bool found = false;
    quat gamma1 = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, Bp1};
    quat gamma2 = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, Bp2};
    while(!found) {
        NTL::ZZ cofac;
        if (p1!= 3067 && p2!=3067) {
            cofac = 3067;
        }
        else {
            cofac = 3319;
        }
        NTL::ZZ target = cofac * ell;
        while (target < 1000 * p1 || target < 1000 * p2) {
            target = cofac * target;
        }
        quat gamma_tmp1 = RepresentInteger(Bp1, target);
        quat testing1 = {{gamma_tmp1[0], gamma_tmp1[1], gamma_tmp1[2], gamma_tmp1[3], NTL::ZZ(2)*gamma_tmp1[4]}, Bp1};
        found = !O01.contains(testing1);

        quat gamma_tmp2 = RepresentInteger(Bp2, target);
        quat testing2 = {{gamma_tmp2[0], gamma_tmp2[1], gamma_tmp2[2], gamma_tmp2[3], NTL::ZZ(2)*gamma_tmp2[4]}, Bp2};
        found = found && !O02.contains(testing2);

        gamma1[0] = gamma_tmp1[0];
        gamma1[1] = gamma_tmp1[1];
        gamma1[2] = gamma_tmp1[2];
        gamma1[3] = gamma_tmp1[3];
        gamma1[4] = gamma_tmp1[4];
        gamma2[0] = gamma_tmp2[0];
        gamma2[1] = gamma_tmp2[1];
        gamma2[2] = gamma_tmp2[2];
        gamma2[3] = gamma_tmp2[3];
        gamma2[4] = gamma_tmp2[4];
    }
    assert(O01.contains(gamma1));assert(O02.contains(gamma2));
    quatlat I1 = O01 * gamma1 + O01 * ell;
    assert(std::get<0>(I1.norm())/std::get<1>(I1.norm()) == ell);
    quatlat I2 = O02 * gamma2 + O02 * ell;
    assert(std::get<0>(I2.norm())/std::get<1>(I2.norm()) == ell);


    std::list<int> ell_list = {NTL::conv<int>(ell)};
    // generating the iterating quaternion
    quat beta1 = find_quaternion_iterator(ell_list, I1 , O01, Bp1);
    quat beta2 = find_quaternion_iterator(ell_list, I2 , O02, Bp2);
    // quat betabar = beta.conjugate();

    std::list<quatlat> ell_id_list1 = O01.left_ideals_of_prime_norm(ell, gamma1, beta1);
    std::list<quatlat> ell_id_list2 = O02.left_ideals_of_prime_norm(ell, gamma2, beta2);

    auto [num_fp21, reorder_list1] = Reordering_list( id_list1, order_jinv_map1 );
    auto [num_fp22, reorder_list2] = Reordering_list( id_list2, order_jinv_map2 );


    long long_ell = NTL::conv<long>(ell);

    clock_t tot_time = 0;



    while (input_list.size() != (size_t)long_ell + 2) {

        // first we pop the first element
        auto [id1,new_j1] = reorder_list1.front();
        auto [id2,new_j2] = reorder_list2.front();


        reorder_list1.pop_front();
        reorder_list2.pop_front();

        std::vector <Fp2> ell_isog_j_list = {};
        std::vector <Fp2> ell_isog_j_list1 = {};
        std::vector <Fp2> ell_isog_j_list2 = {};

        // computing the ell-isogenous j-invariants
        Fp::init(p1_mod);
        Fp2::init(f1);
        ell_isog_j_list1 = SSEvalJinv(id1, order_jinv_map1, O01, ell_id_list1, ell);

        Fp::init(p2_mod);
        Fp2::init(f2);
        ell_isog_j_list2 = SSEvalJinv(id2, order_jinv_map2, O02, ell_id_list2, ell);

        clock_t t = clock();

        assert(ell_isog_j_list1.size() == ell+1);
        assert(ell_isog_j_list2.size() == ell+1);

        // now we combine the two
        Fp::init(m_mod);
        Fp2::init(fm);
        auto new_j = const_CRT_polynomials({new_j1, new_j2}, modulus, 2);
        for (int i=0; i < ell+1; i++) {
            ell_isog_j_list.push_back(const_CRT_polynomials({ell_isog_j_list1[i], ell_isog_j_list2[i]}, modulus ,2));
        }

        FpX sumX;
        Fp2 sum;
        NTL::SetCoeff( sumX, 0, Fp(0) );
        NTL::SetCoeff( sumX, 1, Fp(0) );
        sum = NTL::conv<Fp2>(sumX);
        Fp2X phi_ell_new_j;
        FastInterpolateFromRoots(phi_ell_new_j, ell_isog_j_list);
        assert(NTL::deg(phi_ell_new_j) == ell + 1);
        for (int i = 0; i <= ell + 1; i++) {
            sum += NTL::coeff(phi_ell_new_j, i) * eval_points[i];
        }

        input_list.push_back(new_j);
        output_list.push_back(sum);
        // now if the j_inv belong to Fp2 we can also eval for the conjugate
        if (!is_Fp(new_j1) && !is_Fp(new_j2) && input_list.size() < ell + 2 ) {
            input_list.push_back(Frob(new_j));
            output_list.push_back(Frob(sum));
        }
        tot_time += (clock() - t);


    }


    std::cout << "\n \n poly time = " << (double) (tot_time)/CLOCKS_PER_SEC << "\n \n    ";


    assert(input_list.size() == (size_t) long_ell + 2 );
    // out of the loop now we can interpolate the final result
    auto res_poly = FastInterpolate(input_list, output_list);
    assert(NTL::deg(res_poly) ==  long_ell + 1 );
    auto result = FpX(long_ell + 1, 1);

    for (long i = 0; i < long_ell + 1; i++) {
        assert(is_Fp(coeff(res_poly, i)));
        NTL::SetCoeff(result, i, NTL::coeff(rep(coeff(res_poly,i)),0));
    }

    return result;

}

// CRT variant of SpecialSupersingularEvaluationWeberCRT 
// p is the list of prime modulus
// m is the product of those primes
// num_crt_prime is the number 
// the eval_points are already computed modulo m
FpX SpecialSupersingularEvaluationWeberCRT(const std::vector<FastInteger> &prime_list, const FastInteger &m, int num_crt_prime, const Integer &ell, const std::vector<Fp_elem> eval_points)
{
    //////////////////////////////////////////////////////////////////////////////
    /// CRT variant of SpecialSupersingularEvaluationWeberCRT. Since the primes we use are quite small, they are far below the bound for long integers 
    /// so to gain time in the polynomial interpolation part, we use CRT to combine the elements from several primes into one as 1 operation mod m = p1 p2 
    // is more efficient that 1 operation mod p1 and one operation mod p2 (and same for more primes)
    //////////////////////////////////////////////////////////////////////////////

    std::cout << "SS Eval WeberCRT of size " << num_crt_prime << " with: \n";
    // setting up the polynomials
    std::vector<FpX> mod_polys(num_crt_prime);
    FpX fm;
    for (int i = 0; i < num_crt_prime; i++) {
        Fp_push push(prime_list[i]);
        SetCoeff(mod_polys[i], 2);
        mod_polys[i][0] = Fp(1);
        std::cout << prime_list[i] << " ";
    }
    std::cout << "\n";
    {
        Fp_push push(m);
        SetCoeff(fm, 2);
        fm[1] = Fp(1);
    }
    

    clock_t total_time = clock();
    clock_t tot_interpolation_time = 0;
    clock_t tot_quaternion_time = 0;

    FastInteger fast_ell = convert(ell);
    FastInteger ellmod24 = fast_ell % 24;

    // the final interpolation input and output
    std::vector <Fp2> input_list = {};
    std::vector <Fp2> output_list = {};

    // structure in which is going to be gathered all the values modulo the CRT primes 
    std::vector<std::array<std::vector<ffp2>, 3>> CRT_ell_isog_web_list(num_crt_prime, {std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)}), std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)}), std::vector<ffp2>(fast_ell + 1, {Fp(0), Fp(0)})});


    std::pair<FpX, FpX> phi_ell_new_w;
    phi_ell_new_w.first.SetMaxLength(fast_ell + 1);
    phi_ell_new_w.second.SetMaxLength(fast_ell + 1);

    std::vector<quatalg> Bp(num_crt_prime);
    

    // initializing the quaternion algebras 
    for (int i = 0; i < num_crt_prime; i++) {
        Bp[i] = {Integer(prime_list[i]), Integer(1)};
    }
    std::vector<FastQuatAlg> fast_Bp(num_crt_prime, FastQuatAlg(Bp[0]));
    for (int i = 1; i < num_crt_prime; i++) {
        fast_Bp[i] = FastQuatAlg(Bp[i]);
    }

    std::vector<std::vector<FastQuatLat>> fast_ell_id_list(num_crt_prime);
    std::vector<std::vector <std::pair<FastQuatLat, Key>>> fast_id_list(num_crt_prime);
    std::vector<std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<FastInteger, std::pair<std::pair<FastQuat,FastQuat>,std::pair<FastInteger,FastQuat>>>, weber_full_data>>, KeyHash, KeyEqual>> fast_order_jinv_map(num_crt_prime);
    std::vector<std::vector<std::pair<VerySmallMat,VerySmallMat>>> matrices(num_crt_prime);
    std::vector<std::unordered_set<Jinv, JinvHash, JinvEqual>> input_j_inv_list(num_crt_prime);
    std::vector<std::unordered_set<Jinv, JinvHash, JinvEqual>> weber_tracking_list(num_crt_prime);

    // some precomputation about the modulus
    std::vector<FastInteger> prod_mod = {prime_list[0]};
    std::vector<FastInteger> inverse_mod = {};

    for (int i = 1; i < num_crt_prime; i++) {
        prod_mod.push_back(prod_mod[i - 1] * prime_list[i]);
        SignedBarrettReducer red(prime_list[i]);
        inverse_mod.push_back(InvMod(prod_mod[i - 1], red));
    }

    // first we gather the necessary precomputation
    for (int prime_index = 0; prime_index < num_crt_prime; prime_index++) {
        
        Integer p = Integer(prime_list[prime_index]);

        // init finite field and stuff
        Fp_integer p_mod;
        NTL::conv(p_mod, p);
        Fp::init(p_mod);
        FpX f;
        SetCoeff(f, 2);
        f[0] = Fp(1);
        Fp2::init(f);
        auto start = starting_curve(Bp[prime_index], false);
        quatlat O0 = start.second;

        // generating finite field extensions
        unsigned k_bound = 20; //Some minimum
        // std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
        //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
        std::map<unsigned,Fp2k> Fexts;
        {
            for (unsigned k = 1; k <= k_bound; ++k) {
                // std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
                Fexts.emplace(std::make_pair(k, Fp2k {k}));
            }
            // std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
        }

        // list of ell ideals
        std::vector<quatlat> ell_id_list = left_ideals_of_prime_norm_O0(ell, Bp[prime_index]);
        // transforming the ell ideal list to faster format 
        // TODO could be sped up by directly generating to the fast format, but due to a weird bug it was postponed (this is quite negligible anyway)
        {
            for (auto & id : ell_id_list) {
                fast_ell_id_list[prime_index].push_back( FastQuatLat(id, fast_Bp[prime_index]) );
            }
            SignedBarrettReducer redl(fast_ell);
            for (auto & id : fast_ell_id_list[prime_index]) {
                if (id.good_type) {
                    assert(id.basis[0][2] == 0 && id.basis[0][3] == 0 && id.basis[1][2] == 0 && id.basis[1][3] == 0);
                    id.basis[0][3] = redl.mod((- id.basis[2][1]) * InvMod(id.basis[2][0], redl));
                    id.basis[0][2] = redl.modsqr(InvModSqr( redl.modsqr(id.basis[3][0] * id.basis[3][0] + id.basis[3][1] * id.basis[3][1] - fast_Bp[prime_index].p), redl) << 1);
                    id.basis[1][3] = redl.modsqr(id.basis[0][2] * id.basis[3][1]);
                    id.basis[1][2] = redl.modsqr( - id.basis[0][2] * id.basis[3][0]);
                }        
            }
        }        
        // computation of the list of orders and a list of small ideals we are going to use to collect all the necessary data
        auto [bas0,mat0] = fast_order_to_weber_inv_full_list(fast_order_jinv_map[prime_index], fast_id_list[prime_index], p, Bp[prime_index], fast_Bp[prime_index], Fexts, ell);
        matrices[prime_index] = mat0;

    }

    Fp_integer m_mod;
    NTL::conv(m_mod, m);
    Fp::init(m_mod);
    FpX f;
    SetCoeff(f, 2);
    f[0] = Fp(1);
    Fp2::init(f);

    // std::cout << "total data precomputation (include data conversion to fast format) = " << (double)(clock() - total_time) / (CLOCKS_PER_SEC) << "\n";


    // now we're ready for the enumeration
    // of the weber invariants 
    std::vector<long> id_list_index(num_crt_prime, 0); 
    while (input_list.size() < (size_t) fast_ell + 2) {

        bool is_all_fp2_j_inv = true;
        std::vector<weber_inv_list> list_of_webers(num_crt_prime);

        // computing the ell-isogenous weber structure for all CRT primes
        for (int prime_index = 0; prime_index < num_crt_prime; prime_index++) {

            FastInteger p = prime_list[prime_index];
            Fp_push push(prime_list[prime_index]);
            FpE_push push2(mod_polys[prime_index]);

            clock_t tt = clock();

            // first we select a new ideal that was not already processed
            auto find_weber_precomp = fast_order_jinv_map[prime_index].find(fast_id_list[prime_index][id_list_index[prime_index]].second)->second;
            {
                // this is to decide if we have already handled the galois conjugate 
                // TODO: does this ever happen? because in principle the list is built to avoid that
                Jinv jj_loc= JToJinv(find_weber_precomp.first);
                Jinv jjp_loc = JToJinv(Frob(find_weber_precomp.first));
                auto a = input_j_inv_list[prime_index].insert(jj_loc);
                auto ap = input_j_inv_list[prime_index].insert(jjp_loc);
                bool proceed = a.second && (is_Fp(find_weber_precomp.first) || ap.second) && find_weber_precomp.first != Fp(0) && find_weber_precomp.first != Fp(1728);
                while (!proceed) {
                    id_list_index[prime_index]++;
                    if (id_list_index[prime_index] >= (int) fast_id_list[prime_index].size()) {
                        std::cout << "fast_id_list too small ! \n";
                    }
                    find_weber_precomp = fast_order_jinv_map[prime_index].find(fast_id_list[prime_index][id_list_index[prime_index]].second)->second;
                    jj_loc= JToJinv(find_weber_precomp.first);
                    jjp_loc = JToJinv(Frob(find_weber_precomp.first));
                    a = input_j_inv_list[prime_index].insert(jj_loc);
                    ap = input_j_inv_list[prime_index].insert(jjp_loc);
                    proceed = a.second && (is_Fp(find_weber_precomp.first) || ap.second) && find_weber_precomp.first != Fp(0) && find_weber_precomp.first != Fp(1728);
                }

                // checking if all the j-invariant are Fp2 elements
                is_all_fp2_j_inv = is_all_fp2_j_inv && !is_Fp(find_weber_precomp.first);
            }
            
            // fetching the required precomputed information
            list_of_webers[prime_index] = find_weber_precomp.second.second.inv_list;
            auto weber_enumerator = find_weber_precomp.second.second.enumerator;

            // computing ell-isogenous weber structure
            FastSSEvalWeber(CRT_ell_isog_web_list[prime_index], fast_id_list[prime_index][id_list_index[prime_index]].first, fast_order_jinv_map[prime_index], p, matrices[prime_index], fast_ell_id_list[prime_index], fast_ell, weber_enumerator);
            
            tot_quaternion_time += (clock() - tt);

        }

        // and now exploiting the weber information to compute the interpolation points 
        for (int i = 0; i < 3; i++) {

            // computing the base point 
            ffp2 web_inv_i_CRT;
            std::vector<ffp2> web_inv_i_comp(num_crt_prime);
            for (int index_prime = 0; index_prime < num_crt_prime; index_prime ++) {
                web_inv_i_comp[index_prime] = list_of_webers[index_prime][24 * i]; 
            }
            web_inv_i_CRT = CRT_ffp2(web_inv_i_comp, prime_list, prod_mod, inverse_mod, num_crt_prime);  

            clock_t t_pol = clock();

            std::vector<ffp2> CRT_web_list(fast_ell + 1);
            CRT_ffp2_array_vector(CRT_web_list, CRT_ell_isog_web_list, i, prime_list, prod_mod, inverse_mod, num_crt_prime, fast_ell + 1);

            // interpolation from roots
            FastInterpolateFromRootsKaratsubaPlusTrick(phi_ell_new_w, CRT_web_list);

            // now we compute the evaluation sums
            ffp2 sum;
            std::array<ffp2, 24> sub_sums;
            for (int i =0; i < 24; i++) {
                sub_sums[i] = {Fp(0), Fp(0)};
            }
            // now we compute the sums
            ffp2 coeff;
            for (int k = 0; k <= ell + 1; k++) {
                if (k <= deg(phi_ell_new_w.second)) {
                    coeff = {phi_ell_new_w.first[k], phi_ell_new_w.second[k]};
                }
                else {
                    coeff = {phi_ell_new_w.first[k], Fp(0)};
                }
                
                mul(coeff.first, coeff.first, eval_points[k]);
                mul(coeff.second, coeff.second, eval_points[k]);
                fast_add(sub_sums[k % 24], sub_sums[k % 24], coeff);
            }     

            // compututation of the CRT of the weber invariants + a data telling if all the components are defined over Fp2
            std::vector<ffp2> weber_inv_CRT(24);
            std::vector<bool> is_all_fp2(24, true);
            std::vector<bool> proceed(24, true);

            for (int ii = 0; ii < 24; ii++) {
                std::vector<ffp2> web_inv_comp(num_crt_prime);

                for (int index_prime = 0; index_prime < num_crt_prime; index_prime ++) {
                    auto test = weber_tracking_list[index_prime].insert(JToJinv(Fp2_cast(list_of_webers[index_prime][24 * i + ii]))).second;
                    auto testp = is_Fp(Fp2_cast(list_of_webers[index_prime][24 * i + ii])) || weber_tracking_list[index_prime].insert(JToJinv(Frob(Fp2_cast(list_of_webers[index_prime][24 * i + ii])))).second;
                    assert(test == testp);
                    (void) testp;
                    proceed[ii] = proceed[ii] && test;
                    if (!test) {
                        std::cout << "trying to insert a weber invariant that we already tried \n";
                    }
                    web_inv_comp[index_prime] = list_of_webers[index_prime][24 * i + ii]; 
                    is_all_fp2[ii] = is_all_fp2[ii] && !is_Fp_fast(web_inv_comp[index_prime]);

                }

                weber_inv_CRT[ii] = CRT_ffp2(web_inv_comp, prime_list, prod_mod, inverse_mod, num_crt_prime);  

            }
            
            // local list of things we already saw
            std::unordered_set<Jinv, JinvHash, JinvEqual> input_local_web_inv_list;

            // and now we're ready to compute the output points
            for (int ii = 0; ii < 24 && input_list.size() < ell + 2; ii++) {
                
                ffp2 inverse_web = web_inv_i_CRT;
                fast_inv(inverse_web, inverse_web);
                fast_mul(inverse_web, inverse_web, weber_inv_CRT[ii]);
                std::vector<ffp2> pow(25); get_powers(pow, inverse_web, 24);
                sum = {Fp(0), Fp(0)};
                for (int i = 0; i < 24; i++ ) {
                    fast_mul(coeff, pow[(ellmod24 * ( 25 - i) + 1) % 24], sub_sums[i]);
                    fast_add(sum, sum, coeff);
                }
                // checking that we didn't already added this element 
                // should only happen when on the component curve is defined over Fp
                Jinv w_loc= JToJinv(Fp2_cast(weber_inv_CRT[ii]));
                auto a = input_local_web_inv_list.insert(w_loc);
                bool proceed = a.second; 
                if (proceed) {
                    // and now we push to the input/output lists. Since we are dealing with all Fp2 points we can add directly both the points and its conjugate 
                    // we just need to check 
                    // }
                    input_list.push_back(Fp2_cast(weber_inv_CRT[ii]));
                    output_list.push_back(Fp2_cast(sum));
                    // now if the w_inv belong to Fp2 we can also eval for the conjugate
                    if ((input_list.size() < ell + 2) && is_all_fp2_j_inv) {
                        assert(is_all_fp2[ii]);
                        input_list.push_back(Frob(Fp2_cast(weber_inv_CRT[ii])));
                        output_list.push_back(Frob(Fp2_cast(sum)));
                    }
                    
                }   
            
            }

            tot_interpolation_time += clock() - t_pol;

        }
    }

    // std::cout << "poly time = " << (double) (tot_interpolation_time)/CLOCKS_PER_SEC << "\n";
    // std::cout << "quat time = " << (double)(tot_quaternion_time) / (CLOCKS_PER_SEC) << "\n";    

    long long_ell = fast_ell;

    assert(input_list.size() == (size_t) long_ell + 2 );

    // out of the loop now we can interpolate the final result
    auto res_poly = FastInterpolate(input_list, output_list);
    assert(NTL::deg(res_poly) ==  long_ell + 1 );
    auto result = FpX(long_ell + 1, 1);
    for (long i = 0; i < long_ell + 1; i++) {

#ifndef NDEBUG
        if (!is_Fp(coeff(res_poly, i))){
            std::cout << "not in Fp i = " << i << " \n";
        }
#endif
        assert(is_Fp(coeff(res_poly, i)));
        NTL::SetCoeff(result, i, NTL::coeff(rep(coeff(res_poly,i)),0));
    }

    // std::cout << "total computation time = " << (double)(clock() - total_time) / (CLOCKS_PER_SEC) << "\n";

    return result;

}


FpX_big_elem ModEvalBigCharacteristicWeber(NTL::ZZ p, Fp_big_elem const j, NTL::ZZ l)
{
    //////////////////////////////////////////////////////////////////////////////
    /// Weber variant of ModEvalBigCharacteristic
    ///     Implementation of the ideas in Section 3.5 for the BIG CHAR variant
    //////////////////////////////////////////////////////////////////////////////

    // Initialise mod poly F
    NTL::ZZ_pX F;

    // // Get the bound for prime search
    NTL::ZZ B;
    B = GetBoundWeberBigChar(p, l);

    // Number of coeffs
    Integer j_int = NTL::conv<Integer>(j);
    int l_int = NTL::conv<int>(l);
    int Ncoeffs;

    Ncoeffs = l_int + 2;

    // Constructing js = [j^i mod p for i \in [0, ell+1]]
    std::vector<Integer> js(Ncoeffs);
    js[0] = Integer(1);
    for(int i = 1; i <= l_int+1; i++){
        js[i] = (js[i-1] * j_int) % p;
    }

    // Computing the set of primes
    std::vector<NTL::ZZ> Pl;
    std::cout << "Finding primes..." << std::endl;
    GetPrimesBigCharWeber(Pl, B, l, p);
    std::cout << "Done!\n" << std::endl;

    std::cout << "Set of primes =" << "\n";
    for (auto pp : Pl) {
        std::cout << pp << " ";
    }
    std::cout << "\n";

    int Nprimes = Pl.size();

    // Compute the F mod q and update crt coeffs
    std::cout << "We are working with " << Nprimes << " primes.\n" << std::endl;
    
    // note: we process the primes in reverse order since large primes will probably take longer
    std::atomic<size_t> next_idx = 0;
    std::mutex mtx;
    crt_info crt;

    // non CRT version
    // {
    //     // // Initialise crt structure
        
    //     std::cout << "Initialising CRT coeffs..." << std::endl;
    //     crt_init(crt, Pl, Nprimes, Ncoeffs, p);
    //     std::cout << "Done!" << std::endl;

    //     auto const fun = [&]() {

    //         while (true) {
    //             size_t qidx = ++next_idx;
    //             if (qidx > Pl.size())
    //                 break;

    //             NTL::ZZ const &q = Pl[qidx-1];

    //             std::cerr << "Starting with a new prime: " << q << " this is prime number " << qidx << "/" << Nprimes << std::endl;

    //             // In this function we set this to be in ZZX to not work with two moduli in a function
    //             std::vector<Integer> Fq_coeffs(crt.k);

    //             // Not confusing at all that p is named q, q is named qq... Makes me qqq
    //             Fp_integer q_mod;
    //             NTL::conv(q_mod, q);
    //             Fp::init(q_mod);
    //             FpX f;
    //             SetCoeff(f, 2);
    //             auto qqs = _avail_qs(q, l);
    //             Fp qq;
    //             NTL::conv(qq, qqs.at(0));
    //             f[0] = Fp(qq);
    //             Fp2::init(f);

    //             // Fp tau = NTL::conv<Fp>(j_int);

    //             // computation of the eval list
    //             std::vector<Fp> eval_points = {};
    //             eval_points.push_back(Fp(1));
    //             for (size_t i = 1; i<=l+1; i++) {
    //                 eval_points.push_back(NTL::conv<Fp>(js[i] % q));
    //             }
    //             // eval_points.push_back(tau);
    //             // Fp pow = tau;
    //             assert(eval_points.size() == l + 2);

    //             FpX Fq = SpecialSupersingularEvaluationWeber(q, l, eval_points);
    //             std::cerr << "Done with the polynomial computation!" << std::endl;

    //             // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
    //             for(int i = 0; i <= deg(Fq); i++){
    //                 Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::coeff(Fq, i)); // Don't know if this is the best way to convert between ZZ_pE to ZZ when it lies in ZZ_p
    //             }

    //             //Update CRT sums
    //             {
    //                 std::lock_guard lock(mtx);
    //                 std::cerr << "Updating CRT coeffs for q=" << q << "... " << std::flush;
    //                 crt_update(crt, qidx - 1, Fq_coeffs, crt.k);
    //                 std::cerr << "Done! \n\n" << std::endl;
    //             }
    //         }
    //     };

    //     {
    //         std::vector<std::thread> ts;
    //         for (size_t i = 0; i < num_threads; ++i)
    //             ts.emplace_back(fun);
    //         for (auto &t: ts)
    //             t.join();
    //     }
    // }
    
    int progress = 0;

    {
         // Initialise crt structure
        std::vector<Integer> modulos = {};
        std::vector<std::vector<FastInteger>> prime_list = {};
        std::cout << "Initialising CRT coeffs..." << std::endl;
        batched_crt_init(crt, modulos, prime_list, Pl, Nprimes, Ncoeffs, p);
        std::cout << "Done!\n" << std::endl;

        auto const fun_CRT = [&]() {

            while (true) {
                size_t qidx = ++next_idx;
                if (qidx > modulos.size())
                    break;
                // we treat one CRT 
                NTL::ZZ q;
                q = modulos[qidx-1];

                // In this function we set this to be in ZZX to not work with two moduli in a function
                std::vector<Integer> Fq_coeffs(crt.k);

                // Not confusing at all that p is named q, q is named qq... Makes me qqq
                Fp_integer q_mod;
                NTL::conv(q_mod, q);
                Fp::init(q_mod);
                FpX f;
                SetCoeff(f, 2);
                f[0] = Fp(1);
                Fp2::init(f);

                // computation of the eval list
                std::vector<Fp> eval_points = {};
                eval_points.push_back(Fp(1));
                for (size_t i = 1; i<=l+1; i++) {
                    eval_points.push_back(NTL::conv<Fp>(js[i] % q));
                }
                // eval_points.push_back(tau);
                // Fp pow = tau;
                assert(eval_points.size() == l + 2);
                FpX Fq;
                assert(q < NTL_SP_BOUND);

                progress += (prime_list[qidx - 1]).size();
                std::cout << "progress = " << progress << "/" << Nprimes << "\n";

                Fq = SpecialSupersingularEvaluationWeberCRT(prime_list[qidx - 1], convert(q), (prime_list[qidx - 1]).size(), l, eval_points);
            

                // std::cerr << "Done with the polynomial computation!" << std::endl;

                // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
                for(int i = 0; i <= deg(Fq); i++){
                    Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::coeff(Fq, i)); // Don't know if this is the best way to convert between ZZ_pE to ZZ when it lies in ZZ_p
                }

                //Update CRT sums
                {
                    std::lock_guard lock(mtx);
                    // std::cerr << "Updating CRT coeffs for q=" << q << "... " << std::flush;
                    crt_update(crt, qidx - 1, Fq_coeffs, crt.k);
                    // std::cerr << "Done!" << std::endl;
                    std::cout << "\n" << std::endl;
                }
            }
        };

        // fun_CRT();
        {
            std::vector<std::thread> ts;
            for (size_t i = 0; i < num_threads; ++i)
                ts.emplace_back(fun_CRT);
            for (auto &t: ts)
                t.join();
        }

    
    }

    std::cout << "Finalising CRT coeffs..." << std::endl;
    crt_finalise(crt);
    for(int i = 0; i < crt.k; i++){
        NTL::SetCoeff(F, i, NTL::conv<NTL::ZZ_p>(crt.Cdata[i]));
    }
    std::cout << "Done!\n\n" << std::endl;

    return F;
}
