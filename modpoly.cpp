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
#include "modpoly.hpp"
#include "interpolation.hpp"
#include <unordered_set>


#ifdef NDEBUG
size_t const num_threads = 1 + std::thread::hardware_concurrency();
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
    
    Integer prime = NTL::NextPrime(ell/6 + 250);
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

            FpEX_elem Fq = FastInterpolateFromRoots(j_invariants);
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
    auto webdata = EnumerateAllWeberFast(web, Fexts);
    
    //std::cout << "Getting compatible level structure on E0..." << std::endl;
    bool found_weber = 0;

    for (auto const &webdata_i : webdata){
        auto coeff = webdata_i.second;
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

            FpEX_elem Fq = FastInterpolateFromRoots(invariants);
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


std::pair<mat_Fp,mat_Fp> WeberBasApplicationRemoteEndo(const std::vector<std::pair<mat_Fp,mat_Fp>> &mat0, const quat &gamma, const quatlat &O0) {    
    ///////////////////////////////////////////////////////////////////////////////
    ////// Function applies the matrix of gamma (wrt bas0) to bas
    //////  It assumes that the P16 and Q16 points of bas0 are 
    //////  actually of order 32 so that we can perform the division by two.
    //////  If q = 3, we also assume that P3, Q3 are of order 9
    ///////////////////////////////////////////////////////////////////////////////

    // bool extra2 = (gamma[4] % 2 == 0);
    // bool extra3 = (gamma[4] % gamma.alg.q == 0 && 3 == gamma.alg.q);

    // quat tmp = gamma;
    // // first we apply the endomorphism to the order 16 part of bas0
    // Integer inv_elem = gamma[4]/(extra2 ? (Integer(2)): Integer(1));
    // assert(inv_elem%2 != 0);
    // for (int i = 0; i < 4; i++) {            
    //         tmp[i] = tmp[i] * NTL::InvMod(inv_elem % 16, 16); 
    // }
    // tmp[4] = Integer(1);
    // ecp gammaP016 = evalEndo( tmp, (!extra2 ? (Integer(2)): Integer(1)) * bas0.P16,  Integer(32) / (!extra2 ? (Integer(2)): Integer(1)));
    // ecp gammaQ016 = evalEndo( tmp, (!extra2 ? (Integer(2)): Integer(1)) * bas0.Q16,  Integer(32) / (!extra2 ? (Integer(2)): Integer(1)));
    // assert((16*gammaP016).is_identity());
    // assert(!(8*gammaP016).is_identity());
    // assert((16*gammaQ016).is_identity());
    // assert(!(8*gammaQ016).is_identity());
    // assert(!(8*gammaQ016 - 8*gammaP016).is_identity());

    // // now we do the same for the order 3 part. 
    // inv_elem = gamma[4]/(extra3 ? (Integer(3)): Integer(1));
    // assert(inv_elem%3 != 0);
    // for (int i = 0; i < 4; i++) {            
    //         tmp[i] = gamma[i] * NTL::InvMod(inv_elem % 3, 3); 
    // }
    // tmp[4] = Integer(1);
    // ecp gammaP03 = evalEndo( tmp, bas0.P3,  Integer(9) / (!extra3 ? (Integer(3)): Integer(1)));
    // ecp gammaQ03 = evalEndo( tmp, bas0.Q3,  Integer(9) / (!extra3 ? (Integer(3)): Integer(1)));
    // assert((3*gammaP03).is_identity());
    // assert(!(gammaP03).is_identity());
    // assert((3*gammaQ03).is_identity());
    // assert(!(gammaQ03).is_identity());
    // assert(!(gammaQ03 - gammaP03).is_identity());

    // // now we compute the matrices of order 16 and 3
    // auto M16 = change_of_basis16( {2*bas0.P16, 2*bas0.Q16}, {gammaP016, gammaQ016} );
    // auto M3 = change_of_basis3( {(extra3 ? (Integer(3)): Integer(1)) * bas0.P3, (extra3 ? (Integer(3)): Integer(1)) * bas0.Q3}, {gammaP03, gammaQ03} );
    

    mat_Fp alt_M16;
    mat_Fp alt_M3;
    alt_M16.SetDims(2,2);
    alt_M3.SetDims(2,2);

    // Other way to compute the matrices 
    {
        NTL::mat_ZZ M = O0.basis;
        NTL::ZZ du = O0.denom;
        NTL::ZZ remain,det;
        NTL::vec_ZZ alpha_vec,solve_check;
        alpha_vec.SetLength(4);
        bool res = true;

        for (int i=0; i<4; i++) {
            alpha_vec[i] = du * gamma[i];
            NTL::DivRem(alpha_vec[i],  remain, alpha_vec[i], gamma[4]);
            res = res && (remain==0);
        }
        if (res) {
            solve1(det, solve_check, M, alpha_vec);
            res = res && det == 1;
            assert(!res || solve_check * M == alpha_vec);
        }
        for (int i = 0; i < 4; i++) {

            {
                Fp_push push((Fp_integer(16)));
                alt_M16 = alt_M16 + (NTL::conv<Fp>(solve_check[i]) * mat0[i].first);
            }
            {
                Fp_push push((Fp_integer(3)));
                alt_M3 = alt_M3 + (NTL::conv<Fp>(solve_check[i]) * mat0[i].second);
            }
            
        }
        // assert(alt_M16 == M16);
        // assert(alt_M3 == M3);
    }

    return {alt_M16, alt_M3};
}

mat_Fp invert_mat_16( mat_Fp mat16) {
    Fp det = mat16[0][0] * mat16[1][1] - mat16[1][0] * mat16[0][1];
    det = 1/det;
    mat_Fp inv;
    inv.SetDims(2,2);
    inv[0][0] = mat16[1][1]*det;
    inv[1][1] = mat16[0][0]*det;
    inv[0][1] = - mat16[1][0] * det;
    inv[1][0] = - mat16[0][1] * det;
    return inv;
} 

std::pair<mat_Fp,mat_Fp> invert_mat_pair(const std::pair<mat_Fp,mat_Fp> &mats) {

    mat_Fp M16,M3;
    M16.SetDims(2,2);
    M3.SetDims(2,2);
    {
        Fp_push push((Fp_integer(3)));
        M3 = NTL::inv(mats.second);
    }
    
    {
        Fp_push push((Fp_integer(16)));
        // std::cout << mats.first << " " << Fp::modulus() << "\n";
        M16 = invert_mat_16(mats.first);
    }
    return {M16, M3};
}


quat commutatorfind(const quat &beta1, const quat &beta2) {
    //////////////////////////////////////////////////////////
    //// Finds alpha such that alpha * beta1 = beta2 * alpha
    //// Assumes that beta1, beta2 have denominator 2
    //////////////////////////////////////////////////////////

    // std::cout << beta1.norm().first/beta1.norm().second << " " << beta2.norm().first/beta2.norm().second << "\n"; 
    assert( beta1.norm().first * beta2.norm().second == beta2.norm().first * beta1.norm().second );
    
    NTL::mat_ZZ M;
    M.SetDims(4,4);
    NTL::vec_ZZ v, solve_check;
    v.SetLength(4);
    Integer det;

    auto dx = beta1[0] - beta2[0];
    auto dy = beta1[1] - beta2[1];
    auto dz = beta1[2] - beta2[2];
    auto dt = beta1[3] - beta2[3];

    auto sy = beta1[1] + beta2[1];
    auto sz = beta1[2] + beta2[2];
    auto st = beta1[3] + beta2[3];
    // 
    M[0][0] = dx;
    M[0][1] = -dy;
    M[0][2] = -beta1.alg.p * dz;
    M[0][3] = -beta1.alg.p * dt;

    M[1][0] = dy;
    M[1][1] = dx;
    M[1][2] = beta1.alg.p * st;
    M[1][3] = - beta1.alg.p * sz;

    M[2][0] = dz;
    M[2][2] = dx;
    M[2][1] = - st;
    M[2][3] = sy;

    M[3][0] = dt; 
    M[3][3] = dx;
    M[3][1] = sz;
    M[3][2] = -sy;

    NTL::mat_ZZ newM = NTL::transpose(M);

    // solve1(det, solve_check, M, v);

    NTL::mat_ZZ U;
    U.SetDims(4,4);
    // std::cout << M ;
    auto rank = NTL::image(det, newM, U);
    (void) rank;
    // std::cout << M ;
    // std::cout << U << "\n";
    // std::cout << U[0]*M << "\n";
    assert(rank < 4); 
    // std::cout << M << "\n";

    quat result = {{U[1][0], U[1][1], U[1][2], U[1][3], Integer(1)}, beta1.alg };
    // std::cout << U[0][0] * dt + U[1][0] * sz - U[2][0] * sy + dx * U[3][0] << "\n";
    // std::cout << M*U[0] << "\n"; 
    std::cout << result * beta1 + beta2 * result * Integer(-1) << "\n";
    // assert(result * beta1 == beta2 * result);
    return result;

}



std::vector<std::pair<FpE_elem,std::tuple<bool, Key, std::pair<mat_Fp, mat_Fp>>>> OLDSSEvalWeber(const quatlat &id, const std::unordered_map<Key, std::pair<FpE_elem,std::pair<std::pair<quatlat,quat>, weber_full_data>>, KeyHash, KeyEqual> &order_jinv_map, const quatlat &O0, const weber_bas &bas0, const std::vector<std::pair<mat_Fp,mat_Fp>> &mat0, const std::list<quatlat> &ell_id_list, const Integer &ell)  
{
    ///////////////////////////////////////////////////////////////////
    //// Older version of function SSEvalWeber (slower but bug-free)
    //// For now we use this version
    ///////////////////////////////////////////////////////////////////


    (void) bas0; (void) ell; //FIXME parameters bas0 and ell are unused

    std::vector <std::pair<FpE_elem,std::tuple<bool, Key, std::pair<mat_Fp, mat_Fp>>>> ell_isog_j_list = {}; 
    auto norm = id.norm().first/id.norm().second;
    assert(NTL::GCD(ell,norm)==1);
  
    for (auto ellI : ell_id_list) {
        quatlat I = ellI.copy();
        I._fast_intersect(id);
        quat gamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        
        auto [O,gram] = I.fast_right_order_and_gram(); 
        if (gram[0][0] == 0) {
            gram = compute_gram_order(O);
        }
        Key K = order_invariant_computation_from_gram(O, gram, &gamma); 
    
        
        quat test = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        bool tt1,tt2,tt3,tt4,worked,a;   

        
        // {
        //     quat gamma_I = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        //     quat gamma_J = order_jinv_map.find(K)->second.second.first.second;
        //     assert(order_jinv_map.find(K)->second.second.first.first.contains(gamma_J));
        //     quat iii = {{Integer(0), Integer(1), Integer(0), Integer(0), Integer(1)}, O0.alg};
        //     quat jjj = {{Integer(0), Integer(0), Integer(1), Integer(0), Integer(1)}, O0.alg};
        //     auto small_ngI = get_smallest_element(&gamma_I, I);
        //     assert(I.contains(gamma_I));
        //     std::vector<std::vector<Integer>> tests(4);
        
        //     std::cout << "gamma_I := " << gamma_I << ";\n";
        //     std::cout << "gamma_J := " << gamma_J << ";\n"; 


        //     auto mod = gamma_I[4] * gamma_J[4] * small_ngI;
        //     assert(gamma_I.alg.q == 1);
        //     auto aA = gamma_J[0] * gamma_I[0];
        //     auto aB = gamma_J[0] * gamma_I[1];
        //     auto aC = gamma_J[0] * gamma_I[2];
        //     auto aD = gamma_J[0] * gamma_I[3];
        //     auto bA = - gamma_J[1] * gamma_I[0];
        //     auto bB = - gamma_J[1] * gamma_I[1];
        //     auto bC = - gamma_J[1] * gamma_I[2];
        //     auto bD = - gamma_J[1] * gamma_I[3];
        //     auto cA = - gamma_J[2] * gamma_I[0];
        //     auto cB = - gamma_J[2] * gamma_I[1];
        //     auto cC = - gamma_J[2] * gamma_I[2];
        //     auto cD = - gamma_J[2] * gamma_I[3];
        //     auto dA = - gamma_J[3] * gamma_I[0];
        //     auto dB = - gamma_J[3] * gamma_I[1];
        //     auto dC = - gamma_J[3] * gamma_I[2];
        //     auto dD = - gamma_J[3] * gamma_I[3];  

        //     auto aBpbA = aB + bA;
        //     auto cDmdC = cD - dC;
        //     auto aCpcA = aC + cA;
        //     auto dBmbD = - bD  + dB;
        //     auto aDpdA = aD + dA;
        //     auto bCmcB = bC - cB;
        //     tests[0] = {aA - bB - O0.alg.p* (cC + dD)};
        //     tests[1] = {- aBpbA + O0.alg.p* (cDmdC) };
        //     tests[2] =  { O0.alg.p * (- aCpcA + dBmbD ) };
        //     tests[3] = { O0.alg.p * ( aDpdA - bCmcB ) };
        //     tt1 = (2 * tests[0][0]) % mod == 0;
        //     tt2 = (2 * tests[1][0]) % mod == 0;
        //     tt3 = (2 * tests[2][0]) % mod == 0; 
        //     tt4 = (2 * tests[3][0]) % mod == 0;
        //     bool is_only_one = ((int) tt1 + (int) tt2 + (int) tt3 + (int) tt4) == 1;  
        //         //   
        //     if (is_only_one) {
        //         if (tt1) {
        //             test = {{ tests[0][0], tests[1][0] + 2 *( aBpbA ), aCpcA + dBmbD, aDpdA + bCmcB, mod}, O0.alg};
        //         }
        //         else if (tt2) {
        //             test = {{ tests[1][0], tests[0][0] + 2 * O0.alg.p * (cC + dD) , -aD - cB - bC + dA, aC - dB - bD - cA, mod }, O0.alg}; 
        //         }
        //         else if (tt3) {
        //             test = {{ tests[2][0], O0.alg.p *(aD - bC - cB - dA), tests[0][0] + 2 * ( bB + O0.alg.p * dD) ,  -aB + bA + O0.alg.p * (-dC - cD), mod }, O0.alg}; 
        //         }
        //         else {
        //             test = {{ tests[3][0], O0.alg.p * ( aC + bD - cA + dB),  - aB + bA + O0.alg.p * (cD + dC),  - (tests[0][0] + 2 * (bB + O0.alg.p * cC)), mod }, O0.alg}; 
        //             assert(tt4);
        //         }
        //         worked = I.contains(test) && order_jinv_map.find(K)->second.second.first.first.contains(test.conjugate());
        //     }
        //     else {
        //         // std::cout << "went else" << tt1 << tt2 << tt3 << tt4 << "\n";
        //         test = {{ tests[0][0], tests[1][0] + 2 *( aBpbA ), aCpcA + dBmbD, aDpdA + bCmcB, mod}, O0.alg};
        //         // test = gamma_J.conjugate() * gamma_I;
        //         // test[4] *= small_ngI;
        //         // test.normalize();

        //         // std::cout << test << "\n";
        //         tt1 = tt1 && I.contains(test);
        //         if (true) {
        //         // if (!tt1) {
        //             test = {{ tests[1][0], tests[0][0] + 2 * O0.alg.p * (cC + dD) , -aD - cB - bC + dA, aC - dB - bD - cA, mod }, O0.alg}; 
        //             // test = gamma_J.conjugate() * iii * gamma_I;
        //             // test[4] *= small_ngI;
        //             // test.normalize();
        //         }
        //         // std::cout << test << "\n";
        //         tt2 = tt2 && I.contains(test);
        //         // tt2 = !tt1 && tt2 && I.contains(test);
        //         if (true) {
        //         // if (! (tt1 || tt2) ) {
        //             test = {{ tests[2][0], O0.alg.p *(aD - bC - cB - dA), tests[0][0] + 2 * ( bB + O0.alg.p * dD) ,  -aB + bA + O0.alg.p * (-dC - cD), mod }, O0.alg}; 
        //             // test = gamma_J.conjugate() * jjj * gamma_I;
        //             // test[4] *= small_ngI;
        //             // test.normalize();
        //         }
        //         // std::cout << test << "\n";
        //         tt3 = tt3 && I.contains(test);
        //         // tt3 = !tt1 && !tt2 && tt3 && I.contains(test);
        //         if (true) {
        //         // if (! (tt1 || tt2 || tt3) ) {
        //             test = {{ tests[3][0], O0.alg.p * ( aC + bD - cA + dB),  - aB + bA + O0.alg.p * (cD + dC),  - (tests[0][0] + 2 * (bB + O0.alg.p * cC)), mod }, O0.alg}; 
        //             // test = gamma_J.conjugate() * jjj * iii * gamma_I;
        //             // test[4] *= small_ngI;
        //             // test.normalize();
        //         }
        //         // std::cout << test << "\n";
        //         tt4 = tt4 && I.contains(test);
        //         // tt4 = !tt1 && !tt2 && !tt3 && tt4 && I.contains(test);
        //     }
        //     // in some weird case where all the smallest elements in the ideals have the same norm, it is possible that we were wrong. 
        //     // then we will just use the old method.
        //     worked = (((int) tt1 + (int) tt2 + (int) tt3 + (int) tt4) == 1) && (tt1 || tt2 || tt3 ||tt4) && I.contains(test) && order_jinv_map.find(K)->second.second.first.first.contains(test.conjugate());
        // }  
        (void) worked;
        (void) tt1; (void) tt2; (void) tt3; (void) tt4;      

        auto j_ell_it = order_jinv_map.find(K);
        if (j_ell_it == order_jinv_map.end()) {
            print_key(K);std::cout << "\n";
            std::cout << id << "\n";
            std::cout << ellI << "\n";
            std::cout << O << "\n";
            std::cout << gamma << "\n"; 
        }
        assert(j_ell_it != order_jinv_map.end());
        auto j_ell = j_ell_it->second.first;
        
        // // the curve is defined over Fp
        // if (worked && is_Fp(j_ell) && !(tt1 || tt2)) 
        // {   
        //     std::cout << "rectification ! \n";
        //     quat jj = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, O0.alg};
        //     tt1 = true;
        //     // std::cout << test << "\n";
        //     test = jj * test;
        //     test[4] *= O0.alg.p;
        //     test.normalize();
        //     worked = I.contains(test) && order_jinv_map.find(K)->second.second.first.first.contains(test.conjugate());
        // }


        assert(order_jinv_map.find(K) != order_jinv_map.end());
        I._conjugate();
        quatlat J = I * (order_jinv_map.find(K)->second.second.first.first);
        // std::cout << J << "\n";
        a = isPrincipal_Compute(&gamma, J);
        
        

        if (a) {
            
            std::pair<mat_Fp,mat_Fp> new_web = WeberBasApplicationRemoteEndo(mat0, gamma.conjugate(), O0);
            ell_isog_j_list.push_back({j_ell, { false, K, new_web}});
        }
        else {
            quat jj = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, O0.alg};
            assert(O0.contains(jj));
            J._intersect(O0*jj); 
            bool ap = isPrincipal_Compute(&gamma, J);
            assert(ap); (void) ap;
            std::pair<mat_Fp,mat_Fp> new_web = WeberBasApplicationRemoteEndo(mat0, gamma.conjugate(), O0);
            Fp_integer some_mod;
            NTL::conv(some_mod, O0.alg.p);
            Fp_push push(some_mod);
            {
                assert(rep(j_ell)[0].modulus() == O0.alg.p);
                ell_isog_j_list.push_back(
                    {Frob(j_ell),{true, K, new_web}}
                    );
            }
            
        } 
        // if (worked) {
        //     assert(a == (tt1 || tt2));
        //     if (!(gamma.conjugate() +test).is_zero() && !(gamma.conjugate() + test * Integer(-1)).is_zero()) {
        //         std::cout << "gamma := " << gamma.conjugate() << ";\n";
        //         std::cout << "test := " << test << ";\n";
        //         std::cout << I.basis << "\n";
        //         std::cout << order_jinv_map.find(K)->second.second.first.first.basis << "\n";
        //     }
        //     assert( (gamma.conjugate() +test).is_zero() || (gamma.conjugate() + test * Integer(-1)).is_zero());
        // }
    }

    return ell_isog_j_list;
    
}


/// @brief 
/// @param id 
/// @param order_jinv_map 
/// @param O0 
/// @param bas0 
/// @param mat0 
/// @param ell_id_list 
/// @param ell 
/// @return 
std::vector<std::pair<FpE_elem,std::tuple<bool, Key, std::pair<mat_Fp, mat_Fp>>>> SSEvalWeber(const quatlat &id, const std::unordered_map<Key, std::pair<FpE_elem,std::pair<std::pair<quatlat,quat>, weber_full_data>>, KeyHash, KeyEqual> &order_jinv_map, const quatlat &O0, const weber_bas &bas0, const std::vector<std::pair<mat_Fp,mat_Fp>> &mat0, const std::list<quatlat> &ell_id_list, const Integer &ell)  
{
    ///////////////////////////////////////////////////////////////////
    //// Newer version of function OLDSSEvalWeber 
    //// WARNING: currently unused as it is not bug-free
    ////          To be fixed in the future
    ///////////////////////////////////////////////////////////////////


    (void) bas0; (void) ell; //FIXME parameters bas0 and ell are unused 

    std::vector <std::pair<FpE_elem,std::tuple<bool, Key, std::pair<mat_Fp, mat_Fp>>>> ell_isog_j_list = {}; 
    auto norm = id.norm().first/id.norm().second;
    assert(NTL::GCD(ell,norm)==1);
    // clock_t ecp_time = 0;
    // clock_t find_time = 0;
    // clock_t gamma_time
    // clock_t t;
    quat ii = {{Integer(0), Integer(1), Integer(0), Integer(0), Integer(1)}, O0.alg};
    quat jj = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, O0.alg};
    assert(O0.contains(jj));


    for (auto ellI : ell_id_list) {
        // auto t = tic();
        
        // std::cout << "1st inter " << (tic() - t) << "\n";
        
        // auto t = tic();
        quatlat I = ellI.copy();
        // I._fast_intersect(id);

// #ifndef NDEBUG
//         quatlat I_test = ellI.copy();
//         I_test._intersect(id);
//         for (int i =0; i<4; i++) {
//                     quat bas_el = quat({{I.basis[i][0], I.basis[i][1], I.basis[i][2], I.basis[i][3], I.denom}, I.alg});
//                     if (!I_test.contains(bas_el)) {

//                         std::cout << id.basis << "\n"; 
//                         std::cout << I.basis << "\n";
//                         std::cout << I_test.HNF_basis() << "\n";
//                         std::cout << i << "-th vector not in I_test \n";
//                     }
//                     assert(I_test.contains(bas_el));
//                     auto mtest = I_test.HNF_basis();
//                     quat bas_el_test = quat({{mtest[i][0], mtest[i][1], mtest[i][2], mtest[i][3], I_test.denom}, I.alg});
//                     if (!I.contains(bas_el_test)) {
//                         std::cout << I.basis << "\n";
//                         std::cout << I_test.HNF_basis() << "\n";
//                         std::cout << bas_el_test << "\n";
//                         std::cout << "i=" << i << " is not contained \n";
//                     }
//                 }
        

// #endif         

        auto [O,gram] = I.fast_right_order_and_gram();

     

        quat gamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        quat test = gamma;  
        gram[0][0] = Integer(0);

        if (gram[0][0] == 0) {
            O = I.new_right_order();
            gram = compute_gram_order(O);
        }

//         #ifndef NDEBUG 
//         auto O_test = I.new_right_order();
//         for (int i =0; i<4; i++) {
//                     quat bas_el = quat({{O.basis[i][0], O.basis[i][1], O.basis[i][2], O.basis[i][3], O.denom}, O.alg});
//                     if (!O_test.contains(bas_el)) {

//                         std::cout << I << "\n"; 
//                         std::cout << O << "\n";
//                         std::cout << O_test << "\n";
//                         std::cout << i << "-th vector not in O_test \n";
//                         std::cout << bas_el << "\n";
//                     }
//                     assert(O_test.contains(bas_el));
//                     auto mtest = O_test.basis;
//                     quat bas_el_test = quat({{mtest[i][0], mtest[i][1], mtest[i][2], mtest[i][3], O_test.denom}, I.alg});
//                     if (!I.contains(bas_el_test)) {
//                         std::cout << O << "\n";
//                         std::cout << O_test << "\n";
//                         std::cout << bas_el_test << "\n";
//                         std::cout << "i=" << i << " is not contained in order \n";
//                     }
//                 }
// #endif    

        quat small_OI = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        Key K = order_invariant_computation_from_gram(O, gram, &small_OI);
        auto j_ell = order_jinv_map.find(K)->second.first;
        auto [Jprim,gamma_J] = order_jinv_map.find(K)->second.second.first;
        
        
        // Key K = order_invariant_computation(O); 
        // assert(is_key_equal(K,Ktest));
        // std::cout << "rigo + invariant " << tic() - t << "\n";

        // now we need to decide if we take this one or the frobenius conjugate
        
        // quat gamma_I = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
        
        // auto small_ngI = get_smallest_element(&gamma_I, I);
        // assert(I.contains(gamma_I));


        // std::vector<std::vector<Integer>> tests(4);
        // bool tt1,tt2,tt3,tt4;        

        // auto mod = gamma_I[4] * gamma_J[4] * small_ngI;
        // assert(gamma_I.alg.q == 1);
        // auto aA = gamma_J[0] * gamma_I[0];
        // auto aB = gamma_J[0] * gamma_I[1];
        // auto aC = gamma_J[0] * gamma_I[2];
        // auto aD = gamma_J[0] * gamma_I[3];
        // auto bA = - gamma_J[1] * gamma_I[0];
        // auto bB = - gamma_J[1] * gamma_I[1];
        // auto bC = - gamma_J[1] * gamma_I[2];
        // auto bD = - gamma_J[1] * gamma_I[3];
        // auto cA = - gamma_J[2] * gamma_I[0];
        // auto cB = - gamma_J[2] * gamma_I[1];
        // auto cC = - gamma_J[2] * gamma_I[2];
        // auto cD = - gamma_J[2] * gamma_I[3];
        // auto dA = - gamma_J[3] * gamma_I[0];
        // auto dB = - gamma_J[3] * gamma_I[1];
        // auto dC = - gamma_J[3] * gamma_I[2];
        // auto dD = - gamma_J[3] * gamma_I[3];
        // auto aBpbA = aB + bA;
        // auto cDmdC = cD - dC;
        // auto aCpcA = aC + cA;
        // auto dBmbD = - bD  + dB;
        // auto aDpdA = aD + dA;
        // auto bCmcB = bC - cB;
        // tests[0] = {aA - bB - O0.alg.p* (cC + dD)};
        // tests[1] = {- aBpbA + O0.alg.p* (cDmdC) };
        // tests[2] =  { O0.alg.p * (- aCpcA + dBmbD ) };
        // tests[3] = { O0.alg.p * ( aDpdA - bCmcB ) };
        // tt1 = (2 * tests[0][0]) % mod == 0;
        // tt2 = (2 * tests[1][0]) % mod == 0;
        // tt3 = (2 * tests[2][0]) % mod == 0; 
        // tt4 = (2 * tests[3][0]) % mod == 0;
        // bool is_only_one = ((int) tt1 + (int) tt2 + (int) tt3 + (int) tt4) == 1;  
              
        // if (is_only_one) {
        //     if (tt1) {
        //         test = {{ tests[0][0], tests[1][0] + 2 *( aBpbA ), aCpcA + dBmbD, aDpdA + bCmcB, mod}, O0.alg};
        //     }
        //     else if (tt2) {
        //         test = {{ tests[1][0], tests[0][0] + 2 * O0.alg.p * (cC + dD) , -aD - cB - bC + dA, aC - dB - bD - cA, mod }, O0.alg}; 
        //     }
        //     else if (tt3) {
        //         test = {{ tests[2][0], O0.alg.p *(aD - bC - cB - dA), tests[0][0] + 2 * ( bB + O0.alg.p * dD) ,  -aB + bA + O0.alg.p * (-dC - cD), mod }, O0.alg}; 
        //     }
        //     else {
        //         test = {{ tests[3][0], O0.alg.p * ( aC + bD - cA + dB),  - aB + bA + O0.alg.p * (cD + dC),  - (tests[0][0] + 2 * (bB + O0.alg.p * cC)), mod }, O0.alg}; 
        //         assert(tt4);
        //     }


        // }
        // else {
        //     test = {{ tests[0][0], tests[1][0] + 2 *( aBpbA ), aCpcA + dBmbD, aDpdA + bCmcB, mod}, O0.alg};
        //     tt1 = tt1 &&  I.contains(test);
        //     if (!tt1) {
        //         test = {{ tests[1][0], tests[0][0] + 2 * O0.alg.p * (cC + dD) , -aD - cB - bC + dA, aC - dB - bD - cA, mod }, O0.alg}; 
        //     }

        //     tt2 = !tt1 && tt2 && I.contains(test);
        //     if (! (tt1 || tt2) ) {
        //         test = {{ tests[2][0], O0.alg.p *(aD - bC - cB - dA), tests[0][0] + 2 * ( bB + O0.alg.p * dD) ,  -aB + bA + O0.alg.p * (-dC - cD), mod }, O0.alg};
        //     }

        //     tt3 = !tt1 && !tt2 && tt3 && I.contains(test);
        //     if (! (tt1 || tt2 || tt3) ) {
        //         test = {{ tests[3][0], O0.alg.p * ( aC + bD - cA + dB),  - aB + bA + O0.alg.p * (cD + dC),  - (tests[0][0] + 2 * (bB + O0.alg.p * cC)), mod }, O0.alg}; 
        //     }

        //     tt4 = !tt1 && !tt2 && !tt3 && tt4 && I.contains(test);
        // }
        // in some weird case where all the smallest elements in the ideals have the same norm, it is possible that we were wrong. 
        // then we will just use the old method.
        // bool worked = (tt1 || tt2 || tt3 ||tt4);
        bool worked = false;
        bool a;


        if (worked) {
            // a = tt1 || tt2;
        }
        else {
            quatlat I = ellI.copy();
            I._intersect(id);
            O = I.new_right_order();
            K = order_invariant_computation(O,&gamma); 
            Jprim = order_jinv_map.find(K)->second.second.first.first;
            I._conjugate();
            Jprim = I * (Jprim);
            a = isPrincipal_Compute(&gamma, Jprim);
        }
        
        // t = tic();
        // quat commut = commutatorfind(order_jinv_map.find(K)->second.second.first.second * n1 * n2, small_OI * n1 * n2 );
        // std::cout << "commutator " << tic() - t << "\n";
        // assert(nJ.first/nJ.second * ngI.first/ngI.second == ngJ.first/ngJ.second * nI.first/nI.second);
    
        std::pair<mat_Fp,mat_Fp> new_web;
        if (a) {
            if (worked) {
                new_web = WeberBasApplicationRemoteEndo(mat0, test, O0);
            }
            else {
                new_web = WeberBasApplicationRemoteEndo(mat0, gamma, O0);
            }
            ell_isog_j_list.push_back({j_ell, { false, K, new_web}});
        }
        else {
            if (worked) {
                new_web = WeberBasApplicationRemoteEndo(mat0, test, O0);
            }
            else {
                Jprim._intersect(O0*jj); 
                bool ap = isPrincipal_Compute(&gamma, Jprim);
                assert(ap); (void) ap;
                new_web = WeberBasApplicationRemoteEndo(mat0, gamma.conjugate(), O0);
                NTL::ZZ_pPush push(O0.alg.p);
                new_web = WeberBasApplicationRemoteEndo(mat0, gamma.conjugate(), O0);
            }
            
            NTL::ZZ_pPush push(O0.alg.p);
            {
                assert(rep(j_ell)[0].modulus() == O0.alg.p);
                ell_isog_j_list.push_back(
                    {Frob(j_ell),{true, K, new_web}}
                    );
            }
            
        } 
    }

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

    // std::cout << "before the loop \n";

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
            Fp2X phi_ell_new_j = FastInterpolateFromRoots(ell_isog_j_list);
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
    // std::cout << "q = " << q << "\n";
    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;



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

    // std::list<quatlat> ell_id_list = O0.left_ideals_of_prime_norm(ell, gamma, beta);
    std::list<quatlat> ell_id_list = left_ideals_of_prime_norm_O0(ell, gamma, beta);

    // for (auto idl : ell_id_list) {
    //     std::cout << idl.HNF_basis() << "\n" << idl.denom << "\n";
    //     quatlat O = idl.right_order() ;
    //     std::cout << O.HNF_basis() << "\n" << O.denom << "\n\n";
    // }


    std::unordered_map<Key, std::pair<FpE_elem,std::pair<std::pair<quatlat,quat>, weber_full_data>>, KeyHash, KeyEqual> order_jinv_map;
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


    clock_t total_time = clock();

    auto [bas0,mat0] = order_to_weber_inv_full_list(order_jinv_map, id_list, p, Bp, Fexts, ell);

    std::cout << "total data_enum = " << (double)(clock() - total_time) / (CLOCKS_PER_SEC) << "\n";

    // first we enumerate through the set of ell-ideals 
    std::vector <Fp2> input_list = {};
    std::vector <Fp2> output_list = {};

    

    clock_t tot_interpolation_time = 0;
    clock_t tot_quaternion_time = 0;
    clock_t tot_weber_fetch_time = 0; 

    std::vector<std::vector<std::vector<std::pair<Integer,Integer>>>> coeff_list = EnumerateAllWeberCoeff();

    std::unordered_set<Jinv, JinvHash, JinvEqual> input_weber_list = {}; 
    clock_t t_loop = clock();
    
    long ellmod24 = NTL::conv<long>(ell) %24;

    int count_curve = 0;

    for (auto [id, order_pair] : id_list) {
    
        if (input_list.size() < ell + 2) {
            // std::cout << "treating j = " << new_j << "\n"; 
            clock_t tt = clock();
            count_curve++;
            // auto [new_j, new_web_id] = (order_jinv_map.find(order_pair.second))->second;
            // auto [new_id, new_web] = new_web_id;
            // assert((16*new_web.basis.P16).is_identity());
            // assert((3*new_web.basis.P3).is_identity());

            std::vector<std::pair<Fp2,std::tuple<bool, Key, std::pair<mat_Fp,mat_Fp>>>> ell_isog_j_list = {}; 
            
            ell_isog_j_list = OLDSSEvalWeber(id, order_jinv_map, O0, bas0, mat0, ell_id_list, ell); 
            
            auto weber_list = (order_jinv_map.find(order_pair.second))->second.second.second.enumerator;
            tot_quaternion_time += (clock() - tt);
            
            // auto weber_list = EnumerateAllWeberFast(new_web.basis, Fexts);
            // auto weber_list2 = EnumerateAllWeber(new_web, coeff_list, Fexts);
            
            Fp2 w_ell;


            for (int i = 0; i < 3; i++) {
                Fp2 web_inv = weber_list[24*i].first;
                Jinv ww= JToJinv(web_inv);
                Jinv wwp = JToJinv(Frob(web_inv));
                auto a = input_weber_list.insert(ww);
                auto ap = input_weber_list.insert(wwp);
                bool proceed = a.second && (is_Fp(web_inv) || ap.second) && input_list.size() < ell + 2;
                if (proceed) {
                    std::vector<Fp2> ell_weber_list = {};
                    
                    for (unsigned j = 0; j<ell_isog_j_list.size(); j++) {
                        // auto [j_ell, webbool_ell]: ell_isog_j_list
                        clock_t ttt = clock();
                        // auto [web_ell, need_frob, enumm, mat_pair] = webbool_ell;
                        w_ell = WeberGetFromEnum(weber_list[24*i].second, order_jinv_map.find(std::get<1>(ell_isog_j_list[j].second))->second.second.second.enumerator, std::get<2>(ell_isog_j_list[j].second));
                        tot_weber_fetch_time += (clock() - ttt);
                        bool check = true;
                        // Fp2 w_test;
                        // check = BasCoeffToWeber(&w_test, web_ell, coeffs, Fexts);
                        // assert((power(w_ell,24) - 16)/power(w_ell,8) == (power(w_test,24) - 16)/power(w_test,8));
                        // assert(power(w_ell,24) == power(w_test,24));
                        // assert(w_ell == w_test);
                        // std::cout << power(w_test,24) << " " << power(w_ell,24) << "   " << (power(w_test,24) - 16)/power(w_test,8) << " " << (power(w_ell,24) - 16)/power(w_ell,8) << "    " << power(w_test,8) << " " << power(w_ell,8) << "\n";
                        // std::cout << " j_ell = " << j_ell << " need frob = " << need_frob << "\n";
                        if (check) {
                            if (std::get<0>(ell_isog_j_list[j].second)) {
                                w_ell = Frob(w_ell);
                            }
                            ell_weber_list.push_back(w_ell);   
                        }
                    }
                    if (ell_weber_list.size() == ell + 1 ) {
                        clock_t t = clock();

                        Fp2X phi_ell_new_w = FastInterpolateFromRoots(ell_weber_list);
                        assert(NTL::deg(phi_ell_new_w) == ell + 1);

                        // now we compute the sums
                        FpX sumX;
                        Fp2 sum;
                        

                        for (long ii = 0; ii <24 && input_list.size() < ell + 2; ii++) {
                            
                            Fp2 loc_web = weber_list[24*i + ii].first;
                            Jinv ww_loc= JToJinv(loc_web);
                            Jinv wwp_loc = JToJinv(Frob(loc_web));
                            auto b = input_weber_list.insert(ww_loc);
                            auto bp = input_weber_list.insert(wwp_loc);
                            bool proceed_bis = b.second && (is_Fp(loc_web) || bp.second) && input_list.size() < ell + 2;
                            if (proceed_bis) {
                                Fp2 inverse_web = loc_web/web_inv; 
                                sum = Fp2(0);
                                std::vector<Fp2> pow = get_powers(inverse_web, 24);
                                for (long k = 0; k <= ell + 1; k++) {
                                    long k_pow = (48000000 + ellmod24 * ( 1 - k ) + 1)%24;
                                    sum += NTL::coeff(phi_ell_new_w, k) * pow[k_pow] * eval_points[k];
                                }
                                // std::cout << "\n";
                                input_list.push_back(loc_web);
                                output_list.push_back(sum);
                                // now if the w_inv belong to Fp2 we can also eval for the conjugate
                                if (!is_Fp(loc_web) && input_list.size() < ell + 2 ) {
                                    input_list.push_back(Frob(loc_web));
                                    output_list.push_back(Frob(sum));
                                }
                            }

                            
                        }
                        tot_interpolation_time += (clock() -t);
                    }
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
    // std::cout << "amortized quat time = " << (double)(tot_quaternion_time) / (CLOCKS_PER_SEC * count_curve) << "\n";
    std::cout << "fetch weber time = " << (double)(tot_weber_fetch_time) / (CLOCKS_PER_SEC) << "\n";
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
            for (size_t j = 0; j < input_list.size() ; j++) {
                std::cout << input_list[j] << " " << output_list[j] << "\n";
            }
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
        Fp2X phi_ell_new_j = FastInterpolateFromRoots(ell_isog_j_list);
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
    std::cout << "Done!" << std::endl;

    int Nprimes = Pl.size();

    // Initialise crt structure
    crt_info crt;
    std::cout << "Initialising CRT coeffs..." << std::endl;
    crt_init(crt, Pl, Nprimes, Ncoeffs, p);
    std::cout << "Done!" << std::endl;


    // Compute the F mod q and update crt coeffs
    std::cout << "We are working with " << Nprimes << " primes." << std::endl;

    // note: we process the primes in reverse order since large primes will probably take longer
    std::atomic<size_t> next_idx = 0;
    std::mutex mtx;

    auto const fun = [&]() {

        while (true) {
            size_t qidx = ++next_idx;
            if (qidx > Pl.size())
                break;

            NTL::ZZ const &q = Pl[qidx-1];
            
            std::cerr << "Starting with a new prime: " << q << std::endl;

            // In this function we set this to be in ZZX to not work with two moduli in a function
            std::vector<Integer> Fq_coeffs(crt.k); 

            // Not confusing at all that p is named q, q is named qq... Makes me qqq
            Fp_integer q_mod;
            NTL::conv(q_mod, q);
            Fp::init(q_mod);
            FpX f;
            SetCoeff(f, 2);
            auto qqs = _avail_qs(q, l);
            Fp qq;
            NTL::conv(qq, qqs.at(0));
            f[0] = Fp(qq);
            Fp2::init(f);

            // Fp tau = NTL::conv<Fp>(j_int);

            // computation of the eval list
            std::vector<Fp> eval_points = {};
            eval_points.push_back(Fp(1));
            for (size_t i = 1; i<=l+1; i++) {
                eval_points.push_back(NTL::conv<Fp>(js[i] % q));
            }
            // eval_points.push_back(tau);
            // Fp pow = tau;
            assert(eval_points.size() == l + 2);

            FpX Fq = SpecialSupersingularEvaluationWeber(q, l, eval_points);
            std::cerr << "Done with the polynomial computation!" << std::endl;

            // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
            for(int i = 0; i <= deg(Fq); i++){
                Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::coeff(Fq, i)); // Don't know if this is the best way to convert between ZZ_pE to ZZ when it lies in ZZ_p
            }

            //Update CRT sums
            {
                std::lock_guard lock(mtx);
                std::cerr << "Updating CRT coeffs for q=" << q << "... " << std::flush;
                crt_update(crt, qidx - 1, Fq_coeffs, crt.k);
                std::cerr << "Done! \n" << std::endl;
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
    std::cout << "Done!\n\n" << std::endl;

    return F;
}
