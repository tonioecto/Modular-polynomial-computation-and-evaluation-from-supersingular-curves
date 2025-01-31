///////////////////////////////////////////////////////////////////////////////////////////////
////   General utility functions. Includes:
////      - function to load modular polynomial of level 2
////      - Sutherland's algorithm for supersingularity testing
////      - DLP algorithm
////      - Sieving for primes
/////////////////////////////////////////////////////////////////////////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>


#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "utils.hpp"

FpEX_elem modPoly2(FpE_elem const &j)
{
    ////////////////////////////////////////////////////////////////
    // Construct the modular polynomial of level 2 evaluated at j
    ////////////////////////////////////////////////////////////////
    FpEX_elem Phi;
    NTL::SetCoeff(Phi, 3);

    auto j2 = j*j;
    auto j3 = j2*j;

    Phi[0] = j3 - 162000*j2 + 8748000000*j - 157464000000000;
    Phi[1] = 1488*j2 + 40773375*j + 8748000000;
    Phi[2] = -j2 + 1488*j - 162000;

    return Phi;
}



static size_t log2floor(Fp_integer const &n)
{
    return NTL::NumBits(n) - 1;
}



bool sutherland(long int const &jj, Fp_integer const &p)
{
    ////////////////////////////////////////////////////////////
    // Sutherland's algorithm for supersingularity testing
    // Extremely fast when j is NOT supersingular
    // Optimized following eprint 2022/880
    ////////////////////////////////////////////////////////////

    //fast algorithm not working for very small p
    if (p < 100){
        return sutherland_slow(jj, p);
    }

    Fp_elem::init(p);
    FpX_elem f;
    NTL::BuildIrred(f, 2);
    FpE_elem::init(f);

    FpE_elem j;
    NTL::conv(j, jj);

    if (j == 1728 && Fp_elem::modulus() % 4 == 3){
        return true;
    }

    FpE_elem j_old[3] {j,j,j}, j_new[3];
    FpE_elem j_ext_old, j_ext_new;
    bool found_j_ext = false;

    FpEX_elem Phi_2 = modPoly2(j);

    // First take a step in each direction
    if (myroots(j_new, Phi_2) < 3) {
        return false;
    }

    //check if there is a j-inv not defined over F_p
    for (int i = 0; i < 3; i++) {
        if (NTL::deg(NTL::rep(j_new[i])) == 1){
            j_ext_old = j_old[i];
            j_ext_new = j_new[i];
            found_j_ext = true;
        }
    }

    //if not, take one more step in each path
    if (!found_j_ext){
        for (int i = 0; i < 3; i++) {

            if (found_j_ext){
                continue;
            }

            Phi_2 = modPoly2(j_new[i]);

            // divide out by the previous step
            FpEX_elem prev;
            NTL::SetCoeff(prev, 1);
            prev[0] = -j_old[i];
            Phi_2 /= prev;

            // Check if we're at the floor
            FpE_elem roots[2];
            if (myroots(roots, Phi_2) < 2) {
                return false;
            }

            // break if we found a j-inv over F_p^2
            for (int ii = 0; ii < 2; ii++) {
                if (NTL::deg(NTL::rep(roots[ii])) == 1){
                    j_ext_old = j_new[i];
                    j_ext_new = roots[ii];
                    found_j_ext = true;
                }
            }    
        }
    }

    //must have found a j-inv over F_p^2 in the supersingular case at this point 
    if (!found_j_ext){
        return false;
    }

    // Then keep going along this path, either down the 2-isog vulcano, or just around the supersingular graph
    size_t bound = log2floor(Fp_elem::modulus())/2;
    for (size_t s = 0; s < bound + 1; s++) {  

        Phi_2 = modPoly2(j_ext_new);

        // divide out by the previous step
        FpEX_elem prev;
        NTL::SetCoeff(prev, 1);
        prev[0] = -j_ext_old;
        Phi_2 /= prev;

        // Check if we're at the floor
        FpE_elem roots[2];
        if (myroots(roots, Phi_2) < 2) {
            return false;
        }

        // update j-lists
        j_ext_old = j_ext_new;
        j_ext_new = roots[0];
    }

    return true;
}


bool sutherland_slow(long int const &jj, Fp_integer const &p)
{
    ////////////////////////////////////////////////////////////
    // Sutherland's algorithm for supersingular testing
    // Extremely fast when j is NOT supersingular
    ////////////////////////////////////////////////////////////

    Fp_elem::init(p);
    FpX_elem f;
    NTL::BuildIrred(f, 2);
    FpE_elem::init(f);

    FpE_elem j;
    NTL::conv(j, jj);

    FpE_elem j_old[3] {j,j,j}, j_new[3];

    FpEX_elem Phi_2 = modPoly2(j);

    // First take a step in each direction
    if (myroots(j_new, Phi_2) < 3) {
        return false;
    }

    // Then keep going, either down the 2-isog vulcano, or just around the SS graph
    size_t bound = log2floor(Fp_elem::modulus());
    for (size_t s = 0; s < bound; s++) {

        for (int i = 0; i < 3; i++) {

            Phi_2 = modPoly2(j_new[i]);

            // divide out by the previous step
            FpEX_elem prev;
            NTL::SetCoeff(prev, 1);
            prev[0] = -j_old[i];
            Phi_2 /= prev;

            // Check if we're at the floor
            FpE_elem roots[2];
            if (myroots(roots, Phi_2) < 2) {
                return false;
            }

            // update j-lists
            j_old[i] = j_new[i];
            j_new[i] = roots[0];
        }
    }

    return true;
}


long _BSGS(ecp const &Q, ecp const &base, long base_order) {

    assert (Q.curve() == base.curve());
    assert (&Q.field() == &base.field());
    FpE_push push(Q.field().F);

    long m = NTL::SqrRoot(base_order);
    m += m * m < base_order;
    assert(m * m >= base_order);

    std::unordered_map<ecp, long> baby_steps;

    ecp G = 0*base;
    for (long i = 0; i < m; i++) {
        baby_steps[G] = i;
        G += base;
    }
    assert(G == m * base);

    ecp min_G = -G;
    assert(min_G == -m * base);

    ecp gamma = Q;

    for (long i = 0; i < m; i++) {
        auto hit = baby_steps.find(gamma);

        if (hit != baby_steps.end()) {
            long sol = (hit->second + m*i) % base_order;
            assert(sol >= 0);
            assert(sol * base == Q);
            return sol;
        }
        gamma += min_G;
    }

    return -1;
}

NTL::ZZ DLP(ecp const &Q, ecp const &base, long ell, long e) {
    ecp gamma = NTL::power(NTL::ZZ(ell), e-1)*base;

    NTL::ZZ ell_to_k(1);
    NTL::ZZ x(0);
    for (long k = 0; k < e; k++) {
        ecp H = NTL::power(NTL::ZZ(ell), e-1-k)*(Q - x*base);
        long d = _BSGS(H, gamma, ell);
        if (d == -1) {
            return NTL::ZZ(-1);
        }
        x += ell_to_k*d;
        ell_to_k *= ell;
    }
    assert(x * base == Q);
    return x;
}

std::unordered_map<int, int> factor(Integer const &N) {
    // Trial division for now
    int p = 1;

    std::unordered_map<int, int> factored;
    auto Ni = N;
    while (Ni > 1) {
        p = NTL::NextPrime(p+1);
        bool first = true;
        while ((Ni % p) == 0) {
            Ni /= p;
            if (first) {
                factored[p] = 1;
                first = false;
            } else {
                factored[p] += 1;
            }
        }
    }

    return factored;
}


bool init_sieve(const std::string& filename, std::vector<long int>& sieve_primes) {
    ////////////////////////////////////////////////////////////
    // Load set of primes for sieving from file
    ////////////////////////////////////////////////////////////

    //open file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "\nError: opening file\n\n";
        return false;
    }

    //read by line
    std::string line;
    while (std::getline(file, line)){
        std::stringstream ss(line);
        std::string ell;
        //separate primes by comma
        while (std::getline(ss, ell, ',')) {
            sieve_primes.push_back(std::stol(ell));
        }
    }

    file.close();
    return true;
}


std::vector<long int> sieve_interval(long int &L, long int len_sieve_interval, std::vector<long int> sieve_primes) {
    ////////////////////////////////////////////////////////////////////////////////
    // Sieve for primes in the interval [L,R] with R = L + len_sieve_interval - 1
    ////////////////////////////////////////////////////////////////////////////////

    NTL::ZZ ZZlen_sieve_interval = NTL::ZZ(len_sieve_interval);
    //sieve up to primes of size sqrt(R)
    NTL::ZZ sqrtRZZ = NTL::SqrRoot(L + ZZlen_sieve_interval);
    long int sqrtR = NTL::conv<long int>(sqrtRZZ);
    //check that we have enough primes for sieving
    assert(sieve_primes.back() >= sqrtR);

    std::vector<bool> numbers(len_sieve_interval, true);

    int i = 0;
    long int ell = sieve_primes[i];
    while (ell <= sqrtR) {
        
        //get sieving offset
        long int k = -(L % ell);
        if (k < 0) {
            k += ell;
        }

        //sieve interval
        while (k < len_sieve_interval) {
            numbers[k] = false;
            k += ell;
        }

        i += 1;
        ell = sieve_primes[i];
    }

    //get primes from true entries
    std::vector<long int> new_primes;
    for (long int i = 0; i < len_sieve_interval; i++) {
        if (numbers[i] == true) {
            new_primes.push_back(i + L);
        } 
    }

    return new_primes;
}
