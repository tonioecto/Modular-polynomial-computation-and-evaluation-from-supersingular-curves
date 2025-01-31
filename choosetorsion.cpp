///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////// The main function in this file is choose_torsion, which allows us to create a 
/////// 'dictionary' of prime powers \ell^e and the field extension F_{p^{2k}} over which 
/////// the ell^e-torsion is defined.
/////// To build this 'dictionary' we want to choose the prime powers that we work with 
/////// somewhat optimally so that the cost of EC arithmetic is minimized.
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <algorithm>

#include "choosetorsion.hpp"
#include "costmodel.hpp"


int in_vector(factor_list const &vec, int const &ell, int const &exp) {
    ////////////////////////////////////////////////////////////////////////////////////
    // This function checks if a vector has an entry of the form [ell, exp, ....]
    // Returns index of respective entry or -1 if no corresponding entry
    ////////////////////////////////////////////////////////////////////////////////////

    int found = 0;

    for (const auto& tup : vec) {
        if (std::get<0>(tup) == ell && std::get<1>(tup) == exp){
            return found;
        }
        found += 1;
    }
    return -1;
}

std::vector<std::tuple<int, int>> smooth_part(NTL::ZZ const &N, NTL::ZZ const &B) {
    //////////////////////////////////////////////////////
    // This function performs trial division up to B
    //////////////////////////////////////////////////////


    int p = 1;
    int pow;
    NTL::ZZ le;

    std::vector<std::tuple<int, int>> factored;

    while (p <= B) {
        p = NTL::NextPrime(p+1);
        pow = 0;
        le = p;
        while ((N % le) == 0) {
            pow += 1;
            le *= p;
        }
        if (pow > 0){
            std::tuple<int, int> tup;
            std::get<0>(tup) = NTL::conv<int>(p);
            std::get<1>(tup) = pow;
            factored.push_back(tup);
        }
    }
    return factored;
}

bool sort_factors(ell_tuple const &tup1, ell_tuple const &tup2)
{
    //////////////////////////////////////////////////////
    // This function sorts by last element of tuples
    //////////////////////////////////////////////////////
    return std::get<3>(tup1) < std::get<3>(tup2);
}

void is_prime_power(NTL::ZZ &ell, int &exp, NTL::ZZ const &N)
// WARNING: this funtion assumes that N is not a prime!
{
    ell = NTL::ZZ(2);
    NTL::ZZ bound = SqrRoot(N) + 1;

    while (ell <= bound){
        if (N % ell == 0){
            exp = 1;
            while (N % power(ell, exp+1) == 0){
                exp += 1;
            }

            if (N != power(ell, exp)){
                //return ell = 0 to indicate N is not a prime power
                ell = NTL::ZZ(0);
            }
            break;
        }
        ell = NTL::NextPrime(ell+1);
    }
}


factor_list choose_torsion_naive(NTL::ZZ const &p, NTL::ZZ const &tors_bound, NTL::ZZ const &coprime)
{
    factor_list factors;
    model_function model = cost_model(p);
    NTL::ZZ le = NTL::ZZ(1); 
    NTL::ZZ prod = NTL::ZZ(1);
    NTL::ZZ ell;
    int exp;

    while (prod < tors_bound * 2) {

        le += 1;

        if (NTL::GCD(le, coprime*p) > 1) {
            continue;
        }

        //only use le if it is prime or a prime power
        if (!NTL::ProbPrime(le)) {
            is_prime_power(ell, exp, le);
            if (ell == 0) {
                continue;
            }
        } else {
            ell = le;
            exp = 1;
        }

        //find smallest extension to get le = ell^exp torsion
        int k = 1;
        //Note: first arg of PowerMod must be smaller than third

        /* Related to twists, TODO: implement this this
        while (NTL::PowerMod(p%le, NTL::ZZ(k), le) != 1){
            k += 1;
        }

        
        if (k % 2 == 0 && NTL::PowerMod(p%le, k/2, le) - le == -1) {
            k /= 2;
        }
        Naive way below
        */

        // TODO: Remove for commented out above when twists are working
        NTL::ZZ ord_k(p + 1);
        while (ord_k % le != 0){
            k += 1;
            ord_k = NTL::power(p, long(k)) - NTL::power_long(long(-1), long(k % 2));
        }

        NTL::ZZ temp = NTL::GCD(prod, le);
        prod *= le/temp;

        ell_tuple tup;
        int ellint = NTL::conv<int>(ell);
        std::get<0>(tup) = ellint;
        std::get<1>(tup) = exp;
        std::get<2>(tup) = k;
        std::get<3>(tup) = adj_cost(model(ellint, exp, k), ellint);
        factors.push_back(tup);
    }

    return factors;
}


factor_list choose_torsion(NTL::ZZ const &p, NTL::ZZ const &tors_bound, NTL::ZZ const &coprime)
    /////////////////////////////////////////////////////////////////////////////////////////
    //////// Given a prime p, a bound tors_bound = B, it creates a dictionary 
    ////////    ell^e -> k, which means that E[ell^e] is defined over Fp2k
    ////////    We choose the ell^e in this dictionary so that their product is bigger than B 
    ////////    but the cost is as small as possible. 
    ////////
    ////////    The input 'coprime' allows us to specify that we only want torsion
    ////////    that is coprime to a certain number.
    ////////
    /////////////////////////////////////////////////////////////////////////////////////////
{
    factor_list factors;
    model_function model = cost_model(p);
    NTL::ZZ T;
    NTL::RR max_cost;
    int k = 0;
    
    // Might need to ensure q is part of torsion...
    factors = choose_torsion_naive(p, tors_bound, coprime);

    while (true){
        std::sort(factors.begin(), factors.end(), sort_factors);
        size_t it = 0;
        T = NTL::ZZ(1);
        //get T
        for (const auto& tup : factors) {
            if (T/2 < tors_bound){
                it += 1;
                NTL::ZZ le = NTL::power(NTL::ZZ(std::get<0>(tup)), std::get<1>(tup));
                NTL::ZZ temp = NTL::GCD(T, le);
                T *= le/temp;
                max_cost = std::get<3>(tup);
            }
        }
        //T /= 2; //This will be taken care of somewhere else.
        if (it < factors.size()) {
            factors.erase(factors.begin() + it, factors.end());
        }

        k += 1;
        // search max prime ell we want to seach for with this k
        int maxellbits = 0;
        while (adj_cost(model(NTL::NextPrime(1 << maxellbits), 1, k), NTL::NextPrime(1 << maxellbits)) < max_cost) {
            maxellbits += 1;
        }
        int maxell = 0;
        for (int i = maxellbits; i >= 0; --i) {
            if (adj_cost(model(NTL::NextPrime(maxell | (1 << i)), 1, k), NTL::NextPrime(maxell | (1 << i))) < max_cost) {
                maxell |= 1 << i;
            }
        } 
        if (maxell < 1 || k > 80) {
            break;
        }

        //find small factors
        std::vector<std::tuple<int, int>> on_curve = smooth_part(NTL::power(p, k) - NTL::power(NTL::ZZ(-1), k), NTL::ZZ(maxell));
        //std::vector<std::tuple<int, int>> on_twist = smooth_part(NTL::power(p, k) + NTL::power(NTL::ZZ(-1), k), NTL::ZZ(maxell));
        
        for (const auto& tup : on_curve) {
            for (int i = 1; i <= std::get<1>(tup); ++i) {
                int j = in_vector(factors, std::get<0>(tup), i);
                if (NTL::GCD(NTL::ZZ(std::get<0>(tup)), coprime) > 1) {
                    continue;
                }
                if (j >= 0){
                    std::get<2>(factors[j]) = std::min(k, std::get<2>(factors[j]));
                    std::get<3>(factors[j]) = adj_cost(model(std::get<0>(tup), i, std::get<2>(factors[j])), std::get<0>(tup));
                } else{
                    ell_tuple newtup;
                    std::get<0>(newtup) = std::get<0>(tup);
                    std::get<1>(newtup) = i;
                    std::get<2>(newtup) = k;
                    std::get<3>(newtup) = adj_cost(model(std::get<0>(tup), i, k), std::get<0>(tup));
                    factors.push_back(newtup);
                }
            }
        }
        
        /* Skip on twist for now (twists not yet implemented):
        for (const auto& tup : on_twist) {
            for (int i = 1; i <= std::get<1>(tup); ++i) {
                int j = in_vector(factors, std::get<0>(tup), i);
                if (j >= 0){
                    std::get<2>(factors[j]) = std::min(k, std::get<2>(factors[j]));
                    std::get<3>(factors[j]) = adj_cost(model(std::get<0>(tup), i, std::get<2>(factors[j])), std::get<0>(tup));
                } else{
                    ell_tuple newtup;
                    std::get<0>(newtup) = std::get<0>(tup);
                    std::get<1>(newtup) = i;
                    std::get<2>(newtup) = k;
                    std::get<3>(newtup) = adj_cost(model(std::get<0>(tup), i, k), std::get<0>(tup));
                    factors.push_back(newtup);
                }
            }
        }*/
        
    }   
    
    //remove unnecessary factors if T is too large
    for (int i = factors.size() - 1; i >= 0; --i) {
        if (T/std::get<0>(factors[i]) > tors_bound) {
            T /= std::get<0>(factors[i]);
            factors.erase(factors.begin() + i);
        }
    }

    return factors;
}
