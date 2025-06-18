///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////// The main function in this file is choose_torsion, which allows us to create a
/////// 'dictionary' of prime powers \ell^e and the field extension F_{p^{2k}} over which
///////  the ell^e-torsion is defined.
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <tuple>
#include <vector>

#include "costmodel.hpp"


//need this for using std::pair<int,int> as first element of unordered_map
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& pair) const {
        auto h1 = std::hash<T1>{}(pair.first);
        auto h2 = std::hash<T2>{}(pair.second);

        // combine hash values
        return h1 ^ h2;
    }
};


//store factors in following way:
//<<ell, e>, k>, where ell is prime, e the exponent, k the field extension
using ell_tuple = std::tuple<int, int, int, NTL::RR>;
//using factor_ext = std::pair<std::pair<int,int>, int>;
using factor_list = std::vector<ell_tuple>;


void is_prime_power(NTL::ZZ &ell, int &exp, NTL::ZZ const &N);
factor_list choose_torsion_naive(NTL::ZZ const &p, NTL::ZZ const &tors_bound, NTL::ZZ const &coprime = NTL::ZZ(1));
factor_list choose_torsion(NTL::ZZ const &p, NTL::ZZ const &tors_bound, NTL::ZZ const &coprime = NTL::ZZ(1));
bool sort_factors(ell_tuple const &tup1, ell_tuple const &tup2);
std::vector<int, int> smnooth_part(NTL::ZZ const &N, NTL::ZZ const &B);
int in_vector(factor_list const &vec, int const &ell, int const &exp);


