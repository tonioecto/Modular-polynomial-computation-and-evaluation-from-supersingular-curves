
///////////////////////////////////////////////////////////////////////////////////////////////
////   This header file constructs the structures needed to implement hashmap.cpp
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once
#include <optional>
#include <cassert>
#include <functional>

#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/LLL.h>

#include <gmp.h>
#include "quaternions.hpp"

#include <fplll/fplll.h>

#include "Fp2k.hpp"
#include "quatlatenum.hpp"
#include "utils.hpp"

typedef Fp_elem Fp;
typedef FpE_elem Fp2;
typedef FpX_elem FpX;
typedef FpEX_elem Fp2X; 

static const int LenNumberBytes = 8;

struct Key
{
    unsigned char IntList[3][LenNumberBytes];
};

struct KeyHash
{
    std::size_t operator()( const Key& k) const 
    {

        std::string str1(reinterpret_cast<char*>(const_cast<unsigned char*>(k.IntList[0])));
        std::string str2(reinterpret_cast<char*>(const_cast<unsigned char*>(k.IntList[1])));
        std::string str3(reinterpret_cast<char*>(const_cast<unsigned char*>(k.IntList[2])));

        return std::hash<std::string>()(str1) ^ 
            (std::hash<std::string>()(str2)) ^
            (std::hash<std::string>()(str3));
    }
};

struct KeyEqual
{
    bool operator()(const Key& lhs, const Key& rhs) const
    {
        NTL::ZZ lhs_IntList[3],rhs_IntList[3];
        for (int i=0; i<3 ; i++ ) {
            NTL::ZZFromBytes(lhs_IntList[i],lhs.IntList[i],LenNumberBytes);
            NTL::ZZFromBytes(rhs_IntList[i],rhs.IntList[i],LenNumberBytes);
        }
        return lhs_IntList[0] == rhs_IntList[0] && lhs_IntList[1] == rhs_IntList[1] && lhs_IntList[2] == rhs_IntList[2];
    }
};


void print_key(const Key& k);

inline bool is_key_equal(const Key& lhs, const Key& rhs)
    {
        NTL::ZZ lhs_IntList[3],rhs_IntList[3];
        for (int i=0; i<3 ; i++ ) {
            NTL::ZZFromBytes(lhs_IntList[i],lhs.IntList[i],LenNumberBytes);
            NTL::ZZFromBytes(rhs_IntList[i],rhs.IntList[i],LenNumberBytes);
        }
        return lhs_IntList[0] == rhs_IntList[0] && lhs_IntList[1] == rhs_IntList[1] && lhs_IntList[2] == rhs_IntList[2];
    }

// void print_order_jinv_map(std::unordered_map<Key, FpE_elem, KeyHash, KeyEqual> const &m);

// int map_contains(std::unordered_map<Key, FpE_elem, KeyHash, KeyEqual> const &m, const Key& k);

using mat_t = fplll::ZZ_mat<mpz_t>;
using gso_t = fplll::MatGSOGram<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>>;


// works correctly only if the lattice in input is an order
Key order_invariant_computation(quatlat const &order, quat *small);


// works correctly only if the lattice in input is an order
Key order_invariant_computation_from_gram(quatlat order, NTL::mat_ZZ gram, quat *small);


void inner_order_invariant_computation(NTL::ZZ *IntList, quatlat const &order, quat *small);

Integer get_smallest_element(quat* gamma, quatlat &I);

NTL::mat_ZZ compute_gram_order(quatlat const &order);