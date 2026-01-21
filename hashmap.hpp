
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
#include "fast_quaternions.hpp"

#include <fplll/fplll.h>

#include "Fp2k.hpp"
#include "quatlatenum.hpp"
#include "utils.hpp"

typedef Fp_elem Fp;
typedef FpE_elem Fp2;
typedef FpX_elem FpX;
typedef FpEX_elem Fp2X;

static const int LenNumberBytes = 4;

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

        return std::hash<std::string>()(str1) ^ (std::hash<std::string>()(str2)) ^ (std::hash<std::string>()(str3));
    }
};

struct KeyEqual
{
    bool operator()(const Key& lhs, const Key& rhs) const
    {
        return (!memcmp(lhs.IntList[0], rhs.IntList[0], LenNumberBytes)) && (!memcmp(lhs.IntList[1], rhs.IntList[1], LenNumberBytes)) && (!memcmp(lhs.IntList[2], rhs.IntList[2], LenNumberBytes));
    }
};


void print_key(const Key& k);

inline bool is_key_equal(const Key& lhs, const Key& rhs)
    {
        return (!memcmp(lhs.IntList[0], rhs.IntList[0], LenNumberBytes)) && (!memcmp(lhs.IntList[1], rhs.IntList[1], LenNumberBytes)) && (!memcmp(lhs.IntList[2], rhs.IntList[2], LenNumberBytes));
    }


using mat_t = fplll::ZZ_mat<mpz_t>;
using gso_t = fplll::MatGSOGram<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>>;


// works correctly only if the lattice in input is an order
Key order_invariant_computation(quatlat const &order, quat *small);

// works correctly only if the lattice in input is an order
Key order_invariant_computation_from_gram(quatlat order, NTL::mat_ZZ gram, std::pair<quat,quat> *small);




void inner_order_invariant_computation(NTL::ZZ *IntList, quatlat const &order, quat *small);

Integer get_smallest_element(quat* gamma, quatlat &I);

void GreedyReduction3(NTL::mat_ZZ &gram);

void FastGreedyReduction3(FastMat3 &gram, FastMat3 &Coords);

Key fast_order_invariant_computation_from_gram(const FastQuatLat &order, FastMat3 &gram, std::pair<FastQuat,FastQuat> &small);