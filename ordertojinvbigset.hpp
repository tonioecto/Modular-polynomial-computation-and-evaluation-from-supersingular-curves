///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implements the algorithm OrdersTojInvariantBigSet from
////      https://eprint.iacr.org/2023/064
///////////////////////////////////////////////////////////////////////////////////////////////

#pragma once
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/mat_ZZ_p.h>

#include "hashmap.hpp"
#include "endring.hpp"
#include "quaternions.hpp"
#include "choosetorsion.hpp"
#include "klpt.hpp"
#include "Fp2k.hpp"
#include <list>
#include "getweber.hpp"

struct Jinv
{
    unsigned char IntList[2][LenNumberBytes];
};

struct JinvHash
{
    std::size_t operator()( const Jinv& k) const
    {

        std::string str1(reinterpret_cast<char*>(const_cast<unsigned char*>(k.IntList[0])));
        std::string str2(reinterpret_cast<char*>(const_cast<unsigned char*>(k.IntList[1])));

        return std::hash<std::string>()(str1) ^ (std::hash<std::string>()(str2));
    }
};

struct JinvEqual
{
    bool operator()(const Jinv& lhs, const Jinv& rhs) const
    {
        NTL::ZZ lhs_IntList[2],rhs_IntList[2];
        for (int i=0; i<2 ; i++ ) {
            NTL::ZZFromBytes(lhs_IntList[i],lhs.IntList[i],LenNumberBytes);
            NTL::ZZFromBytes(rhs_IntList[i],rhs.IntList[i],LenNumberBytes);
        }
        return lhs_IntList[0] == rhs_IntList[0] && lhs_IntList[1] == rhs_IntList[1];
    }
};




Jinv JToJinv(const FpE_elem &j);

bool is_Fp(const Fp2 &j);
Fp2 Frob(const Fp2 &j);

quat find_quaternion_iterator(std::list<int>& prime_list, const quatlat& I, const quatlat& O0, const quatalg &Bp);

void order_to_jinv_full_list(std::unordered_map<Key, std::pair<FpE_elem, quatlat>, KeyHash, KeyEqual> &m, std::vector<std::pair<quatlat, std::pair<quatlat, Key>>> &ideal_list, const NTL::ZZ &p, const quatalg &Bp, const std::map<unsigned,Fp2k> &Fexts, const Integer &coprime = Integer(1));

std::pair<weber_bas,std::vector<std::pair<SmallMatFp,SmallMatFp>>> order_to_weber_inv_full_list(std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<quatlat,quat>, weber_full_data>>, KeyHash, KeyEqual> &m, std::vector<std::pair<quatlat, std::pair<quatlat, Key>>> &ideal_list, const NTL::ZZ &p, const quatalg &Bp, const std::map<unsigned,Fp2k> &Fexts, const Integer &coprime = Integer(1));
