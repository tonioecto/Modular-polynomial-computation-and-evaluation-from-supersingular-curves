///////////////////////////////////////////////
////   General utility functions
///////////////////////////////////////////////


#pragma once

#include <cassert>
#include <vector>
#include <unordered_map>

#include <NTL/ZZ.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pEX.h>

#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"

#include <time.h>
#include <locale.h>

inline uint64_t
cpucycles(void)
{
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (uint64_t)time.tv_sec * 1000000000 + time.tv_nsec;

}


__inline__ uint64_t
tic(void)
{
    return (uint64_t)cpucycles();
}


typedef NTL::ZZ Integer; //NTL::ZZ is probably overkill, optimise later

NTL::ZZ DLP(ecp const &Q, ecp const &base, long ell, long e);
std::unordered_map<int, int> factor(Integer const &N);
bool sutherland(long int const &jj, Fp_integer const &p); // For supersingularity testing
bool sutherland_slow(long int const &jj, Fp_integer const &p); // For supersingularity testing
bool init_sieve(const std::string& filename, std::vector<long int>& sieve_primes);
std::vector<long int> sieve_interval(long int &L, long int len_sieve_interval, std::vector<long int> sieve_primes);


// size_t DLP(ecp Q, ecp base, int ell) {return DLP(Q, base, ell, 1);};

// Product tree: inspired from SCALLOP code (https://github.com/isogeny-scallop/scallop)
template <class I>
I product_tree(std::vector<I> const &leaves) {
            if (leaves.empty())
                throw std::logic_error("no leaves");
            auto prev = leaves; // <- copies leaves to not modify the original list
            while (prev.size() > 1) {
                std::vector<I> next;
                {
                    for (size_t i = 0; i < prev.size()-1; i += 2)
                        next.push_back(prev[i] * prev[i+1]);
                    if (prev.size() % 2)
                        next.push_back(prev.back());
                }
                prev = next;
            }
            return prev[0];
        };


