
#include <NTL/ZZ.h>
#include "crt.hpp"

struct modpoly_crt_struct{
    std::vector<NTL::ZZ> ps; // primes 
    int Nprimes;
    int Ncoeffs;
    NTL::ZZ P;
    crt_context crt;
    int index;
};

void modpoly_crt_precomp(modpoly_crt_struct poly_crt, std::vector<NTL::ZZ> primes, int Nprimes, int Ncoeffs, NTL::ZZ P);

