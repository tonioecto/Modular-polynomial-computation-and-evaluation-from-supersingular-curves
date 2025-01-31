// Adapted from Sutherland's code for constructing Hilbert Class Polynomials (see packages on his website)

// All CRT cases treated as "small" in the sense of section 6 of "Computing Hilbert class polynomials with the Chinese Remainder Theorem")

#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <time.h>
#include "crt.hpp"
#include "modpoly_crt.hpp"


void modpoly_crt_precomp(modpoly_crt_struct poly_crt, std::vector<NTL::ZZ> primes, int Nprimes, int Ncoeffs, NTL::ZZ P){
    
    clock_t t;
    t = clock()

    poly_crt.ps = primes
    poly_crt.Nprimes = Nprimes;
    poly_crt.Ncoeffs = Ncoeffs;
    poly_crt.P = P
    
    crt_init(poly_crt.crt, primes, Nprimes, Ncoeffs, P);

    t = clock() - t;

    printf("CRT precomputation completed in %f seconds.\n",((float)t)/CLOCKS_PER_SEC);

    poly_crt.index = -1;
}

// void modpoly_crt_postcomp(){
    
// }