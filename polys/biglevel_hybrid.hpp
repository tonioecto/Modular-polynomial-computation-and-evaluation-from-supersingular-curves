
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
// #include "crt.hpp"

// struct modpoly_crt{
//     std::vector<NTL::ZZ> ps; // primes 
//     int Nprimes;
//     int Ncoeffs;
//     NTL::ZZ P;
//     crt_info crt;
//     int index;
// };

// void modpoly_crt_precomp(modpoly_crt &poly_crt, std::vector<NTL::ZZ> const primes, int Nprimes, int Ncoeffs, NTL::ZZ const P);
void get_primes(std::vector<NTL::ZZ> &Pl, NTL::ZZ const B, NTL::ZZ const ell);
NTL::ZZ_pX jInv_interpolate(std::vector<NTL::ZZ_p> const js, std::vector<std::vector<NTL::ZZ_p>> const iso_js, std::vector<NTL::ZZ_p> const xs, NTL::ZZ ell);
void special_ss_eval(std::vector<NTL::ZZ> &Fq_coeffs, NTL::ZZ const q, std::vector<NTL::ZZ> const xs, NTL::ZZ const ell);
void biglevel(NTL::ZZ_pX &F, NTL::ZZ const p, NTL::ZZ_p const j, NTL::ZZ const ell);
void for_testing_js(std::vector<NTL::ZZX> &js_polys, std::vector<std::vector<NTL::ZZX>> &iso_js_polys, NTL::ZZ const q);
void for_testing_primes(std::vector<NTL::ZZ> &Pl);