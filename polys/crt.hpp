#pragma once

#include <iostream>
#include <vector>
#include <NTL/ZZ.h>


struct crt_info{
    std::vector<NTL::ZZ> m; // n moduli m_1,..., m_n, with M=prod m_i
    std::vector<NTL::ZZ> a; // n values a_i = 1/M_i mod m_i where M_i = M/m_i = prod_{j!=i} m_j
    NTL::ZZ P;							// output modulus
    NTL::ZZ MP;						// prod m_i mod P
    NTL::ZZ X,Y, Z;						// work variables
	int n;							// number of moduli
	int k;							// number of coefficients
	// int j;							// next coefficient to enumerate (after finalizing)
    std::vector<NTL::ZZ> Cdata; // Coefficients (TODO: Sutherland uses "limbs" for this?)
    std::vector<NTL::ZZ> sdata; // sj's in crt_update algo (TODO: Sutherland uses "limbs" for this?)
    int delta;  
};

void crt_coeff(std::vector<NTL::ZZ> &a, std::vector<NTL::ZZ> const m, int n);
void crt_init(crt_info &crt, std::vector<NTL::ZZ> const m, int n, int k, NTL::ZZ const P);
void crt_update(crt_info &crt, int i, std::vector<NTL::ZZ> const c, int k);
void crt_finalize_coeff(crt_info &crt, int j, NTL::ZZ tf);
void crt_finalise(crt_info &crt);