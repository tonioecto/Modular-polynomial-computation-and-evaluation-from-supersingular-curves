///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing the CRT method
////   Adapted from Sutherland's software package classpoly:
//// 							https://math.mit.edu/~drew/classpoly.html
////   which implements the algorithms in the papers:
////
////   Andrew V. Sutherland, Computing Hilbert class polynomials with the
////   Chinese Remainder Theorem, Mathematics of Computation 80 (2011), 501-538.
////
////   Andreas Enge and Andrew V. Sutherland, Class invariants by the CRT method,
////   Ninth Algorithmic Number Theory Symposium (ANTS IX), Lecture Notes in
////   Computer Science 6197 (2010), 142-156.
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include <cassert>

#include "crt.hpp"



void crt_coeff(std::vector<NTL::ZZ> &a, std::vector<NTL::ZZ> const m, int n){
	//////////////////////////////////////////////////////////////////////////////////////
	/// Given n relatively prime moduli m[i], computes a[i] = \prod_{j!=i}m[j] mod m[i].
	//////////////////////////////////////////////////////////////////////////////////////

	std::size_t nbits = NTL::NumBits(n);
	std::vector<NTL::ZZ> e(n+nbits);
	NTL::ZZ w0,w1,w2;
	std::vector<std::vector<NTL::ZZ>> lm(nbits);
	std::vector<NTL::ZZ> lcnt(nbits);
	int j,k;

	// Base cases
	if (n == 1){
		a[0] = NTL::ZZ(1);
		return;
	}
	if (n == 2){
		a[0] = m[1] % m[0];
		a[1] = m[0] % m[1];
		return;
	}
	
	// Build product tree level 1
	lcnt[0] = n; lm[1] = e; k = 1;
	for(j = 0; j+1 < lcnt[k-1]; j+= 2){
		NTL::mul(lm[k][j>>1], m[j], m[j+1]);
	}
	if(j < lcnt[k-1]){
		lm[k][j>>1] = m[j];
		j+= 2;
	}
	
	lcnt[k] = j>>1;
	// shift e by lcnt[k]
	e.erase(e.begin(), e.begin() + NTL::conv<int>(lcnt[k]));

	// Build product tree levels 2 and above
	for(k = 2; lcnt[k-1] > 2 ; k++) {
		lm[k] = e;
		for(j = 0; j+1 < lcnt[k-1]; j+= 2){
			NTL::mul(lm[k][j>>1], lm[k-1][j], lm[k-1][j+1]);
		}
		if(j < lcnt[k-1]){
			lm[k][j>>1] = lm[k-1][j];
			j+=2;
		}
		lcnt[k] = j>>1;
		e.erase(e.begin(), e.begin() + NTL::conv<int>(lcnt[k]));
	}
	
	if ( lcnt[--k] != 2 ) {
		std::cout << "Error, " << lcnt[k] << " nodes at top of product tree, expected exactly 2" << std::endl;
		throw;
	}

	// Set nodes at top level to contain reduced complements, we must have k>=1
	w0 = lm[k][1]%lm[k][0];
	lm[k][1] = lm[k][0]%lm[k][1];
	lm[k][0] = w0;

	

	// Reduce complements down the tree
	for(k--; k > 0; k--){
		for (j = 0; j+1 < lcnt[k]; j+=2){
			// Reduce parent complement
			w1 = lm[k+1][j>>1]%lm[k][j];							
			// Multiply right sibling into left and reduce (enlarges complement to include sibling)
			NTL::MulMod(w2, w1,lm[k][j+1],lm[k][j]);
			// Reduce parent complement
			w1 = lm[k+1][j>>1]%lm[k][j+1];						
			// Multiply left sibling into right and reduce
			NTL::MulMod(lm[k][j+1], w1,lm[k][j],lm[k][j+1]);
			
			lm[k][j] = w2;
		}
		// If no sibling, just copy parent value, its already reduced
		if(j < lcnt[k]){
			lm[k][j] = lm[k+1][j>>1];
		}	
	}
	
		// Reduce complements to level 0, which is the output a[]
		for (j = 0; j+1 < lcnt[k]; j+= 2){
			// Reduce parent complement
			w1 = lm[k+1][j>>1]%m[j];
			
			// Multiply right sibling into left and reduce (enlarges complement to include sibling)
			NTL::MulMod(w0, w1, m[j+1], m[j]);
			a[j] = w0;		
			// Reduce parent complement
			
			w1 = lm[k+1][j>>1]%m[j+1];
			// Multiply left sibling into right and reduce (ditto)
			NTL::MulMod(w0,w1,m[j], m[j+1]);
			a[j+1] = w0;	
		}
		if(j < lcnt[k]){
			a[j] = lm[k+1][j>>1];
		}
		
}

void crt_init(crt_info &crt, std::vector<NTL::ZZ> const m, int n, int k, NTL::ZZ const P){
	/////////////////////////////////////////////////////
	// Algorithm 2.3 of Hilbert CRT paper by Sutherland
	/////////////////////////////////////////////////////

	int i;
	NTL::ZZ X;
	// Compute a[i] = M_i mod m[i] where M_i = prod_{j\ne i} m_i
	
	crt.a.resize(n);
	crt_coeff(crt.a, m, n);
	
	// Compute M mod P iteratively, don't bother with a tree (unless P is huge, it wouldn't make much of a difference)
	X = m[0];
	for(i = 1; i < n; i++){
		NTL::MulMod(X, X, m[i], P);
	}

	crt.MP = X;
	
	crt.P = P;

	// _ecrt_init in Sutherland's code (minus parallelisation and "limbs" stuff)
	crt.n = n;
	crt.k = k;
	// crt.j = -1;


	crt.Cdata.resize(k);
	crt.sdata.resize(k);
	crt.delta = NTL::NumBits(n)+1;
	crt.m.resize(m.size());
	for(i = 0; i < n; i++){
		crt.m[i] = m[i];
	}

	// a holds 1/M_i mod m[i]
	for(i = 0; i < n; i++){
		NTL::InvMod(crt.a[i], crt.a[i], m[i]);	
	}
}

void crt_update(crt_info &crt, int i, std::vector<NTL::ZZ> const c, int k){
	////////////////////////////////////////////////
	// Algorithm 2.4 of Hilbert CRT paper

	////////////////////////////////////////////////
	NTL::ZZ a;
	NTL::ZZ m;

	assert(k == crt.k);
	assert(i <= crt.n);
	
	a = crt.a[i];
	m = crt.m[i];

	// compute d[i] = a_iM_i = a(M/m) mod P
	// compute M_i mod P using division mod P
	// TODO: add batched inversion
	NTL::InvMod(crt.X, m%crt.P, crt.P);
	NTL::MulMod(crt.X, crt.X, crt.MP, crt.P);
	// Y = d[i] = a_iM_i mod P
	NTL::MulMod(crt.Y, crt.X, a, crt.P);
	
	
	for(int j = 0; j < k; j++){
		NTL::mul(crt.X, crt.Y, c[j]);
		// Add to Cj
		NTL::add(crt.Cdata[j], crt.Cdata[j], crt.X);
		crt.Cdata[j] %= crt.P;
		
		// compute floor((2^delta*c[j]*a[i])/m[i])
		NTL::mul(crt.X, c[j], a);
		NTL::power2(crt.Z, crt.delta);
		NTL::mul(crt.X, crt.X, crt.Z);
		NTL::div(crt.X, crt.X, m);

		// Add to sj
		NTL::add(crt.sdata[j], crt.sdata[j], crt.X);
	}
	}

void crt_finalize_coeff(crt_info &crt, int j, NTL::ZZ tf){
	/////////////////////////////////////////////////////
	// Subroutine of Algorithm 2.5 of Hilbert CRT paper
	/////////////////////////////////////////////////////

	NTL::add(crt.X, crt.sdata[j], tf);
	NTL::power2(crt.Y, crt.delta);
	NTL::div(crt.X, crt.X, crt.Y);
	NTL::MulMod(crt.X, crt.X, crt.MP, crt.P);
	NTL::sub(crt.Cdata[j], crt.Cdata[j], crt.X);
	crt.Cdata[j] %= crt.P;
}

void crt_finalise(crt_info &crt){
	////////////////////////////////////////
	// Algorithm 2.5 of Hilbert CRT paper
	////////////////////////////////////////

	int j;
	NTL::ZZ tf;
	tf = NTL::ZZ(3) << (crt.delta-2);

	for(j = 0; j < crt.k; j++){
		crt_finalize_coeff(crt, j, tf);
	}
}
