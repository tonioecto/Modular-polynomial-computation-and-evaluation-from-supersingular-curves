// Adapted from Sutherland's code for constructing Hilbert Class Polynomials (see packages on his website)

#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include <cassert>
#include "crt.hpp"

// TODO: use arrays if the vector has fixed length, otherwise use std::vector with push_back or resize
// std::vector<std::vector<int>> best way to create matrix?

/// Given n relatively prime moduli m[i], computes a[i] = \prod_{j!=i}m[j] mod m[i].
void crt_coeff(std::vector<NTL::ZZ> &a, std::vector<NTL::ZZ> const m, int n){

	std::size_t nbits = NTL::NumBits(n);
	std::vector<NTL::ZZ> e(n+nbits); // These "workspace values" in Sutherland code
	NTL::ZZ w0,w1,w2; // These "workspace values" in Sutherland code
	std::vector<std::vector<NTL::ZZ>> lm(nbits);
	std::vector<NTL::ZZ> lcnt(nbits); 
	int j,k;

	// PRINTING STUFF
	// std::cout << "crt_coeff inputs :" << std::endl;

	// for ( j = 0 ; j < n ; j++ ){
	// 	std::cout << m[j] << std::endl;
	// }

	// Base cases are the same as above
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
	// shift e by lcnt[k] -- is there a better way to do this?
	e.erase(e.begin(), e.begin() + NTL::conv<int>(lcnt[k]));

	// PRINTING STUFF
	// std::cout << "Level " << k << " nodes:" << std::endl;
	// std::cout << lcnt[k] << std::endl;

	// for ( j = 0 ; j < lcnt[k] ; j++ ){
	// 	std::cout << lcnt[k] << std::endl;
	// 	std::cout << lm[k][j] << std::endl;
	// }

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

		// PRINTING STUFF
		// std::cout << "Level " << k << " nodes:" << std::endl;
		// std::cout << lcnt[k] << std::endl;

		// for ( j = 0 ; j < lcnt[k] ; j++ ){
		// 	std::cout << lcnt[k] << std::endl;
		// 	std::cout << lm[k][j] << std::endl;
		// }
	}
	
	// PRINTING STUFF
		// std::cout << "lcnt[--k]: " << lcnt[--k] << std::endl;

	if ( lcnt[--k] != 2 ) { 
		std::cout << "Error, " << lcnt[k] << " nodes at top of product tree, expected exactly 2" << std::endl;
		abort();
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

// Agorithm 2.3 of Hilbert CRT paper by Sutherland
void crt_init(crt_info &crt, std::vector<NTL::ZZ> const m, int n, int k, NTL::ZZ const P){
	// TODO: add paralelisation as in Sutherland's code

	int i,j;
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

// Algorithm 2.4 of Hilbert CRT paper
// TODO: support batching
void crt_update(crt_info &crt, int i, std::vector<NTL::ZZ> const c, int k){
	
	int n;
	NTL::ZZ a;
	NTL::ZZ m;

	assert(k == crt.k);
	// assert(crt.j == -1); // otherwise call to crt_update while coefficient enumeration is in progress
	assert(i <= crt.n); // otherwise call to crt_update with modulus index > modulus count
	
	a = crt.a[i];
	m = crt.m[i];

	// compute d[i] = a_iM_i = a(M/m) mod P
	// compute M_i mod P using division mod P
	NTL::InvMod(crt.X, m, crt.P);
	NTL::MulMod(crt.X, crt.X, crt.MP, crt.P); 
	// Y = d[i] = a_iM_i mod P
	NTL::MulMod(crt.Y, crt.X, a, crt.P);
	
	// TODO: batching stuff
	
	for(int j = 0; j < k; j++){
		NTL::mul(crt.X, crt.Y, c[j]); // don't reduce mod P to save time
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

// Subroutine of Algorithm 2.5 of Hilbert CRT paper
void crt_finalize_coeff(crt_info &crt, int j, NTL::ZZ tf){

	std::cout << "sj: " << crt.sdata[j] << std::endl;
	NTL::add(crt.X, crt.sdata[j], tf); // add tf to sj
	std::cout << "add tf to sj: " << crt.X << std::endl;
	NTL::power2(crt.Y, crt.delta); // compute Y = 2^delta
	std::cout << "compute Y = 2^delta: " << crt.Y << std::endl;
	NTL::div(crt.X, crt.X, crt.Y); // divide by Y = 2^delta
	std::cout << "divide X by Y = 2^delta: " << crt.X << std::endl;
	NTL::MulMod(crt.X, crt.X, crt.MP, crt.P); // multiply by M
	std::cout << "multiply by MP: " << crt.X << std::endl;
	std::cout << "crt.Cdata[j] before: " << crt.Cdata[j] << std::endl;
	NTL::sub(crt.Cdata[j], crt.Cdata[j], crt.X);
	crt.Cdata[j] %= crt.P;
	std::cout << "crt.Cdata[j] after subtracting by " << crt.X << "mod P: " << crt.Cdata[j] << std::endl;
	std::cout << "" << std::endl;
}

// Algorithm 2.5 of Hilbert CRT paper
void crt_finalise(crt_info &crt){

	int j;

	// tf = (2^delta*3)/4
	
	NTL::ZZ tf;
	tf = NTL::ZZ(3) << (crt.delta-2);

	//TODO: add batching
	for(j = 0; j < crt.k; j++){
		crt_finalize_coeff(crt, j, tf);
	}
	// crt.j = 0;
}