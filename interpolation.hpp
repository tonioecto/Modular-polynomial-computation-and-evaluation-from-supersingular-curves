/////////////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing interpolation of univariate polynomials over finite fields (header file)
/////////////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <iostream>
#include <vector>
#include "Fp2k.hpp"
#include "fast_ff.hpp"


std::vector<std::vector<FpEX_elem>> BuildTree(std::vector<FpEX_elem> const &m, unsigned k, size_t pow2k);
std::vector<FpE_elem> GoingDownTree(FpEX_elem f, std::vector<FpE_elem> inputs, std::vector<std::vector<FpEX_elem>> tree, int k, size_t pow2k);
FpEX_elem LinearCombo(std::vector<FpEX_elem> const c, std::vector<std::vector<FpEX_elem>> const tree, int k, size_t const pow2k);
FpEX_elem FastInterpolate(std::vector<FpE_elem> const inputs, std::vector<FpE_elem> const outputs);
// interpolate poly from its roots
// use native NTL native Fp2X multiplication which seems to be 6 times slower than FpX multiplication 
void FastInterpolateFromRoots(FpEX_elem &poly, std::vector<FpE_elem> const &roots);
// same as above
// performs the Fp2X multiplication from 3 FpX multiplication with Karatsuba (we know that Fp2 is Fp[X] / (X^2 + 1))
// due to memory management issues, this is more than 3 times slower than FpX multiplication, but it is still much better than the version based on NTL native Fp2X mul
void FastInterpolateFromRootsKaratsuba(std::pair<FpX_elem, FpX_elem> &poly, std::vector<ffp2> const & roots);

// same as above
// performs the Fp2X multiplication from 3 FpX multiplication with Karatsuba (we know that Fp2 is Fp[X] / (X^2 + 1))
// due to memory management issues, this is more than 3 times slower than FpX multiplication, but it is still much better than the version based on NTL native Fp2X mul
void FastInterpolateFromRootsKaratsubaPlusTrick(std::pair<FpX_elem, FpX_elem> &poly, std::vector<ffp2> const & roots);