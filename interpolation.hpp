/////////////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing interpolation of univariate polynomials over finite fields (header file)
/////////////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <iostream>
#include <vector>
#include "Fp2k.hpp"

std::vector<std::vector<FpEX_elem>> BuildTree(std::vector<FpEX_elem> const &m, unsigned k, size_t pow2k);
std::vector<FpE_elem> GoingDownTree(FpEX_elem f, std::vector<FpE_elem> inputs, std::vector<std::vector<FpEX_elem>> tree, int k, size_t pow2k);
FpEX_elem LinearCombo(std::vector<FpEX_elem> const c, std::vector<std::vector<FpEX_elem>> const tree, int k, size_t const pow2k);
FpEX_elem FastInterpolate(std::vector<FpE_elem> const inputs, std::vector<FpE_elem> const outputs);
FpEX_elem FastInterpolateFromRoots(std::vector<FpE_elem> const roots);