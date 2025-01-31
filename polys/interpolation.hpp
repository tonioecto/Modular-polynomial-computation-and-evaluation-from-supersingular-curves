#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>

std::vector<std::vector<NTL::ZZ_pEX>> build_tree(std::vector<NTL::ZZ_pEX> const m, int k, size_t const pow2k);
std::vector<NTL::ZZ_pE> going_down_tree(NTL::ZZ_pEX f, std::vector<NTL::ZZ_pE> inputs, std::vector<std::vector<NTL::ZZ_pEX>> tree, int k, size_t const pow2k);
NTL::ZZ_pEX linear_combination(std::vector<NTL::ZZ_pE> const inputs, std::vector<NTL::ZZ_pE> const c, std::vector<std::vector<NTL::ZZ_pEX>> const tree, int k, size_t const pow2k);
NTL::ZZ_pEX fast_interpolate(std::vector<NTL::ZZ_pE> const inputs, std::vector<NTL::ZZ_pE> const outputs);
NTL::ZZ_pEX fast_interpolate_from_roots(std::vector<NTL::ZZ_pE> const roots);
