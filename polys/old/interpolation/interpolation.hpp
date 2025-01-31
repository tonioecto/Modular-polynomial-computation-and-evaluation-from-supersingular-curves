#include <NTL/ZZ_pX.h>
#include <iostream>
#include <vector>

std::vector<std::vector<NTL::ZZ_pX>> build_tree(std::vector<NTL::ZZ_pX> const m);
std::vector<NTL::ZZ_p> going_down_tree(NTL::ZZ_pX f, std::vector<NTL::ZZ_p> inputs, std::vector<std::vector<NTL::ZZ_pX>> tree);
NTL::ZZ_pX linear_combination(std::vector<NTL::ZZ_p> const inputs, std::vector<NTL::ZZ_p> const c, std::vector<std::vector<NTL::ZZ_pX>> const tree);
NTL::ZZ_pX fast_interpolate(std::vector<NTL::ZZ_p> const inputs, std::vector<NTL::ZZ_p> const outputs);
