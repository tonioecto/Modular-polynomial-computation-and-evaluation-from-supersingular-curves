
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cmath> 
#include "bivariates.hpp" 


std::vector<std::vector<ZZ_pEXY>> BuildTreeBivar(std::vector<ZZ_pEXY> const m, unsigned k, size_t pow2k);
ZZ_pEXY FastInterpolateBivar(std::vector<NTL::ZZ_pEX> const inputs, std::vector<NTL::ZZ_pEX> const outputs);