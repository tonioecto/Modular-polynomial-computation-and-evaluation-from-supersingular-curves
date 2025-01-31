///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements algorithms to perform ideal-to-isogeny translation  (header)
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <map>

#include <NTL/ZZ_pE.h>

#include "quaternions.hpp"
#include "isog.hpp"
#include "utils.hpp"
#include "Fp2k.hpp"

int torsionToFieldDegree(NTL::ZZ ell_e);

ecp evalEndo(quat const &alpha, ecp P, NTL::ZZ order_P);

std::vector<std::pair<ecp,std::pair<int, int>>> idealToKernel(quat const &alpha, std::unordered_map<int,int>& facN, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases);
isog_chain idealToIsogeny(quat const &alpha, NTL::ZZ const &N, ec const &E0, std::map<unsigned,Fp2k> &FieldExtensions, std::map<NTL::ZZ,std::pair<ecp,ecp>> &TorsionBases);
