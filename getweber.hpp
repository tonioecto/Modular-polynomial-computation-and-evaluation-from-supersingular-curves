
///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements the method to enable the Weber variants (header file)
///////////////////////////////////////////////////////////////////////////////////////////////


#pragma once

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "choosetorsion.hpp"
#include "hashmap.hpp"
#include "smallint.hpp"
#include "getresultants.hpp"

typedef std::pair<Integer,Integer> IntegerPair;

typedef std::vector<std::pair< Fp2, std::vector<std::vector<SmallIntegerPair>>>> weber_enum;


struct weber_bas
{
    ecp P16;
    ecp Q16;
    ecp P3;
    ecp Q3;
};

struct weber_full_data
{
    weber_bas basis;
    weber_enum enumerator;
};


struct weber_enum_poly_precomp {

    ZZ_pEXY R1;
    ZZ_pEXYZ R2;
    FpEX_elem G0;
    FpEX_elem H0;
    FpEX_elem I0;
    FpEX_elem J0;
    ZZ_pEXY Xs16;
    ZZ_pEXY D1;
    ZZ_pEXY D2;
    ZZ_pEXY N1;
    ZZ_pEXY N2;

};

weber_enum_poly_precomp SetWeberPrecomp();




Fp2 WeberGetFromEnum(const std::vector<std::vector<SmallIntegerPair>> &coeffs, const weber_enum &webdat, const std::pair<SmallMatFp, SmallMatFp> change_mats);
std::vector<Fp2> get_powers(Fp2 a, unsigned k);

FpE_elem CommonRootTwoResultants(std::vector<FpEX_elem> rs);
FpE_elem GetCommonRoot(std::vector<FpE_elem> cs);
FpE_elem _getWeberThirdPowerFromRoot(FpE_elem x, FpE_elem y);
FpE_elem _getWeberThirdPower(bool *check, std::vector<FpE_elem> cs);
FpE_elem _getWeberEighthPower(FpE_elem gamma2, FpE_elem t_inv);
std::vector<std::vector<ecp>> GetLevelStructureFromWeber(ec E, FpE_elem w, const std::map<unsigned,Fp2k> &Fexts);
bool GetWeberOfLevelStruct(Fp2 *w, const ec E, const std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts);
bool GetWeberOfLevelStruct_j0(Fp2 *w, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts);
bool GetWeberOfImage(Fp2 *w, ec E, isog phi, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts);
bool GetWeberOfImage_chain(Fp2* w, isog_chain phi, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts, bool retry = false);
FpE_elem GetWeberDomain(FpE_elem j);
NTL::ZZ_pE GetWeberDomainBig(NTL::ZZ_pE j);
//Assures the result is defined over Fp
FpE_elem GetWeberDomainFp(FpE_elem j);
NTL::ZZ_pE GetWeberDomainFpBig(NTL::ZZ_pE j);
std::vector<FpE_elem> GetWeberDomainAll(FpE_elem const &j);
std::vector<std::vector<ecp>> BasCoeffToLevelStructure(const weber_bas &web, const std::vector<std::vector<SmallIntegerPair>> coeff);
bool BasCoeffToWeber(Fp2 *w, const weber_bas &web, const std::vector<std::vector<SmallIntegerPair>> coeff, const std::map<unsigned,Fp2k> &Fexts);

weber_enum EnumerateAllWeberFast(const weber_bas &web, const std::map<unsigned,Fp2k> &Fexts, const weber_enum_poly_precomp *precomp);
Fp2 WeberGetFromEnum(const std::vector<std::vector<SmallIntegerPair>> &coeffs, const weber_enum &webdat, const std::pair<NTL::mat_ZZ_p, NTL::mat_ZZ_p> change_mats);


Fp2X eval_phi11_weber( Fp2 w );

inline FpE_elem w_to_j(FpE_elem const &w) {
    FpE_elem w_24 = NTL::power(w, 24);
    return NTL::power(w_24 - FpE_elem(16), 3)/w_24;
}

inline NTL::ZZ_pE w_to_j_BIG(NTL::ZZ_pE const &w) {
    NTL::ZZ_pE w_24 = NTL::power(w, 24);
    return NTL::power(w_24 - NTL::ZZ_pE(16), 3)/w_24;
}

inline NTL::ZZ_p w_to_j_prime(NTL::ZZ_p const &w) {
    NTL::ZZ_p w_24 = NTL::power(w, 24);
    return NTL::power(w_24 - NTL::ZZ_p(16), 3)/w_24;
}
