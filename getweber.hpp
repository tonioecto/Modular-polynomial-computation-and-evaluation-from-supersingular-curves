
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
#include "mont.hpp"

typedef std::pair<Integer,Integer> IntegerPair;

typedef std::array<std::pair< ffp2, std::vector<std::vector<VerySmallIntegerPair>>>, 3> weber_enum;
typedef std::array<ffp2, 72> weber_inv_list;


struct weber_bas
{
    ecp P16;
    ecp Q16;
    ecp P3;
    ecp Q3;
};

struct weber_full_data
{
    bool check;
    weber_inv_list inv_list;
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

struct fast_weber_enum_poly_precomp {

    ffp2XY R1;
    ffp2XYZ R2;
    FpEX_elem G0;
    FpEX_elem H0;
    FpEX_elem I0;
    FpEX_elem J0;
    ffp2XY Xs16;
    ffp2XY D1;
    ffp2XY D2;
    ffp2XY N1;
    ffp2XY N2;
};

fast_weber_enum_poly_precomp Set(const weber_enum_poly_precomp *t);





void get_powers(std::vector<Fp2> &powers,const Fp2 &a, unsigned k);


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

// weber_full_data EnumerateAllWeberFast(const weber_bas &web, const std::map<unsigned,Fp2k> &Fexts, const weber_enum_poly_precomp *precomp);
weber_full_data EnumerateAllWeberFastFast(const weber_bas &web, const std::map<unsigned,Fp2k> &Fexts, const fast_weber_enum_poly_precomp *precomp);

std::array<std::vector<std::vector<VerySmallIntegerPair>>, 72>  EnumerateAllWeberCoeff();

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

unsigned char WeberGetFromEnum(const std::vector<std::vector<VerySmallIntegerPair>> &coeffs, const std::pair<VerySmallMat, VerySmallMat> &change_mats);

std::pair<VerySmallMat,VerySmallMat> WeberBasApplicationRemoteEndo(const std::vector<std::pair<VerySmallMat,VerySmallMat>> &mat0, FastQuat &gamma);