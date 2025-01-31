///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements the method described in Secion 3.5 of the accompanying paper to 
////   enable the Weber variants of the ModEvalBigLevel and ModEvalBigChar algorithms
////   (in modpolys.cpp)
////
////   We use ideas from  https://github.com/mariascrs/SplitSearcher
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <iostream>
#include <vector>
#include "getresultants.hpp"
#include "ecp.hpp"
#include "isog.hpp"
#include "id2iso.hpp"
#include "hashmap.hpp"
#include <unordered_set>
#include "getweber.hpp"

///////////////////////////////////////////////////////////////////////////////
// Implementing elimination theory algos needed to obtain the weber invariant
//   Inspired by https://github.com/mariascrs/SplitSearcher
///////////////////////////////////////////////////////////////////////////////


std::vector<FpE_elem> get_powers(FpE_elem a, unsigned k){
    
    std::vector<FpE_elem> powers(k+1);

    FpE_elem A = FpE_elem(1);
    powers[0] = A;

    for(unsigned i = 1; i <= k; i++){
        A = A*a;
        powers[i] = A;
    }

    return powers;
}

std::vector<FpEX_elem> EvaluateResultants(ZZ_pEXY R1, ZZ_pEXYZ R2, std::vector<FpE_elem> cs)
{

    FpE_elem c0 = cs[0];
    FpE_elem c1 = cs[1];

    std::vector<FpEX_elem> rs(2);

    rs[0] = EvaluateBivariate(R1,c0);
    rs[1] = EvaluateTrivariate(R2,c0,c1);

    return rs;
}


FpE_elem CommonRootTwoResultants(std::vector<FpEX_elem> rs)
{
    FpEX_elem r0 = rs[0];
    FpEX_elem r1 = rs[1];

    // TODO: implement fast GCDs?
    FpEX_elem g = GCD(r0,r1);

    return FpE_elem(-g[0]);
}

FpE_elem GetCommonRoot(std::vector<FpE_elem> cs){
    
    auto R1 = GetBivariateResultant();
    auto R2 = GetTrivariateResultant();

    auto rs = EvaluateResultants(R1, R2, cs);
    // std::cout << rs[0] << "\n" << rs[1] << "\n";
    return CommonRootTwoResultants(rs);   


};

////////////////////////////////////////////////////////
// Algorithms to get powers of the Weber invariant
////////////////////////////////////////////////////////


FpE_elem _getWeberThirdPowerFromRoot(FpE_elem x, FpE_elem y){

    FpEX_elem n1,n2,n3,n4,n5;
    NTL::SetCoeff(n1,2,FpE_elem(1));
    NTL::SetCoeff(n2,6,FpE_elem(-1));
    NTL::SetCoeff(n3,2,FpE_elem(20));
    NTL::SetCoeff(n4,6,FpE_elem(-16));
    NTL::SetCoeff(n5,2,FpE_elem(64));
    
    NTL::Vec<FpEX_elem> cs_N1;
    cs_N1.SetLength(17);
    cs_N1[0] = n1;
    cs_N1[4] = n2;
    cs_N1[8] = n3;
    cs_N1[12] = n4;
    cs_N1[16] = n5;

    NTL::Vec<FpEX_elem> cs_D1;
    cs_D1.SetLength(16);
    cs_D1[7] = FpEX_elem(4);
    cs_D1[15] = FpEX_elem(32);

    ZZ_pEXY N1 = ZZ_pEXY(cs_N1);
    ZZ_pEXY D1 = ZZ_pEXY(cs_D1);

    auto N1_eval = EvaluateBivariate(N1, y);
    auto D1_eval = EvaluateBivariate(D1, y);

    return NTL::eval(N1_eval,x)/NTL::eval(D1_eval,x);

};

FpE_elem _getWeberThirdPower(bool *check, std::vector<FpE_elem> cs){

    auto y = GetCommonRoot(cs);

    // Evaluate Xs16p at Y = y (TODO: more efficient at something other than Xs16p?)

    FpEX_elem f1,f2,f3,f4;
    NTL::SetCoeff(f1,8,FpE_elem(16));
    NTL::SetCoeff(f1,0,FpE_elem(-16));
    NTL::SetCoeff(f2,4,FpE_elem(-16));
    NTL::SetCoeff(f3,8,FpE_elem(1));
    NTL::SetCoeff(f4,4,FpE_elem(-1));

    NTL::Vec<FpEX_elem> cs_Xs16;
    cs_Xs16.SetLength(13);
    cs_Xs16[0] = f4;
    cs_Xs16[4] = f3;
    cs_Xs16[8] = f2;
    cs_Xs16[12] = f1;

    ZZ_pEXY Xs16 = ZZ_pEXY(cs_Xs16);

    FpEX_elem n11,n12,n13,n14,n15,n16,n17,n18;
    SetCoeff(n11,4,FpE_elem(1));
    SetCoeff(n12,5,FpE_elem(-1));
    SetCoeff(n13,6,FpE_elem(1));
    SetCoeff(n14,7,FpE_elem(-1));
    SetCoeff(n15,4,FpE_elem(16));
    SetCoeff(n16,5,FpE_elem(-16));
    SetCoeff(n17,6,FpE_elem(16));
    SetCoeff(n18,7,FpE_elem(-16));

    FpEX_elem n21,n22,n23,n24,n25,n26,n27,n28,n29;
    SetCoeff(n21,3,FpE_elem(-1));
    SetCoeff(n22,7,FpE_elem(1));
    SetCoeff(n22,3,FpE_elem(14));
    SetCoeff(n22,2,FpE_elem(-4));
    SetCoeff(n23,1,FpE_elem(-16));
    SetCoeff(n23,2,FpE_elem(48));
    SetCoeff(n23,3,FpE_elem(-120));
    SetCoeff(n23,6,FpE_elem(4));
    SetCoeff(n23,7,FpE_elem(-14));
    SetCoeff(n24,0,FpE_elem(-64));
    SetCoeff(n24,1,FpE_elem(160));
    SetCoeff(n24,2,FpE_elem(-368));
    SetCoeff(n24,3,FpE_elem(728));
    SetCoeff(n24,5,FpE_elem(16));
    SetCoeff(n24,6,FpE_elem(-48));
    SetCoeff(n24,7,FpE_elem(120));
    SetCoeff(n25,7,FpE_elem(-736));
    SetCoeff(n25,6,FpE_elem(384));
    SetCoeff(n25,5,FpE_elem(-192));
    SetCoeff(n25,4,FpE_elem(128));
    SetCoeff(n25,3,FpE_elem(-3328));
    SetCoeff(n25,2,FpE_elem(1984));
    SetCoeff(n25,1,FpE_elem(-1088));
    SetCoeff(n25,0,FpE_elem(512));
    SetCoeff(n26,7,FpE_elem(3200));
    SetCoeff(n26,6,FpE_elem(-1792));
    SetCoeff(n26,5,FpE_elem(768));
    SetCoeff(n26,3,FpE_elem(11648));
    SetCoeff(n26,2,FpE_elem(-7936));
    SetCoeff(n26,1,FpE_elem(5120));
    SetCoeff(n26,0,FpE_elem(-3072));
    SetCoeff(n27,7,FpE_elem(-12288));
    SetCoeff(n27,6,FpE_elem(9216));
    SetCoeff(n27,5,FpE_elem(-7168));
    SetCoeff(n27,4,FpE_elem(6144));
    SetCoeff(n27,3,FpE_elem(-30720));
    SetCoeff(n27,2,FpE_elem(23552));
    SetCoeff(n27,1,FpE_elem(-17408));
    SetCoeff(n27,0,FpE_elem(12288));
    SetCoeff(n28,7,FpE_elem(24576));
    SetCoeff(n28,6,FpE_elem(-16384));
    SetCoeff(n28,5,FpE_elem(8192));
    SetCoeff(n28,3,FpE_elem(57344));
    SetCoeff(n28,2,FpE_elem(-49152));
    SetCoeff(n28,1,FpE_elem(40960));
    SetCoeff(n28,0,FpE_elem(-32768));
    SetCoeff(n29,7,FpE_elem(-65536));
    SetCoeff(n29,6,FpE_elem(65536));
    SetCoeff(n29,5,FpE_elem(-65536));
    SetCoeff(n29,4,FpE_elem(65536));
    SetCoeff(n29,3,FpE_elem(-65536));
    SetCoeff(n29,2,FpE_elem(65536));
    SetCoeff(n29,1,FpE_elem(-65536));
    SetCoeff(n29,0,FpE_elem(65536));

    NTL::Vec<FpEX_elem> cs_N1, cs_N2;
    cs_N1.SetLength(12);
    cs_N1[0] = n11;
    cs_N1[1] = n12;
    cs_N1[2] = n13;
    cs_N1[3] = n14;
    cs_N1[8] = n15;
    cs_N1[9] = n16;
    cs_N1[10] = n17;
    cs_N1[11] = n18;

    cs_N2.SetLength(33);
    cs_N2[0] = n21;
    cs_N2[4] = n22;
    cs_N2[8] = n23;
    cs_N2[12] = n24;
    cs_N2[16] = n25;
    cs_N2[20] = n26;
    cs_N2[24] = n27;
    cs_N2[28] = n28;
    cs_N2[32] = n29;


    NTL::Vec<FpEX_elem> cs_D1, cs_D2;
    cs_D1.SetLength(13);
    cs_D1[12] = FpEX_elem(4);
    cs_D2.SetLength(17);
    cs_D2[16] = FpEX_elem(8);
    cs_D2[8] = FpEX_elem(1);

    ZZ_pEXY D1 = ZZ_pEXY(cs_D1);
    ZZ_pEXY D2 = ZZ_pEXY(cs_D2);
    ZZ_pEXY N1 = ZZ_pEXY(cs_N1);
    ZZ_pEXY N2 = ZZ_pEXY(cs_N2);

    auto d1 = EvaluateBivariate(D1,y);
    auto d2 = EvaluateBivariate(D2,y);
    auto n1 = EvaluateBivariate(N1,y);
    auto n2 = EvaluateBivariate(N2,y);

    FpEX_elem F1 = n1 - cs[0]*d1;
    FpEX_elem F2 = n2 - cs[1]*d2;

    Fp2X g = GCD(F1,F2);
        // assert(NTL::deg(g) > 0);
    if (NTL::deg(g) == 0) {
        (*check) = false;
        FpE_elem a;
        NTL::conv(a, 0);
        return a;
    }
        
    (*check) = true;
    // seems to be true most of the time but not always
    if (NTL::deg(g) == 1) {
        return _getWeberThirdPowerFromRoot(-g[0],y);
    }
    else {
            Fp2X F = EvaluateBivariate(Xs16,y);
            F /= LeadCoeff(F);
            Fp2X gg = GCD(g,F);
            if (NTL::deg(gg) != 1) {
            }
            assert(NTL::deg(gg)==1);
             return _getWeberThirdPowerFromRoot(-gg[0],y);
    }

};



FpE_elem _getWeberEighthPower(FpE_elem gamma2, FpE_elem t_inv){

    FpEX_elem f1, f2;

    NTL::SetCoeff(f1, 3, FpE_elem(1));
    NTL::SetCoeff(f1, 0, FpE_elem(t_inv));
    NTL::SetCoeff(f2, 3, FpE_elem(1));
    NTL::SetCoeff(f2, 1, FpE_elem(-gamma2));
    NTL::SetCoeff(f2, 0, FpE_elem(-16));

    // std::cout << "polys:" << std::endl;
    // std::cout << f1 << std::endl;
    // std::cout << f2 << std::endl;
    auto g = GCD(f1,f2);
    // std::cout << g << std::endl;
    
    assert(NTL::deg(g) == 1);

    return FpE_elem(-g[0]); 

};



FpE_elem _getGammaTwoInvariant(std::vector<ecp> tors3, ec E){
    
    assert((3*tors3[0]).is_identity());
    assert((3*tors3[1]).is_identity());
    assert((3*tors3[2]).is_identity());
    assert((3*tors3[3]).is_identity());
    // Ordering should be fixed from the first step
    std::vector<FpE_elem> tors3_image = {tors3[0].aff_x(), 
                                           tors3[1].aff_x(),
                                           tors3[2].aff_x(), 
                                           tors3[3].aff_x()};



    assert((tors3_image[0] != tors3_image[1]) && (tors3_image[0] != tors3_image[3]) && (tors3_image[0] != tors3_image[2]) && (tors3_image[2] != tors3_image[1]) && (tors3_image[3] != tors3_image[1]) && (tors3_image[3] != tors3_image[2]));                                       

    Fp2 A_image = (E).a();
    Fp2 c3 = 2*A_image;
    Fp2 c4 = -24*c3;

    Fp2 c;
    
    if (NTL::deg(rep(tors3_image[0])) > 1 || NTL::deg(rep(tors3_image[1])) > 1 || NTL::deg(rep(tors3_image[2])) > 1 || NTL::deg(rep(tors3_image[2])) > 1) {
        // the extension degree might be 4
        // assert(NTL::deg(rep(tors3_image[0])) == 3);
        Fp2 c_FpE;
        {
            FpE_push push(tors3[0].field().F);
            c_FpE = tors3_image[0]*tors3_image[1] + tors3_image[2]*tors3_image[3];
        }

        c = coerce(c_FpE, tors3[0].field());
    }
    else {
        c = tors3_image[0]*tors3_image[1] + tors3_image[2]*tors3_image[3];
    }

    

    return c4/(c3 - 3*c);

};



FpE_elem _getTInvariant(ecp P2, ec E){
    // P2 is a two torsion point 
    assert(!P2.is_identity());
    assert((2*(P2)).is_identity());
    isog phi2 = isogeny(P2, 2);

    auto E2 = phi2.get_codomain();

    return E2.disc()/E.disc();
};


////////////////////////////////////////////////////////////////////////////////
// For fetching the polynomials needed (as described in the accompanying paper)
////////////////////////////////////////////////////////////////////////////////


FpEX_elem _getI08(){
    FpEX_elem I0;

    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("4096"), NTL::conv<NTL::ZZ>("98304"), NTL::conv<NTL::ZZ>("847872"), NTL::conv<NTL::ZZ>("3092480"), NTL::conv<NTL::ZZ>("4436736"), NTL::conv<NTL::ZZ>("3379200"), NTL::conv<NTL::ZZ>("1564160"), NTL::conv<NTL::ZZ>("468480"), NTL::conv<NTL::ZZ>("92976"), NTL::conv<NTL::ZZ>("12160"), NTL::conv<NTL::ZZ>("1008"), NTL::conv<NTL::ZZ>("48"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(I0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return I0;
};

FpEX_elem _getJ08(){
    FpEX_elem J0;
    
    std::vector<NTL::ZZ> cs = {NTL::conv<NTL::ZZ>("128"), NTL::conv<NTL::ZZ>("80"), NTL::conv<NTL::ZZ>("16"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {1, 2, 3, 4};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(J0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return J0;
};

FpEX_elem _getG04(){
    FpEX_elem G0;

    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("68719476736"), NTL::conv<NTL::ZZ>("12884901888"), NTL::conv<NTL::ZZ>("855638016"), NTL::conv<NTL::ZZ>("23068672"), NTL::conv<NTL::ZZ>("208896"), NTL::conv<NTL::ZZ>("768"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3, 4, 5, 6};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(G0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return G0;
};

FpEX_elem _getH04(){
    FpEX_elem H0;
    
    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("16"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {4, 5};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(H0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return H0;
};

FpEX_elem _getI04(){
    FpEX_elem I0;

    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("4096"), NTL::conv<NTL::ZZ>("12288"), NTL::conv<NTL::ZZ>("13056"), NTL::conv<NTL::ZZ>("5632"), NTL::conv<NTL::ZZ>("816"), NTL::conv<NTL::ZZ>("48"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3, 4, 5, 6};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(I0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return I0;
};

FpEX_elem _getJ04(){
    FpEX_elem J0;
    
    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("16"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {1, 2};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(J0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return J0;
};

FpEX_elem _getI02(){
    FpEX_elem I0;

    std::vector<NTL::ZZ> cs = {  NTL::conv<NTL::ZZ>("4096"), NTL::conv<NTL::ZZ>("768"), NTL::conv<NTL::ZZ>("48"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(I0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return I0;
};

FpEX_elem _getJ02(){
    FpEX_elem J0;
    
    std::vector<NTL::ZZ> cs = { NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {1};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(J0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return J0;
};




FpEX_elem _getG0(){
    FpEX_elem G0;

    std::vector<NTL::ZZ> cs = {NTL::conv<NTL::ZZ>("1152921504606846976"), NTL::conv<NTL::ZZ>("6917529027641081856"), NTL::conv<NTL::ZZ>("19887895954468110336"), NTL::conv<NTL::ZZ>("36461142583191535616"), NTL::conv<NTL::ZZ>("47841738841556779008"), NTL::conv<NTL::ZZ>("47787695646028333056"), NTL::conv<NTL::ZZ>("37724965228622381056"), NTL::conv<NTL::ZZ>("24114242729778610176"), NTL::conv<NTL::ZZ>("12682896313210109952"), NTL::conv<NTL::ZZ>("5546125766502121472"), NTL::conv<NTL::ZZ>("2028993128165277696"), NTL::conv<NTL::ZZ>("622699889175822336"), NTL::conv<NTL::ZZ>("160256972254347264"), NTL::conv<NTL::ZZ>("34465746750799872"), NTL::conv<NTL::ZZ>("6152208865296384"), NTL::conv<NTL::ZZ>("901713888280576"), NTL::conv<NTL::ZZ>("106811460943872"), NTL::conv<NTL::ZZ>("9993958981632"), NTL::conv<NTL::ZZ>("714316447744"), NTL::conv<NTL::ZZ>("37065719808"), NTL::conv<NTL::ZZ>("1285091328"), NTL::conv<NTL::ZZ>("25587712"), NTL::conv<NTL::ZZ>("213504"), NTL::conv<NTL::ZZ>("768"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(G0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return G0;
};

FpEX_elem _getH0(){
    FpEX_elem H0;
    
    std::vector<NTL::ZZ> cs = {NTL::conv<NTL::ZZ>("4096"), NTL::conv<NTL::ZZ>("8192"), NTL::conv<NTL::ZZ>("7168"), NTL::conv<NTL::ZZ>("3584"), NTL::conv<NTL::ZZ>("1104"), NTL::conv<NTL::ZZ>("208"), NTL::conv<NTL::ZZ>("22"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {16, 17, 18, 19, 20, 21, 22, 23};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(H0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return H0;
};

FpEX_elem _getI0(){
    FpEX_elem I0;

    std::vector<NTL::ZZ> cs = {NTL::conv<NTL::ZZ>("4096"), NTL::conv<NTL::ZZ>("393216"), NTL::conv<NTL::ZZ>("13664256"), NTL::conv<NTL::ZZ>("204701696"), NTL::conv<NTL::ZZ>("1285091328"), NTL::conv<NTL::ZZ>("4633214976"), NTL::conv<NTL::ZZ>("11161194496"), NTL::conv<NTL::ZZ>("19519451136"), NTL::conv<NTL::ZZ>("26077016832"), NTL::conv<NTL::ZZ>("27518124032"), NTL::conv<NTL::ZZ>("23468814336"), NTL::conv<NTL::ZZ>("16434548736"), NTL::conv<NTL::ZZ>("9552059904"), NTL::conv<NTL::ZZ>("4639475712"), NTL::conv<NTL::ZZ>("1889647104"), NTL::conv<NTL::ZZ>("645654016"), NTL::conv<NTL::ZZ>("184560432"), NTL::conv<NTL::ZZ>("43863552"), NTL::conv<NTL::ZZ>("8577664"), NTL::conv<NTL::ZZ>("1358208"), NTL::conv<NTL::ZZ>("169968"), NTL::conv<NTL::ZZ>("16192"), NTL::conv<NTL::ZZ>("1104"), NTL::conv<NTL::ZZ>("48"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(I0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return I0;
};

FpEX_elem _getJ0(){
    FpEX_elem J0;

    std::vector<NTL::ZZ> cs = {NTL::conv<NTL::ZZ>("512"), NTL::conv<NTL::ZZ>("1408"), NTL::conv<NTL::ZZ>("1664"), NTL::conv<NTL::ZZ>("1104"), NTL::conv<NTL::ZZ>("448"), NTL::conv<NTL::ZZ>("112"), NTL::conv<NTL::ZZ>("16"), NTL::conv<NTL::ZZ>("1")};
    std::vector<int> expon = {1, 2, 3, 4, 5, 6, 7, 8};

    assert(expon.size() == cs.size());

    for(size_t i = 0; i < cs.size(); i++){
        NTL::SetCoeff(J0, expon[i], NTL::conv<FpE_elem>(cs[i]));
    };

    return J0;
};


bool _getCoeffs(std::vector<FpE_elem> &output, std::vector<ecp> tors16, ec E){
    // tors16 contains two 16-torsion points on E where jE = j(E)

    ecp P1 = tors16[0];
    ecp P2 = tors16[1];
    assert((16*P1).is_identity());
    assert((16*P2).is_identity());
    assert(!((8*P1).is_identity()));
    assert(!((8*P2).is_identity()));
    assert(!((8*(P1-P2)).is_identity()));
    // ec E = phi.get_codomain();
    auto jE = E.j_invariant();

    std::vector<std::pair<ecp,std::pair<int, int>>> kerGens1, kerGens2;
    kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P1, std::pair<int, int>(2,4)));
    kerGens2.push_back(std::pair<ecp,std::pair<int, int>>(P2, std::pair<int, int>(2,4)));
    isog_chain phi_P1(kerGens1);
    isog_chain phi_P2(kerGens2);
    auto E1 = phi_P1.get_codomain();
    auto jE1 = E1.j_invariant();
    auto E2 = phi_P2.get_codomain();
    assert(E == phi_P1.get_domain());
    assert(E == phi_P2.get_domain());
    auto jE2 = E2.j_invariant();


    auto G0 = _getG0();
    auto H0 = _getH0();
    auto I0 = _getI0();
    auto J0 = _getJ0();

    FpEX_elem G, F1, F2;
    G = G0 - jE*H0;
    F1 = I0 - jE1*J0;
    F2 = I0 - jE2*J0;

    FpEX_elem f1 = GCD(G, F1);
    FpEX_elem f2 = GCD(G, F2);

    if (NTL::deg(f1) == 0 || NTL::deg(f2) == 0) { 
        return false;
    }

    
    std::vector<Fp2> r1 = {};
    std::vector<Fp2> r2 = {};
    if (NTL::deg(f1)==1) {
        r1.push_back(-f1[0]);
    }
    if (NTL::deg(f2)==1) {
        r2.push_back(-f2[0]);
    }


    if((NTL::deg(f1) > 1) ){

        // std::cout << "weird case 1 \n";

        auto v1 = FindRoots(f1);
        
        // there are several possibilies we need to find which one is the correct
        ec EE2 = phi_P1.get_middle_codomain(0);
        ec EE4 = phi_P1.get_middle_codomain(1);
        ec EE8 = phi_P1.get_middle_codomain(2);
        Fp2 t1 = EE2.disc()/E.disc();
        Fp2 t2 = EE4.disc()/EE2.disc();
        Fp2 t3 = EE8.disc()/EE4.disc();
        Fp2 t4 = E1.disc()/EE8.disc();

        bool found = false; 
        int i1 = 0;
        while (i1 < NTL::deg(f1)) {

            Fp2 w016 = v1[i1];
            Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
            Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
            Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
            Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
            Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
            Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
            found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 


            if (found) {r1.push_back(v1[i1]);}
            // if (!found){
                i1++;
            // }
            // 
        }
    }
    if (NTL::deg(f2)>1) 
    {

        auto v2 = FindRoots(f2);
        int i2 = 0;
        bool found = false;
        ec EE2 = phi_P2.get_middle_codomain(0);
        ec EE4 = phi_P2.get_middle_codomain(1);
        ec EE8 = phi_P2.get_middle_codomain(2);
        // std::cout << jE << " " << EE2.j_invariant() << " " << EE4.j_invariant() << " " << EE8.j_invariant() << " " << jE2 << "\n";
        Fp2 t1 = EE2.disc()/E.disc();
        Fp2 t2 = EE4.disc()/EE2.disc();
        Fp2 t3 = EE8.disc()/EE4.disc();
        Fp2 t4 = E2.disc()/EE8.disc();
        // std::cout << t1 << " " << t2 << " " << t3 << " " << t4 << "\n";
        while (i2 < NTL::deg(f2)) {
            Fp2 w016 = v2[i2];
            Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
            Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
            Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
            Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
            Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
            Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
            found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
            if (found) {r2.push_back(v2[i2]);}
            // if (!found){
                i2++;    
            // }
            // 
            
       }
        // assert(found);
        // r2 = v2[i2];    

    }
    if (r1.size() == 1 && r2.size() == 1) {
        output.push_back(r1[0]);
        output.push_back(r2[0]);
        return true;
    }    
        
    return false;
};


////////////////////////////////////////////////////////////////////////////
// Obtaining level structure on the domain using the known weber invariant 
////////////////////////////////////////////////////////////////////////////

std::vector<ecp> _getConsistentThreeTorsion(ec E, FpE_elem A, FpE_elem w, const Fp2k &Fext){
    
    auto tors3 = E.allTorsionPoints(Fext, 3, 1);
    assert((3*tors3[0]).is_identity());
    assert((3*tors3[1]).is_identity());
    assert((3*tors3[2]).is_identity());
    assert((3*tors3[3]).is_identity());

    // Eventually we want this to be x-only and so we give just the x-coords not the points in full
    
    FpE_elem x1 = tors3[0].aff_x();
    ecp P1 = tors3[0];
    int i = 1;
    FpE_elem x2 = tors3[1].aff_x();
    ecp P2 = tors3[1];
    while (x2 == x1) {
        i = i+1;
        x2 = tors3[i].aff_x();
        P2 = tors3[i];
    }
    i = i+1;
    FpE_elem x3 = tors3[i].aff_x();
    ecp P3 = tors3[i];
    while (x3 == x1 || x3==x2) {
        i = i+1;
        x3 = tors3[i].aff_x();
        P3 = tors3[i];
    }
    i = i+1;
    FpE_elem x4 = tors3[i].aff_x();
    ecp P4 = tors3[i];
    while (x4 == x1 || x4==x2 || x4 == x3 ) {
        i = i+1;
        x4 = tors3[i].aff_x();
        P4 = tors3[i];
    }
    std::vector<ecp> new_tors = {P1, P2, P3, P4};
    // std::vector<ecp> new_tors = {P1,P2,P3,P4};

    FpE_elem c1 = 2*A;
    FpE_elem c2 = -24*c1;
    FpE_elem gamma21 = c2/(c1 - 3*(x1*x2 + x3*x4));
    FpE_elem gamma22 = c2/(c1 - 3*(x1*x3 + x2*x4));
    // Can delete after we've fully debuged:
    FpE_elem gamma23 = c2/(c1 - 3*(x1*x4 + x3*x2));

    auto gamma2_up = power(w, 24) - 16;

    std::vector<int> ordering(4);


    if (gamma21*power(w, 8) == gamma2_up) {
        ordering = {0,1,2,3};
    } else if (gamma22*power(w, 8) == gamma2_up){
        ordering = {0,2,1,3};
    } else{
        assert(gamma23*power(w, 8) == gamma2_up);
        ordering = {0,3,2,1};
    };


    return {new_tors[ordering[0]], new_tors[ordering[1]], new_tors[ordering[2]], new_tors[ordering[3]]};
    
};


ecp _getConsistentTwoTorsion(ec E, FpE_elem w, const Fp2k &Fext){

    size_t i = 0;
    auto w24 = NTL::power(w,24);
    auto Edisc = 1/E.disc();
    auto tors2 = E.allTorsionPoints(Fext, 2, 1);

    while(true){
        isog phi2 = isogeny(tors2[i], 2);
        auto E2 = phi2.get_codomain();

        if(-E2.disc()*Edisc == w24){
            return tors2[i];
        }
        i += 1;
    }

};

std::vector<ecp> allCyclicGenerators(ec E, const Fp2k &Fext, int ell, int e) {

    std::vector<ecp> result = {};
    
    // generate the basis 
    auto Basis = E.torsionBasis(Fext, ell, e);
    assert((16*Basis.first).is_identity());
    assert(!(8*Basis.first).is_identity());
    assert((16*Basis.second).is_identity());
    assert(!(8*Basis.second).is_identity());
   

    // we generate the set of coeffs
    std::list<std::pair<Integer,Integer>> coeff_list = {};
    Integer ell_e = NTL::power_ZZ(ell,e);
    Integer iterate = Integer(0);
    while (iterate < ell_e) {
        coeff_list.push_back({Integer(1),iterate});
        if (NTL::GCD(iterate,Integer(ell))!=1) {
            coeff_list.push_back({iterate,Integer(1)});
        }
        iterate++;
    }
    for (auto coeff: coeff_list) {
        result.push_back(coeff.first* Basis.first + coeff.second * Basis.second);
        assert((16*(coeff.first* Basis.first + coeff.second * Basis.second)).is_identity());
        assert(!(8*(coeff.first* Basis.first + coeff.second * Basis.second)).is_identity());
    }

    return result;
}

std::vector<ecp> _getConsistentSixteenTorsion(ec E, FpE_elem w, const Fp2k &Fext){
    //std::cout << "in 16 thing" << std::endl;
    auto tors16 = allCyclicGenerators(E, Fext, 2, 4); 
    //std::cout << "allCyclicGens done" << std::endl;
    // E.allTorsionPoints(Fext, 2, 4);

    // for(size_t i = 0; i < tors16.size(); i++){
    //     std::cout << tors16[i] << std::endl;
    // }

    std::cout << "In 16 thing" << std::endl;
    auto w3 = NTL::power(w,3);
    auto jE = E.j_invariant();
    
    std::vector<std::pair<ecp,std::pair<int, int>>> kerGens1, kerGens2;
    FpEX_elem G;
    auto G0 = _getG0();
    auto H0 = _getH0();
    auto I0 = _getI0();
    auto J0 = _getJ0();
    G = G0 - jE*H0;

    std::cout << "Got the polys, entering loop..." << std::endl;

    for(size_t i = 1; i < tors16.size()-1; i++){
        // std::cout << "i: " << i << std::endl;
        ecp P1 = tors16[i];
        assert((16*P1).is_identity());
        assert(!(8*P1).is_identity());
        kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P1, std::pair<int, int>(2,4)));
        assert(kerGens1.size() == 1);
        isog_chain phi_P1(kerGens1);
        kerGens1.pop_back();
        auto E1 = phi_P1.get_codomain();
        auto jE1 = E1.j_invariant();
        FpEX_elem F1;
        F1 = I0 - jE1*J0;
        FpEX_elem f1 = GCD(G, F1);
        assert(NTL::deg(f1) == 1);
        
        for(size_t j = i+1; j < tors16.size(); j++){
            // std::cout << "j: " << j << std::endl;
            //std::cout << "looping i = " << i << ", j = " << j << std::endl;
            ecp P2 = tors16[j];
            assert( (16*P2).is_identity());
            assert( !(8*P2).is_identity() );
            //std::cout << P1.normalized() << std::endl;
            // std::cout << P2.normalized() << std::endl;
            if (!(8 * (P1-P2)).is_identity()) {
                assert( !(8 * (P1-P2)).is_identity());
                kerGens2.push_back(std::pair<ecp,std::pair<int, int>>(P2, std::pair<int, int>(2,4)));
                assert(kerGens2.size() == 1);
                isog_chain phi_P2(kerGens2);
                kerGens2.pop_back();
                auto E2 = phi_P2.get_codomain();
                auto jE2 = E2.j_invariant();
                FpEX_elem F2;
                F2 = I0 - jE2*J0;
                FpEX_elem f2 = GCD(G, F2);
                assert(NTL::deg(f2) == 1);
                std::vector<FpE_elem> cs = {FpE_elem(-f1[0]), FpE_elem(-f2[0])};
                bool thirdpow_check;
                auto w_check = _getWeberThirdPower(&thirdpow_check, cs);

                if(w_check == w3 && thirdpow_check){
                    return {tors16[i], tors16[j]};
                }
            }
        }
    }
    // Return two identity points to signal that it failed
    return {ecp(std::make_shared<const ec>(E), Fext), ecp(std::make_shared<const ec>(E), Fext)}; 
}


std::vector<std::vector<ecp>> GetLevelStructureFromWeber(ec E, FpE_elem w, const std::map<unsigned,Fp2k> &Fexts){

    // std::cout << "      Getting Torsion..." << std::endl;
    // auto tors = GetTorsionForLevelStructure(E, Fext);

    unsigned k2 = torsionToFieldDegree(Integer(2));
    auto jt2 = Fexts.find(k2);
    assert(jt2 != Fexts.end());  
    assert(jt2->first == jt2->second.k);
    unsigned k3 = torsionToFieldDegree(Integer(3));
    auto jt3 = Fexts.find(k3);
    assert(jt3 != Fexts.end());  
    assert(jt3->first == jt3->second.k);
    unsigned k16 = torsionToFieldDegree(Integer(16));
    auto jt16 = Fexts.find(k16);
    assert(jt16 != Fexts.end());  
    assert(jt16->first == jt16->second.k);

    std::cout << "      Getting 2-Torsion..." << std::endl;
    ecp P2 = _getConsistentTwoTorsion(E, w, jt2->second);

    std::cout << "      Getting 3-Torsion..." << std::endl;
    auto tors3 = _getConsistentThreeTorsion(E, E.a(), w, jt3->second);

    //  Get 16-torsion for w^3
    std::cout << "      Getting 16-Torsion..." << std::endl;
    auto tors16 = _getConsistentSixteenTorsion(E, w, jt16->second);

    return {{P2}, tors3, tors16};
};

FpE_elem GetWeberDomain(FpE_elem j){
    FpEX_elem W,f,g;
    NTL::SetCoeff(f,24,FpE_elem(1));
    NTL::SetCoeff(f,0,FpE_elem(-16));
    NTL::SetCoeff(g, 24, FpE_elem(-j));
    W = NTL::power(f,3) + g;
    return NTL::FindRoot(W);
}

FpE_elem GetWeberDomainFp(FpE_elem j){
    FpX_elem W,f,g;
    NTL::SetCoeff(f,24,Fp_elem(1));
    NTL::SetCoeff(f,0,Fp_elem(-16));
    NTL::SetCoeff(g, 24, Fp_elem(-NTL::ConstTerm(NTL::rep(j))));
    W = NTL::power(f,3) + g;
    Fp_elem w_Fp = NTL::FindRoot(W);
    FpE_elem w;
    NTL::conv(w, w_Fp);
    return w;
}

NTL::ZZ_pE GetWeberDomainFpBig(NTL::ZZ_pE j){
    NTL::ZZ_pX W,f,g;
    NTL::SetCoeff(f,24,NTL::ZZ_p(1));
    NTL::SetCoeff(f,0,NTL::ZZ_p(-16));
    NTL::SetCoeff(g, 24, NTL::ZZ_p(-NTL::ConstTerm(NTL::rep(j))));
    W = NTL::power(f,3) + g;
    //std::cout << W << std::endl;
    NTL::ZZ_p w_Fp = NTL::FindRoot(W);
    NTL::ZZ_pE w;
    NTL::conv(w, w_Fp);
    return w;
}

NTL::ZZ_pE GetWeberDomainBig(NTL::ZZ_pE j){
    NTL::ZZ_pEX W,f,g;
    NTL::SetCoeff(f,24,NTL::ZZ_pE(1));
    NTL::SetCoeff(f,0,NTL::ZZ_pE(-16));
    NTL::SetCoeff(g, 24, NTL::ZZ_pE(-j));
    W = NTL::power(f,3) + g;
    //std::cout << W << std::endl;
    NTL::FindRoots(W);
    return NTL::FindRoot(W);
}

std::vector<FpE_elem> GetWeberDomainAll(FpE_elem const &j){
    
    FpEX_elem W,f,g;
    NTL::SetCoeff(f,24,FpE_elem(1));
    NTL::SetCoeff(f,0,FpE_elem(-16));
    NTL::SetCoeff(g, 24, FpE_elem(-j));
    W = NTL::power(f,3) + g;
    std::vector<FpE_elem> out;

    for (auto const &r: NTL::FindRoots(W)) {
        out.push_back(r);
    }

    return out;
}


////////////////////////////////////////////////////
///// Getting weber invariant of image curve 
////////////////////////////////////////////////////

bool IsEqual(const std::pair<Integer,Integer> &pair1, const std::pair<Integer,Integer> &pair2) {
    return (pair1.first == pair2.first) && (pair1.second == pair2.second);
}

bool IsMinusMod3(const std::pair<Integer,Integer> &pair1, const std::pair<Integer,Integer> &pair2) {
    return ( pair1.first + pair2.first == 0 || pair1.first+pair2.first == 3) && ( pair1.second + pair2.second == 0 || pair1.second+pair2.second == 3);  
}

bool IsEqualOrMinusMod3(const std::pair<Integer,Integer> &pair1, const std::pair<Integer,Integer> &pair2) {
    return (IsEqual(pair1,pair2) || IsMinusMod3(pair1,pair2));
}

bool IsProjectiveEquivalent(const std::pair<Integer,Integer> &pair1, const std::pair<Integer,Integer> &pair2, const Integer &modulus) {
    return ((pair1.first * pair2.second - pair2.first * pair1.second)%modulus == 0);
}

bool IsVecProjectiveEquivalent(const std::vector<std::pair<Integer,Integer>> &vec1, const std::vector<std::pair<Integer,Integer>> &vec2, const Integer &modulus) {
    return ( (IsProjectiveEquivalent(vec1[0], vec2[0], modulus) && IsProjectiveEquivalent(vec1[1], vec2[1], modulus)) 
    ||  (IsProjectiveEquivalent(vec1[0], vec2[1], modulus) && IsProjectiveEquivalent(vec1[1], vec2[0], modulus)));
}

std::vector<std::pair<Integer,Integer>> PairApplication(const std::vector<std::pair<Integer,Integer>> &vec1, const std::vector<std::pair<Integer,Integer>> &vec2, const Integer &modulus) {
    return { {(vec1[0].first * vec2[0].first + vec1[1].first * vec2[0].second)%modulus,
              (vec1[0].second * vec2[0].first + vec1[1].second * vec2[0].second)%modulus}, 
             {(vec1[0].first * vec2[1].first + vec1[1].first * vec2[1].second)%modulus,
              (vec1[0].second * vec2[1].first + vec1[1].second * vec2[1].second)%modulus} };
} 

std::pair<Integer,Integer> Normalize(const std::pair<Integer,Integer> bas) {
    Fp_push push((Fp_integer(16)));
    if (bas.first%2 == 1) {
        return {Integer(1), NTL::conv<Integer>( NTL::conv<Fp>(bas.second)/NTL::conv<Fp>(bas.first) ) };
    }
    else {
        return {NTL::conv<Integer>( NTL::conv<Fp>(bas.first)/NTL::conv<Fp>(bas.second) ), Integer(1)  };
    }
}



bool IsWeberCoeffsEquivalent(const std::vector<std::vector<std::pair<Integer,Integer>>> coeffs1, const std::vector<std::vector<std::pair<Integer,Integer>>> coeffs2) {
    assert(IsEqual( coeffs1[0][0], coeffs2[0][0]));
    assert(
                    ( IsEqualOrMinusMod3(coeffs1[1][0],coeffs2[1][0]) && IsEqualOrMinusMod3(coeffs1[1][1],coeffs2[1][1]) )
                    || ( IsEqualOrMinusMod3(coeffs1[1][0],coeffs2[1][1]) && IsEqualOrMinusMod3(coeffs1[1][1],coeffs2[1][0]) )
                    || ( !IsEqualOrMinusMod3(coeffs1[1][0],coeffs2[1][0]) && !IsEqualOrMinusMod3(coeffs1[1][0],coeffs2[1][1]) && !IsEqualOrMinusMod3(coeffs1[1][1],coeffs2[1][1]) && !IsEqualOrMinusMod3(coeffs1[1][1],coeffs2[1][0]) )
                    );
    // assert(check3);


    std::pair<Integer,Integer> coeffs1161 = coeffs1[2][0];
    std::pair<Integer,Integer> coeffs1162 = coeffs1[2][1];
    std::pair<Integer,Integer> coeffs2161 = coeffs1[2][0];
    std::pair<Integer,Integer> coeffs2162 = coeffs1[2][1];
        
    bool check16 = false;    
    std::vector<std::vector<std::pair<Integer,Integer>>> alt_coeff = {};
    alt_coeff.push_back( { {Integer(1), Integer(0)}, {Integer(0), Integer(1)} } );
    alt_coeff.push_back( { {Integer(2), Integer(1)}, {Integer(1), Integer(10)} } );
    alt_coeff.push_back( { {Integer(1), Integer(2)}, {Integer(10), Integer(1)} } );
    alt_coeff.push_back( { {Integer(1), Integer(4)}, {Integer(4), Integer(1)} } );
    alt_coeff.push_back( { {Integer(1), Integer(8)}, {Integer(8), Integer(1)} } );
    alt_coeff.push_back( { {Integer(1), Integer(12)}, {Integer(12), Integer(1)} } );
    alt_coeff.push_back( { {Integer(1), Integer(6)}, {Integer(14), Integer(1)} } );
    alt_coeff.push_back( { {Integer(6), Integer(1)}, {Integer(1), Integer(14)} } );
    for (auto aa : alt_coeff) {
        auto bb = PairApplication(coeffs2[2], aa, Integer(16));
        // std::cout << bb[0].first << " " << bb[0].second << " " << bb[1].first << " " << bb[1].second;    
        check16 = check16 || IsVecProjectiveEquivalent(coeffs1[2],bb , Integer(16)); 
        // bb = {Normalize(bb[0]), Normalize(bb[1]) };
        // std::cout << "      " << bb[0].first << " " << bb[0].second << " " << bb[1].first << " " << bb[1].second << " " << check16 << " \n";
    } 


    return check16;
}


Fp2 WeberGetFromEnum(const std::vector<std::vector<std::pair<Integer,Integer>>> &coeffs, const weber_enum &webdat, const std::pair<mat_Fp, mat_Fp> change_mats) {

    std::vector<std::pair<Integer,Integer>> order_2 = {};
    std::vector<std::pair<Integer,Integer>> order_3 = {};
    std::vector<std::pair<Integer,Integer>> order_16 = {};
    
    assert(coeffs.size() == 3);

    // order 2
    {
        Fp_push push((Fp_integer(2)));
        order_2.push_back( Zp_mat_application(change_mats.first, coeffs[0][0]) );
    }
    // order 3
    {
        Fp_push push((Fp_integer(3)));
        order_3.push_back( Zp_mat_application(change_mats.second, coeffs[1][0])); 
        order_3.push_back( Zp_mat_application(change_mats.second, coeffs[1][1]));
    }
    // order 16
    {
        Fp_push push((Fp_integer(16)));
        order_16.push_back( Zp_mat_application(change_mats.first, coeffs[2][0]));
        order_16.push_back( Zp_mat_application(change_mats.first, coeffs[2][1]));
    }

    int j2;
    if (order_2[0].first == 0) {j2 = 0;}
    else if (order_2[0].second == 0) {j2 = 1;}
    else {j2 = 2;} 

    int j3;
    if (order_3[0].first == 0) {
        assert(order_3[1].first != 0);
        if (order_3[1].second == 0) {j3 = 0;}
        else if (order_3[1].first == order_3[1].second) {j3 = 2;}
        else {j3 = 1;}
    }
    else if (order_3[1].first == 0) {
        if (order_3[0].second == 0) {j3 = 0;}
        else if (order_3[0].first == order_3[0].second) {j3 = 2;}
        else {j3 = 1;}
    }
    else {
        assert(order_3[0].first !=0 && order_3[1].first != 0);
        if (order_3[0].second == 0) {
            if (order_3[1].first == order_3[1].second) { j3 = 1;}
            else {j3 = 2;}
        }
        else if (order_3[1].second == 0) {
            if (order_3[0].first == order_3[0].second) { j3 = 1;}
            else {j3 = 2;}
        }
        else {j3 = 0;}
    }

    std::vector<std::vector<std::pair<Integer,Integer>>> new_coeffs = {order_2, order_3, order_16};

    Fp_push push((Fp_integer(16)));
    long web_index;
    // alternate way to find the correct valus 
    // the treament first depends on the values mod 2
    if (order_2[0].first == 0 && order_2[0].second  == 1) {
        
        std::vector<long> list_of_coeffs = { 1, 7, 13, 3, 9, 15, 5, 11 };
        std::vector<long> list_of_indices = {0, 3, 6, 1, 4, 7, 2, 5};

        // checking if we need to swap 
        if (order_16[0].second%2 == 1) {
            assert(order_16[1].second%2 == 0);
            Integer tmp = order_16[0].first;
            order_16[0].first = order_16[1].first;
            order_16[1].first = tmp;
            tmp = order_16[0].second;
            order_16[0].second = order_16[1].second;
            order_16[1].second = tmp;

        }
        assert(order_16[0].second%2 == 0);
        assert(order_16[1].second%2 == 1);
        order_16[0].second = NTL::conv<Integer>(NTL::conv<Fp>(order_16[0].second)/NTL::conv<Fp>(order_16[0].first));
        order_16[0].first = 1;
        order_16[1].second = NTL::conv<Integer>(NTL::conv<Fp>(order_16[1].second)/NTL::conv<Fp>(order_16[1].first));
        order_16[1].first = 1;

        long index1 = NTL::conv<long>(order_16[0].second) / 2;
        long index2 = list_of_indices[NTL::conv<long>(order_16[1].second - 1) / 2];
        web_index = (list_of_coeffs[(8 + index2 - index1)%8] - 1)/2;
        // std::cout << "alt index = " << final_index << "\n";
    }
    else if (order_2[0].first == 1 && order_2[0].second  == 0) {
        
        std::vector<long> list_of_coeffs = { 1, 7, 5, 11, 9, 15, 13, 3 };
        std::vector<long> list_of_indices = {0, 7, 2, 1, 4, 3, 6, 5};

        // checking if we need to swap 
        if (order_16[0].first%2 == 1) {
            assert(order_16[1].first%2 == 0);
            Integer tmp = order_16[0].first;
            order_16[0].first = order_16[1].first;
            order_16[1].first = tmp;
            tmp = order_16[0].second;
            order_16[0].second = order_16[1].second;
            order_16[1].second = tmp;

        }
        assert(order_16[0].first%2 == 0);
        assert(order_16[1].first%2 == 1);
        order_16[0].first = NTL::conv<Integer>(NTL::conv<Fp>(order_16[0].first)/NTL::conv<Fp>(order_16[0].second));
        order_16[0].second = 1;
        order_16[1].second = NTL::conv<Integer>(NTL::conv<Fp>(order_16[1].second)/NTL::conv<Fp>(order_16[1].first));
        order_16[1].first = 1;

        long index1 = NTL::conv<long>(order_16[0].first) / 2;
        long index2 = list_of_indices[NTL::conv<long>(order_16[1].second - 1) / 2];
        web_index = (list_of_coeffs[(8 + index2 - index1)%8] - 1)/2;
        // std::cout << "alt index = " << final_index << "\n";
    }
    else {
        std::vector<long> list_of_coeffs;
        std::vector<long> list_of_indices; 
        // checking if we need to swap 
        if (order_16[1].first%2 == 1) {
            assert(order_16[1].first%2 == 1);
            Integer tmp = order_16[0].first;
            order_16[0].first = order_16[1].first;
            order_16[1].first = tmp;
            tmp = order_16[0].second;
            order_16[0].second = order_16[1].second;
            order_16[1].second = tmp;

        }
        assert(order_16[0].first%2 == 1);
        assert(order_16[0].second%2 == 0);
        assert(order_16[1].first%2 == 0);
        assert(order_16[1].second%2 == 1);
        order_16[0].second = NTL::conv<Integer>(NTL::conv<Fp>(order_16[0].second)/NTL::conv<Fp>(order_16[0].first));
        order_16[0].first = 1;
        order_16[1].first = NTL::conv<Integer>(NTL::conv<Fp>(order_16[1].first)/NTL::conv<Fp>(order_16[1].second));
        order_16[1].second = 1;
        // the list depends on the value mod 4
        if (( 16 +  order_16[0].second - order_16[1].first) %4 == 2) {
            list_of_coeffs = {0, 2, 4, 6, 8, 10, 12, 14 };
            list_of_indices = {0, 1, 2, 3, 4, 5, 6, 7};
        }
        else {
            list_of_coeffs = {0, 10, 4, 14, 8, 2, 12, 6};
            list_of_indices = {0, 5, 2, 7, 4, 1, 6, 3};
        }
        long index1 = NTL::conv<long>(order_16[0].second) / 2;
        long index2 = list_of_indices[NTL::conv<long>(order_16[1].first) / 2];
        web_index = (list_of_coeffs[(8 + index2 - index1)%8])/2;
    }

    // for (long i = 0; i < 8 ; i++) {
    //     if (IsWeberCoeffsEquivalent(new_coeffs, webdat[i + 8 * j3 + 24 * j2].second)) {
    //         assert(i == web_index);
    //         return webdat[i + 8 * j3 +  24 * j2].first;
    //     } 
    // }
    // assert(0);
    assert(IsWeberCoeffsEquivalent(new_coeffs, webdat[web_index + 8 * j3 + 24 * j2].second));
    assert( web_index + 8* j3 + 24 * j2 < 72 );
    return webdat[ web_index + 8* j3 + 24 * j2 ].first;
}

bool GetWeberOfLevelStruct_j0(Fp2 *w, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts) {
    weber_bas w0 = { levelstruc[2][0], levelstruc[2][1], levelstruc[1][0], levelstruc[1][1]};
    weber_enum wb0 = EnumerateAllWeberFast(w0, Fexts);
    *w = wb0[48].first;
    return true;
}


bool GetWeberOfLevelStruct(Fp2 *w, ec E, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts) {


    if (E.j_invariant() == 0 || E.j_invariant() == FpE_elem(1728)){
        //std::cout << "Hopefully this is called?" << std::endl;
        return GetWeberOfLevelStruct_j0(w, levelstruc, Fexts);
    }

    ecp P2 = levelstruc[0][0];
    std::vector<ecp> tors3 = levelstruc[1];
    std::vector<ecp> tors16 = levelstruc[2];

    // std::cout << "  Getting Weber^8..." << std::endl;
    // std::cout << "      Computing the t-invariant" << std::endl;
    auto t_invariant = _getTInvariant(P2, E);
    // std::cout << "      done: " << t_invariant << std::endl;
    // std::cout << "      Computing Gamma_2" << std::endl;
    auto gamma2 = _getGammaTwoInvariant(tors3, E); 
    // std::cout << "      done: " << gamma2 << std::endl;
    auto w8 = _getWeberEighthPower(gamma2, t_invariant);
    // std::cout << "Got Weber^8: " << w8 << std::endl;

    // std::cout << "Getting Weber^3..." << std::endl;
    // std::cout << "      Getting necessary constants" << std::endl;
    std::vector<Fp2> cs = {};
    auto cs_check = _getCoeffs(cs, tors16, E);
    if (!cs_check) {
        return false;
    }
    // std::cout << "      done: " << cs[0] << ", " << cs[1] << std::endl;
    bool thirdpow_check;
    auto w3 = _getWeberThirdPower(&thirdpow_check, cs);
    //call the special one
    if (!(thirdpow_check)) {
        return GetWeberOfLevelStruct_j0(w, levelstruc, Fexts);
    }
    // std::cout << "Got Weber^3: " << w3 << std::endl;

    // std::cout << "Computing Weber invariant..." << std::endl;

    if (w3 == NTL::conv<Fp2>(0)) {
        return false;
    }

    FpEX_elem W1, W2;
    NTL::SetCoeff(W1, 3, FpE_elem(1));
    NTL::SetCoeff(W1, 0, FpE_elem(-w3));
    NTL::SetCoeff(W2, 8, FpE_elem(1));
    NTL::SetCoeff(W2, 0, FpE_elem(-w8));

    auto g = NTL::GCD(W1,W2);
    assert(NTL::deg(g) <= 1);
    if (NTL::deg(g)==0) {
        // std::cout << "problem in weber for curve " << E.j_invariant()    << " with w3 = " << w3 << " w8 = " << w8 << "\n"; 
        return false;
    }
    else {
        assert(NTL::deg(g) == 1);
        *w = Fp2(-g[0]);
        return true;
    }

    // return FpE_elem(-g[0]); 
}

bool GetWeberOfImage_chain(Fp2 *w, isog_chain phi, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts, bool retry){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gets Weber invariant of E' where phi: E --> E' with E: y^2 = x^3 + Ax + B, and the weber
    //      invariant of E is w_domain
    // We need a two torsion point P2, two 16-torsion points {P16, Q16} and three torsion points 
    // tors3 on E
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    ecp P2 = levelstruc[0][0];
    std::vector<ecp> tors3 = levelstruc[1];
    std::vector<ecp> tors16 = levelstruc[2];

    //std::cout << "Pushing points" << std::endl;
    std::vector<std::vector<ecp>> new_levelstruct = { {phi(P2)}, { phi(tors3[0]), phi(tors3[1]), phi(tors3[2]), phi(tors3[3]) }, { phi(tors16[0]), phi(tors16[1]) } };
    //std::cout << "Finished pushing points, now time for weber stuff" << std::endl;
    bool check = GetWeberOfLevelStruct(w, phi.get_codomain(), new_levelstruct, Fexts); 
    if (!(check) && retry) {
        check = GetWeberOfLevelStruct_j0(w, new_levelstruct, Fexts); 
    }
    return check;
};

bool GetWeberOfImage(Fp2 *w, ec E, isog phi, std::vector<std::vector<ecp>> levelstruc, const std::map<unsigned,Fp2k> &Fexts){
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Gets Weber invariant of E' where phi: E --> E' with E: y^2 = x^3 + Ax + B, and the weber
    //      invariant of E is w_domain
    // We need a two torsion point P2, two 16-torsion points {P16, Q16} and three torsion points 
    // tors3 on E
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    FpE_elem jE = E.j_invariant();

    ecp P2 = levelstruc[0][0];
    std::vector<ecp> tors3 = levelstruc[1];
    std::vector<ecp> tors16 = levelstruc[2];

    std::vector<std::vector<ecp>> new_levelstruct = { {phi(P2)}, { phi(tors3[0]), phi(tors3[1]), phi(tors3[2]), phi(tors3[3]) }, { phi(tors16[0]), phi(tors16[1]) } };

    bool check = GetWeberOfLevelStruct( w ,phi.get_codomain(), new_levelstruct, Fexts); 
    assert(check);
    return check;
    
};

////////////////////////////////////////////////////
///// For enumerating all Weber invariants
////////////////////////////////////////////////////


Fp2X eval_phi5_weber( Fp2 w ) {

    Fp2X res;

    NTL::SetCoeff(res, 6, NTL::conv<Fp2>(1));
    NTL::SetCoeff(res, 5, -NTL::power(w,5));
    NTL::SetCoeff(res, 4, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 3, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 2, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 1, 4 * w);
    NTL::SetCoeff(res, 0, NTL::power(w,6));
    return res;

}

Fp2X eval_phi7_weber( Fp2 w ) {

    Fp2X res;

    NTL::SetCoeff(res, 8, NTL::conv<Fp2>(1));
    NTL::SetCoeff(res, 7, -NTL::power(w,7));
    NTL::SetCoeff(res, 4, 7*NTL::power(w,4));
    NTL::SetCoeff(res, 5, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 6, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 3, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 2, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 1, -8 * w);
    NTL::SetCoeff(res, 0, NTL::power(w,8));
    return res;

}

Fp2X eval_phi11_weber( Fp2 w ) {
    Fp2X res;
    NTL::SetCoeff(res, 12, NTL::conv<Fp2>(1));
    NTL::SetCoeff(res, 11, -NTL::power(w,11));
    NTL::SetCoeff(res, 9, 11*NTL::power(w,9));
    NTL::SetCoeff(res, 7, -44*NTL::power(w,7));
    NTL::SetCoeff(res, 5, 88*NTL::power(w,5));
    NTL::SetCoeff(res, 3, -88*NTL::power(w,3));
    NTL::SetCoeff(res, 10, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 8, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 6, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 4, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 2, NTL::conv<Fp2>(0));
    NTL::SetCoeff(res, 1, 32 * w);
    NTL::SetCoeff(res, 0, NTL::power(w,12));
    return res;
}

std::vector<std::vector<ecp>> BasCoeffToLevelStructure(const weber_bas &web, const std::vector<std::vector<std::pair<Integer,Integer>>> coeff) {
    auto c2 = coeff[0][0];
    auto c31 = coeff[1][0];
    auto c32 = coeff[1][1]; 
    ecp P3 = c31.first * web.P3 + c31.second * web.Q3;
    ecp Q3 = c32.first * web.P3 + c32.second * web.Q3;
    auto c161 = coeff[2][0];
    auto c162 = coeff[2][1];
    ecp P2 = (c2.first*8)*web.P16 + (c2.second*8)*web.Q16;
    assert((2*P2).is_identity());
    assert(!(P2).is_identity());
    assert((3*P3).is_identity());
    assert(!(P3).is_identity());
    std::vector<ecp> tors3 = {P3, Q3, P3+Q3, P3+2*Q3};
    std::vector<ecp> tors16 = {c161.first*web.P16+c161.second*web.Q16, c162.first*web.P16+c162.second*web.Q16};
    std::vector<std::vector<ecp>> levelstruct = { {P2}, tors3, tors16};
    return levelstruct;
}

bool BasCoeffToWeber(Fp2 *w, const weber_bas &web, const std::vector<std::vector<std::pair<Integer,Integer>>> coeff, const std::map<unsigned,Fp2k> &Fexts) {    
    auto c2 = coeff[0][0];
    auto c31 = coeff[1][0];
    auto c32 = coeff[1][1]; 
    ecp P3 = c31.first * web.P3 + c31.second * web.Q3;
    ecp Q3 = c32.first * web.P3 + c32.second * web.Q3;
    auto c161 = coeff[2][0];
    auto c162 = coeff[2][1];
    ecp P2 = (c2.first*8)*web.P16 + (c2.second*8)*web.Q16;
    assert((2*P2).is_identity());
    assert(!(P2).is_identity());
    assert((3*P3).is_identity());
    assert(!(P3).is_identity());
    std::vector<ecp> tors3 = {P3, Q3, P3+Q3, P3+2*Q3};
    std::vector<ecp> tors16 = {c161.first*web.P16+c161.second*web.Q16, c162.first*web.P16+c162.second*web.Q16};
    std::vector<std::vector<ecp>> levelstruct = { {P2}, tors3, tors16};
    // clock_t t = clock();
    bool check = GetWeberOfLevelStruct(w, P2.curve(), levelstruct, Fexts);
    // std::cout << "weber time = " << (double) (clock() - t)/CLOCKS_PER_SEC << "\n";

    // If the first computation has failed, possibly due to some weird special case with some power of 2 small endomorphisms
    // We try alternate 16 torsion points that should give the same Weber invariant but that hopefully are not in the same bad special case
    if (!check) {
        std::cout << "need alt points \n";
        std::vector<std::pair<std::pair<Integer,Integer>,std::pair<Integer,Integer>>> alt_coeff = {};
        alt_coeff.push_back( { {Integer(2), Integer(1)}, {Integer(1), Integer(10)} } );
        alt_coeff.push_back( { {Integer(1), Integer(4)}, {Integer(4), Integer(1)} } );
        alt_coeff.push_back( { {Integer(1), Integer(6)}, {Integer(14), Integer(1)} } );
        alt_coeff.push_back( { {Integer(6), Integer(1)}, {Integer(1), Integer(14)} } );
        alt_coeff.push_back( { {Integer(1), Integer(2)}, {Integer(10), Integer(1)} } );
        for (auto cc : alt_coeff) {
                std::vector<std::vector<ecp>> alt_levelstruct = { {P2}, tors3, {cc.first.first*tors16[0] + cc.first.second * tors16[1], cc.second.first*tors16[0] + cc.second.second*tors16[1]}};
                check = GetWeberOfLevelStruct(w, P2.curve(), alt_levelstruct, Fexts);
                if (check) {
                    return check;
                }

        }
    }

    // Stil doesn't work
    // We can compute 5 and 7 isogenous coefficients and do a GCD
    if (!check) {
        std::cout << "last resort !\n";
        Fp2 w5,w7;
    
        ec E = P2.curve();

        // first the 5 torsion 
        unsigned k = torsionToFieldDegree(Integer(5));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());  
        assert(jt->first == jt->second.k);
        auto bas5 = E.torsionBasis(jt->second, int(5), 1);
        k = torsionToFieldDegree(Integer(7));
        jt = Fexts.find(k);
        assert(jt != Fexts.end());  
        assert(jt->first == jt->second.k);
        auto bas7 = E.torsionBasis(jt->second, int(7), 1);

        std::vector<std::pair<ecp,std::pair<int, int>>> ker5, ker7;

        std::vector<std::pair<Integer,Integer>> c5 = {};
        c5.push_back({Integer(0), Integer(1)});
        c5.push_back({Integer(1), Integer(0)});
        c5.push_back({Integer(1), Integer(2)});
        c5.push_back({Integer(1), Integer(3)});
        c5.push_back({Integer(1), Integer(4)});

        bool check5 = false;
        for (auto cc5 : c5) {
            ker5.push_back(std::pair<ecp,std::pair<int, int>>(cc5.first*bas5.first + cc5.second * bas5.second, std::pair<int, int>(5,1)));
            isog_chain phi5(ker5);
            ker5.pop_back();
            check5 = GetWeberOfImage_chain(&w5, phi5, levelstruct, Fexts);
            if (check5) {
                break;
            }
        }
        assert(check5);

        std::vector<std::pair<Integer,Integer>> c7 = {};  
        c7.push_back( { Integer(0), Integer(1)});
        c7.push_back( { Integer(1), Integer(0)});
        c7.push_back( { Integer(1), Integer(2)});
        c7.push_back( { Integer(1), Integer(3)});
        c7.push_back( { Integer(1), Integer(4)});
        c7.push_back( { Integer(1), Integer(1)});
        c7.push_back( { Integer(1), Integer(5)});
        c7.push_back( { Integer(1), Integer(6)});

        bool check7 = false;
        for (auto cc7: c7) {
            ker7.push_back(std::pair<ecp,std::pair<int, int>>(cc7.first*bas7.first + cc7.second * bas7.second, std::pair<int, int>(7,1)));
            isog_chain phi7(ker7);
            ker7.pop_back();
            check7 = GetWeberOfImage_chain(&w7, phi7, levelstruct, Fexts);
            if (check7) {
                break;
            }
        }
        assert(check7);     
        
         
        if (check5 && check7) {
            Fp2X p5 = eval_phi5_weber(w5);
            Fp2X p7 = eval_phi7_weber(w7);
            Fp2X g = NTL::GCD(p5, p7);
            check = (NTL::deg(g)==1);
            if (check) {
                *w = -g[0];
            } 
        }
        
    }

    return check;
    

}

std::vector<std::vector<std::vector<std::pair<Integer,Integer>>>> EnumerateAllWeberCoeff() {

    std::vector<std::vector<std::vector<std::pair<Integer,Integer>>>> result = {};
    // order 2 
    std::vector<std::pair<Integer,Integer>> order_2 = { {Integer(0),Integer(1)}, {Integer(1),Integer(0)}, {Integer(1), Integer(1)} };
    // order 3
    std::vector<std::vector<std::pair<Integer,Integer>>> order_3 = { { {Integer(1), Integer(0)}, {Integer(0), Integer(1)} }, { {Integer(1),Integer(0)}, {Integer(1), Integer(1)} }, { {Integer(1), Integer(0)}, {Integer(1), Integer(2)} } };

    // // pairs of relevant order_16 generators
    std::vector<std::vector<std::pair<Integer,Integer>>> order_16_pairs = {};
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(3)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(5)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(7)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(9)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(11)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(13)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(15)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(3)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(5)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(7)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(9)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(11)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(13)} } );
    order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(15)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(0),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(2),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(4),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(6),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(8),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(10),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(12),Integer(1)} } );
    order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(14),Integer(1)} } );

    // now we combine everything togeter 
    for (size_t i = 0; i < order_2.size(); i++) {
        for (size_t j = 0; j < order_3.size(); j++) {
            for (size_t k =0; k < order_16_pairs.size(); k++) {
                if (!(((order_2[i].first == (order_16_pairs[k][0].first % 2)) && (order_2[i].second == (order_16_pairs[k][0].second % 2))) || ((order_2[i].first == (order_16_pairs[k][1].first % 2)) && (order_2[i].second == (order_16_pairs[k][1].second % 2))))) {
                    result.push_back( { {order_2[i]}, order_3[j], order_16_pairs[k] } );
                }
                
            }
        }
    }
    return result;
}



std::vector<std::pair< Fp2, std::vector<std::vector<std::pair<Integer,Integer>>>>> EnumerateAllWeber(const weber_bas &web, const std::vector<std::vector<std::vector<std::pair<Integer,Integer>>>> &coeff_list, const std::map<unsigned,Fp2k> &Fexts) {

    std::cout << "curve = " << web.P3.curve().j_invariant() << "\n";
    assert(web.P3.curve() == web.P16.curve());
    std::vector<std::pair<Fp2, std::vector<std::vector<std::pair<Integer,Integer>>>>> result = {};

    std::cout << "number of levelstruct = " << coeff_list.size() << "\n";

    clock_t t = clock();
    for (size_t i = 0; i< coeff_list.size(); i++) {
        if (result.size() == 72) {
            goto endloop;
        }

        Fp2 w;
        bool check = BasCoeffToWeber(&w, web, coeff_list[i], Fexts);
        // assert(check);


        if (check) {
            result.push_back( {w, coeff_list[i]} );
        }

        

    }
    endloop: 

    std::cout << "number of weber invariants = " << result.size() << "\n";
    std::cout << "time to compute them :" << (double) (clock() - t)/CLOCKS_PER_SEC << "\n";
    return result;
}

NTL::Vec<Fp2X> add_bivariate(const NTL::Vec<Fp2X> &v1, const NTL::Vec<Fp2X> &v2) {
    NTL::Vec<Fp2X> result;
    int max = std::max(v1.length(), v2.length());
    int min = std::min(v1.length(), v2.length());
    result.SetLength(max);
    for (int i = 0; i < max ; i++) {
        if (i < min) {
            result[i] = v1[i] + v2[i];
        }
        else {
            if (min == v1.length()) {
                result[i] = v2[i];
            }
            else {
                result[i] = v1[i];
            }
        }
    }
    return result;
}


weber_enum SpecialEnumFastWeber(const weber_bas web, const std::map<unsigned,Fp2k> &Fexts) {
        // std::cout << "special enum ! \n";

        weber_enum result = {};
        ec E = web.P3.curve();

        // std::cout << "special case for j=" << E.j_invariant() << "\n";

        // first the 5 torsion 
        unsigned k = torsionToFieldDegree(Integer(5));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());  
        assert(jt->first == jt->second.k);
        auto bas5 = E.torsionBasis(jt->second, int(5), 1);
        k = torsionToFieldDegree(Integer(7));
        jt = Fexts.find(k);
        assert(jt != Fexts.end());  
        assert(jt->first == jt->second.k);
        auto bas7 = E.torsionBasis(jt->second, int(7), 1);
        k = torsionToFieldDegree(Integer(11));
        jt = Fexts.find(k);
        assert(jt != Fexts.end());  
        assert(jt->first == jt->second.k);
        auto bas11 = E.torsionBasis(jt->second, int(11), 1);

        std::vector<std::pair<ecp,std::pair<int, int>>> ker5, ker7, ker11;

        std::vector<std::pair<Integer,Integer>> c5 = {};
        weber_enum w5 = {};
        weber_enum w7 = {};
        weber_enum w11 = {};
        c5.push_back({Integer(0), Integer(1)});
        c5.push_back({Integer(1), Integer(0)});
        c5.push_back({Integer(1), Integer(2)});
        c5.push_back({Integer(1), Integer(3)});
        c5.push_back({Integer(1), Integer(4)});

        std::vector<std::pair<Integer,Integer>> c7 = {};  
        c7.push_back( { Integer(0), Integer(1)});
        c7.push_back( { Integer(1), Integer(0)});
        c7.push_back( { Integer(1), Integer(2)});
        c7.push_back( { Integer(1), Integer(3)});
        c7.push_back( { Integer(1), Integer(4)});
        c7.push_back( { Integer(1), Integer(1)});
        c7.push_back( { Integer(1), Integer(5)});
        c7.push_back( { Integer(1), Integer(6)});

        std::vector<std::pair<Integer,Integer>> c11 = {};  
        c11.push_back( { Integer(0), Integer(1)});
        c11.push_back( { Integer(1), Integer(0)});
        c11.push_back( { Integer(1), Integer(2)});
        c11.push_back( { Integer(1), Integer(3)});
        c11.push_back( { Integer(1), Integer(4)});
        c11.push_back( { Integer(1), Integer(1)});
        c11.push_back( { Integer(1), Integer(5)});
        c11.push_back( { Integer(1), Integer(6)});
        c11.push_back( { Integer(1), Integer(7)});
        c11.push_back( { Integer(1), Integer(8)});
        c11.push_back( { Integer(1), Integer(9)});
        c11.push_back( { Integer(1), Integer(10)});

        bool check5 = false;
        bool check = false;
        for (size_t i5 = 0; i5 < c5.size(); i5++) {
            auto cc5 = c5[i5];

            ker5.push_back(std::pair<ecp,std::pair<int, int>>(cc5.first*bas5.first + cc5.second * bas5.second, std::pair<int, int>(5,1)));
            isog_chain phi5(ker5);
            check5 = (phi5.get_codomain().j_invariant() != Fp2(1728) && phi5.get_codomain().j_invariant() != Fp2(287496) && phi5.get_codomain().j_invariant() != Fp2(0));
            ker5.pop_back();
            if (check5) {
                w5 = EnumerateAllWeberFast( {phi5(web.P16), phi5(web.Q16), phi5(web.P3), phi5(web.Q3)}, Fexts);
                check5 = (w5.size() == 72);
            }
            if (!check5) {
                continue;
            }
            bool check7 = false;
            for (size_t i7 = 0; i7 < c7.size(); i7++) {
                auto cc7 = c7[i7];
                ker7.push_back(std::pair<ecp,std::pair<int, int>>(cc7.first*bas7.first + cc7.second * bas7.second, std::pair<int, int>(7,1)));
                isog_chain phi7(ker7);
                check7 = (phi7.get_codomain().j_invariant() != Fp2(1728) && phi7.get_codomain().j_invariant() != Fp2(287496) && phi7.get_codomain().j_invariant() != Fp2(0));

                ker7.pop_back();
                if (check7) {
                    w7 = EnumerateAllWeberFast( {phi7(web.P16), phi7(web.Q16), phi7(web.P3), phi7(web.Q3)}, Fexts);
                    check7 = (w7.size() == 72);
                }
                if (!check7)
                    continue;

                bool check11 = false;
                for (size_t i11 = 0; i11 < c11.size(); i11++) {
                    auto cc11 = c11[i11];

                    //std::cout << i5 << " " << i7 <<  " " << i11 << "\n";
                    ker11.push_back(std::pair<ecp,std::pair<int, int>>(cc11.first*bas11.first + cc11.second * bas11.second, std::pair<int, int>(11,1)));
                    isog_chain phi11(ker11);

                    check11 = (phi11.get_codomain().j_invariant() != Fp2(1728) && phi11.get_codomain().j_invariant() != Fp2(287496) && phi11.get_codomain().j_invariant() != Fp2(0));
                    ker11.pop_back();
                    if (check11) {
                        w11 = EnumerateAllWeberFast( {phi11(web.P16), phi11(web.Q16), phi11(web.P3), phi11(web.Q3)}, Fexts);
                        check11 = (w11.size() == 72);
                    }
                    if (!check11)
                        continue;

                    if (check5 && check7 && check11) {
                        // std::cout << "tried a full round \n"; 
                        assert(w5.size() == w7.size());
                        assert(w11.size() == w7.size());
                        check = true;
                        for (unsigned i = 0; i<w5.size(); i++) 
                        {
                            Fp2X p5 = eval_phi5_weber(w5[i].first);
                            Fp2X p7 = eval_phi7_weber(w7[i].first);
                            Fp2X p11 = eval_phi11_weber(w11[i].first);
                            Fp2X g = NTL::GCD(p5, p7);
                            g = NTL::GCD(g, p11);
                            check = check && (NTL::deg(g)==1);
                            if (check) {
                                // *w = -g[0];
                                result.push_back( {-g[0], w5[i].second} );

                            } 
                            else {
                                result = {};
                            }
                        }
                    }  
                    if (result.size() == 72) {
                        goto endloop;
                    }
                }
            }
            
        }
    
        assert(result.size()==72);

        endloop:
        assert(result.size()==72);    
        return result;
        
}


weber_enum OldEnumerateAllWeberFast(const weber_bas &web, const std::map<unsigned,Fp2k> &Fexts) {


    clock_t tot = clock(); (void) tot;
    clock_t t = clock();
    ec E = web.P3.curve();
    auto jE = E.j_invariant();
    

    // first check if we are not in some special case
    if (jE ==Fp2(1728) || jE == Fp2(287496)) {
        return SpecialEnumFastWeber(web, Fexts);
    }

    // first we start by enumerating the X016 points that because this where things can fail most dramatically 
    std::vector<std::pair<Fp2, std::vector<std::vector<std::pair<Integer,Integer>>>>> result = {};
    // std::unordered_set<Jinv, JinvHash, JinvEqual> weber_list = {}; 

    // std::cout << "curve = " << web.P3.curve().j_invariant() << "\n";
    assert(web.P3.curve() == web.P16.curve());

    
    // Initalizing all needed polynomials
    auto R1 = GetBivariateResultant();
    auto R2 = GetTrivariateResultant();
    auto G0 = _getG0();
    auto H0 = _getH0();
    auto I0 = _getI0();
    auto J0 = _getJ0();
    Fp2X G, F1, F2;
    
    G = G0 - jE*H0;
    FpEX_elem f1,f2,f3,f4;
    NTL::Vec<FpEX_elem> cs_D1, cs_D2;
    NTL::Vec<FpEX_elem> cs_Xs16;
    NTL::Vec<FpEX_elem> cs_N1, cs_N2;
    FpEX_elem n11,n12,n13,n14,n15,n16,n17,n18;
    FpEX_elem n21,n22,n23,n24,n25,n26,n27,n28,n29;
    {
        NTL::SetCoeff(f1,8,FpE_elem(16));
        NTL::SetCoeff(f1,0,FpE_elem(-16));
        NTL::SetCoeff(f2,4,FpE_elem(-16));
        NTL::SetCoeff(f3,8,FpE_elem(1));
        NTL::SetCoeff(f4,4,FpE_elem(-1));
        cs_Xs16.SetLength(13);
        cs_Xs16[0] = f4;
        cs_Xs16[4] = f3;
        cs_Xs16[8] = f2;
        cs_Xs16[12] = f1;

    
        SetCoeff(n11,4,FpE_elem(1));
        SetCoeff(n12,5,FpE_elem(-1));
        SetCoeff(n13,6,FpE_elem(1));
        SetCoeff(n14,7,FpE_elem(-1));
        SetCoeff(n15,4,FpE_elem(16));
        SetCoeff(n16,5,FpE_elem(-16));
        SetCoeff(n17,6,FpE_elem(16));
        SetCoeff(n18,7,FpE_elem(-16));

        SetCoeff(n21,3,FpE_elem(-1));
        SetCoeff(n22,7,FpE_elem(1));
        SetCoeff(n22,3,FpE_elem(14));
        SetCoeff(n22,2,FpE_elem(-4));
        SetCoeff(n23,1,FpE_elem(-16));
        SetCoeff(n23,2,FpE_elem(48));
        SetCoeff(n23,3,FpE_elem(-120));
        SetCoeff(n23,6,FpE_elem(4));
        SetCoeff(n23,7,FpE_elem(-14));
        SetCoeff(n24,0,FpE_elem(-64));
        SetCoeff(n24,1,FpE_elem(160));
        SetCoeff(n24,2,FpE_elem(-368));
        SetCoeff(n24,3,FpE_elem(728));
        SetCoeff(n24,5,FpE_elem(16));
        SetCoeff(n24,6,FpE_elem(-48));
        SetCoeff(n24,7,FpE_elem(120));
        SetCoeff(n25,7,FpE_elem(-736));
        SetCoeff(n25,6,FpE_elem(384));
        SetCoeff(n25,5,FpE_elem(-192));
        SetCoeff(n25,4,FpE_elem(128));
        SetCoeff(n25,3,FpE_elem(-3328));
        SetCoeff(n25,2,FpE_elem(1984));
        SetCoeff(n25,1,FpE_elem(-1088));
        SetCoeff(n25,0,FpE_elem(512));
        SetCoeff(n26,7,FpE_elem(3200));
        SetCoeff(n26,6,FpE_elem(-1792));
        SetCoeff(n26,5,FpE_elem(768));
        SetCoeff(n26,3,FpE_elem(11648));
        SetCoeff(n26,2,FpE_elem(-7936));
        SetCoeff(n26,1,FpE_elem(5120));
        SetCoeff(n26,0,FpE_elem(-3072));
        SetCoeff(n27,7,FpE_elem(-12288));
        SetCoeff(n27,6,FpE_elem(9216));
        SetCoeff(n27,5,FpE_elem(-7168));
        SetCoeff(n27,4,FpE_elem(6144));
        SetCoeff(n27,3,FpE_elem(-30720));
        SetCoeff(n27,2,FpE_elem(23552));
        SetCoeff(n27,1,FpE_elem(-17408));
        SetCoeff(n27,0,FpE_elem(12288));
        SetCoeff(n28,7,FpE_elem(24576));
        SetCoeff(n28,6,FpE_elem(-16384));
        SetCoeff(n28,5,FpE_elem(8192));
        SetCoeff(n28,3,FpE_elem(57344));
        SetCoeff(n28,2,FpE_elem(-49152));
        SetCoeff(n28,1,FpE_elem(40960));
        SetCoeff(n28,0,FpE_elem(-32768));
        SetCoeff(n29,7,FpE_elem(-65536));
        SetCoeff(n29,6,FpE_elem(65536));
        SetCoeff(n29,5,FpE_elem(-65536));
        SetCoeff(n29,4,FpE_elem(65536));
        SetCoeff(n29,3,FpE_elem(-65536));
        SetCoeff(n29,2,FpE_elem(65536));
        SetCoeff(n29,1,FpE_elem(-65536));
        SetCoeff(n29,0,FpE_elem(65536));
    
        cs_N1.SetLength(12);
        cs_N1[0] = n11;
        cs_N1[1] = n12;
        cs_N1[2] = n13;
        cs_N1[3] = n14;
        cs_N1[8] = n15;
        cs_N1[9] = n16;
        cs_N1[10] = n17;
        cs_N1[11] = n18;

        cs_N2.SetLength(33);
        cs_N2[0] = n21;
        cs_N2[4] = n22;
        cs_N2[8] = n23;
        cs_N2[12] = n24;
        cs_N2[16] = n25;
        cs_N2[20] = n26;
        cs_N2[24] = n27;
        cs_N2[28] = n28;
        cs_N2[32] = n29;
        cs_D1.SetLength(13);
        cs_D1[12] = FpEX_elem(4);
        cs_D2.SetLength(17);
        cs_D2[16] = FpEX_elem(8);
        cs_D2[8] = FpEX_elem(1);
    }
    ZZ_pEXY Xs16 = ZZ_pEXY(cs_Xs16);
    ZZ_pEXY D1 = ZZ_pEXY(cs_D1);
    ZZ_pEXY D2 = ZZ_pEXY(cs_D2);
    ZZ_pEXY N1 = ZZ_pEXY(cs_N1);
    ZZ_pEXY N2 = ZZ_pEXY(cs_N2);

    // the two main 16 subgroup
    // we compute their point on X0(16)
    ecp P1 = web.P16;
    ecp P2 = web.Q16;
    assert((16*P1).is_identity());
    assert((16*P2).is_identity());
    assert(!((8*P1).is_identity()));
    assert(!((8*P2).is_identity()));
    assert(!((8*(P1-P2)).is_identity()));
    

    t = clock();
    clock_t isog_time = clock() - t;
    std::vector<std::pair<ecp,std::pair<int, int>>> kerGens1, kerGens2;
    // computing the w016 invariants for the 2 main subgroups
    Fp2 w0161,w0162;
    {
        kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P1, std::pair<int, int>(2,4)));
        kerGens2.push_back(std::pair<ecp,std::pair<int, int>>(P2, std::pair<int, int>(2,4)));
        isog_chain phi_P1(kerGens1);
        isog_chain phi_P2(kerGens2);
        auto E1 = phi_P1.get_codomain();
        auto jE1 = E1.j_invariant();
        auto E2 = phi_P2.get_codomain();
        assert(E == phi_P1.get_domain());
        assert(E == phi_P2.get_domain());
        auto jE2 = E2.j_invariant();
        
        F1 = I0 - jE1*J0;
        F2 = I0 - jE2*J0;

        FpEX_elem f1 = GCD(G, F1);
        FpEX_elem f2 = GCD(G, F2);

        if (NTL::deg(f1) == 0 || NTL::deg(f2) == 0) { 
            assert(0);
        }
    
        std::vector<Fp2> r1 = {};
        std::vector<Fp2> r2 = {};
        if (NTL::deg(f1)==1) {
            r1.push_back(-f1[0]);
        } 
        else if (NTL::deg(f1) > 1) {
            ec EE2 = phi_P1.get_middle_codomain(0);
            ec EE4 = phi_P1.get_middle_codomain(1);
            ec EE8 = phi_P1.get_middle_codomain(2);
            // alternate try
            // computing the t-invariants
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E1.disc()/EE8.disc();
            auto v1 = FindRoots(f1);
            bool found = false; 
            int i1 = 0;
            while (i1 < NTL::deg(f1)) {
                Fp2 w016 = v1[i1];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r1.push_back(v1[i1]);  }
                i1++;
            }
        }
        if (NTL::deg(f2)==1) {
            r2.push_back(-f2[0]);
        }
        else if (NTL::deg(f2) > 1) {

            ec EE2 = phi_P2.get_middle_codomain(0);
            ec EE4 = phi_P2.get_middle_codomain(1);
            ec EE8 = phi_P2.get_middle_codomain(2);
            // alternate try
            // computing the t-invariants
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E2.disc()/EE8.disc();

            bool found = false; 
            int i2 = 0;
            auto v2 = FindRoots(f2);
            while (i2 < NTL::deg(f2)) {
                // we compute the X0(4) invariants associated 
                Fp2 w016 = v2[i2];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r2.push_back(v2[i2]);  }
                i2++;
            }
        }

        if (r1.size() == 1) {
                w0161 = r1[0];
        }
        else {return SpecialEnumFastWeber(web, Fexts);}
        if (r2.size() == 1) {
                w0162 = r2[0]; 
        }
        else {return SpecialEnumFastWeber(web, Fexts);}

    }
    
    // now we compute the 16 others
    std::vector<std::pair<Integer,Integer>> coeff_list_16 = {};
    {
        coeff_list_16.push_back( {Integer(1),Integer(1)}  );
        coeff_list_16.push_back( {Integer(1),Integer(3)}  );
        coeff_list_16.push_back( {Integer(1),Integer(5)}  );
        coeff_list_16.push_back( {Integer(1),Integer(7)}  );
        coeff_list_16.push_back( {Integer(1),Integer(9)}  );
        coeff_list_16.push_back( {Integer(1),Integer(11)}  );
        coeff_list_16.push_back( {Integer(1),Integer(13)}  );
        coeff_list_16.push_back( {Integer(1),Integer(15)}  );
        coeff_list_16.push_back( {Integer(0),Integer(1)}  );
        coeff_list_16.push_back( {Integer(2),Integer(1)}  );
        coeff_list_16.push_back( {Integer(4),Integer(1)}  );
        coeff_list_16.push_back( {Integer(6),Integer(1)}  );
        coeff_list_16.push_back( {Integer(8),Integer(1)}  );
        coeff_list_16.push_back( {Integer(10),Integer(1)}  );
        coeff_list_16.push_back( {Integer(12),Integer(1)}  );
        coeff_list_16.push_back( {Integer(14),Integer(1)}  );

    }
    std::vector<Fp2> w016_list = {}; 
    // computing the w016 invariants for the remaing subgroups
    for (auto coeff : coeff_list_16 ) {
        clock_t t_loc = clock();
        ecp P = coeff.first * P1 + coeff.second * P2;
        kerGens1.pop_back();
        kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P, std::pair<int, int>(2,4)));
        isog_chain phi_P1(kerGens1);
        auto E1 = phi_P1.get_codomain();
        auto jE1 = E1.j_invariant();
        assert(E == phi_P1.get_domain());
        isog_time += (clock() - t_loc);
        F1 = I0 - jE1*J0;
        FpEX_elem f1 = GCD(G, F1);
        if (NTL::deg(f1) == 1) {
            w016_list.push_back(-f1[0]);
            Fp2 w016 = -f1[0];
        }
        else {
            auto v1 = FindRoots(f1);

            // there are several possibilies we need to find which one is the correct
            ec EE2 = phi_P1.get_middle_codomain(0);
            ec EE4 = phi_P1.get_middle_codomain(1);
            ec EE8 = phi_P1.get_middle_codomain(2);
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E1.disc()/EE8.disc();

            bool found = false; 
            int i1 = 0;
            std::vector<Fp2> r1 = {};
            while (i1 < NTL::deg(f1)) {
                // we compute the X0(4) invariants associated 
                Fp2 w016 = v1[i1];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r1.push_back(v1[i1]);  }
                i1++;
            }
            if (r1.size() == 1) {
                w016_list.push_back(r1[0]);
            }
            else {
                // std::cout << r1.size() << " " << r1[0] << " " << r1[1] << "\n";
                // assert(0);
                return SpecialEnumFastWeber(web, Fexts);
            }
        }
    }


    t = clock();
    // Now that we have computed all relevant w016 invariants we can proceed to the computation of the xs16 invariants
    // For that we precompute some informations related to the two main ones.    
    // First precomputing the powers
    std::vector<Fp2> pow1 = get_powers(w0161, 34);
    std::vector<Fp2> pow2 = get_powers(w0162, 34);

    // Doing the bivariate part
    Fp2X biv1 = Fp2X(0); 
    Fp2X biv2 = Fp2X(0);
    for (unsigned int i = 0 ; i <= R1.dY ; i++ ) {
        biv1 += pow1[i]*R1.coeffs[i];
        biv2 += pow2[i]*R1.coeffs[i];
    }

    // Doing the trivariate ones
    std::vector<std::vector<Fp2>> resultant_ws16(24);


    
    NTL::Vec<Fp2X> triv1_vec,triv2_vec;
    triv1_vec.SetLength(R2.dZ + 1);
    triv2_vec.SetLength(R2.dZ + 1);
    for (unsigned int i = 0; i <= R2.dZ ; i++) {
        triv1_vec[i] = Fp2X(0); 
        triv2_vec[i] = Fp2X(0); 
        for (unsigned int j = 0; j <= R2.coeffs[i].dY; j++ ) {
            triv1_vec[i] += pow1[j] * R2.coeffs[i].coeffs[j];
            triv2_vec[i] += pow2[j] * R2.coeffs[i].coeffs[j];
        }
    }

    for (int id = 0; id < 16; id++) {
        std::vector<Fp2> powi = get_powers(w016_list[id], R2.dZ+1);
        Fp2X triv1j = Fp2X(0);
        Fp2X triv2j = Fp2X(0);
        
        for (unsigned int i = 0; i <= R2.dZ ; i++) {
            triv1j += powi[i] * triv1_vec[i]; 
            if (id < 8) {
                triv2j += powi[i] * triv2_vec[i];
            }
        }
        if (id < 8) {
            resultant_ws16[id] = {w0161,w016_list[id],CommonRootTwoResultants({biv1, triv1j})};
            resultant_ws16[id+8] = {w0162, w016_list[id],CommonRootTwoResultants({biv2, triv2j})};
        } 
        if (id >=8) {
            resultant_ws16[id+8] = {w0161,w016_list[id],CommonRootTwoResultants({biv1, triv1j})};
        }
    }

    t = clock();
    std::vector<Fp2> w3_list = {};
    // Computing all third roots
    for (auto yy : resultant_ws16) {
        Fp2 y = yy[2];

        auto d1 = EvaluateBivariate(D1,y);
        auto d2 = EvaluateBivariate(D2,y);
        auto n1 = EvaluateBivariate(N1,y);
        auto n2 = EvaluateBivariate(N2,y);

        Fp2X F1 = n1 - yy[0]*d1;
        Fp2X F2 = n2 - yy[1]*d2;
        Fp2X g = GCD(F1,F2);
        
        // Seems to be true most of the time but not always
        if (NTL::deg(g) == 1) {
            w3_list.push_back(_getWeberThirdPowerFromRoot(-g[0],y));
        }
        else {
            Fp2X F = EvaluateBivariate(Xs16,y);
            F /= LeadCoeff(F);
            Fp2X gg = GCD(g,F);
            if (NTL::deg(gg) != 1) {
                std::cout << "problematic curve = " << jE << " deg = " << NTL::deg(gg) << " " << yy[0] << " " << yy[1] << "\n";

            }
        }
    }


    // computing all t-invariants
    std::vector<Fp2> tinv_list = { _getTInvariant(8*web.Q16, E), _getTInvariant(8*web.P16, E), _getTInvariant(8*(web.P16+web.Q16), E) };

    // computing gamma2-invariants
    std::vector<std::pair<ecp,ecp>> ord3_bases = {{web.P3, web.Q3}, {web.P3, web.Q3 + web.P3}, {web.P3, 2*web.Q3 + web.P3}};
    std::vector<Fp2> gamma2_list = {};
    for (auto bas3 : ord3_bases) {
        gamma2_list.push_back( _getGammaTwoInvariant( {bas3.first, bas3.second, bas3.first + bas3.second, bas3.first + 2*bas3.second}, E ));
    }



    // constructing the coefficients
    // order 2 
    std::vector<std::pair<Integer,Integer>> order_2 = { {Integer(0),Integer(1)}, {Integer(1),Integer(0)}, {Integer(1), Integer(1)} };
    // order 3
    std::vector<std::vector<std::pair<Integer,Integer>>> order_3 = { { {Integer(1), Integer(0)}, {Integer(0), Integer(1)} }, { {Integer(1),Integer(0)}, {Integer(1), Integer(1)} }, { {Integer(1), Integer(0)}, {Integer(1), Integer(2)} } };

    // // pairs of relevant order_16 generators
    std::vector<std::vector<std::pair<Integer,Integer>>> order_16_pairs = {};
    {
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(3)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(5)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(7)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(9)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(11)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(13)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(15)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(3)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(5)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(7)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(9)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(11)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(13)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(15)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(0),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(2),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(4),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(6),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(8),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(10),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(12),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(14),Integer(1)} } );
    }

    // now we combine everything together 
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) { 
            // first we combine t-inv and gamma2 together
            Fp2 w8 = _getWeberEighthPower(gamma2_list[j], tinv_list[i]);
            Fp2X W1, W2;
            NTL::SetCoeff(W2, 8, FpE_elem(1));
            NTL::SetCoeff(W2, 0, FpE_elem(-w8));
            
            for (size_t k =0; k < 8; k++) {
                // if (!(((order_2[i].first == (order_16_pairs[k][0].first % 2)) && (order_2[i].second == (order_16_pairs[k][0].second % 2))) || ((order_2[i].first == (order_16_pairs[k][1].first % 2)) && (order_2[i].second == (order_16_pairs[k][1].second % 2))))) {
                //     result.push_back( { {order_2[i]}, order_3[j], order_16_pairs[k] } );
                // }
                NTL::SetCoeff(W1, 3, FpE_elem(1));
                NTL::SetCoeff(W1, 0, FpE_elem(-w3_list[8*i + k]));
                
            
                auto g = NTL::GCD(W1,W2);
                assert(NTL::deg(g) <= 1);
                if (NTL::deg(g)==0) {
                    // std::cout << "problem in weber for curve " << E.j_invariant()    << " with w3 = " << w3 << " w8 = " << w8 << "\n"; 
                    assert(0);
                }
                else {
                    assert(NTL::deg(g) == 1);
                    Fp2 w = Fp2(-g[0]);
                    result.push_back( {w, {{order_2[i]}, order_3[j], order_16_pairs[8*i + k]}} );
                }
                
            }
        }
    }
    // std::cout << "total time for fast weber enum " << (double) (clock() - tot)/CLOCKS_PER_SEC << "\n";
    assert( power(w3_list[0]/w3_list[1],8) == 1 );
    assert(result[0].first/result[1].first == power(w3_list[0]/w3_list[1],3));

    assert(power(result[0].first/result[1].first,2) == -result[0].first/result[2].first);
    assert(result[0].first/result[2].first == result[1].first/result[3].first);
    assert(result[0].first/result[1].first == -result[48].first/result[49].first);
    assert(result[0].first/result[2].first == -result[48].first/result[50].first);

    return result;

}

weber_enum EnumerateAllWeberFast(const weber_bas &web, const std::map<unsigned,Fp2k> &Fexts) {


    clock_t tot = clock(); (void) tot;
    clock_t t = clock();
    ec E = web.P3.curve();
    auto jE = E.j_invariant();
    

    // first check if we are not in some special case
    if (jE ==Fp2(1728) || jE == Fp2(287496) || jE == Fp2(0)) {
        return SpecialEnumFastWeber(web, Fexts);
    }

    // First we start by enumerating the X016 points that because this where things can fail most dramatically 
    std::vector<std::pair<Fp2, std::vector<std::vector<std::pair<Integer,Integer>>>>> result(72);

    assert(web.P3.curve() == web.P16.curve());

    
    // Initalizing all needed polynomials
    auto R1 = GetBivariateResultant();
    auto R2 = GetTrivariateResultant();
    auto G0 = _getG0();
    auto H0 = _getH0();
    auto I0 = _getI0();
    auto J0 = _getJ0();
    Fp2X G, F1, F2;
    
    G = G0 - jE*H0;
    FpEX_elem f1,f2,f3,f4;
    NTL::Vec<FpEX_elem> cs_D1, cs_D2;
    NTL::Vec<FpEX_elem> cs_Xs16;
    NTL::Vec<FpEX_elem> cs_N1, cs_N2;
    FpEX_elem n11,n12,n13,n14,n15,n16,n17,n18;
    FpEX_elem n21,n22,n23,n24,n25,n26,n27,n28,n29;
    {
        NTL::SetCoeff(f1,8,FpE_elem(16));
        NTL::SetCoeff(f1,0,FpE_elem(-16));
        NTL::SetCoeff(f2,4,FpE_elem(-16));
        NTL::SetCoeff(f3,8,FpE_elem(1));
        NTL::SetCoeff(f4,4,FpE_elem(-1));
        cs_Xs16.SetLength(13);
        cs_Xs16[0] = f4;
        cs_Xs16[4] = f3;
        cs_Xs16[8] = f2;
        cs_Xs16[12] = f1;

    
        SetCoeff(n11,4,FpE_elem(1));
        SetCoeff(n12,5,FpE_elem(-1));
        SetCoeff(n13,6,FpE_elem(1));
        SetCoeff(n14,7,FpE_elem(-1));
        SetCoeff(n15,4,FpE_elem(16));
        SetCoeff(n16,5,FpE_elem(-16));
        SetCoeff(n17,6,FpE_elem(16));
        SetCoeff(n18,7,FpE_elem(-16));

        SetCoeff(n21,3,FpE_elem(-1));
        SetCoeff(n22,7,FpE_elem(1));
        SetCoeff(n22,3,FpE_elem(14));
        SetCoeff(n22,2,FpE_elem(-4));
        SetCoeff(n23,1,FpE_elem(-16));
        SetCoeff(n23,2,FpE_elem(48));
        SetCoeff(n23,3,FpE_elem(-120));
        SetCoeff(n23,6,FpE_elem(4));
        SetCoeff(n23,7,FpE_elem(-14));
        SetCoeff(n24,0,FpE_elem(-64));
        SetCoeff(n24,1,FpE_elem(160));
        SetCoeff(n24,2,FpE_elem(-368));
        SetCoeff(n24,3,FpE_elem(728));
        SetCoeff(n24,5,FpE_elem(16));
        SetCoeff(n24,6,FpE_elem(-48));
        SetCoeff(n24,7,FpE_elem(120));
        SetCoeff(n25,7,FpE_elem(-736));
        SetCoeff(n25,6,FpE_elem(384));
        SetCoeff(n25,5,FpE_elem(-192));
        SetCoeff(n25,4,FpE_elem(128));
        SetCoeff(n25,3,FpE_elem(-3328));
        SetCoeff(n25,2,FpE_elem(1984));
        SetCoeff(n25,1,FpE_elem(-1088));
        SetCoeff(n25,0,FpE_elem(512));
        SetCoeff(n26,7,FpE_elem(3200));
        SetCoeff(n26,6,FpE_elem(-1792));
        SetCoeff(n26,5,FpE_elem(768));
        SetCoeff(n26,3,FpE_elem(11648));
        SetCoeff(n26,2,FpE_elem(-7936));
        SetCoeff(n26,1,FpE_elem(5120));
        SetCoeff(n26,0,FpE_elem(-3072));
        SetCoeff(n27,7,FpE_elem(-12288));
        SetCoeff(n27,6,FpE_elem(9216));
        SetCoeff(n27,5,FpE_elem(-7168));
        SetCoeff(n27,4,FpE_elem(6144));
        SetCoeff(n27,3,FpE_elem(-30720));
        SetCoeff(n27,2,FpE_elem(23552));
        SetCoeff(n27,1,FpE_elem(-17408));
        SetCoeff(n27,0,FpE_elem(12288));
        SetCoeff(n28,7,FpE_elem(24576));
        SetCoeff(n28,6,FpE_elem(-16384));
        SetCoeff(n28,5,FpE_elem(8192));
        SetCoeff(n28,3,FpE_elem(57344));
        SetCoeff(n28,2,FpE_elem(-49152));
        SetCoeff(n28,1,FpE_elem(40960));
        SetCoeff(n28,0,FpE_elem(-32768));
        SetCoeff(n29,7,FpE_elem(-65536));
        SetCoeff(n29,6,FpE_elem(65536));
        SetCoeff(n29,5,FpE_elem(-65536));
        SetCoeff(n29,4,FpE_elem(65536));
        SetCoeff(n29,3,FpE_elem(-65536));
        SetCoeff(n29,2,FpE_elem(65536));
        SetCoeff(n29,1,FpE_elem(-65536));
        SetCoeff(n29,0,FpE_elem(65536));
    
        cs_N1.SetLength(12);
        cs_N1[0] = n11;
        cs_N1[1] = n12;
        cs_N1[2] = n13;
        cs_N1[3] = n14;
        cs_N1[8] = n15;
        cs_N1[9] = n16;
        cs_N1[10] = n17;
        cs_N1[11] = n18;

        cs_N2.SetLength(33);
        cs_N2[0] = n21;
        cs_N2[4] = n22;
        cs_N2[8] = n23;
        cs_N2[12] = n24;
        cs_N2[16] = n25;
        cs_N2[20] = n26;
        cs_N2[24] = n27;
        cs_N2[28] = n28;
        cs_N2[32] = n29;
        cs_D1.SetLength(13);
        cs_D1[12] = FpEX_elem(4);
        cs_D2.SetLength(17);
        cs_D2[16] = FpEX_elem(8);
        cs_D2[8] = FpEX_elem(1);
    }
    ZZ_pEXY Xs16 = ZZ_pEXY(cs_Xs16);
    ZZ_pEXY D1 = ZZ_pEXY(cs_D1);
    ZZ_pEXY D2 = ZZ_pEXY(cs_D2);
    ZZ_pEXY N1 = ZZ_pEXY(cs_N1);
    ZZ_pEXY N2 = ZZ_pEXY(cs_N2);

    // the two main 16 subgroup
    // we compute their point on X0(16)
    ecp P1 = web.P16;
    ecp P2 = web.Q16;
    assert((16*P1).is_identity());
    assert((16*P2).is_identity());
    assert(!((8*P1).is_identity()));
    assert(!((8*P2).is_identity()));
    assert(!((8*(P1-P2)).is_identity()));
    

    t = clock();
    clock_t isog_time = clock() - t;
    std::vector<std::pair<ecp,std::pair<int, int>>> kerGens1, kerGens2;
    // computing the w016 invariants for the 2 main subgroups
    Fp2 w0161,w0162;
    {
        kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P1, std::pair<int, int>(2,4)));
        kerGens2.push_back(std::pair<ecp,std::pair<int, int>>(P2, std::pair<int, int>(2,4)));
        isog_chain phi_P1(kerGens1);
        isog_chain phi_P2(kerGens2);
        auto E1 = phi_P1.get_codomain();
        auto jE1 = E1.j_invariant();
        auto E2 = phi_P2.get_codomain();
        assert(E == phi_P1.get_domain());
        assert(E == phi_P2.get_domain());
        auto jE2 = E2.j_invariant();
        
        F1 = I0 - jE1*J0;
        F2 = I0 - jE2*J0;

        FpEX_elem f1 = GCD(G, F1);
        FpEX_elem f2 = GCD(G, F2);

        if (NTL::deg(f1) == 0 || NTL::deg(f2) == 0) { 
            assert(0);
        }
    
        std::vector<Fp2> r1 = {};
        std::vector<Fp2> r2 = {};
        if (NTL::deg(f1)==1) {
            r1.push_back(-f1[0]);
        } 
        else if (NTL::deg(f1) > 1) {
            ec EE2 = phi_P1.get_middle_codomain(0);
            ec EE4 = phi_P1.get_middle_codomain(1);
            ec EE8 = phi_P1.get_middle_codomain(2);
            // alternate try
            // computing the t-invariants
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E1.disc()/EE8.disc();
            auto v1 = FindRoots(f1);
            bool found = false; 
            int i1 = 0;
            while (i1 < NTL::deg(f1)) {
                Fp2 w016 = v1[i1];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r1.push_back(v1[i1]);  }
                i1++;
            }
        }
        if (NTL::deg(f2)==1) {
            r2.push_back(-f2[0]);
        }
        else if (NTL::deg(f2) > 1) {

            ec EE2 = phi_P2.get_middle_codomain(0);
            ec EE4 = phi_P2.get_middle_codomain(1);
            ec EE8 = phi_P2.get_middle_codomain(2);
            // alternate try
            // computing the t-invariants
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E2.disc()/EE8.disc();

            bool found = false; 
            int i2 = 0;
            auto v2 = FindRoots(f2);
            while (i2 < NTL::deg(f2)) {
                // we compute the X0(4) invariants associated 
                Fp2 w016 = v2[i2];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r2.push_back(v2[i2]);  }
                i2++;
            }
        }

        if (r1.size() == 1) {
                w0161 = r1[0];
        }
        else {return SpecialEnumFastWeber(web, Fexts);}
        if (r2.size() == 1) {
                w0162 = r2[0]; 
        }
        else {return SpecialEnumFastWeber(web, Fexts);}

    }
    
    // now we compute the 16 others
    std::vector<std::pair<Integer,Integer>> coeff_list_16 = {};
    {
        coeff_list_16.push_back( {Integer(1),Integer(1)}  );
        coeff_list_16.push_back( {Integer(1),Integer(3)}  );
    }
    std::vector<Fp2> w016_list = {}; 
    // computing the w016 invariants for the remaing subgroups
    for (auto coeff : coeff_list_16 ) {
        clock_t t_loc = clock();
        ecp P = coeff.first * P1 + coeff.second * P2;
        kerGens1.pop_back();
        kerGens1.push_back(std::pair<ecp,std::pair<int, int>>(P, std::pair<int, int>(2,4)));
        isog_chain phi_P1(kerGens1);
        auto E1 = phi_P1.get_codomain();
        auto jE1 = E1.j_invariant();
        assert(E == phi_P1.get_domain());
        isog_time += (clock() - t_loc);
        F1 = I0 - jE1*J0;
        FpEX_elem f1 = GCD(G, F1);
        if (NTL::deg(f1) == 1) {
            w016_list.push_back(-f1[0]);
            Fp2 w016 = -f1[0];
        }
        else {
            auto v1 = FindRoots(f1);
            // There are several possibilies we need to find which one is the correct
            ec EE2 = phi_P1.get_middle_codomain(0);
            ec EE4 = phi_P1.get_middle_codomain(1);
            ec EE8 = phi_P1.get_middle_codomain(2);
            Fp2 t1 = EE2.disc()/E.disc();
            Fp2 t2 = EE4.disc()/EE2.disc();
            Fp2 t3 = EE8.disc()/EE4.disc();
            Fp2 t4 = E1.disc()/EE8.disc();

            bool found = false; 
            int i1 = 0;
            std::vector<Fp2> r1 = {};
            while (i1 < NTL::deg(f1)) {
                // we compute the X0(4) invariants associated 
                Fp2 w016 = v1[i1];
                Fp2 w041 = NTL::power(w016,4) / ( NTL::power(w016,3) + 6 * NTL::power(w016,2) + 16 * w016 + 16 );
                Fp2 w042 = ( NTL::power(w016,4) + 8* NTL::power(w016,3) + 24 * NTL::power(w016,2) + 32 * w016);
                Fp2 w021 = NTL::power(w041,2) / (w041 + 16) ;
                Fp2 w022 = (NTL::power(w041,2) + 16 * w041);
                Fp2 w023 = NTL::power(w042,2) / (w042 + 16) ;
                Fp2 w024 = (NTL::power(w042,2) + 16 * w042);
                found = t1 * w021 == 4096 &&
                    t2 * w022 == 4096 &&
                    t3 * w023 == 4096 && 
                    t4 * w024 == 4096; 
                if (found) {r1.push_back(v1[i1]);  }
                i1++;
            }
            if (r1.size() == 1) {
                w016_list.push_back(r1[0]);
            }
            else {
                return SpecialEnumFastWeber(web, Fexts);
            }
        }
    }

    t = clock();
    // Now that we have computed all relevant w016 invariants we can proceed to the computation of the xs16 invariants
    // For that we precompute some informations related to the two main ones.    
    // First precomputing the powers
    std::vector<Fp2> pow1 = get_powers(w0161, 34);
    std::vector<Fp2> pow2 = get_powers(w0162, 34);

    // Doing the bivariate part
    Fp2X biv1 = Fp2X(0); 
    Fp2X biv2 = Fp2X(0);
    for (unsigned int i = 0 ; i <= R1.dY ; i++ ) {
        biv1 += pow1[i]*R1.coeffs[i];
        biv2 += pow2[i]*R1.coeffs[i];
    }

    // Doing the trivariate ones
    std::vector<std::vector<Fp2>> resultant_ws16(4);


    
    NTL::Vec<Fp2X> triv1_vec,triv2_vec;
    triv1_vec.SetLength(R2.dZ + 1);
    triv2_vec.SetLength(R2.dZ + 1);
    for (unsigned int i = 0; i <= R2.dZ ; i++) {
        triv1_vec[i] = Fp2X(0); 
        triv2_vec[i] = Fp2X(0); 
        for (unsigned int j = 0; j <= R2.coeffs[i].dY; j++ ) {
            triv1_vec[i] += pow1[j] * R2.coeffs[i].coeffs[j];
            triv2_vec[i] += pow2[j] * R2.coeffs[i].coeffs[j];
        }
    }

    w016_list.push_back(w0162);

    for (int id = 0; id < 3; id++) {
        std::vector<Fp2> powi = get_powers(w016_list[id], R2.dZ+1);
        Fp2X triv1j = Fp2X(0);
        Fp2X triv2j = Fp2X(0);
        
        for (unsigned int i = 0; i <= R2.dZ ; i++) {
            triv1j += powi[i] * triv1_vec[i];
            if (id == 0) {
                triv2j += powi[i] * triv2_vec[i];
            }
        }
        if (id == 0) {
            resultant_ws16[0] = {w0161,w016_list[id],CommonRootTwoResultants({biv1, triv1j})};
            resultant_ws16[1] = {w0162, w016_list[id],CommonRootTwoResultants({biv2, triv2j})};
        }
        else if (id == 1) {
            resultant_ws16[3] = {w0161,w016_list[id],CommonRootTwoResultants({biv1, triv1j})};
        } 
        else {
            assert(id==2);
            resultant_ws16[2] = {w0161, w016_list[id],CommonRootTwoResultants({biv1, triv1j})};
        }
    }


    t = clock();
    std::vector<Fp2> w3_list = {};
    // Computing all third roots
    for (auto yy : resultant_ws16) {
        Fp2 y = yy[2];
        auto d1 = EvaluateBivariate(D1,y);
        auto d2 = EvaluateBivariate(D2,y);
        auto n1 = EvaluateBivariate(N1,y);
        auto n2 = EvaluateBivariate(N2,y);

        Fp2X F1 = n1 - yy[0]*d1;
        Fp2X F2 = n2 - yy[1]*d2;
        Fp2X g = GCD(F1,F2);
        // assert(NTL::deg(g) > 0);
        if (NTL::deg(g) == 0) {
            // std::cout << "g was empty \n";
            return SpecialEnumFastWeber(web, Fexts);
        }
        
        // seems to be true most of the time but not always
        if (NTL::deg(g) == 1) {
            w3_list.push_back(_getWeberThirdPowerFromRoot(-g[0],y));
        }
        else {
            Fp2X F = EvaluateBivariate(Xs16,y);
            F /= LeadCoeff(F);
            Fp2X gg = GCD(g,F);
            if (NTL::deg(gg) != 1) {
                // std::cout << "problematic curve = " << jE << " deg = " << NTL::deg(g) << NTL::deg(gg) << " " << yy[0] << " " << yy[1] << "\n";
            }
            assert(NTL::deg(gg)==1);
            w3_list.push_back(_getWeberThirdPowerFromRoot(-gg[0],y));
        }
    }
    Fp2 zeta_80 = w3_list[0]/w3_list[3];
    Fp2 zeta_8 = power(zeta_80,5);
    Fp2 zeta_8bis = power(zeta_80,7);
    assert(power(zeta_8, 8) == Fp2(1));
    Fp2 zeta_4 = -power(zeta_8,2);
    std::vector<Fp2> z8_list = { zeta_8, zeta_8bis, - zeta_8 };
    std::vector<Fp2> z4_list = { zeta_4, zeta_4, - zeta_4 };

    // computing all t-invariants
    std::vector<Fp2> tinv_list = { _getTInvariant(8*web.Q16, E), _getTInvariant(8*web.P16, E), _getTInvariant(8*(web.P16+web.Q16), E) };

    // computing gamma2-invariants
    std::vector<std::pair<ecp,ecp>> ord3_bases = {{web.P3, web.Q3}, {web.P3, web.Q3 + web.P3}, {web.P3, 2*web.Q3 + web.P3}};
    std::vector<Fp2> gamma2_list = {};
    for (auto bas3 : ord3_bases) {
        gamma2_list.push_back( _getGammaTwoInvariant( {bas3.first, bas3.second, bas3.first + bas3.second, bas3.first + 2*bas3.second}, E ));
    }



    // Constructing the coefficients
    // order 2 
    std::vector<std::pair<Integer,Integer>> order_2 = { {Integer(0),Integer(1)}, {Integer(1),Integer(0)}, {Integer(1), Integer(1)} };
    // order 3
    std::vector<std::vector<std::pair<Integer,Integer>>> order_3 = { { {Integer(1), Integer(0)}, {Integer(0), Integer(1)} }, { {Integer(1),Integer(0)}, {Integer(1), Integer(1)} }, { {Integer(1), Integer(0)}, {Integer(1), Integer(2)} } };

    // pairs of relevant order_16 generators
    std::vector<std::vector<std::pair<Integer,Integer>>> order_16_pairs = {};
    {
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(3)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(5)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(7)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(9)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(11)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(13)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(1),Integer(15)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(3)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(5)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(7)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(9)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(11)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(13)} } );
        order_16_pairs.push_back( { {Integer(0),Integer(1)}, {Integer(1),Integer(15)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(0),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(2),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(4),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(6),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(8),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(10),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(12),Integer(1)} } );
        order_16_pairs.push_back( { {Integer(1),Integer(0)}, {Integer(14),Integer(1)} } );
    }

    // now we combine everything together 
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < 3; j++) { 
            Fp2X W1,W2;
            // first we combine t-inv and gamma2 together
            Fp2 w8 = _getWeberEighthPower(gamma2_list[j], tinv_list[i]);
            
            NTL::SetCoeff(W2, 8, FpE_elem(1));
            NTL::SetCoeff(W2, 0, FpE_elem(-w8));
            NTL::SetCoeff(W1, 3, FpE_elem(1));
            NTL::SetCoeff(W1, 0, FpE_elem(-w3_list[i])); 
            auto g = NTL::GCD(W1,W2);
            assert(NTL::deg(g) <= 1);
            if (NTL::deg(g)==0) {
                // std::cout << "problem in weber for curve " << E.j_invariant()    << " with w3 = " << w3 << " w8 = " << w8 << "\n"; 
                assert(0);
            }
            else {
                assert(NTL::deg(g) == 1);
                Fp2 w = Fp2(-g[0]);
                result[24*i + 8*j] = {w, {{order_2[i]}, order_3[j], order_16_pairs[8*i]}};
                result[24*i + 8*j + 1] = {z8_list[i] * w, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 1]}};
                result[24*i + 8*j + 2] = {z4_list[i] * w, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 2]}};
                result[24*i + 8*j + 3] = {z4_list[i] * result[24*i + 8*j + 1].first, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 3]}};
                result[24*i + 8*j + 4] = {- w, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 4]}};
                result[24*i + 8*j + 5] = {- result[24*i + 8*j + 1].first, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 5]}};
                result[24*i + 8*j + 6] = {- result[24*i + 8*j + 2].first, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 6]}};
                result[24*i + 8*j + 7] = {- result[24*i + 8*j + 3].first, {{order_2[i]}, order_3[j], order_16_pairs[8*i + 7]}};
            }   


        }
    }

    return result;
}


