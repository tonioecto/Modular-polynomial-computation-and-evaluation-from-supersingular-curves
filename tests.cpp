#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_p.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pE.h>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <list>


#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"
#include "utils.hpp"
#include "quaternions.hpp"
#include "id2iso.hpp"
#include "quatlatenum.hpp"
#include "endring.hpp"
#include "modpoly.hpp"
#include "klpt.hpp"
#include "hashmap.hpp"
#include "multivariates.hpp"
#include "ordertojinvbigset.hpp"

#include <cassert>
#include <ctime>

void test_SpecialSupersingularEvaluation() {

    // checking p = 3 mod 4;

    // choice of prime
    Integer p(18131);
    // Integer p(1051);

    // init of field
    Fp_integer p_mod;
    NTL::conv(p_mod, p);
    Fp::init(p_mod);
    FpX f;
    SetCoeff(f, 2);
    // f[0] = Fp(7);
    // f[1] = Fp(18127);
    // f[2] = Fp(1);
    f[0] = 1;
    Fp2::init(f);

    // j-invariant
    size_t ell = 5;
    Integer ZZell(5);
    Fp j = Fp(12);

    // computation of the eval list
    std::vector<Fp> eval_points = {};
    eval_points.push_back(Fp(1));
    eval_points.push_back(j);
    Fp pow = j;
    for (size_t i = 2; i<=ell+1; i++) {
        pow = j * pow;
        eval_points.push_back(pow);
    }
    assert(eval_points.size() == ell + 2);
    auto Phi = SpecialSupersingularEvaluation( p, ZZell, eval_points);



    std::cout << " \n \nOutput of SpecialSupersingularEval for p = 3 mod 4 : \n";
    std::cout << Phi << "\n \n";

    std::cout << " \n \nOutput of Supersingular Evaluation : \n";
    std::cout << "[[11831] [17900] [13484] [13108] [9422] [4136] [1]]" << "\n \n";


    p=71909;
    NTL::conv(p_mod, p);
    Fp::init(p_mod);
    SetCoeff(f, 2);
    f[0] = Fp(3);
    Fp2::init(f);

    ell = 7;
    ZZell= 7;

    std::vector<Fp> new_eval_points = {};
    new_eval_points.push_back(Fp(1));
    new_eval_points.push_back(j);
    pow = j;
    for (size_t i = 2; i<=ell+1; i++) {
        pow = j * pow;
        new_eval_points.push_back(pow);
    }
    assert(new_eval_points.size() == ell + 2);
    auto new_Phi = SpecialSupersingularEvaluation( p, ZZell, new_eval_points);




    // checking p = 1 mod 4;
    std::cout << " \n \nOutput of SpecialSupersingularEval for p = 1 mod 4 : \n";
    std::cout << new_Phi << "\n \n";

    std::cout << " \n \nOutput of Supersingular Evaluation : \n";
    std::cout << "[[65003] [35779] [70505] [10525] [57521] [33681] [47728] [32793] [1]] \n \n";

}


void test_SpecialSupersingularEvaluationWeber() {

    // checking p = 3 mod 4;

    // choice of prime
    // Integer p(18131);
    Integer p(13003);
    // Integer p(823);
    // Integer p(479);
    // Integer p(337);
    // Integer p(277);
    // Integer p(59);

    // init of field
    Fp_integer p_mod;
    NTL::conv(p_mod, p);
    Fp::init(p_mod);
    FpX f;
    SetCoeff(f, 2);
    // f[0] = Fp(7);
    // f[1] = Fp(18127);
    // f[2] = Fp(1);
    f[0] = 1;
    Fp2::init(f);

    // j-invariant
    size_t ell = 8011;
    // size_t ell = 229;
    // size_t ell = 811;
    // size_t ell = 2741;
    // size_t ell = 5701;
    Integer ZZell(ell);
    Fp tau = Fp(2);

    // computation of the eval list
    std::vector<Fp> eval_points = {};
    eval_points.push_back(Fp(1));
    eval_points.push_back(tau);
    Fp pow = tau;
    for (size_t i = 2; i<=ell+1; i++) {
        pow = tau * pow;
        eval_points.push_back(pow);
    }
    assert(eval_points.size() == ell + 2);
    auto Phi = SpecialSupersingularEvaluationWeber( p, ZZell, eval_points);



    std::cout << " \n \nOutput of SpecialSupersingularEval for p = 3 mod 4 : \n";
    std::cout << Phi << "\n \n";

    std::cout << "\n\nSupposed output \n";
    std::cout << "[21 585 494 177 402 976 460 1 150 577 843 748 153 357 280 431 454 369 210 188 449 383 698 735 122 462 521 953 698 895 780 419 18 462 1048 739 339 967 283 835 71 1043 289 662 121 331 18 142 788 646 567 809 388 571 624 693 982 533 446 263 846 643 705 145 74 374 917 77 519 153 902 799 610 1034 313 19 393 169 794 1 212 675 722 49 1030 562 417 380 1031 428 438 38 269 458 310 446 252 209 338 373 383 449 496 971 806 889 441 510 941 155 678 541 461 983 611 1046 1026 966 969 709 951 646 854 284 992 800 519 973 264 939 1020 661 694 658 351 1044 92 106 148 851 712 249 903 37 91 912 351 865 72 616 443 142 957 1032 725 453 367 32 731 877 397 903 587 519 532 817 663 242 1025 120 455 342 1034 177 172 1015 878 124 380 704 233 712 140 988 162 588 865 628 572 474 49 250 325 411 674 508 282 515 370 513 1004 1030 466 556 240 512 548 1000 398 613 905 413 592 798 182 271 177 603 913 677 888 944 543 163 1] \n";

    std::cout << "\n\n";

}

void test_ModEvalBigCharWeber(){
    std::cout << "Entering the test \n";

    NTL::ZZ p = NTL::NextPrime( NTL::power(Integer(2),29) );
    while (p % 4 == 1) {
        p = NTL::NextPrime(p + 1);
    }
    std::cout << "choice of big prime = " << p << "\n";

    NTL::ZZ_p::init(p);

    std::cout << "Testing Modular Evaluation for Big Char with Webers for p = " << p << std::endl;

    // NTL::ZZ l = NTL::ZZ(100003);
    NTL::ZZ l(101); // Level
    // Integer l(11);
    Fp_big_elem j = Fp_big_elem(2);

    //bool weber = false;
    //auto F = ModEvalBigLevel(p, j, l, weber);
    auto t = clock();
    auto F = ModEvalBigCharacteristicWeber(p, j, l);

    std::cout << "ModEvalBigCharWeber took " << (double) (clock()-t)/CLOCKS_PER_SEC << " seconds \n";

    std::cout << "Phi_ell = " << F << "\n\n";

    //std::cout << "Supposed output = " << "[207162232 138521550 305968275 174073185 352015099 149742940 480403223 238119311 14899027 364363891 277800298 315628770 96156811 142525374 344290142 128318772 420678799 471939641 246068059 129077453 352914432 129097656 63254404 532228408 377251769 312751983 145850017 393236377 23112267 523304411 209164381 13190603 169055942 42505548 378862041 466964402 502910090 271518085 171309075 50275568 232602181 502822696 47567917 497107244 478027451 95166378 158341651 391590456 179208684 5617259 25635802 86392412 343977418 255371736 31037424 498824295 165294885 131486149 25667923 132728838 158688818 169933810 396635670 279339869 188468679 420152301 471592472 535512029 302811917 73981049 131630619 280419556 523176516 41804911 201149311 333554320 464857155 197611339 149119561 273529863 129752714 307454458 365386927 486997018 292163959 357979101 275267307 60040085 434689172 93632411 84499215 188868089 401438770 179997411 361268344 397579196 4556521 333455524 376543509 324617194 345311022 107536323 294078293 499049044 115202482 455343940 181976285 97197050 268007110 97435749 513158051 157873558 324024599 199560808 272719965 56826139 440406102 63532039 500519701 334554038 295642232 386702835 286210788 267814135 263998620 212049372 198333256 8964 400987618 152843391 57533987 184804816 217317417 140145335 10147839 361483581 98768261 170757092 304308734 441230068 262548023 122072700 427405641 478882023 264661919 234565130 3339966 294775017 31374358 447895621 170624259 472043832 135703667 87601709 12973065 125263253 453545630 362076793 97228431 295704883 503844703 117115742 501749256 149214560 111474691 9011699 424650242 159231814 335326892 421303812 283225981 112125963 61879022 92450551 398852596 427838987 383132656 304574031 337054649 524243724 203214245 310400932 300684402 345853760 479829667 444935645 282730784 490282124 239111130 300797859 482631412 225068742 372076914 42777929 371698258 435598075 234234854 238424531 316739488 291161245 347813289 319384801 494044850 107961502 229690811 39577104 25589141 272933366 68499476 188598662 199557272 128322509 395091679 345702077 314485697 37870335 402035282 194809441 15652007 124296820 467860264 457225235 266478735 107026518 1]\n";
}

void test_SupersingularEvaluation() {
    //NTL::ZZ p(62207);
    //NTL::ZZ p(1000003);
    Fp_integer p;
    NTL::conv(p, 1093);
    Fp_elem::init(p);
    FpX_elem f;
    SetCoeff(f, 2);
    f[0] = Fp_elem(19);
    FpE_elem::init(f);

    //FpE_elem j = FpE_elem(352946);
    FpE_elem j = FpE_elem(697);
    //ec E(FpE_elem(1010), FpE_elem(11161));
    //FpE_elem j = E.j_invariant();

    //auto roots = SupersingularEvaluation(p, j, NTL::ZZ(101), true);
    auto roots = SupersingularEvaluation(p, j, 11, true);

    FpEX_elem Phi_l(1);
    FpEX_elem X; NTL::SetX(X);
    for (auto j : roots) {
        std::cout << j << std::endl;
        Phi_l *= X - j;
    }

    std::cout << "\n\n\n~~~~~~ recovered poly: ~~~~~~" << std::endl;
    std::cout << Phi_l << std::endl;
}

void test_ModEvalBigLevel(){

    NTL::ZZ p = NTL::NextPrime( NTL::power(Integer(2),29) );
    while (p % 4 == 1) {
        p = NTL::NextPrime(p + 1);
    }

    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    SetCoeff(f, 2);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);

    std::cout << "Testing Modular Evaluation for Big Level for p = " << p << std::endl;

    //NTL::ZZ l = NTL::ZZ(100003);
    long l = 5;
    //NTL::ZZ l = NTL::ZZ(7); // Level
    NTL::ZZ_pE j = NTL::ZZ_pE(5);

    //bool weber = false;
    //auto F = ModEvalBigLevel(p, j, l, weber);

    auto F = ModEvalBigLevel(p, j, l);

    for(int i = 0; i <= deg(F); i++){
        std::cout << "Coefficient of x^" << i << " is: " << NTL::coeff(F,i) << std::endl;
    }
}

void test_SupersingularEvaluationWeber() {
    //NTL::ZZ p(62207);
    Fp_integer p;
    //NTL::conv(p, 18131);
    //NTL::conv(p, 2791);
    //NTL::conv(p, 743);
    //NTL::conv(p, 257);
    NTL::conv(p, 1073741909);
    //NTL::ZZ p(1000003);
    Fp_elem::init(p);
    FpX_elem f;
    SetCoeff(f, 2);
    f[0] = Fp_elem(3);
    FpE_elem::init(f);

    //FpE_elem j = FpE_elem(352946);
    FpE_elem j = FpE_elem(601543879);
    //ec E(FpE_elem(1010), FpE_elem(11161));
    //FpE_elem j = E.j_invariant();

    FpE_elem w = GetWeberDomainFp(j);
    //FpE_elem w = FpE_elem(384);
    //FpE_elem w = FpE_elem(123);
    //FpE_elem w = FpE_elem(2);
    long l = 307;
    auto roots = SupersingularEvaluationWeber(p, w, l);

    FpEX_elem Phi_l(1);
    FpEX_elem Phi_l_classic(1);
    FpEX_elem X; NTL::SetX(X);
    for (auto w : roots) {
        //std::cout << w << std::endl;

        FpE_elem j = NTL::power((NTL::power(w,24)-16),3)/NTL::power(w,24);
        Phi_l *= X - w;
        Phi_l_classic *= X - j;
    }

    std::cout << "\n\n\n~~~~~~ recovered weber poly: ~~~~~~" << std::endl;
    std::cout << Phi_l << std::endl;
    std::cout << "\n~~~~~~ classic poly: ~~~~~~" << std::endl;
    std::cout << Phi_l_classic << std::endl;
}

void test_ModEvalBigLevelWeber(){
    //assert(false);
    NTL::ZZ p = NTL::NextPrime( NTL::power(Integer(2),29) );
    while (p % 4 == 1) {
        p = NTL::NextPrime(p + 1);
    }

    NTL::ZZ_p::init(p);

    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    SetCoeff(f, 2);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);


    // NTL::ZZ l = NTL::ZZ(100003);
    long l = 101; // Level
    FpE_big_elem w = FpE_big_elem(12);

    std::cout << "Testing Modular Evaluation for Big Level with Webers for p = " << p << ", ell = " << l << std::endl;

    auto F = ModEvalBigLevelWeber(p, w, l);

    for(int i = 0; i <= deg(F); i++){
        std::cout << "Coefficient of x^" << i << " is: " << NTL::coeff(F,i) << std::endl;
    }

    // std::cout << "supposed output = " << eval_phi11_weber(j) << "\n";
}

int main()
{
    //test_SupersingularEvaluation();
    //test_SupersingularEvaluationWeber();
    test_ModEvalBigCharWeber();
    test_ModEvalBigLevelWeber();
    test_ModEvalBigLevel();
    return 0;
}
