#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <iostream>
#include <vector>
#include "crt.hpp"
#include "biglevel_hybrid.hpp"
#include "interpolation.hpp"

#include <cassert>
#include <chrono>


void test_fast_interpolate(int Ntests, int Ninputs){

    // This may fail because we construct the polys in an adhoc way, just rerun the tests if so...

    std::cout << "Testing fast_interpolate..." << std::endl;

    NTL::ZZ p;
    int bits = 40;
    NTL::RandomPrime(p, bits);
    std::cout << "Using prime " << p << " which is " << bits << " bits long." << std::endl;
    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    SetCoeff(f, 1);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);

    std::cout << "Running " << Ntests << " tests, and number of inputs is " << Ninputs << std::endl;


    for(int i = 0; i < Ntests; i++){
        NTL::ZZ_pEX F_real;
        NTL::random(F_real, Ninputs-1);
        NTL::vec_ZZ_pE inputs, outputs;
        NTL::random(inputs, Ninputs);
        
        NTL::eval(outputs, F_real, inputs);

        std::vector<NTL::ZZ_pE> inputs1(Ninputs);
        std::vector<NTL::ZZ_pE> outputs1(Ninputs);
        for(int i = 0; i < Ninputs; i++){
            if(outputs[i]*inputs[i] != 0){
                // To try and minimise the bug where we can't invert
                inputs1[i] = inputs[i];
                outputs1[i] = outputs[i];
                }
        }
        NTL::ZZ_pEX F = fast_interpolate(inputs1, outputs1);

        assert(F == F_real);

    }

    std::cout << "Fast interpolation function is working!" << std::endl;

    //TIMING
    
    NTL::ZZ_pEX F_real;
    NTL::random(F_real, Ninputs-1);
    NTL::vec_ZZ_pE inputs, outputs;
    NTL::random(inputs, Ninputs);
    
    NTL::eval(outputs, F_real, inputs);

    std::vector<NTL::ZZ_pE> inputs1(Ninputs);
    std::vector<NTL::ZZ_pE> outputs1(Ninputs);
    for(int i = 0; i < Ninputs; i++){
        inputs1[i] = inputs[i];
        outputs1[i] = outputs[i];
    }
    NTL::ZZ_pEX F;

    auto start = std::chrono::high_resolution_clock::now();   
    F = fast_interpolate(inputs1, outputs1);
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration_fast = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    NTL::ZZ_pEX F_inbuilt;
    auto start1 = std::chrono::high_resolution_clock::now();
    NTL::interpolate(F_inbuilt, inputs, outputs);
    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration_inbuilt = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    
    assert(F_inbuilt == F);

    std::cout << "It takes " << duration_fast.count() << " microseconds compared to " << duration_inbuilt.count() << " microseconds." << std::endl;
    std::cout << "" << std::endl;
    
    
}



void test_fast_interpolate_from_roots(int Ntests, int Ninputs){

    // This may fail because we construct the polys in an adhoc way, just rerun the tests if so...

    std::cout << "Testing fast_interpolate_from_roots..." << std::endl;

    NTL::ZZ p;
    int bits = 40;
    NTL::RandomPrime(p, bits);
    std::cout << "Using prime " << p << " which is " << bits << " bits long." << std::endl;
    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    SetCoeff(f, 1);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);

    std::cout << "Running " << Ntests << " tests, and number of inputs is " << Ninputs << std::endl;


    for(int i = 0; i < Ntests; i++){
        NTL::vec_ZZ_pE inputs;
        NTL::random(inputs, Ninputs);
        NTL::ZZ_pEX F_real = NTL::BuildFromRoots(inputs);
        

        std::vector<NTL::ZZ_pE> inputs1(Ninputs);
        for(int i = 0; i < Ninputs; i++){
            if(inputs[i] != 0){
                // To try and minimise the bug where we can't invert
                inputs1[i] = inputs[i];
                }
        }
        auto F = fast_interpolate_from_roots(inputs1);

        assert(F == F_real);

    }

    std::cout << "Fast interpolation from roots function is working!" << std::endl;

    //TIMING
    NTL::vec_ZZ_pE inputs, outputs;
    NTL::random(inputs, Ninputs);

    NTL::ZZ_pEX F_real = NTL::BuildFromRoots(inputs);
    

    std::vector<NTL::ZZ_pE> inputs1(Ninputs);
    for(int i = 0; i < Ninputs; i++){
        inputs1[i] = inputs[i];
    }
    NTL::ZZ_pEX F;

    auto start = std::chrono::high_resolution_clock::now();   
    F = fast_interpolate_from_roots(inputs1);
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration_fast = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    NTL::ZZ_pEX F_inbuilt;
    auto start1 = std::chrono::high_resolution_clock::now();
    F_inbuilt = NTL::BuildFromRoots(inputs);
    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration_inbuilt = std::chrono::duration_cast<std::chrono::microseconds>(stop1 - start1);
    
    assert(F_inbuilt == F);

    std::cout << "It takes " << duration_fast.count() << " microseconds compared to " << duration_inbuilt.count() << " microseconds." << std::endl;
    std::cout << "" << std::endl;
    
    
}

// void test_special_ss_eval(){
//     std::cout << "Testing special_ss_eval..." << std::endl;

//     NTL::ZZ p = NTL::ZZ(1048583);
//     NTL::ZZ q = NTL::ZZ(12889);
//     NTL::ZZ j = NTL::ZZ(1728);
//     NTL::ZZ ell = NTL::ZZ(2);

//     // Want to check if the modulus switching is working correctly
//     NTL::ZZ_p::init(p);
//     NTL::ZZ_pX f;
//     SetCoeff(f, 2);
//     f[0] = NTL::ZZ_p(1);
//     NTL::ZZ_pE::init(f);

//     NTL::ZZX Fq;

//     int ell_int = NTL::conv<int>(ell);
//     std::vector<NTL::ZZ> xs(ell_int+1);
//     for(int i = 0; i <= ell_int; i++){
//         NTL::power(xs[i], j, i+1);
//     }    
    
//     NTL::ZZX tmp1, tmp2, tmp3;
//     SetCoeff(tmp1, 1);
//     SetCoeff(tmp2, 1);
//     SetCoeff(tmp3, 1);
//     tmp1[0] = 202735;
//     tmp1[1] = 189820;
//     tmp2[0] = 5805;
//     tmp2[1] = 259948;
//     tmp3[0] = 45523;
//     tmp3[1] = 241713;
//     std::vector<NTL::ZZX> js_polys = {tmp1, tmp2, tmp3};
    
//     NTL::ZZX tmp11, tmp12, tmp13;
//     SetCoeff(tmp11, 1);
//     SetCoeff(tmp12, 1);
//     SetCoeff(tmp13, 1);
//     tmp11[0] = 5805;
//     tmp11[1] = 259948;
//     tmp12[0] = 203959;
//     tmp12[1] = 99584;
//     tmp13[0] = 107009;
//     tmp13[1] = 43540;
//     NTL::ZZX tmp21, tmp22, tmp23;
//     SetCoeff(tmp21, 1);
//     SetCoeff(tmp22, 1);
//     SetCoeff(tmp23, 1);
//     tmp21[0] = 45523;
//     tmp21[1] = 241713;
//     tmp22[0] = 50040;
//     tmp22[1] = 211646;
//     tmp23[0] = 202735;
//     tmp23[1] = 189820;
//     NTL::ZZX tmp31, tmp32, tmp33;
//     SetCoeff(tmp31, 1);
//     SetCoeff(tmp32, 1);
//     SetCoeff(tmp33, 1);
//     tmp31[0] = 5805;
//     tmp31[1] = 259948;
//     tmp32[0] = 88943;
//     tmp32[1] = 169810;
//     tmp33[0] = 155122;
//     tmp33[1] = 58473;

//     std::vector<std::vector<NTL::ZZX>> iso_js_polys = {{tmp11, tmp12, tmp13}, {tmp21, tmp22, tmp23}, {tmp31, tmp32, tmp33}};

//     // js = [202735*z2 + 189820, 5805*z2 + 259948, 45523*z2 + 241713]
//     // iso_js = [[5805*z2 + 259948, 203959*z2 + 99584, 107009*z2 + 43540],
//     // [45523*z2 + 241713, 50040*z2 + 211646, 202735*z2 + 189820],
//     // [5805*z2 + 259948, 88943*z2 + 169810, 155122*z2 + 58473]]

//     std::cout << "Initialised everything..." << std::endl;

//     special_ss_eval(Fq, q, xs, ell);

//     // X^3 + 209721*X^2 + 98550*X + 179157
//     NTL::ZZX Fq_real;
//     SetCoeff(Fq_real, 3);
//     SetCoeff(Fq_real, 2, 209721);
//     SetCoeff(Fq_real, 1, 98550);
//     SetCoeff(Fq_real, 0, 179157);

//     assert(Fq == Fq_real);

//     std::cout << "special_ss_eval function is working!" << std::endl;
//     std::cout << "" << std::endl;

// }


void test_biglevel(){

    NTL::ZZ p = NTL::ZZ(8647);
    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    SetCoeff(f, 2);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);

    NTL::ZZ_p j = NTL::ZZ_p(1728);
    NTL::ZZ ell = NTL::ZZ(2);

    NTL::ZZ_pX F;
    biglevel(F, p, j, ell);
    for(int i = 0; i <= deg(F); i++){
        std::cout << "Coefficient of x^" << i << " is: " << NTL::coeff(F,i) << std::endl;
    }
}

int main(){

    int Ntests, Ninputs;
    // std::cout << "Number of tests: " << std::endl;
    // std::cin >> Ntests;
    // std::cout << "Number of inputs: " << std::endl;
    // std::cin >> Ninputs;

    Ntests = 10;
    Ninputs = pow(2,8)-13;
    test_fast_interpolate_from_roots(Ntests, Ninputs);
    test_fast_interpolate(Ntests, Ninputs);
    

    // test_special_ss_eval();

    // test_biglevel();

    return 0;
}
