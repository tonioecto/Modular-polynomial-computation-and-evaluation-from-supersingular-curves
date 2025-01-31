#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include "interpolation.hpp"

#include <cassert>


void test_buildtree(){
    NTL::ZZ_p::init(NTL::ZZ(1019));

    std::vector<std::vector<NTL::ZZ_pX>> M;

    NTL::ZZ_pX m1, m2, m3;
    NTL::SetCoeff(m1, 0, NTL::ZZ_p(1));
    NTL::SetCoeff(m1, 1);
    NTL::SetCoeff(m2, 0, NTL::ZZ_p(2));
    NTL::SetCoeff(m2, 1);
    NTL::SetCoeff(m3, 0, NTL::ZZ_p(3));
    NTL::SetCoeff(m3, 1);

    std::vector<NTL::ZZ_pX> m = {m1, m2, m3};

    std::cout << "Testing buildtree..." << std::endl;
    
    M = build_tree(m);
    for(int j = 0; j <= 2; j++){
        for(int i = 0; i < M[j].size(); i++){
        std::cout << "M[" << j << "][" << i << "] = " << M[j][i] << std::endl;
        }
    }

}

void test_going_down_tree(){

    NTL::ZZ_p::init(NTL::ZZ(1019));
    NTL::ZZ_pX f;
    NTL::SetCoeff(f, 0);
    NTL::SetCoeff(f, 2);

    std::vector<NTL::ZZ_p> inputs = {NTL::ZZ_p(1), NTL::ZZ_p(2), NTL::ZZ_p(3), NTL::ZZ_p(4), NTL::ZZ_p(5), NTL::ZZ_p(6), NTL::ZZ_p(7), NTL::ZZ_p(8)};

    std::vector<NTL::ZZ_pX> m(inputs.size());
    for(int i = 0; i < inputs.size(); i++){
        NTL::ZZ_pX tmp;
        NTL::SetCoeff(tmp, 1);
        NTL::SetCoeff(tmp, 0, NTL::ZZ_p(-inputs[i]));
        m[i] = tmp;
    }

    std::vector<std::vector<NTL::ZZ_pX>> tree = build_tree(m);    

    std::vector<NTL::ZZ_p> evals = going_down_tree(f, inputs, tree);
    // for(int i = 0; i < evals.size(); i++){
    //     std::cout << "evals[" << i << "] = " << evals[i] << std::endl;
    // }

    std::vector<NTL::ZZ_p> evals_real = {NTL::ZZ_p(2), NTL::ZZ_p(5), NTL::ZZ_p(10), NTL::ZZ_p(17), NTL::ZZ_p(26), NTL::ZZ_p(37), NTL::ZZ_p(50), NTL::ZZ_p(65)};

    
    if(evals_real == evals){
        std::cout << "Function is all good!" << std::endl;
    }
    else{
        std::cout << "There's a problem." << std::endl;
    }

}


void test_linear_combo(){

    NTL::ZZ_p::init(NTL::ZZ(1019));
    std::vector<NTL::ZZ_p> inputs = {NTL::ZZ_p(1), NTL::ZZ_p(2), NTL::ZZ_p(3), NTL::ZZ_p(4), NTL::ZZ_p(5), NTL::ZZ_p(6), NTL::ZZ_p(7)};
    std::vector<NTL::ZZ_p> c = {NTL::ZZ_p(2), NTL::ZZ_p(3), NTL::ZZ_p(6), NTL::ZZ_p(8), NTL::ZZ_p(11), NTL::ZZ_p(14), NTL::ZZ_p(23)};

    std::vector<NTL::ZZ_pX> m(inputs.size());
    for(int i = 0; i < inputs.size(); i++){
        NTL::ZZ_pX tmp;
        NTL::SetCoeff(tmp, 1);
        NTL::SetCoeff(tmp, 0, NTL::ZZ_p(-inputs[i]));
        m[i] = tmp;
    }

    std::vector<std::vector<NTL::ZZ_pX>> tree = build_tree(m); 

    std::cout << "Testing Linear combination function...";

    NTL::ZZ_pX r = linear_combination(inputs, c, tree);

    for(int i = 0; i <= NTL::deg(r); i++){
        std::cout << "Coefficient of X^" << i << ": " << coeff(r, i) << std::endl;
    }

}

void test_fast_interpolate(int Ntests, int Ninputs){


    NTL::ZZ_p::init(NTL::ZZ(12889));

    

    std::cout << "Running " << Ntests << " tests, and number of inputs is " << Ninputs << std::endl;

    for(int i = 0; i < Ntests; i++){
        NTL::ZZ_pX F_real;
        NTL::random(F_real, Ninputs-1);
        NTL::vec_ZZ_p inputs, outputs;
        NTL::random(inputs, Ninputs);
        
        NTL::eval(outputs, F_real, inputs);

        std::vector<NTL::ZZ_p> inputs1(Ninputs);
        std::vector<NTL::ZZ_p> outputs1(Ninputs);
        for(int i = 0; i < Ninputs; i++){
            inputs1[i] = inputs[i];
            outputs1[i] = inputs[i];
        }

        NTL::ZZ_pX F = fast_interpolate(inputs1, outputs1);

        assert(F == F_real);

    }

    std::cout << "Function working!" << std::endl;
    
    
}

int main(){

    // test_buildtree();
    // test_going_down_tree();
    // test_linear_combo();
    int Ntests = 5;
    int Ninputs = 12;
    test_fast_interpolate(Ntests, Ninputs);
    return 0;
}
