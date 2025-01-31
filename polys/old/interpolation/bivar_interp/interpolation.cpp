#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <NTL/vec_ZZ.h>
#include <iostream>
#include <vector>
#include "bivariates.hpp"





void interpolate(NTL::ZZ_pX &f, const std::vector<NTL::ZZ_p> outputs, const std::vector<NTL::ZZ_p> inputs){
    

}


void reconstruct(bivar_poly &f, const std::vector<NTL::ZZ_p> eval_js, const std::vector<NTL::vec_ZZ_p> iso_js, int Noutputs){
    // iso_js[i] are the j-invariants l-isogeneous to eval_js[i]
    // i.e., they are the roots of \Phi_{l}(X, eval_js[i])


    // Construct the output \Phi_{l}(X, eval_js[i]) from the roots iso_js[i] 
    std::vector<NTL::ZZ_pX> outputs(Noutputs);
    for(int i = 0; i < Noutputs; i++){
        NTL::BuildFromRoots(outputs[i], iso_js[i]);
    }


}