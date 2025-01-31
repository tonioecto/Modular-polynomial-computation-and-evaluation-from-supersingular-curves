
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cmath> 
#include "bivariates.hpp" 


//// Interpolation for Bivariates

std::vector<std::vector<ZZ_pEXY>> BuildTreeBivar(std::vector<ZZ_pEXY> const m, unsigned k, size_t pow2k){
    // k is such that 2^(k-1) < m.size() <= 2^k
    std::vector<std::vector<ZZ_pEXY>> tree(k+1); 
    tree[0] = m;

    size_t level = pow2k;

    // Construct all levels of tree
    for(unsigned i = 1; i <= k; i++){
        level >>= 1;
        for(unsigned j = 0; j < level; j++){
            tree[i].push_back(tree[i-1][2*j]*tree[i-1][2*j+1]);
        }
    }

    return tree;
}



// ZZ_pEXY FastInterpolateBivar(std::vector<NTL::ZZ_pEX> const inputs, std::vector<NTL::ZZ_pEX> const outputs){

//     size_t n = inputs.size();
//     int k = ceil(log2(n)); 

//     // Pad out inputs and outputs to be size a power of 2 once so we dont have to do it in the rest of the functions
//     size_t pow2k = pow(2,k);
//     std::vector<NTL::ZZ_pEX> inputs2k = inputs;
//     std::vector<NTL::ZZ_pEX> outputs2k = outputs;
//     if(pow2k != n){
//         std::size_t diff = pow2k-n;
//         std::vector<NTL::ZZ_pEX> extra(diff, NTL::ZZ_pEX(1));
//         inputs2k.insert(inputs2k.end(), extra.begin(), extra.end());
//         outputs2k.insert(outputs2k.end(), extra.begin(), extra.end());
//     }


//     std::vector<ZZ_pEXY> m(pow2k);
//     for(unsigned i = 0; i < pow2k; i++){
//         NTL::Vec<NTL::ZZ_pEX> cs;
//         cs.SetLength(2);
//         cs[0] = NTL::ZZ_pEX(1);
//         cs[1] = NTL::ZZ_pEX(-inputs2k[i]);
//     }

//     // Build tree
//     std::vector<std::vector<ZZ_pEXY>> tree = BuildTreeBivar(m, k, pow2k);
    
//     // // m' is the derivative of prod(m), cheaper than div
//     // NTL::ZZ_pEX mp;
//     // NTL::diff(mp, tree[k][0]); 
    
//     // std::vector<NTL::ZZ_pE> mp_eval = GoingDownTree(mp, inputs2k, tree, k, pow2k);
    
//     // // mp_eval.resize(n);
//     // std::vector<NTL::ZZ_pEX> c(pow2k);
//     // // TODO: batch inverting mp_eval may be faster
//     // for(unsigned i = 0; i < n; i++){
//     //     NTL::ZZ_pE tmp;
//     //     NTL::div(tmp, outputs2k[i], mp_eval[i]);
//     //     NTL::SetCoeff(c[i], 0, tmp);
//     // }
//     // for(unsigned i = n; i < pow2k; i++){
//     //     NTL::SetCoeff(c[i], 0, 0);
//     // }
    
//     // return m[0];

// }