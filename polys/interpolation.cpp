
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <iostream>
#include <vector>
#include <cmath> 

std::vector<std::vector<NTL::ZZ_pEX>> build_tree(std::vector<NTL::ZZ_pEX> const m, int k, size_t pow2k){
    // k is such that 2^(k-1) < m.size() <= 2^k
    std::vector<std::vector<NTL::ZZ_pEX>> tree(k+1); 
    tree[0] = m;

    size_t level = pow2k;

    // Construct all levels of tree
    for(int i = 1; i <= k; i++){
        level >>= 1;
        for(int j = 0; j < level; j++){
            tree[i].push_back(tree[i-1][2*j]*tree[i-1][2*j+1]);
        }
    }

    return tree;
}


std::vector<NTL::ZZ_pE> going_down_tree(NTL::ZZ_pEX f, std::vector<NTL::ZZ_pE> inputs, std::vector<std::vector<NTL::ZZ_pEX>> tree, int k, size_t pow2k){
    
    size_t level = pow2k;
    level >>= 1;

    // BASE CASE
    std::vector<NTL::ZZ_pEX> r0(level), r1(level), tmp0(level), tmp1(level);
    
    r0[0] = f%tree[k-1][0];
    r1[0] = f%tree[k-1][1];

    // Going down rest of tree (except last step)
    std::size_t pow2ki = 2;
    for(int i = k-2; i > 0; i--){
        for(int j = 0; j < (pow2ki >> 1); j++){
            int ind = 2*j;
            tmp0[ind] = r0[j]%tree[i][ind];
            tmp0[ind+1] = r0[j]%tree[i][ind+1];
            tmp1[ind] = r1[j]%tree[i][pow2ki+ind];
            tmp1[ind+1] = r1[j]%tree[i][pow2ki+ind+1];
        }
        r0 = tmp0;
        r1 = tmp1;
        pow2ki <<= 1;
    }

    //  Last step
    std::vector<NTL::ZZ_pE> r(level << 1);
    for(int j = 0; j < (pow2ki >> 1); j++){
        // These should actually be in NTL::ZZ_p now
        int ind = 2*j;
        r[ind] = NTL::LeadCoeff(r0[j]%tree[0][ind]);
        r[ind+1] = NTL::LeadCoeff(r0[j]%tree[0][ind + 1]);
        r[pow2ki + ind] = NTL::LeadCoeff(r1[j]%tree[0][(pow2ki) + ind]);
        r[pow2ki + ind+1] = NTL::LeadCoeff(r1[j]%tree[0][(pow2ki) + ind + 1]);
    }

    return r;
}


NTL::ZZ_pEX linear_combination(std::vector<NTL::ZZ_pE> const inputs, std::vector<NTL::ZZ_pEX> const c, std::vector<std::vector<NTL::ZZ_pEX>> const tree, int k, size_t const pow2k){

    
    std::vector<NTL::ZZ_pEX> R = c;
    size_t level = pow2k;

    // compute the linear combinations
    std::vector<NTL::ZZ_pEX> tmp(level);
    for(int i = 1; i <= k; i++){
        level >>= 1;
        for(int j = 0; j < level; j++){
            int ind = 2*j;
            tmp[j] = R[ind]*tree[i-1][ind+1]+R[ind+1]*tree[i-1][ind];
        }
        R = tmp;
    }

    return R[0];

}



NTL::ZZ_pEX fast_interpolate(std::vector<NTL::ZZ_pE> const inputs, std::vector<NTL::ZZ_pE> const outputs){

    size_t n = inputs.size();
    int k = ceil(log2(n)); 

    // Pad out inputs and outputs to be size a power of 2 once so we dont have to do it in the rest of the functions
    size_t pow2k = pow(2,k);
    std::vector<NTL::ZZ_pE> inputs2k = inputs;
    std::vector<NTL::ZZ_pE> outputs2k = outputs;
    if(pow2k != n){
        std::size_t diff = pow2k-n;
        std::vector<NTL::ZZ_pE> extra(diff, NTL::ZZ_pE(1));
        inputs2k.insert(inputs2k.end(), extra.begin(), extra.end());
        outputs2k.insert(outputs2k.end(), extra.begin(), extra.end());
    }


    std::vector<NTL::ZZ_pEX> m(pow2k);
    for(int i = 0; i < n; i++){
        NTL::SetCoeff(m[i], 1);
        NTL::SetCoeff(m[i], 0, NTL::ZZ_pE(-inputs[i]));
    }
    for(int i = n; i < pow2k; i++){
        NTL::SetCoeff(m[i], 0, NTL::ZZ_pE(1));
    }

    // Build tree
    std::vector<std::vector<NTL::ZZ_pEX>> tree = build_tree(m, k, pow2k);
    
    // m' is the derivative of prod(m), cheaper than div
    NTL::ZZ_pEX mp;
    NTL::diff(mp, tree[k][0]); 
    
    std::vector<NTL::ZZ_pE> mp_eval = going_down_tree(mp, inputs2k, tree, k, pow2k);
    
    // mp_eval.resize(n);
    std::vector<NTL::ZZ_pEX> c(pow2k);
    // Maybe batch inverting mp_eval may be faster
    for(int i = 0; i < n; i++){
        NTL::ZZ_pE tmp;
        NTL::div(tmp, outputs2k[i], mp_eval[i]);
        NTL::SetCoeff(c[i], 0, tmp);
    }
    for(int i = n; i < pow2k; i++){
        NTL::SetCoeff(c[i], 0, 0);
    }
    
    return linear_combination(inputs2k, c, tree, k, pow2k);

}




NTL::ZZ_pEX fast_interpolate_from_roots(std::vector<NTL::ZZ_pE> const roots){

    size_t n = roots.size();
    int k = ceil(log2(n)); 

    // Pad out inputs and outputs to be size a power of 2 once so we dont have to do it in the rest of the functions
    size_t pow2k = pow(2,k);
    std::vector<NTL::ZZ_pE> inputs = roots;
    if(pow2k != n){
        std::size_t diff = pow2k-n;
        std::vector<NTL::ZZ_pE> extra(diff, NTL::ZZ_pE(1));
        inputs.insert(inputs.end(), extra.begin(), extra.end());
    }

    // std::cout << "inputs:" << std::endl;
    // for(int i = 0; i < inputs.size(); i++){
    //     std::cout << "i = " << inputs[i] << std::endl;
    // }


    std::vector<NTL::ZZ_pEX> m(pow2k);
    for(int i = 0; i < n; i++){
        NTL::SetCoeff(m[i], 1);
        NTL::SetCoeff(m[i], 0, NTL::ZZ_pE(-roots[i]));
    }
    for(int i = n; i < pow2k; i++){
        NTL::SetCoeff(m[i], 0, NTL::ZZ_pE(1));
    }

    // Build tree
    std::vector<std::vector<NTL::ZZ_pEX>> tree = build_tree(m, k, pow2k);
    
    // return the top of the tree which is prod(m)
    return tree[k][0];
}