#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_p.h>
#include <iostream>
#include <vector>
#include <cmath> 
#include <chrono>

// TODO: use arrays if the vector has fixed length, otherwise use std::vector with push_back or resize

std::vector<std::vector<NTL::ZZ_pX>> build_tree(std::vector<NTL::ZZ_pX> const m){
    // k is such that 2^(k-1) < m.size() <= 2^k
    int k = ceil(log2(m.size())); 

    std::vector<std::vector<NTL::ZZ_pX>> tree(k+1); 
    // Should we use an array of pointers?


    std::size_t pow2k = pow(2,k);

    // Setting up M
    // // adjusting if m.size() is not a power of 2 (quite messy rn)
    if(m.size() < pow2k){
        std::vector<NTL::ZZ_pX> tmp(pow2k);
        for(int i = 0; i < m.size(); i++){
            tmp[i] = m[i];
        }
        for(int i = m.size(); i < pow2k; i++){
            NTL::SetCoeff(tmp[i],0);
        }
        tree[0] = tmp;
    }
    else{
        tree[0] = m;
    }

    // Construct all levels of tree
    for(int i = 1; i <= k; i++){
        pow2k = pow2k >> 1;
        tree[i].resize(pow2k);
        for(int j = 0; j < pow2k; j++){
            tree[i][j] = tree[i-1][2*j]*tree[i-1][2*j+1];
        }
    }

    return tree;
}


std::vector<NTL::ZZ_p> going_down_tree(NTL::ZZ_pX f, std::vector<NTL::ZZ_p> inputs, std::vector<std::vector<NTL::ZZ_pX>> tree){
    
    int k = ceil(log2(inputs.size())); 
    std::size_t pow2k = pow(2,k);

    // BASE CASE
    std::vector<NTL::ZZ_pX> r0(1), r1(1);
    rem(r0[0], f, tree[k-1][0]); 
    rem(r1[0], f, tree[k-1][1]); 

    // Going down rest of tree (except last step)
    for(int i = k-2; i > 0; i--){
        std::size_t pow2ki = (pow2k >> i);
        std::vector<NTL::ZZ_pX> tmp0(pow2ki >> 1);
        std::vector<NTL::ZZ_pX> tmp1(pow2ki >> 1);
        for(int j = 0; j < (pow2ki >> 2); j++){
            rem(tmp0[2*j], r0[j], tree[i][2*j]);
            rem(tmp0[2*j+1], r0[j], tree[i][2*j + 1]);
            rem(tmp1[2*j], r1[j], tree[i][(pow2ki >> 1) + 2*j]);
            rem(tmp1[2*j+1], r1[j], tree[i][(pow2ki >> 1) + 2*j + 1]);
        }
        r0.resize(pow2ki >> 1);
        r0 = tmp0;
        r1.resize(pow2ki >> 1);
        r1 = tmp1;
    }

    //  Last step
    std::vector<NTL::ZZ_pX> tmp0(pow2k >> 1);
    std::vector<NTL::ZZ_pX> tmp1(pow2k >> 1);
    for(int j = 0; j < pow2k >> 2; j++){
        // These should actually be in NTL::ZZ_p now
        rem(tmp0[2*j], r0[j], tree[0][2*j]);
        rem(tmp0[2*j+1], r0[j], tree[0][2*j + 1]);
        rem(tmp1[2*j], r1[j], tree[0][(pow2k >> 1) + 2*j]);
        rem(tmp1[2*j+1], r1[j], tree[0][(pow2k >> 1) + 2*j + 1]);
    }
    std::vector<NTL::ZZ_p> r(pow2k);
    for(int i = 0; i < (pow2k >> 1); i++){
        r[i] = NTL::LeadCoeff(tmp0[i]);
        r[(pow2k >> 1) + i] = NTL::LeadCoeff(tmp1[i]);
    }

    return r;
}


NTL::ZZ_pX linear_combination(std::vector<NTL::ZZ_p> const inputs, std::vector<NTL::ZZ_p> const c, std::vector<std::vector<NTL::ZZ_pX>> const tree){
    
    int k = ceil(log2(inputs.size())); 
    std::size_t pow2k = pow(2,k);
    std::vector<NTL::ZZ_pX> R(pow2k);


    // adjusting if number of inputs is not a power of 2
    std::vector<NTL::ZZ_pX> tmp(pow2k);
    for(int i = 0; i < c.size(); i++){
        NTL::SetCoeff(tmp[i],0, c[i]);
    }
    if(c.size() < pow2k){
        for(int i = c.size(); i < pow2k; i++){
            NTL::SetCoeff(tmp[i],0);
        }
    }
    R = tmp;

    // compute the linear combinations
    for(int i = 1; i <= k; i++){
        pow2k = pow2k >> 1;
        std::vector<NTL::ZZ_pX> tmp(pow2k);
        for(int j = 0; j < pow2k; j++){
            tmp[j] = R[2*j]*tree[i-1][2*j+1] + R[2*j+1]*tree[i-1][2*j];
        }
        R.resize(pow2k);
        R = tmp;
    }
    return R[0];
}


NTL::ZZ_pX fast_interpolate(std::vector<NTL::ZZ_p> const inputs, std::vector<NTL::ZZ_p> const outputs){

    std::size_t n = inputs.size();
    int k = ceil(log2(n)); 

    std::vector<NTL::ZZ_pX> m(n);
    for(int i = 0; i < n; i++){
        NTL::SetCoeff(m[i], 1);
        NTL::SetCoeff(m[i], 0, NTL::ZZ_p(-inputs[i]));
    }

    // Build tree
    std::vector<std::vector<NTL::ZZ_pX>> tree = build_tree(m);

    // Construct the mi' = prod(m)/(x-ui) - (x-u1)*...*(x-u(i-1))*(x-u(i+1))*..(x-un)
    
    // NTL::SetCoeff(mp, 0, 0);
    // for(int j = 0; j < n; j++){
    //     NTL::ZZ_pX tmp;
    //     NTL::div(tmp, tree[k][0], m[j]);
    //     mp += tmp;
    // }

    NTL::ZZ_pX mp = NTL::diff(tree[k][0]); // m' is the derivative of prod(m), should be cheaper than div
    
    std::vector<NTL::ZZ_p> mp_eval = going_down_tree(mp, inputs, tree);

    std::vector<NTL::ZZ_p> c(n);
    for(int i = 0; i < n; i++){
        NTL::ZZ_p tmp;
        NTL::inv(tmp, mp_eval[i]);
        c[i] = outputs[i]*tmp;
    }

    NTL::ZZ_pX F = linear_combination(inputs, c, tree);

    return F;
}




