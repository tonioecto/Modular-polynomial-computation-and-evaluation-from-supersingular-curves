///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing interpolation of univariate polynomials over finite fields
////   We follow §10 from:
////            Joachim von zur Gathen and Jürgen Gerhard. Modern Computer Algebra. 
////            Cambridge University Press, 1999.    
////
////   Note: for polynomials of large-ish degree, this outperforms NTL in-built.
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <cmath> 
#include "interpolation.hpp"


std::vector<std::vector<FpEX_elem>> BuildTree(std::vector<FpEX_elem> const &m, unsigned k, size_t pow2k){
    // k is such that 2^(k-1) < m.size() <= 2^k
    std::vector<std::vector<FpEX_elem>> tree(k+1); 
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


std::vector<FpE_elem> GoingDownTree(FpEX_elem f, std::vector<FpE_elem> inputs, std::vector<std::vector<FpEX_elem>> tree, int k, size_t pow2k)
{

    size_t level = pow2k;
    level >>= 1;

    // BASE CASE
    std::vector<FpEX_elem> r0(level), r1(level), tmp0(level), tmp1(level);
    
    r0[0] = f%tree[k-1][0];
    r1[0] = f%tree[k-1][1];

    // Going down rest of tree (except last step)
    std::size_t pow2ki = 2;
    for(unsigned i = k-2; i > 0; i--){
        for(unsigned j = 0; j < (pow2ki >> 1); j++){
            unsigned ind = 2*j;
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
    std::vector<FpE_elem> r(level << 1);
    for(unsigned j = 0; j < (pow2ki >> 1); j++){
        // These should actually be in NTL::ZZ_p now
        unsigned ind = 2*j;
        r[ind] = NTL::LeadCoeff(r0[j]%tree[0][ind]);
        r[ind+1] = NTL::LeadCoeff(r0[j]%tree[0][ind + 1]);
        r[pow2ki + ind] = NTL::LeadCoeff(r1[j]%tree[0][(pow2ki) + ind]);
        r[pow2ki + ind+1] = NTL::LeadCoeff(r1[j]%tree[0][(pow2ki) + ind + 1]);
    }

    return r;
}


FpEX_elem LinearCombo(std::vector<FpEX_elem> const c, std::vector<std::vector<FpEX_elem>> const tree, int k, size_t const pow2k)
{

    std::vector<FpEX_elem> R = c;
    size_t level = pow2k;

    // compute the linear combinations
    std::vector<FpEX_elem> tmp(level);
    for(int i = 1; i <= k; i++){
        level >>= 1;
        for(unsigned j = 0; j < level; j++){
            unsigned ind = 2*j;
            tmp[j] = R[ind]*tree[i-1][ind+1]+R[ind+1]*tree[i-1][ind];
        }
        R = tmp;
    }

    return R[0];

}


FpEX_elem FastInterpolate(std::vector<FpE_elem> const inputs, std::vector<FpE_elem> const outputs){

    size_t n = inputs.size();
    int k = ceil(log2(n)); 

    // Pad out inputs and outputs to be size a power of 2 once so we don't have to do it in the rest of the functions
    size_t pow2k = pow(2,k);
    std::vector<FpE_elem> inputs2k = inputs;
    std::vector<FpE_elem> outputs2k = outputs;
    if(pow2k != n){
        std::size_t diff = pow2k-n;
        std::vector<FpE_elem> extra(diff, FpE_elem(1));
        inputs2k.insert(inputs2k.end(), extra.begin(), extra.end());
        outputs2k.insert(outputs2k.end(), extra.begin(), extra.end());
    }


    std::vector<FpEX_elem> m(pow2k);
    for(unsigned i = 0; i < n; i++){
        NTL::SetCoeff(m[i], 1);
        NTL::SetCoeff(m[i], 0, FpE_elem(-inputs[i]));
    }
    for(unsigned i = n; i < pow2k; i++){
        NTL::SetCoeff(m[i], 0, FpE_elem(1));
    }

    // Build tree
    std::vector<std::vector<FpEX_elem>> tree = BuildTree(m, k, pow2k);
    
    // m' is the derivative of prod(m), cheaper than div
    FpEX_elem mp;
    NTL::diff(mp, tree[k][0]); 
    
    std::vector<FpE_elem> mp_eval = GoingDownTree(mp, inputs2k, tree, k, pow2k);
    
    std::vector<FpEX_elem> c(pow2k);
    // TODO: batch inverting mp_eval may be faster
    for(unsigned i = 0; i < n; i++){
        FpE_elem tmp;
        NTL::div(tmp, outputs2k[i], mp_eval[i]);
        NTL::SetCoeff(c[i], 0, tmp);
    }
    for(unsigned i = n; i < pow2k; i++){
        NTL::SetCoeff(c[i], 0, 0);
    }
    
    return LinearCombo(c, tree, k, pow2k);

}


FpEX_elem FastInterpolateFromRoots(std::vector<FpE_elem> const roots){

    size_t n = roots.size();
    int k = ceil(log2(n)); 

    // Pad out inputs and outputs to be size a power of 2 once so we dont have to do it in the rest of the functions
    size_t pow2k = pow(2,k);
    std::vector<FpE_elem> inputs = roots;
    if(pow2k != n){
        std::size_t diff = pow2k-n;
        std::vector<FpE_elem> extra(diff, FpE_elem(1));
        inputs.insert(inputs.end(), extra.begin(), extra.end());
    }

    std::vector<FpEX_elem> m(pow2k);
    for(unsigned i = 0; i < n; i++){
        NTL::SetCoeff(m[i], 1);
        NTL::SetCoeff(m[i], 0, FpE_elem(-inputs[i]));
    }
    for(unsigned i = n; i < pow2k; i++){
        NTL::SetCoeff(m[i], 0, FpE_elem(1));
    }

    // Build tree
    std::vector<std::vector<FpEX_elem>> tree = BuildTree(m, k, pow2k);

    // return the top of the tree which is prod(m)
    return tree[k][0];
}
