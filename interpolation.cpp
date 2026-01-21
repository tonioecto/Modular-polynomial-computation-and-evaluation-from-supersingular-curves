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
#include <cassert>
#include "interpolation.hpp"
#include "utils.hpp"



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


std::vector<FpE_elem> GoingDownTree(FpEX_elem f, std::vector<std::vector<FpEX_elem>> tree, int k, size_t pow2k)
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

    std::vector<FpE_elem> mp_eval = GoingDownTree(mp, tree, k, pow2k);

    std::vector<FpEX_elem> c(pow2k);
    // TODO: batch inverting mp_eval may be faster
    for(unsigned i = 0; i < n; i++){
        FpE_elem tmp;
        // NTL::div(tmp, outputs2k[i], mp_eval[i]);
        fast_inv(tmp, mp_eval[i]);
        NTL::mul(tmp, tmp, outputs2k[i]);
        NTL::SetCoeff(c[i], 0, tmp);
    }
    for(unsigned i = n; i < pow2k; i++){
        NTL::SetCoeff(c[i], 0, 0);
    }

    return LinearCombo(c, tree, k, pow2k);

}





FpEX_elem FastInterpolateFromRootsRec(std::vector<FpE_elem> const &roots, int level, int index, int max_depth, int len) {

    if (level == max_depth) {
        FpEX_elem m;
        if (index < len) {
            m.SetLength(2);
            NTL::SetCoeff(m, 1, 1);
            NTL::SetCoeff(m, 0, -roots[index]);
            // std::cout << "level = " << level << " index = " << index << " " << m << "\n";
            return m;
        }
        else {
            m.SetLength(1);
            NTL::SetCoeff(m, 0, 1);
            return m;
        }
    }
    else {
        auto res = FastInterpolateFromRootsRec(roots, level + 1, index << 1, max_depth, len) * FastInterpolateFromRootsRec(roots, level + 1, 1 + (index << 1), max_depth, len);
        // std::cout << "level = " << level << " index = " << index << " " << res << "\n";
        return res;
        // FpEX_elem temp;
        // temp.SetLength(1 << (max_depth - level));
        
    }

}

void FastInterpolateFromRoots(FpEX_elem &poly, std::vector<FpE_elem> const &roots) {

    size_t n = roots.size();
    size_t k = ceil(log2(n));

    poly = FastInterpolateFromRootsRec(roots, 0, 0, k, n);

}

std::pair<FpX_elem,FpX_elem> custom_fp2X_mul(std::pair<FpX_elem,FpX_elem> &p1, std::pair<FpX_elem,FpX_elem> &p2) {
    
        FpX_elem temp;
        NTL::add(temp, p1.first, p1.second);
        NTL::mul(temp, temp, (p2.first + p2.second));
        NTL::mul(p2.first, p1.first, p2.first);
        NTL::mul(p2.second, p1.second, p2.second); 

        NTL::sub(p1.first, p2.first, p2.second);
        NTL::add(p1.second, p2.first, p2.second);
        NTL::sub(p1.second, temp, p1.second);
        return p1;

}

void explicit_FFT_custom_fp2X_mul(std::pair<FpX_elem, FpX_elem> &res, std::pair<FpX_elem,FpX_elem>& p1, std::pair<FpX_elem,FpX_elem>& p2) {
        // TODO : set this constant properly

        

        if (deg(p1.first) > 45) {
            long d = (deg(p1.first) + deg(p2.first));
            long k = NTL::NextPowerOfTwo(d + 1);
            // std::cout << deg(p1.first) << " " << deg(p2.first) << " " << d << " " << k << "\n";

            // assert(deg(p1.first) == deg(p1.second) + 1);
            
            // long d = (deg(p1.first) << 1);
            

            // clock_t t = tic();
            NTL::fftRep R1(NTL::INIT_SIZE, k), R2(NTL::INIT_SIZE, k), R3(NTL::INIT_SIZE, k), R4(NTL::INIT_SIZE, k), Rt1(NTL::INIT_SIZE, k), Rt2(NTL::INIT_SIZE, k);
            TofftRep_trunc(R1, p1.first, k, d+1);
            TofftRep_trunc(R2, p1.second, k, d+1);
            TofftRep_trunc(R3, p2.first, k, d+1);
            TofftRep_trunc(R4, p2.second, k, d+1);

            // std::cout << "conv " << tic() - t << "\n";
            // t = tic();

            mul(Rt1, R1, R3);
            mul(Rt2, R2, R4);
            add(R1, R1, R2);
            add(R3, R3, R4);
            mul(R1, R1, R3);
            add(R2, Rt1, Rt2);
            sub(R3, R1, R2);

            sub(R2, Rt1, Rt2);

            // std::cout << "ops  " << tic() -t << "\n";

            // t = tic();
            FromfftRep(res.second, R3, 0, d - 1);
            FromfftRep(res.first, R2, 0, d);

            // std::cout << "back " << tic() -t << "\n";

            

            // Fp temp = Fp(1);
            // for (int i = 0; i <= deg(res.second); i++ ) {
            //     temp = temp * res.first[i];
            //     temp = temp * res.second[i];
            //     temp = temp * res.first[i];
            // }
            

            // return p1;
        }
        else {
            FpX_elem temp;
            NTL::add(temp, p1.first, p1.second);
            NTL::mul(temp, temp, (p2.first + p2.second));
            NTL::mul(p2.first, p1.first, p2.first);
            NTL::mul(p2.second, p1.second, p2.second); 

            NTL::sub(res.first, p2.first, p2.second);
            NTL::add(res.second, p2.first, p2.second);
            NTL::sub(res.second, temp, res.second);
            // return p1;
        }
}



// TODO we could be something like 10% faster by changing a bit the way the roots are processed. Indeed, when the degree is exactly a power of two,
// FFT is a bit inefficient as it needs 2^k + 1 roots, so it would be more efficient to use polynomials are degree smaller than 2^k. When roots.size() is not too close to a power of 2, we could quite easily do this change by working on polynomials of degree 63 (instead of 64) (64 is the first power of 2 above the FFT threshold). 
// From some experiments it appears that this would allow us to gain 10% time.  
void FastInterpolateFromRootsKaratsubaRec(std::pair<FpX_elem, FpX_elem>& m, std::vector<ffp2> const &roots, int level, int index, int max_depth, int len) {

    if (level == max_depth - 1) {
        int new_index = index << 1;
        if (new_index < (len - 1)) {
            // two roots
            m.first.SetLength(3);
            m.second.SetLength(2);
            // NTL::SetCoeff(m.first, 2);
            m.first[2] = Fp(1);
            
            // NTL::SetCoeff(m.first, 1, - (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0)));
            m.first[1] = - (roots[new_index].first + roots[new_index + 1].first);  
            // assert( (- (roots[new_index]._zz_pE__rep[0] + roots[new_index + 1]._zz_pE__rep[0])) ==  (- (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0))));
            // NTL::SetCoeff(m.second, 1, - (NTL::coeff(rep(roots[new_index]), 1) + NTL::coeff(rep(roots[new_index + 1]), 1)));
            m.second[1] = - (roots[new_index].second + roots[new_index + 1].second);  

            // karatsuba mul
            // Fp temp = (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index]), 1)) * (NTL::coeff(rep(roots[new_index + 1]), 0) + NTL::coeff(rep(roots[new_index + 1]), 1));
            Fp temp = (roots[new_index].first + roots[new_index].second) * (roots[new_index + 1].first + roots[new_index + 1].second);

            // Fp t1 =  (NTL::coeff(rep(roots[new_index]), 0) * NTL::coeff(rep(roots[new_index + 1]), 0));
            // Fp t2 =  (NTL::coeff(rep(roots[new_index]), 1) * NTL::coeff(rep(roots[new_index + 1]), 1));
            Fp t1 = roots[new_index].first * roots[new_index + 1].first;
            Fp t2 = roots[new_index].second * roots[new_index + 1].second;
            // NTL::SetCoeff(m.first, 0, t1 - t2);
            // NTL::SetCoeff(m.second, 0, temp - (t1 + t2));
            m.first[0] = t1 - t2;
            m.second[0] = temp - (t1 + t2);
            
        }
        else if (new_index == len - 1) {
            // only one root
                m.first.SetLength(2); // real part
                m.second.SetLength(1); // imaginary part

                m.first[1] = Fp(1);
                m.first[0] = - roots[new_index].first;
                m.second[0] = - roots[new_index].second;
                // NTL::SetCoeff(m.first, 1);
                // NTL::SetCoeff(m.first, 0, - NTL::coeff(rep(roots[new_index]),0));
                // NTL::SetCoeff(m.second, 0, - NTL::coeff(rep(roots[new_index]),1));
        }
        else {

            m.first.SetLength(1);
            m.second.SetLength(1);
            // NTL::SetCoeff(m.first, 0, 1);
            // NTL::set(m.first);
            m.first[0] = Fp(1);
            // NTL::SetCoeff(m.second, 0, 0);
            m.second = m.second.zero();
            // return m;
        }
    }
    // if (level == max_depth) {
    //     if (index < len) {
    //         // only one root
    //         m.first.SetLength(2); // real part
    //         m.second.SetLength(1); // imaginary part

    //         NTL::SetCoeff(m.first, 1, 1);
    //         NTL::SetCoeff(m.first, 0, - NTL::coeff(rep(roots[ index ]), 0));
    //         NTL::SetCoeff(m.second, 0, - NTL::coeff(rep(roots[ index ]), 1));
    //     }
    //     else {
    //         m.first.SetLength(1);
    //         m.second.SetLength(1);
    //         NTL::set(m.first);
    //         // m.second = m.second.zero();
    //         NTL::SetCoeff(m.second, 0, 0);
    //         // return m;
    //     }
    // }
    else {

        // we do karatsuba multiplication by hand
        // int size = 1 << (max_depth - level + 2);
        // return custom_fp2X_mul(FastInterpolateFromRootsKaratsubaRec(roots, level + 1, 1 + (index << 1), max_depth, len), FastInterpolateFromRootsKaratsubaRec(roots, level + 1, index << 1, max_depth, len));
        FastInterpolateFromRootsKaratsubaRec(m, roots, level + 1, 1 + (index << 1), max_depth, len);
        std::pair<FpX, FpX> m2;
        FastInterpolateFromRootsKaratsubaRec(m2, roots, level + 1, index << 1, max_depth, len);
        explicit_FFT_custom_fp2X_mul(m, m, m2);

    }

}


void FastInterpolateFromRootsKaratsubaRec2(std::pair<FpX_elem, FpX_elem>& m, std::vector<ffp2> const &roots, int level, int index, int max_depth, int len) {

    if (level == max_depth - 1) {
        int new_index = index % 63;
        if (index < (len - 1) && new_index != 62) {
            // two roots
            m.first.SetLength(3);
            m.second.SetLength(2);
            // NTL::SetCoeff(m.first, 2);
            m.first[2] = Fp(1);
            
            // NTL::SetCoeff(m.first, 1, - (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0)));
            m.first[1] = - (roots[index].first + roots[index + 1].first);  
            // assert( (- (roots[new_index]._zz_pE__rep[0] + roots[new_index + 1]._zz_pE__rep[0])) ==  (- (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0))));
            // NTL::SetCoeff(m.second, 1, - (NTL::coeff(rep(roots[new_index]), 1) + NTL::coeff(rep(roots[new_index + 1]), 1)));
            m.second[1] = - (roots[index].second + roots[index + 1].second);  

            // karatsuba mul
            // Fp temp = (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index]), 1)) * (NTL::coeff(rep(roots[new_index + 1]), 0) + NTL::coeff(rep(roots[new_index + 1]), 1));
            Fp temp = (roots[index].first + roots[index].second) * (roots[index + 1].first + roots[index + 1].second);

            // Fp t1 =  (NTL::coeff(rep(roots[new_index]), 0) * NTL::coeff(rep(roots[new_index + 1]), 0));
            // Fp t2 =  (NTL::coeff(rep(roots[new_index]), 1) * NTL::coeff(rep(roots[new_index + 1]), 1));
            Fp t1 = roots[index].first * roots[index + 1].first;
            Fp t2 = roots[index].second * roots[index + 1].second;
            // NTL::SetCoeff(m.first, 0, t1 - t2);
            // NTL::SetCoeff(m.second, 0, temp - (t1 + t2));
            m.first[0] = t1 - t2;
            m.second[0] = temp - (t1 + t2);
            
        }
        else if (index == len - 1 || (new_index == 62 && index < len - 1)) {
            // only one root
                m.first.SetLength(2); // real part
                m.second.SetLength(1); // imaginary part

                m.first[1] = Fp(1);
                m.first[0] = - roots[index].first;
                m.second[0] = - roots[index].second;
                // NTL::SetCoeff(m.first, 1);
                // NTL::SetCoeff(m.first, 0, - NTL::coeff(rep(roots[new_index]),0));
                // NTL::SetCoeff(m.second, 0, - NTL::coeff(rep(roots[new_index]),1));
        }
        else {

            m.first.SetLength(1);
            m.second.SetLength(1);
            // NTL::SetCoeff(m.first, 0, 1);
            // NTL::set(m.first);
            m.first[0] = Fp(1);
            // NTL::SetCoeff(m.second, 0, 0);
            m.second = m.second.zero();
            // return m;
        }
    }
    // if (level == max_depth) {
    //     if (index < len) {
    //         // only one root
    //         m.first.SetLength(2); // real part
    //         m.second.SetLength(1); // imaginary part

    //         NTL::SetCoeff(m.first, 1, 1);
    //         NTL::SetCoeff(m.first, 0, - NTL::coeff(rep(roots[ index ]), 0));
    //         NTL::SetCoeff(m.second, 0, - NTL::coeff(rep(roots[ index ]), 1));
    //     }
    //     else {
    //         m.first.SetLength(1);
    //         m.second.SetLength(1);
    //         NTL::set(m.first);
    //         // m.second = m.second.zero();
    //         NTL::SetCoeff(m.second, 0, 0);
    //         // return m;
    //     }
    // }
    else {

    int step = max_depth - level;

    if (step <= 6) { 

        step = 1 << (step - 1);

        FastInterpolateFromRootsKaratsubaRec2(m, roots, level + 1, index + step, max_depth, len);
        std::pair<FpX, FpX> m2;
        FastInterpolateFromRootsKaratsubaRec2(m2, roots, level + 1, index, max_depth, len);
        explicit_FFT_custom_fp2X_mul(m, m, m2);
    }
    else {
        step = 63 * (1 << (step - 6 - 1) );

        FastInterpolateFromRootsKaratsubaRec2(m, roots, level + 1, index + step, max_depth, len);
        std::pair<FpX, FpX> m2;
        FastInterpolateFromRootsKaratsubaRec2(m2, roots, level + 1, index, max_depth, len);
        explicit_FFT_custom_fp2X_mul(m, m, m2);
    }
        

    }

}



void FastInterpolateFromRootsKaratsuba(std::pair<FpX_elem, FpX_elem> &poly, std::vector<ffp2> const & roots) {
    size_t n = roots.size();

    size_t k = ceil(log2(n));
    FastInterpolateFromRootsKaratsubaRec(poly, roots, 0, 0, k, n);
    
    
}

// // multiply p1 by (X - root)
// void fast_1_to_n_mul(std::pair<FpX_elem, FpX_elem> &res, std::pair<FpX_elem, FpX_elem> &p1, const ffp2 &root) {
    
//     // increase the degree by 1
//     int d = deg(p1.second);
//     // assert(deg(p1.first) == d + 1);
//     res.first.SetLength(deg(p1.first) + 2);
//     res.second.SetLength(d + 2);

//     // first coefficient
//     ffp2 new_co0,new_co1;
//     fast_mul(new_co0, root, {p1.first[0], p1.second[0]});
//     negate(new_co0, new_co0);

//     for (int i = 0; i < d; i++) {
//         fast_mul(new_co1, root, {p1.first[i + 1], p1.second[i + 1]});
//         sub(new_co1.first, new_co1.first, p1.first[i]);
//         sub(new_co1.second, new_co1.second, p1.second[i]);
        
//         res.first[i] = new_co0.first;
//         res.second[i] = new_co0.second;

//         negate(new_co0, new_co1);

//     }
//     for (int i = d; i < deg(p1.first); i++) {
//         // finishing
//         mul(new_co1.first, root.first, p1.first[i + 1]);
//         mul(new_co1.second, root.second, p1.first[i + 1]);
//         sub(new_co1.first, new_co1.first, p1.first[i]);
//         sub(new_co1.second, new_co1.second, p1.second[i]);

//         res.first[i] = new_co0.first;
//         res.second[i] = new_co0.second;

//         negate(new_co0, new_co1);
//     }

    
//     res.first[deg(p1.first) + 1] = p1.first[deg(p1.first)];
//     res.first[deg(p1.first)] = new_co0.first;
//     res.second[deg(p1.first)] = new_co0.second;

    

// }

// ffp2 FastInterpolateFromRootsKaratsubaRecPlusTrick(std::pair<FpX_elem, FpX_elem>& m, std::vector<ffp2> const &roots, int level, int index, int max_depth, int len) {

//     if (level == max_depth - 1) {
//         int new_index = index << 1;
//         int indmod64 = new_index % 64;
//         if (new_index < (len - 1) && indmod64 != 62) {
//             // two roots
//             m.first.SetLength(3);
//             m.second.SetLength(2);
//             // NTL::SetCoeff(m.first, 2);
//             m.first[2] = Fp(1);
            
//             // NTL::SetCoeff(m.first, 1, - (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0)));
//             m.first[1] = - (roots[new_index].first + roots[new_index + 1].first);  
//             // assert( (- (roots[new_index]._zz_pE__rep[0] + roots[new_index + 1]._zz_pE__rep[0])) ==  (- (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index + 1]), 0))));
//             // NTL::SetCoeff(m.second, 1, - (NTL::coeff(rep(roots[new_index]), 1) + NTL::coeff(rep(roots[new_index + 1]), 1)));
//             m.second[1] = - (roots[new_index].second + roots[new_index + 1].second);  

//             // karatsuba mul
//             // Fp temp = (NTL::coeff(rep(roots[new_index]), 0) + NTL::coeff(rep(roots[new_index]), 1)) * (NTL::coeff(rep(roots[new_index + 1]), 0) + NTL::coeff(rep(roots[new_index + 1]), 1));
//             Fp temp = (roots[new_index].first + roots[new_index].second) * (roots[new_index + 1].first + roots[new_index + 1].second);

//             // Fp t1 =  (NTL::coeff(rep(roots[new_index]), 0) * NTL::coeff(rep(roots[new_index + 1]), 0));
//             // Fp t2 =  (NTL::coeff(rep(roots[new_index]), 1) * NTL::coeff(rep(roots[new_index + 1]), 1));
//             Fp t1 = roots[new_index].first * roots[new_index + 1].first;
//             Fp t2 = roots[new_index].second * roots[new_index + 1].second;
//             // NTL::SetCoeff(m.first, 0, t1 - t2);
//             // NTL::SetCoeff(m.second, 0, temp - (t1 + t2));
//             m.first[0] = t1 - t2;
//             m.second[0] = temp - (t1 + t2);

//             return {Fp(0), Fp(0)};
            
//         }
//         else if (new_index < (len - 1) && indmod64 == 62) {
//             m.first.SetLength(2); // real part
//             m.second.SetLength(1); // imaginary part

//             m.first[1] = Fp(1);
//             m.first[0] = - roots[new_index].first;
//             m.second[0] = - roots[new_index].second;

//             return roots[new_index + 1];
//         }
//         else if (new_index == len - 1) {
//             // only one root
//                 m.first.SetLength(2); // real part
//                 m.second.SetLength(1); // imaginary part

//                 m.first[1] = Fp(1);
//                 m.first[0] = - roots[new_index].first;
//                 m.second[0] = - roots[new_index].second;
//                 // NTL::SetCoeff(m.first, 1);
//                 // NTL::SetCoeff(m.first, 0, - NTL::coeff(rep(roots[new_index]),0));
//                 // NTL::SetCoeff(m.second, 0, - NTL::coeff(rep(roots[new_index]),1));
//                 return {Fp(0), Fp{0}};
//         }
//         else {

//             m.first.SetLength(1);
//             m.second.SetLength(1);
//             // NTL::SetCoeff(m.first, 0, 1);
//             // NTL::set(m.first);
//             m.first[0] = Fp(1);
//             // NTL::SetCoeff(m.second, 0, 0);
//             m.second = m.second.zero();
//             // return m;
//             return {Fp(0), Fp{0}};
//         }
//     }
//     else {

//         // we do karatsuba multiplication by hand
//         // int size = 1 << (max_depth - level + 2);
//         // return custom_fp2X_mul(FastInterpolateFromRootsKaratsubaRec(roots, level + 1, 1 + (index << 1), max_depth, len), FastInterpolateFromRootsKaratsubaRec(roots, level + 1, index << 1, max_depth, len));
//         ffp2 r1 = FastInterpolateFromRootsKaratsubaRecPlusTrick(m, roots, level + 1, 1 + (index << 1), max_depth, len);
//         // once we're above the FFT threshold we need to multiply by (X - r1)
        
//         std::pair<FpX, FpX> m2;
//         ffp2 r2 = FastInterpolateFromRootsKaratsubaRecPlusTrick(m2, roots, level + 1, index << 1, max_depth, len);

//         if (deg(m2.first) >= 63 && deg(m2.first) == ((1 <<  (max_depth - level - 1) ) - 1)) {

//         // clock_t t = tic();
        
//             fast_1_to_n_mul(m2, m2, r2);


//         // std::cout << "1 dim mul " << tic() - t << "\n";
//         //    std::cout << "deg after " << deg(m2.first) << "\n";
//         }
//         explicit_FFT_custom_fp2X_mul(m, m, m2);
//         return r1;

//     }

// }

void FastInterpolateFromRootsKaratsubaPlusTrick(std::pair<FpX_elem, FpX_elem> &poly, std::vector<ffp2> const & roots) {
    size_t n = roots.size();

    size_t k1 = ceil(log2((n / 63) + 1));
    size_t k2 = ceil(log2(n));
    if (k1 + 6 > k2) {
        FastInterpolateFromRootsKaratsubaRec(poly, roots, 0, 0, k2, n);
    }
    else {
        FastInterpolateFromRootsKaratsubaRec2(poly, roots, 0, 0, k1 + 6, n);
    }
    
}

