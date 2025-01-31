// Computing Algorithm 4 and 5 in eval.tex

// All CRT cases treated as "small" in the sense of section 6 of "Computing Hilbert class polynomials with the Chinese Remainder Theorem")

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZX.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <time.h>
#include "crt.hpp"
#include "interpolation.hpp"

// TODO: use arrays if the vector has fixed length, otherwise use std::vector with push_back or resize


// just adding this function so that we can run tests without having the underlying elliptic curve function stuff -- we can remove it after

void for_testing_primes(std::vector<NTL::ZZ> &Pl){

    Pl.push_back(NTL::ZZ(1031));
    Pl.push_back(NTL::ZZ(2063));
    Pl.push_back(NTL::ZZ(8627));

}


void for_testing_js(std::vector<NTL::ZZX> &js_polys, std::vector<std::vector<NTL::ZZX>> &iso_js_polys, NTL::ZZ const q){
    
    // Need to have 4 for a degree 3 poly
    NTL::ZZX f1, f2, f3, f4, f5, f6, f7, f8, f9, g1, g2, g3;
    SetCoeff(f1, 1); SetCoeff(f2, 1); SetCoeff(f3, 1);
    SetCoeff(f4, 1); SetCoeff(f5, 1); SetCoeff(f6, 1);
    SetCoeff(f7, 1); SetCoeff(f8, 1); SetCoeff(f9, 1);
    SetCoeff(g1, 1); SetCoeff(g2, 1); SetCoeff(g3, 1);
    
    f1[1] = 603;
    f1[0] = 89;
    f2[1] = 719;
    f2[0] = 976;
    f3[1] = 0;
    f3[0] = 960;
    g1[1] = 0;
    g1[0] = 615;

    f4[1] = 1094;
    f4[0] = 206;
    f5[1] = 0;
    f5[0] = 46;
    f6[1] = 0;
    f6[0] = 358;
    g2[1] = 0;
    g2[0] = 1606;

    f7[1] = 5375;
    f7[0] = 3219;
    f8[1] = 1823;
    f8[0] = 3933;
    f9[1] = 7661;
    f9[0] = 5722;
    g3[1] = 3292;
    g3[0] = 6294;

    if(q == NTL::ZZ(1031)){
        js_polys.push_back(f1); 
        js_polys.push_back(f2); 
        js_polys.push_back(f3);
        js_polys.push_back(g1);
    }
    else if(q == NTL::ZZ(2063)){
        js_polys.push_back(f4); 
        js_polys.push_back(f5); 
        js_polys.push_back(f6);
        js_polys.push_back(g2);
    }
    else{
        js_polys.push_back(f7); 
        js_polys.push_back(f8); 
        js_polys.push_back(f9);
        js_polys.push_back(g3);
    }

    NTL::ZZX f11, f12, f13, f21,f22,f23, f31,f32,f33,f41,f42,f43,f51,f52,f53;
    NTL::ZZX f61,f62,f63,f71,f72,f73,f81,f82,f83,f91,f92,f93;
    NTL::ZZX g11, g12, g13, g21, g22, g23, g31, g32, g33;

    SetCoeff(f11, 1); SetCoeff(f12, 1); SetCoeff(f13, 1);
    SetCoeff(f21, 1); SetCoeff(f22, 1); SetCoeff(f23, 1);
    SetCoeff(f31, 1); SetCoeff(f32, 1); SetCoeff(f33, 1);
    SetCoeff(f41, 1); SetCoeff(f42, 1); SetCoeff(f43, 1);
    SetCoeff(f51, 1); SetCoeff(f52, 1); SetCoeff(f53, 1);
    SetCoeff(f61, 1); SetCoeff(f62, 1); SetCoeff(f63, 1);
    SetCoeff(f71, 1); SetCoeff(f72, 1); SetCoeff(f73, 1);
    SetCoeff(f81, 1); SetCoeff(f82, 1); SetCoeff(f83, 1);
    SetCoeff(f91, 1); SetCoeff(f92, 1); SetCoeff(f93, 1);
    SetCoeff(g11, 1); SetCoeff(g12, 1); SetCoeff(g13, 1);
    SetCoeff(g21, 1); SetCoeff(g22, 1); SetCoeff(g23, 1);
    SetCoeff(g31, 1); SetCoeff(g32, 1); SetCoeff(g33, 1);


    f11[1] = 719;
    f11[0] = 976;
    f12[1] = 119;
    f12[0] = 537;
    f13[1] = 9;
    f13[0] = 225;

    f21[1] = 0;
    f21[0] = 960;
    f22[1] = 0;
    f22[0] = 8;
    f23[1] = 603;
    f23[0] = 89;

    f31[1] = 0;
    f31[0] = 615;
    f32[1] = 719;
    f32[0] = 976;
    f33[1] = 312;
    f33[0] = 976;

    g11[1] = 0;
    g11[0] = 960;
    g12[1] = 0;
    g12[0] = 783;
    g13[1] = 0;
    g13[0] = 708;



    f41[1] = 0;
    f41[0] = 46;
    f42[1] = 1656;
    f42[0] = 1844;
    f43[1] = 1236;
    f43[0] = 935;

    f51[1] = 0;
    f51[0] = 358;
    f52[1] = 1094;
    f52[0] = 206;
    f53[1] = 969;
    f53[0] = 206;

    f61[1] = 0;
    f61[0] = 1606;
    f62[1] = 0;
    f62[0] = 1305;
    f63[1] = 0;
    f63[0] = 46;

    g21[1] = 0;
    g21[0] = 1635;
    g22[1] = 0;
    g22[0] = 868;
    g23[1] = 0;
    g23[0] = 358;





    f71[1] = 1823;
    f71[0] = 3933;
    f72[1] = 5203;
    f72[0] = 1998;
    f73[1] = 2183;
    f73[0] = 1028;

    f81[1] = 7661;
    f81[0] = 5722;
    f82[1] = 5375;
    f82[0] = 3219;
    f83[1] = 2116;
    f83[0] = 1546;

    f91[1] = 3292;
    f91[0] = 6294;
    f92[1] = 5121;
    f92[0] = 6006;
    f93[1] = 1823;
    f93[0] = 3933;

    g31[1] = 2352;
    g31[0] = 7906;
    g32[1] = 7661;
    g32[0] = 5722;
    g33[1] = 4469;
    g33[0] = 2663;

    if(q == NTL::ZZ(1031)){
        iso_js_polys.push_back({f11,f12,f13}); 
        iso_js_polys.push_back({f21,f22,f23}); 
        iso_js_polys.push_back({f31,f32,f33});
        iso_js_polys.push_back({g11,g12,g13});
    }
    else if(q == NTL::ZZ(2063)){
        iso_js_polys.push_back({f41,f42,f43}); 
        iso_js_polys.push_back({f51,f52,f53}); 
        iso_js_polys.push_back({f61,f62,f63});
        iso_js_polys.push_back({g21,g22,g23});
    }
    else{
        iso_js_polys.push_back({f71,f72,f73}); 
        iso_js_polys.push_back({f81,f82,f83}); 
        iso_js_polys.push_back({f91,f92,f93});
        iso_js_polys.push_back({g31,g32,g33});
    }

}

// TODO: replace with michaels stuff 
void get_primes(std::vector<NTL::ZZ> &Pl, NTL::ZZ const B, NTL::ZZ const ell){
    NTL::ZZ tmp = NTL::ZZ(1);
    NTL::ZZ q = NTL::NextPrime(ell+1);
    while(tmp < B){
        tmp *= q;
        Pl.push_back(q);
        q = NTL::NextPrime(q+1);
    }
}

NTL::ZZ_pX jInv_interpolate(std::vector<NTL::ZZ_pE> const js, std::vector<std::vector<NTL::ZZ_pE>> const iso_js, std::vector<NTL::ZZ_pE> const xs, NTL::ZZ ell){
    // Input the j-invariants js and the isogeneous j-invariants iso_js, where iso_js[i] are the j-invariants
    // that are ell-isogenous to js[i],  x0,..,xl values in Fp and the level ell of the modular polynomial
    int ell_int = NTL::conv<int>(ell);

    std::vector<NTL::ZZ_pE> outputs(ell_int+2);

    for(int i = 0; i <= ell_int+1; i++){
        // Kind of annoying: buildfromroots needs NTL::vec_ZZ_p and so lets do some adhoc conversion... Maybe we need to find a better way to do this
        NTL::vec_ZZ_pE rts;
        for(int j = 0; j <= ell; j++){
            rts.append(iso_js[i][j]);
        }
        
        NTL::ZZ_pEX F = NTL::BuildFromRoots(rts);

        outputs[i] = NTL::ZZ_pE(0);
        for(int j = 0; j <= deg(F); j++){
            NTL::ZZ_pE tmp;
            tmp = NTL::coeff(F, j);
            tmp *= xs[j];
            outputs[i] += tmp;
        }
    }

    NTL::ZZ_pEX Ptmp = fast_interpolate(js, outputs);

    NTL::ZZ_pX Phi;
    NTL::SetCoeff(Phi, deg(Ptmp));
    for(int i = 0; i <= deg(Ptmp); i++){
        NTL::ZZ_pX f = NTL::conv<NTL::ZZ_pX>(NTL::coeff(Ptmp, i));
        NTL::SetCoeff(Phi, i, NTL::LeadCoeff(f));
    }
    return Phi;
}




// lets input the js and iso_js for now to test -- TODO: replace with jonathan+lorenz's stuff

void special_ss_eval(std::vector<NTL::ZZ> &Fq_coeffs, NTL::ZZ const q, std::vector<NTL::ZZ> const xs, NTL::ZZ const ell){
    // Argh we need to be able to handle different moduli to use BuildFromRoots in interpolate function

    int ell_int = NTL::conv<int>(ell);

    // This is just here while we can't use the real function to compute the js and iso_js
    NTL::ZZ_pPush push(q);
    NTL::ZZ_pX f;
    SetCoeff(f, 2);
    f[0] = NTL::ZZ_p(1);
    NTL::ZZ_pE::init(f);

    std::vector<NTL::ZZX> js_polys;
    std::vector<std::vector<NTL::ZZX>> iso_js_polys;
    for_testing_js(js_polys, iso_js_polys, q);


    std::vector<NTL::ZZ_pE> js(ell_int+2);
    std::vector<std::vector<NTL::ZZ_pE>> iso_js(ell_int+2);
    for(int i = 0; i < ell_int+2; i++){
        js[i] = NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(js_polys[i]));
        for(int j = 0; j < iso_js_polys[0].size(); j++){
            iso_js[i].push_back(NTL::conv<NTL::ZZ_pE>(NTL::conv<NTL::ZZ_pX>(iso_js_polys[i][j])));
        }
    }     


    std::vector<NTL::ZZ_pE> xqs(ell_int+2);
    for(int i = 0; i <= ell_int+1; i++){
        xqs[i] = NTL::conv<NTL::ZZ_pE>(xs[i]);
    }
    
    NTL::ZZ_pX Ftmp = jInv_interpolate(js, iso_js, xqs, ell);


    for(int i = 0; i <= deg(Ftmp); i++){

        Fq_coeffs[i] = NTL::conv<NTL::ZZ>(NTL::coeff(Ftmp, i));
        // We view the coeffs as being in NTL::ZZ as we want to go back to working with modulus p
    
    }

}   


// Not sure we can test this until we have jonathan and lorenz's functions so let's postpone until thins

void biglevel(NTL::ZZ_pX &F, NTL::ZZ const p, NTL::ZZ_p const j, NTL::ZZ const ell){

    // Set bound as in paper Alg 3
    NTL::ZZ b = 6*ell*ceil(log(ell))+18*ell+(ell+1)*ceil(log(p))+ceil(log(ell+2));
    NTL::ZZ B;
    NTL::power(B, NTL::ZZ(2), NTL::conv<int>(b)); 


    // Number of coeffs
    int ell_int = NTL::conv<int>(ell);
    int Ncoeffs = ell_int + 2;

    // Constructing js = [j^i mod p for i \in [1, ell+1]]
    std::vector<NTL::ZZ> js(ell_int+2); //Want this to be over Fp?
    js[0] = NTL::ZZ(1);
    for(int i = 1; i <= ell_int+1; i++){
        js[i] = js[i-1]*NTL::conv<NTL::ZZ>(j);
    } 

    // Computing the set of primes 
    std::vector<NTL::ZZ> Pl; 
    // get_primes(Pl, B, ell);

    // for testing:
    for_testing_primes(Pl);

    int Nprimes = Pl.size();

    // Initialise crt structure
    crt_info crt;
    crt_init(crt, Pl, Nprimes, Ncoeffs, p);

    // Compute the F mod q and update crt coeffs
    for(int i = 0; i < Nprimes; i++){

        NTL::ZZ q = Pl[i];
        // In this function we set this to be in ZZX to not work with two moduli in a function
        std::vector<NTL::ZZ> Fq_coeffs(crt.k); 

        
        // Compute Fq 
        special_ss_eval(Fq_coeffs, q, js, ell);
        
        //correct output up to here
    

        //Update CRT sums
        crt_update(crt, i, Fq_coeffs, crt.k);
    }

    crt_finalise(crt);
    for(int i = 0; i < crt.k; i++){
        NTL::SetCoeff(F, i, NTL::conv<NTL::ZZ_p>(crt.Cdata[i])); 
    } 

}