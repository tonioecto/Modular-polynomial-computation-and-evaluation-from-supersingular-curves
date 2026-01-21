#include "modpoly.hpp"

#include <NTL/ZZ.h>
#include <NTL/vec_lzz_pE.h>

#include "interpolation.hpp"
#include <unordered_set>
#include "fast_quaternions.hpp"

#include <cassert>
#include <ctime>




int main(int argc, char **argv)
{
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i)
        args.push_back(argv[i]);

    if (args.size() < 5) {
        bad_args:
        std::cerr << "usage:\n";
        std::cerr << "  " << argv[0] << " bigl j [p] [l] [j]\n";
        std::cerr << "  " << argv[0] << " bigl w [p] [l] [w]\n";
        std::cerr << "  " << argv[0] << " bigc w [p] [l] [w]\n";
        std::cerr << "  " << argv[0] << " ssse w [p] [l] [w]\n";
        return 1;
    }
    if (args[0] != "bigc" and args[0] != "bigl" and args[0] != "ssse" and args[0] !="bench")
        goto bad_args;
    if (args[1] != "j" and args[1] != "w")
        goto bad_args;


    NTL::ZZ p, l, w;
    NTL::ZZ_p inv;
    NTL::ZZ_pX F;
    if (!(std::istringstream(args[2]) >> p))
        goto bad_args;

    NTL::ZZ_p::init(p);
    NTL::ZZ_pX f;
    BuildIrred(f, 2);
    NTL::ZZ_pE::init(f);


    if (!(std::istringstream(args[3]) >> l))
        goto bad_args;
    if (!(std::istringstream(args[4]) >> w))
        goto bad_args;

    inv = NTL::conv<NTL::ZZ_p>(w);



    if (args[0] == "bigl" && args[1] == "j")
        F = ModEvalBigLevel(p, NTL::ZZ_pE(inv), NTL::conv<long>(l));
    else if (args[0] == "bigl" && args[1] == "w")
        F = ModEvalBigLevelWeber(p, NTL::ZZ_pE(inv), NTL::conv<long>(l));
    else if (args[0] == "bigc" && args[1] == "w") {
        clock_t t = clock();
        F = ModEvalBigCharacteristicWeber(p, inv, l);
        std::cout << "Total ModEvalBigCW time = " << (double) (clock() - t) / (CLOCKS_PER_SEC) << "\n";
    }
        
    else if (args[0] == "ssse" && args[1] == "w") {
        Fp_elem::init(NTL::conv<long>(p));
        Fp_elem tau;
        NTL::conv(tau, w);

        // computation of the eval list
        std::vector<Fp> eval_points = {};
        eval_points.push_back(Fp(1));
        eval_points.push_back(tau);
        Fp pow = tau;
        for (size_t i = 2; i<=l+Integer(1); i++) {
            pow = tau * pow;
            eval_points.push_back(pow);
        }


        auto Phi = SpecialSupersingularEvaluationWeber( p, l, eval_points);
        // std::cout << Phi << "\n";
        long ell = NTL::conv<long>(l);
        std::cout << "Phi = "<< Phi[0] << " " << Phi[1] << " "<< Phi[2] << " ........ " << Phi[ell - 1] << " " << Phi[ell] << " " << Phi[ell + 1] << "\n";
    }
    else if (args[0] == "bench") 
    {

    
        // this benches various stuff about finite field arithmetic 
        // a part of this is deprecated
        FastInteger p1 = 3011;
        // FastInteger p3 = 3019;
        // FastInteger p2 = 24407;
        FastInteger ell = 10;


        // finite field extension arithmetic
        {

            // NTL context for p1
            Fp::init(p1);
            FpX ff1;
            SetCoeff(ff1, 2);
            ff1[0] = Fp(1);
            Fp2::init(ff1);
            unsigned k_start = 1;
            unsigned k_bound = 20; //Some minimum
            std::cerr << "Generating field exts up to: " << k_bound << "..." << std::flush;
            //TODO eventually we should generate these field extensions on the fly and cache them, similar to how we now do it for torsion bases. at the moment this fails because references to the Fp2k object are stored in all kinds of other objects and those references are invalidated when the std::map is modified. possible solution: use std::shared_ptr for Fp2k references, just like we do for ec references.
            std::map<unsigned, Fp2k> Fexts;
            {
                for (unsigned k = k_start; k <= k_bound; ++k) {
                    std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... " << k << std::flush;
                    Fexts.emplace(std::make_pair(k, Fp2k {k}));
                }
                std::cerr << "\r\x1b[KGenerating field exts up to: " << k_bound << "... done" << std::endl;
            }


            std::map<unsigned, ffp2X> ffexts;
            {
                for (unsigned k = k_start; k <= k_bound; ++k) {
                    auto fpk = Fexts.find(k)->second;
                    Fp2X mod_poly_k = fpk.mod_Fp2;
                    std::cout <<" mod poly = " << mod_poly_k << "\n";
                    ffp2X ffp2_mod_poly = ffp2X(mod_poly_k);
                    ffexts.emplace( std::make_pair(k, ffp2_mod_poly) );
                }
            }


            for (unsigned k = k_start; k <= k_bound; ++k) {
                auto fpk = Fexts.find(k)->second;
                
                FpE_elem j;
                {   
                    
                    FpE_push Push(fpk.F);
                    j = NTL::random_zz_pE();
                    std::cout << "mod = " << j.modulus() << "\n";
                    FpE_elem pow = j;
                    clock_t tim = tic();
                    for (int i = 0; i < ell; i++) {
                        pow = j * pow;
                        sqrt(pow);

                    }
                    std::cout << "k = " << k << " NTL      time = " << tic() - tim << "  " << pow << "\n";

                }

                {
                    FpE_push Push(fpk.F);
                    FpE_elem pow = j;
                    clock_t tim = tic();
                    for (int i = 0; i < ell; i++) {
                        SpecialMul(pow, pow, j);
                        auto sq = fast_sqrt(fpk, pow);
#ifndef NDEBUG 
                        auto sqq = sqrt(fpk, pow);
                        if (sq) {
                            // std::cout << "is sq \n";
                            assert(sqq);
                            assert(*sq == *sqq || *sq == -(*sqq));
                        }
                        else {
                            // std::cout << "is not sq \n";
                            assert(!sqq);
                        }
#endif
                    }
                    std::cout << "k = " << k << " NTL spec time = " << tic() - tim << "  " << pow << "\n";
                }

                // {
                //     auto mod_poly = ffexts.find(k)->second;

                //     ffp2k j_rand = random_ffp2k(&mod_poly);
                //     j_rand.normalize();
                //     ffp2k pow = j_rand;
                //     clock_t tim = tic();
                //     for (int i = 0; i < ell; i++) {
                //         mul(pow , j_rand,  pow);
                //     }
                //     std::cout << "k = " << k << " custom   time = " << tic() - tim << "  " << pow.def_poly << "\n";

                // }

            }
        }


        // NTL context for p1
        // Fp::init(p1);
        // FpX ff1;
        // SetCoeff(ff1, 2);
        // ff1[0] = Fp(1);
        // Fp2::init(ff1);

        // Fp2 j1,j2;

        // std::vector<Fp2> roots1(ell + 1);
        // std::vector<ffp2> roots1_ffp2(ell + 1);
        // std::vector<ffp2> roots2(ell + 1);

        // // sampling random j1
        // j1 = NTL::random_zz_pE();

        // FpX t;
        // t.SetLength(2);

        // // Fp2 pow = j1;
        // ffp2 j1f = {rep(j1)[0], rep(j1)[1]};
        // ffp2 powf = j1f;
        // clock_t t2 = tic();
        // for (int i = 0; i < ell + 1; i++) {
        //     fast_mul(powf, j1f, powf);
        //     fast_add(powf, powf, j1f);
        //     fast_mul(powf, j1f, powf);
            
        // }
        // std::cout << "pair Fp mul        " << tic() - t2 << " [" << powf.first << " " <<powf.second << "]\n";
        // powf = j1f;
        // NTL::mulmod_t pinv = NTL::PrepMulMod(p1);
        // t2 = tic();
        // for (int i = 0; i < ell + 1; i++) {
        //     fast_fast_mul(powf, j1f, powf, p1, pinv);
        //     fast_add(powf, powf, j1f);
        //     fast_fast_mul(powf, j1f, powf, p1, pinv);
        // }
        // std::cout << "fast pair Fp mul   " << tic() - t2 << " [" << powf.first << " " <<powf.second << "]\n";

        // Fp2 pow = j1;
        // t2 = tic();
        // for (int i = 0; i < ell + 1; i++) {
        //     faster_mul(pow, j1, pow, t);
        //     fast_add(pow, pow, j1, t);
        //     faster_mul(pow, j1, pow, t);
        // }
        // std::cout << "custom Fp2 mul     " << tic() - t2 << " " << pow << "\n";

        // pow = j1;
        // clock_t t1 = tic();
        // for (int i = 0; i < ell + 1; i++) {
        //     NTL::mul(pow, j1, pow);
        //     add(pow, pow, j1);
        //     mul(pow, pow, j1);
        // }
        // std::cout << "native NTL Fp2 mul " << tic() - t1 << " " << pow << "\n";

        // for (int i = 0; i < 1; i++) {
        //     j1 = NTL::random_zz_pE();
        //     pow = j1;

        //     // creation of roots
        //     for (int i = 0; i < ell + 1; i++) {

        //         // NTL
        //         roots1[i] = pow;
        //         roots1_ffp2[i] = from_Fp2(pow);
        //         pow = pow * j1;

        //     }

        //     t2 = tic();
        //     std::pair<FpX_elem, FpX_elem> test1;
        //     FastInterpolateFromRootsKaratsuba(test1, roots1_ffp2);
        //     t2 = tic() - t2;

        //     std::cout << "\n\n";

        //     clock_t t3 = tic();
        //     std::pair<FpX_elem, FpX_elem> test11;
        //     FastInterpolateFromRootsKaratsubaPlusTrick(test11, roots1_ffp2);
        //     t3 = tic() - t3;
        
        //     // std::cout << "basic NTL res for p1 = " << f1[0] << " ... " << f1[ell-1] << ", " << f1[ell] << "\nin time = " << (double) t1 / CLOCKS_PER_SEC << "\n";
        //     std::cout << "hand NTL res for p1         = " << test1.first[0] << "," << test1.second[0] << " ... " << test1.first[ell-1] << "," << test1.second[ell-1] << ", " << test1.first[ell] <<","<< test1.second[ell]<< "\nin time = " << (double) t2 / CLOCKS_PER_SEC << "\n";
        //     std::cout << "hand NTL res for p1 + trick = " << test11.first[0] << "," << test11.second[0] << " ... " << test11.first[ell-1] << "," << test11.second[ell-1] << ", " << test11.first[ell] <<","<< test11.second[ell]<< "\nin time = " << (double) t3 / CLOCKS_PER_SEC << "\n";
        
        // }

        // // CRT interpolation from roots
        // {
        //     Fp::init(p2);
        //     FpX ff2;
        //     SetCoeff(ff2, 2);
        //     ff2[0] = Fp(1);
        //     Fp2::init(ff2);

        //     j2 = NTL::random_zz_pE();

        
        //     pow = j2;
        //     // creation of roots
        //     for (int i = 0; i < ell + 1; i++) {

        //         // NTL
        //         roots2[i] = from_Fp2(pow);
        //         pow = pow * j2;

        //     }

        //     t1 = tic();
        //     std::pair<FpX_elem, FpX_elem> f2;
        //     FastInterpolateFromRootsKaratsuba(f2, roots2);
        //     t1 = tic() - t1;

        //     std::cout << "hand NTL res for p2 in time = " << (double) t1 / CLOCKS_PER_SEC << "\n";

        //     Fp::init(p3);
        //     FpX ff3;
        //     SetCoeff(ff3, 2);
        //     ff3[0] = Fp(1);
        //     Fp2::init(ff3);
        //     std::vector<ffp2> roots3(ell + 1);

        //     Fp2 j3 = NTL::random_zz_pE();
        //     pow = j3;
        //     // creation of roots
        //     for (int i = 0; i < ell + 1; i++) {

        //         // NTL
        //         roots3[i] = from_Fp2(pow);
        //         pow = pow * j3;
        //     }


        //     t2 = tic();
        //     std::pair<FpX_elem, FpX_elem> f3;
        //     FastInterpolateFromRootsKaratsuba(f3, roots3);
        //     t2 = tic() - t2;

        //     std::cout << "hand NTL res for p3 = in time = " << (double) t2 / CLOCKS_PER_SEC << "\n";


        //     auto m = p1 * p2 * p3;
        //     Fp_integer m_mod;
        //     NTL::conv(m_mod, m);
        //     Fp::init(m_mod);
        //     Fp::init(m);
        //     FpX fm;
        //     SetCoeff(fm,2);
        //     fm[0] = Fp(1);
        //     Fp2::init(fm);

        //     std::vector<ffp2> rootsm(ell + 1);

        //     t2 = tic();

        //     // constituting 
        //     std::vector<FastInteger> modulus = {p1, p2, p3};
        //     std::vector<FastInteger> mod_prod = {p1, p1 * p2, p1 * p2 * p3};
        //     std::vector<FastInteger> inv(2);
        //     auto red2 = SignedBarrettReducer(p2);
        //     auto red3 = SignedBarrettReducer(p3);
        //     inv[0] = InvMod(p1, red2);
        //     inv[1] = InvMod(p1 * p2, red3);

        //     // creation of roots
        //     for (int i = 0; i < ell + 1; i++) {

        //         // NTL
        //         rootsm[i] = CRT_ffp2( {roots1_ffp2[i], roots2[i], roots3[i]}, modulus, mod_prod, inv, 3);

        //     }

        //     std::cout << p1 << "  " << p2 << " " << p3 << "\n";
        //     std::cout << roots1_ffp2[0] << " " << roots2[0] << " " << roots3[0] << " " << rootsm[0] << "\n";

        //     std::pair<FpX_elem, FpX_elem> testm;
        //     testm.first.SetLength(ell + 1);
        //     testm.first.SetLength(ell + 1);
        //     FastInterpolateFromRootsKaratsubaPlusTrick(testm, rootsm);
        //     t2 = tic() - t2;

        //     std::cout << "CRT + hand NTL res mod m = in time = " << (double) t2 / CLOCKS_PER_SEC << "\n";

        // }

        





    }
    else
        goto bad_args;

    for(int i = 0; i <= NTL::deg(F); i++) {
        if (i < 3 || i >= NTL::deg(F) - 3) {
            std::cout << "Coefficient of x^" << i << " is: " << NTL::coeff(F,i) << "\n";
        }
    }
        
    std::cout << std::flush;

    return 0;
}
