
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// This file includes functions for extension field arithemetic via the structure Fp2k
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>

#include <NTL/lzz_p.h>
#include <NTL/lzz_pE.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pEXFactoring.h>

#include <NTL/vec_lzz_p.h>
#include <NTL/mat_lzz_p.h>

#include <optional>


/*typedef NTL::ZZ_p Fp_elem;
typedef NTL::ZZ_pE FpE_elem;
typedef NTL::ZZ_pX FpX_elem;
typedef NTL::ZZ_pEX FpEX_elem;
typedef NTL::ZZ_pEContext FpE_context;
typedef NTL::ZZ_pEXModulus FpEX_modulus;
typedef NTL::ZZ_pPush Fp_push;
typedef NTL::ZZ_pEPush FpE_push;

typedef NTL::vec_ZZ_p vec_Fp;
typedef NTL::mat_ZZ_p mat_Fp;

typedef NTL::ZZ Fp_integer; //Should be used for when passing the modulus to Fp_elem

inline Fp_elem random_Fp_elem() { return NTL::random_ZZ_p();}
inline FpE_elem random_FpE_elem() { return NTL::random_ZZ_pE();}
*/

// We chose to use zz_p rather that ZZ_p, etc. (for efficiency reasons)
typedef NTL::zz_p Fp_elem;
typedef NTL::zz_pE FpE_elem;
typedef NTL::zz_pX FpX_elem;
typedef NTL::zz_pEX FpEX_elem;
typedef NTL::zz_pEContext FpE_context;
typedef NTL::zz_pEXModulus FpEX_modulus;
typedef NTL::zz_pPush Fp_push;
typedef NTL::zz_pEPush FpE_push;

typedef NTL::vec_zz_p vec_Fp;
typedef NTL::mat_zz_p mat_Fp;

typedef long Fp_integer; //Should be used for when passing the modulus to Fp_elem

inline Fp_elem random_Fp_elem() { return NTL::random_zz_p();}
inline FpE_elem random_FpE_elem() { return NTL::random_zz_pE();}


size_t myroots(FpE_elem *roots, FpEX_elem const &f);

//Some functions in utils are hardcoded for this...
size_t myroots_ZZ(NTL::ZZ_pE *roots, NTL::ZZ_pEX const &f);

std::optional<FpE_elem> sqrt(FpE_elem const &alpha);

struct Fp2k
{
    FpE_context F;        // The field extension
    FpE_elem Fp2_gen;         // The (image of the) generator of Fp2
    unsigned k;                 // The degree as an extension of Fp2
    size_t Fp2_gen_nonzero;
    FpEX_elem mod_Fp2;
    mat_Fp frob_action;    // For large fields this uses quite a bit of memory, but probably worth it
    Fp2k(unsigned degree);
    Fp2k() = delete;

    FpE_elem frob(FpE_elem alpha) const;

    // Compute the iota map when generating the fields, as a reduction of curve with CM
    // Could consider doing this only for Fp2 (currently this essentially precomputes the lifting for each field extension, using alot of memory)
    FpE_elem starting_a;
    FpE_elem starting_b;

    FpEX_elem iota_x_num;
    FpEX_elem iota_x_denom;
    FpEX_elem iota_y_num;
    FpEX_elem iota_y_denom;

    bool maximal; // whether to use iota or 2iota + 1

    size_t iota_degree;
};

std::optional<FpE_elem> sqrt(Fp2k Fext, FpE_elem const &alpha);
FpE_elem lift(FpE_elem const &alpha, Fp2k const &Fext);
FpE_elem coerce(FpE_elem const &alpha, Fp2k const &Fext);
FpEX_elem lift(FpEX_elem const &f, Fp2k const &Fext);
FpEX_elem coerce(FpEX_elem const &f, Fp2k const &Fext);
FpEX_elem MinPoly(FpE_elem const &alpha, Fp2k const &Fext);

