///////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// The code in this file implements the main algorithms in the accompanying paper (header)
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "utils.hpp"
#include "endring.hpp"
#include "choosetorsion.hpp"
#include "quatlatenum.hpp"
#include "crt.hpp"
#include "klpt.hpp"
#include "ordertojinvbigset.hpp"

#include <list>

typedef NTL::ZZ Integer;
typedef NTL::ZZ_p Fp_big_elem;
typedef NTL::ZZ_pE FpE_big_elem;
typedef NTL::ZZ_pX FpX_big_elem;
typedef NTL::ZZ_pEX FpEX_big_elem;

std::vector<NTL::ZZ> _avail_qs(NTL::ZZ const &p, NTL::ZZ const &ell);

///// Big Level ////
std::vector<FpE_elem> SupersingularEvaluation(Fp_integer p, FpE_elem const &j, long l, bool useKLPT = true);
FpX_big_elem ModEvalBigLevel(NTL::ZZ p, NTL::ZZ_pE const j, long l);
std::vector<FpE_elem> SupersingularEvaluationWeber(Fp_integer p, FpE_elem w, long l);
FpX_big_elem ModEvalBigLevelWeber(NTL::ZZ p, NTL::ZZ_pE const w, long l);

//// Big Char /////
FpX SpecialSupersingularEvaluation(const Integer &p, const Integer &ell, const std::vector<Fp_elem> eval_points);
FpX SpecialSupersingularEvaluationWeber(const Integer &p, const Integer &ell, const std::vector<Fp_elem> eval_points);

Fp2 const_CRT_polynomials( const std::vector<Fp2> &values, const std::vector<Integer> &modulus, const long len);

FpX SpecialSupersingularEvaluationCRT(const Integer &p1, const Integer &p2, const Integer ell, const std::vector<Fp_elem> eval_points);
FpX_big_elem ModEvalBigCharacteristicWeber(NTL::ZZ p, Fp_big_elem const j, NTL::ZZ l);


