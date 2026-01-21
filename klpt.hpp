///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements the KLPT algorithm (and variants) -- header
///////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <optional>

#include "quaternions.hpp"
#include "quatlatenum.hpp"


//Bunch of subroutines required in klpt

// enumerates small solutions (x,y) to A*x + B*y = C (mod N)
void small_solutions(NTL::ZZ const &N, NTL::ZZ const &C, NTL::ZZ const &A, NTL::ZZ const &B,
                     std::function<bool(NTL::ZZ const &x, NTL::ZZ const &y)> const &fun, double maxdist);

quat RepresentInteger(quatalg const &B, NTL::ZZ const &M);
quat DetRepresentInteger(quatalg const &B, NTL::ZZ const &M);
std::optional<std::pair<NTL::ZZ, NTL::ZZ>> Cornacchia(const NTL::ZZ d,const NTL::ZZ &M);
std::pair<NTL::ZZ, NTL::ZZ> IdealModConstraint(quatlat const &I, quat const &gamma);
std::pair<quat, NTL::ZZ> KLPT(quatlat const &J, factor_list const &fac_list, NTL::ZZ const &coprime = NTL::ZZ(1));
std::pair<quat, NTL::ZZ> KLPT_conj(quat *conj, quatlat const &J, factor_list const &fac_list, NTL::ZZ const &coprime = NTL::ZZ(1));
