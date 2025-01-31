///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////// The code in this file implements the conversion between NTL and GMP data types
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <cassert>
#include <gmp.h>
#include <NTL/ZZ.h>

void ntl2gmp(mpz_t out, NTL::ZZ const &num)
{
    thread_local std::vector<unsigned char> bs;

    if (NTL::IsZero(num)) {
        mpz_set_ui(out, 0);
        return;
    }

    size_t l = NTL::NumBytes(num);
    if (l > bs.size())
        bs.resize(l);

    int sgn = NTL::sign(num);
    assert(sgn == +1 || sgn == -1);

    NTL::BytesFromZZ(bs.data(), num, l);
    mpz_import(out, l, -1, 1, 0, 0, bs.data());

    if (sgn < 0)
        mpz_neg(out, out);
}

void gmp2ntl(NTL::ZZ &out, mpz_t const num)
{
    thread_local std::vector<unsigned char> bs;

    if (!mpz_cmp_ui(num, 0)) {
        NTL::clear(out);
        return;
    }

    size_t l = (mpz_sizeinbase(num, 2) + 7) / 8;
    if (l > bs.size())
        bs.resize(l);

    int sgn = mpz_sgn(num);
    assert(sgn == +1 || sgn == -1);

    size_t ll;
    mpz_export(bs.data(), &ll, -1, 1, 0, 0, num);
    assert(ll <= l);

    NTL::ZZFromBytes(out, bs.data(), ll);

    if (sgn < 0)
        out = -out;
}

