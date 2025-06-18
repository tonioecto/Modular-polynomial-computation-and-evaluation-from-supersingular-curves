//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// This file includes functions for extension field arithemetic via the structure Fp2k
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <optional>
#include "Fp2k.hpp"

size_t myroots(FpE_elem *roots, FpEX_elem const &f)
{
    // WARNING: This is only marginally faster than NTL::CanZass().
    assert(NTL::IsOne(NTL::LeadCoeff(f)));
    FpEX_elem x; NTL::SetX(x);
    auto sfd = NTL::SquareFreeDecomp(f);
    size_t idx = 0;
    for (auto const &[g,m]: sfd) {
        FpEX_modulus F;
        NTL::build(F, g);
        auto h = NTL::PowerXMod(FpE_elem::cardinality(), F);
        auto g1 = NTL::GCD(g, h - x);
        for (auto const &r: NTL::FindRoots(g1))
            for (ssize_t i = 0; i < m; ++i)
                roots[idx++] = r;
    }
    return idx;
}

size_t myroots_ZZ(NTL::ZZ_pE *roots, NTL::ZZ_pEX const &f)
{
    // WARNING: This is only marginally faster than NTL::CanZass().
    assert(NTL::IsOne(NTL::LeadCoeff(f)));
    NTL::ZZ_pEX x; NTL::SetX(x);
    auto sfd = NTL::SquareFreeDecomp(f);
    size_t idx = 0;
    for (auto const &[g,m]: sfd) {
        NTL::ZZ_pEXModulus F;
        NTL::build(F, g);
        auto h = NTL::PowerXMod(NTL::ZZ_pE::cardinality(), F);
        auto g1 = NTL::GCD(g, h - x);
        for (auto const &r: NTL::FindRoots(g1))
            for (ssize_t i = 0; i < m; ++i)
                roots[idx++] = r;
    }
    return idx;
}

std::optional<FpE_elem> sqrt(FpE_elem const &alpha)
    // Take EXTREME caution to only call this when ZZ_pE agrees with alpha
    // This is the same as with all arithmetic in field extensions...
{
    FpEX_elem f;
    FpE_elem roots[2];
    SetCoeff(f, 2);
    f[0] = -alpha;
    if (myroots(roots, f) == 0) {
        return {};
    }
    return roots[0];
}

std::optional<FpE_elem> sqrt(Fp2k Fext, FpE_elem const &alpha)
{
    FpE_push Push(Fext.F);
    return sqrt(alpha);
}

FpX_elem _construct_f_from_f2(FpEX_elem const &f2) {
    FpEX_elem f2_conj;
    size_t deg = NTL::deg(f2);
    auto p = Fp_elem::modulus();

    NTL::SetCoeff(f2_conj, deg);
    for (size_t i = 0 ; i < deg ; i++) {
        f2_conj[i] = NTL::power(f2[i], p);
    }

    FpEX_elem f_Fp2 = f2*f2_conj;
    FpX_elem f;
    size_t deg_f = NTL::deg(f_Fp2);

    NTL::SetCoeff(f, deg_f);
    for (size_t i = 0 ; i < deg_f ; i++) {
        if (f_Fp2[i] != 0) {
            f[i] = NTL::rep(f_Fp2[i])[0];
        }
    }

    return f;
}

mat_Fp _compute_frobenius_matrix(unsigned degree) {
    mat_Fp frob_action;
    NTL::ZZ p;
    NTL::conv(p, Fp_elem::modulus());
    frob_action.SetDims(degree, degree);
    for (unsigned i = 0 ; i < degree ; i++) {
        FpX_elem i_th_basis_poly;
        NTL::SetCoeff(i_th_basis_poly, i);
        FpE_elem i_th_basis;
        NTL::conv(i_th_basis, i_th_basis_poly);
        FpE_elem i_th_image = NTL::power(i_th_basis, p);
        FpX_elem image_poly = NTL::rep(i_th_image);
        for (unsigned j = 0 ; j < degree ; j++) {
            if (NTL::deg(image_poly) >= (long) j) {
                frob_action[i][j] = NTL::rep(i_th_image)[j];
            } else {
                frob_action[i][j] = 0;
            }
        }
    }
    return frob_action;
}

// ZZ_pE must be Fp2
Fp2k::Fp2k(unsigned degree)
{
    assert(degree > 0);
    assert(FpE_elem::degree() == 2);

    this->k = degree;

    FpE_context F;
    if (degree == 1) {
        F.save();
        FpX_elem Fp2_gen_poly;
        FpE_elem Fp2_gen;
        SetCoeff(Fp2_gen_poly, 1);
        conv(Fp2_gen, Fp2_gen_poly);
        this->F = std::move(F);
        this->Fp2_gen = Fp2_gen;
        this->Fp2_gen_nonzero = 1;
        this->frob_action = _compute_frobenius_matrix(2);
    }

    if (degree > 1) {
        FpEX_elem f_2;
        FpX_elem f;
        bool found = false;
        for (size_t i = 0 ; i < 1000 ; i++) {
            NTL::BuildIrred(f_2, degree);
            f = _construct_f_from_f2(f_2);

            if (NTL::SquareFreeDecomp(f)[0].b == 1) { //<- check if irreducible
                this->mod_Fp2 = f_2;
                found = true;
                break;
            }
        }
        if (!(found)) {
            throw;
        }

        auto Fp2_modulus = FpE_elem::modulus();
        assert(NTL::deg(Fp2_modulus) == 2);

        FpE_elem roots[2];
        { //init Fp2k temp
            FpE_push push(f);
            F.save();
            this->frob_action = _compute_frobenius_matrix(2*k);

            FpEX_elem mod;

            // Convert Fp2 modulus to a polynomial over Fp2k
            NTL::conv(mod, Fp2_modulus);

            // Solve to find a generator
            myroots(roots, mod);

            this->F = std::move(F);

            // Have to pick the right one w.r.t the modulus we chose earlier
            auto r1 = roots[0];
        }

        FpEX_elem r0;
        FpX_elem i_poly;
        NTL::SetCoeff(i_poly, 1);
        FpE_elem i;
        NTL::conv(i, i_poly);

        NTL::conv(r0, NTL::rep(roots[0]));
        FpEX_elem r0_minpoly = NTL::MinPolyMod(r0 % this->mod_Fp2, this->mod_Fp2);
        if (r0_minpoly[0] == i) {
            this->Fp2_gen = roots[1]; // something got turned upside down, but w/e it works now
        } else {
            this->Fp2_gen = roots[0];
        }
        for (Fp2_gen_nonzero = 1; NTL::IsZero(NTL::rep(Fp2_gen)[Fp2_gen_nonzero]); ++Fp2_gen_nonzero);
    }

    std::pair<FpE_elem, FpE_elem> starting_curve_coeffs;
    FpEX_elem xmap_up;
    FpEX_elem xmap_down;
    FpEX_elem ymap_up;
    FpEX_elem ymap_down;
    {
        NTL::ZZ p;
        NTL::conv(p, Fp_elem::modulus());
        FpE_push push(this->F);

        FpE_elem w = this->Fp2_gen;
        FpE_elem q = -(w*w);
        if (q == 1) {
            this->iota_degree = 1;
            this->starting_a = 1;
            this->starting_b = 0;
            assert (w*w == -1);

            NTL::SetCoeff(xmap_up, 1);
            NTL::SetCoeff(xmap_down, 0);
            NTL::SetCoeff(ymap_up, 0);
            NTL::SetCoeff(ymap_down, 0);

            xmap_up[1] = -1;
            ymap_up[0] = w;

            this->maximal = true;
        } else if (q == 3) {
            this->iota_degree = 1;
            this->starting_a = 0;
            this->starting_b = 1;
            assert (w*w == -3);

            NTL::SetCoeff(xmap_up, 1);
            NTL::SetCoeff(xmap_down, 0);
            NTL::SetCoeff(ymap_up, 0);
            NTL::SetCoeff(ymap_down, 0);

            xmap_up[1] = (w-1)/2;

            this->maximal = false;
        } else if (q == 7) {
            this->iota_degree = 2;
            this->starting_a = -35;
            this->starting_b = 98;

            assert (w*w == -7);
            NTL::SetCoeff(xmap_up, 2);
            NTL::SetCoeff(xmap_down, 1);
            NTL::SetCoeff(ymap_up, 2);
            NTL::SetCoeff(ymap_down, 2);
            this->maximal = false;

            xmap_up[0] = 35*w - 63;
            xmap_up[0] /= 8;
            xmap_up[1] = w + 7;
            xmap_up[1] /= 4;
            xmap_up[2] = -w - 3;
            xmap_up[2] /= 8;
            xmap_down[0] = w - 7;
            xmap_down[0] /= 2;
            xmap_down[1] = 1;
            ymap_up[0] = -21*w - 119;
            ymap_up[0] /= 16;
            ymap_up[1] = -3*w + 7;
            ymap_up[1] /= 4;
            ymap_up[2] = w - 5;
            ymap_up[2] /= 16;
            ymap_down[0] = -7*w + 21;
            ymap_down[0] /= -2;
            ymap_down[1] = -w + 7;
            ymap_down[2] = -1;
        } else if (q == 11) {
            this->iota_degree = 3;
            this->starting_a = -1056;
            this->starting_b = 13552;

            assert (w*w == -11);
            NTL::SetCoeff(xmap_up, 3);
            NTL::SetCoeff(xmap_down, 2);
            NTL::SetCoeff(ymap_up, 3);
            NTL::SetCoeff(ymap_down, 3);
            this->maximal = false;

            xmap_up[0] = 20768*w + 73568;
            xmap_up[0] /= 9;
            xmap_up[1] = -352*w - 1936;
            xmap_up[1] /= 3;
            xmap_up[2] = -4*w + 44;
            xmap_up[2] /= 3;
            xmap_up[3] = w - 5;
            xmap_up[3] /= 18;
            xmap_down[0] = 88*w + 440;
            xmap_down[1] = -4*w - 44;
            xmap_down[2] = 1;
            ymap_up[0] = -24640*w + 15488;
            ymap_up[0] /= 27;
            ymap_up[1] = 88*w - 2024;
            ymap_up[1] /= 9;
            ymap_up[2] = 10*w + 22;
            ymap_up[2] /= 3;
            ymap_up[3] = -w - 4;
            ymap_up[3] /= 27;
            ymap_down[0] = 2816*w + 7744;
            ymap_down[1] = -264*w - 1320;
            ymap_down[2] = -6*w - 66;
            ymap_down[2] /= -1;
            ymap_down[3] = -1;

        } else if (q == 19)  {
            this->iota_degree = 6; //For some reason degree 6 in bottom y map
            this->starting_a = -152;
            this->starting_b = 722;

            assert (w*w == -19);
            NTL::SetCoeff(xmap_up, 5);
            NTL::SetCoeff(xmap_down, 4);
            NTL::SetCoeff(ymap_up, 6);
            NTL::SetCoeff(ymap_down, 6);
            this->maximal = false;

            xmap_up[0] = 46208*w + 219488;
            xmap_up[0] /= 5;
            xmap_up[1] = -17328*w - 103968;
            xmap_up[1] /= 5;
            xmap_up[2] = 2128*w + 17328;
            xmap_up[2] /= 5;
            xmap_up[3] = -76*w - 1216;
            xmap_up[3] /= 5;
            xmap_up[4] = -2*w + 38;
            xmap_up[4] /= 5;
            xmap_up[5] = w - 9;
            xmap_up[5] /= 50;
            xmap_down[0] = 31768*w + 147288;
            xmap_down[0] /= 25;
            xmap_down[1] = -456*w - 2888;
            xmap_down[1] /= 1;
            xmap_down[2] = 266*w + 2546;
            xmap_down[2] /= 5;
            xmap_down[3] = -2*w - 38;
            xmap_down[3] /= 1;
            xmap_down[4] = 1;
            xmap_down[4] /= 1;
            ymap_up[0] = 877952*w + 3072832;
            ymap_up[0] /= 125;
            ymap_up[1] = -600704*w - 658464;
            ymap_up[1] /= 125;
            ymap_up[2] = 112632*w - 147288;
            ymap_up[2] /= 125;
            ymap_up[3] = -608*w + 54872;
            ymap_up[3] /= 125;
            ymap_up[4] = -1672*w - 5852;
            ymap_up[4] /= 125;
            ymap_up[5] = 27*w + 57;
            ymap_up[5] /= 25;
            ymap_up[6] = -2*w - 7;
            ymap_up[6] /= 125;
            ymap_down[0] = 18875968*w + 38629888;
            ymap_down[0] /= -125;
            ymap_down[1] = -450528*w - 1316928;
            ymap_down[1] /= -5;
            ymap_down[2] = 528504*w + 2174664;
            ymap_down[2] /= -25;
            ymap_down[3] = -2432*w - 14440;
            ymap_down[3] /= -1;
            ymap_down[4] = 684*w + 6384;
            ymap_down[4] /= -5;
            ymap_down[5] = -3*w - 57;
            ymap_down[5] /= -1;
            ymap_down[6] = -1;
        } else {
            throw; //Have to be extremely unlucky
        }
    }

    this->iota_x_num = xmap_up;
    this->iota_x_denom = xmap_down;
    this->iota_y_num = ymap_up;
    this->iota_y_denom = ymap_down;
    // std::cout << "Done with one!" << std::endl;
}

FpE_elem lift(FpE_elem const &alpha, Fp2k const &Fext)
//////////////////////////////////////////////////////////////
///// Takes an element alpha of Fp2, and an extension Fp2k,
///// Returns the of alpha in Fp2k
//////////////////////////////////////////////////////////////

{

    if (alpha == 0 || Fext.k == 1 || NTL::deg(NTL::rep(alpha)) != 1) { //last case is if alpha is over Fp
        return alpha;
    }

    FpE_elem alpha_bar;

    Fp_elem a = NTL::rep(alpha)[0], b = NTL::rep(alpha)[1];
    {
        FpE_push push(Fext.F);
        alpha_bar = a + Fext.Fp2_gen*b;
    }
    return alpha_bar;
}



FpE_elem coerce(FpE_elem const &alpha, Fp2k const &Fext) //super confusing that the first ZZ_pE here refers to Fp2, and the second refers to Fp2k...
//////////////////////////////////////////////////////////////
///// Takes an element alpha of Fp2k, and an extension Fp2k,
///// Returns the of alpha in Fp2k
//////////////////////////////////////////////////////////////

{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return alpha;
    }

    FpE_elem alpha_low;
    {
        auto const idx = Fext.Fp2_gen_nonzero;
        if (NTL::IsZero(NTL::coeff(NTL::rep(alpha), idx))) {
            // If alpha is over Fp, this is needed to avoid representation of the form [x_0, 0]
            alpha_low = alpha;
        } else {
            FpE_push push(Fext.F);
            FpX_elem alpha_poly;
            SetCoeff(alpha_poly, 1);
            alpha_poly[0] = trace(alpha) / (2*Fext.k);
            alpha_poly[1] = NTL::coeff(NTL::rep(alpha), idx) / NTL::rep(Fext.Fp2_gen)[idx];
            conv(alpha_low, alpha_poly);
        }
    }
    //assert (lift(alpha_low, Fext) == alpha);
    return alpha_low;
}

FpEX_elem lift(FpEX_elem const &f, Fp2k const &Fext)
//////////////////////////////////////////////////////////////
///// Takes a polynomial over Fp2, and an extension Fp2k,
///// Returns the polynomial in Fp2k.
//////////////////////////////////////////////////////////////

{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return f;
    }
    if (f == 0){ // 0 is zero
        return f;
    }
    FpE_push push(Fext.F);

    FpEX_elem f_high;
    size_t d = NTL::deg(f);

    NTL::SetCoeff(f_high, d);

    f_high[0] = lift(coeff(f, 0), Fext);
    for (size_t i = 1; i <= d; ++i){
        f_high[i] = lift(coeff(f, i), Fext);
    }

    return f_high;
}


FpEX_elem coerce(FpEX_elem const &f, Fp2k const &Fext)
//////////////////////////////////////////////////////////////
///// Takes a polynomial over Fp2k, and an extension Fp2k,
///// Returns the polynomial in Fp2. Assumes it exists.
//////////////////////////////////////////////////////////////

{
    if (Fext.k == 1){ // Nothing to be done in this case.
        return f;
    }

    FpEX_elem f_low;
    size_t d = NTL::deg(f);

    NTL::SetCoeff(f_low, d);

    f_low[0] = coerce(coeff(f, 0), Fext);
    for (size_t i = 1; i <= d; ++i){
        f_low[i] = coerce(coeff(f, i), Fext);
    }

    return f_low;
}

FpE_elem Fp2k::frob(FpE_elem alpha) const
{
    vec_Fp alpha_vec;
    NTL::conv(alpha_vec, NTL::rep(alpha));
    alpha_vec.SetLength(this->k*2);

    vec_Fp image_alpha_vec = alpha_vec*this->frob_action;

    FpE_elem image_alpha;
    FpX_elem image_alpha_poly;
    NTL::conv(image_alpha_poly, image_alpha_vec);
    NTL::conv(image_alpha, image_alpha_poly);

    return image_alpha;
}

FpEX_elem MinPoly(FpE_elem const &alpha, Fp2k const &Fext)
//////////////////////////////////////////////////////////////
////// Computes the min poly of alpha \in Fext over Fp2
//////////////////////////////////////////////////////////////
{
    if (Fext.k == 1) {
        FpEX_elem f;
        NTL::SetCoeff(f, 1);
        f[0] = -alpha;
        return f;
    }

    FpEX_elem a_min_poly;

    FpEX_elem alpha_over_Fp2;
    NTL::conv(alpha_over_Fp2, NTL::rep(alpha));

    a_min_poly = NTL::MinPolyMod(alpha_over_Fp2 % Fext.mod_Fp2, Fext.mod_Fp2);

    return a_min_poly;
}
