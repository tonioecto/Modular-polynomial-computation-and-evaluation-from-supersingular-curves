///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   This code implements algorithms to for isogeny computations between ellipitc curves 
////     over finite fields
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <cassert>
#include <optional>

#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include "Fp2k.hpp"
#include "ec.hpp"
#include "ecp.hpp"
#include "isog.hpp"


FpEX_elem _derivative(FpEX_elem const &F)
{
    FpEX_elem f;
    if (NTL::IsZero(F)) {
        f = 0;
        return f;
    }
    size_t deg = NTL::deg(F);
    if (NTL::deg(F) == 0) {
        f = 0;
        return f;
    }
    NTL::SetCoeff(f, deg-1);

    for (size_t i = 1; i <= deg; i++) {
        FpE_elem coeff = NTL::coeff(F, i);
        f[i-1] = i*coeff;
    }

    return f;
}

isog isogeny(ecp const &K, int degree) {

    FpEX_elem h_K;
    if (degree > 2) {
        h_K = kernel_polynomial(K, degree);
    } else {
        NTL::SetCoeff(h_K, 1);
        h_K[0] = -coerce(K.aff_x(), K.field());
    }
    if (degree > 2) {
        FpEX_elem psi_der = _derivative(h_K);
        FpEX_elem psi_der_der = _derivative(psi_der);
        FpEX_elem psi_der_der_der = _derivative(psi_der_der);
        int n = (degree - 1)/2;
        assert (n == NTL::deg(h_K));

        FpE_elem s1 = -NTL::coeff(h_K, n-1), s2 = NTL::coeff(h_K, n-2), s3 = -NTL::coeff(h_K, n-3);
        FpE_elem a = K.curve().a(), b = K.curve().b();

        auto s1_sq = s1*s1;
        FpE_elem t = 6*(s1_sq - 2*s2) + n*2*a;
        FpE_elem w = 10*(s1_sq*s1 - 3*s1*s2 + 3*s3) + 6*a*s1 + n*4*b;

        ec E1(a - 5*t, b - 7*w);

        return isog(K.curve_ptr(), std::make_shared<const ec>(E1), degree, h_K, psi_der, psi_der_der, psi_der_der_der);
    } else {
        FpE_elem a = K.curve().a(), b = K.curve().b();

    
        FpE_elem x0 = -coeff(h_K, 0);
        FpE_elem t = 3*x0*x0 + a;
        FpE_elem w = x0*t;

        ec E1(a - 5*t, b - 7*w);
        return isog(K.curve_ptr(), std::make_shared<const ec>(E1), degree, h_K, t);
    }
}

std::pair<FpE_elem, FpE_elem> _prod_derivative_pair(std::pair<FpE_elem, FpE_elem> const &f1, std::pair<FpE_elem, FpE_elem> const &f2) {
    std::pair<FpE_elem, FpE_elem> out(f1.first*f2.first, f1.first*f2.second + f1.second*f2.first);
    return out;
}

ecp isog::operator()(ecp const &P) const
{
    assert (P.curve() == this->get_domain());

    if (P.is_identity()) {
        return this->get_codomain()(P.field());
    }

    if (this->_degree == 2) {
        return this->_even_evaluation(P);
    }

    FpE_push push(P.field().F);

    FpE_elem x = P.aff_x();

    FpEX_elem psi = lift(this->psi, P.field());
    FpEX_elem psi_der = lift(this->psi_der, P.field());
    FpEX_elem psi_der_der = lift(this->psi_der_der, P.field());
    FpEX_elem psi_der_der_der = lift(this->psi_der_der_der, P.field());

    FpE_elem psi_x = NTL::coeff(psi, 0);
    FpE_elem psi_der_x = NTL::coeff(psi_der, 0);
    FpE_elem psi_der_der_x = NTL::coeff(psi_der_der, 0);
    FpE_elem psi_der_der_der_x = NTL::coeff(psi_der_der_der, 0);

    size_t n = NTL::deg(psi);
    FpE_elem xi(1);

    
    for (size_t i = 1; i <= n; i++) { //NTL::coeff returns zero when i > degree
        xi *= x;
        psi_x += xi*NTL::coeff(psi, i);
        psi_der_x += xi*NTL::coeff(psi_der, i);
        psi_der_der_x += xi*NTL::coeff(psi_der_der, i);
        psi_der_der_der_x += xi*NTL::coeff(psi_der_der_der, i);
    }

    // Evaluating the derivative without expanding the polynomial Pogchamp
    std::pair<FpE_elem, FpE_elem> psi_x_pair(psi_x, psi_der_x);
    std::pair<FpE_elem, FpE_elem> psi_der_x_pair(psi_der_x, psi_der_der_x);
    std::pair<FpE_elem, FpE_elem> psi_der_der_x_pair(psi_der_der_x, psi_der_der_der_x);

    FpE_elem b4 = lift(this->get_domain().a(), P.field())*2, b6 = lift(this->get_domain().b(), P.field())*4;
    std::pair<FpE_elem, FpE_elem> aux_1_pair((4*x*x*x + 2*b4*x + b6), (12*x*x + 2*b4));
    std::pair<FpE_elem, FpE_elem> aux_2_pair((6*x*x + b4), (12*x));
    std::pair<FpE_elem, FpE_elem> aux_3_pair((this->_degree*x + 2*NTL::coeff(psi, n-1)), this->_degree);

    //begin the evaluation
    auto psi_x_sq_pair = _prod_derivative_pair(psi_x_pair, psi_x_pair);
    auto psi_der_x_sq_pair = _prod_derivative_pair(psi_der_x_pair, psi_der_x_pair);
    auto psi_x_psi_der_der_x_pair = _prod_derivative_pair(psi_x_pair, psi_der_der_x_pair);

    std::pair<FpE_elem, FpE_elem> temp(psi_der_x_sq_pair.first - psi_x_psi_der_der_x_pair.first, psi_der_x_sq_pair.second - psi_x_psi_der_der_x_pair.second);
    auto summand_1 = _prod_derivative_pair(aux_1_pair, temp);

    auto psi_x_psi_der_x_pair = _prod_derivative_pair(psi_der_x_pair, psi_x_pair);
    auto summand_2 = _prod_derivative_pair(psi_x_psi_der_x_pair, aux_2_pair);

    auto summand_3 = _prod_derivative_pair(psi_x_sq_pair, aux_3_pair);
    
    std::pair<FpE_elem, FpE_elem> N_x_pair(summand_1.first - summand_2.first + summand_3.first, summand_1.second - summand_2.second + summand_3.second);

    FpE_elem one(1);
    std::pair<FpE_elem, FpE_elem> psi_x_sq_inv_pair(one/psi_x_sq_pair.first, -psi_x_sq_pair.second/(psi_x_sq_pair.first*psi_x_sq_pair.first));
    auto eval_res = _prod_derivative_pair(N_x_pair, psi_x_sq_inv_pair);

    auto new_x = eval_res.first;
    auto new_y = P.aff_y()*eval_res.second;

    auto phi_P = this->get_codomain()(P.field(), new_x, new_y);

    return phi_P;
}

ecp isog::_even_evaluation(ecp const &P) const
{
    FpE_push push(P.field().F);
    FpE_elem x0 = lift(-NTL::coeff(this->psi, 0), P.field()), x = P.aff_x(), y = P.aff_y();
    FpE_elem tlift = lift(t, P.field());
    FpE_elem xmx0 = x-x0;
    FpE_elem new_x = x + tlift/xmx0;
    FpE_elem new_y = y*(1 - tlift/(xmx0*xmx0));
 
    ecp phi_P(this->get_codomain_ptr(), P.field(), new_x, new_y);
    return phi_P;
}

int _least_primitive(unsigned l) {
    NTL::zz_pPush push(l);
    
    int a0 = 2;
    NTL::zz_p a;

    auto fac_order = factor(NTL::ZZ(l-1));
    long radical_order(1);
    for (auto ell_e : fac_order) {
        radical_order *= ell_e.first;
    }

    auto cofac = (l-1)/radical_order;

    for (size_t i = 0; i < l ; i++) {
        bool gen = true;
        a = NTL::power(NTL::zz_p(a0), cofac);
        for (auto ell_e : fac_order) {
            if (NTL::power(a, radical_order/ell_e.first) == 1) {
                gen = false;
                break;
            }
        }
        if (gen) {
            return a0;
        };
        a0 += 1;
    }
    throw;
}

FpEX_elem kernel_polynomial(ecp const&K, int degree) {
    ///////////////////////////////////////
    /////// K point of order degree, 
    /////// Returns the kernel poly of <K>
    ///////////////////////////////////////

    FpEX_elem Psi;

    FpE_elem x = K.aff_x();
    if (degree == 3) {
        NTL::SetCoeff(Psi, 1);
        Psi[0] = -coerce(x, K.field());  
        return Psi; 
    }

    FpEX_elem f0 = MinPoly(x, K.field());
    int k = NTL::deg(f0);

    size_t m = (degree-1)/(2*k);

    std::vector<FpEX_elem> fs;
    fs.reserve(m);
    
    fs.push_back(f0);

    int a = _least_primitive(degree);
    ecp Ki = K;

    for (size_t i = 1; i < m; i++) {
        Ki = a*Ki;
        fs.push_back(MinPoly(Ki.aff_x(), K.field()));
    }

    Psi = product_tree<FpEX_elem>(fs);

    return Psi;
}

isog_chain::isog_chain(std::vector<std::pair<ecp,std::pair<int, int>>> kerGens) //Somehow this horrible structure feels like bad practice
{
    auto start_id2ker = std::chrono::steady_clock::now();
    std::chrono::duration<long, std::milli> duration_compisog = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());
    std::chrono::duration<long, std::milli> duration_pushpoints = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()-std::chrono::steady_clock::now());
    
    while (kerGens.size() > 0) {

        std::pair<ecp, std::pair<int, int>> Ker = kerGens.back();
        kerGens.pop_back();

        ecp K = Ker.first; std::pair<int, int> deg = Ker.second;
        int ell = deg.first; int e = deg.second;
                
        //std::cout << "ISOG CHAIN: Computing isogeny over Fp2k, k = " << K.field().k << "of degree ell = " << ell << std::endl;
        auto start_isog = std::chrono::steady_clock::now();
        isog_parts.push_back(isogeny(NTL::power(NTL::ZZ(ell), e-1)*K, ell));
        duration_compisog = std::chrono::duration_cast<std::chrono::milliseconds>(duration_compisog + std::chrono::steady_clock::now() - start_isog);

        //std::cout << "ISOG CHAIN: Done! Now pushing points!" << std::endl;
        auto start_push = std::chrono::steady_clock::now();
        std::vector<std::pair<ecp,std::pair<int, int>>> newKerGens;
        for (std::pair<ecp, std::pair<int, int>> KerPart : kerGens) {
            ecp Kpart = KerPart.first; std::pair<int, int> deg = KerPart.second;
            //std::cout << "ISOG CHAIN: deg = " << deg.first << "^" << deg.second << std::endl;
            newKerGens.push_back(std::pair<ecp, std::pair<int, int>>(isog_parts.back()(Kpart), deg));
        }

        if (e > 1) {
            //std::cout << "ISOG CHAIN: special, point itself: deg = " << deg.first << "^" << deg.second << std::endl;
            std::pair<int, int> deg(ell, e-1);
            newKerGens.push_back(std::pair<ecp, std::pair<int, int>>(isog_parts.back()(K), deg));
        }

        duration_pushpoints = std::chrono::duration_cast<std::chrono::milliseconds>(duration_pushpoints + std::chrono::steady_clock::now() - start_push);
        //std::cout << "Done!" << std::endl;
        kerGens = newKerGens;
    }
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_id2ker);
    // std::cout << ">>> kernel to isog took: " << duration.count() << " milliseconds" << std::endl;
    // std::cout << ">>> On computing isogs: " << duration_compisog.count() << " milliseconds" << std::endl;
    // std::cout << ">>> On pushing points: " << duration_pushpoints.count() << " milliseconds" << std::endl;
}
