///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implements the algorithm OrdersTojInvariantBigSet from:
////
////        Antonin Leroux. Computation of Hibert Class polynomials and modular
////        polynomials from supersingular elliptic curves.
////        https://eprint.iacr.org/2023/064
////
////   (c.f. Algorithm 1)
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include "ordertojinvbigset.hpp"

#include <algorithm>
#include <unordered_set>
#include <ctime>



void print_faclist(std::unordered_map<int, int> const &m)
{
    for (auto it = m.cbegin(); it != m.cend(); ++it) {
        std::cout << "{" << (*it).first << ": " << (*it).second << "}\n";
    }
}


//
std::map<NTL::ZZ,std::pair<ecp,ecp>> BasesPreProcessing(quat const &alpha, quat const &beta, std::unordered_map<int, int>& facN, ec const &E0, const std::map<unsigned,Fp2k> &FieldExtensions) {
    // I = O0<alpha, N>
    // Field extensions should maybe be something like an unordered map instead?

    std::map<NTL::ZZ,std::pair<ecp,ecp>> result_map;
    std::map<NTL::ZZ,std::pair<ecp,ecp>> TorsionBases;

    for (const auto& [ell,e] : facN) {
        auto ell_e_m_1 = NTL::power(NTL::ZZ(ell), e-1);
        auto ell_e = ell_e_m_1 * ell;
        auto ell_e_p_1 = ell_e * ell;
        assert((2*alpha.alg.q) % alpha[4] == 0);
        bool extra1 = (alpha[4] % 2 == 0 && ell == 2) || (alpha[4] % alpha.alg.q == 0 && ell == alpha.alg.q);
        bool extra2 = (beta[4] % 2 == 0 && ell == 2) ||  (beta[4] % alpha.alg.q == 0 && ell == alpha.alg.q);
        auto ell_ext = ell_e * (extra1 ? (ell) : 1) * (extra2 ? (ell) : 1);
        quat endo = alpha.conjugate() * alpha[4];

        auto it = TorsionBases.find(ell_ext);
        if (it == TorsionBases.end()) {
            unsigned k = torsionToFieldDegree(ell_ext);
            auto jt = FieldExtensions.find(k);
            assert(jt != FieldExtensions.end());  //TODO
            assert(jt->first == jt->second.k);
            auto bas = E0.torsionBasisDet(jt->second, ell, e + extra1 + extra2);
            it = TorsionBases.emplace(std::make_pair(ell_ext, bas)).first;
        }
        auto const &Basis = it->second;
        assert(!(ell_ext * Basis.first) && ell_ext/ell * Basis.first);
        assert(!(ell_ext * Basis.second) && ell_ext/ell * Basis.second);

        auto const project = [&](ecp const &pt) {
            return evalEndo(endo, pt, ell_ext);
        };

        ecp K = project(Basis.first);
        assert(!(ell_e * (extra2 ? (ell) : 1) * K));
        if (!(ell_e_m_1 * (extra2 ? (ell) : 1) * K))
            K = project(Basis.second);
        quat sec_endo = beta;

        sec_endo[4] = 1;
        Integer inv_elem = beta[4]/(extra2 ? (ell) : 1);
        for (int i = 0; i < 4; i++) {
            sec_endo[i] = sec_endo[i] * NTL::InvMod(inv_elem % ell_e, ell_e);
        }

        ecp K2 = evalEndo(sec_endo,K,ell_ext/(extra1 ? (ell) : 1));

        K = (extra2 ? (ell) : 1) * K;

        assert(!(ell_e * K));
        assert((ell_e_m_1) * K);
        assert(!(ell_e * K2));
        assert((ell_e_m_1) * K2);

        assert(!((K-K2).is_identity()));
        assert(!((ell_e_m_1 * K- ell_e_m_1 * K2).is_identity()));
        result_map.emplace(std::make_pair(ell_e, std::make_pair(K,K2)));
    }

    return result_map;
}



Jinv JToJinv(const FpE_elem &j) {
        Jinv J;
        NTL::ZZ IntList[2];
        IntList[0] = rep( NTL::coeff(rep(j),0));
        IntList[1] = rep( NTL::coeff(rep(j),1));
        NTL::BytesFromZZ(J.IntList[0], IntList[0], LenNumberBytes);
        NTL::BytesFromZZ(J.IntList[1], IntList[1], LenNumberBytes);
        return J;
}


bool is_Fp(const Fp2 &j) {
    return NTL::coeff(rep(j),1) == 0;
}

Fp2 Frob(const Fp2 &j) {
    // std::cout << "\n p= " << NTL::coeff(rep(j),1).modulus() << "\n";
    Fp2 jp;
    FpX jpp;
    // std::cout << jp[0];
    Fp temp;
    // std::cout << "a " << NTL::coeff(rep(j),1)<< " \n";
    NTL::negate(temp, NTL::coeff(rep(j),1));
    // std::cout << "b " << temp << " \n";
    NTL::SetCoeff(jpp, 1, temp);
    NTL::SetCoeff(jpp, 0, NTL::coeff(rep(j),0));
    NTL::conv(jp, jpp);
    //assert( !ProbPrime(rep(j)[0].modulus()) || jp == power(j, rep(j)[0].modulus()));
    return jp;
}

void print_key(const Key& k) {
    NTL::ZZ IntList[3];
    for (int i=0; i<3 ; i++ ) {
            NTL::ZZFromBytes(IntList[i],k.IntList[i],LenNumberBytes);
    }
    std::cout << "(" << IntList[0] << "," << IntList[1] << "," << IntList[2] << ")";
}

void order_to_jinv_full_list(std::unordered_map<Key, std::pair<FpE_elem, quatlat>, KeyHash, KeyEqual> &m, std::vector<std::pair<quatlat,std::pair<quatlat, Key>>> &id_list, const  NTL::ZZ &p, const quatalg &Bp, const std::map<unsigned,Fp2k> &Fexts, const Integer &coprime) {

    NTL::ZZ q = Bp.q;

    // we start by creating the map and the quaternion algebra

    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;

    FpE_elem j0 = start.first.j_invariant();

    quat dum = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};

    Key K0 = order_invariant_computation(O0, &dum);
    m.insert({K0, {j0, O0}});
    id_list.push_back( {O0, {O0, K0}} );

    // torsion selected
    factor_list tors_list = choose_torsion(p, 2*p/3, coprime);

    NTL::ZZ S(1);
    std::list<int> prime_list;
    std::unordered_map<int,int> fac_list;
    unsigned k_bound = Fexts.size(); //this is the number of extensions that have been computed.
    for (const auto& tup : tors_list) {
        int fac = std::get<0>(tup);
        S *= fac; // ell occurs once for each exponent....
        // ;
        auto is_in = std::find(prime_list.begin(), prime_list.end(), fac);
        if (is_in == prime_list.end()){
            prime_list.push_back(fac);
            fac_list[fac] = 1;
        }
        else {
            fac_list[fac] = fac_list[fac]+1;
        }
        if (std::get<2>(tup) > (ssize_t) k_bound) {
            // k_bound = std::get<2>(tup);
            // we need more extensions that have been computed !!
            assert(0);
        }
    }
    prime_list.sort();
    prime_list.reverse();





    // now we generate one ideal of order S
    // This is a bit empirical
    bool found = false;
    quat gamma = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, Bp};
    while(!found) {
        NTL::ZZ cofac;
        if (p != 3067) {
            cofac = 3067;
        }
        else {
            cofac = 3319;
        }
        NTL::ZZ target = cofac * S;
        while (target < 1000* p) {
            target = cofac * target;
        }
        quat gamma_tmp = DetRepresentInteger(Bp, target);
        quat testing = {{gamma_tmp[0], gamma_tmp[1], gamma_tmp[2], gamma_tmp[3], NTL::ZZ(2)*gamma_tmp[4]}, Bp};
        found = !O0.contains(testing);
        gamma[0] = gamma_tmp[0];
        gamma[1] = gamma_tmp[1];
        gamma[2] = gamma_tmp[2];
        gamma[3] = gamma_tmp[3];
        gamma[4] = gamma_tmp[4];
    }

    assert(O0.contains(gamma));

    quatlat I = O0 * gamma + O0 * S;
    assert(std::get<0>(I.norm())/std::get<1>(I.norm()) == S);
    // std::cout << S << "\n";
    // print_faclist(fac_list);



    // generating the iterating quaternion
    quat beta = find_quaternion_iterator(prime_list, I , O0, Bp);
    quat betabar = beta.conjugate();

    // now we generate the bases
    std::map<NTL::ZZ,std::pair<ecp,ecp>> TorsionBases = BasesPreProcessing(gamma, beta, fac_list, start.first, Fexts);

    std::list<std::pair<ecp,ecp>> torsion_list = {};

    for (auto l : prime_list) {
        NTL::ZZ ell(l);
        NTL::ZZ elle = power(ell,fac_list[l]);
        auto it = TorsionBases.find(elle);
        torsion_list.push_back(it->second);
        assert(!(elle*it->second.first));
        assert(!(elle*it->second.second));
        assert(!((it->second.second-it->second.first).is_identity()));
    }

    std::list<std::tuple<quatlat,ec,std::list<std::pair<ecp,ecp>>>> global_list = {{O0,start.first,torsion_list}};
    clock_t enum_time = clock(); (void) enum_time;

    NTL::ZZ target_num;
    if (p%12 == 1) {
        target_num = p/12;
    }
    else if (p%12 ==7  || p%12 == 5) {
        target_num = p/12 + 1;
    }
    else {
        assert(p%12 == 11);
        target_num = p/12 + 2;
    }

    std::unordered_map<Jinv, std::pair<quatlat,Key> ,JinvHash, JinvEqual> j_inv_list = {};

    // iterating through the set of S-ideals, we go factor by factor
    int order_count = 1;
    int count_mistake = 0;
    // std::cout << "target = " << target_num <<"\n";

    clock_t quaternion_time = 0;
    clock_t rigo_time = 0;
    clock_t isogeny_time = 0;
    clock_t reduc_time = 0;
    clock_t t;

    int isog_count = 0;
    int quat_count = 0;

    for (auto l : prime_list) {
        // std::cout << "round for" << l << " \n";
        std::list<std::tuple<quatlat,ec,std::list<std::pair<ecp,ecp>>>>  new_list = {};
        std::list<std::tuple<NTL::ZZ,NTL::ZZ>> coeff_list = {};
        // first we create the list of coefficients
        NTL::ZZ ell(l);
        NTL::ZZ ell_e = power(ell,fac_list[l]);
        NTL::ZZ ell_e_m_1 = power(ell,fac_list[l]-1);
        NTL::ZZ iterate = NTL::ZZ(0);
        while (iterate < ell_e) {
            coeff_list.push_back({NTL::ZZ(1),iterate});
            if (NTL::GCD(iterate,NTL::ZZ(ell))!=1) {
                coeff_list.push_back({iterate,NTL::ZZ(1)});
            }
            iterate++;
        }
        int count = 0;
        std::unordered_map<Key, FpE_elem, KeyHash, KeyEqual> local_map;
        for (auto tup : global_list) {
            quatlat K = std::get<0>(tup);
            ec E = std::get<1>(tup);
            std::list<std::pair<ecp,ecp>> TorBas = std::get<2>(tup);
            std::pair<ecp,ecp> Basis = TorBas.front();
            ecp P = Basis.first;
            ecp Q = Basis.second;
            assert(!((P-Q).is_identity()));
            assert(!(ell_e*P));
            assert(!(ell_e*Q));
            TorBas.pop_front();
            for (auto coeff : coeff_list) {
                // std::cout << std::get<0>(coeff) << " : " << std::get<1>(coeff) << "\n";
                count++;
                t = clock();
                quat gen = (gamma * (quat({{std::get<0>(coeff),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(1)},Bp}) + betabar * std::get<1>(coeff)  ));
                quatlat J = create_from_generator_O0( gen, ell_e);
                // quatlat J = O0 * gen + O0 * ell_e; :
                J._intersect(K);
                // std::cout << J.norm().first << " " << J.norm().second << "\n";
                // std::cout << K.norm().first << " " << K.norm().second << "\n";
                quaternion_time += clock() - t;
                t = clock();
                quatlat rigo = J.right_order();
                rigo_time += clock() - t;
                t = clock();
                Key K = order_invariant_computation(rigo, &dum);
                reduc_time += clock() - t;
                quat_count++;
                auto a_loc = local_map.insert({K,j0});
                if (a_loc.second) {
                    t = clock();
                    isog_count++;
                    ecp ker = std::get<0>(coeff) * P + std::get<1>(coeff) * Q;
                    assert((ell_e_m_1 * ker));
                    assert(!(ell_e * ker));
                    // now to the isogeny computation
                    std::vector<std::pair<ecp,std::pair<int, int>>> kerGen = {{ker, {l,fac_list[l]}} };
                    // computing the isogeny
                    isog_chain phi = isog_chain(kerGen);
                    // pushing the other bases
                    std::list<std::pair<ecp,ecp>> new_torsion_list = {};
                    for (auto bas : TorBas) {
                        new_torsion_list.push_back( {phi(bas.first), phi(bas.second)});
                    }
                    FpE_elem j = phi.get_codomain().j_invariant();
                    isogeny_time += clock() - t;
                    new_list.push_back({J,phi.get_codomain(),new_torsion_list});

                    auto a = m.insert({K, {j, J}});
                    if (a.second) {
                        // check if the element is not already contained in the list
                        // TODO this is only used for debug purpose
                        Jinv jtemp = JToJinv(j);
                        auto search = j_inv_list.find(jtemp);
                        if (search != j_inv_list.end()) {
                            count_mistake ++;
                            quatlat test = J.conjugate() * search->second.first;
                            isPrincipal(test);
                            print_key(order_invariant_computation(rigo,&dum));
                            print_key(K);
                            print_key(search->second.second);
                            print_key(order_invariant_computation(search->second.first.right_order(),&dum));
                            std::cout << "\n";



                        }

                        order_count ++;
                        Jinv jnew = JToJinv(j);
                        j_inv_list.insert({jnew,{J,K}});
                        id_list.push_back( {J, {rigo, K}} );
                        // we add one more if the curve is not defined over Fp
                        if (!is_Fp(j)) {
                            order_count ++;
                            jnew = JToJinv(Frob(j));
                            quat jj = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(1)}, Bp};
                            J._intersect(O0*jj);
                            j_inv_list.insert({jnew,{J,K}});
                        }
                    }
                }
                if (order_count >= target_num) {
                    goto endloop;
                }


            }

        }
        // std::cout << "number of elements in the new list " << count << " temp count = " << order_count << " \n";
        global_list = new_list;

    }
    std::cerr << "exhausted the list and did not find all the curves" << std::endl;
    assert(0);

endloop:
    assert(count_mistake == 0);
    std::cout << "number of curves :" << order_count << "\n";
    std::cout << "size of the list of curves : " << m.size() << "\n";
    std::cout << "number of mistakes :" << count_mistake << "\n";
    // std::cout << "time to compute  : " << (double) (clock() - enum_time)/CLOCKS_PER_SEC << "\n";
    // std::cout << "time for quaternion operations (amortized per 100) : " << (double) (100 * quaternion_time)/(CLOCKS_PER_SEC * quat_count) << "\n";
    // std::cout << "time for right-ordr operations (amortized per 100) : " << (double) (100 * rigo_time)/(CLOCKS_PER_SEC * quat_count) << "\n";
    // std::cout << "time for reduction operations  (amortized per 100) : " << (double) (100 * reduc_time)/(CLOCKS_PER_SEC * quat_count) << "\n";
    // std::cout << "time for isogeny operations    (amortized per 100) : " << (double) (100 * isogeny_time)/(CLOCKS_PER_SEC * isog_count) << "\n";
    // std::cout << "total time for quaternion operations : " << (double) (quaternion_time)/(CLOCKS_PER_SEC) << "\n";
    // std::cout << "total time for right-ordr operations : " << (double) (rigo_time)/(CLOCKS_PER_SEC) << "\n";
    // std::cout << "total time for reduction operations  : " << (double) (reduc_time)/(CLOCKS_PER_SEC) << "\n";
    // std::cout << "total time for isogeny operations    : " << (double) (isogeny_time)/(CLOCKS_PER_SEC) << "\n";

}



std::pair<weber_bas,std::vector<std::pair<VerySmallMat,VerySmallMat>>> fast_order_to_weber_inv_full_list(std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<FastInteger, std::pair<std::pair<FastQuat,FastQuat>,std::pair<FastInteger,FastQuat>>>, weber_full_data>>, KeyHash, KeyEqual> &m, std::vector<std::pair<FastQuatLat, Key>> &id_list, const Integer &p, const quatalg &Bp, const FastQuatAlg &fast_Bp, const std::map<unsigned,Fp2k> &Fexts, const Integer &coprime) {

    // we restrict to p = 3 mod 4
    assert(p % 4 == 3);    

    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;

    // FastQuat init 
    FastQuatLat fast_O0 = FastQuatLat(O0, fast_Bp);

    
    Fp2 j0 = Fp2(1728);
    Key K0;
    FastIntToBytes(4, K0.IntList[0],LenNumberBytes);
    FastIntToBytes(convert(p), K0.IntList[1],LenNumberBytes);
    FastIntToBytes(convert(p + 1),K0.IntList[2],LenNumberBytes);

    // torsion selected
    factor_list tors_list = choose_torsion(p, 10 * p, 6 * coprime);

    NTL::ZZ S(1);
    std::list<int> prime_list;
    std::unordered_map<int,int> fac_list;
    unsigned k_bound = Fexts.size(); //Some minimum
    for (const auto& tup : tors_list) {
        int fac = std::get<0>(tup);
        S *= fac; // ell occurs once for each exponent....
        // ;
        auto is_in = std::find(prime_list.begin(), prime_list.end(), fac);
        if (is_in == prime_list.end()){
            prime_list.push_back(fac);
            fac_list[fac] = 1;
        }
        else {
            fac_list[fac] = fac_list[fac]+1;
        }
        if (std::get<2>(tup) > (ssize_t) k_bound) {
            k_bound = std::get<2>(tup);
        }
    }
    prime_list.sort();
    prime_list.reverse();

    // now we generate one ideal of order S
    // This is a bit empirical
    bool found = false;
    quat gamma = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)}, Bp};
    NTL::ZZ cofac;
    {  
        if (p != 2099) {
            cofac = 2099;
        }
        else {
            cofac = 977;
        }
        NTL::ZZ target = cofac * S;
        while (target < 10000 * p) {
            target = cofac * target;
        }
        while(!found) {
            // std::cout << "target = " << target << "\n";
            quat gamma_tmp = DetRepresentInteger(Bp, target);
            quat testing = {{gamma_tmp[0], gamma_tmp[1], gamma_tmp[2], gamma_tmp[3], NTL::ZZ(2)*gamma_tmp[4]}, Bp};
            found = !O0.contains(testing);
            gamma[0] = gamma_tmp[0];
            gamma[1] = gamma_tmp[1];
            gamma[2] = gamma_tmp[2];
            gamma[3] = gamma_tmp[3];
            gamma[4] = gamma_tmp[4];
            target = cofac * target;
        }
    }

    assert(O0.contains(gamma));

    quatlat I = O0 * gamma + O0 * S;
    assert(std::get<0>(I.norm())/std::get<1>(I.norm()) == S);


    // generating the iterating quaternion
    quat beta = find_quaternion_iterator(prime_list, I , O0, Bp);
    quat betabar = beta.conjugate();

    FastQuat fast_betabar = FastQuat(betabar, fast_Bp);
    FastQuat fast_gamma = FastQuat(gamma, fast_Bp);

    // now we generate the bases
    // for the list of torsion points we are going to need
    std::map<Integer, std::pair<ecp, ecp>> TorsionBases = BasesPreProcessing(gamma, beta, fac_list, start.first, Fexts);
    std::list<std::pair<ecp, ecp>> torsion_list = {};
    for (auto l : prime_list) {
        NTL::ZZ ell(l);
        NTL::ZZ elle = power(ell,fac_list[l]);
        auto it = TorsionBases.find(elle);
        torsion_list.push_back(it->second);
        assert(!(elle*it->second.first));
        assert(!(elle*it->second.second));
        assert(!((it->second.second-it->second.first).is_identity()));
    }

    // now we generate 3 and 16 torsion bases for Weber application
    std::pair<ecp,ecp> bas2 = torsion_list.front();
    std::pair<ecp,ecp> bas3 = torsion_list.front();
    bool extra3 = Bp.q == 3;
    {
        unsigned k = torsionToFieldDegree(Integer(3) * (extra3 ? (Integer(3)): Integer(1)));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());
        assert(jt->first == jt->second.k);
        bas3 = start.first.torsionBasisDet(jt->second, int(3), 1 + (extra3 ? (1): 0));
    }
    {
        unsigned k = torsionToFieldDegree(Integer(32));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());
        assert(jt->first == jt->second.k);
        bas2 = start.first.torsionBasisDet(jt->second, int(2), 5);
    }
    weber_bas bas0 = {bas2.first, bas2.second, bas3.first, bas3.second};
    bas2.first = 2 * bas2.first;
    bas2.second = 2 * bas2.second;
    bas3.first =  (extra3 ? (Integer(3)) : Integer(1)) * bas3.first;
    bas3.second =  (extra3 ? (Integer(3)) : Integer(1)) * bas3.second;
    torsion_list.push_back(bas3);
    torsion_list.push_back(bas2);
    assert((16*bas2.first).is_identity());
    assert((3*bas3.first).is_identity());


    // precomputing the action of the endo ring on the weber_basis
    std::vector<std::pair<VerySmallMat,VerySmallMat>> mat0 = {};
    bool extra2 = (O0.denom % 2 == 0);

    for (int j = 0; j<4; j++) {
        Integer inv_elem = O0.denom/(extra2 ? (Integer(2)): Integer(1));
        quat tmp = { { O0.basis[j][0], O0.basis[j][1], O0.basis[j][2],O0.basis[j][3], Integer(1) }, O0.alg };
        assert(inv_elem%2 != 0);
        for (int i = 0; i < 4; i++) {
                tmp[i] = tmp[i] * NTL::InvMod(inv_elem % 16, 16);
        }
        ecp gammaP016 = evalEndo( tmp, (!extra2 ? (Integer(2)): Integer(1)) * bas0.P16,  Integer(32) / (!extra2 ? (Integer(2)): Integer(1)));
        ecp gammaQ016 = evalEndo( tmp, (!extra2 ? (Integer(2)): Integer(1)) * bas0.Q16,  Integer(32) / (!extra2 ? (Integer(2)): Integer(1)));
        assert((16*gammaP016).is_identity());
        assert((16*gammaQ016).is_identity());

        // now we do the same for the order 3 part.
        inv_elem = O0.denom/(extra3 ? (Integer(3)): Integer(1));
        assert(inv_elem%3 != 0);
        for (int i = 0; i < 4; i++) {
                tmp[i] = O0.basis[j][i] * NTL::InvMod(inv_elem % 3, 3);
        }
        tmp[4] = Integer(1);
        ecp gammaP03 = evalEndo( tmp, bas0.P3,  Integer(9) / (!extra3 ? (Integer(3)): Integer(1)));
        ecp gammaQ03 = evalEndo( tmp, bas0.Q3,  Integer(9) / (!extra3 ? (Integer(3)): Integer(1)));
        assert((3*gammaP03).is_identity());
        assert((3*gammaQ03).is_identity());

        // now we compute the matrices of order 16 and 3
        auto temp_M16 = change_of_basis16( {2*bas0.P16, 2*bas0.Q16}, {gammaP016, gammaQ016} );
        auto temp_M3 = change_of_basis3( {(extra3 ? (Integer(3)): Integer(1)) * bas0.P3, (extra3 ? (Integer(3)): Integer(1)) * bas0.Q3}, {gammaP03, gammaQ03} );
       
        

        VerySmallMat M16, M3;
        M16[0][0] = (int) temp_M16.mat[0][0];
        M16[0][1] = (int) temp_M16.mat[0][1];
        M16[1][0] = (int) temp_M16.mat[1][0];
        M16[1][1] = (int) temp_M16.mat[1][1];

        // std::cout << temp_M16 << " " << (int) M16[0][0] << " " << (int) M16[1][0] << " " << (int) M16[0][1] << " " << (int) M16[1][1] << "\n";

        M3[0][0] = (int) temp_M3.mat[0][0];
        M3[0][1] = (int) temp_M3.mat[0][1];
        M3[1][0] = (int) temp_M3.mat[1][0];
        M3[1][1] = (int) temp_M3.mat[1][1];

        mat0.push_back({M16,M3});

    }
    weber_bas w0 = { bas2.first, bas2.second, bas3.first, bas3.second};
    w0.P3.normalize();
    w0.Q3.normalize();
    w0.P16.normalize();
    w0.Q16.normalize();

    // initiliazing the list of element for the enumeration 
    std::list<std::tuple<FastQuatLat,ec,std::list<std::pair<ecp,ecp>>>> global_list = {{fast_O0, start.first, torsion_list}};

    // computing the number of supersingular curves
    Integer target_num;
    if (p%12 == 1) {
        target_num = p/12;
    }
    else if (p%12 ==7  || p%12 == 5) {
        target_num = p/12 + 1;
    }
    else {
        assert(p%12 == 11);
        target_num = p/12 + 2;
    }
    
    // some precomputation for fast enumeration of all the weber invariants
    weber_enum_poly_precomp precomp = SetWeberPrecomp();
    fast_weber_enum_poly_precomp fast_precomp = Set(&precomp);


    weber_full_data w0_data = EnumerateAllWeberFastFast(w0, Fexts, &fast_precomp);
    FastQuat ii = { {0, 1, 0, 0, 1}, fast_O0.alg};
    FastQuat kk = { {0, 0, 0, 1, 1}, fast_O0.alg};
    std::pair<std::pair<FastInteger,std::pair<std::pair<FastQuat, FastQuat>, std::pair<FastInteger,FastQuat>>>, weber_full_data> bos = { {1, {{ii * 2, kk},{1, ii * kk}}}, w0_data };
    m.insert({K0, {j0, bos }});

    std::unordered_map<Jinv, std::pair<FastQuatLat,Key> ,JinvHash, JinvEqual> j_inv_list = {};

    // iterating through the set of S-ideals, we go factor by factor
    int order_count = 1;

    clock_t quaternion_time = 0;
    clock_t isogeny_time = 0;
    clock_t weber_time = 0;
    clock_t t;

    int isog_count = 0;
    int quat_count = 0;

    FastInteger Prod = 1;
    int num_prime = 0;

    for (auto l : prime_list) {

        std::list<std::tuple<FastQuatLat,ec,std::list<std::pair<ecp,ecp>>>>  new_list = {};
        std::list<IntegerPair> coeff_list = {};
        
        // first we create the list of coefficients
        Integer ell(l);
        Integer ell_e = power(ell,fac_list[l]);
        Integer ell_e_m_1 = power(ell,fac_list[l]-1);
        Integer iterate = NTL::ZZ(0);
        while (iterate < ell_e) {
            coeff_list.push_back({NTL::ZZ(1),iterate});
            if (NTL::GCD(iterate,NTL::ZZ(ell))!=1) {
                coeff_list.push_back({iterate,NTL::ZZ(1)});
            }
            iterate++;
        }
        int count = 0;
        std::unordered_map<Key, FpE_elem, KeyHash, KeyEqual> local_map;

        // precomputing all relevant information for fast inversion and so on
        SignedBarrettReducer red_elle = SignedBarrettReducer(convert(ell_e));
        SignedBarrettReducer red_elle_S = SignedBarrettReducer(convert(ell_e) * Prod);
        SignedBarrettReducer red_S = SignedBarrettReducer(Prod);

        FastInteger inv = InvMod2(convert(ell_e), red_S);
        FastInteger inv2sqr = InvMod2Sqr(convert(ell_e * ell_e), red_S);

        Prod = Prod * convert(ell_e);

        for (auto tup : global_list) {

            FastQuatLat K_id = std::get<0>(tup);
            
            ec E = std::get<1>(tup);
            std::list<std::pair<ecp,ecp>> TorBas = std::get<2>(tup);
            std::pair<ecp,ecp> Basis = TorBas.front();
            ecp P = Basis.first;
            ecp Q = Basis.second;
            assert(!((P-Q).is_identity()));
            assert(!(ell_e*P));
            assert(!(ell_e*Q));
            TorBas.pop_front();

            for (auto coeff : coeff_list) {

                count++;
                t = clock();
                FastQuat gen = (fast_gamma * (FastQuat({{ convert(std::get<0>(coeff)), 0, 0, 0, 1}, fast_Bp}) + fast_betabar * convert(std::get<1>(coeff))  ));
                // TODO : we redo this computation several time for nothing
                // it would make more sense to iterate in the other way (first iterate over coefficient, and the iterate over the global_list)
                // 
                FastQuatLat J = create_from_generator_O0( gen, convert(ell_e));
                J.good_type = false; // TODO handle good_type = true
                assert(J.basis[2][0] >=0 );

                // and now we precompute some value used to speed-up the computation of fast_right_order_and_gram

                // if the ideal is split, we know it is equivalent to O0 and so we have it already
                bool is_split = J.basis[1][0] != 0;
                if (is_split) {
                    continue;
                }

                if (num_prime != 0) {
                    // precompute some value to help ultimately with the maximal order computation
                    if (J.good_type) {
                        assert(J.basis[0][2] == 0 && J.basis[0][3] == 0 && J.basis[1][2] == 0 && J.basis[1][3] == 0);
                        J.basis[0][3] = red_elle.mod((- J.basis[2][1]) * InvMod(J.basis[2][0], red_elle));
                        J.basis[0][2] = red_elle.modsqr(InvModSqr( red_elle.modsqr(J.basis[3][0] * J.basis[3][0] + J.basis[3][1] * J.basis[3][1] - fast_Bp.p), red_elle) << 1);
                        J.basis[1][3] = red_elle.modsqr(J.basis[0][2] * J.basis[3][1]);
                        J.basis[1][2] = red_elle.modsqr( - J.basis[0][2] * J.basis[3][0]);
                    } 
                    J._fast_intersect(K_id, red_elle_S, red_S, inv, inv2sqr);
                } 
                else {
                    // no need to intersect because K_id is O0
                    // we still need to make the precomputation
                    if (J.good_type) {
                        assert(J.basis[0][2] == 0 && J.basis[0][3] == 0 && J.basis[1][2] == 0 && J.basis[1][3] == 0);
                        J.basis[0][3] = red_elle.mod((- J.basis[2][1]) * InvMod(J.basis[2][0], red_elle));
                        J.basis[0][2] = red_elle.modsqr(InvModSqr( red_elle.mod2sqr(J.basis[3][0] * J.basis[3][0] + J.basis[3][1] * J.basis[3][1] - fast_Bp.p) >> 1 , red_elle));
                        J.basis[1][3] = red_elle.mod2sqr(J.basis[0][2] * J.basis[3][1]);
                        if (IsOdd(J.basis[0][2])) {
                            J.basis[1][2] = red_elle.mod2sqr( red_elle.mod2sqr(- J.basis[0][2]) * J.basis[3][0]);
                        }
                        else {
                            J.basis[1][2] = red_elle.mod2sqr( red_elle.mod2sqr(- J.basis[0][2]) * J.basis[3][0] + red_elle.bsqr);
                        }
                    }

                }
                // std::cout << "J = " << J << "\n";


                auto [rigo,gram] = J.fast_right_order_and_gram(red_elle_S);
          
                // std::cout << rigo << "\n";
                if (gram[0][0] == 0) {
                    // this shouldn't happen
                    assert(0);
                    // quatlat rigo_test = J.right_order();
                    // rigo = rigo_test;
                    // gram = compute_HNF_gram_order(&rigo);
                    
                }
                FastQuat small_endo = {{0, 0, 0, 0, 1}, fast_O0.alg};
                std::pair <FastQuat,FastQuat> small_endo_pair = {small_endo, small_endo};

                Key K = fast_order_invariant_computation_from_gram(rigo, gram, small_endo_pair);
                // print_key(K); std::cout << "\n";

#ifndef NDEBUG 
                quat gen_test = (gamma * (quat({{ std::get<0>(coeff), Integer(0), Integer(0), Integer(0), Integer(1)}, Bp}) + betabar * std::get<1>(coeff)  ));
                auto Jtest = create_from_generator_O0( gen_test, ell_e);
                quatlat K_test = K_id.FastQuatLat_to_quatlat(Bp);
                Jtest._fast_intersect(K_test);
                
                auto [rigotest,gramtest] = Jtest.fast_right_order_and_gram();
                
                quat ggamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
                std::pair<quat,quat> gamma_pair = {ggamma, ggamma};
                Key Ktest = order_invariant_computation_from_gram(rigotest, gramtest, &gamma_pair);


                if (!is_key_equal(K,Ktest)) {
                    std::cout << "Keys are not equal !\n";
                    std::cout << "\n\n" << gen << "\n" << gen_test << "\n" << J << "\n" << Jtest << "\n";
                    std::cout << rigo.basis << "\n";
                    std::cout << rigotest.basis << "\n";
                        // 
                    // std::cout << gramtest << "\n";
                    // std::cout << gram << "\n";   
                    print_key(K);
                    print_key(Ktest);
                    std::cout << "\n";
                    assert(0);
                }
                
                // std::cout << "rigotest = " << rigotest << "\n";
#endif 

                
                // Key K = order_invariant_computation_from_gram(rigo, gram, &small_endo_pair);

                quaternion_time += clock() - t;
                quat_count++;
                auto a_loc = local_map.insert({K, j0});
                // if this order was not computed already
                if (a_loc.second) {

                    t = clock();
                    isog_count++;
                    ecp ker = std::get<0>(coeff) * P + std::get<1>(coeff) * Q;
                    assert((ell_e_m_1 * ker));
                    assert(!(ell_e * ker));
                    // now to the isogeny computation
                    std::vector<std::pair<ecp,std::pair<int, int>>> kerGen = {{ker, {l,fac_list[l]}} };
                    // computing the isogeny
                    isog_chain phi = isog_chain(kerGen);
                    // pushing the other bases
                    std::list<std::pair<ecp,ecp>> new_torsion_list = {};
                    for (auto bas : TorBas) {
                        new_torsion_list.push_back( {phi(bas.first), phi(bas.second)});
                    }
                    FpE_elem j = phi.get_codomain().j_invariant();
                    isogeny_time += clock() - t;
                    new_list.push_back({J, phi.get_codomain(), new_torsion_list});

                    // onto the weber computation
                    bas2 = new_torsion_list.back();
                    new_torsion_list.pop_back();
                    bas3 = new_torsion_list.back();
                    new_torsion_list.push_back(bas2);

                    weber_bas web = { bas2.first, bas2.second, bas3.first, bas3.second};

                    assert((16*web.P16).is_identity());
                    assert((16*web.Q16).is_identity());
                    assert((3*web.P3).is_identity());
                    assert((3*web.Q3).is_identity());
                    assert(J.basis[3][0] >= 0);

                    FastQuat prod_quat = small_endo_pair.first * small_endo_pair.second;
                    prod_quat[0] = 0;
                    prod_quat[4] *= 2;
                    prod_quat.normalize();
                    
                    auto nn = prod_quat.integral_norm() / small_endo_pair.first.alg.p;
                    assert(prod_quat.integral_norm() % small_endo_pair.first.alg.p == 0);

                    auto nJ = J.norm().first/J.norm().second;

                    auto a = m.insert({K, {j, {{nJ, {small_endo_pair, {nn, prod_quat}}}, w0_data}}});

                    // if it was not already inserted
                    if (a.second) {
                        t = clock();
                        assert(J.basis[3][0] >= 0);
                        assert(J.basis[2][0] >= 0);
                        auto it = m.find(K);
                        assert(it!=m.end());

                        // computation of the weber invariants
                        it->second.second.second = EnumerateAllWeberFastFast(web, Fexts, &fast_precomp);
                        weber_time += clock() - t;
                        if (J.basis[3][0] < 0) {
                            J.basis[3][0] += 2 * J.norm().first/J.norm().second;
                        }
                        if (J.basis[2][0] < 0) {
                            J.basis[2][0] += 2 * J.norm().first/J.norm().second;
                        }

                        order_count++;
                        if (!is_Fp(j)) {
                            order_count++;
                        }

                    }
                }
                if (order_count >= target_num) {
                    goto endloop;
                }


            }

        }
        global_list = new_list;
        num_prime++;

    }
    std::cerr << "exhausted the list and did not find all the curves for p =" << O0.alg.p <<std::endl;
    // shouldn't happen
    assert(0);

endloop:

    // now is time to compute a list of very small ideals that we are going to use for our enumeration
    std::vector<std::pair<FastQuatLat, Key>> small_id_list;
    std::unordered_set<Key,KeyHash,KeyEqual> Keys;

    Keys.insert(K0);

    Integer small_prime_norm = Integer(5);

    auto weber_coeff_list = EnumerateAllWeberCoeff();


    SignedBarrettReducer redp = SignedBarrettReducer(fast_Bp.p);

    // this should be the number of orders we need
    // TODO a possibility to reduce the norm we need would be to consider ideal whose norm is product of small primes as well
    while ((int) Keys.size() < std::min(5 + (NTL::conv<int>(coprime) / 100), (int) m.size()) ) {

        std::vector<FastQuatLat> norm_id_list = left_ideals_of_prime_norm_O0( convert(small_prime_norm), Bp, fast_Bp);

        SignedBarrettReducer redn = SignedBarrettReducer(convert(small_prime_norm));

        for (auto & id : norm_id_list) {

            if (id.basis[3][0] < 0) {
                id.basis[3][0] += 2 * convert(small_prime_norm);
            }
            if (id.basis[2][0] < 0) {
                id.basis[2][0] += 2 * convert(small_prime_norm);
            }
                
            if (id.basis[3][1] == 0) {
                id.good_type = false;
            }

            if (id.good_type) {
                assert(id.basis[0][2] == 0 && id.basis[0][3] == 0 && id.basis[1][2] == 0 && id.basis[1][3] == 0);
                id.basis[0][3] = redn.mod((- id.basis[2][1]) * InvMod(id.basis[2][0], redn));
                id.basis[0][2] = redn.modsqr(InvModSqr( redn.mod2sqr(id.basis[3][0] * id.basis[3][0] + id.basis[3][1] * id.basis[3][1] - fast_Bp.p) >> 1 , redn));
                id.basis[1][3] = redn.mod2sqr(id.basis[0][2] * id.basis[3][1]);
                if (IsOdd(id.basis[0][2])) {
                    id.basis[1][2] = redn.mod2sqr( redn.mod2sqr(- id.basis[0][2]) * id.basis[3][0]);
                }
                else {
                    id.basis[1][2] = redn.mod2sqr( redn.mod2sqr(- id.basis[0][2]) * id.basis[3][0] + redn.bsqr);
                }
            }
                

            auto [rigo,gram] = id.fast_right_order_and_gram(redn);
            
            if (gram[0][0] != 0 ) {


                FastQuat small_endo = {{0, 0, 0, 0, 1}, fast_Bp};
                std::pair <FastQuat, FastQuat> small_endo_pair = {small_endo, small_endo};
                Key K = fast_order_invariant_computation_from_gram(rigo, gram, small_endo_pair);
                

#ifndef NDEBUG 
                auto Jtest = id.FastQuatLat_to_quatlat(Bp);
                auto [rigotest,gramtest] = Jtest.fast_right_order_and_gram();
                quat ggamma = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
                std::pair<quat,quat> gamma_pair = {ggamma, ggamma};
                Key Ktest = order_invariant_computation_from_gram(rigotest, gramtest, &gamma_pair);


                if (!is_key_equal(K,Ktest)) {
                    // std::cout << "Keys are not equal !\n";
                    // std::cout << "\n\n" << gen << "\n" << gen_test << "\n" << id << "\n" << Jtest << "\n";
                    std::cout << rigo.basis << "\n";
                    std::cout << rigotest.basis << "\n";
                        // 
                    // std::cout << gramtest << "\n";
                    // std::cout << gram << "\n";   
                    print_key(K);
                    print_key(Ktest);
                    std::cout << "\n";
                    assert(0);
                }
#endif 
                
                auto a = Keys.insert(K);
                
                // the key is not contained already we can add the ideal
                if (a.second) {
                    
                                
                    // and now we need to modify m to contain the smaller ideal 
                    auto j_it = m.find(K)->second;

                    Fp2 j_inv = j_it.first;

                    // use in priority fp2 curves as they are more efficient due to the symmetry coming from Galois conjugates
                    // This is particularly important for the CRT variant 
                    if (is_Fp(j_inv)) {
                        id_list.push_back({id, K});
                    }
                    else {
                        id_list.insert(id_list.begin(), {id, K});
                    }

                    

                    // the first step is to compute the isomorphism to go from the ideal currently stored in m to the new small ideal
                    auto normJ = j_it.second.first.first;
                    auto is_equiv = commutatorfind(small_endo_pair, j_it.second.first.second.first, {convert(small_prime_norm), normJ}, j_it.second.first.second.second.second, j_it.second.first.second.second.first, is_Fp(j_it.first), small_endo, redp);
                    // std::cout << "alphaJ " << j_it.second.first.second.first.first << " small endo = " << small_endo << "\n";
                    // std::cout << "weber enum" << j_it.second.second.enumerator[0].first << " " << j_it.second.second.enumerator[1].first << "\n";

                    // and now we compute the change of basis matrix
                    std::pair<VerySmallMat,VerySmallMat> new_web;
                    // small_endo = small_endo.conjugate();

                    new_web = WeberBasApplicationRemoteEndo(mat0, small_endo); 

                    // now we apply the change of matrix to the weber enumerator and inv_list
                    weber_enum new_web_enum;
                    int count = 0;
                    for (auto [w, coeffs] : j_it.second.second.enumerator) {

                        unsigned index = WeberGetFromEnum(coeffs, new_web);
                        auto new_w = j_it.second.second.inv_list[index];

                        if (is_equiv) {
                            new_web_enum[count] = {new_w, coeffs};
                        }
                        else {
                            Fp_integer some_mod;
                            NTL::conv(some_mod, Bp.p);
                            Fp_push push(some_mod);
                            new_web_enum[count] = {fast_Frob(new_w), coeffs};
                        }
                        count++;
                    }
                    
                    weber_inv_list new_web_list;
                    count = 0;
                    for (auto coeffs : weber_coeff_list) {

                        unsigned index = WeberGetFromEnum(coeffs, new_web);
                        auto new_w = j_it.second.second.inv_list[index];

                        if (is_equiv) {
                            new_web_list[count] = new_w;
                        }
                        else {
                            Fp_integer some_mod;
                            NTL::conv(some_mod, O0.alg.p);
                            Fp_push push(some_mod);
                            new_web_list[count] = fast_Frob(new_w);
                        }
                        count++;
                    }

                    // and now we can modify j_it 

                    // applying the Frobenius if needed
                    if (!is_equiv) {
                        Fp_integer some_mod;
                        NTL::conv(some_mod, O0.alg.p);
                        Fp_push push(some_mod);
                        j_it.first = Frob(j_it.first);
                    }

                    // replacing the norm of the ideal
                    j_it.second.first.first = convert(small_prime_norm);
                    // normally the smallest element is not used anymore 
                    // replacing the pair of small elements
                    j_it.second.first.second.first = small_endo_pair;

                    // computing the product as small_endo and replacing it
                    {
                        small_endo = small_endo_pair.first * small_endo_pair.second;
                        small_endo[0] = 0;
                        small_endo[4] *= 2;
                        small_endo.normalize();
                        auto nn = small_endo.integral_norm() / small_endo_pair.first.alg.p;
                        assert(small_endo.integral_norm() % small_endo_pair.first.alg.p == 0);

                        

                        // and replacing it
                        j_it.second.first.second.second = {nn, small_endo};
                    }

                    // replacing the weber enumerator
                    j_it.second.second.enumerator = new_web_enum;
                    j_it.second.second.inv_list = new_web_list;
                    
                    
                    // and making the change in m
                    m.insert_or_assign(K, j_it);                    
                    
                }
            }
            
        }


        small_prime_norm = NTL::NextPrime(small_prime_norm + 1);
    }

    // std::cout << "number of curves :" << order_count;
    // std::cout << "total time for quaternion operations : " << (double) (quaternion_time)/(CLOCKS_PER_SEC) << "\n";
    // std::cout << "total time for isogeny operations    : " << (double) (isogeny_time)/(CLOCKS_PER_SEC) << "\n";
    // std::cout << "total time for weber operations in the enumeration  : " << (double) (weber_time)/(CLOCKS_PER_SEC) << "\n";

    return {bas0, mat0};

}
