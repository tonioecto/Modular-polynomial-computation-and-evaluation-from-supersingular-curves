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

quat find_quaternion_iterator(std::list<int>& prime_list, const quatlat& I, const quatlat& O0, const quatalg &Bp) {

    // std::cout << O0 << "\n";

    int is_first = 1;

    quat result = { {NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(1), NTL::ZZ(1), NTL::ZZ(1)} , Bp};

    NTL::ZZ accumulated_multiple(1);
    for (auto l : prime_list) {
        NTL::ZZ ell(l);
        // generating the iterating quaternion mod ell
        // we need to find a quaternion element whose norm is coprime to ell that is not contained in ZZ + (I+ O_0 ell)
        quatlat I_ell = I.copy();
        I_ell.reduce_norm_cyclic(O0, ell);

        assert(std::get<0>(I_ell.norm())/std::get<1>(I_ell.norm()) == ell);
        quatlat Eich_ell = I_ell._compute_order(true);
        Eich_ell._intersect(I_ell.left_order());
        bool found = false;
        for (int i1 = 0; i1 < ell && !found; i1++) {
            for (int i2 = 0; i2 < ell && !found; i2++) {
                for (int i3 = 0; i3 < ell && !found; i3++) {
                    for (int i4 = 0; i4 < ell && !found; i4++) {
                        quat gen = { { NTL::ZZ(i1), NTL::ZZ(i2), NTL::ZZ(i3), NTL::ZZ(i4), NTL::ZZ(1)} , Bp};

                        if ((NTL::GCD(ell,O0.denom) == ell) || ell == Bp.q) {
                            gen[4] = ell;
                            NTL::vec_ZZ vec;
                            vec.SetLength(4);
                            vec[0] = gen[0];
                            vec[1] = gen[1];
                            vec[2] = gen[2];
                            vec[3] = gen[3];
                            vec = vec * O0.basis;
                            gen[0] = vec[0];
                            gen[1] = vec[1];
                            gen[2] = vec[2];
                            gen[3] = vec[3];
                        }
                        else {
                            NTL::vec_ZZ vec;
                            vec.SetLength(4);
                            vec[0] = gen[0];
                            vec[1] = gen[1];
                            vec[2] = gen[2];
                            vec[3] = gen[3];
                            vec = vec * O0.basis;
                            gen[0] = vec[0];
                            gen[1] = vec[1];
                            gen[2] = vec[2];
                            gen[3] = vec[3];
                        }



                        auto pair_norm = gen.norm();
                        NTL::ZZ norm = std::get<0>(pair_norm)/std::get<1>(pair_norm);
                        assert(O0.contains(gen));
                        if (NTL::GCD(norm,ell) == 1) {
                            // now we check the second condition
                            if (!Eich_ell.contains(gen)) {
                                // for (int scal=0; scal < ell; scal++) {
                                //     quat check_quat = { {NTL::ZZ(scal), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1)} , Bp};
                                //     check_quat = check_quat + gen;
                                //     assert(!I_ell.contains(check_quat));
                                // }
                                for (int i=0; i<4; i++) {
                                    NTL::ZZ acc = accumulated_multiple;
                                    if (is_first) {
                                        result[i] = gen[i];
                                    }
                                    else {
                                        NTL::ZZ m1, m2;
                                        NTL::mul(m1, result[i], gen[4]);
                                        NTL::mul(m2, gen[i], result[4]);
                                        NTL::CRT(m1, acc, m2, ell);
                                        result[i] = m1;
                                    }
                                }
                                accumulated_multiple *= ell;
                                result[4] *= gen[4];
                                if (is_first) {
                                    is_first = 0;
                                }
                                found = true;
                                auto pn = result.norm();
                            }

                        }

                    }
                }
            }
        }
        if (!found) {
            std::cout << "quat iter loop failed \n";
        }
    }



    return result;

}

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

    // std::cout << beta << "\n";

    std::map<NTL::ZZ,std::pair<ecp,ecp>> result_map;
    std::map<NTL::ZZ,std::pair<ecp,ecp>> TorsionBases;
    for (const auto& [ell,e] : facN) {
        // std::cout << ell << "\n";
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
            //std::cout << "Constructing torsion basis of E[" << ell << "^" << e << "] over Fp2k, k = " << k << std::endl;
            auto jt = FieldExtensions.find(k);
            assert(jt != FieldExtensions.end());  //TODO
//            if (jt == FieldExtensions.end())
//                jt = FieldExtensions.emplace(std::make_pair(k, Fp2k {k})).first;
            assert(jt->first == jt->second.k);
            auto bas = E0.torsionBasis(jt->second, ell, e + extra1 + extra2);
            it = TorsionBases.emplace(std::make_pair(ell_ext, bas)).first;
        }
        //std::cout << "Finding kernel generator of order" << ell << "^" << e << std::endl;
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
        //  * (extra2 ? ((int) beta[4]) : 1);
        //  * (!extra2 ? (NTL::InvMod(beta[4], ell_e)));
        ecp K2 = evalEndo(sec_endo,K,ell_ext/(extra1 ? (ell) : 1));

        K = (extra2 ? (ell) : 1) * K;

        assert(!(ell_e * K));
        assert((ell_e_m_1) * K);
        assert(!(ell_e * K2));
        assert((ell_e_m_1) * K2);

        assert(!((K-K2).is_identity()));
        assert(!((ell_e_m_1 * K- ell_e_m_1 * K2).is_identity()));
        result_map.emplace(std::make_pair(ell_e, std::make_pair(K,K2)));


        //std::cout << "Done!" << std::endl;
        // kerGens.push_back(std::pair<ecp,std::pair<int, int>>(K, {ell,e}));
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
        quat gamma_tmp = RepresentInteger(Bp, target);
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


std::pair<weber_bas,std::vector<std::pair<SmallMatFp,SmallMatFp>>> order_to_weber_inv_full_list(std::unordered_map<Key, std::pair<FpE_elem, std::pair<std::pair<quatlat,quat>, weber_full_data>>, KeyHash, KeyEqual> &m, std::vector<std::pair<quatlat,std::pair<quatlat, Key>>> &id_list, const  NTL::ZZ &p, const quatalg &Bp, const std::map<unsigned,Fp2k> &Fexts, const Integer &coprime) {

    NTL::ZZ q = Bp.q;

    // we start by creating the map and the quaternion algebra
    // std::cout << Bp.p << " " << Bp.q << "\n";
    auto start = starting_curve(Bp, false);
    quatlat O0 = start.second;



    Fp2 j0 = start.first.j_invariant();
    quat dum = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
    Key K0 = order_invariant_computation(O0, &dum);

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
    // std::cout << S << "\n";
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
        quat gamma_tmp = RepresentInteger(Bp, target);
        quat testing = {{gamma_tmp[0], gamma_tmp[1], gamma_tmp[2], gamma_tmp[3], NTL::ZZ(2)*gamma_tmp[4]}, Bp};
        found = !O0.contains(testing);
        gamma[0] = gamma_tmp[0];
        gamma[1] = gamma_tmp[1];
        gamma[2] = gamma_tmp[2];
        gamma[3] = gamma_tmp[3];
        gamma[4] = gamma_tmp[4];
        target = cofac * target;
    }

    assert(O0.contains(gamma));

    quatlat I = O0 * gamma + O0 * S;
    assert(std::get<0>(I.norm())/std::get<1>(I.norm()) == S);


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
    // now we generate 3 and 16 torsion bases
    std::pair<ecp,ecp> bas2 = torsion_list.front();
    std::pair<ecp,ecp> bas3 = torsion_list.front();
    bool extra3 = Bp.q == 3;
    {
        unsigned k = torsionToFieldDegree(Integer(3) * (extra3 ? (Integer(3)): Integer(1)));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());
        assert(jt->first == jt->second.k);
        bas3 = start.first.torsionBasis(jt->second, int(3), 1 + (extra3 ? (1): 0));
    }
    {
        unsigned k = torsionToFieldDegree(Integer(32));
        auto jt = Fexts.find(k);
        assert(jt != Fexts.end());
        assert(jt->first == jt->second.k);
        bas2 = start.first.torsionBasis(jt->second, int(2), 5);
    }
    // std::pair<ecp,ecp> bas = {bas2.first + bas3.first, bas2.second + bas3.second};


    weber_bas bas0 = {bas2.first, bas2.second, bas3.first, bas3.second};
    bas2.first = 2* bas2.first;
    bas2.second = 2* bas2.second;
    bas3.first =  (extra3 ? (Integer(3)) : Integer(1)) * bas3.first;
    bas3.second =  (extra3 ? (Integer(3)) : Integer(1)) * bas3.second;
    torsion_list.push_back(bas3);
    torsion_list.push_back(bas2);
    assert((16*bas2.first).is_identity());
    assert((3*bas3.first).is_identity());


    // precomputing the action of the endo ring on the weber_basis
    std::vector<std::pair<SmallMatFp,SmallMatFp>> mat0 = {};
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
        // assert(!(8*gammaP016).is_identity());
        assert((16*gammaQ016).is_identity());
        // assert(!(8*gammaQ016).is_identity());
        // assert(!(8*gammaQ016 - 8*gammaP016).is_identity());

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
        // assert(!(gammaP03).is_identity());
        assert((3*gammaQ03).is_identity());
        // assert(!(gammaQ03).is_identity());
        // assert(!(gammaQ03 - gammaP03).is_identity());

        // now we compute the matrices of order 16 and 3
        auto M16 = change_of_basis16( {2*bas0.P16, 2*bas0.Q16}, {gammaP016, gammaQ016} );
        auto M3 = change_of_basis3( {(extra3 ? (Integer(3)): Integer(1)) * bas0.P3, (extra3 ? (Integer(3)): Integer(1)) * bas0.Q3}, {gammaP03, gammaQ03} );
        mat0.push_back({M16,M3});

    }



    std::list<std::tuple<quatlat,ec,std::list<std::pair<ecp,ecp>>>> global_list = {{O0,start.first,torsion_list}};

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
    weber_bas w0 = { bas2.first, bas2.second, bas3.first, bas3.second};

    weber_enum_poly_precomp precomp = SetWeberPrecomp();

    weber_enum wb0 = EnumerateAllWeberFast(w0, Fexts, &precomp);
    weber_full_data w0_data = {w0,wb0};
    quat ii = { {Integer(0),Integer(1),Integer(0),Integer(0), Integer(1)}, O0.alg};
    std::pair<std::pair<quatlat,quat>,weber_full_data> bos = { {O0, ii}, w0_data };
    m.insert({K0, {j0, bos }});





    std::unordered_map<Jinv, std::pair<quatlat,Key> ,JinvHash, JinvEqual> j_inv_list = {};

    // iterating through the set of S-ideals, we go factor by factor
    int order_count = 1;
    int count_mistake = 0;
    // std::cout << "target = " << target_num <<"\n";

    clock_t quaternion_time = 0;
    clock_t isogeny_time = 0;
    clock_t weber_time = 0;
    clock_t t;

    int isog_count = 0;
    int quat_count = 0;

        for (auto l : prime_list) {
        // std::cout << "round for" << l << " \n";
        std::list<std::tuple<quatlat,ec,std::list<std::pair<ecp,ecp>>>>  new_list = {};
        std::list<IntegerPair> coeff_list = {};
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
            quatlat K_id = std::get<0>(tup);
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
                assert(J.basis[2][0]>=0);

                // if the ideal is split, we know it is equivalent to O0 and so we have already
                bool is_split = J.basis[1][0] != 0;
                if (is_split) {
                    continue;
                }

                J._fast_intersect(K_id);

#ifndef NDEBUG
                // testing
                quatlat J_test = create_from_generator( gen, ell_e, O0);
                J_test._intersect(K_id);
                for (int i =0; i<4; i++) {
                    quat bas_el = quat({{J.basis[i][0], J.basis[i][1], J.basis[i][2], J.basis[i][3], J.denom}, J.alg});
                    if (!J_test.contains(bas_el)) {

                        std::cout << create_from_generator_O0(gen, ell_e).basis << "\n";
                        std::cout << K_id.basis << "\n";
                        std::cout << J.basis << "\n";
                        std::cout << J_test.HNF_basis() << "\n";
                        std::cout << i << "-th vector not in J_test \n";
                    }
                    assert(J_test.contains(bas_el));
                    auto mtest = J_test.HNF_basis();
                    quat bas_el_test = quat({{mtest[i][0], mtest[i][1], mtest[i][2], mtest[i][3], J_test.denom}, J.alg});
                    if (!J.contains(bas_el_test)) {
                        std::cout << J.basis << "\n";
                        std::cout << J_test.HNF_basis() << "\n";
                        std::cout << bas_el_test << "\n";
                        std::cout << "i=" << i << " is not contained \n";
                    }
                }


#endif
                auto [rigo,gram] = J.fast_right_order_and_gram();
                gram[0][0] = 0;
                if (gram[0][0] == 0) {
                    quatlat rigo_test = J.right_order();
                    rigo = rigo_test;
                    gram = compute_gram_order(rigo);
                }
                quat small_endo = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};
                // Key Ktest = order_invariant_computation(rigo, &small_endo);
                Key K = order_invariant_computation_from_gram(rigo, gram, &small_endo);
                // quat small_endo = {{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, O0.alg};

                quaternion_time += clock() - t;
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
                    // auto bas = new_torsion_list.back();
                    new_list.push_back({J,phi.get_codomain(),new_torsion_list});
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
                    auto a = m.insert({K, {j, {{J,small_endo}, {web, {}} }}});

                    if (a.second) {
                        t = clock();
                        assert(J.basis[3][0] >= 0);
                        assert(J.basis[2][0] >= 0);
                        auto it = m.find(K);
                        assert(it!=m.end());
                        // it->second = nK;
                        it->second.second.second.enumerator = EnumerateAllWeberFast(web, Fexts, &precomp);
                        assert(it->second.second.second.enumerator.size() == 72);
                        for (size_t i = 0; i < 72; i++) {
                            if (it->second.second.second.enumerator[i].second.size() != 3) {
                                std::cout << i << " " << it->second.second.second.enumerator[i].second.size() << "\n";
                            }
                            assert(it->second.second.second.enumerator[i].second.size() == 3);
                        }
                        weber_time += clock() - t;
                        auto nnn = get_smallest_element(&(it->second.second.first.second),J);
                        if (J.basis[3][0] < 0) {
                            J.basis[3][0] += 2 * J.norm().first/J.norm().second;
                        }
                        if (J.basis[2][0] < 0) {
                            J.basis[2][0] += 2 * J.norm().first/J.norm().second;
                        }
                        // check if the element is not already contained in the list
                        // TODO this is only used for debug purpose
                        Jinv jtemp = JToJinv(j);
                        auto search = j_inv_list.find(jtemp);
                        assert(J.basis[3][0] >= 0);
                        if (search != j_inv_list.end()) {
                            count_mistake ++;
                            std::cout << "\nMISTAKE \n";
                            quatlat test = J.conjugate() * search->second.first;
                            isPrincipal(test);
                            print_key(order_invariant_computation(rigo, &dum));
                            print_key(K);
                            print_key(search->second.second);
                            print_key(order_invariant_computation(search->second.first.right_order(),&dum));
                            std::cout << "\n";

                        }
                        assert(J.basis[3][0] >= 0);
                        order_count ++;
                        Jinv jnew = JToJinv(j);
                        j_inv_list.insert({jnew,{J,K}});
                        assert(J.basis[3][0] >= 0);
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
    std::cerr << "exhausted the list and did not find all the curves for p =" << O0.alg.p <<std::endl;
    assert(0);

endloop:

    id_list.push_back( {O0, {O0, K0}} );

    assert(count_mistake == 0);
    // std::cout << "number of curves :" << order_count;
    std::cout << "total time for quaternion operations : " << (double) (quaternion_time)/(CLOCKS_PER_SEC) << "\n";
    std::cout << "total time for isogeny operations    : " << (double) (isogeny_time)/(CLOCKS_PER_SEC) << "\n";
    std::cout << "total time for weber operations in the enumeration  : " << (double) (weber_time)/(CLOCKS_PER_SEC) << "\n";

    return {bas0, mat0};

}
