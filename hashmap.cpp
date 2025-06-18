///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   The purpose of this file is to build code that computes an invariant corresponding to
////     a maximal order, for the purpose of using in a hash table
////
////   The invariant is computed by viewing the order as a quaternion lattice, and computing
////    its successive minima
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////

#include "hashmap.hpp"

NTL::mat_ZZ compute_gram_order(quatlat const &order) {
    NTL::mat_ZZ gram;
    gram.SetDims(3,3);
    NTL::mat_ZZ HNF;
    if (order.basis[0][0] == order.denom && order.basis[0][1]==0 && order.basis[0][2]==0 && order.basis[0][3]==0){
        HNF = order.basis;
    }
    else {
        HNF = order.HNF_basis();
    }
    std::array<quat, 3> elts {{
        {{NTL::ZZ(0), 2*(HNF[1][1]), 2*(HNF[1][2]), 2*(HNF[1][3]), order.denom}, order.alg},
        {{NTL::ZZ(0), 2*(HNF[2][1]), 2*(HNF[2][2]), 2*(HNF[2][3]), order.denom}, order.alg},
        {{NTL::ZZ(0), 2*(HNF[3][1]), 2*(HNF[3][2]), 2*(HNF[3][3]), order.denom}, order.alg},
    }};


    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            auto pair = elts[i] * elts[j].conjugate();
            assert(pair[0] %(pair[4]) ==0);
            gram[i][j] = pair[0]/(pair[4]);
        }
    }
    return gram;

}

void inner_order_invariant_computation(NTL::ZZ *IntList
,quatlat const &order, quat *small
) {

    (void) small;
    // taking the HNF
    NTL::mat_ZZ HNF;
    if (order.basis[0][0] == order.denom && order.basis[0][1]==0 && order.basis[0][2]==0 && order.basis[0][3]==0){
        // std::cout << "no HNF \n";
        HNF = order.basis;
    }
    else {
        // refac = order.denom;
        HNF = order.HNF_basis();
        // std::cout << HNF << " \n";
        // std::cout << order.denom << "\n";
    }

    // now we select the three basis elements that we are interested in
    // it suffices to take the 2nd,3rd and 4th basis elements and apply x -> 2x - tr(x)
    std::array<quat, 3> elts {{
        {{NTL::ZZ(0), 2*(HNF[1][1]), 2*(HNF[1][2]), 2*(HNF[1][3]), order.denom}, order.alg},
        {{NTL::ZZ(0), 2*(HNF[2][1]), 2*(HNF[2][2]), 2*(HNF[2][3]), order.denom}, order.alg},
        {{NTL::ZZ(0), 2*(HNF[3][1]), 2*(HNF[3][2]), 2*(HNF[3][3]), order.denom}, order.alg},
    }};

    mat_t Gram;
    mat_t Red;
    Gram.resize(3, 3);
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {

            auto pair = elts[i] * elts[j].conjugate();

            // if (pair[0] %(pair[4]) !=0) {
            //     std::cout << i << " " << j << "\n";
            // }
            assert(pair[0] %(pair[4]) ==0);
            auto val = pair[0]/(pair[4]);

            ntl2gmp(Gram[i][j].get_data(), val);
        }
    }

    // std::cout << "Gram time = " << (tic() - t) << "\n";
    Red.gen_identity(3);
    mat_t Redinv;
    gso_t gso(Gram, Red, Redinv, fplll::GSO_INT_GRAM);
    // std::cout << "delta = " << fplll::LLL_DEF_DELTA << " eta = " << fplll::LLL_DEF_ETA << " \n";

    fplll::LLLReduction<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> lllobj(gso,
    fplll::LLL_DEF_DELTA,
    fplll::LLL_DEF_ETA,
    // 0.9999999,
    // 0.5000001,
    0);
    // std::vector<fplll::Strategy> strat = fplll::load_strategies_json(fplll::strategy_full_path(fplll::default_strategy()));
    // int blocksize = 3;
    // fplll::BKZParam bkz_param(blocksize, strat);
    // fplll::BKZReduction bkz_red(gso, lllobj, bkz_param);
    // bkz_red.bkz();

    // std::cout << fplll::default_strategy();

    lllobj.lll();
    // std::cout << "LLL time = " << (tic() - t) << "\n";

    {
        NTL::mat_ZZ gram_test;
        gram_test.SetDims(3,3);
        for (unsigned j = 0; j < 3; ++j) {
            auto &c = Gram[j][j].get_data();
            // assert(mpz_divisible_p(c, d));
            // mpz_divexact(c, c, d);
            gmp2ntl(IntList[j], c);
            for (int i=0; i<3; i++) {
                auto &c = Gram[j][i].get_data();
                gmp2ntl(gram_test[j][i], c);
            }
        }


        // reordering if needed
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
        }

        int bound = 1;

        // gram_test
        for (int i1 = -bound; i1 < bound + 1; i1++) {
            for (int i2 = -bound; i2 < bound + 1; i2++) {
                for (int i3 = -bound; i3 < bound + 1; i3++) {
                    NTL::mat_ZZ v1,v2;
                    v1.SetDims(1,3);
                    v2.SetDims(3,1);
                    v1[0][0] = NTL::ZZ(i1);
                    v2[0][0] = NTL::ZZ(i1);
                    v1[0][1] = NTL::ZZ(i2);
                    v2[1][0] = NTL::ZZ(i2);
                    v1[0][2] = NTL::ZZ(i3);
                    v2[2][0] = NTL::ZZ(i3);
                    NTL::mat_ZZ res = v1 * gram_test * v2;
                    auto norm = res[0][0];
                    assert(norm >= 0);
                    if (i3==0 && i2!=0 && norm < IntList[1]) {
                        // std::cout << "gram_test fails !! \n";
                        IntList[1] = norm;
                    }
                    if (i3!=0 && norm < IntList[2]) {
                        IntList[2] = norm;
                        // std::cout << "gram_test fails !! \n";
                    }
                }
            }
        }

        // reordering if needed
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
        }



    }


}

void inner_order_invariant_computation_from_gram(NTL::ZZ *IntList, NTL::mat_ZZ gram_input
// , quatlat order, quat *small
) {
    // (void) small;
    // (void) order;
    mat_t Gram;
    mat_t Red;
    Gram.resize(3, 3);
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            ntl2gmp(Gram[i][j].get_data(), gram_input[i][j]);
        }
    }

    // std::cout << "Gram time = " << (tic() - t) << "\n";
    // t = tic();
    Red.gen_identity(3);
    mat_t Redinv;
    gso_t gso(Gram, Red, Redinv, fplll::GSO_INT_GRAM);
    // std::cout << "delta = " << fplll::LLL_DEF_DELTA << " eta = " << fplll::LLL_DEF_ETA << " \n";

    fplll::LLLReduction<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> lllobj(gso,
    fplll::LLL_DEF_DELTA,
    fplll::LLL_DEF_ETA,
    // 0.9999999,
    // 0.5000001,
    0);

    lllobj.lll();

    std::array<std::array<FastInteger,3>,3> Red_Mat{{{0,0,0},{0,0,0},{0,0,0}}};
    for (int i = 0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // NTL::ZZ c;
            // gmp2ntl(Red_Mat[i][j], Red[i][j].get_data());
            Red_Mat[i][j] = mpz_get_si(Red[i][j].get_data());
        }
    }


    // std::cout << "Gram + LLL time = " << (tic() - t) << "\n";

    {
        std::array<std::array<FastInteger,3>,3> gram_test{{{0,0,0},{0,0,0},{0,0,0}}};
        for (unsigned j = 0; j < 3; ++j) {
            auto &c = Gram[j][j].get_data();
            // assert(mpz_divisible_p(c, d));
            // mpz_divexact(c, c, d);
            gmp2ntl(IntList[j], c);
            for (int i=0; i<3; i++) {
                auto &c = Gram[j][i].get_data();
                // gmp2ntl(NTL_gram_test[j][i], c);
                gram_test[i][j] = mpz_get_si(c);
            }
        }


        // reordering if needed
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            // std::swap(Red_Mat[1][0], Red_Mat[2][0]);
            // std::swap(Red_Mat[1][1], Red_Mat[2][1]);
            // std::swap(Red_Mat[1][2], Red_Mat[2][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        // std::cout << "read and swap = " << (tic() - t) << "\n";
        int bound = 1;
        // gram_test
        for (int i1 = -bound; i1 < bound + 1; i1++) {
            for (int i2 = -bound; i2 < bound + 1; i2++) {
                for (int i3 = 0; i3 < bound + 1; i3++) {

                    // NTL::mat_ZZ v1,v2;
                    // v1.SetDims(1,3);
                    // v2.SetDims(3,1);
                    auto norm = i1*i1 * gram_test[0][0] + i2*i2 * gram_test[1][1] + i3 * i3 * gram_test[2][2] +
                    2 * (i1 * (i2 * gram_test[0][1] + i3 * gram_test[0][2]) + i2*i3 * gram_test[1][2] );
                    // v1[0][0] = NTL::ZZ(i1);
                    // v2[0][0] = NTL::ZZ(i1);
                    // v1[0][1] = NTL::ZZ(i2);
                    // v2[1][0] = NTL::ZZ(i2);
                    // v1[0][2] = NTL::ZZ(i3);
                    // v2[2][0] = NTL::ZZ(i3);
                    // NTL::mat_ZZ res = v1 * NTL_gram_test * v2;
                    // auto norm = res[0][0];
                    // assert(norm == res[0][0]);
                    assert(norm >= 0);

                    if (i3==0 && i2!=0 && norm < IntList[1]) {
                        IntList[1] = Integer(norm);
                        // v1[0][0] = Integer(1);
                        // v1[0][1] = Integer(0);
                        // v1[0][2] = Integer(0);
                        // auto NTL_new01 = (v1 * NTL_gram_test * v2)[0][0];
                        auto new01 = i1 * gram_test[0][0] + i2 * gram_test[0][1] + i3 * gram_test[0][2];
                        // assert(new01 == NTL_new01);
                        // v1[0][0] = Integer(0);
                        // v1[0][2] = Integer(1);
                        // auto NTL_new12 = (v1 * NTL_gram_test * v2)[0][0];
                        auto new12 = i1 * gram_test[2][0] + i2 * gram_test[2][1] + i3 * gram_test[2][2];
                        // assert(new12 == NTL_new12);
                        gram_test[1][1] = norm;
                        gram_test[0][1] = new01;
                        gram_test[1][0] = new01;
                        gram_test[2][1] = new12;
                        gram_test[1][2] = new12;
                        // Red_Mat[1][0] = i1 * Red_Mat[0][0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
                    }
                    if (i3!=0 && norm < IntList[2]) {
                        // v1[0][0] = Integer(1);
                        // v1[0][1] = Integer(0);
                        // v1[0][2] = Integer(0);
                        // auto NTL_new02 = (v1 * NTL_gram_test * v2)[0][0];
                        auto new02 = i1 * gram_test[0][0] + i2 * gram_test[0][1] + i3 * gram_test[0][2];
                        // v1[0][0] = Integer(0);
                        // v1[0][1] = Integer(1);
                        // auto NTL_new12 = (v1 * NTL_gram_test * v2)[0][0];
                        auto new12 = i1 * gram_test[1][0] + i2 * gram_test[1][1] + i3 * gram_test[1][2];
                        gram_test[2][2] = norm;
                        gram_test[0][2] = new02;
                        gram_test[2][0] = new02;
                        gram_test[2][1] = new12;
                        gram_test[1][2] = new12;
                        IntList[2] = Integer(norm);
                        // Red_Mat[2] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
                    }
                }
            }
        }
        // std::cout << "combi time = " << (tic() - t) << "\n";
        // reordering if needed
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            // std::swap(Red_Mat[2][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[2][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[2][2], Red_Mat[1][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        for (int i1 = -bound; i1 < bound + 1; i1++) {
            for (int i2 = -bound; i2 < bound + 1; i2++) {
                for (int i3 = 0; i3 < bound + 1; i3++) {

                    // NTL::mat_ZZ v1,v2;
                    // v1.SetDims(1,3);
                    // v2.SetDims(3,1);
                    auto norm = i1*i1 * gram_test[0][0] + i2*i2 * gram_test[1][1] + i3 * i3 * gram_test[2][2] +
                    2 * (i1 * (i2 * gram_test[0][1] + i3 * gram_test[0][2]) + i2*i3 * gram_test[1][2] );
                    // v1[0][0] = NTL::ZZ(i1);
                    // v2[0][0] = NTL::ZZ(i1);
                    // v1[0][1] = NTL::ZZ(i2);
                    // v2[1][0] = NTL::ZZ(i2);
                    // v1[0][2] = NTL::ZZ(i3);
                    // v2[2][0] = NTL::ZZ(i3);
                    // NTL::mat_ZZ res = v1 * gram_test * v2;
                    // auto norm = res[0][0];
                    assert(norm >= 0);

                    if (i3==0 && i2!=0 && norm < IntList[1]) {
                        IntList[1] = Integer(norm);
                        // v1[0][0] = Integer(1);
                        // v1[0][1] = Integer(0);
                        // v1[0][2] = Integer(0);
                        // auto new01 = (v1 * gram_test * v2)[0][0];
                        auto new01 = i1 * gram_test[0][0] + i2 * gram_test[0][1] + i3 * gram_test[0][2];
                        // v1[0][0] = Integer(0);
                        // v1[0][2] = Integer(1);
                        // auto new12 = (v1 * gram_test * v2)[0][0];
                        auto new12 = i1 * gram_test[2][0] + i2 * gram_test[2][1] + i3 * gram_test[2][2];
                        gram_test[1][1] = norm;
                        gram_test[0][1] = new01;
                        gram_test[1][0] = new01;
                        gram_test[2][1] = new12;
                        gram_test[1][2] = new12;
                        // Red_Mat[1][0] = i1 * Red_Mat[0][0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
                    }
                    if (i3!=0 && norm < IntList[2]) {
                        // v1[0][0] = Integer(1);
                        // v1[0][1] = Integer(0);
                        // v1[0][2] = Integer(0);
                        // auto new02 = (v1 * gram_test * v2)[0][0];
                        auto new02 = i1 * gram_test[0][0] + i2 * gram_test[0][1] + i3 * gram_test[0][2];
                        // v1[0][0] = Integer(0);
                        // v1[0][1] = Integer(1);
                        // auto new12 = (v1 * gram_test * v2)[0][0];
                        auto new12 = i1 * gram_test[1][0] + i2 * gram_test[1][1] + i3 * gram_test[1][2];
                        gram_test[2][2] = norm;
                        gram_test[0][2] = new02;
                        gram_test[2][0] = new02;
                        gram_test[2][1] = new12;
                        gram_test[1][2] = new12;
                        IntList[2] = Integer(norm);
                        // Red_Mat[2] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
                    }
                }
            }
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            // std::swap(Red_Mat[2][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[2][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[2][2], Red_Mat[1][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            // std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            // std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            // std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }

    // NTL::Vec<NTL::ZZ> row;
    // row.SetLength(4);
    // for (unsigned j = 0; j < 3; ++j) {
    //     auto c = Red_Mat[2][j];
    //     row += 2 * c * order.basis[j+1];
    //     row[0] -= 2 * c* order.basis[j+1][0];
    // }

    // (*small)[0] = row[0];
    // (*small)[1] = row[1];
    // (*small)[2] = row[2];
    // (*small)[3] = row[3];
    // (*small)[4] = order.denom;
    // assert(order.contains(*small));
    // if (!(small->norm().first/small->norm().second == gram_test[0][0]))
    // {
    //     std::cout << small->norm().first/small->norm().second << " " << gram_test[0][0] << "\n";
    // }
    // assert(small->norm().first/small->norm().second == IntList[0] );
    }



}

Key order_invariant_computation(quatlat const &order, quat *small) {
        //////////////////////////////////////////////////////////////////
        /// Computes the invariant corresponding to the given order
        //////////////////////////////////////////////////////////////////

        Key K;
        NTL::ZZ IntList[3];

        // std::cout << clock() << "\n";
        // clock_t t = clock();
        // std::cout << "t =" << t << "\n";
        inner_order_invariant_computation(IntList, order, small);
        // clock_t t2 = clock();
        // std::cout << "t2 = " << t2 << "\n";
        // std::cout << "t2 - t =" << t2 - t << "\n";
        // std::cout << IntList[0] << IntList[1] << IntList[2] << "\n";
        NTL::BytesFromZZ(K.IntList[0], IntList[0], LenNumberBytes);
        NTL::BytesFromZZ(K.IntList[1], IntList[1], LenNumberBytes);
        NTL::BytesFromZZ(K.IntList[2], IntList[2], LenNumberBytes);
        return K;
}

Key order_invariant_computation_from_gram(quatlat order, NTL::mat_ZZ Gram, quat *small) {
        Key K;
        NTL::ZZ IntList[3];
        (void) order;
        (void) small;

        // std::cout << clock() << "\n";
        // clock_t t = clock();
        // std::cout << "t =" << t << "\n";
        // auto t = tic();
        inner_order_invariant_computation_from_gram(IntList, Gram
        // , order, small
        );
        // clock_t t2 = clock();
        // std::cout << "t2 = " << t2 << "\n";
        // std::cout << "t2 - t =" << t2 - t << "\n";
        // std::cout << IntList[0] << IntList[1] << IntList[2] << "\n";
        NTL::BytesFromZZ(K.IntList[0], IntList[0], LenNumberBytes);
        NTL::BytesFromZZ(K.IntList[1], IntList[1], LenNumberBytes);
        NTL::BytesFromZZ(K.IntList[2], IntList[2], LenNumberBytes);
        return K;
}

Integer get_smallest_element(quat* gamma, quatlat &I) {

    auto norm = I.norm().first/I.norm().second;
    // Gram computation
    NTL::mat_ZZ gram;
    gram.SetDims(4,4);
    if (I.alg.p%4==3) {

        if (I.basis[1][0] == 0) {
            // assert( I.basis[3][0] + I.basis[2][1] == 2*norm );
            I.basis[3][0] = - I.basis[2][1];

            // gram computation
            gram[0][0] = 2 * norm;
            gram[0][1] = Integer(0);
            gram[1][0] = gram[0][1];
            gram[0][2] = I.basis[2][0];
            gram[2][0] = gram[0][2];
            gram[0][3] = I.basis[3][0];
            gram[3][0] = gram[0][3];
            gram[1][1] = gram[0][0];
            // could do better here ?
            // gram[1][2] = I.basis[2][1];
            gram[1][2] = - gram[0][3];

            gram[2][1] = gram[1][2];
            gram[3][1] = gram[0][2];
            gram[1][3] = gram[3][1];

            // could do better here ?
            auto sq = I.alg.p + I.basis[2][0] * I.basis[2][0];
            gram[2][2] = (sq + I.basis[2][1] * I.basis[2][1])/(2*norm);

            // could do better here ?
            // gram[3][2] = (I .basis[2][0] * I.basis[3][0] + I.basis[2][1] * I.basis[3][1])/(2*norm);
            gram[3][2] = Integer(0);
            gram[2][3] = gram[3][2];

            // could do better here ?
            // gram[3][3] = (sq + I.basis[3][0] * I .basis[3][0])/(2*norm);
            gram[3][3] = gram[2][2];
        }
        else {
            //TODO
            std::array<quat, 4> elts {{
            {{I.basis[0][0], I.basis[0][1], I.basis[0][2], I.basis[0][3], NTL::ZZ(1)}, I.alg},
            {{I.basis[1][0], I.basis[1][1], I.basis[1][2], I.basis[1][3], NTL::ZZ(1)}, I.alg},
            {{I.basis[2][0], I.basis[2][1], I.basis[2][2], I.basis[2][3], NTL::ZZ(1)}, I.alg},
            {{I.basis[3][0], I.basis[3][1], I.basis[3][2], I.basis[3][3], NTL::ZZ(1)}, I.alg},
             }};
            // NTL::mat_ZZ gram_check;
            // gram_check.SetDims(4,4);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j) {
                    gram[i][j] = (elts[i] * elts[j].conjugate())[0] / (2*norm);
            }
            // std::cout << gram_check << "\n";
            // std::cout << gram << "\n";

            // assert(0);}
        }
    }
    else {
        assert(0);
    }
    mat_t Gram;
    mat_t Red;
    Gram.resize(4, 4);
    for (unsigned i = 0; i < 4; ++i) {
        for (unsigned j = 0; j < 4; ++j) {

                auto val = gram[i][j];
            // std::cout << val << " ";
            ntl2gmp(Gram[i][j].get_data(), val);
        }
    }
    Red.gen_identity(4);
    mat_t Redinv;
    gso_t gso(Gram, Red, Redinv, fplll::GSO_INT_GRAM);
    fplll::LLLReduction<fplll::Z_NR<mpz_t>, fplll::FP_NR<double>> lllobj(gso,
    fplll::LLL_DEF_DELTA,
    fplll::LLL_DEF_ETA,
    // 0.9999999,
    // 0.5000001,
    0);
    lllobj.lll();
    NTL::Vec<NTL::ZZ> row;
    row.SetLength(4);
    for (unsigned j = 0; j < 4; ++j) {
        NTL::ZZ c;
        gmp2ntl(c, Red[0][j].get_data());
        row += c * I.basis[j];
    }

    (*gamma)[0] = row[0];
    (*gamma)[1] = row[1];
    (*gamma)[2] = row[2];
    (*gamma)[3] = row[3];
    (*gamma)[4] = I.denom;
    // Integer small_norm;
    // gmp2ntl(small_norm, Gram[0][0].get_data());
    // std::cout << small_norm << "\n";
    assert(I.contains(*gamma));

    Integer small_norm;
    gmp2ntl(small_norm, Gram[0][0].get_data());
    small_norm = small_norm/2;

   if ((gamma->norm().first/gamma->norm().second) != small_norm * (I.norm().first/I.norm().second)) {
    std::cout << gamma->norm().first/gamma->norm().second << " " << small_norm * (I.norm().first/I.norm().second) << "\n";
}

    assert((gamma->norm().first/gamma->norm().second) == small_norm * (I.norm().first/I.norm().second));

    return small_norm;
}
