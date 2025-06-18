///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
////   Code implementing the necessary basic quaternion arithmetic
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include "quaternions.hpp"
#include "utils.hpp"
#include "endring.hpp"


quatlat::quatlat(NTL::mat_ZZ const &gens, NTL::ZZ const &denom_, quatalg const &alg_) : alg{alg_}, basis{gens}, denom{denom_}, gen{{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, alg_ }
{
    assert(basis.NumCols() == 4);
    normalize();
    the_norm.first = the_norm.second = 0;
}


quatlat::quatlat(NTL::mat_ZZ const &gens, std::pair<NTL::ZZ,NTL::ZZ> const &norm, NTL::ZZ const &denom_, quatalg const &alg_) : alg{alg_}, basis{gens}, denom{denom_}, gen{{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, alg_ }
{
    assert(basis.NumCols() == 4);
    the_norm = norm;
    normalize();
}

quatlat::quatlat(NTL::mat_ZZ const &gens, NTL::ZZ const &denom_, quatalg const &alg_, bool do_normalize) : alg{alg_}, basis{gens}, denom{denom_}, gen{{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, alg_ }
{
    assert(basis.NumCols() == 4);
    if (do_normalize) {
        normalize();
    }
    the_norm.first = the_norm.second = 0;
}


quatlat::quatlat(NTL::mat_ZZ const &gens, std::pair<NTL::ZZ,NTL::ZZ> const &norm, NTL::ZZ const &denom_, quatalg const &alg_, bool do_normalize) : alg{alg_}, basis{gens}, denom{denom_}, gen{{Integer(0), Integer(0), Integer(0), Integer(0), Integer(1)}, alg_ }
{
    // assert(basis.NumCols() == 4);
    the_norm = norm;
    if (do_normalize) {
        normalize();
    }

}

bool quatlat::is_order() const { auto [n,d] = norm(); assert(!NTL::IsZero(d)); return n == d; }

quatlat quatalg::maximal_order(bool const &surface) const
{
    NTL::mat_ZZ mat;
    NTL::ZZ denom;
    mat.SetDims(4, 4);
    NTL::ZZ q(this->q);
    NTL::ZZ two_q = q*2;

    if (q == 1) {
        denom = 2;
        //Special case with q == 1 (p % 4 == 3), where we need surface/floor things
        mat[0][0] = mat[1][1] = 2;
        mat[2][1^surface] = mat[2][2] = 1;
        mat[3][0^surface] = mat[3][3] = 1;
    }

    else if (q == 2) {
        denom = 4;
        mat[0][0] = 2; mat[0][1] = 0; mat[0][2] = 2; mat[0][3] = 2;
        mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 2; mat[1][3] = 1;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 4; mat[2][3] = 0;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 4;
    }

    else {
        NTL::ZZ p(this->p);
        assert (q % 4 == 3);
        NTL::ZZ a(0);
        while (!((a*a*p + 1) % q == 0)) {
            a += 1;
        }
        denom = q*2;
        mat[0][0] = q; mat[0][1] = q; mat[0][2] = 0; mat[0][3] = 0;
        mat[1][0] = 0; mat[1][1] = 2; mat[1][2] = 0; mat[1][3] = a*2;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = q; mat[2][3] = q;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = denom;
    }

    quatlat O(mat, denom, *this, false);
    assert (O.is_order());
    return O;
};

quatlat quatalg::maximal_order_with_quat_for_check(quat *alpha, bool const &surface) const
{
    NTL::mat_ZZ mat;
    NTL::ZZ denom;
    mat.SetDims(4, 4);
    NTL::ZZ q(this->q);
    NTL::ZZ two_q = q*2;

    if (q == 1) {
        denom = 2;
        //Special case with q == 1 (p % 4 == 3), where we need surface/floor stuff
        mat[0][0] = mat[1][1] = 2;
        mat[2][1^surface] = mat[2][2] = 1;
        mat[3][0^surface] = mat[3][3] = 1;
    }

    else if (q == 2) {
        denom = 4;
        mat[0][0] = 2; mat[0][1] = 0; mat[0][2] = 2; mat[0][3] = 2;
        mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 2; mat[1][3] = 1;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 4; mat[2][3] = 0;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 4;
    } else {
        NTL::ZZ p(this->p);
        assert (q % 4 == 3);
        NTL::ZZ a(0);
        while (!((a*a*p + 1) % q == 0)) {
            a += 1;
        }
        denom = q*2;
        mat[0][0] = q; mat[0][1] = q; mat[0][2] = 0; mat[0][3] = 0;
        mat[1][0] = 0; mat[1][1] = 2; mat[1][2] = 0; mat[1][3] = a*2;
        mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = q; mat[2][3] = q;
        mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = denom;

        (*alpha)[0] = NTL::ZZ(0);
        (*alpha)[1] = NTL::ZZ(1);
        (*alpha)[2] = NTL::ZZ(0);
        (*alpha)[3] = NTL::ZZ(a);
        (*alpha)[4] = NTL::ZZ(1);
    }

    quatlat O(mat, denom, *this, false);
    assert (O.is_order());
    return O;
};


void quatlat::normalize(bool safe) {
        // NTL::ZZ det;
        size_t numrows = basis.NumRows();
        // std::cout << basis << "\n" ;
        double delta = 0.5;
        // clock_t t =clock();
        // for (int i = 0; i < numrows; i++) {
        //     basis[i][2] *=32;
        //     basis[i][3] *=32;
        // }


        size_t rank;
        if ((!(safe)) && NTL::NumBits(basis[0][0]) < 240 && NTL::NumBits(basis[0][1]) < 240) {
            // std::cout << "fast_LLL \n";
            rank = NTL::LLL_FP(basis, delta, 0 , 0, 0);
        }
        else {
            // std::cout << "slow_LLL \n";
            rank = NTL::LLL_XD(basis, delta, 0 , 0, 0);
        }
        // std::cout << "LLL time " << (double) (clock() - t) / CLOCKS_PER_SEC << "\n";

        if (rank != numrows)
            for (size_t i = 0; i < rank; ++i)
                std::swap(basis[i], basis[numrows-rank+i]);
        basis.SetDims(rank, 4);

        assert(basis.NumRows() == 4);
        assert(basis.NumCols() == 4);

        // for (int i = 0; i < rank; i++) {
        //     basis[i][2] /=32;
        //     basis[i][3] /=32;
        // }

        assert(denom != 0);
        if (denom < 0) {
            basis = -basis;
            denom = -denom;
        }

        NTL::ZZ g = denom;
        for (size_t i = 0; i < (size_t) basis.NumRows(); ++i) {
            if (NTL::IsOne(g)) break;
            for (size_t j = 0; j < (size_t) basis.NumCols(); ++j) {
                if (NTL::IsOne(g)) break;
                g = NTL::GCD(g, basis[i][j]);
            }
        }

        assert(g > 0);
        if (!NTL::IsOne(g)) {
            denom /= g;
            for (size_t i = 0; i < (size_t) basis.NumRows(); ++i)
                for (size_t j = 0; j < (size_t) basis.NumCols(); ++j)
                    basis[i][j] /= g;
        }
        // basis = this->HNF_basis();

    }


quatlat quatlat::conjugate() const
       {
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4, 4);
        for (unsigned i = 0; i < 4; ++i) {
            newbasis[i][0] =  basis[i][0];
            newbasis[i][1] = -basis[i][1];
            newbasis[i][2] = -basis[i][2];
            newbasis[i][3] = -basis[i][3];
        }
        if(NTL::IsZero(the_norm.second)) {
            return {newbasis, denom, alg};
        }
        return {newbasis, the_norm, denom, alg, false};
    }

void quatlat::_conjugate()
{
    for (unsigned i = 0; i < 4; ++i) {
            basis[i][1] = -basis[i][1];
            basis[i][2] = -basis[i][2];
            basis[i][3] = -basis[i][3];
        }
}

quatlat quatlat::HNF_conjugate() const
{
   NTL::mat_ZZ newbasis;
        newbasis.SetDims(4, 4);
        for (unsigned i = 1; i < 4; ++i) {
            newbasis[i][0] = -basis[i][0];
            newbasis[i][1] = basis[i][1];
            newbasis[i][2] = basis[i][2];
            newbasis[i][3] = basis[i][3];
        }
        newbasis[0][0] = basis[0][0];
        if(NTL::IsZero(the_norm.second)) {
            return {newbasis, denom, alg};
        }
        return {newbasis, the_norm, denom, alg, false};
}

quatlat quatlat::inverse() const
{
    auto conj = conjugate();
    auto nrm = norm();
    return {conj.basis*nrm.second, {nrm.second, nrm.first}, denom*nrm.first, alg, false};
}

quatlat quatlat::operator+(quatlat const &other) const
{
        size_t n = this->basis.NumRows();
        size_t m = other.basis.NumRows();
        NTL::mat_ZZ newbasis = this->basis * other.denom;
        newbasis.SetDims(n+m, 4);
        for (size_t i = 0; i < m; ++i)
            newbasis[n+i] = other.basis[i] * this->denom;
        auto newdenom = this->denom * other.denom;

        return {newbasis, newdenom, alg};

}

quatlat quatlat::operator*(quatlat const &other) const
{
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(16, 4);
        quat a{alg}, b{alg};
        a[4] = this->denom;
        b[4] = other.denom;
        auto newdenom = this->denom * other.denom;
        for (size_t i = 0; i < 4; ++i) {
            for (size_t k = 0; k < 4; ++k)
                a[k] = this->basis[i][k];
            for (size_t j = 0; j < 4; ++j) {
                for (size_t k = 0; k < 4; ++k)
                    b[k] = other.basis[j][k];

                auto c = a * b;
                // std::cerr << "(" << a << ") * (" << b << ") = " << c << std::endl;
                for (unsigned k = 0; k < 4; ++k)
                    newbasis[4*i+j][k] = c[k] * newdenom / c[4];
            }
        }
        if (NTL::IsZero(the_norm.second) || NTL::IsZero(other.norm_no_comput().second)) {
            return {newbasis, newdenom, alg};
        }
        // if we deal with integral ideals we can improve
        auto nrm1 = norm();
        auto nrm2 = other.norm();
        return {newbasis, {nrm1.first * nrm2.first, nrm1.second * nrm2.second} , newdenom, alg};

}


quatlat quatlat::operator*(quat const &other) const
{
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4, 4);
        quat a{alg};
        a[4] = this->denom;
        auto newdenom = this->denom * other[4];
        for (size_t i = 0; i < 4; ++i) {
            for (size_t k = 0; k < 4; ++k) {
                a[k] = this->basis[i][k];

                auto c = a * other;
                for (unsigned k = 0; k < 4; ++k)
                    newbasis[i][k] = c[k] * newdenom / c[4];
            }
        }

        if (NTL::IsZero(the_norm.second)) {
            return {newbasis, newdenom, alg};
        }
        auto nrm1 = norm();
        auto nrm2 = other.norm();
        return {newbasis, {nrm1.first * nrm2.first, nrm1.second * nrm2.second} , newdenom, alg, false};
    }


quatlat quatlat::operator*(NTL::ZZ const &other) const
{
    NTL::mat_ZZ newbasis = basis;
    newbasis *= other;
    if (NTL::IsZero(the_norm.second)) {
        return {newbasis, denom, alg};
    }
    auto nrm = norm();
    return {newbasis, {nrm.first * other * other, nrm.second} ,denom, alg, false};
}

NTL::mat_ZZ quatlat::normalized_invariant() const {
    //quatlat this_copy = (*this);
    //this_copy.normalize(true);
    auto basis_mat = this->HNF_basis();
    for (unsigned i = 0; i < 4; ++i) {
        bool swap_sign = false; //first non-zero coeff should be positive
        for (unsigned j = 0; j < 4; ++j) {
            if (basis_mat[i][j] > 0) {
                break;
            }
            if (basis_mat[i][j] < 0) {
                swap_sign = true;
                break;
            }
        }
        if (swap_sign) {
            for (unsigned j = 0; j < 4; ++j) {
                basis_mat[i][j] = -basis_mat[i][j];
            }
        }
    }
    //std::cout << "INVARIANT COMPUTATION" << std::endl;
    //std::cout << basis_mat << std::endl;
    //std::cout << this_copy.denom << std::endl;
    return basis_mat;
}

bool quatlat::operator==(quatlat const &other) const {
    quatlat this_copy = (*this);
    quatlat other_copy = other;
    return &alg == &other.alg && this->normalized_invariant() == other.normalized_invariant();
}

quatlat quatlat::conjugated_by_j() const {
    NTL::mat_ZZ newbasis = basis;
    for (unsigned i = 0; i < 4; ++i) {
        newbasis[i][1] = -newbasis[i][1];
        newbasis[i][3] = -newbasis[i][3];
    }

    return {newbasis, denom, alg, false};
}

std::pair<NTL::ZZ,NTL::ZZ> quatlat::norm() const
 {
        if (NTL::IsZero(the_norm.second)) {
            auto L = (*this) * this->conjugate();
            NTL::vec_ZZ v, x;
            v.SetLength(4);
            NTL::set(v[0]);
            NTL::ZZ num;
            NTL::solve1(num, x, L.basis, v);
            auto g = NTL::GCD(num, L.denom);
            the_norm = std::make_pair(num/g, L.denom/g);
        }
        return the_norm;
}


void quatlat::_intersect(quatlat const &other, bool safe)
{
        auto bas1 = NTL::transpose(basis * other.denom);
        auto bas2 = NTL::transpose(other.basis * denom);
        NTL::ZZ det1, det2;
        NTL::mat_ZZ ker1, ker2;
        NTL::inv(det1, ker1, bas1);
        NTL::inv(det2, ker2, bas2);
        assert(ker1.NumRows() == 4 && ker2.NumRows() == 4);

        NTL::mat_ZZ ker;
        ker.SetDims(8, 4);
        for (unsigned j = 0, i = 0; j < 4; ++j) {
            for (i = 0; i < ker1.NumRows(); ++i)
                ker[i][j] = ker1[i][j] * det2;
            for (i = 0; i < ker2.NumRows(); ++i)
                ker[4+i][j] = ker2[i][j] * det1;
        }

        NTL::ZZ det;
        size_t rank = NTL::image(det, ker);
        assert(rank == 4);
        for (size_t i = 0; i < rank; ++i)
            std::swap(ker[i], ker[4+i]);
        ker.SetDims(4, 4);

        NTL::transpose(ker, ker);

        NTL::inv(det, basis, ker);
        denom *= other.denom * det;
        assert(denom % (det1 * det2) == 0);
        denom /= det1 * det2;


        normalize(safe);
        auto nrm1 = the_norm;
        auto nrm2 = other.norm_no_comput();
        if ( nrm1.second != 0 && nrm2.second !=0 && nrm1.first % nrm1.second ==0 && nrm2.first % nrm2.second == 0 && NTL::GCD(nrm1.first/nrm1.second, nrm2.first/nrm2.second) == 1 ) {
            the_norm.first = (nrm1.first * nrm2.first) / (nrm2.second * nrm1.second);
            the_norm.second = NTL::ZZ(1);
        }
        else {
            the_norm.first = the_norm.second = 0;
        }

}

void quatlat::_fast_intersect(quatlat const &other) {

    if (this->alg.p%4 == 3) {
        // requires both basis of ideals to be in HNF
        assert( basis[0][1] == 0 && basis[2][3] == 0 && other.basis[0][1] == 0 && other.basis[2][3] == 0 );
        assert( denom == Integer(2) && other.denom == 2);
        assert( other.norm().second == 1 && this->the_norm.second == 1);
        auto n1 = other.norm().first;
        auto n = this->the_norm.first;
        auto newnorm = n * n1;
        // std::cout << *this << "\n";
        // std::cout << other << "\n";

        if (n1 > 1) {

            basis[0][0] *= n1;
            if (basis[1][0] == 0 && other.basis[1][0] == 0) {
                basis[1][1] *= n1;
                auto before = basis[2][0];
                assert(basis[2][0] >= 0);
                assert(other.basis[2][0] >=0);
                NTL::CRT(basis[2][0], n, other.basis[2][0], 2*n1);
                n = this->the_norm.first;

                assert(basis[2][0] % (2*n1) == other.basis[2][0] % (2*n1));

                NTL::CRT(basis[2][1], n, other.basis[2][1], 2*n1);
                assert(basis[2][1] % (2*n1) == other.basis[2][1]);
                basis[3][1] = basis[2][0];
                basis[3][0] = - basis[2][1];
            }
            else {
                NTL::CRT(basis[1][0], n, other.basis[1][0], 2*n1);
                n = this->the_norm.first;
                NTL::CRT(basis[1][1], n, other.basis[1][1], 2*n1);
                n = this->the_norm.first;

                auto mod = newnorm;

                // reducing
                auto a = NTL::GCD(mod,basis[1][1]);
                auto b = basis[1][1]/(2*a);
                basis[1][1] = 2*a;
                while (NTL::GCD(b,mod) != 1) {
                    b += mod/a;
                }
                assert(NTL::GCD(b,mod) == 1);
                auto c = NTL::InvMod(b % mod,mod);


                basis[1][0] = 2*((basis[1][0]/2 * c) % mod);

                NTL::CRT(basis[2][0], n, other.basis[2][0], 2*n1);
                n = this->the_norm.first;
                NTL::CRT(basis[2][1], n, other.basis[2][1], 2*n1);
                n = this->the_norm.first;
                NTL::CRT(basis[2][2], n, other.basis[2][2], 2*n1);
                n = this->the_norm.first;

                // reducing
                basis[2][1] = (basis[2][1] - basis[2][2])/2;
                a = NTL::GCD(mod,basis[2][2]);
                b = basis[2][2]/(a);
                while (NTL::GCD(b,mod) != 1) {
                    b += 2*mod/a;
                }
                basis[2][2] = a;
                c = NTL::InvMod(b % mod, mod);
                basis[2][0] = 2*((basis[2][0]/2 * c) % mod);
                basis[2][1] = basis[2][2] + 2*(((basis[2][1]) * c) % mod);

                // and reducing by previous vector
                Integer quotient = (basis[2][1]/basis[1][1] );
                basis[2][1] = basis[2][1] - basis[1][1]*quotient;
                basis[2][0] = (basis[2][0] - basis[1][0]*quotient) % (2*mod);

                NTL::CRT(basis[3][0], n, other.basis[3][0], 2*n1);
                n = this->the_norm.first;
                assert(basis[3][0] % (2*n1) == other.basis[3][0]);
                // std::cout << basis[3][1] << " " << n << " " << other.basis[3][1] + 2*n1 << " " << 2*n1 << "\n";
                NTL::CRT(basis[3][1], n, (other.basis[3][1] + 2*n1), 2*n1);
                // std::cout << basis[3][1] << " " << n << "\n";
                n = this->the_norm.first;
                assert(basis[3][1] % (2*n1) == other.basis[3][1] %(2*n1));
                NTL::CRT(basis[3][2], n, other.basis[3][2], 2*n1);
                n = this->the_norm.first;
                NTL::CRT(basis[3][3], n, other.basis[3][3], 2*n1);
                // n = this->the_norm.first;

                // std::cout << basis << "\n";

                // reducing
                basis[3][0] = (basis[3][0] - basis[3][3])/2;
                a = NTL::GCD(mod,basis[3][3]);
                b = basis[3][3]/(a);
                basis[3][3] = a;
                c = NTL::InvMod(b % mod, mod);
                basis[3][0] = basis[3][3] + 2*((basis[3][0] * c) % mod);
                basis[3][1] = 2*((((basis[3][1] -basis[3][2])/2) * c) % mod);
                basis[3][2] = (((basis[3][2]) * c) % mod);
                basis[3][1] += basis[3][2];

                // and reducing by previous vector
                quotient = (basis[3][2]/basis[2][2] );
                basis[3][2] = basis[3][2] - basis[2][2]*quotient;
                basis[3][1] = (basis[3][1] - basis[2][1]*quotient) % (2*mod);
                basis[3][0] = (basis[3][0] - basis[2][0]*quotient) % (2*mod);
                // and by previous vector
                quotient = (basis[3][1]/basis[1][1] );
                basis[3][1] = (basis[3][1] - basis[1][1]*quotient) % (2*mod);
                basis[3][0] = (basis[3][0] - basis[1][0]*quotient) % (2*mod);



                // std::cout << "final intersect basis = \n" <<  basis << "\n";
            }

            if (basis[1][0] < 0) {
                basis[1][0] += 2 * newnorm;
            }
            if (basis[2][0] < 0) {
                basis[2][0] += 2 * newnorm;
                }
            if (basis[2][1] < 0) {
                basis[2][1] += 2 * newnorm;
            }
            if (basis[3][0] < 0) {
                basis[3][0] += 2 * newnorm;
                }
            if (basis[3][1] < 0) {
                basis[3][1] += 2 * newnorm;
            }
            if (basis[3][2] < 0) {
                basis[3][2] += 2 * newnorm;
            }
        }
        // else we have nothing to do
        this->the_norm = {newnorm, Integer(1)};
    }
    else
    {
        this->_intersect(other);
    }
    // std::cout << basis << "\n";


}

quatlat quatlat::_compute_order(bool right_order=false) const
{
        reset_norm();
        auto nrm = norm();

        quatlat L(basis*nrm.second, denom*nrm.first, alg);

        for (unsigned k = 0; k < 4; ++k) {

            quat a{alg};
            for (unsigned j = 0; j < 4; ++j)
                a[j] = basis[k][j];
            a[4] = denom;
            a.invert();

            NTL::mat_ZZ mat;
            mat.SetDims(4, 4);
            NTL::ZZ d(1);
            for (unsigned i = 0; i < 4; ++i) {
                quat b{alg};
                for (unsigned j = 0; j < 4; ++j)
                    b[j] = basis[i][j];
                b[4] = denom;
                auto c = right_order ? a * b : b * a;
                mat *= c[4];
                for (unsigned j = 0; j < 4; ++j)
                    mat[i][j] = c[j] * d;
                d *= c[4];
            }

            quatlat N(mat, d, alg);


            // std::cerr << "N:\n" << N << "\n";
            L._intersect(N, true);
            // std::cerr << "L:\n" << L << "\n";
        }

        L.reset_norm();
        assert(L.is_order());
        return L;
}

quatlat quatlat::new_right_order() const {

        // compute the norm
        auto nrm = norm();
        // this will be the result
        // quatlat L(basis, nrm ,denom, alg);

        // quatlat L = inverse();

        // now we multiply
        // L = L * (*this);
        // return L;
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(6, 4);

        quat a{alg}, b{alg};
        a[4] = this->denom*nrm.first;
        b[4] = this->denom;
        auto newdenom = this->denom * a[4];
        // newbasis[6][0] = newdenom;
        // newbasis[6][1] = 0;
        // newbasis[6][2] = 0;
        // newbasis[6][3] = 0;

        int count = 0;
        for (size_t i = 0; i < 3; ++i) {
            a[0] = this->basis[i][0];
            for (size_t k = 1; k < 4; ++k)
                a[k] = -this->basis[i][k];
            for (size_t j = i+1; j < 4; ++j) {
                for (size_t k = 0; k < 4; ++k)
                    b[k] = this->basis[j][k];

                auto c = a * b;

                //std::cerr << "(" << a << ") * (" << b << ") = " << c << std::endl;
                for (unsigned k = 0; k < 4; ++k)
                    newbasis[count][k] = c[k] * newdenom * nrm.second / c[4];

                count ++;

            }
        }

         // normalizing
         double delta = 0.5;
         // NTL::LLLCheckFct check = 0;
         size_t rank = NTL::LLL_FP(newbasis, delta, 0 , 0, 0);
         size_t numrows = 6;
         assert(rank==3 || rank == 4);

        // for (int i = 0; i < numrows; i++) {
        //     newbasis[i][2] *=32;
        //     newbasis[i][3] *=32;
        // }

         if (rank != numrows)
             for (size_t i = 0; i < rank; ++i)
                 std::swap(newbasis[i+(4 - rank)], newbasis[numrows-rank+i]);

         newbasis.SetDims(4, 4);
         assert(newbasis.NumRows() == 4);
         assert(newbasis.NumCols() == 4);



         if (rank == 3) {
             newbasis[0][0] = newdenom;
             newbasis[0][1] = 0;
             newbasis[0][2] = 0;
             newbasis[0][3] = 0;
         }

        assert(denom != 0);
        if (denom < 0) {
            newbasis = -newbasis;
            newdenom = -newdenom;
        }

         NTL::ZZ g = newdenom;
        for (size_t i = 0; i < (size_t) newbasis.NumRows(); ++i) {
            if (NTL::IsOne(g)) break;
            for (size_t j = 0; j < (size_t) newbasis.NumCols(); ++j) {
                if (NTL::IsOne(g)) break;
                g = NTL::GCD(g, newbasis[i][j]);
            }
        }

         assert(g > 0);
         if (!NTL::IsOne(g)) {
             newdenom /= g;
             for (size_t i = 0; i < (size_t) newbasis.NumRows(); ++i)
                 for (size_t j = 0; j < (size_t) newbasis.NumCols(); ++j)
                     newbasis[i][j] /= g;
        }

        return { newbasis, {NTL::ZZ(1), NTL::ZZ(1)}, newdenom, alg , false};

}

std::pair<quatlat,NTL::mat_ZZ> quatlat::fast_right_order_and_gram() {

    NTL::mat_ZZ gram;
    gram.SetDims(3,3);
    if (alg.p %4 == 3) {
        Integer norm = the_norm.first/the_norm.second;
        // testing if we are in some weird case
        Integer a = basis[3][0] * basis[3][0] + basis[3][1] * basis[3][1] - alg.p * ( basis[3][2] * basis[3][2] + basis[3][3] * basis[3][3]);
        assert(a%2 ==0);
        a = a/2;
        if (NTL::GCD(a , norm) != 1) {
            // std::cout << "weird case!!!!! \n";
            // return { { basis, the_norm, denom, alg, false }, gram };
            quatlat O = this->new_right_order();
            return { O, gram};
        }

        Integer normsqr = norm*norm;
        Integer newdenom = 2*norm;
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4,4);
        // first element is 1
        newbasis[0][0] = newdenom;
        newbasis[0][1] = Integer(0);
        newbasis[0][2] = Integer(0);
        newbasis[0][3] = Integer(0);

        Integer c = basis[3][1] * basis[2][2] - basis[2][1] * basis[3][2];

        // third element is of the form A * j + x * k
        // where A = GCD(c, norm)
        Integer A = NTL::GCD(c,norm);
        newbasis[2][0] = Integer(0);
        newbasis[2][1] = Integer(0);
        newbasis[2][2] = newdenom * A;

        // last element is (1 + (N/A) k)/2
        newbasis[3][0] = norm;
        newbasis[3][1] = Integer(0);
        newbasis[3][2] = Integer(0);
        newbasis[3][3] = normsqr/A;

        // computing the value of x
        Integer x;
        if (c == 0) {
            x = 0;
        }
        else {

            c = c / A;
            auto old_c = c;
            while (NTL::GCD(c,norm) != 1) {
                basis[2][0] += 2*norm;
                basis[3][1] += 2*norm;
                c = basis[3][1] * basis[2][2] - basis[2][1] * basis[3][2];
                assert(A == NTL::GCD(c,norm));
                c = c/A;
                a = basis[3][0] * basis[3][0] + basis[3][1] * basis[3][1] - alg.p * ( basis[3][2] * basis[3][2] + basis[3][3] * basis[3][3]);
                assert(a%2 ==0);
                a = a/2;
            }
            if (NTL::GCD(c, norm) != 1) {
                std::cout << basis << "\n";
                std::cout << old_c << " " << c << " " << A << " " << norm << "\n";
            }
            assert( NTL::GCD(c, norm) ==1 );
            c = NTL::InvMod(c % norm, norm);
            x =  (((basis[3][1] * basis[2][3] - basis[2][1] * basis[3][3]) % norm) * c) % norm ;
        }
        newbasis[2][3] = x * newdenom;


        // 2nd element is of the form (i + y * j + z * k) / newdenom
        newbasis[1][0] = Integer(0);
        newbasis[1][1] = Integer(1);
        // let gen be the last basis element
        // a priori this is derived from Conjugate(gen) * i * gen / norm = (a * i + b * j + c * k) / newdenom
        // first we compute b,c (a was computed previously)
        Integer b = basis[3][1] * basis[3][2] - basis[3][0] * basis[3][3];
        c = basis[3][1] * basis[3][3] - basis[3][0] * basis[3][2] ;

        // and new we reduce it
        assert(NTL::GCD(a , norm) == 1);
        Integer d = NTL::InvMod(a % normsqr,normsqr);
        newbasis[1][3] = ((c * d) % (2*normsqr));
        Integer is_odd = Integer(1);
        if (a % 2 == 1 && d %2 == 1) {
            is_odd = Integer(0);
        }

        newbasis[1][2] = ((b * d + normsqr*is_odd  ) % (2*normsqr));

        // computing the gram matrix
        gram[0][0] = (Integer(1) + this->alg.p * ( newbasis[1][2] * newbasis[1][2] +  newbasis[1][3] * newbasis[1][3]))/normsqr;
        gram[1][0] = Integer(2) * this->alg.p * (newbasis[1][2] * A +  newbasis[1][3] * x)/norm;
        gram[0][1] = gram[1][0];
        gram[0][2] = this->alg.p * newbasis[1][3] / A;
        gram[2][0] = gram[0][2];
        auto quot = norm / A;
        gram[1][2] = 2 * this->alg.p * quot * x;
        gram[2][1] = gram[1][2];
        gram[1][1] = 4 * this->alg.p * ( A*A + x*x );
        gram[2][2] = this->alg.p * quot * quot;


        return {{newbasis, {Integer(1), Integer(1)}, newdenom, alg, false },gram};
    }
    else {

        return {this->new_right_order(), gram};
    }




}

NTL::mat_ZZ quatlat::HNF_basis() const
{
        NTL::mat_ZZ M = this->basis;
        NTL::mat_ZZ H;
        NTL::HNF(H, M, NTL::determinant(M));
        return H;
}

NTL::mat_ZZ quatlat::LLL_basis() const
{
        NTL::mat_ZZ M = this->basis;
        NTL::ZZ scal[4];
        {
            NTL::ZZ s(1);
            s <<= NTL::NumBits(alg.p);  //TODO: should this be smaller or bigger?
            NTL::RR big = NTL::conv<NTL::RR>(s);
            scal[0] = NTL::CeilToZZ(NTL::conv<NTL::RR>(big));
            scal[1] = NTL::CeilToZZ(NTL::conv<NTL::RR>(big) * sqrt(NTL::conv<NTL::RR>(alg.q)));
            scal[2] = NTL::CeilToZZ(NTL::conv<NTL::RR>(big) * sqrt(NTL::conv<NTL::RR>(alg.p)));
            scal[3] = NTL::CeilToZZ(NTL::conv<NTL::RR>(big) * sqrt(NTL::conv<NTL::RR>(alg.p * alg.q)));
        }
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                M[i][j] *= scal[j];
        NTL::ZZ det;
        NTL::LLL(det, M);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j) {
                assert(M[i][j] % scal[j] == 0);
                M[i][j] /= scal[j];
            }
        return M;
}

void quatlat::enumerate_shortish(unsigned long bnd, std::function<bool(quat const &)> const &fun) const
{
    std::vector<quat> basis;
    {
        NTL::mat_ZZ mat = LLL_basis();
        for (unsigned i = 0; i < 4; ++i)
            basis.emplace_back(std::array<NTL::ZZ,5>{mat[i][0], mat[i][1], mat[i][2], mat[i][3], denom}, alg);
    }

    /*
     * TODO this enumeration can be done better; see for instance here:
     * https://github.com/sagemath/sage/blob/10.5/src/sage/modules/free_module.py#L2479-L2505
     */
    for (long max = 1; (unsigned long) max <= bnd; ++max) {
        long vec[4];
        for (vec[0] = 0; vec[0] <= max; ++vec[0])
            for (vec[1] = vec[0] ? -max : 0; vec[1] <= max; ++vec[1])
                for (vec[2] = vec[0] || vec[1] ? -max : 0; vec[2] <= max; ++vec[2])
                    for (vec[3] = -max; vec[3] <= max; ) {

                        quat const elt = basis[0]*vec[0] + basis[1]*vec[1] + basis[2]*vec[2] + basis[3]*vec[3];
                        if (fun(elt)) {
                        // std::cerr << "FOUND: " << max << " ~> " << vec[0] << " " << vec[1] << " " << vec[2] << " " << vec[3] << std::endl;
                            return;
                        }

                        if (abs(vec[0]) == max || abs(vec[1]) == max || abs(vec[2]) == max)
                            ++vec[3];
                        else if (vec[3] == -max)
                            vec[3] = +max;
                        else
                            break;
                    }
        }
}

quat quatlat::compute_gen_cyclic() const {
    this->reset_norm();
    assert (this->norm().second == 1);
    auto N = this->norm().first;
    quat alpha(alg);
    this->enumerate_shortish(100, //TODO <- good value here?
        [&](quat const &el) -> bool {
        assert (el.norm().second == 1);
        NTL::ZZ N = el.norm().first;
        if ((NTL::GCD(N*N, el.norm().first) == N)) { //Test that el is a generator.
            alpha[0] = el[0];
            alpha[1] = el[1];
            alpha[2] = el[2];
            alpha[3] = el[3];
            alpha[4] = el[4];
            return true;
        }
        return false;
    });
    return alpha;
}

quatlat quatlat::copy() const {
        return {basis, the_norm, denom, alg, false};
}


bool quatlat::contains(quat const &alpha) const
{
        NTL::mat_ZZ M = basis;
        NTL::ZZ du = denom;
        NTL::ZZ remain,det;
        NTL::vec_ZZ alpha_vec,solve_check;
        alpha_vec.SetLength(4);
        bool res = true;

        for (int i=0; i<4; i++) {
            alpha_vec[i] = denom* alpha[i];
            NTL::DivRem(alpha_vec[i],  remain, alpha_vec[i], alpha[4]);
            res = res && (remain==0);
        }
        if (res) {
            solve1(det, solve_check, M, alpha_vec);
            res = res && det == 1;
            assert(!res || solve_check * M == alpha_vec);
        }

        return res;

}

bool quatlat::fast_contains(quat const &alpha, Integer const &norm) const
{
    auto modnorm = alpha[4] * norm;
    // if we are in usual HNF format
    if (basis[1][0] == 0 && (basis[2][2] == 1 || basis[2][2] == 1) && basis[3][3] == 1) {
        if (alpha[4] != 2 && alpha[4] != 1) {
            return false;
        }
        assert(alpha[4] == 2 || alpha[4] == 1);
        return (alpha[0] % modnorm == ((alpha[2] * basis[2][0] + alpha[3] * basis[3][0])) % modnorm)
        && (alpha[1] % modnorm == ((alpha[2] * basis[2][1] + alpha[3] * basis[3][1] )) % modnorm);
    }
    else {
        return this->contains(alpha);
    }
}


void quatlat::reduce_norm_cyclic(const quatlat &left_order, const NTL::ZZ &norm)
{
        quatlat J = left_order * norm;
        *this = this->special_add(J);
        the_norm.first = norm;
        the_norm.second = NTL::ZZ(1);
}

quatlat quatlat::special_add(const quatlat &other) const
{
        size_t n = this->basis.NumRows();
        size_t m = other.basis.NumRows();
        NTL::mat_ZZ newbasis = this->basis * other.denom;
        newbasis.SetDims(n+m, 4);
        for (size_t i = 0; i < m; ++i)
            newbasis[n+i] = other.basis[i] * this->denom;
        auto newdenom = this->denom * other.denom;

        // if we deal with integral ideals we can improve
        if (NTL::IsZero(the_norm.second) || NTL::IsZero(other.norm_no_comput().second)) {
            return {newbasis, newdenom, alg};
        }
        else {
            auto nrm1 = norm();
            auto nrm2 = other.norm();
            if ( (nrm1.first % nrm1.second ==0) && ((nrm2.first % nrm2.second ==0))) {
                return {newbasis, {NTL::GCD(nrm1.first/nrm1.second, nrm2.first/nrm2.second), NTL::ZZ(1)}, newdenom, alg};
            }
        }

        return {newbasis, newdenom, alg};
}

void quatlat::set_generator(quat const &new_gen) {
    assert(new_gen.alg.p == this->alg.p && new_gen.alg.q == this->alg.q);
    this->gen[0] = new_gen[0];
    this->gen[1] = new_gen[1];
    this->gen[2] = new_gen[2];
    this->gen[3] = new_gen[3];
    this->gen[4] = new_gen[4];
}

quat quatlat::get_generator() const {
    if (!this->gen.is_zero()) {
        return this->gen;
    }
    else {
        // it could be computed
        return this->gen;
    }
}


quatlat create_from_generator(const quat &gen, const NTL::ZZ &norm, const quatlat &left_order) {

    clock_t t=clock(); (void) t;
    quatlat I = left_order * gen;
    // std::cout << "O * gen time " << (double) (clock() - t)/CLOCKS_PER_SEC << "\n";
    I.reduce_norm_cyclic(left_order, norm);
    I.set_generator(gen);
    return I;
}

quatlat create_from_generator_O0(const quat &gen, const NTL::ZZ &norm) {

    if (gen.alg.p%4!=3)
    {
        assert(0);
    }
    else {
        Integer newdenom = Integer(2);
        NTL::mat_ZZ newbasis;
        newbasis.SetDims(4,4);
        newbasis[0][0] = 2*norm;
        newbasis[0][1] = Integer(0);
        newbasis[0][2] = Integer(0);
        newbasis[0][3] = Integer(0);
        newbasis[1][2] = Integer(0);
        newbasis[1][3] = Integer(0);
        Integer c = (gen[2]*gen[2] + gen[3]*gen[3])%norm;
        if (NTL::GCD(norm, c) == 1) {
            newbasis[1][0] = Integer(0);
            newbasis[1][1] = newbasis[0][0];
            newbasis[2][2] = Integer(1);
            newbasis[2][3] = Integer(0);
            newbasis[3][3] = Integer(1);
            newbasis[3][2] = Integer(0);
            Integer Inv = NTL::InvMod((2*c)%norm, norm);
            newbasis[2][0] = 2*((2*(gen[2] * gen[0] + gen[1] * gen[3]) * Inv * NTL::InvMod(Integer(2),norm)) % norm);
            if (newbasis[2][0] < 0) {
                newbasis[2][0] += 2*norm;
            }
            newbasis[3][1] = newbasis[2][0];
            newbasis[2][1] = 1 + 2*(((gen[2] * gen[1] - gen[0] * gen[3] - c) * Inv ) % norm);
            if (newbasis[2][1] < 0) {
                newbasis[2][1] += 2*norm;
            }
            newbasis[3][0] = 2*norm - newbasis[2][1];
            if (newbasis[3][0] < 0) {
                newbasis[3][0] += 2*norm;
            }


            return {newbasis, {norm, Integer(1)}, newdenom, gen.alg, false};
        }
        else
        {

            // we are in one of the two split ideals
            // if the norm is not a prime we should be within the order to jinv enumeration and we are going to skip
            // this ideal anyway
            if (!NTL::ProbPrime(norm)) {
                newbasis[1][0] = 1;
                return {newbasis, {norm, Integer(1)}, newdenom, gen.alg, false};
            }

            // we compute one quaternion element in ZZ[omega] of norm divisible by norm
            // std::cout << norm << " " << norm  - gen.alg.q << "split ideal... \n";
            assert(NTL::Jacobi(Integer(norm  - gen.alg.q),norm));
            Integer a;

            // weird cases where NTL modular sqrt just times out
            if (norm == 25) {
                // weird failure case
                a = Integer(7);
            }
            else if (norm == 169) {
                a = Integer(99);
            }
            else  {
                a = NTL::SqrRootMod(Integer(norm  - gen.alg.q),norm);
            }
            // auto bb = NTL::SqrRootMod(Integer(48), Integer(49));
            // else {}

            // std::cout << "out of sqrt \n";
            quat special_gen = quat({ {a, Integer(1), Integer(0), Integer(0), Integer(1)}, gen.alg });
            quat prod = (gen)*(special_gen);
            // if the result is on NO0 then the generator is the conjugate of special_gen
            if ( prod[0] % norm == 0 && prod[1] % norm == 0 && prod[2] % norm == 0, prod[3] % norm == 0 ) {
                newbasis[1][0] = 2*norm - 2*a;
                newbasis[1][1] = Integer(2);
                newbasis[2][0] = Integer(0);
                newbasis[2][1] = norm;
                newbasis[2][2] = norm;
                newbasis[2][3] = Integer(0);
                newbasis[3][0] = Integer(1);
                newbasis[3][1] = a;
                newbasis[3][2] = a;
                newbasis[3][3] = Integer(1);
                return {newbasis, {norm, Integer(1)}, newdenom, gen.alg, false};
            }
            else {
                newbasis[1][0] = 2*a;
                newbasis[1][1] = Integer(2);
                newbasis[2][0] = Integer(0);
                newbasis[2][1] = norm;
                newbasis[2][2] = norm;
                newbasis[2][3] = Integer(0);
                newbasis[3][0] = Integer(1);
                newbasis[3][1] = norm - a;
                newbasis[3][2] = norm - a;
                newbasis[3][3] = Integer(1);
                return {newbasis, {norm, Integer(1)}, newdenom, gen.alg, false};
            }
        }
    }

    auto O0 = starting_curve(gen.alg, false).second;
    return create_from_generator(gen, norm, O0);
    // return *

}


quatlat connecting_ideal(const quatlat &left_order, const quatlat &right_order) {
    // WARNING: there may be uncaught bugs in this function
    auto I_conn_frac = left_order*right_order;
    I_conn_frac.reset_norm();
    return I_conn_frac.norm().second * I_conn_frac;
}

void quatlat::right_ideals_of_norm(NTL::ZZ const &ell, std::function<void(quatlat const &)> const &fun)
{
        if (!is_order())
            throw std::logic_error("not an order");

        NTL::ZZ_pPush push(ell);

        auto mat1 = NTL::ident_mat_ZZ_p(2);
        auto [mati, matj] = alg.splitting(ell);
        auto matk = mati * matj;
        NTL::mat_ZZ_p const *mats[4] = {&mat1, &mati, &matj, &matk};

        // can ignore denominator for now because it's all projective anyway

        NTL::mat_ZZ mat;
        mat.SetDims(4, 4);
        for (unsigned k = 0; k < 4; ++k) {
            NTL::mat_ZZ cur;
            cur.SetDims(2,2);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                    cur[j/2][j%2] += basis[k][i] * NTL::rep((*mats[i])[j/2][j%2]);
                    // std::cerr << cur << std::endl;
            for (unsigned j = 0; j < 4; ++j)
                mat[k][j] = cur[j/2][j%2];
                // std::cerr << mat[k] << std::endl;
        }

        // We do need to be careful when the denominator is divisible by ell
        //     currently it fails in this case; TODO: fix

        NTL::mat_ZZ_p matmod, matinv;
        matmod.SetDims(4,4);
        for (unsigned i = 0; i < 4; ++i)
            for (unsigned j = 0; j < 4; ++j)
                NTL::conv(matmod[i][j], mat[i][j]);
        NTL::inv(matinv, matmod);

        NTL::vec_ZZ_p rhs;
        rhs.SetLength(4);
        NTL::ZZ_p &x = rhs[0], &y = rhs[2];

        NTL::set(x);
        for (NTL::ZZ i; i <= ell; ++i) {

            auto sol = rhs * matinv;
            assert(sol * matmod == rhs);

            quat elt{alg};
            NTL::set(elt[4]);
            for (unsigned i = 0; i < 4; ++i)
                for (unsigned j = 0; j < 4; ++j)
                    elt[j] += NTL::rep(sol[i]) * basis[i][j];
            assert(NTL::IsOne(elt.norm().second));
            assert(NTL::IsZero(elt.norm().first % ell));

            // std::cerr << rhs << std::endl;

            auto I = ell*(*this) + elt*(*this);
            assert(I.norm().first == ell && NTL::IsOne(I.norm().second));
            fun(I);

            if (NTL::IsZero(i)) {
                NTL::clear(x);
                NTL::set(y);
            }
            else
                ++x;
        }
}


std::list<quatlat> quatlat::left_ideals_of_prime_norm(NTL::ZZ const &ell, const quat &first_gen, const quat &iter) {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Enumerates through the set of ell-ideals as ideals of the form * first_gen * (C + D *iter) + this * ell
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    assert(this->is_order());
    assert(NTL::ProbPrime(ell));

    std::list<std::tuple<NTL::ZZ,NTL::ZZ>> coeff_list = {};
    std::list<quatlat> id_list = {};

    // create the list of coeff
    NTL::ZZ iterate = NTL::ZZ(0);
    while (iterate < ell) {
        coeff_list.push_back({NTL::ZZ(1),iterate});
        if (NTL::GCD(iterate,NTL::ZZ(ell))!=1) {
            coeff_list.push_back({iterate,NTL::ZZ(1)});
        }
        iterate++;
    }

    for (auto coeff : coeff_list) {
        quat gen = (first_gen * (quat({{std::get<0>(coeff),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(1)},this->alg}) + iter * std::get<1>(coeff)  ));
        quatlat J = create_from_generator( gen, ell, *this);
        id_list.push_back(J);
    }

    return id_list;
}

// enumerate throught the set of ell-ideals as ideals of the form * first_gen * (C + D *iter) + this * ell
std::list<quatlat> left_ideals_of_prime_norm_O0(NTL::ZZ const &ell, const quat &first_gen, const quat &iter) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Enumerates through the set of ell-ideals of O_0 as ideals of the form * first_gen * (C + D *iter) + this * ell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    assert(NTL::ProbPrime(ell));

    std::list<std::tuple<NTL::ZZ,NTL::ZZ>> coeff_list = {};
    std::list<quatlat> id_list = {};

    // create the list of coeff
    NTL::ZZ iterate = NTL::ZZ(0);
    while (iterate < ell) {
        coeff_list.push_back({NTL::ZZ(1),iterate});
        if (NTL::GCD(iterate,NTL::ZZ(ell))!=1) {
            coeff_list.push_back({iterate,NTL::ZZ(1)});
        }
        iterate++;
    }

    for (auto coeff : coeff_list) {
        quat gen = (first_gen * (quat({{std::get<0>(coeff),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(0),NTL::ZZ(1)},first_gen.alg}) + iter * std::get<1>(coeff)  ));
        quatlat J = create_from_generator_O0( gen, ell);
        id_list.push_back(J);
    }

    return id_list;
}


