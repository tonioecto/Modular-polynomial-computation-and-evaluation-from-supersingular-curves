#include "smallint.hpp"

static const std::vector<SmallInteger> invmod16 = {1, 11, 13, 7, 9, 3, 5, 15};

SmallInteger InvModSpecial(SmallInteger &a, const SmallInteger &m) {

    if (m == 16) {
        return invmod16[NTL::conv<long>(a)/2];
    }

    return NTL::conv<SmallInteger>(1/NTL::conv<Fp>(a));
}


SmallMatFp SmallMatFp::operator+(const SmallMatFp &other) const {
    assert (this->mod == other.mod);
    std::array<std::array<SmallInteger,2>,2> matrix{{ {(this->mat[0][0] + other.mat[0][0])%this->mod, (this->mat[0][1] + other.mat[0][1])%this->mod}, { (this->mat[1][0] + other.mat[1][0])%this->mod, (this->mat[1][1] + other.mat[1][1])%this->mod}  }};
    SmallMatFp M(this->mod, matrix);
    return M;
}

SmallMatFp SmallMatFp::operator*(const SmallInteger &other) const {
    return {this->mod, {{ { (this->mat[0][0] * other)%this->mod, (this->mat[0][1] * other)%this->mod}, { (this->mat[1][0] * other)%this->mod, (this->mat[1][1] * other)%this->mod}  }}};
}

SmallIntegerPair SmallMatFp::operator*(const SmallIntegerPair &x) const {
    return { ((this->mat[0][0] * x.first)%this->mod + (this->mat[1][0] * x.second)%this->mod)%this->mod, ((this->mat[0][1] * x.first)%this->mod + (this->mat[1][1] * x.second)%this->mod) % this->mod };
}

void SmallMatFp::normalize() {
    if(this->mat[0][0] < 0) {
        this->mat[0][0]+=this->mod;
    }
    if(this->mat[0][1] < 0) {
        this->mat[0][1]+=this->mod;
    }
    if(this->mat[1][0] < 0) {
        this->mat[1][0]+=this->mod;
    }
    if(this->mat[1][1] < 0) {
        this->mat[1][1]+=this->mod;
    }
}

std::pair<SmallInteger,SmallInteger> bi_dlp_3(const ecp &P, const std::pair<ecp,ecp> &bas) {
    ecp P3 = P;
    ecp R3 = bas.first;
    ecp S3 = bas.second;
    SmallInteger a1 = 0;
    SmallInteger a2 = 0;

    if (P3 == R3) {
        a1 = 1; a2 = 0;
    }
    else if (P3 == -R3) {
        a1 = 2; a2 = 0;
    }
    else if (P3 == S3) {
        a1 = 0; a2 = 1;
    }
    else if (P3 == -S3) {
        a1 = 0; a2 = 2;
    }
    else if (P3 == S3+R3) {
        a1 = 1; a2 = 1;
    }else if (P3 == S3-R3) {
        a1 = 2; a2 = 1;
    }
    else if (P3 == -S3+R3) {
        a1 = 1; a2 = 2;
    }else if (P3 == -S3-R3) {
        a1 = 2; a2 = 2;
    }

    assert( P3 == a1* R3 + a2*S3 );
    return {a1,a2};

}

std::pair<SmallInteger,SmallInteger> bi_dlp_16(const ecp &P, const std::pair<ecp,ecp> &bas) {
    SmallInteger b1 = 0;
    SmallInteger b2 = 0;
    ecp Ptemp = P;
    ecp R16 = bas.first;
    ecp S16 = bas.second;
    ecp R2 = 8*R16;
    ecp S2 = 8*S16;
    for (int i =0; i< 4; i++) {
        auto powni = NTL::power_long(2,3-i);
        auto powi = NTL::power_long(2,i);
        ecp Qtemp = powni*Ptemp;
        if (Qtemp == R2) {
            b1 += powi;
            Ptemp -= powi * R16;
        }
        else if (Qtemp == S2) {
            b2 += powi;
            Ptemp -= powi * S16;
        }
        else if (Qtemp == R2+S2) {
            b1 += powi;
            b2 += powi;
            Ptemp -= powi *(R16 + S16);
        }
        else {
            assert(Qtemp.is_identity());
        }
        // if ()
    }


    assert( P == b1* bas.first + b2* bas.second);
    return {b1,b2};
}

SmallMatFp change_of_basis3(const std::pair<ecp, ecp> &bas1, const std::pair<ecp, ecp> &bas2) {

    auto [a,b] = bi_dlp_3(bas2.first, bas1);
    auto [c,d] = bi_dlp_3(bas2.second, bas1);
    SmallMatFp M(3);
    M.mat[0][0] = a;
    M.mat[0][1] = b;
    M.mat[1][0] = c;
    M.mat[1][1] = d;

    return M;
}

SmallMatFp change_of_basis16(const std::pair<ecp, ecp> &bas1, const std::pair<ecp, ecp> &bas2) {

    auto [a,b] = bi_dlp_16(bas2.first, bas1);
    auto [c,d] = bi_dlp_16(bas2.second, bas1);
    SmallMatFp M(16);
    M.mat[0][0] = NTL::conv<SmallInteger>(a);
    M.mat[0][1] = NTL::conv<SmallInteger>(b);
    M.mat[1][0] = NTL::conv<SmallInteger>(c);
    M.mat[1][1] = NTL::conv<SmallInteger>(d);

    return M;
}
