#include "fast_ff.hpp"
#include <cassert>
#include <gperftools/profiler.h>

void to_Fp2 (Fp2 &t, const ffp2 &a) {
    if (is_zero(a)) {
        t = Fp2(0);
    }
    else if (NTL::IsZero(a.second)) {
        t._zz_pE__rep.SetLength(1);
        t._zz_pE__rep[0] = a.first;
    }
    else {
        t._zz_pE__rep.SetLength(2);
        t._zz_pE__rep[0] = a.first;
        t._zz_pE__rep[1] = a.second;
    }
}

Fp2 Fp2_cast(const ffp2 &a) {
    Fp2 t;
    to_Fp2(t, a);
    return t;
}

void negate(ffp2 &a, const ffp2 &b) {
    negate(a.first, b.first);
    negate(a.second, b.second);
}

ffp2 fast_Frob(const ffp2 &x) {
    Fp a;
    negate(a, x.second);
    return {x.first, a};
}

ffp2 from_Fp2(const Fp2 &a) {
    if (NTL::IsZero(a)) {
        return {Fp(0), Fp(0)};
    }
    else if (NTL::deg(a._zz_pE__rep) == 0) {
        return {a._zz_pE__rep[0], Fp(0)};
    }
    else {
        return {a._zz_pE__rep[0], a._zz_pE__rep[1]};
    }
}

bool is_zero(const ffp2 &a) {
    return NTL::IsZero(a.first) && NTL::IsZero(a.second); 
}

bool is_one(const ffp2 &a) {
    return NTL::IsOne(a.first) && NTL::IsZero(a.second); 
}


// ffp2
void fast_mul(ffp2 &c, const ffp2 &a, const ffp2 &b) {
    Fp temp1;
    add(temp1, a.first, a.second);
    mul(temp1, temp1, b.first + b.second);
    mul(c.first, a.first, b.first);
    mul(c.second, a.second, b.second);
    sub(temp1, temp1, c.first + c.second);
    sub(c.first, c.first, c.second);
    c.second._zz_p__rep = temp1._zz_p__rep;
} 



// void fast_mul_with_precomp(ffp2 &c, const ffp2 &a, const ffp2 &b, const std::vector<NTL::mulmod_precon_t> &prep) {
//     Fp temp1;
//     add(temp1, a.first, a.second);
//     mul(temp1, temp1, b.first + b.second);
//     mul(c.first, a.first, b.first);
//     mul(c.second, a.second, b.second);
//     sub(temp1, temp1, c.first + c.second);
//     sub(c.first, c.first, c.second);
//     c.second._zz_p__rep = temp1._zz_p__rep;
// }



void fast_sqr(ffp2 &s, const ffp2 &a) {
    Fp temp;
    mul(temp, (a.first + a.second), (a.first - a.second));
    mul(s.second, a.first, a.second);
    add(s.second, s.second, s.second);
    s.first._zz_p__rep = temp._zz_p__rep;
}

    
void fast_add(ffp2 &c, const ffp2 &a, const ffp2 &b) {
    add(c.first, a.first, b.first);
    add(c.second, a.second, b.second);
}

void fast_inv(ffp2 &y, const ffp2 &x) {
     
    // compute N = n(x) mod p
    Fp N = x.second * x.second + x.first * x.first;
    assert(N != 0);
    NTL::inv(N, N);
    // NTL::SetCoeff(y._zz_pE__rep, 0, N * NTL::coeff(x._zz_pE__rep, 0));
    // NTL::SetCoeff(y._zz_pE__rep, 1, -N * NTL::coeff(x._zz_pE__rep, 1));
    mul(y.first, N, x.first);
    mul(y.second, N, -x.second);
}

void get_powers(std::vector<ffp2> &powers, const ffp2 &a, int k) {
    // assumes the vector has been properly initialized
    assert(powers.size() >= (size_t) k + 1);
    powers[0] = {Fp(1), Fp(0)};
    for (int i = 1; i <=k; i++) {
        fast_mul(powers[i], powers[i - 1], a);
    }
}

void fast_pow(ffp2 &y, const ffp2 &x, long e) {

    ffp2 t = x;
    y = {Fp(1), Fp(0)};
    while (e > 0) {
        if (e & 1) {
            fast_mul(y, y, t);
        }
        fast_sqr(t, t);
        e >>= 1;
    }

}


// ffp2X
ffp2X::ffp2X(const std::vector<ffp2> &co, int degree) {
    this->deg = degree;
    assert(co.size() == (size_t) degree + 1);
    for (int i = 0; i <= degree; i++) {
        this->coeffs.push_back(co[i]);
    }
}

ffp2X::ffp2X(const Fp2X &poly) {
    this->deg = NTL::deg(poly);
    for (int i = 0; i <= this->deg; i++) {
        this->coeffs.push_back(from_Fp2(poly[i]));
    }
}

ffp2X::ffp2X(int degree, long coeff) {
    this->deg = degree;
    for (int i = 0; i < degree; i++) {
        this->coeffs.push_back({Fp(0), Fp(0)});
    }
    this->coeffs.push_back({Fp(coeff), Fp(0)});
}

bool ffp2X::is_zero_poly() const {

    // starts at the degree
    return this->deg < 0 || is_zero(this->coeffs[this->deg]);
}



void ffp2X::normalize() {
    while (is_zero(this->coeffs[this->deg])) {
        this->deg--;
    }
    if (this->deg >= 0) {
        this->coeffs.resize(this->deg + 1);
    }
    
}

void ffp2X::resize(int k) {
    this->coeffs.resize(k);
} 

ffp2 ffp2X::lead_coeff() {
    return this->coeffs[this->deg];
}

void swap(ffp2X &a, ffp2X &b) {
    std::swap(a.deg, b.deg);
    std::swap(a.coeffs, b.coeffs);
}

void ffp2X::mul(const ffp2 &x) {
    for (int i = 0; i <= this->deg; i++) {
        fast_mul(this->coeffs[i], this->coeffs[i], x);
    }
}

void sub(ffp2X &c, const ffp2X &a, const ffp2X &b) {

    
    if (a.deg >= b.deg ) {
        ffp2 t;
        c = ffp2X(a.coeffs, a.deg);
        for (int i = 0; i <= b.deg; i++) {
            negate(t, b.coeffs[i]);
            fast_add(c.coeffs[i], c.coeffs[i], t);
        }
    } 
    else {
        c = ffp2X( b.coeffs, b.deg);
        for (int i = 0; i <= a.deg; i++) {
            negate(c.coeffs[i], c.coeffs[i]);
            fast_add(c.coeffs[i], c.coeffs[i], a.coeffs[i]);
        }
        for (int i = a.deg + 1; i <= b.deg; i++) {
            negate(c.coeffs[i], c.coeffs[i]);
        }
    }

    if (a.deg == b.deg) {
        c.normalize();
    }

}


bool ffp2X::operator==(const ffp2X &other) {
    if (this->deg != other.deg) {
        return false;
    }
    else {
        if (this->deg < 0) {
            return true;
        }

        for (int i = 0; i <= this->deg; i++) {
            if (this->coeffs[i] != other.coeffs[i]) {
                return false;
            }
        }

    }
    return true;
}

void add(ffp2X &c, const ffp2X &a, const ffp2X &b) {

    
    if (a.deg >= b.deg ) {
        ffp2 t;
        c = ffp2X(a.coeffs, a.deg);
        for (int i = 0; i <= b.deg; i++) {
            fast_add(c.coeffs[i], a.coeffs[i], b.coeffs[i]);
        }
        // for (int i = b.deg + 1; i <= a.deg ) {
        //     c.coeffs[i] 
        // }
    } 
    else {
        c = ffp2X(b.coeffs, b.deg);
        for (int i = 0; i <= a.deg; i++) {
            fast_add(c.coeffs[i], b.coeffs[i], a.coeffs[i]);
        }
    }

    if (a.deg == b.deg) {
        c.normalize();
    }

}


// naive polynomial multiplication
void mul(ffp2X &c, const ffp2X &a, const ffp2X &b) {

    std::vector<ffp2> coefficients(a.deg + b.deg + 1, {Fp(0), Fp(0)});
    // 

    ffp2 temp;

    for (int i = 0; i <= a.deg; i++) {
        for (int j = 0; j <= b.deg; j++) {
            fast_mul(temp, a.coeffs[j], b.coeffs[j] );
            fast_add(coefficients[i + j], coefficients[i + j], temp);
        }
        
    } 
    c.deg = a.deg + b.deg;
    c.resize(c.deg + 1);
    for (int i = 0; i <= c.deg; i++) {
        c.coeffs[i] = coefficients[i];
    }
    // c = ffp2X(coefficients, a.deg + b.deg);
    
}


void ffp2X::Evaluate(ffp2 &x, const ffp2 &a) const {
    std::vector<ffp2> powers(this->deg + 1);
    get_powers(powers, a, this->deg);
    x = {Fp(0), Fp(0)};
    ffp2 t;
    for (int i = 0; i <= this->deg; i++) {
        fast_mul(t, this->coeffs[i], powers[i]);
        fast_add(x, x, t);
    }
}

void ffp2X::fast_PlainRem(const ffp2X& b) {
    long db, dq, da;
    long i, j, LCIsOne;

    ffp2 LCInv, t;
    ffp2 s;

    da = this->deg;
    db = b.deg;
    
    assert(db >= 0);

    // we do nothing
    if (da < db) {
        return;
    }

    {
        if (is_one(b.coeffs[db]))
            LCIsOne = 1;
        else {
            LCIsOne = 0;
            fast_inv(LCInv, b.coeffs[db]);
        }

    }

    
    dq = da - db;

    for (i = dq; i >= 0; i--) {
        if (!LCIsOne) {
            fast_mul(t, this->coeffs[i + db], LCInv);
            negate(t, t);
        }
        else {
            negate(t, this->coeffs[i + db]);
        }
            
        
        for (j = db-1; j >= 0; j--) {

            fast_mul(s, t, b.coeffs[j]);
            fast_add(this->coeffs[i + j], this->coeffs[i + j], s);

        }
    }
    this->deg = db - 1;
    this->normalize();

}

// exploits a known special structure of b to be much faster (b = X^k + aX + b)
void special_fast_PlainRem(ffp2X &a, const ffp2X& b) {
    long db, dq, da;
    long i, LCIsOne;

    ffp2 LCInv, t;
    ffp2 s;

    da = a.deg;
    db = b.deg;
    
    assert(db >= 0);

    // we do nothing
    if (da < db) {
        return;
    }

    {
        if (is_one(b.coeffs[db]))
            LCIsOne = 1;
        else {
            LCIsOne = 0;
            fast_inv(LCInv, b.coeffs[db]);
        }

    }

    
    dq = da - db;

    for (i = dq; i >= 0; i--) {
        if (!LCIsOne) {
            fast_mul(t, a.coeffs[i + db], LCInv);
            negate(t, t);
        }
        else {
            negate(t, a.coeffs[i + db]);
        }
            
        fast_mul(s, t, b.coeffs[1]);
        fast_add(a.coeffs[i + 1], a.coeffs[i + 1], s);
        fast_mul(s, t, b.coeffs[0]);
        fast_add(a.coeffs[i], a.coeffs[i], s);
    }
    a.deg = db - 1;
    a.normalize();

}



void fast_PlainGCD_2(ffp2X& x, const ffp2X& a, const ffp2X& b)
{
   ffp2 t;

   if (b.is_zero_poly())
      x = ffp2X(a.coeffs, a.deg);
   else if (a.is_zero_poly())
      x = ffp2X(b.coeffs, b.deg);
   else {
      ffp2X u, v;

      u = ffp2X(a.coeffs, a.deg);
      v = ffp2X(b.coeffs, b.deg);;
      do {
        u.fast_PlainRem(v);
        swap(u, v);
      } while (!v.is_zero_poly());

      x = ffp2X(u.coeffs, u.deg);
   }

   if (x.is_zero_poly()) return;
   if (is_one(x.lead_coeff())) return;

   /* make gcd monic */
   fast_inv(t, x.lead_coeff());
   x.mul(t);
}



// ffp2XY
ffp2XY::ffp2XY(const std::vector<ffp2X> &cs) {

    int l = cs.size();
    this->degY = l-1;

    // std::cout << "size = " << l << "\n";

    int DX = 0;
    ffp2X cos;
    for (int i = 0; i < l; i++){

        // std::cout << " i = " << i << "deg = " << cs[i].deg << cs[i] << "\n";

        auto D = cs[i].deg;
        cos = ffp2X(cs[i].coeffs, D)    ;
        this->coeffs.push_back(cs[i]);
        if(D > DX){
            DX = D;
        }
    }
    this->degX = DX;
}

ffp2XY::ffp2XY(const ZZ_pEXY &F) {
    this->degX = F.dX;
    this->degY = F.dY;

    ffp2X cos;
    
    for (int i = 0; i <= degY; i++) {
        cos = ffp2X(F.coeffs[i]);
        this->coeffs.push_back(cos);
    }
}

void ffp2XY::Evaluate(ffp2X &x, const ffp2 &a) const {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Evaluate bivariate in variables X,Y at Y = a to get a polynomial in X.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<ffp2> powers(this->degY + 1);
    get_powers(powers, a, this->degY);

    std::vector<ffp2> new_cs(this->degX + 1, {Fp(0), Fp(0)});
    
    ffp2 t;

    for(int i = 0; i <= this->degY; i++) {
        for (int j = 0; j <= this->coeffs[i].deg; j++) {
            fast_mul(t, powers[i], this->coeffs[i].coeffs[j]);
            fast_add(new_cs[j], new_cs[j], t); 
        }
    }

    x = ffp2X(new_cs, this->degX);
}

void ffp2XY::Evaluate(ffp2X &x, const std::vector<ffp2> &powers) const {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Evaluate bivariate in variables X,Y at Y = a to get a polynomial in X.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // std::cout << "in evaluate \n";
    std::vector<ffp2> new_cs(this->degX + 1, {Fp(0), Fp(0)});
    
    ffp2 t;

    for(int i = 0; i <= this->degY; i++) {
        for (int j = 0; j <= this->coeffs[i].deg; j++) {
            fast_mul(t, powers[i], this->coeffs[i].coeffs[j]);
            fast_add(new_cs[j], new_cs[j], t); 
        }
    }
    // std::cout << "before create \n";
    x = ffp2X(new_cs, this->degX);
    // std::cout << "after create" << "\n";
}

void ffp2XY::Evaluate(ffp2 &t, const ffp2 &x, const ffp2 &y) const {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Evaluate bivariate in variables X,Y at Y = a to get a polynomial in X.
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<ffp2> powersy(this->degY);
    get_powers(powersy, y, this->degY);

    ffp2X evaly;
    this->Evaluate(evaly, powersy);
    evaly.Evaluate(t, x);
}



// ffp2XYZ 
ffp2XYZ::ffp2XYZ(const ZZ_pEXYZ &F) {
    this->degZ = F.dZ;

    for (int i = 0; i <= this->degZ; i++) {
        this->coeffs.push_back(ffp2XY(F.coeffs[i]));
    }
}


/// Karatsuba multiplication (recursive)
void karatsuba_mul(ffp2X &C, const ffp2X &A, const ffp2X &B) {

    // std::cout << "Karatsuba mul with " << A.deg << " " << B.deg << "\n";
    
    if (A.deg < 0 || B.deg < 0) {
        ffp2X R;
        C = R;
    }

    assert(A.deg >=0 && B.deg >= 0);
    int n = std::max(A.deg, B.deg) + 1;

    // Threshold: below this, use naive multiplication
    if (A.deg < 2 || B.deg < 2 || n <= 4) { // tweak depending on your platform
        return mul(C, A, B);
    }

    int m = n / 2;

    // Split A into A0 + A1*x^m
    ffp2X A0(std::min((int) A.deg, m - 1));
    ffp2X A1(std::max(0, (int) A.deg - m));
    A0.coeffs.assign(A.coeffs.begin(), A.coeffs.begin() + std::min<int>(A.deg + 1, m));
    if (A.deg >= m)
        A1.coeffs.assign(A.coeffs.begin() + m, A.coeffs.end());
    // normalize(A0);
    // normalize(A1);
    A0.normalize();
    A1.normalize();

    // Split B into B0 + B1*x^m
    ffp2X B0(std::min((int) B.deg, m - 1));
    ffp2X B1(std::max(0, (int) B.deg - m));
    B0.coeffs.assign(B.coeffs.begin(), B.coeffs.begin() + std::min<int>(B.deg + 1, m));
    if (B.deg >= m)
        B1.coeffs.assign(B.coeffs.begin() + m, B.coeffs.end());
    B0.normalize();
    B1.normalize();

    // Recursive calls:
    ffp2X Z0, Z2, Z1, A0A1, B0B1;
    // std::cout << "1st sub call to " << A0.deg << " " << B0.deg << "\n";
    karatsuba_mul(Z0, A0, B0);
    // std::cout << "1st call to " << A0.deg << " " << B0.deg << " done \n";             // A0 * B0
    // std::cout << "2nd sub call to " << A1.deg << " " << B1.deg << "\n";
    karatsuba_mul(Z2, A1, B1);
    // std::cout << "2nd call to " << A1.deg << " " << B1.deg << " done \n";              // A1 * B1
    add(A0A1, A0, A1);                // A0 + A1
    add(B0B1, B0, B1);                // B0 + B1
    // std::cout << "3rd sub call to " << A0A1.deg << " " << B0B1.deg << "\n";
    karatsuba_mul(Z1, A0A1, B0B1);         // (A0 + A1)(B0 + B1)
    // std::cout << "3rd call to " << A0A1.deg << " " << B0B1.deg << " done \n";     
    


    // Z1 = Z1 - Z0 - Z2
    sub(Z1, Z1, Z0);
    sub(Z1, Z1, Z2);

    // std::cout << Z0 << "\n" << Z1 << "\n" <<  Z2 << "\n";
    // // Combine result:
    // // R = Z0 + Z1*x^m + Z2*x^(2m)
    ffp2X R(Z2.deg + 2*m);
    // copy Z0
    for (int i = 0; i <= Z0.deg; i++) {
        R.coeffs[i] = Z0.coeffs[i];
    }
    // add Z1 shifted by m
    for (int i = 0; i <= Z1.deg; i++) {
        fast_add(R.coeffs[i + m], R.coeffs[i + m], Z1.coeffs[i]);
    }
    // add Z2 shifted by 2m
    for (int i = 0; i <= Z2.deg; i++) {
        fast_add(R.coeffs[i + 2*m], R.coeffs[i + 2*m], Z2.coeffs[i]);
    }

    

    R.normalize();
    // std::cout << "R= " << R << "\n";
    C = R;
}

/// Public wrapper
void poly_mul_karatsuba(ffp2X &C, const ffp2X &A, const ffp2X &B) {
    karatsuba_mul(C, A, B);
}


// ffp2k 
void ffp2k::normalize() {
    special_fast_PlainRem(this->def_poly, *this->mod_poly);
    this->def_poly.normalize();
}

void mul(ffp2k &c, const ffp2k &b, ffp2k &a) {
    // std::cout << "mul with " << a.def_poly.deg << " " << b.def_poly.deg << "\n";

    // ffp2X test;
    // poly_mul_karatsuba(test, a.def_poly, b.def_poly);
    // test.normalize();

    // assumes that a,b,c have the same mod_poly
    mul(c.def_poly, a.def_poly, b.def_poly);
    // std::cout << test << "\n" <<c.def_poly << "\n"; 

     c.normalize();
}

ffp2k random_ffp2k(const ffp2X *mod_poly) 
{
    std::vector<ffp2> coefficients(mod_poly->deg);
    for (int i = 0; i < mod_poly->deg; i++) {
        coefficients[i] = {NTL::random_zz_p(), NTL::random_zz_p()};
    }

    ffp2X def_pol; 
    def_pol = ffp2X(coefficients, mod_poly->deg - 1);
    // special_fast_PlainRem(def_pol, *mod_poly);

    return ffp2k(def_pol, mod_poly);

}




// Fast Fp2 handling 
// removing all NTL annoying stuff

void fast_set(Fp2 &c, const FpX &t) {
    c._zz_pE__rep.SetLength(NTL::deg(t) + 1);
    for (int i=0; i <= NTL::deg(t); i++) {
        c._zz_pE__rep[i]._zz_p__rep = t[i]._zz_p__rep;
    }
}

void fast_mul(Fp2 &c, const Fp2 &a, const Fp2 &b) {
    Fp temp;
    mul(temp, NTL::coeff(a._zz_pE__rep, 0) + NTL::coeff(a._zz_pE__rep, 1), NTL::coeff(b._zz_pE__rep, 0) + NTL::coeff(b._zz_pE__rep, 1));
    // mul(c._zz_pE__rep[0], a._zz_pE__rep[0], b._zz_pE__rep[0]);
    NTL::SetCoeff(c._zz_pE__rep, 0, NTL::coeff(a._zz_pE__rep, 0) * NTL::coeff(b._zz_pE__rep, 0));
    // mul(c._zz_pE__rep[1], a._zz_pE__rep[1], b._zz_pE__rep[1]);
    NTL::SetCoeff(c._zz_pE__rep, 1, NTL::coeff(a._zz_pE__rep, 1) * NTL::coeff(b._zz_pE__rep, 1));
    sub(temp, temp, NTL::coeff(c._zz_pE__rep, 0) + NTL::coeff(c._zz_pE__rep, 1));
    // sub(c._zz_pE__rep[0], c._zz_pE__rep[0], c._zz_pE__rep[1]);
    NTL::SetCoeff(c._zz_pE__rep, 0, NTL::coeff(c._zz_pE__rep, 0) - NTL::coeff(c._zz_pE__rep, 1));
    // c._zz_pE__rep[1] = temp;
    NTL::SetCoeff(c._zz_pE__rep, 1, temp);

}


void fast_mul_rep(FpX &c, const FpX &a, const FpX &b) {
    
    if (deg(a) == 0 && deg(b) == 0) {
        // c[0] = a[0] * b[0];
        // c[1] = Fp(0);
        mul(c[0], a[0], b[0]);
        c[1]._zz_p__rep = 0;
    }
    else if (deg(a) == 0) {
        // c[0] = a[0] * b[0];
        // c[1] = a[0] * b[1];
        mul(c[0], a[0], b[0]);
        mul(c[1], a[0], b[1]);
    }
    else if (deg(b) == 0) {
        // c[0] = b[0] * a[0];
        // c[1] = b[0] * a[1];
        mul(c[0], b[0], a[0]);
        mul(c[1], b[0], a[1]);
    }
    else {
        Fp temp;
        mul(temp, a[0] + a[1], b[0] + b[1]);
        mul(c[0], a[0], b[0]);
        mul(c[1], a[1], b[1]);
        sub(temp, temp, c[0] + c[1]);
        sub(c[0], c[0], c[1]);
        // c[1] = temp;
        c[1]._zz_p__rep = temp._zz_p__rep;
    }

    
}



void faster_mul(Fp2 &c, const Fp2 &a, const Fp2 &b, FpX &t) {
#ifndef NDEBUG
    Fp2 reca = a;
    Fp2 recb = b;
#endif
    if (deg(a._zz_pE__rep) < 0 || deg(b._zz_pE__rep) < 0) {
        c._zz_pE__rep = FpX(0);
    }
    else {
        // c._zz_pE_rep.SetLength(2);
        fast_mul_rep(t, a._zz_pE__rep, b._zz_pE__rep);
        // conv(c._zz_pE__rep, t);
        fast_set(c, t);
        
    }
#ifndef NDEBUG 
    Fp2 test;
    fast_mul(test, reca, recb);
    if (deg(rep(test)) == 0) {
        assert( test == c || (deg(rep(c)) == 1 && rep(test)[0] == rep(c)[0] && rep(c)[1] == 0));
    }
    else if (NTL::IsZero(test)) {
        assert(NTL::IsZero(c) || (rep(c)[0] == 0 && rep(c)[1] == 0));
    }
    else {
        assert(test == c);
    }
    // std::cout << rep(test) << " " << rep(c) << "\n";
    // assert(test == c);
#endif

}




void fast_sqr(Fp2 &s, const Fp2 &a) {
#ifndef NDEBUG
    Fp2 rec = a;
#endif
    
    Fp temp;
    temp = (NTL::coeff(a._zz_pE__rep, 0) + NTL::coeff(a._zz_pE__rep, 1)) * (NTL::coeff(a._zz_pE__rep, 0) - NTL::coeff(a._zz_pE__rep, 1));
    NTL::SetCoeff(s._zz_pE__rep, 1, NTL::coeff(a._zz_pE__rep, 0) * NTL::coeff(a._zz_pE__rep, 1));
    NTL::SetCoeff(s._zz_pE__rep, 1, NTL::coeff(s._zz_pE__rep, 1) + NTL::coeff(s._zz_pE__rep, 1));
    NTL::SetCoeff(s._zz_pE__rep, 0, temp);
#ifndef NDEBUG 
    assert(Fp2(rec * rec) == s);
#endif
}


void fast_pow(Fp2 &y, const Fp2 &x, long e) {
    FpX temp;
    temp.SetLength(2);
#ifndef NDEBUG 
    Fp2 rec = x;
    long rec_e = e;
#endif
    Fp2 t = x;
    y = Fp2(1);
    while (e > 0) {
        if (e & 1) {
            faster_mul(y, y, t, temp);
        }
        fast_sqr(t, t);
        e >>= 1;
    }
#ifndef NDEBUG
    assert(y == NTL::power(rec, rec_e));
#endif 

}


bool fast_sqrt(Fp2 &s, const Fp2 &a) {
    
#ifndef NDEBUG 
    Fp2 rec = a;
#endif
    if (a == 0) {
        s = 0;
        return true;
    }

    // Fp2 test;
    // fast_pow(test, a, (a.modulus().val()[0].modulus() * a.modulus().val()[0].modulus()  -1) >> 1);
    if (1) {
        
        Fp x0 = NTL::power( NTL::power(NTL::coeff(a._zz_pE__rep, 0), 2) + NTL::power(NTL::coeff(a._zz_pE__rep, 1), 2), (a.modulus().val()[0].modulus() + 1) >> 2 );
        x0 = NTL::coeff(a._zz_pE__rep, 0) + x0;
        Fp t0 = x0 + x0;
        Fp x1 = NTL::power(t0, (a.modulus().val()[0].modulus() - 3) >> 2);
        x0 = x1 * x0;
        x1 = x1 * NTL::coeff(a._zz_pE__rep, 1);
        Fp t1 = NTL::power(2 * x0, 2);
        if (t0 == t1) {
            NTL::SetCoeff(s._zz_pE__rep, 0, x0);
            NTL::SetCoeff(s._zz_pE__rep, 1, x1);
        }
        else {
            NTL::SetCoeff(s._zz_pE__rep, 0, x1);
            NTL::SetCoeff(s._zz_pE__rep, 1, -x0);
        }
        if ( s==0 ) {
            return false;
        }
        #ifndef NDEBUG
            Fp2 new_test = s;
            fast_sqr(new_test, new_test);
            assert(new_test == rec);
        #endif
        return true;
    }
    else {
        assert(0);
        return false;
    }


}


void fast_inv(Fp2 &y, const Fp2 &x) {
     
    // compute N = n(x) mod p
    Fp N = NTL::coeff(x._zz_pE__rep, 1) * NTL::coeff(x._zz_pE__rep, 1) + NTL::coeff(x._zz_pE__rep, 0) * NTL::coeff(x._zz_pE__rep, 0);
    assert(N != 0);
    NTL::inv(N, N);
    NTL::SetCoeff(y._zz_pE__rep, 0, N * NTL::coeff(x._zz_pE__rep, 0));
    NTL::SetCoeff(y._zz_pE__rep, 1, -N * NTL::coeff(x._zz_pE__rep, 1));

}

// 
bool correc_is_zero(FpEX_elem& x) {
    for (int i = deg(x); i >= 0; i--) {
        // std::cout << (int) NTL::IsZero(x[i]) << "\n";
        // std::cout << (int) (deg(x[i]._zz_pE__rep) == 1) << (int) (x[i]._zz_pE__rep[0] == 0) << (int) (x[i]._zz_pE__rep[1] == 0) << "\n"; 
        // std::cout << (int) (deg(x[i]._zz_pE__rep) == 1 && x[i]._zz_pE__rep[0] == 0 && x[i]._zz_pE__rep[1] == 0) << "\n";
        if (!NTL::IsZero(x[i]) && !(deg(x[i]._zz_pE__rep) == 1 && x[i]._zz_pE__rep[0] == 0 && x[i]._zz_pE__rep[1] == 0) ) {
            x.normalize();
            return false;
        }
        x[i]._zz_pE__rep.normalize();
    }
    x.normalize();
    return true;
}

void fast_add(Fp2 &c, const Fp2 &a, const Fp2 &b, FpX &temp) {
    if (NTL::deg(a._zz_pE__rep) < 0) {
        c = b;
    }
    else if (NTL::deg(b._zz_pE__rep) < 0) {
        c = a;
    }
    else if (NTL::deg(b._zz_pE__rep) == 0) {
        temp = a._zz_pE__rep;
        add(temp[0], temp[0], b._zz_pE__rep[0]);
        fast_set(c, temp);
    }
    else if (NTL::deg(a._zz_pE__rep) == 0) {
        temp = b._zz_pE__rep;
        add(temp[0], temp[0], a._zz_pE__rep[0]);
        fast_set(c, temp);
        // conv(c._zz_pE__rep, temp);
    }
    else {
        add(temp[0], a._zz_pE__rep[0], b._zz_pE__rep[0]);
        add(temp[1], a._zz_pE__rep[1], b._zz_pE__rep[1]);
        // conv(c._zz_pE__rep, temp);
        fast_set(c, temp);
    }   
    
}

void fast_PlainRem(FpEX_elem& r, const FpEX_elem& a, const FpEX_elem& b)
{
    long db, dq;
    long i, j, LCIsOne;
    const Fp2 *bp;
    
    FpX temp;
    temp.SetLength(2);

    bp = b.rep.elts();

    Fp2 LCInv, t;
    Fp2 s;

    long da = deg(a);
    
    
    // Fp2 *xp = new Fp2(1);
    std::vector<Fp2> xp(da + 1);
    db = deg(b);

    assert(db >= 0);

    if (da < db) {
        r = a;
        return;
    }

    if (NTL::IsOne(bp[db]))
        LCIsOne = 1;
    else {
        LCIsOne = 0;
        fast_inv(LCInv, bp[db]);
    }
    dq = da - db;


    for (i = 0; i <= da; i++) {
        xp[i]._zz_pE__rep = a[i]._zz_pE__rep;
    }

    for (i = dq; i >= 0; i--) {
        t =  xp[i + db];
        if (!LCIsOne) {
            faster_mul(t, t, LCInv, temp);
        }
            
        negate(t, t);

        for (j = db-1; j >= 0; j--) {

            faster_mul(s, t, bp[j], temp);
            fast_add(xp[i + j], xp[i + j], s, temp);
            
            // if (i==0 && j==0 ){
                // std::cout << " 0 0 " << bp[j] << t << s << xp[i + j] << "\n";
            // }
        }
    }

    // std::cout << xp[0] << "\n";
    r.SetLength(db);
    for (i = 0; i < db; i++) {
        r[i]._zz_pE__rep = xp[i]._zz_pE__rep;
    }

    r.normalize();

}



void fast_PlainGCD(FpEX_elem& x, const FpEX_elem& a, const FpEX_elem& b)
{
   Fp2 t;

//    std::cout << "starting GCD with A = " << a <<"\nB = " << B << "\n";

   if (NTL::IsZero(b))
      x = a;
   else if (NTL::IsZero(a))
      x = b;
   else {
      long n = std::max(deg(a),deg(b)) + 1;
      FpEX_elem u, v;
      u.SetLength(n);
      v.SetLength(n);

      u = a;
      v = b;
      do {
        fast_PlainRem(u, u, v);
        swap(u, v);
      } while (!correc_is_zero(v));

      x = u;
   }

   if (correc_is_zero(x)) return;
   if (NTL::IsOne(NTL::LeadCoeff(x))) return;

   /* make gcd monic */
   fast_inv(t, NTL::LeadCoeff(x));
//    NTL::inv(t, NTL::LeadCoeff(x)); 
   NTL::mul(x, x, t); 
}


void SpecialPlainReminder(FpX_elem& r, const FpX_elem& a, const FpX_elem& b)
{
   long da, db, dq, i, LCIsOne;
   const Fp *bp;
   Fp *xp;


   Fp LCInv, t;
   Fp s;

   da = deg(a);
   db = deg(b);

   if (da < db) {
      r = a;
      return;
   }

   bp = b.rep.elts();

   if (IsOne(bp[db]))
      LCIsOne = 1;
   else {
      LCIsOne = 0;
      inv(LCInv, bp[db]);
   }

   NTL::vec_zz_p x;

   if (&r == &a)
      xp = r.rep.elts();
   else {
      x = a.rep;
      xp = x.elts();
   }

   dq = da - db;

   long p = Fp::modulus();
   NTL::mulmod_t pinv = Fp::ModulusInverse();
   std::vector<long> indices = {0, 1, 2, db /2 ,(db/2) +1};

   for (i = dq; i >= 0; i--) {
      t = xp[i+db];
      if (!LCIsOne)
         mul(t, t, LCInv);
      negate(t, t);

      long T = rep(t);
      NTL::mulmod_precon_t Tpinv = NTL::PrepMulModPrecon(T, p, pinv); 

    
      for (long j : indices) {
            long S = NTL::MulModPrecon(rep(bp[j]), T, p, Tpinv);
            S = NTL::AddMod(S, rep(xp[i+j]), p);
            xp[i+j].LoopHole() = S;
         
      }
   }

   r.rep.SetLength(db);
   if (&r != &a) {
      for (i = 0; i < db; i++)
         r.rep[i] = xp[i];
   }
   r.normalize();
}

// 
void SpecialMul(FpE_elem& x, const FpE_elem& a, const FpE_elem& b)
{

#ifndef NDEBUG 
    // this uses a special function to reduce to the modulus assuming the modulus polynomial has a specific shape
    auto PP = FpE_elem::modulus();
    // std::cout << PP << "\n";
    for (int i = 0; i < deg(PP); i++) {
        if (i != 0 && i != 1 && i != 2 && i != deg(PP)/2 && i != (deg(PP)/2 + 1)) {
            // std::cout << i << "/" << deg(PP) << " " << NTL::coeff(PP,i) << "\n";
            assert(NTL::coeff(PP,i) == 0);
        }
        
    }
    auto testa = a;
    auto testb = b;
#endif 

   FpX_elem t;

   mul(t, a._zz_pE__rep, b._zz_pE__rep);
   if (deg(b._zz_pE__rep) <= 3) {
    rem(x._zz_pE__rep, t, FpE_elem::modulus());
   }
   else {
    SpecialPlainReminder(x._zz_pE__rep, t, FpE_elem::modulus());
   }
   
   
#ifndef NDEBUG 
   assert(x == testa * testb);
#endif
}
