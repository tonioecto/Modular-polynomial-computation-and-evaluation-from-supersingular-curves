#include "fast_quaternions.hpp"
#include "utils.hpp"
#include <numeric>

// --------------------------  Integer Functions ---------------------------- //

FastInteger convert(const Integer &x) {
    Integer pow2 = Integer( (long)2<<29);
    Integer y = x % pow2;
    long low_bits = NTL::conv<long>( y );
    y = (x - y)/pow2;
    long high_bits = NTL::conv<long>( y );

    FastInteger result = (FastInteger) high_bits;
    result = result << 30;
    result = result + (FastInteger) low_bits;
    return result; 

}

// not used anymore
FastInteger Euclid(const FastInteger &a, const FastInteger &b) {
    if (b == 0) {
        return a;
    }
    FastInteger gcd = Euclid(b, a % b);

    return gcd;
} 


// Function to implement
// Stein's Algorithm
FastInteger bin_gcd(FastInteger a, FastInteger b)
{
    /* GCD(0, b) == b; GCD(a, 0) == a,
       GCD(0, 0) == 0 */
    if (a == 0)
        return b;
    if (b == 0)
        return a;

    /*Finding K, where K is the
      greatest power of 2
      that divides both a and b. */
    unsigned char k;
    for (k = 0; ((a | b) & 1) == 0; ++k) 
    {
        a >>= 1;
        b >>= 1;
    }

    /* Dividing a by 2 until a becomes odd */
    while ((a & 1) == 0)
        a >>= 1;

    /* From here on, 'a' is always odd. */
    do
    {
        /* If b is even, remove all factor of 2 in b */
        while ((b & 1) == 0)
            b >>= 1;

        /* Now a and b are both odd.
           Swap if necessary so a <= b,
           then set b = b - a (which is even).*/
        if (a > b)
            std::swap(a, b); // Swap u and v.

        b = (b - a);
    }while (b != 0);

    /* restore common factors of 2 */
    return a << k;
}


// TODO there some Thomas pornin optimizations that we could potentially use?
FastInteger binary_extended_gcd(FastInteger a, FastInteger b, FastInteger &x, FastInteger &y) {
    // base case: if a or b is 0, return the other number as gcd, with appropriate coefficients
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    if (b == 0){
        x = 1;
        y = 0;
        return a;
    } 

    /*Finding K, where K is the
      greatest power of 2
      that divides both a and b. */
    unsigned char k;
    for (k = 0; ((a | b) & 1) == 0; ++k) 
    {
        a >>= 1;
        b >>= 1;
    }

    // Step 2: Initialize variables for the extended GCD
    FastInteger u = a;
    FastInteger v = b;             
    x = 1; y = 0;           
    FastInteger x2 = 0; FastInteger y2 = 1;            

    //  Step 3: Main loop to compute GCD while updating Bézout coefficients
    while (u != v) {
        if (!IsOdd(u)) {
            // If u is even, divide by 2
            u >>= 1;
            // Adjust (x1, y1) to keep Bézout identity valid
            if (((x | y) & 1) == 0) {
                x >>= 1;
                y >>= 1;
            } 
            else {
                // If odd, we adjust using the original a, b to keep x1, y1 as integers
                x = (x + b) >> 1;
                y = (y - a) >> 1;
            }
        }
        else if (!IsOdd(v)) {
            //  Same logic for v if it is even
            v >>= 1;
            if (((x2 | y2) & 1) == 0) {
                x2 >>= 1;
                y2 >>= 1;
            } 
            else {
                // If odd, we adjust using the original a, b to keep x1, y1 as integers
                x2 = (x2 + b) >> 1;
                y2 = (y2 - a) >> 1;
            }
        }
        else if (u > v) {
            //  Reduce u by subtracting v
            u -= v;
            x -= x2;
            y -= y2; 
        }
        else {
            // Reduce v by subtracting u
            v -= u;
            x2 -= x;
            y2 -= y;
        }
    }

    // Step 4: Multiply back the common factor of 2 (if any)
    // At this point, u == v == gcd(a, b)
    return u << k;  
}


FastInteger extendedEuclid(const FastInteger &a, const FastInteger &b, FastInteger &x, FastInteger &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }

    FastInteger x1, y1;
    FastInteger gcd = extendedEuclid(b, a % b, x1, y1);

    x = y1;
    y = x1 - (a / b) * y1;

    return gcd;
}


FastInteger GCD(const FastInteger &a1, const FastInteger &a2) {

    return bin_gcd(abs(a1), abs(a2));

}

// Function to compute the inverse of a modulo m
// Result is -1 if the result does not exist
// red was init with m so m = red.b
FastInteger InvMod(const FastInteger &a, const SignedBarrettReducer &red) {
    FastInteger x, y;
    FastInteger g = extendedEuclid(a, red.b, x, y);

    if (g != 1) {
        // L'inverse n'existe pas si a et m ne sont pas premiers entre eux
        return -1;
    } else {
        // On s'assure que l'inverse est positif
        return red.mod(x);
    }
}

// same as above but now we want the inverse of a mod m^2 
FastInteger InvModSqr(const FastInteger &a, const SignedBarrettReducer &red) {
    FastInteger x, y;
    FastInteger g = extendedEuclid(a, red.bsqr, x, y);

    if (g != 1) {
        // L'inverse n'existe pas si a et m ne sont pas premiers entre eux
        return -1;
    } else {
        // On s'assure que l'inverse est positif
        return red.modsqr(x);
    }
}

// same as above but now we want the inverse of a mod 2 m^2 
FastInteger InvMod2Sqr(const FastInteger &a, const SignedBarrettReducer &red) {
    FastInteger x, y;
    FastInteger g = extendedEuclid(a, red.bsqr << 1, x, y);

    if (g != 1) {
        // L'inverse n'existe pas si a et m ne sont pas premiers entre eux
        return -1;
    } else {
        // On s'assure que l'inverse est positif
        return red.mod2sqr(x);
    }
}


// same as above but now we want the inverse of a mod 2m 
FastInteger InvMod2(const FastInteger &a, const SignedBarrettReducer &red) {
    FastInteger x, y;
    FastInteger g = extendedEuclid(a, red.b << 1, x, y);

    if (g != 1) {
        // L'inverse n'existe pas si a et m ne sont pas premiers entre eux
        return -1;
    } else {
        // On s'assure que l'inverse est positif
        return red.mod2(x);
    }
}

// tried to improve InvModSqr with this but for now it does not work (due to sizes constraints), when using it on a case where normsqr is too big it overflows and it is a mess
// there's probably a good way of doing it but it is really not clear if it will be that much faster than the old method
FastInteger NewInvMod(const FastInteger &y, const FastInteger &m, const SignedBarrettReducer &red, FastInteger precompInv, unsigned char iter_bound) {
    FastInteger a = y;
    FastInteger b = m;
    unsigned char r = 0;
    FastInteger f0 = 1;
    FastInteger g0 = 0;
    FastInteger f1 = 0;
    FastInteger g1 = 1;

    while (a != 0) {
        r++;
        if (!IsOdd(a)) {
            a >>= 1;
        }
        else {
            if (a < b) {
                std::swap(a,b);
                // std::swap(u,v);
                std::swap(f0,f1);
                std::swap(g0,g1);
            }
            a = (a - b) >> 1;
            f0 = (f0 - f1);
            g0 = (g0 - g1);
            // if (IsOdd(u) == IsOdd(v)) {
                // u = (u -v) >> 1;
            // }
            // else {
                // u = red.modsqr((u - v) * k);
            // }
            
        }
        f1 <<= 1; 
        g1 <<= 1;
    }
    f1 = red.modsqr(f1);
    f1 <<= (iter_bound - r);
    f1 = red.modsqr(f1);

    return red.mulmodsqr(f1, precompInv);
}


// CRT mod (p1 x 2(p2)) where p1, p2 are two odd numbers
// redM is init with p1 p2 
// red2 is init with p2 
// inv is p1 mod (2 p2)
FastInteger CRT(const FastInteger &a1, const FastInteger &p1, const FastInteger &a2, const SignedBarrettReducer &redM, const SignedBarrettReducer &red2, const FastInteger &inv) {

    // Formule : x ≡ a1 + (a2 - a1) * (p1^{-1} mod (2 p2)) × p1  (mod M)
    FastInteger result = redM.mod2(a1 + red2.mod2((a2 - a1) * inv) * p1);
    return result;
    // return result;
}

// same as above but we want the result mod (p1')^2 p2^2 
// p2 is red2.b 
// p1 = (p1')^2 
// inv is p1 mod p2^2
FastInteger CRTsqr(const FastInteger &a1, const FastInteger &p1, const FastInteger &a2, const SignedBarrettReducer &redM, const SignedBarrettReducer &red2, const FastInteger &inv) {

    // Formule : x ≡ a1 + (a2 - a1) * p1^{-1} mod p2 × p1  (mod M)
    FastInteger result = redM.modsqr(a1 + red2.modsqr((a2 - a1) * inv) * p1);
    return result;
    // return result;
}

// same as above but we want the result mod 2 (p1')^2 p2^2 
// p2 is red2.b 
// p1 = (p1')^2 
// inv is p1 mod p2^2
FastInteger CRT2sqr(const FastInteger &a1, const FastInteger &p1, const FastInteger &a2, const SignedBarrettReducer &redM, const SignedBarrettReducer &red2, const FastInteger &inv) {

    // Formule : x ≡ a1 + (a2 - a1) * p1^{-1} mod p2 × p1  (mod M)
    FastInteger result = redM.mod2sqr(a1 + red2.mod2sqr((a2 - a1) * inv) * p1);
    return result;
    // return result;
}

FastFloat IntToFloat(FastInteger a1, FastInteger a2) {
    return ((FastFloat) a1)/((FastFloat) a2);
}

// (a1/a2)^2
FastFloat IntToFloatSqr(FastInteger a1, FastInteger a2) {
    auto a = ((FastFloat) a1)/((FastFloat) a2);
    return a * a;
}

FastInteger Rounding(FastFloat f) {
    return llround(f);  
}

bool IsOdd(const FastInteger a) {
    return (a & (FastInteger) 1);
}

FastInteger SqrRoot(const FastInteger a) {
    return (FastInteger) sqrt((FastFloat) a);
}

bool IsIntegralSquare(const FastInteger a) {
    auto s = SqrRoot(a);
    return (s*s == a);
}

// Convert long long to byte array
void FastIntToBytes(FastInteger value, unsigned char* buffer, int BUFFER_SIZE)
{
    for (int i = 0; i < BUFFER_SIZE; i++)
    {
        buffer[i] = ((value >> (8 * i)) & 0XFF);
    }
}

// Convert byte array to long long
FastInteger BytesToFastInt(unsigned char* buffer, int BUFFER_SIZE)
{
    FastInteger recoveredValue = 0;
    for (int i = 0; i < BUFFER_SIZE; i++)
    {
        auto byteVal = (((FastInteger)buffer[i]) << (8 * i));
        recoveredValue |= byteVal;
    }
    return recoveredValue;
}

// mod b 
FastInteger SignedBarrettReducer::mod(FastInteger a) const {
    FastInteger abs_a = a > 0 ? a : -a;
    abs_a = abs_a - ((FastInteger) (((__int128_t)abs_a * m) >> 61 )) * b;
    return a >= 0 ? abs_a : -abs_a + b;
}

    // mod bsqr = b²
FastInteger SignedBarrettReducer::modsqr(FastInteger a) const {
    FastInteger abs_a = a > 0 ? a : -a;
    abs_a = abs_a - ((FastInteger) (((__int128_t)abs_a * msqr) >> 61 )) * bsqr;
    return a >= 0 ? abs_a : -abs_a + bsqr;
}

// mod 2*b
FastInteger SignedBarrettReducer::mod2(FastInteger a) const {
    FastInteger abs_a = a > 0 ? a : -a;
    abs_a = abs_a - ((FastInteger) (((__int128_t)abs_a * m2) >> 61 )) * (b << 1);
    return a >= 0 ? abs_a : -abs_a + (b << 1);
}

// mod 2*b²
FastInteger SignedBarrettReducer::mod2sqr(FastInteger a) const {
    FastInteger abs_a = a > 0 ? a : -a;
    abs_a = abs_a - ((FastInteger) (((__int128_t)abs_a * m2sqr) >> 61 )) * (bsqr << 1);
    return a >= 0 ? abs_a : -abs_a + (bsqr << 1);
}

FastInteger SignedBarrettReducer::mulmodsqr(FastInteger a, FastInteger c) const {
    __int128_t mul = (__int128_t) a * c ;
    __int128_t mul_abs = mul > 0 ? mul : -mul;
    mul_abs = (FastInteger) ( mul_abs - (( mul_abs * (__int128_t) msqr) >> 61 ) * (__int128_t) (bsqr) ); 
    return mul >= 0 ? mul_abs : - mul_abs + bsqr;
}


// -------------------------- Mat functions 

bool FastMat3::operator==(NTL::mat_ZZ &other) {

    return mat[0][0] == other[0][0] &&
           mat[0][1] == other[0][1] &&
           mat[0][2] == other[0][2] &&
           mat[1][0] == other[1][0] &&
           mat[1][1] == other[1][1] &&
           mat[1][2] == other[1][2] &&
           mat[2][0] == other[2][0] &&
           mat[2][1] == other[2][1] &&
           mat[2][2] == other[2][2]; 
}

bool FastMat3::operator==(FastMat3 &other) {

    return mat[0][0] == other[0][0] &&
           mat[0][1] == other[0][1] &&
           mat[0][2] == other[0][2] &&
           mat[1][0] == other[1][0] &&
           mat[1][1] == other[1][1] &&
           mat[1][2] == other[1][2] &&
           mat[2][0] == other[2][0] &&
           mat[2][1] == other[2][1] &&
           mat[2][2] == other[2][2]; 
}


// ------------------------- FastQuat Functions 

void FastQuat::normalize() {

        if (coeffs[4] != 1 && coeffs[4] != 2) {
            // usually those are the smallest coefficients
            FastInteger g = GCD(coeffs[2], coeffs[3]);
            if (g != 1) {
                g = GCD(g, coeffs[0]);
            }
            if (g != 1) {
                g = GCD(g, coeffs[1]);
            }
            if (g != 1) {
                g = GCD(g, coeffs[4]);
            }
    
            if (g != 1)
                for (auto &c: coeffs)
                    c /= g;
            } 
        
    }

// normalize by 2 or 4 if possible 
void FastQuat::normalize2() {
    if (!(IsOdd(coeffs[0]) || IsOdd(coeffs[1]) || IsOdd(coeffs[2]) || IsOdd(coeffs[3])) && !IsOdd(coeffs[4])) {
        coeffs[0] >>= 1;
        coeffs[1] >>= 1;
        coeffs[2] >>= 1;
        coeffs[3] >>= 1;
        coeffs[4] >>= 1;
        if (!(IsOdd(coeffs[0]) || IsOdd(coeffs[1]) || IsOdd(coeffs[2]) || IsOdd(coeffs[3])) && !IsOdd(coeffs[4])) {
            coeffs[0] >>= 1;
            coeffs[1] >>= 1;
            coeffs[2] >>= 1;
            coeffs[3] >>= 1;
            coeffs[4] >>= 1;
        }
    }
}

FastQuat::FastQuat(quat const &a, FastQuatAlg const &_alg) : alg{_alg} {

    this->coeffs[0] = convert(a[0]);
    this->coeffs[1] = convert(a[1]);
    this->coeffs[2] = convert(a[2]);
    this->coeffs[3] = convert(a[3]);
    this->coeffs[4] = convert(a[4]);
}

FastInteger FastQuat::integral_norm() const {
    return (coeffs[0] * coeffs[0] + alg.q * coeffs[1] * coeffs[1] + alg.p * ( coeffs[2] * coeffs[2] + alg.q * coeffs[3] * coeffs[3])) / (coeffs[4] * coeffs[4]);
}

FastInteger FastQuat::scalar_remove() {
    FastInteger g = 1;
    int index = 3;
    while (index >= 0 && coeffs[index] == 0) {
        index--;
    }
    if (index != -1) {
        g = coeffs[index];
        for (int i = 1; i < 4 && g != 1; i++) {
            if (coeffs[ (4 + index - i) & 3 ] != 0) {
                g = GCD(g, coeffs[(4 + index - i) & 3]);
            }
        }
    
        if (g != 1) {
            coeffs[0] = coeffs[0] / g;
            coeffs[1] = coeffs[1] / g;
            coeffs[2] = coeffs[2] / g;
            coeffs[3] = coeffs[3] / g;
        }
        
    
    }
    return g;
    
}

void FastQuat::scalar_mul(const FastInteger N) {
    this->coeffs[0] *= N;
    this->coeffs[1] *= N;
    this->coeffs[2] *= N;
    this->coeffs[3] *= N;
}

bool FastQuat::is_zero() const
    {
        auto const &[a,b,c,d,e] = coeffs;
        assert(e != 0);
        return a == 0 && b == 0 && c == 0 && d == 0;
    }

// if alpha can be divided by 2 in O0, do it
void divide_by_two_O0(FastQuat *quat_sol) {
    if ((IsOdd((*quat_sol)[0]) == IsOdd((*quat_sol)[3]) ) && ( IsOdd((*quat_sol)[1]) == IsOdd((*quat_sol)[2]) )) {
        (*quat_sol)[4] = 2;
    }
}

FastQuat FastQuat::operator*(FastQuat const &other) const
    {
        assert(alg.q == 1);
        assert (&other.alg == &alg); 

        auto const &[a,b,c,d,e] = this->coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;

        FastInteger r,s,t,u;
        // r = a*A - b*B - alg.p * (c*C + d*D);
        // s = a*B + A*b + alg.p * ( c*D - C*d);
        // t = a*C + A*c - b*D + B*d;
        // u = a*D + A*d + b*C - B*c;

        // this is 12 multiplication instead of 18
        FastInteger temp;
        FastInteger m1 = (a + b) * ( A + B);
        temp = (a - b) * (A - B);
        FastInteger m2 = (m1 - temp) >> 1;
        m1 = -(m1 + temp) >> 1;
        s = - (c + d) * (C + D);
        temp = - (c - d) * (C - D);
        r = (s + temp) >> 1;
        s = (s - temp) >> 1; 

        temp = (a + b + c + d) * (A + B + C + D);
        u = (a - b + c - d) * (A - B + C - D);

        t = ((u + temp) >> 1) + r + m1 - (b * D << 1);
        u = ((temp - u) >> 1) - m2 + s - (B * c << 1); 

        r = m1 + (a * A << 1) + alg.p * r;
        s = m2 + alg.p * (s + (c * D << 1));
        

        return {{r,s,t,u, e*E}, alg};
    }

void FastQuat::mul_left(FastQuat const &other) {
        auto const &[a,b,c,d,e] = other.coeffs;
        auto const [A,B,C,D,E] = this->coeffs;

        // this is 12 multiplication instead of 18
        FastInteger temp;
        FastInteger m1 = (a + b) * ( A + B);
        temp = (a - b) * (A - B);
        FastInteger m2 = (m1 - temp) >> 1;
        m1 = -(m1 + temp) >> 1;
        this->coeffs[1] = - (c + d) * (C + D);
        temp = - (c - d) * (C - D);
        this->coeffs[0] = (this->coeffs[1] + temp) >> 1;
        this->coeffs[1] = (this->coeffs[1] - temp) >> 1; 

        temp = (a + b + c + d) * (A + B + C + D);
        this->coeffs[3] = (a - b + c - d) * (A - B + C - D);

        this->coeffs[2] = ((this->coeffs[3] + temp) >> 1) + this->coeffs[0] + m1 - (b * D << 1);
        this->coeffs[3] = ((temp - this->coeffs[3]) >> 1) - m2 + this->coeffs[1] - (B * c << 1); 

        this->coeffs[0] = m1 + (a * A << 1) + alg.p * this->coeffs[0];
        this->coeffs[1] = m2 + alg.p * (this->coeffs[1] + (c * D << 1));

        this->coeffs[4] = e * E;
        

        // return {{r,s,t,u, e*E}, alg};
}


FastQuat FastQuat::operator+(FastQuat const &other) const
{
        if (&other.alg != &alg) // object identity
            throw;
        auto const &[a,b,c,d,e] = coeffs;
        auto const &[A,B,C,D,E] = other.coeffs;
        return {{a*E+A*e, b*E+B*e, c*E+C*e, d*E+D*e, e*E}, alg};
}

// multiplication of quaternions whose coefficients are a bit too big, so we use floating point the value by fac
std::array<FastFloat,5> float_mul(FastQuat alpha, FastQuat beta, FastInteger const fac)
    {
        assert(alpha.alg.q == 1);
        if (&alpha.alg != &beta.alg) // object identity
            throw;
        // auto const &[a,b,c,d,e] = alpha.coeffs;
        // auto const &[A,B,C,D,E] = beta.coeffs;
        FastInteger A = beta[0];
        FastInteger B = beta[1];
        FastInteger C = beta[2];
        FastInteger D = beta[3];
        FastInteger E = beta[4];
        FastInteger e = alpha[4]; 
        FastFloat fa = IntToFloat(alpha[0], fac);
        FastFloat fb = IntToFloat(alpha[1], fac);
        FastFloat fc = IntToFloat(alpha[2], fac);
        FastFloat fd = IntToFloat(alpha[3], fac);

        FastFloat r = fa*A - alpha.alg.q*fb*B - alpha.alg.p* (fc*C + alpha.alg.q*fd*D);
        FastFloat s = fa*B + A*fb + alpha.alg.p* ( fc*D - C*fd);
        FastFloat t = fa*C + A*fc - alpha.alg.q*fb*D + alpha.alg.q*B*fd;
        FastFloat u = fa*D + A*fd + fb*C - B*fc;
        // std::cout << e << " " << E << "\n";
        return {r,s,t,u, (FastFloat) E * e};
    }

quat FastQuat_to_quat(FastQuat alpha, const quatalg &Bp) {
    return {{Integer(alpha[0]), Integer(alpha[1]), Integer(alpha[2]), Integer(alpha[3]), Integer(alpha[4])}, Bp};
}


// ------------------------- FastQuatLat functions



FastQuatAlg::FastQuatAlg(const quatalg &A) {
    this->p = convert(A.p);
    this->q = convert(A.q);
}

FastQuatLat::FastQuatLat(FastMat4 const &gens, std::pair<FastInteger, FastInteger> const &norm, FastInteger const &denom_, FastQuatAlg const &alg_, bool do_normalize, bool good_type_) : alg{alg_}, basis{gens}, denom{denom_}, good_type{good_type_}
{

    if (do_normalize) {
        // this should not happen
        assert(0);
    }
    the_norm.first = norm.first;
    the_norm.second = norm.second;
}

FastQuatLat::FastQuatLat(quatlat const L, FastQuatAlg const &_alg) : alg{_alg} {

    this->basis[0][0] = convert(L.basis[0][0]);
    this->basis[0][1] = convert(L.basis[0][1]);
    this->basis[0][2] = convert(L.basis[0][2]);
    this->basis[0][3] = convert(L.basis[0][3]);
    this->basis[1][0] = convert(L.basis[1][0]);
    this->basis[1][1] = convert(L.basis[1][1]);
    this->basis[1][2] = convert(L.basis[1][2]);
    this->basis[1][3] = convert(L.basis[1][3]);
    this->basis[2][0] = convert(L.basis[2][0]);
    this->basis[2][1] = convert(L.basis[2][1]);
    this->basis[2][2] = convert(L.basis[2][2]);
    this->basis[2][3] = convert(L.basis[2][3]);
    this->basis[3][0] = convert(L.basis[3][0]);
    this->basis[3][1] = convert(L.basis[3][1]);
    this->basis[3][2] = convert(L.basis[3][2]);
    this->basis[3][3] = convert(L.basis[3][3]);
    this->denom = convert(L.denom);
    auto [a1,a2] = L.norm();
    this->the_norm = {convert(a1), convert(a2)};
    auto norm = this->the_norm.first / this->the_norm.second;
    if (norm != 1 && this->basis[1][0] == 0 && GCD(norm, this->basis[3][1]) == 1) {
        this->good_type = true;
    }
    else {
        this->good_type = false;
    }

}

bool FastQuatLat::operator==(FastQuatLat &other) {
    return this->alg.p == other.alg.p && 
    this->denom * other.basis[0][0] == other.denom * this->basis[0][0] &&
    this->denom * other.basis[0][1] == other.denom * this->basis[0][1] &&
    this->denom * other.basis[0][2] == other.denom * this->basis[0][2] &&
    this->denom * other.basis[0][3] == other.denom * this->basis[0][3] && 
    this->denom * other.basis[1][0] == other.denom * this->basis[1][0] &&
    this->denom * other.basis[1][1] == other.denom * this->basis[1][1] &&
    this->denom * other.basis[1][2] == other.denom * this->basis[1][2] &&
    this->denom * other.basis[1][3] == other.denom * this->basis[1][3] && 
    this->denom * other.basis[2][0] == other.denom * this->basis[2][0] &&
    this->denom * other.basis[2][1] == other.denom * this->basis[2][1] &&
    this->denom * other.basis[2][2] == other.denom * this->basis[2][2] &&
    this->denom * other.basis[2][3] == other.denom * this->basis[2][3] && 
    this->denom * other.basis[3][0] == other.denom * this->basis[3][0] &&
    this->denom * other.basis[3][1] == other.denom * this->basis[3][1] &&
    this->denom * other.basis[3][2] == other.denom * this->basis[3][2] &&
    this->denom * other.basis[3][3] == other.denom * this->basis[3][3];  
}


std::pair<FastInteger, FastInteger> FastQuatLat::norm() const
 {
        if ( the_norm.second == 0) {
            // this shouldn't happen
            assert(0);
        }
        return the_norm;
}


void FastQuatLat::_fast_intersect(FastQuatLat &other, const SignedBarrettReducer &red, const SignedBarrettReducer &red1, const FastInteger &inv, const FastInteger &inv2sqr) {

    // only work for p=3 mod 4
    assert(this->alg.p%4 == 3);

    // requires both basis of ideals to be in HNF
    assert( basis[0][1] == 0 && basis[2][3] == 0 && other.basis[0][1] == 0 && other.basis[2][3] == 0 );
    assert( denom == 2 && other.denom == 2);
    assert( other.norm().second == 1 && this->the_norm.second == 1);
    auto n1 = other.norm().first;
    auto n = this->the_norm.first;
    auto newnorm = red.b;

    if (n1 > 1) {
        basis[0][0] *= n1;
        if (basis[1][0] == 0 && other.basis[1][0] == 0) {
            basis[1][1] *= n1;
            assert(basis[2][0] >= 0);
            assert(other.basis[2][0] >=0);
            basis[2][0] = CRT(basis[2][0], n, other.basis[2][0], red, red1, inv);
            basis[2][1] = CRT(basis[2][1], n, other.basis[2][1], red ,red1, inv);


            if (good_type && other.good_type) {    
                good_type = true;
                basis[0][3] = CRT(basis[0][3], n, other.basis[0][3], red, red1, inv);
                // these ones are mod 2 n^2 n1^2
                basis[1][3] = CRT2sqr(basis[1][3], n * n, other.basis[1][3], red, red1, inv2sqr);
                basis[1][2] = CRT2sqr(basis[1][2], n * n, other.basis[1][2], red, red1, inv2sqr);
            }
            else {
                good_type = false;
                basis[3][1] = basis[2][0];
                basis[3][0] = - basis[2][1];
                basis[0][3] = 0;
                basis[0][2] = 0;
                basis[1][2] = 0;
                basis[1][3] = 0;
            }
            
            
        }
        else {
            basis[0][3] = 0;
            basis[0][2] = 0;
            basis[1][2] = 0;
            basis[1][3] = 0;
            basis[1][0] = CRT(basis[1][0], n, other.basis[1][0], red, red1, inv);
            basis[1][1] = CRT(basis[1][1], n, other.basis[1][1], red, red1, inv);
            auto mod = newnorm;
            // reducing
            auto a = GCD(mod,basis[1][1]);
            auto b = basis[1][1]/(2*a);
            basis[1][1] = 2*a;
            while (GCD(b,mod) != 1) {
                b += mod/a;
            }
            assert(GCD(b,mod) == 1);
            auto c = InvMod(red.mod(b), red);

            basis[1][0] = 2*((basis[1][0]/2 * c) % mod);
            basis[2][0] = CRT(basis[2][0], n, other.basis[2][0], red, red1, inv);
            basis[2][1] = CRT(basis[2][1], n, other.basis[2][1], red, red1, inv);
            basis[2][2] = CRT(basis[2][2], n, other.basis[2][2], red, red1, inv);

            // reducing
            basis[2][1] = (basis[2][1] - basis[2][2])/2;
            a = GCD(mod,basis[2][2]);
            b = basis[2][2]/(a);
            while (GCD(b,mod) != 1) {
                b += 2 * mod/a;
            }
            basis[2][2] = a;
            c = InvMod(red.mod(b), red);
            basis[2][0] = 2*((basis[2][0]/2 * c) % mod);
            basis[2][1] = basis[2][2] + 2*(((basis[2][1]) * c) % mod);
            // and reducing by previous vector
            FastInteger quotient = (basis[2][1] / basis[1][1] );
            basis[2][1] = basis[2][1] - basis[1][1] * quotient;
            basis[2][0] = (basis[2][0] - basis[1][0] * quotient) % (2 * mod);
            basis[3][0] = CRT(basis[3][0], n, other.basis[3][0], red, red1, inv);
            // assert(basis[3][0] % (2 * n1) == other.basis[3][0]);
            basis[3][1] = CRT(basis[3][1], n, (other.basis[3][1] + 2 * n1), red, red1, inv);
            // assert(basis[3][1] % (2 * n1) == other.basis[3][1] %(2 * n1));
            basis[3][2] = CRT(basis[3][2], n, other.basis[3][2], red, red1, inv);
            basis[3][3] = CRT(basis[3][3], n, other.basis[3][3], red, red1, inv);

            basis[3][0] = (basis[3][0] - basis[3][3])/2;
            a = GCD(mod,basis[3][3]);
            b = basis[3][3]/(a);
            basis[3][3] = a;
            c = InvMod(b % mod, red);
            basis[3][0] = basis[3][3] + 2*((basis[3][0] * c) % mod);
            basis[3][1] = 2*((((basis[3][1] -basis[3][2])/2) * c) % mod);
            basis[3][2] = (((basis[3][2]) * c) % mod);
            basis[3][1] += basis[3][2];
            // and reducing by previous vector
            quotient = (basis[3][2]/basis[2][2] );
            basis[3][2] = basis[3][2] - basis[2][2]*quotient;
            basis[3][1] = (basis[3][1] - basis[2][1]*quotient) % (2 * mod);
            basis[3][0] = (basis[3][0] - basis[2][0]*quotient) % (2 * mod);
            // and by previous vector
            quotient = (basis[3][1]/basis[1][1] );
            basis[3][1] = (basis[3][1] - basis[1][1]*quotient) % (2 * mod);
            basis[3][0] = (basis[3][0] - basis[1][0]*quotient) % (2 * mod);
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
    this->the_norm = {newnorm, 1};
}

FastQuatLat FastQuatLat::new_right_order() {
    
    quatalg Bp {Integer(this->alg.p), Integer(this->alg.q)};
    quatlat copy = this->FastQuatLat_to_quatlat(Bp); 
    return FastQuatLat(copy.new_right_order(), this->alg);
}



std::pair<FastQuatLat, FastMat3> FastQuatLat::fast_right_order_and_gram(const SignedBarrettReducer &red) {
    
    FastMat3 gram; 
    gram[0][0] = 0;

    FastInteger norm = red.b;

    assert(alg.p %4 == 3);
    assert(the_norm.first % the_norm.second == 0);

    FastInteger normsqr = red.bsqr;
    FastInteger newdenom = norm << 1 ;
    FastMat4 newbasis;

    // first element is 1
    newbasis[0][0] = newdenom;
    newbasis[0][1] = 0;
    newbasis[0][2] = 0;
    newbasis[0][3] = 0;
    newbasis[3][0] = norm;
    newbasis[3][1] = 0;
    newbasis[3][2] = 0;
    newbasis[2][0] = 0;
    newbasis[2][1] = 0;    
    newbasis[1][0] = 0;
    newbasis[1][1] = 1;

    // we assume some specific shape of the basis that allows to simplify the computation
    if (good_type) {
        
        assert(basis[3][2] == 0 && basis[3][3] == 1 && basis[2][2] == 1 && basis[2][3] == 0);
    
        // third element is of the form 1 * j + x * k
        newbasis[2][2] = newdenom;
        // last element is (1 + (N) k)/2
        newbasis[3][3] = normsqr;
    
        // FastInteger c = basis[2][0];
        FastInteger x;
        x = basis[0][3]; // precomputed value 
        if (x > norm) {
            x -= norm;
        }
        newbasis[2][3] = x * newdenom;

        newbasis[1][3] = red.mod2sqr(basis[1][3]);
        newbasis[1][2] = red.mod2sqr(basis[1][2]);

        gram[0][0] = Rounding( IntToFloat(this->alg.p, normsqr) * ((FastFloat) newbasis[1][2] * newbasis[1][2] + (FastFloat) newbasis[1][3] * newbasis[1][3] ) );
        gram[1][0] = Rounding(2 * IntToFloat(this->alg.p, norm) * (FastFloat) (newbasis[1][2] +  newbasis[1][3] * x));
        gram[0][1] = gram[1][0];
        gram[0][2] = this->alg.p * (newbasis[1][3]);
        gram[2][0] = gram[0][2];
        gram[1][2] = (this->alg.p * norm * x) << 1;
        gram[2][1] = gram[1][2];
        gram[1][1] = (this->alg.p * ( 1 + x * x )) << 2;
        gram[2][2] = this->alg.p * normsqr;


    }
    else {
        // testing if we are in some weird case
        FastInteger a = basis[3][0] * basis[3][0] + basis[3][1] * basis[3][1] - alg.p * (basis[3][2] * basis[3][2] + basis[3][3] * basis[3][3]);
        assert(a%2 ==0);
        a = a >> 1;
    
        // weird rare case for which we use NTL 
        if (GCD(a , norm) != 1) {
            quatalg Bp {Integer(this->alg.p), Integer(this->alg.q)};
            FastQuatLat O = this->new_right_order();
            quatlat slowO = O.FastQuatLat_to_quatlat(Bp);
            auto slowGram = compute_HNF_gram_order(&slowO);
            FastMat3 gram;
            gram[0][0] = convert(slowGram[0][0]);
            gram[0][1] = convert(slowGram[0][1]);
            gram[0][2] = convert(slowGram[0][2]);
            gram[1][0] = convert(slowGram[1][0]);
            gram[1][1] = convert(slowGram[1][1]);
            gram[1][2] = convert(slowGram[1][2]);
            gram[2][0] = convert(slowGram[2][0]);
            gram[2][1] = convert(slowGram[2][1]);
            gram[2][2] = convert(slowGram[2][2]);
            return {FastQuatLat(slowO, O.alg), gram};
        }

    
        FastInteger c = basis[3][1] * basis[2][2] - basis[2][1] * basis[3][2];
        // third element is of the form A * j + x * k
        // where A = GCD(c, norm)
        FastInteger A = GCD(c ,norm);
        newbasis[2][2] = newdenom * A;

        // last element is (1 + (N/A) k)/2
        newbasis[3][3] = normsqr / A;

        // computing the value of x
        FastInteger x;
        if (c == 0) {
            x = 0;
        }
        else {
            c = c / A;
            while (GCD(c,norm) != 1) {
                basis[2][0] += 2 * norm;
                basis[3][1] += 2 * norm;
                c = basis[3][1] * basis[2][2] - basis[2][1] * basis[3][2];
                assert(A == GCD(c,norm));
                c = c / A;
                a = basis[3][0] * basis[3][0] + basis[3][1] * basis[3][1] - alg.p * ( basis[3][2] * basis[3][2] + basis[3][3] * basis[3][3]);
                assert(a % 2 ==0);
                a = a / 2;
            }
            assert( GCD(c, norm) ==1 );
            c = InvMod(c, red);
            x = red.mod((basis[3][1] * basis[2][3] - basis[2][1] * basis[3][3]) * c);
        }

        newbasis[2][3] = x * newdenom;

        // 2nd element is of the form (i + y * j + z * k) / newdenom
        
        
        // let gen be the last basis element
        // a priori this is derived from Conjugate(gen) * i * gen / norm = (a * i + b * j + c * k) / newdenom
        // first we compute b,c (a was computed previously)
        newbasis[1][2] = (basis[3][1] - basis[3][0]) * (basis[3][2] + basis[3][3]);
        c = (basis[3][1] + basis[3][0]) * (basis[3][2] - basis[3][3]);
        FastInteger b = (newbasis[1][2] + c) >> 1;
        c = (newbasis[1][2] - c) >> 1;

        // and new we reduce it
        assert(GCD(a , norm) == 1);

        FastInteger d = InvModSqr(red.modsqr(a), red);

        newbasis[1][3] = red.mod2sqr(c * d);
        if (IsOdd(a) && IsOdd(d)) {
            newbasis[1][2] = (red.mod2sqr( b * d ));
        }
        else {
            newbasis[1][2] = (red.mod2sqr( b * d + normsqr ));
        }
        gram[0][0] = Rounding( IntToFloat(this->alg.p, normsqr) * ((FastFloat) newbasis[1][2] * newbasis[1][2] + (FastFloat) newbasis[1][3] * newbasis[1][3] ) );
        gram[1][0] = Rounding(2 * IntToFloat(this->alg.p, norm) * (FastFloat) (newbasis[1][2] * A +  newbasis[1][3] * x));
        gram[0][1] = gram[1][0];
        gram[0][2] = this->alg.p * newbasis[1][3] / A;
        gram[2][0] = gram[0][2];
        c = norm / A;
        gram[1][2] = 2 * this->alg.p * c * x;
        gram[2][1] = gram[1][2];
        gram[1][1] = 4 * this->alg.p * ( A * A + x * x );
        gram[2][2] = this->alg.p * c * c;
        
    }

    return {{newbasis, {1, 1}, newdenom, alg, false, false}, gram};
}


quatlat FastQuatLat::FastQuatLat_to_quatlat(const quatalg &Bp) {

    NTL::mat_ZZ mat;
    mat.SetDims(4,4);

    mat[0][0] = Integer(basis[0][0]);
    mat[0][1] = Integer(basis[0][1]);
    mat[0][2] = Integer(basis[0][2]);
    mat[0][3] = Integer(basis[0][3]);
    mat[1][0] = Integer(basis[1][0]);
    mat[1][1] = Integer(basis[1][1]);
    mat[1][2] = Integer(basis[1][2]);
    mat[1][3] = Integer(basis[1][3]);
    mat[2][0] = Integer(basis[2][0]);
    mat[2][1] = Integer(basis[2][1]);
    mat[2][2] = Integer(basis[2][2]);
    mat[2][3] = Integer(basis[2][3]);
    mat[3][0] = Integer(basis[3][0]);
    mat[3][1] = Integer(basis[3][1]);
    mat[3][2] = Integer(basis[3][2]);
    mat[3][3] = Integer(basis[3][3]);
    Integer denom = Integer(this->denom);
    std::pair<FastInteger, FastInteger> a = this->norm();
    std::pair<Integer,Integer> the_norm = {Integer(a.first), Integer(a.second)};

    return {mat, the_norm, denom, Bp, false};

}

FastQuatLat FastQuatLat::copy() const {
        return {basis, the_norm, denom, alg, false, good_type};
}


FastQuatLat create_from_generator_O0(const FastQuat &gen, const FastInteger &norm) {

        assert(gen.alg.p % 4 == 3);

        FastInteger newdenom = 2;
        FastMat4 newbasis;
        newbasis[0][0] = 2 * norm;
        newbasis[0][1] = 0;
        newbasis[0][2] = 0;
        newbasis[0][3] = 0;
        newbasis[1][2] = 0;
        newbasis[1][3] = 0;

        auto gen2 = gen[2] % norm;
        auto gen3 = gen[3] % norm;

        FastInteger c = (gen2 * gen2 + gen3 * gen3) % norm;

        bool good_type;
        
        if (GCD(norm, c) == 1) {

            newbasis[1][0] = 0;
            newbasis[1][1] = newbasis[0][0];
            newbasis[2][2] = 1;
            newbasis[2][3] = 0;
            newbasis[3][3] = 1;
            newbasis[3][2] = 0;

            FastInteger Inv = NTL::InvMod((2 * c) % norm, norm);

            newbasis[2][0] = 2 * ((2*(( (gen2 * gen[0]) % norm + (gen[1] * gen3) % norm) * Inv) % norm * NTL::InvMod(2, norm)) % norm);
            if (newbasis[2][0] < 0) {
                newbasis[2][0] += 2*norm;
            }
            
            newbasis[3][1] = newbasis[2][0];
            newbasis[2][1] = 1 + 2 * ((( (gen2 * gen[1]) % norm - (gen[0] * gen3)% norm - c) * Inv ) % norm);
            
            if (newbasis[2][1] < 0) {
                newbasis[2][1] += 2*norm;
            }
            
            newbasis[3][0] = 2*norm - newbasis[2][1];
            
            if (newbasis[3][0] < 0) {
                newbasis[3][0] += 2*norm;
            }

            if (newbasis[3][1] < 0) {
                newbasis[3][1] += 2*norm;
            }


            good_type = true;
        }
        else
        {

            // we are in one of the two split ideals
            // if the norm is not a prime we should be within the order to jinv enumeration and we are going to skip
            // this ideal anyway
            if (!NTL::ProbPrime(norm)) {
                newbasis[1][0] = 1;
                return {newbasis, {norm, 1}, newdenom, gen.alg, false, false};
            }

            // we compute one quaternion element in ZZ[omega] of norm divisible by norm
            // std::cout << norm << " " << norm  - gen.alg.q << "split ideal... \n";
            assert(NTL::Jacobi(Integer(norm  - gen.alg.q), Integer(norm)));
            
            FastInteger a;
            a = convert(NTL::SqrRootMod(Integer(norm  - gen.alg.q), Integer(norm)));

            
            FastQuat special_gen = FastQuat({ {a, 1, 0, 0, 0}, gen.alg });
            FastQuat prod = gen * special_gen;

            // if the result is on NO0 then the generator is the conjugate of special_gen
            if ( prod[0] % norm == 0 && prod[1] % norm == 0 && prod[2] % norm == 0, prod[3] % norm == 0 ) {
                newbasis[1][0] = 2*norm - 2*a;
                newbasis[1][1] = 2;
                newbasis[2][0] = 0;
                newbasis[2][1] = norm;
                newbasis[2][2] = norm;
                newbasis[2][3] = 0;
                newbasis[3][0] = 1;
                newbasis[3][1] = a;
                newbasis[3][2] = a;
                newbasis[3][3] = 1;
                good_type = false;
            }
            else {
                newbasis[1][0] = 2*a;
                newbasis[2][0] = 0;
                newbasis[2][1] = norm;
                newbasis[2][2] = norm;
                newbasis[2][3] = 0;
                newbasis[3][0] = 1;
                newbasis[3][1] = norm - a;
                newbasis[3][2] = norm - a;
                newbasis[3][3] = 1;
                good_type = false;
            }
        }

        return {newbasis, {norm, 1}, newdenom, gen.alg, false, good_type};

}


// enumerate throught the set of ell-ideals as ideals of the form * first_gen * (C + D *iter) + O0 * ell
// ell is prime and >= 5
std::vector<FastQuatLat> left_ideals_of_prime_norm_O0(FastInteger const &ell, const quatalg &Bp, const FastQuatAlg &fast_Bp) {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Enumerates through the set of ell-ideals of O_0 as ideals of the form * first_gen * (C + D *iter) + O0 * ell
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    assert(ell != fast_Bp.p);

    // first, we find the two elements we need, first gen and iter
    // performed with NTL
    bool found = false;
    quat first_gen = {{NTL::ZZ(0), NTL::ZZ(0), NTL::ZZ(1), NTL::ZZ(0), NTL::ZZ(1)}, Bp};
    while(!found) {
        first_gen[0] = (-NTL::power(first_gen[1],2) - Bp.p) % Integer(ell);
        assert(first_gen[0] >= 0);
        long legendre = NTL::Jacobi( first_gen[0], Integer(ell));
        if (legendre == 1) {
            found = true;
            first_gen[0] = NTL::SqrRootMod( first_gen[0], Integer(ell)); 
        }
        else {
            first_gen[1]++;
        }
        assert(first_gen[1] < ell); // we should ALWAYS find a solution before reaching that point
    }

    quatlat O0 = Bp.maximal_order(false);

    quatlat I = create_from_generator_O0( first_gen, Integer(ell));
    assert(std::get<0>(I.norm())/std::get<1>(I.norm()) == ell);

    std::list<int> ell_list = {(int) ell};
    // generating the iterating quaternion
    quat iter = find_quaternion_iterator(ell_list, I, O0, Bp);

    assert(NTL::ProbPrime(ell));

    std::vector<std::tuple<FastInteger,FastInteger>> coeff_list = {};
    std::vector<FastQuatLat> id_list = {};

    // create the list of coeff
    FastInteger iterate = 0;
    while (iterate < ell) {
        coeff_list.push_back({1, iterate});
        if (GCD(iterate, ell) != 1) {
            coeff_list.push_back({iterate, 1});
        }
        iterate++;
    }

    FastQuat fast_first_gen = FastQuat(first_gen, fast_Bp);
    FastQuat fast_iter = FastQuat(iter, fast_Bp);


    for (auto coeff : coeff_list) {
        FastQuat gen = (fast_first_gen * (FastQuat({{std::get<0>(coeff), 0, 0, 0, 1}, fast_first_gen.alg}) + fast_iter * std::get<1>(coeff) ));
        FastQuatLat J = create_from_generator_O0( gen, ell);
        id_list.push_back(J);
    }

    return id_list;
}


// ---------------------------------  Other functions 

bool commutatorfind(const std::pair<FastQuat,FastQuat> &beta1, const std::pair<FastQuat,FastQuat> &beta2, const std::pair<FastInteger,FastInteger> &n, const FastQuat &prod_quat, const FastInteger &prod_quat_norm, bool is_FP, FastQuat &quat_sol, const SignedBarrettReducer &redp) {
    //////////////////////////////////////////////////////////
    //// Finds alpha such that alpha * beta1.first = beta2.first * alpha
    //// and alpha * beta1.second = beta2.second * alpha
    //// Assumes that beta1, beta2 have denominator 2
    //////////////////////////////////////////////////////////
    assert( beta1.first.integral_norm() == beta2.first.integral_norm());
    assert( beta1.second.integral_norm() == beta2.second.integral_norm() );


    FastInteger nprod = n.first * n.second;
    FastInteger N21;

    assert(n.first == beta1.first[4]);
    assert(n.first == beta1.second[4]);
    assert(n.second == beta2.first[4]);
    assert(n.second == beta2.second[4]);

    FastInteger M00, M01, M02, M03, M10, M11, M12, M13, M20, M21, M22, M23;
    M11 = beta1.first[2] * n.second;
    M10 = beta1.first[3] * n.second;
    M13 =  beta2.first[1] * n.first; // tempy
    M03 = beta1.first[1] * n.second + M13;
    M02 =  beta2.first[3] * n.first; // tempt
    M22 = beta2.first[2] * n.first; // tempz
    M00 = M11 - M22 ;
    M01 = - (M10 + M02);
    M10 = M10 - M02;
    M11 = M11 + M22;

    M20 = beta1.second[2] * n.second;
    M22 = beta2.second[2] * n.first; // tempzeta
    M21 =  beta1.second[3] * n.second;
    M02 = beta2.second[3] * n.first; // temptau
    M20 = M20 - M22;
    M21 = - (M21 + M02);
    M12 = beta2.second[1] * n.first; // temprho
    M23 = M12  + beta1.second[1] * n.second;

    // FastFloat norm_quot = IntToFloat(n.second, n.first);
    // FastFloat inv_norm_quot = 1 / norm_quot;

    // compared to the basic version we will divide the ones coming from the first vector by fac1 and the one coming from the second vectors by fac2
    // TODO we could probably do better by doing the division during the matrix computation to avoid some unnecessary computations 
    FastFloat fM00, fM01, fM02, fM03, fM10, fM11, fM12, fM13, fM20, fM21, fM22, fM23;
    // fM22 = (FastFloat) beta2.first[2];
    // fM11 = (FastFloat) beta1.first[2] * norm_quot;
    // fM00 = fM11 - fM22;
    // fM11 = fM11 + fM22; 
    // fM10 = (FastFloat) beta1.first[3] * norm_quot;
    // fM02 = (FastFloat)  beta2.first[3]; // tempt
    // fM01 = - (fM10 + fM02);
    // fM10 = fM10 - fM02;
    // fM13 =  (FastFloat) beta2.first[1]; // tempy
    // fM03 = (FastFloat) beta1.first[1] * norm_quot + fM13;

    fM11 = IntToFloat(M11, n.first);
    fM10 = IntToFloat(M10, n.first);
    fM13 = IntToFloat(M13, n.first);
    fM03 = IntToFloat(M03, n.first);
    fM00 = IntToFloat(M00, n.first);
    fM01 = IntToFloat(M01, n.first);

    // fM20 = (FastFloat) beta1.second[2];
    // fM22 = (FastFloat) beta2.second[2] * inv_norm_quot; // tempzeta
    // fM21 = (FastFloat)  beta1.second[3];
    // fM02 = (FastFloat) beta2.second[3] * inv_norm_quot; // temptau
    // fM20 = fM20 - fM22;
    // fM21 = - (fM21 + fM02);
    // fM12 = (FastFloat) beta2.second[1] * inv_norm_quot; // temprho
    // fM23 = fM12 + (FastFloat) beta1.second[1];
    
    
    fM20 = IntToFloat(M20, n.second);
    fM21 = IntToFloat(M21, n.second);
    fM12 = IntToFloat(M12, n.second);
    fM23 = IntToFloat(M23, n.second); 
    fM22 = IntToFloat(M22, n.second);
    fM02 = IntToFloat(M02, n.second);

    // in most cases the solution is equal to  alpha2 * (j * beta1 - beta2 * j) + (j * beta1 - beta2 * j) * alpha1 / p
    // but it is possible that the sign choices were bad and in that case the solution is - alpha2 * (j * beta1 + beta2 * j) + (j * beta1 + beta2 * j) * alpha1 / p
    quat_sol[0] = Rounding(fM03 * fM21 - fM23 * fM01);
    if (is_FP && redp.mod(quat_sol[0]) == 0 ) {
        // quat_sol[0] =  ( (M03 - 2 * M13) * (M21 + 2 * M02) - ( (M23 - 2 * M12) * (-M10) ) ) / nprod ; 
        quat_sol[0] = Rounding( (fM03 - 2 * fM13) * (fM21 + 2 * fM02) - ( (fM23 - 2 * fM12) * (-fM10) ) ); 
        // quat_sol[1] =  (M11 * (M23 - 2 * M12) - (M03 - 2 * M13) * (M20 + 2 * M22) ) / nprod ;
        quat_sol[1] =  Rounding(fM11 * (fM23 - 2 * fM12) - (fM03 - 2 * fM13) * (fM20 + 2 * fM22) );
        // if we were very unlucky and in fact both first coordinates is divisible by p, we may have made a mistake and so we switch back to the other formula 
        if (redp.mod(quat_sol[0]) == 0 && redp.mod(quat_sol[1]) == 0){
            quat_sol[0] =  Rounding(fM03 * fM21 - fM23 * fM01); 
            // quat_sol[1] =  (M00 * M23 - M03 * M20) / nprod ;
            quat_sol[1] = Rounding(fM00 * fM23 - fM03 * fM20);
            // quat_sol[2] =  - (((M03 - 2 * M13 ) * M23) / beta1.first.alg.p +  M11 * M20 - M10 * M21) /nprod ;
            quat_sol[2] = Rounding( - (((fM03 - 2 * fM13 ) * fM23) / beta1.first.alg.p +  fM11 * fM20 - fM10 * fM21));
            // quat_sol[3] =  -(M21 * M00 -  M20 * M01) / nprod ;
            quat_sol[3] = Rounding(-(fM21 * fM00 -  fM20 * fM01));
        }
        else {
            // quat_sol[2] =  - (((M03) * (M23 - 2 * M12)) / beta1.first.alg.p +  M00 * (M20 + 2 * M22) + M01 * (M21 + 2 * M02) ) /nprod ;
            quat_sol[2] = Rounding(- (((fM03) * (fM23 - 2 * fM12)) / beta1.first.alg.p + fM00 * (fM20 + 2 * fM22) + fM01 * (fM21 + 2 * fM02) ));
            // quat_sol[3] =  -((M21 + 2* M02) * M11 -  (M20 + 2 * M22) * ((-M10))) / nprod ;
            quat_sol[3] = Rounding(-((fM21 + 2* fM02) * fM11 -  (fM20 + 2 * fM22) * ((-fM10))));
        }
    }
    else {
        // quat_sol[1] =  (M00 * M23 - M03 * M20) / nprod ;
        quat_sol[1] = Rounding(fM00 * fM23 - fM03 * fM20);
        // quat_sol[2] =  - (((M03 - 2 * M13 ) * M23) / beta1.first.alg.p +  M11 * M20 - M10 * M21) /nprod;
        quat_sol[2] =  Rounding(- (((fM03 - 2 * fM13 ) * fM23) / beta1.first.alg.p +  fM11 * fM20 - fM10 * fM21));
        // quat_sol[3] =  -(M21 * M00 -  M20 * M01) / nprod ;
        quat_sol[3] = Rounding(-(fM21 * fM00 -  fM20 * fM01));
    }
    N21 = quat_sol.scalar_remove(); // this contains the removed factor

    N21 = GCD(nprod, N21); // this is the factor that was potentially removed 

    bool boo = !quat_sol.is_zero();
    // in some unlucky case we may have chosen the wrong sign
    if (!boo) {
        // we try by changing the sign of the elements of the second basis
        if (is_FP) {
            // quat_sol[0] =  ( M03 * M21 - M23 * M01 ) / nprod ;
            // quat_sol[1] =  (M00 * M23 - M03 * M20) / nprod ;
            // quat_sol[2] =  - (((M03 - 2 * M13 ) * M23) / beta1.first.alg.p +  M11 * M20 - M10 * M21) /nprod ;
            // quat_sol[3] =  -(M21 * M00 -  M20 * M01) / nprod ;
            quat_sol[0] = Rounding(fM03 * fM21 - fM23 * fM01);
            quat_sol[1] = Rounding(fM00 * fM23 - fM03 * fM20);
            quat_sol[2] =  Rounding(- (((fM03 - 2 * fM13 ) * fM23) / beta1.first.alg.p +  fM11 * fM20 - fM10 * fM21));
            quat_sol[3] = Rounding(-(fM21 * fM00 -  fM20 * fM01));
        }
        else {
            // quat_sol[0] =  ( (M03 - 2 * M13) * (M21 + 2 * M02) - ( (M23 - 2 * M12) * (-M10) ) ) / nprod ; // 
            // quat_sol[1] =  (M11 * (M23 - 2 * M12) - (M03 - 2 * M13) * (M20 + 2 * M22) ) / nprod ;
            // quat_sol[2] =  - (((M03) * (M23 - 2 * M12)) / beta1.first.alg.p +  M00 * (M20 + 2 * M22) + M01 * (M21 + 2 * M02) ) /nprod ;
            // quat_sol[3] =  -((M21 + 2* M02) * M11 -  (M20 + 2 * M22) * ((-M10))) / nprod ;
            quat_sol[0] = Rounding( (fM03 - 2 * fM13) * (fM21 + 2 * fM02) - ( (fM23 - 2 * fM12) * (-fM10) ) ); 
            quat_sol[1] =  Rounding(fM11 * (fM23 - 2 * fM12) - (fM03 - 2 * fM13) * (fM20 + 2 * fM22) );
            quat_sol[2] = Rounding(- (((fM03) * (fM23 - 2 * fM12)) / beta1.first.alg.p + fM00 * (fM20 + 2 * fM22) + fM01 * (fM21 + 2 * fM02) ));
            quat_sol[3] = Rounding(-((fM21 + 2* fM02) * fM11 -  (fM20 + 2 * fM22) * ((-fM10))));
        }
        
        N21 = quat_sol.scalar_remove(); // this contains the removed factor
        N21 = GCD(nprod, N21); // this is the factor that was potentially removed 
        boo = !quat_sol.is_zero();

        // if it is sill zero, we try alpha2 * (k * beta1 - beta2 * k) + (k * beta1 - beta2 * k) * alpha1 / p
        if (!boo) {
            // // TODO this could be greatly improved, but at the same time, it is very very rare
            quatalg Bp {Integer(quat_sol.alg.p), Integer(quat_sol.alg.q)};
            quat kkk = {{Integer(0), Integer(0), Integer(0), Integer(1), Integer(1)}, Bp };
            quat a1 = FastQuat_to_quat(beta1.first, Bp);
            quat b1 = FastQuat_to_quat(beta1.second, Bp);
            quat a2 = FastQuat_to_quat(beta2.first, Bp);
            quat b2 = FastQuat_to_quat(beta2.second, Bp);

            quat slow_quat_sol = a2 * (kkk * b1 + b2 * kkk * Integer(-1)) + (kkk * b1 + b2 * kkk * Integer(-1)) * a1;
            slow_quat_sol[4] = 1;
            N21 = convert(slow_quat_sol.scalar_remove()); 
            quat_sol = FastQuat(slow_quat_sol, quat_sol.alg);
            N21 = GCD(nprod, N21); 
            boo = !quat_sol.is_zero();

            // if it is sill zero, we try alpha2 * (k * beta1 + beta2 * k) - (k * beta1 + beta2 * k) * alpha1 / p
            if (!boo) {
                // TODO this could be greatly improved, but at the same time, it is very very rare
                quatalg Bp {Integer(quat_sol.alg.p), Integer(quat_sol.alg.q)};
                quat kkk = {{Integer(0), Integer(0), Integer(0), Integer(1), Integer(1)}, Bp };
                quat a1 = FastQuat_to_quat(beta1.first, Bp);
                quat b1 = FastQuat_to_quat(beta1.second, Bp);
                quat a2 = FastQuat_to_quat(beta2.first, Bp);
                quat b2 = FastQuat_to_quat(beta2.second, Bp);
            
                quat slow_quat_sol = a2 * (kkk * b1 + b2 * kkk ) * (-1) + (kkk * b1 + b2 * kkk ) * a1;
                slow_quat_sol[4] = 1;
                N21 = convert(slow_quat_sol.scalar_remove()); // this contains the removed factor
                quat_sol = FastQuat(slow_quat_sol, quat_sol.alg);
            
                N21 = GCD(nprod, N21); // this is the factor that was potentially removed 
                boo = !quat_sol.is_zero();

                if (!boo) {
                    // shouldn't happen
                    std::cout << is_FP << "\n";
                    std::cout << "alpha1 = " << beta1.first << "\n";
                    std::cout << "beta1  = " << beta1.second << "\n";
                    std::cout << "alpha2 = " << beta2.first << "\n";
                    std::cout << "beta2  = " << beta2.second << "\n";
                    assert(0);
                }
            }
        }
    }

    // divide by two if we can
    divide_by_two_O0(&quat_sol);

    FastInteger norm = quat_sol.integral_norm();
    norm = (nprod * beta1.first.alg.p * prod_quat_norm) / norm;
    // we have nprod * p * prod_quat_norm / norm = p^bp * (prod_quat_norm)^b * ell^2 where bp,b are in {0,1} and ell divides N21 
    // we need to find the values of bp and b
    // there is a potential problem to find the value of the bit b when prod_quat_norm is a square having a common factor with N21 

    bool bp = (redp.mod(norm) == 0); // if bp is true, then p did not divide the norm of quat_sol
    if (bp) {
        norm = norm / beta1.first.alg.p;
    }
    bool isSquareProd = IsIntegralSquare(prod_quat_norm);
    bool isSquare = IsIntegralSquare(norm);

    // if prod_quat_norm does not divide norm, then b must be 0
    // if prod_quat_norm is not a square and norm / p^bp is a square, then b must be zero
    if ( GCD(prod_quat_norm, norm) == prod_quat_norm && (isSquareProd || !isSquare) ) {
        
        // we know that norm is a square
        // if prod_quat_norm is one, then b must be 1 (this is the Fp case)
        // if prod_quat_norm does not divide N21, prod_quat_norm cannot divide ell, and so b must be 1
        // if norm is not a square, then b must be 1
        if (prod_quat_norm != 1 && N21 % prod_quat_norm == 0 && isSquare) {

            // we still don't know the value of b, we must try to multiply by prod_quat to see what happens            
            // FastQuat test_sol = prod_quat * quat_sol;
            FastQuat test_sol = quat_sol;

            if (bp) {
                test_sol[4] *= prod_quat_norm; 
            } else {
                test_sol[4] *= prod_quat_norm * beta1.first.alg.p;
            }
            test_sol.normalize();

            if (test_sol[4] == 1 or test_sol[4] == 2) {
                // in that case we are going to assume that b = 0 
                quat_sol = test_sol;
                // it remains to divide by ell if we can

                quat_sol.scalar_mul(SqrRoot(norm));
            
                return !bp;
            }
            else {
                // in that case we know b must be 1 
                // it only remains to multiply by ell
                quat_sol.scalar_mul(SqrRoot(norm / (prod_quat_norm)));
                return bp;
            }

        }
        else {
            // b must be 1 this means that quat_sol is the solution we were looking for 
            // and we only need to renormalize by ell
            // unless the j-invariant is in Fp, then we are in a special case, and if bp is 0, then we need to make a correction 
            if (is_FP && !bp) {
                quat_sol.mul_left( prod_quat );
                nprod = SqrRoot(norm / (prod_quat_norm));
                quat_sol[0] *= (nprod);
                quat_sol[1] *= (nprod);
                quat_sol[2] *= (nprod);
                quat_sol[3] *= (nprod);
                nprod = (IsOdd(prod_quat[4]) ? prod_quat[4] : (prod_quat[4] >> 1));
                quat_sol[4] /= (nprod);
                nprod *= quat_sol.alg.p;
                quat_sol[0] /= (nprod);
                quat_sol[1] /= (nprod);
                quat_sol[2] /= (nprod);
                quat_sol[3] /= (nprod);
                quat_sol.normalize2();
                
            }
            else {
                quat_sol.scalar_mul(SqrRoot(norm / (prod_quat_norm)));
            }

            return bp || is_FP;
        }  
        
    }
    else {
        // b must be 0 

        // old way
        // quat_sol = prod_quat * quat_sol;
        // if (bp) {
        //     quat_sol[4] *= prod_quat_norm;
        // } else {
        //     quat_sol[4] *= prod_quat_norm * beta1.first.alg.p;
        // }
        // quat_sol.scalar_mul(SqrRoot(norm));

        quat_sol.mul_left( prod_quat );
        nprod = SqrRoot(norm);
        quat_sol[0] *= (nprod);
        quat_sol[1] *= (nprod);
        quat_sol[2] *= (nprod);
        quat_sol[3] *= (nprod);
        nprod = (IsOdd(prod_quat[4]) ? prod_quat[4] : (prod_quat[4] >> 1));
        quat_sol[4] /= (nprod);
        if (bp) {
            nprod *= prod_quat_norm;
        }
        else {
            nprod *= prod_quat_norm * quat_sol.alg.p;
        }
        
        quat_sol[0] /= (nprod);
        quat_sol[1] /= (nprod);
        quat_sol[2] /= (nprod);
        quat_sol[3] /= (nprod);
        quat_sol.normalize2();

        return !bp;
    }

    return bp;

}
