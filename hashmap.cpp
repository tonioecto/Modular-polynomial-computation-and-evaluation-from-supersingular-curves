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
#include <bit>
#include <algorithm>
#include <cstdint>
#include <stdexcept>

// SLOW LLL for low dimension (from Low-Dimensional Lattice Basis Reduction Revisited paper by Stehlé and Nguyen)

Integer DivRound(Integer a, Integer b) {
    Integer x;
    NTL::div(x,a + b/2,b);
    
    return x;     
}

void CVP2(NTL::mat_ZZ &gram, const int i, const int j, int k) { 
    Integer x = DivRound(gram[i][j],gram[i][i]);

    Integer y = gram[i][j] - x * gram[i][i];

    gram[j][j] = gram[j][j] - x * (y + gram[i][j]);

    gram[i][j] = y;
    gram[j][i] = y;
    gram[j][k] = gram[j][k] - x * gram[i][k];
    gram[k][j] = gram[j][k];

}

Integer NormSubGram(Integer x1, Integer x2, NTL::mat_ZZ gram, int i, int j, int k) {
    return gram[k][k] + x1 * x1 * gram[i][i] + x2 * x2 * gram[j][j] + 2 * x1 * x2 * gram[i][j] - 2 * ( x1 * gram[i][k] + x2 * gram[j][k] ); 
}


void CVP3(NTL::mat_ZZ &gram, const int i, const int j, const int k) {

    // std::cout << "i " << i << " j " << j << " k " << k << "\n";

    Integer t1 = gram[i][i] * gram[j][j] - gram[i][j] * gram[i][j];
    assert(t1 >= 0);
    // std::cout << "     t1 = " << t1 << "\n"; 
    // std::cout << "x1 * t1 = " << gram[j][j] * gram[i][k] - gram[i][j] * gram[j][k] << "\n";
    Integer x1 = DivRound( (gram[j][j] * gram[i][k] - gram[i][j] * gram[j][k]), t1);
    Integer x2 = DivRound( (gram[i][i] * gram[k][j] - gram[i][k] * gram[i][j]), t1);


    // std::cout << "x12 = " << x1 << " " << x2 << " \n";

    Integer smallest = NormSubGram(x1, x2, gram, i, j ,k);
    Integer temp = smallest;
    // std::cout << "first guess = " << smallest << "\n";
    Integer new_x1 = x1;
    Integer new_x2 = x2;
    
    // std::vector< std::pair<Integer, Integer>> candidates = { {Integer(-1), Integer(-1)}, {Integer(-1), Integer(0)}, {Integer(-1), Integer(1)}, {Integer(0), Integer(-1)}, {Integer(0), Integer(1)}, {Integer(1), Integer(-1)}, {Integer(1), Integer(0)}, {Integer(1), Integer(1)} };

    // Integer new_x1 = x1;
    // Integer new_x2 = x2;

    // for (auto [eps1, eps2] : candidates) {
    //     Integer new_norm = NormSubGram(x1 + eps1, x2 + eps2, gram, i, j, k);
    //     assert(new_norm > 0);
    //     if (new_norm < smallest) {
    //         std::cout<< eps1 << " " << eps2 << "\n";
    //         smallest = new_norm;
    //         new_x1 = x1 + eps1;
    //         new_x2 = x2 + eps2;
    //     }
    // } 

    // apparently the only admissible corrections are only (0,1), (0,-1), (1,0) and (-1,0) 
    // we can determine the sign by the cross product and so we have only two values to test
    Integer s1 = x1 * gram[i][i] + x2 * gram[i][j] - gram[i][k];
    Integer s2 = x2 * gram[j][j] + x1 * gram[i][j] - gram[j][k];
    if (s1 >=0) {
        s1 = temp + gram[i][i] - 2 * s1;
        if (s1 < smallest) {
            smallest = s1;
            new_x1 = x1 + Integer(-1);
        }
    }
    else {
        s1 = temp + gram[i][i] + 2 * s1;
        if (s1 < smallest) {
            smallest = s1;
            new_x1 = x1 + Integer(1);
        }
    }
    if (s2 >=0) {
        s2 = temp + gram[j][j] - 2 * s2;
        if (s2 < smallest) {
            smallest = s2;
            new_x2 = x2 + Integer(-1);
            new_x1 = x1;
        }
    }
    else {
        s2 = temp + gram[j][j] + 2 * s2;
        if (s2 < smallest) {
            smallest = s2;
            new_x2 = x2 + Integer(1);
            new_x1 = x1;
        }
    }

    // std::cout << "2nd guess " << smallest << "\n";

    gram[k][k] = smallest;
    gram[i][k] = gram[i][k] - new_x1 * gram[i][i] - new_x2 * gram[i][j];
    gram[k][i] = gram[i][k];
    gram[k][j] = gram[k][j] - new_x1 * gram[i][j] - new_x2 * gram[j][j];
    gram[j][k] = gram[k][j];    
}

void SortGram(NTL::mat_ZZ &gram) {
    if (gram[0][0] > gram[1][1]) {
            std::swap(gram[0][0], gram[1][1]);
            std::swap(gram[2][1], gram[2][0]);
            std::swap(gram[0][2], gram[1][2]);
        }
        if (gram[1][1] > gram[2][2]) {
            std::swap(gram[2][2], gram[1][1]);
            std::swap(gram[0][1], gram[0][2]);
            std::swap(gram[1][0], gram[2][0]);
        }
        if (gram[0][0] > gram[1][1]) {
            std::swap(gram[0][0], gram[1][1]);
            std::swap(gram[2][1], gram[2][0]);
            std::swap(gram[0][2], gram[1][2]);
        }
}

void GreedyReduction3(NTL::mat_ZZ &gram) {

    // std::cout << "Greedy \n";
    // std::cout << "gram= " << gram << "\n";

    

    // first we sort by increasing order
    SortGram(gram);

    int k = 1;
    int count = 0;
    while (k < 3 && count <= 100) {
        
        count++;
        // if (count  < 5) {
            // std::cout << count << " k = " << k << " " << gram << "\n";
        // }
        // std::cout << "k = " << k << "\n";
        // std::cout << gram << "\n";

        assert(gram[0][0] * gram[1][1] - gram[0][1] * gram[0][1] >= 0);
        assert(gram[0][0] * gram[2][2] - gram[0][2] * gram[0][2] >= 0);
        assert(gram[2][2] * gram[1][1] - gram[2][1] * gram[2][1] >= 0);
        
        if (k == 1) {
            Integer old = gram[1][1];
            // std::cout << "before CVP2" << gram << "\n";
            CVP2(gram, 0, 1, 2);
            // std::cout << "after CVP2" << gram << "\n";
// #ifndef NDEBUG 
            // if (old >= gram[1][])
// #endif 
            assert(old >= gram[1][1]);

            if (gram[1][1] >= gram[0][0]) {
                k = k + 1;      
            }
            else {
                // we swap the first and second vector
                std::swap(gram[0][0], gram[1][1]);
                std::swap(gram[2][1], gram[2][0]);
                std::swap(gram[0][2], gram[1][2]);
            }
        }
        else if (k == 2) {
            Integer old = gram[2][2];
            CVP3(gram, 0, 1, 2);
            assert(gram[2][2] <= old);

            if (gram[2][2] >= gram[1][1]) {
                k = k + 1;
            }
            else {
                // maybe this is a bit unefficient as we know the first and second vectors are already sorted
                SortGram(gram);
                k = 1;
            }

        } 


    }

} 

FastInteger div_round(FastInteger a, FastInteger b) { return (a^b)<0 ? (a - (b >> 1) )/b : (a + (b >> 1))/b; }


void FastCVP2(FastMat3 &gram, const size_t i, const size_t j, const size_t k, FastMat3 &Coords) { 
    
    FastInteger x = div_round(gram[i][j], gram[i][i]);
    gram[i][j] -= x * gram[i][i];

    // the j-th vector b[j] <- b[j] - x * b[i]
    gram[j][j] = gram[j][j] - x * (gram[i][j] + gram[j][i]);

    gram[j][i] = gram[i][j];
    gram[j][k] = gram[j][k] - x * gram[i][k];
    gram[k][j] = gram[j][k];

    Coords[j][0] -= x * Coords[i][0];
    Coords[j][1] -= x * Coords[i][1];
    Coords[j][2] -= x * Coords[i][2];
}


// this is the number of zeros, so the real bit number is #sizeof(FastInteger) - bit_length(x) 
// assumes the input is positive
size_t bit_length(FastInteger x) {
    // if (x >= 0) {
        return std::countl_zero((unsigned long long) x);
    // }
    // else {
        // return bit_length(-x - 1);
    // }
}

FastInteger absolute_value(FastInteger x) {
    uint64_t temp = x >> 63;     // make a mask of the sign bit
    x ^= temp;                   // toggle the bits if value is negative
    x += temp & 1;  
    return x;
}


void binLagrange(FastMat3 &gram, size_t &i, size_t &j, size_t &k, FastMat3 &Coords) {
    int bitlen = bit_length(gram[i][i]);

    while (true) {

        if (gram[j][j] < gram[i][i]) {
            gram[j][i] = gram[i][j];
            std::swap(i,j); 
            bitlen = bit_length(gram[i][i]);
        }

        if (gram[i][j] >= 0) {

            if ((gram[i][j] << 1) <= gram[i][i]) {
                gram[j][i]  = gram[i][j];
                gram[k][j] = gram[j][k];
                gram[k][i] = gram[i][k];
                return ;
            }
            size_t s = std::max(0, bitlen - (int) bit_length(gram[i][j]) );

            // b[j] <- b[j] - 2^s b[i]
            Coords[j][0] -= Coords[i][0] << s;
            Coords[j][1] -= Coords[i][1] << s;
            Coords[j][2] -= Coords[i][2] << s;
            // TODO there may be better to do, in particular, some part of this computation is already done to check the return condition
            gram[j][j] += (gram[i][i] << (s << 1)) - (gram[i][j] << (s + 1));
            gram[i][j] -= (gram[i][i] << s);
            gram[j][k] -= gram[i][k] << s;
        }
        else {


            if ((-gram[i][j] << 1) <= gram[i][i]) {
                gram[j][i]  = gram[i][j];
                gram[k][j] = gram[j][k];
                gram[k][i] = gram[i][k];
                return ;
            }
            size_t s = std::max(0, bitlen - (int) bit_length(- gram[i][j] - 1) );


           // b[j] <- b[j] + 2^s b[i]
            Coords[j][0] += Coords[i][0] << s;
            Coords[j][1] += Coords[i][1] << s;
            Coords[j][2] += Coords[i][2] << s;

            gram[j][j] += (gram[i][i] << (s << 1)) + (gram[i][j] << (s + 1));
            gram[i][j] += (gram[i][i] << s);
            gram[j][k] += gram[i][k] << s;
            
        }

    }
    
}

void FastLagrange2(FastMat3 &gram, size_t &i, size_t &j, size_t &k, FastMat3 &Coords) {

    FastCVP2(gram, i, j, k, Coords);

    while (gram[j][j] < gram[i][i]) {
        std::swap(i, j);
        // if ((gram[i][j] << 1) > gram[i][i]) {
            FastCVP2(gram, i, j, k, Coords);
        // }
        // else {return;}
        
    }

}

void FastCVP3(FastMat3 &gram, const size_t i, const size_t j, const size_t k, FastMat3 &Coords) {

    // FastInteger t1 = gram[i][i] * gram[j][j] - gram[i][j] * gram[i][j];
    // FastInteger x1 = FastDivRound( (gram[j][j] * gram[i][k] - gram[i][j] * gram[j][k]), t1);
    // FastInteger x2 = FastDivRound( (gram[i][i] * gram[k][j] - gram[i][k] * gram[i][j]), t1);

    FastInteger x1,x2;
    // std::cout << gram << "\n";

    if (gram[i][j] == 0 ) {
        assert(gram[i][i] != 0);
        x1 = div_round(gram[i][k], gram[i][i]);
        x2 = div_round(gram[k][j], gram[j][j]);
    }
    else {
        FastFloat t1 = ((FastFloat) gram[i][i]) * IntToFloat(gram[j][j], gram[i][j]) - (FastFloat) gram[i][j];
        x1 = Rounding( (((FastFloat) gram[j][j]) * IntToFloat(gram[i][k], gram[i][j]) - (FastFloat) gram[k][j]) / t1 );

        // the idea below doesn't seem to work
        // auto tx1 = round_div_branchless( round_div_branchless(gram[i][k]<< 32,gram[i][j]) - round_div_branchless(gram[k][j] << 32,gram[j][j]) ,  round_div_branchless(gram[i][i] << 32,gram[i][j]) - round_div_branchless(gram[i][j] << 32,gram[j][j]) );

        x2 = Rounding( (((FastFloat) gram[i][i]) * IntToFloat(gram[j][k], gram[i][j]) - (FastFloat) gram[i][k]) / t1 );  
    }

    // some common computations
    gram[k][i] = x1 * gram[i][i];
    gram[k][j] = x2 * gram[j][j];
    gram[j][i] = x2 * gram[i][j];
    gram[k][k] += x1 * gram[k][i] + x2 * gram[k][j] + ((x1 * gram[j][i]) << 1) - (( x1 * gram[i][k] + x2 * gram[j][k] ) << 1); 

    // apparently the only admissible corrections are only (0,1), (0,-1), (1,0) and (-1,0) 
    // we can determine the sign by the cross product and so we have only two values to test

    // we can already precompute those values
    gram[i][k] -= gram[k][i] + gram[j][i];
    gram[j][k] -= gram[k][j] + x1 * gram[i][j];
    gram[j][i] = gram[i][j];

    if (gram[i][k] < 0) {
        gram[k][i] = gram[i][i] + (gram[i][k] << 1);
        if (gram[k][i] < 0) {
            if (gram[j][k] < 0) {
                gram[k][j] = gram[j][j] + (gram[j][k] << 1);
                if (gram[k][j] < gram[k][i]) {
                    gram[k][k] += gram[k][j];
                    x2--;
                    gram[i][k] += gram[i][j];
                    gram[j][k] += gram[j][j];
                }
                else {
                    gram[k][k] += gram[k][i];
                    x1--;
                    gram[i][k] += gram[i][i];
                    gram[j][k] += gram[i][j];
                }
            }
            else {
                gram[k][j] = gram[j][j] - (gram[j][k] << 1);
                if (gram[k][j] < gram[k][i]) {
                    gram[k][k] += gram[k][j];
                    x2++;
                    gram[i][k] -= gram[i][j];
                    gram[j][k] -= gram[j][j];
                }
                else {
                    gram[k][k] += gram[k][i];
                    x1--;
                    gram[i][k] += gram[i][i];
                    gram[j][k] += gram[i][j];
                }
            }
        }
        else {
            if (gram[j][k] < 0) {
                gram[k][j] = gram[j][j] + (gram[j][k] << 1);
                if (gram[k][j] < 0) {
                    gram[k][k] += gram[k][j];
                    x2--;
                    gram[i][k] += gram[i][j];
                    gram[j][k] += gram[j][j];
                }
            }
            else {
                gram[k][j] = gram[j][j] - (gram[j][k] << 1);
                if (gram[k][j] < 0) {
                    gram[k][k] += gram[k][j];
                    x2++;
                    gram[i][k] -= gram[i][j];
                    gram[j][k] -= gram[j][j];
                }
            }
        }
        
    }
    else {
        gram[k][i] = gram[i][i] - (gram[i][k] << 1);
        if (gram[k][i] < 0) {
            if (gram[j][k] < 0) {
                gram[k][j] = gram[j][j] + (gram[j][k] << 1);
                if (gram[k][j] < gram[k][i]) {
                    gram[k][k] += gram[k][j];
                    x2--;
                    gram[i][k] += gram[i][j];
                    gram[j][k] += gram[j][j];
                }
                else {
                    gram[k][k] += gram[k][i]; 
                    x1++;
                    gram[i][k] -= gram[i][i];
                    gram[j][k] -= gram[i][j];
                }
            }
            else {
                gram[k][j] = gram[j][j] - (gram[j][k] << 1);
                if (gram[k][j] < gram[k][i]) {
                    gram[k][k] += gram[k][j];
                    x2++;
                    gram[i][k] -= gram[i][j];
                    gram[j][k] -= gram[j][j];
                }
                else {
                    gram[k][k] += gram[k][i];
                    x1++;
                    gram[i][k] -= gram[i][i];
                    gram[j][k] -= gram[i][j];
                }
            }
        }
        else {
            if (gram[j][k] < 0) {
                gram[k][j] = gram[j][j] + (gram[j][k] << 1);
                if (gram[k][j] < 0) {
                    gram[k][k] += gram[k][j];
                    x2--;
                    gram[i][k] += gram[i][j];
                    gram[j][k] += gram[j][j];
                }
            }
            else {
                gram[k][j] = gram[j][j] - (gram[j][k] << 1);
                if (gram[k][j] < 0) {
                    gram[k][k] += gram[k][j];
                    x2++;
                    gram[i][k] -= gram[i][j];
                    gram[j][k] -= gram[j][j];
                }
            }
        }
        
    }
    
    gram[k][i] = gram[i][k];
    gram[k][j] = gram[j][k];    

    // b[k] <- b[k] - new_x1 * b[i] - new_x2 * b[j];
    Coords[k][0] -= x1 * Coords[i][0] + x2 * Coords[j][0];
    Coords[k][1] -= x1 * Coords[i][1] + x2 * Coords[j][1];
    Coords[k][2] -= x1 * Coords[i][2] + x2 * Coords[j][2];
}

void FastSortGram(FastMat3 &gram, FastMat3 &Coords) {
    if (gram[0][0] > gram[1][1]) {
            std::swap(gram[0][0], gram[1][1]);
            std::swap(gram[2][1], gram[2][0]);
            std::swap(gram[0][2], gram[1][2]);
            std::swap(Coords[0][0], Coords[1][0]);
            std::swap(Coords[0][1], Coords[1][1]);
            std::swap(Coords[0][2], Coords[1][2]);

        }
        if (gram[1][1] > gram[2][2]) {
            std::swap(gram[2][2], gram[1][1]);
            std::swap(gram[0][1], gram[0][2]);
            std::swap(gram[1][0], gram[2][0]);
            std::swap(Coords[2][0], Coords[1][0]);
            std::swap(Coords[2][1], Coords[1][1]);
            std::swap(Coords[2][2], Coords[1][2]);
        }
        if (gram[0][0] > gram[1][1]) {
            std::swap(gram[0][0], gram[1][1]);
            std::swap(gram[2][1], gram[2][0]);
            std::swap(gram[0][2], gram[1][2]);
            std::swap(Coords[0][0], Coords[1][0]);
            std::swap(Coords[0][1], Coords[1][1]);
            std::swap(Coords[0][2], Coords[1][2]);
        }
}

void FastSortIndex(FastMat3 &gram, size_t &i, size_t &j, size_t &k) {
    if (gram[i][i] > gram[j][j]) {
        std::swap(i,j);
    }
    if (gram[j][j] > gram[k][k]) {
        std::swap(j,k);
    }
    if (gram[i][i] > gram[j][j]) {
        std::swap(i,j);
    }
}


// this is also known as Semaev algorithm
void FastGreedyReduction3(FastMat3 &gram, FastMat3 &Coords) {

    // std::cout << "\n\n new_greedy \n";
    // FastSortGram(gram, Coords);

    int current_dim = 1;
    int count = 0;

    size_t i,j,k;
    i = 0;
    j = 1; 
    k = 2;

    FastSortIndex(gram, i, j, k);

    while (current_dim < 3 && count <= 1000) {
        count++;

        if (current_dim == 1) {
            // auto gg = gram;
            // auto cc = Coords;
            // clock_t t = tic();
            FastLagrange2(gram, i, j, k, Coords);
            // t = tic() - t;
            // std::cout << "FL2 " << t << "\n"; 
            // // it maybe slightly faster (not 100% clear)
            // t = tic();
            // binLagrange(gg, i, j, k, cc);
            // t = tic() - t; 
            // std::cout << "BL2 " << t << "\n"; 
            current_dim++;
        }
        else if (current_dim == 2) {
#ifndef NDEBUG 
            FastInteger old = gram[k][k];
#endif
            FastCVP3(gram, i, j, k, Coords);
            
            assert(gram[k][k] <= old);

            if (gram[k][k] >= gram[j][j]) {
                current_dim++;
            }
            else {
                // maybe this is a bit unefficient as we know the first and second vectors are already sorted
                FastSortIndex(gram, i, j, k);
                // FastSortGram(gram, Coords);
                current_dim = 1;
            }

        } 

    }

    FastSortGram(gram, Coords);
} 

// we apply BinaryLagrange Reduction on all pairs 
// and then finish off with Semaev
void FastLagrange3(FastMat3 &gram, FastMat3 &Coords) {
    size_t i,j,k;
    i = 0;
    j = 1; 
    k = 2;

    FastSortIndex(gram, i, j, k);

    // std::cout << "input =" << gram << "\n";

    bool not_reduced = true;
    FastInteger old;
    int count = 0;
    while (not_reduced && count < 1000) {
        // std::cout << "loop !" << "\n";
        count++;
        not_reduced = false;
        
        // (i,j)
        old = gram[i][i];
        FastLagrange2(gram, i, j, k, Coords);
        not_reduced = gram[i][i] < old;

        // (i,k)
        old = gram[i][i];
        FastLagrange2(gram, i, k, j, Coords);
        not_reduced = not_reduced || (gram[i][i] < old);

        // (j,k)
        old = gram[j][j];
        FastLagrange2(gram, j, k, i, Coords);
        not_reduced = not_reduced || (gram[j][j] < old);

        if (!not_reduced) {
            FastCVP3(gram, i, j, k, Coords);
            not_reduced = gram[k][k] < gram[j][j];
            if (not_reduced) {
                FastSortIndex(gram, i, j, k);
            }
            
        }

    }

    // std::cout << "input = " << gram << "\n";

    // FastGreedyReduction3(gram, Coords);
    // std::cout << "output = " << gram << "\n\n\n";
    // FastCVP3(gram, i, j, k, Coords);
    FastSortGram(gram, Coords);

    

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
, quatlat order, std::pair<quat,quat> *small
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

    // std::array<std::array<FastInteger,3>,3> Red_Mat{{{0,0,0},{0,0,0},{0,0,0}}};
    NTL::mat_ZZ Red_Mat;
    Red_Mat.SetDims(3,3);
    for (int i = 0; i<3; i++) {
        for (int j=0; j<3; j++) {
            // NTL::ZZ c;
            gmp2ntl(Red_Mat[i][j], Red[i][j].get_data());
            // Red_Mat[i][j] = mpz_get_si(Red[i][j].get_data());
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
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            std::swap(Red_Mat[1][0], Red_Mat[2][0]);
            std::swap(Red_Mat[1][1], Red_Mat[2][1]);
            std::swap(Red_Mat[1][2], Red_Mat[2][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
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
                        Red_Mat[1] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
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
                        Red_Mat[2] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
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
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            std::swap(Red_Mat[2][0], Red_Mat[1][0]);
            std::swap(Red_Mat[2][1], Red_Mat[1][1]);
            std::swap(Red_Mat[2][2], Red_Mat[1][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
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
                        Red_Mat[1] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
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
                        Red_Mat[2] = i1 * Red_Mat[0] + i2 * Red_Mat[1] + i3 * Red_Mat[2];
                    }
                }
            }
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }
        if (IntList[1] > IntList[2]) {
            std::swap(IntList[1], IntList[2]);
            std::swap(gram_test[2][2], gram_test[1][1]);
            std::swap(gram_test[0][1], gram_test[0][2]);
            std::swap(gram_test[1][0], gram_test[2][0]);
            std::swap(Red_Mat[2][0], Red_Mat[1][0]);
            std::swap(Red_Mat[2][1], Red_Mat[1][1]);
            std::swap(Red_Mat[2][2], Red_Mat[1][2]);
        }
        if (IntList[0] > IntList[1]) {
            std::swap(IntList[0], IntList[1]);
            std::swap(gram_test[0][0], gram_test[1][1]);
            std::swap(gram_test[2][1], gram_test[2][0]);
            std::swap(gram_test[0][2], gram_test[1][2]);
            std::swap(Red_Mat[0][0], Red_Mat[1][0]);
            std::swap(Red_Mat[0][1], Red_Mat[1][1]);
            std::swap(Red_Mat[0][2], Red_Mat[1][2]);
        }

    // otherwise what follows doesn't work

    NTL::Vec<NTL::ZZ> row;
    row.SetLength(4);
    for (unsigned j = 0; j < 3; ++j) {
        for (unsigned i = 1; i <= 3; i++) {
            row[i] += Red_Mat[0][j]*order.basis[j+1][i]; 
        }
    }

    (small->first)[0] = row[0];
    (small->first)[1] = row[1];
    (small->first)[2] = row[2];
    (small->first)[3] = row[3];
    (small->first)[4] = order.denom;

    row[0] = Integer(0);
    row[1] = Integer(0);
    row[2] = Integer(0);
    row[3] = Integer(0);
    for (unsigned j = 0; j < 3; ++j) {
        for (unsigned i = 1; i <= 3; i++) {
            row[i] += Red_Mat[1][j]*order.basis[j+1][i]; 
        }
    }

    if (gram_test[0][1] <0) {
        row = -row;
    }

    (small->second)[0] = row[0];
    (small->second)[1] = row[1];
    (small->second)[2] = row[2];
    (small->second)[3] = row[3];
    (small->second)[4] = order.denom;

    small->first[4] >>= 1;
    small->second[4] >>= 1;

    // assert(order.contains(small->first));
    // assert(order.contains(small->second));
    if (!(small->first.norm().first/small->first.norm().second == gram_test[0][0]))
    {
        std::cout << small->first.norm().first/small->first.norm().second << " " << gram_test[0][0] << " " << gram_test[1][1] << "\n";
    }
    assert(small->first.norm().first/small->first.norm().second == IntList[0] );
    assert(small->second.norm().first/small->second.norm().second == IntList[1] );
    // small->first.normalize();
    // small->second.normalize();
    

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

Key order_invariant_computation_from_gram(quatlat order, NTL::mat_ZZ Gram, std::pair<quat,quat> *small) {
        Key K;
        NTL::ZZ IntList[3];
        // (void) order;
        // (void) small;

        // std::cout << clock() << "\n";
        // clock_t t = clock();
        // std::cout << "t =" << t << "\n";
        // auto t = tic();
        inner_order_invariant_computation_from_gram(IntList, Gram
        , order, small
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



Key fast_order_invariant_computation_from_gram(const FastQuatLat &FastO, FastMat3 &FastGram, std::pair<FastQuat,FastQuat> &fast_gamma_pair) {


        
        FastMat3 Coords;
        Coords[0][0] = 1;
        Coords[1][1] = 1;
        Coords[2][2] = 1;
        Coords[0][1] = 0;
        Coords[0][2] = 0;
        Coords[1][0] = 0;
        Coords[1][2] = 0;
        Coords[2][1] = 0;
        Coords[2][0] = 0;

        FastLagrange3(FastGram, Coords);

        fast_gamma_pair.first[4] = FastO.denom;
        fast_gamma_pair.second[4] = FastO.denom;

        for (unsigned i = 1; i <= 3; i++) {
            for (unsigned j = 0; j < 3; ++j) {
                fast_gamma_pair.first[i]  += Coords[0][j] * FastO.basis[j+1][i];
                fast_gamma_pair.second[i] += Coords[1][j] * FastO.basis[j+1][i]; 
    
            }
            if (FastGram[0][1] < 0) {
                fast_gamma_pair.second[i] = - fast_gamma_pair.second[i]; 
            }
        }
        fast_gamma_pair.first[4] >>= 1;
        fast_gamma_pair.second[4] >>= 1;

        Key K;

        FastIntToBytes(FastGram[0][0],K.IntList[0],LenNumberBytes);
        FastIntToBytes(FastGram[1][1],K.IntList[1],LenNumberBytes);
        FastIntToBytes(FastGram[2][2],K.IntList[2],LenNumberBytes);

        return K;
}