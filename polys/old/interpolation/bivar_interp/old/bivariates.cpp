#include <NTL/ZZ_pX.h>
#include <iostream>
#include <vector>
#include "bivariates.hpp"

int degree(const bivar_poly F){
    return F.cs.length() - 1;
}


NTL::ZZ_pX reciprocal_kronecker_sub(const bivar_poly F, const int Ly, const int Lx, const int N){
    // This can be change to make multiplication more efficient https://arxiv.org/abs/0712.4046
    NTL::ZZ_pX f;
    f.SetLength((Ly-1)*N+Lx-1);

    int degF = Ly - 1;

    for(int i = degF; i > -1; i--){
        NTL::ZZ_pX Fi = F.cs[degF-i];
        for(int j = deg(Fi); j > -1; j--){
            int ind = i*N+j;
            NTL::ZZ_p new_c = coeff(f,ind) + coeff(Fi, j);
            NTL::SetCoeff(f, ind, new_c);
        };
    }

    return f;

}

void bivariate_mult(bivar_poly &FG, const bivar_poly F, const bivar_poly G, const int Ly, const int Lx){
    
    int N = 2*Lx-1;

    NTL::ZZ_pX f = reciprocal_kronecker_sub(F, Ly, Lx, N);
    NTL::ZZ_pX g = reciprocal_kronecker_sub(G, Ly, Lx, N);
    NTL::ZZ_pX fg; 
    NTL::mul(fg, f, g);

    int degFG = degree(F)*degree(G);
    FG.cs.SetLength(degFG+1);

    for(int i = 0; i < degFG+1; i++){
        FG.cs[i].SetLength(N);
        for(int j = N-1; j > -1; j--){
            FG.cs[i][j] = coeff(fg, j+i*N);
        }
        FG.cs[i].normalize();
    }

    
}



// int main()
// {
//     NTL::ZZ_p::init(NTL::ZZ(1019));
//     NTL::ZZ_pX f1,f2,f3,f4,f5;
//     NTL::SetCoeff(f1, 2, NTL::ZZ_p(1));
//     NTL::SetCoeff(f1, 1, NTL::ZZ_p(2));
//     NTL::SetCoeff(f1, 0, NTL::ZZ_p(3));
//     NTL::SetCoeff(f2, 3, NTL::ZZ_p(5));
//     NTL::SetCoeff(f2, 0, NTL::ZZ_p(6));
//     NTL::SetCoeff(f3, 2, NTL::ZZ_p(1));
//     NTL::SetCoeff(f3, 1, NTL::ZZ_p(2));
//     NTL::SetCoeff(f3, 0, NTL::ZZ_p(3));
//     NTL::SetCoeff(f4, 2, NTL::ZZ_p(1));
//     NTL::SetCoeff(f4, 1, NTL::ZZ_p(2));
//     NTL::SetCoeff(f4, 0, NTL::ZZ_p(3));
//     NTL::SetCoeff(f5, 2, NTL::ZZ_p(1));
//     NTL::SetCoeff(f5, 1, NTL::ZZ_p(2));
//     NTL::SetCoeff(f5, 0, NTL::ZZ_p(3));
    
//     int Ly = 3;
//     int Lx = 4;
    

//     bivar_poly F;
//     F.cs.SetLength(Ly);
//     F.cs[0] = f1;
//     F.cs[1] = f2;
//     F.cs[2] = f3;

//     bivar_poly G;
//     G.cs.SetLength(Ly);
//     G.cs[1] = f5;
//     G.cs[0] = f4;
    

//     bivar_poly FG;
//     bivariate_mult(FG, F, G, Ly, Lx);

//     std::cout << "FG: " << FG.cs << std::endl;

//     return 0;
// }

