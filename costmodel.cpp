///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
/////// The code in this file builds the cost model with which we choose the 'optimal' 
///////  torsion we work with in choosetorsion.cpp
///////
/////// The code is heavily inspired by:
///////         https://github.com/friends-of-quaternions/deuring/blob/main/cost.py
///////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZ.h>

#include "costmodel.hpp"


NTL::RR karatsuba(int const &k, NTL::RR const &logp, NTL::RR const &loglogp)
{    
    NTL::RR log3 = NTL::log(NTL::RR(3))/NTL::log(NTL::RR(2));
    NTL::RR kpowlog3;
    NTL::conv(kpowlog3, k);
    kpowlog3 = NTL::exp(log3 * NTL::log(kpowlog3));
    if (logp < 32) {
        NTL::RR temp = NTL::RR(1)/30 + logp/1000 + NTL::power(logp, 2)/500000;
        if (logp >= 30) {
            temp += NTL::RR(1)/50;
            if (logp >= 31) {
                temp += NTL::RR(1)/30;
            }
        } 
        return kpowlog3 * temp;
    } else if (logp <= 64) {
        return kpowlog3 * (NTL::RR(1)/30 + logp/190);
    } else {
        return kpowlog3 * (NTL::RR(1)/30 + loglogp/50);
    }
}


model_function cost_model(NTL::ZZ const &p)
{
    NTL::RR logp = NTL::log(p)/NTL::log(NTL::RR(2));
    NTL::RR loglogp = NTL::log(logp)/NTL::log(NTL::RR(2));
    NTL::RR log3 = NTL::log(NTL::RR(3))/NTL::log(NTL::RR(2));
    int cutoff;
    
    //TODO: update constants
    NTL::RR c1 = NTL::RR(0.01);
    //NTL::RR c2 = NTL::RR(2.65);
    NTL::RR c2 = NTL::RR(0.65); //Trying penalising higher extension fields less
    NTL::RR c3 = NTL::RR(1.25);
    NTL::RR c4 = NTL::RR(0.022);

    //TODO: measure actual values here
    //for now copied values from DFTP
    if (logp < 32){
        cutoff = 79;
    }
    else if (logp <= 64){
        cutoff = 44;
    }
    else {
        cutoff = 55;
    }

    auto quasilinear = [](int const &k, NTL::RR const &logp, NTL::RR const &loglogp, int const &cutoff) -> NTL::RR {
        auto fun0 = [logp](int x) -> NTL::RR {
            NTL::RR xRR;
            NTL::conv(xRR, x);
            return xRR * NTL::log(xRR)/NTL::log(NTL::RR(2)) * (NTL::RR(1)/10 + logp/200);
        };

        NTL::RR off = karatsuba(cutoff, logp, loglogp) - fun0(cutoff);
        NTL::RR temp = fun0(k);
        return off + temp;
    };

    auto field_model = [&quasilinear](int const &kk, NTL::RR const &logp, NTL::RR const &loglogp, int const &cutoff) -> NTL::RR {
        NTL::RR oneop;

        if (kk < cutoff)
        {
            oneop = karatsuba(kk, logp, loglogp);
        }
        else
        {
            oneop = quasilinear(kk, logp, loglogp, cutoff);
        }
        return oneop;
    };

    return [c1, c2, c3, c4, logp, loglogp, cutoff, &field_model](int const &ell,  int const &exp, int const &k) -> NTL::RR {
        NTL::RR kRR;
        NTL::conv(kRR, k);
        NTL::RR logell;
        NTL::conv(logell, ell);
        logell = NTL::log(logell)/NTL::log(NTL::RR(2));
        NTL::RR oneop = field_model(k, logp, loglogp, cutoff);
        return exp*c2*exp*oneop*logell + c3*exp*ell*(kRR+c4*logell*logell) + c1*oneop*ell*logell;
    };
}


NTL::RR adj_cost(NTL::RR const &cost, int const &ell)
{
    NTL::RR temp;

    temp = NTL::log(NTL::RR(ell))/NTL::log(NTL::RR(2));

    if (temp < 1){
        temp = NTL::RR(1);
    }

    return cost/temp;
}
