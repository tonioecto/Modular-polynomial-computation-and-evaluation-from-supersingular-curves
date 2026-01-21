#pragma once

#include <NTL/ZZ.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <array>
#include <string>
#include <sstream>
#include <functional>
#include "choosetorsion.hpp"
#include "hashmap.hpp"

typedef int SmallInteger;
typedef std::array<std::array<unsigned char, 2>, 2> VerySmallMat; 
typedef std::pair<unsigned char, unsigned char> VerySmallIntegerPair;

typedef std::pair<SmallInteger,SmallInteger> SmallIntegerPair;

// unsigned char MulMod16(const unsigned char &a, const unsigned char &b, ) {

// } 

unsigned char InvModSpecial(const unsigned char &a);

class SmallMatFp {
    public:
        SmallInteger mod;
        std::array<std::array<SmallInteger,2>,2> mat;

        SmallMatFp(const SmallInteger &m) : mod{m}, mat{{{0,0},{0,0}}} {};
        SmallMatFp(const SmallInteger &m, const std::array<std::array<SmallInteger,2>,2> _mat) : mod{m}, mat{_mat} {};

        SmallMatFp operator+(const SmallMatFp &other) const;
        SmallMatFp operator*(const SmallInteger &other) const;
        SmallIntegerPair operator*(const SmallIntegerPair &other) const;
        friend std::ostream& operator<<(std::ostream& o, SmallMatFp const &M)
        {
                return o
                << "["
                << M.mat[0][0] << ", " << M.mat[0][1] << "\n " << M.mat[1][0] << ", " << M.mat[1][1] <<  "] mod " << M.mod << "\n";
        }
        void normalize();
};



// Find the matrix M to such that bas2 = M * bas1 for points of order 48
SmallMatFp change_of_basis16(const std::pair<ecp, ecp> &bas1, const std::pair<ecp, ecp> &bas2);
SmallMatFp change_of_basis3(const std::pair<ecp, ecp> &bas1, const std::pair<ecp, ecp> &bas2);
