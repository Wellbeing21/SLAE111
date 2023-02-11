//
// Created by gorilla on 10.02.23.
//

#ifndef SLAE_TRIDIAGONAL_H
#define SLAE_TRIDIAGONAL_H
#include <vector>
#include <array>
struct line {
    double f;
    double s;
    double t;
};
class TriMat {
    //private field
    std::vector <line> data;
public: //pity
    TriMat(const std::vector<line>& data1);
    //overdosing operator
    std::array<double, 3> operator[] (int i) const;


};




#endif //SLAE_TRIDIAGONAL_H
