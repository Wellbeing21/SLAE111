//
// Created by gorilla on 10.02.23.
//

#include "tridiagonal.h"

TriMat::TriMat(const std::vector<line>& data1) {
 data = data1;

}

std::array<double, 3> TriMat::operator[](int i) const{

    static std::array<double, 3> arr;
    arr[0] = data[i].f;
    arr[1] = data[i].s;
    arr[2] = data[i].t;
    //initialization
    return arr;

}
