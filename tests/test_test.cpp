#include <gtest/gtest.h>
#include <iostream>
#include "../solvers/solve_tridiagonal.hpp"
TEST(A, B){
    std::vector<line> Mat = {{0, 2, 1},{1,10,-5},{1,-5,2},{1,4,0}};
    std::vector<double> free1 = {-5, -18, -40, -27};
    TriMat matr1(Mat);
    std::vector<double> solut = findsolution(Mat, free1);
    for (auto i: solut) {
        std::cout << i << " ";


    }
}