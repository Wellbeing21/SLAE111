#include <gtest/gtest.h>
#include <iostream>
#include "../Src/solvers/solve_tridiagonal.hpp"
#include "../Src/CSR/csrmatrix.h"
TEST(A, B){
    std::vector<line> Mat = {{0, 2, 1},{1,10,-5},{1,-5,2},{1,4,0}};
    std::vector<double> free1 = {-5, -18, -40, -27};
    TriMat matr1;
    matr1.apply(Mat);
    std::vector<double> solut = solver::findsolution(matr1, free1);
    for (auto i: solut) {
        std::cout << i << " ";
    }
}
TEST (B, C){
    std::vector<double> b;
    std::vector<DOK> a;
    a.resize(4);
    //there's massive: 1    5
    //                 0    2
    //                 0    0
    //                 0    1
    a[0] = {0, 0, 1};
    a[1] = {0, 1, 5};
    a[2] ={1, 1, 2};
    a[3] = {3, 1, 1};
    CSR_space::Csr_matrix matrix(a);
    std::cout << matrix(3,1) << "\n";
    //1-stt element
    std::cout << matrix(1,1) << "\n";
    //second
    std::cout << matrix(2,0) << "\n";
    //third
    b.resize(2);
    b = {1, 5};
    b = matrix*b;
    for (auto el : b) {
        std::cout << el << " ";
    }

}