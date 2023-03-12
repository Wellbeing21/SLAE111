#include <gtest/gtest.h>
#include "../Src/solvers/solve_tridiagonal.hpp"
#include "../Src/CSR/csrmatrix.h"
using namespace CSR_space;
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
TEST (D, E){
    double e = 0.000001;
    //test for yakobi, its working
    std::vector<DOK> A {{0,0,0.7},{0,1,0.2}, {1,1,0.2}, {2,0,0.1}, {2,2,0.3} };
    Csr_matrix M(A);
    M.out(M);
    std::vector<double> x0{1, 1, 1}, b{1,1,1};
    //taking vector, now we need to use yakobi method
    x0 = M.ReshYak(M,x0,b,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x0[i] << "\n";
    }

}
TEST (F, G){
    double e = 0.000001;
    //it's test for my G-Z method
    std::vector<DOK> A1 {{0,0,4},{0,1,1},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix M1(A1);
    M1.out(M1);
    std::vector<double> x01{1, 1, 3}, b1{1,1,1};
    x01 = M1.ReshGZ(M1,x01,b1,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x01[i] << "\n";
    }

}
TEST (H, I){
    double e = 0.000001;
    double t = 0.1;
    std::vector<DOK> A2 {{0,0,3},{0,1,0.5},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix M2(A2);
    M2.out(M2);
    std::vector<double> x02{1, 1, 1}, b2{1,1,1};
    x02= M2.ReshIT(M2,x02,b2,t,e);
    for (int i = 0; i < 3; i++) {
        std::cout << x02[i] << "\n";
    }
    //it's working
}
TEST (IT, Smart){
//new work
    double e = 0.000001;
    std::vector<DOK> AA {{0,0,3},{0,1,0.5},{1,0,1},{1,1,2}, {2,2,1} };
    Csr_matrix M3(AA);
    M3.out(M3);
    std::vector<double> x03{1, 1, 1}, b3{1,1,1};

x03= M3.ReshITS(M3,x03,b3,3,e,2.,8.);
    for (int i = 0; i < 3; i++) {
        std::cout << x03[i] << "\n";
    }


}
