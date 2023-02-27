#include <gtest/gtest.h>
#include "../Src/HouseHolder/househ.h"

TEST(QR, B) {

    std::vector<double> a = { 3, 2, 4, 4};
    int n = sqrt(a.size());
    Matr Q(n), R(n);
    Matr A(a);
    A.out();
    Q = A.QRH().first;
    Matr B(a);
    R = B.QRH().second;
    Q.out();
    R.out();
    //Householder method

    Matr C(a);
    Q = C.QRG().first;
    Q.out();
    Matr D(a);
    R = D.QRG().second;
    R.out();

}
TEST(QR, T) {
    std::vector<double> a = { 3, 1,4,4, 5, 3, 1, 2, 3};
    Matr A(a);
    A.out();
    A.QRH().first.out();
    Matr B(a);
    B.QRH().second.out();
//householder working

    Matr C(a);
    C.out();
    C.QRG().first.out();
    Matr D(a);
    D.QRG().second.out();

    //All is working;)
}