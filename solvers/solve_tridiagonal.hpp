//
// Created by gorilla on 10.02.23.
//
#include "../Tridiagonal_matrix/tridiagonal.h"
std::vector <double> findsolution(const TriMat& Matr, std::vector<double> free) {
    int n = free.size();
    std::vector<double> p_(n);
    std::vector<double> q_(n);
    std::vector<double> sol(n);
    p_[0] = - Matr[0][2] / Matr[0][1];
    q_[0] = - free[0] / Matr[0][1];
    for (int i = 1; i < n; i++) {
        p_[i] = - Matr[i - 1][2] / (Matr[i - 1][0] * p_[i-1] + Matr[i-1][1]);
        q_[i] = (free[i-1] - Matr[i-1][0] * q_[i-1])/(Matr[i-1][0]* p_[i-1] + Matr[i-1][1]);
}
    sol[n-1] = (free[n-1] - Matr[n-1][0] * q_[n-1])/(Matr[n-1][0] * p_[n-1] + Matr[n-1][1]);
    for (int i = n-2; i >= 0; i--) {
    sol[i] = p_[i+1] * sol[i+1] + q_[i+1];

    }
    return sol;
}



