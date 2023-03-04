#ifndef SLAE_CSRMATRIX_H_GLEB
#define SLAE_CSRMATRIX_H_GLEB
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
struct DOK {
    int i;
    int j;
    double val;
};

namespace CSR_space {
    class Csr_matrix {
    private:


    public:
        std::vector<double> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> sup_row;

        double operator()(std::size_t i, std::size_t j);
        void out(Csr_matrix &A);

        std::vector<double> operator*(std::vector<double> fre);
        std::vector<double> ReshIT(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e);
        std::vector<double> ReshGZ(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        double MulGZU(std::vector<double> fre, int i);
        double MulGZD(std::vector<double> fre, int i);
        std::vector<double> ReshYak(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);
        std::vector<double> MulYa(std::vector<double> fre);
        Csr_matrix(std::vector<DOK> vec);

    };

    void out(Csr_matrix &A);

    std::vector<double> ReshIT(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e);

    std::vector<double> ReshGZ(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e);

    double MulGZU(std::vector<double> fre, int i);
};




#endif