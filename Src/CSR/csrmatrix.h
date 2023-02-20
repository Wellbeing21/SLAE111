#ifndef SLAE_CSRMATRIX_H_GLEB
#define SLAE_CSRMATRIX_H_GLEB

#include <vector>
struct DOK {
    int i;
    int j;
    double val;
};

namespace CSR_space {
    class Csr_matrix {
    private:
        std::vector<double> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> sup_row;

    public:


        double operator()(std::size_t i, std::size_t j);

        std::vector<double> operator*(std::vector<double> fre);

        Csr_matrix(std::vector<DOK> vec);

    };
};




#endif