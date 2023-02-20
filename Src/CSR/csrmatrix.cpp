//
// Created by gorilla on 18.02.23.
//

#include "csrmatrix.h"

        double CSR_space::Csr_matrix::operator()(std::size_t i, std::size_t j) {
            int c = 0;
            for (std::size_t it = sup_row[i]; it < sup_row[i+1]; it++) {
                if (col_ind[it] == j) {
                    return data[it];
                    c = 1;
                }
                else if (c == 0) {
                    return 0;
                }
            }
        }

        std::vector<double> CSR_space::Csr_matrix::operator*(std::vector<double> fre) {
            std::vector<double> a;
            std::size_t n = sup_row.size() - 1;
            a.resize(n);
            for (int it = 0; it < n; it++) {
                for (int jt = sup_row[it]; jt < sup_row[it+1]; jt++) {
                    a[it] += data[jt] * fre[col_ind[jt]];
                }
            }
            return a;
        }

        CSR_space::Csr_matrix::Csr_matrix(std::vector<DOK> vec) {
            data.resize(vec.size());
            sup_row.resize(vec[vec.size() - 1].i + 2);
            col_ind.resize(vec.size());
            sup_row[0] = 0;

            std::size_t j = 1;
            for (std::size_t i = 0; i < vec.size(); i++) {
                data[i] = vec[i].val;
                col_ind[i] = vec[i].j;

                if (i > 0 && vec[i].i - vec[i - 1].i != 0) {
                    for (int it = 0; it < vec[i].i - vec[i - 1].i; it++) {
                        sup_row[j] = i;
                        j++;
                    }
                }
            }

            sup_row[j] = vec.size();
        };