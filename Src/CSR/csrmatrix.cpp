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
void CSR_space::Csr_matrix::out(Csr_matrix &A) {
    int n = (A.sup_row.size()-1);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout  << std::setw(6) << A(i,j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}
std::vector<double> CSR_space::Csr_matrix::ReshIT(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e) {
    int n = A.sup_row.size() - 1;
    std::vector<double> mul, bf;
    double r = 1;
    while (not (fabs(r) < e)) {
        mul = A*x0;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - t * mul[i] + t * b[i];
        }
        bf = A*x0;
        r = fabs(fabs(bf[0]) - fabs(b[0]));
        for (int i = 0; i < n; i++) {
            if (fabs(fabs(bf[i]) - fabs(b[i])) > r) {
                r = (fabs(fabs(bf[i]) - fabs(b[i])));
            }
        }

    }
    return x0;

}
//same RESH with simple iterations, but with smart polinoms
//n = degree of our polynom
std::vector<double> CSR_space::Csr_matrix::ReshITS(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, int rd, double e, double lam1, double lam2) {
    int n = A.sup_row.size() - 1;
    int num = pow(2, rd);
    std::vector<double> mul, bf, t(num), ord(num), x(num);
    //first, we need to make our t-vector with cycle with num roots
    double sina = sin(4 * atan(1) / num), cosa = cos(4 * atan(1) / num), sinb = sin(4 * atan(1) / (2 * num));
    x[0] = cos(4 * atan(1) / (2 * num));
    t[0] = 1 / ((lam1 + lam2) / 2 + x[0] * (lam2 - lam1) / 2);
    for (int i = 1; i < num; i++) {
        x[i] = x[i - 1] * cosa - sinb * sina;
        t[i] = 1 / ((lam1 + lam2) / 2 + x[i] * (lam2 - lam1) / 2);
        sinb = sinb * cosa + x[i - 1] * sina;
    }
    //ok we found xk and tk
    //now we need order here to take tk with it's number
    double j = num / 2;
    ord[0] = 0;
    ord[j] = 1;
    for (int i = 1; i < n; i++) {
        j = j/2;
        for (int jt = j; jt < num; jt += 2*j) {
            ord[jt] = pow(2,i+1) - ord[jt - j] - 1;
        }
    }
    //and we found our order so, let's start iterations
    double r = 1;
    int por = 0;
    int iter = 0;
    while (not (fabs(r) < e)) {
        iter++;
        mul = A*x0;
        for (int i = 0; i < n; i++) {
            x0[i] = x0[i] - t[ord[por]] * mul[i] + t[ord[por]] * b[i];
        }
        por++;
        if (por >= num-1) {
            por = 0;
        }
        //we increase our t_por
        bf = A*x0;
        r = fabs(fabs(bf[0]) - fabs(b[0]));
        for (int i = 0; i < n; i++) {
            if (fabs(fabs(bf[i]) - fabs(b[i])) > r) {
                r = (fabs(fabs(bf[i]) - fabs(b[i])));
            }
        }

    }
    std::cout << iter <<"- this is the number of iterations" << "\n";
    return x0;

}
std::vector<double> CSR_space::Csr_matrix::ReshGZ(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e) {
    int n = A.sup_row.size() -1;
    double r = 1;
    std::vector<double> bf;
    while (not fabs(r) < e) {
        for (int i = 0; i < n; i++) {
            x0[i] = b[i] / A(i,i) - MulGZU(x0,i)/A(i,i) - MulGZD(x0,i)/A(i,i);
        }
        bf = A*x0;
        r = fabs(fabs(bf[0]) - fabs(b[0]));
        for (int i = 0; i < n; i++) {
            if (fabs(fabs(bf[i]) - fabs(b[i])) > r) {
                r = (fabs(fabs(bf[i]) - fabs(b[i])) > r);
            }
        }
        //finded supremum


    }
    return x0;
}
//multiplication wihh only UPTriangle matrix U
double CSR_space::Csr_matrix::MulGZU(std::vector<double> fre, int i) {
    double a;
    std::size_t n = sup_row.size() - 1;
    for (int jt = sup_row[i]; jt < sup_row[i+1]; jt++) {
        if (col_ind[jt] <= i) { continue;}
        a += data[jt] * fre[col_ind[jt]];
    }
    return a;
};
double CSR_space::Csr_matrix::MulGZD(std::vector<double> fre, int i) {
    double a = 0;
    std::size_t n = sup_row.size() - 1;
    for (int jt = sup_row[i]; jt < sup_row[i+1]; jt++) {
        if (col_ind[jt] >= i) { continue;}
        a += data[jt] * fre[col_ind[jt]];
    }
    return a;
}



std::vector<double> CSR_space::Csr_matrix::ReshYak(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e) {
    int n = (A.sup_row.size() - 1);
    std::vector<double> ans;
    double r = 1;
    //multiplied without diogonale
    while (not fabs(r) < e) {
        x0 =  A.MulYa(x0);
        for (int i = 0; i < n; i++) {
            x0[i] = b[i] / A(i, i) - x0[i] / A(i,i);
        }
        ans = A * x0;
        r = fabs(ans[0]) - fabs(b[0]);
        for (int i = 1; i < n; i++) {
            if (fabs(fabs(ans[i]) - fabs(b[i])) > r) {
                r = fabs(ans[i]) - fabs(b[i]);
            }
        }
    }
    //count our neuvyazky po supremumu
    return x0;

}
//multiplication without diogonale elements
std::vector<double> CSR_space::Csr_matrix::MulYa(std::vector<double> fre) {
    std::vector<double> a;
    std::size_t n = sup_row.size() - 1;
    a.resize(n);
    for (int it = 0; it < n; it++) {
        for (int jt = sup_row[it]; jt < sup_row[it+1]; jt++) {
            if (col_ind[jt] == it) { continue;}
            a[it] += data[jt] * fre[col_ind[jt]];
        }
    }
    return a;
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