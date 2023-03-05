#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
struct DOK {
    int i;
    int j;
    double val;
};

    class Csr_matrix {
    private:
        std::vector<double> data;
        std::vector<std::size_t> col_ind;
        std::vector<std::size_t> sup_row;

    public:
        double operator()(std::size_t i, std::size_t j) {
            int c = 0;
            for (std::size_t it = sup_row[i]; it < sup_row[i+1]; it++) {
                if (col_ind[it] == j) {
                    return data[it];
                    c = 1;
                }
            }
            if (c == 0) {
                return 0;
            }
        }
        void out(Csr_matrix &A) {
            int n = (A.sup_row.size()-1);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    std::cout  << std::setw(6) << A(i,j) << " ";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }

        double ReshITc(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e, size_t c) {
            int n = A.sup_row.size() - 1;
            c = 0;
            std::vector<double> mul, bf;
            double r = 1;
            while (not (fabs(r) < e)) {
                c+=1;
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
            return c;

        }


        std::vector<double> ReshIT(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double t, double e) {
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

        std::vector<double> ReshGZ(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e) {
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
        double MulGZU(std::vector<double> fre, int i) {
            double a;
            std::size_t n = sup_row.size() - 1;
                for (int jt = sup_row[i]; jt < sup_row[i+1]; jt++) {
                    if (col_ind[jt] <= i) { continue;}
                    a += data[jt] * fre[col_ind[jt]];
                }
            return a;
        }
        double MulGZD(std::vector<double> fre, int i) {
            double a = 0;
            std::size_t n = sup_row.size() - 1;
            for (int jt = sup_row[i]; jt < sup_row[i+1]; jt++) {
                if (col_ind[jt] >= i) { continue;}
                a += data[jt] * fre[col_ind[jt]];
            }
            return a;
        }



        std::vector<double> ReshYak(Csr_matrix &A, std::vector<double> x0, std::vector<double> b, double e) {
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
        std::vector<double> MulYa(std::vector<double> fre) {
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

        std::vector<double> operator*(std::vector<double> fre) {
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

        Csr_matrix(std::vector<DOK> vec) {
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

    };


int main() {
//    double e = 0.000001;
//    //test for yakobi, its working
//    std::vector<DOK> A {{0,0,0.7},{0,1,0.2}, {1,1,0.2}, {2,0,0.1}, {2,2,0.3} };
//    Csr_matrix M(A);
//    M.out(M);
//    std::vector<double> x0{1, 1, 1}, b{1,1,1};
//    //taking vector, now we need to use yakobi method
//    x0 = M.ReshYak(M,x0,b,e);
//    for (int i = 0; i < 3; i++) {
//        std::cout << x0[i] << "\n";
//    }
//
////it's test for my G-Z method
//    std::vector<DOK> A1 {{0,0,4},{0,1,1},{1,0,1},{1,1,2}, {2,2,1} };
//    Csr_matrix M1(A1);
//    M1.out(M1);
//    std::vector<double> x01{1, 1, 3}, b1{1,1,1};
//    x01 = M1.ReshGZ(M1,x01,b1,e);
//    for (int i = 0; i < 3; i++) {
//        std::cout << x01[i] << "\n";
//    }
//
////it's test for my iterations method
//    double t = 0.1;
//    std::vector<DOK> A2 {{0,0,3},{0,1,0.5},{1,0,1},{1,1,2}, {2,2,1} };
//    Csr_matrix M2(A2);
//    M2.out(M2);
//    std::vector<double> x02{1, 1, 1}, b2{1,1,1};
//    x02= M2.ReshIT(M2,x02,b2,t,e);
//    for (int i = 0; i < 3; i++) {
//        std::cout << x02[i] << "\n";
//    }


    //These were tests for methods
    //3-rd
    double t3 = 1e-3, e3 = 1e-12, c=0;
    std::vector<DOK> a {
            {0, 0, 10},
            {0, 1, 1},
            {1, 0, 1},
            {1, 1, 7},
            {2, 1, 0.1},
            {2, 2, 1}
    };
    Csr_matrix M3(a);
    M3.out(M3);
    std::vector<double> b3{20, 30, 1}, x3{0,0,0};
    for (int i = 0; i < 100; i++) {
        t3 += 1e-4;
        c = M3.ReshITc(M3,x3,b3,t3,e3,c);
        std::cout << t3 << " " << c  <<"\n";
    }


    return 0;
}
