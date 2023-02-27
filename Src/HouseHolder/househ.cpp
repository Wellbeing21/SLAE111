#include "househ.h"
    Matr::Matr(std::vector<std::vector<double>> &a) {
        //with vector of vectors
        int n = a.size();
        vv.resize(n * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vv[i * n + j] = a[i][j];
            }
        }
    }
    Matr::Matr(int n) {
        //with zeros
        vv.resize(n * n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                vv[i * n + j] = 0;
            }
        }
    }
    Matr::Matr (std::vector <double> &a) {
        int n = a.size();
        vv.resize(n);
        for (int i = 0; i < n; i++) {
            vv[i] = a[i];
        }
    }
    void Matr::out() {
        int n = sqrt( vv.size());
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                std::cout << std::setw(10) << vv[i * n + j] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }
    double Matr::operator()(int i, int j) const{
        int n = sqrt(vv.size());
        return vv[i * n + j];
    }
//with const if we don't need to change something
    double Matr::operator()(int i, int j) {
        int n = sqrt(vv.size());
        return vv[i * n + j];
    }
//without(to change smth)
    Matr Matr::operator* (Matr &M) {
        int n = sqrt(M.vv.size());
        Matr m(n);
//created empty matrix
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    m.vv[i * n + j] += vv[i * n + k] * M.vv[k * n + j];
                }
            }
        }
        return m;
    }
    std::pair<Matr, Matr> Matr::QRH () {
        int n = sqrt(vv.size());
        Matr Q(n);
        //try to make Q matrix
        double modm = 0;
        std::vector<double> x(n);
        for (int i = 0; i < n; i++){
            x[i] = vv[i * n];
            modm += vv[i * n] * vv[i * n];
        }
        modm = sqrt(modm);
        x[0] += modm;
        modm = 0;
        //made 1-st vector normal
        for (int i = 0; i < n; i++) {
            Q.vv[i * n + i] += 1;
            modm += x[i] * x[i];
        }
        //made I and |V|^2
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Q.vv[i * n + j] -= 2 / modm * x[i] * x[j];
            }
        }
        //YEAH i made P1

        for (int ht = 0; ht < n - 1; ht++) {
            //householder times

            std::vector<double> nor(n);
            double mod = 0;
            for (int i = ht; i < n; i++) {
                mod += vv[i * n + ht] * vv[i * n + ht];
            }
            mod = sqrt(mod);
            //take our |x|
            for (int i = ht; i < n; i++) {
                nor[i] = vv[i * n + ht];
            }
            nor[ht] += mod;
            //created normal vector
            mod = 0;
            for (int i = ht; i < n; i++) {
                mod += nor[i] * nor[i];
            }
            //created |V|^2
            ////was for (int vn = ht; vn < n; vn++)
            for (int vn = 0; vn < n; vn++) {
                //vector number
                double sc = 0;
                //scalar multiplication
                for (int i = ht; i < n; i++) {
                    sc += vv[i * n + vn] * nor[i];
                }
                //finded our scalar multiplication for every colomn vector
                for (int i = ht; i < n; i++) {
                    vv[i * n + vn] -= 2 * sc / mod * nor[i];
                    // was Q.vv[i  * n + vn] -= 2 * sc / mod * nor[i];
                }
                sc = 0;
                if (ht > 0) {
                    for (int i = ht; i < n; i++) {
                        sc += Q.vv[vn * n + i] * nor[i];
                    }

                    for (int i = ht; i < n; i++) {
                        Q.vv[vn  * n + i] -= 2 * sc / mod * nor[i];
                        //but i can't do it like this, because scalar multiplication isn't the same
                    }
                }

            }
        }
        //and we need to transpose
        Matr R(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ( i > j && fabs(R.vv[i * n + j]) < 0.00000000001) {
                    R.vv[i * n + j] = 0;
                }
                else {R.vv[i * n + j] =  vv[i * n + j];}
            }
        }
        return std::make_pair(Q, R);


    }
    std::pair<Matr, Matr> Matr::QRG () {
        int n = sqrt(vv.size());
        //here we need to check if [0,0] ne 0, so
        if (vv[0] == 0) {
            bool stop = false;
            int i = 0;
            while (vv[i * n] == 0) {
                if (stop) {
                    break;
                }
                i++;
                if (vv[i * n] != 0) {
                    stop = true;
                    //take givence rotation with this element
                    for (int j = 0; j < n; j++) {
                        double temp;
                        temp = vv[0 * n + j];
                        vv[0 * n + j] = vv[i * n + j];
                        vv[i * n + j] = -temp;

                    }
                }
            }

        }
        ////its working: if we have zeros in first colomn, they will rotate to value
        Matr Q(n);
        for (int i = 0; i < n; i++) {
            Q.vv[i * n + i] = 1;
        }
        //created P1
        double c = vv[0] / sqrt(vv[0] * vv[0] + vv[1 * n] * vv[1 * n]), s = vv[1 * n] / sqrt(vv[0] * vv[0] + vv[1 * n] * vv[1 * n]);
        Q.vv[0] = c;
        Q.vv[1 * n + 1] = c;
        Q.vv[0 * n + 1] = s;
        Q.vv[1 * n + 0] = -s;
        ////P1's working
        for (int gt = 0; gt < n; gt++) {
            for (int sn = gt + 1; sn < n; sn++) {

                //here we have Givence time and string number
                double sinn = vv[sn * n + gt] / sqrt(vv[sn * n + gt] * vv[sn * n + gt] + vv[gt * n + gt] * vv[gt * n + gt]);
                double coss = vv[gt * n + gt] / sqrt(vv[sn * n + gt] * vv[sn * n + gt] + vv[gt * n + gt] * vv[gt * n + gt]);
                //for every rotation we need to find our sinn and coss(for every string)
                for (int j = 0; j < n; j++) {
                    double base = 0, string = 0, baseq = 0, stringq = 0;
                    base += vv[gt * n + j] * coss + vv[sn * n + j] * sinn;
                    string += vv[sn * n + j] * coss - vv[gt * n + j] * sinn;
                    //temporarily to multiply
                    //i - nomer coordinatu
                    //starts with next
                    vv[gt * n + j] = base;
                    vv[sn * n + j] = string;
                    if ((sn > gt + 1) or (gt > 0)) {
                        baseq += Q.vv[gt * n + j] * coss + Q.vv[sn * n + j] * sinn;
                        stringq += Q.vv[sn * n + j] * coss - Q.vv[gt * n + j] * sinn;
                        Q.vv[gt * n + j] = baseq;
                        Q.vv[sn * n + j] = stringq;
                    }

                    //transform our sn string here
                }

                //we multiplied here our sn string on matrix, but also first string
////its working also


            }
        }
        Matr R(vv);
        return std::make_pair(Q, R);
    }