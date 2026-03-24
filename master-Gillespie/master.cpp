#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <string>
#include "master.h"

double U(){
    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(generator);
}

void evolution(int x1_0, int x2_0, int x3_0,
            double k1, double k2, double k3, double k4,
            int t_max, int P_max, std::string filename){
    
    
    int N = 50;
    double delta_time = static_cast<double>(t_max)/N;

    std::vector<double> h1(N, 0.0);
    std::vector<double> h2(N, 0.0);

    std::ofstream output_x1x2x3(filename + "_x1x2x3.txt");

    for (int p = 1; p <= P_max; p++){
        std::vector<double> h0(N, 0.0);
        std::vector<int> ncount(N, 0);

        double t = 0.0;
        double x1 = static_cast<double>(x1_0);
        double x2 = static_cast<double>(x2_0);
        double x3 = static_cast<double>(x3_0);

        while (t < t_max){
            double Gamma1 = k1;
            double Gamma2 = k2;
            double Gamma3 = k3*x1*x2;
            double Gamma4 = k4*x3;
            double Gamma_max = Gamma1 + Gamma2 + Gamma3 + Gamma4;

            double U1 = U();
            double dt = -std::log(U1)/Gamma_max;
            double U2 = U();

            if (U2 <= Gamma1/Gamma_max) { // zdarzenie 1
                x1 += 1.0;
            } else if (U2 <= (Gamma1 + Gamma2)/Gamma_max) { // zdarzenie 2
                x2 += 1.0;
            } else if (U2 <= (Gamma1 + Gamma2 + Gamma3)/Gamma_max) { // zdarzenie 3
                x1 -= 1.0;
                x2 -= 1.0;
                x3 += 1.0;
            } else { // zdarzenie 4
                x3 -= 1.0; 
            }

            t += dt;

            output_x1x2x3 << t << " " << x1 << " " << x2 << " " << x3 << "\n";

            int l = static_cast<int>(std::floor(t/delta_time));
            if (l < N) {
                h0[l] += x3;
                ncount[l]++;
            }
        }

        for (int i = 0; i < N; i++) {
            if (ncount[i] > 0) {
                double x3t = h0[i]/ncount[i];
                h1[i] += x3t;
                h2[i] += x3t*x3t;
            }
        }
    }

    output_x1x2x3.close();
    std::ofstream output(filename + ".txt");
    for(int l = 0; l < N; l++){
        double x1_3t = h1[l]/P_max;
        double x2_3t = h2[l]/P_max;
        double sigma_x3t = std::sqrt((x2_3t - x1_3t*x1_3t)/(P_max));
        double t = (l + 0.5)*delta_time;

        output << t << " " << x1_3t << " " << sigma_x3t << "\n";
    }
    output.close();
;
}

int main(){

double k1 = 1.0;
double k2 = 1.0;
double k3 = 0.001;
double k4 = 0.01;
int x1_0 = 120;
int x2_0 = 80;
int x3_0 = 1;
double tmax = 200.0;
int Pmax = 5;

//evolution(x1_0, x2_0, x3_0, k1, k2, k3, k4, tmax, Pmax, "histogram_wyniki");

Pmax = 100;
evolution(x1_0, x2_0, x3_0, k1, k2, k3, k4, tmax, Pmax, "histogram_wyniki_100");

Pmax = 1;
evolution(x1_0, x2_0, x3_0, k1, k2, k3, k4, tmax, Pmax, "histogram_wyniki_1");

}