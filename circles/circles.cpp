#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "circles.h"

using vector = std::vector<double>;
using matrix = std::vector<vector>;

vector random_position(double scale, double x0, double y0){

    // losuję punkt wewnątrz koła o R = scale o środku w punkcie
    // (x0, y0) z rozkładem jednorodnym

    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    double U1 = dist(generator);
    double U2 = dist(generator);

    double x = std::sqrt(-2.0*std::log(1-U1))*std::cos(2.0*PI*U2);
    double y = std::sqrt(-2.0*std::log(1-U1))*std::sin(2.0*PI*U2);
    double norm = std::sqrt(x*x + y*y);
    x /= norm; y /= norm;

    double U3 = std::sqrt(dist(generator));
    x = U3*x*scale + x0;
    y = U3*y*scale + y0;

    //double norm_2 = 1/PI/scale/scale;

    return {x, y};
}

int circle(double x_0, double y_0, double R, double x, double y){
    if ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) <= R*R) return 1;
    else return 0;
}

matrix common_part(double x0_1, double y0_1, double R_1,
                double x0_2, double y0_2, double R_2,
                std::initializer_list<int> N){

        // CAŁKOWANIE PRZEBIEGA PO KOLE 1
        // LICZONA JEST CZ. WSPÓLNA Z KOŁEM 2
        // dla n w N - liczona jest w. średnia i odch. std.

        int N_max = *std::max_element(N.begin(), N.end());

        double mu = 0;
        double mu2, mu_local;
        vector pos;
        int inside_1, inside_2;
        vector std_dev, n, mu_vec;
        
        for(int i = 1; i<N_max+1; i++){

            pos = random_position(R_1, x0_1, y0_1);
            inside_1 = circle(x0_1, y0_1, R_1, pos[0], pos[1]);
            inside_2 = circle(x0_2, y0_2, R_2, pos[0], pos[1]);

            if(inside_1 && inside_2){
                mu+=1.0;
            }

            if (std::find(N.begin(), N.end(), i) != N.end()){
                n.push_back(i);
                mu_local = R_1*R_1*PI*mu/(double)i;
                mu_vec.push_back(mu_local);
                mu2 = PI*R_1*R_1*mu_local;
                std_dev.push_back(std::sqrt((mu2 - mu_local*mu_local)/(double)i));
            }
        }
        
    return {n, mu_vec, std_dev};
    }

void plot_circles(std::initializer_list<matrix> matrices){

    std::ofstream plik("temp.dat");
    if(!plik.is_open()){ return; }
    vector pos;
    for(matrix m : matrices){
        for(auto element : m){
            plik << element[0] << " ";
        }
        plik << "\n";
        for(auto element : m){
            plik << element[1] << " ";
        }
        plik << "\n";
    }
    plik.close();
    std::string command = "python3 lab4_circles.py";
    int wynik = system(command.c_str());
    remove("temp.dat");
}

void plot_plots(std::initializer_list<matrix> n_mu_stddev){
    
    std::ofstream plik("temp.dat");
    if(!plik.is_open()){ return; }
    for(matrix m : n_mu_stddev){
        for(size_t i = 0; i < m[0].size(); i++){
            plik << m[0][i] << " ";
        }
        plik << "\n";
        for(size_t i = 0; i < m[1].size(); i++){
            plik << m[1][i] << " ";
        }
        plik << "\n";
        for(size_t i = 0; i < m[2].size(); i++){
            plik << m[2][i] << " ";
        }
        plik << "\n";
    }
    plik.close();
    std::string command = "python3 lab4_wykresy.py";
    int wynik = system(command.c_str());
    remove("temp.dat");
}

int main(){

    // parametry kół A i B
    double R_A = 1.0;
    double R_B = std::sqrt(2.0)*R_A;
    double x0_A = 1.0; double y0_A = 0.0;
    double x0_B = 0.0; double y0_B = 0.0;
    
    // ZADANIE 1
    // liczba punktów losowanych w każdym z kół
    int N = pow(10.0, 4.0);

    double x0_A_changed = R_A + R_B; // podmianka
    matrix K_A, K_B;
    for (int i = 0; i<N; i++){
        K_A.push_back(random_position(R_A, x0_A_changed, y0_A));
        K_B.push_back(random_position(R_B, x0_B, y0_B));
    }
    plot_circles({K_A, K_B});
    
    // ZADANIE 2
    N = pow(10.0, 6.0);

    // a)
    x0_A = R_B + 0.5*R_A;
    matrix nmd_a = common_part(x0_A, y0_A, R_A, x0_B, y0_B, R_B,
                                {100, 1000, 10000, 100000, 1000000});
    
    // b)
    x0_A = 0.0;
    matrix nmd_b = common_part(x0_A, y0_A, R_A, x0_B, y0_B, R_B,
                                {100, 1000, 10000, 100000, 1000000});

    // c)
    x0_A = R_B + 0.5*R_A;
    matrix nmd_c = common_part(x0_B, y0_B, R_B, x0_A, y0_A, R_A,
                                {100, 1000, 10000, 100000, 1000000});

    // d)
    x0_A = 0.0;
    matrix nmd_d = common_part(x0_B, y0_B, R_B, x0_A, y0_A, R_A,
                                {100, 1000, 10000, 100000, 1000000});

    
    plot_plots({nmd_a, nmd_b, nmd_c, nmd_d});
}