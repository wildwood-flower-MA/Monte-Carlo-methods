#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include "elements.h"

namespace cst {
    inline constexpr double pi = 3.14159265358979323846;
    inline const std::vector<double> e_x = {1.0, 0.0};
    inline const std::vector<double> e_y = {0.0, 1.0};
}

// 1.1 rozkład sferycznie konturowany- normalny

inline double random_uniform_01(){
    /** Mersenne Twister **/

    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(generator);
}

std::vector<double> random_2D_BoxMuller(){

    double rand_1 = random_uniform_01();
    double rand_2 = random_uniform_01();
    std::vector<double> xy = {std::sqrt(-2.0*std::log(1-rand_1))*std::cos(2.0*cst::pi*rand_2),
                            std::sqrt(-2.0*std::log(1-rand_1))*std::sin(2.0*cst::pi*rand_2)};
    
    return xy;
}

// 1.2 Rozkład jednorodny w kole K^2(0,1)

inline std::vector<double> normalise(const std::vector<double>& pair){

    if(pair[0] == 0.0 && pair[1] == 0.0){
        return pair;
    }

    double norm = (double)0.0;
    for(auto vi = pair.begin(); vi != pair.end(); vi++){
        norm += (*vi)*(*vi); 
    }
    std::vector<double> xy = {pair[0]/std::sqrt(norm), pair[1]/std::sqrt(norm)};

    return xy;
}

std::vector<double> random_circle() {

    std::vector<double> norm_vec = normalise(random_2D_BoxMuller());
    double R = std::sqrt(random_uniform_01());

    return {R*norm_vec[0], R*norm_vec[1]};
}

// 1.3 Transformacja afiniczna: koło → elipsa
// 1.3.1 Wybór osi i skalowanie

std::vector<double> operator*(const std::vector<std::vector<double>>& A,
    const std::vector<double>& v) {
    return { A[0][0]*v[0] + A[0][1]*v[1], A[1][0]*v[0] + A[1][1]*v[1] };
}

std::vector<double> operator*(double c, const std::vector<double>& v) {
    // skalowanie wektora stałą c

    std::vector<double> xy = {c*v[0], c*v[1]};
return xy;
}

std::vector<double> operator+(const std::vector<double>& v, const std::vector<double>& u) {
    // dodawanie wektorów v+u
return {v[0] + u[0], v[1] + u[1]};
}

std::vector<std::vector<double>> create_R(double alpha){
    // tworzenie macierzy rotacji w konwencji
    return {{std::cos(alpha), std::sin(alpha)}, {-std::sin(alpha), std::cos(alpha)}};
}

// 1.4 Wyznaczanie macierzy kowariancji

std::pair<double,double> avg_x_y(const std::vector<std::vector<double>>& pairs){
    int size = pairs.size();
    double sum_x = 0.0;
    double sum_y = 0.0;
    for(int i = 0; i<size; i++){
        sum_x += pairs[i][0];
        sum_y += pairs[i][1];
    }
    return {sum_x/(double)size, sum_y/(double)size};
}

std::pair<double,double> avg_x2_y2(const std::vector<std::vector<double>>& pairs){
    int size = pairs.size();
    double sum_x2 = 0.0;
    double sum_y2 = 0.0;
    for(int i = 0; i<size; i++){
        sum_x2 += pow(pairs[i][0],2.0);
        sum_y2 += pow(pairs[i][1],2.0);
    }
    return {sum_x2/(double)size, sum_y2/(double)size};
}

double avg_xy(const std::vector<std::vector<double>>& pairs){
    int size = pairs.size();
    double sum_xy = 0.0;
    for(int i = 0; i<size; i++){
        sum_xy += pairs[i][0]*pairs[i][1];
    }
    return sum_xy/(double)size;
}

std::vector<std::vector<double>> cov_matrix(const std::vector<std::vector<double>>& pairs){

    std::pair<double,double> avg = avg_x_y(pairs);
    std::pair<double,double> avg2 = avg_x2_y2(pairs);
    double avgxy = avg_xy(pairs);
    double c00 = avg2.first - avg.first*avg.first;
    double c = avgxy - avg.first*avg.second;
    double c11 = avg2.second - avg.second*avg.second;

    return {{c00, c}, {c ,c11}};
}

double corr(const std::vector<std::vector<double>>& pairs){

    std::vector<std::vector<double>> cor_mat = cov_matrix(pairs);
    return cor_mat[1][0]/std::sqrt(cor_mat[0][0]*cor_mat[1][1]);
    }

// 1.5 Transformacja afiniczna a macierz kowariancji dla rozkładu gaussowskiego

void transform(const std::vector<std::vector<double>> A,
    std::vector<std::vector<double>>& pairs,
    double a1, double a2,
    double b1, double b2,
    double c1, double c2){
    /*
    konwencja jest taka:
    a1, a2 - skalowanie po transformacji A,
    b1, b2 - skalowanie przed transformacją A,
    c1, c2 - przesunięcie po transformacji A
    */

    for (std::vector<double>& pair : pairs) {

        pair = {b1*pair[0], b2*pair[1]};
        pair = A*pair;
        pair = {a1*pair[0], a2*pair[1]};
        pair = pair + (c1*cst::e_x) + (c2*cst::e_y);

    }
}

// OGÓLNE

void zapisz(std::ostream& my_out, const std::vector<std::vector<double>>& pairs){
    for(auto iter = pairs.begin(); iter != pairs.end(); iter++){
        my_out << (*iter)[0] << " " << (*iter)[1];
        if (std::next(iter) != pairs.end()) {
            my_out << "\n";
        }
    }
}

void rysuj(const std::vector<std::vector<double>>& pairs, std::string name){
    std::ofstream plik("lab3.dat");
    if(!plik.is_open()){
        return;
    }
    zapisz(plik, pairs);
    plik.close();
    std::string command = "python3 lab3.py "+name;
    int wynik = system(command.c_str());
    //std::cout<<wynik;
    remove("lab3.dat");
}

int main(){

    int N_normalny = pow(10,4);
    std::vector<std::vector<double>> Normalny_2D;
    Normalny_2D.reserve(N_normalny);
    for(int i=0; i<N_normalny; i++){
        Normalny_2D.emplace_back(random_2D_BoxMuller());
    }
    rysuj(Normalny_2D, "rozklad_normalny_N(0,1)");

    int N_sferycznie_konturowany = pow(10,4);
    std::vector<std::vector<double>> SK_2D;
    SK_2D.reserve(N_sferycznie_konturowany);
    for(int i=0; i<N_sferycznie_konturowany; i++){
        SK_2D.emplace_back(random_circle());
    }
    rysuj(SK_2D, "rozklad_sferycznie_konturowany");

    double alpha = cst::pi/(double)4.0;
    transform(create_R(-alpha), SK_2D, 1.0, 1.0, 1.0, 0.2, 0.0, 0.0);
    rysuj(SK_2D, "rozklad_jednorodny_w_kole");

    std::vector<std::vector<double>> cm = cov_matrix(SK_2D);
    std::cout<<"Macierz_kowariancji_(jednorodny):"<<"\n";
    std::cout<<cm[0][0]<<" "<<cm[0][1]<<"\n"<<cm[1][0]<<" "<<cm[1][1]<<"\n";

    double cr = corr(SK_2D);
    std::cout<<"wsp. korelacji: "<<cr<<"\n";

    std::cout<<"\n";
    transform(create_R(-alpha), Normalny_2D, 1.0, 1.0, 1.0, 0.2, 0.0, 0.0);
    rysuj(Normalny_2D, "rozklad_normalny");

    std::vector<std::vector<double>> cm_normalny = cov_matrix(Normalny_2D);
    std::cout<<"Macierz_kowariancji_(normalny):"<<"\n";
    std::cout<<cm_normalny[0][0]<<" "<<cm_normalny[0][1]<<"\n"<<cm_normalny[1][0]<<" "<<cm_normalny[1][1]<<"\n";

    double cr_normalny = corr(Normalny_2D);
    std::cout<<"wsp. korelacji: "<<cr_normalny<<"\n";

}