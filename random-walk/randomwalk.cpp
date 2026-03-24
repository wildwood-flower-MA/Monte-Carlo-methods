#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "randomwalk.h"

inline double c1(double x){ // -3 -> 3: 6
    return 1.0 + std::tanh(x);
}

inline double c2(double x){ // 0 -> 10: 1.4711276
    return 1.0/(1.0 + x*x);
}

inline double c3(double x){ // 0 -> 1: 0.24609375
    return pow(std::cos(PI*x), 10.0);
}

double u(){
    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(generator);
}

std::vector<int> histogramuj(vector& sequence, int N){
    // N bins

    std::vector<int> histogram = std::vector<int>(N, 0);
    int bin;
    for(double element : sequence){
        bin = static_cast<int>(N*element)%N;
        histogram[bin]++;
    }

    return histogram;
}

// METODA PODSTAWOWA
result podstawowa(int N, double a, double b, double (*prob_dist)(double)){
    // returns average value (integral) and variance of a given
    // probability distribution prob_dist
    // from a sample of length N, in a range a -> b 

    double avg = 0.0;
    double avg2 = 0.0;
    double U;
    vector sequence;
    sequence.reserve(N);

    vector us;
    us.reserve(N);

    for(int i = 0; i < N; i++){
        U = u();
        us.push_back(U);
        sequence.push_back(prob_dist(a + (b - a)*U));
        avg += sequence.back();    
        avg2 += sequence.back()*sequence.back();
    }
    avg *= (b-a)/N;
    avg2 *= (b-a)*(b-a)/N;

    return {avg, std::sqrt((avg2 - avg*avg)/N)/avg, histogramuj(us, 10), std::vector<int>(1, N)};
}

// METODA LOSOWANIA SYSTEMATYCZNEGO (nieoptymalne)
result systematyczne_nopt(int N, int M, double a, double b,
    double (*prob_dist)(double)){
    // N - liczba losowan; M - liczba podprzedzialow;
    // a -> b - przedzial calkowania; prob_dist - funkcja calkowana
    
    int N_m;
    double x_m, x_m1, x;
    double delta_x = (b-a)/static_cast<double>(M);
    vector averages = std::vector<double>(M, 0.0);
    vector averages2 = std::vector<double>(M, 0.0);
    vector variances = std::vector<double>(M, 0.0);
    double g_bar = 0.0;
    double var = 0.0;
    vector sequence;
    sequence.reserve(N);

    for(int m = 1; m <= M; m++){

        x_m = a + delta_x*((double)m-1.0);
        x_m1 = x_m + delta_x;
        N_m = static_cast<int>(static_cast<double>(N)/M);
        for(int i = 0; i < N_m; i++){

            x = x_m + (x_m1 - x_m)*u();
            sequence.push_back((b-a)*prob_dist(x));
            averages[m-1] += sequence.back()/(double)N_m;
            averages2[m-1] += sequence.back()*sequence.back()/(double)N_m;
            variances[m-1] = averages2[m-1] - averages[m-1]*averages[m-1];
        }

        g_bar += averages[m-1]/(double)M;
        var += variances[m-1]/(double)M/(double)M/(double)N_m;
    }

    return {g_bar, std::sqrt(var)/g_bar, histogramuj(sequence, 10), std::vector<int>(M, N_m)};
}

// METODA LOSOWANIA SYSTEMATYCZNEGO (optymalne)
result systematyczne_opt(int N, int N_pre, int M, double a, double b,
    double (*prob_dist)(double)){
    // N - liczba losowan; N_pre - liczba losowan w celu okreslenia sigma_m;
    // M - liczba podprzedzialow; a -> b - przedzial calkowania; prob_dist - funkcja calkowana
    
    int N_m;
    double x_m, x_m1, x;
    double delta_x = (b-a)/M;
    vector averages = std::vector<double>(M, 0.0);
    vector averages2 = std::vector<double>(M, 0.0);
    vector variances = std::vector<double>(M, 0.0);
    double g_bar = 0.0;
    double var = 0.0;
    vector sigma_ms = vector(M, 0.0);
    double sigma_ms_suma = 0.0;
    std::tuple<double, double, std::vector<int>, std::vector<int>> avg_var;
    vector sequence;
    sequence.reserve(N);
    std::vector<int> Ns(M, 0);
    
    for(int m = 0; m < M; m++){

        x_m = a + delta_x*m;
        x_m1 = x_m + delta_x;
        avg_var = systematyczne_nopt(N_pre, 10, x_m, x_m1, prob_dist);
        sigma_ms[m] = std::get<1>(avg_var) * std::get<0>(avg_var);
        sigma_ms_suma += sigma_ms[m];
    }

    for(int m = 1; m <= M; m++){
        x_m = a + delta_x*(m-1);
        x_m1 = x_m + delta_x;
        N_m = sigma_ms[m-1]/sigma_ms_suma*(double)N;
        if(N_m == 0) {
            N_m = 1;
        }

        Ns[m-1] += N_m;
        for(int i = 0; i < (int)N_m; i++){

            x = x_m + (x_m1 - x_m)*u();
            sequence.push_back((b-a)*prob_dist(x));
            averages[m-1] += sequence.back()/N_m;
            averages2[m-1] += sequence.back()*sequence.back()/N_m;
        }

        variances[m-1] = averages2[m-1] - averages[m-1]*averages[m-1];
        g_bar += averages[m-1]/(double)M;
        var += variances[m-1]/(double)M/(double)M/N_m;
    }

    return {g_bar, std::sqrt(var)/g_bar, histogramuj(sequence, 10), Ns};
}

void wypisz_wyniki(result& wyniki){

    double integral = std::get<0>(wyniki);
    double rel_error = std::get<1>(wyniki);
    std::vector<int> hist = std::get<2>(wyniki);
    std::vector<int> Ns = std::get<3>(wyniki);

    std::cout << "wartosc calki = " << integral << " | blad wzgledny = " << rel_error << "\n";
    std::cout << "histogram: ";
    for(int bin : hist){
        std::cout << bin << ", ";
    }
    std::cout << "\nliczba losowan w binie: ";
    for(int n : Ns){
        std::cout << n << ", ";
    }
}

int main(){

std::vector<int> K = {2, 3, 4, 5};
Result C1_podst, C1_syst, C1_syst_opt;
Result C2_podst, C2_syst, C2_syst_opt;
Result C3_podst, C3_syst, C3_syst_opt;

for(int k : K){
    C1_podst.push_back(podstawowa(pow(10,k), -3.0, 3.0, c1));
    C1_syst.push_back(systematyczne_nopt(pow(10,k), 10, -3.0, 3.0, c1));
    if(k>=1000){
        C1_syst_opt.push_back(systematyczne_opt(pow(10,k), 1000, 10, -3.0, 3.0, c1));
    } else {
        C1_syst_opt.push_back(systematyczne_opt(pow(10,k), 100, 10, -3.0, 3.0, c1));
    }
}


// calkowanie c2
for(int k : K){
    C2_podst.push_back(podstawowa(pow(10,k), 0.0, 10.0, c2));
    C2_syst.push_back(systematyczne_nopt(pow(10,k), 10, 0.0, 10.0, c2));
    if(k>=3){
        C2_syst_opt.push_back(systematyczne_opt(pow(10,k), 1000, 10, 0.0, 10.0, c2)); // zadanie 3.
    } else {
        C2_syst_opt.push_back(systematyczne_opt(pow(10,k), 100, 10, 0.0, 10.0, c2)); // zadanie 3.
    }
}

// calkowanie c3
for(int k : K){
    C3_podst.push_back(podstawowa(pow(10,k), 0.0, 1.0, c3)); // zadanie 1.
    C3_syst.push_back(systematyczne_nopt(pow(10,k), 10, 0.0, 1.0, c3)); // zadanie 2.
    if(k>=1000){
        C3_syst_opt.push_back(systematyczne_opt(pow(10,k), 1000, 10, 0.0, 1.0, c3)); // zadanie 3.
    } else {
        C3_syst_opt.push_back(systematyczne_opt(pow(10,k), 100, 10, 0.0, 1.0, c3)); // zadanie 3.
    }
}


std::cout << "Metoda podstawowa\n"; // METODA PODSTAWOWA

std::cout << "C1:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C1_podst[i]);
    std::cout << "\n";
}

std::cout << "C2:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C2_podst[i]);
    std::cout << "\n";
}
std::cout << "C3:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C3_podst[i]);
    std::cout << "\n";
}

std::cout << "C1:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C1_syst[i]);
    std::cout << "\n";
}
std::cout << "C2:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C2_syst[i]);
    std::cout << "\n";
}
std::cout << "C3:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C3_syst[i]);
    std::cout << "\n";
}

std::cout << "C1:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C1_syst_opt[i]);
    std::cout << "\n";
}
std::cout << "C2:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C2_syst_opt[i]);
    std::cout << "\n";
}
std::cout << "C3:\n";
for(int i = 0; i<K.size(); i++){
    std::cout << "k = "<< K[i] << "\n";
    wypisz_wyniki(C3_syst_opt[i]);
    std::cout << "\n";
}


}