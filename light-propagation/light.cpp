#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <string>
#include "light.h"

int main(){

    DYFUZJA_FOTONOW_2D ph_diff;

    // WARSTWY
    ph_diff.liczba_warstw = 3;

    // ------( legenda )------
    // 0 - mu_a (absorpcja)
    // 1 - mu_s (rozpraszanie)
    // 2 - d (grubosc)
    // 3 - g (anizotropia)
    // 4 - n (wsp. zalamania)
    // -----------------------

    // bazowy zestaw pm
    /*
    ph_diff.xmax = 0.2;
    ph_diff.x_zrodla = 0.1;
    ph_diff.dx_zrodla = 0.0;
    ph_diff.x_detekcji = 0.15;
    ph_diff.dx_detekcji = 0.01;
    ph_diff.nx = 100;
    ph_diff.ny = 100;
    ph_diff.rx0 = 0.0;
    ph_diff.ry0 = 1.0;

    ph_diff.dane_warstw[1][0] = 1.0;
    ph_diff.dane_warstw[1][1] = 10.0;
    ph_diff.dane_warstw[1][2] = 0.02;
    ph_diff.dane_warstw[1][3] = 0.75;
    ph_diff.dane_warstw[1][4] = 1.3;

    ph_diff.dane_warstw[2][0] = 1.0;
    ph_diff.dane_warstw[2][1] = 190.0;
    ph_diff.dane_warstw[2][2] = 0.02;
    ph_diff.dane_warstw[2][3] = 0.075;
    ph_diff.dane_warstw[2][4] = 1.0;

    ph_diff.dane_warstw[3][0] = 10.0;
    ph_diff.dane_warstw[3][1] = 90.0;
    ph_diff.dane_warstw[3][2] = 0.02;
    ph_diff.dane_warstw[3][3] = 0.95;
    ph_diff.dane_warstw[3][4] = 1.0;
    */

    // g++ -O3 -march=native -flto light.cpp -o light.exe

    ph_diff.xmax = 0.2;
    ph_diff.x_zrodla = 0.1;
    ph_diff.dx_zrodla = 0.0;
    ph_diff.x_detekcji = 0.15;
    ph_diff.dx_detekcji = 0.01;
    ph_diff.nx = 100;
    ph_diff.ny = 100;
    ph_diff.rx0 = 0.0;
    ph_diff.ry0 = 1.0;

    ph_diff.dane_warstw[1][0] = 1.0;
    ph_diff.dane_warstw[1][1] = 10.0;
    ph_diff.dane_warstw[1][2] = 0.02;
    ph_diff.dane_warstw[1][3] = 0.75;
    ph_diff.dane_warstw[1][4] = 1.0;

    ph_diff.dane_warstw[2][0] = 10.0;
    ph_diff.dane_warstw[2][1] = 210.0;
    ph_diff.dane_warstw[2][2] = 0.02;
    ph_diff.dane_warstw[2][3] = 0.75;
    ph_diff.dane_warstw[2][4] = 1.5;

    ph_diff.dane_warstw[3][0] = 10.0;
    ph_diff.dane_warstw[3][1] = 90.0;
    ph_diff.dane_warstw[3][2] = 0.02;
    ph_diff.dane_warstw[3][3] = 0.95;
    ph_diff.dane_warstw[3][4] = 1.0;

    ph_diff.inicjalizuj();
    int N = 200000; // liczba wiązek fotonowych

    ph_diff.zapisz_wszystkie_sciezki = 0;
    ph_diff.zapisz_sciezki_detekcji = 0;

    for(int k=0; k<N; k++){
        ph_diff.pojedyncza_sciezka();
    }

    std::ofstream plik_absorpcja("absorption.dat");
    for(auto& row : ph_diff.absorpcja){
        for(auto element : row){
            plik_absorpcja << element << " ";
        }
        plik_absorpcja << std::endl;
    }

    std::ofstream plik_refl_transm("reflt_tramsm.dat");
    for(auto& el : ph_diff.odbicie){
            plik_refl_transm << el << " ";
    }
    plik_refl_transm << std::endl;
    for(auto& el : ph_diff.transmisja){
            plik_refl_transm << el << " ";
    }

    std::ofstream plik_pm("pm.dat");
    plik_pm << ph_diff.rx0 << " ";
    plik_pm << ph_diff.ry0 << " ";

    plik_pm << ph_diff.dane_warstw[1][4] << " "; // n
    plik_pm << ph_diff.dane_warstw[1][0] << " "; // mu_a
    plik_pm << ph_diff.dane_warstw[1][1] << " "; // mu_s

    plik_pm << ph_diff.dane_warstw[2][4] << " ";
    plik_pm << ph_diff.dane_warstw[2][0] << " ";
    plik_pm << ph_diff.dane_warstw[2][1] << " ";

    plik_pm << ph_diff.dane_warstw[3][4] << " ";
    plik_pm << ph_diff.dane_warstw[3][0] << " ";
    plik_pm << ph_diff.dane_warstw[3][1] << " ";

    plik_pm << ph_diff.dane_warstw[2][3] << " "; // g(2)

    return 0;
}