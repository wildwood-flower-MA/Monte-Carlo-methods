#ifndef DFSION_H
#define DFSION_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <random>

#define PI 3.141592653589793
using vector = std::vector<double>;
using matrix = std::vector<vector>;

struct Parametry {

    // PUNKT KREACJI
    double x_birth; // położenie w "x" punktu kreacji
    double y_birth; // położenie w "y" punktu kreacji
    double dn; // liczba dodanych cząstek na krok czasowy
    int N_max; // maksymalna liczba cząstek w obszarze
    int N_0; // startowa liczba cząstek w obszarze

    // OBSZAR DYFUZJI
    double x_diff; // położenie w "x" obszaru dyfuzji
    double y_diff; // położenie w "y" obszaru dyfuzji
    double R; // promień obszaru

    // OBSZAR ZNIKANIA
    double x_death; // położenie w "x" środka obszaru znikania
    double y_death; // położenie w "y" środka obszaru znikania
    double r_death; // promień obszaru znikania

    // DYFUZJA
    double D; // stała dyfuzji
    double dt; // krok czasowy
    int N; // liczba Kraków czasowych
};

void vector_rotation(double, double, double, double&, double&, double&, double);
double particle_translation(const std::vector<double>&, std::vector<double>&, const Parametry&);

matrix make_sand(const Parametry&);
matrix make_wyniki(const Parametry&);
void spill_sand_to_a_file(std::ofstream&, const matrix&, const Parametry&);
void spill_wyniki_to_a_file(std::ofstream&, const matrix&, const Parametry&);
inline std::vector<double> delta_r(const Parametry&);
void roam(std::vector<double>&, const Parametry&);
void dyfuzja(matrix&, matrix&, const Parametry&);

#endif