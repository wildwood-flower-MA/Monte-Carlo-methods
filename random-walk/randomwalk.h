#ifndef RANDOMWALK_H
#define RANDOMWALK_H
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
using result = std::tuple<double, double, std::vector<int>, std::vector<int>>;
using Result = std::vector<std::tuple<double, double, std::vector<int>, std::vector<int>>>;

inline double c1(double);
inline double c2(double);
inline double c3(double);
double u();
std::vector<int> histogramuj(vector&, int);
result podstawowa(int, double, double, double (*)(double));
result systematyczne_nopt(int, int, double, double, double (*)(double));
result systematyczne_opt(int, int, int, double, double, double (*)(double));
void wypisz_wyniki(std::tuple<double, double, std::vector<int>>&);

#endif
