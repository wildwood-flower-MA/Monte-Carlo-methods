#ifndef CIRCLES_H
#define CIRCLES_H
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

vector random_position(double, double, double);
int circle(double, double, double, double, double);
matrix common_part(double, double, double,
    double, double, double,
    std::initializer_list<int>);

void plot_circles(std::initializer_list<matrix>);
void plot_plots(std::initializer_list<matrix>);

#endif
