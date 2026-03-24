#ifndef MASTER_H
#define MASTER_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <random>

#define PI 3.141592653589793
using vector = std::vector<double>;
using matrix = std::vector<std::vector<double>>;

double U();
void evolution(int x1_0, int x2_0, int x3_0,
    double k1, double k2, double k3, double k4,
    int t_max, std::string filename);

#endif