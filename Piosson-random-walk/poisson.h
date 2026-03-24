#ifndef POISSON_H
#define POISSON_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <random>

double U();

void save_to_file(const std::string& filename, const std::vector<std::vector<double>>& matrix);

double _rho(double x, double y, 
            double rho_max,
            double x_max, double y_max,
            double sigma_rho);

std::vector<std::vector<double>> metoda_relaksacyjna(int nx, int ny, double delta,
                                                    double epsilon, double omega,
                                                    double tol, int itmax,
                                                    double rho_max, double sigma_rho,
                                                    double VL, double VB, double VT);

std::vector<std::vector<std::vector<double>>> metoda_MonteCarlo(int nx, int ny, double delta,
                                                                double epsilon, int Nchains,
                                                                int nlength, int block_B,
                                                                double rho_max, double sigma_rho,
                                                                double VL, double VB, double VT);

#endif