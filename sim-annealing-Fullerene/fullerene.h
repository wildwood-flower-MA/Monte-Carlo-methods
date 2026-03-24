#ifndef FULLERENE_H
#define FULLERENE_H

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <random>

using vector = std::vector<float>;
using matrix = std::vector<vector>;

struct Parametry{

    float R0;
    float R1;
    float R2;
    float De;
    float S;
    float lambda;
    float delta;
    float a0;
    float c0;
    float d0;

    float M;

    float w_r;
    float w_phi;
    float w_theta;

    float beta;
    float beta_max;
    float beta_min;

    float W_all;
    float p;
    float it_max;

    int n;
    float r_0;

    int version;

    std::string file_vectors;
    std::string file_energy;
    std::string file_pcf;

};

inline float f_cut(const float r, const Parametry& pm);
float B(int i, int j, matrix& vectors, const Parametry& pm);

float V_i(int i, matrix& vectors, const Parametry& pm, float vRj_constant, float vAj_constant);
float V_tot(matrix& vectors, const Parametry& pm, float vRj_constant, float vAj_constant);
vector PCF(matrix& vectors, const Parametry& pm);

inline void update_vec_cartesian(vector& vec);
inline void update_vec_radial(vector& vec);

void update_particle(matrix& vectors, const Parametry& pm, int i, float beta, float vRj_constant, float vAj_constant);
void update_global_radius(matrix& vectors, const Parametry& pm, float beta, float& V);
void SA_algorithm(matrix& vectors, const Parametry& pm);
matrix starting_positions(const Parametry& pm);

#endif