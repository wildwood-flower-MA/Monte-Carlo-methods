#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <string>
#include <chrono>
#include "fullerene.h"

inline float dis(){
    
    static std::mt19937_64 gen{ std::random_device{}() };
    static std::uniform_real_distribution<float> dist{ 0.0, 1.0 };
    return dist(gen);
}

inline float f_cut(const float r, const Parametry& pm){

    if(r <= pm.R1){
        return 1.0;
    } else if(r <= pm.R2){
        return 0.5*(1.0 + std::cos(M_PI*(r - pm.R1)/(pm.R2 - pm.R1)));
    } else {
        return 0.0;
    }
}

float B(int i, int j, matrix& vectors, const Parametry& pm){
    
    auto compute_B = [&](int a, int b){

        vector r_ab = {vectors[b][0] - vectors[a][0],
                    vectors[b][1] - vectors[a][1],
                    vectors[b][2] - vectors[a][2]};

        float r_ab_len_2 = r_ab[0]*r_ab[0] + r_ab[1]*r_ab[1] + r_ab[2]*r_ab[2];

        //if (r_ab_len_2 > pm.R2*pm.R2){
        //    return 1.0;
        //}

        float r_ab_len = std::sqrt(r_ab_len_2);

        float zeta = 0.0;
        for (int k = 0; k < pm.n; ++k){

            if (k == a || k == b){
                continue;
            }

            vector r_ak = {vectors[k][0] - vectors[a][0],
                        vectors[k][1] - vectors[a][1],
                        vectors[k][2] - vectors[a][2]};
            
            float r_ak_len = std::sqrt(r_ak[0]*r_ak[0] + r_ak[1]*r_ak[1] + r_ak[2]*r_ak[2]);
            if (r_ak_len > pm.R2){
                continue;
            }
            float cos_theta = (r_ab[0]*r_ak[0] + r_ab[1]*r_ak[1] + r_ab[2]*r_ak[2])/(r_ab_len*r_ak_len);

            float g_theta = pm.a0*(1.0 + pm.c0*pm.c0/pm.d0/pm.d0
                - pm.c0*pm.c0/(pm.d0*pm.d0 + (1.0 + cos_theta)*(1.0 + cos_theta)));
            
            if (pm.version == 1){

                    if(cos_theta > 0){ 
                        zeta = 10.0; // UNIEMOZLIWIA 4 WIAZANIA 
                    } else {
                        zeta += f_cut(r_ak_len, pm)*g_theta;
                    }

                } else {

                    zeta += f_cut(r_ak_len, pm)*g_theta;

                }

        }

        return std::pow(1.0 + zeta, -pm.delta);
    };

    float B_ij = compute_B(i, j);
    float B_ji = compute_B(j, i);

    return 0.5*(B_ij + B_ji);
}

float V_i(int i, matrix& vectors, const Parametry& pm,
            float vRj_constant, float vAj_constant){

    float v_i = 0.0;
    for(int j = 0; j < pm.n; j++){

        if(i == j){
            continue;
        }

        vector r_ij = {vectors[j][0] - vectors[i][0],
                    vectors[j][1] - vectors[i][1],
                    vectors[j][2] - vectors[i][2]};

        float r2 = (r_ij[0]*r_ij[0]
                    + r_ij[1]*r_ij[1]
                    + r_ij[2]*r_ij[2]);

        if(r2 > pm.R2*pm.R2){
            continue;
        }

        float r = std::sqrt(r2);
        float v_R_j = pm.De*std::exp(vRj_constant*(r - pm.R0))/(pm.S - 1.0);
        float v_A_j = pm.De*pm.S*std::exp(vAj_constant*(r - pm.R0))/(pm.S - 1.0);
        v_i += f_cut(r, pm)*(v_R_j -  B(i, j, vectors, pm)*v_A_j); 

    }

    return v_i;
}

float V_tot(matrix& vectors, const Parametry& pm, float vRj_constant, float vAj_constant){

    float v_tot = 0.0;
    for(int i = 0; i < pm.n; i++){
        v_tot += 0.5*V_i(i, vectors, pm, vRj_constant, vAj_constant);
    }

    return v_tot;
}

vector PCF(matrix& vectors, const Parametry& pm){

    vector pcf(pm.M, 0.0);

    float r_av = 0.0;
    for(vector v : vectors){
        r_av += std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    r_av /= pm.n;
    float dr = 2.5*r_av/pm.M;

    for(int i = 0; i < pm.n; i++){        
        for(int j = i + 1; j < pm.n; j++){

            float dx = vectors[j][0] - vectors[i][0];
            float dy = vectors[j][1] - vectors[i][1];
            float dz = vectors[j][2] - vectors[i][2];
            float r = std::sqrt(dx*dx + dy*dy + dz*dz);

            int m = static_cast<int>(std::floor(r/dr));
            
            if(m < pm.M){                
                pcf[m] += 4.0*r_av*r_av/pm.n/pm.n/r/dr;
            }
        }
    }
    
    return pcf;
}

// od teraz wektory w postaci [x, y, z, r, phi, theta]
inline void update_vec_cartesian(vector& vec){

    vec[3] = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    vec[4] = std::atan2(vec[1], vec[0]);
    vec[5] = std::acos(vec[2]/vec[3]);
    
}

inline void update_vec_radial(vector& vec){

    vec[0] = vec[3]*std::sin(vec[5])*std::cos(vec[4]); 
    vec[1] = vec[3]*std::sin(vec[5])*std::sin(vec[4]);
    vec[2] = vec[3]*std::cos(vec[5]);

}

void update_particle(matrix& vectors, const Parametry& pm, int i, float beta, float vRj_constant, float vAj_constant){

    int size = pm.n;

    float V_old = V_i(i, vectors, pm, vRj_constant, vAj_constant);
    vector vector_old = vectors[i];

    vectors[i][3] += vectors[i][3]*(2.0*dis() - 1.0)*pm.w_r;
    vectors[i][4] += vectors[i][4]*(2.0*dis() - 1.0)*pm.w_phi;
    vectors[i][5] += vectors[i][5]*(2.0*dis() - 1.0)*pm.w_theta;

    if(vectors[i][4] > 2.0*M_PI){
        vectors[i][4] -= 2.0*M_PI;
    }
    if(vectors[i][4] < 0.0){
        vectors[i][4] += 2.0*M_PI;
    }

    if(vectors[i][5] < 0.0 || vectors[i][5] > M_PI){
        vectors[i][5] = vector_old[5];
    }

    update_vec_radial(vectors[i]);

    float p_acc = std::min(1.0f, std::exp(-beta*(V_i(i, vectors, pm, vRj_constant, vAj_constant) - V_old)));
    if (dis() > p_acc){
        vectors[i] = vector_old;
    }

}

void update_global_radius(matrix& vectors, const Parametry& pm, float beta, float& V,
                        float vRj_constant, float vAj_constant){

    vector old_radiuses(pm.n);
    for(int i = 0; i < pm.n; i++){
        old_radiuses[i] = vectors[i][3];
    }
    //float V_here = V_tot(vectors, pm, vRj_constant, vAj_constant);

    float U = dis();
    for(int i = 0; i < pm.n; i++){
        vectors[i][3] = vectors[i][3]*(1.0 + pm.W_all*(2.0*U - 1.0));
        update_vec_radial(vectors[i]);
    }

    float V_current = V_tot(vectors, pm, vRj_constant, vAj_constant);
    float p_acc = std::min(1.0f, std::exp(-beta*(V_current - V)));

    if(dis() > p_acc){
        for(int i = 0; i < pm.n; i++){
            vectors[i][3] = old_radiuses[i];
            update_vec_radial(vectors[i]);
        }
    } else {    
        V = V_current;
    }

}

void SA_algorithm(matrix& vectors, const Parametry& pm){

    float vRj_constant = -std::sqrt(2.0*pm.S)*pm.lambda;
    float vAj_constant = -std::sqrt(2.0/pm.S)*pm.lambda;

    std::ofstream file_energy(pm.file_energy);
    std::ofstream file_vectors(pm.file_vectors);
    std::ofstream file_pcf(pm.file_pcf);

    float V = V_tot(vectors, pm, vRj_constant, vAj_constant);

    for(int it = 0; it < pm.it_max; it++){

        std::cout << "\riteracja: " << it << "/" << pm.it_max << std::flush;

        float beta = pm.beta_min + (pm.beta_max - pm.beta_min)*std::pow(static_cast<float>(it + 1)/pm.it_max, pm.p);

        for(int i = 0; i < pm.n; i++){
            update_particle(vectors, pm, i, beta, vRj_constant, vAj_constant);
        }
        update_global_radius(vectors, pm, beta, V, vRj_constant, vAj_constant);

        if(it % 100 == 0){
            file_energy << V << std::endl;
            for(auto& vector : vectors){
                file_vectors << vector[0] << " " << vector[1] << " " << vector[2] << std::endl;
            }
            vector pcf = PCF(vectors, pm);
            for(float bin : pcf){
                file_pcf << bin << " "; 
            }
            file_pcf << std::endl;
        }
    }

    file_energy.close();
    file_vectors.close();
    file_pcf.close();
}

matrix starting_positions(const Parametry& pm){

    matrix positions(pm.n, vector(6, 0.0f));
    for (int i = 0; i < pm.n; i++){

        positions[i][3] = pm.r_0;


        //positions[i][4] = 2.0f*M_PI*dis();
        //positions[i][5] = std::acos(1.0f - 2.0f*dis());

        positions[i][4] = 2.0f*M_PI*dis();
        positions[i][5] = M_PI*dis();
        update_vec_radial(positions[i]);
    }

    return positions;
}

int main(){

    Parametry pm;

    pm.R0 = 1.315;
    pm.R1 = 1.7;
    pm.R2 = 2.0;
    pm.De = 6.325;
    pm.S = 1.29;
    pm.lambda = 1.5;
    pm.delta = 0.80469;
    pm.a0 = 0.011304;
    pm.c0 = 19.0;
    pm.d0 = 2.5;

    pm.M = 100;
    pm.w_r = 0.0001;
    pm.w_phi = 0.05;
    pm.w_theta = 0.05;
    pm.beta_max = 100.0;
    pm.beta_min = 1.0;
    pm.W_all = 0.0001;
    pm.p = 2.0;
    pm.n = 60;
    pm.r_0 = 3.50;
    pm.it_max = 100000;
    pm.version = 0;

    // ZADANIE 1.
    float vRj_constant = -std::sqrt(2.0*pm.S)*pm.lambda;
    float vAj_constant = -std::sqrt(2.0/pm.S)*pm.lambda;

    // Zadanie 1.
    matrix atoms_positions_c60;
    std::ifstream input_file("atoms_positions_c60.dat");
    float x, y, z;
    while(input_file >> x >> y >> z){
        atoms_positions_c60.push_back({x, y, z});
    }
    input_file.close();
    float E_total = V_tot(atoms_positions_c60, pm, vRj_constant, vAj_constant);
    std::cout << E_total << std::endl;

    
    /*
    // Zadanie 2. // ZMIENIĆ WARUNEK PRZY LICZENIU B
    pm.version = 0;
    pm.file_energy = "zad2_energy.dat";
    pm.file_vectors = "zad2_vectors.dat";
    pm.file_pcf = "zad2_pcf.dat";
    auto start_1 = std::chrono::high_resolution_clock::now();

    matrix atoms_positions_zad2 = starting_positions(pm);
    SA_algorithm(atoms_positions_zad2, pm);

    auto end_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_1 = end_1 - start_1;
    std::cout << "\n Runtime zad. 2.: " << elapsed_1.count() << " s" << std::endl;
    
    // ZADANIE 3.
    pm.version = 1;
    pm.file_energy = "zad3_energy.dat";
    pm.file_vectors = "zad3_vectors.dat";
    pm.file_pcf = "zad3_pcf.dat";
    auto start_2 = std::chrono::high_resolution_clock::now();

    matrix atoms_positions_zad3 = starting_positions(pm);
    SA_algorithm(atoms_positions_zad3, pm);

    auto end_2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_2 = end_2 - start_2;
    std::cout << "\n Runtime zad. 3.: " << elapsed_2.count() << " s" << std::endl;
    

    // ZADANIE 4.
    pm.version = 1;
    pm.r_0 = 2.50;
    pm.file_energy = "zad4_energy.dat";
    pm.file_vectors = "zad4_vectors.dat";
    pm.file_pcf = "zad4_pcf.dat";
    auto start_3 = std::chrono::high_resolution_clock::now();

    matrix atoms_positions_zad3 = starting_positions(pm);
    SA_algorithm(atoms_positions_zad3, pm);

    auto end_3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_3 = end_3 - start_3;
    std::cout << "\n Runtime zad. 4.: " << elapsed_3.count() << " s" << std::endl;
    */
    pm.version = 1;
    pm.r_0 = 2.50;

    // ZAD 6.1
    pm.beta_max = 50.0; // ZMIENIŁEM -zad6
    pm.beta_min = 0.5; // ZMIENIŁEM -zad6
    pm.p = 1.5; // ZMIENIŁEM -zad6

    pm.file_energy = "zad6_1_energy.dat";
    pm.file_vectors = "zad6_1_vectors.dat";
    pm.file_pcf = "zad6_1_pcf.dat";
    auto start_6_1 = std::chrono::high_resolution_clock::now();

    matrix atoms_positions_zad6_1 = starting_positions(pm);
    SA_algorithm(atoms_positions_zad6_1, pm);

    auto end_6_1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_6_1 = end_6_1 - start_6_1;
    std::cout << "\n Runtime zad. 6.1.: " << elapsed_6_1.count() << " s" << std::endl;


    // ZAD 6.2
    pm.beta_max = 25.0; // ZMIENIŁEM -zad6
    pm.beta_min = 0.2; // ZMIENIŁEM -zad6
    pm.p = 1; // ZMIENIŁEM -zad6

    pm.file_energy = "zad6_2_energy.dat";
    pm.file_vectors = "zad6_2_vectors.dat";
    pm.file_pcf = "zad6_2_pcf.dat";
    auto start_6_2 = std::chrono::high_resolution_clock::now();

    matrix atoms_positions_zad6_2 = starting_positions(pm);
    SA_algorithm(atoms_positions_zad6_2, pm);

    auto end_6_2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_6_2 = end_6_2 - start_6_2;
    std::cout << "\n Runtime zad. 6.2.: " << elapsed_6_2.count() << " s" << std::endl;

}