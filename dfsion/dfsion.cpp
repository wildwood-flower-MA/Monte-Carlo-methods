#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "dfsion.h"

namespace cst {
    inline constexpr double pi = 3.14159265358979323846;
}

/////////////////////////////////////////////////////////////////
void vector_rotation(double tx, double ty, double tz, double& x, double& y, double& z, double teta){

    double p0, p1, p2, p3;
    double q0, q1, q2, q3;
    double r0, r1, r2, r3;

    double t = std::sqrt(tx*tx + ty*ty + tz*tz);
    double cs = std::cos(teta/2.0);
    double sn = std::sin(teta/2.0);

    p0 = cs;
    p1 = tx/t*sn;
    p2 = ty/t*sn;
    p3 = tz/t*sn;

    r0 = 0.0;
    r1 = x;
    r2 = y;
    r3 = z;

    q0 = r0*p0 - r1*(-p1) - r2*(-p2) - r3*(-p3);
    q1 = r1*p0 + r0*(-p1) - r3*(-p2) + r2*(-p3);
    q2 = r2*p0 + r3*(-p1) + r0*(-p2) - r1*(-p3);
    q3 = r3*p0 - r2*(-p1) + r1*(-p2) + r0*(-p3);

    r0 = p0*q0 - p1*q1 - p2*q2 - p3*q3;
    r1 = p1*q0 + p0*q1 - p3*q2 + p2*q3;
    r2 = p2*q0 + p3*q1 + p0*q2 - p1*q3;
    r3 = p3*q0 - p2*q1 + p1*q2 + p0*q3;

    x = r1;
    y = r2;
    z = r3;
}

double particle_translation(const std::vector<double>& grain_of_sand_before, std::vector<double>& grain_of_sand, const Parametry& pm){

    double xr = pm.x_diff;
    double yr = pm.y_diff;
    double Rr = pm.R;
    double xa = pm.x_death;
    double ya = pm.y_death;
    double Ra = pm.r_death;

    double x1 = grain_of_sand_before[0];
    double y1 = grain_of_sand_before[1];
    double x2 = grain_of_sand[0];
    double y2 = grain_of_sand[1];
    double& exist = grain_of_sand[2];
    double a, b, c, delta, beta1, beta2, alfa1, alfa2, alfa;

    a = std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2);
    b = 2*((x2 - x1)*(x1 - xa) + (y2 - y1)*(y1 - ya));
    c = std::pow(x1 - xa, 2) + std::pow(y1 - ya, 2) - std::pow(Ra, 2);
    delta = b*b - 4*a*c;
    if(delta >= 0){
        beta1 = (-b - std::sqrt(delta))/(2*a);
        beta2 = (-b + std::sqrt(delta))/(2*a);
        if((beta1 >= 0 && beta1 <= 1) || (beta2 >= 0 && beta2 <= 1)){
            exist = 0.0;
            return 0.0;
        }
    }

    alfa = -1.0;
    a = std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2);
    b = 2*((x2 - x1)*(x1 - xr) + (y2 - y1)*(y1 - yr));
    c = std::pow(x1 - xr, 2) + std::pow(y1 - yr, 2) - std::pow(Rr, 2);
    delta = b*b - 4*a*c;
    if(delta >= 0){
        alfa1 = (-b - std::sqrt(delta))/(2*a);
        alfa2 = (-b + std::sqrt(delta))/(2*a);
        if(alfa1 >= 0 && alfa1 <= 1){
            alfa = alfa1;
        }else if(alfa2 >= 0 && alfa2 <= 1){
            alfa = alfa2;
        }
    }

    if(alfa < 0){
        grain_of_sand[0] = x2;
        grain_of_sand[1] = y2;
        return 0.0;
    }else{
        double x3 = x1 + alfa*(x2 - x1);
        double y3 = y1 + alfa*(y2 - y1);
        double z3 = 0.0;

        double x13 = x1 - x3;
        double y13 = y1 - y3;
        double z13 = 0.0;
        double norm = std::sqrt(std::pow(x13, 2) + std::pow(y13, 2) + std::pow(z13, 2));
        x13 /= norm;
        y13 /= norm;
        z13 /= norm;

        double tx = 0.0 - x3;
        double ty = 0.0 - y3;
        double tz = 0.0 - z3;

        vector_rotation(tx, ty, tz, x13, y13, z13, cst::pi);

        double length = std::sqrt(std::pow(x2 - x3, 2) + std::pow(y2 - y3, 2));
        grain_of_sand[0] = x3 + x13*1.0E-6;
        grain_of_sand[1] = y3 + y13*1.0E-6;
        grain_of_sand[0] = x3 + length*x13;
        grain_of_sand[1] = y3 + length*y13;

        return length;
    }
}
/////////////////////////////////////////////////////////////////

matrix make_sand(const Parametry& pm){
    matrix every_grain_of_sand(pm.N_max, vector{pm.x_birth, pm.y_birth, 0.0});
    for(int i = 0; i < pm.N_0; i++){
        every_grain_of_sand[i][2] = 1.0;
    }
    return every_grain_of_sand;
}

matrix make_wyniki(const Parametry& pm){
    matrix wyniki(pm.N, vector{0.0, 0.0, 0.0, 0.0});
    return wyniki;
}

void spill_sand_to_a_file(std::ofstream& file, const matrix& every_grain_of_sand, const Parametry& pm){
    file << pm.x_birth << " " << pm.y_birth << " " << pm.dn << " " << pm.N_max << " "
        << pm.x_diff << " " << pm.y_diff << " " << pm.R << " "
        << pm.x_death << " " << pm.y_death << " " << pm.r_death << " "
        << pm.D << " " << pm.dt << " " << pm.N << "\n";
    for(const vector& grain_of_sand : every_grain_of_sand){
        if(grain_of_sand[2] == 1.0){
            file << grain_of_sand[0] << " " << grain_of_sand[1] << "\n";
        }
    }
}

void spill_wyniki_to_a_file(std::ofstream& file, const matrix& wyniki, const Parametry& pm){
    file << pm.x_birth << " " << pm.y_birth << " " << pm.dn << " " << pm.N_max << " "
        << pm.x_diff << " " << pm.y_diff << " " << pm.R << " "
        << pm.x_death << " " << pm.y_death << " " << pm.r_death << " "
        << pm.D << " " << pm.dt << " " << pm.N << "\n";
    for(const vector& wynik : wyniki){
        file << wynik[0] << " " << wynik[1] << " " << wynik[2] << " " << wynik[3] << "\n";
    }
}

inline vector delta_r(const Parametry& pm){
    double sigma = std::sqrt(2.0*pm.D*pm.dt);

    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

    double rand_1 = dist(generator);
    double rand_2 = dist(generator);
    vector dr = {sigma*std::sqrt(-2.0*std::log(1-rand_1))*std::cos(2.0*cst::pi*rand_2),
                sigma*std::sqrt(-2.0*std::log(1-rand_1))*std::sin(2.0*cst::pi*rand_2)};
    return dr;
}

void roam(vector& grain_of_sand, const Parametry& pm){
    vector grain_of_sand_copy = grain_of_sand;
    vector dr = delta_r(pm);
    grain_of_sand[0] += dr[0];
    grain_of_sand[1] += dr[1];
    double length = 1.0;
    int safety = 0;

    // Zmiana: Zabezpieczenie przed uwiezieniem czastki w nieskonczonej petli odbic
    while(length > 1e-6 && safety++ < 10){
        length = particle_translation(grain_of_sand_copy, grain_of_sand, pm);
    }
}

void dyfuzja(matrix& every_grain_of_sand, matrix& wyniki, const Parametry& pm){
    for(int t = 0; t < pm.N; t++){
        double av_x = 0.0;
        double av_y = 0.0;
        double av_xy = 0.0;
        double av_x2 = 0.0;
        double av_y2 = 0.0;

        int i = 0;
        
        int dn_int = static_cast<int>(pm.dn);
        for(vector& grain_of_sand : every_grain_of_sand){
            if(grain_of_sand[2] == 0.0){
                grain_of_sand[0] = pm.x_birth;
                grain_of_sand[1] = pm.y_birth;
                grain_of_sand[2] = 1.0;

                if(++i >= dn_int){
                    break;
                }
            }
        }

        for(vector& grain_of_sand : every_grain_of_sand){
            if(grain_of_sand[2] == 1.0){
                roam(grain_of_sand, pm);
            }
        }
        
        for(vector& grain_of_sand : every_grain_of_sand){
            wyniki[t][0] += grain_of_sand[2];
        }

        
        if(wyniki[t][0] > 0.0){
            for(vector& grain_of_sand : every_grain_of_sand){
                av_x += grain_of_sand[2]*grain_of_sand[0];
                av_y += grain_of_sand[2]*grain_of_sand[1];
                av_xy += grain_of_sand[2]*grain_of_sand[0]*grain_of_sand[1];
                av_x2 += grain_of_sand[2]*grain_of_sand[0]*grain_of_sand[0];
                av_y2 += grain_of_sand[2]*grain_of_sand[1]*grain_of_sand[1];
            }
            
            av_x /= wyniki[t][0];
            av_y /= wyniki[t][0];
            av_xy /= wyniki[t][0];
            av_x2 /= wyniki[t][0];
            av_y2 /= wyniki[t][0];

            double Dxx = 0.5*(av_x2 - av_x*av_x)/(t+1)/pm.dt;
            double Dyy = 0.5*(av_y2 - av_y*av_y)/(t+1)/pm.dt;
            double Dxy = 0.5*(av_xy - av_x*av_y)/(t+1)/pm.dt;

            wyniki[t][1] = Dxx;
            wyniki[t][2] = Dyy;
            wyniki[t][3] = Dxy;
        }
    }
}

int main(){
    
    Parametry pm1;

    pm1.x_birth = -4.5;
    pm1.y_birth = 0.0;
    pm1.dn = 50*0.1;
    pm1.N_max = 10000;
    pm1.N_0 = 0;

    pm1.x_diff = 0.0;
    pm1.y_diff = 0.0;
    pm1.R = 5.0;

    pm1.x_death = 3.0;
    pm1.y_death = 0.0;
    pm1.r_death = 0.5;

    pm1.D = 1.0;
    pm1.dt = 0.1;
    pm1.N = 5/pm1.dt;

    matrix every_grain_of_sand_1 = make_sand(pm1);
    matrix wyniki_1 = make_wyniki(pm1);
    dyfuzja(every_grain_of_sand_1, wyniki_1, pm1);

    std::ofstream zadanie1("zadanie1.dat");
    spill_sand_to_a_file(zadanie1, every_grain_of_sand_1, pm1);
    zadanie1.close();

    std::ofstream zadanie1_wyniki("zadanie1_wyniki.dat");
    spill_wyniki_to_a_file(zadanie1_wyniki, wyniki_1, pm1);
    zadanie1_wyniki.close();

    return 0;
}