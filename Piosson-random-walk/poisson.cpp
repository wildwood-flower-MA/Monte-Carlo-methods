#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <string>
#include "poisson.h"

// OGÓLNE

double U(){
    static std::random_device random;
    static std::mt19937 generator(random());
    static std::uniform_real_distribution<double> dist(0.0, 1.0);

    return dist(generator);
}

void save_to_file(const std::string& filename, const std::vector<std::vector<double>>& matrix){

    std::ofstream outfile(filename);
    if (!outfile) return;

    for(const auto& v : matrix){
        for(double element : v){
            outfile << element << " ";
        }
        outfile << std::endl;
    }
}

double _rho(double x, double y, 
            double rho_max,
            double x_max, double y_max,
            double sigma_rho){

    double exponented = -0.5*((x - x_max/2)*(x - x_max/2) + (y - y_max/2)*(y - y_max/2))/sigma_rho/sigma_rho;

    return std::exp(exponented);
}


// METODA RELAKSACYJNA

std::vector<std::vector<double>> metoda_relaksacyjna(int nx, int ny, double delta,
                                                    double epsilon, double omega,
                                                    double tol, int itmax,
                                                    double rho_max, double sigma_rho,
                                                    double VL, double VB, double VT){

    double Fold = 0;
    double Fnew = 0;

    std::vector<std::vector<double>> V(nx + 1, std::vector<double>(ny + 1, 0.0));
    std::vector<std::vector<double>> rho(nx + 1, std::vector<double>(ny + 1, 0.0));
    
    double x_max = delta*nx;
    double y_max = delta*ny;

    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){

            double x = delta*i;
            double y = delta*j;

            rho[i][j] = _rho(x, y, rho_max, x_max, y_max, sigma_rho);
            
        }
    }

    for(int j=0; j <= ny;j++){

        V[0][j] = VL*std::sin(M_PI*j/ny);

        for(int i = 0; i <= nx; i++){

            V[i][0] = VB*std::sin(M_PI*i/nx);
            V[i][ny] = VT*std::sin(M_PI*i/nx);

        }
    }

    for(int it = 1; it < itmax; it++){

        for(int i = 1; i < nx; i++){
            for(int j = 1; j < ny; j++){

                V[i][j] = (1 - omega)*V[i][j] + 0.25*omega*(  V[i + 1][j]
                                                            + V[i - 1][j]
                                                            + V[i][j + 1]
                                                            + V[i][j - 1]
                                                            + rho[i][j]*delta*delta/epsilon
                                                            );

            }
        } 
    
        // WB: Neumanna
        for(int j = 0; j < ny; j++){
            V[nx][j] = V[nx - 1][j];
        }

        Fold = Fnew;
        Fnew = 0;

        for(int i = 1; i < nx; i++){
            for(int j = 0; j < ny; j++){
                double Ex = (V[i + 1][j] - V[i - 1][j])/delta/2.0;
                double Ey = (V[i][j + 1] - V[i][j - 1])/delta/2.0;
                Fnew += (Ex*Ex + Ey*Ey)/2.0 - rho[i][j]*V[i][j];

            }
        }

        if(std::abs((Fnew - Fold)/Fnew) < tol){
            break;
        }

    }

    return V;
}

// METODA MONTE CARLO

std::vector<std::vector<std::vector<double>>> metoda_MonteCarlo(int nx, int ny, double delta,
                                                                double epsilon, int Nchains,
                                                                int nlength, int block_B,
                                                                double rho_max, double sigma_rho,
                                                                double VL, double VB, double VT){

    std::vector<std::vector<double>> V(nx + 1, std::vector<double>(ny + 1, 0.0));
    std::vector<std::vector<double>> sigma_V(nx + 1, std::vector<double>(ny + 1, 0.0));
    std::vector<std::vector<double>> B(nx + 1, std::vector<double>(ny + 1, 0.0));
    std::vector<std::vector<double>> S(nx + 1, std::vector<double>(ny + 1, 0.0));

    std::vector<std::vector<double>> rho(nx + 1, std::vector<double>(ny + 1, 0.0));
    double x_max = delta*nx;
    double y_max = delta*ny;

    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){

            double x = delta*i;
            double y = delta*j;

            rho[i][j] = _rho(x, y, rho_max, x_max, y_max, sigma_rho);
            
        }
    }

    for(int j=0; j <= ny;j++){

        V[0][j] = VL*std::sin(M_PI*j/ny);
        B[0][j] = 1.0;

        for(int i = 0; i <= nx; i++){

            V[i][0] = VB*std::sin(M_PI*i/nx);
            V[i][ny] = VT*std::sin(M_PI*i/nx);
            B[i][0] = 1.0;
            B[i][ny] = 1.0;

        }
    }

    for(int i0 = 1; i0 < nx; i0++){
        for(int j0 = 1; j0 < ny; j0++){

            double sum_V1 = 0.0;
            double sum_V2 = 0.0;
            int kchains = 0;

            // dla wezla generujemy N lancuchow
            for(int N = 1; N <= Nchains; N++){
                
                int i = i0;
                int j = j0;
                double g = 0.0;

                for(int n = 1; n <= nlength; n++){

                    // bladzenie
                    double U1 = U();
                    int m = static_cast<int>(std::floor(4*U1));
                    if(m == 0){
                        i--;
                    } else if(m == 1){
                        i++;
                    } else if(m == 2){
                        j--;
                    } else if(m == 3){
                        j++;
                    }
                    
                    // odbicie na prawym brzegu
                    if(i == nx + 1){
                        i = nx - 1;
                    }

                    // absorpcja - WB Dirichleta - koniec dobrego
                    if (B[i][j] == 1){

                        double dV = V[i][j] + g; // wzor (17)
                        sum_V1 += dV; // wzor (18)
                        sum_V2 += dV*dV; // wzor (19)

                        kchains++;
                        break;
                    }

                    g += 0.25*rho[i][j]*delta*delta/epsilon;
                }
            }

            double V1 = sum_V1/kchains;
            double V2 = sum_V2/kchains;
            V[i0][j0] = V1;
            sigma_V[i0][j0] = std::sqrt((V2 - V1*V1)/kchains);
            B[i0][j0] = block_B;
            S[i0][j0] = static_cast<double>(kchains)/static_cast<double>(Nchains);
        
        }
    }

    return {V, sigma_V, S};
}

int main(){

    int nx = 30;
    int ny = 30;
    double delta = 0.1;
    double epsilon = 1.0;
    double omega = 1.8;
    double tol = 1e-6;
    int itmax = 10000;
    double xmax = delta*nx;
    double ymax = delta*ny;
    double rho_max = 1.0;
    double sigma_rho = xmax/10.0;
    double VL = 1.0;
    double VB = -1.0;
    double VT = -1.0;

    std::vector<std::vector<double>> potencjal_dokl = metoda_relaksacyjna(nx, ny, delta,
                                                                        epsilon, omega,
                                                                        tol, itmax,
                                                                        rho_max, sigma_rho,
                                                                        VL, VB, VT);

    save_to_file("V_dokladnie.dat", potencjal_dokl);
    std::cout << "0. zrobione" << std::endl;

    int Nchains;
    int nlength;
    int block_B;

    Nchains = 100;
    nlength = 100;
    block_B = 0;
    std::vector<std::vector<std::vector<double>>> Vsigma_VS1 = metoda_MonteCarlo(nx, ny, delta,
                                                                                epsilon, Nchains,
                                                                                nlength, block_B,
                                                                                rho_max, sigma_rho,
                                                                                VL, VB, VT);
    save_to_file("V_1.dat", Vsigma_VS1[0]);
    save_to_file("sigma_V_1.dat", Vsigma_VS1[1]);
    save_to_file("S_1.dat", Vsigma_VS1[2]);
    std::cout << "1. zrobione" << std::endl;

    Nchains = 100;
    nlength = 100;
    block_B = 1;
    std::vector<std::vector<std::vector<double>>> Vsigma_VS2 = metoda_MonteCarlo(nx, ny, delta,
                                                                                epsilon, Nchains,
                                                                                nlength, block_B,
                                                                                rho_max, sigma_rho,
                                                                                VL, VB, VT);
    save_to_file("V_2.dat", Vsigma_VS2[0]);
    save_to_file("sigma_V_2.dat", Vsigma_VS2[1]);
    save_to_file("S_2.dat", Vsigma_VS2[2]);
    std::cout << "2. zrobione" << std::endl;

    Nchains = 300;
    nlength = 300;
    block_B = 1;                                               
    std::vector<std::vector<std::vector<double>>> Vsigma_VS3 = metoda_MonteCarlo(nx, ny, delta,
                                                                                epsilon, Nchains,
                                                                                nlength, block_B,
                                                                                rho_max, sigma_rho,
                                                                                VL, VB, VT);
                                                                                
    save_to_file("V_3.dat", Vsigma_VS3[0]);
    save_to_file("sigma_V_3.dat", Vsigma_VS3[1]);
    save_to_file("S_3.dat", Vsigma_VS3[2]);
    std::cout << "3. zrobione" << std::endl;
}
