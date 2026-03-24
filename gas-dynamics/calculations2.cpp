#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <omp.h>
#include <time.h>
#include <vector>
#include <iomanip> 
#include <string>
#include <cstdlib>
#include <filesystem>
#include "gasdynamics.h"

using namespace std;

int main(){
    
    // Lewy
    GasDynamics ob_left;
    ob_left.read("i_left.dat");
    ob_left.init();
    ob_left.write_position_velocity("rv_left.dat");

    // Probierz
    GasDynamics ob_right;
    ob_right.read("i_right.dat");
    ob_right.init();
    ob_right.write_position_velocity("rv_right.dat");

    // Reprezentacja
    std::ofstream file("pos_vel_start.dat");
    std::ifstream rv_left("rv_left.dat");
    file << rv_left.rdbuf();
    rv_left.close();
    std::ifstream rv_right("rv_right.dat");
    file << rv_right.rdbuf();
    rv_right.close();
    file.close();

    GasDynamics ob;
    ob.read("i.dat");
    ob.init();
    ob.nthreads = 1; //obliczenia na jednym rdzeniu
    ob.icol = 1; //cząstki zderzają się

    ob.evolution(0.0, 2000);

    return 0;

}