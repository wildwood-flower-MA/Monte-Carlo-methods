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
    
    GasDynamics ob;
    ob.read("i.dat"); // wczytujemy dane zpliku wejściowego
    ob.init(); // automatyczna inicjalizacja położeń i prędkości
    
    //ob.write_position_velocity("rv.dat"); //zapis ustawień początkowych
    
    ob.nthreads = 1; //obliczenia na jednym rdzeniu
    ob.icol = 1 ; //cząstki zderzają się
    
    ob.evolution(0.0, 0); // wykonujemy 20 tysięcy kroków (tmax - nieznany)
    
    //ob.hist_velocity_all("hist2.dat", 5.0 ,50); //zapis histogramu prędkosci do pliku
    //ob.write_position_velocity("rv.dat"); //zapis położeń i prędkości końcowych do pliku

    return 0;

}