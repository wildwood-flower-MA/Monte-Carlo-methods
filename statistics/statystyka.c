#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "statystyka.h"

double f(double x){
    return 0.8*(1.0 + x - x*x*x);
}

double F(double x){
    return 0.8*x + (2.0*x*x - x*x*x*x)/5.0;
}

double rand_unif(void){
    return (double)rand()/(double)(RAND_MAX);
}

// ROZKLAD ZLOZONY
double rand_zlozony(double g1, double g2){

    double random_number_1 = rand_unif();
    double random_number_2 = rand_unif();

    if(random_number_1 <= g1){
        return random_number_2;
    } else {
        return sqrt(1.0-sqrt(1-random_number_2));
    }
}

// LANCUCH MARKOWA
double p_acc(double Xi, double Xi1, double (*f)(double)){    
    // f jest gestoscia p-stwa zmiennej X

    double some_value = f(Xi1)/f(Xi);

    if(some_value >= 1.0){
        return (double)1.0;
    } else {
        return some_value;
    }
}

double rand_Markow(double x0, double delta){

    double x1 = x0 + (2*rand_unif() - 1)*delta;
    if(x1 <= 1 && x1 >= 0 && rand_unif() <= p_acc(x0, x1, f)){
        return x1;
    } else {
        return x0;
    }
}

// METODA ELIMINACJI
double rekurencja_przy_eliminacji(double* random_number, double* generator, int* bezpiecznik){

    if(*bezpiecznik > 100000){
        return (double)0;
    }
    if(*generator > f(*random_number)){

        *random_number = rand_unif();
        *generator = 1.15*rand_unif();
        (*bezpiecznik)++;

        return rekurencja_przy_eliminacji(random_number, generator, bezpiecznik);
    } else {
        *bezpiecznik = 0;
        return *random_number;
    }
}

double rand_eliminacja(void){
    int bezpiecznik = 0;
    double random_number = rand_unif();
    double generator = 1.15*rand_unif();

    return rekurencja_przy_eliminacji(&random_number, &generator, &bezpiecznik);
}

// HISTOGRAM
double* histogram_nowy(double* array, size_t size, int N){

    double* bins = (double*)calloc(N, sizeof(double));
    int j;
    for(int n = 0; n < size; n++){
        j = (int)floor(N*array[n]);
        *(bins+j) += (double)N/(double)size;
    }
    return bins;
}

void histogramuj_nowy(int N, double* random, size_t size, FILE* plik){
    // N - liczba binow histogramu brana pod uwagę
    // random - tablica z wylosowanymi liczbami
    // size - dlugosc tej tablicy
    // plik - plik do zapisu

    printf("Przedział: %f --> %f\n", 0.0, 1.0);

    double* hist = histogram_nowy(random, size, N);
    if (hist== NULL) {
        return;
    }
    for (int n = 0; n < N; n++) {
        fprintf(plik, "%f ", hist[n]);
    }
    fprintf(plik, "%f %f", 0.0, 1.0);
    for (int n = 0; n < N; n++) {
        printf("bin nr %f.: %f\n", n+1, hist[n]);
    }
    free(hist);
}

// CHI^2
double Chi2_nowy(double* random, size_t size, int N, double (*F)(double)){
    // random - tablica z wylosowanymi liczbami
    // size - dlugosc tej tablicy
    // N - liczba binow histogramu brana pod uwagę
    // F - dystrybuanta

    double chi2 = 0.0;
    double* hist = histogram_nowy(random, size, N);

    double fn = 0;
    for(int n = 0; n < N; n++){
        fn = -F(0.1*(double)n)+F(0.1*((double)n+1.0));
        chi2 += pow(0.1*hist[n]*size-fn*size, 2)/fn/size;
    }

    free(hist);
    //free(bin_separating_value);
    //free(F_values);

    return chi2;
}

double statystyka_C(double* random1, double* random2, size_t size){

    // POPRAWKA 1: Inicjalizacja zerami
    double s1 = 0.0, s2 = 0.0, mu1 = 0.0, mu2 = 0.0;
    
    for(int n = 0; n < size; n++){
        mu1 += random1[n];
        mu2 += random2[n];
    }
    mu1 /= size;
    mu2 /= size;

    for(int n = 0; n < size; n++){
        s1 += pow(random1[n]-mu1, 2);
        s2 += pow(random2[n]-mu2, 2);
    }
    s1 /= (size-1);
    s2 /= (size-1);

    return (mu1-mu2)/sqrt(s1/size+s2/size);
}

/*double Chi2(double* random, size_t size, int N, double (*F)(double)){
    // random - tablica z wylosowanymi liczbami
    // size - dlugosc tej tablicy
    // N - liczba binow histogramu brana pod uwagę
    // F - dystrybuanta

    double chi2 = 0.0;
    int* hist = histogram(random, size, N);

    double* bin_separating_value = _bin_separating_value(random, size, N);
    double* F_values = (double*)malloc((N+1)*sizeof(double));
    for(int n = 0; n <= N; n++){
        *(F_values + n) = F(*(bin_separating_value + n));
    }

    double fn = 0;
    for(int n = 0; n < N; n++){
        fn = (F_values[n+1]-F_values[n])*size;
        chi2 += pow(hist[n]-fn, 2)/fn;
    }

    free(hist);
    free(bin_separating_value);
    free(F_values);

    return chi2;
}*/

/*int* histogram(double* array, size_t size, int N){
    // array - tablica
    // size - jej rozmiar
    // N - liczba binow histogramu

    double* bin_separating_value = _bin_separating_value(array, size, N);

    int* bins = (int*)calloc(N, sizeof(int));
    for(int n = 0; n < size; n++)
        for(int m = 0; m < N; m++){
            int b1 = bin_separating_value[m]<array[n];
            int b2 = bin_separating_value[m+1]>=array[n];
            if(b1 && b2){
                bins[m]++;
            }
        }
    free(bin_separating_value);
    return bins;
}

void histogramuj(int N, double* random, size_t size, FILE* plik){
    // N - liczba binow histogramu brana pod uwagę
    // random - tablica z wylosowanymi liczbami
    // size - dlugosc tej tablicy
    // plik - plik do zapisu

    double* minmax = min_max(random, size);
    printf("Przedział: %f --> %f\n", minmax[0], minmax[1]);
    free(minmax);

    int* hist = histogram(random, size, N);
    if (hist== NULL) {
        return;
    }
    for (int n = 0; n < N; n++) {
        fprintf(plik, "%d ", hist[n]);
    }
    fprintf(plik, "%f %f", minmax[0], minmax[1]);
    for (int n = 0; n < N; n++) {
        printf("bin nr %d.: %d\n", n+1, hist[n]);
    }
    free(hist);
}*/

/*double* min_max(double* array, size_t size){

    double* minmax = calloc(2,sizeof(double));

    double max_value = *array;
    double min_value = *array;
    double current_value = 0.0;

    for(int n = 0; n < size; n++){
        current_value = *(array+n);
        if(max_value <= current_value){
            max_value = current_value;
        }
        if(min_value > current_value){
            min_value = current_value;
        }
    }
    minmax[0] = min_value;
    minmax[1] = max_value;
    return minmax;
}

double* _bin_separating_value(double* array, size_t size, int N){

    double* minmax = min_max(array, size);
    double range = minmax[1] - minmax[0];
    double* bin_separating_value = malloc((N+1) * sizeof(double));

    for(int n = 0; n <= N; n++){
        bin_separating_value[n] = minmax[0] + (double)n*range/(double)N;
    }
    free(minmax);
    return bin_separating_value;
}*/

int main(void){

    srand(time(NULL));
    FILE* plik;

    plik = fopen("lab2.dat", "w");
    if (plik == NULL) {
        perror("Kaszana");
        return 1;
    }

    FILE* plik2;
    plik2 = fopen("lab2_chain.dat", "w");
    if (plik2 == NULL) {
        perror("Kaszana");
        return 1;
    }

    int n_losowanch = pow(10,6);

    double* random_zlozony = malloc(n_losowanch*sizeof(double));
    double* random_Markow_05 = malloc(n_losowanch*sizeof(double));
    double* random_Markow_005 = malloc(n_losowanch*sizeof(double));
    double* random_eliminacja = malloc(n_losowanch*sizeof(double));

    // LOSOWANIE
    for(int i = 0; i<n_losowanch; i++){
        random_zlozony[i] = rand_zlozony(0.8, 0.2);
        random_eliminacja[i] = rand_eliminacja();
    }

    // Poczatek lancuchow Markowa (warunek poczatkowy z metody eliminacji)
    *random_Markow_05 = rand_eliminacja();
    *random_Markow_005 = rand_eliminacja();
    
    // Generowanie lancuchow Markowa
    for(int i = 1; i<n_losowanch; i++){
        random_Markow_05[i] = rand_Markow(random_Markow_05[i-1], 0.5);
        random_Markow_005[i] = rand_Markow(random_Markow_005[i-1], 0.05);
    }

    //LICZBA BINOW
    int N = 10;

    // HISTOGRAM ZLOZONY
    printf("\nHistogram zlozony, biny:\n");
    histogramuj_nowy(N, random_zlozony, n_losowanch, plik);
    fprintf(plik, "\n");

    // HISTOGRAMY MARKOW 
    printf("\nHistogram Markow 0.5, biny:\n");
    histogramuj_nowy(N, random_Markow_05, n_losowanch, plik);
    fprintf(plik, "\n");

    printf("\nHistogram Markow 0.05, biny:\n");
    histogramuj_nowy(N, random_Markow_005, n_losowanch, plik);
    fprintf(plik, "\n");

    // HISTOGRAM ELIMINACJA
    printf("\nHistogram eliminacja, biny:\n");
    histogramuj_nowy(N, random_eliminacja, n_losowanch, plik);
    fprintf(plik, "\n");

    // CHI^2
    double chi2_zlozony = Chi2_nowy(random_zlozony, n_losowanch, N, F);
    double chi2_Markow_05 = Chi2_nowy(random_Markow_05, n_losowanch, N, F);
    double chi2_Markow_005 = Chi2_nowy(random_Markow_005, n_losowanch, N, F);
    double chi2_eliminacja = Chi2_nowy(random_eliminacja, n_losowanch, N, F);

    printf("\nzlozony, chi2 = %f", chi2_zlozony);
    printf("\nM05, chi2 = %f", chi2_Markow_05);
    printf("\nM005, chi2 = %f", chi2_Markow_005);
    printf("\nelimin, chi2 =%f", chi2_eliminacja);

    // STATYSTYKA C
    double statystyka = statystyka_C(random_Markow_05, random_Markow_005, n_losowanch);
    printf("\nStatystyka C: %f\n", statystyka);

    free(random_zlozony);
    free(random_Markow_05);
    free(random_Markow_005);
    free(random_eliminacja);
    fflush(plik);
    fclose(plik);

    n_losowanch = 1000;
    double* random_Markow_05_test = (double*)malloc(n_losowanch * sizeof(double));
    double* random_Markow_005_test = (double*)malloc(n_losowanch * sizeof(double));
    
    *random_Markow_05_test = rand_eliminacja();
    *random_Markow_005_test = rand_eliminacja();
    
    for(int i = 1; i<n_losowanch; i++){
        random_Markow_05_test[i] = rand_Markow(random_Markow_05_test[i-1], 0.5);
        random_Markow_005_test[i] = rand_Markow(random_Markow_005_test[i-1], 0.05);
    }

    for(int i = 0; i < n_losowanch; i++){
        fprintf(plik2, "%f ", random_Markow_05_test[i]);
        fprintf(plik2, "%f\n", random_Markow_005_test[i]);
    }
    fflush(plik2);
    fclose(plik2);
    free(random_Markow_05_test);
    free(random_Markow_005_test);

    return 0;
}