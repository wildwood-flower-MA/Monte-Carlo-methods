#ifndef STATYSTYKA_H
#define STATYSTYKA_H

double f(double); // f. gestosci p-stwa
double F(double); // dystrybuanta
double rand_unif(void); // losowanie z plaskiego [0,1]

double rand_zlozony(double, double);

double p_acc(double, double, double (*)(double));
double rand_Markow(double, double);

double rekurencja_przy_eliminacji(double*, double*, int*);
double rand_eliminacja(void);

//double* min_max(double*, size_t);
//double* _bin_separating_value(double*, size_t, int);
//int* histogram(double*, size_t, int);
//void histogramuj(int, double*, size_t, FILE*);

double* histogram_nowy(double*, size_t, int);
void histogramuj_nowy(int, double*, size_t, FILE*);
double Chi2(double*, size_t, int, double (*)(double));

double Chi2_nowy(double*, size_t, int, double (*)(double));

double statystyka_C(double*, double*, size_t);


#endif