#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(void){

    srand(time(NULL));
    FILE* plik;

    plik = fopen("lab1.dat", "w");
    if (plik == NULL) {
        perror("Kaszana");
        return 1;
    }

    double p = 0.0;
    int power_of_ten = 7;

    int* array = malloc(pow(10, power_of_ten)*sizeof(int));
    double random_number;

    double sum_X1 = 0;
    double sum_X2 = 0;

    double X, X2, err_X, var_num, var_teo, err_var;

    double* vals = calloc(3, sizeof(double));
    vals[0] = 0.1; vals[1] = 0.5; vals[2] = 0.9;
    int k = 2;
    for(int o = 0; o < 3; o++){
        p = vals[o];
        k = 2;
        for(int i = 1; i <= pow(10, power_of_ten); i++){
            random_number = (double)rand()/(double)(RAND_MAX);
            if(random_number<p){
                array[i-1] = 1;
            } else {
                array[i-1] = 0;
            }
            sum_X1 += array[i-1];
            sum_X2 += array[i-1]*array[i-1]; // = array[i-1], ale formalnie ok

            if(i == pow(10, k)){

                X = sum_X1/(double)i;
                X2 = sum_X2/(double)i;
                err_X = fabs(X-p)/p;
                var_num = (X2-X*X)/(double)i;
                var_teo = (p - p*p)/(double)i;
                err_var = fabs((var_num-var_teo)/var_teo);

                fprintf(plik, "%.17lf", X); fprintf(plik, " ");
                fprintf(plik, "%.17lf", X2); fprintf(plik, " ");
                fprintf(plik, "%.17lf", err_X); fprintf(plik, " ");
                fprintf(plik, "%.17lf", var_num); fprintf(plik, " ");
                fprintf(plik, "%.17lf", var_teo); fprintf(plik, " ");
                fprintf(plik, "%.17lf", err_var); fprintf(plik, " ");
                fprintf(plik, "%.17lf", pow(10, k)); fprintf(plik, "\n");
                
                k++;
            }
        }
        fprintf(plik, "\n");
        sum_X1 = 0;
        sum_X2 = 0;
    }
    fflush(plik);
    fclose(plik);
    free(array);
    free(vals);

return 0;
}
