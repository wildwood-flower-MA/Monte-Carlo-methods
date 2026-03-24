#pragma once
#include <vector>
#include <random>
#include <string>

class GasDynamics {
private:
    class CZASTKA {
    public:
        double x, y, vx, vy, v, mc, rc;
        double x0, y0, vx0, vy0;
        double vxt = 0.0;
        double vyt = 0.0;
        double cv;
        int ic, ix, iy, irun;
        int indx = -1;
        int ncol = 0;
        int nbound_col = 0;
        double path = 0.0;
    };

public:
    double kb = 1.38E-23;
    double xmin, xmax, ymin, ymax;
    int nx, ny, nxy, n_mieszanina;
    double temp;
    double tempi[20];
    std::vector<std::vector<double>> temp_komorki;
    std::vector<std::vector<double>> gestosc_komorki;
    std::vector<std::vector<double>> vx_komorki;
    std::vector<std::vector<double>> vy_komorki;
    std::vector<std::vector<double>> cisnienie_komorki;
    std::vector<std::vector<std::vector<double>>> tensor_cisnienia;
    std::vector<double> predkosc_x;
    std::vector<double> predkosc_y;
    std::vector<double> cisnienie_x;
    std::vector<double> cisnienie_y;
    std::vector<int> indx;
    std::vector<std::vector<int>> indx0;
    std::vector<CZASTKA> czastki;
    int nc[5];
    int ntot;
    double mc[5], rc[5];
    int init_dist;
    double dt, suma_czasu;
    double delta_x, delta_y;
    int wezly_zewn = 4;
    int wezly;
    double krawedz_zewn[100][2];
    double krawedz[100][2];
    int icol = 1;
    int liczba_watkow = 1;
    double droga_swobodna_num, droga_swobodna_teo, droga_swobodna_enskog, gamma_col;
    double wsp_cv = 0.0;
    double wsp_dyfuzji = 0.0;
    unsigned ziarno = 12345;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> losuj_U = std::uniform_real_distribution<double>(0.0, 1.0);
    std::normal_distribution<double> losuj_N = std::normal_distribution<double>(0.0, 1.0);

    void czytaj(const char*);
    void inicjalizuj();
    void hist_predkosci_wszystkie(const char*, double, int);
    void zapisz_nptv(const char*, int);
    void ewolucja(double, int);
    void zapisz_polozenie_predkosc(const char*);
    void zapisz_granice_komorek(const char*);
    double temperatura_efektywna();
    void oblicz_temperature_komorki();
    void oblicz_gestosc_komorki();
    void oblicz_predkosc_komorki();
    void oblicz_srednia_droge_swobodna();
    void oblicz_autokorelacje();
    void komorka_ix_iy(double, double, int*, int*);
    void oblicz_cisnienie();
    int globalny_indeks_komorki(int, int);
    void krok();
    void sortuj();
    void rozklad_poczatkowy();
    int przeciecie(double, double, double, double, double, double, double, double, double*, double*, double*, double*, int*);
    int odbicie_od_bariery(CZASTKA&, double*);
    void swobodny_lot(CZASTKA&, double*);
    int zderzenia_czastek(int, int);
    void rozpraszanie_czastek_kule(CZASTKA&, CZASTKA&, int);
    void rozklad_maxwella(double*, double*, double, double);
    void oblicz_tensor_cisnienia_oddzialywanie(CZASTKA&, CZASTKA&, CZASTKA&, CZASTKA&);
    void oblicz_tensor_cisnienia_kinetyczny();
    void test();
};