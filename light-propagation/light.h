#ifndef LIGHT_H
#define LIGHT_H

#include<cmath>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <random>
#include <cstdio>

using vector = std::vector<double>;
using matrix = std::vector<std::vector<double>>;

class DYFUZJA_FOTONOW_2D {
    private: 
        class WIAZKA {
            public:
                double w;
                double x, y;
                double x_nowe, y_nowe;
                double rx, ry;
                int warstwa;
                bool zywa = true;
                int dlugosc;
                vector sciezka;
        };
    
    public:     
        int liczba_warstw;
        int nx, ny;
        double xmax, ymax, dx, dy;
        double x_zrodla, dx_zrodla;
        double x_detekcji, dx_detekcji;
        double p_min;
        double w_min;
        double rx0, ry0;
        
        int zapisz_wszystkie_sciezki;
        int zapisz_sciezki_detekcji;
        
        WIAZKA wiazka;
        double abs_lustrzana;
        
        matrix absorpcja;
        vector odbicie;
        vector transmisja;
        matrix dane_warstw;
        
        DYFUZJA_FOTONOW_2D();
        void inicjalizuj();
        double losuj_jednorodnie();
        void pojedyncza_sciezka();
        void ruletka();
        void oblicz_nowe_polozenie();
        void rozprosz_w_warstwie(); 
        void rozprosz_na_brzegach_gora_dol();
        double znak(double);
        void przeciecie_odcinkow(double, double, double, double, double, double, double, double, double&, double&, int&);
        void zapisz_sciezki_do_pliku();
};

DYFUZJA_FOTONOW_2D::DYFUZJA_FOTONOW_2D() {
    zapisz_wszystkie_sciezki = 0;
    zapisz_sciezki_detekcji = 0;
    liczba_warstw = 0;
    nx = 0;
    ny = 0;
    xmax = 0.0;
    ymax = 0.0;
    dx = 0.0;
    dy = 0.0;
    p_min = 0.1;
    w_min = 1.0E-4;
    abs_lustrzana = 0.0;
    
    rx0 = 0.0;
    ry0 = 1.0;
    
    int M = 20;
    dane_warstw.resize(M, vector(M, 0.0));     
    for(int i = 0; i < M; i++){
        dane_warstw[i][3] = 0.0; 
        dane_warstw[i][4] = 1.0; 
    }
}       

void DYFUZJA_FOTONOW_2D::inicjalizuj() {
    absorpcja.resize(nx + 1, vector(ny + 1, 0.0));
    odbicie.resize(nx + 1, 0.0);
    transmisja.resize(nx + 1, 0.0);
    
    dx = xmax/nx;
    ymax = 0.0;
    for(int i = 1; i <= liczba_warstw; i++) ymax += dane_warstw[i][2];
    dy = ymax/ny;
    
    for(int i = 1; i <= liczba_warstw + 1; i++){
        dane_warstw[i][5] = dane_warstw[i-1][6];  
        dane_warstw[i][6] += dane_warstw[i][5] + dane_warstw[i][2];   
    }
    
    if((x_zrodla - dx_zrodla/2.0) < 0 || (x_zrodla + dx_zrodla/2.0) > xmax){
        printf("zrodlo poza obszarem - STOP\n");
        exit(0);
    }
    
    FILE *fp;
    fp = fopen("sciezki_detekcji_zrodla.dat", "w");
    fclose(fp);
    fp = fopen("wszystkie_sciezki.dat", "w");
    fclose(fp);
}

double DYFUZJA_FOTONOW_2D::losuj_jednorodnie() {
    return (double)rand()/RAND_MAX;
}

void DYFUZJA_FOTONOW_2D::pojedyncza_sciezka() {
    if(dx_zrodla >= 0.0  && x_zrodla >= 0.0 && x_zrodla <= xmax){
        wiazka.x = x_zrodla + dx_zrodla/2.0*(2.0*losuj_jednorodnie() - 1.0);
        wiazka.y = 1.0E-10;
        wiazka.rx = rx0;
        wiazka.ry = ry0;
        wiazka.w = 1.0;
        wiazka.zywa = true;
        wiazka.warstwa = 1;
        
        wiazka.sciezka.clear();
        wiazka.sciezka.push_back(wiazka.x);
        wiazka.sciezka.push_back(wiazka.y);
        wiazka.sciezka.push_back(wiazka.w);
    }

    int l = 0;
    while(wiazka.zywa == true){
        oblicz_nowe_polozenie(); 
        rozprosz_w_warstwie(); 
        rozprosz_na_brzegach_gora_dol();
        ruletka();

        wiazka.sciezka.push_back(wiazka.x);
        wiazka.sciezka.push_back(wiazka.y);
        wiazka.sciezka.push_back(wiazka.w);
        l++;
    }
    
    zapisz_sciezki_do_pliku();
}

void DYFUZJA_FOTONOW_2D::zapisz_sciezki_do_pliku() {
    if(zapisz_sciezki_detekcji == 1 && std::abs(wiazka.x - x_detekcji) <= dx_detekcji/2.0){
        FILE *fp;
        fp = fopen("sciezki_detekcji_zrodla.dat", "a");
        fprintf(fp, "\n");
        
        for(int i = 0; i < (int)wiazka.sciezka.size(); i += 3){
            double x = wiazka.sciezka[i];
            double y = wiazka.sciezka[i+1];
            double w = wiazka.sciezka[i+2];
            fprintf(fp, "%15.5E  %15.5E   %15.5E \n", x, y, w);
        }
        fclose(fp);
    }
    
    if(zapisz_wszystkie_sciezki == 1){
        FILE *fp;
        fp = fopen("wszystkie_sciezki.dat", "a");
        fprintf(fp, "\n");
        
        for(int i = 0; i < (int)wiazka.sciezka.size(); i += 3){
            double x = wiazka.sciezka[i];
            double y = wiazka.sciezka[i+1];
            double w = wiazka.sciezka[i+2];
            fprintf(fp, "%15.5E  %15.5E   %15.5E \n", x, y, w);
        }
        fclose(fp);
    }
}

double DYFUZJA_FOTONOW_2D::znak(double x) {
    if(x < 0.0) return -1.0;
    else return 1.0;
}

void DYFUZJA_FOTONOW_2D::przeciecie_odcinkow(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double& x_przeciecia, double& y_przeciecia, int& czy_przecina) {
    double a, b, c, d, b1, b2, alfa, beta, det;
    a = x2 - x1;
    b = -(x4 - x3);
    c = y2 - y1;
    d = -(y4 - y3);
    b1 = x3 - x1;
    b2 = y3 - y1;
    
    det = a*d - b*c;
    czy_przecina = 0;
    if(std::abs(det) > 1.0E-50){
        alfa = (d*b1 - b*b2)/det;
        beta = (-c*b1 + a*b2)/det;
        if(alfa >= 0.0 && alfa <= 1.0 && beta >= 0.0 && beta <= 1.0){
            x_przeciecia = x1 + alfa*(x2 - x1);
            y_przeciecia = y1 + alfa*(y2 - y1);
            czy_przecina = 1;
        }
    }
}

void DYFUZJA_FOTONOW_2D::rozprosz_na_brzegach_gora_dol() {
    if(wiazka.dlugosc == 1){
        double y_dol = dane_warstw[wiazka.warstwa][5];
        double y_gora = dane_warstw[wiazka.warstwa][6];
        double x_przeciecia, y_przeciecia;
        int i, j;
        
        double x1, y1, x2, y2, x3, y3, x4, y4;
        int czy_przecina, ktory;
        
        x1 = wiazka.x;
        y1 = wiazka.y;
        x2 = wiazka.x_nowe;
        y2 = wiazka.y_nowe;      
        ktory = 0; 
        
        if(ktory == 0){
            x3 = 0.0;
            y3 = y_dol;
            x4 = 0.0;
            y4 = y_gora;
            przeciecie_odcinkow(x1, y1, x2, y2, x3, y3, x4, y4, x_przeciecia, y_przeciecia, czy_przecina);   
            if(czy_przecina == 1) ktory = 1;
        }
        if(ktory == 0){
            x3 = xmax;
            y3 = y_dol;
            x4 = xmax;
            y4 = y_gora;
            przeciecie_odcinkow(x1, y1, x2, y2, x3, y3, x4, y4, x_przeciecia, y_przeciecia, czy_przecina);              
            if(czy_przecina == 1) ktory = 2;
        }
        if(ktory == 0){
            x3 = 0.0;
            y3 = y_gora;
            x4 = xmax;
            y4 = y_gora;
            przeciecie_odcinkow(x1, y1, x2, y2, x3, y3, x4, y4, x_przeciecia, y_przeciecia, czy_przecina);              
            if(czy_przecina == 1) ktory = 3;
        }
        if(ktory == 0){
            x3 = 0.0;
            y3 = y_dol;
            x4 = xmax;
            y4 = y_dol;
            przeciecie_odcinkow(x1, y1, x2, y2, x3, y3, x4, y4, x_przeciecia, y_przeciecia, czy_przecina);              
            if(czy_przecina == 1) ktory = 4;
        }
        
        i = std::round(x_przeciecia/dx);
        j = std::round(y_przeciecia/dy);
        
        if(ktory == 1 || ktory == 2){
            absorpcja[i][j] += wiazka.w;   
            wiazka.x = x_przeciecia;
            wiazka.y = y_przeciecia;
            wiazka.dlugosc = 0;
            wiazka.zywa = false;
            return;
        }
                
        if(ktory == 3 || ktory == 4){
            double no, nn, sin_alfa_kryt, sin_alfa;
            
            no = dane_warstw[wiazka.warstwa][4];
            if(ktory == 3) nn = dane_warstw[wiazka.warstwa + 1][4];
            if(ktory == 4) nn = dane_warstw[wiazka.warstwa - 1][4];
            
            sin_alfa_kryt = nn/no;
            sin_alfa = std::abs(wiazka.rx);
            
            if(sin_alfa > sin_alfa_kryt){ 
                wiazka.ry = -wiazka.ry;
                wiazka.x = x_przeciecia;
                wiazka.y = y_przeciecia + wiazka.ry*1.0E-10; 
                wiazka.dlugosc = 0;
                return;
            } else {   
                double alfa_o = std::abs(std::asin(wiazka.rx));
                double alfa_n = std::abs(std::asin(no/nn*wiazka.rx));
                double odbicie_wsp = (std::pow(std::sin(alfa_o - alfa_n), 2)/std::pow(std::sin(alfa_o + alfa_n), 2) + std::pow(std::tan(alfa_o - alfa_n), 2)/std::pow(std::tan(alfa_o + alfa_n), 2))/2.0;
                double u1 = losuj_jednorodnie();
                if(u1 < odbicie_wsp){ 
                    wiazka.ry = -wiazka.ry;
                    wiazka.x = x_przeciecia;
                    wiazka.y = y_przeciecia + wiazka.ry*1.0E-10; 
                    wiazka.dlugosc = 0;
                    return;
                } else {  
                    double rx = wiazka.rx; 
                    double ry = wiazka.ry; 
                    wiazka.rx = no/nn*rx; 
                    wiazka.ry = znak(ry)*std::sqrt(1.0 - std::pow(no/nn*rx, 2));  
                    wiazka.x = x_przeciecia;
                    wiazka.y = y_przeciecia + wiazka.ry*1.0E-10; 
                    wiazka.dlugosc = 0;
                    
                    if(ktory == 3) wiazka.warstwa++;
                    else if(ktory == 4) wiazka.warstwa--;
                    
                    if(wiazka.warstwa == (liczba_warstw + 1)){
                        transmisja[i] += wiazka.w;
                        wiazka.w = 0.0;
                        wiazka.zywa = false;
                    }
                    else if(wiazka.warstwa == 0){
                        odbicie[i] += wiazka.w;
                        wiazka.w = 0.0;
                        wiazka.zywa = false;
                    }
                    return;
                }
            }
        }           
    }
}

void DYFUZJA_FOTONOW_2D::rozprosz_w_warstwie() {
    double y1 = dane_warstw[wiazka.warstwa][5];
    double y2 = dane_warstw[wiazka.warstwa][6];
    
    if(wiazka.x_nowe > 0.0 && wiazka.x_nowe < xmax && wiazka.y_nowe > y1 && wiazka.y_nowe < y2 && wiazka.dlugosc == 1){
        
        double g, u1, u2, cos_teta, znak_wart, sin_teta, rx, ry;
        u1 = losuj_jednorodnie();
        g = dane_warstw[wiazka.warstwa][3];
        if(g > 1.0E-3){
            cos_teta = (1.0 + g*g - std::pow((1.0 - g*g)/(1.0 - g + 2.0*g*u1), 2))/2.0/g;
        } else {
            cos_teta = 2.0*u1 - 1.0;
        }
        u2 = losuj_jednorodnie();
        if(u2 <= 0.5){
            znak_wart = -1.0;
        } else {
            znak_wart = 1.0;
        }
        
        sin_teta = znak_wart*std::sqrt(1.0 - cos_teta*cos_teta);
        rx = wiazka.rx;
        ry = wiazka.ry;
        wiazka.rx = cos_teta*rx - sin_teta*ry;
        wiazka.ry = sin_teta*rx + cos_teta*ry;

        double dw = 0.0;
        double wsp_abs = dane_warstw[wiazka.warstwa][0];
        double wsp_rozpraszania = dane_warstw[wiazka.warstwa][1];
        dw = wsp_abs/(wsp_abs + wsp_rozpraszania)*wiazka.w;
        wiazka.w = wiazka.w - dw;
        int i, j;
        i = std::round(wiazka.x_nowe/dx);
        j = std::round(wiazka.y_nowe/dy);
        absorpcja[i][j] += dw;
        wiazka.x = wiazka.x_nowe;
        wiazka.y = wiazka.y_nowe;
        wiazka.dlugosc = 0;
    }
}

void DYFUZJA_FOTONOW_2D::oblicz_nowe_polozenie() {
    double wsp_abs = dane_warstw[wiazka.warstwa][0];
    double wsp_rozpraszania = dane_warstw[wiazka.warstwa][1];
    double wsp = wsp_abs + wsp_rozpraszania;
    double s = -std::log(losuj_jednorodnie())/wsp;
    wiazka.x_nowe = wiazka.x + wiazka.rx*s;    
    wiazka.y_nowe = wiazka.y + wiazka.ry*s;
    wiazka.dlugosc = 1;
}

void DYFUZJA_FOTONOW_2D::ruletka() {
    if(wiazka.w < w_min && wiazka.zywa == true){
        double p = losuj_jednorodnie();
        if(p <= p_min){
            wiazka.w = wiazka.w/p;
        } else {
            wiazka.w = 0.0;
            wiazka.zywa = false;
        }
    }
}