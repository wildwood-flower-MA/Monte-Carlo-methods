#ifndef ELEMENTS_H
#define ELEMENTS_H

// 1.1 rozkład sferycznie konturowany- normalny
inline double random_uniform_01();
std::vector<double> random_2D_BoxMuller();

// 1.2 Rozkład jednorodny w kole K^2(0,1)
inline std::vector<double> normalise(const std::vector<double>&);
std::vector<double> random_circle(const std::vector<double>&);
std::vector<double> random_circle();

// 1.3 Transformacja afiniczna: koło → elipsa
// 1.3.1 Wybór osi i skalowanie
std::vector<double> operator*(const std::vector<std::vector<double>>&, const std::vector<double>&);
std::vector<double> operator*(double c, const std::vector<double>&);
std::vector<double> operator+(const std::vector<double>&, const std::vector<double>&);
std::vector<std::vector<double>> create_R(double);

// 1.4 Wyznaczanie macierzy kowariancji

std::pair<double,double> avg_x_y(const std::vector<std::vector<double>>&);
std::pair<double,double> avg_x2_y2(const std::vector<std::vector<double>>&);
double avg_xy(const std::vector<std::vector<double>>&);
std::vector<std::vector<double>> cov_matrix(const std::vector<std::vector<double>>&);
double corr(const std::vector<std::vector<double>>& pairs);
void transform(const std::vector<std::vector<double>>,std::vector<std::vector<double>>&,double, double,double,double,double, double);

// 1.5 Transformacja afiniczna a macierz kowariancji dla rozkładu gaussowskiego

// OGÓLNE

void zapisz(std::ostream&, const std::vector<std::vector<double>>&);
void rysuj(const std::vector<std::vector<double>>&, std::string);

#endif 