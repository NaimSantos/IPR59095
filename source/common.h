#include <iostream>
#include <iomanip>  //std::setw
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iterator> //std::istream_iterator


std::vector<double> tdma2(const std::vector<std::vector<double>>& Mat, const std::vector<double>& X);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> void print_array_1D(const std::vector<T> A);


struct Fluid{
	double rho;      // massa específica
	double mu;       // viscosidade
	double c;        // compressibilidade
};
struct Rock{
	double phi_ref;  // porosidade de referencia
	double k_x;      // permeabilidade
	double k_y;      // permeabilidade
	double k_z;      // permeabilidade
	double c_phi;
};

