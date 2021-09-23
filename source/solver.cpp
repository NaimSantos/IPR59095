#include <iostream>
#include <iomanip>  //std::setw
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iterator> //std::istream_iterator

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

void gauss_siedel_solver(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double>& X);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> void print_array_1D(const std::vector<T> A);


int main(int argc, char* argv[]){

	constexpr double x_len {1000.0};
	constexpr double y_len {100.0};
	constexpr double z_len {10.0};
	constexpr double rho {1000.0};       // massa específica
	constexpr double phi {1.0};          // porosidade
	constexpr double phi_ref {1.0};      // porosidade
	constexpr double mu {1.0};           // viscosidade
	constexpr double k_x {1.0};          // permeabilidade
	constexpr double k_y {1.0};          // permeabilidade
	constexpr double k_z {1.0};          // permeabilidade
	constexpr double c_phi {1.0};        // compressibilidade
	constexpr double B_ref {1.0};        // B de referência
	constexpr double P_ref {3000.0};     // pressão de referência
	constexpr double FVF {0.5};          // fator volume formação

	// ---- Exemplo com o TDMA ---- //
	std::vector<double> a1 = {0, 1, 2, 3};      //diagonal inferior
	std::vector<double> b1 = {2, 3, 5, 8};      //diagonal principal
	std::vector<double> c1 = {1, 2, 1, 0};      //diagonal superior
	std::vector<double> d1 = {7, 19, 31, 52};   // b

	tdma_solver(a1, b1, c1, d1);
	std::cout << "Solucao (TDMA): ";
	print_array_1D<double>(d1);
}

double evaluate_B(const double p){
	return B_ref / (1 + c*(p - P_ref));
}

void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d){
	auto n = static_cast<int>(d.size()-1);

	c[0] /= b[0];
	d[0] /= b[0];

	for (int i = 1; i < n; i++){
		c[i] = (c[i] ) / (b[i] - a[i]*c[i-1]);
		d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
	}

	d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

	for (int i = n; i-- > 0;){
		d[i] = d[i] - (c[i]*d[i+1]);
	}
}

template <typename T>
void print_array_2D(const std::vector<std::vector<T>> A){
	auto nrow = A.size();
	auto ncol = A[0].size();
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << A[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

template <typename T>
void print_array_1D(const std::vector<T> A){
	auto nrow = A.size();
	for(int i = 0; i < nrow; i++){
		std::cout << A[i] << ' ';
	}
	std::cout << '\n';
}
