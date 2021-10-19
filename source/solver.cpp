#include "common.h"

constexpr double k_x {10e-15};       // permeabilidade
constexpr double phi_ref {0.25};     // porosidade
constexpr double P_ini {4.5e7};     // porosidade
constexpr double Lx {500.0};         // dimensão em x
constexpr double Ly {40.0};          // dimensão em y
constexpr double Lz {10.0};          // dimensão em z
constexpr double B_ref {1.5};        // B de referência
constexpr double rho_o {820.0};      // massa específica do óleo
constexpr double mu_o {1.2e-3};         // viscosidade
constexpr double c_phi {6.0e-10};        // compressibilidade
constexpr double FVF {1.0};          // fator volume formação

int main(int argc, char* argv[]){
	std::vector<std::vector<double>> VX = {{2, 1, 0, 0}, {1, 3, 2, 0}, {0, 2, 5, 1}, {0, 0, 3, 8}};
	std::vector<double> X0 = {7, 19, 31, 52};
	print_array_2D(VX);
	auto res = tdma2(VX, X0);
	std::cout << "O resultado eh: "<< std::endl;
	print_array_1D(res);
}

std::vector<double> tdma2(const std::vector<std::vector<double>>& Mat, const std::vector<double>& X){
	auto n_row = Mat.size();
	auto n_col = Mat[0].size();
	if (n_row != n_col || n_row == 0)
		std::cout << "Matriz nao-quadrada detectada" << std::endl;

	std::vector<double> A (n_row, 0.0);
	std::vector<double> B (n_row, 0.0);
	std::vector<double> C (n_row, 0.0) ;
	std::vector<double> D = X;

	for (size_t i = 0; i < (n_row-1); i++){
		size_t j = i;
		B[i] = Mat[i][j];       // diagonal principal
		A[i+1] = Mat[i+1][j];   // diagonal inferior
		C[i] = Mat[i][j+1];     // diagonal superior
	}
	B[n_row-1] = Mat[n_row-1][n_col-1];

	C[0] /= B[0];
	D[0] /= B[0];
	int n = D.size()-1;
	for (int i = 1; i < n; i++){
		C[i] = (C[i] ) / (B[i] - A[i]*C[i-1]);
		D[i] = (D[i] - A[i]*D[i-1]) / (B[i] - A[i]*C[i-1]);
	}

	D[n] = (D[n] - A[n]*D[n-1]) / (B[n] - A[n]*C[n-1]);

	for (int i = n; i-- > 0;){
		D[i] = D[i] - (C[i]*D[i+1]);
	}

	return D;
}

double evaluate_B(const double p){
	return B_ref / (1 + c_phi*(p - P_ini));
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
