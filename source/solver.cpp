#include <iostream>
#include <iomanip>  //std::setw
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator> //std::istream_iterator

double evaluate_B(const double p);
std::vector<double> tdma2(const std::vector<std::vector<double>>& Mat, const std::vector<double>& X);
template <typename T> void print_array_1D(const std::vector<T> A);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> std::vector<T> linspace(const double xi, const double xf, int Num);
void save_data(const std::vector<double>& X, const std::vector<double>& Y);

constexpr double k_x {10e-15};       // permeabilidade
constexpr double phi_ref {0.25};     // porosidade
constexpr double P_ini {4.5e7};      // porosidade
constexpr double Lx {10.0};          // dimensão em x
constexpr double Ly {40.0};          // dimensão em y
constexpr double Lz {10.0};          // dimensão em z
constexpr double B0 {1.5};           // B de referência
constexpr double mu {1.2e-3};        // viscosidade
constexpr double c_ref {6.0e-10};    // compressibilidade
constexpr double Vb {Lx*Ly*Lz};      // volume
constexpr double A_x {Ly*Lz};        // área
constexpr int N {11};                // número de células
constexpr double dx = Lx/N;
constexpr double ti {0.0};
constexpr double tf {2.0};
constexpr double dt {0.5};
constexpr auto nsteps = static_cast<int>((tf - ti)/dt);

int main(int argc, char* argv[]){
	std::vector<std::string> args(argv, argv + argc);
	std::cout << "Argumentos obtidos via linha de comando:" << std::endl;
	for (const auto& arg : args){
		std::cout << arg << std::endl;
	}
	auto X = linspace<double>(0.0, Lx, N);
	// Vetores:
	std::vector<double> P (N, P_ini);                                    // vetor com as pressões
	std::vector<double> P_old (N, P_ini);                                // vetor com as pressões anteriores
	std::vector<double> B (N, B0);                                       // vetor B;
	std::vector<std::vector<double>> T (N, std::vector<double>(N, 0.0)); // matriz de transmissibilidades
	
	
	// Variáveis utilizadas no processo iterativo:
	double Ei = 0.0;                        // Transmissibilidade à direita da célula i
	double Wi = 0.0;                        // Transmissibilidade à esquerda da célula i
	double Bi = 0.0;                        // B(p) na célula i
	double Bi_prev = 0.0;                   // B(p) na célula i - 1
	double Bi_next = 0.0;                   // B(p) na célula i + 1
	double Bh_prev = 0.0;                   // média harmônica i - 1/2
	double Bh_next = 0.0;                   // média harmônica i + 1/2
	
	// Iteraçao no tempo:
	for (size_t n = 0; n < nsteps; n++){

		// Iteração no domínio espacial:
		
		for (size_t i = 0; i < N; i++){

			Bi = evaluate_B(P[i]);
			// Contorno esquerdo:
			if (i == 0){                            
				T[i][i] = 1.0;                      // CORRIGIR
				P[i] = P[i+1];                      // CORRIGIR
			}
			// Contono direito:
			else if (i == N-1){
				T[i][i] = 1.0;                     // CORRIGIR
				P[i] = P_ini;                      // CORRIGIR
			}
			// Células internas:
			else {
				Bi_prev = evaluate_B(P[i-1]);
				Bi_next = evaluate_B(P[i+1]);
				
				Bh_prev = 1.0/(0.5 * (1.0/Bi_prev + 1.0/Bi));
				Bh_next = 1.0/(0.5 * (1.0/Bi + 1.0/Bi_next));

				Wi = (A_x * k_x)/(dx * mu * (Bh_prev));
				Ei = (A_x * k_x)/(dx * mu * (Bh_next));

				T[i][i-1] = Wi;                                         // termo anterior
				T[i][i] =  - ((Vb*phi_ref*c_ref)/(Bi*dt) + Wi + Ei ) ;  // termo na diagonal
				T[i][i+1] = Ei ;                                        // termo posterior
	
				P[i] = - P[i] * (Vb*phi_ref*c_ref)/(Bi*dt);

			}
		}
		std::cout << "Iteracao " << n << std::endl;
		print_array_2D(T);
		//std::cout << "Pressao:" << std::endl;
		//print_array_1D(P);
		P = tdma2(T, P);
	}
	save_data(X, P);
	
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
	print_array_1D(D);
	return D;
}

double evaluate_B(const double p){
	return B0 / (1.0 + c_ref*(p - P_ini));
}

template <typename T>
void print_array_1D(const std::vector<T> A){
	auto nrow = A.size();
	for(int i = 0; i < nrow; i++){
		std::cout << A[i] << ' ';
	}
	std::cout << '\n';
}

template <typename T>
void print_array_2D(const std::vector<std::vector<T>> A){
	auto nrow = A.size();
	auto ncol = A[0].size();
	for(int i = 0; i < nrow; i++){
		for(int j = 0; j < ncol; j++){
			std::cout << std::setw(12) << A[i][j] << ' ';
		}
		std::cout << '\n';
	}
}

void save_data(const std::vector<double>& X, const std::vector<double>& Y){
	std::fstream saver{"output_data.txt", std::ios::out|std::ios::trunc};
	
	const auto N = X.size();
	const auto M = Y.size();
	if ( N != M)
		return;
	for (int i = 0; i < N; i++){
		saver << std::setw(10) << X[i] << " " << std::setw(10) << Y[i] << std::endl;
	}
}

template <typename T>
std::vector<T> linspace(const double xi, const double xf, int Num){
	if (Num == 0 || Num == 1)
		Num++;

	std::vector<T> V (Num, 0.0);

	auto h = (xf - xi) / (Num-1);
	auto n = static_cast<int>(V.size());
	for (int i = 0; i < n; i++){
		V[i] = xi + i*h;
	}

	return V;
}

/*
	if (n==1 && i==2){
		std::cout << B[i-1] << '\t' << Bh_prev << '\t' << Bi << '\t' << Bh_next << '\t' << B[i+1] << std::endl;
	}

*/