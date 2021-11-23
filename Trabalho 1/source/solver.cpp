/*
	Descrição: Solução para a pressão em um escoamento monofásico
	Autor: Naim J. S. Carvalho (njscarvalho@iprj.uerj.br)
	Data: 20 de Novembro de 2021
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

void evaluate_pressure(std::vector<std::vector<double>>& Trans, std::vector<double>& P);
double evaluate_B(const double p);
std::vector<double> solve_by_tdma(const std::vector<std::vector<double>>& Mat, const std::vector<double>& X);
template <typename T> void print_array_1D(const std::vector<T> A);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> std::vector<T> linspace(const double xi, const double xf, int Num);
void save_pressure_data(const std::vector<double>& X, const std::vector<double>& Y);
void save_pressure_evolution(std::fstream& saver, double Dia, double Press1, double Press2 = 0.0, double Press3 = 0.0, double Press4 = 0.0);
double evaluate_trans(const double B1, const double B2);
void print_time_info();

//Fatores de conversão:
constexpr double factor_L {0.3048};         // metro <-> pé
constexpr double factor_mu {0.001};         // Pa.s  <-> cp
constexpr double factor_kx {0.9869233};     // micrometro quadrado <-> darcy
constexpr double factor_P {6.894757};       // kPascal <-> psia
constexpr double factor_c {0.1450377};      // kPascal <-> psia
constexpr double factor_q {0.1589873};      // std m^3/d <-> STD/D
constexpr double factor_t {86400};          // segundos <-> dia
constexpr double betac {86.4e-06};          // transmissibilidades
constexpr double alphac {5.614583};         // volumes

// Variáveis do problema e da simulação:
constexpr double k_x {0.01};               // permeabilidade, em para mili m^2
constexpr double phi_ref {0.25};           // porosidade
constexpr double P_ini {45000};            // Pressão inicial
constexpr double Lx {5000.0};              // dimensão em x
constexpr double Ly {40.0};                // dimensão em y
constexpr double Lz {10.0};                // dimensão em z
constexpr double B0 {1.5};                 // B de referência
constexpr double mu {1.2e-3};              // viscosidade
constexpr double c_ref {6.0e-7};           // compressibilidade
constexpr double Vb {Lx*Ly*Lz};            // volume
constexpr double A_x {Ly*Lz};              // área
constexpr double D {30.0};                 // vazão no lado esquerdo
constexpr int N {32};                      // número de células
constexpr double dx = Lx/N;
constexpr double ti {0.0};
constexpr double tf {365*3};
constexpr double dt {0.5};
constexpr auto nsteps = static_cast<int>((tf - ti)/dt);


int main(int argc, char* argv[]){
	print_time_info();

	std::vector<double> Pos = linspace<double>(0.0, Lx, N);              // vetor para plotar P por x
	std::vector<double> Pressure (N, P_ini);                             // vetor com as pressões
	std::vector<std::vector<double>> T (N, std::vector<double>(N, 0.0)); // matriz de transmissibilidades

	evaluate_pressure(T, Pressure);
	save_pressure_data(Pos, Pressure);

	std::cout << "\nFim da execucao!" << std::endl;
}

void evaluate_pressure(std::vector<std::vector<double>>& Trans, std::vector<double>& P){
	// Posições de registro da evolução temporal (x=0, x=0.1Lx, x=0.2Lx e x=0.9Lx):
	int pos0_0 = 0;
	auto pos0_1 = static_cast<int>(0.1*N);
	auto pos0_2 = static_cast<int>(0.2*N);
	auto pos0_9 = static_cast<int>(0.9*N);


	// Variáveis utilizadas no processo iterativo:
	double Ei = 0.0;                        // Transmissibilidade à direita da célula i
	double Wi = 0.0;                        // Transmissibilidade à esquerda da célula i
	double Bi = 0.0;                        // B(p) na célula i
	double Bi_prev = 0.0;                   // B(p) na célula i - 1
	double Bi_next = 0.0;                   // B(p) na célula i + 1
	double gamma = Vb*phi_ref*c_ref/(B0);

	// Para salvar a evolução no tempo:
	std::string filename {"time_evolution_data.txt"};
	std::fstream time_file{filename, std::ios::out|std::ios::trunc};
	time_file << std::setw(10) << "Dias " << std::setw(10) << "P_0.0Lx P_0.1Lx P_0.2Lx P_0.9Lx " << std::endl;
	std::fstream time_data{filename, std::ios::out|std::ios::app};

	// Iteraçao no tempo:
	for (size_t n = 1; n <= nsteps; n++){

		// Iteração no domínio espacial:
		for (size_t i = 0; i < N; i++){

			Bi = evaluate_B(P[i]);

			// Contorno esquerdo:
			if (i == 0){
				Bi_next = evaluate_B(P[i+1]);
				Bi_prev = evaluate_B(P[i] - D*dx);
				Ei = betac*evaluate_trans(Bi, Bi_next);
				Wi = betac*evaluate_trans(Bi, Bi_prev);

				Trans[i][i] = - (Ei + gamma/dt);
				Trans[i][i+1] = Ei;
				P[i] = -P[i]*(gamma/dt) + Wi*(D*dx);
			}

			// Contorno direito:
			else if (i == N-1){
				Bi_prev = evaluate_B(P[i-1]);
				Wi = betac*evaluate_trans(Bi_prev, Bi);

				Trans[i][i-1] = Wi;
				Trans[i][i] = - (Wi + gamma/dt);
				P[i] = -P[i]*(gamma/dt);
			}

			// Células internas:
			else {
				Bi_prev = evaluate_B(P[i-1]);
				Bi_next = evaluate_B(P[i+1]);

				Wi = betac*evaluate_trans(Bi_prev, Bi);
				Ei = betac*evaluate_trans(Bi, Bi_next);

				Trans[i][i-1] = Wi;                           // termo à esquerda
				Trans[i][i] =  - (Wi + (gamma/dt) + Ei );     // termo na diagonal
				Trans[i][i+1] = Ei;                           // termo à direita

				P[i] = - P[i] * (gamma/dt);
			}

		}
		P = solve_by_tdma(Trans, P);

		// Registrar a evolução no tempo:
		save_pressure_evolution(time_data, n*dt, P[pos0_0], P[pos0_1], P[pos0_2], P[pos0_9]);
	}
}

double evaluate_B(const double p){
	return B0 / (1.0 + c_ref*(p - P_ini));
}

double evaluate_trans(const double B1, const double B2){
	double res1 = (A_x*k_x)/(dx*mu*B1);
	double res2 = (A_x*k_x)/(dx*mu*B2);
	return 0.5*(res1 + res2);
}

void save_pressure_data(const std::vector<double>& X, const std::vector<double>& Y){
	static int a {1};
	std::string init {"pressure_data_"};
	std::string filename = init + std::to_string(a) + ".txt";
	std::fstream saver{filename, std::ios::out|std::ios::trunc};

	saver << std::setw(10) << "x (m) " << std::setw(10) << "Pressao (kPa)" << std::endl;
	const auto N = X.size();
	const auto M = Y.size();
	if ( N != M)
		return;
	for (int i = 0; i < N; i++){
		saver << std::setw(10) << X[i] << " " << std::setw(10) << std::setprecision(10) << Y[i] << std::endl;
	}
	a++;
}

void save_pressure_evolution(std::fstream& saver, double Dia, double Press1, double Press2, double Press3, double Press4){
	saver <<  std::setw(5) << Dia << " " << std::setprecision(9) << std::setw(11) << Press1 << " " << std::setw(11) << Press2 << " " << std::setw(11) << Press3 << " " << std::setw(11) << Press4 << std::endl;
}

std::vector<double> solve_by_tdma(const std::vector<std::vector<double>>& Mat, const std::vector<double>& X){

	auto n_row = Mat.size();
	auto n_col = Mat[0].size();
	if (n_row != n_col || n_row == 0)
		std::cout << "Matriz nao-quadrada detectada" << std::endl;

	std::vector<double> A (n_row, 0.0);
	std::vector<double> B (n_row, 0.0);
	std::vector<double> C (n_row, 0.0) ;
	std::vector<double> D = X;

	// Preenchimento:
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
		C[i] = C[i] / (B[i] - A[i]*C[i-1]);
		D[i] = (D[i] - A[i]*D[i-1]) / (B[i] - A[i]*C[i-1]);
	}

	D[n] = (D[n] - A[n]*D[n-1]) / (B[n] - A[n]*C[n-1]);

	for (int i = n; i-- > 0;){
		D[i] = D[i] - (C[i]*D[i+1]);
	}

	return D;
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
			std::cout << std::setw(6) << A[i][j] << ' ';
		}
		std::cout << '\n';
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

void print_time_info(){
	std::cout << "Tempo total: "  << tf << " dias" << std::endl;
	std::cout << "Passo de tempo: " << dt << " dia(s)" << std::endl;
	std::cout << "Numero de passos de tempo a calcular: " << nsteps  << std::endl;
}