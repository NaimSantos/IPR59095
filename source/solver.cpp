#include <iostream>
#include <iomanip>  //std::setw
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
void save_data(const std::vector<double>& X, const std::vector<double>& Y);
double media_harmonica(const double a, const double b);

constexpr double k_x {10e-15};       // permeabilidade
constexpr double phi_ref {0.25};     // porosidade
constexpr double P_ini {4.5e7};      // Pressão inicial
constexpr double Lx {5000.0};        // dimensão em x
constexpr double Ly {40.0};          // dimensão em y
constexpr double Lz {10.0};          // dimensão em z
constexpr double B0 {1.5};           // B de referência
constexpr double mu {1.2e-3};        // viscosidade
constexpr double c_ref {6.0e-10};    // compressibilidade
constexpr double Vb {Lx*Ly*Lz};      // volume
constexpr double A_x {Ly*Lz};        // área
constexpr int N {50};                // número de células
constexpr double dx = Lx/N;
constexpr double ti {0.0};
constexpr double tf {1000000.0};
constexpr double dt {10};
constexpr auto nsteps = static_cast<int>((tf - ti)/dt);

int main(int argc, char* argv[]){
	std::vector<std::string> args(argv, argv + argc);
	std::cout << "Argumentos obtidos via linha de comando:" << std::endl;
	for (const auto& arg : args){
		std::cout << arg << std::endl;
	}

	auto Pos = linspace<double>(0.0, Lx, N);                             // vetor para plotar P por x
	std::vector<double> Pressure (N, P_ini);                             // vetor com as pressões
	std::vector<std::vector<double>> T (N, std::vector<double>(N, 0.0)); // matriz de transmissibilidades

	evaluate_pressure(T, Pressure);
	save_data(Pos, Pressure);

	/* TODO:
	 função para salvar a evolução da pressao com o tempo em pontos dados.
	*/
}

void evaluate_pressure(std::vector<std::vector<double>>& Trans, std::vector<double>& P){

	// Variáveis utilizadas no processo iterativo:
	double Ei = 0.0;                        // Transmissibilidade à direita da célula i
	double Wi = 0.0;                        // Transmissibilidade à esquerda da célula i
	double Bi = 0.0;                        // B(p) na célula i
	double Bi_prev = 0.0;                   // B(p) na célula i - 1
	double Bi_next = 0.0;                   // B(p) na célula i + 1
	double Bh_prev = 0.0;                   // média harmônica i - 1/2
	double Bh_next = 0.0;                   // média harmônica i + 1/2
	double gamma = 0.0;
	double D = 3.0;                         // vazão no lado esquerdo

	// Iteraçao no tempo:
	for (size_t n = 1; n <= nsteps; n++){

		// Iteração no domínio espacial:
		for (size_t i = 0; i < N; i++){

			Bi = evaluate_B(P[i]);
			gamma = Vb*phi_ref*c_ref/Bi;

			// Contorno esquerdo:
			if (i == 0){
				Bi_next = evaluate_B(P[i+1]);
				Bh_next = media_harmonica(Bi, Bi_next);
				Ei = (A_x*k_x)/(dx*mu*Bh_next);

				Trans[i][i] = - (Ei + gamma/dt);
				Trans[i][i+1] = Ei;
				P[i] = -P[i]*(gamma/dt) + dx*D;
			}
			// Contono direito:
			else if (i == N-1){
				Bi_prev = evaluate_B(P[i-1]);
				Bh_prev = media_harmonica(Bi_prev, Bi);
				Wi = (A_x*k_x)/(dx*mu*Bh_prev);

				Trans[i][i-1] = Wi;
				Trans[i][i] = - (Wi + gamma/dt);
				P[i] = -P[i]*(gamma/dt);
			}
			// Células internas:
			else {
				Bi_prev = evaluate_B(P[i-1]);
				Bi_next = evaluate_B(P[i+1]);

				Bh_prev = media_harmonica(Bi_prev, Bi);
				Bh_next = media_harmonica(Bi, Bi_next);

				Wi = (A_x*k_x)/(dx*mu*Bh_prev);
				Ei = (A_x*k_x)/(dx*mu*Bh_next);

				Trans[i][i-1] = Wi;                                      // termo à esquerda
				Trans[i][i] =  - (Wi + (gamma/dt) + Ei ) ;               // termo na diagonal
				Trans[i][i+1] = Ei ;                                     // termo à direita

				P[i] = - P[i] * (gamma/dt);
			}
		}
		if (n%1000 == 0)
			std::cout << "Passo de tempo " << n << std::endl;

		P = solve_by_tdma(Trans, P);
	}
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
	static int a {1};
	std::string init {"output_data_"};
	std::string filename = init + std::to_string(a) + ".txt";
	std::fstream saver{filename, std::ios::out|std::ios::trunc};

	saver << std::setw(10) << "x (m)" << std::setw(10) << "Pressao (kPa)" << std::endl;
	const auto N = X.size();
	const auto M = Y.size();
	if ( N != M)
		return;
	for (int i = 0; i < N; i++){
		saver << std::setw(10) << X[i] << " " << std::setw(10) << Y[i]/1000 << std::endl;
	}
	a++;
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
double media_harmonica(const double a, const double b){
	return 1.0/(0.5 * (1.0/a + 1.0/b));
}
/*
	if (n==1 && i==2){
		std::cout << B[i-1] << '\t' << Bh_prev << '\t' << Bi << '\t' << Bh_next << '\t' << B[i+1] << std::endl;
	}

*/