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
void save_pressure_data(const std::vector<double>& X, const std::vector<double>& Y);
void save_pressure_evolution(std::fstream& saver, const double Dia, const double Press);
double media_harmonica(const double a, const double b);
void teste_tdma();

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
constexpr int N {100};               // número de células
constexpr double dx = Lx/N;
constexpr int factor {86400};        // fator de conversão segundos/dia
constexpr double D {4.0/factor};     // vazão no lado esquerdo
constexpr double ti {0.0};
constexpr double tf {100*factor};
constexpr double dt {0.5*factor};
constexpr auto nsteps = static_cast<int>((tf - ti)/dt);

int main(int argc, char* argv[]){
	/*
	std::vector<std::string> args(argv, argv + argc);
	std::cout << "Argumentos obtidos via linha de comando:" << std::endl;
	for (const auto& arg : args){
		std::cout << arg << std::endl;
	}
	*/
	if (argc > 1)
		teste_tdma();

	auto Pos = linspace<double>(0.0, Lx, N);                             // vetor para plotar P por x
	std::vector<double> Pressure (N, P_ini);                             // vetor com as pressões
	std::vector<std::vector<double>> T (N, std::vector<double>(N, 0.0)); // matriz de transmissibilidades

	evaluate_pressure(T, Pressure);
	save_pressure_data(Pos, Pressure);
	
}

void evaluate_pressure(std::vector<std::vector<double>>& Trans, std::vector<double>& P){

	// Variáveis utilizadas no processo iterativo:
	double Ei = 0.0;                        // Transmissibilidade à direita da célula i
	double Wi = 0.0;                        // Transmissibilidade à esquerda da célula i
	double Bi = 0.0;                        // B(p) na célula i
	double Bi_prev = 0.0;                   // B(p) na célula i - 1
	double Bi_next = 0.0;                   // B(p) na célula i + 1
	double gamma = Vb*phi_ref*c_ref/B0;     // Vb*phi_ref*c_ref/Bi;

	// Para salvar a evolução no tempo:
	std::string filename {"time_evolution_data.txt"};
	std::fstream time_file{filename, std::ios::out|std::ios::trunc};
	time_file << std::setw(10) << "Dias " << std::setw(10) << "Pressao (kPa)" << std::endl;
	std::fstream time_data{filename, std::ios::out|std::ios::app};


	// Iteraçao no tempo:
	for (size_t n = 1; n <= nsteps; n++){

		// Iteração no domínio espacial:
		for (size_t i = 0; i < N; i++){

			Bi = evaluate_B(P[i]);

			// Contorno esquerdo:
			if (i == 0){
				Bi_next = evaluate_B(P[i+1]);
				Ei = (A_x*k_x)/(dx*mu*media_harmonica(Bi, Bi_next));

				Trans[i][i] = - (Ei + gamma/dt);
				Trans[i][i+1] = Ei;
				P[i] = -P[i]*(gamma/dt) + D;
			}

			// Contorno direito:
			else if (i == N-1){
				Bi_prev = evaluate_B(P[i-1]);
				Wi = (A_x*k_x)/(dx*mu*media_harmonica(Bi_prev, Bi));

				Trans[i][i-1] = Wi;
				Trans[i][i] = - (Wi + gamma/dt);
				P[i] = -P[i]*(gamma/dt);
			}

			// Células internas:
			else {
				Bi_prev = evaluate_B(P[i-1]);
				Bi_next = evaluate_B(P[i+1]);

				Wi = (A_x*k_x)/(dx*mu*media_harmonica(Bi_prev, Bi));
				Ei = (A_x*k_x)/(dx*mu*media_harmonica(Bi, Bi_next));

				Trans[i][i-1] = Wi;                                      // termo à esquerda
				Trans[i][i] =  - (Wi + (gamma/dt) + Ei ) ;               // termo na diagonal
				Trans[i][i+1] = Ei ;                                     // termo à direita

				P[i] = - P[i] * (gamma/dt);
			}
		}

		P = solve_by_tdma(Trans, P);

		// Registrar a evolução no tempo:
		/*
		double total_steps = nsteps*1.0;
		if (n%5000 == 0){
			double prg = 100*n/total_steps;
			std::cout << "Progresso: " << std::setprecision(2) << std::fixed <<  prg << " \%" << std::endl;
		}
		*/
		save_pressure_evolution(time_data, n*dt, P[0]);
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
			std::cout << std::setw(6) << A[i][j] << ' ';
		}
		std::cout << '\n';
	}
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
		saver << std::setw(10) << X[i] << " " << std::setw(10) << Y[i]/1000 << std::endl;
	}
	a++;
}

void save_pressure_evolution(std::fstream& saver, const double Dia, const double Press){
	saver << std::setw(10) << Dia/factor << " " << std::setw(10) << Press/1000 << std::endl;
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

void teste_tdma(){

	std::vector<std::vector<double>> A {
										{4.0, 1.0, 0.0, 0.0},
										{1.0, 3.0, 2.0, 0.0},
										{0.0, 2.0, 7.0, 6.0},
										{0.0, 0.0, -4.0, 3.0}
										};
	std::vector<double> B {11, 14, 31.5, 1.5};
	std::vector<double> Res = solve_by_tdma(A, B);
	std::cout << "----Teste com o Algoritmo de Thomas----" << std::endl;
	std::cout << "Matriz A: " << std::endl;
	print_array_2D<double>(A);
	std::cout << "Vetor B: " << std::endl;
	print_array_1D<double>(B);
	std::cout << "\nSolucao via TDMA: ";
	for (const auto& x : Res)
		std::cout << x << " ";
	std::cout << std::endl;
}