#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
//#include <string>


void fill_initial_cond(std::vector<double>& V);
double function_s(double x);
double function_phi(double x);
void solve_via_upwind(std::vector<std::vector<double>>& Q_up);
std::vector<double> linspace(double start, double end, int num);
void solve_exat(std::vector<double>& V, double tf);
void save_data(const std::vector<double>& X, const std::vector<double>& E, const std::vector<double>& A);

// Variáveis do domínio do problema e da simulação:
constexpr double C {0.8};                         // Número de Courant
constexpr double u {0.5};                         // Velocidade de propagação
constexpr int N {501};                            // Número de 'nós' na malha
constexpr double L{1.0};                          // Domínio espacial
constexpr auto dx = L/N;                          // Refinamento da discretização
constexpr auto dt = C*dx/u;                       // Passo de tempo calculado
constexpr double t0 {0.0};                        // Início da simulação
constexpr double tf1 {1.0};                       // Tempo de interesse 1
constexpr double tf2 {5.0};                       // Tempo de interesse 2
constexpr auto nsteps = static_cast<int>(tf1/dt); // Número de passos de tempo
constexpr auto CFL = u*dt/dx;

int main (int argc, char* argv[]){
	auto X = linspace(0.0, L, N);
	auto Ext = linspace(0.0, L, N);
	std::vector<std::vector<double>> Q_up (nsteps, std::vector<double>(N, 0.0));
	
	solve_via_upwind(Q_up);
	solve_exat(Ext, tf1);
	/*
	for (const auto& p : Q_up[nsteps-1])
		std::cout << p << std::endl;
	*/
	save_data(X, Ext, Q_up[nsteps-1]);
	std::cout << "\nExecution reached the end" << std::endl;
}

double function_s(double x){
    return (x>=0.6 && x<=0.8) ? 1.0 : 0.0;
}
double function_phi(double x){
    return std::exp(-200*(std::pow(x-0.3, 2))) + function_s(x);
}

void solve_via_upwind(std::vector<std::vector<double>>& Q_up){
	std::cout << "Upwind solver called..." << std::endl;
	fill_initial_cond(Q_up[0]);
    // Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
        // Iteração nas células espaciais:
		size_t i_previous = 0;
		for (int i = 0; i < N; i++){
            // Periodicidade:
            i == 0 ? i_previous = N-1 : i_previous = i-1;
            // Calcula os fluxos:
            Q_up[n][i] = Q_up[n-1][i] - C*(Q_up[n-1][i] - Q_up[n-1][i_previous]);
		}
	}
	std::cout << "Upwind solver finished." << std::endl;
}
void fill_initial_cond(std::vector<double>& V){
	std::cout << "Filling initial condition..." << std::endl;
	auto k = V.size();
	if (k == 0)
		return;
	for (int i = 0; i < k; i++){
		V[i] = function_phi(i*dx);
	}
}

void solve_exat(std::vector<double>& V, double tf){
	auto p = V.size();
	for (size_t i = 0; i < p; i++)
		V[i] =  function_phi(i*dx - u*(tf-t0));
}

std::vector<double> linspace(double start, double end, int num){
	std::vector<double> Res (num, 0.0);
	auto h = (end - start) / (num-1);
	for (int i = 0; i < num; i++){
		Res[i] = start + i*h;
	}
	return Res;
}
void save_data(const std::vector<double>& X, const std::vector<double>& E, const std::vector<double>& A){
	auto k = X.size();
	if (k != E.size() || k != A.size()){
		std::cout << "\nAviso: vetores com dimensoes incompatíveis." << std::endl;
		return;
	}
	std::fstream printer {"data.txt", std::ios::out|std::ios::trunc};
	printer << "Pos Exato Numerico" << std::endl;
	for (size_t i = 0; i < k; i++){
		printer << X[i] << ' ' << E[i] << ' ' << A[i] << std::endl;
	}
	std::cout << "\nData saved to \"data.txt\"" << std::endl;
}
