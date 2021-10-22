#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>  //std::max, std::min
#include <functional>  //std::function


double function_s(double x);
double function_phi(double x);
double phi_koren(double teta);
double phi_ospre(double teta);
double phi_albada(double teta);
void solve_via_upwind(std::vector<std::vector<double>>& Q);
void solve_via_lax(std::vector<std::vector<double>>& Q);
void solve_via_beam_warming(std::vector<std::vector<double>>& Q);
void solve_via_fromm(std::vector<std::vector<double>>& Q);
void solve_via_highresolution(std::vector<std::vector<double>>& Q, std::function<double (double)> func);
void fill_initial_cond(std::vector<double>& V);
void solve_exat(std::vector<double>& V, double tf);
void save_data(const std::vector<double>& X, const std::vector<double>& E, const std::vector<double>& A);
std::vector<double> linspace(double start, double end, int num);

// Variáveis do domínio do problema e da simulação:
constexpr double C {0.8};                         // Número de Courant
constexpr double u {0.5};                         // Velocidade de propagação da onda
constexpr int N {50};                             // Número de 'nós' na malha
constexpr double L{1.0};                          // Domínio espacial
constexpr auto dx = L/N;                          // Refinamento da discretização
constexpr auto dt = C*dx/u;                       // Passo de tempo calculado
constexpr double t0 {0.0};                        // Início da simulação
constexpr double tf1 {10.0};                       // Tempo de interesse 1
constexpr auto nsteps = static_cast<int>(tf1/dt); // Número de passos de tempo
constexpr double eta = 1e-8;                      // Constante para previnir divisão por zero


int main (int argc, char* argv[]){
	std::vector<std::string> args(argv, argv + argc);

	auto X = linspace(0.0, L, N);
	auto Ext = linspace(0.0, L, N);

	// Exemplo com método de primeira ordem:
	// std::vector<std::vector<double>> Q_up (nsteps, std::vector<double>(N, 0.0));
	// solve_via_upwind(Q_up);
	// fill_initial_cond(Ext);
	// save_data(X, Ext, Q_up[nsteps-1]);

	std::cout << "\nArgumentos obtidos via linha de comando:" << std::endl;
	for (const auto& e : args){
		std::cout << e << std::endl;
	}

	// Exemplo com método de alta ordem:
	std::vector<std::vector<double>> Q_minmod (nsteps, std::vector<double>(N, 0.0));
	solve_via_highresolution(Q_minmod, phi_albada);
	fill_initial_cond(Ext);
	save_data(X, Ext, Q_minmod[nsteps-1]);
}

double function_s(double x){
	return (x>=0.6 && x<=0.8) ? 1.0 : 0.0;
}

double function_phi(double x){
	return std::exp(-200*(std::pow(x-0.3, 2))) + function_s(x);
}

// Função limitadora de fluxo do método TVD Koren: minmod
double phi_koren(double teta){
	return std::max<double>(0, std::min<double>(2*teta, std::min<double>((1 + 2*teta)/3, 2)));
}

// Função limitadora de fluxo do método TVD Ospre: superbee
double phi_ospre(double teta){
	return 1.5*(teta*teta + teta)/(teta*teta + teta + 1);
}

// Função limitadora de fluxo do método TVD van Albada: SC
double phi_albada(double teta){
	return (teta*teta + teta)/(teta*teta + 1);
}

void solve_via_upwind(std::vector<std::vector<double>>& Q){
	std::cout << "Upwind solver called..." << std::endl;
	fill_initial_cond(Q[0]);
	// Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
		// Iteração nas células espaciais:
		size_t i_prev = 0;
		for (int i = 0; i < N; i++){
			// Periodicidade:
			(i == 0) ? i_prev = N-1 : i_prev = i-1;
			// Calcula os fluxos:
			Q[n][i] = Q[n-1][i] - C*(Q[n-1][i] - Q[n-1][i_prev]);
		}
	}
	std::cout << "Upwind solver finished." << std::endl;
}

void solve_via_lax(std::vector<std::vector<double>>& Q){
	std::cout << "Lax solver called..." << std::endl;
	fill_initial_cond(Q[0]);
	//Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
		// Iteração nas células espaciais:
		size_t i_prev = 0;
		size_t i_next = 0;
		for (int i = 0; i < N; i++){
			// Periodicidade:
			(i == 0) ? i_prev = N-1 : i_prev = i-1;
			(i == N-1) ? i_next = 0 : i_next = i+1;
			// Calcula os fluxos:
			Q[n][i] = Q[n-1][i] - 0.5*C*((Q[n-1][i_next] - Q[n-1][i_prev]) - C*(Q[n-1][i_prev] - 2*Q[n-1][i] + Q[n-1][i_next]));
		}
	}
	std::cout << "Upwind solver finished." << std::endl;
}

void solve_via_beam_warming(std::vector<std::vector<double>>& Q){
	std::cout << "Beam-Warming solver called..." << std::endl;
	fill_initial_cond(Q[0]);
	//Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
		// Iteração nas células espaciais:
		size_t i_prev = 0;
		size_t i_prev2 = 0;
		for (int i = 0; i < N; i++){
			// Periodicidade:
			(i == 0) ? i_prev = N-1 : i_prev = i-1;
			(i_prev == 0) ? i_prev2 = N-1 : i_prev2 = i_prev-1;
			// Calcula os fluxos:
			Q[n][i] = Q[n-1][i] - C*(Q[n-1][i] - Q[n-1][i_prev]) - 0.5*C*(1-C)*(Q[n-1][i] - 2*Q[n-1][i_prev] + Q[n-1][i_prev2]);
		}
	}
}

void solve_via_fromm(std::vector<std::vector<double>>& Q){
	std::cout << "Fromm solver called..." << std::endl;
	fill_initial_cond(Q[0]);
	//Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
		// Iteração nas células espaciais:
		size_t i_next = 0;
		size_t i_prev = 0;
		size_t i_prev2 = 0;
		for (int i = 0; i < N; i++){
			// Periodicidade:
			(i == N-1) ? i_next = 0 : i_next = i+1;
			(i == 0) ? i_prev = N-1 : i_prev = i-1;
			(i_prev == 0) ? i_prev2 = N-1 : i_prev2 = i_prev-1;
			// Calcula os fluxos:
			Q[n][i] = Q[n-1][i] - 0.25*C*((Q[n-1][i_next] + 3*Q[n-1][i] - 5*Q[n-1][i_prev] + Q[n-1][i_prev2]) - C*(Q[n-1][i_next] - Q[n-1][i] - Q[n-1][i_prev] + Q[n-1][i_prev2]));
		}
	}
}

void solve_via_highresolution(std::vector<std::vector<double>>& Q, std::function<double (double)> func){
	std::cout << "High Resolution solver called..." << std::endl;
	fill_initial_cond(Q[0]);
	// Iteração no tempo:
	for (int n = 1; n < nsteps; n++){
		// Iteração nas células espaciais:
		size_t i_next = 0;
		size_t i_next2 = 0;
		size_t i_prev = 0;
		size_t i_prev2 = 0;
		for (int i = 0; i < N; i++){
			// Periodicidade:
			(i == 0) ? i_prev = N-1 : i_prev = i-1;
			(i_prev == 0) ? i_prev2 = N-1 : i_prev2 = i_prev-1;
			(i == N-1) ? i_next = 0 : i_next = i+1;
			(i_next == N-1) ? i_next2 = 0 : i_next2 = i_next+1;
			
			double teta_prev = 0.0;
			double teta_next = 0.0;
			if (u > 0){
				teta_prev = (Q[n-1][i_prev] - Q[n-1][i_prev2])/(Q[n-1][i] - Q[n-1][i_prev] + eta);
				teta_next = (Q[n-1][i] - Q[n-1][i_prev])/(Q[n-1][i_next] - Q[n-1][i] + eta);
				Q[n][i] = Q[n-1][i] - C*(Q[n-1][i] - Q[n-1][i_prev]) - 0.5*C*(1 - C)*(func(teta_next)*(Q[n-1][i_next] - Q[n-1][i]) - func(teta_prev)*(Q[n-1][i] - Q[n-1][i_prev]));
			}
			else{
				teta_prev = (Q[n-1][i_next] - Q[n-1][i])/(Q[n-1][i] - Q[n-1][i_prev] + eta);
				teta_next = (Q[n-1][i_next2] - Q[n-1][i_next])/(Q[n-1][i_next] - Q[n-1][i] + eta);
				Q[n][i] = Q[n-1][i] - C*(Q[n-1][i_next] - Q[n-1][i]) + 0.5*C*(1 + C)*(func(teta_next)*(Q[n-1][i_next] - Q[n-1][i]) - func(teta_prev)*(Q[n-1][i] - Q[n-1][i_prev]));
			}
		}
	}
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
		V[i] = function_phi(i*dx - u*(tf-t0));
}

void save_data(const std::vector<double>& X, const std::vector<double>& E, const std::vector<double>& A){
	auto k = X.size();
	if (k != E.size() || k != A.size()){
		std::cout << "\nAviso: vetores com dimensoes incompativeis." << std::endl;
		return;
	}
	std::fstream printer {"data.txt", std::ios::out|std::ios::trunc};
	printer << "Posicao Exato Numerico" << std::endl;
	printer << std::fixed;
	for (size_t i = 0; i < k; i++){
		printer << std::setprecision(8) << X[i] << ' ' << E[i] << ' ' << A[i] << std::endl;
	}
	std::cout << "\nData saved to \"data.txt\"" << std::endl;
}

std::vector<double> linspace(double start, double end, int num){
	std::vector<double> Res (num, 0.0);
	auto h = (end - start) / (num-1);
	for (int i = 0; i < num; i++){
		Res[i] = start + i*h;
	}
	return Res;
}