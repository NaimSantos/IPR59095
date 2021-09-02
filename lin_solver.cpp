#include <iostream>
#include <cmath>
#include <vector>


//Protótipos das funções:
template<typename T> bool is_diagonal_dom(const std::vector<std::vector<T>>& A);
template<typename T> void gauss_siedel_solver(const std::vector<std::vector<T>>& A, const std::vector<T>& B, std::vector<T>& X);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> void print_array_1D(const std::vector<T> A);

int main(int argc, char* argv[]){

	// ---- Exemplo com Gauss-Siedel -----//
	std::vector<std::vector<double>> A { {2, 1, 0, 0},
                            {1, 3, 2, 0},
                            {0, 2, 5, 1},
                            {0, 0, 3, 8} };
	std::vector<double> B {7, 19, 31, 52};
	auto n = B.size();
	std::vector<double> X(n, 0.0); // estimativa inicial: um vetor com n zeros

	std::cout << "Matriz de Coeficientes:\n";
	print_array_2D(A);
	std::cout << "Termos independentes: ";
	print_array_1D(B);
	std::cout << "Estimativas iniciais: " ;
	print_array_1D(X);


	// ---- Exemplo com o TDMA ---- //
	std::vector<double> a1 = {0, 1, 2, 3};      //diagonal inferior
	std::vector<double> b1 = {2, 3, 5, 8};      //diagonal principal
	std::vector<double> c1 = {1, 2, 1, 0};      //diagonal superior
	std::vector<double> d1 = {7, 19, 31, 52};   // b

	tdma_solver(a1, b1, c1, d1);
	gauss_siedel_solver(A, B, X);


	std::cout << "\nSolucao (Gauss-Siedel): ";
	print_array_1D (X);
	std::cout << "Solucao (TDMA): ";
	print_array_1D<double>(d1);

}


// Resolve um sistema Ax = B pelo método de Gauss-Siedel, onde X um vetor com as estimativas iniciais:
template<typename T>
void gauss_siedel_solver(const std::vector<std::vector<T>>& A, const std::vector<T>& B, std::vector<T>& X){
	// Caso não seja diagonal dominante, a convergência não é garantida
	/*
	if (!is_diagonal_dom(A))
		return;
	*/
	// Dimensões não são compatíveis:
	if ((A.size() != A[0].size()) || A.size() != B.size())
		return;

	auto Y = B;
	auto E = X;

	int counter {0};
	bool teste {false};
	const double eps {0.000001};

	int m = A.size();
	int n = A[0].size();
	while(!teste && counter<40){
		teste = true;
		for (int i = 0; i < m; i++){
			Y[i] = (B[i] / A[i][i]);
			for (int j = 0; j < n; j++){
				if (j == i)
					continue;
				Y[i] = Y[i] - ((A[i][j] / A[i][i]) * X[j]);
				X[i] = Y[i]; // Escreve em X a estimativa encontrada
			}
			auto res = std::fabs(((X[i] - E[i]) / X[i])) <= eps;
			teste = teste & res;
			E[i] = X[i];
		}
		counter++;
	}
}

/* Resolve um sistema Ax = B pelo algoritmo de Thomas, onde A é uma matriz tridiagonal
	Parâmetros:
	a: objeto std::vetor com a diagonal inferior
	b: objeto std::vetor com a diagonal principal
	c: objeto std::vetor com a diagonal superior
	d: objeto std::vetor com os termos independentes
	c não é preservado nesta versão
	Ao final, d contém a solução.
*/
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


/////---------Outros utilitários-------/////

// Verifica se um objeto std::vector<std::vector> (vetor de vetores), de elementos do tipo T, é uma matriz diagonal dominante:
template<typename T>
bool is_diagonal_dom(const std::vector<std::vector<T>>& A){
	auto n_rows = A.size();
	auto n_col = A[0].size();
	if (n_rows != n_col)
		return false;

	for (int i = 0; i < n_rows; i++){
		double sum = 0.0;
		double diag = std::fabs(A[i][i]);
		for (int j = 0; j < n_col; j++){
			if (i == j)
				continue;
			sum += std::fabs(A[i][j]);
		}
		if (sum > diag)
			return false;
	}
	return true;
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
