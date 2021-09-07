#pragma once

#include <vector>
//Protótipos das funções:

void gauss_siedel_solver(const std::vector<std::vector<double>>& A, const std::vector<double>& B, std::vector<double>& X);
void tdma_solver(const std::vector<double>& a, const std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
template <typename T> void print_array_2D(const std::vector<std::vector<T>> A);
template <typename T> void print_array_1D(const std::vector<T> A);
void read_data();
void read_2();
