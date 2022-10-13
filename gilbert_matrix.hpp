#pragma once
#include "matrix.hpp"

template<std::floating_point T>
BandMatrix<T> create_Gilbert_matrix(int size) {

	BandMatrix<T> band_matrix(size, 2 * (size - 1) + 1);

	int half_band = size - 1;

	for (int i = 0; i < size; ++i) {
		band_matrix.di[i] = 1.0 / (i + i + 1);
	}

	for (int i = 0; i < size; ++i) {

		int i_1 = i + 1;
		int l_j = std::min(half_band - 1, half_band - i);
		for (int j = 0; j < i && l_j < half_band; ++j, ++l_j) {
			band_matrix.al[i][l_j] = 1.0 / (i_1 + j);
		}

		int u_j = 0;
		for (int j = i + 1; j < size && u_j < half_band; ++j, ++u_j) {
			band_matrix.au[i][u_j] = 1.0 / (i_1 + j);
		}
	}

	return band_matrix;
}

template<std::floating_point T>
vector<T> create_vector_x_star(int size) {
	vector<T> vec(size, 0.0);
	std::iota(vec.begin(), vec.end(), 1.0);
	return vec;
}