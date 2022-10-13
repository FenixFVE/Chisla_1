#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <concepts>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>

using std::vector;

template<std::floating_point T>
class BandMatrix {
public:
	int matrix_size;
	int band_size;
	vector<vector<T>> al;
	vector<vector<T>> au;
	vector<T> di;
	T eps;

	BandMatrix(int matrix_size, int band_size) {
		this->matrix_size = matrix_size;
		this->band_size = band_size;
		int half_band = band_size / 2;
		al = vector<vector<T>>(matrix_size, vector<T>(half_band, 0.0));
		au = vector<vector<T>>(matrix_size, vector<T>(half_band, 0.0));
		di = vector<T>(matrix_size, 0.0);
		eps = (typeid(T) == typeid(float)) ? FLT_EPSILON : DBL_EPSILON;
	}

	// Ux=y, y=?
	void multiply_upper_matrix_on_vector(const vector<T>& input, vector<T>& output) const {
		if (matrix_size != input.size() || input.size() != output.size())
			throw std::exception("error: upper matrix on vector multiplication incompatible sizes");

		// diagonal filled with 1
		for (int i = 0; i < input.size(); ++i) {
			output[i] = input[i];
		}
		

		int half_band = band_size / 2;
		for (int i = 0; i < matrix_size; ++i) {
			int i_1 = i + 1;
			int max = std::min(matrix_size - i_1, half_band);
			for (int j = 0; j < max; ++j) {
				output[i] += au[i][j] * input[j + i_1];
			}
		}
		
	}

	// Lx=y, y=?
	void multiply_lower_matrix_on_vector(const vector<T>& input, vector<T>& output) const {
		if (matrix_size != input.size() || input.size() != output.size())
			throw std::exception("error: lower matrix on vector multiplication incompatible sizes");
		
		// diagonal
		for (int i = 0; i < input.size(); ++i) {
			output[i] = di[i] * input[i];
		}

		int half_band = band_size / 2;
		for (int i = 1; i < matrix_size; ++i) {
			int i_1 = i - half_band;
			int start = std::max(0, -i_1);
			for (int j = start; j < half_band; ++j) {
				output[i] += al[i][j] * input[j + i_1];
			}
		}

	}

	// Ax=y, y=?
	void multiply_matrix_on_vector(const vector<T>& input, vector<T>& output) const {
		if (matrix_size != input.size() || input.size() != output.size())
			throw std::exception("error: matrix on vector multiplication incompatible sizes");

		// diagonal
		for (int i = 0; i < input.size(); ++i) {
			output[i] = di[i] * input[i];
		}

		int half_band = band_size / 2;

		for (int i = 1; i < matrix_size; ++i) {
			int i_1 = i - half_band;
			int start = std::max(0, -i_1);
			for (int j = start; j < half_band; ++j) {
				output[i] += al[i][j] * input[j + i_1];
			}
		}

		for (int i = 0; i < matrix_size; ++i) {
			int i_1 = i + 1;
			int max = std::min(matrix_size - i_1, half_band);
			for (int j = 0; j < max; ++j) {
				output[i] += au[i][j] * input[j + i_1];
			}
		}
	}

	// A=LU, L=?, U=?
	template<std::floating_point F>
	void LU_decompose_matrix() {

		for (T elem : di) {
			if (std::abs(elem) <= eps) {
				throw std::exception("matrix is not decomposable");
			}
		}

		int half_band = band_size / 2;
		for (int i = 0; i < matrix_size; ++i) {
			
			// lower
			int start = std::max(0, half_band - i);
			for (int j = start; j < half_band; ++j) {
				F sum = 0.0;
				int l_j = j - 1;
				int u_i = i + j - half_band - 1;
				int u_j = 0;
				while (l_j >= 0 && u_i >= 0 && u_j < half_band) {
					sum += al[i][l_j--] * au[u_i--][u_j++];
				}
				al[i][j] -= sum;
			}
			
			// center
			{
				F sum = 0.0;
				int l_j = half_band - 1;
				int u_i = i - 1;
				int u_j = 0;
				while (l_j >= 0 && u_i >= 0 && u_j < half_band) {
					sum += al[i][l_j--] * au[u_i--][u_j++];
				}
				di[i] -=  sum;
				if (std::abs(di[i]) <= eps) {
					throw std::exception("matrix is not decomposable");
				}
			}
			
			// upper
			int end = std::min(half_band - 1, band_size - i);
			for (int j = 0; j <= end; ++j) {
				F sum = 0.0;
				int l_j = half_band - 1;
				int u_i = i - 1;
				int u_j = j + 1;
				while (u_i >= 0 && u_j < half_band) {
					sum += al[i][l_j--] * au[u_i--][u_j++];
				}
				au[i][j] = (au[i][j] - sum) / di[i];
			}
			
		}

	}

	// Lx=y, x=?
	void compute_SLE_by_forward_move(const vector<T>& input, vector<T>& output) const {
		if (matrix_size != input.size() || input.size() != output.size())
			throw std::exception("error: SLE computation by forward move incompatible size");

		for (int i = 0; i < input.size(); ++i) {
			output[i] = input[i];
		}

		int half_band = band_size / 2;
		for (int k = 0; k < matrix_size; ++k) {
			output[k] /= di[k];
			int i = k + 1;
			for (int j = half_band - 1; j >= 0 && i < matrix_size; --j, ++i) {
				output[k + half_band - j] -= output[k] * al[i][j];
			}
		}
	}

	// Ux=y, x=?
	void compute_SLE_by_backward_move(const vector<T>& input, vector<T>& output) const {
		if (matrix_size != input.size() || input.size() != output.size())
			throw std::exception("error: SLE computation by backward move incompatible size");

		for (int i = 0; i < input.size(); ++i) {
			output[i] = input[i];
		}

		int half_band = band_size / 2;
		for (int k = matrix_size - 1; k >= 0; --k) {
			int i = k - 1;
			for (int j = 0; i >= 0 && j < half_band; --i, ++j) {
				output[k - 1 - j] -= output[k] * au[i][j];
			}
		}

	}

	// Ax=y, x=?
	template<std::floating_point F>
	void compute_SLE(const vector<T>& input, vector<T>& output, vector<T>& buffer) {
		if (matrix_size != input.size() || input.size() != output.size() || output.size() != buffer.size())
			throw std::exception("error: SLE computation incompatible size");

		this->LU_decompose_matrix<F>();
		this->compute_SLE_by_forward_move(input, buffer);
		this->compute_SLE_by_backward_move(buffer, output);

	}
};

template<std::floating_point T>
BandMatrix<T> read_matrix(std::string matrix_file) {

	std::ifstream stream(matrix_file);
	if (!stream.good()) throw std::exception("invalid matrix file");

	int matrix_size = 0, band_size = 0;
	stream >> matrix_size >> band_size;
	BandMatrix<T> matrix(matrix_size, band_size);

	int half_band = matrix.band_size / 2;
	for (int i = 0; i < matrix_size; ++i) {
		for (int j = 0; j < half_band; ++j) {
			stream >> matrix.al[i][j];
		}
	}
	for (int i = 0; i < matrix_size; ++i) {
		stream >> matrix.di[i];
	}
	for (int i = 0; i < matrix_size; ++i) {
		for (int j = 0; j < half_band; ++j) {
			stream >> matrix.au[i][j];
		}
	}
	
	stream.close();
	return matrix;
}

template<std::floating_point T>
vector<T> read_vector(std::string vector_file) {

	std::ifstream stream(vector_file);
	if (!stream.good()) throw std::exception("invalid vector file");

	int vector_size = 0;
	stream >> vector_size;
	vector<T> vec(vector_size);

	for (int i = 0; i < vector_size; ++i) {
		stream >> vec[i];
	}

	stream.close();
	return vec;
}

template<std::floating_point T>
void print_vector(const vector<T>& vec, std::string vector_file, int precision) {

	std::ofstream stream(vector_file);
	if (!stream.good()) throw std::exception("invalid vector file");

	stream << std::fixed << std::setprecision(precision + 1);
	for (T x : vec) {
		stream << x << '\n';
	}

}

