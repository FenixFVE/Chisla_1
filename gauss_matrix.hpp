#pragma once
#include "matrix.hpp"


template<std::floating_point T>
class GaussMatrix {
public:
	int size;
	vector<vector<T>> matrix;
	T eps;

	GaussMatrix(int size) {
		this->size = size;
		matrix = vector<vector<T>>(size, vector<T>(size, 0.0));
		eps = (typeid(T) == typeid(float)) ? FLT_EPSILON : DBL_EPSILON;
	}

	// Ax=y
	void convert_to_stepped_view(vector<T>& y) {
		T max, coefficient, temporal;
		int row_number;
		for (int i = 0; i < size - 1; ++i) {
			row_number = i;
			max = std::abs(matrix[i][i]);
			for (int j = i + 1; j < size; ++j) {
				temporal = std::abs(matrix[j][i]);
				if (max < temporal) {
					max = temporal;
					row_number = j;
				}
			}
			if (i != row_number) {
				std::swap(matrix[i], matrix[row_number]);
				std::swap(y[i], y[row_number]);
			}
			if (std::abs(matrix[i][i]) > eps) {
				for (int j = i + 1; j < size; ++j) {
					coefficient = matrix[j][i] / matrix[i][i];
					matrix[j][i] = 0.0;
					y[j] -= coefficient * y[i];
					for (int k = i + 1; k < size; ++k) {
						matrix[j][k] -= coefficient * matrix[i][k];
					}
				}
			}
		}
	}

	// Ax=y, y=?
	void multiply_matrix_on_vector(const vector<T>& input, vector<T>& output) {
		if (size != input.size() || input.size() != output.size())
			throw std::exception("error: matrix on vector multiplication incompatible sizes");
		for (int i = 0; i < size; ++i) {
			T sum = 0.0;
			for (int j = 0; j < size; ++j) {
				sum += matrix[i][j] * input[j];
			}
			output[i] = sum;
		}
	}

	// Ax=y, x=?
	void compute_SLE(const vector<T>& input, vector<T>& output) const {
		if (size != input.size() || input.size() != output.size())
			throw std::exception("error: SLE computation incompatible size");

		for (int i = size - 1; i >= 0; --i) {
			T sum = 0;
			for (int j = i + 1; j < size; ++j) {
				sum += matrix[i][j] * output[j];
			}
			if (std::abs(matrix[i][i]) <= eps) {
				throw std::exception("matrix is not decomposable");
			}
			output[i] = (input[i] - sum) / matrix[i][i];
		}
	}
	
};

template<std::floating_point T>
GaussMatrix<T> read_gauss_matrix(std::string matrix_file) {

	std::ifstream stream(matrix_file);
	if (!stream.good()) throw std::exception("invalid gauss matrix file");

	int size;
	stream >> size;
	GaussMatrix<T> matrix(size);

	for (int x = 0; x < size; ++x) {
		for (int y = 0; y < size; ++y) {
			stream >> matrix.matrix[x][y];
		}
	}

	return matrix;
}