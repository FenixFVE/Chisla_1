#pragma once

#include "matrix.hpp"
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using std::vector;
using std::string;

void vector_reading_test() {
	vector<double> actual = read_vector<double>("test_vector_1.txt");
	vector<double> expected = { 1, 2, 3, 4, 5 };

	if (actual != expected) {
		throw std::exception("vector reading test failed");
	}
}

// - - -

void matrix_reading_test() {
	BandMatrix<double> actual = read_matrix<double>("test_matrix_1.txt");

	int expected_matrix_size = 5;
	int expected_band_size = 7;
	vector<vector<double>> expected_al = {
		{ 0, 0, 0 },
		{ 0, 0, 5 },
		{ 0, 10, 11 },
		{ 15, 16, 17 },
		{ 20, 21, 22 }
	};
	vector<double> expected_di = { 1, 6, 12, 18, 23 };
	vector<vector<double>> expected_au = {
		{ 2, 3, 4 },
		{ 7, 8, 9 },
		{ 13, 14, 0 },
		{ 19, 0, 0 },
		{ 0, 0, 0}
	};

	auto exception = std::exception("matrix reading test failed");

	if (expected_matrix_size != actual.matrix_size || expected_band_size != actual.band_size) {
		throw exception;
	}

	if (expected_di != actual.di) {
		throw exception;
	}

	for (int i = 0; i < expected_al.size(); ++i) {
		if (expected_al[i] != actual.al[i] || expected_au[i] != actual.au[i]) {
			throw exception;
		}
	}
	
}

void matrix_reading_test_2() {
	BandMatrix<double> actual = read_matrix<double>("test_matrix_2.txt");

	int expected_matrix_size = 5;
	int expected_band_size = 3;
	vector<vector<double>> expected_al = { {0}, {3}, {6}, {9}, {12} };
	vector<double> expected_di = { 1, 4, 7, 10, 13 };
	vector<vector<double>> expected_au = { {2}, {5}, {8}, {11}, {0} };

	auto exception = std::exception("matrix reading test_2 failed");
	
	if (expected_matrix_size != actual.matrix_size || expected_band_size != actual.band_size) {
		throw exception;
	}

	if (expected_di != actual.di) {
		throw exception;
	}

	for (int i = 0; i < expected_al.size(); ++i) {
		if (expected_al[i] != actual.al[i] || expected_au[i] != actual.au[i]) {
			throw exception;
		}
	}

}

// - - -

void multiply_upper_matrix_on_vector_test() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_1.txt");
	vector<double> vec = read_vector<double>("test_vector_1.txt");
	vector<double> actual(vec.size());

	vector<double> expected = { 30, 100, 125, 99, 5 };
	
	matrix.multiply_upper_matrix_on_vector(vec, actual);

	if (actual != expected) {
		throw std::exception("upper matrix multiplication on vector test failed");
	}
}

void multiply_upper_matrix_on_vector_test_2() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> vec = read_vector<double>("test_vector_1.txt");
	vector<double> actual(vec.size());

	vector<double> expected = { 5, 17, 35, 59, 5 };

	matrix.multiply_upper_matrix_on_vector(vec, actual);

	if (actual != expected) {
		throw std::exception("upper matrix multiplication on vector test_2 failed");
	}
}

// - - - 

void multiply_lower_matrix_on_vector_test() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_1.txt");
	vector<double> vec = read_vector<double>("test_vector_1.txt");
	vector<double> actual(vec.size());

	vector<double> expected = { 1, 17, 68, 170, 306 };

	matrix.multiply_lower_matrix_on_vector(vec, actual);

	if (actual != expected) {
		throw std::exception("lower matrix multiplication on vector test failed");
	}
}

void multiply_lower_matrix_on_vector_test_2() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> vec = read_vector<double>("test_vector_1.txt");
	vector<double> actual(vec.size());

	vector<double> expected = { 1, 11, 33, 67, 113 };

	matrix.multiply_lower_matrix_on_vector(vec, actual);

	if (actual != expected) {
		throw std::exception("lower matrix multiplication on vector test_2 failed");
	}
}

// - - -

void compute_SLE_by_backward_move_test() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_1.txt");
	vector<double> vec = { 30, 100, 125, 99, 5 };
	vector<double> actual(5);

	vector<double> expected = { 1, 2, 3, 4, 5 };

	matrix.compute_SLE_by_backward_move(vec, actual);

	if (actual != expected) {
		throw std::exception("backward SLE computition test failed");
	}
}

void compute_SLE_by_backward_move_test_2() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> vec = { 5, 17, 35, 59, 5 };
	vector<double> actual(5);

	vector<double> expected = { 1, 2, 3, 4, 5 };

	matrix.compute_SLE_by_backward_move(vec, actual);
	
	if (actual != expected) {
		throw std::exception("backward SLE computition test_2 failed");
	}
}

// - - -

void compute_SLE_by_forward_move_test() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_1.txt");
	vector<double> vec = { 1, 17, 68, 170, 306 };
	vector<double> actual(5);

	vector<double> expected = { 1, 2, 3, 4, 5 };

	matrix.compute_SLE_by_forward_move(vec, actual);
	
	if (actual != expected) {
		throw std::exception("forward SLE computition test failed");
	}
}

void compute_SLE_by_forward_move_test_2() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> vec = { 1, 11, 33, 67, 113 };
	vector<double> actual(5);

	vector<double> expected = { 1, 2, 3, 4, 5 };

	matrix.compute_SLE_by_forward_move(vec, actual);
	
	if (actual != expected) {
		throw std::exception("forward SLE computition test_2 failed");
	}
}

// - - -

void LU_decompose_matrix_test() {
	// Ax = y
	// LUx = L(Ux) = y
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> x = { 1, 2, 3, 4, 5 }; // x
	vector<double> y(5); //
	matrix.multiply_matrix_on_vector(x, y);

	matrix.LU_decompose_matrix<double>();
	vector<double> t1(5);
	matrix.multiply_upper_matrix_on_vector(x, t1);
	vector<double> t2(5);
	matrix.multiply_lower_matrix_on_vector(t1, t2);

	double eps = 0.000'000'000'1;
	for (int i = 0; i < t2.size(); ++i) {
		if (abs(t2[i] - y[i]) > eps) {
			throw std::exception("LU decomposition test failed");
		}
	}
}

// - - -

void compute_SLE_test() {
	BandMatrix<double> matrix = read_matrix<double>("test_matrix_2.txt");
	vector<double> x = { 1, 2, 3, 4, 5 };
	vector<double> y(5);
	matrix.multiply_matrix_on_vector(x, y);

	vector<double> buffer(5);
	vector<double> x_actual(5);

	matrix.compute_SLE<double>(y, x_actual, buffer);

	double eps = 0.000'000'000'1;
	for (int i = 0; i < x.size(); ++i) {
		if (abs(x[i] - x_actual[i]) > eps) {
			throw std::exception("SLE computition test failed");
		}
	}
}

// - - -

void testing() {
	try {

		vector_reading_test();
		std::cout << "vector reading tested\n";

		matrix_reading_test();
		matrix_reading_test_2();
		std::cout << "matrix reading tested\n";

		multiply_upper_matrix_on_vector_test();
		multiply_upper_matrix_on_vector_test_2();
		std::cout << "upper matrix multiplication on vector tested\n";

		multiply_lower_matrix_on_vector_test();
		multiply_lower_matrix_on_vector_test_2();
		std::cout << "lower matrix multiplication on vector tested\n";

		compute_SLE_by_backward_move_test();
		compute_SLE_by_backward_move_test_2();
		std::cout << "backward SLE computition tested\n";

		compute_SLE_by_forward_move_test();
		compute_SLE_by_forward_move_test_2();
		std::cout << "forward SLE computition tested\n";

		LU_decompose_matrix_test();
		std::cout << "LU decomposition tested\n";

		compute_SLE_test();
		std::cout << "SLE computition tested\n";

		std::cout << "tests result: everything works\n";
	}
	catch (std::exception e) {
		std::cout << e.what() << '\n';
	}
}