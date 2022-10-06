
#include "matrix.hpp"
#include "testing.hpp"
#include <iostream>
#include <cmath>
#include <exception>

template<std::floating_point T, std::floating_point K>
void research_4(int k) {
	try {
		vector<T> x_star = read_vector<T>("vector_x_star.txt");
		BandMatrix<T> A_k = read_matrix<T>("matrix_A.txt");
		vector<T> F_k(x_star.size());
		vector<T> x_k(x_star.size());
		vector<T> buffer(x_star.size());

		A_k.di[0] += pow(10.0, -k);
		
		A_k.multiply_matrix_on_vector(x_star, F_k);
		
		A_k.compute_SLE<K>(F_k, x_k, buffer);

		int precision = (typeid(T) == typeid(float)) ? 7 : 15;
		print_vector<T>(x_k, "vector_x_k.txt", precision);
	}
	catch (std::exception e) {
		std::cout << e.what() << '\n';
	}
}

template<std::floating_point T, std::floating_point K>
void reserch_5(int k) {
	try {

	}
	catch (std::exception e) {
		std::cout << e.what() << '\n';
	}
}



template<class T>
void _main() {
	//testing();
	//research_4<float, float>(5);
	auto t = create_Gilbert_matrix(5);
	int half_band = 4;
	for (int i = 0; i < t.matrix_size; ++i) {
		for (int j = 0; j < half_band; ++j)
			std::cout  << t.al[i][j] << ' ';
		std::cout << "| " << t.di[i] << " | ";
		for (int j = 0; j < half_band; ++j)
			std::cout << t.au[i][j] << ' ';
		std::cout << '\n';
	}
}

int main() {

	bool _float = true;
	if (_float) {
		_main<float>();
	}
	else {
		_main<double>();
	}


	return 0;
}