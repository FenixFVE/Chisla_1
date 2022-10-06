
#include "matrix.hpp"
#include "testing.hpp"
#include <iostream>
#include <cmath>
#include <exception>

template<std::floating_point T, std::floating_point K>
void research_4(int k, int precision) {
	try {
		vector<T> x_star = read_vector<T>("vector_x_star.txt");
		BandMatrix<T> A_k = read_matrix<T>("matrix_A.txt");
		vector<T> F_k(x_star.size());
		vector<T> x_k(x_star.size());
		vector<T> buffer(x_star.size());

		A_k.di[0] += pow(10.0, -k);
		
		A_k.multiply_matrix_on_vector(x_star, F_k);
		A_k.compute_SLE<K>(F_k, x_k, buffer);

		print_vector<T>(x_k, "vector_x_k.txt", precision);

		for (int i = 0; i < x_k.size(); ++i) {
			x_star[i] -= x_k[i];
		}

		print_vector<T>(x_star, "vector_of_dif.txt", precision);
	}
	catch (std::exception e) {
		std::cout << e.what() << '\n';
	}
}

template<class T>
void _main() {
	//testing();
	research_4<float, float>(0, 7);

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