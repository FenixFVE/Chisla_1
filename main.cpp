
#include "matrix.hpp"
#include "gauss_matrix.hpp"
#include "gilbert_matrix.hpp"
#include "testing.hpp"
#include <iostream>
#include <cmath>
#include <exception>

template<std::floating_point T, std::floating_point K>
void research_4(int k) {
	vector<T> x_star = read_vector<T>("vector_x_star.txt");
	BandMatrix<T> A_k = read_matrix<T>("matrix_A.txt");
	vector<T> F_k(x_star.size(), 0);
	vector<T> x_k(x_star.size(), 0);
	vector<T> buffer(x_star.size(), 0);

	A_k.di[0] += pow(10.0, -k);
	
	A_k.multiply_matrix_on_vector(x_star, F_k);
	
	A_k.compute_SLE<K>(F_k, x_k, buffer);

	int precision = (typeid(T) == typeid(float)) ? 7 : 15;
	print_vector<T>(x_k, "vector_x_k.txt", precision);

}

template<std::floating_point T, std::floating_point K>
void research_5(int size) {
	vector<T> x_star = create_vector_x_star<T>(size);
	BandMatrix<T> A_k = create_Gilbert_matrix<T>(size);
	vector<T> F_k(size, 0);
	vector<T> x_k(size, 0);
	vector<T> buffer(size, 0);

	A_k.multiply_matrix_on_vector(x_star, F_k);

	A_k.compute_SLE<K>(F_k, x_k, buffer);

	int precision = (typeid(T) == typeid(float)) ? 7 : 15;
	print_vector<T>(x_k, "vector_x_k.txt", precision);
}

template<std::floating_point T>
void research_7(int k) {
	vector<T> x_star = read_vector<T>("vector_x_star.txt");
	GaussMatrix<T> A_k = read_gauss_matrix<T>("gauss_matrix_A.txt");
	vector<T> F_k(x_star.size(), 0);
	vector<T> x_k(x_star.size(), 0);

	A_k.multiply_matrix_on_vector(x_star, F_k);

	A_k.matrix[0][0] += pow(10.0, -k);

	A_k.convert_to_stepped_view(F_k);
	A_k.compute_SLE(F_k, x_k);

	int precision = (typeid(T) == typeid(float)) ? 7 : 15;
	print_vector<T>(x_k, "vector_x_k.txt", precision);
}


int main() {

	try {
		std::ifstream stream("config.txt");
		if (!stream.good()) throw std::exception("Cant open config file");
		std::string task, type, sum_type;
		int k;

		stream >> task >> type >> sum_type >> k;

		if ((type != "float" && type != "double") || (sum_type != "float" && sum_type != "double")) {
			throw std::exception("Invalid type");
		}

		if (task == "testing") {
			testing();
		}
		else if (task == "research_4") {
			if (type == "float") {
				if (sum_type == "float") {
					research_4<float, float>(k);
					std::cout << "research_4 float float " << k << '\n';
				}
				else {
					research_4<float, double>(k);
					std::cout << "research_4 float double " << k << '\n';
				}
			}
			else {
				if (sum_type == "float") {
					research_4<double, float>(k);
					std::cout << "research_4 double float " << k << '\n';
				}
				else {
					research_4<double, double>(k);
					std::cout << "research_4 double double " << k << '\n';
				}
			}
		}
		else if (task == "research_5") {
			if (type == "float") {
				if (sum_type == "float") {
					research_5<float, float>(k);
					std::cout << "research_5 float float " << k << '\n';
				}
				else {
					research_5<float, double>(k);
					std::cout << "research_5 float double " << k << '\n';
				}
			}
			else {
				if (sum_type == "float") {
					research_5<double, float>(k);
					std::cout << "research_5 double float " << k << '\n';
				}
				else {
					research_5<double, double>(k);
					std::cout << "research_5 double double " << k << '\n';
				}
			}
		}
		else if (task == "research_7") {
			if (type == "float") {
				research_7<float>(k);
				std::cout << "research_7 float " << k << '\n';
			}
			else {
				research_7<double>(k);
				std::cout << "research_7 double " << k << '\n';
			}
		}
		else {
			throw std::exception("Invalid task");
		}
	}
	catch (std::exception e) {
		std::cout << "Exception: " << e.what() << '\n';
	}


	return 0;
}