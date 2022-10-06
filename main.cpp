
#include "matrix.hpp"
#include "testing.hpp"
#include <iostream>

template<class T>
void _main() {
	//testing();
	
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