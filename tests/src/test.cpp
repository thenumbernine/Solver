#include "Solvers/GMRes.h"
#include "Solvers/DenseInverse.h"
#include "Solvers/Newton.h"
#include <iostream>
int main() {
	size_t n = 2;
	double x[2] = {-1,-1};
	double b[2] = {1,1};
	double a[2*2] = {
		2,3,
		4,5
	};

	Solvers::HouseholderQR<double>().solveLinear(n, x, a, b);
	std::cout << "QR:" << std::endl;
	std::cout << x[0] << ", " << x[1] << std::endl;

	x[0] = -1;
	x[1] = -1;
	Solvers::GMRes<double>(n, x, b, [&](double *y, const double* x) {
		for (int i = 0; i < n; ++i) {
			double sum = 0;
			for (int j = 0; j < n; ++j) {
				sum = sum + a[i+n*j] * x[j];
			}
			y[i] = sum;
		}
	}, 1e-7, 20, 2).solve();
	std::cout << "GMRes:" << std::endl;
	std::cout << x[0] << ", " << x[1] << std::endl;

	x[0] = -1;
	x[1] = -1;
	Solvers::Newton<double>(n, x, [&](double *y, const double* x) {
		for (int i = 0; i < n; ++i) {
			double sum = 0; 
			for (int j = 0; j < n; ++j) {
				sum = sum + a[i+n*j] * x[j];
			}
			y[i] = sum - b[i];
		}
	}, 1e-7, 100, 1e-7, 20, 2).solve();
	std::cout << "JFNK:" << std::endl;
	std::cout << x[0] << ", " << x[1] << std::endl;
}
