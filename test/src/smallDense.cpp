#include "Solvers/GMRes.h"
#include "Solvers/DenseInverse.h"
#include "Solvers/JFNK.h"
#include <iostream>

void test_smallDense() {
	size_t n = 3;
	double x[3] = {-1,-1,-1};
	double b[3] = {0,0,1};
	double a[3*3] = {
		1,0,0,
		0,1,0,
		0,0,1
	};

	Solvers::HouseholderQR<double>().solveLinear(n, x, a, b);
	std::cout << "QR:" << std::endl;
	std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;

	x[0] = -1;
	x[1] = -1;
	x[2] = -1;
	Solvers::GMRes<double>(n, x, b, [&](double *y, const double* x) {
		for (int i = 0; i < (int)n; ++i) {
			double sum = 0;
			for (int j = 0; j < (int)n; ++j) {
				sum = sum + a[i+n*j] * x[j];
			}
			y[i] = sum;
		}
	}, 1e-7, 20, 2).solve();
	std::cout << "GMRes:" << std::endl;
	std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;

	x[0] = -1;
	x[1] = -1;
	x[2] = -1;
	Solvers::JFNK<double> jfnk(n, x, [&](double* y, const double* x) {
		y[0] = 0;
		y[1] = -x[2];
		y[2] = x[1] - 1;
	}, 1e-10, 100, [&](size_t n, double* x, double* b, Solvers::JFNK<double>::Func linearFunc) -> std::shared_ptr<Solvers::Krylov<double>> {
		return std::make_shared<Solvers::GMRes<double>>(
			n,
			x,
			b,
			linearFunc,
			1e-20,
			10*n,
			n);
	});
	jfnk.lineSearch = &Solvers::JFNK<double>::lineSearch_bisect;
	jfnk.lineSearchMaxIter = 100;
	jfnk.stopCallback = [&]()->bool{
		printf("jfnk %d\t%.16f\t%f\n", jfnk.getIter(), jfnk.getResidual(), jfnk.getAlpha());
		return false;
	};
	std::shared_ptr<Solvers::GMRes<double>> gmres = std::dynamic_pointer_cast<Solvers::GMRes<double>>(jfnk.getLinearSolver());
	gmres->stopCallback = [&]()->bool{
		printf("gmres %d\t%d\t%.16f\n", jfnk.getIter(), gmres->getIter(), gmres->getResidual());
		return false;
	};


	jfnk.solve();
	
	std::cout << "JFNK:" << std::endl;
	std::cout << x[0] << ", " << x[1] << ", " << x[2] << std::endl;
}

