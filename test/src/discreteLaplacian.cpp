#include "Solvers/ConjGrad.h"
#include "Solvers/ConjRes.h"
#include "Solvers/GMRes.h"
#include "Solvers/JFNK.h"
#include <memory.h>
#include <vector>
#include <algorithm>

#ifndef M_PI	//thanks, MSVC
#define M_PI	3.14159265358979323846264338327950288
#endif

void test_discreteLaplacian() {
	size_t n = 50;
	std::vector<double> rho(n * n);
	std::vector<double> phi(n * n);
	double h2 = 1;

	//phi,xx + phi,yy = rho
	//solve for phi
	for (int i = 0; i < (int)n; ++i) {
		for (int j = 0; j < (int)n; ++j) {
			const double r = 15;
			double dx = (i+.5) - n/2;
			double dy = (j+.5) - n/2;
			rho[i + n * j] = (fabs(dx) < r && fabs(dy) < r)
				? ((dx < 0 ? 1 : -1) * (dy < 0 ? 1 : -1))
				: 0;
			phi[i + n * j] = rho[i + n * j];
		}
	}

	Solvers::Krylov<double>::Func A = [&](double* y, const double* x) {
		for (int i = 0; i < (int)n; ++i) {
			int ip = std::min<int>(i+1, n-1);
			int im = std::max<int>(i-1, 0);
			for (int j = 0; j < (int)n; ++j) {
				int jp = std::min<int>(j+1, n-1);
				int jm = std::max<int>(j-1, 0);
//				if (i == 0 || j == 0 || i == (int)(n-1) || j == (int)(n-1)) {
//					y[i + n * j] = 0;
//				} else if (rho[i + n * j] != 0) {
//					y[i + n * j] = -rho[i + n * j];
//				} else 
				{
					y[i + n * j] =
						(x[ip + n * j]
						+ x[im + n * j]
						+ x[i + n * jp]
						+ x[i + n * jm]
						- 4. * x[i + n * j]) / (h2 * 4. * M_PI);
				}
			}
		}
	};
	
	FILE* solverFile = fopen("solver.txt", "w");
	fprintf(solverFile, "#iter residual alpha\n");

#if 1	//using linear solvers

#if 0	//has a memory access error
	Solvers::ConjGrad<double> solver(n * n, phi.data(), rho.data(), A, 1e-7, n * n * 10);
#endif

#if 0	//works, but has poor convergence
	Solvers::ConjRes<double> solver(n * n, phi.data(), rho.data(), A, 1e-20, -1);
#endif

#if 0	//using gmres with restart proportional to gridsize ... works!
	Solvers::GMRes<double> solver(n * n, phi.data(), rho.data(), A, 1e-7, n * n * 10, n * n);
#endif

#if 1	//using gmres with restart proportional to constant ... works! but slower, of course ... O(exp(x)) instead of O(exp(x^2))
	Solvers::GMRes<double> solver(n * n, phi.data(), rho.data(), A, 1e-7, n * n * 10, 10);
#endif
	
	solver.stopCallback = [&]()->bool{
		fprintf(solverFile, "%d\t%.16f\n", solver.getIter(), solver.getResidual());
		return false;
	};
#endif

#if 0	//using JFNK nonlinear solver
	FILE* gmresFile = fopen("gmres.txt", "w");
	fprintf(gmresFile, "#jfnk_iter gmres_iter residual\n");
	Solvers::JFNK<double> solver(
		n * n,		//n 
		phi.data(),	//x
		[&](double* y, const double* phi) {
			A(y, phi);
			for (int i = 0; i < (int)(n * n); ++i) {
				y[i] -= rho[i];
				//y[i] *= y[i];
			}
		},			//F(x) to minimize
		1e-7,		//stop epsilon
		n * n,		//max iter
		[&](size_t n, double* x, double* b, Solvers::JFNK<double>::Func linearFunc) -> std::shared_ptr<Solvers::Krylov<double>> {
			return std::make_shared<Solvers::GMRes<double>>(
				n,			//n
				x,			//initial x
				b,			//b for A(x) = b
				linearFunc,	//linearFunc(y,x) <=> y = A(x)
				1e-7,		//epsilon
				n * 10,	//gmres max iteration
				n);		//restart iteration
		}			//function for creating Krylov linear function
	);
	solver.lineSearch = &Solvers::JFNK<double>::lineSearch_none;
	//solver.lineSearch = &Solvers::JFNK<double>::lineSearch_bisect;
	solver.lineSearchMaxIter = 20;	//250
	solver.stopCallback = [&]()->bool{
		fprintf(solverFile, "%d\t%.16f\t%f\n", solver.getIter(), solver.getResidual(), solver.getAlpha());
		fflush(solverFile);
		fprintf(gmresFile, "\n");
		fflush(gmresFile);
		return false;
	};
	std::shared_ptr<Solvers::GMRes<double>> gmres = std::dynamic_pointer_cast<Solvers::GMRes<double>>(solver.getLinearSolver());
	gmres->stopCallback = [&]()->bool{
		fprintf(gmresFile, "%d\t%d\t%.16f\n", solver.getIter(), gmres->getIter(), gmres->getResidual());
		fflush(gmresFile);
		return false;
	};
#endif
	
	solver.solve();

	printf("#x y phi\n");
	for (int i = 0; i < (int)n; ++i) {
		for (int j = 0; j < (int)n; ++j) {
			printf("%d %d %.16f\n", i, j, phi[i+n*j]);
		}
		printf("\n");
	}

	fclose(solverFile);
#if 0	//if we're using the jfnk
	fclose(gmresFile);
#endif
}
