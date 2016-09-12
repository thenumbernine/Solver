#include "Solvers/JFNK.h"
void test_discreteLaplacian() {
	size_t n = 50;
	double rho[n * n];
	double phi[n * n];
	double h2 = 1;

	//phi,xx + phi,yy = rho
	//solve for phi
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			const double r = 15;
			double dx = (i+.5) - n/2;
			double dy = (j+.5) - n/2;
			rho[i + n * j] = (fabs(dx) < r && fabs(dy) < r)
				? ((dx < 0 ? 1 : -1) * (dy < 0 ? 1 : -1))
				: 0;
		}
	}
	memcpy(phi, rho, sizeof(phi));

	Solvers::Krylov<double>::Func A = [&](double* y, const double* x) {
		for (int i = 0; i < n; ++i) {
			int ip = std::min<int>(i+1, n-1);
			int im = std::max<int>(i-1, 0);
			for (int j = 0; j < n; ++j) {
				int jp = std::min<int>(j+1, n-1);
				int jm = std::max<int>(j-1, 0);
				if (i == 0 || j == 0 || i == n-1 || j == n-1) {
					y[i + n * j] = 0;
				} else if (rho[i + n * j] != 0) {
					y[i + n * j] = rho[i + n * j];
				} else {
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

#if 1	//works!
	Solvers::GMRes<double>(n * n, phi, rho, A, 1e-7, n * n * 10, n * n).solve();
#endif

#if 0
	Solvers::JFNK<double> jfnk(
		n * n,		//n 
		phi,		//initial x
		[&](double* y, const double* phi) {
			A(y, phi);
			for (int i = 0; i < n * n; ++i) {
				y[i] -= rho[i];
			}
		},			//F(x) to minimize
		1e-7,		//stop epsilon
		2, 			//max iter
		[&](size_t n, double* x, double* b, Solvers::JFNK<double>::Func linearFunc) -> std::shared_ptr<Solvers::Krylov<double>> {
			return std::make_shared<Solvers::GMRes<double>>(
				n,			//n
				x,			//initial x
				b,			//b for A(x) = b
				linearFunc,	//linearFunc(y,x) <=> y = A(x)
				1e-7,		//epsilon
				n * n * 10,	//gmres max iteration
				n * n);		//restart iteration
		}			//function for creating Krylov linear function
	);
	jfnk.lineSearch = &Solvers::JFNK<double>::lineSearch_none;
	/*jfnk.stopCallback = [&]()->bool{
		printf("jfnk iter %d alpha %f residual %.16f\n", jfnk.iter, jfnk.alpha, jfnk.residual);
		return false;
	};*/
	/*jfnk.gmres.stopCallback = [&]()->bool{
		printf("gmres iter %d residual %.16f\n", jfnk.gmres.iter, jfnk.gmres.residual);
		return false;
	};*/
	jfnk.solve();
#endif

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			printf("%d %d %.16f\n", i, j, phi[i+n*j]);
		}
		printf("\n");
	}
}
