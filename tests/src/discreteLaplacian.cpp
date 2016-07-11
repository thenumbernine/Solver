#include "Solvers/JFNK.h"
void test_discreteLaplacian() {
	size_t n = 100;
	double rho[n * n];
	double phi[n * n];
	double h2 = 1;

	//phi,xx + phi,yy = rho
	//solve for phi
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			rho[i + n * j] = (i == n/2 && j == n/2) ? 1 : 0;
		}
	}

	std::function<void(double*, const double*)> A = [&](double* rho, const double* phi) {
		for (int i = 0; i < n; ++i) {
			int ip = std::min<int>(i+1, n-1);
			int im = std::max<int>(i-1, 0);
			for (int j = 0; j < n; ++j) {
				int jp = std::min<int>(j+1, n-1);
				int jm = std::max<int>(j-1, 0);
				if (i == 0 || j == 0 || i == n-1 || j == n-1) {
					rho[i + n * j] = 0;
				} else {
					rho[i + n * j] =
						(phi[ip + n * j]
						+ phi[im + n * j]
						+ phi[i + n * jp]
						+ phi[i + n * jm]
						- 4 * phi[i + n * j]) / (h2 * 4 * M_PI);
				}
			}
		}
	};

#if 0	//works!
	Solvers::GMRes<double>(n * n, phi, rho, A, 1e-7, n * n * 10, n * n).solve();
#endif

#if 1
	Solvers::JFNK<double> jfnk(n * n, phi, [&](double* y, const double* phi) {
		A(y, phi);
		for (int i = 0; i < n * n; ++i) {
			y[i] -= rho[i];
		}
	}, 1e-7, 2, 1e-7, n * n * 10, n * n);
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
