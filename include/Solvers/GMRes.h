#pragma once

#include "Solvers/Krylov.h"
#include "Solvers/DenseInverse.h"	//optional DenseInverse solver parameter
#include <stdlib.h>	//size_t

namespace Solvers {

/*
source:
Saad, Schultz (1986). "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems." SIAM Journal of Statistical Computations vol. 7 no. 3 July 1986
*/
template<typename real>
struct GMRes : public Krylov<real> {
	typedef Krylov<real> Super;

	typedef typename Super::LinearFunc LinearFunc;

	GMRes(
		size_t n,
		real* x,
		const real* b,
		LinearFunc A,
		real epsilon = 1e-7,
		int maxiter = -1,
		int restart = -1);
	virtual ~GMRes();
	
	size_t restart;				//how many iterations to restart.
	
	virtual void solve();

	//n = problem size, m = restart
	real* r;	//[n] residual
	real* v;	//[n,m+1] linear projection
	real* h;	//[m+1,m] lower dimensional space mapping - upper triangular matrix
	real* cs;	//[m] cosine of Givens rotations
	real* sn;	//[m] sine of Givens rotations
	real* y;	//[m+1] back-substitution of h from s
	real* s;	//[m+1] progressively solved
	real* w;	//[n] vHat in the paper, solved with h via elimination

	void updateX(size_t m, size_t n, real* x, real* h, real* s, real* v, real* y, int i);
	void genrot(real* cs, real* sn, real a, real b);
	void rotate(real* dx, real* dy, real cs, real sn);
};

}
