#pragma once

#include "Solver/Krylov.h"
#include "Solver/DenseInverse.h"	//optional DenseInverse solver parameter
#include <stdlib.h>	//size_t

namespace Solver {

/*
source:
Saad, Schultz (1986). "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems." SIAM Journal of Statistical Computations vol. 7 no. 3 July 1986

note that the MInv inherited from Krylov typically doesn't allow in-place operations,
but my GMRES always uses MInv for in-place operations (i.e. the output and input vectors are the same)
*/
template<typename real>
struct GMRES : public Krylov<real> {
	typedef Krylov<real> Super;

	typedef typename Super::Func Func;

	GMRES(
		size_t n,
		real* x,
		const real* b,
		Func A,
		real epsilon = 1e-7,
		int maxiter = -1,
		int restart = -1);
	virtual ~GMRES();
	
	virtual void solve();

protected:
	size_t restart;				//how many iterations to restart.
	
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
