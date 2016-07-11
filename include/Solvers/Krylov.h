#pragma once

#include <functional>

namespace Solvers {

template<typename real>
struct Krylov {
	/*
	solves the system y = A x
	accepts x as constant, A as parameter structure
	stores result in y
	*/
	typedef std::function<void(real* y, const real* x)> LinearFunc;
	
	Krylov(size_t n, real* x, const real* b, LinearFunc A, real epsilon_ = 1e-7, int maxiter = -1);
	virtual ~Krylov();
	
	//user-provided / initially populated
	size_t n;
	real* x;								//initial guess
	const real* b;								//solution
	LinearFunc A;					//linear function
	LinearFunc MInv;				//optional.  linear function of inverse of the preconditioner.  currently must be able to operate with the input and output the same memory

	std::function<bool()> stopCallback;

	real epsilon;							//optional.  default 1e-10
	int maxiter;							//optional.  default 'n'
	
	//member variables
	int iter;								//current iteration
	real residual;						//current residual

	/*
	returns the residual scalar value
	r = residual
	b = solution vector
	*/
	virtual real calcResidual(real rNormL2, real bNormL2, const real* r);
	
	/*
	determines whether to stop
	krylov = krylov args
	returns false to keep going, true to stop
	*/
	virtual bool stop();

	virtual void solve() = 0;
};

}
