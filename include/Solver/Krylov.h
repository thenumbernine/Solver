#pragma once

#include <functional>

namespace Solver {

template<typename real>
struct Krylov {
	/*
	solves the system y = A x
	accepts x as constant, A as parameter structure
	stores result in y
	*/
	typedef std::function<void(real* y, const real* x)> Func;
	
	Krylov(size_t n, real* x, const real* b, Func A, real epsilon_ = 1e-7, int maxiter = -1);
	virtual ~Krylov();
	
	virtual void solve() = 0;

protected:
	//user-provided / initially populated
	size_t n;
	real* x;								//initial guess
	const real* b;								//solution
public:
	Func A;					//linear function
	Func MInv;				//optional.  linear function of inverse of the preconditioner.  currently must be able to operate with the input and output the same memory

	std::function<bool()> stopCallback;

	real epsilon;							//optional.  default 1e-10
	int maxiter;							//optional.  default 'n'

	int getIter() const { return iter; }
	real getResidual() const { return residual; }

public:
	typedef enum {
		NOT_STOPPED,
		STOP_CALLBACK,
		STOP_RESIDUAL_NOT_FINITE,
		STOP_RESIDUAL_WITHIN_EPSILON,
		STOP_REACHED_MAXITER,
	} stopReason_t;
	stopReason_t stopReason;

protected:	
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
};

}


#include <cmath>	//isfinite

namespace Solver {

/*
initialized the object with default arguments
assigns default callbacks
n = problem size

after krylov_init, the caller is still expected to provide x, b, A, and override any other paramters
*/
template<typename real>
Krylov<real>::Krylov(size_t n_, real* x_, const real* b_, Func A_, real epsilon_, int maxiter_)
: n(n_)
, x(x_)
, b(b_)
, A(A_)
, epsilon(epsilon_)
, maxiter(maxiter_)
, stopReason(NOT_STOPPED)
{
	if (maxiter == -1) maxiter = n;
}

template<typename real>
Krylov<real>::~Krylov() {}


template<typename real>
real Krylov<real>::calcResidual(real rNormL2, real bNormL2, const real* r) {
	return rNormL2;
	//most implementations I see rely on L2 norms
	//return bNormL2 == 0 ? rNormL2 : rNormL2 / bNormL2;
	//but an AI trick is that -- with incredibly high dimension space -- the L2 norm fails, and a (weighted) L1 norm works better
	//return vec_normL1(r, n) / fmax(1., vec_normL1(b, n));
	//...but I'm getting oob values... even in the case of re-evaluating the L2 norm ...
	//return vec_normL2(r, n) / fmax(1., vec_normL2(b, n));
}

template<typename real>
bool Krylov<real>::stop() {
	if (stopCallback && stopCallback()) {
		stopReason = STOP_CALLBACK;
		return true;
	}
	if (!std::isfinite(residual)) {
		stopReason = STOP_RESIDUAL_NOT_FINITE;
		return true;
	}
	if (residual < epsilon) {
		stopReason = STOP_RESIDUAL_WITHIN_EPSILON;
		return true;
	}
	if (iter >= maxiter) {
		stopReason = STOP_REACHED_MAXITER;
		return true;
	}
	return false;
}

}
