#include "Solver/Krylov.h"
#include <math.h>	//isfinite

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
	//most implementations I see rely on L2 norms
	return bNormL2 == 0 ? rNormL2 : rNormL2 / bNormL2;
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
	if (!isfinite(residual)) {
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

template struct Krylov<float>;
template struct Krylov<double>;

}
