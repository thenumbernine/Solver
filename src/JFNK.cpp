#include "Solvers/JFNK.h"
#include "Solvers/Vector.h"
#include <string.h>	//memcpy
#include <limits>
#include <math.h>	//isfinite
#include <assert.h>

namespace Solvers {

template<typename real>
JFNK<real>::JFNK(
	size_t n_,
	real* x_,
	Func F_,
	double stopEpsilon_,
	int maxiter_,
	std::function<std::shared_ptr<Krylov<real>>(size_t n, real* F_of_x, real* dx, Func linearFunc)> createLinearSolver)
: n(n_)
, x(x_)
, F(F_)
, lineSearch(&JFNK::lineSearch_bisect)
, maxAlpha(1)
, lineSearchMaxIter(20)
, jacobianEpsilon(1e-6)
, stopEpsilon(stopEpsilon_)
, maxiter(maxiter_)
, dx(new real[n])
, F_of_x(new real[n])
, x_plus_dx(new real[n])
, F_of_x_plus_dx(new real[n])
, x_minus_dx(new real[n])
, F_of_x_minus_dx(new real[n])
, residual(0)
, alpha(0)
, iter(0)
, linearSolver(createLinearSolver(n, F_of_x, dx, [&](real* y, const real* x) {
	return this->krylovLinearFunc(y,x);
}))
{
	//assume x has the initial content
	//use x as the initial dx
	memcpy(dx, x, sizeof(real) * n);
}

template<typename real>
JFNK<real>::~JFNK() {
	delete[] dx;
	delete[] F_of_x;
	delete[] x_plus_dx;
	delete[] F_of_x_plus_dx;
	delete[] x_minus_dx;
	delete[] F_of_x_minus_dx;
}

//solve dF(x[n])/dx[n] x = F(x[n]) for x
template<typename real>
void JFNK<real>::krylovLinearFunc(real* y, const real* dx) {
#if 0
	// https://en.wikipedia.org/wiki/Machine_epsilon
	// machine epsilon for double precision is 2^-53 ~ 1.11e-16
	// sqrt machine epsilon is ~ 1e-8
	//but the Knoll, Keyes 2003 paper says they use 1e-6
	// the paper doesn't say the norm for this particular equation, but in all other equations they use a L2 norm
	real sqrtMachineEpsilon = 1e-6;
	real xNorm = vec_normL2(x);
	real dxNorm = vec_normL2(dx);
	real epsilon = sqrt( (1 + xNorm) / dxNorm ) * sqrtMachineEpsilon;
#else
	//looks like in my config I'm using 1e-6, which is the sqrt(machine epsilon) that the paper describes
	real epsilon = jacobianEpsilon;
#endif

	for (int i = 0; i < n; ++i) {
		x_plus_dx[i] = x[i] + dx[i] * epsilon;
		x_minus_dx[i] = x[i] - dx[i] * epsilon;
	}
	
	F(F_of_x_plus_dx, x_plus_dx);	//F(x + dx * epsilon)
	F(F_of_x_minus_dx, x_minus_dx);	//F(x - dx * epsilon)

	/*
	Knoll, Keyes "Jacobian-Free JFNK-Krylov Methods" 2003 
	 shows "(f(x+epsilon v) - f(x)) / epsilon" for first order
	 and "(f(x+epsilon v) - f(x-epxilon v)) / epsilon" for second order
	  shouldn't the latter have "2 epsilon" on the bottom?
	 shouldn't they both have "|v|" on the bottom?
	
	using epsilon converges in 30 iterations, using 2*epsilon converges in 10
	*/
	real denom = 2. * epsilon;	//jacobianEpsilon;// * vec_normL2((vec_t){.v=dx, .n=n});
	
	//TODO shouldn't this be divided by epsilon times |dx| ?
	//(F(x + dx * epsilon) - F(x - dx * epsilon)) / (2 * |dx| * epsilon)
	for (int i = 0; i < n; ++i) {
		y[i] = (F_of_x_plus_dx[i] - F_of_x_minus_dx[i]) / denom;		//F(x + dx * epsilon) - F(x - dx * epsilon)
	}
}

template<typename real>
real JFNK<real>::residualAtAlpha(real alpha) {
	
	//advance by fraction along dx
	for (int i = 0; i < n; ++i) {
		x_plus_dx[i] = x[i] - dx[i] * alpha;
	}
	
	//calculate residual at x
	F(F_of_x_plus_dx, x_plus_dx);
	
	//divide by n to normalize, so errors remain the same despite vector size
	real stepResidual = Vector<real>::normL2(n, F_of_x_plus_dx) / (real)n;
	
	//for comparison's sake, convert nans to flt_max's
	//this will still fail the isfinite() conditions, *and* it will correctly compare when searching for minimas
	if (stepResidual != stepResidual) stepResidual = std::numeric_limits<real>::max();

	return stepResidual;
}

template<typename real>
real JFNK<real>::lineSearch_none() {
	residual = residualAtAlpha(maxAlpha);
	return maxAlpha;
}

template<typename real>
real JFNK<real>::lineSearch_linear() {
	real alpha = 0;
	residual = std::numeric_limits<real>::max();
	
	for (int i = 0; i <= lineSearchMaxIter; ++i) {
		real stepAlpha = maxAlpha * (real)i / (real)lineSearchMaxIter;
		real stepResidual = residualAtAlpha(stepAlpha);
		if (stepResidual < residual) {
			residual = stepResidual;
			alpha = stepAlpha;
		}
	}

	return alpha;
}

template<typename real>
real JFNK<real>::lineSearch_bisect() {
	real alphaL = 0;
	real alphaR = maxAlpha;
	real residualL = residualAtAlpha(alphaL);
	real residualR = residualAtAlpha(alphaR);

	for (int i = 0; i < lineSearchMaxIter; ++i) {
		real alphaMid = .5 * (alphaL + alphaR);
		real residualMid = residualAtAlpha(alphaMid);
		if (residualMid > residualL && residualMid > residualR) break;
		if (residualMid < residualL && residualMid < residualR) {	//better than both?  see which edge has the lead
			if (residualL <= residualR) {	//<= to prefer sticking closer to the origin in case of equality
				alphaR = alphaMid;
				residualR = residualMid;
			} else {
				alphaL = alphaMid;
				residualL = residualMid;
			}
		} else if (residualMid < residualL) {
			alphaL = alphaMid;
			residualL = residualMid;
		} else {
			alphaR = alphaMid;
			residualR = residualMid;
		}
	}
	
	residual = fmin(residualL, residualR);
	return residualL < residualR ? alphaL : alphaR;
}

/*
performs update of iteration x[n+1] = x[n] - ||dF/dx||^-1 F(x[n])
*/
template<typename real>
void JFNK<real>::update() {	

	//first calc F(x[n])
	F(F_of_x, x);	

	//solve dF(x[n])/dx[n] dx[n] = F(x[n]) for dx[n]
	//treating dF(x[n])/dx[n] = I gives us the (working) explicit version
	linearSolver->solve();

//the next step in matching the implicit to the explicit (whose results are good) is making sure the line search is going the correct distance 
	//update x[n] = x[n] - alpha * dx[n] for some alpha
	alpha = (this->*lineSearch)();

	if (!alpha) {
		//fail code? one will be set in the sim_t calling function at least.
	} else if (!isfinite(residual)) {
		//fail code as well?  likewise, one will be set in the caller.
	} else {

		//don't error here -- instead let the outer loop handle this 
		//if (private->alpha == 0) errorStr("stuck"); 

		//set x[n+1] = x[n] - alpha * dx[n]
		for (int i = 0; i < n; ++i) {
			x[i] -= dx[i] * alpha;
		}
	}
}

template<typename real>
void JFNK<real>::solve() {
	
	for (; iter < maxiter; ++iter) {
		update();
		if (stopCallback && stopCallback()) break;
		if (residual < stopEpsilon) break;
	}
}

template struct JFNK<float>;
template struct JFNK<double>;

}
