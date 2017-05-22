#pragma once

#include "Solver/GMRES.h"
#include "Solver/Vector.h"
#include <memory>

namespace Solver {

/*
source:
Knoll, Keyes "Jacobian-Free Newton-Krylov Methods" 2003

LinearSolver is constructed with the createLinearSolver lambda
and must have a .solve() routine to solve for a single iteration
*/
template<typename real>
struct JFNK {

	typedef std::function<void(real* y, const real* x)> Func;

	JFNK(
		size_t n,
		real* x,
		Func F,
		double stopEpsilon,
		int maxiter,
		std::function<std::shared_ptr<Krylov<real>>(size_t n, real* dx, real* F_of_x, Func linearFunc)> createLinearSolver
		= [](size_t n, real* dx, real* F_of_x, Func linearFunc) -> std::shared_ptr<Krylov<real>> {
			return std::make_shared<GMRES<real>>(n, dx, F_of_x, linearFunc, 1e-20, 10 * n, n);
		});
	virtual ~JFNK();

	/*
	perform a single newton iteration
	newton = newton structure
	maxAlpha = scale applied to solved dx update step size
	lineSearchMaxIter = number of divisions to break maxAlpha * dx into when line searching
	*/
	void update();

	/*
	run all iterations until maxiter is reached or until stopEpsilon is reached
	*/
	void solve();

protected:
	size_t n;
	
	//external buffers for the caller to provide

	//current state / at which we are converging about
	//this is *not* the residual that the JFNK minimizes.  That is the private->dx vector.
	//size. equal to gmres.krylov.n. I'm keeping it separate for when gmres is disabled.
	real* x;

public:
	//function which we're minimizing wrt
	Func F;

protected:
	void krylovLinearFunc(real* y, const real* x);

public:
	/*
	don't do any extra searching -- just take the full step
	*/
	real lineSearch_none();

	/*
	assumes current step is in newton->private->dx
	finds best alpha along line, using no more than 'lineSearchMaxIter' iterations 
	reads and writes x, writes to x_plus_dx and F_of_x_plus_dx 
	*/
	real lineSearch_linear();

	/*
	successively subdivide line search to find best alpha
	*/	
	real lineSearch_bisect();

	//line search method
	real (JFNK::*lineSearch)();

	//line search scalar
	real maxAlpha;

	//line search max iter
	int lineSearchMaxIter;

	//epsilon for computing jacobian
	real jacobianEpsilon;

	//stop epsilon
	real stopEpsilon;

	//stop max iter
	int maxiter;

protected:
	virtual real calcResidual(const real* x, real alpha) const;
	
	real residualAtAlpha(real alpha);
	
	//step to solve (df/du)^-1 * du via GMRES
	real* dx;

	//function value at that point
	real* F_of_x;
	
	//temporary buffers
	real* x_plus_dx;
	real* F_of_x_plus_dx;

	real* x_minus_dx;
	real* F_of_x_minus_dx;

public:
	real getResidual() const { return residual; }
	real getAlpha() const { return alpha; }
	int getIter() const { return iter; }
	std::shared_ptr<Krylov<real>> getLinearSolver() const { return linearSolver; }
protected:
	//residual of best solution along the line search
	real residual;

	//alpha of best solution along the line search
	real alpha;

	//current iteration
	int iter;

	std::shared_ptr<Krylov<real>> linearSolver;

public:
	std::function<bool()> stopCallback;
};

}
