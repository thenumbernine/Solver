#pragma once

#include "Common/Exception.h"
#include <stddef.h>	//size_t

namespace Solver {

template<typename real>
struct DenseInverse {
	
	/*
	solves for x the equation a x = b
	a is the matrix, sized n * n, memory stored as a[i + n * j] <- column-major
	b is the given solution, size n
	x is the input to solve for, size n
	*/
	virtual void solveLinear(size_t n, real* x, const real* a, const real* b) {
		throw Common::Exception() << "not implemented";
	}

	/*
	solves x for a x = b with a upper-triangular using back-substitution
	a is size m * n stored column major, for m >= n
	b is size m
	x is size n
	x and b can be the same memory

	can operate in-place for x and b using the same memory

	while m is the height of a, only the first n rows are used.
	the only functional reason m is provided is for the column stride 
	*/
	void backSubstituteUpperTriangular(size_t m, size_t n, real* x, const real* a, const real* b);
};

template<typename real>
struct HouseholderQR : public DenseInverse<real> {
	typedef DenseInverse<real> Super;

	/*
	helper function for Householder QR
	a is a m * jmax matrix
	v is scratch vector size m (portions k through m are used)
	*/
	void applyQ(real* a, int m, int k, int jmin, int jmax, real* v);
	
	/*
	Algorithm 10.1 from Trefethen and Bau "Numerical Linear Algebra"
	solves qt^T r = a for unitary q and upper-triangular r
	results:
	qt[m][m] is the transpose of the unitary matrix solving q r = a
	a[m][n] is the linear system on input, and r on output
	m, n are the sizes, for m >= n
	*/
	void householderQR(size_t m, size_t n, real* qt, real* a);

	/*
	solves A x = b for x using Householder QR method
	a is the matrix, sized n * n, memory stored as a[i + n * j] <- column-major
	b is the given solution, size n
	x is the input to solve for, size n
	*/
	virtual void solveLinear(size_t n, real* x, const real* a, const real* b);

	/*
	solves A x = b for x
	a is the matrix, sized m * n, memory stored column major as a[i + m * j]
	for n = b.n, m = x.n
	*/
	void solveLinear_leastSquares(size_t m, size_t n, real* x, const real* a, const real* b);

	/*
	using HouseholderQR, construct an inverse matrix
	a size n * n
	ainv size n * n
	a and ainv can have the same memory
	*/
	void matrixInverse(size_t n, real* ainv, const real* a);
};

}
