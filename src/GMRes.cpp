#include "Solvers/GMRes.h"
#include "Solvers/DenseInverse.h"
#include "Solvers/Vector.h"
#include <math.h>
#include <memory.h>
#include <assert.h>

namespace Solvers {

template<typename real>
GMRes<real>::GMRes(size_t n, real* x, const real* b, LinearFunc A, real epsilon, int maxiter, int restart_)
: Super(n, x, b, A, epsilon, maxiter)
, restart(restart_)
{
	if (restart == -1) restart = n;
	r = new real[n];
	v = new real[n * (restart + 1)];
	h = new real[(restart + 1) * restart];
	cs = new real[restart];
	sn = new real[restart];
	y = new real[restart + 1];
	s = new real[restart + 1];
	w = new real[n];
}

template<typename real>
GMRes<real>::~GMRes() {
	delete[] w;
	delete[] s;
	delete[] y;
	delete[] sn;
	delete[] cs;
	delete[] h;
	delete[] v;
	delete[] r;
}

/*
update x
using the gmres steps:
	y = h(1:i, 1:i) \ s(1:i)
	x = x + v(:, 1:i) * y
h is size i * i is the system to invert
s is size i
v is size n * i is the basis
y size up to (restart+1) is a scratch vector
i runs the range from 1 to restart+1
m is restart size / h is sized (m+1) * m - m is used for determining h's element's addresses when back-substituting
n is the size of v - used for linear combinations of y and v to adjust x
*/
template<typename real>
void GMRes<real>::updateX(size_t m, size_t n, real* x, real* h, real* s, real* v, real* y, int i) {
	//y = h(1:i, 1:i) \ s(1:i)
	DenseInverse<real>().backSubstituteUpperTriangular(m+1, i, y, h, s);
	//x = x + v(:, 1:i) * y
	for (int j = 0; j < i; ++j) {
		for (int k = 0; k < n; ++k) {
			x[k] += v[k + n * j] * y[j];
		}
	}
}

template<typename real>
void GMRes<real>::genrot(real* cs, real* sn, real a, real b) {
	if (b == 0) {
		*cs = 1;
		*sn = 0;
	} else if (fabs(b) > fabs(a)) {
		real tmp = a / b;
		*sn = 1. / sqrt(1. + tmp * tmp);
		*cs = tmp * *sn;
	} else {
		real tmp = b / a;
		*cs = 1. / sqrt(1. + tmp * tmp);
		*sn = tmp * *cs;
	}
}

template<typename real>
void GMRes<real>::rotate(real* dx, real* dy, real cs, real sn) {
	real tmp = cs * *dx + sn * *dy;
	*dy = -sn * *dx + cs * *dy;
	*dx = tmp;
}

/*
sources:
Saad, Schultz. 1986. "GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems"
http://www.netlib.org/templates/cpp/gmres.h
http://www.netlib.org/templates/matlab/gmres.m
*/
template<typename real>
void GMRes<real>::solve() {
	size_t n = this->n;
	int m = restart;

	memset(v, 0, sizeof(real) * (m + 1) * n);
	memset(h, 0, sizeof(real) * (m + 1) * m);
	memset(cs, 0, sizeof(real) * m);
	memset(sn, 0, sizeof(real) * m);
	memset(s, 0, sizeof(real) * (m + 1));

	this->iter = 0;

	real bNormL2 = Vector<real>::normL2(n, this->b);

	//r = MInv(b - A(x))
	this->A(r, this->x);
	for (int i = 0; i < n; ++i) {
		r[i] = this->b[i] - r[i];
	}
	if (this->MInv) this->MInv(r, r);
	real rNormL2 = Vector<real>::normL2(n, r);

	this->residual = this->calcResidual(rNormL2, bNormL2, r);
	if (this->stop()) {
	} else {
		int done = 0;
		for (this->iter = 1; this->iter <= this->maxiter && !done;) {
			//v[0] = r/|r|
			for (int i = 0; i < n; ++i) {
				v[i] = r[i] / rNormL2;
			}

			//s = |r|*e1
			memset(s + 1, 0, sizeof(real) * m);
			s[0] = rNormL2;

			//construct orthonormal basis using Gram-Schmidt
			int i = 0;
			for (; i < m; ++i, ++this->iter) {
				//w = MInv(A(v[i]))
				this->A(w, v + n * i);
				if (this->MInv) this->MInv(w, w);
				for (int k = 0; k <= i; ++k) {
					h[k + (m + 1) * i] = Vector<real>::dot(n, w, v + n * k);
					//w = w - h[k][i] * v[k]
					for (int l = 0; l < n; ++l) {
						w[l] -= v[l + n * k] * h[k + (m + 1) * i];
					}
				}
				//h[i+1][i] = |w|
				real wNormL2 = Vector<real>::normL2(n, w);
				//if |w| = 0 then we get a '"lucky" breakdown' according to the GMRES paper
				if (wNormL2 == 0) {
					++i;
					break;
				}
				h[(i+1) + (m+1)*i] = wNormL2;
				//v[i+1] = w / h[i+1][i] = w/|w|
				for (int k = 0; k < n; ++k) {
					v[k + n * (i+1)] = w[k] / h[(i+1) + (m+1)*i];
				}
				//apply Givens rotation
				for (int k = 0; k < i; ++k) {
					rotate(&h[k+(m+1)*i], &h[k+1+(m+1)*i], cs[k], sn[k]);
				}
				//generate plane rotation from h[i][i], h[i+1][i]
				genrot(&cs[i], &sn[i], h[i+(m+1)*i], h[i+1+(m+1)*i]);
				//apply plane rotation to s[i], s[i+1]
#if 0
				//the matlab implementation assumed sn[i] * s[i+1] == 0
				//while the C code relied on the same unoptimized Givens rotation
				//however does that assumption hold when the iteration is restarted?
				//only if (since s is initialized to beta * e1) s is initialized within the restarted loop
				//(which is the case for all my code examples, but not for the original pseudocode, right?)
				rotate(&s[i], &s[i+1], cs[i], sn[i]);
#else
				{
					real tmp = cs[i] * s[i];
					s[i+1] = -sn[i] * s[i];
					s[i] = tmp;
				}
#endif

				//apply plan rotation to h[i][i], h[i+1][i]
				//here's another place where the matlab optimizes but the c doesn't
#if 0
				rotate(&h[i+(m+1)*i], &h[i+1+(m+1)*i], cs[i], sn[i]);
				if (fabs(h[i+1+(m+1)*i]) >= 1e-10) {
					this->residualStr("expected %.16f to be within 1e-10", fabs(h[i+1+(m+1)*i]));
				}
#else
				h[i+(m+1)*i] = cs[i] * h[i+(m+1)*i] + sn[i] * h[i+1+(m+1)*i];
				h[i+1+(m+1)*i] = 0;
#endif

				this->residual = this->calcResidual(fabs(s[i+1]), bNormL2, r);
				if (this->stop()) {
					//updateX(this->x, h, s, v, y, i+1, m, n);
#if 0
this->A(r, this->x);
vec_sub(r, this->b, r);
if (this->MInv) this->MInv(r, r);
#endif
					++i;
					done = 1;
					break;
				}
			}

			//if (done) break;
			updateX(m, n, this->x, h, s, v, y, i);
			if (done) break;

			//r = MInv(b - A(x))
			this->A(r, this->x);
			for (int k = 0; k < n; ++k) {
				r[k] = this->b[k] - r[k];
			}
			if (this->MInv) this->MInv(r, r);
			rNormL2 = Vector<real>::normL2(n, r);
			this->residual = this->calcResidual(rNormL2, bNormL2, r);
			if (this->stop()) {
				break;
			}
		}
	}

}

template struct GMRes<float>;
template struct GMRes<double>;

}
