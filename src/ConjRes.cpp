#include "Solvers/ConjRes.h"
#include "Solvers/Vector.h"
#include <string.h>	//memcpy

namespace Solvers {

template<typename real>
void ConjRes<real>::solve() {
	real* r = new real[this->n];
	real* p = new real[this->n];
	real* Ap = new real[this->n];
	real* Ar = new real[this->n];
	real* MInvAp = !this->MInv ? Ap : new real[this->n];
	
	real bNormL2 = Vector<real>::normL2(this->n, this->b);

	//r = this->MInv(this->b - this->A(this->x))
	this->A(r, this->x);
	for (int i = 0; i < this->n; ++i) {
		r[i] = this->b[i] - r[i];
	}
	if (this->MInv) this->MInv(r, r);

	real rNormL2 = Vector<real>::normL2(this->n, r);
	this->residual = this->calcResidual(rNormL2, bNormL2, r);

	do {
		if (this->stop()) break;

		this->A(Ar, r);
		real rAr = Vector<real>::dot(this->n, r, Ar);
		memcpy(p, r, sizeof(real) * this->n);
		this->A(Ap, p);
		for (this->iter = 1; this->iter <= this->maxiter; ++this->iter) {
			//alpha = dot(r, this->A(r)) / dot(this->A(p), this->MInv(this->A(p)))
			if (this->MInv) this->MInv(MInvAp, Ap);
			real alpha = rAr / Vector<real>::dot(this->n, Ap, MInvAp);
			
			for (int i = 0; i < this->n; ++i) {
				this->x[i] += p[i] * alpha;
				r[i] -= MInvAp[i] * alpha;
			}
			
			rNormL2 = Vector<real>::normL2(this->n, r);
			this->residual = this->calcResidual(rNormL2, bNormL2, r);
			if (this->stop()) break;
		
			this->A(Ar, r);
			real nrAr = Vector<real>::dot(this->n, r, Ar);
			real beta = nrAr / rAr;

			rAr = nrAr;

			for (int i = 0; i < this->n; ++i) {
				p[i] *= beta;
				p[i] += r[i];
				Ap[i] *= beta;
				Ap[i] += Ar[i];
			}
		}
	} while (0);

	delete[] r;
	delete[] p;
	delete[] Ap;
	delete[] Ar;
	if (this->MInv) delete[] MInvAp;
}

template struct ConjRes<float>;
template struct ConjRes<double>;

}
