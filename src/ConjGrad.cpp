#include "Solvers/ConjGrad.h"
#include "Solvers/Vector.h"
#include <string.h>	//memcpy

namespace Solvers {

template<typename real>
void ConjGrad<real>::solve() {
	real* r = new real[this->n];
	real* p = new real[this->n];
	real* Ap = new real[this->n];
	real* MInvR = !this->MInv ? r : new real[this->n];
	
	real bNormL2 = Vector<real>::normL2(this->n, this->b);

	//r = this->b - this->A(this->x)
	this->A(r, this->x);
	for (int i = 0; i < (int)this->n; ++i) {
		r[i] = this->b[i] - r[i];
	}

	//MInvR = this->MInv(r)
	if (this->MInv) this->MInv(MInvR, r);	//else MInvR is already r ...
	
	real rDotMInvR = Vector<real>::dot(this->n, r, MInvR);
	real rNormL2 = Vector<real>::normL2(this->n, r);
	this->residual = this->calcResidual(rNormL2, bNormL2, r);
	do {
		if (this->stop()) break;
		memcpy(p, MInvR, sizeof(real) * this->n);
		for (this->iter = 1; this->iter <= this->maxiter; ++this->iter) {
			//alpha = dot(r, this->MInv(r)) / dot(p, this->A(p))
			this->A(Ap, p);
			real alpha = rDotMInvR / Vector<real>::dot(this->n, p, Ap);
			
			for (int i = 0; i < (int)this->n; ++i) {
				this->x[i] += p[i] * alpha;
				r[i] -= Ap[i] * alpha;
			}
			
			rNormL2 = Vector<real>::normL2(this->n, r);
			this->residual = this->calcResidual(rNormL2, bNormL2, r);
			if (this->stop()) break;
			
			if (this->MInv) this->MInv(MInvR, r);
			real nRDotMInvR = Vector<real>::dot(this->n, r, MInvR);
			real beta = nRDotMInvR / rDotMInvR;
	
			for (int i = 0; i < (int)this->n; ++i) {
				p[i] *= beta;
				p[i] += MInvR[i];
			}
			rDotMInvR = nRDotMInvR;
		}
	} while (0);
	
	delete[] r;
	delete[] p;
	delete[] Ap;
	if (this->MInv) delete[] MInvR;
}

template struct ConjGrad<float>;
template struct ConjGrad<double>;

}
