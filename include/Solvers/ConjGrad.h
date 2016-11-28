#pragma once

#include "Solvers/Krylov.h"

namespace Solvers {

template<typename real>
struct ConjGrad : public Krylov<real> {
	typedef Krylov<real> Super;
	using Super::Super; 
	virtual void solve();
};

}
