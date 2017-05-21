#pragma once

#include "Solver/Krylov.h"

namespace Solver {

template<typename real>
struct ConjGrad : public Krylov<real> {
	typedef Krylov<real> Super;
	using Super::Super; 
	virtual void solve();
};

}
