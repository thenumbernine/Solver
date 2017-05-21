#pragma once

#include "Solver/Krylov.h"

namespace Solver {

template<typename real>
struct ConjRes : public Krylov<real> {
	typedef Krylov<real> Super;
	using Super::Super;
	virtual void solve();
};

}
