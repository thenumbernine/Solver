#pragma once

#include <cmath>
#include <stdlib.h>	//size_t

namespace Solver {

template<typename real>
struct Vector {
	static real dot(size_t n, const real* a, const real* b) {
		real s = 0;
		for (int i = 0; i < (int)n; ++i) {
			s += a[i] * b[i];
		}
		return s;
	}
	
	static real normL2(size_t n, const real* v) {
		return sqrt(dot(n,v,v));
	}
};

}
