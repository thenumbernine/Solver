#include "Solver/DenseInverse.h"

namespace Solver {

template struct DenseInverse<float>;
template struct DenseInverse<double>;

template struct HouseholderQR<float>;
template struct HouseholderQR<double>;

}
