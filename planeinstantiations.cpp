#include "libmygl/plane.h"

#include "libdnstd/Double.h"
#include "libdnstd/Complex.h"
#include "libdnstd/DStack.h"

template class Plane<Double>;
template class Plane<math::Complex>;
template class Plane<utils::DStack<DString> >;
template class Plane<utils::DStack<Double> >;
