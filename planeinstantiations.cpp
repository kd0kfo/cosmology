#include "plane.cpp"
#include "Double.h"
#include "Complex.h"
#include "DStack.h"

template class Plane<Double>;
template class Plane<math::Complex>;
template class Plane<utils::DStack<DString> >;
template class Plane<utils::DStack<Double> >;
