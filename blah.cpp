#include "libmygl/plane.h"
#include <map>
#include <iostream>
#include <string>

typedef Plane<Double> plane_t;
int main()
{

  using namespace std;

  try{
  map<string,plane_t> blah;
  blah["3"] = plane_t(3,3,42);
  std::cout << "value: " << blah["3"].getValue(21,1) << std::endl;

  }
  catch(DavidException de)
    {
      de.stdOut();
      return de.getCode();
    }
  return 0;
}
