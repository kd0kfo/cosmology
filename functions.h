#ifndef FUNCTIONS_H
#include <string>
#include <map>
#include "libmygl/plane.h"

extern double ans;
typedef Plane<Double> plane_t;
std::map<std::string,plane_t> plane_dict;

bool create_plane(const char* name, double n, double m)
{
  plane_dict[name] = plane_t((int)n,(int)m,0.0);
  return true;
}

double print_plane(const char* name)
{
  if(plane_dict.find(name) == plane_dict.end())
    {
      printf("No such plane: %s\n",name);
      return ans;
    }
  const plane_t& plane = plane_dict[name];
  
  printf("dim: %d by %d\n",plane.numberOfRows(), plane.numberOfColumns());
  return ans;
}

#endif

