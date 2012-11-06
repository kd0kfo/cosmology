#include "libmygl/plane.h"
#include <cstdio>

void usage(){printf("./plane2netcdf <input> <output>\n");}

int main(int argc, char **argv)
{
  Plane<Double> *original = NULL;
  if(argc < 3)
    {
      usage();
      return 1;
    }
  
  original = Plane<Double>::readPlane(argv[1]);
  
  original->writeCDF(argv[2],true);
  
  return 0;
}
