#ifndef FUNCTIONS_H
#include <string>
#include <map>
#include <iostream>
#include "libmygl/plane.h"
#include "libdnstd/Double.h"
#include "symrec.h"

extern double ans[];
typedef Plane<math::Complex> plane_t;

plane_t *create_plane(double n, double m)
{
  return (new plane_t((int)n,(int)m,0.0));
}

symrec* print_plane(symrec** vars, size_t size)
{
  using namespace std;
  if(vars == NULL || *vars == NULL || !vars[0]->isPlane)
    return NULL;
  const plane_t* plane = vars[0]->value.planeptr;
  printf("dim: %d by %d\n",plane->numberOfRows(), plane->numberOfColumns());
  cout << "max: " << plane->getMaxValue() << "   min: " << plane->getMinValue() << "   sum: " << plane->getTotalValue() << endl;
  return NULL;
}

symrec* clear_plane(symrec** vars,size_t size)
{
  if(vars == NULL || vars[0] == NULL)
    return NULL;;
  if(vars[0]->isPlane)
    {
      delete vars[0]->value.planeptr;
      vars[0]->value.var[0] = vars[0]->value.var[1] = 0;
      vars[0]->isPlane = false;
    }
  return NULL;
}

symrec* add_planes(symrec** vars,size_t size)
{
  if(size < 2 || vars == NULL || vars[0] == NULL || vars[1] == NULL)
    return NULL;
  if(!vars[0]->isPlane || !vars[1]->isPlane)
    return NULL;
  symrec* ptr = (symrec*)malloc(sizeof(symrec));
  ptr->name = (char *) calloc (1,sizeof(char));
  ptr->type = VAR;
  ptr->isPlane = true;
  ptr->value.planeptr = plane_t::addPlanes(vars[0]->value.planeptr,vars[1]->value.planeptr);
  
  return ptr;
}

class DavidException;
symrec* open_plane(symrec** vars,size_t size)
{
  if(vars == NULL || vars[0] == NULL || vars[1] == NULL || vars[1]->name == NULL || size != 2)
    return NULL;
  try{
    vars[0]->value.planeptr = plane_t::readPlane(vars[1]->name);
    vars[0]->isPlane = true;
  }
  catch(DavidException de)
    {
      de.stdErr();
    }
  return NULL;
}

symrec* save_plane(symrec** vars,size_t size)
{
  if(vars == NULL || vars[0] == NULL || !vars[0]->isPlane || vars[1] == NULL || vars[0]->name == NULL || size != 2)
    return NULL;
  try{
    vars[0]->value.planeptr->savePlane(vars[1]->name);
  }
  catch(DavidException de)
    {
      de.stdErr();
    }
  return NULL;
}

#endif

