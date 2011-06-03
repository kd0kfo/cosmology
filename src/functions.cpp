#include <string>
#include <map>
#include <iostream>
#include <algorithm>
#include <fftw3.h>
#include "functions.h"
#include "libmygl/plane.h"
#include "libdnstd/Double.h"
#include "symrec.h"

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
  plane_t* product = plane_t::addPlanes(vars[0]->value.planeptr,vars[1]->value.planeptr);
  delete vars[0]->value.planeptr;
  vars[0]->value.planeptr = product;
  return NULL;
}

symrec* subtract_planes(symrec** vars,size_t size)
{
  if(size < 2 || vars == NULL || vars[0] == NULL || vars[1] == NULL)
    return NULL;
  if(!vars[0]->isPlane || !vars[1]->isPlane)
    return NULL;
  plane_t* product = plane_t::subtractPlanes(vars[0]->value.planeptr,vars[1]->value.planeptr);
  delete vars[0]->value.planeptr;
  vars[0]->value.planeptr = product;
  return NULL;
}

symrec* multiply_planes(symrec** vars,size_t size)
{
  using std::min;
  if(size < 2 || vars == NULL || vars[0] == NULL || vars[1] == NULL)
    return NULL;
  if(!vars[0]->isPlane || !vars[1]->isPlane)
    return NULL;
  int nn[2];
  plane_t *lhs,*rhs;
  lhs = vars[0]->value.planeptr;
  rhs = vars[1]->value.planeptr;
  nn[0] = min(lhs->numberOfRows(),rhs->numberOfRows());
  nn[1] = min(lhs->numberOfColumns(),rhs->numberOfColumns());
  for(int i = 0;i<nn[0];i++)
    for(int j = 0;j<nn[1];j++)
      lhs->setValue(i,j, lhs->getValue(i,j)*rhs->getValue(i,j) );
  return NULL;
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

symrec* copy_plane(symrec** vars, size_t size)
{
  if(vars == NULL || vars[0] == NULL || !vars[0]->isPlane || vars[1] == NULL || vars[0]->name == NULL || size != 2)
    return NULL;
  
  symrec* source = vars[0];
  symrec* dest = vars[1];

  if(dest->isPlane)
    delete dest->value.planeptr;
  else
    dest->isPlane = true;
  dest->value.planeptr = new plane_t(*source->value.planeptr);
  return NULL;
}

symrec* fourier_plane(symrec** vars, size_t size)
{
  if(vars == NULL || *vars == NULL || !vars[0]->isPlane)
    return NULL;
  plane_t* plane = vars[0]->value.planeptr;
  fftw_complex* data;
  fftw_plan plan ;
  int* nn = plane->getDimensions();
  data = (fftw_complex*)malloc(nn[0]*nn[1]*sizeof(fftw_complex));
  plan = fftw_plan_dft(2,nn ,data,data,FFTW_FORWARD,FFTW_ESTIMATE);
  for(int i = 0;i<nn[0];i++)
    for(int j = 0;j<nn[1];j++)
      {
	const plane_t::data_type& val = plane->getValue(i,j);
	data[ j+i*nn[1] ][0] = val.getRealPart();
	data[ j+i*nn[1] ][1] = val.getImaginaryPart();
      }

  fftw_execute(plan);

 for(int i = 0;i<nn[0];i++)
    for(int j = 0;j<nn[1];j++)
      {
	plane_t::data_type val(data[ j+i*nn[1] ][0],data[ j+i*nn[1] ][1]);
	plane->setValue(i,j,val) ;
      }
 
  fftw_destroy_plan(plan);
  fftw_free(data);
  delete nn;
  return NULL;
}

symrec* ifourier_plane(symrec** vars, size_t size)
{
  if(vars == NULL || *vars == NULL || !vars[0]->isPlane)
    return NULL;
  plane_t* plane = vars[0]->value.planeptr;
  fftw_complex* data;
  fftw_plan plan ;
  int* nn = plane->getDimensions(), norm;
  norm = nn[0]*nn[1];
  data = (fftw_complex*)malloc(nn[0]*nn[1]*sizeof(fftw_complex));
  plan = fftw_plan_dft(2,nn ,data,data,FFTW_BACKWARD,FFTW_ESTIMATE);
  for(int i = 0;i<nn[0];i++)
    for(int j = 0;j<nn[1];j++)
      {
	const plane_t::data_type& val = plane->getValue(i,j);
	data[ j+i*nn[1] ][0] = val.getRealPart();
	data[ j+i*nn[1] ][1] = val.getImaginaryPart();
      }

  fftw_execute(plan);

 for(int i = 0;i<nn[0];i++)
    for(int j = 0;j<nn[1];j++)
      {
	plane_t::data_type val(data[ j+i*nn[1] ][0]/norm,data[ j+i*nn[1] ][1]/norm);
	plane->setValue(i,j,val) ;
      }
 
  fftw_destroy_plan(plan);
  fftw_free(data);
  delete nn;
  return NULL;
}
