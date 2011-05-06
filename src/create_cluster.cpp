/**
 * 
 * This file is part of makecluster, a program that creates a mass
 * distribution data structure (2-D).
 *
 * Copyright 2007, 2010 David Coss, PhD
 *
 * makecluster is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * makecluster is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with makecluster.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>

#include "libdnstd/DavidException.h"

#include "libmygl/Cosmology.h"

#include "create_cluster.h"

Create_Cluster::Create_Cluster()
{
  randy = new utils::DRandom();
}
Create_Cluster::~Create_Cluster()
{
  delete randy;
  randy = 0;
}

long Create_Cluster::resetRandom(long seed) throw (DavidException)
{
  if(randy == 0)
    throw DavidException("For some reason, the random number generator has not be initialized", DavidException::DEFAULT_ERROR_CODE);

  return randy->changeSeed(seed);
}

Double Create_Cluster::calculatePosition(const int clusterType,double const * const parameters)const
{
  Double returnMe = 0;
  double * position = new double[3];
  position[0] = randy->random();
  while(position[0] == 0)
    {
      position[0] = randy->random();
    }

  position[1] = randy->random()*2*D_PI;
  position[2] = randy->random()*D_PI;

  switch(clusterType)
    {
    case Create_Cluster::NFW:
      {
	returnMe = Double(getNFWDistance(position,parameters),position[1],position[2]);
	break;
      }
    case Create_Cluster::NSIS:
      {
	returnMe = getNSISDistance(position,parameters);
	returnMe.setValue(1,position[1]);
	returnMe.setValue(2,D_PI/2);
	break;
      }
    case Create_Cluster::UNIFORM_ELLIPSE:
      {
	while(returnMe == 0)
	  {
	    //position[0] = (randy->random()*2-1)*parameters[1];
	    //position[1] = (randy->random()*2-1)*parameters[1];
	    returnMe = getUniformEllipseDistance(position,parameters);
	    position[0] = randy->random();
	    position[1] = randy->random()*2*D_PI;
	  }
	break;
      }
    default:
      break;
    }
  delete [] position;
  return returnMe;
}


double Create_Cluster::getNFWDistance(double const * const position, double const * const parameters)const 
{

  if(parameters == 0)
    throw DavidException("No Parameters given",DavidException::INVALID_ARGUMENT_ERROR_CODE);

  double r_ell, phi, theta,  q,  s, maximumDistance;
  r_ell = position[0];
  phi = position[1];
  theta = position[2];
  q = parameters[0];
  s = parameters[1];
  maximumDistance = parameters[2];

  double r,x,y,z;

  x =  cos(phi)*sin(theta);
  y =  sin(phi)*sin(theta);
  z =  cos(theta);

  double omega_sqrd = x*x + y*y/(q*q) + z*z/(s*s);
  r = r_ell*r_ell/omega_sqrd;
  r = sqrt(r);

  double a = 1;

  double norm = log(1+a)-a/(a+1);
  //double norm = log(a);

  double returnMe = log(1+r)-r/(1+r);  
  
  return maximumDistance*returnMe/norm;
}

double Create_Cluster::getNSISDistance(double const * const position, double const * const parameters)const
{

  if(parameters == 0)
    throw DavidException("No Parameters given",DavidException::INVALID_ARGUMENT_ERROR_CODE);

  double r,  coreRadius,  maximumDistance,coeff;
  r = position[0];
  coreRadius = parameters[0];
  maximumDistance = parameters[1];
  coeff = log(maximumDistance+sqrt(maximumDistance*maximumDistance+coreRadius*coreRadius))-log(coreRadius);

  return coreRadius*sinh(r*coeff);
}

Double Create_Cluster::getUniformEllipseDistance(double const * const position,double const * const parameters)const
{
  double axialRatio = parameters[0];
  double semiMajorAxis = parameters[1];
  double r = sqrt(position[0]);
  double x = r*cos(position[1]);
  double y = r*sin(position[1]);

  double distance = pow(x,2)+pow(y/axialRatio,2);
  if(distance > 1)
    return 0.0;

  distance = semiMajorAxis*sqrt(x*x+y*y);
  return Double(distance,position[1],position[2]);
  
}



