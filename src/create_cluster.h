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
#ifndef CREATE_CLUSTER_CPP
#define CREATE_CLUSTER_CPP

#include "libdnstd/Double.h"
#include "libdnstd/DRandom.h"

class Create_Cluster
{

 public:
  Create_Cluster();
  ~Create_Cluster();
  
  enum DISTRIBUTION_TYPES{NFW/* NFW ellipse (triaxial)*/,
			  NSIS/*Non singular isothermal sphere*/,
			  UNIFORM_ELLIPSE/* Uniform ellipse*/
  };

  /*
   * Generates a new point based on a distribution
   * Positions are in arcseconds.
   * 
   * @param maxDistance int optional max size. Default = 500
   * @param parameters double array of numbers used to build the distribution
   * @return Double 3-D position
   */
  Double calculatePosition(const int clusterType,double const * const parameters)const;
    
  /**
   * Resets the built-in random number generator. Requires a positive LONG
   *
   * @param seed long (positive)
   * @return long last seed used
   * @throw DavidException
   */
  long resetRandom(long seed) throw (DavidException);

  

 private:
  /**
   * Parameterizes distance used to build NFW
   *
   * Array variables:
   * r_ell = position[0];
   * phi = position[1];
   * theta = position[2];
   * q = parameters[0];
   * s = parameters[1];
   * maximumDistance = parameters[2];
   */
  double getNFWDistance(double const * const position, double const * const parameters)const;
  
  /**
   * Parameterizes distance used to build NSIS
   *
   * position = position[0];
   * coreRadius = parameters[0];
   * maximumDistance = parameters[1];
   */
  double getNSISDistance(double const * const position, double const * const parameters)const;

  /**
   * Parameterizes distance used to build Uniform Ellipse
   *
   * position[0] = random distance as a fraction of the elliptical boundary at the random angles, position[1] and position[2].
   * 0 <= position[1] <= 2*pi;
   * 0 <= position[2] <= pi;
   *
   * parameters[0] = axial ratio;
   * parameters[1] = semi-major axis length
   *
   */
  Double getUniformEllipseDistance(double const * const position,double const * const parameters)const;

  utils::DRandom * randy;// Random number generator.

};

#endif

