#ifndef CREATE_CLUSTER_CPP
#define CREATE_CLUSTER_CPP

#include <math.h>
#include <time.h>
#include <iostream>
#include <fstream>

#include "libdnstd/Double.h"
#include "libdnstd/DRandom.h"
#include "libdnstd/DavidException.h"

#include "libdnstd/Cosmology.h"


#ifdef __VERBOSE__
#define VERBOSE_PRINT(X) std::cout << X << std::endl;
#else
#define VERBOSE_PRINT(X) 
#endif

/**
 * Create a cluster based on a given distribution
 * Author: David Coss
 * Copyright 2009
 */
class Create_Cluster
{

 public:
  Create_Cluster();
  ~Create_Cluster();

  /*
   * Generates a new point based on a distribution
   * Positions are in arcseconds.
   * @param maxDistance int optional max size. Default = 500
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

  
  static const int NFW = 1;///< NFW ellipse (triaxial)
  static const int NSIS = 2;///<Non singular isothermal sphere
  static const int UNIFORM_ELLIPSE = 3;///<Uniform ellipse
  

 private:
  /**
   * r_ell = position[0];
   * phi = position[1];
   * theta = position[2];
   * q = parameters[0];
   * s = parameters[1];
   * maximumDistance = parameters[2];
   */
  double getNFWDistance(double const * const position, double const * const parameters)const;
  
  /**
   * position = position[0];
   * coreRadius = parameters[0];
   * maximumDistance = parameters[1];
   */
  double Create_Cluster::getNSISDistance(double const * const position, double const * const parameters)const;

  /**
   * Uniform Ellipse
   * position[0] = random distance as a fraction of the elliptical boundary at the random angles, position[1] and position[2].
   * 0 <= position[1] <= 2*pi;
   * 0 <= position[2] <= pi;
   *
   * parameters[0] = axial ratio;
   * parameters[1] = semi-major axis length
   *
   */
  Double Create_Cluster::getUniformEllipseDistance(double const * const position,double const * const parameters)const;

  utils::DRandom * randy;

};

#endif

