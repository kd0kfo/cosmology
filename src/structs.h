/**
 * 
 * This file is part of ray_trace_ellipse, a program which calculates
 * gravitational lensing light deflection angles and ray traces
 * background images.
 *
 * Copyright 2007, 2010 David Coss, PhD
 *
 * ray_trace_elipse is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ray_trace_elipse is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with ray_trace_elipse.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef COSMOLOGY_STRUCTS_H
#define COSMOLOGY_STRUCTS_H 1

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_MPI

typedef struct {
  int rank;
  int num_ranks;
  std::string hostname;
}MPIData;

#endif


struct ray_trace_arguments {

  std::string lensMassDeflectionPlane, sourceBMPFilename;
  std::string savesourcelocations;
  std::string fileNamePrefix;///< prefix of files created in the simulation
  std::string makeMassDensity;// flag to indicate whether the mass should be made before the lensing.
  std::string parameter_name;// filename for parameters

  bool runExistingDeflection, createDeflection, runSim;
  bool useTimeStamp, drawRemovedArea, verbose, useRandom;
  bool includeCricalCurveAndCaustic;///< Whether Critical Curves and Caustics should be used.

  double *offset;
  int gridSpace;
  
  Plane<Double> *lensMassDensity;
};

#endif

