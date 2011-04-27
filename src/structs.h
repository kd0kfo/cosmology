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

  std::string lensMassDeflectionPlane, sourceBMPFilename,mainPrefix;
  std::string savesourcelocations;
  std::string fileNamePrefix;///< prefix of files created in the simulation
  std::string makeMassDensity;// flag to indicate whether the mass should be made before the lensing.
  std::string parameter_name;// filename for parameters

  bool runExistingDeflection, createDeflection, runSim, runAsDaemon;
  bool useTimeStamp, drawRemovedArea, verbose, useRandom;
  bool includeCricalCurveAndCaustic;///< Whether Critical Curves and Caustics should be used.

  double *offset;
  int gridSpace;
  
  Plane<Double> *lensMassDensity;
};

#endif

