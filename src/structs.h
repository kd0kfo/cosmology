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

#endif

