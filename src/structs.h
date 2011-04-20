#ifdef COSMOLOGY_STRUCTS_H
#define COSMOLOGY_STRUCTS_H 1


#ifdef USE_MPI
#include "mpi_utils.h"

typedef struct{
  int rank;
  int num_ranks;
  std::string hostname;
}MPIData;

#endif

#endif

