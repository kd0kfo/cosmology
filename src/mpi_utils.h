#ifndef COSMOLOGY_MPI_UTILS_H
#define COSMOLOGY_MPI_UTILS_H 1

#include <mpi.h>
#include <string>
#include <sstream>

#include "libdnstd/DavidException.h"

#include "libmygl/glellipse.h"

#include "structs.h"

extern MPIData mpi_data;

namespace utils{

  static const int MASTER_RANK = 0;

int init_mpi(int *argc, char ***argv, int *mpi_rank, int *mpi_size);
 void mpi_recombine(GLAlgorithm& gls, MPI_Comm plane_creators);
 void mpi_recombine(Plane<math::Complex> *plane, MPI_Comm plane_creators);
 void mpi_adjust_glellipsebounds(int *glellipseBounds, size_t N);
};//end namespace utils

#endif//COSMOLOGY_MPI_UTILS_H

