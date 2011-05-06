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
#ifndef COSMOLOGY_MPI_UTILS_H
#define COSMOLOGY_MPI_UTILS_H 1

#include <mpi.h>
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

