/**
 * 
 * This file is part of physcalc, an interactive calculator utility
 * designed to carry out lensing calculations and to manipulate plane
 * data structures.
 *
 * Copyright 2007, 2010 David Coss, PhD
 *
 * physcalc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * physcalc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with physcalc.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FUNCTIONS_H

#include "physcalc.h"
#include "symrec.h"
typedef Plane<math::Complex> plane_t;


plane_t *create_plane(double n, double m);
symrec* print_plane(symrec** vars, size_t size);
symrec* clear_plane(symrec** vars,size_t size);
symrec* add_planes(symrec** vars,size_t size);
symrec* subtract_planes(symrec** vars,size_t size);
symrec* multiply_planes(symrec** vars,size_t size);
symrec* open_plane(symrec** vars,size_t size);
symrec* copy_plane(symrec** vars, size_t size);
symrec* save_plane(symrec** vars,size_t size);
symrec* fourier_plane(symrec** vars, size_t size);
symrec* ifourier_plane(symrec** vars, size_t size);

#endif

