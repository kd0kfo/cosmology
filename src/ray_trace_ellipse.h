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
#ifndef RAY_TRACE_CPP
#define RAY_TRACE_CPP

#include <getopt.h>

#include "libdnstd/DRandom.h"

#include "libmygl/densityprofile.h"
#include "libmygl/plane.h"


enum option_keys{ SOURCE_BMP = 1,DEFLECTION_MAP,FORCE_RUN,
		  STOP_RUN, FILE_PREFIX, USE_GRID, RBG_COLOR, XOFFSET, YOFFSET,
		  SUBTRACT_PLANES, ADD_PLANES, SQUARE_PLANE,
		  SAVESOURCE_LOC,DRAW_REMOVED_AREA};

static const struct option ray_trace_options[] = 
  {
    {"newlens", required_argument, NULL,'n'},
    {"createsurfacemassdensity",required_argument,NULL,'c'},
    {"usesurfacemassdensity",required_argument,NULL,'m'},
    {"sourcebmp",required_argument,NULL,SOURCE_BMP},
    {"deflectionmap",required_argument,NULL,DEFLECTION_MAP},
    {"lens",required_argument,NULL,'l'},
    {"forcerun",no_argument,NULL,FORCE_RUN},
    {"stoprun",no_argument,NULL,STOP_RUN},
    {"prefix",required_argument,NULL,FILE_PREFIX},
    {"usegrid",required_argument,NULL,USE_GRID},
    {"bgcolor",required_argument,NULL,RBG_COLOR},
    {"xoffset",required_argument,NULL,XOFFSET},
    {"yoffset",required_argument,NULL,YOFFSET},
    {"parameters",required_argument,NULL,'p'},
    {"glellipsebounds",required_argument,NULL,'g'},
    {"timestamp",no_argument,NULL,'t'},
    {"addplanes",required_argument,NULL,ADD_PLANES},
    {"subtractplanes",required_argument,NULL,SUBTRACT_PLANES},
    {"squareplane",required_argument,NULL,SQUARE_PLANE},
    {"savesourcelocations",no_argument,NULL,SAVESOURCE_LOC},
    {"daemon",no_argument,NULL,'d'},
    {"drawremovedarea",no_argument,NULL,DRAW_REMOVED_AREA},
    {"help",no_argument,NULL,'h'},
    {"verbose",no_argument,NULL,'v'},
    {0,0,0,0}
  };


/** \mainpage Documentation of the Main Executables.
 * Copyright 2007, 2010 David Coss, PhD
 */

template <class T> void buildSource(Plane<T> *, double * sourceParams,T);///< constructs the source plane.
template <class T> void constructLens(Plane<T> * source, Plane<T> * lens, T valueOfSource, T lensValue);///< constructs the lens plane.

/**
 * Main start method.
 * Before this method is called arguments are parsed and then sent to sub_main
 */
int sub_main(struct ray_trace_arguments *args);

/**
 * Intermediate Main Method.
 * This method is called before sub_main. It acts as a way to route information
 * whether it comes from main(int,char**) or BOINC
 */
int super_main(int argc, char** argv);

/**
 * Writes Parameter info.
 * Writes all of the parameters to two text files: 
 * one that the simulator can re-read and one 
 * that is better formatted for people.
 *
 * @param params General Parameters.
 * @param lensParams Lens Parameters.
 * @param sourceParams Array of Parameters for each source
 */
bool writeParameters(const struct ray_trace_arguments *args,
		     const struct general_parameters& params,const struct lens_parameters& lensParams,const struct source_parameters& sourceParams);

/**
 * The main simulation routine.
 * This is the meat of the program. The actually does the ray tracing.
 *
 * @param lens Plane<Double> representing the lens plane
 * @param sources Plane<Double> containing the source background images
 * @param massDensity DensityProfile containing the 2D mass density
 * @throw DavidException Exception thrown upon error (if you're lucky)
 */
int simulation(struct ray_trace_arguments *args,Plane<Double> **lens, Plane<Double> **sources, DensityProfile **massDensity) throw (DavidException);

/**
 * Creates parameters for a lens to be used in the simulation
 *
 * @param lensParams double array of parameters to use in the simulation
 * @param specificLens int specific lens to use
 */
bool createLensParams(struct lens_parameters *lensParams, int specificLens);
int run_simulation(struct ray_trace_arguments *args);

/**
 * Fills the argument structure will default data
 */
void default_arguments(struct ray_trace_arguments *args);

/**
 * Parses arguments sent from command line.
 */
int parseArgs(int argc, char** argv,struct ray_trace_arguments *args);

/**
 * Reads the Parameters from the Parameters file.
 *
 * @param params general parameters.
 * @param lensParams Lens parameters.
 * @param sourceParams Array of Parameters for each Source.
 */
void loadParameters(struct ray_trace_arguments *args,struct general_parameters * params,struct lens_parameters * lensParams,struct source_parameters *sourceParams);

/**
 * Creates a 2 Dimensional Mass Array from lens parameters.
 */
void createSurfaceMassDensity(const std::string& fileName, struct ray_trace_arguments& args);

/**
 * Returns a std::string containing the current time.
 * Format: Wed May 16 16:44:42 2007
 */
std::string getTime();

/**
 * Takes the largest dimension of a Plane and makes it square.
 *
 * @param fileName const char* name of file
 */
void squarePlane(const char * fileName);

utils::DRandom * randy;///< Random number Generator.
Plane<math::Complex> * deflectionPlane;
int overallParamNumber;
Double bgColor;
void drawEllipse(const struct ray_trace_arguments *args,const char * parameters);
void verbosePrint(const struct ray_trace_arguments *args, const char * frmt, ...);
int * glellipseBounds;
#endif
