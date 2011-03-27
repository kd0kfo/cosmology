#ifndef RAY_TRACE_CPP
#define RAY_TRACE_CPP

#if __VERBOSE__ 
#define VERBOSE_PRINT(x) std::cout << x << std::endl; 
#else 
#define VERBOSE_PRINT(x)
#endif

#ifdef __DEBUG__
#define DEBUG_PRINT(x) std::cout << x << std::endl; 
#else 
#define DEBUG_PRINT
#endif

#include <iostream>
#include <fstream>
#include <math.h>

#ifdef _WIN32
#include <time.h>
#endif

#include "libdnstd/utils.h"
#include "libdnstd/Double.h"
#include "libdnstd/DRandom.h"
#include "libdnstd/Complex.h"
#include "libdnstd/StringTokenizer.h"

#include "libmygl/plane.h"
#include "libmygl/glellipse.h"
#include "libmygl/densityprofile.h"
#include "libmygl/planecreator.h"

#ifndef _WIN32
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#ifndef __USE_BOINC__
#include "mydaemon.cpp"
#define DAEMON_NAME "ray_trace_ellipse"
#endif

template class Plane<math::Complex>;


/** \mainpage Documentation of the Main Executables.
 * Created By David Coss, 2007
 */
template <class T> void buildSource(Plane<T> *, double * sourceParams,T);///< constructs the source plane.
template <class T> void constructLens(Plane<T> * source, Plane<T> * lens, T valueOfSource, T lensValue);///< constructs the lens plane.
bool writeSourceInfo(double ** sourceParams, int specificSource, int numberOfSources, utils::DRandom * randy, bool useRandom, int numberOfPixels);///< writes the source info(eg location) to a text file.

/**
 * Main start method.
 * Before this method is called arguments are parsed and then sent to sub_main
 */
int sub_main(int argc, char** argv);

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
 * @param params General Parameters.
 * @param lensParams Lens Parameters.
 * @param sourceParams Array of Parameters for each source
 * @param numberOfSources.
 * @return bool true if writing was successful.
 */
bool writeParameters(double * params,double * lensParams,double ** sourceParams, int numberOfSources);

/**
 * The main simulation routine.
 * This is the meat of the program. The actually does the ray tracing.
 *
 * @param lens Plane<Double> representing the lens plane
 * @param sources Plane<Double> containing the source background images
 * @param massDensity DensityProfile containing the 2D mass density
 * @param numberOfSources int number of source planes to be lensed
 * @param specificSource int source to be run if only one source is used
 * @param specificLens int specific lens to be used if only one is used
 * @param parameter1 int option parameter (as of now, not used give it anything)
 * @return int zero if successful
 * @throw DavidException Exception thrown upon error (if you're lucky)
 */
int simulation(Plane<Double> * lens, Plane<Double> * sources, DensityProfile * massDensity, int numberOfSources, int specificSource, int specificLens, int parameter1) throw (DavidException);

/**
 * Creates parameters for a lens to be used in the simulation
 *
 * @param lensParams double array of parameters to use in the simulation
 * @param specificLens int specific lens to use
 * @return bool true (always, don't ask why)
 */
bool createLensParams(double * lensParams, int specificLens);///< Sets the lens parameters.
int simulationSetup();

/**
 * Parses arguments sent from command line.
 */
int parseArgs(int argc, char** argv);

/**
 * Parses the Parameter std::string Array.
 * @param filename
 * @return std::string*
 * @see loadParameters(double * params,double * lensParams,double ** sourceParams,int numParams, int numLensParams, int numberOfSources)
 */
std::string * paramParser(const char * fileName);//Size will be set in method. size is the size of the std::string vector

/**
 * Reads the Parameters from the Parameters file.
 * @param array of general parameters.
 * @param lensParams Lens parameters.
 * @param sourceParams Array of Parameters for each Source.
 * @param numParams Number of General Lens parameters.
 * @param numLensParams Number of Lens parameters.
 * @param numberOfSources Number of Sources.
 */
void loadParameters(double * params,double * lensParams,double ** sourceParams,int numParams, int numLensParams, int numberOfSources);

/**
 * Creates a 2 Dimensional Mass Array from lens parameters.
 */
void createSurfaceMassDensity(std::string parameters, std::string fileName);

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


std::string fileNamePrefix;///< prefix of files created in the simulation
utils::DRandom * randy;///< Random number Generator.
bool includeCricalCurveAndCaustic;///< Whether Critical Curves and Caustics should be used.
Plane<math::Complex> * deflectionPlane;
char * lensMassDeflectionPlane;
char * sourceBMPFilename;
bool runExistingDeflection;
bool createDeflection;
bool runSim;
std::string mainPrefix;
bool runAsDaemon;
int gridSpace;
int overallParamNumber;
std::string * parameterArray;
std::string * savesourcelocations;
Double bgColor;
void drawEllipse(const char * parameters);
Plane<Double> * lensMassDensity;
void verbosePrint(const char * string);
void verbosePrint(std::string string){return verbosePrint(string.c_str());}
void verbosePrint(double val){return verbosePrint(Double(val).str().c_str());}
int * glellipseBounds;
bool useTimeStamp;
bool drawRemovedArea;
double * offset;//Offset of the computational grid (in pixels)
#endif
