#ifndef FLATTEN_CPP
#define FLATTEN_CPP

#include <iostream>
#include <fstream>
#include <time.h>
#include <cmath>
#include <deque>

#include "libdnstd/DavidException.h"
#include "libdnstd/StringTokenizer.h"
#include "libdnstd/Double.h"
#include "libdnstd/utils.h"

#include "libmygl/planecreator.h"

#ifndef __USE_BOINC__
#include "libmygl/EasyBMP/EasyBMP.h"
#include "mydaemon.cpp"
#define DAEMON_NAME "flatten"
#endif

#if __DEBUG__
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x) 1;
#endif

class Flatten{

 public:

 Flatten(bool _runAsDaemon = false,
	    bool _useTimeStamp = false,
	    bool _scaleImage = false,
	    bool _drawImage = false,
	    int _lineOfSight = 2,
	   int _timeInSeconds = 0);

 ~Flatten();


/**
 * Parse Command Line Arguments
 * @param argc int number of arguments
 * @param argv char* array of arguments argv[0] is the executable name.
 */
int parseArgs(int argc, char** argv);

 void verbosePrint(const std::string bean);
void verbosePrint(const double d){return verbosePrint(Double(d).str());}
void verbosePrint(const int i){return verbosePrint((double) i);}

/**
 * Finds the lines with the largest and smalles value of the specified column.
 * The values are stored in an array where the [0] value is the largest line and the smallest is in [1] in the array.
 * The values [2] and [3] are the highest and lowest values respectively
 * If not column is specified the left most is used.
 * Numbering starts from the top and left and begins with zero. So the second column and the second row would both be labeled with a 1.
 *
 * @param fileName conat char * file to read.
 * @param columnToSearch int number of the column to search through.
 * @param lineToStartAt int line to start searching at. Default: top most column.
 * @returns Line number (first line = 0).
 */
double *  findBounds(const char * fileName,int columnToSearch = 0,int lineToStartAt = 0);

/**
 * Returns the current time as a std::string.
 */
std::string getTime();

/**
 * Reads the cluster file.
 * The Cluster file is read and pointer to an array which contains
 * particle line of sight locations (INSERT UNITS HERE)
 * @param fileName const char* name of the file to be read.
 * @param xWidth int increment of width of x values in the Plane
 * @param yWidth int increment of width of y values in the Plane
 * @return Plane<Double> pointer.
 */
Plane<std::deque<Double> > * parseClusterFile(const char* fileName, Double  resolution);


void drawdequePlane(Plane<std::deque<Double> > * planeToPrint,const char * fileName, int width, int height);

void writedequePlane(const char * fileName, Plane<std::deque<Double> > * plane);

/**
 * Prints the help information
 */
void stdout_help();

int createSurfaceMassDensity(const char* fileName, Double resolution, const char* outLENS);

bool runAsDaemon,useTimeStamp,scaleImage,drawImage;
std::string fileID;
 double numberOfParticles,timeInSeconds, totalMass;
double lineOfSight;//0=x,1=y,2=z
double * xBounds;
double * yBounds;
std::string clusterFilename,clusterOutBMP,clusterOutLENS;
Double clusterResolution;
 double unitConversion;
};

#endif
