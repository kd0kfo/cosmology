#ifndef FLATTEN_CPP
#define FLATTEN_CPP

#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include "DavidException.h"
#include "DString.h"
#include "StringTokenizer.h"
#include "Double.h"
#include "planecreator.h"
#include "DStack.h"

#ifndef __USE_BOINC__
#include "EasyBMP.h"
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

void verbosePrint(const DString bean);
void verbosePrint(const char * bean){return verbosePrint(DString(bean));}
void verbosePrint(const double d){return verbosePrint(Double(d).toDString());}
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
 * Returns the current time as a DString.
 */
DString getTime();

/**
 * Reads the cluster file.
 * The Cluster file is read and pointer to an array which contains
 * particle line of sight locations (INSERT UNITS HERE)
 * @param fileName const char* name of the file to be read.
 * @param xWidth int increment of width of x values in the Plane
 * @param yWidth int increment of width of y values in the Plane
 * @return Plane<Double> pointer.
 */
Plane<utils::DStack<Double> > * parseClusterFile(const char* fileName, Double  resolution);
Plane<utils::DStack<Double> > * parseClusterFile(DString fileName, Double resolution){return parseClusterFile(fileName.getString(),resolution);}


void drawDStackPlane(Plane<utils::DStack<Double> > * planeToPrint,const char * fileName, int width, int height);

void writeDStackPlane(const char * fileName, Plane<utils::DStack<Double> > * plane);

/**
 * Prints the help information
 */
void stdout_help();

int createSurfaceMassDensity(char* fileName, Double resolution, char* outLENS);

bool runAsDaemon,useTimeStamp,scaleImage,drawImage;
DString fileID;
 double numberOfParticles,timeInSeconds, totalMass;
double lineOfSight;//0=x,1=y,2=z
double * xBounds;
double * yBounds;
DString clusterFilename,clusterOutBMP,clusterOutLENS;
Double clusterResolution;
 double unitConversion;
};

#endif
