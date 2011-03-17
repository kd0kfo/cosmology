#ifndef UTILITIES_CPP
#define UTILITIES_CPP

#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>

#include "libdnstd/DavidException.h"
#include "libdnstd/DString.h"
#include "libdnstd/StringTokenizer.h"
#include "libdnstd/Double.h"

#include "libmygl/plane.h"
#include "libmygl/Cosmology.h"
#include "libmygl/shearimage.h"

#include "Functions.h"
#include "Rainbow.h"

#ifndef __USE_BOINC__
#ifndef MYDAEMON_CPP
#include "mydaemon.cpp"
#define DAEMON_NAME "flatten"
#endif
#endif

#if __DEBUG__
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x) 1;
#endif

class Utilities{

 public:
  Utilities();
  ~Utilities();

/**
 * Acts as an intermediary between int main() and the main processes.
 */
int sub_main(int argc,char** argv);


/**
 * Draws a Double Plane as a bitmap
 *
 * The plane draw may have red grid lines drawn on top of it.
 */
bool drawPlane(Plane<Double> * planeToDraw,DString fileName,bool useGrid = false, int gridSpace = 1000);


/**
 * Draws a Histogram of a Double Plane.
 *
 * @see bool drawPlane
 */
static bool drawShearMap(Plane<Double> * planeToDraw, DString fileName, bool useGrid = false, int gridSpace = 1000);

/**
 * Prints the given String to Standard Output
 */
void verbosePrint(const DString bean);
void verbosePrint(const char * bean){return verbosePrint(DString(bean));}
void verbosePrint(const double d){return verbosePrint(Double(d).toDString());}
void verbosePrint(const int i){return verbosePrint((double) i);}


/**
 * Normalizes the given plane to a given value.
 */
static void normalize(Plane<Double> * normMe, double normalizationConstant);


/**
 * Outputs the given Double plane in a format which can be used in IDL
 * Optional parameters allows for the use of a formatting that is better for Supermongo
 */
void writePlotableData(Plane<Double> * plane, DString fileName, bool useSM = false);

/**
 * Draws an ellipse as a plane ASCII file.
 *
 * Values within the ellipse are proportional to the distance from the center.
 *
 * @param coefficient double proportionality constant for the values within the ellipse
 * @param exponent double value of the exponent of the distance dependence of the values within the ellipse
 * @param semimajorLength double length of the semimajor axis
 * @param eccentricity double value of the eccentricity of the ellipse (e=1-b/a)
 * @param semimajorAxisAngle double angle in radius w.r.t. the "horizontal" axis
 * @param centerX integer x coordinate of the center of the ellipse
 * @param centerY integer y coordinate of the center of the ellipse
 * @param objectInPlane depreicated parameter
 * @param width integer width of the plane
 * @param height integer height of the plane
 *
 */
 void drawEllipse(DString fileName, double coefficient, double exponent, double semimajorLength, double eccentricity, double semimajorAxisAngle, int centerX, int centerY, Double objectInPlane, int width, int height) const;

/**
 * Creates a plane with ellipses
 */
 Plane<Double> * drawEllipse(double coefficient, double exponent, double semimajorLength, double eccentricity, double semimajorAxisAngle, int centerX, int centerY,Double objectInPlane, int width, int height) const;

/**
 * Creates a crude histogram.
 *
 * @param filename DString
 * @param plane Plane<Double>*
 * @param numberOfBins int optional parameter to set the number of bins (default = 10)
 */
 static  Plane<Double> * createHistogram(Plane<Double> * plane, int numberOfBins, int scale, int gridSize, int * colorBarDimensions);
 
 static void drawHistogram(DString filename, Plane<Double> * plane, int numberOfBins = 10, int scale = 0, int gridSize = -1, int * colorBarDimensions = 0);


 /**
  * Draws a plane of the given dimensions containing a rainbow with the given scale. This can be used as a lengend for a histogram.
  *
  * @return Plane<Double>*
  */
 static Plane<Double> * createColorBar(int rows, int column, int scale);

 /**
  * Adds the given colorbar to the plane.
  * Note: The colobar must be smaller than or equal in dimension of the plane.
  * Otherwise an exception is thrown.
  * 
  * @param colorBar Plane<Double>*
  * @param plane Plane<Double>*
  * @throw DavidException
  */
 static void addColorBar(Plane<Double> * colorBar, Plane<Double> * plane) throw (DavidException);

 /**
  * Creates a Double array containing the rainbow divided spread out into the
  * given number of bins.
  *
  * @param numberOfBins int
  * @return Double* size = numberOfBins
  */
 static Double * createColorBins(int numberOfBins);

 /**
  * Createds a distorted image from shear and convergence maps
  *
  */
 static Plane<Double> * shearImage(Plane<Double> * original, Plane<Double>& convergence, Plane<math::Complex>& shearMap);

/**
 * Returns the current timeas a DString
 */
 static DString getTime();

 /**
  * Sends the help and usage information to stdout
  */
 void stdOutHelp();

 static const int SCALE_LOG = 1;///<natural log scale
 static const int SCALE_EXP = 2;///<Exponential (e) scale
 static const int SCALE_INVERSE = 3;///<Inverse (1/r) scale

 private:
 
 /**
  * Returns 3-dim int array of the color for the bin of value.
  * The 3 elements are a value from 0 to 255, representing RGB respectively.
  *
  * Scale factors:
  * 0 = no scaling
  * 1 = log
  * 2 = exponential
  *
  * @param value double current value
  * @param max double maximum value
  * @param min double minimum value
  * @param numberOfBins int number of bins
  * @param scale int
  * @return int[3]
  */
 static int getBin(double value, double min, double max, int numberOfBins,int scale);

 bool useTimeStamp,runAsDaemon;
 double h;//< Hubble Constant Parameter (default 0.7)
 double omegaM;//< Omega Mass Value
 
};//end Utilities

#endif
