/**
 * 
 * This file is part of flatten, a program that takes a 3-D particle
 * distribution and produces a two-dimentional projection.
 *
 * Copyright 2007, 2010 David Coss, PhD
 *
 * flatten is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * flatten is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with flatten.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef FLATTEN_CPP
#define FLATTEN_CPP

#include <deque>

#include "libdnstd/Double.h"

#include "libmygl/plane.h"

class Flatten{

 public:

  Flatten(   bool _useTimeStamp = false,
	     bool _scaleImage = false,
	     bool _drawImage = false,
	     int _lineOfSight = 2,
	     int _timeInSeconds = 0);

  ~Flatten();


  /**
   * Parse Command Line Arguments
   *
   * @param argc int number of arguments
   * @param argv char* array of arguments argv[0] is the executable name.
   */
  int parseArgs(int argc, char** argv);

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
   * particle line of sight locations
   *
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

  bool useTimeStamp,scaleImage,drawImage;
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
