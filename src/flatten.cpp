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
#include "flatten.h"
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
#endif

Flatten::Flatten(bool _useTimeStamp,
	    bool _scaleImage,
	    bool _drawImage,
	    int _lineOfSight,
	    int _timeInSeconds)
{
  useTimeStamp = _useTimeStamp;
  scaleImage = _scaleImage;
  drawImage = _drawImage;
  lineOfSight = _lineOfSight;
  timeInSeconds = _timeInSeconds;
  
  totalMass = 0;
  numberOfParticles = 0;
  xBounds = yBounds = 0;
  unitConversion = 1;

}

Flatten::~Flatten()
{
  delete [] xBounds;
  delete [] yBounds;
  xBounds = yBounds = 0;
}

int Flatten::parseArgs(int argc, char** argv)
{
  
  if(argc == 1)
    {
      stdout_help();
      return 0;
    }
  
	bool parsedParams = false;
	std::string bean;
	for(int i = 1;i<argc;i++)
	{
	  bean = utils::lower_case(argv[i]);
		int j = i;
		if(bean == "--createsurfacemassdensity")
		  {
		    clusterFilename = argv[++j];
		    std::string resDString(argv[++j]);
		    clusterResolution = Double(resDString);
		    clusterOutLENS = argv[++j];
		  }
		else if(bean == "--lineofsight")
		  {
		    lineOfSight = (int) Double(argv[i+1]).doubleValue();
		  }
		else if(bean == "--findbounds")
		  {

		    std::string fileName = argv[i+1];
		    Double searchColumn(argv[i+2]);
		    Double startingPoint(argv[i+3]);
		    
		    double * bounds = findBounds(fileName.c_str(),searchColumn.doubleValue(),startingPoint.doubleValue());

		    printf("Maximum Value is %g  located at line %d\n",bounds[2],(size_t)bounds[0]);
		    printf("Minimum Value is %g  located at line %d\n",bounds[3],(size_t)bounds[1]);
		    delete [] bounds;
		    return 0;
		  }		    
		else if(bean == "--scaledimage")
		  {
		    scaleImage = true;
		  }
		else if(bean == "--drawimage")
		  {
		    drawImage = true;
		    clusterOutBMP = argv[i+1];
		  }
		else if(bean == "--unit")
		  {
		    unitConversion = Double(argv[i+1]).doubleValue();
		  }
		else if(bean == "--timeevolve")
		  {
		    timeInSeconds = Double(argv[i+1]).doubleValue();
		    timeInSeconds *= 365*24*3600;
		    DEBUG_PRINT("time is: ");
		    DEBUG_PRINT(timeInSeconds);
		  }
		else if(bean == "--help")
		{
		  stdout_help();
		  return 0;
		}

	}

	return -42;

}

void Flatten::stdout_help()
{
  printf("Options are:\n");
  printf("--createsurfacemassdensity <input> <majoraxislength> <2dplanename>\n");
  printf("Creates 2-D massdensity (double values) based on the 3-D box of particles <input>. <majoraxislength> represents the number of pixels in the longer of x and y. Then the cluster is then written as an ascii data <outLENS>.\n");
  printf("--findbounds <filename> <column> <startingline>\n");
  printf("Gives the Maximum and minimum values and the lines they are at. Note: Counting is ZERO INDEX\n");
  printf("--timeevolve <time>\n");
  printf("Evolves the cluster data by the given number of years, based on the cluster particle velocities\n");
  printf("--timestamp\n");
  printf("Places a timestamp at the beginning of every STDOUT message.\n");
  printf("--lineofsight\n");
  printf("Sets the line of sight. X = 0, Y = 1, Z = 2 (default)\n");
  printf("--unit <unit factor>\n");
  printf("Multiplies results by <unit factor>\n");
  printf("--daemon\n");
  printf("Runs the process as a daemon\n");
}

std::string Flatten::getTime()
{
  time_t timey;
  time ( &timey );
  return ctime (&timey);
}  


double*  Flatten::findBounds(const char * fileName,int columnToSearch,int lineToStartAt)
{

  using namespace std;
  using utils::StringTokenizer;

  fstream readMe(fileName,ios::in);

  char buffer[256];
  bool keepGoing = true;

  double * returnMe = new double[4];
  returnMe[2] = returnMe[3] = returnMe[1] = returnMe[0] = 0;

  int currentRow = -1;

  while(readMe.getline(buffer,256))
    {
      if(buffer[0] == '#')
	continue;

      currentRow++;
      
      if(currentRow >= lineToStartAt)
	{
	  std::string temp(buffer);
	  
	  StringTokenizer tokie(temp," ");
	  
	  for(int i = 0;i<columnToSearch;i++)
	    tokie.nextToken();
	  

	  try
	    {
	      Double value(tokie.nextToken());
	      if(currentRow == lineToStartAt)
		{
		  returnMe[2] = returnMe[3] = value.doubleValue();
		  returnMe[0] = returnMe[1] = currentRow;
		}

	      if(value.doubleValue() > returnMe[2])
		{
		  returnMe[2] = value.doubleValue();
		  returnMe[0] = currentRow;
		}
	      
	      if(value.doubleValue() < returnMe[3])
		{
		  returnMe[3] = value.doubleValue();
		  returnMe[1] = currentRow;
		}
	      
	    }
	  catch(DavidException de)
	    {
	      if(de.getCode() == DavidException::STRING_TOKENIZER_ERROR_CODE)
		{
		  DEBUG_PRINT(std::string("Column ") + Double(columnToSearch).str() + std::string(" does not exist"));
		}
	      de.stdErr();
	      throw de;
	    }
	}
    }

  return returnMe;
  
}

Plane<std::deque<Double> > * Flatten::parseClusterFile(const char* fileName, Double resolution)
{
  using namespace std;
  using utils::StringTokenizer;
  using std::deque;

  double xInc, yInc;

  delete [] xBounds;
  delete [] yBounds;

  //xBounds and yBounds are the highest and lowest values of x and y
  //coordinates, in the 2-D plane. So x and y mean, 2-D coordinates 
  //not 1st and 2nd column in 3-D particle data
  if(xBounds == NULL)
    {
      if(lineOfSight == 2)
	xBounds = findBounds(fileName,0,1);
      else if(lineOfSight == 1)
	xBounds = findBounds(fileName,0,1);
      else
	xBounds = findBounds(fileName,1,1);
    }
  if(yBounds == NULL)
    {
      if(lineOfSight == 2)
	yBounds = findBounds(fileName,1,1);
      else if(lineOfSight == 1)
	yBounds = findBounds(fileName,2,1);
      else
	yBounds = findBounds(fileName,2,1);
    }

  //The larger dimension is picked out, x=1,y=2 (2-D coordinates)
  //then that dimension is used to determine the pixel size
  double pixelSize = max(yBounds[2]-yBounds[3],xBounds[2]-xBounds[3])/resolution.doubleValue();
  
  //Find center of the distribution
  //Offset the smaller dimension to center the grid
  double xCenter = (xBounds[2]-xBounds[3])/2+xBounds[3];
  double yCenter = (yBounds[2]-yBounds[3])/2+yBounds[3];
  double xOffset = 0;double yOffset = 0;
  if(xBounds[2]-xBounds[3] > yBounds[2]-yBounds[3])
    yOffset = resolution.doubleValue()*pixelSize/2 + yBounds[3] - yCenter;
  else
    xOffset = resolution.doubleValue()*pixelSize/2 + xBounds[3] - xCenter;
  
  //open cluster data
  fstream readMe(fileName,ios::in);
  int particleCounter = 0;
  int outOfBox = 0;
  if(!readMe.is_open())
    throw DavidException(std::string("Could not open: ") + fileName,DavidException::IO_ERROR_CODE);

  
  char buffer[256];
  bool keepGoing = true;
  
  printf("Reading %s\n", fileName);

  readMe.getline(buffer,256);
  printf("%s\n",buffer);
  StringTokenizer tokie(buffer," ");
  printf("Tokenizer started\n");
  fileID = tokie.nextToken();//fileID

  totalMass = Double(tokie.nextToken()).doubleValue();//total mass of cluster
  printf("Total mass: %g Solar Masses\n",totalMass);
  tokie.nextToken();//virial radius

  tokie.nextToken();//
  tokie.nextToken();// CENTER OF MASS in *Mpc*
  tokie.nextToken();//
  
  numberOfParticles = Double(tokie.nextToken()).doubleValue();//total particle count
  printf("There are %f total particles\n",numberOfParticles);

  deque<Double> stacky;//particles in a column along the line of sight are stored in a stack. Then this stack is popped to sum the total particles, with weight if desired.

  Plane<deque<Double> > * returnMe = new Plane<deque<Double> >(resolution.toInt(),resolution.toInt(),stacky);
  std::string newHeader = std::string("#total mass = ") + Double(totalMass).str() +std::string(" solar masses. Dimensions ") + Double(resolution*pixelSize).str() + std::string("x") + Double(resolution*pixelSize).str() + " kpc";
  printf("Parsing file: %s\n",fileName);
  
  
  while(readMe.getline(buffer,256))
    {
      if(buffer[0] == '#')
	continue;

      tokie = StringTokenizer(buffer, " ");

      Double X(tokie.nextToken());//
      Double Y(tokie.nextToken());//position in *kpc*
      Double Z(tokie.nextToken());//

      //move the particles if a time other than 0 was given 
      //or if no time was given at all
      Double Vx,Vy,Vz;
      Vx = Vy = Vz = 0.0;
      if(tokie.hasMoreTokens())
	Vx = Double(tokie.nextToken());
      if(tokie.hasMoreTokens())
	Vy = Double(tokie.nextToken());
      if(tokie.hasMoreTokens())
	Vz = Double(tokie.nextToken());

      if(timeInSeconds != 0)
	{
	  X += Vx*timeInSeconds*unitConversion;
	  Y += Vy*timeInSeconds*unitConversion;
	  Z += Vz*timeInSeconds*unitConversion;
	  
	}
 
      //assign the chosen spacial coordinate for the 2-D project coordinates
      if(lineOfSight == 0)
	{
	  Double temp = Z;
	  Z = X;
	  X = Y;
	  Y = temp;

	  temp = Vz;
	  Vz = Vx;
	  Vx = Vy;
	  Vy = temp;

	}
      else if(lineOfSight == 1)
	{
	  Double temp = Z;
	  //Z = Y;
	  //Y = X;
	  //X = temp;
	  Z = Y;
	  Y = temp;

	  temp = Vz;
	  Vz = Vy;
	  Vy = temp;

	}
      
      //calculate pixel location of particle in the grid
      //Find distance from center of distribution. Calculate distance from
      //edge of grid in pixels.
      int x = Double((X.doubleValue()-xCenter + xOffset)/pixelSize+resolution.doubleValue()/2 ).toInt();
      int y = Double((Y.doubleValue()-yCenter + yOffset)/pixelSize+resolution.doubleValue()/2 ).toInt();
      particleCounter++;

      //If particle is outside of grid, warn. Otherwise, add it to the grid stack
      if(x >= returnMe->numberOfRows() || y >= returnMe->numberOfColumns() || x < 0 || y < 0)
	{
	  //if the particle moves out of the box, we're going to ignore it(for now)
	  outOfBox++;
	  printf("Particle moved out of box\n");
	}
      else
	{
	  deque<Double> previousStack = returnMe->getValue(x,y);
	  previousStack.push_back(Z.doubleValue());
	  returnMe->setValue(x,y,previousStack);
	}
    }  

  printf("%d by %d plane created: %d total particles used with %d particles out of the box\n",
	 returnMe->numberOfRows(), returnMe->numberOfColumns(),particleCounter,outOfBox);

  newHeader += std::string("\n#Made with " ) + Double(particleCounter - outOfBox).str() + " total particles";
  
  //Set Header
  if(newHeader != "")
    returnMe->setHeader(newHeader);
  
  return returnMe;
}

void Flatten::drawdequePlane(Plane<std::deque<Double> > * planeToPrint,const char * fileName, int width, int height)
{

#ifndef __USE_BOINC__
  BMP image;

  image.SetSize(width,height);
  image.SetBitDepth(24);

  using std::deque;

  for(int i = 0;i< width - 1;i++)
    for(int j =0;j<height - 1;j++)
      {
	deque<Double> stacky = planeToPrint->getValue(i,j);
	
	if(stacky.size())
	  {
	    image(i+1,j+1)->Red = 255;
	    image(i+1,j+1)->Green = 0;
	    image(i+1,j+1)->Blue = 0;
	  }
	else
	  {
	    image(i+1,j+1)->Red = 0;
	    image(i+1,j+1)->Green = 0;
	    image(i+1,j+1)->Blue = 0;
	  }
      }

  image.WriteToFile(fileName);
  printf("%s has been written\n",fileName);
  #endif
}

void Flatten::writedequePlane(const char * fileName, Plane<std::deque<Double> > * plane)
{

  using std::deque;
  int width = plane->numberOfRows();
  int height = plane->numberOfColumns();

  using namespace std;
  fstream outfile(fileName,ios::out);

  if(!outfile.is_open()){
    throw DavidException(std::string("Could not open: ")+fileName,DavidException::IO_ERROR_CODE);
  }

  outfile << plane->getHeader() << std::endl;
  outfile << width << std::endl;
  outfile << height << std::endl;

  double normalize = totalMass/numberOfParticles;
  double curr = 0;
  deque<Double> stacky;
  int counter = 0;
  for(int i = 0;i < width; i++)
    for(int j = 0;j< height;j++)
      {
	
	stacky = plane->getValue(i,j);
	curr = 0;

	while(stacky.size())
	  {
	    stacky.pop_back();
	    curr += normalize;
	  }

	outfile << curr << std::endl;
      }
  outfile << width << std::endl;
  outfile << height << std::endl;

  outfile.close();
}

int Flatten::createSurfaceMassDensity(const char* fileName, Double resolution, const char* outLENS)
{
		    using std::deque;
		    
		    Plane<std::deque<Double> > * cluster = parseClusterFile(fileName,resolution);
		    
#ifndef __USE_BOINC__
		    if(drawImage)
		      drawdequePlane(cluster,clusterOutBMP.c_str(),cluster->numberOfRows(),cluster->numberOfColumns());
#endif
		    printf("writing plane\n");
		    writedequePlane(outLENS,cluster);

		    delete cluster;
		    delete [] xBounds;
		    delete [] yBounds;
		    xBounds = yBounds = 0;		
      
		    return 0;

}
