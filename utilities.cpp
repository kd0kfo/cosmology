#include "utilities.h"


Utilities::Utilities()
{
    useTimeStamp = false;
    runAsDaemon = false;
    h = 0.73;
    omegaM = 0.27;
}

Utilities::~Utilities()
{

}

int Utilities::sub_main(int argc,char** argv)
{


  
  bool useSM = false;//use supermongo

  int returnMe = 0;

  if(argc == 1)
    {
      stdOutHelp();
      return 0;
    }

  for(int i = 1;i<argc;i++)
    {

      DString bean(argv[i]);
      bean.toLowerCase();
      if(bean.equals("--drawplane"))
	{
	  int gridSize = 0;
	  if(i+3 < argc)
	    {
	      gridSize = (int) Double(argv[i+3]).doubleValue();
	    }

	  Plane<Double> * drawMe = Plane<Double>::readPlane(DString(argv[i+1]));

	  
	  drawPlane(drawMe,argv[i+2],(gridSize != 0),(gridSize == 0) ? 1000 : gridSize);
	  
	  delete drawMe;
	  drawMe = 0;

	  returnMe = 1;
	}
      else if(bean.equals("--shearimage"))
	{

	  Plane<Double> * shearMe = Plane<Double>::readPlane(DString(argv[i+1]));
	  Plane<Double> * convergence = Plane<Double>::readPlane(DString(argv[i+2]));

	  Plane<math::Complex> * shearMap = Plane<math::Complex>::readPlane(DString(argv[i+3]));
	  Plane<Double> * sheared = Utilities::shearImage(shearMe,*convergence,*shearMap);

	  sheared->savePlane(DString(argv[i+4]));
	  
	  delete shearMe;
	  delete convergence;
	  delete shearMap;
	  delete sheared;
	  shearMap = 0;
	  shearMe =  sheared = convergence =  0;
	  return 1;
	  
	}
      else if(bean.equals("--omegam"))
	{
	  omegaM = Double(argv[i+1]).doubleValue();
	}
      else if(bean.equals("--h"))
	{
	  h = Double(argv[i+1]).doubleValue();
	}
      else if(bean.equals("--drawellipse"))
	{
	  int j = i+1;
	  DString fileName(argv[j++]);
	  DEBUG_PRINT(fileName);
	  Double coefficient(argv[j++]);
	  Double exponent(argv[j++]);
	  int WIDTH = (int) Double(argv[j++]).doubleValue();
	  int HEIGHT = (int) Double(argv[j++]).doubleValue();
	  Double A(argv[j++]);
	  Double E(argv[j++]);
	  Double angle(argv[j++]);
	  int X = (int) Double(argv[j++]).doubleValue();
	  int Y = (int) Double(argv[j++]).doubleValue();
	  Double object;

	  drawEllipse(fileName, coefficient.doubleValue(), exponent.doubleValue(), A.doubleValue(), E.doubleValue(),angle.doubleValue(),X, Y, object, WIDTH, HEIGHT);

	  return 0;
	}
      else if(bean.equals("--normalize"))
	{
	  DString plane(argv[i+1]);
	  Double normalizationConstant(argv[i+2]);
	  

	  Plane<Double>  * normMe = Plane<Double>::readPlane(plane);
	  normalize(normMe,normalizationConstant.doubleValue());
	  
	  normMe->savePlane(plane);	  
	}		
      else if(bean.equals("--drawshearmap"))
	{
	  DString plane(argv[i+1]);
	  DString bmp(argv[i+2]);
	  Plane<Double> * shearMap = Plane<Double>::readPlane(plane);
	  drawShearMap(shearMap,bmp);
	  delete [] shearMap; shearMap = 0;
	}
      else if(bean.equals("--redshifttodistance"))
	{
	  Double redshift(argv[i+1]);
	  
	  double distance = Cosmology::redshiftToDistance(0,redshift.doubleValue(),omegaM,h*100);
	  verbosePrint(Double(distance*h).toDString() + " h^-1 Mpc");
	  returnMe = 1;
	}
      else if(bean.equals("--timestamp"))
	{
	  useTimeStamp = true;
 }
      else if(bean.equals("--writeplottabledata"))
	{
	  DString plane(argv[i+1]);
	  DString output(argv[i+2]);

	  Plane<Double> * rewriteMe = Plane<Double>::readPlane(plane);
	  Functions::writePlotableData(rewriteMe, output, useSM);
	  delete rewriteMe;
	  rewriteMe = 0;
	}
      else if(bean.equals("--sm"))
	{
	  useSM = true;
	}
      else if(bean.equals("--help"))
	{
	  stdOutHelp();
	  return 0;
	  
	}

    }

}



bool Utilities::drawPlane(Plane<Double> * planeToDraw,DString fileName,bool useGrid,int gridSpace)
{

  int rows = planeToDraw->numberOfRows();
  int columns = planeToDraw->numberOfColumns();

  Double zero(0.0);
  Double color(255.0);
  Plane<Double> drawMe(rows,columns,zero);

  DEBUG_PRINT("for loop");

  for(int i = 0; i < rows;i++)
    for(int j = 0;j< columns; j++)
      {
	//DEBUG_PRINT(planeToDraw->getValue(i,j));

	if(planeToDraw->getValue(i,j) != 0)
	  {
	    drawMe.setValue(i,j,color);
	  }
      }

  DEBUG_PRINT("done with loop");

  drawMe.draw(fileName, false, useGrid,gridSpace);

  return true;

}

void Utilities::verbosePrint(const DString bean)
{
#if __VERBOSE__
  DString printString = bean;
  DString timeMe = getTime();
  if(useTimeStamp || runAsDaemon)
    printString = timeMe.substring(0,timeMe.length() - 1) + DString(": ") + printString;

  std::cout << printString << std::endl;
#endif
}

DString Utilities::getTime()
{

  time_t timey;

  time ( &timey );

  return DString(ctime (&timey));
  

}  


bool Utilities::drawShearMap(Plane<Double> * planeToDraw, DString fileName, bool useGrid, int gridSpace)
{

  double max = planeToDraw->getMaxValue().doubleValue();

  int rows = planeToDraw->numberOfRows();
  int columns = planeToDraw->numberOfColumns();

  for(int i = 0 ;i< rows;i++)
    for(int j = 0; j< columns;j++)
      {
	Double curr = planeToDraw->getValue(i,j).doubleValue()/max;

	if(curr >= 0.66)
	  {
	    curr = Double((1-curr.doubleValue())*255/0.33,0.0,0.0);
	  }
	else if(curr >= 0.33 && curr < 0.66)
	  {
	    curr = Double(0.0,(0.66-curr.doubleValue())*255/0.33,0.0);
	  }
	else
	  {
	    curr = Double(0.0,0.0,(0.33-curr.doubleValue())*255/0.33);
	  }

	planeToDraw->setValue(i,j,curr);
      }
  planeToDraw->draw(fileName,false,useGrid,gridSpace);
  planeToDraw->savePlane(fileName+".txt");
}

void Utilities::normalize(Plane<Double> * normMe, double normalizationConstant)
{
  int rows = normMe->numberOfRows();
  int columns = normMe->numberOfColumns();

  for(int i = 0; i<rows;i++)
    for(int j=0;j<columns;j++)
      {
	Double curr(normMe->getValue(i,j).doubleValue()/normalizationConstant);
	normMe->setValue(i,j,curr);
      }
}






void Utilities::drawEllipse(DString fileName, double coefficient, double exponent, double semimajorLength, double eccentricity, double semimajorAxisAngle, int centerX, int centerY,Double objectInPlane, int width, int height) const
{
  Plane<Double> * drawMe = drawEllipse( coefficient, exponent,  semimajorLength,  eccentricity, semimajorAxisAngle, centerX, centerY, objectInPlane, width, height);
  drawMe->savePlane(fileName);
  
  delete drawMe;
  drawMe = 0;
}

Plane<Double> * Utilities::drawEllipse(double coefficient, double exponent, double semimajorLength, double eccentricity, double semimajorAxisAngle, int centerX, int centerY,Double objectInPlane, int width, int height) const
{
  bool haveWarned = false;
  Double background(0.0);
  Plane<Double> * mainPlane = new Plane<Double>(width,height,background);
  double a = semimajorLength;
  double e = eccentricity;
  double b = a*(1-e);
  
  double angle = -1*semimajorAxisAngle;//planes are backwards, planning fallacy :(
  
  int leftBound = -1-((int) a)+centerX;
  leftBound = -width/2;
  int rightBound = 1+((int) a)+centerX;
  rightBound = width/2; 
  int lowerBound = -1-((int) b)+centerY;
  lowerBound = -height/2;
  int upperBound = 1+((int) b)+centerY;
  upperBound = height/2;
  
  if(!haveWarned && abs(leftBound) > width/2 || abs(rightBound) > width/2 || abs(lowerBound) > height/2 || abs(upperBound) > height/2)
    {
      //throw DavidException("Ellipse out of Plane Bounds",DavidException::INVALID_ARGUMENT_ERROR_CODE);
      
      DEBUG_PRINT("Ellipse went out of plane bounds");
      haveWarned = true;
    }

  for(int i = leftBound ; i<rightBound;i++)
    for(int j = lowerBound ; j<upperBound;j++)
      {
	double x = i-centerX;
	double y = j-centerY;
	double ellipse = pow(x*cos(angle)+y*sin(angle),2.0)/pow(a,2.0)+pow(y*cos(angle)-x*sin(angle),2.0)/pow(b,2.0);
	if(ellipse <= 1)
	  {
	    double r =  pow(i-centerX,2.0)+pow(j-centerY,2.0);
	    Double curr(coefficient*pow(r,exponent));
	    if(r == 0 && exponent == 0)
	      curr = coefficient;
	    if(r == 0 && exponent < 0)
	      curr = 0.0;
	    
	    
	    Double X(i+width/2);
	    Double Y(j+height/2);
	    if(X.intValue() > 0 && Y.intValue() > 0 && X.intValue()  < mainPlane->numberOfRows() && Y.intValue()  < mainPlane->numberOfColumns()) 
		mainPlane->setValue(X.intValue(),Y.intValue(),curr);

	  }
      }
  return mainPlane;
  
}

Plane<Double> * Utilities::createHistogram(Plane<Double> * plane, int numberOfBins, int scale, int gridSize, int * colorBarDimensions)
{
  
   double * bounds = new double[2];
  bounds[1] = plane->getMaxValue().doubleValue();
  bounds[0] = plane->getMinValue().doubleValue();

  int rows, columns;

  Plane<Double> * drawMe = new Plane<Double>(rows = plane->numberOfRows(),columns = plane->numberOfColumns(),0.0);

  Plane<Double> * colorBar;

  if(colorBarDimensions != 0)
    colorBar = createColorBar(colorBarDimensions[0],colorBarDimensions[1],scale);
  else
    colorBar = 0;

  Double * RGBBins = createColorBins(numberOfBins);

  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {
	double curr = plane->getValue(i,j).doubleValue();
	int bin = getBin(curr, bounds[0], bounds[1], numberOfBins, scale) - 1;
	Double newVal(0.0);
	if(bin >= 0 && bin < numberOfBins)
	  newVal = RGBBins[bin];
	drawMe->setValue(i,j,newVal);
      }

  if(colorBar != 0)
    addColorBar(colorBar,drawMe);

  delete [] RGBBins;
  delete [] bounds;
  delete colorBar;
  colorBar = 0;
  bounds = 0;

  return drawMe;
}

void Utilities::drawHistogram(DString filename, Plane<Double> * plane, int numberOfBins, int scale, int gridSize, int * colorBarDimensions)
{

  Plane<Double> * drawMe = createHistogram(plane, numberOfBins, scale, gridSize, colorBarDimensions);

  if(gridSize <= 0)
    drawMe->draw(filename);
  else
    drawMe->draw(filename,false, true, gridSize);

  delete drawMe;
  drawMe = 0;
}

int Utilities::getBin(double value, double min, double max, int numberOfBins, int scale)
{
  double  bin; 

  double tmp = min;

  switch(scale)
    {
    case SCALE_LOG://log
      {
	value = log(value-min+1);
	min = 0;
	max = log(max-min+1);
	break;
      }
    case SCALE_EXP://exp
      {
	value = exp(value);
	min = exp(min);
	max = exp(max);
	break;
      }
    case SCALE_INVERSE://1/r
      {
	value = 1/value;
	min = 1/max;
	max = 1/tmp;
      }
    default:
      break;
    }

  bin =  (value-min)/(max-min);

  double * intPart = new double;
  modf(bin*numberOfBins,intPart);
  bin = *intPart;

  delete intPart;
  intPart = 0;

  if(bin >= 0)
    return (int) bin;
  else
    return (int) -1*bin;
}


Plane<Double> * Utilities::createColorBar(int rows, int columns, int scale)
{
  Plane<Double> * returnMe = new Plane<Double>(rows,columns,0.0);

  int numberOfBins;
  Double * RGBBins = createColorBins(numberOfBins = 500);
  Double blackPixel(0.0,0.0,0.0);
  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {

	int bin = getBin(i+1.0,1.0,rows, numberOfBins,scale);
	if(bin < numberOfBins && bin >= 0)
	  returnMe->setValue(i,j,RGBBins[bin]);
	else
	  returnMe->setValue(i,j,blackPixel);
      }

  delete [] RGBBins;
  RGBBins = 0;
  return returnMe;
}

void Utilities::addColorBar(Plane<Double> * colorBar, Plane<Double> * plane) throw (DavidException)
{

  if(colorBar->numberOfRows() > plane->numberOfRows() || colorBar->numberOfColumns() > plane->numberOfColumns())
    throw DavidException("The color bar extends beyond the dimensions of the plane", DavidException::PLANE_OUT_OF_BOUNDS_ERROR_CODE);

  for(int i = 0;i<colorBar->numberOfRows();i++)
    for(int j = 0;j<colorBar->numberOfColumns();j++)
      {
	Double curr = colorBar->getValue(i,j);
	plane->setValue(i,j,curr);
      }

}

Double * Utilities::createColorBins(int numberOfBins)
{
  Double * RGBBins = new Double[numberOfBins];

  for(int i = 0;i<numberOfBins;i++)
    {
      double specificBin = ((double) i)/numberOfBins;
      //specificBin *= (double) 2/3;
      unsigned char r,g,b;
      Rainbow::getRainbow(specificBin,r,g,b);
      RGBBins[i] = Double(r,g,b);
    }

  return RGBBins;
}


Plane<Double> * Utilities::shearImage(Plane<Double> * original, Plane<Double>& convergence, Plane<math::Complex>& shearMap)
{
  ShearImage si(convergence, shearMap);

  Plane<Double> * returnMe = si.createDistortedImage(original);
  
  return returnMe;
}

void Utilities::stdOutHelp()
{
  verbosePrint("--drawplane planeName bmpname <gridsize>");
  verbosePrint("Draws the given plane, planeName, as a bitmap, bmpname, with the optional grid of size, gridsize.");
  verbosePrint("--omegam #");
  verbosePrint("Optional flag to manually set a value of Omega(mass).");
  verbosePrint("--h #");
  verbosePrint("Optional flag to manually set a value of the Hubble Parameter.");
  verbosePrint("--redshifttodistance redshift");
  verbosePrint("Gives the distance for the given redshift. Note: this should be the last flag.");
	  
  verbosePrint("--writeplottabledata input output");
  verbosePrint("Converts a plant to data that can be read by IDL");
  verbosePrint("--shearimage <ellipse info> <convergence> <shear> <outputfile>");
  verbosePrint("Takes a plane of ellipse information and shears the ellipses based on the given convergence and shear maps. convergence plane is real and the shear map is a Complex plane");
  verbosePrint("--drawshearmap inputfile outputfile");
  verbosePrint("Makes a Histogram bitmap from a given plane.");
  verbosePrint("--normalize planefile normalization");
  verbosePrint("Normalizes a plane to the given constant.");
  verbosePrint("--timestamp");
  verbosePrint("Adds a timestamp to every line printed.");
  verbosePrint("--drawellipse filename coefficient exponent semimajoraxisangle width height a e x y");
  verbosePrint("Draws an ellipse with values of coefficient*distance^exponent. The ellipse is centered at (x,y) and has a semimajor axis length of a and eccentricity e at an angle in radius w.r.t. the horizontal axis. Width and height are the dimensions of the plane.");
  verbosePrint("--sm");
  verbosePrint("Formats files for use with Supermongo");
}
