/**
 * Created By David Coss, 2007
 */
#include "ray_trace_ellipse.h"


#ifndef __USE_BOINC__
int main(int argc, char** argv)
{
  int retval = 0;
  try
    {
#ifdef USE_MPI
      retval = utils::init_mpi(&argc,&argv,&mpi_data.rank,&mpi_data.num_ranks);
      if(retval)
	{
	  char msg[MPI_MAX_ERROR_STRING];
	  int size;
	  MPI_Error_string(retval,msg,&size);
	  std::cerr << "ERROR - MPI initialization failed." << std::endl
		    << "Reason: " << msg << " Error Code " << retval;
	}
      else
	{
	  mpi_data.hostname = utils::get_hostname();
	  std::cout << "MPI started on " << mpi_data.hostname 
		    << " with rank " << mpi_data.rank << std::endl;
	}
#endif
      if(retval == 0)
	retval = super_main(argc,argv);
    }
  catch(DavidException de)
    {
#ifdef USE_MPI
      std::cerr << "Process of Rank " << mpi_data.rank  << "on host " << mpi_data.hostname << std::endl;
#endif
      std::cerr << "*************" << std::endl << "Program Ended due to " << de.getCause() << std::endl;
      std::cerr << "The following message was provided:" << std::endl;
      de.stdErr();

      retval = de.getCode();
    }

#ifdef USE_MPI
  MPI_Finalize();
#endif

  return retval;
}
#endif
//endif boinc



int super_main(int argc, char** argv)
{
  runAsDaemon = false;

  mainPrefix = "run-";

  sourceBMPFilename = "source.bmp";

  bgColor = 0.0;

  savesourcelocations = 0;
	
  lensMassDensity = 0;

  lensMassDeflectionPlane = "lens";

  offset = 0;
	
  gridSpace = -1000;

  drawRemovedArea = false;

  glellipseBounds = new int[4];
  for(int i = 0;i<4;i++)
    glellipseBounds[i] = -1;

  useTimeStamp = false;
  int retval = parseArgs(argc,argv);
  if( retval != -42)
    {
      return retval;
    }

  if(!runAsDaemon)
    {
      retval = sub_main(argc, argv);
      return retval;
    }

  int returnValue = 0;
	
#ifndef __USE_BOINC__	
  //Running as daemon if you get this far
  // Initialize the logging interface 
  openlog( DAEMON_NAME, LOG_PID, LOG_LOCAL5 );
  syslog( LOG_INFO, "starting" );
	
  // One may wish to process command line arguments here 
  // Daemonize 
  daemonize( "/var/lock/subsys/" DAEMON_NAME );
	
  // Now we are a daemon -- do the work for which we were paid 
	
  returnValue = sub_main(argc, argv);
	
  // Finish up 
  syslog( LOG_NOTICE, "terminated" );
  closelog();
#endif
      
  return returnValue;
}

int sub_main(int argc, char** argv)
{
  int returnMe = 0;

  try{
    std::string startMessage = "ray_trace_ellipse started";
    if(runAsDaemon)
      startMessage += " as daemon process.";
    verbosePrint(startMessage);
    returnMe = simulationSetup();
  }
  catch(DavidException de)
    {

#ifdef __DEBUG__
      DEBUG_PRINT(std::string("Type: ") + de.getType());
      DEBUG_PRINT(std::string("Code: ") + Double(de.getCode()).str());
      DEBUG_PRINT(std::string("Message: "));
#endif


      de.stdOut();

      return de.getCode();

    }

  delete savesourcelocations;
  delete [] parameterArray;

  return returnMe;
}

int simulationSetup()
{
  using utils::DRandom;
  DRandom * randy = 0;

  int N,height,width,numberOfSources,numberOfLensVariations;


  if(parameterArray == 0)
    {
      N = 480;
      height = N;
      width = N;
      numberOfSources = 1;
      numberOfLensVariations = 1;
	    
    }
  else
    {

      N = (int) Double(parameterArray[3]).doubleValue();
      height = N;
      width = N;
      numberOfSources = (int) Double(parameterArray[2]).doubleValue();
      numberOfLensVariations = 1;
    }

#ifdef USE_MPI
  utils::mpi_adjust_glellipsebounds(::glellipseBounds,N);
  DEBUG_PRINT("Process " << mpi_data.rank << ": Making plane dim: " << glellipseBounds[0] << ", " << glellipseBounds[1] << ", " << glellipseBounds[2] << ", " << glellipseBounds[3]);

#endif

  Double zeroDouble(0.0);
  Plane<Double> * lens = new Plane<Double>(width,height, bgColor);

  Plane<Double> * sources = 0;
#ifndef __USE_BOINC__
  if(access(sourceBMPFilename.c_str(),R_OK) == 0)
    sources = Plane<Double>::bmpToPlane(sourceBMPFilename);
#endif
  Plane<Double> * lensMassPlane = new Plane<Double>(width,height, zeroDouble);
  Plane<math::Complex> * deflectionPlane = 0;

  DensityProfile * massDensity = 0;

  //Run once per source
  for(int i = 0; i<numberOfSources;i++)
    {
      using std::string;
      verbosePrint(string("Sim ") + Double(i+1).str() + string(" of ") + Double(numberOfSources).str() + ": ");

      fileNamePrefix = mainPrefix + Double((double) i).str();
      std::string fullLensFilename = fileNamePrefix  + "lensedimage.bmp";
      simulation(lens,sources, massDensity, 1,i,3,i);


#ifndef __USE_BOINC__
      lens->draw(fullLensFilename,false,gridSpace > 0,gridSpace);
		
      if(sources != NULL)
	{
	  std::string fullSourceFilename = fileNamePrefix + "sources.bmp";
	  DEBUG_PRINT("Drawing sources");
	  DEBUG_PRINT(sources);
	  sources->draw(fullSourceFilename,false,gridSpace > 0,gridSpace);
	}
#endif
      DEBUG_PRINT("DELETING");
      delete lens;
      DEBUG_PRINT("deleted lens");
      delete sources;
      DEBUG_PRINT("delted sources");
      Double cSucksALot(0.0);
      lens = new Plane<Double>(width,height, bgColor);
      sources = new Plane<Double>(width,height, bgColor);
    }


  //making a combined image if multiple sources are used.
  if(numberOfSources > 1)
    {
      verbosePrint("Combined Sim");
      fileNamePrefix = mainPrefix + Double(numberOfSources).str();
      simulation(lens,sources, massDensity, numberOfSources,-1,0,numberOfSources);
#ifndef __USE_BOINC__
      lens->draw(fileNamePrefix+"lensedimage.bmp",false,gridSpace < 0,gridSpace);
      sources->draw(fileNamePrefix+"sources.bmp",false,gridSpace < 0,gridSpace);
#endif
    }



  DEBUG_PRINT("DELETING for the last time");
  delete randy;
  delete lens;
  DEBUG_PRINT("BUT NOT HERE");
  delete sources;
  DEBUG_PRINT("BUT NOT HERE Either");
  if(massDensity != 0)
    {
      DEBUG_PRINT("Clearing all fields");
      massDensity->clearAllFields();
      DEBUG_PRINT("all fields cleared");
      delete massDensity;
      DEBUG_PRINT("mass density deleted");
    }
	
  massDensity = 0;
  randy = 0;
  lens = sources = lensMassPlane = 0;
	
  return 0;
}



int simulation(Plane<Double> * lens, Plane<Double> * sources, DensityProfile * massDensity, int numberOfSources, int specificSource, int specificLens, int parameter1) throw (DavidException)
{
	
	
  using namespace std;
  using utils::DRandom;
	
  int N = lens->numberOfRows();//Plane is always square
	
  includeCricalCurveAndCaustic = false;
  bool useRandom = false;//source placement
  if(useRandom)
    randy = new DRandom();
	
  int numberOfObsParameters = 6;
  double * lensParams = new double[10];
  double * params = new double[numberOfObsParameters];
  double ** sourceParams = new double*[numberOfSources];

  if(parameterArray == 0)
    {
      params[0] = N; //N number of pixels, area = 
      params[1] = .38; //arc length of each pixel in arcseconds
      params[2] = 6.67*pow((double) 10,(double) -8);//G, newtonian grav. constant
      params[3] = 2*pow((double) 10,(double) 33);//Solar mass, grams
      params[4] = 3.09*pow((double) 10, (double) 24);//distance scale Mpc
      params[5] = 3*pow(10.0,10.0);//speed of light in cm/s
		
      createLensParams(lensParams,specificLens);

      writeSourceInfo(sourceParams,specificSource, numberOfSources, randy, useRandom, N);
      //UNCOMMENT ABOVE LINE AFTER TEST!!!
    }
  else
    {
      DEBUG_PRINT("Loading parameters");
      loadParameters(params,lensParams,sourceParams,numberOfObsParameters,10,numberOfSources);

    }
  delete massDensity;//No memory leaks here!

  DEBUG_PRINT("Creating Mass density");
  massDensity = new DensityProfile(lensMassDensity, lensParams, params);
  GLAlgorithm * gls = new GLAlgorithm[numberOfSources];

  for(int i = 0;i<numberOfSources;i++)
    {
      DEBUG_PRINT(lensMassDeflectionPlane);	    
      if(runExistingDeflection)
	{
	  // USE EXISTING DEFLECTION PLANE
	  std::string deflectionFilename(lensMassDeflectionPlane);
	  
	  verbosePrint(std::string("Parsing Deflection Plane: ")+deflectionFilename);
	  deflectionPlane = Plane<math::Complex>::readPlane(deflectionFilename);
			
	  if(deflectionPlane == 0)
	    DEBUG_PRINT("deflectionPlane is null");
	  gls[i] = GLAlgorithm(deflectionPlane,params,sourceParams[0],sources);
	}
      if(createDeflection)
	{

	  // CREATE NEW DEFLECTION PLANE 
	  verbosePrint(std::string("Creating Deflection Plane: ") + lensMassDeflectionPlane);
	  gls[i] = GLAlgorithm(massDensity,sourceParams[0],sources,glellipseBounds,offset);

#ifdef USE_MPI
		  MPI_Comm lens_group;
		  int has_lens = (gls[i].getDeflectionPlane() != 0) ? 1 : 0;
		  MPI_Comm_split(MPI_COMM_WORLD,has_lens,mpi_data.rank,&lens_group);
#endif

	  if(gls[i].getDeflectionPlane() != 0)
	    {
#ifdef USE_MPI
	      utils::mpi_recombine(gls[i],lens_group);
	      if(mpi_data.rank == utils::MASTER_RANK)
	    	  gls[i].getDeflectionPlane()->savePlane(lensMassDeflectionPlane);

#else
	      gls[i].getDeflectionPlane()->savePlane(lensMassDeflectionPlane);
#endif
	    }
#ifdef USE_MPI
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
      if(!createDeflection && !runExistingDeflection)
	{
	  verbosePrint("Running based on parameters");
	  gls[i] = GLAlgorithm(deflectionPlane,params,sourceParams[0],sources);
	}
		  
    }

#ifndef __USE_BOINC__
  verbosePrint("Running Sim");
  if(runSim)//Run ray trace, else no sim
    {
      Plane<math::Complex> * sourceLocations = 0;
      if(savesourcelocations != 0)
	sourceLocations = new Plane<math::Complex>(N,N,math::Complex(0.0,0.0));


      
      int percentFinished = 0;
      Double sourceValue;
      for(int i = 0;i<N;i++)
	{
	  for(int j = 0;j<N;j++)
	    {
	      for(int k = 0;k<numberOfSources;k++)
		{
		  try
		    {
		      sourceValue = lens->getValue(i,j) + gls[k].mapsToSource(i,j);
		      if(sourceLocations != 0 && savesourcelocations != 0)
			{
			  double * loc = gls[k].newLocation(i,j);
			  math::Complex curr = math::Complex(loc[0],loc[1]);
			  sourceLocations->setValue(i,j,curr);
			  delete [] loc;
			}
		      lens->setValue(i,j,sourceValue);
		    }
		  catch(DavidException de)
		    {
		      DEBUG_PRINT(de.getMessage());
		    }
		}//end for(k...)
	    }
	  if((i*100/N) >= (percentFinished+5))
	    {
	      percentFinished = (int) i*100/N;
	      verbosePrint("Percent finished: ");
	      VERBOSE_PRINT(percentFinished);
	    }
	}

      verbosePrint("Simulation Complete");

      if(savesourcelocations != 0 && sourceLocations != 0)
	{
	  verbosePrint(std::string("Writing source plane to ") + *savesourcelocations);
#ifdef USE_MPI
	  utils::mpi_recombine(sourceLocations,MPI_COMM_WORLD);
#endif
	  sourceLocations->savePlane(*savesourcelocations);
	}
    }//end if(runsim)
  else
    {
      DEBUG_PRINT("No sim");
    }
#endif

	
#ifndef __USE_BOINC__
  writeParameters(params,lensParams,sourceParams,numberOfSources);
#endif

  DEBUG_PRINT("MADE IT THIS FAR!");
  delete [] gls;

  for(int i = 0; i< numberOfSources;i++)
    delete [] sourceParams[i];
  delete [] sourceParams;

  gls = 0;
  params = 0;
  lensParams = 0;//don't delete here. Contained in massDensity which is kept as a pointer.
  sourceParams = 0;

  return 0;
}/**/


template <class T> void buildSource(Plane<T> * source, double * sourceParams, T newValue)
{
  
  double xCenter = sourceParams[0];
  double yCenter = sourceParams[1];
  double a = sourceParams[3];
  double e = sourceParams[4];
	
  for(int i = 0; i<source->numberOfColumns();i++)
    {
      for(int j = 0;j<source->numberOfRows();j++)
	{
	  double x = i - (source->numberOfColumns()/2);
	  double y = j - (source->numberOfRows()/2);
	  double test = pow( x-xCenter,2.0)/pow(a,2.0)+pow( y-yCenter,2)/pow(a*(1-e),2);
	  if(test <= 1)
	    source->setValue(x+(source->numberOfColumns()/2),y+(source->numberOfRows()/2),newValue);
	}
    }


}


template <class T> void constructLens(Plane<T> * source,Plane<T> * lens,T valueOfSource,T lensValue)
{

  int rows = lens->numberOfRows();
  int columns = lens->numberOfColumns();

  for(int i = 0;i<columns;i++)
    {
      for(int j = 0; j<rows; j++)
	{
	  if(source->getValue(i,j) == valueOfSource)
	    {
	      double x = i - (columns/2);
	      double y = j - (rows/2);
	      double l = sqrt(x*x + y*y);
	      double angle = (x != 0) ? atan((y/x)) : 3.14159/2;//note trig functions use radians
	      double * newl = new double[2];
	      newl[0] = .5*(l+sqrt(l*l+4));
	      newl[1] = .5*(l-sqrt(l*l+4));

	      for(int q = 0;q<2;q++)
		{
		  double newX = (newl[q])*cos(angle);
		  double newY = (newl[q])*sin(angle);
		  if(newX < rows/2 && newX > -1*rows/2 && newY < rows/2 && newY > -1*columns/2)
		    lens->setValue(newX+(rows/2),newY+(columns/2),lensValue);
		}
	      delete [] newl;
	      newl = 0;
	    }
	}
    }
	



}


//If a specific sources is wanted specify it in specificSource, else specificSource < 0
bool writeSourceInfo(double ** sourceParams, int specificSource, int numberOfSources, utils::DRandom * randy, bool useRandom, int numberOfPixels)
{
  std::string fullFileName = fileNamePrefix+"sourceinfo.txt";
  std::ofstream sourceLOC(fullFileName.c_str());

  int N = numberOfPixels;

  using namespace std;

  for(int i = 0;i<numberOfSources;i++)
    {
      sourceParams[i] = new double[7];
				
      //sourceParams[i][0] x coordinate of center
      //sourceParams[i][1] y coordinate of center
      sourceParams[i][2] = 2000;//source redshift
      sourceParams[i][3] = 10;//semimajor axis length in arcseconds
      sourceParams[i][4] = 0;//source ellipticity, e=1-b/a
      sourceParams[i][5] = 0;//orientation of source semimajor axis in degrees
      sourceParams[i][6] = 255;//Source flux
				
		
      if(useRandom)
	{
	  throw DavidException("FIX RANDOM SOURCE GENERATOR: SORRY ABOUT THIS :(");
	  /*
	    sourceParams[i][0] = N*randy->random()/200;
	    if((randy->random()/3) % 2 == 0)
	    sourceParams[i][0] *= -1;
	    sourceParams[i][1] = N*randy->random()/200;
	    if((randy->random()/3) % 2 == 0)
	    sourceParams[i][1] *= -1;
	    if(sourceLOC.is_open())
	    sourceLOC << sourceParams[i][0] << ", " << sourceParams[i][1] << std::endl;
	    /**/
	}
      else if(specificSource >= 0)
	{
	  switch(specificSource)
	    {
	    case 1:
	      sourceParams[i][0] = 50;//x coordinate of source center;
	      sourceParams[i][1] = 50;//y coordinate of source center;
	      break;
	    case 2:
	      sourceParams[i][0] = 100;//x coordinate of source center;
	      sourceParams[i][1] = 100;//y coordinate of source center;
	      break;
	    case 3:
	      sourceParams[i][0] = 150;//x coordinate of source center;
	      sourceParams[i][1] = 150;//y coordinate of source center;
	      break;
	    case 4:
	      sourceParams[i][0] = 200;//x coordinate of source center;
	      sourceParams[i][1] = 200;//y coordinate of source center;
	      break;
	    case 5:
	      sourceParams[i][0] = -68;//x coordinate of source center;
	      sourceParams[i][1] = -150;//y coordinate of source center;
	      break;
	    default:
	      sourceParams[i][0] = 0;//x coordinate of source center;
	      sourceParams[i][1] = 0;//y coordinate of source center;
	      break;
	    }
	  if(sourceLOC.is_open())
	    sourceLOC << sourceParams[i][0] << ", " << sourceParams[i][1] << std::endl;
	}
      else//For multiple, non-random sources
	{
	  switch(i)
	    {
	    case 1:
	      sourceParams[i][0] = 50;//x coordinate of source center;
	      sourceParams[i][1] = 50;//y coordinate of source center;
	      break;
	    case 2:
	      sourceParams[i][0] = 100;//x coordinate of source center;
	      sourceParams[i][1] = 100;//y coordinate of source center;
	      break;
	    case 3:
	      sourceParams[i][0] = 150;//x coordinate of source center;
	      sourceParams[i][1] = 150;//y coordinate of source center;
	      break;
	    case 4:
	      sourceParams[i][0] = 200;//x coordinate of source center;
	      sourceParams[i][1] = 200;//y coordinate of source center;
	      break;
	    case 5:
	      sourceParams[i][0] = -68;//x coordinate of source center;
	      sourceParams[i][1] = -150;//y coordinate of source center;
	      break;
	    default:
	      sourceParams[i][0] = 0;//x coordinate of source center;
	      sourceParams[i][1] = 0;//y coordinate of source center;
	      break;
	    }
	  if(sourceLOC.is_open())
	    sourceLOC << sourceParams[i][0] << ", " << sourceParams[i][1] << std::endl;
	}

		

    }
  sourceLOC.close();

  return true;
}

/**
 * writeParameters writes the parameters to 2 files: "human readable" file and a file which can be parsed by ray_trace
 */
bool writeParameters(double * params,double * lensParams,double ** sourceParams, int numberOfSources)
{
  std::string fullFileName = fileNamePrefix+"human-parameters.txt";
  std::string fullFileName2 = fileNamePrefix+"parameters.txt";
  std::ofstream outfile(fullFileName.c_str());
  std::ofstream outfile2(fullFileName2.c_str());
  using namespace std;

  int numParams = 6;
  int numLensParams = 10;
	
  if(outfile.is_open() && outfile2.is_open())
    {

      outfile2 << numParams << ";number of parameters" << endl;
      outfile2 << numLensParams <<";number of lens parameters" << endl;
      outfile2 << numberOfSources <<";Number of Sources" << endl;

      outfile << "General Parameters:" << endl;
      for(int i = 0;i<numParams;i++)
	{	
	  switch(i)
	    {
	    case 0:
	      outfile << "Length and Width in pixels: " << params[i] << endl;
	      outfile2 << params[i] << ";Length and Width in pixels " << endl;
	      break;
	    case 1:
	      outfile << "Arc length per pixel: " << params[i] << endl;
	      outfile2 << params[i] << ";Arc length per pixel" <<  endl;
	      break;
	    case 2:
	      outfile << "Value of Newton's Constant: " << params[i] << endl;
	      outfile2 <<  params[i] << ";Value of Newton's Constant" << endl;
	      break;
	    case 3:
	      outfile << "Solar mass: " << params[i] << endl;
	      outfile2 << params[i] << ";Solar mass" << endl;
	      break;
	    case 4:
	      outfile << "Distance Scale: " << params[i] << endl;
	      outfile2 << params[i] << ";Distance Scale" << endl;
	      break;
	    case 5:
	      outfile << "Speed of Light" << params[i] << endl;
	      outfile2 << params[i] << ";Speed of Light" << endl;
	    default:
	      break;
	    }
	}
      outfile << "Critical Curves included: ";
      if(includeCricalCurveAndCaustic)
	{
	  outfile << "Yes" << endl;
	  outfile2 << 1 ;
	}
      else
	{
	  outfile << "No" << endl;
	  outfile2 << 0 ;
	}
      outfile << endl;
      outfile2 << ";Critical Curves included" << endl;
    }
  else
    return false;

  if(outfile.is_open())
    {
      outfile << "Lens Parameters:" << endl;
      for(int i = 0;i<numLensParams;i++)
	{	
	  switch(i)
	    {
	    case 0:
	      outfile << "x coordinate of lens center: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";x coordinate of lens center"  << endl;
	      break;
	    case 1:
	      outfile << "y coordinate of lens center: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";y coordinate of lens center"  << endl;
	      break;
	    case 2:
	      outfile << "lens redshift or distance: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";lens redshift or distance" << lensParams[i] << endl;
	      break;
	    case 5:
	      outfile << "lens velocity of dispersion: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";lens velocity of dispersion" << lensParams[i] << endl;
	      break;
	    case 6:
	      outfile << "lens mass density (in solar masses per square arcsecond): " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";lens mass density (in solar masses per square arcsecond)" << endl;
	      break;
	    case 7:
	      outfile << "Shear: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";Shear" << endl;
	      break;
	    case 8:
	      outfile << "Semi-major axis length: " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";Semi-major axis length" << endl;
	      break;
	    case 9:
	      outfile << "Eccentricity(e=1-b/a): " << lensParams[i] << endl;
	      outfile2 << lensParams[i] << ";Eccentricity(e=1-b/a)" << endl;
	      break;
	    default:
	      outfile2 << "42;useless array parameter PLEASE use me next" << endl;
	      break;
	    }
	}
      outfile << endl;
    }
  else
    return false;

  if(outfile.is_open())
    {
      outfile << "Source Parameters:" << endl;
      for(int sourceNumber = 0;sourceNumber<numberOfSources;sourceNumber++)
	{	
	  outfile << "Source " << sourceNumber + 1 << ": " << endl;
	  for(int i = 0;i<7;i++)
	    {
	      switch(i)
		{
		case 0:
		  outfile << "x coordinate of source center: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";x coordinate of source center" << endl;
		  break;
		case 1:
		  outfile << "y coordinate of source center: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";y coordinate of source center" << endl;
		  break;
		case 2:
		  outfile << "source redshift or distance: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";source redshift or distance" << endl;
		  break;
		case 3:
		  outfile << "semimajor axis length in arcseconds: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";semimajor axis length in arcseconds" << endl;
		  break;
		case 4:
		  outfile << "source eccentricity: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";source eccentricity" << endl;
		  break;
		case 5:
		  outfile << "orientation of source semimajor axis in degrees: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";orientation of source semimajor axis in degrees" << endl;
		  break;
		case 6:
		  outfile << "Source flux: " << sourceParams[sourceNumber][i] << endl;
		  outfile2 << sourceParams[sourceNumber][i] << ";Source flux" << endl;
		  break;
		default:
		  break;
		}
	    }
	  outfile << endl;
	}
      outfile << endl;
    }
  else
    return false;

  outfile.close();

  return true;

}


bool createLensParams(double * lensParams, int specificLens)
{
  lensParams[0] = 0;//x coordinate of lens center
  lensParams[1] = 0;//y coordinate of lens center
  lensParams[2] = 1000;//lens distance in Mpc (typically 1000)
  lensParams[3] = 600;//lens velocity of dispersion
  lensParams[4] = 2.2*pow(10.0,7.0);//dispersion velocity in glnonsingular in cm/s
  lensParams[5] = 2.2*pow(10.0,7.0);//dispersion velocity in glSIS in cm/s
  lensParams[6] = pow((double) 10, (double) 16); //in Solar Masses, density if applicable, originally 10^12
  //lensParams[6] = 0.000285883;//Saturn in solar masses
  lensParams[7] = 0;//Shear in glSISwithShear
  lensParams[8] = .38;//a, semi-major axis length of lens, in arcseconds
  lensParams[9] = 0;//e, ellipticity e= 1-b/a, b is the semi-minor axis length

  return true;
}


int parseArgs(int argc, char** argv)
{
	
  for(int i = 0;i<argc;i++)
    {
      DEBUG_PRINT(argv[i]);
    }

  runExistingDeflection = true;
  createDeflection = !(runExistingDeflection);
  runSim = true;

  bool parsedParams = false;

  for(int i = 0;i<argc;i++)
    {
      std::string bean = utils::lower_case(argv[i]);
      if(bean == "--newlens")
	{
	  DEBUG_PRINT("newlens");
	  lensMassDeflectionPlane = argv[i+1];
	  runExistingDeflection = false;
	  createDeflection = !(runExistingDeflection);
	  runSim = false;
	}
      else if(bean == "--sourcebmp")
	{
	  sourceBMPFilename = argv[i+1];
	  verbosePrint(argv[i+1]);
	  verbosePrint("source:");
	  verbosePrint(sourceBMPFilename);
	
	}
      else if(bean == "--bgcolor")
	{
	  VERBOSE_PRINT("--bgcolor has been disabled");
	  //for(int j = 0;j<3;j++)
	  //  bgColor.setValue(j,Double(argv[i+j+1]).doubleValue());
	}
      else if(bean == "--offset")
	{
	  if(offset != 0)
	    delete [] offset;
	  
	  offset = new double[2];
	  offset[0] = Double(argv[i+1]).doubleValue();
	  offset[1] = Double(argv[i+2]).doubleValue();
	}
      else if(bean == "--prefix")
	{
	  mainPrefix = argv[i+1];
	}
      else if(bean == "--deflectionmap" || bean == "--lens" )
	{
	  lensMassDeflectionPlane = argv[i+1];
	  runExistingDeflection = true;
	  createDeflection = !(runExistingDeflection);
	  runSim = true;
	}
      else if(bean == "--forcerun")
	{
	  runSim = true;
	}
      else if(bean == "--stoprun")
	{
	  runSim = false;
	}
      else if(bean == "--usegrid")
	{
	  gridSpace = (int) Double(argv[i+1]).doubleValue();
	  DEBUG_PRINT("grid size");
	  DEBUG_PRINT(gridSpace);
	}
      else if(bean == "--drawsource")
	{
	  drawEllipse(argv[i+1]);
	  return 0;
	}
      else if(bean == "--parameters" || bean == "--parameter")
	{
	  parameterArray = paramParser(argv[i+1]);
	  parsedParams = true;
	}
      else if(bean == "--createsurfacemassdensity")
	{
	  std::string parameters = argv[i+1];
	  std::string fileName = argv[i+2];
	  createSurfaceMassDensity(parameters,fileName);

	  return 0;
	}
      else if(bean == "--usesurfacemassdensity")
	{
	  std::string fileName = argv[i+1];
	  lensMassDensity = Plane<Double>::readPlane(fileName);
	}
      else if(bean == "--daemon")
	{
	  runAsDaemon = true;
	}
      else if(bean == "--timestamp")
	{
	  useTimeStamp = true;
	}
      else if(bean == "--addplanes")
	{
	  if(argc < i+4)
	    throw DavidException("I need to two planes to add and a file name for the sum.","Insufficient Arguments",DavidException::INVALID_ARGUMENT_ERROR_CODE);

	  std::string left,right,out;
	  left = argv[i+1];
	  right = argv[i+2];
	  out = argv[i+3];
	  using utils::StringTokenizer;
	  using math::Complex;
	  Complex zero(0,0);
	  Plane<Complex> * leftPlane = new Plane<Complex>(0,0,zero);
	  Plane<Complex> * rightPlane = new Plane<Complex>(0,0,zero);
	  leftPlane = Plane<Complex>::readPlane(left);
	  rightPlane = Plane<Complex>::readPlane(right);
	  verbosePrint(std::string("Adding ")+left + std::string(" to ") + right);
	  DEBUG_PRINT(leftPlane);
	  Plane<Complex> * outPlane = Plane<Complex>::addPlanes(leftPlane,rightPlane);
		    
	  verbosePrint(std::string("Saving sum as ")+out);
	  outPlane->savePlane(out);

	  delete leftPlane;
	  delete rightPlane;
	  delete outPlane;
	  leftPlane = rightPlane = outPlane = 0;

	  return 0;
	}
      else if(bean == "--subtractplanes")
	{
	  using utils::StringTokenizer;
	  using math::Complex;
	  Complex zero(0,0);
	  Plane<Complex> * leftPlane = new Plane<Complex>(0,0,zero);
	  Plane<Complex> * rightPlane = new Plane<Complex>(0,0,zero);
	  std::string left,right,out;
		    
	  leftPlane = Plane<Complex>::readPlane(left = argv[i+1]);
	  rightPlane = Plane<Complex>::readPlane(right = argv[i+2]);
		    
	  verbosePrint(std::string("Subtracting ") + left + std::string(" and ") + right);
	  Plane<Complex> * outPlane = Plane<Complex>::subtractPlanes(leftPlane,rightPlane);

	  verbosePrint(std::string("Saving difference as ")+(out = argv[i+3]));
	  outPlane->savePlane(out);

	  delete leftPlane;
	  delete rightPlane;
	  delete outPlane;
	  leftPlane = rightPlane = outPlane = 0;
		    
	  return 0;
	}
      else if(bean == "--glellipsebounds" || bean == "--glellipse")
	{
	  using utils::StringTokenizer;
	  StringTokenizer tokie = StringTokenizer(argv[i+1],",");

	  for(int i = 0;i<4;i++)
	    glellipseBounds[i] = (int) Double(tokie.nextToken()).doubleValue();
			  
	}
      else if(bean == "--savesourcelocations")
	{
	  savesourcelocations = new std::string(argv[i+1]);
	}
      else if(bean == "--drawremovedarea")
	{
	  drawRemovedArea = true;
	}
      else if(bean == "--squareplane")
	{
	  squarePlane(argv[i+1]);
	  return 0;
	}
      else if(bean == "--help")
	{
#ifndef __USE_BOINC__
	  verbosePrint("Options are:");
	  verbosePrint("--newLens <lens>");
	  verbosePrint("Creates a new lens saved as <lens>");
	  verbosePrint("--createsurfacemassdensity <parameterfile> <outfile>");
	  verbosePrint("Creates 2-D massdensity (double values) based on the given parameters.");
	  verbosePrint("--usesurfacemassdensity <infile>");
	  verbosePrint("Uses an already created 2-D massdensity (double values).");

	  verbosePrint("--sourceBMP <source.bmp>");
	  verbosePrint("Uses <source.bmp> as the source plane");
	  verbosePrint("--deflectionMap <deflectionMap.txt>");
	  verbosePrint("Uses <deflectionMap.txt> as the deflection plane");
	  verbosePrint("--lens <deflectionMap.txt>");
	  verbosePrint("Uses <deflectionMap.txt> as the deflection plane. Same as --deflectionMap  (typing lens was a habit)");
	  verbosePrint("--forceRun");
	  verbosePrint("Forces the Simulation to run, whether or not a deflection map was created or read");
	  verbosePrint("--stopRun");
	  verbosePrint("Forces the Simulation NOT to run, whether or not a deflection map was created or read");
	  verbosePrint("--prefix");
	  verbosePrint("Sets the output file(s) prefix ");
	  verbosePrint("--useGrid (size)");
	  verbosePrint("Draws a grid on output images. Size is the pixel spacing of the grid.");
	  verbosePrint("--bgcolor R G B");
	  verbosePrint("Sets the background color of the images to the RGB value provided");
	  verbosePrint("--offset x y");
	  verbosePrint("Set a pixel offset for the calculation.");
	  verbosePrint("--parameters <file name>");
	  verbosePrint("Parse file provided file for parameters.");
	  verbosePrint("--drawsource a,e,width,height,x,y,R,G,B,filename");
	  verbosePrint("Draws an RGB source of eccentricity e, semimajor axis length a, at (x,y), saved as width by height sized a bitmap called filename. Program exits when completed.");
	  verbosePrint("--glellipsebounds left,right,upper,lower");
	  verbosePrint("Sets bounds on calculation of deflection grid. This allows embarrisingly parallel processing.");
	  verbosePrint("--glellipse left,right,upper,lower");
	  verbosePrint("Same as glellipsebounds. Added to keep typos from killing processing time.");
	  verbosePrint("--timestamp");
	  verbosePrint("Places a timestamp at the beginning of every STDOUT message.");
	  verbosePrint("--addplanes leftplane rightplane outplane");
	  verbosePrint("Adds leftplane and rightplane and saves the sum as outplane");
	  verbosePrint("--subtractplanes leftplane rightplane outplane");
	  verbosePrint("Subtracts leftplane and rightplane and saves the difference as outplane");
	  verbosePrint("--squareplane filename");
	  verbosePrint("Takes the plane, filename, and makes it a square by buffing the shorter dimension with zeros");
#endif			

			
	  return 0;
	}


    }

  if(!parsedParams)
    parameterArray = 0;

  return -42;

}

/**
 * This Method Parses Parameters from a text file and returns them as
 * an entry in a std::string Array. The first three entries are the Number of General Parameters,
 * Number of Lens Parameters and Number of Sources, respectively. So the total array size is
 * (Number of General Parameters) + (Number of Lens Parameters) + 7*(Number of Sources)
 */
std::string * paramParser(const char * fileName)
{
  std::ifstream infile(fileName);
  using namespace std;
		
  int stringSize = 150;
  int totalParams, numParams, numLensParams, sourceNumber;
	
  std::string * paramsVector;
  int vectorCount = 0;
  size_t npos = std::string::npos;

  if(infile.is_open())
    {
      Double tmp;

      char * curr = new char[stringSize];
      infile.getline(curr,stringSize);
      std::string currDString(curr);
      if(currDString.find(";") != npos)
		currDString = currDString.substr(0,currDString.find(";"));
      tmp = Double::parseString(currDString);
      numParams = (int) tmp.doubleValue();

      curr = new char[stringSize];
      infile.getline(curr,stringSize);
      currDString = curr;
      if(currDString.find(';') != npos)
	currDString = currDString.substr(0,currDString.find(";"));
      numLensParams = (int) Double::parseString(currDString).doubleValue();

      curr = new char[stringSize];
      infile.getline(curr,stringSize);
      currDString = curr;
      if(currDString.find(';') != npos)
	currDString = currDString.substr(0,currDString.find(";"));
      sourceNumber = (int) Double::parseString(currDString).doubleValue();
		
      totalParams = numParams + numLensParams + sourceNumber*7 + 3 +1;//+3 to include the total number of parameters(eg numParams) +1 for critical curve inclusion (0 = false, 1 = true)
		
      paramsVector = new std::string[totalParams];
      paramsVector[vectorCount++] = Double((double) numParams).str();
      paramsVector[vectorCount++] = Double((double) numLensParams).str();
      paramsVector[vectorCount++] = Double((double) sourceNumber).str();
		
      for(int i = 0;i< totalParams - 3;i++)
	{
	  curr = new char[stringSize];
	  infile.getline(curr,stringSize);
	  currDString = curr;
	  if(currDString.find(';') != npos)
	    currDString = currDString.substr(0,currDString.find(";"));
	  paramsVector[vectorCount++] = currDString;
	}

    }

	
  overallParamNumber = totalParams;
  return paramsVector;
}

void loadParameters(double * params,double * lensParams,double ** sourceParams,int numParams, int numLensParams, int numberOfSources)
{
  int vectorCount = 3;//This tells the parser how much to skip in the parameter file header which contains the number of parameters in each category.

  for(int i = 0;i<numParams;i++)
    params[i] = Double(parameterArray[vectorCount++]).doubleValue();


  //Ignoreing critical curve stuff for now
  vectorCount++;

  //"Lens Parameters:"
  for(int i = 0;i<numLensParams;i++)
    {			
      
      lensParams[i] =Double(parameterArray[vectorCount++]).doubleValue();
    }
  //"Source Parameters:"
  for(int sourceNumber = 0;sourceNumber<numberOfSources;sourceNumber++)
    {	
      sourceParams[sourceNumber] = new double[7];
      for(int i = 0;i<7;i++)
	{
	  sourceParams[sourceNumber][i] = Double(parameterArray[vectorCount++]).doubleValue();
	  
	}
      
    }
  
  //Sets the glellipseBounds to the whole array if the bounds are not yet set.
		
		
  if(glellipseBounds[1] < 0)
    glellipseBounds[1] = params[0];
  if(glellipseBounds[3] < 0)
    glellipseBounds[3] = params[0];
  if(glellipseBounds[0] < 0)
    glellipseBounds[0] = 0;
  if(glellipseBounds[2] < 0)
    glellipseBounds[2] = 0;
  
}



void drawEllipse(const char * parameters)
{
  using utils::StringTokenizer;
  using utils::PlaneCreator;

  std::string params(parameters);

  std::string a,e,width,height,x,y,R,G,B,filename;

  StringTokenizer tokie(params,",");

  a = tokie.nextToken();
  e = tokie.nextToken();
  width = tokie.nextToken();
  height = tokie.nextToken();
  x = tokie.nextToken();
  y = tokie.nextToken();
  R = tokie.nextToken();
  G = tokie.nextToken();
  B = tokie.nextToken();
  filename = tokie.nextToken();

  int r = (int) Double(R).doubleValue();
  int g = (int) Double(G).doubleValue();
  int b = (int) Double(B).doubleValue();
  double A = Double(a).doubleValue();
  double E = Double(e).doubleValue();
  int WIDTH = (int) Double(width).doubleValue();
  int HEIGHT = (int) Double(height).doubleValue();
  int X = (int) Double(x).doubleValue();
  int Y = (int) Double(y).doubleValue();
  Double object(r,g,b);
  Double background(0,0,0);

  PlaneCreator<Double> * newPlane = new PlaneCreator<Double>(WIDTH, HEIGHT, background);
  newPlane->addEllipse(A, E, X, Y, object);
  Plane<Double> * planizzle = newPlane->getPlane();
	
  verbosePrint(std::string("Drawing ")+filename);
  planizzle->draw(filename);
	
  delete newPlane;
  delete planizzle;

}

void createSurfaceMassDensity(std::string parameters, std::string fileName)
{
  parameterArray = paramParser(parameters.c_str());

  int numberOfSources = 1;
  double * params = new double[6];
  double * lensParams = new double[10];
  double ** sourceParams = new double*[numberOfSources];
  loadParameters(params,lensParams,sourceParams,6,10,numberOfSources);
  
  Plane<Double> * dummy = new Plane<Double>(params[0],params[0],42);
  DensityProfile * profile = new DensityProfile(dummy,lensParams,params);

  profile->drawPlane();
  profile->getPlane()->savePlane(fileName);

  delete [] params;
  delete [] lensParams;
  for(int i = 0 ; i<numberOfSources;i++)
    delete [] sourceParams[i];

  delete [] sourceParams;

  delete profile;

}

std::string getTime()
{

  time_t timey;

  time ( &timey );

  return std::string(ctime (&timey));
  

}


void verbosePrint(const char * string)
{

#if __VERBOSE__
  std::string printString(string);

  if(runAsDaemon || useTimeStamp)
    {
      std::string timeMe = getTime();
      printString = timeMe.substr(0,timeMe.length() - 1) + std::string(": ") + printString;
    }

  std::cout << printString << std::endl;
#endif

}


void squarePlane(const char * fileName)
{
  
  using math::Complex;
  
  Complex cZero(0.0,0.0);

  Plane<Complex> * oldPlane = new Plane<Complex>(0,0,cZero);
  oldPlane = Plane<Complex>::readPlane(fileName);

  int width,height;
  DEBUG_PRINT("here");
  DEBUG_PRINT(oldPlane);

  if(oldPlane->numberOfRows() > oldPlane->numberOfColumns())
    width = height = oldPlane->numberOfRows();
  else
    width = height = oldPlane->numberOfColumns();

  Plane<math::Complex> * newPlane = new Plane<math::Complex>(width,height,cZero);



  for(int i = 0;i < oldPlane->numberOfRows();i++)
    for(int j = 0; j < oldPlane->numberOfColumns();j++)
      {
	Complex curr = oldPlane->getValue(i,j);
	newPlane->setValue(i,j,curr);
      }

  newPlane->savePlane(fileName);

  delete oldPlane;
  delete newPlane;
}
