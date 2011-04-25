/**
 * Created By David Coss, 2007
 */
#include "ray_trace_ellipse.h"

enum option_keys{ CREATE_MASS = 1,SOURCE_BMP,DEFLECTION_MAP,FORCE_RUN,
		  STOP_RUN, FILE_PREFIX, USE_GRID, RBG_COLOR, XOFFSET, YOFFSET,
		  SUBTRACT_PLANES, ADD_PLANES, SQUARE_PLANE,
		  SAVESOURCE_LOC,DRAW_REMOVED_AREA};

static const struct option ray_trace_options[] = 
  {
    {"newlens", required_argument, NULL,'n'},
    {"createsurfacemassdensity",required_argument,NULL,CREATE_MASS},
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

#ifndef __USE_BOINC__
int main(int argc, char** argv)
{
  int retval = 0;
  glellipseBounds = NULL;
  parameterArray = NULL;
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
  struct ray_trace_arguments args;

  default_arguments(&args);


  glellipseBounds = new int[4];
  for(int i = 0;i<4;i++)
    glellipseBounds[i] = -1;

  int retval = parseArgs(argc,argv,&args);
  
  if(args.makeMassDensity.size() != 0)
    {
      createSurfaceMassDensity(args.makeMassDensity);
      if(args.offset != NULL)
	delete [] args.offset;
      return 0;
    }

  if( retval != -42)
    {
      return retval;
    }

  if(!args.runAsDaemon)
    {
      retval = sub_main(&args);
      if(args.offset != NULL)
	delete [] args.offset;
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
	
  returnValue = sub_main(&args);
	
  // Finish up 
  syslog( LOG_NOTICE, "terminated" );
  closelog();
#endif

  if(args.offset != NULL)
    delete [] args.offset;
  
  return returnValue;
}

int sub_main(struct ray_trace_arguments *args)
{
  int returnMe = 0;

  try{
    std::string startMessage = "ray_trace_ellipse started";
    if(args->runAsDaemon)
      startMessage += " as daemon process.";
    verbosePrint(args,startMessage);
    returnMe = simulationSetup(args);
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

  delete [] parameterArray;

  return returnMe;
}

int simulationSetup(struct ray_trace_arguments *args)
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
  utils::mpi_adjust_glellipsebounds(glellipseBounds,N);
  DEBUG_PRINT("Process " << mpi_data.rank << ": Making plane dim: " << glellipseBounds[0] << ", " << glellipseBounds[1] << ", " << glellipseBounds[2] << ", " << glellipseBounds[3]);

#endif

  Double zeroDouble(0.0);
  Plane<Double> * lens = new Plane<Double>(width,height, bgColor);

  Plane<Double> * sources = 0;
#ifndef __USE_BOINC__
  if(access(args->sourceBMPFilename.c_str(),R_OK) == 0)
    sources = Plane<Double>::bmpToPlane(args->sourceBMPFilename);
#endif
  Plane<Double> * lensMassPlane = new Plane<Double>(width,height, zeroDouble);
  Plane<math::Complex> * deflectionPlane = 0;

  DensityProfile * massDensity = 0;

  //Run once per source
  for(int i = 0; i<numberOfSources;i++)
    {
      using std::string;
      verbosePrint(args,string("Sim ") + Double(i+1).str() + string(" of ") + Double(numberOfSources).str() + ": ");

      args->fileNamePrefix = args->mainPrefix + Double((double) i).str();
      std::string fullLensFilename = args->fileNamePrefix  + "lensedimage.bmp";
      simulation(args,lens,sources, massDensity, 1,i,3);


#ifndef __USE_BOINC__
      lens->draw(fullLensFilename,false,args->gridSpace > 0,args->gridSpace);
		
      if(sources != NULL)
	{
	  std::string fullSourceFilename = args->fileNamePrefix + "sources.bmp";
	  DEBUG_PRINT("Drawing sources");
	  DEBUG_PRINT(sources);
	  sources->draw(fullSourceFilename,false,args->gridSpace > 0,args->gridSpace);
	}
#endif
      DEBUG_PRINT("DELETING");
      delete lens;
      DEBUG_PRINT("deleted lens");
      delete sources;
      DEBUG_PRINT("delted sources");
      lens = new Plane<Double>(width,height, bgColor);
      sources = new Plane<Double>(width,height, bgColor);
    }


  //making a combined image if multiple sources are used.
  if(numberOfSources > 1)
    {
      verbosePrint(args,"Combined Sim");
      args->fileNamePrefix = args->mainPrefix + Double(numberOfSources).str();
      simulation(args,lens,sources, massDensity, numberOfSources,-1,0);
#ifndef __USE_BOINC__
      lens->draw(args->fileNamePrefix+"lensedimage.bmp",false,args->gridSpace < 0,args->gridSpace);
      sources->draw(args->fileNamePrefix+"sources.bmp",false,args->gridSpace < 0,args->gridSpace);
#endif
    }

  delete randy;
  delete lens;
  delete sources;
  if(massDensity != 0)
    {
      massDensity->clearAllFields();
      delete massDensity;
    }
	
  massDensity = 0;
  randy = 0;
  lens = sources = lensMassPlane = 0;
	
  return 0;
}



int simulation(struct ray_trace_arguments *args,Plane<Double> * lens, Plane<Double> * sources, DensityProfile * massDensity, int numberOfSources, int specificSource, int specificLens) throw (DavidException)
{
	
	
  using namespace std;
  using utils::DRandom;
	
  int N = lens->numberOfRows();//Plane is always square
	
  if(args == NULL)
    throw DavidException("simulation: Need arguments. Got null pointer.",DavidException::NULL_POINTER);
  args->includeCricalCurveAndCaustic = false;
#if 0// deprecated??
  bool useRandom = false;//source placement
  if(useRandom)
    randy = new DRandom();
#endif

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

      writeSourceInfo(args,sourceParams,specificSource, numberOfSources, randy, args->useRandom, N);
    }
  else
    {
      verbosePrint(args,"Loading parameters");
      loadParameters(params,lensParams,sourceParams,numberOfObsParameters,10,numberOfSources);

    }
  delete massDensity;//No memory leaks here!

  verbosePrint(args,"Creating Mass density");
  massDensity = new DensityProfile(args->lensMassDensity, lensParams, params);
  GLAlgorithm * gls = new GLAlgorithm[numberOfSources];

  for(int i = 0;i<numberOfSources;i++)
    {
      if(args->runExistingDeflection)
	{
	  // USE EXISTING DEFLECTION PLANE
	  const std::string& deflectionFilename = args->lensMassDeflectionPlane;
	  
	  verbosePrint(args,std::string("Parsing Deflection Plane: ")+deflectionFilename);
	  deflectionPlane = Plane<math::Complex>::readPlane(deflectionFilename);
	  gls[i] = GLAlgorithm(deflectionPlane,params,sourceParams[0],sources);
	}
      if(args->createDeflection)
	{

	  // CREATE NEW DEFLECTION PLANE 
	  verbosePrint(args,std::string("Creating Deflection Plane: ") + args->lensMassDeflectionPlane);
	  gls[i] = GLAlgorithm(massDensity,sourceParams[0],sources,glellipseBounds,args->offset);

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
	    	  gls[i].getDeflectionPlane()->savePlane(args->lensMassDeflectionPlane);

#else
	      gls[i].getDeflectionPlane()->savePlane(args->lensMassDeflectionPlane);
#endif
	    }
#ifdef USE_MPI
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
      if(!args->createDeflection && !args->runExistingDeflection)
	{
	  verbosePrint(args,"Running based on parameters");
	  gls[i] = GLAlgorithm(deflectionPlane,params,sourceParams[0],sources);
	}
		  
    }

#ifndef __USE_BOINC__
  verbosePrint(args,"Running Sim");
  if(args->runSim)//Run ray trace, else no sim
    {
      size_t X,Y,X0,Y0;
      Plane<math::Complex> * sourceLocations = 0;
      if(args->savesourcelocations.size() != 0)
	sourceLocations = new Plane<math::Complex>(N,N,math::Complex(0.0,0.0));

      X = lens->numberOfRows();
      Y = lens->numberOfColumns();
      X0 = 0;
      Y0 = 0;
      if(glellipseBounds[0] > 0)
	X0 = glellipseBounds[0];
      if(glellipseBounds[1] < N)
	X = glellipseBounds[1];
      if(glellipseBounds[2] > 0)
	Y0 = glellipseBounds[2];
      if(glellipseBounds[3] > 0)
	Y = glellipseBounds[3];
      
#ifdef USE_MPI
  DEBUG_PRINT("Process " << mpi_data.rank << ": Making plane dim: " << glellipseBounds[0] << ", " << glellipseBounds[1] << ", " << glellipseBounds[2] << ", " << glellipseBounds[3]);

#endif

      int percentFinished = 0;
      Double sourceValue;
      for(size_t i = X0;i<X;i++)
	{
	  for(size_t j = Y0;j<Y;j++)
	    {
	      for(size_t k = 0;k<numberOfSources;k++)
		{
		  try
		    {
		      sourceValue = lens->getValue(i,j) + gls[k].mapsToSource(i,j);
		      if(sourceLocations != 0 && args->savesourcelocations.size() != 0)
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
	      verbosePrint(args,"Percent finished: ");
	      verbosePrint(args,percentFinished);
	    }
	}

      verbosePrint(args,"Simulation Complete");

      if(args->savesourcelocations.size() && sourceLocations != 0)
	{
	  verbosePrint(args,std::string("Writing source plane to ") + args->savesourcelocations);
#ifdef USE_MPI
	  utils::mpi_recombine(sourceLocations,MPI_COMM_WORLD);
	  if(mpi_data.rank == utils::MASTER_RANK)
	    sourceLocations->savePlane(args->savesourcelocations);
#else
	  sourceLocations->savePlane(args->savesourcelocations);
#endif
	}
    }//end if(runsim)
  else
    {
      verbosePrint(args,"No sim");
    }
#endif

	
#ifndef __USE_BOINC__
  writeParameters(args,params,lensParams,sourceParams,numberOfSources);
#endif

  delete [] gls;
  for(int i = 0; i< numberOfSources;i++)
    delete [] sourceParams[i];
  delete [] sourceParams;

  gls = 0;
  params = 0;
  lensParams = 0;//don't delete here. Contained in massDensity which is kept as a pointer.
  sourceParams = 0;

  return 0;
}


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
bool writeSourceInfo(struct ray_trace_arguments *args,
		     double ** sourceParams, 
		     int specificSource, int numberOfSources, 
		     utils::DRandom * randy, bool useRandom, 
		     int numberOfPixels)
{
  using namespace std;

  string fullFileName;

  if(args == NULL)
    throw DavidException("writeSourceInfo: Need arguments. Got null pointer",DavidException::NULL_POINTER);

  fullFileName = args->fileNamePrefix + "sourceinfo.txt";

  std::ofstream sourceLOC(fullFileName.c_str());

  int N = numberOfPixels;

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
				
		
      if(args->useRandom)
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
bool writeParameters(struct ray_trace_arguments *args,
		     double * params,double * lensParams,double ** sourceParams, 
		     int numberOfSources)
{
  if(args == NULL)
    throw DavidException("writeParameters: Need arguments. Got null pointer",DavidException::NULL_POINTER);

  std::string fullFileName = args->fileNamePrefix+"human-parameters.txt";
  std::string fullFileName2 = args->fileNamePrefix+"parameters.txt";
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
      if(args->includeCricalCurveAndCaustic)
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


void default_arguments(struct ray_trace_arguments *args)
{
  if(args == NULL)
    return;

  args->runExistingDeflection = true;
  args->createDeflection = false;
  args->runSim = true;
  args->offset = NULL;
  args->lensMassDensity = NULL;
  args->gridSpace = -1000;
  args->mainPrefix = "run-";
  args->sourceBMPFilename = "source.bmp";
  args->lensMassDeflectionPlane = "lens";
  args->drawRemovedArea = false;
  args->useTimeStamp = false;
  args->runAsDaemon = false;
  args->verbose = false;
  args->useRandom = false;
}

int planes_op(char *csv_names,const char op)
{
  char *token, *saveptr;
  std::string lhs_name, rhs_name, result_name;
  token = strtok_r(csv_names,",",&saveptr);
  if(token == NULL)
    throw DavidException("planes_op: Usage <left plane>,<right plane>,<result name>",DavidException::DATA_FORMAT_ERROR);
  lhs_name = token;

  token = strtok_r(NULL,",",&saveptr);
  if(token == NULL)
    throw DavidException("planes_op: Usage <left plane>,<right plane>,<result name>",DavidException::DATA_FORMAT_ERROR);
  rhs_name = token;
  
  token = strtok_r(NULL,",",&saveptr);
  if(token == NULL)
    throw DavidException("planes_op: Usage <left plane>,<right plane>,<result name>",DavidException::DATA_FORMAT_ERROR);
  result_name = token;
  
  if(!access(lhs_name.c_str(),R_OK))
    throw DavidException("planes_op: Could not open left plane.",DavidException::IO_ERROR_CODE);
  if(!access(rhs_name.c_str(),R_OK))
    throw DavidException("planes_op: Could not open right plane.",DavidException::IO_ERROR_CODE);
  if(!access(result_name.c_str(),W_OK))
    throw DavidException("planes_op: Cannot write result file.",DavidException::IO_ERROR_CODE);
  
  Plane<Double> *lhs = NULL, *rhs = NULL, *result = NULL;
  try{
    lhs = Plane<Double>::readPlane(lhs_name);
    rhs = Plane<Double>::readPlane(rhs_name);
    switch(op)
      {
      case '+':
	result = Plane<Double>::addPlanes(lhs,rhs);
	break;
      case '-':
	result = Plane<Double>::subtractPlanes(lhs,rhs);
	break;
      default:
	{
	  std::string error = "planes_op: Unknown operator: " + op;
	  throw DavidException(error,DavidException::FORMAT_ERROR_CODE);
	}
      }
    if(result != NULL)
      result->savePlane(result_name);
    else
      throw DavidException("planes_op: Result was NULL.",DavidException::DATA_FORMAT_ERROR);
  }
  catch(DavidException de)
    {
      if(lhs != NULL)
	delete lhs;
      if(rhs != NULL)
	delete rhs;
      if(lhs != NULL)
	delete result;
      throw de;
    }

  delete lhs;
  delete rhs;
  delete result;
  
  return 0;
}

void parse_bounds(char *csv_bounds)
{
  char *token, *saveptr;

  if(glellipseBounds == NULL)
    glellipseBounds = new int[4];

    
  token = strtok_r(csv_bounds,",",&saveptr);
  for(size_t i = 0;i<4;i++)
    {
      if(token == NULL)
	throw DavidException("parse_bounds: Usage <left plane>,<right plane>,<result name>",DavidException::DATA_FORMAT_ERROR);
      if(sscanf(token,"%i",&glellipseBounds[i]) != 1)
	throw DavidException("parse_bounds: Could not read the bounds.",DavidException::FORMAT_ERROR_CODE);
      token = strtok_r(NULL,",",&saveptr);
    }
}

int parseArgs(int argc, char** argv, struct ray_trace_arguments *args)
{
  static const char short_opts[] = "dg:hl:m:n:p:tv";
  int option_index = 0, c;
  void print_help();

  while(true)
    {
      c = getopt_long(argc,argv,short_opts,ray_trace_options,&option_index);
      if(c == -1)
	break;

      switch(c)
	{
	case 'n':
	  {
	    args->lensMassDeflectionPlane = optarg;
	    args->runExistingDeflection = false;
	    args->createDeflection = true;
	    args->runSim = false;
	    break;
	  }
	case SOURCE_BMP:
	  args->sourceBMPFilename = optarg;
	  break;
	case RBG_COLOR:
	  std::cout << "--bgcolor has been disabled" << std::endl;
	  break;
	case XOFFSET:
	  {
	    if(args->offset == NULL)
	      args->offset = new double[2];
	    sscanf(optarg,"%f",&args->offset[0]);
	    break;
	  }
	case YOFFSET:
	  {
	    if(args->offset == NULL)
	      args->offset = new double[2];
	    sscanf(optarg,"%f",&args->offset[1]);
	    break;
	  }
	case FILE_PREFIX:
	  args->mainPrefix = optarg;
	  break;
	case 'l':case DEFLECTION_MAP:
	  args->lensMassDeflectionPlane = optarg;
	  args->runExistingDeflection = true;
	  args->createDeflection = false;
	  args->runSim = true;
	  break;
	case FORCE_RUN:
	  args->runSim = false;
	  break;
	case STOP_RUN:
	  args->runSim = true;
	  break;
	case USE_GRID:
	  args->gridSpace = (int) Double(optarg).doubleValue();
	  break;
	case 'p':
	  parameterArray = paramParser(optarg);
	  break;
	case CREATE_MASS:
	  args->makeMassDensity = optarg;
	  break;
	case 'm':
	  args->lensMassDensity = Plane<Double>::readPlane(optarg);
	  break;
	case 'd':
	  args->runAsDaemon = true;
	  break;
	case 't':
	  args->useTimeStamp = true;
	  break;
	case ADD_PLANES:
	  return planes_op(optarg,'+');
	case SUBTRACT_PLANES:
	  return planes_op(optarg,'-');
	case 'g':
	  parse_bounds(optarg);
	  break;
	case SAVESOURCE_LOC:
	  args->savesourcelocations = optarg;
	  break;
	case DRAW_REMOVED_AREA:
	  args->drawRemovedArea = true;
	  break;
	case SQUARE_PLANE:
	  squarePlane(optarg);
	  return 0;
	case 'v':
	  args->verbose = true;
	case 'h':
	default:
	  print_help();
	  return 0;
	}//end of switch(c)
    }

  return -42;
}


void print_help()
{
  struct ray_trace_arguments blah;
  default_arguments(&blah);
  blah.verbose = true;
  verbosePrint(&blah, "Options are:");
  verbosePrint(&blah, "--newlens <lens>");
  verbosePrint(&blah, "Creates a new lens saved as <lens>. Causes simulation to no happen. Override this with --forcerun");
  verbosePrint(&blah, "--createsurfacemassdensity <parameterfile> <outfile>");
  verbosePrint(&blah, "Creates 2-D massdensity (double values) based on the given parameters.");
  verbosePrint(&blah, "--usesurfacemassdensity <infile>");
  verbosePrint(&blah, "Uses an already created 2-D massdensity (double values).");

  verbosePrint(&blah, "--sourceBMP <source.bmp>");
  verbosePrint(&blah, "Uses <source.bmp> as the source plane");
  verbosePrint(&blah, "--deflectionMap <deflectionMap.txt>");
  verbosePrint(&blah, "Uses <deflectionMap.txt> as the deflection plane");
  verbosePrint(&blah, "--lens <deflectionMap.txt>");
  verbosePrint(&blah, "Uses <deflectionMap.txt> as the deflection plane. Same as --deflectionMap  (typing lens was a habit)");
  verbosePrint(&blah, "--forceRun");
  verbosePrint(&blah, "Forces the Simulation to run, whether or not a deflection map was created or read");
  verbosePrint(&blah, "--stopRun");
  verbosePrint(&blah, "Forces the Simulation NOT to run, whether or not a deflection map was created or read");
  verbosePrint(&blah, "--prefix");
  verbosePrint(&blah, "Sets the output file(s) prefix ");
  verbosePrint(&blah, "--useGrid (size)");
  verbosePrint(&blah, "Draws a grid on output images. Size is the pixel spacing of the grid.");
  verbosePrint(&blah, "--bgcolor R G B");
  verbosePrint(&blah, "Sets the background color of the images to the RGB value provided");
  verbosePrint(&blah, "--offset x y");
  verbosePrint(&blah, "Set a pixel offset for the calculation.");
  verbosePrint(&blah, "--parameters <file name>");
  verbosePrint(&blah, "Parse file provided file for parameters.");
  verbosePrint(&blah, "--drawsource a,e,width,height,x,y,R,G,B,filename");
  verbosePrint(&blah, "Draws an RGB source of eccentricity e, semimajor axis length a, at (x,y), saved as width by height sized a bitmap called filename. Program exits when completed.");
  verbosePrint(&blah, "--glellipsebounds left,right,upper,lower");
  verbosePrint(&blah, "Sets bounds on calculation of deflection grid. This allows embarrisingly parallel processing.");
  verbosePrint(&blah, "--glellipse left,right,upper,lower");
  verbosePrint(&blah, "Same as glellipsebounds. Added to keep typos from killing processing time.");
  verbosePrint(&blah, "--timestamp");
  verbosePrint(&blah, "Places a timestamp at the beginning of every STDOUT message.");
  verbosePrint(&blah, "--addplanes leftplane rightplane outplane");
  verbosePrint(&blah, "Adds leftplane and rightplane and saves the sum as outplane");
  verbosePrint(&blah, "--subtractplanes leftplane rightplane outplane");
  verbosePrint(&blah, "Subtracts leftplane and rightplane and saves the difference as outplane");
  verbosePrint(&blah, "--squareplane filename");
  verbosePrint(&blah, "Takes the plane, filename, and makes it a square by buffing the shorter dimension with zeros");
  verbosePrint(&blah, "--savesourcelocations filename");
  verbosePrint(&blah, "Runs raytracing and saves the location of source points in the lens plane.");
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

  if(fileName == NULL)
    return NULL;

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
	  sourceParams[sourceNumber][i] = Double(parameterArray[vectorCount++]).doubleValue();
    }
  
  //Sets the glellipseBounds to the whole array if the bounds are not yet set.
		
		
  if(glellipseBounds[1] < 0)
    glellipseBounds[1] = (int)floor(params[0]);
  if(glellipseBounds[3] < 0)
    glellipseBounds[3] = (int)floor(params[0]);
  if(glellipseBounds[0] < 0)
    glellipseBounds[0] = 0;
  if(glellipseBounds[2] < 0)
    glellipseBounds[2] = 0;
  
}



void drawEllipse(const struct ray_trace_arguments *args, const char * parameters)
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
	
  verbosePrint(args,std::string("Drawing ")+filename);
  planizzle->draw(filename);
	
  delete newPlane;
  delete planizzle;

}

void createSurfaceMassDensity(const std::string& fileName)
{
  if(parameterArray==NULL)
    throw DavidException("createSurfaceMassDensity: Cannot create a mass density map without paramters.",DavidException::NULL_POINTER);

  int numberOfSources = 1;
  double * params = new double[6];
  double * lensParams = new double[10];
  double ** sourceParams = new double*[numberOfSources];
  loadParameters(params,lensParams,sourceParams,6,10,numberOfSources);
  
  Plane<Double> * dummy = new Plane<Double>(params[0],params[0],42.0);
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


void verbosePrint(const struct ray_trace_arguments *args, const char * string)
{
  if(args->verbose == 0)
    return;

  if(args->runAsDaemon || args->useTimeStamp)
    std::cout << getTime() << ": ";
  std::cout << string << std::endl;
}

void squarePlane(const char * fileName)
{
  
  using math::Complex;
  
  static const Complex cZero(0.0,0.0);

  Plane<Complex> * oldPlane = new Plane<Complex>(0,0,cZero);
  oldPlane = Plane<Complex>::readPlane(fileName);

  int width,height;

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
