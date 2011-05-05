/**
 * Created By David Coss, 2007
 */
#include "ray_trace_ellipse.h"

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


void DefaultParameters(struct general_parameters *params)
{
  if(params == NULL)
    return;

  params->ndim = 480; //N number of pixels, area = 
  params->arcsec_per_pixel = .38; //arc length of each pixel in arcseconds
  params->G = 6.67*pow((double) 10,(double) -8);//G, newtonian grav. constant
  params->solar_mass = 2*pow((double) 10,(double) 33);//Solar mass, grams
  params->distance_scale = 3.09*pow((double) 10, (double) 24);//distance scale Mpc
  params->c = 3*pow(10.0,10.0);//speed of light in cm/s
}

int main(int argc, char** argv)
{
  int retval = 0;
  time_t start_time = time(NULL),end_time;
  
  printf("ray_trace_ellipse started on %s\n",ctime(&start_time));
 
  glellipseBounds = NULL;
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
      std::cerr << "Process of Rank " << mpi_data.rank  << " on host " << mpi_data.hostname << std::endl;
#endif
      std::cerr << "*************" << std::endl << "Program Ended due to " << de.getCause() << std::endl;
      std::cerr << "The following message was provided:" << std::endl;
      de.stdErr();

      retval = de.getCode();
    }

#ifdef USE_MPI
  MPI_Finalize();
#endif

  end_time = time(NULL);
  printf("Calculation duration: %f second\n",difftime(end_time,start_time));
  printf("ray_trace_ellispe ended on %s\n",ctime(&end_time));
  return retval;
}



int super_main(int argc, char** argv)
{
  struct ray_trace_arguments args;

  default_arguments(&args);


  glellipseBounds = new int[4];
  for(int i = 0;i<4;i++)
    glellipseBounds[i] = -1;

  int retval = parseArgs(argc,argv,&args);
  if(args.runAsDaemon)
    printf("DAEMON PROCESS\n");

  if(args.makeMassDensity.size() != 0)
      createSurfaceMassDensity(args.makeMassDensity,args);

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
    returnMe = run_simulation(args);
  }
  catch(DavidException de)
    {
      std::cerr << "Exception Thrown!" << std::endl;
      std::cerr << "Type: " << de.getType() << std::endl;
      std::cerr << "Code: " << de.getCode() << std::endl;
      std::cerr << "Message: " << de.getMessage() << std::endl;
      returnMe = de.getCode();
    }
  return returnMe;
}

int run_simulation(struct ray_trace_arguments *args)
{
  using utils::DRandom;
  DRandom * randy = NULL;

  Plane<Double> * lens = NULL;
  Plane<Double> * sources = NULL;
  Plane<Double> * lensMassPlane = NULL;
  Plane<math::Complex> * deflectionPlane = NULL;
#ifndef __USE_BOINC__
  if(access(args->sourceBMPFilename.c_str(),R_OK) == 0)
    sources = Plane<Double>::bmpToPlane(args->sourceBMPFilename);
#endif
  DensityProfile * massDensity = NULL;

  using std::string;
  
  std::string fullLensFilename = args->fileNamePrefix  + "lensedimage.bmp";
  simulation(args,&lens,&sources, &massDensity);


#ifndef __USE_BOINC__
  if(lens != NULL)
    lens->draw(fullLensFilename,false,args->gridSpace > 0,args->gridSpace);
		
  if(sources != NULL)
    {
      std::string fullSourceFilename = args->fileNamePrefix + "sources.bmp";
      verbosePrint(args,"Drawing sources");
      sources->draw(fullSourceFilename,false,args->gridSpace > 0,args->gridSpace);
    }
#endif

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



int simulation(struct ray_trace_arguments *args,Plane<Double> **lens_addr, Plane<Double> **sources_addr, DensityProfile **massDensity_addr) throw (DavidException)
{
	
	
  using namespace std;
  using utils::DRandom;
  Plane<Double> *lens, *sources;
  DensityProfile *massDensity;
  int N = 0;
	
  if(args == NULL || lens_addr == NULL || sources_addr == NULL || massDensity_addr == NULL)
    throw DavidException("simulation: Need arguments. Got null pointer.",DavidException::NULL_POINTER);
  
  massDensity = *massDensity_addr;
  args->includeCricalCurveAndCaustic = false;

  struct lens_parameters lensParams;
  struct general_parameters params;
  struct source_parameters sourceParams;

  if(args == NULL || args->parameter_name.size() == 0)
    {
      createLensParams(&lensParams,0);
      DefaultParameters(&params);
    }
  else
    {
      verbosePrint(args,"Loading parameters");
      loadParameters(args,&params,&lensParams,&sourceParams);

    }

  N = params.ndim;
  if(*lens_addr == NULL)
    *lens_addr = new Plane<Double>(N,N,0.0);
  lens = *lens_addr;

  if(*sources_addr == NULL)
    *sources_addr = new Plane<Double>(N,N,0.0);
  sources = *sources_addr;


  
#ifdef USE_MPI
  utils::mpi_adjust_glellipsebounds(glellipseBounds,N);
  DEBUG_PRINT("Process " << mpi_data.rank << ": Making plane dim: " << glellipseBounds[0] << ", " << glellipseBounds[1] << ", " << glellipseBounds[2] << ", " << glellipseBounds[3]);

#endif

  // setup planes
  delete massDensity;//No memory leaks here!

  verbosePrint(args,"Creating Mass density");
  massDensity = new DensityProfile(&args->lensMassDensity, &lensParams, &params);  *massDensity_addr = massDensity;
  
  GLAlgorithm gl;

  if(args->runExistingDeflection)
    {
      // USE EXISTING DEFLECTION PLANE
      const std::string& deflectionFilename = args->lensMassDeflectionPlane;
	  
      verbosePrint(args,std::string("Parsing Deflection Plane: ")+deflectionFilename);
      deflectionPlane = Plane<math::Complex>::readPlane(deflectionFilename);
      gl = GLAlgorithm(deflectionPlane,&params,&sourceParams,sources);
    }
  if(args->createDeflection)
    {

      // CREATE NEW DEFLECTION PLANE 
      verbosePrint(args,std::string("Creating Deflection Plane: ") + args->lensMassDeflectionPlane);
      gl = GLAlgorithm(massDensity,&sourceParams,sources,glellipseBounds,args->offset);

#ifdef USE_MPI
      MPI_Comm lens_group;
      int has_lens = (gl.getDeflectionPlane() != 0) ? 1 : 0;
      MPI_Comm_split(MPI_COMM_WORLD,has_lens,mpi_data.rank,&lens_group);
#endif

      if(gl.getDeflectionPlane() != 0)
	{
#ifdef USE_MPI
	  utils::mpi_recombine(gl,lens_group);
	  if(mpi_data.rank == utils::MASTER_RANK)
	    gl.getDeflectionPlane()->savePlane(args->lensMassDeflectionPlane);

#else
	  gl.getDeflectionPlane()->savePlane(args->lensMassDeflectionPlane);
#endif
	}
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
  if(!args->createDeflection && !args->runExistingDeflection)
    {
      verbosePrint(args,"Running based on parameters");
      gl = GLAlgorithm(deflectionPlane,&params,&sourceParams,sources);
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
	      try
		{
		  sourceValue = lens->getValue(i,j) + gl.mapsToSource(i,j);
		  if(sourceLocations != 0 && args->savesourcelocations.size() != 0)
		    {
		      double * loc = gl.newLocation(i,j);
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
  writeParameters(args,params,lensParams,sourceParams);
#endif

  *massDensity_addr = massDensity;
  return 0;
}


template <class T> void buildSource(Plane<T> * source, const struct source_parameters& sourceParams, T newValue)
{
  
  double xCenter = sourceParams.xcenter;
  double yCenter = sourceParams.ycenter;
  double a = sourceParams.semimajor_length;
  double e = sourceParams.eccentricity;
	
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



/**
 * writeParameters writes the parameters to 2 files: "human readable" file and a file which can be parsed by ray_trace
 */
bool writeParameters(const struct ray_trace_arguments *args,
		     const struct general_parameters& params,const struct lens_parameters& lensParams,const struct source_parameters& sourceParams)
{
  using namespace std;
  using namespace utils;
  if(args == NULL)
    throw DavidException("writeParameters: Need arguments. Got null pointer",DavidException::NULL_POINTER);

  static const string fullFileName = args->fileNamePrefix+"human-parameters.txt";
  static const string fullFileName2 = args->fileNamePrefix+"parameters.txt";
  ofstream outfile(fullFileName.c_str());
  XMLParser xmldoc;
  XMLNode *xml = new XMLNode("parameters");
  XMLNode *xmlparams = new XMLNode("general_parameters");
  XMLNode *xmllens = new XMLNode("lens_parameters");
  ostringstream buff;
  
  xml->insertChild(xmlparams);
  xml->insertChild(xmllens);

  if(outfile.is_open())
    {
      
      outfile << "General Parameters:" << endl;
      outfile << "Length and Width in pixels: " << params.ndim << endl;
      buff.clear();buff << params.ndim;
      xmlparams->insertChild("ndim",buff.str());
      outfile << "Arc length per pixel: " << params.arcsec_per_pixel << endl;
      buff.clear();buff << params.arcsec_per_pixel ;
      xmlparams->insertChild("arcsec_per_pixel",buff.str());
      outfile << "Value of Newton's Constant: " << params.G << endl;
      buff.clear();buff << params.G ;
      xmlparams->insertChild("G",buff.str());
      outfile << "Solar mass: " << params.G << endl;
      buff.clear();buff << params.G ;
      xmlparams->insertChild("G",buff.str());
      outfile << "Distance Scale: " << params.distance_scale << endl;
      buff.clear();buff << params.distance_scale ;
      xmlparams->insertChild("distance_scale",buff.str());
      outfile << "Speed of Light" << params.c << endl;
      buff.clear();buff << params.c ;
      xmlparams->insertChild("c",buff.str());
      outfile << "Critical Curves included: ";
      if(params.showCriticalCurves)
	{
	  outfile << "Yes" << endl;
	  xmlparams->insertChild("c","1");
	}
      else
	{
	  outfile << "No" << endl;
	  xmlparams->insertChild("c","0");
	}
      outfile << endl;
    }
  else
    return false;

  if(outfile.is_open())
    {
      outfile << "Lens Parameters:" << endl;
      outfile << "x coordinate of lens center: " << lensParams.xcenter << endl;
      buff.clear();buff << lensParams.xcenter ;
      xmllens->insertChild("xcenter",buff.str());
      outfile << "y coordinate of lens center: " << lensParams.ycenter << endl;
      buff.clear();buff << lensParams.ycenter ;
      xmllens->insertChild("ycenter",buff.str());
      outfile << "lens redshift or distance: " << lensParams.redshift << endl;
      buff.clear();buff << lensParams.redshift ;
      xmllens->insertChild("redshift",buff.str());
      outfile << "lens velocity of dispersion: " << lensParams.velocity_dispersion << endl;
      buff.clear();buff << lensParams.velocity_dispersion ;
      xmllens->insertChild("velocity_dispersion",buff.str());
      outfile << "mass density (optional): " << lensParams.mass_density << endl;
      buff.clear();buff << lensParams.mass_density;
      xmllens->insertChild("mass_density",buff.str());
      outfile << "Shear: " << lensParams.shear << endl;
      buff.clear();buff << lensParams.shear;
      xmllens->insertChild("shear",buff.str());
      outfile << "Semi-major axis length: " << lensParams.semimajor_length << endl;
      buff.clear();buff << lensParams.semimajor_length;
      xmllens->insertChild("semimajor_length",buff.str());
      outfile << "Eccentricity(e=1-b/a): " << lensParams.eccentricity << endl;
      buff.clear();buff << lensParams.eccentricity;
      xmllens->insertChild("eccentricity",buff.str());
      outfile << endl;
    }
  else
    return false;

  if(outfile.is_open())
    {
      XMLNode *xmlsourceParams = new XMLNode("source_parameters");
      outfile << "Source Parameters:" << endl;
      xml->insertChild(xmlsourceParams);
      outfile << "x coordinate of source center: " << sourceParams.xcenter << endl;
      buff.clear();buff << sourceParams.xcenter;
      xmlsourceParams->insertChild("xcenter",buff.str());
      outfile << "y coordinate of source center: " << sourceParams.ycenter << endl;
      buff.clear();buff << sourceParams.ycenter;
      xmlsourceParams->insertChild("ycenter",buff.str());
      outfile << "source redshift or distance: " << sourceParams.redshift << endl;
      buff.clear();buff << sourceParams.redshift;
      xmlsourceParams->insertChild("redshift",buff.str());
      outfile << "semimajor axis length in arcseconds: " << sourceParams.semimajor_length << endl;
      buff.clear();buff << sourceParams.semimajor_length;
      xmlsourceParams->insertChild("semimajor_length",buff.str());
      outfile << "source eccentricity: " << sourceParams.eccentricity << endl;
      buff.clear();buff << sourceParams.eccentricity;
      xmlsourceParams->insertChild("eccentricity",buff.str());
      outfile << "orientation of source semimajor axis in radians: " << sourceParams.orientation << endl;
      buff.clear();buff << sourceParams.orientation;
      xmlsourceParams->insertChild("orientation",buff.str());
      outfile << "Source flux: " << sourceParams.flux << endl;
      buff.clear();buff << sourceParams.flux;
      xmlsourceParams->insertChild("flux",buff.str());
      outfile << endl;
    }
  else
    return false;

  outfile.close();

  xmldoc.setHead(xml);
  try{
    xmldoc.write(fullFileName2);
  }
  catch(DavidException de)
    {
      de.stdErr();
      return false;
    }

  //xmlparser will delete the nodes...
  return true;

}


bool createLensParams(struct lens_parameters *lensParams, int specificLens)
{
  // currently, there is only 1 lens type. specificLens is ignored
  lensParams->xcenter = 0;//x coordinate of lens center
  lensParams->ycenter = 0;//y coordinate of lens center
  lensParams->redshift = 4.2;//lens distance in Mpc (typically 1000)
  lensParams->velocity_dispersion= 600;//lens velocity of dispersion
  lensParams->mass_density = 1e16; //in Solar Masses, density if applicable, originally 10^12
  lensParams->shear = 0;//Shear in glSISwithShear
  lensParams->semimajor_length = .38;//a, semi-major axis length of lens, in arcseconds
  lensParams->eccentricity = 0;//e, ellipticity e= 1-b/a, b is the semi-minor axis length

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
  args->fileNamePrefix = "run-";
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
  static const char short_opts[] = "c:dg:hl:m:n:p:tv";
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
	  args->fileNamePrefix = optarg;
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
	  args->parameter_name = optarg;
	  break;
	case 'c':
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
	  break;
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

void loadGeneralParameters(struct general_parameters * params, const utils::XMLNode *node)
{
  utils::XMLNode *curr;

  if(params == NULL || node == NULL)
    return;

  std::istringstream buff;
  double val;
  for(curr = node->children;curr != NULL;curr = curr->siblings)
    {
      const std::string& name = curr->getName();
      buff.clear();buff.str(curr->getText());
      buff >> val;
      if(buff.fail())
	{
	  std::ostringstream error;
	  error << "loadGeneralParameters: Invalid value of " << curr->getText() << " for " << name;
	  throw DavidException(error,DavidException::FORMAT_ERROR_CODE);
	}
	 
      if(name == "ndim")
	params->ndim = (int)ceil(val);
      else if(name == "arcsec_per_pixel")
	params->arcsec_per_pixel = val;
      else if(name == "G")
	params->G = val;
      else if(name == "solar_mass")
	params->solar_mass = val;
      else if(name == "distance_scale")
	params->distance_scale = val;
      else if(name == "showcriticalcurves")
	params->showCriticalCurves;
      else if(name == "c")
	params->c = val;
      else
	std::cerr << "WARNING - Unknown General Parameter: " << name << std::endl;
    }
}

void loadLensParameters(struct lens_parameters * params, const utils::XMLNode *node)
{
  utils::XMLNode *curr;

  if(params == NULL || node == NULL)
    return;

  std::istringstream buff;
  double val;
  for(curr = node->children;curr != NULL;curr = curr->siblings)
    {
      const std::string& name = curr->getName();
      buff.clear();buff.str(curr->getText());
      buff >> val;
      if(buff.fail())
	{
	  std::ostringstream error;
	  error << "loadLensParameters: Invalid value of " << curr->getText() << " for " << name;
	  throw DavidException(error,DavidException::FORMAT_ERROR_CODE);
	}
      
      if(name == "xcenter")
	params->xcenter = (int)ceil(val);
      else if(name == "ycenter")
	params->ycenter = (int)ceil(val);
      else if(name == "redshift")
	params->redshift = val;
      else if (name == "velocity_dispersion")
	params->velocity_dispersion = val;
      else if (name == "mass_density")
	params->mass_density = val;
      else if (name == "shear")
	params->shear = val;
      else if (name == "semimajor_length")
	params->semimajor_length = val;
      else if (name == "eccentricity")
	params->eccentricity = val;
      else
	std::cerr << "WARNING - Unknown Lens Parameter: " << name << std::endl;
    }
}

struct source_parameters readSourceParameters(const utils::XMLNode *node)
{
  utils::XMLNode *curr;

  if(node == NULL)
    throw DavidException("readGeneralParameters: NULL POINTER provided for the XML data. Could not load paramters.",DavidException::NULL_POINTER);

  struct source_parameters params;
  std::istringstream buff;
  double val;
  for(curr = node->children;curr != NULL;curr = curr->siblings)
    {
      const std::string& name = curr->getName();
      buff.clear();buff.str(curr->getText());
      buff >> val;
      if(buff.fail())
	{
	  std::ostringstream error;
	  error << "loadSourceParameters: Invalid value of " << curr->getText() << " for " << name;
	  throw DavidException(error,DavidException::FORMAT_ERROR_CODE);
	}

      if(name == "xcenter")
	params.xcenter = (int)ceil(val);
      else if(name == "ycenter")
	params.ycenter = (int)ceil(val);
      else if(name == "redshift")
	params.redshift = val;
      else if(name == "semimajor_length")
	params.semimajor_length = val;
      else if(name == "eccentricity")
	params.eccentricity = val;
      else if(name == "orientation")
	params.orientation = val;
      else if(name == "flux")
	params.flux = val;
      else
	std::cerr << "WARNING - Unknown Source Parameter: " << name << std::endl;	
    }      
  return params;
}


void loadParameters(struct ray_trace_arguments *args,struct general_parameters * params,struct lens_parameters * lensParams,struct source_parameters *sourceParams)
{
  using namespace utils;

  if(args == NULL || params == NULL || lensParams == NULL || sourceParams == NULL)
    {
      std::cerr << "WARNING - Null pointer was provided to loadParameters. Therefore, no parameters could be loaded." << std::endl;
      return;
    }
    XMLParser parser;
  const XMLNode *curr;
  parser.parse(args->parameter_name);
  for(curr = parser.getHead()->children;curr != NULL;curr = curr->siblings)
    {
      if(curr->getName().find("general") != std::string::npos)
	loadGeneralParameters(params,curr);
      else if(curr->getName().find("lens") != std::string::npos)
	loadLensParameters(lensParams,curr);
      else if(curr->getName().find("source") != std::string::npos)
	*sourceParams = readSourceParameters(curr);
      else
	throw DavidException("loadParameters: Invalid parameter tag: " + curr->getName(),DavidException::FORMAT_ERROR_CODE);
    }
  

  // Default values of glellipseBounds if not set
  if(glellipseBounds == NULL)
    {
      glellipseBounds = new int[4];
      memset(glellipseBounds,-1,4);
    }
  if(glellipseBounds[1] < 0)
    glellipseBounds[1] = (int)floor(params->ndim);
  if(glellipseBounds[3] < 0)
    glellipseBounds[3] = (int)floor(params->ndim);
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

void createSurfaceMassDensity(const std::string& fileName, struct ray_trace_arguments& args)
{
  int numberOfSources = 1;
  struct general_parameters params;
  struct lens_parameters lensParams;
  struct source_parameters sourceParams;
  loadParameters(&args,&params,&lensParams,&sourceParams);
  
  Plane<Double> * dummy = NULL;
  DensityProfile * profile = new DensityProfile(&dummy,&lensParams,&params);

  profile->drawPlane();
  profile->getPlane()->savePlane(fileName);

  delete profile;
  delete dummy;
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
