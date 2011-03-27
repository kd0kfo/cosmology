#include "flatten.h"

int wrappedmain(int argc, char** argv)
{

  bool runAsDaemon = false;
  bool useTimeStamp = false;
  bool scaleImage = false;
  bool drawImage = false;
  int lineOfSight = 2;
  int timeInSeconds = 0;
  
  Flatten f(runAsDaemon, useTimeStamp, scaleImage,drawImage,lineOfSight,timeInSeconds);

  int returnMe = 0;
    /**/
  

  int argReturned = f.parseArgs(argc,argv);
  if(argReturned != -42)
    return argReturned;


  if(!f.runAsDaemon)
    {
      DEBUG_PRINT("runasdaemon is false");
      returnMe = f.createSurfaceMassDensity(f.clusterFilename.c_str(), f.clusterResolution,f.clusterOutLENS.c_str());
    }
  else
    {
      DEBUG_PRINT("running as daemon");
#ifndef __USE_BOINC__	
	//Running as daemon if you get this far
	// Initialize the logging interface 
	openlog( DAEMON_NAME, LOG_PID, LOG_LOCAL5 );
	syslog( LOG_INFO, "starting" );
	
	// One may wish to process command line arguments here 
	// Daemonize 
	daemonize( "/var/lock/subsys/" DAEMON_NAME );
	
	// Now we are a daemon -- do the work for which we were paid 
	returnMe = f.createSurfaceMassDensity(f.clusterFilename.c_str(), f.clusterResolution,f.clusterOutLENS.c_str());
	
	// Finish up 
	syslog( LOG_NOTICE, "terminated" );
	closelog();
#else
	DEBUG_PRINT("not running daemon");
	returnMe = f.createSurfaceMassDensity(f.clusterFilename.c_str(), f.clusterResolution,f.clusterOutLENS.c_str());
#endif
    }

  //delete [] xBounds;
  //delete [] yBounds;
  //xBounds = yBounds = 0;



  return returnMe;

}

int main(int argc, char** argv)
{
  
  try{
    return wrappedmain(argc, argv);
      }
  catch(DavidException de)
    {
#ifdef __DEBUG__
      DEBUG_PRINT("Type: " << de.getType());
      DEBUG_PRINT("Code: " << de.getCode());
      DEBUG_PRINT("Message: ");
#endif

      de.stdOut();
      
      return de.getCode();
     }
}
