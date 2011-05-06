#include "flatten.h"

int wrappedmain(int argc, char** argv)
{
  bool useTimeStamp = false;
  bool scaleImage = false;
  bool drawImage = false;
  int lineOfSight = 2;
  int timeInSeconds = 0;
  
  Flatten f(useTimeStamp, scaleImage,drawImage,lineOfSight,timeInSeconds);

  int argReturned = f.parseArgs(argc,argv);
  if(argReturned != -42)
    return argReturned;

  printf("running as daemon\n");
  return f.createSurfaceMassDensity(f.clusterFilename.c_str(), f.clusterResolution,f.clusterOutLENS.c_str());
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
