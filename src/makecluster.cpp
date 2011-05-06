#include <fstream>

#include "libmygl/Cosmology.h"

#include "libdnstd/utils.h"

#include "create_cluster.h"

int uebermain(int argc, char** argv);

int main(int argc, char** argv)
{

  try
    {
      return uebermain(argc, argv);
    }
  catch(DavidException de)
    {
      de.stdErr();
      return de.getCode();
    }
  return -1;
}

int uebermain(int argc, char** argv)
{

  static const char* errorString = "Usage: makecluster filename [total mass/convergence] [number of Particles] [redshift] [cluster type] [parameters] \nCluster types and parameters:\n\n nfw:[semimajor axis length] [ratio to next biggest] [ratio to smallest]\n nsis:[maximum distance] [core radius]\n ellipse: [axial ratio] [semimajor axis length]";

  if(argc < 6)
    {
      std::cout << errorString << std::endl;
      return 1;
    }
  
  int counter = 0;
  


  std::string fileName(argv[1]);
  double totalMass = Double(argv[2]).doubleValue();
  double numberOfPoints = Double(argv[3]).doubleValue();
  double redshift = Double(argv[4]).doubleValue();
  int typeNumber = -1;
  double * parameters;
  std::string clusterType = utils::lower_case(argv[5]);

  double arcToKPC = 3.24077649e-22*Cosmology::arcsecondsToCentimeters(1,redshift);
  if(clusterType == "nfw")
    {
      if(argc < 9)
	throw DavidException(errorString);

      parameters = new double[3];
      parameters[2] = Double(argv[6]).doubleValue();
      parameters[0] = Double(argv[7]).doubleValue();
      parameters[1] = Double(argv[8]).doubleValue();
      typeNumber = Create_Cluster::NFW;
    }
  else if(clusterType == "nsis")
    {
      if(argc < 8)
	throw DavidException(errorString);

      parameters = new double[2];
      parameters[1] = Double(argv[6]).doubleValue();
      parameters[0] = Double(argv[7]).doubleValue();
      typeNumber = Create_Cluster::NSIS;
    }
  else if(clusterType == "ellipse")
    {
      if(argc < 8)
	throw DavidException(errorString);
      
      parameters = new double[2];
      parameters[0] = Double(argv[6]).doubleValue();
      parameters[1] = Double(argv[7]).doubleValue();
      typeNumber = Create_Cluster::UNIFORM_ELLIPSE;
    }
  else
    {
      parameters = 0;
    }
  
  if(fileName.at(0) == '-')
    throw DavidException(errorString);
  
  std::ofstream myfile(fileName.c_str(), std::ios::out);

  if(!myfile.is_open())
    throw DavidException(std::string("Could not open: ")+fileName,DavidException::IO_ERROR_CODE);

  try{
    
    Create_Cluster cc;
    
    cc.resetRandom(time(NULL));
    
    float percentFinished = 0.;

    myfile << "#" << time(NULL) << " " << totalMass << " " << 0 << " " << 0.0 << " " << 0.0 << " " << 0.0 << " " << numberOfPoints << std::endl;
    
    for(int i = 0; i<numberOfPoints;i++)
      {
	counter++;
	Double newCoords = cc.calculatePosition(typeNumber,parameters);

	double x, y, z;
	x = newCoords.getValue(0) * cos(newCoords.getValue(1))*sin(newCoords.getValue(2));
	y = newCoords.getValue(0) * sin(newCoords.getValue(1))*sin(newCoords.getValue(2));
	z = newCoords.getValue(0) * cos(newCoords.getValue(2));
	
	myfile << x << "    " << y << "    " << z << std::endl;
	   
	if((i*100/numberOfPoints) >= (percentFinished+5))
	  {
	    percentFinished = i*100./numberOfPoints;
	    printf("Percent finished: %02.02f\n",percentFinished);
	  }
	
      }
    
  }
  catch(DavidException de)
    {

      de.stdOut();
      return de.getCode();
    }
  catch(...)
    {
      std::cout << "some other error" << std::endl;
      return 1;
    }

  std::cout << "#: " << counter << std::endl;

  return 0;
  
}
