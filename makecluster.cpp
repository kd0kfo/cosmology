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

  DString errorString("Usage: makecluster filename [total mass/convergence] [number of Particles] [redshift] [cluster type] [parameters] \nCluster types and parameters:\n\n nfw:[semimajor axis length] [ratio to next biggest] [ratio to smallest]\n nsis:[maximum distance] [core radius]\n ellipse: [axial ratio] [semimajor axis length]");

  if(argc < 6)
    {
      std::cout << errorString << std::endl;
      return 1;
    }
  
  int counter = 0;
  


  DString fileName(argv[1]);
  double totalMass = Double(argv[2]).doubleValue();
  double numberOfPoints = Double(argv[3]).doubleValue();
  double redshift = Double(argv[4]).doubleValue();
  DString clusterType(argv[5]);int typeNumber = -1;
  double * parameters;

  double arcToKPC = 3.24077649e-22*Cosmology::arcsecondsToCentimeters(1,redshift);
  clusterType.toLowerCase();
  if(clusterType.equals("nfw"))
    {
      if(argc < 9)
	throw DavidException(errorString);

      parameters = new double[3];
      parameters[2] = Double(argv[6]).doubleValue();
      parameters[0] = Double(argv[7]).doubleValue();
      parameters[1] = Double(argv[8]).doubleValue();
      typeNumber = Create_Cluster::NFW;
    }
  else if(clusterType.equals("nsis"))
    {
      if(argc < 8)
	throw DavidException(errorString);

      parameters = new double[2];
      parameters[1] = Double(argv[6]).doubleValue();
      parameters[0] = Double(argv[7]).doubleValue();
      typeNumber = Create_Cluster::NSIS;
    }
  else if(clusterType.equals("ellipse"))
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
  
  if(fileName.charAt(0) == '-')
    throw DavidException(errorString);
  
  std::ofstream myfile(fileName.toCharArray(), std::ios::out);

  if(!myfile.is_open())
    throw DavidException(DString("Could not open: ")+fileName,DavidException::IO_ERROR_CODE);

  try{
    
    Create_Cluster cc;
    
    cc.resetRandom(time(NULL));
    
    int percentFinished = 0;

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
	    percentFinished = (int) i*100/numberOfPoints;
	    VERBOSE_PRINT("Percent finished: ");
	    VERBOSE_PRINT(percentFinished);
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
