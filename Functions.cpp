#include "Functions.h"


std::set<std::string> * Functions::getFunctions()
{
  using namespace std;
  FunctionHelp fh;
  const map<string, string>* functions = fh.getHelp();
  set<string> * returnMe = new std::set<string>;
  map<string,string>::const_iterator function = functions->begin();
  for(;function != functions->end();function++)
    returnMe->insert(function->first);
  return returnMe;
}




Functions::Functions(const DString& functionType)
{
  this->functionType = new DString(""); *(this->functionType) = functionType;
}




bool Functions::isFunction(const DString& bean)
{
  
  using namespace std;

  bool returnMe;
  FunctionHelp fh;
  const map<string, string>* functions = fh.getHelp();

  returnMe = (functions->find(bean) != functions->end());

  return returnMe;
}


/**
 * Returns the result of a function
 * Object[] uses the follow structure: os[0] = Object[LTree[]], os[1] = LinAl, os[2] = misc
 *
 * @param bean String name of function
 * @param os Object[] parameters of function(see above)
     * @return LTree result
     * @throws DavidException
     */
DString Functions::doFunction(const DString& bean_, DString * parameters, int numberOfParameters,double * gparameters, int numberOfGparameters, std::map<std::string,plane_t> * storedStuff, DString * currentDirectory)
{
  DString returnMe("Result: ");
  DString bean = bean_;
  bean.toLowerCase();

  if(!isFunction(bean))
    {
      throw new DavidException(DString("Error: ")+bean+" is not a defined function.");
    }
  else if(bean.equals("r2d"))
    {

      if(numberOfParameters < 1)
	throw DavidException("I need to know the redshift.");

      Double z1,z2;
      if(numberOfParameters == 2)
	{
	  z1 = Double(parameters[0]);
	  z2 = Double(parameters[1]);
	}
      else
	{
	  z1 = 0.0;
	  z2 = Double(parameters[0]);
	}

      double result = Cosmology::redshiftToDistance(z1.doubleValue(),z2.doubleValue());
      
      returnMe += Double(result).toDString() + DString(" Mpc, based on WMAP 5-year data");
    }
  else if(bean.equals("shear"))
    {
      if(storedStuff->size() < 2 || numberOfParameters < 3)
	throw DavidException("I need a convergence, a kernel and an output name.");
      DString * currentParameters = new DString[4];
      currentParameters[0] = parameters[0];currentParameters[1] = "1";currentParameters[2] = parameters[0] + "-transform";currentParameters[3] = "fourier";
      returnMe = doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);//transform convergence
      returnMe += "\n";
      
      currentParameters[0] = parameters[1];currentParameters[2] = parameters[1]+"-transform";currentParameters[3] = "fourier";//Transform kernel
      returnMe += doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);
      returnMe += "\n";

      for(int newParam = 0;newParam < 3;newParam++)
	currentParameters[newParam] = parameters[newParam] + "-transform";
      currentParameters[3] = "directproduct";//multiply
      returnMe += doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);
      returnMe += "\n";
      
      currentParameters[0] = parameters[2]+"-transform";currentParameters[1] = "-1";currentParameters[2] = parameters[2];currentParameters[3] = "fourier";//transform back to shear in real space
      returnMe += doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);
      returnMe += "\n";
      
      currentParameters[0] = parameters[2];currentParameters[1] = "";currentParameters[3] = "fixdiag";//properly arrange 2-d FFT array
      returnMe += doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);
      returnMe += "\n";
      
      currentParameters[0] = parameters[2];currentParameters[1] = "3.14159";currentParameters[3] = "normalize";//normalize shear integral (1/pi)
      returnMe += doFunction(currentParameters[3],currentParameters,3,0,0,storedStuff,currentDirectory);
      
      for(int newParam = 0;newParam < 3;newParam++)
	storedStuff->erase(parameters[newParam] + "-transform");//remove transform. Maybe add away to keep them if needed, eg keep transform boolean boolean
      
      //done calculating shear
      delete [] currentParameters;currentParameters = 0;

    }
  else if(bean.equals("fixdiag"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 1)
	throw DavidException("One Plane and two Parameters are needed");

      DString plane = parameters[0];
      DString newPlane = plane;
      if(numberOfParameters >= 2 && parameters[1] != "")
	newPlane = parameters[1];

      if(storedStuff->find(plane) == storedStuff->end())
	throw DavidException(plane + " was not found.", DavidException::INVALID_ARGUMENT_ERROR_CODE);
      
      int oldRows, oldColumns;

      plane_t * oldPlane = new plane_t((*storedStuff)[plane]);
      plane_t * newGuy = new plane_t(oldRows = oldPlane->numberOfRows(),oldColumns = oldPlane->numberOfColumns(),0.0);
      
      for(int i = 0;i<oldRows;i++)
	for(int j = 0;j<oldColumns;j++)
	  {
	    Double curr = oldPlane->getValue(i,j);
	    int x = (i+oldRows/2) % oldRows;
	    int y = (j+oldColumns/2) % oldColumns;
	    newGuy->setValue(x,y,curr);
	  }

      

      (*storedStuff)[newPlane] = *newGuy;
      
      delete oldPlane;
      delete newGuy;
      oldPlane = newGuy = 0;

      returnMe = plane + " was fixed";
      if(plane != newPlane)
	returnMe += DString(" and saved as ")+newPlane;

      returnMe += ".";
      
      
    } 
  else if(bean.equals("crit_density"))
    {
      if(numberOfParameters < 2)
	throw DavidException("I need to know the redshift values.");
      
      Double z1(parameters[0]);
      Double z2(parameters[1]);
      Double sigma_c = Cosmology::criticalDensity(z1.doubleValue(),z2.doubleValue());
      returnMe += sigma_c.toDString() + " is the critical density.";
    }
  else if(bean.equals("new"))
    {
      if(numberOfParameters < 3)
	throw DavidException("I need to know the name of the new plane and its dimensions.");

      Double value(0.0);

      if(numberOfParameters >= 4)
	value = Double(parameters[3]);

      Double rows(parameters[1]);
      Double columns(parameters[2]);

      plane_t newGuy(rows.toInt(),columns.toInt(),value);
      
      (*storedStuff)[parameters[0]] = newGuy;
      
      returnMe = DString("A new plane was saved as ")+parameters[0];
       
    }
  else if(bean.equals("drawellipses"))
    {
      if(numberOfParameters < 2)
	throw DavidException("I need to know the name of the plane and file");
      
      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0]+" was not found.");
      
      plane_t original = (*storedStuff)[parameters[0]];
      
      plane_t drawMe(original.numberOfRows(),original.numberOfColumns(),0.0);
      
      Double color(255);

      int columns = original.numberOfColumns();
      int rows = original.numberOfRows();
      for(int i = 0;i<original.numberOfRows();i++)
	for(int j = 0;j<original.numberOfColumns();j++)
	  {
	    Double curr = original.getValue(i,j);
	    double a = curr.getValue(1);
	    double e = curr.getValue(0);
	    double angle = curr.getValue(2);
	    if(a != 0 || e != 0)
	      {
		Utilities u;
		plane_t * buffer = u.drawEllipse(1,0,a,e,angle,i-rows/2,j-columns/2,color,rows,columns);
		plane_t * otherBuffer = plane_t::addPlanes(&drawMe,buffer);
		drawMe = *otherBuffer;
		
		delete buffer;
		delete otherBuffer;
		otherBuffer = buffer = 0;
	      }
	  }
      
					(*storedStuff)[parameters[1]] = drawMe;

      returnMe = DString("The images of ") + parameters[0] + DString(" was saved as ") + parameters[1];

    }
  else if(bean.equals("flipmatrix"))
    {
      if(numberOfParameters < 2)
	throw DavidException("I need to know the name of the matrix and the operation to perform");
      
      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found");
      
      if(!parameters[1].equals("transpose") && !parameters[1].equals("horizontal") && !parameters[1].equals("vertical"))
	throw DavidException(DString("What does ") + parameters[1] + " mean?");
	 
      DString toDo = parameters[1];

      plane_t original = (*storedStuff)[parameters[0]];
      plane_t * newPlane = 0;
      int rows = original.numberOfRows();
      int columns = original.numberOfColumns();
      Double curr;

      if(toDo.equals("transpose"))
	newPlane = original.transpose();
      else if(toDo.equals("horizontal"))
	{
	  newPlane = new plane_t(rows,columns,0.0);
	  for(int i = 0;i<rows;i++)
	    for(int j = 0;j<columns;j++)
	      {
		curr = original.getValue(i,j);
		newPlane->setValue(rows-1-i,j,curr);
	      }
	}
      else if(toDo.equals("vertical"))
	{
	  newPlane = new plane_t(rows,columns,0.0);
	  for(int i = 0;i<rows;i++)
	    for(int j = 0;j<columns;j++)
	      {
		curr = original.getValue(i,j);
		newPlane->setValue(i,columns-1-j,curr);
	      }
	}
      else 
	newPlane = 0;

      if(newPlane == 0)
	throw DavidException("Some error occured. Sorry, that's all I know.");
      
      DString newPlaneName = parameters[0];
      if(numberOfParameters >= 3)
	newPlaneName = parameters[2];
      
					(*storedStuff)[newPlaneName] = *newPlane;
      
      delete newPlane;
      
      returnMe = DString("The ") + parameters[1] + DString(" flipmatrix result was stored as ") + newPlaneName;
       
    }
  else if(bean.equals("addellipseinfo"))
    {
      if(numberOfParameters < 4)
	throw DavidException("I need to know the name of the plane and ellipse information.");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found.");
      
      Double row(parameters[1]);
      Double column(parameters[2]);
      Double a(parameters[3]);
      Double e(0.0);

      if(numberOfParameters >= 5)
	e = Double(parameters[4]);
      

      plane_t original = (*storedStuff)[parameters[0]];
      
      if(row < 0 || row >= original.numberOfRows())
	throw DavidException("That point is out of bounds.");
      if(column < 0 || column >= original.numberOfColumns())
	throw DavidException("That point is out of bounds.");
      
      Double newValue(e.doubleValue(),a.doubleValue(),0.0);

      original.setValue(row.intValue(),column.intValue(), newValue);
      
      (*storedStuff)[parameters[0]] = original;
      
      returnMe = DString("An ellipse was added to ") + parameters[0] + DString(" at (") + row.toDString() + DString(", ") + column.toDString() + ")";

    }
  else if(bean.equals("fourier"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 3)
	throw DavidException("One Plane and two Parameters are needed");

      DString data = parameters[0];
      Double fourierdirection(parameters[1]);//1 for fourier, -1 for inverse
      DString result = parameters[2];

      plane_t * dPlane = &(*storedStuff)[data];
      
      std::complex<double> * junx = PlaneToComplexArray(dPlane);

      int * nn = new int[2];
      nn[0] = dPlane->numberOfRows();
      nn[1] = dPlane->numberOfColumns();

      if((nn[0] & (nn[0]  - 1)) != 0 || (nn[1] & (nn[1] - 1)) != 0)
	{
	  delete junx;
	  delete [] nn;
	  
	  throw DavidException("The dimensions supplied to the FFT must be a power of 2.");
	}

      unsigned fftDirection = FFTW_FORWARD;
      
      if(fourierdirection.doubleValue() == -1)
	fftDirection = FFTW_BACKWARD;
      
      math::FFT::ndim_discrete(junx,nn,2, fftDirection);

      double normconst = 1;
      if(fourierdirection == -1)
	normconst = nn[0]*nn[1];

      dPlane = ComplexArrayToPlane(junx,nn[0],nn[1], normconst);
      
      (*storedStuff)[result] = *dPlane;

      delete dPlane;
      delete [] junx;
      delete [] nn;
      junx = 0;
      dPlane = 0;
      nn = 0;
      DString resultString = "";
      if(fourierdirection == 1.0)
	resultString = "The fourier transform of ";
      else
	resultString = "The inverse fourier transform of ";
	
      returnMe =  resultString + parameters[0] + DString(" is stored as ") + parameters[2];
    }
  else if(bean.equals("header"))
    {
      if(storedStuff->size() == 0 || numberOfParameters < 1)
	throw DavidException("I need a plane");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(DString("I cannot find ") + parameters[0]);

      plane_t& editedPlane = (*storedStuff)[parameters[0]];
      
      if(numberOfParameters > 1)
	{
	  editedPlane.setHeader(parameters[1]);
	  returnMe = DString("The header for ") + parameters[0] + DString(" now is \"") + parameters[1] + "\"";      
	}
      else
	{
	  DString * currentHeader = editedPlane.getHeader();
	  if(currentHeader == 0)
	    returnMe = DString("There currently is no header for ") + parameters[0];
	  else
	    returnMe = DString("The current header for ") + parameters[0] + DString(" is ") + *currentHeader;
	}
    }
  else if(bean.equals("convolve"))
    {
      returnMe = "Verify that this works (function convolve)";
      /*
      if(storedStuff->getNumberOfKeys() < 2 || numberOfParameters < 3)
	throw DavidException("Two Planes and three Parameters are needed");

      if(!storedStuff->containsKey(parameters[0]))
	throw DavidException(parameters[0] + " could not be found.");

      if(!storedStuff->containsKey(parameters[1]))
	throw DavidException(parameters[1] + " could not be found.");

      plane * left = new plane(storedStuff->get(parameters[0]));
      plane * right = new plane(storedStuff->get(parameters[1]));

      if(left->numberOfRows() != right->numberOfRows() && left->numberOfColumns() != right->numberOfColumns())
	{
	  delete left;
	  delete right;
	  throw DavidException(DString("The dimensions of ") + parameters[0] + DString(" must equal the dimensions of ") + parameters[1]);
	}

	unsigned long * nn = new unsigned long[2];
      
      nn[0] = left->numberOfRows();
      nn[1] = left->numberOfColumns();

      if(nn[0] % 2 != 0 || nn[1] % 2 != 0 )
	{
	  delete left;
	  delete right;
	  delete [] nn;
	  throw DavidException("All dimensions must be a power of 2.");
	}

      std::complex<double> * leftDoubles = PlaneToComplexArray(left);
      std::complex<double> * rightDoubles = PlaneToComplexArray(right);

      FFT::ndim_discrete(leftDoubles,nn,2,1);
      FFT::ndim_discrete(rightDoubles,nn,2,1);

      std::complex<double> * convolve = new complex<double>[nn[0]*nn[1]];

      for(int i = 0;i<nn[0]*nn[1];i++)
	convolve[i] = leftDoubles[i] * rightDoubles[i];

      FFT::ndim_discrete(convolve-1,nn-1,2,-1);

      plane_t * result = ComplexArrayToPlane(convolve,nn[0],nn[1]);
      int norm = nn[0] * nn[1];
      for(int i = 0;i<nn[0];i++)
	for(int j =0;j<nn[1];j++)
	  {
	    Double curr;
	    for(int k = 0;k<3;k++)
	      curr.setValue(k,result->getValue(i,j).getValue(k)/norm);
	    result->setValue(i,j,curr);
	  }

      delete left;
      delete right;
      delete [] nn;
      
      if(storedStuff->containsKey(parameters[2]))
	storedStuff->remove(parameters[2]);
      
      storedStuff->put(parameters[2],*result);
      
      delete result;
      delete [] leftDoubles;
      delete [] rightDoubles;
      delete [] convolve;

      returnMe = DString("The convolution of ") + parameters[0] + DString(" and ") + parameters[1] + DString(" was saved as ") + parameters[2];
      /**/
    }
  else if(bean.equals("sum"))
    {
      if(storedStuff->size() == 0 || numberOfParameters < 1)
	throw DavidException("I need to know what to sum :-/ ");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " does not exist.");

      const plane_t& sumMe = (*storedStuff)[parameters[0]];
      Double sum = 0;
      int * ndims = sumMe.getDimensions();
      
      for(int i = 0;i < ndims[0];i++)
	for(int j = 0;j < ndims[1];j++)
	  sum += sumMe.getValue(i,j);
      
      delete [] ndims;
      
      returnMe = DString("The sum of all elements in ") + parameters[0] + DString(" is ") + sum.toDString();
      
    }
  else if(bean.equals("directproduct"))
    {
      if(storedStuff->size() < 2 || numberOfParameters < 3)
	throw DavidException("I need 2 planes and the name to be given to the product.");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found");

      if(storedStuff->find(parameters[1]) == storedStuff->end())
	throw DavidException(parameters[1] + " was not found");

      const plane_t * left = &(*storedStuff)[parameters[0]];
      const plane_t * right = &(*storedStuff)[parameters[1]];
      
      plane_t * product = 0;
      try
	{
	  std::cerr << "functions.cpp:519 inefficient. should replace with operating on a plane rather than returning one, which is then copied." << std::endl;
	  product = directProduct(left,right);
	}
      catch(DavidException de)
	{
	  if(product != 0)
	    delete product;
	  throw de;//dying gracefully
	}
       
      (*storedStuff)[parameters[2]] = *product;
      if(product)
	delete product;

      returnMe = DString("The direct product of ") + parameters[0] + DString(" and ") + parameters[1] + DString(" was saved as ") + parameters[2];
    }
  else if(bean.equals("multiply"))
    {
      if(storedStuff->size() < 2 || numberOfParameters < 3)
	throw DavidException("I need 2 planes and the name to be given to the product.");
      
      DString A = parameters[0];
      DString B = parameters[1];
      DString Product = parameters[2];

      if(storedStuff->find(A) == storedStuff->end())
	throw DavidException(A + " could not be found.");
      
      if(storedStuff->find(B) == storedStuff->end())
	throw DavidException(B + " could not be found.");


      DEBUG_PRINT("Loading " << A);
      const plane_t& a = (*storedStuff)[A];
      DEBUG_PRINT("Loading " << B);
      const plane_t& b = (*storedStuff)[B];

      if(b.numberOfRows() != a.numberOfColumns())
	{
	  throw DavidException(DString("The number of columns in ") + A + DString(" must equal the number of rows in ") + B);
	}

      DEBUG_PRINT("Multiplying...");
      std::vector<Double> c = Functions::matrixMultiply(a.getPlaneArray(),b.getPlaneArray(),a.numberOfRows(),b.numberOfRows(), b.numberOfColumns());

      plane_t product(c,a.numberOfRows(),b.numberOfColumns());


      (*storedStuff)[Product] = product;

      returnMe = DString("The product of ") + A + DString(" and ") + B + DString(" has been saved as ") + Product;


    }
  else if(bean.equals("add"))
    {
      if(storedStuff->size() < 2 || numberOfParameters < 3)
	throw DavidException("I need 2 planes and the name to be given to the sum.");
      
      DString A = parameters[0];
      DString B = parameters[1];
      DString sumString = parameters[2];

      if(storedStuff->find(A) == storedStuff->end())
	throw DavidException(A + " could not be found.");
      
      if(storedStuff->find(B) == storedStuff->end())
	throw DavidException(B + " could not be found.");


      DEBUG_PRINT("Loading " << A);
      const plane_t * a = &(*storedStuff)[A];
      DEBUG_PRINT("Loading " << B);
      const plane_t * b = &(*storedStuff)[B];

      if(b->numberOfRows() != a->numberOfRows() && b->numberOfColumns() != a->numberOfColumns())
	{
	  throw DavidException(DString("The dimensions of ") + A + DString(" must equal the dimensions of ") + B);
	}

      VERBOSE_PRINT("Adding...");
      std::vector<Double> c = Functions::addMatrices(a->getPlaneArray(), b->getPlaneArray(),a->numberOfRows(),a->numberOfColumns());

      DEBUG_PRINT("Functions.cpp:602 inefficient: don't create it on the heap, then copy it somewhere else :-/");
      plane_t * sum  = new plane_t(c,a->numberOfRows(),b->numberOfColumns());
      
      (*storedStuff)[sumString] = *sum;
      delete sum;
      a = b = sum = 0;

      returnMe = DString("The sum of ") + A + DString(" and ") + B + DString(" has been saved as ") + sumString;


    }
  else if(bean.equals("modulus"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 1)
	throw DavidException("I need a plane.");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found.");

      DString saveName = parameters[0];

      if(numberOfParameters > 1)
	saveName = parameters[1];

      using math::Complex;
      Plane<Complex> * hmmm = convertDoublePlaneToComplexPlane(&(*storedStuff)[parameters[0]]);

      for(int i = 0;i<hmmm->numberOfRows();i++)
	for(int j = 0;j<hmmm->numberOfColumns();j++)
	  {
	    Complex curr(hmmm->getValue(i,j).modulus(),0.0);
	    hmmm->setValue(i,j,curr);
	  }

      plane_t * gimmy = convertComplexPlaneToDoublePlane(hmmm);

      (*storedStuff)[saveName] = *gimmy;

      delete gimmy;
      delete hmmm;
      gimmy = 0;
      hmmm = 0;

      returnMe = DString("The modulus of ") + parameters[0] + DString(" is saved as ") + saveName;

    }
  else if(bean.equals("extract"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 3)
	throw DavidException("I need a plane and I need to know what to do with the plane.");

      DString plane = parameters[1];
      DString saveAs = parameters[numberOfParameters-1];
      
      if(storedStuff->find(plane) == storedStuff->end())
	throw DavidException(plane + " was not found.");

      plane_t verytmpPlane = (*storedStuff)[plane];
      Plane<math::Complex> * tmpPlane = convertDoublePlaneToComplexPlane(&verytmpPlane);
      Plane<math::Complex> * finalResult;
      
      std::vector<DString> doThis;
      
      for(int i = 0;i<numberOfParameters-1;i++)
	if(i != 1)
	  doThis.push_back(parameters[i]);

      try
	{
	  finalResult = extractPortion(tmpPlane,doThis);
	}
      catch(DavidException de)
	{
	  delete tmpPlane;
	  tmpPlane = 0;
	  throw de;
	}

      plane_t * veryfinalResult = convertComplexPlaneToDoublePlane(finalResult);
      (*storedStuff)[saveAs] = *veryfinalResult;
      
      delete veryfinalResult;
      delete finalResult;
      delete tmpPlane;
      tmpPlane = 0;
      finalResult = 0;
      veryfinalResult = 0;
      
      returnMe = doThis.at(0) + DString(" was obtained from " ) + plane + DString(" and saved as ") + saveAs;


    }
  else if(bean.equals("bitmap"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 2)
	throw DavidException("I need a plane to save as a bitmap and the name of the bitmap file.");

      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found.");

      bool useBinary = false;
      int useGrid = 0;
      bool useWhiteBackground = false;

      for(int i = 2;i<numberOfParameters;i++)
	{      
	  if(parameters[i] == "--binarymode")
	    useBinary = true;
	  if(parameters[i] == "--whitebackground")
	    useWhiteBackground = true;
	  if(parameters[i].length() > 8 && parameters[i].substring(0,7) == "--grid=")
	    {
	      DString gridSize =  parameters[i].substring(7);
	      
	      DEBUG_PRINT(DString("Grid size: ") + gridSize);

	      if(Double::isDouble(gridSize))
		{
		  Double val(gridSize);
		  if(val.doubleValue() < 0)
		    throw DavidException("The gridsize must be a positive number!", DavidException::INVALID_ARGUMENT_ERROR_CODE);
		  
		  useGrid = val.toInt();
		}
	    }//end grid process
	}

      if(!useBinary)
	{
	  if(!useWhiteBackground)
	    (*storedStuff)[parameters[0]].draw(parameters[1]);
	  else
	    (*storedStuff)[parameters[0]].draw(parameters[1],false,false,1000,255,255,255);
	}
      else
	{
	  plane_t& poop = (*storedStuff)[parameters[0]];
	  Double red(255,0.0,0.0);
	  for(int i = 0;i<poop.numberOfRows();i++)
	    for(int j = 0;j<poop.numberOfColumns();j++)
	      {
		Double curr = poop.getValue(i,j);
		for(int k = 0;k<3;k++)
		  if(curr.getValue(k) != 0)
		    {
		      poop.setValue(i,j,red);
		    }
	      }

	  //drawing
	  double bgRed,bgGreen,bgBlue;
	  if(useWhiteBackground)
	    bgRed = bgGreen = bgBlue = 255;
	  else
	    bgRed = bgGreen = bgBlue = 0;
	  
	  if(useGrid > 0)
	    poop.draw(fullFilename(parameters[1],currentDirectory),false,true,useGrid,bgRed,bgGreen,bgBlue);
	  else
	    poop.draw(fullFilename(parameters[1],currentDirectory),false,false,1000,bgRed,bgGreen,bgBlue);
	  
	}
      
      returnMe = parameters[0] + DString(" was saved as ") + parameters[1];
      if(useBinary)
	returnMe += ", where every point which is non-zero is a red pixel";
      if(useGrid > 0)
	returnMe += DString(". A grid was used with a space of ")+Double(useGrid).toDString() + " pixels";
      returnMe += ".";
	  
    }
  else if(bean.equals("count"))
    {
       if(storedStuff->size() < 1 || numberOfParameters < 1)
	 throw DavidException("I need a plane.");
      
       plane_t& stuff = (*storedStuff)[parameters[0]];
       

       double totalCount = 0;
       double dZero = 0;
       for(int i = 0;i < stuff.numberOfRows(); i++)
	 for(int j = 0;j < stuff.numberOfColumns();j++)
	   if(stuff.getValue(i,j) != dZero)
	     totalCount++;
  
       returnMe  = DString("The total number of non-zero values in ") + parameters[0] + DString(" is ") + Double(totalCount).toDString();
       
    }
  else if(bean.equals("max") || bean.equals("min"))
    {
       if(storedStuff->size() < 1 || numberOfParameters < 1)
	 throw DavidException("I need a plane.");
       
       bool findMax = bean.equals("max");
       
       double value = 0.0;
       double value2 = 0.0;
       double * coords = new double[2];
       double * coords2 = new double[2];
       coords[0] = coords[1] = 0.0;

       plane_t& stuff = (*storedStuff)[parameters[0]];
       int rows = stuff.numberOfRows();
       int columns = stuff.numberOfColumns();
       for(int i = 0;i<rows;i++)
	 for(int j = 0;j<columns;j++)
	   {

	     double curr = stuff.getValue(i,j).doubleValue();
	     double curr2 = stuff.getValue(i,j).getValue(1);

	     if(i == 0 && j == 0)
	       {
		 value = curr;
		 value2 = curr2;
	       }

	     if(findMax && curr > value)
	       {
		 value = curr;
		 coords[0] = i;
		 coords[1] = j;
	       }
	     if(!findMax && curr < value)
	       {
		 value = curr;
		 coords[0] = i;
		 coords[1] = j;
	       }
	     if(findMax && curr2 > value2)
	       {
		 value2 = curr2;
		 coords2[0] = i;
		 coords2[1] = j;
	       }
	     if(!findMax && curr2 < value2)
	       {
		 value2 = curr2;
		 coords2[0] = i;
		 coords2[1] = j;
	       }
	       

	   }

       returnMe += Double(value).toDString() + DString(" @ (") + Double(coords[0]).toDString()+ DString(",")+Double(coords[1]).toDString()+")";  
       if(value2 != 0.0)
	 returnMe += DString("\n and ") + Double(value2).toDString() + DString("i @ (") + Double(coords2[0]).toDString()+ DString(",")+Double(coords2[1]).toDString()+")";


       delete [] coords;
       delete [] coords2;
       coords = coords2 = 0;
    }
  else if(bean.equals("normalize"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 2)
	throw DavidException("I need a plane and a number");
      
      plane_t& plane = (*storedStuff)[parameters[0]];
      
      double norm = Double(parameters[1]).doubleValue();
      int rows = plane.numberOfRows();
      int columns = plane.numberOfColumns();
      
      for(int i = 0;i<rows;i++)
	for(int j = 0; j<columns;j++)
	  {
	    Double curr;
	    for(int k = 0;k<3;k++)
	      curr.setValue(k,plane.getValue(i,j).getValue(k)/norm);
	    plane.setValue(i,j,curr);
	  }
      
      (*storedStuff)[parameters[0]] = plane;
      
      returnMe = DString("Every value in ") + parameters[0] + DString(" was divided by ") + parameters[1];
    }
  else if(bean.equals("squareplane"))
    {
      if(storedStuff->size() < 1 || numberOfParameters < 1)
	throw DavidException("I need a plane.");
      
      plane_t * plane = &(*storedStuff)[parameters[0]];
      
      int size = -1;
      if(numberOfParameters >= 2)
	size = (int) Double(parameters[1]).doubleValue();
      plane_t * newGuy = Functions::squarePlane(plane,size);

      (*storedStuff)[parameters[0]] = *newGuy;

      returnMe = parameters[0] + DString(" now has the dimensions ") + Double(newGuy->numberOfRows()).toDString() + DString(" by ") + Double(newGuy->numberOfColumns()).toDString() + ".";

      delete newGuy;
      plane = newGuy = 0;
      
    }
  else if(bean.equals("flatten"))
    {
      using utils::DStack;
      if( numberOfParameters < 4)
	throw DavidException("I need to know the file of the cube, resolution, line of sight (0-2) and the variable name (time in seconds is an optional parameter).");
      
      DString filename = parameters[0];
      Double resolution = Double(parameters[1]);
      int los = (int) Double(parameters[2]).doubleValue();
      DString name = parameters[3];
      
      double timeInSeconds = 0;
      if(numberOfParameters >= 5)
	timeInSeconds = Double(parameters[4]).doubleValue();

      DEBUG_PRINT("this los");
      DEBUG_PRINT(los);
      
      Flatten f(false, false, false,false,los,timeInSeconds);

      Plane<utils::DStack<Double> > * plane = f.parseClusterFile(filename, resolution);

      plane_t Dplane(plane->numberOfRows(),plane->numberOfColumns(),0.0);

      int rows = Dplane.numberOfRows();
      int columns = Dplane.numberOfColumns();
      
      double beginningZValue = 9495.09;
      double normalize = 3.16*pow(10.0,8.0);
      double sigma = 1000;//kpc

      if(los == 0)
	beginningZValue = 3565.26;
      else if(los == 1)
	beginningZValue = 55899.2;

      for(int i =0;i<rows;i++)
	for(int j = 0;j<columns;j++)
	  {
	    DStack<Double> stacky = plane->getValue(i,j);
	    Double curr = 0.0;
	    while(!stacky.isEmpty())
	      {
		double value = stacky.pop().doubleValue() - beginningZValue;
		curr += normalize*exp(-value*value/(sigma*sigma));
	      }

	    Dplane.setValue(i,j,curr);
	  }

      (*storedStuff)[name] = Dplane;
            
      delete plane;
      plane = 0;

    }
  else if(bean.equals("histogram"))
    {
      if(numberOfParameters < 2 || storedStuff->size() < 1)
	throw DavidException("I need at least 1 variable and a filename.");

      const plane_t * histoMe = &(*storedStuff)[parameters[0]];

      double bins = 10;
      int scale = 0;

      int gridSize = -1;

      int * colorBarDimensions = 0;
      
      if(numberOfParameters >= 3)
	{
	  for(int i = 2;i<numberOfParameters;i++)
	    {
	      if(Double::isDouble(parameters[i]))
		{
		  bins = Double(parameters[i]).doubleValue();
		}
	      else if(parameters[i] == "--log")
		scale = 1;
	      else if(parameters[i] == "--exp")
		scale = 2;
	      else if(parameters[i] == "--inverse")
		scale = 3;
	      if(parameters[i].length() >= 6)
		{
		  if(parameters[i].substring(0,6) == "--grid")
		    {
		      if(parameters[i].contains('='))
			{
			  DString number = parameters[i].substring(parameters[i].indexOf('=')+1);
			  if(Double::isDouble(number))
			    gridSize = Double(number).intValue();
			}
		    }
		}
	      if(parameters[i].length() >= 14)
		{

		  if(parameters[i].substring(0,11) == "--colorbar=")
		    {
		      if(parameters[i].contains(',') && parameters[i].indexOf(',') < parameters[i].length() - 1)
			{
			  DString widthString = parameters[i].substring(parameters[i].indexOf('=')+1,parameters[i].indexOf(','));
			  DString heightString = parameters[i].substring(parameters[i].indexOf(',') + 1);

			  if(Double::isDouble(widthString) && Double::isDouble(heightString))
			    {
			      delete [] colorBarDimensions;
			      colorBarDimensions = new int[2];
			      colorBarDimensions[0] = Double(widthString).intValue();
			      colorBarDimensions[1] = Double(heightString).intValue();
			      std::cout << ("making colorbar") << std::endl;
			    }
			}
		    }
		}//end colorbar's else if
	    }//end for(i...)
	}//end if(numberOfParameters...)
      
		      

      Utilities::drawHistogram(fullFilename(parameters[1],currentDirectory),histoMe,bins,scale,gridSize, colorBarDimensions);

      returnMe = parameters[0] + DString(" was drawn as ") + parameters[1] + DString(" with ") + Double(bins).toDString() + " bins";

      if(gridSize >= 0)
	returnMe += DString(" using a ") + Double((double) gridSize).toDString() + " pixel grid";
      
      returnMe += ".";

      switch(scale)
	{
	case 1:
	  returnMe += " A log scale was used.";
	  break;
	case 2:
	  returnMe += " An exponential scale was used.";
	  break;
	case 3:
	  returnMe += " A 1/distance scale was used.";
	  break;
	default:
	  break;
	}

      if(colorBarDimensions != 0)
	returnMe += " A colorbar lengend was added.";

      delete [] colorBarDimensions;
      colorBarDimensions = 0;
    }
  else if(bean.equals("plotable"))
    {
      if(numberOfParameters < 2 || storedStuff->size() < 1)
	throw DavidException("I need at least 2 variables, one of which is a stored plane.");

      
      const plane_t * printMe = &(*storedStuff)[parameters[0]];
      bool writeParam = Functions::PARAMETER_USE_SM;
      
      if(numberOfParameters >= 3)
	if(parameters[2].equals("--flip_coordinates"))
	  writeParam = Functions::PARAMETER_FLIP_COORDINATES;

      writePlotableData(printMe,fullFilename(parameters[1],currentDirectory),writeParam);

      returnMe = parameters[0] + DString(" has been written as ") + parameters[1];
      if(writeParam == Functions::PARAMETER_FLIP_COORDINATES)
	returnMe += " and the coordinates were flipped.";
      
      
    }
  else if(bean.equals("floor"))
    {
      if(numberOfParameters < 2 || storedStuff->size() < 1)
	throw DavidException("I need at least a plane and a number");
      
      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found.");
      
      const plane_t * plane = &(*storedStuff)[parameters[0]];
      Double DFloor(parameters[1]);

      bool inverseFloor = false;
      bool useAbs = false;
      if(numberOfParameters >= 3)
	{
	  useAbs = useAbs || parameters[2].equals("--abs");
	  inverseFloor = inverseFloor || parameters[2].equals("--inverse");
	}
      if(numberOfParameters >= 4)
	{
	  useAbs = useAbs || parameters[3].equals("--abs");
	  inverseFloor = inverseFloor || parameters[3].equals("--inverse");
	}

      plane_t * newplane = floorMe(plane,DFloor.doubleValue(),useAbs,inverseFloor);
      
      (*storedStuff)[parameters[0]] = *newplane;
      
      delete newplane;
      newplane = 0;

      returnMe += "Floor done (new comment someday?)";
      
    }
  else if(bean.equals("conjugate"))
    {
      if(numberOfParameters < 1 || storedStuff->size() < 1)
	throw DavidException("I need the name of the matrix to take the conjuage of");
      
      if(storedStuff->find(parameters[0]) == storedStuff->end())
	throw DavidException(parameters[0] + " was not found.");

      DString newGuy = parameters[0];

      if(numberOfParameters > 1)
	newGuy = parameters[1];

      using math::Complex;
      const plane_t& stuffedStuff = (*storedStuff)[parameters[0]];
      Plane<Complex> * original = convertDoublePlaneToComplexPlane(&stuffedStuff);

      Plane<Complex> * newPlane = complexConjugate(original);
      plane_t * saveMe = convertComplexPlaneToDoublePlane(newPlane);

      (*storedStuff)[newGuy] = *saveMe;
      

      delete original;
      delete newPlane;
      delete saveMe;
      original = newPlane = 0;
      saveMe = 0;

      returnMe = DString("The complex conjuage of ") + parameters[0] + " was saved";
      if(newGuy != parameters[0])
	returnMe += DString(" as ") + newGuy;

      returnMe += ".";

      
    }
  else if(bean.equals("dir"))
    {
      DString dir;
      
      if(currentDirectory == 0)
	currentDirectory = new DString("./");
      
      dir = *currentDirectory;

      if(numberOfParameters >= 1)
	dir = parameters[0];

      DIR * currDir;
      dirent * theEntry;
      
      currDir = opendir(dir.toCharArray());

      if(currDir == 0)
	throw DavidException("Error in viewing directory",DavidException::IO_ERROR_CODE);

      returnMe = DString("Contents of ")+dir;

      while(theEntry = readdir(currDir))
	{
	  returnMe += DString(" ") + theEntry->d_name;
	}
      
      closedir(currDir);
      

    }
  else if(bean.equals("chdir"))
    {
      

      if(currentDirectory == 0)
	currentDirectory = new DString("./");

      if(numberOfParameters < 1)
	{
	  DString curr = *currentDirectory;
	  return DString("Current Directory: ") + curr;
	}
      
      DString lastDir = *currentDirectory;
      
      *currentDirectory = parameters[0];
      
      if(currentDirectory->charAt(currentDirectory->length() -1) != '/')
	*currentDirectory += "/";
      
      returnMe = DString("Current Directory: " )+ *currentDirectory + DString("\nPrevious directory: ") + lastDir;
     
    }
 

  //end of do function else business


return returnMe;
  
}/**/
	
DString Functions::listFunctions()    
{
  using namespace std;
  
  FunctionHelp fh;
  map<string,string> * functions = fh.getHelp();
  map<string,string>::const_iterator function = functions->begin();

  DString returnMe = "";
  DString space = " ";
  for(;function != functions->end();function++)
    returnMe += space+function->first;

  return returnMe;
}

DString Functions::getHelp(const DString& bean)
{
  using namespace std;
  if(bean == "list")
    {
      DString returnMe("List of functions: ");
      returnMe += Functions::listFunctions();
      return returnMe;
    }
  if(!isFunction(bean))
    return "No such function.";
  FunctionHelp fh;
  map<string,string> * hm = fh.getHelp();//do not delete will be deleted by ~fh
  if(hm->find(bean) == hm->end())
    return DString("No help for ")+bean;

  return (*hm)[bean];
}


std::complex<double> * Functions::PlaneToComplexArray(plane_t * dPlane)
{

  int rows = dPlane->numberOfRows();
  int columns = dPlane->numberOfColumns();
  
  std::complex<double> * cplane = new std::complex<double>[rows*columns];
  int counter = 0;
  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {
	cplane[counter] = std::complex<double>(dPlane->getValue(i,j).doubleValue(),dPlane->getValue(i,j).getValue(1));
	counter++;
      }

  return cplane;
}


  plane_t *  Functions::ComplexArrayToPlane(std::complex<double> * cPlane, int rows, int columns, double norm)
{

  plane_t * dPlane = new plane_t(rows,columns, 0.0);
  int counter = 0;
  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {
	  Double curr(cPlane[counter].real()/norm,cPlane[counter].imag()/norm,0);
	  dPlane->setValue(i,j,curr);
	  counter++;
      }

  return dPlane;

}


std::vector<Double> Functions::matrixMultiply(std::vector<Double> a, std::vector<Double> b, int I, int J, int K)
{
  std::vector<Double> returnMe(I*J,0.0);
  Double tmp;
  for(int i = 0;i<I;i++)
    {
      for(int j = 0;j<J;j++)
	{
	  tmp = returnMe[j+i*J];
	  
	  for(int k = 0;k<K;k++)
	    {
	      tmp += a[k+i*K]*b[j+k*J];
	    }
	}
    }

  return returnMe;
  

}

plane_t * Functions::directProduct(const plane_t * a, const plane_t * b)
{
  
  using math::Complex;
  
  int I = a->numberOfRows();
  int J = a->numberOfColumns();

  if(I != b->numberOfRows() && J != b->numberOfColumns())
    throw DavidException("The dimensions of both planes must be equal.", DavidException::PLANE_OUT_OF_BOUNDS_ERROR_CODE);

  plane_t * returnMe = new plane_t(I,J,0.0);

  for(int i = 0;i<I;i++)
    for(int j = 0;j<J;j++)
      {
	Complex left(a->getValue(i,j));
	Complex right(b->getValue(i,j));
	Double product = (Double) (left * right);
	returnMe->setValue(i,j,product);
      }

  return returnMe;
}


plane_t *  Functions::squarePlane(const plane_t * oldPlane, const int& newWidth)
{
  
  Double cZero(0.0);

  int bounds;
  
  if(oldPlane->numberOfRows() > oldPlane->numberOfColumns())
    bounds = oldPlane->numberOfRows();
  else
    bounds = oldPlane->numberOfColumns();

  if(newWidth > 0)
    bounds = newWidth;

  plane_t * newPlane = new plane_t(bounds,bounds,cZero);

  for(int i = 0;i < oldPlane->numberOfRows();i++)
    for(int j = 0; j < oldPlane->numberOfColumns();j++)
      {
	Double curr = oldPlane->getValue(i,j);
	newPlane->setValue(i,j,curr);
      }

  return newPlane;
}


void Functions::writePlotableData(plane_t * plane, DString fileName, int writeParam)
{

  int rows,columns;

  switch(writeParam)
    {
    case Functions::PARAMETER_USE_SM:
      {
	rows = plane->numberOfRows();
	columns = plane->numberOfColumns();
	break;
      }
    case Functions::PARAMETER_FLIP_COORDINATES: default:
      {
	rows = plane->numberOfColumns();
	columns = plane->numberOfRows();
	break;
      }
    }  


  using namespace std;

  fstream file(fileName.toCharArray(),ios::out);

  if(!file.is_open())
    throw DavidException(DString("Cannot open") + fileName, DavidException::IO_ERROR_CODE);

  for(int i = 0 ; i < rows ; i++)
    {
      double x = (i-rows/2);
      for(int j = 0; j<columns;j++)
	{
	
	  double y =(j-columns/2);
	  
	  switch(writeParam)
	    {
	    case Functions::PARAMETER_USE_SM:
	      {
		y = -1*y-1;
		break;
	      }
	    case Functions::PARAMETER_FLIP_COORDINATES: default:
	      break;
	    }
	
	  file << "     " << x << "     " << y << "     ";
	  
	  switch(writeParam)
	    {
	    case Functions::PARAMETER_USE_SM:
	      {
		file << math::Complex(plane->getValue(i,j)).toDString() << endl;
		break;
	      }
	    case Functions::PARAMETER_FLIP_COORDINATES: default:
	      {
		file << math::Complex(plane->getValue(j,i)).toDString() << endl;
		break;
	      }
	    }
	}
    }

  file.close();
  
}


plane_t * Functions::convertComplexPlaneToDoublePlane(const Plane<math::Complex> * oldPlane)
{

  int rows = oldPlane->numberOfRows();
  int columns = oldPlane->numberOfColumns();

  plane_t * returnMe = new plane_t(rows,columns, 0.0);

  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {
	math::Complex cCurr = oldPlane->getValue(i,j);
	Double curr(cCurr.getRealPart(),cCurr.getImaginaryPart(),0.0);
	returnMe->setValue(i,j,curr);
      }

  if(oldPlane->getHeader() != 0)
    returnMe->setHeader(*(oldPlane->getHeader()));

  return returnMe;

}

Plane<math::Complex> * Functions::convertDoublePlaneToComplexPlane(const plane_t * oldPlane)
{

  int rows = oldPlane->numberOfRows();
  int columns = oldPlane->numberOfColumns();

  using math::Complex;

  Complex CZero(0.0,0.0);
  
  Plane<Complex> * returnMe = new Plane<Complex>(rows,columns,CZero);

  for(int i = 0;i<rows;i++)
    for(int j = 0;j<columns;j++)
      {
	Double cCurr = oldPlane->getValue(i,j);
	Complex curr(cCurr.getValue(0),cCurr.getValue(1));
	returnMe->setValue(i,j,curr);
      }
  
  if(oldPlane->getHeader() != 0)
    returnMe->setHeader(*(oldPlane->getHeader()));
  
  return returnMe;
}

plane_t * Functions::floorMe(plane_t * original, double floor, bool useAbs, bool inverseFloor)
{


  int rows = original->numberOfRows();
  int columns = original->numberOfColumns();
  
  plane_t * returnMe = new plane_t(rows,columns,0.0);

  for(int i = 0;i<rows;i++)
    for(int j = 0; j<columns;j++)
      {
	Double curr = original->getValue(i,j);
	
	for(int k = 0;k<3;k++)
	  {
	    if(useAbs)
	      {
		if(abs(curr.getValue(k)) > floor && !inverseFloor)
		  curr.setValue(k,0.0);
		else if(abs(curr.getValue(k)) < floor && inverseFloor)
		  curr.setValue(k,0.0);
	      }
	    else
	      {
		if(curr.getValue(k) > floor && !inverseFloor)
		  curr.setValue(k,0.0);
		else if(curr.getValue(k) < floor && inverseFloor)
		  curr.setValue(k,0.0);
	      }
	  }
	
	returnMe->setValue(i,j,curr);
      }

  return returnMe;

}

Plane<math::Complex> * Functions::complexConjugate(Plane<math::Complex> * original)
{

  int rows = original->numberOfRows();
  int columns = original->numberOfColumns();
  
  using math::Complex;
  Complex c(0.0,0.0);
  
  Plane<math::Complex> * returnMe = new Plane<math::Complex>(rows,columns,c);
  
  for(int i = 0;i<rows;i++)
    for(int j =0;j<columns;j++)
      {
	Complex curr = original->getValue(i,j).getConjugate();
	returnMe->setValue(i,j,curr);
      }

  return returnMe;

}


std::vector<Double> Functions::addMatrices(std::vector<Double> A,std::vector<Double> B,int rows,int columns) throw (DavidException)
{

  std::vector<Double> sum(rows*columns,0.0);

  for(int i = 0;i<rows*columns;i++)
    sum[i] = A[i] + B[i];

  return sum;

}


Plane<math::Complex> * Functions::getRealPart(const Plane<math::Complex> * original)
{
  using math::Complex;
  Complex cZero(0.0,0.0);
  Plane<Complex> * returnMe = new Plane<Complex>(original->numberOfRows(),original->numberOfColumns(),cZero);
  
  for(int i = 0;i<returnMe->numberOfRows();i++)
    for(int j = 0;j<returnMe->numberOfColumns();j++)
      {
	Complex tmp(original->getValue(i,j).getRealPart(),0.0);
	returnMe->setValue(i,j,tmp);
      }

  return returnMe;

}

Plane<math::Complex> * Functions::getImaginaryPart(const Plane<math::Complex> * original)
{
  using math::Complex;
  Complex cZero(0.0,0.0);
  Plane<Complex> * returnMe = new Plane<Complex>(original->numberOfRows(),original->numberOfColumns(),cZero);
  
  for(int i = 0;i<returnMe->numberOfRows();i++)
    for(int j = 0;j<returnMe->numberOfColumns();j++)
      {
	Complex tmp(original->getValue(i,j).getImaginaryPart(),0.0);
	returnMe->setValue(i,j,tmp);
      }

  return returnMe;

}

Plane<math::Complex> * Functions::extractPortion(const Plane<math::Complex> * original, const std::vector<DString>& parameters) throw (DavidException)
{

  if(parameters.size() == 0)
    throw DavidException("Look I need to know what to do here.", DavidException::INVALID_ARGUMENT_ERROR_CODE);
  
  Plane<math::Complex> * finalResult = 0;

  DString whatToDo = parameters.at(0);

  if(whatToDo == "real")
    {
      finalResult = getRealPart(original);
    }
  else if(whatToDo == "imaginary")
    {
      finalResult = getImaginaryPart(original);
    }
  else if(whatToDo == "subset")
    {
      if(parameters.size() < 5)
	throw DavidException("I need to know the first and last row and columns.",DavidException::INVALID_ARGUMENT_ERROR_CODE);
      
      int * bounds = new int[4];

      
      
      for(int i = 0;i<4;i++)
	bounds[i] = Double(parameters.at(i+1)).intValue();

      
      if(bounds[0] < 0 || bounds[2] < 0)
	throw DavidException("The lowest row or column number is 0.",DavidException::INVALID_ARGUMENT_ERROR_CODE);

      if(bounds[1] > original->numberOfRows())
	throw DavidException(DString("The lowest row is ") + Double(1.0*original->numberOfRows()),DavidException::INVALID_ARGUMENT_ERROR_CODE);
      
      if(bounds[3] > original->numberOfColumns())
	throw DavidException(DString("The lowest column number is ") + Double(1.0*original->numberOfColumns()),DavidException::INVALID_ARGUMENT_ERROR_CODE);
      

      math::Complex cZero(0.0,0.0);
      finalResult = new Plane<math::Complex>(bounds[1]-bounds[0]+1,bounds[3]-bounds[2]+1,cZero);

      for(int i = bounds[0];i<bounds[1]+1;i++)
	for(int j = bounds[2];j<bounds[3]+1;j++)
	  {
	    math::Complex curr = original->getValue(i,j);
	    finalResult->setValue(i-bounds[0],j-bounds[2],curr);
	  }
      delete [] bounds;
    }
  else
    throw DavidException(whatToDo + " is not a valid thing to do.");

  return finalResult;
}


DString Functions::fullFilename(DString fileName, DString * currentDirectory)
{
  if(fileName.charAt(0) == '/')
    return fileName;

  if(currentDirectory == 0)
    currentDirectory = new DString("./");
  
  return *currentDirectory + fileName;
}
