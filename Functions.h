#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include <dirent.h>//Directory entries
#include <vector>
#include <complex>

#include "libdnstd/DHashMap.h"
#include "libdnstd/Double.h"
#include "libdnstd/DArray.h"
#include "libdnstd/StringTokenizer.h"
#include "libdnstd/DStack.h"
#include "libdnstd/Complex.h"
#include "libdnstd/fft.h"

#include "libmygl/Cosmology.h"
#include "libmygl/planecreator.h"


#include "flatten.h"
#include "utilities.h"

	class Functions
	{
	  
	public:
	  Functions(const DString& bean);
		static bool isFunction(const DString& bean);
		static DString doFunction(const DString& bean, DString * parameters,int numberOfParameters,double * gparameters, int numberOfGparameters, DHashMap<Plane<Double> > * storedStuff, DString * currentDirectory);
		static DString listFunctions();
		static DString getHelp(const DString&);
		static std::complex<double> * PlaneToComplexArray(Plane<Double> * dPlane);
		static Plane<Double> *  ComplexArrayToPlane(std::complex<double> * cPlane, int rows, int columns, double normalizationConst = 1);
		static std::vector<Double> matrixMultiply(std::vector<Double> a, std::vector<Double> b, int I, int J, int K);
		static Plane<Double> * directProduct(Plane<Double> * a, Plane<Double> * b);
		static std::vector<Double> addMatrices(std::vector<Double> A,std::vector<Double> B,int rows,int columns) throw (DavidException);
		static Plane<Double> *  squarePlane(Plane<Double> * oldPlane,int newWidth = -1);
		static void writePlotableData(Plane<Double> * plane, DString fileName, int writeParameter);
		static Plane<Double> * convertComplexPlaneToDoublePlane(Plane<math::Complex> * oldPlane);
		static Plane<math::Complex> * convertDoublePlaneToComplexPlane(Plane<Double> * oldPlane);
		static utils::DArray<DString> * getFunctions();
		static Plane<Double> * floorMe(Plane<Double> * original, double floor, bool useAbs = false, bool inverseFloor = false);
		static Plane<math::Complex> * complexConjugate(Plane<math::Complex> * original);
		static Plane<math::Complex> * getRealPart(const Plane<math::Complex> * original);
		static Plane<math::Complex> * getImaginaryPart(const Plane<math::Complex> * original);
		static Plane<math::Complex> * extractPortion(const Plane<math::Complex> * original, utils::DArray<DString>& parameters) throw (DavidException);
		
		/**
		 * Creates an absolute file name.
		 */
		static DString fullFilename(DString fileName, DString * currentDirectory);


		static const int PARAMETER_USE_SM = 1;
		static const int PARAMETER_FLIP_COORDINATES = 2;
	private:
		DString * functionType;
	};
	

class FunctionHelp
{
      
	private :
		DHashMap<DString> * help;

	public:
		FunctionHelp(){
		  help = new DHashMap<DString>;
		  help->put("fixdiag","Fixes the diagonal problem with fourier transform of shear");
		  help->put("r2d",DString("Gives the angular diameter distance to a given redshift"));
		  help->put("fourier",DString("Fourier Transforms the given matrix"));
		  help->put("shear",DString("Usage: shear <convergence> <kernel> <output name>\n     Calculates the shear for a given convergence."));
		  help->put("multiply",DString("Multiplies two matrices and saves the result under the given variable."));
		  help->put("add",DString("Adds two matrices and saves the result under the given variable."));
		  help->put("bitmap",DString("Creates a bitmap from the given plane. The optional \"--binarymode\" argument puts a red pixel everywhere the value in the matrix is nonzero.\nThe optional \"--whitebackground\" draws the bitmap with a white background (meaning 0 values are white pixels)."));
		  help->put("max",DString("Gives the maximum element in a plane."));
		  help->put("count",DString("Gives the count of non-zero elements in a plane"));
		  help->put("min",DString("Gives the minimum element in a plane."));
		  help->put("crit_density",DString("Calculates the Critical Density (g/cm^2) for given redshifts\nUsage:\ncrit_density z1 z2"));
		  help->put("extract","Gets the specified part of the plane (ie Real, Imaginary, etc)\nUsage: extract <what to do> <plane> [optional parameters]   <result plane>\nThings to be done: real, imaginary, subset ");
		  help->put("normalize",DString("Divides every element in a plane by the given number.\n Usage: normalize plane number"));
		  help->put("squareplane",DString("Turns the plane into a square plane by adding zero to the smaller dimension.\nUsage: squareplane <plane name> <new width/height>"));
		  help->put("flatten",DString("Flattens a cube to a 2D plane along the given line of sight"));
		  help->put("plotable",DString("Writes the data in a format that can be plotted in another program. An additional flag \"--flip_coordinates\" can be added to  flip the coordinates.\n Usage: plotable <variable name> <outputfile> [options]"));
		  help->put("floor",DString("Sets all values beyond the given value equal to zero. Adding the --abs tag will use the absolute value to determine whether or not to drop the value. Adding the --inverse tag will take all values below the given number instead of above."));
		  help->put("convolve",DString("Convolves 2 matrices"));
		  help->put("directproduct",DString("Multiplies 2 matrices element by element (as opposed to \"multiply\" which uses standard matrix multiplication"));
		  help->put("conjugate",DString("Takes the complex conjuage of the elements in the given matrix"));
		  help->put("modulus",DString("Calculates the modulus of each element in a matrix."));
		  help->put("histogram",DString("Creates a crude histogram. There are flags for optional scales (ie log) as well.\nUsage: histogram <variable> <filename> [#bins] [--log] [--exp]"));
		  help->put("dir",DString("List the contents of the given directory. If no directory is provided, the current director is listed."));
		  help->put("addellipseinfo","Adds ellipse information to the plane the given pixel. Note: if no ellipticity is given, a circle is drawn.\nUsage: addellipseinfo <plane> <x> <y> <semimajoraxislength> [ellipticity = 0] ");
		  help->put("new","Creates a new empty plane.\nUsage: new <name> <rows> <columns> [value = 0]");
		  help->put("drawellipses","Draws an ellipse diagram\nUsage: <plane> <outputplane>");
		  help->put("chdir","Changes the current directory. If no directory is give, then the current directory is displayed");
		  help->put("sum","Gives the sum of the elements of a matrix.");
		  help->put("flipmatrix","Exhanges matrix elements based on the given method.\nUsage:flipmatrix <matrix> <option> [new matrix name]\n<matrix> -- Name of Matrix to use.\n<option>: transpose, horizontal, vertical\n[new matrix new] -- optional matrix to store result in. If no name is provided, the original is used.");
		  help->put("header","Sets the value of the plane header. If no header is provided, the current header is displayed.\nUsage: header <plane> [new header]");
		}

		~FunctionHelp(){delete help;help = 0;}
        
        DHashMap<DString> * getHelp()
        {
            return help;
        }
    };



#endif

