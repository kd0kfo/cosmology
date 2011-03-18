#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP

#include <dirent.h>//Directory entries
#include <string>
#include <vector>
#include <set>
#include "libmygl/plane.h"
#include <map>
#include <complex>

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

#ifndef plane_t
typedef Plane<Double> plane_t;
#endif

	class Functions
	{
	  
	public:
	  Functions(const DString& bean);
		static bool isFunction(const DString& bean);
		static DString doFunction(const DString& bean, DString * parameters,int numberOfParameters,double * gparameters, int numberOfGparameters, std::map<std::string, plane_t> * storedStuff, DString * currentDirectory);
		static DString listFunctions();
		static DString getHelp(const DString&);
		static std::complex<double> * PlaneToComplexArray(plane_t * dPlane);
		static plane_t *  ComplexArrayToPlane(std::complex<double> * cPlane, int rows, int columns, double normalizationConst = 1);
		static std::vector<Double> matrixMultiply(std::vector<Double> a, std::vector<Double> b, int I, int J, int K);
		static plane_t * directProduct(const plane_t * a, const plane_t * b);
		static std::vector<Double> addMatrices(std::vector<Double> A,std::vector<Double> B,int rows,int columns) throw (DavidException);
		static plane_t *  squarePlane(const plane_t * oldPlane,const int& newWidth = -1);
		static void writePlotableData(plane_t * plane, DString fileName, int writeParameter);
		static plane_t * convertComplexPlaneToDoublePlane(const Plane<math::Complex> * oldPlane);
		static Plane<math::Complex> * convertDoublePlaneToComplexPlane(const plane_t * oldPlane);
		static std::set<std::string> * getFunctions();
		static plane_t * floorMe(plane_t * original, double floor, bool useAbs = false, bool inverseFloor = false);
		static Plane<math::Complex> * complexConjugate(Plane<math::Complex> * original);
		static Plane<math::Complex> * getRealPart(const Plane<math::Complex> * original);
		static Plane<math::Complex> * getImaginaryPart(const Plane<math::Complex> * original);
		static Plane<math::Complex> * extractPortion(const Plane<math::Complex> * original, const std::vector<DString>& parameters) throw (DavidException);
		
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
  std::map<std::string,std::string> * help;

 public:
		FunctionHelp(){
		  help = new std::map<std::string,std::string>;
		  (*help)["fixdiag"] = "Fixes the diagonal problem with fourier transform of shear";
		  (*help)["r2d"] = "Gives the angular diameter distance to a given redshift";
		  (*help)["fourier"] = "Fourier Transforms the given matrix";
		  (*help)["shear"] = "Usage: shear <convergence> <kernel> <output name>\n     Calculates the shear for a given convergence.";
		  (*help)["multiply"] = "Multiplies two matrices and saves the result under the given variable.";
		  (*help)["add"] = "Adds two matrices and saves the result under the given variable.";
		  (*help)["bitmap"] = "Creates a bitmap from the given plane. The optional \"--binarymode\" argument puts a red pixel everywhere the value in the matrix is nonzero.\nThe optional \"--whitebackground\" draws the bitmap with a white background (meaning 0 values are white pixels).";
		  (*help)["max"] = "Gives the maximum element in a plane.";
		  (*help)["count"] = "Gives the count of non-zero elements in a plane";
		  (*help)["min"] = "Gives the minimum element in a plane.";
		  (*help)["crit_density"] = "Calculates the Critical Density (g/cm^2) for given redshifts\nUsage:\ncrit_density z1 z2";
		  (*help)["extract"] = "Gets the specified part of the plane (ie Real, Imaginary, etc)\nUsage: extract <what to do> <plane> [optional parameters]   <result plane>\nThings to be done: real, imaginary, subset ";
		  (*help)["normalize"] = "Divides every element in a plane by the given number.\n Usage: normalize plane number";
		  (*help)["squareplane"] = "Turns the plane into a square plane by adding zero to the smaller dimension.\nUsage: squareplane <plane name> <new width/height>";
		  (*help)["flatten"] = "Flattens a cube to a 2D plane along the given line of sight";
		  (*help)["plotable"] = "Writes the data in a format that can be plotted in another program. An additional flag \"--flip_coordinates\" can be added to  flip the coordinates.\n Usage: plotable <variable name> <outputfile> [options]";
		  (*help)["floor"] = "Sets all values beyond the given value equal to zero. Adding the --abs tag will use the absolute value to determine whether or not to drop the value. Adding the --inverse tag will take all values below the given number instead of above.";
		  (*help)["convolve"] = "Convolves 2 matrices";
		  (*help)["directproduct"] = "Multiplies 2 matrices element by element (as opposed to \"multiply\" which uses standard matrix multiplication";
		  (*help)["conjugate"] = "Takes the complex conjuage of the elements in the given matrix";
		  (*help)["modulus"] = "Calculates the modulus of each element in a matrix.";
		  (*help)["histogram"] = "Creates a crude histogram. There are flags for optional scales (ie log) as well.\nUsage: histogram <variable> <filename> [#bins] [--log] [--exp]";
		  (*help)["dir"] = "List the contents of the given directory. If no directory is provided, the current director is listed.";
		  (*help)["addellipseinfo"] = "Adds ellipse information to the plane the given pixel. Note: if no ellipticity is given, a circle is drawn.\nUsage: addellipseinfo <plane> <x> <y> <semimajoraxislength> [ellipticity = 0] ";
		  (*help)["new"] = "Creates a new empty plane.\nUsage: new <name> <rows> <columns> [value = 0]";
		  (*help)["drawellipses"] = "Draws an ellipse diagram\nUsage: <plane> <outputplane>";
		  (*help)["chdir"] = "Changes the current directory. If no directory is give, then the current directory is displayed";
		  (*help)["sum"] = "Gives the sum of the elements of a matrix.";
		  (*help)["flipmatrix"] = "Exhanges matrix elements based on the given method.\nUsage:flipmatrix <matrix> <option> [new matrix name]\n<matrix> -- Name of Matrix to use.\n<option>: transpose, horizontal, vertical\n[new matrix new] -- optional matrix to store result in. If no name is provided, the original is used.";
		  (*help)["header"] = "Sets the value of the plane header. If no header is provided, the current header is displayed.\nUsage: header <plane> [new header]";
		}

		~FunctionHelp(){delete help;help = 0;}
        
		const std::map<std::string,std::string> * getHelp()
		  {
		    return help;
		  }
};



#endif

