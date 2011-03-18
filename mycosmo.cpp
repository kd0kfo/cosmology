/** \Mainpage Assembler for the compilation of programs
 * for the PIC microcontroller.
 *
 * Created By David Coss, 2010
 * Free Software under the terms of the GNU Public License version 3.
 */
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>

#include "globals.h"

//this prevents the GUI from being used.
#include "libdnstd/StringTokenizer.h"
#include "Parser.h"

#include "Build.h"
#include "args_parser.cpp"

using namespace std;

class DavidException;

void init(Object *obj)
{
  if(obj == 0)
    return;
  
  obj->_string = 0;
  obj->_int = 0;
  obj->_float = 0;
}


ios::fmtflags radix;
map<arg_t,int> equs;
map<arg_t,Object> free_store;//maps variable names to the number of elements in its array. addr_t = 1 for a single variable (instead of array); this is the default.
vector<string> precompiledCode;

std::string print_function(const Command& arg)
{
  std::cout << arg.getCommandWord() << std::endl;
  return arg.getCommandWord();
}

Help makeHelp(){
  Help returnMe = mycosmoHelp();
  returnMe["malloc"] = "Allocates memory for the specified type. Can optionally supply a number of variable spaces to allocate. Default type: ushort\n Possible Types: plane<type>(NxM) (NxM Plane<type>), float (2bytes), char(1byte), ushort(1byte)";
  return returnMe;
}//for tab complete
Help readline_help_function(){return mycosmoHelp();}

string formatHex(const int& unformatted,ios::fmtflags currRadix)
{
  ostringstream os;
  os.precision(2);
  os.setf(currRadix,ios::basefield);
  os << unformatted;
  string returnMe = os.str();
  while(returnMe.size() < 2)
    returnMe = "0" + returnMe;
  return returnMe;
}

string formatHex(const string& unformatted)
{
  string semiformatted = unformatted;
  istringstream is;
  ios::fmtflags currRadix = radix;

  //check if format of the string overrides the default radix.
  if(semiformatted.substr(0,2) == "0x" || semiformatted.substr(0,2) == "0X")
    {
      currRadix = ios::hex;//because we've already formatted it into hex. So the "decimals" are literally in hex.
      semiformatted.erase(0,2);
    }
  else
    {
      switch(semiformatted.at(0))
	{
	case '.': case 'd': case 'D':
	  {
	    currRadix = ios::dec;
	    semiformatted.erase(0,1);
	    break;
	  }
	case 'o': case 'O':
	  {
	    currRadix = ios::oct;
	    semiformatted.erase(0,1);
	    break;
	  }
	default:
	  break;
	}
    }
  is.clear();is.str(semiformatted);is.setf(currRadix,ios::basefield);
  int val = -1;
  is >> val;
  
  if(is.fail())
    throw DavidException("Could not compile: " + unformatted);
  
  return formatHex(val,ios::hex);
}


template <class T> T strToType(const string& arg)
{
  istringstream buff(arg);
  T val;
  buff >> val;
  if(buff.fail())
    throw DavidException("strToType: Could not string convert to requested type: " + arg);
  return val;
}
template int strToType(const string&);
template float strToType(const string&);



int main(int argc, char **argv)
{
  radix = ios::hex;
  //so that there aren't multiple main's
  #define __HAVE_MAIN__ 1

  struct assembler_arguments args;
  args.source_filename = "";
  args.output_filename = "a.mycosmo";

  parse_args(argc,argv,args);
  
  if(argc > 1)
	{
		try{
		  Parser p;
		  p.program_name = argv[0];
                  p.setHelp(makeHelp());
		  p.main(argc,argv);
		  return 0;
		}
		catch(DavidException de)
		{
		  de.stdOut();
		  return de.getCode();
		}
	}
	else{
		try{
			Parser p;
			p.program_name = argv[0];
			p.setHelp(mycosmoHelp());
			p.main();
			return 0;
		}
		catch(DavidException de)
		{
		  de.stdOut();
		  return de.getCode();
		}
	}

}

