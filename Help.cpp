#ifndef HELP_CPP
#define HELP_CPP

#include "DHashMap.h"

	class Help
	{
		
		public: 
			Help();
			~Help();
			
			static DString getVersion(){return "1.1";}
			DString getHelp(const DString&);
			DString getHelp(const char * bean){return getHelp(DString(bean));}
			

		private: 
			DHashMap<DString> * contents;
	};
    
Help::Help()
{
  contents = new DHashMap<DString>;
  contents->put("list", DString("Usage: list\nLists all variables defined."));
  contents->put("define",DString("Usage: defines x y\n Defines x to be the expression y."));
  contents->put("det", DString("Usage: det A\nGives the determinate of A."));
  contents->put("load", DString("Usage: load <filename> x\nLoads the saved Variable as x."));
  contents->put("print",DString("Usage: print x\nPrints x. \n The optional argument of \"--real\" uses double values to save and is much faster. This should only be used if the array elements are REAL."));
  contents->put("save",DString("Usage: save x <filename>\nSaves x as <filename>.\n The optional argument of \"--real\" uses double values to save and is much faster. This should only be used if the array elements are REAL."));
  contents->put("remove",DString("Usage: remove x\nRemoves x."));
  contents->put("rename",DString("Usage: rename <old> <new>\n     Renames the plane."));
  contents->put("edit",DString("Usage: edit A m n a,b,...\nEdits the A matrix of size mxn as elements a,b,c,.."));
  contents->put("help",DString("Usage: help <topic>\nOffers help on the topic"));
  contents->put("last",DString("Usage: last\nPrints the last expression used. The last expression is represented by the variable $last$"));
  contents->put("dir",DString("Usage: dir [directory]\nlists the contents of the directory."));
  contents->put("quit",DString("Usage: quit\nExits linal."));
  contents->put("exit", contents->get("quit"));
  contents->put("physics",DString("Usage: physics <equation> parameter\nGives the output of <equations with the give parameters."));
  contents->put("import",DString("Usage: import\nimports constants."));
  contents->put("set",DString("Usage: set parameter <value>\nSets the parameter to the value <value>."));
  contents->put("run",DString("Usage: run scriptname\n Runs the script of commands which must end with a semicolon."));
  contents->put("functions", DString("List of Functions: "));
  contents->put("about", DString("Gives a description of the given variable."));
  contents->put("history",DString("Give the commands that have been executed."));
}

DString Help::getHelp(const DString& keyWord)
{
  DString returnMe = DString("No help for ")+keyWord+".";

  if(contents->containsKey(keyWord))
    returnMe =  contents->get(keyWord);

  return returnMe;
}

Help::~Help()
{
  
  delete contents;
  contents = 0;
}


#endif


