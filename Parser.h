#ifndef PARSER_CPP
#define PARSER_CPP

#include "Command.h"
#include "CommandWords.h"
#include "Help.cpp"
#include "Functions.h"
#include "DArray.h"

#include "planecreator.h"

class David{
 public:
  static void say(const char * bean){ std::cout << bean << std::endl;}
  static void say(const DString& bean){say(bean.toCharArray());}
};

class Parser
{
 public:
#ifdef __DEBUG__
  int DEBUG_COUNTER;
#endif

  Parser();
  ~Parser();
  int main(utils::DArray<DString> *);
  int main(int argc, char **args);
  int main();
  bool processCommand(Command * command, DString * whatIsSaid);
 private:
  CommandWords * commands;  // holds all valid command words
  Help * helper;
  utils::DArray<DString> * ans;
  double * parameters;
  int numberOfParameters;
  DString * autoSave;
  DString * currentDirectory;///<The current working directory
  DHashMap<Plane<Double> > * storedStuff;
  int commandCount;
			
  Command * getCommand();
  bool processCommand(Command * command);///<Returns True if the program should quit
  void printHelp();
  void start();
  bool openFile(DString fileName);
  DString getHistoryValue(Command * command);
  DString getSwapValue(Command * command);

};
            
#ifndef MASS_TO_SHEAR_CPP
static void verbosePrint(DString bean){std::cout << bean << std::endl;}
static void verbosePrint(int bean){std::cout << bean << std::endl;}
#endif

#endif

