#ifndef COMMANDWORDS_CPP
#define COMMANDWORDS_CPP

#include "libdnstd/DArray.h"
#include "libdnstd/DString.h"


class CommandWords
{
 public:
  CommandWords();
  bool isCommand(DString&);
  DString showAll();
  //DArray getCommandWords();
  
 private: 
  utils::DArray<DString> * validCommands;
};


#endif
