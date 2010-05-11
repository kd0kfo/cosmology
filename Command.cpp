#include "Command.h"


	/**
	* Creates a new Command with the specified words. 
	*
	* @param words String command supplied.
	*/
	Command::Command(const DString& _words)
	{
	  using utils::DArray;
	  using utils::StringTokenizer;

	  args = new DArray<DString>;
	  commandWord = new DString("");
	  secondWord = new DString("");
	  this->words = new DString;
	  *(this->words) = _words;
	  
	  DString test("");
	  bool emptyCommand = words->equals(test);
	  if(!emptyCommand)
	    {
	      DString test2(*words);
	      StringTokenizer * toke = new StringTokenizer(test2);
	      while(toke->hasMoreTokens())
		{
		  args->put(toke->nextToken());
			}
	      *commandWord = (DString) args->get(0);
	      if(args->size() >= 2)
		*secondWord = (DString) args->get(1);
			delete toke;
			toke = 0;
	    }
	  else
	    {
	      *words += "help me";
	      *commandWord += "help";
	    }
	  
	  
	}
  
  DString Command::getCommandWord()
  {
    return *commandWord;
  }

    /**
     * Returns the second word of this command.
     */
	DString Command::getSecondWord()
    {
		return (secondWord == 0) ? DString("") : *secondWord;
    }

    /**
     * Return true if this command was not a known command.
     */
	bool Command::isUnknown()
    {
        return !(new CommandWords())->isCommand(*commandWord);
    }

    /**
     * Determines whether or not a command has a second word, which is usually 
     * the expression or parameter of the command
     */
	bool Command::hasSecondWord()
    {
        return (secondWord != 0);
    }
    
    /**
     * Returns an ArrayList of the words of the supplied command (including the command itself).
     */
utils::DArray<DString>& Command::getWords() const
    {
    	return *args;
    }
    
/**
     * Prints a variable in the LinAl object
     *
     * @param command Command Command to be Processed
     * @param la LinAl Stores variables
     */
DString Command::print(const Command& command, double * parameters, int numberOfParameters, DHashMap<Plane<Double> > * storedStuff, DString * currentDirectory)
  {
    utils::DArray<DString> args = command.getWords();
    DString toDo = args.get(0);
    DString returnMe;
    if(!Functions::isFunction(toDo))
      throw DavidException(toDo + " is not a valid function.","Invalid Function",DavidException::INVALID_ARGUMENT_ERROR_CODE);

    returnMe = Functions::doFunction(toDo,args.getArray()+1,args.size()-1,parameters,numberOfParameters, storedStuff, currentDirectory);    
    

    return returnMe;
  }//end "print"
  
    /**
     * Loads a file as a variable in the LinAl object
     *
     * @param command Command Command to be Processed
     * @param la LinAl Stores variables
     */
  utils::DArray<DString> * Command::load(DString fileName)
    {

      utils::DArray<DString> * returnMe = new utils::DArray<DString>;

      using namespace std;
      ifstream file (fileName.toCharArray());
      int buffersize = 200;
      char str[buffersize];
      if (file.is_open())
	{
	  while(file.getline(str,buffersize))
	    {
	      DString stickMeIn = str;
	      if(stickMeIn.contains('#'))
		stickMeIn = stickMeIn.substring(0,stickMeIn.indexOf('#'));
	      if(stickMeIn != "")
		returnMe->put(stickMeIn);
	    }
	}
      else
	throw DavidException(DString("Could not open: ") + fileName);


      return returnMe;
      
		
    }/**/

  void Command::save(const Command& command, DString& fileName, utils::DArray<DString> * output)
  {
    using namespace std;
    // verify that the file length is correct (it wasn't under Win95)
    ofstream myfile(fileName.toCharArray(), ios::out);
    DEBUG_PRINT(fileName);
    
    if(myfile.is_open())
      {
	
	for(int i = 0;i<output->size();i++)
	  myfile << output->get(i) << endl;
	myfile.close();
 }
    else
      {
	DEBUG_PRINT(DString("Could not open: ") + fileName);
      }
  }
