#include "Command.h"

	/**
	* Creates a new Command with the specified words. 
	*
	* @param words String command supplied.
	*/
	Command::Command(const arg_t& _words)
	{
	  commandWord = "";
	  secondWord = "";
	  this->words = _words;
	  
	  if(words.size() != 0)
	    {
	      utils::StringTokenizer toke(_words);
	      	while(toke.hasMoreTokens())
			{
			  std::string currString = toke.nextToken();
			  args.push_back(currString);
			}
                commandWord = *args.begin();
			if(args.size() >= 2)
                            secondWord = args.at(1);
	    }
	  else
	    {
	      words += "help me";
	      commandWord += "help";
	    }
	}
  
    /**
     * Return true if this command was not a known command.
     */
    bool Command::isUnknown()const
    {
        CommandWords blah;
      return !(blah.isCommand(commandWord));
    }

    /**
     * Determines whether or not a command has a second word, which is usually 
     * the expression or parameter of the command
     */
	bool Command::hasSecondWord()const
    {
        return (secondWord.size() > 0);
    }
    
    
    /**
     * Loads a file as a variable in the LinAl object
     *
     * @param command Command Command to be Processed
     * @param la LinAl Stores variables
     */
  args_t Command::load(const arg_t& fileName)
    {

      using namespace std;
      args_t returnMe;

      ifstream file (fileName.c_str());
      string line;
      if (file.is_open())
	{
          line = "";
	  while(!file.eof())
	    {
	      file >> line;
              returnMe.push_back(line);
	    }
	}
      else
	throw DavidException("Could not open: " + fileName);

      return returnMe;
      
		
    }/**/

  
  arg_t Command::print(const Command& command)
  {
      return Command::print(command,print_function);
  }

  arg_t Command::print(const Command& command,std::string the_function(const Command&))
  {
      return the_function(command.getWholeCommandString());
  }

