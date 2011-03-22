#include "Command.h"

	/**
	* Creates a new Command with the specified words. 
	*
	* @param words String command supplied.
	*/
Command::Command(const arg_t& _words)
	{
	  using namespace grmr;
	  static size_t num_keys = num_keywords();
	  this->words = _words;

	  if(words.size() != 0)
	    {
	      utils::StringTokenizer toke(_words);
	      	while(toke.hasMoreTokens())
			{
			  token_t newToken;
			  newToken.precedence = 0;
			  newToken._str = toke.nextToken();
			  if(newToken._str.size() == 0)
			    continue;
			  newToken.type = token_t::VARIABLE;
			  for(size_t i = 0;i<num_keys;i++)
			    if(newToken._str == keywords[i])
			      newToken.type = token_t::KEYWORD;
			  args.push_back(newToken);
			}
                commandWord = *args.begin();
		if(args.size() >= 2)
		  secondWord = args.at(1);
	    }
	  else
	    {
	      words += "help me";
	      commandWord.type = grmr::token_t::KEYWORD;commandWord._str = "help";
	    }
	}
  
    /**
     * Return true if this command was not a known command.
     */
    bool Command::isUnknown()const
    {
        CommandWords blah;
      return !(blah.isCommand(commandWord._str));
    }

    /**
     * Determines whether or not a command has a second word, which is usually 
     * the expression or parameter of the command
     */
	bool Command::hasSecondWord()const
    {
        return (secondWord._str.size() > 0);
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

  
arg_t Command::print()
  {
    return print_function(*this);
  }

