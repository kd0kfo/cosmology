#include "CommandWords.h"

	CommandWords::CommandWords()
	{
	  validCommands = new utils::DArray<DString>;
	  validCommands->put("about");
	  validCommands->put("copy");
	  validCommands->put("rename");
	  validCommands->put("history");
	  validCommands->put("last");
	  validCommands->put("list");
	  validCommands->put("load");
	  validCommands->put("print");
	  validCommands->put("remove");
	  validCommands->put("save");
	  validCommands->put("help");
	  validCommands->put("exit");
	  validCommands->put("quit");
	  validCommands->put("run");
	  validCommands->put("autosave");
	}

  bool CommandWords::isCommand(DString& aString)
    {
        for(ushort i = 0; i < validCommands->size(); i++) 
            if(validCommands->get(i) == aString)
                return true;
		return false;
    }

  DString CommandWords::showAll() 
    {
    	DString bean = "";
	
	DString * junks = validCommands->getArray();
	
	DString::alphabetize(junks,validCommands->size());
	
        for(int i = 0; i < validCommands->size(); i++) {
            bean += DString(" ") + junks[i];
        }

	delete [] junks;
	junks = 0;
        return bean + " and functions";
    }

