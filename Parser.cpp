#include "Parser.h"

#ifndef __USE_MAIN__
int main(int argc,char** argv)
{

  Parser * P = new Parser();
  int returnMe = 0;
  try{
    returnMe = P->main(argc,argv);
  }
  catch(DavidException de)
    {
      #ifdef __DEBUG__
            DEBUG_PRINT(DString("Type: ") + de.getType());
	    DEBUG_PRINT(DString("Code: ") + Double(de.getCode()).toDString());
	    DEBUG_PRINT(DString("Message: "));
      #endif
	    
	    returnMe = de.getCode();
    }
  delete P;
  return returnMe;
  
}
#endif __USE_MAIN__

Parser::Parser() 
{
#ifdef __DEBUG__
  int _DEBUG_COUNTER = 0;
#endif
  commands = new CommandWords();
  helper = new Help();
  ans = new utils::DArray<DString>;
  currentDirectory = new DString("./");
  autoSave = 0;
  commandCount = 0;

  parameters = new double[2];
  parameters[0] = 0.7;//Hubble Parameter, H_0 = h*100 km/s/Mpc WMAP 3-yr values
  parameters[1] = 0.26;//Mass Density Parameter, WMAP 3-yr values

  storedStuff = new DHashMap<Plane<Double> >;
}

Parser::~Parser() 
{
  delete commands;
  delete helper;
  delete autoSave;
  delete currentDirectory;
  delete ans;
  delete parameters;
  delete storedStuff;
  helper = 0;
  currentDirectory = autoSave = 0;
  ans = 0;
  parameters = 0;
  storedStuff = 0;
}

int Parser::main()
{ 
  try
    {
      this->start();
      return 1;
    }
  catch(DavidException de)
    {
      #ifdef __DEBUG__
            DEBUG_PRINT(DString("Type: ") + de.getType());
	    DEBUG_PRINT(DString("Code: ") + Double(de.getCode()).toDString());
	    DEBUG_PRINT(DString("Message: "));
      #endif
	    
	    return de.getCode();
    }
  
}

int Parser::main(int argc,char **args)
{
  using utils::DArray;
	
  int returnMe = 0;
  
  DArray<DString> * arrargs = new DArray<DString>();
  for(int i=1;i<argc;i++)
    {
      arrargs->put(DString(args[i]));
    }
  returnMe = main(arrargs);
  delete arrargs;
  
  return returnMe;
}

int Parser::main(utils::DArray<DString> * args)
{
  int returnMe = 0;
  if(args->length() >= 1 && (args->get(0) == "--help" || (args->get(0)).equals("--usage")))
    {
      David::say("Cosmo Environment");
      David::say("Usage: ./mycosmo.exe [options]");
      David::say("Options:");
      David::say("--console, Start the program within a console.");
      David::say("--help, this");
      return 0;
    }
  else if(args->length() >=1 && ((args->get(0)).equals("--gui") || args->get(0).equals("-x")) )
    {
    }
  else if(args->length() >= 1 && (args->get(0).equals("--console")))
    {
      try
	{
	  start();
	}
      catch(DavidException de)
	{
#ifdef __DEBUG__
            DEBUG_PRINT(DString("Type: ") + de.getType());
	    DEBUG_PRINT(DString("Code: ") + Double(de.getCode()).toDString());
	    DEBUG_PRINT(DString("Message: "));
#endif
	    return de.getCode();
	}
    }
  else if(args->length() >= 2 && (args->get(0).equals("--load")))
    {
      DString fileName = args->get(1);
      DString variableName = args->get(2);
      DavidException * gde = 0;
      try
	{
	  Plane<math::Complex> * cPlane = Plane<math::Complex>::readPlane(fileName);
	  Plane<Double> * newGuy = Functions::convertComplexPlaneToDoublePlane(cPlane);
	      
	  storedStuff->put(variableName,*newGuy);
	      
	  ans->put(DString("load ") + fileName + DString(" ") + variableName);

	  delete newGuy;
	  delete cPlane;
	  cPlane = 0;
	  newGuy = 0;
	}
      catch(DavidException de)
	{
	  gde = new DavidException(de);
	}
	  
      if(gde != 0)
	{
	  gde->stdErr();
	  returnMe = gde->getCode();
	  delete gde;
	  gde = 0;
	  return returnMe;
	}
      else
	start();
	  
    }
  else if(args->length() >= 2 && (args->get(0).equals("--run")))
    {
      DString fileName = args->get(1);
      utils::DArray<DString> * toBeRun;
      DavidException * gDe = 0;
      int lineNumber = 0;
      bool shouldQuit = false;
      try{
	toBeRun = Command::load(fileName);
	Command tmpCommand("blah");
	for(int i = 0;i<toBeRun->size();i++)
	  {
	    lineNumber++;
	    tmpCommand = Command(toBeRun->get(i));
	    if(processCommand(&tmpCommand))
	      shouldQuit = true;
	  }
      }
      catch(DavidException de)
	{
	  gDe = new DavidException(de);
	}

      delete ans;
      ans = toBeRun;
      if(gDe != 0)
	{
	  std::cout << "Error on line " << lineNumber << " of " << fileName << ": " << gDe->getMessage() << std::cout;
	  returnMe = gDe->getCode();
	  delete gDe;
	  shouldQuit = true;
	}
      if(!shouldQuit)
	start();
    }
  else 
    {
      start();//GUI.main(new String[]{});
    }
    
  return returnMe;
} 

Command * Parser::getCommand() 
{
  using namespace std;
  char inputLine[256];   // will hold the full input line
  cout << "> ";     // print prompt
		
  try{
    cin.get(inputLine,256);
  }
  catch(...)
    {
      throw DavidException("Input error");
    }
  cout.flush();
  cin.ignore(255,'\n');
						
  if(strlen(inputLine) == 0 || DString(inputLine) == "")
    {
      return new Command(DString("help"));
    }
  else
    {
      return new Command(DString(inputLine));
    }
  /*else
    return new Command(null);*/
}

/**
 * Print out a list of valid command words.
 */
void Parser::printHelp()
{
  David::say(commands->showAll().toCharArray());
}

void Parser::start()
{
  David::say("Welcome!\nFor help with the program, enter help");
  // Enter the main command loop.  Here we repeatedly read commands and
  // execute them until the game is over.

  bool finished = false;
  while (! finished) {
    try{
      Command * command = getCommand();
      if(command->getCommandWord().length() >= 2)
	if(command->getCommandWord().charAt(0) == '!')
	  {
	    DString newCommand = getHistoryValue(command);
	    David::say(newCommand);
	    *command = Command(newCommand);
	  }
	else if(command->getCommandWord().charAt(0) == '^')
	  {
	    DString newCommand = getSwapValue(command);
	    David::say(newCommand);
	    *command = Command(newCommand);
	  }
	    
      utils::StringTokenizer subcommands(command->getWholeCommandString(),";");

      while(subcommands.hasMoreTokens())
	{
	  *command = Command(subcommands.nextToken());
	  finished = processCommand(command);
	}
      delete command;
      command = 0;
    }
    catch(DavidException e)
      {
	David::say((DString) "Process Error: "+e.getMessage());
      }/**/
  }
  std::cout << "Thank you.  Good bye." << std::endl;
}


bool Parser::processCommand(Command * command)
{
  DString * dummyString = new DString;
  bool returnMe = processCommand(command,dummyString);
  delete dummyString;
  dummyString = 0;
  return returnMe;
}

bool Parser::processCommand(Command * command, DString * whatIsSaid)
{
	
  commandCount++;

  if(commandCount % 10 == 0 && autoSave != 0)
    Command::save(*command,*autoSave,ans);


  bool wantToQuit = false;
		
  DString ansBuffer = command->getWords().get(0);
  for(int i = 1;i<command->getWords().size();i++)
    ansBuffer += DString(" ") + command->getWords().get(i);

  if(command->getCommandWord() != "save" && command->getCommandWord() != "run" && command->getCommandWord() != "autosave")
    ans->put(ansBuffer);
	  


  if(command->isUnknown()) 
    {
      //LTree tmp = AlPars::infixLTree(command->getCommandWord(),*la);
      DString printString = Command::print(*command, parameters,numberOfParameters,storedStuff, currentDirectory);
      David::say(printString);
      *whatIsSaid = printString;
			
      return false;
    }

	

  DString commandWord = command->getCommandWord();
  if (commandWord.equals("help") && command->getSecondWord() == "")
    {
      printHelp();
    }
  else if(commandWord.equals("help") && command->getSecondWord() != "")
    {
				
      if(command->getSecondWord() == "functions")
	{
	  DString extraPart = "list";
	  if(command->getWords().size() >= 3)
	    extraPart = command->getWords().get(2);
	  *whatIsSaid = Functions::getHelp(extraPart);
	  David::say(*whatIsSaid);
	}
      else if(Functions::isFunction(command->getSecondWord()))
	{
	  *whatIsSaid = Functions::getHelp(command->getSecondWord());
	  David::say(*whatIsSaid);
	}
      else
	{
	  *whatIsSaid = helper->getHelp(command->getSecondWord());
	  David::say(*whatIsSaid);
	}
    }
  else if(commandWord.equals("autosave"))
    {
      if(command->getWords().size() == 1)
	{
	  if(autoSave == 0)
	    David::say("autosave is current not set");
	  else
	    David::say(DString("autosaving to ")+*autoSave);
	  return false;
	}

      if(autoSave != 0)
	delete autoSave;
	    
      if(command->getWords().get(1) == "off")
	{
	  David::say("autosave is off.");
	  autoSave = 0;
	  return false;
	}

      autoSave = new DString(command->getWords().get(1));
      David::say(DString("autosaving to ") + *autoSave);
	    
      return false;

    }
  else if(commandWord.equals("about"))
    {
	    
      if(command->getWords().size() == 1 || command->getSecondWord().equals("all"))
	{
	  if(storedStuff->getNumberOfKeys() == 0)
	    {
	      David::say("There are no stored variables.");
	      return false;
	    }
		
	  DIterator<Plane<Double> > * dude = new DIterator<Plane<Double> >(storedStuff);
	  DString * currentKey = new DString();
	  while(dude->hasNext())
	    {
	      const Plane<Double>& junk = dude->next(currentKey);
	      *whatIsSaid = *currentKey + ": ";
	      *whatIsSaid += Double(junk.numberOfRows()).toDString() + DString(" by ") + Double(junk.numberOfColumns()).toDString() + " array.";
	      David::say(*whatIsSaid);
	    }
	  return false;
	}
	    
      if(!storedStuff->isKey(command->getSecondWord()))
	{
	  David::say(command->getSecondWord() + " was not found.");
	  return false;
	}
      Plane<Double> * junk = new Plane<Double>(storedStuff->get(command->getSecondWord()));
      *whatIsSaid = command->getSecondWord() + ": ";
      *whatIsSaid += Double(junk->numberOfRows()).toDString() + DString(" by ") + Double(junk->numberOfColumns()).toDString() + " array.";
      David::say(*whatIsSaid);

      delete junk;
    }
  else if(commandWord.equals("last"))
    {
      David::say(ans->get(ans->size()-1));
    }
  else if(commandWord.equals("print"))
    {
      if(command->getWords().size() < 4)
	{
	  David::say("I need the name of the plane and the row and column numbers in order to print the value.");
	  return false;
	}
      DString * words = command->getWords().getArray();

      if(!storedStuff->isKey(words[1]))
	{
	  David::say(words[1] + " was not found.");
	  return false;
	}

      int i = (int) Double(words[2]).doubleValue();
      int j = (int) Double(words[3]).doubleValue();

      DString stuff;
      try
	{
	  Plane<Double> poop(storedStuff->get(words[1]));
	  Double junk = poop.getValue(i,j);
	  stuff = junk.toDString();
	  if(junk.getValue(1) != 0 || junk.getValue(2) != 0)
	    {
	      stuff = DString("(") + stuff + DString(", ")+Double(junk.getValue(1)).toDString();
	      if(junk.getValue(2) != 0)
		stuff += DString(", ")+Double(junk.getValue(2));
	      stuff += ")";
	    }
		
	}
      catch(DavidException de)
	{
	  de.stdOut();
	  return false;
	}
	    
      *whatIsSaid = words[1] + DString("@(") + words[2] + DString(",") + words[3] + DString(") equals ") + stuff;
      David::say(*whatIsSaid);
    }
  else if(commandWord.equals("remove"))
    {
	 
      if(command->getSecondWord().equals("all"))
	{
	  storedStuff->removeAll();
	  David::say("All variables have been removed.");
	  return false;
	}

      utils::DArray<DString> deleteUs = command->getWords();

      for(int i = 1; i< deleteUs.size();i++)
	{

	  DString curr = deleteUs.get(i);
	  storedStuff->removeKey(curr);
	  David::say(curr + " has been removed.");
	}
    }
  else if(commandWord.equals("list"))
    {
      std::vector<DString> keys = storedStuff->getKeys();
      int numberOfKeys = keys.size();

      if(numberOfKeys == 0)
	{
	  David::say("There are no stored planes");
	  return false;
	}

      for(int i = 0;i<numberOfKeys;i++)
	*whatIsSaid += keys[i] + " ";
      David::say(*whatIsSaid);
    }
  else if(commandWord.equals("copy"))
    {

      if(command->getWords().size() < 3)
	{
	  David::say("Name of original and name new plane are both required.");
	  return false;
	}

      DString old = command->getWords().get(1);
      DString _new = command->getWords().get(2);

      if(!storedStuff->isKey(old))
	{
	  David::say(old + " was not found.");
	  return false;
	}
	  
      if(storedStuff->isKey(_new))
	storedStuff->removeKey(_new);

      storedStuff->put(_new,storedStuff->get(old));
      David::say(old + DString(" has been copied to ") + _new);
    }
  else if(commandWord.equals("rename"))
    {

      if(command->getWords().size() < 3)
	{
	  David::say("New and old name of the plane are required.");
	  return false;
	}

      DString old = command->getWords().get(1);
      DString _new = command->getWords().get(2);

      if(!storedStuff->isKey(old))
	{
	  David::say(old + " was not found.");
	  return false;
	}
	  
      if(storedStuff->isKey(_new))
	storedStuff->removeKey(_new);

      storedStuff->put(_new,storedStuff->get(old));
      storedStuff->removeKey(old);
      David::say(old + DString(" has been renamed ") + _new);
    }
  else if(commandWord.equals("history"))
    {
      for(int i = 0;i<ans->size();i++)
	David::say(Double((double) i).toDString() + DString(" ")+ ans->get(i));

      return false;
    }
  else if(commandWord.equals("run") || commandWord.equals("open"))
    {
      if(!command->hasSecondWord())
	{
	  David::say("Uh, what do you want me to open?");
	  return false;
	}
      utils::DArray<DString> * toRun = 0;
      Command * tmpguy = new Command("blah");
      try
	{
	  DString loadMe = command->getSecondWord();
	  if(loadMe.charAt(0) != '/')
	    {
	      if(currentDirectory == 0)
		currentDirectory = new DString("./");
	      
	      loadMe = *currentDirectory + loadMe;
	    }
	    
	  toRun = Command::load(loadMe);
	  for(int i = 0;i<toRun->size();i++)
	    {
	      *tmpguy = Command(toRun->get(i));
	      processCommand(tmpguy);
	    }
	}
      catch(DavidException de)
	{
	  David::say(de.getMessage());
	}
      delete tmpguy;
      delete toRun;
      return false;
    }

  else if(commandWord.equals("load"))
    {
	  
      if(command->getWords().size() < 3)
	{
	  David::say("Name of file and variable name are both required.");
	  return false;
	}

      if(command->getWords().get(2).equals("all"))
	{
	  David::say("Cannot store a plane as \"all\".");
	  return false;
	}

      bool isComplexPlane = true;
      for(int i = 0;i<command->getWords().size();i++)
	{
	  if(command->getWords().get(i).equals("--real"))
	    isComplexPlane = false;
	}

      DString variable = command->getWords().get(2);
      if(storedStuff->isKey(variable))
	{
	  David::say(variable + " already existed, but now it's being replaced.");
	  storedStuff->removeKey(variable);
	}
	  
      if(currentDirectory == 0)
	currentDirectory = new DString("./");
	  
      David::say(DString("Loading ") + *currentDirectory + command->getSecondWord());
      Plane<Double> * newGuy;
      if(!isComplexPlane)
	newGuy = Plane<Double>::readPlane(*currentDirectory + command->getSecondWord());
      else
	{
	  Plane<math::Complex> * tmpCPlane = Plane<math::Complex>::readPlane(*currentDirectory + command->getSecondWord());
	  newGuy = Functions::convertComplexPlaneToDoublePlane(tmpCPlane);
	  delete tmpCPlane;
	  tmpCPlane = 0;
	}
	  
      storedStuff->put(command->getWords().get(2),*newGuy);
      delete newGuy;
      newGuy = 0;
      *whatIsSaid = command->getWords().get(2) + " has been loaded.";
      David::say(*whatIsSaid);
    }//end load part
  else if(commandWord.equals("save"))
    {

      if(currentDirectory == 0)
	currentDirectory = new DString("./");

      if(command->getWords().size() == 1)
	{
	  if(autoSave == 0)
	    {
	      David::say("I need to know the filename to which I am saving.");
	    }
	  else
	    {
	      DString saveMe = *currentDirectory + *autoSave;
	      Command::save(*command,saveMe,ans);
	    }
	  return false;
	}

      if(command->getWords().size() == 2)
	{
	  DString fileName = command->getSecondWord();
	  DString saveMe = *currentDirectory + fileName;
	  Command::save(*command,saveMe,ans);
	  return false;
	}

      if(!storedStuff->isKey(command->getSecondWord()))
	{
	  David::say(DString("Sorry, ") + command->getSecondWord() + " is not a stored plane.");
	  return false;
	}

      bool isComplexPlane = true;
      for(int i = 0;i<command->getWords().size();i++)
	{
	  if(command->getWords().get(i).equals("--real"))
	    isComplexPlane = false;
	}

      if(!isComplexPlane)
	storedStuff->get(command->getSecondWord()).savePlane(*currentDirectory + command->getWords().get(2));
      else
	{
	  Plane<Double> tmpDoublePlane = storedStuff->get(command->getSecondWord());
	  Plane<math::Complex> * printMe = Functions::convertDoublePlaneToComplexPlane(&tmpDoublePlane);
	  printMe->savePlane(command->getWords().get(2));
		
	  delete printMe;
	  printMe = 0;
	}

      David::say(command->getSecondWord() + DString(" saved as ") + command->getWords().get(2));
    }
  else if (commandWord.equals("quit") || commandWord.equals("exit")) 
    {
      wantToQuit = true;
      if(autoSave != 0)
	Command::save(*command,*autoSave,ans);

    }

  return wantToQuit;/**/

	   
}
	
DString Parser::getHistoryValue(Command * command)
{

  int historyCount = ans->size() - 1;
  if(!command->getCommandWord().substring(0,2).equals("!!"))
    historyCount = (int) Double(command->getCommandWord().substring(1)).doubleValue();
  
  if(historyCount < 0)
    throw DavidException("I need a positive number for the place in history");
  
  if(historyCount >= ans->size())
    throw DavidException("Ummm, that is in the future.");

  DString newCommand = ans->get(historyCount);
  utils::DArray<DString> words = command->getWords();
  for(int i = 1; i< words.size();i++)
    newCommand += DString(" ") + words.get(i);
  return newCommand;
	 
}

DString Parser::getSwapValue(Command * command)
{


  using utils::StringTokenizer;
  StringTokenizer tokie(command->getWholeCommandString(),"^");

  DString lhs,rhs;
  DavidException hopefullyNot("There need to be two ^ each with a string after them.");

  if(!tokie.hasMoreTokens())
    throw hopefullyNot;

  lhs = tokie.nextToken();

  if(!tokie.hasMoreTokens())
    throw hopefullyNot;

  rhs = tokie.nextToken();


  DString lastCommand;

  if(ans->size() >= 1)
    lastCommand = ans->get(ans->size() - 1);
  else
    throw DavidException("There are no previous commands");



  StringTokenizer newTokie(lastCommand);

  DString newCommand;

  while(newTokie.hasMoreTokens())
    {
      DString curr = newTokie.nextToken();
      if(curr == lhs)
	curr = rhs;

      newCommand += DString(" ") + curr;
    }


  return newCommand.trim();

}
