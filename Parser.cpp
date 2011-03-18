#include "Parser.h"

#ifdef USE_READLINE
#include "readline_functions.cpp"
#endif

Parser::Parser() {
#ifdef USE_READLINE
  init_readline();
#endif
    commands = new CommandWords();
    h = new Help();
    inputBuffer = new arg_t("");
    history = new args_t();
    program_name = "";
}

Parser::~Parser() {
    delete commands;
    delete h;
    delete history;
    h = 0;
    history = 0;
}

void Parser::main() {
    this->start();
}

void Parser::main(int argc, char **args) {
    args_t arrargs;
    this->program_name = args[0];
    for (int i = 1; i < argc; i++) {
        arrargs.push_back(args[i]);
    }
    main(arrargs);
}

void Parser::main(args_t& args) {
    if (args.size() >= 1 && args.at(0) == "--version") {
        std::cout << "Version " + Build::getVersion() + " built on " + Build::getBuild() << std::endl;
        return;
    }
    if (args.size() >= 1 && args.at(0) == "--build") {
        return;
    }

    if (args.size() >= 1 && (args.at(0) == "--help" || args.at(0) == "--usage")) {
        std::cout << this->program_name <<  Build::getVersion() << std::endl;
        std::cout << "Last built on: " << Build::getBuild() << std::endl;
	args_usage();
    } else if (args.size() >= 1 && args.at(0) == "--console") {
        start();
    } else {
        start();
    }

}

Command * Parser::getCommand() {
    using namespace std;

    try {
#ifndef USE_READLINE
    char inputLine[256]; // will hold the full input line
        cout << "> "; // print prompt
        cin.get(inputLine, 256);
        cin.ignore(255, '\n');
	cout.flush();
	return new Command(inputLine);
#else
	char *inputLine = readline("> ");
        add_history(inputLine);
	return new Command(inputLine);
#endif
    } catch (...) {
        throw DavidException("Input error");
    }
    

}

/**
 * Print out a list of valid command words.
 */
void Parser::printHelp(arg_t& whatIsSaid) {
    if(whatIsSaid.size() != 0)
      {
	whatIsSaid = h->getHelp(whatIsSaid);
	return;
      }
    whatIsSaid = "For help, enter help <command>\nCommands:\n";
    for(Help::const_iterator it = h->begin();it != h->end();it++)
      whatIsSaid += it->first + " ";
}

void Parser::start() {
  std::cout << "Welcome to " << this->program_name << ", version:" << getVersion() << "\nFor help with the program, enter help" << std::endl;
    // Enter the main command loop.  Here we repeatedly read commands and
    // execute them until the game is over.

    bool finished = false;
    while (!finished) {
        try {
            Command * command = getCommand();
            finished = processCommand(*command);
            delete command;
            command = 0;
        } catch (DavidException e) {
            std::cout << "Process Error: " << e.what() << std::endl;
        }
    }
    std::cout << "Thank you.  Good bye." << std::endl;
}

bool Parser::processCommand(Command& command) {
    std::string dummyString = "";
    bool returnMe = processCommand(command, dummyString);
    return returnMe;
}

bool Parser::processCommand(Command& command, std::string& whatIsSaid) {
    bool wantToQuit = false;
    
    if (command.getCommandWord() != "save"  && command.getCommandWord() != "exit" && command.getCommandWord() != "quit" && command.getCommandWord() != "history")
        history->push_back(command.getWholeCommandString());

    std::string commandWord = command.getCommandWord();
    if (commandWord== "help"){
      whatIsSaid = command.getSecondWord();
        printHelp(whatIsSaid);
        std::cout << whatIsSaid << std::endl;
    } else if (commandWord == "version") {
        std::cout << "Version: " << Build::getVersion() << std::endl;
    } else if (commandWord == "history") {
        if (history->size() == 0) {
            std::cout << "No commands have been entered yet." << std::endl;
            return false;
        }
        std::cout << "Previous Commands:" << std::endl;
        int counter = 0;
        for (args_t::iterator it = history->begin();
                it != history->end(); it++)
                {
                std::cout << counter++ << ": " << *it << std::endl;
                }
    }else if (commandWord == "quit" || commandWord == "exit") {
        wantToQuit = true;
    }
    else
    {
      whatIsSaid = Command::print(command.getWholeCommandString());
      std::cout << whatIsSaid << std::endl;
    }

    return wantToQuit; /**/
}

arg_t Parser::getHistoryValue(Command& command) const {

    size_t historyCount = history->size() - 1;
    std::istringstream int_buff(command.getCommandWord().substr(1));
    
    if (command.getCommandWord().substr(0, 2) != "!!")
    {
	    int_buff >> historyCount;
    	if(int_buff.fail() || historyCount < 0)
        	throw DavidException("I need a positive number for the place in history");
    }

    if (historyCount >= history->size())
        throw DavidException("Ummm, that is in the future.");

    std::string newCommand = history->at(historyCount);
    args_t words = command.getWords();
    for (size_t i = 1; i < words.size(); i++)
        newCommand += " " + words.at(i);
    return newCommand;

}

void Parser::setHelp(const Help& newHelp)
{
  delete h;
  h = new Help(newHelp);
}

