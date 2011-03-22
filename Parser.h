#ifndef PARSER_CPP
#define PARSER_CPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <sstream>

#include "Build.h"
#include "Command.h"
#include "CommandWords.h"
#include "Help.h"

#ifdef USE_READLINE
#ifdef _READLINE_H_
#  include "posixstat.h"
#  include "readline.h"
#  include "history.h"
#else
#  include <sys/stat.h>
#  include <readline/readline.h>
#  include <readline/history.h>
#endif
#endif /**USE_READLINE**/

extern void args_usage();
class DavidException;

  class Parser
  {
  public:
#ifdef __DEBUG__
    int DEBUG_COUNTER;
#endif

    Parser();
    ~Parser();
    static std::string getVersion(){return Build::getVersion();}
    /**
     * Starts the parser.
     */
    void start();

    void setHelp(const Help& newHelp);
    void main(args_t& arguments);
    void main(int argc, char **args);

    /**
     * Same as start()
     */
    void main();

    bool processCommand(Command& command, std::string& whatIsSaid);
    arg_t getHistoryValue(Command& command) const;

    std::string program_name;
    

  private:
    CommandWords * commands;  // holds all valid command words
    Help*  h;
    arg_t* inputBuffer;
    args_t* history;


    Command * getCommand();
    bool processCommand(Command& command);
    void printHelp(arg_t& whatIsSaid);
    bool openFile(const std::string& fileName);


  };

#endif
