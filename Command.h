#ifndef COMMAND_CPP
#define COMMAND_CPP

#include <iostream>
#include <fstream>
#include <vector>

#include "libdnstd/StringTokenizer.h"
#include "CommandWords.h"
#include "Help.h"

#ifndef DEBUG_PRINT
#ifdef __DEBUG__
#define DEBUG_PRINT(x) std::cout << x << std::endl;
#else
#define DEBUG_PRINT(x) ;
#endif
#endif

class DavidException;

	class Command
	{
		public:
			Command(const arg_t&);

			static arg_t getVersion(){return "2.0";}

			const arg_t& getCommandWord()const{return commandWord;}
			const arg_t& getSecondWord()const{return secondWord;}
			bool isUnknown()const;
			bool hasSecondWord()const;
			const args_t& getWords() const{return args;}
			const arg_t& getWholeCommandString() const{return words;}

			//Calculation methods
			static arg_t print(const Command& command);
                        static arg_t set(const Command& command);
                        static args_t load(const arg_t& fileName);

		private:
                    static arg_t print(const Command& command,std::string (*)(const Command&));
			arg_t words;
			arg_t commandWord;
			arg_t secondWord;
                        args_t args;
	};

        extern std::string print_function(const Command& arg);

#endif
