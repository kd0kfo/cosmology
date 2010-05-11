#ifndef COMMAND_CPP
#define COMMAND_CPP

#include "DString.h"
#include "DArray.h"
#include "StringTokenizer.h"
#include "CommandWords.h"
#include "Functions.h"

	class Command
	{
		public:
			Command(const DString&);

			static DString getVersion(){return (DString) "1.1";}

			DString getCommandWord();
			DString getSecondWord();
			bool isUnknown();
			bool hasSecondWord();
			utils::DArray<DString>& getWords() const ;
			DString getWholeCommandString() const{return *words;}
			static utils::DArray<DString> * Command::load(DString fileName);
			static void save(const Command& command, DString& fileName, utils::DArray<DString> * output);

			//Calculation methods
			/*
			 * Does a calculation and prints the result
			 */
			static DString print(const Command& command,double * parameters, int numberOfParameters, DHashMap<Plane<Double> > * storedStuff,DString * currentDirectory);
            
		private:
			DString * words;
			DString * commandWord;
			DString * secondWord;
			utils::DArray<DString> * args;
	};

#endif 

