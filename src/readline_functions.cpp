#include <readline/readline.h>
#include <readline/history.h>
#include <vector>

#include "Help.h"
extern Help readline_help_function();
char* command_generator PARAMS((const char*, int));
char** parser_completion PARAMS((const char*, int, int));

std::vector<std::string> readline_matches;
int readline_current_match;

void init_readline()
{
  rl_readline_name = "Parser";
  readline_current_match = 0;
  rl_attempted_completion_function = parser_completion;
  rl_bind_key('\t', rl_possible_completions);
}


char** parser_completion(const char* text, int start, int end)
{
  char** matches;
  
  matches = (char**)NULL;

  if(start == 0)
    {
      readline_current_match = 0;
      matches = rl_completion_matches(text,command_generator);
    }
  
  return matches;
  }/**/

char* command_generator(const char* text, int state)
{
  Help commands = readline_help_function();
  Help::const_iterator it = commands.begin();
  int length,match_counter;
  
  match_counter = 0;
  if(readline_current_match >= commands.size())
    return (char*)NULL;

  length = strlen(text);
  for(;it != commands.end();it++)
    {
      if(!strncmp(it->first.c_str(),text,length) && match_counter >= readline_current_match)
	{
	  char* returnMe = (char*) calloc(it->first.size()+1,sizeof(char));
	  strcat(returnMe,it->first.c_str());
	  readline_current_match++;
	  return returnMe;
	}
      match_counter++;
    }

  return (char*) NULL;
  }
