/**
 * Provides functions to compile subroutines.
 */

#include "Help.h"

#define PROGRAM_NAME "mycosmo"
#define GUI_CPP 1
#define GENERAL_ERROR 1
#define IO_ERROR 2
#define SYNTAX_ERROR 3

const int mycosmoVersion[] = {0,0,0};

typedef struct {
  std::string* _string;
  int* _int;
  float* _float;
}Object;


Help mycosmoHelp()
{
    using namespace std;
    Help returnMe;
    returnMe["nothing"] = "to do";
    return returnMe;
}
