#include "utilities.h"
#include "DavidException.h"

int main(int argc, char** argv)
{
#ifdef __DEBUG__
  try
#endif

  {
    Utilities u;
    return u.sub_main(argc,argv);
  }
#ifdef __DEBUG__
  catch(DavidException de)
    {
      de.stdOut();
    }
#endif

}
