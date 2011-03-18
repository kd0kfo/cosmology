#ifndef BUILD_CPP
#define BUILD_CPP

#include <string>

class Build
{
 public:
  static std::string getBuild(){return "9/28/2010 20:00";}
  static std::string getVersion(){return getMajorVersion() + "." + getMinorVersion() +"."+getRevision();}
  static std::string getMajorVersion(){return "0";}
  static std::string getMinorVersion(){return "1";}
  static std::string getRevision(){return "2";}
  static std::string todo;
};

#endif
