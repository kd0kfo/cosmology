#ifndef __RAINBOW_CPP__
#define __RAINBOW_CPP__

//---------------------------------------------------------------------------
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
//---------------------------------------------------------------------------
//From http://www.richelbilderbeek.nl/CppRainbow.htm

 class Rainbow
{
public:
  static void getRainbow(const double, unsigned char&, unsigned char&, unsigned char&);

  static const unsigned char GetRed(const double x);
  static const unsigned char GetGreen(const double x);
  static const unsigned char GetBlue(const double x);
  static char * yuvToRGB(double y, double u, double v);
  static char * hsvToRGB(double h, double s, double v);

};

#endif
