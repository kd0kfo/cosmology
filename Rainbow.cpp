#include "Rainbow.h"


void Rainbow::getRainbow(
  const double x,
  unsigned char& r,
  unsigned char& g,
  unsigned char& b)
{
  /** ORiginal /**/
  const int r0 = GetRed(x);
  const int g0 = GetGreen(x);
  const int b0 = GetBlue(x);
  const int max = std::max(r0, std::max(g0,b0));
  assert(max!=0);

  r = 255.0 * static_cast<double>(r0) / static_cast<double>(max);
  g = 255.0 * static_cast<double>(g0) / static_cast<double>(max);
  b = 255.0 * static_cast<double>(b0) / static_cast<double>(max);
  /**/


  /** Mine /*
  char* rgb = hsvToRGB(x*350,1,1);
  
  r = rgb[0];
  g = rgb[1];
  b = rgb[2];
  
  delete [] rgb;
  rgb = 0;/**/
}
//---------------------------------------------------------------------------
//From http://www.richelbilderbeek.nl/CppRainbow.htm
const unsigned char Rainbow::GetRed(const double x)
{
  assert( x >= 0.0 && x < 1.0);
  /** Original Function /*
  const double f = std::max(0.0,
    (x < 0.5
    ?  std::cos(x * 1.5 * M_PI)
    : -std::sin(x * 1.5 * M_PI)
    ) );/**/

  /** My Change /*
  const double f = std::max(0.0,
			    (x < 0.5
			     ?  std::cos(x * M_PI)
			     : 0));/**/
  
  const double f = std::max(0.0,-2*x+1);
			
  assert( f >= 0.0);
  assert( f <= 1.0);
  const double y = 255.0 * f;
  assert( static_cast<int>(y) < 256 );
  return static_cast<unsigned char>(y);
}
//---------------------------------------------------------------------------
//From http://www.richelbilderbeek.nl/CppRainbow.htm
const unsigned char Rainbow::GetGreen(const double x)
{
  assert( x >= 0.0 && x < 1.0);

  /** Original /*
      const double f = std::max(0.0, std::sin( x * 1.5 * M_PI ) );/**/

  /** My changes/*
  const double f = std::max(0.0, std::sin( 2*x*M_PI-0.5*M_PI ));/**/

  const double f = (x < 0.5) ? std::max(0.0,3*x-0.5) : std::max(0.0,-3*x+2.5);

  assert( f >= 0.0);
  assert( f <= 1.0);
  const double y = 255.0 * f;
  assert( static_cast<int>(y) < 256 );
  return static_cast<unsigned char>(y);
}
//---------------------------------------------------------------------------
//From http://www.richelbilderbeek.nl/CppRainbow.htm
const unsigned char Rainbow::GetBlue(const double x)
{
  assert( x >= 0.0 && x < 1.0);

  /**Original /*
  const double f = std::max(0.0, -std::cos( x * 1.5 * M_PI ) );/**/

  /** MY changes /*
  const double f = std::max(0.0, -std::cos( x *  M_PI ) );/**/

  const double f = std::max(0.0,2*x-1);

  assert( f >= 0.0);
  assert( f <= 1.0);
  const double y = 255.0 * f;
  assert( static_cast<int>(y) < 256 );
  return static_cast<unsigned char>(y);
}
//---------------------------------------------------------------------------

char * Rainbow::yuvToRGB(double y, double u, double v)
{
  char * returnMe = new char[3];

  double r = y+1.13983*v;
  double g = y-0.39465*u-0.58060*v;
  double b = y+2.03211*u;
  
  assert( static_cast<int>(255.0*r) < 256 );
  returnMe[0] = static_cast<unsigned char>(255.0*r);

  assert( static_cast<int>(255.0*g) < 256 );
  returnMe[1] = static_cast<unsigned char>(255.0*g);

  assert( static_cast<int>(255.0*b) < 256 );
  returnMe[2] = static_cast<unsigned char>(255.0*b);

  return returnMe;
}


char * Rainbow::hsvToRGB(double h, double s, double v)
{
  int h_i = static_cast<int>(floor(h/60)) % 6;
  double f = h/60 - floor(h/60);

  double p,q,t;
  p = v*(1-s);
  q = v*(1-f*s);
  t = v*(1-s+f*s);

  double * rgb = new double[3];
  
  switch(h_i)
    {
    case 0:
      {
	rgb[0] = v;
	rgb[1] = t;
	rgb[2] = p;
	break;
      }
    case 1:
      {
	rgb[0] = q;
	rgb[1] = v;
	rgb[2] = p;
	break;
      }
    case 2:
      {
	rgb[0] = p;
	rgb[1] = v;
	rgb[2] = t;
	break;
      }
    case 3:
      {
	rgb[0] = p;
	rgb[1] = q;
	rgb[2] = v;
	break;
      }
    case 4:
      {
	rgb[0] = t;
	rgb[1] = p;
	rgb[2] = v;
	break;
      }
    case 5: default:
      {
	rgb[0] = v;
	rgb[1] = p;
	rgb[2] = q;
	break;
      }
    }

  
  char * returnMe = new char[3];
  
  for(int i = 0;i<3;i++)
    {
      assert( static_cast<int>(255.0*rgb[i]) < 256 );
      returnMe[i] = static_cast<unsigned char>(255.0*rgb[i]);
    }
  
  delete [] rgb;
  rgb = 0;
  return returnMe;
}
