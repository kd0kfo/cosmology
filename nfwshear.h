#ifndef __NFWSHEAR__
#define __NFWSHEAR__

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "libdnstd/DavidException.h"
#include "libdnstd/Complex.h"
#include "libdnstd/Cosmology.h"
#include "libdnstd/DArray.cpp"

template class utils::DArray<DString>;

#ifdef __DEBUG__
#define DEBUG_PRINT(X) std::cout << X << std::endl;
#else
#define DEBUG_PRINT(X) 1;
#endif

int calculateShear(int argc, char** argv);

double ax(double x, double y);
double daxdx(double x,double y);///<Derviative of Alpha_x wrt x at the point (x,y)
double daxdy(double x,double y);///<Derviative of Alpha_x wrt y at the point (x,y)
double daydx(double x,double y);///<Derviative of Alpha_y wrt x at the point (x,y)
double daydy(double x,double y);///<Derviative of Alpha_y wrt y at the point (x,y)
double totalmass(double x,double y);///<Calculates the total mass within an ellipse that contains the point (x,y)
double sigma(double x, double y);///<Calculates the mass density at the point (x,y)
double ax(double x, double y);///<Calculates Alpha_x at the point (x,y)
double ay(double x, double y);///<Calculates Alpha_y at the point (x,y)

//Functions in the integrand(s)
double daxdx(double x, double y, double m);
double daxdy(double x, double y, double m);
double daydx(double x, double y, double m);
double daydy(double x, double y, double m);
double dAdx(double x, double y, double m);
double dAdy(double x, double y, double m);
double dBdx(double x, double y, double m);
double dBdy(double x, double y, double m);
double sigmaval(double x, double y, double m);
double daprimedx(double x, double y, double m);
double daprimedy(double x, double y, double m);
double dbprimedx(double x, double y, double m);
double dbprimedy(double x, double y, double m);
double dlambdadx(double x, double y, double m);
double dlambdady(double x, double y, double m);
double getq2d(double s, double q, double theta, double phi);
double getf(double s, double q, double theta, double phi);
double getlambda(double x , double y, double m);
double ax(double x, double y, double m);
double ay(double x, double y, double m);
double Acoefficient(double x, double y, double m);
double Bcoefficient(double x, double y, double m);
//integration routines
double integral(double x, double y, double m0, double m, int numberOfSteps,double integrand(double x, double y, double m) );
double simpsonRule(double x, double y, double a, double b, double f(double x, double y, double m) );
double compSimpsonsRule(double x, double y, double a, double b, double f(double x, double y, double m) );
double romberg(double func(double z), double a, double b);

utils::DArray<DString> parseParameterFile(DString fileName);

//Constants
double a;
double b;
double q2d;
double rhos;
double Rs;
double sqrtF;//function of above constants to be set @ beginning of program

#endif

