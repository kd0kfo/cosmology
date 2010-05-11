#include "nfwshear.h"

int main(int argc, char** argv)
{
  try{
    return calculateShear(argc,argv);
  }
  catch(DavidException de)
    {
      de.stdOut();
      return de.getCode();
    }
}

int calculateShear(int argc, char** argv)
{

  using namespace std;

  if(argc == 1)
    throw DavidException("Need a file name");

  if(argc == 2)
    throw DavidException("Need the parameter filename");
  DString parameterFile(argv[2]);
  int parametersPosition = 0;
  utils::DArray<DString> parameters = parseParameterFile(parameterFile);

  bool calculateShear = true;
  bool calculateConvergence = false;
  if(argc >= 4)
    {
      if(DString(argv[3]) == "--deflectionangle")
	{
	  std::cout << "calculating deflection angle" << std::endl;
	  calculateShear = false;
	}
      else if(DString(argv[3]) == "--convergence")
	{
	  std::cout << "calculating convergence" << std::endl;
	  calculateShear = false;
	  calculateConvergence =true;
	}
    }

  ofstream myfile(argv[1], ios::out);
  DString space("  ");
  if(!myfile.is_open())
    throw DavidException(DString("Cannot open")+DString(argv[1]));

  int size = Double(parameters.get(parametersPosition++)).toInt();
  DEBUG_PRINT(size);
  double * shear1 = new double[size+1];
  double * shear2 = new double[size+1];

  double x0 = Double(parameters.get(parametersPosition++)).doubleValue();
  double y0 = Double(parameters.get(parametersPosition++)).doubleValue();
  double redshiftL = Double(parameters.get(parametersPosition++)).doubleValue();
  double redshiftS = Double(parameters.get(parametersPosition++)).doubleValue();  
  double sigma_crit = Cosmology::criticalDensity(redshiftL,redshiftS);
  double s,q,phi,theta;//s is the smallest to largest, q is middle to largest
  q = Double(parameters.get(parametersPosition++)).doubleValue();
  s = Double(parameters.get(parametersPosition++)).doubleValue();
  phi = Double(parameters.get(parametersPosition++)).doubleValue();
  theta = Double(parameters.get(parametersPosition++)).doubleValue();
  q2d = getq2d(s,q,phi,theta);//triaxial
  //q2d = getq2d(s=1,q=1,phi = 0.25*D_PI,theta = 0.25*D_PI);//sphere
  sqrtF = sqrt(getf(s,q,phi,theta));
  b=q2d;
  a = 1;
  rhos = 1;
  Rs = Double(parameters.get(parametersPosition++)).doubleValue();//in pixels
  double coefficient = 4*rhos*Rs*q2d/(sigma_crit*sqrtF);
  
  double x,y,averageCounter;
  shear1[0] = shear2[0] = 0;
  std::cout << "Calculating Shear" << std::endl;
  myfile << "#Rs = " << Rs << ", s = " << s << ", q = " << q << ", q2d = " << q2d << ", f = " << sqrtF*sqrtF << std::endl;
  myfile << "#Critical Density = " << sigma_crit << ", coefficient = " << coefficient << ", angles in radians = (" << phi << ", " << theta << ")" << std::endl;

  double * angles = new double[20];
  for(int i = 0;i<20;i++)
    angles[i] = i*(2*D_PI)/20;

  double gamma1,gamma2;
  for(int i = 1;i<= size;i++)
    {
      double radius = i/Rs;//integrand (i.e. sigmaval) is in terms of Rs
      averageCounter = 0;
      shear1[i] = 0;shear2[i] = 0;
      for(int j = 0;j < 20 ;j++)
	{
	  averageCounter++;
	  x = radius*cos(angles[j]);
	  y = radius*sin(angles[j]);

	  if(calculateConvergence)
	    {
	      shear1[i] += sigma(x,y);
	    }
	  else if(calculateShear)
	    {
	      gamma1 = 0.5*(daydy(x,y)-daxdx(x,y));;
	      gamma2 = daxdy(x,y);
	      shear1[i] += sqrt(gamma1*gamma1+gamma2*gamma2);
	    }
	  else
	    {
	      shear1[i] += ax(x,y);
	      shear2[i] += ay(x,y);
	    }

	}
	  shear1[i] = shear1[i]/averageCounter;//average of values over a quarter circle
	  shear2[i] = shear2[i]/averageCounter;
	  myfile << shear1[i] << std::endl;
	  DEBUG_PRINT(100.0*i/size << "%");
    }
  DEBUG_PRINT("Saving as " << argv[1]);

  //math::Complex curr;
  //  for(int i = 1;i< size+1;i++)
  // {
  //   curr = math::Complex(shear1[i],shear2[i]);
  //   myfile << curr.toDString() << std::endl;
  // }
  delete [] shear1;delete [] shear2;delete [] angles;shear1 = shear2 = angles = 0;
  return 0;
 }

 double getq2d(double s, double q, double theta, double phi)
 {
   double A = pow(cos(theta)/s,2)*(pow(sin(phi),2)+pow(cos(phi)/q,2))+pow(sin(theta)/q,2);
   double B = cos(theta)*sin(2*phi)*(1-1/(q*q))/(s*s);
   double C = (pow(sin(phi)/q,2)+pow(cos(phi),2))/(s*s);

   double radical = (A-C)*(A-C)+B*B;
   radical = sqrt(radical);

   double returnMe = A+C-radical;
   returnMe = returnMe/(A+C+radical);
   return sqrt(returnMe);
 }

 double getf(double s, double q, double theta, double phi)
 {
   return pow(sin(theta),2)*(pow(cos(phi),2)+pow(sin(phi)/q,2))+pow(cos(theta)/s,2);
 }

 double getlambda(double x , double y, double m)
 {
   double F = m*m*(a*a+b*b)-x*x-y*y;
   double G = m*m*(m*m*a*a*b*b-x*x*b*b-y*y*a*a);

   if(F*F < 4*G)
     throw DavidException("lambda is complex: F^2 < 4*G");


   if( (-F+sqrt(F*F-4*G))/2 > 0)
     return (-F+sqrt(F*F-4*G))/2;
   else
     throw DavidException("Lambda cannot be negative");// return (-F-sqrt(F*F-4*G))/2;
   
 }

 double simpsonRule(double x,double y, double a, double b, double f(double x, double y, double m))
 {
   throw DavidException("Don't use this.");
   double returnMe = f(x,y,a) + 4*f(x,y,(a+b)/2)+f(x,y,b);
   return returnMe*(b-a)/6;

 }

 double compSimpsonsRule(double x, double y, double a, double b, double f(double x, double y, double m))
 {
   double returnMe = f(x,y,a) + f(x,y,b);
   double n = 10000;
   double step = (b-a)/n;

   for(double i = 2;i<n;i += 2)
     {
       returnMe += 2*f(x,y,a+i*step);
     }

   for(double i = 1;i<n;i += 2)
     {
       returnMe += 4*f(x,y,a+i*step);
     }

   return returnMe*step/3;
 }

 double integral(double x, double y, double m0, double m, int numberOfSteps,double integrand(double x, double y, double m) )
 {
   double returnMe = 0;
   double stepSize = (m-m0)/numberOfSteps;
   for(double i = m0;i<m;i = i + stepSize)
     {
       returnMe += compSimpsonsRule(x,y,i,i+stepSize,integrand);
       //returnMe += simpsonRule(x, y, i,i+stepSize,integrand);
       //returnMe += integrand(i)*stepSize;
     }

   return returnMe;
 }

 double sigma(double x, double y)
 {
   return sigmaval(x,y,sqrt(x*x+y*y/(q2d*q2d)));
 }

double totalmass(double x, double y)
{
  return compSimpsonsRule(x,y,0.000001,sqrt(x*x+y*y/(q2d*q2d)),sigmaval);
}

double ax(double x,double y)
{
  return compSimpsonsRule(x,y,0.000001,sqrt(x*x+y*y/(q2d*q2d))-0.00001,ax);
}

double ax(double x, double y, double m)
{
  return Acoefficient(x,y,m)*sigmaval(x,y,m)*m;
}

double ay(double x,double y)
{
  return compSimpsonsRule(x,y,0.000001,sqrt(x*x+y*y/(q2d*q2d))-0.00001,ay);
}

double ay(double x, double y, double m)
{
  return Bcoefficient(x,y,m)*sigmaval(x,y,m)*m;
}

 double daxdx(double x,double y)
 {
   double epsilon = 0.000001;
   double m = sqrt(x*x+y*y/(q2d*q2d))-epsilon;
   double limitDerivative = x/m;
   limitDerivative *= Acoefficient(x,y,m)*sigmaval(x,y,m)*m;
   
   return compSimpsonsRule(x,y,epsilon,m,daxdx) + limitDerivative;
 }

 double daxdx(double x, double y, double m)
 {
   double sigma = sigmaval(x,y,m);
   return dAdx(x,y,m)*sigma*m;
 }

 double daxdy(double x, double y)
 {
   double epsilon = 0.000001;
   double m = sqrt(x*x+y*y/(q2d*q2d))-epsilon;
   double limitDerivative = (y/(q2d*q2d*m))*Acoefficient(x,y,m)*sigmaval(x,y,m)*m;

   return compSimpsonsRule(x, y, epsilon, m,daxdy) + limitDerivative;
 }

 double daxdy(double x, double y, double m)
 {
   double sigma = sigmaval(x,y,m);
   return dAdy(x,y,m)*sigma*m;
 }

 double daydy(double x, double y)
 {
   double epsilon = 0.000001;
   double m = sqrt(x*x+y*y/(q2d*q2d))-epsilon;
   double limitDerivative = (y/(q2d*q2d*m))*Bcoefficient(x,y,m)*sigmaval(x,y,m)*m;

   return compSimpsonsRule(x, y, epsilon, m, daydy) + limitDerivative;
 }

 double daydy(double x, double y, double m)
 {
   double sigma = sigmaval(x,y,m);
   return dBdy(x,y,m)*sigma*m;
 }

 double daydx(double x, double y)
 {
   double epsilon = 0.000001;
   double m = sqrt(x*x+y*y/(q2d*q2d))-epsilon;
   double limitDerivative = (x/m)*Bcoefficient(x,y,m)*sigmaval(x,y,m)*m;
   return compSimpsonsRule(x, y, epsilon, m, daydx ) + limitDerivative;
 }
 double daydx(double x, double y, double m)
 {
   double sigma = sigmaval(x,y,m);
   return dBdx(x,y,m)*sigma*m;
 }


 double sigmaval(double x, double y, double m)
 {
   double stuff = (m<1) ? sqrt(1-m*m) : sqrt(m*m -1);
   double returnMe;
   if(m < 1)
     returnMe = atanh(stuff)-stuff;//factor of 2: mistake in flores' notes?
   else
     returnMe = stuff-atan(stuff);

   return returnMe/pow(stuff,3);

 }

 double dAdx(double x, double y, double m)
 {

   double * ap = new double[4];
   double * bp = new double[4];
   ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
   bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));

   for(int i = 1;i<4;i++)
     {
       ap[i] = ap[i-1]*ap[0];
       bp[i] = bp[i-1]*bp[0];
     }

   double l = pow(x*bp[1],2)+pow(y*ap[1],2);

   double aprimeTerm = bp[2]*x/l - 4*ap[3]*bp[2]*x*y*y/(l*l);
   aprimeTerm *= daprimedx(x,y,m);

   double bprimeTerm = 3*ap[0]*bp[1]*x/l - 4*ap[0]*bp[3]*bp[1]*pow(x,3)/(l*l);
   bprimeTerm *= dbprimedx(x,y,m);

   double returnMe =  aprimeTerm + bprimeTerm + ap[0]*bp[2]/l - ap[0]*bp[2]*2*x*x*bp[3]/(l*l);
   
   delete [] ap;delete [] bp;
   return returnMe;
 }

 double dAdy(double x, double y, double m)
 {
   double * ap = new double[4];
   double * bp = new double[4];
   ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
   bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));
   for(int i = 1;i<4;i++)
     {
       ap[i] = ap[i-1]*ap[0];
       bp[i] = bp[i-1]*bp[0];
     }

   double l = pow(x*bp[1],2)+pow(y*ap[1],2);

   double aprimeTerm = bp[2]*x/l - 4*ap[3]*bp[2]*x*y*y/(l*l);
   aprimeTerm *= daprimedy(x,y,m);

   double bprimeTerm = 3*ap[0]*bp[1]*x/l - 4*ap[0]*bp[3]*bp[1]*pow(x,3)/(l*l);
   bprimeTerm *= dbprimedy(x,y,m);

   double returnMe = aprimeTerm + bprimeTerm - 2*x*y*ap[3]*ap[0]*bp[2]/(l*l);

   delete [] ap;delete [] bp;
   return returnMe;  
 }



 double dBdx(double x, double y, double m)
 {
   double * ap = new double[4];
   double * bp = new double[4];

   ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
   bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));

   for(int i = 1;i<4;i++)
     {
       ap[i] = ap[i-1]*ap[0];
       bp[i] = bp[i-1]*bp[0];
     }

   double l = pow(x*bp[1],2)+pow(y*ap[1],2);

   double aprimeTerm = 3*ap[1]*bp[0]*y/l - 4*bp[0]*ap[3]*ap[1]*pow(y,3)/(l*l);
   aprimeTerm *= daprimedx(x,y,m);

   double bprimeTerm = ap[2]*y/l - 4*ap[2]*bp[3]*x*x*y/(l*l);
   bprimeTerm *= dbprimedx(x,y,m);

   double returnMe = aprimeTerm + bprimeTerm - 2*ap[2]*bp[3]*bp[0]*x*y/(l*l);

   delete [] ap;delete [] bp;
   return returnMe;
 }


 double dBdy(double x, double y, double m)
 {
   double * ap = new double[4];
   double * bp = new double[4];

   ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
   bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));

   for(int i = 1;i<4;i++)
     {
       ap[i] = ap[i-1]*ap[0];
       bp[i] = bp[i-1]*bp[0];
     }

   double l = pow(x*bp[1],2)+pow(y*ap[1],2);

   double aprimeTerm = 3*ap[1]*bp[0]*y/l - 4*bp[0]*ap[3]*ap[1]*pow(y,3)/(l*l);
   aprimeTerm *= daprimedy(x,y,m);

   double bprimeTerm = ap[2]*y/l - 4*ap[2]*bp[3]*x*x*y/(l*l);
   bprimeTerm *= dbprimedy(x,y,m);

   double returnMe = aprimeTerm + bprimeTerm + ap[2]*bp[0]*(1-2*ap[3]*y*y/l)/l;
  

   delete [] ap; delete [] bp;
   return returnMe;
}

double daprimedx(double x, double y, double m)
{
  return dlambdadx(x,y,m)/(2*sqrt(m*m*a*a+getlambda(x,y,m)));
}

double daprimedy(double x, double y, double m)
{
  return dlambdady(x,y,m)/(2*sqrt(m*m*a*a+getlambda(x,y,m)));
}

double dbprimedx(double x, double y, double m)
{
  return dlambdadx(x,y,m)/(2*sqrt(m*m*b*b+getlambda(x,y,m)));
}

double dbprimedy(double x, double y, double m)
{
  return dlambdady(x,y,m)/(2*sqrt(m*m*b*b+getlambda(x,y,m)));
}

double dlambdadx(double x, double y, double m)
{
  //This is what I found by simplifying the radical in lambda
  /**
  double denominator = pow(m*m*(a*a-b*b),2)+pow(x*x*+y*y,2)-2*m*m*(a*a-b*b)*(x*x-y*y);
  denominator = sqrt(denominator);
  double returnMe = x*(x*x+y*y)-x*m*m*(a*a-b*b)*(x*x-y*y);
  returnMe = returnMe/denominator;
  return returnMe + x;/**/

  //if(x*x+y*y/(q2d*q2d) < 1)
  //return 0;

  //double F = (a*a+b*b)-x*x-y*y;
  double F = m*m*(a*a+b*b)-x*x-y*y;
  //double G = a*a*b*b-x*x*b*b-y*y*a*a;
  double G = m*m*(m*m*a*a*b*b-x*x*b*b-y*y*a*a);
  
  double dFdx = -2*x;
  double dGdx = -2*m*m*b*b*x;
  //double dGdx = -2*b*b*x;
  
  if(F*F > 4*G)
    return -0.5*dFdx+0.5*(F*dFdx-2*dGdx)/sqrt(F*F-4*G);
  else
    throw DavidException("Imaginary Lambda in dlambdadx");
}

double dlambdady(double x, double y, double m)
{
  //This is what I found by simplifying the radical in lambda
  /**
  double denominator = pow(m*m*(a*a-b*b),2)+pow(x*x*+y*y,2)-2*m*m*(a*a-b*b)*(x*x-y*y);
  denominator = sqrt(denominator);
  double returnMe = y*(x*x+y*y)+y*m*m*(a*a-b*b)*(x*x-y*y);
  returnMe = returnMe/denominator;
  return returnMe + y;
  /**/

  //if(x*x+y*y/(q2d*q2d) < 1)
  //  return 0;

  //double F = (a*a+b*b)-x*x-y*y;
  double F = m*m*(a*a+b*b)-x*x-y*y;
  //double G = a*a*b*b-x*x*b*b-y*y*a*a;
  double G = m*m*(m*m*a*a*b*b-x*x*b*b-y*y*a*a);

  double dFdy = -2*y;
  //double dGdy = -2*a*a*y;
  double dGdy = -2*m*m*a*a*y;
  
  if(F*F > 4*G)
    return -0.5*dFdy+0.5*(F*dFdy-2*dGdy)/sqrt(F*F-4*G);
  else
    throw DavidException("Imaginary Lambda in dlambdady");
}

double Acoefficient(double x, double y, double m)
{
  double * ap = new double[4];
  double * bp = new double[4];

   ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
   bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));

   for(int i = 1;i<4;i++)
     {
       ap[i] = ap[i-1]*ap[0];
       bp[i] = bp[i-1]*bp[0];
     }

   double l = pow(x*bp[1],2)+pow(y*ap[1],2);
   
   double returnMe = ap[0]*bp[2]*x/l;
   delete [] ap;delete [] bp;
   return returnMe;
}

double Bcoefficient(double x, double y, double m)
{
  double * ap = new double[4];
  double * bp = new double[4];

  ap[0] = sqrt(m*m*a*a+getlambda(x,y,m));
  bp[0] = sqrt(m*m*b*b+getlambda(x,y,m));

  for(int i = 1;i<4;i++)
    {
      ap[i] = ap[i-1]*ap[0];
      bp[i] = bp[i-1]*bp[0];
    }

  double l = pow(x*bp[1],2)+pow(y*ap[1],2);
   
  double returnMe = ap[2]*bp[0]*y/l;
  delete [] ap;delete [] bp;
  return returnMe;
}


double romberg(double func(double z), double a, double b)
{
  double h = b - a;     // coarsest panel size
  double dR;			      // convergence
  int np = 1;           // Current number of panels
  const int N = 25;     // maximum iterations
  double prec = 1e-8;	  // desired precision
  double * R = new double[N*N];

  for(int i = 0;i<N*N;i++)
    R[i] = 0;

  // Compute the first term R(1,1)
  R[0] = h/2 * (func(a) + func(b));

  // Loop over the desired number of rows, i = 2,...,N
  int i,j,k;
  for( i=1; i<N; ++i )
  {
    // Compute the summation in the recursive trapezoidal rule
    h /= 2.0;          // Use panels half the previous size
    np *= 2;           // Use twice as many panels
    double sumT = 0.0;
    for( k=1; k<=(np-1); k+=2 ) 
      sumT += func( a + k*h);

    // Compute Romberg table entries R(i,1), R(i,2), ..., R(i,i)
    R[N*i] = 0.5 * R[N*(i-1)] + h * sumT;   
    int m = 1;
    for( j=1; j<i; ++j )
    {
      m *= 4;
      R[N*i+j] = R[N*i+j-1] + (R[N*i+j-1] - R[N*(i-1)+j-1]) / (m-1);
    }
    dR = (j > 1) ? R[N*i+j-1] - R[N*(i-1)+(j-2)] : R[0];
    if (fabs(dR) < prec)
      return R[N*i+j-1];
  }
  return R[N*N-1];
}


utils::DArray<DString> parseParameterFile(DString fileName)
{
  std::ifstream infile(fileName.toCharArray());
  using namespace std;
		
  int stringSize = 150;
  int totalParams;
	
  utils::DArray<DString> paramsVector;
  int vectorCount = 0;

  if(infile.is_open())
    {
      Double tmp;

      char * curr = new char[stringSize];
      infile.getline(curr,stringSize);
      DString currDString(curr);
      if(currDString.contains(';'))
		currDString = currDString.substring(0,currDString.indexOf(";"));
      tmp = Double::parseDString(currDString);
      totalParams = (int) tmp.doubleValue();

      for(int i = 0;i< totalParams;i++)
	{
	  curr = new char[stringSize];
	  infile.getline(curr,stringSize);
	  currDString = curr;
	  if(currDString.contains(';'))
	    currDString = currDString.substring(0,currDString.indexOf(";"));
	  paramsVector.put(currDString);
	}

    }
  else
    throw DavidException(DString("Could not open") + fileName,DavidException::IO_ERROR_CODE);

  return paramsVector;
}
