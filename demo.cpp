#include <iostream>
#include <boost/format.hpp>
#include <fstream>
#include "demo.h"

double quad(double x)
{
	return 3.*x*x + 1;
}
double dquad_dx(double x)
{
	return 6.*x;
}


double fun(double x,double y)
{
	return exp(-x*x+3.*y);
}

double dfundx(double x,double y)
{
	return -2.*x*exp(-x*x+3.*y);
}

double dfundy(double x,double y)
{
	return 3.*exp(-x*x+3.*y);
}

using namespace std;
int main()
{
	Py_Initialize();
    initscipy_tools();
    
    int size_x = 170;
    double x[size_x];
    double y[size_x];
    for(int i=0; i<size_x; i++)
    {
    	x[i] = 0.+ ((double) i)/((double) size_x-1)*8.;
    	y[i] = quad(x[i]); 
    }
	
	scipy_spline_1d sp(&x[0], &y[0], size_x, 0);
	
	ofstream f("spline.dat");
	f << "x   y   y_sp   dydx   dydx_sp" << endl;
	for (int i=0; i<1000; i++)
	{
		double xval = x[0] + ((double) i)/((double) 1000-1)*8.;
		f << boost::format("%.6e   %.6e   %.6e   %.6e   %.6e\n") 
		%xval %quad(xval) %sp(xval,0) %dquad_dx(xval) %sp(xval,1);
               
	}
	
	f.close();
	/**
	Now test 2d splines
	**/
	int size_x2d = 700;
	int size_y2d = 2000;
	double *x2d;
	double *y2d;
	double *z;
	x2d = new double[size_x2d];
	y2d = new double[size_y2d];
	z = new double[size_x2d*size_y2d];
	
	for (int i=0; i<size_x2d; i++)
	{
		x2d[i] = -5. + ((double)i)/((double) size_x2d -1)*5.;
	}
	
	for (int i=0; i<size_y2d; i++)
	{
		y2d[i] = -6. + ((double)i)/((double) size_y2d -1)*7.;
	}
	for (int i=0; i<size_x2d; i++)
	{
		for (int j=0; j<size_y2d; j++)
		{
			z[size_y2d*i+j] = fun(x2d[i],y2d[j]);
		}
	}
	
	scipy_spline_2d sp2d(x2d, y2d,z,size_x2d,size_y2d,0);
	
	ofstream f2d("2dspline.dat");
	f2d << "x   y   z   z_sp   dzdx   dzdx_sp   dzdy   dzdy_sp" << endl;
	for (int i=0; i<1000; i++)
	{
		double xval = x2d[0] + ((double) i)/((double) 1000-1)*x2d[size_x2d-1];
		double yval = y2d[0] + ((double) i)/((double) 1000-1)*y2d[size_y2d-1];
		
		f2d << boost::format("%.6e   %.6e   %.6e   %.6e   %.6e   %.6e   %.6e   %.6e\n") 
		%xval %yval %fun(xval,yval) %sp2d(xval,yval,0,0) 
		%dfundx(xval,yval) %sp2d(xval,yval,1,0)
		%dfundy(xval,yval) %sp2d(xval,yval,0,1);
               
	}
	f2d.close();
	
	delete [] x2d; delete [] y2d; delete [] z;
	Py_Finalize(); 
	return 0;
}
