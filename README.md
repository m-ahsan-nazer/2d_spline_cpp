# 2d_spline_cpp

The gsl library does not have any 2d interpolation routines. Here I am  have provided a c++ interface to the scipy 2d interpolation routine RectBivariateSpline via cython. Scipy uses fitpack by P. Dierckx see http://www.netlib.org/dierckx/ . The two fortran subroutines I have wrapped are splev and bispeu. I could not find bispeu at netlib and the fortran code here has been taken from the scipy source files at https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack. 

The 1d and 2d splines are setup at python/numpy/cython speed but evaluated at c++/fortran speed. Personnaly I use the 2d splines with differential equation solvers which require repeated spline evaluation at a single point but it is possibly to evaluate the splines at (x,y) where x and y are arrays with a single call to the fortran subroutines. You will have to modify the classes scipy_spline_1d and scipy_spline_2d. 

During the fortran file compilation a warning is issued "Actual argument contains too few elements for dummy argument 'h'". If you wish to get rid of the warning then change the size of array h. I have left it the way it was in the scipy source files.  
