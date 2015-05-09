#!python2.7
#cython: boundscheck=False, wraparound=False
from scipy.interpolate import UnivariateSpline as US
from scipy.interpolate import RectBivariateSpline as RBS
import numpy as np
cimport numpy as np
from cython.operator cimport dereference as deref
from libc.string cimport memcpy
from cython.view cimport array as cvarray

cdef extern from "numpy/arrayobject.h":
	void * PyArray_DATA(np.ndarray arr) #(PyArrayObject *arr)


cdef extern from "splev_bispeu.h":
	cdef cppclass spline_1d_param:
		int n
		double *t
		double *c
		int k
		void operator()(int)
	cdef cppclass spline_2d_param:
		int nx
		double *tx
		int xy
		double *ty
		int nc
		double *c
		int kx
		int ky
		#(int nx_, int ny_, int nc_, int kx_, int ky_)
		void operator()(int, int, int, int, int)
#
cdef public void set_1d_spline(spline_1d_param * sp_param,
                               spline_1d_param * sp_deriv_param, 
                              double *x_, double *y_, int size,int smoothing):
	cdef np.ndarray[dtype=np.double_t,ndim=1] x = np.ascontiguousarray(<np.double_t[:size]> x_)
	cdef np.ndarray[dtype=np.double_t,ndim=1] y = np.ascontiguousarray(<np.double_t[:size]> y_)
	cdef int s, nu
	s = smoothing
	nu = 1
	sp = US(x,y,s=smoothing)
	sp_deriv = US(x,sp(x,nu=nu),s=s)
	#n,t,c,k,ier = data[7],data[8],data[9],data[5],data[-1]
	deref(sp_param)(sp._data[7])
	sp_param.n = sp._data[7]
	sp_param.k = sp._data[5]
	deref(sp_deriv_param)(sp_deriv._data[7])
	sp_deriv_param.n = sp_deriv._data[7]
	sp_deriv_param.k = sp_deriv._data[5]
	
	cdef int i
	for i in xrange(sp_param.n):
		sp_param.t[i] = sp._data[8][i]
		sp_param.c[i] = sp._data[9][i]
		sp_deriv_param.t[i] = sp_deriv._data[8][i]
		sp_deriv_param.c[i] = sp_deriv._data[9][i]
	"""
	#memcpy is faster
	memcpy(sp_param.t,<double *> PyArray_DATA(sp._data[8]),sizeof(double)*sp._data[8].size)
	memcpy(sp_param.c,<double *> PyArray_DATA(sp._data[9]),sizeof(double)*sp._data[9].size)
	memcpy(sp_deriv_param.t,<double *> PyArray_DATA(sp_deriv._data[8]),
	       sizeof(double)*sp_deriv._data[8].size)
	memcpy(sp_deriv_param.c,<double *> PyArray_DATA(sp_deriv._data[9]),
	       sizeof(double)*sp._deriv_data[9].size)
	"""




cdef public void set_2d_spline(spline_2d_param * sp_param,
                               spline_2d_param * spdx_param,
                               spline_2d_param * spdy_param, 
                               double *x_, double *y_,double *z_, 
                               int size_x,int size_y,int smoothing):
	cdef np.ndarray[dtype=np.double_t,ndim=1] x = np.ascontiguousarray(<np.double_t[:size_x]> x_)
	cdef np.ndarray[dtype=np.double_t,ndim=1] y = np.ascontiguousarray(<np.double_t[:size_y]> y_)
	cdef np.ndarray[dtype=np.double_t,ndim=1] z = \
	                      np.ascontiguousarray(<np.double_t[:size_x*size_y]> z_)
	
	cdef int s, dx, dy
	s=smoothing
	sp = RBS(x,y,z.reshape((size_x,size_y)),s=s)
	dx = 1
	dy = 0
	spdx = RBS(x,y,sp(x,y,dx=dx,dy=dy,grid=True),s=s)
	dx = 0
	dy = 1
	spdy = RBS(x,y,sp(x,y,dx=dx,dy=dy,grid=True),s=s)
	##(int nx_, int ny_, int nc_, int kx_, int ky_)
	cdef int nx_, ny_, nc_, kx_, ky_
	nx_ = sp.tck[0].size
	ny_ = sp.tck[1].size
	nc_ = sp.tck[2].size
	kx_ , ky_ = sp.degrees 
	deref(sp_param)(nx_,ny_,nc_,kx_,ky_)
	deref(spdx_param)(nx_,ny_,nc_,kx_,ky_)
	deref(spdy_param)(nx_,ny_,nc_,kx_,ky_)
	"""
	cdef int i
	for i in xrange(nx_):
		sp_param.tx[i]   = sp.tck[0][i]
		spdx_param.tx[i] = spdx.tck[0][i]
		spdy_param.tx[i] = spdy.tck[0][i]

	for i in xrange(ny_):
		sp_param.ty[i]   = sp.tck[1][i]
		spdx_param.ty[i] = spdx.tck[1][i]
		spdy_param.ty[i] = spdy.tck[1][i]

	for i in xrange(nc_):
		sp_param.c[i]   = sp.tck[2][i]
		spdx_param.c[i] = spdx.tck[2][i]
		spdy_param.c[i] = spdy.tck[2][i]
	"""
	#memcpy is faster
	memcpy(sp_param.tx,<double *> PyArray_DATA(sp.tck[0]),sizeof(double)*sp.tck[0].size)
	memcpy(spdx_param.tx,<double *> PyArray_DATA(spdx.tck[0]),sizeof(double)*spdx.tck[0].size)
	memcpy(spdy_param.tx,<double *> PyArray_DATA(spdy.tck[0]),sizeof(double)*spdy.tck[0].size)
	
	memcpy(sp_param.ty,<double *> PyArray_DATA(sp.tck[1]),sizeof(double)*sp.tck[1].size)
	memcpy(spdx_param.ty,<double *> PyArray_DATA(spdx.tck[1]),sizeof(double)*spdx.tck[1].size)
	memcpy(spdy_param.ty,<double *> PyArray_DATA(spdy.tck[1]),sizeof(double)*spdy.tck[1].size)

	memcpy(sp_param.c,<double *> PyArray_DATA(sp.tck[2]),sizeof(double)*sp.tck[2].size)
	memcpy(spdx_param.c,<double *> PyArray_DATA(spdx.tck[2]),sizeof(double)*spdx.tck[2].size)
	memcpy(spdy_param.c,<double *> PyArray_DATA(spdy.tck[2]),sizeof(double)*spdy.tck[2].size)

	
