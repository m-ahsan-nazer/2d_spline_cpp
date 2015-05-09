#ifndef __PYX_HAVE__scipy_tools
#define __PYX_HAVE__scipy_tools


#ifndef __PYX_HAVE_API__scipy_tools

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

__PYX_EXTERN_C DL_IMPORT(void) set_1d_spline(spline_1d_param *, spline_1d_param *, double *, double *, int, int);
__PYX_EXTERN_C DL_IMPORT(void) set_2d_spline(spline_2d_param *, spline_2d_param *, spline_2d_param *, double *, double *, double *, int, int, int);

#endif /* !__PYX_HAVE_API__scipy_tools */

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initscipy_tools(void);
#else
PyMODINIT_FUNC PyInit_scipy_tools(void);
#endif

#endif /* !__PYX_HAVE__scipy_tools */
