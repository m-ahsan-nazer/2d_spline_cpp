from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

extension_modules =  Extension("scipy_tools", 
              sources=["scipy_tools.pyx"],
              language="c++",
              )
setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [extension_modules],
    include_dirs=[np.get_numarray_include(),np.get_include()]
)

