# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
from .types cimport *

cdef extern from "gsl/gsl_sf_dawson.h":

  double  gsl_sf_dawson(double x)

  int  gsl_sf_dawson_e(double x, gsl_sf_result * result)

