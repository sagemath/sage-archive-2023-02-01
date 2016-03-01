# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
from .types cimport *

cdef extern from "gsl/gsl_sf_elementary.h":

  int  gsl_sf_multiply_e(double x, double y, gsl_sf_result * result)

  int  gsl_sf_multiply_err_e(double x, double dx, double y, double dy, gsl_sf_result * result)

