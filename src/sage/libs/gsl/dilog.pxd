# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
from .types cimport *

cdef extern from "gsl/gsl_sf_dilog.h":

  double  gsl_sf_dilog(double x)

  int  gsl_sf_dilog_e(double x, gsl_sf_result * result)

  int  gsl_sf_complex_dilog_e(double r, double theta, gsl_sf_result * result_re, gsl_sf_result * result_im)

