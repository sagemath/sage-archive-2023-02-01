# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
from .types cimport *

cdef extern from "gsl/gsl_diff.h":
  int gsl_diff_central ( gsl_function *f, double x,
                        double *result, double *abserr)

  int gsl_diff_backward ( gsl_function *f, double x,
                         double *result, double *abserr)

  int gsl_diff_forward ( gsl_function *f, double x,
                        double *result, double *abserr)

