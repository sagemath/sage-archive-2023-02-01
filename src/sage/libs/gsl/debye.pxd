# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
from .types cimport *

cdef extern from "gsl/gsl_sf_debye.h":

  double  gsl_sf_debye_1(double x)

  int  gsl_sf_debye_1_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_2(double x)

  int  gsl_sf_debye_2_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_3(double x)

  int  gsl_sf_debye_3_e(double x, gsl_sf_result * result)

  double  gsl_sf_debye_4(double x)

  int  gsl_sf_debye_4_e(double x, gsl_sf_result * result)

