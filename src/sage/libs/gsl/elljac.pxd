# distutils: libraries = GSL_LIBRARIES
# distutils: library_dirs = gsl_library_dirs
# distutils: include_dirs = gsl_include_dirs
cdef extern from "gsl/gsl_sf_elljac.h":

  int  gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn)

