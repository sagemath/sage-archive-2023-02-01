
cdef extern from "gsl/gsl_sf_result.h":
  ctypedef struct gsl_sf_result:
    double val
    double err


  ctypedef struct gsl_sf_result_e10:
    double val
    double err
    int    e10

