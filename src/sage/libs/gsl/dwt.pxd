from sage.libs.gsl.wavelet cimport *
from .array cimport GSLDoubleArray

cdef class DiscreteWaveletTransform(GSLDoubleArray):
    cdef gsl_wavelet* wavelet
    cdef gsl_wavelet_workspace* workspace
