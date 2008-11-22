
include "python.pxi"

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void free(void *ptr)
    void *malloc(size_t size)
    void *realloc(void *ptr, size_t size)
    size_t strlen(char *s)
    char *strcpy(char *dest, char *src)

cdef extern from "string.h":
    void *memset(void *dest, int c, size_t n)
    void *memcpy(void *dest, void *src, size_t n)

cdef extern from "stdio.h":
    ctypedef struct FILE
    int printf(char *format, ...)
    int fprintf(FILE *stream, char *format, ...)
    int sprintf(char *str, char *format, ...)
    FILE *fopen(char *path, char *mode)
    int fclose(FILE *strea)
    cdef FILE *stdout
    int scanf(char *format, ...)

cdef extern from "math.h":
    double sqrt(double x)
    float roundf(float x)    # linux-ish and non-standard; avoid!
    double ldexp(double x, int exp)
    double frexp(double x, int *exp)

from sage.libs.gmp.all cimport *
cdef extern from "gmp.h":
    pass # cython bug sometimes includes this in the wrong place

##########################################################################
# stdsage.pxi declares the macros, etc., that got used a lot in SAGE.
##########################################################################

