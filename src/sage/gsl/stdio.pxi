cdef extern from "stdio.h":
  ctypedef struct FILE
  FILE *fopen(char *path, char *mode)
  int fclose(FILE *strea)
  cdef FILE *stdout
  int scanf(char *format, ...)

