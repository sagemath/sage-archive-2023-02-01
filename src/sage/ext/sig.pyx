include 'python_exc.pxi'
include 'stdsage.pxi'

## cdef extern from "Python.h":
##     PyObject* PyExc_KeyboardInterrupt

## cdef extern from "signal.h":
##     int SIGSEGV
##     int SIGBUS
##     int SIGFPE
##     ctypedef void (*sighandler_t)(int)
##     void signal(int, sighandler_t)

## cdef extern from "stdio.h":
##     struct FILE
##     int fprintf(FILE *stream, char *format, ...)
##     FILE* stderr

## cdef extern from "stdlib.h":
##     void exit(int)

## cdef void msg(char* s):
##     fprintf(stderr, '\n\n------------------------------------------------------------\n')
##     fprintf(stderr, s)
##     fprintf(stderr, "This probably occured because a *compiled* component\n")
##     fprintf(stderr, "of SAGE has a bug in it (typically accessing invalid memory)\n")
##     fprintf(stderr, "or is not properly wrapped with _sig_on, _sig_off.\n")
##     fprintf(stderr, "You might want to run SAGE under gdb with 'sage -gdb' to debug this.\n")
##     fprintf(stderr, "SAGE will now terminate (sorry).\n")
##     fprintf(stderr, '------------------------------------------------------------\n\n')

## cdef void sig_handle_sigsegv(int n):
##     msg("Unhandled SIGSEGV: A segmentation fault occured in SAGE.\n")
##     PyErr_SetString(<object>PyExc_KeyboardInterrupt, "")
##     exit(1)

## cdef void sig_handle_sigbus(int n):
##     msg("Unhandled SIGBUS: A bus error occured in SAGE.\n")
##     PyErr_SetString(<object>PyExc_KeyboardInterrupt, "")
##     exit(1)

## cdef void sig_handle_sigfpe(int n):
##     msg("Unhandled SIGFPE: An unhandled floating point exception occured in SAGE.\n")
##     PyErr_SetString(<object>PyExc_KeyboardInterrupt, "")
##     exit(1)

def get_bad_sigs():
#    signal(SIGSEGV, sig_handle_sigsegv)
#    signal(SIGBUS, sig_handle_sigbus)
#    signal(SIGFPE, sig_handle_sigfpe)
    init_csage()

