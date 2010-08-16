#include <pari/pari.h>
#include "interrupt.h"


// global catch variable !
// this means that the code is not reentrant -- beware !
// THAT MEANS NO CALLING TO PARI from inside the trap....
// Should replace by a stack, that would work.
static long __catcherr = 0;

#define _pari_raise(errno) { \
        PyErr_SetObject(PyExc_PariError, PyInt_FromLong(errno)); \
    }

#define _pari_endcatch { err_leave(__catcherr); }

/* Careful with pari_retries !
 * It should not be in a register, we flag it as "volatile".
 */
#define _pari_catch { \
        long pari_errno; \
        long volatile pari_retries = 0; \
        jmp_buf __env; \
        __catcherr = 0; \
        if ((pari_errno=setjmp(__env))) { \
            _pari_trap(pari_errno, pari_retries); \
            if(PyErr_Occurred()) { \
                _pari_endcatch; \
                return NULL; \
            } \
            pari_retries++; \
        } \
        __catcherr = err_catch(CATCH_ALL, &__env); \
    }
