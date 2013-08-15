#include <pari/pari.h>
#include "interrupt.h"


// global catch variable !
// this means that the code is not reentrant -- beware !
// THAT MEANS NO CALLING TO PARI from inside the trap....
// Should replace by a stack, that would work.
static volatile long sage_pari_catcherr = 0;

#define _pari_raise(errno) {                                                  \
        PyErr_SetObject(PyExc_PariError, PyInt_FromLong(errno));              \
    }

#define _pari_endcatch {                                                      \
         err_leave(sage_pari_catcherr);                                       \
    }

/* Careful with pari_retries, it must be declared volatile!
 * Also note that "if (pari_errno=setjmp(...))" is not legal C/C++
 */
#define _pari_catch                                                           \
    jmp_buf _pari_catch_env;                                                  \
    {                                                                         \
        volatile long pari_retries = 0;                                       \
        sage_pari_catcherr = 0;                                               \
        long pari_errno = setjmp(_pari_catch_env);                            \
        if (pari_errno) {                                                     \
            _pari_trap(pari_errno, pari_retries);                             \
            if(PyErr_Occurred()) {                                            \
                _pari_endcatch;                                               \
                return NULL;                                                  \
            }                                                                 \
            pari_retries++;                                                   \
        }                                                                     \
        sage_pari_catcherr = err_catch(CATCH_ALL, &_pari_catch_env);          \
    }
