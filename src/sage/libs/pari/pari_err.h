#include <pari/pari.h>
#include "interrupt.h"


/*****************************************
   Interrupts and PARI exception handling
 *****************************************/
#define pari_catch_sig_on() sig_on(); _pari_catch;
#define pari_catch_sig_str(s) sig_str(s); _pari_catch;
#define pari_catch_sig_off() _pari_endcatch; sig_off();


// global catch variable !
// this means that the code is not reentrant -- beware !
// THAT MEANS NO CALLING TO PARI from inside the trap....
// Should replace by a stack, that would work.
static volatile long sage_pari_catcherr = 0;

/* Careful with pari_retries, it must be declared volatile!
 * Also note that "pari_errno = setjmp(...)" is not legal C.
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

#define _pari_endcatch {                                                      \
         err_leave(sage_pari_catcherr);                                       \
    }

