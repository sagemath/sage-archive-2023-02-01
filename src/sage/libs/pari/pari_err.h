#include <pari/pari.h>
#include "interrupt.h"


// global catch variable !
// this means that the code is not reentrant -- beware !
static void *__catcherr = NULL;

#define _pari_raise(errno) { \
        PyErr_SetObject(PyExc_RuntimeError, PyInt_FromLong(errno)); \
    }

#define _pari_endcatch { err_leave(&__catcherr); }

#define _pari_catch { \
        long pari_errno; \
        jmp_buf __env; \
        __catcherr = NULL; \
        if ((pari_errno=setjmp(__env))) { \
            _pari_endcatch; \
            _pari_raise(pari_errno); \
            return NULL; \
        } \
        __catcherr = err_catch(CATCH_ALL, &__env); \
    }

#define _pari_sig_on _sig_on; _pari_catch;
#define _pari_sig_str(s) _sig_str(s); _pari_catch;
#define _pari_sig_off _pari_endcatch; _sig_off;


