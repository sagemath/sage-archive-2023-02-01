#include <pari/pari.h>
#include "interrupt.h"


// global catch variable !
// this means that the code is not reentrant -- beware !
// THAT MEANS NO CALLING TO PARI from inside the trap....
// Should replace by a stack, that would work.
static void *__catcherr = NULL;

#define _pari_raise(errno) { \
        PyErr_SetObject(PyExc_PariError, PyInt_FromLong(errno)); \
    }

#define _pari_endcatch { err_leave(&__catcherr); }

/* Careful with pari_retries !
 * It should not be in a register, we flag it as "volatile".
 */
#define _pari_catch { \
        long pari_errno; \
        long volatile pari_retries = 0; \
        jmp_buf __env; \
        __catcherr = NULL; \
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

#define _pari_sig_on _sig_on; _pari_catch;
#define _pari_sig_str(s) _sig_str(s); _pari_catch;
#define _pari_sig_off _pari_endcatch; _sig_off;


int
gcmp_sage(GEN x, GEN y)
{
  long tx = typ(x), ty = typ(y), f;
  pari_sp av;

  if (is_intreal_t(tx))
    { if (is_intreal_t(ty)) return mpcmp(x,y); }
  else
  {
    if (tx==t_STR)
    {
      if (ty != t_STR) return 1;
      f = strcmp(GSTR(x),GSTR(y));
      return f > 0? 1
                  : f? -1: 0;
    }
    if (tx != t_FRAC)
    {
      if (ty == t_STR) return -1;
      return -1;
      /* pari_err(typeer,"comparison"); */
    }
  }
  if (ty == t_STR) return -1;
  if (!is_intreal_t(ty) && ty != t_FRAC)
    return 1; /* pari_err(typeer,"comparison"); */
  av=avma; y=gneg_i(y); f=gsigne(gadd(x,y)); avma=av; return f;
}
