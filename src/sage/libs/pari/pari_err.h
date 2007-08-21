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

inline int strcmp_to_cmp(int f) {
    if (f > 0) {
      return 1;
    } else if (f) {
      return -1;
    } else {
      return 0;
    }
}

inline int
gcmp_sage(GEN x, GEN y)
{
  long tx = typ(x), ty = typ(y), f;
  GEN tmp;
  pari_sp av;

  if (is_intreal_t(tx) && is_intreal_t(ty)) {
    /* Compare two numbers that can be considered as reals. */
    return mpcmp(x,y);
  }

  /***** comparing strings *****/
  if (tx==t_STR) {
    /* Compare two strings */
    if (ty != t_STR)  {
      return 1;
    }

    return strcmp_to_cmp(strcmp(GSTR(x),GSTR(y)));
  }
  if (ty == t_STR) /* tx is not a string */
     return -1;
  /***** end comparing strings *****/

  /*if (!is_intreal_t(ty) && ty != t_FRAC)  */
  /*     return 1; */
 /* pari_err(typeer,"comparison"); */

  av = avma;
  char *c, *d;
  c = GENtostr(x);
  d = GENtostr(y);
  f = strcmp_to_cmp(strcmp(c, d));
  free(c);
  free(d);
  avma = av;
  return f;

  /*
  av = avma;
  y = gneg_i(y);
  tmp = gadd(x,y);
  switch(typ(tmp)) {
     case t_INT:
     case t_REAL:
        return signe(tmp);
     case t_FRAC:
        return signe(tmp[1]);
  }
  avma = av;
  */
}
