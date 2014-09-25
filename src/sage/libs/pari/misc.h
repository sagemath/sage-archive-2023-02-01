/*****************************************
   PARI misc macros and functions
 *****************************************/
#ifndef SAGE_LIBS_PARI_MISC_H
#define SAGE_LIBS_PARI_MISC_H

#include <pari/pari.h>
#include <string.h>


/*
  Like gcmp(x, y), but returns 2 if comparison fails.
  Adapted from the PARI function gen2.c:gequal_try().
*/
inline int gcmp_try(GEN x, GEN y) {
    int i;
    pari_CATCH(CATCH_ALL) {
        GEN E = pari_err_last();
        switch(err_get_num(E))
        {
            case e_STACK: case e_MEM: case e_ALARM:
            pari_err(0, E); /* rethrow */
        }
        i = 2;
    } pari_TRY {
        i = gcmp(x, y);
    } pari_ENDCATCH;
    return i;
}

/* Compare the string representations of x and y.  */
inline int gcmp_string(GEN x, GEN y) {
    sig_block();
    char *a = GENtostr(x);
    char *b = GENtostr(y);
    int i = strcmp(a, b);
    free(a);
    free(b);
    sig_unblock();
    if (i == 0) return 0;
    return i > 0 ? 1 : -1;
}


int factorint_withproof_sage(GEN* ans, GEN x, GEN cutoff) {
  /*
  Factors and proves that the prime factors are really prime.
  If any aren't an ERROR condition (signal) is raised.

  INPUT:
     x -- a t_INT
     cutoff -- only check for primality of numbers at least this large.
  */

  GEN F = factorint(x, 0);
  *ans = F;

  long i, l;
  if (lg(F) == 1) return F; // x = 1
  F = gel(F,1); l = lg(F);
  for (i = 1; i < l; i++) {
    GEN p = gel(F,i);
    if (mpcmp(p, cutoff) > 0 && !isprime(p)) {
      char *c, *d;
      c = GENtostr(x);
      d = GENtostr(p);
      fprintf(stderr, "***\nPARI's factor(%s): Found composite pseudoprime %s (very rare and exciting -- PLEASE REPORT!!)\n***\n",
                   c, d);
      fprintf(stderr, "Do not worry, SAGE will further factor the number until each factor is proven prime.\n");
      free(c);
      free(d);
      return 1;
    }
  }
  return 0;
}

#endif  /* SAGE_LIBS_PARI_MISC_H */
