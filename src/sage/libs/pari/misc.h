/*****************************************
   PARI misc macros and functions
 *****************************************/
#ifndef SAGE_LIBS_PARI_MISC_H
#define SAGE_LIBS_PARI_MISC_H

#include <pari/pari.h>


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
