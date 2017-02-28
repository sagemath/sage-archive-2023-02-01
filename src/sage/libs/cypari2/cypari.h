/*
 * Additional macros and fixes for the PARI headers. This is meant to
 * be included after including <pari/pari.h>
 */

#undef coeff  /* Conflicts with NTL (which is used by SageMath) */


/* Array element assignment */
#define set_gel(x, n, z)         (gel(x,n) = z)
#define set_gmael(x, i, j, z)    (gmael(x,i,j) = z)
#define set_gcoeff(x, i, j, z)   (gcoeff(x,i,j) = z)
#define set_uel(x, n, z)         (uel(x,n) = z)
