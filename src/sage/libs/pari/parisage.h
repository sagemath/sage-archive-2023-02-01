/* Include the relevant PARI includes in a safe way */
#include <pari/paricfg.h>
#include <pari/pari.h>

#undef coeff  /* Conflicts with NTL */
#undef ulong  /* Conflicts with FLINT */


/* Array element assignment */
#define set_gel(x, n, z)         (gel(x,n) = z)
#define set_gmael(x, i, j, z)    (gmael(x,i,j) = z)
#define set_gcoeff(x, i, j, z)   (gcoeff(x,i,j) = z)
