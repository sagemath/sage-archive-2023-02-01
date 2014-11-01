/* Include the relevant PARI includes in a safe way */
#include <pari/paricfg.h>
#include <pari/pari.h>

#undef coeff  /* Conflicts with NTL */
