from sage.rings.memory import pmem_malloc

import sage.libs.ntl.all  as ntl

#import sage.libs.cf.cf as cf

pmem_malloc()


from sage.libs.pari.all   import pari, pari_gen, allocatemem, PariError

from sage.libs.mwrank.all  import (mwrank_EllipticCurve, mwrank_MordellWeil,
                                   mwrank_initprimes,
                                   get_precision as mwrank_get_precision,
                                   set_precision as mwrank_set_precision)


import symmetrica.all as symmetrica

from cremona.all import CremonaModularSymbols

