"""
Andrew Sutherland's Probabilistic Image of Galois Algorithm

AUTHOR:

   - William Stein, 2010-03 -- wrote the Cython wrapper
   - Sutherland -- wrote the C code and designed the algorithm
"""

#############################################################################
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


include "../../ext/stdsage.pxi"
include "../../ext/cdefs.pxi"
from sage.rings.integer cimport Integer

cdef extern from "galrep.h":
    int galrep_ecdata_load (char *filename)
    int galrep_gl2data_load (char *filename)
    int galrep_ec_modl_image (int ell, mpz_t A, mpz_t B, int errbnd)
    int galrep_ec_modl_images (int *ccs, int min, int max, mpz_t A, mpz_t B, int errbnd)
    int galrep_ec_non_surjective (int min, int max, mpz_t A, mpz_t B, int errbnd)
    int GALREP_MAX_ELL


# There is no need to support pickling for this class, since it is
# only used for providing an API for some functions.

cdef class GalRep:
    """
    Andrew Sutherland's Probabilistic Image of Galois Algorithm
    """
    cdef object primes

    def __init__(self, galrep_ecdata=None, galrep_gl2data=None):
        """
        INPUT:

            - galrep_ecdata -- (default: None); if given, specifies galrep_ecdata.dat filename

            - galrep_gl2data -- (default: None); if given, specifies galrep_gl2data.dat filename

        EXAMPLES::

            sage: galrep.GalRep()
            Andrew Sutherland's Probabilistic Image of Galois Algorithm
            sage: galrep.GalRep(os.path.join(sage.misc.misc.SAGE_DATA,'galrep','galrep_ecdata.dat'))
            Andrew Sutherland's Probabilistic Image of Galois Algorithm
        """
        from sage.misc.misc import SAGE_DATA
        import os
        if galrep_ecdata is None:
            galrep_ecdata = os.path.join(SAGE_DATA, 'galrep', 'galrep_ecdata.dat')
        if galrep_gl2data is None:
            galrep_gl2data = os.path.join(SAGE_DATA, 'galrep', 'galrep_gl2data.dat')
        galrep_ecdata_load(galrep_ecdata)
        galrep_gl2data_load(galrep_gl2data)
        self.primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

    def mod_ell_image(self, Integer A, Integer B, int ell):
        """
        Returns the conjugacy class id of the image of the mod-ell
        Galois representation attached to the elliptic curve
        y^2=x^3+A*x+B.  The id code 0 signifies a surjective
        representation.

        INPUT:

            - A -- integer
            - B -- integer
            - ell -- prime integer (must be <= 59)

        OUTPUT:

            - integer

        EXAMPLES::

            sage: galrep.GalRep().mod_ell_image(-432,8208,3)
            0
            sage: galrep.GalRep().mod_ell_image(-432,8208,2)
            0
            sage: galrep.GalRep().mod_ell_image(-432,8208,5)
            9
            sage: galrep.GalRep().mod_ell_image(-432,8208,0)
            Traceback (most recent call last):
            ...
            ValueError: ell (=0) must be a prime <= 59
            sage: galrep.GalRep().mod_ell_image(-432,8208,67)
            Traceback (most recent call last):
            ...
            ValueError: ell (=67) must be a prime <= 59
        """
        cdef int a = galrep_ec_modl_image(ell, A.value, B.value, 0)
        if a == -1: raise ValueError, "ell (=%s) must be a prime <= %s"%(ell, GALREP_MAX_ELL)
        return a

    def mod_ell_images(self, Integer A, Integer B, int min=2, int max=59):
        """
        Returns list of the conjugacy class id's of the images of the
        mod-ell Galois representation attached to the elliptic curve
        y^2=x^3+A*x+B for ell in the interval [min, max].  The id code
        0 signifies a surjective representation.

        INPUT:

            - A -- integer
            - B -- integer
            - min -- integer (default: 2)
            - max -- integer (default: 59)  (must be <= 59)

        OUTPUT:

            - list of primes

        EXAMPLES::

            sage: galrep.GalRep().mod_ell_images(-432,8208)
            [0, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: galrep.GalRep().mod_ell_images(-432,8208,2,10)
            [0, 0, 9, 0]
        """
        if min >= max:
            raise ValueError, "max must be bigger than min"
        if max > GALREP_MAX_ELL:
            raise ValueError, "max must be <= %s"%GALREP_MAX_ELL
        cdef int* ccs = <int*> sage_malloc(sizeof(int)*len(self.primes))
        cdef int i, a = galrep_ec_modl_images (ccs, min, max, A.value, B.value, 0)
        cdef list v = []
        for i in range(len(self.primes)):
            if self.primes[i] > max: break
            v.append(ccs[i])
        sage_free(ccs)
        return v

    def non_surjective_primes(self, Integer A, Integer B, int min=2, int max=59):
        """
        Returns a list of the primes ell in the interval [min, max]
        (default=[2,59]) such that the mod-ell Galois representation
        attached to y^2=x^3+A*x+B is highly likely to be
        non-surjective.

        INPUT:

            - A -- integer
            - B -- integer
            - min -- integer (default: 2)
            - max -- integer (default: 59)  (must be <= 59)

        OUTPUT:

            - list of primes

        EXAMPLES::

            sage: galrep.GalRep().non_surjective_primes(-432,8208)
            [5]
            sage: galrep.GalRep().non_surjective_primes(-432,8208,7,59)
            []
            sage: galrep.GalRep().non_surjective_primes(-432,8208,1,10)
            [5]
        """
        if min >= max:
            raise ValueError, "max must be bigger than min"
        cdef int a = galrep_ec_non_surjective (min, max, A.value, B.value, 0)
        if a == -1: raise ValueError, "min and max must be <= %s"%GALREP_MAX_ELL
        cdef int i, j=1
        v = []
        for i in range(len(self.primes)):
            if a & j:
                v.append(self.primes[i])
            j = j << 1
        return v

    def __repr__(self):
        """
        EXAMPLES::

            sage: galrep.GalRep().__repr__()
            "Andrew Sutherland's Probabilistic Image of Galois Algorithm"
        """
        return "Andrew Sutherland's Probabilistic Image of Galois Algorithm"


