"""
Singular C function and class declaration

AUTHOR: Martin Albrecht <malb@informatik.uni-bremen.de>
"""

################################################################################
#
################################################################################

###############################################################################
#   SAGE: System for Algebra and Geometry Experimentation
#       Copyright (C) 2005, 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################

include "singular-cdefs.pxi"

from sage.rings.rational_field import RationalField
from sage.rings.finite_field import FiniteField_prime_modn


cdef extern from "stdsage.h":
    ctypedef void PyObject
    object PY_NEW(object t)
    int PY_TYPE_CHECK(object o, object t)
    PyObject** FAST_SEQ_UNSAFE(object o)
    void init_csage()

cdef class Conversion:
    """
    A work around class to export the contained methods/functions
    """

    cdef public Rational si2sa_QQ(self, number *n, ring *_ring):
        """
        Converts a SINGULAR rational number to a SAGE rational number.

        INPUT:
            n -- number
            _ring -- singular ring, used to check type of n

        OUTPUT:
            SAGE rational number matching n
        """

        #TYPECHECK HERE

        cdef number *nom
        cdef number *denom
        cdef mpq_t _z

        cdef mpz_t nom_z, denom_z

        cdef Rational z

        mpq_init(_z)

        ##  Immediate integers handles carry the tag 'SR_INT', i.e. the last bit is 1.
        ##  This distuingishes immediate integers from other handles which  point  to
        ##  structures aligned on 4 byte boundaries and therefor have last bit  zero.
        ##  (The second bit is reserved as tag to allow extensions of  this  scheme.)
        ##  Using immediates as pointers and dereferencing them gives address errors.

        nom = nlGetNom(n, _ring)
        mpz_init(nom_z)

        if (SR_HDL(nom) & SR_INT): mpz_set_si(nom_z, SR_TO_INT(nom))
        else: mpz_set(nom_z,&nom.z)

        mpq_set_num(_z,nom_z)
        n_Delete(&nom,_ring)
        mpz_clear(nom_z)

        denom = nlGetDenom(n, _ring)
        mpz_init(denom_z)

        if (SR_HDL(denom) & SR_INT): mpz_set_si(denom_z, SR_TO_INT(denom))
        else: mpz_set(denom_z,&denom.z)

        mpq_set_den(_z, denom_z)
        n_Delete(&denom,_ring)
        mpz_clear(denom_z)

        z = Rational()
        z.set_from_mpq(_z)
        return z

    cdef public number *sa2si_QQ(self, Rational r, ring *_ring):
        """
        """
        cdef number *z
        return nlInit2gmp( mpq_numref(r.value), mpq_denref(r.value) )

    cdef public object si2sa(self, number *n, ring *_ring, object base):
        if PY_TYPE_CHECK(base, FiniteField_prime_modn):
            return base(nInt(n))

        elif PY_TYPE_CHECK(base, RationalField):
            return self.si2sa_QQ(n,_ring)
        else:
            raise ValueError, "cannot convert SINGULAR number"

    cdef public number *sa2si(self, Element elem, ring * _ring):
        if PY_TYPE_CHECK(elem._parent, FiniteField_prime_modn):
            return n_Init(int(elem),_ring)

        elif PY_TYPE_CHECK(elem._parent, RationalField):
            return self.sa2si_QQ(elem, _ring)
        else:
            raise ValueError, "cannot convert SINGULAR number"
