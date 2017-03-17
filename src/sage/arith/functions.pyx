"""
Fast Arithmetic Functions
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

include "cysignals/signals.pxi"
from sage.libs.gmp.mpz cimport mpz_lcm
from sage.rings.integer cimport Integer
from sage.structure.element cimport coercion_model

cpdef LCM_list(v):
    """
    Return the LCM of an interable ``v``.

    Elements of ``v`` are converted to Sage objects if they aren't already.

    This function is used, e.g., by :func:`~sage.arith.misc.lcm`.

    INPUT:

    -  ``v`` -- an iterable

    OUTPUT: integer

    EXAMPLES::

        sage: from sage.arith.functions import LCM_list
        sage: w = LCM_list([3,9,30]); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    The inputs are converted to Sage integers.

    ::

        sage: w = LCM_list([int(3), int(9), '30']); w
        90
        sage: type(w)
        <type 'sage.rings.integer.Integer'>

    TESTS::

        sage: from sage.structure.sequence import Sequence
        sage: from sage.arith.functions import LCM_list
        sage: l = Sequence(())
        sage: LCM_list(l)
        1

    This is because ``lcm(0,x) = 0`` for all ``x`` (by convention)::

        sage: LCM_list(Sequence(srange(100)))
        0

    So for the lcm of all integers up to 10 you must do this::

        sage: LCM_list(Sequence(srange(1,100)))
        69720375229712477164533808935312303556800

    Note that the following example did not work in QQ[] as of 2.11,
    but does in 3.1.4; the answer is different, though equivalent::

        sage: R.<X> = ZZ[]
        sage: LCM_list(Sequence((2*X+4,2*X^2,2)))
        2*X^3 + 4*X^2
        sage: R.<X> = QQ[]
        sage: LCM_list(Sequence((2*X+4,2*X^2,2)))
        X^3 + 2*X^2
    """
    cdef Integer x, z = Integer(1)

    itr = iter(v)
    for elt in itr:
        sig_check()
        if isinstance(elt, Integer):
            x = (<Integer>elt)
        elif isinstance(elt, (int, long, str)):
            x = Integer(elt)
        else:
            a,b = coercion_model.canonical_coercion(z, elt)
            return LCM_generic(itr, a.lcm(b))
        mpz_lcm(z.value, z.value, x.value)

    return z

cdef LCM_generic(itr, ret):
    for vi in itr:
        sig_check()
        a,b = coercion_model.canonical_coercion(vi, ret)
        ret = a.lcm(b)
        if not ret:
            return ret
    return ret

