# -*- coding: utf-8 -*-
#
#  NB for doctests:   The first time that the c++ library giac is loaded a message appears.
#                     This message is version and arch dependant.
"""
Wrappers for Giac functions

<Paragraph description>

AUTHORS:

- Martin Albrecht (2015-07-01): initial version
- Han Frederic (2015-07-01): initial version

EXAMPLES::

    sage: from sage.libs.giac import groebner_basis_libgiac # optional - giacpy
    sage: P = PolynomialRing(QQ, 6, 'x')
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: B = groebner_basis_libgiac(I.gens()) # optional - giacpy, random
    sage: B # optional - giacpy
    Polynomial Sequence with 45 Polynomials in 6 Variables
"""

#*****************************************************************************
#       Copyright (C) 2013 Frederic Han <frederic.han@imj-prg.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.proof.all import polynomial as proof_polynomial
from sage.rings.polynomial.multi_polynomial_sequence import PolynomialSequence

class GiacSettingsDefaultContext:
    """
    Context preserve libgiac settings.
    """

    def __enter__(self):
        """
        EXAMPLE::

            sage: from sage.libs.giac import GiacSettingsDefaultContext  # optional - giacpy
            sage: from giacpy import giacsettings # optional - giacpy
            sage: giacsettings.proba_epsilon = 1/4 # optional - giacpy, random
            sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1/2 # optional - giacpy, random
            sage: giacsettings.proba_epsilon # optional - giacpy
            0.25

        """
        try:
            from giacpy import giacsettings, libgiac
        except ImportError:
            raise ImportError("""One of the optional packages giac or giacpy is missing""")

        self.proba_epsilon = giacsettings.proba_epsilon
        self.debuginfolevel = libgiac('debug_infolevel()')

    def __exit__(self, typ, value, tb):
        """
        EXAMPLE::

            sage: from sage.libs.giac import GiacSettingsDefaultContext  # optional - giacpy
            sage: from giacpy import giacsettings # optional - giacpy
            sage: giacsettings.proba_epsilon = 1/4 # optional - giacpy, random
            sage: with GiacSettingsDefaultContext(): giacsettings.proba_epsilon = 1/8 # optional - giacpy, random
            sage: giacsettings.proba_epsilon # optional - giacpy
            0.25

        """
        try:
            from giacpy import giacsettings, libgiac
        except ImportError:
            raise ImportError("""One of the optional packages giac or giacpy is missing""")

        libgiac('debug_infolevel')(self.debuginfolevel)
        giacsettings.proba_epsilon = self.proba_epsilon

def local_giacsettings(func):
    """
    Decorator to preserve Giac's epsilon settings.

    EXAMPLE::

    """
    from sage.misc.decorators import sage_wraps

    @sage_wraps(func)
    def wrapper(*args, **kwds):
        """
        Execute function in ``GiacSettingsDefaultContext``.
        """
        with GiacSettingsDefaultContext():
            return func(*args, **kwds)

    return wrapper


@local_giacsettings
def groebner_basis_libgiac(gens, epsilon=None, prot=False, *args, **kwds):
    """
    """
    try:
        from giacpy import libgiac, giacsettings
    except ImportError:
        raise ImportError("""One of the optional packages giac or giacpy is missing""")

    try:
        iter(gens)
    except TypeError:
        gens = gens.gens()

    # get the ring from gens
    P = iter(gens).next().parent()
    K = P.base_ring()
    p = K.characteristic()

    if K.is_prime_field() and p == 0:
        F = libgiac(gens)
    elif K.is_prime_field() and p < 2**31:
        F = (libgiac(gens) % p)
    else:
        raise NotImplementedError("Only prime fields of cardinal < 2^31 are implemented in Giac for Groebner bases.")

    if P.term_order() != "degrevlex":
        raise NotImplementedError("Only degrevlex term orderings are supported in Giac Groebner bases.")

    # proof or probabilistic reconstruction
    if epsilon is None:
        if proof_polynomial():
            giacsettings.proba_epsilon = 0
        else:
            giacsettings.proba_epsilon = 1e-15
    else:
        giacsettings.proba_epsilon = epsilon

    # prot
    if (prot==True):
        libgiac('debug_infolevel(2)')

    # compute de groebner basis with giac
    gb_giac = F.gbasis([P.gens()], "revlex")

    return PolynomialSequence(gb_giac, P, immutable=True)
