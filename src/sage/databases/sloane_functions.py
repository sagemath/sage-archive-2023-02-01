r"""
Functions that compute some of the sequences in Sloane's tables

EXAMPLES:
   Type sloane.[tab] to see a list of the sequences that are defined.
   sage: d = sloane.A000005; d
    The integer sequence tau(n), which is the number of divisors of n.
    sage: d(1)
    1
    sage: d(6)
    4
    sage: d(100)
    9

Type \code{d._eval??} to see how the function that computes an individual
term of the sequence is implemented.

The input must be a positive integer:
    sage: d(0)
    Traceback (most recent call last):
    ...
    ValueError: input n (=0) must be a positive integer
    sage: d(1/3)
    Traceback (most recent call last):
    ...
    TypeError: Unable to coerce rational (=1/3) to an Integer.

You can also change how a sequence prints:
    sage: d = sloane.A000005; d
    The integer sequence tau(n), which is the number of divisors of n.
    sage: d.rename('(..., tau(n), ...)')
    sage: d
    (..., tau(n), ...)
    sage: d.reset_name()
    sage: d
    The integer sequence tau(n), which is the number of divisors of n.
"""

########################################################################
#
# To add your own new sequence here, do the following:
#
# 1. Add a new class to Section II below, which you should
#    do by copying an existing class and modifying it.
#    Make sure to at least define _eval and _repr_.
#    NOTES:  (a) define the _eval method only, which you may
#                assume has as input a *positive* SAGE integer.
#            (b) define the list method if there is a faster
#                way to compute the terms of the sequence than
#                just calling _eval (which is the default definition
#                of list).
#            (c) *AVOID* using gp.method if possible!  Use pari(obj).method()
#            (d) In many cases the function that computes a given integer
#                sequence belongs elsewhere in SAGE.  Put it there and make
#                your class in this file just call it.
#            (e) _eval should always return a SAGE integer.
#
# 2. Add an instance of your class in Section III below.
#
# 3. Type "sage -br" to rebuild SAGE, then fire up the notebook and
#    try out your new sequence.  Click the text button to get a version
#    of your session that you then include as a docstring.
#
# 4. Send a patch using
#      sage: hg_sage.ci()
#      sage: hg_sage.send('patchname')
#    (Email it to sage-dev@groups.google.com or post it online.)
#
########################################################################

########################################################################
# I. Define the generic Sloane sequence class.
########################################################################

# just used for handy .load, .save, etc.
from sage.structure.sage_object import SageObject
from sage.misc.misc import srange

class SloaneSequence(SageObject):
    def _repr_(self):
        raise NotImplementedError

    def __getitem__(self, n):
        return self(n)

    def __call__(self, n):
        m = Integer(n)
        if m <= 0:
            raise ValueError, "input n (=%s) must be a positive integer"%n
        return self._eval(m)

    def _eval(self, n):
        # this is what you implement in the derived class
        # the input n is assumed to be a *SAGE* integer >= 1
        raise NotImplementedError

    def list(self, n):
        return [self._eval(i) for i in srange(Integer(1),n+1)]


########################################################################
# II. Actual implementations of Sloane sequences.
########################################################################

# You may have to import more here when defining new sequences
import sage.rings.arith as arith
from sage.rings.integer import Integer


class A000005(SloaneSequence):
    r"""
    The sequence $tau(n)$, which is the number of divisors of $n$.

    This sequence is also denoted $d(n)$ (also called $\tau(n)$ or
    $\sigma_0(n)$), the number of divisors of n.

    EXAMPLES:
        sage: d = sloane.A000005; d
        The integer sequence tau(n), which is the number of divisors of n.
        sage: d(1)
        1
        sage: d(6)
        4
        sage: d(51)
        4
        sage: d(100)
        9
        sage: d(0)
        Traceback (most recent call last):
        ...
        ValueError: input n (=0) must be a positive integer
        sage: d.list(10)
        [1, 2, 2, 3, 2, 4, 2, 4, 3, 4]

    AUTHOR:
        -- Jaap Spies (2006-12-10)
        -- William Stein (2007-01-08)
    """
    def _repr_(self):
        return "The integer sequence tau(n), which is the number of divisors of n."

    def _eval(self, n):
        return arith.number_of_divisors(n)

    def list(self, n):
        return [self(i) for i in range(1,n+1)]


#############################################################
# Create the sloane object, off of which all the sequence
# objects hang.
#############################################################

class Sloane(SageObject):
    pass
sloane = Sloane()
sloane.A000005 = A000005()

