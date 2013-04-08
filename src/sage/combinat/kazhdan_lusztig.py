"""
Kazhdan-Lusztig Polynomials
"""
#*****************************************************************************
#       Copyright (C) 2008 Daniel Bump <bump at match.stanford.edu>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import is_Polynomial
from sage.functions.other import floor
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_mpair
from sage.structure.sage_object import SageObject

class KazhdanLusztigPolynomial(SageObject):
    def __init__(self, W, q, trace=False):
        """
        A class for Kazhdan-Lusztig polynomials

        INPUT:

            - ``W`` - a Weyl Group.
            - ``q`` - an indeterminate.

        OPTIONAL:

            - trace

        The parent of q may be a PolynomialRing or a LaurentPolynomialRing.

        Set trace=True in order to see intermediate results. This is fun
        and instructive.

        EXAMPLES ::

            sage: W = WeylGroup("B3",prefix="s")
            sage: [s1,s2,s3]=W.simple_reflections()
            sage: R.<q>=LaurentPolynomialRing(QQ)
            sage: KL=KazhdanLusztigPolynomial(W,q)
            sage: KL.P(s2,s3*s2*s3*s1*s2)
            q + 1

        A faster implementation (using the optional package Coxeter 3) is given by::

            sage: W = CoxeterGroup(['B', 3], implementation='coxeter3') # optional - coxeter3
            sage: W.kazhdan_lusztig_polynomial([2], [3,2,3,1,2])        # optional - coxeter3
            q + 1
        """
        self._coxeter_group = W
        self._q = q
        self._trace = trace
        self._one = W.one()
        self._base_ring = q.parent()
        if is_Polynomial(q):
            self._base_ring_type = "polynomial"
        elif isinstance(q, LaurentPolynomial_mpair):
            self._base_ring_type = "laurent"
        else:
            self._base_ring_type = "unknown"

    @cached_method
    def R(self, x, y):
        """
        Returns the Kazhdan-Lusztig R polynomial.

        INPUT:

        - ``x``, ``y`` -- elements of the underlying Coxeter group

        EXAMPLES ::

           sage: R.<q>=QQ[]
           sage: W = WeylGroup("A2", prefix="s")
           sage: [s1,s2]=W.simple_reflections()
           sage: KL = KazhdanLusztigPolynomial(W, q)
           sage: [KL.R(x,s2*s1) for x in [1,s1,s2,s1*s2]]
           [q^2 - 2*q + 1, q - 1, q - 1, 0]

        """
        if x == 1:
            x = self._one
        if y == 1:
            y = self._one
        if x == y:
            return self._base_ring.one()
        if not x.bruhat_le(y):
            return self._base_ring.zero()
        if y.length() == 0:
            if x.length() == 0:
                return self._base_ring.one()
            else:
                return self._base_ring.zero()
        s = self._coxeter_group.simple_reflection(y.first_descent(side="left"))
        if (s*x).length() < x.length():
            ret = self.R(s*x,s*y)
            if self._trace:
                print "  R(%s,%s)=%s"%(x, y, ret)
            return ret
        else:
            ret = (self._q-1)*self.R(s*x,y)+self._q*self.R(s*x,s*y)
            if self._trace:
                print "  R(%s,%s)=%s"%(x, y, ret)
            return ret

    @cached_method
    def P(self, x, y):
        """
        Returns the Kazhdan-Lusztig P polynomial.

        If the rank is large, this runs slowly at first but speeds up
        as you do repeated calculations due to the caching.

        INPUT:

        - ``x``, ``y`` -- elements of the underlying Coxeter group

        .. SEEALSO::

            :mod:`~sage.libs.coxeter3.coxeter_group.CoxeterGroup.kazhdan_lusztig_polynomial`
            for a faster implementation using Fokko Ducloux's Coxeter3 C++ library.

        EXAMPLES ::

            sage: R.<q>=QQ[]
            sage: W = WeylGroup("A3", prefix="s")
            sage: [s1,s2,s3]=W.simple_reflections()
            sage: KL = KazhdanLusztigPolynomial(W, q)
            sage: KL.P(s2,s2*s1*s3*s2)
            q + 1
        """
        if x == 1:
            x = self._one
        if y == 1:
            y = self._one
        if x == y:
            return self._base_ring.one()
        if not x.bruhat_le(y):
            return self._base_ring.zero()
        if y.length() == 0:
            if x.length() == 0:
                return self._base_ring.one()
            else:
                return self._base_ring.zero()
        p = sum(-self.R(x,t)*self.P(t,y) for t in self._coxeter_group.bruhat_interval(x,y) if t != x)
        tr = floor((y.length()-x.length()+1)/2)
        try:
            ret = p.truncate(tr)
        except StandardError:
            ret = laurent_polynomial_truncate(p, tr)
        if self._trace:
            print "    P(%s,%s)=%s"%(x, y, ret)
        return ret

def laurent_polynomial_truncate(p, n):
    """
    Truncates the Laurent polynomial p, returning only terms
    of degree <n, similar to the truncate method for polynomials.

    EXAMPLES ::

        sage: from sage.combinat.kazhdan_lusztig import laurent_polynomial_truncate
        sage: P.<q> = LaurentPolynomialRing(QQ)
        sage: laurent_polynomial_truncate((q+q^-1)^3+q^2*(q+q^-1)^4,3)
        6*q^2 + 3*q + 4 + 3*q^-1 + q^-2 + q^-3
    """
    pdict = p._dict()
    dict = {}
    for k in pdict:
        if k[0] < n:
            dict[k] = pdict[k]
    return p.parent()(dict)
