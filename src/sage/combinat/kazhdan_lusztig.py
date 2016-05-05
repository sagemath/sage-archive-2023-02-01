r"""
Kazhdan-Lusztig Polynomials

AUTHORS:

- Daniel Bump (2008): initial version

- Alan J.X. Guo (2014-03-18): ``R_tilde()`` method.

"""
#*****************************************************************************
#       Copyright (C) 2008 Daniel Bump <bump at match.stanford.edu>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.polynomial.polynomial_element import is_Polynomial
from sage.functions.other import floor
from sage.misc.cachefunc import cached_method
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_generic
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.combinat.root_system.coxeter_group import CoxeterGroup

class KazhdanLusztigPolynomial(UniqueRepresentation, SageObject):
    """
    A Kazhdan-Lusztig polynomial.

    INPUT:

    - ``W`` -- a Weyl Group
    - ``q`` -- an indeterminate

    OPTIONAL:

    - ``trace`` -- if ``True``, then this displays the trace: the intermediate
      results. This is instructive and fun.

    The parent of ``q`` may be a :class:`PolynomialRing` or a
    :class:`LaurentPolynomialRing`.

    REFERENCES:

    .. [KL79] \D. Kazhdan and G. Lusztig. *Representations of Coxeter
       groups and Hecke algebras*. Invent. Math. **53** (1979).
       no. 2, 165--184. :doi:`10.1007/BF01390031` :mathscinet:`MR0560412`

    .. [Dy93] \M. J. Dyer. *Hecke algebras and shellings of Bruhat
       intervals*. Compositio Mathematica, 1993, 89(1): 91-115.

    .. [BB05] \A. Bjorner, F. Brenti. *Combinatorics of Coxeter
       groups*. New York: Springer, 2005.

    EXAMPLES::

        sage: W = WeylGroup("B3",prefix="s")
        sage: [s1,s2,s3] = W.simple_reflections()
        sage: R.<q> = LaurentPolynomialRing(QQ)
        sage: KL = KazhdanLusztigPolynomial(W,q)
        sage: KL.P(s2,s3*s2*s3*s1*s2)
        1 + q

    A faster implementation (using the optional package Coxeter 3) is given by::

        sage: W = CoxeterGroup(['B', 3], implementation='coxeter3') # optional - coxeter3
        sage: W.kazhdan_lusztig_polynomial([2], [3,2,3,1,2])        # optional - coxeter3
        q + 1
    """
    def __init__(self, W, q, trace=False):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: W = WeylGroup("B3",prefix="s")
            sage: R.<q> = LaurentPolynomialRing(QQ)
            sage: KL = KazhdanLusztigPolynomial(W,q)
            sage: TestSuite(KL).run()
        """
        self._coxeter_group = W
        self._q = q
        self._trace = trace
        self._one = W.one()
        self._base_ring = q.parent()
        if is_Polynomial(q):
            self._base_ring_type = "polynomial"
        elif isinstance(q, LaurentPolynomial_generic):
            self._base_ring_type = "laurent"
        else:
            self._base_ring_type = "unknown"

    @cached_method
    def R(self, x, y):
        """
        Return the Kazhdan-Lusztig `R` polynomial.

        INPUT:

        - ``x``, ``y`` -- elements of the underlying Coxeter group

        EXAMPLES::

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
    def R_tilde(self, x, y):
        r"""
        Return the Kazhdan-Lusztig `\tilde{R}` polynomial.

        Information about the `\tilde{R}` polynomials can be found in
        [Dy93]_ and [BB05]_.

        INPUT:

        - ``x``, ``y`` -- elements of the underlying Coxeter group

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: W = WeylGroup("A2", prefix="s")
            sage: [s1,s2] = W.simple_reflections()
            sage: KL = KazhdanLusztigPolynomial(W, q)
            sage: [KL.R_tilde(x,s2*s1) for x in [1,s1,s2,s1*s2]]
            [q^2, q, q, 0]
        """
        if x == 1:
            x = self._one
        if y == 1:
            y = self._one
        if not x.bruhat_le(y):
            return self._base_ring.zero()
        if x == y:
            return self._base_ring.one()
        s = self._coxeter_group.simple_reflection(y.first_descent(side="right"))
        if (x * s).length() < x.length():
            ret = self.R_tilde(x * s, y * s)
            if self._trace:
                print " R_tilde(%s,%s)=%s" % (x, y, ret)
            return ret
        else:
            ret = self.R_tilde(x * s, y * s) + self._q * self.R_tilde(x, y * s)
            if self._trace:
                print " R_tilde(%s,%s)=%s" % (x, y, ret)
            return ret

    @cached_method
    def P(self, x, y):
        """
        Return the Kazhdan-Lusztig `P` polynomial.

        If the rank is large, this runs slowly at first but speeds up
        as you do repeated calculations due to the caching.

        INPUT:

        - ``x``, ``y`` -- elements of the underlying Coxeter group

        .. SEEALSO::

            :mod:`~sage.libs.coxeter3.coxeter_group.CoxeterGroup.kazhdan_lusztig_polynomial`
            for a faster implementation using Fokko Ducloux's Coxeter3 C++ library.

        EXAMPLES::

            sage: R.<q> = QQ[]
            sage: W = WeylGroup("A3", prefix="s")
            sage: [s1,s2,s3] = W.simple_reflections()
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
        ret = p.truncate(tr)
        if self._trace:
            print "    P({},{})={}".format(x, y, ret)
        return ret

