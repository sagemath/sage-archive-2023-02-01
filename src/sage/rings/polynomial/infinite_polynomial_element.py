"""
Elements of Infinite Polynomial Rings

AUTHORS:

- Simon King <simon.king@nuigalway.ie>
- Mike Hansen <mhansen@gmail.com>

An Infinite Polynomial Ring has generators `x_\\ast, y_\\ast,...`, so
that the variables are of the form `x_0, x_1, x_2, ..., y_0, y_1,
y_2,...,...` (see :mod:`~sage.rings.polynomial.infinite_polynomial_ring`).
Using the generators, we can create elements as follows::

    sage: X.<x,y> = InfinitePolynomialRing(QQ)
    sage: a = x[3]
    sage: b = y[4]
    sage: a
    x_3
    sage: b
    y_4
    sage: c = a*b+a^3-2*b^4
    sage: c
    x_3^3 + x_3*y_4 - 2*y_4^4

Any Infinite Polynomial Ring ``X`` is equipped with a monomial ordering.
We only consider monomial orderings in which:

    ``X.gen(i)[m] > X.gen(j)[n]`` `\iff` ``i<j``, or ``i==j`` and ``m>n``

Under this restriction, the monomial ordering can be lexicographic
(default), degree lexicographic, or degree reverse lexicographic.
Here, the ordering is lexicographic, and elements can be compared
as usual::

    sage: X._order
    'lex'
    sage: a > b
    True

Note that, when a method is called that is not directly implemented
for 'InfinitePolynomial', it is tried to call this method for the
underlying *classical* polynomial. This holds, e.g., when applying the
``latex`` function::

    sage: latex(c)
    x_{3}^{3} + x_{3} y_{4} - 2 y_{4}^{4}

There is a permutation action on Infinite Polynomial Rings by
permuting the indices of the variables::

    sage: P = Permutation(((4,5),(2,3)))
    sage: c^P
    x_2^3 + x_2*y_5 - 2*y_5^4

Note that ``P(0)==0``, and thus variables of index zero are invariant
under the permutation action.  More generally, if ``P`` is any
callable object that accepts non-negative integers as input and
returns non-negative integers, then ``c^P`` means to apply ``P`` to
the variable indices occurring in ``c``.

TESTS:

We test whether coercion works, even in complicated cases in which
finite polynomial rings are merged with infinite polynomial rings::

    sage: A.<a> = InfinitePolynomialRing(ZZ,implementation='sparse',order='degrevlex')
    sage: B.<b_2,b_1> = A[]
    sage: C.<b,c> = InfinitePolynomialRing(B,order='degrevlex')
    sage: C
    Infinite polynomial ring in b, c over Infinite polynomial ring in a over Integer Ring
    sage: 1/2*b_1*a[4]+c[3]
    1/2*a_4*b_1 + c_3

"""

#*****************************************************************************
#       Copyright (C) 2009 Simon King <king@mathematik.nuigalway.ie>
#                          and Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer_ring import ZZ
from sage.rings.integer import Integer
from sage.structure.element import RingElement
from sage.misc.cachefunc import cached_method
import copy

def InfinitePolynomial(A, p):
    """
    Create an element of a Polynomial Ring with a Countably Infinite Number of Variables.

    Usually, an InfinitePolynomial is obtained by using the generators
    of an Infinite Polynomial Ring (see :mod:`~sage.rings.polynomial.infinite_polynomial_ring`)
    or by conversion.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring.
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``.

    ASSUMPTIONS:

    In the dense implementation, it must be ensured that the argument
    ``p`` coerces into ``A._P`` by a name preserving conversion map.

    In the sparse implementation, in the direct construction of an
    infinite polynomial, it is *not* tested whether the argument ``p``
    makes sense in ``A``.

    EXAMPLES::

        sage: from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        sage: X.<alpha> = InfinitePolynomialRing(ZZ)
        sage: P.<alpha_1,alpha_2> = ZZ[]

    Currently, ``P`` and ``X._P`` (the underlying polynomial ring of
    ``X``) both have two variables::

        sage: X._P
        Multivariate Polynomial Ring in alpha_1, alpha_0 over Integer Ring

    By default, a coercion from ``P`` to  ``X._P`` would not be name preserving.
    However, this is taken care for; a name preserving conversion is impossible,
    and by consequence an error is raised::

        sage: InfinitePolynomial(X, (alpha_1+alpha_2)^2)
        Traceback (most recent call last):
        ...
        TypeError: Could not find a mapping of the passed element to this ring.

    When extending the underlying polynomial ring, the construction of
    an infinite polynomial works::

        sage: alpha[2]
        alpha_2
        sage: InfinitePolynomial(X, (alpha_1+alpha_2)^2)
        alpha_2^2 + 2*alpha_2*alpha_1 + alpha_1^2

    In the sparse implementation, it is not checked whether the
    polynomial really belongs to the parent::

        sage: Y.<alpha,beta> = InfinitePolynomialRing(GF(2), implementation='sparse')
        sage: a = (alpha_1+alpha_2)^2
        sage: InfinitePolynomial(Y, a)
        alpha_1^2 + 2*alpha_1*alpha_2 + alpha_2^2

    However, it is checked when doing a conversion::

        sage: Y(a)
        alpha_2^2 + alpha_1^2

    """
    from sage.all import parent
    if hasattr(A,'_P'):
        if parent(p) is A._P or (A._P.base_ring().has_coerce_map_from(parent(p))):
            return InfinitePolynomial_dense(A, p)
        # MPolynomialRing_polydict is crab. So, in that case, use sage_eval
        from sage.rings.polynomial.multi_polynomial_ring import MPolynomialRing_polydict
        if isinstance(A._P, MPolynomialRing_polydict):
            from sage.rings.polynomial.infinite_polynomial_ring import GenDictWithBasering
            from sage.misc.sage_eval import sage_eval
            p = sage_eval(repr(p), GenDictWithBasering(A._P,A._P.gens_dict()))
            return InfinitePolynomial_dense(A, p)
        else:
            # Now there remains to fight the oddities and bugs of libsingular.
            PP = p.parent()
            if A._P.has_coerce_map_from(PP):
                if A._P.ngens() == PP.ngens(): # coercion is sometimes by position!
                    f = PP.hom(PP.variable_names(),A._P)
                    try:
                        return InfinitePolynomial_dense(A, f(p))
                    except (ValueError, TypeError):
                        # last desparate attempt: String conversion
                        from sage.misc.sage_eval import sage_eval
                        from sage.rings.polynomial.infinite_polynomial_ring import GenDictWithBasering
                        # the base ring may be a function field, therefore
                        # we need GenDictWithBasering
                        return InfinitePolynomial_dense(A, sage_eval(repr(p), GenDictWithBasering(A._P, A._P.gens_dict())))
                return InfinitePolynomial_dense(A, A._P(p))
            # there is no coercion, so, we set up a name-preserving map.
            SV = set([repr(x) for x in p.variables()])
            f = PP.hom([x if x in SV else 0 for x in PP.variable_names()], A._P)
            try:
                return InfinitePolynomial_dense(A, f(p))
            except (ValueError, TypeError):
                # last desparate attempt: String conversion
                from sage.misc.sage_eval import sage_eval
                from sage.rings.polynomial.infinite_polynomial_ring import GenDictWithBasering
                # the base ring may be a function field, therefore
                # we need GenDictWithBasering
                return InfinitePolynomial_dense(A, sage_eval(repr(p), GenDictWithBasering(A._P, A._P.gens_dict())))
    return InfinitePolynomial_sparse(A, p)

class InfinitePolynomial_sparse(RingElement):
    """
    Element of a sparse Polynomial Ring with a Countably Infinite Number of Variables.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring in sparse implementation
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``.

    Of course, one should not directly invoke this class, but rather
    construct elements of ``A`` in the usual way.

    EXAMPLES::

        sage: A.<a> = QQ[]
        sage: B.<b,c> = InfinitePolynomialRing(A,implementation='sparse')
        sage: p = a*b[100] + 1/2*c[4]
        sage: p
        a*b_100 + 1/2*c_4
        sage: p.parent()
        Infinite polynomial ring in b, c over Univariate Polynomial Ring in a over Rational Field
        sage: p.polynomial().parent()
        Multivariate Polynomial Ring in b_100, b_0, c_4, c_0 over Univariate Polynomial Ring in a over Rational Field

    """
    # Construction and other basic methods
    # We assume that p is good input. Type checking etc. is now done
    # in the _element_constructor_ of the parent.
    def __init__(self, A, p):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[1] + x[2]
            sage: a == loads(dumps(a))
            True

        """
        self._has_footprint = False
        self._footprint = {}
        self._p = p
        RingElement.__init__(self, A)

    def _repr_(self):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: str(x[1] + x[2])  # indirect doctest
            'x_2 + x_1'

        """
        return repr(self._p)

    def _hash_(self):
        """
        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1]
            sage: hash(a) # indirect doctest
            -6172640511012239345   # 64-bit
            -957478897             # 32-bit

        """
        return hash(self._p)

    def polynomial(self):
        """
        Return the underlying polynomial.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(GF(7))
            sage: p=x[2]*y[1]+3*y[0]
            sage: p
            x_2*y_1 + 3*y_0
            sage: p.polynomial()
            x_2*y_1 + 3*y_0
            sage: p.polynomial().parent()
            Multivariate Polynomial Ring in x_2, x_1, x_0, y_2, y_1, y_0 over Finite Field of size 7
            sage: p.parent()
            Infinite polynomial ring in x, y over Finite Field of size 7

        """
        return self._p

    def __call__(self, *args, **kwargs):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ,implementation='sparse')
            sage: a = x[0] + x[1]
            sage: a(x_0=2,x_1=x[1])
            x_1 + 2
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=3)
            x_0 + 3
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=x[100])
            x_100 + x_0

        """
        #Replace any InfinitePolynomials by their underlying polynomials
        if hasattr(self._p,'variables'):
            V = [str(x) for x in self._p.variables()]
        else:
            V = []
        for kw in kwargs:
            value = kwargs[kw]
            if isinstance(value, InfinitePolynomial_sparse):
                kwargs[kw] = value._p
                V.append(kw)
                if hasattr(value._p,'variables'):
                    V.extend([str(x) for x in value._p.variables()])
        args = list(args)
        for i, arg in enumerate(args):
            if isinstance(arg, InfinitePolynomial_sparse):
                args[i] = arg._p
                if hasattr(arg._p,'variables'):
                    V.extend([str(x) for x in arg._p.variables()])
        V=list(set(V))
        V.sort(cmp=self.parent().varname_cmp,reverse=True)
        if V:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(),V,order=self.parent()._order)
        else:
            return self
        res = R(self._p)(*args, **kwargs)
        try:
            from sage.misc.sage_eval import sage_eval
            return sage_eval(repr(res), self.parent()._gens_dict)
        except Exception:
            return res

    def _getAttributeNames(self):
        """
        This method implements tab completion, see :trac:`6854`.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: import sagenb.misc.support as s
            sage: p = x[3]*x[2]
            sage: s.completions('p.co',globals(),system='python') # indirect doctest
            ['p.coefficient', 'p.coefficients', 'p.constant_coefficient', 'p.content']

        """
        return dir(self._p)

    def __dir__(self):
        """
        This method implements tab completion, see :trac:`6854`.

        TESTS::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: import sagenb.misc.support as s
            sage: p = x[3]*x[2]
            sage: s.completions('p.co',globals(),system='python') # indirect doc test
            ['p.coefficient', 'p.coefficients', 'p.constant_coefficient', 'p.content']
            sage: 'constant_coefficient' in dir(p) # indirect doctest
            True
        """
        return dir(self._p)

    def __getattr__(self, s):
        """
        NOTE:

        This method will only be called if an attribute of ``self``
        is requested that is not known to Python. In that case,
        the corresponding attribute of the underlying polynomial
        of ``self`` is returned.

        EXAMPLES:

        Elements of Infinite Polynomial Rings have no genuine
        ``_latex_`` method. But the method inherited from the
        underlying polynomial suffices::

            sage: X.<alpha> = InfinitePolynomialRing(QQ)
            sage: latex(alpha[3]*alpha[2]^2) # indirect doctest
            \alpha_{3} \alpha_{2}^{2}

        Related with tickets :trac:`6854` and :trac:`7580`, the attribute
        ``__methods__`` is treated in a special way, which
        makes introspection and tab completion work::

            sage: import sagenb.misc.support as s
            sage: p = alpha[3]*alpha[2]^2
            sage: s.completions('p.co',globals(),system='python') # indirect doc test
            ['p.coefficient', 'p.coefficients', 'p.constant_coefficient', 'p.content']
            sage: 'constant_coefficient' in dir(p) # indirect doctest
            True

        """
        if s=='__members__':
            return dir(self._p)
        if s=='__methods__':
            return [X for X in dir(self._p) if hasattr(self._p,X) and ('method' in str(type(getattr(self._p,X))))]
        try:
            return getattr(self._p,s)
        except AttributeError:
            raise AttributeError('%s has no attribute %s'%(self.__class__, s))

    def ring(self):
        """
        The ring which ``self`` belongs to.

        This is the same as ``self.parent()``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ,implementation='sparse')
            sage: p = x[100]*y[1]^3*x[1]^2+2*x[10]*y[30]
            sage: p.ring()
            Infinite polynomial ring in x, y over Integer Ring

        """
        return self.parent()

    def is_unit(self):
        r"""
        Answers whether ``self`` is a unit

        EXAMPLES::

            sage: R1.<x,y> = InfinitePolynomialRing(ZZ)
            sage: R2.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 1 + x[2]
            sage: R1.<x,y> = InfinitePolynomialRing(ZZ)
            sage: R2.<a,b> = InfinitePolynomialRing(QQ)
            sage: (1+x[2]).is_unit()
            False
            sage: R1(1).is_unit()
            True
            sage: R1(2).is_unit()
            False
            sage: R2(2).is_unit()
            True
            sage: (1+a[2]).is_unit()
            False

        TESTS::

            sage: R.<x> = InfinitePolynomialRing(ZZ.quotient_ring(8))
            sage: [R(i).is_unit() for i in range(8)]
            [False, True, False, True, False, True, False, True]
        """
        if len(self.variables()) > 0:
            return False
        else:
            return self.base_ring()(self._p).is_unit()

    @cached_method
    def variables(self):
        """
        Return the variables occurring in ``self`` (tuple of elements of some polynomial ring).

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: p = x[1] + x[2] - 2*x[1]*x[3]
            sage: p.variables()
            (x_3, x_2, x_1)
            sage: x[1].variables()
            (x_1,)
            sage: X(1).variables()
            ()

        """
        if hasattr(self._p, 'variables'):
            return tuple(self._p.variables())
        return ()

    @cached_method
    def max_index(self):
        r"""
        Return the maximal index of a variable occurring in ``self``, or -1 if ``self`` is scalar.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p=x[1]^2+y[2]^2+x[1]*x[2]*y[3]+x[1]*y[4]
            sage: p.max_index()
            4
            sage: x[0].max_index()
            0
            sage: X(10).max_index()
            -1
        """
        return max([Integer(str(X).split('_')[1]) for X in self.variables()]+[-1])

    # Basic arithmetics
    def _add_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[1] + x[2] # indirect doctest
            x_2 + x_1

        """
        # One may need a new parent for  self._p and x._p
        try:
            return InfinitePolynomial_sparse(self.parent(),self._p+x._p)
        except Exception:
            pass
        ## We can now assume that self._p and x._p actually are polynomials,
        ## hence, their parent is not simply the underlying ring.
        VarList = list(set(self._p.parent().variable_names()).union(set(x._p.parent().variable_names())))
        VarList.sort(cmp=self.parent().varname_cmp,reverse=True)
        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(),VarList,order=self.parent()._order)
        else:
            R = self._p.base_ring()
        return InfinitePolynomial_sparse(self.parent(),R(self._p) + R(x._p))

    def _mul_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: x[2]*x[1] # indirect doctest
            x_2*x_1

        """
        try:
            return InfinitePolynomial_sparse(self.parent(),self._p*x._p)
        except Exception:
            pass
        ## We can now assume that self._p and x._p actually are polynomials,
        ## hence, their parent is not just the underlying ring.
        VarList = list(set(self._p.parent().variable_names()).union(set(x._p.parent().variable_names())))
        VarList.sort(cmp=self.parent().varname_cmp,reverse=True)
        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(),VarList,order=self.parent()._order)
        else:
            R = self._p.base_ring()
        return InfinitePolynomial_sparse(self.parent(),R(self._p) * R(x._p))

    def gcd(self, x):
        """
        computes the greatest common divisor

        EXAMPLES::

            sage: R.<x>=InfinitePolynomialRing(QQ)
            sage: p1=x[0]+x[1]**2
            sage: gcd(p1,p1+3)
            1
            sage: gcd(p1,p1)==p1
            True
        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_sparse(self.parent(),self._p.gcd(x._p))

    def _rmul_(self, left):
        """
        TEST::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: R.from_base_ring(4)   # indirect doctest
            4

        """
        return InfinitePolynomial_sparse(self.parent(),left*self._p)

    def _lmul_(self, right):
        """
        TEST::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: alpha[3]*4   # indirect doctest
            4*alpha_3

        """
        return InfinitePolynomial_sparse(self.parent(),self._p*right)

    def _div_(self, x):
        r"""
        Division of Infinite Polynomials.

        EXAMPLES:

        Division by a rational over `\QQ`::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: x[0]/2
            1/2*x_0

        Division by an integer over `\ZZ`::

            sage: R.<x> = InfinitePolynomialRing(ZZ, implementation='sparse')
            sage: p = x[3]+x[2]
            sage: q = p/2
            sage: q
            1/2*x_3 + 1/2*x_2
            sage: q.parent()
            Infinite polynomial ring in x over Rational Field

        Division by a non-zero element::

            sage: R.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: 1/x[1]
            1/x_1
            sage: (x[0]/x[0])
            x_0/x_0
            sage: qt = 1/x[2]+2/x[1]; qt
            (2*x_2 + x_1)/(x_2*x_1)
            sage: qt.parent()
            Fraction Field of Infinite polynomial ring in x over Rational Field

            sage: z = 1/(x[2]*(x[1]+x[2]))+1/(x[1]*(x[1]+x[2]))
            sage: z.parent()
            Fraction Field of Infinite polynomial ring in x over Rational Field
            sage: factor(z)
            x_1^-1 * x_2^-1
        """
        if not x.variables():
            p = self.base_ring()(x._p)
            divisor = self.base_ring().one_element()/p  # use induction
            OUTP = self.parent().tensor_with_ring(divisor.base_ring())
            return OUTP(self)*OUTP(divisor)
        else:
            from sage.rings.fraction_field_element import FractionFieldElement
            field = self.parent().fraction_field()
            # there remains a problem in reduction
            return FractionFieldElement(field, self, x, reduce=False)

    def _sub_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[2] - x[1] # indirect doctest
            x_2 - x_1

        """
        try:
            return InfinitePolynomial_sparse(self.parent(),self._p-x._p)
        except Exception:
            pass
        ## We can now assume that self._p and x._p actually are polynomials,
        ## hence, their parent is not just the underlying ring.
        VarList = list(set(self._p.parent().variable_names()).union(x._p.parent().variable_names()))
        VarList.sort(cmp=self.parent().varname_cmp,reverse=True)
        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(),VarList, order=self.parent()._order)
        else:
            R = self._p.base_ring()
        return InfinitePolynomial_sparse(self.parent(),R(self._p) - R(x._p))

    def __pow__(self, n):
        """
        Exponentiation by an integer, or action by a callable object

        NOTE:

        The callable object must accept non-negative integers as input
        and return non-negative integers. Typical use case is a
        permutation, that will result in the corresponding permutation
        of variables.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ,implementation='sparse')
            sage: p = x[10]*y[2]+2*x[1]*y[3]
            sage: P = Permutation(((1,2),(3,4,5)))
            sage: p^P # indirect doctest
            x_10*y_1 + 2*x_2*y_4

        """
        P = self.parent()
        if callable(n):
            if (self._p.parent() == self._p.base_ring()):
                return self
            if not (hasattr(self._p,'variables') and self._p.variables()):
                return self
            if hasattr(n,'to_cycles') and hasattr(n,'__len__'): # duck typing Permutation
                # auxiliary function, necessary since n(m) raises an error if m>len(n)
                l = len(n)
                p = lambda m: n(m) if 0<m<=l else m
            else: # Permutation group element
                p = n
            q = lambda s: s[0]+'_'+str(p(ZZ(s[1])))
            newVars = [q(X.split('_')) for X in self._p.parent().variable_names()]
            if not newVars:
                return self
            copyVars = copy.copy(newVars)
            newVars = list(set(list(self._p.parent().variable_names())+newVars))
            newVars.sort(cmp=self.parent().varname_cmp, reverse=True)
            if newVars == list(self._p.parent().variable_names()):
                newR = self._p.parent()
            else:
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                newR = PolynomialRing(self._p.base_ring(), newVars,order=P._order)
            mapR = self._p.parent().hom(copyVars,newR)
            return InfinitePolynomial_sparse(self.parent(), mapR(self._p))
        return InfinitePolynomial_sparse(self.parent(), self._p**n)

    # Basic tools for Buchberger algorithm:
    # order, leading term/monomial, symmetric cancellation order

    def __cmp__(self, x):
        """
        Comparison of Infinite Polynomials.

        NOTE:

        Let x and y are generators of the parent of self. We only consider
        monomial orderings in which
            x[m] > y[n] iff x appears earlier in the list of generators than y, or
                            x==y and m>n
        Under this restriction, the monomial ordering can be 'lex' (default),
        'degrevlex' or 'deglex'.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a = x[10]^3
            sage: b = x[1] + x[2]
            sage: c = x[1] + x[2]
            sage: d = y[1] + x[2]
            sage: a == a # indirect doctest
            True
            sage: b == c # indirect doctest
            True
            sage: a == b # indirect doctest
            False
            sage: c > d # indirect doctest
            True

        TESTS::

        A classical and an infinite sparse polynomial ring. Note that
        the Sage coercion system allows comparison only if a common
        parent for the two rings can be constructed. This is why we
        have to have the order 'degrevlex'.
        ::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ,order='degrevlex', implementation='sparse')
            sage: Y.<z,x_3,x_1> = QQ[]
            sage: x[3] == x_3 # indirect doctest
            True

        Two infinite polynomial rings in different implementation and
        order::

            sage: Y = InfinitePolynomialRing(QQ,['x','y'],order='deglex',implementation='dense')
            sage: x[2] == Y(x[2]) # indirect doctest
            True

        An example in which a previous version had failed::

            sage: X.<x,y> = InfinitePolynomialRing(GF(3), order='degrevlex', implementation='sparse')
            sage: p = Y('x_3*x_0^2 + x_0*y_3*y_0')
            sage: q = Y('x_1*x_0^2 + x_0*y_1*y_0')
            sage: p < q # indirect doctest
            False

        """
        # We can assume that self.parent() is x.parent(),
        # but of course the underlying polynomial rings
        # may be widely different, and the sage coercion
        # system can't guess what order we want.
        from sage.all import parent
        R1 = parent(self._p)
        R2 = parent(x._p)
        if (hasattr(R1,'has_coerce_map_from') and R1.has_coerce_map_from(R2)) or (hasattr(R2,'has_coerce_map_from') and R2.has_coerce_map_from(R1)):
            return cmp(self._p, x._p)
        VarList = list(set(self._p.parent().variable_names()).union(x._p.parent().variable_names()))
        VarList.sort(cmp=self.parent().varname_cmp,reverse=True)
        if VarList:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(self._p.base_ring(),VarList,order=self.parent()._order)
        else:
            R = self._p.base_ring()
        if (self._p.parent() is self._p.base_ring()) or not self._p.parent().gens():
                fself = self._p.base_ring()
        else:
            fself = self._p.parent().hom(self._p.parent().variable_names(),R)
        if (x._p.parent() is x._p.base_ring()) or not x._p.parent().gens():
                fx = x._p.base_ring()
        else:
            fx = x._p.parent().hom(x._p.parent().variable_names(),R)
        return cmp(fself(self._p), fx(x._p))

    @cached_method
    def lm(self):
        """
        The leading monomial of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30]+x[10]*y[1]^3*x[1]^2
            sage: p.lm()
            x_10*x_1^2*y_1^3

        """
        if hasattr(self._p,'lm'):
            return InfinitePolynomial(self.parent(), self._p.lm())
        if self._p==0:
            return self
        if hasattr(self._p,'variable_name'): # if it is univariate
            return InfinitePolynomial(self.parent(),self._p.parent().gen()**max(self._p.exponents()))
        return self # if it is scalar

    @cached_method
    def lc(self):
        """
        The coefficient of the leading term of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30]+3*x[10]*y[1]^3*x[1]^2
            sage: p.lc()
            3

        """
        if hasattr(self._p,'lc'):
            return self._p.lc()
        if hasattr(self._p,'variable_name'): # univariate case
            return self._p.leading_coefficient()
        # scalar case
        return self._p

    @cached_method
    def lt(self):
        """
        The leading term (= product of coefficient and monomial) of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30]+3*x[10]*y[1]^3*x[1]^2
            sage: p.lt()
            3*x_10*x_1^2*y_1^3

        """
        if hasattr(self._p,'lt'):
            return InfinitePolynomial(self.parent(), self._p.lt())
        if self._p==0:
            return self
        if hasattr(self._p,'variable_name'): # if it is univariate
            return InfinitePolynomial(self.parent(), self._p.leading_coefficient()*self._p.parent().gen()**max(self._p.exponents()))
        return self # if it is scalar

    def tail(self):
        """
        The tail of ``self`` (this is ``self`` minus its leading term).

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = 2*x[10]*y[30]+3*x[10]*y[1]^3*x[1]^2
            sage: p.tail()
            2*x_10*y_30

        """
        return self-self.lt()

    def squeezed(self):
        """
        Reduce the variable indices occurring in ``self``.

        OUTPUT:

        Apply a permutation to ``self`` that does not change the order of
        the variable indices of ``self`` but squeezes them into the range
        1,2,...

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ,implementation='sparse')
            sage: p = x[1]*y[100] + x[50]*y[1000]
            sage: p.squeezed()
            x_2*y_4 + x_1*y_3

        """
        Indices = set([0]+[Integer(str(Y).split('_')[1]) for Y in self.variables()])
        Indices = sorted(Indices)
        P = lambda n: Indices.index(n) if Indices.__contains__(n) else n
        return self**P

    def footprint(self):
        """
        Leading exponents sorted by index and generator.

        OUTPUT:

        ``D`` -- a dictionary whose keys are the occurring variable indices.

        ``D[s]`` is a list ``[i_1,...,i_n]``, where ``i_j`` gives the
        exponent of ``self.parent().gen(j)[s]`` in the leading
        term of ``self``.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = x[30]*y[1]^3*x[1]^2+2*x[10]*y[30]
            sage: sorted(p.footprint().items())
            [(1, [2, 3]), (30, [1, 0])]

        TESTS:

        This is a test whether it also works when the underlying polynomial ring is
        not implemented in libsingular::

            sage: X.<x> = InfinitePolynomialRing(ZZ)
            sage: Y.<y,z> = X[]
            sage: Z.<a> = InfinitePolynomialRing(Y)
            sage: Z
            Infinite polynomial ring in a over Multivariate Polynomial Ring in y, z over Infinite polynomial ring in x over Integer Ring
            sage: type(Z._P)
            <class 'sage.rings.polynomial.multi_polynomial_ring.MPolynomialRing_polydict_with_category'>
            sage: p = a[12]^3*a[2]^7*a[4] + a[4]*a[2]
            sage: sorted(p.footprint().items())
            [(2, [7]), (4, [1]), (12, [3])]

        """
        if not self._has_footprint:
            PARENT = self.parent()
            l = len(self.parent()._names)
            # get the pairs (shift,exponent) of the leading monomial, indexed by the variable names
            Vars = self._p.parent().variable_names()
            from sage.rings.polynomial.multi_polynomial_libsingular import MPolynomial_libsingular
            if isinstance(self._p, MPolynomial_libsingular):
                L = [(Vars[i].split('_'),e) for i,e in enumerate(self._p.lm().exponents(as_ETuples=False)[0]) if e]
            elif hasattr(self._p,'lm'):
                # self._p  is multivariate, but not libsingular, hence,
                # exponents is slow and does not accept the optional argument as_ETuples.
                # Thus, fall back to regular expressions
                L = PARENT._find_varpowers.findall(repr(self.lm()._p))
                L = [((x[0:2]),int(x[2]) if x[2] else 1) for x in L]
            else: # it is a univariate polynomial -- this should never happen, but just in case...
                L = [(Vars[0].split('_'),self._p.degree())]
            for t in L:
                n = t[0][0]       # the variable *n*ame
                s = int(t[0][1])  # the variable *s*hift
                if s not in self._footprint:
                    self._footprint[s] = [0]*l
                self._footprint[s][self.parent()._name_dict[n]] = t[1]   # the exponent
            self._has_footprint = True
        return self._footprint

    def symmetric_cancellation_order(self,other):
        """
        Comparison of leading terms by Symmetric Cancellation Order, `<_{sc}`.

        INPUT:

        self, other -- two Infinite Polynomials

        ASSUMPTION:

        Both Infinite Polynomials are non-zero.

        OUTPUT:

        ``(c, sigma, w)``, where

        * c = -1,0,1, or None if the leading monomial of ``self`` is smaller, equal,
          greater, or incomparable with respect to ``other`` in the monomial
          ordering of the Infinite Polynomial Ring
        * sigma is a permutation witnessing
          ``self`` `<_{sc}` ``other`` (resp. ``self`` `>_{sc}` ``other``)
          or is 1 if ``self.lm()==other.lm()``
        * w is 1 or is a term so that
          ``w*self.lt()^sigma == other.lt()`` if `c\\le 0`, and
          ``w*other.lt()^sigma == self.lt()`` if `c=1`

        THEORY:

        If the Symmetric Cancellation Order is a well-quasi-ordering
        then computation of Groebner bases always terminates. This is
        the case, e.g., if the monomial order is lexicographic. For
        that reason, lexicographic order is our default order.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: (x[2]*x[1]).symmetric_cancellation_order(x[2]^2)
            (None, 1, 1)
            sage: (x[2]*x[1]).symmetric_cancellation_order(x[2]*x[3]*y[1])
            (-1, [2, 3, 1], y_1)
            sage: (x[2]*x[1]*y[1]).symmetric_cancellation_order(x[2]*x[3]*y[1])
            (None, 1, 1)
            sage: (x[2]*x[1]*y[1]).symmetric_cancellation_order(x[2]*x[3]*y[2])
            (-1, [2, 3, 1], 1)

        """
        PARENT = self.parent()
        other = PARENT(other)
        slt = self.lt()
        olt = other.lt()
        rawcmp = cmp(self.lm(),other.lm())
        if rawcmp == 0:
            if olt==0:
                return (0, 1, 1)
            return (0, 1, self.lc()/other.lc())
        if rawcmp == -1:
            Fsmall = dict([[k[0], [e for e in k[1]]] for k in self.footprint().items()])
            Fbig = dict([[k[0], [e for e in k[1]]] for k in other.footprint().items()])
            ltsmall = slt
            ltbig = olt
        else:
            Fbig = dict([[k[0], [e for e in k[1]]] for k in self.footprint().items()])
            Fsmall = dict([[k[0], [e for e in k[1]]] for k in other.footprint().items()])
            ltbig = slt
            ltsmall = olt
        # Case 1: one of the Infintie Polynomials is scalar.
        if not Fsmall:
            return (rawcmp, 1, ltbig/ltsmall)
        # "not Fbig" is now impossible, because we only consider *global* monomial orderings.
        # These are the occurring shifts:
        Lsmall = sorted(Fsmall.keys())
        Lbig   = Fbig.keys()
        Lbig.sort()
        P = range(Lbig[-1]+1)
        gens = xrange(PARENT.ngens())
        if Lsmall[0]==0:
            if 0 not in Fbig:
                return (None,1,1)
            Lsmall.pop(0)
            Lbig.pop(0)
            ExpoSmall = Fsmall[0]
            ExpoBig = Fbig[0]
            for k in gens:
                if ExpoBig[k]<ExpoSmall[k]:
                    return (None,1,1)
                ExpoBig[k]-=ExpoSmall[k]
        lenBig = len(Lbig)
        j = -1              # will have Lbig[j] -- a shift of the bigger polynomial
        for i in Lsmall:   # i is a shift of the smaller polynomial
            j += 1
            ExpoSmall = Fsmall[i]
            while (j<lenBig):
                found = False
                if Lbig[j]>=i:
                    ExpoBigSave = [e for e in Fbig[Lbig[j]]]
                    ExpoBig = Fbig[Lbig[j]]
                    found = True
                    for k in gens:
                        if ExpoBig[k]<ExpoSmall[k]:
                            found = False
                            Fbig[Lbig[j]] = ExpoBigSave
                            break
                        ExpoBig[k]-=ExpoSmall[k]
                if found:
                    break
                j+=1
            if j==lenBig:
                return (None,1,1)  ## no "increasing" permutation transforms
                                    ## the smaller monomial into a factor of
                                    ## the bigger monomial
            tmp = P[i]
            P[i] = Lbig[j]
            P[Lbig[j]]=tmp
        # now, P defines an 'up-shift' permutation, slt^P divides olt, and
        # Fbig contains the exponents for olt/slt^P.
        OUT = PARENT(PARENT._base(ltbig.lc()/ltsmall.lc()))
        for shift, Expo in Fbig.items():
            for g in gens:
                if Expo[g]:
                    OUT *= PARENT.gen(g)[shift].__pow__(Expo[g])
        from sage.combinat.permutation import Permutation
        return (rawcmp, Permutation(P[1:]), OUT)

    def coefficient(self, monomial):
        """
        Returns the coefficient of a monomial in this polynomial.

        INPUT:

        - A monomial (element of the parent of self) or
        - a dictionary that describes a monomial (the keys
          are variables of the parent of self, the values
          are the corresponding exponents)

        EXAMPLES:

        We can get the coefficient in front of monomials::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = 2*x[0]*x[1] + x[1] + x[2]
            sage: a.coefficient(x[0])
            2*x_1
            sage: a.coefficient(x[1])
            2*x_0 + 1
            sage: a.coefficient(x[2])
            1
            sage: a.coefficient(x[0]*x[1])
            2

        We can also pass in a dictionary::

            sage: a.coefficient({x[0]:1, x[1]:1})
            2

        """
        if self._p==0:
            res = 0
        elif isinstance(monomial, self.__class__):
            if not (self.parent().has_coerce_map_from(monomial.parent())):
                res = 0
            else:
                if hasattr(self._p,'variables'):
                    VarList = [str(X) for X in self._p.variables()]
                else:
                    VarList = []
                if hasattr(monomial._p,'variables'):
                    VarList.extend([str(X) for X in monomial._p.variables()])
                VarList = list(set(VarList))
                VarList.sort(cmp=self.parent().varname_cmp,reverse=True)
                from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
                if len(VarList)==1:
                    R = PolynomialRing(self._p.base_ring(),VarList+['xx'],order=self.parent()._order)
                                                                 ## 'xx' is guaranteed to be no variable
                                                                 ## name of monomial, since coercions
                                                                 ## were tested before
                    res = PolynomialRing(self._p.base_ring(),VarList,order=self.parent()._order)(R(self._p).coefficient(R(monomial._p)))
                else:
                    R = PolynomialRing(self._p.base_ring(),VarList,order=self.parent()._order)
                    res = R(self._p).coefficient(R(monomial._p))
        elif isinstance(monomial, dict):
            if monomial:
                I = monomial.iterkeys()
                K = I.next()
                monomial.__delitem__(K)
                res = self.coefficient(K).coefficient(monomial)
            else:
                return self
        else:
            raise TypeError("Objects of type %s have no coefficients in InfinitePolynomials"%(type(monomial)))
        return self.parent()(res)

    ## Essentials for Buchberger
    def reduce(self, I, tailreduce=False, report=None):
        """
        Symmetrical reduction of ``self`` with respect to a symmetric ideal (or list of Infinite Polynomials).

        INPUT:

        - ``I`` -- a :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal` or a list
          of Infinite Polynomials.
        - ``tailreduce`` -- (bool, default ``False``) *Tail reduction* is performed if this
          parameter is ``True``.
        - ``report`` -- (object, default ``None``) If not ``None``, some information on the
          progress of computation is printed, since reduction of huge polynomials may take
          a long time.

        OUTPUT:

        Symmetrical reduction of ``self`` with respect to ``I``, possibly with tail reduction.

        THEORY:

        Reducing an element `p` of an Infinite Polynomial Ring `X` by
        some other element `q` means the following:

        1. Let `M` and `N` be the leading terms of `p` and `q`.
        2. Test whether there is a permutation `P` that does not does not diminish the variable
           indices occurring in `N` and preserves their order, so that there is some term `T\in X`
           with `TN^P = M`. If there is no such permutation, return `p`
        3. Replace `p` by `p-T q^P` and continue with step 1.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: p = y[1]^2*y[3]+y[2]*x[3]^3
            sage: p.reduce([y[2]*x[1]^2])
            x_3^3*y_2 + y_3*y_1^2

        The preceding is correct: If a permutation turns
        ``y[2]*x[1]^2`` into a factor of the leading monomial
        ``y[2]*x[3]^3`` of ``p``, then it interchanges the variable
        indices 1 and 2; this is not allowed in a symmetric
        reduction. However, reduction by ``y[1]*x[2]^2`` works, since
        one can change variable index 1 into 2 and 2 into 3::

            sage: p.reduce([y[1]*x[2]^2])
            y_3*y_1^2

        The next example shows that tail reduction is not done, unless
        it is explicitly advised.  The input can also be a Symmetric
        Ideal::

            sage: I = (y[3])*X
            sage: p.reduce(I)
            x_3^3*y_2 + y_3*y_1^2
            sage: p.reduce(I, tailreduce=True)
            x_3^3*y_2

        Last, we demonstrate the ``report`` option::

            sage: p=x[1]^2+y[2]^2+x[1]*x[2]*y[3]+x[1]*y[4]
            sage: p.reduce(I, tailreduce=True, report=True)
            :T[2]:>
            >
            x_1^2 + y_2^2

        The output ':' means that there was one reduction of the
        leading monomial. 'T[2]' means that a tail reduction was
        performed on a polynomial with two terms. At '>', one round of
        the reduction process is finished (there could only be several
        non-trivial rounds if `I` was generated by more than one
        polynomial).

        """
        from sage.rings.polynomial.symmetric_reduction import SymmetricReductionStrategy
        if hasattr(I,'gens'):
            I = I.gens()
        if (not I):
            return self
        I = list(I)
        S = SymmetricReductionStrategy(self.parent(),I, tailreduce)
        return S.reduce(self, report=report)

    ## Further methods
    def stretch(self, k):
        """
        Stretch ``self`` by a given factor.

        INPUT:

        ``k`` -- an integer.

        OUTPUT:

        Replace `v_n` with `v_{n\\cdot k}` for all generators `v_\\ast` occurring in self.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1] + x[2]
            sage: a.stretch(2)
            x_4 + x_2 + x_0

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1] + y[0]*y[1]; a
            x_1 + x_0 + y_1*y_0
            sage: a.stretch(2)
            x_2 + x_0 + y_2*y_0

        TESTS:

        The following would hardly work in a dense implementation,
        because an underlying polynomial ring with 6001 variables
        would be created. This is avoided in the sparse
        implementation::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: a = x[2] + x[3]
            sage: a.stretch(2000)
            x_6000 + x_4000

        """
        P = lambda n: k*n
        return self.__pow__(P)


class InfinitePolynomial_dense(InfinitePolynomial_sparse):
    """
    Element of a dense Polynomial Ring with a Countably Infinite Number of Variables.

    INPUT:

    - ``A`` -- an Infinite Polynomial Ring in dense implementation
    - ``p`` -- a *classical* polynomial that can be interpreted in ``A``.

    Of course, one should not directly invoke this class, but rather
    construct elements of ``A`` in the usual way.

    This class inherits from
    :class:`~sage.rings.polynomial.infinite_polynomial_element.InfinitePolynomial_sparse`. See
    there for a description of the methods.
    """
    # Construction and other basic methods
##    def __init__(self, A, p): # is inherited from the dense implementation

    def __call__(self, *args, **kwargs):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: a = x[0] + x[1]
            sage: a(x_0=2,x_1=x[1])
            x_1 + 2
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field
            sage: a(x_1=3)
            x_0 + 3
            sage: _.parent()
            Infinite polynomial ring in x over Rational Field

            sage: a(x_1=x[100])
            x_100 + x_0

        """
        #Replace any InfinitePolynomials by their underlying polynomials
        for kw in kwargs:
            value = kwargs[kw]
            if isinstance(value, InfinitePolynomial_sparse):
                kwargs[kw] = value._p
        args = list(args)
        for i, arg in enumerate(args):
            if isinstance(arg, InfinitePolynomial_sparse):
                args[i] = arg._p
        self._p = self.parent().polynomial_ring()(self._p)
        res = self._p(*args, **kwargs)
        try:
            return self.parent()(res)
        except ValueError:
            return res

    def __cmp__(self, x):
        """
        TESTS::

        A classical and an infinite polynomial ring::

            sage: X.<x,y> = InfinitePolynomialRing(ZZ,order='degrevlex')
            sage: Y.<z,x_3,x_1> = QQ[]
            sage: x[3] == x_3
            True

        Two infinite polynomial rings with different order and
        implementation::

            sage: Y = InfinitePolynomialRing(QQ,['x','y'],order='deglex',implementation='sparse')
            sage: x[2] == Y(x[2])
            True

        An example in which a previous version had failed::

            sage: X.<x,y> = InfinitePolynomialRing(GF(3), order='degrevlex', implementation='dense')
            sage: p = Y('x_3*x_0^2 + x_0*y_3*y_0')
            sage: q = Y('x_1*x_0^2 + x_0*y_1*y_0')
            sage: p < q
            False

        """
        # We can assume that self and x belong to the same ring.
        # We can not assume yet that self._p and
        # x._p are already defined over self.parent()._P
        # It won't hurt to change self in place.
        # But, to be on the safe side...
        try:
            self._p = self.parent()._P(self._p)
        except Exception:
            pass
        try:
            x._p = x.parent()._P(x._p)
        except Exception:
            pass
        return cmp(self._p,x._p)

    # Basic arithmetics
    def _add_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[1] + x[2] # indirect doctest
            x_2 + x_1

        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(self.parent(),self._p + x._p)

    def _mul_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[2]*x[1] # indirect doctest
            x_2*x_1

        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(self.parent(),self._p * x._p)

    def _rmul_(self, left):
        """
        TEST::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ)
            sage: R.from_base_ring(4)   # indirect doctest
            4

        """
        return InfinitePolynomial_dense(self.parent(),left*self._p)

    def _lmul_(self, right):
        """
        TEST::

            sage: R.<alpha,beta> = InfinitePolynomialRing(QQ)
            sage: alpha[3]*4   # indirect doctest
            4*alpha_3

        """
        return InfinitePolynomial_dense(self.parent(),self._p*right)

    def _sub_(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[2] - x[1] # indirect doctest
            x_2 - x_1

        """
        P = self.parent()
        self._p = P._P(self._p)
        x._p = P._P(x._p)
        return InfinitePolynomial_dense(self.parent(), self._p - x._p)



    def __pow__(self, n):
        """
        Exponentiation by an integer, or action by a callable object

        NOTE:

        The callable object must accept non-negative integers as input
        and return non-negative integers. Typical use case is a
        permutation, that will result in the corresponding permutation
        of variables.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: x[10]^3
            x_10^3
            sage: p = x[10]*y[2]+2*x[1]*y[3]
            sage: P = Permutation(((1,2),(3,4,5)))
            sage: p^P
            x_10*y_1 + 2*x_2*y_4

        """
        P = self.parent()
        if callable(n):
            if (self._p.parent() == self._p.base_ring()):
                return self
            if not (hasattr(self._p,'variables') and self._p.variables()):
                return self
            if hasattr(n,'to_cycles') and hasattr(n,'__len__'): # duck typing Permutation
                # auxiliary function, necessary since n(m) raises an error if m>len(n)
                l = len(n)
                p = lambda m: n(m) if 0<m<=l else m
            else: # Permutation group element
                p = n

            # determine whether the maximal index must be raised
            oldMax = P._max
            newMax = max([p(X) for X in range(oldMax+1)]+[oldMax])
            if newMax > P._max:
                P.gen()[newMax]
            self._p = P._P(self._p)
            # next, determine the images of variable names
            PP = P._P
            PPgens = PP.gens()

            newVars = []
            sh = PP.ngens()//P.ngens() - 1
            blocklength = sh
            nM = sh+1
            for i in range(P.ngens()):
                newVars.extend([PPgens[sh-p(j)] for j in range(blocklength,-1,-1)])
                sh += nM
            mapR = PP.hom(newVars, PP)
            return InfinitePolynomial_dense(self.parent(), mapR(self._p))

        # else, n is supposed to be an integer
        return InfinitePolynomial_dense(self.parent(), self._p**n)

