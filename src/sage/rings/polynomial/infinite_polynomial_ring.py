"""
Infinite Polynomial Rings

By Infinite Polynomial Rings, we mean polynomial rings in a countably infinite
number of variables. The implementation consists of a wrapper around the
current *finite* polynomial rings in Sage.

AUTHORS:

- Simon King <simon.king@uni-jena.de>
- Mike Hansen <mhansen@gmail.com>

An Infinite Polynomial Ring has finitely many generators `x_\\ast, y_\\ast,...`
and infinitely many variables of the form `x_0, x_1, x_2, ..., y_0, y_1, y_2,...,...`.
We refer to the natural number `n` as the *index* of the variable `x_n`.

INPUT:

- ``R``, the base ring, which must in fact be a *field*
- ``names``, a list of generator names. Generator names must be pairwise different. We only
  allow generator names that are single characters.
- ``order`` (optional string). The default order is ``'lex'`` (lexicographic). ``'deglex'``
  is degree lexicographic, and ``'degrevlex'`` (degree reverse lexicographic) is possible
  but discouraged.

Each generator ``x`` produces an infinite sequence of variables ``x[1], x[2], ...``
which are printed on screen as ``x1, x2, ...`` and are latex typeset as `x_{1}, x_{2}`.
Then, the Infinite Polynomial Ring is formed by polynomials in these variables.

By default, the monomials are ordered lexicographically. Alternatively,
degree (reverse) lexicographic ordering is possible as well. However, we do not
guarantee that the computation of Groebner bases will terminate in this case.

In either case, the variables of a Infinite Polynomial Ring X are ordered
according to the following rule:
  ``X.gen(i)[m] < X.gen(j)[n]`` if and only if ``i<j or (i==j and m<n)``

We provide a 'dense' and a 'sparse' implementation. In the dense implementation,
the Infinite Polynomial Ring carries a finite polynomial ring that comprises *all*
variables up to the maximal index that has been used so far. This is potentially
a very big ring and may also comprise many variables that are not used.

In the sparse implementation, we try to keep the underlying finite polynomial rings
small, using only those variables that are really needed. By default, we use the
dense implementation, since it usually is much faster.

We provide coercion from an Infinite Polynomial Ring ``X`` over a
field ``F_X`` to an Infinite Polynomial ring ``Y`` over a field ``F_Y`` (regardless
of dense or sparse implementation and of monomial ordering) if and only if there is
a coercion from ``F_X`` to ``F_Y``, and the set of generator names of ``X`` is a
subset of the set of generator names of ``Y``. The coercion map is name preserving.
We also allow coercion from a *classical* polynomial ring to ``X`` if base rings
and variable names fit.


EXAMPLES::

    sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
    sage: A.<a,b> = InfinitePolynomialRing(QQ, order='deglex')

    sage: f = x[5] + 2; f
    x5 + 2
    sage: g = 3*y[1]; g
    3*y1
    sage: g._p.parent()
    Univariate Polynomial Ring in y1 over Rational Field

    sage: f2 = a[5] + 2; f2
    a5 + 2
    sage: g2 = 3*b[1]; g2
    3*b1
    sage: A.polynomial_ring()
    Multivariate Polynomial Ring in b5, b4, b3, b2, b1, b0, a5, a4, a3, a2, a1, a0 over Rational Field

Of course, we provide the usual polynomial arithmetic::

    sage: f+g
    3*y1 + x5 + 2
    sage: p = x[10]^2*(f+g); p
    3*y1*x10^2 + x10^2*x5 + 2*x10^2

There is a permutation action on the variables, by permuting positive variable indices::

    sage: P = Permutation(((10,1)))
    sage: p^P
    3*y10*x1^2 + x5*x1^2 + 2*x1^2

Note that `x_0^P = x_0`, since the permutations only change *positive* variable indices.

We also implemented ideals of Infinite Polynomial Rings. Here, it is thoroughly assumed that
the ideals are set-wise invariant under the permutation action. We therefore refer to these
ideals as *Symmetric Ideals*. Symmetric Ideals are finitely generated modulo addition, multiplication
by ring elements and permutation of variables, and (at least in the default case of a lexicographic
order), one can compute Groebner bases::

    sage: I = (x[1]*y[2])*X
    sage: I.groebner_basis()
    [y1*x2, y2*x1]
    sage: J = A*(a[1]*b[2])
    sage: J.groebner_basis()
    [b1*a2, b2*a1]

For more details, see :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal`.

"""
#*****************************************************************************
#       Copyright (C) 2009 Simon King <simon.king@uni-jena.de> and
#                          Mike Hansen <mhansen@gmail.com>,
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
from sage.rings.ring import CommutativeRing
from sage.structure.all import Parent, SageObject
from sage.structure.factory import UniqueFactory
from sage.misc.cachefunc import cached_method
import copy, operator, sys, re

###############################################################
## Ring Factory framework

class InfinitePolynomialRingFactory(UniqueFactory):
    """
    A factory for creating infinite polynomial ring elements.  It
    handles making sure that they are unique as well as handling
    pickling.  For more details, see
    :class:`~sage.structure.factory.UniqueFactory` and
    :mod:`~sage.rings.polynomial.infinite_polynomial_ring`.

    EXAMPLES::

        sage: X.<x> = InfinitePolynomialRing(QQ)
        sage: X2.<x> = InfinitePolynomialRing(QQ)
        sage: X is X2
        True
        sage: X3.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: X is X3
        False

        sage: X is loads(dumps(X))
        True

    """
    def create_key(self, R, names=('x',), order='lex', implementation='dense'):
        """
        Creates a key which uniquely defines the infinite polynomial ring.

        TESTS::

            sage: InfinitePolynomialRing.create_key(QQ)
            (Rational Field, ('x',), 'lex', 'dense')
            sage: InfinitePolynomialRing.create_key(QQ, 'y')
            (Rational Field, 'y', 'lex', 'dense')
            sage: InfinitePolynomialRing.create_key(QQ, names='y', order='deglex', implementation='sparse')
            (Rational Field, 'y', 'deglex', 'sparse')
            sage: InfinitePolynomialRing.create_key(QQ, names=['x','y'], implementation='dense')
            (Rational Field, ('x', 'y'), 'lex', 'dense')
        """
        if isinstance(names, list):
            names = tuple(names)
        return (R, names, order, implementation)

    def create_object(self, version, key):
        """
        Returns the infinite polynomial ring corresponding to the key ``key``.

        TESTS::

            sage: InfinitePolynomialRing.create_object('1.0', (QQ, 'x', 'deglex', 'sparse'))
            Infinite polynomial ring in x over Rational Field

        """
        if key[-1]=='dense':
            return InfinitePolynomialRing_dense(*key[:-1])
        return InfinitePolynomialRing_sparse(*key[:-1])

InfinitePolynomialRing = InfinitePolynomialRingFactory('InfinitePolynomialRing')

##############################################################
##  The sparse implementation

class InfinitePolynomialRing_sparse(CommutativeRing):
    """
    Sparse implementation of Infinite Polynomial Rings.

    An Infinite Polynomial Ring with generators `x_\\ast, y_\\ast, ...` over a field
    `F` is a free commutative `F`-algebra generated by `x_0, x_1, x_2, ...,
    y_0, y_1, y_2, ..., ...` and is equipped with a permutation action
    on the generators, namely `x_n^P = x_{P(n)}, y_{n}^P=y_{P(n)}, ...`
    for any permutation `P` (note that variables of index zero are invariant
    under such permutation).

    It is known that any permutation invariant ideal in an Infinite Polynomial Ring
    is finitely generated modulo the permutation action --
    see :class:`~sage.rings.polynomial.symmetric_ideal.SymmetricIdeal`
    for more details.

    Usually, an instance of this class is created using ``InfinitePolynomialRing``
    with the optional parameter ``implementation='sparse'``. This takes care
    of uniqueness of parent structures. However, a direct construction is possible,
    in principle::

        sage: X.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: Y.<x,y> = InfinitePolynomialRing(QQ, implementation='sparse')
        sage: X is Y
        True
        sage: from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing_sparse
        sage: Z = InfinitePolynomialRing_sparse(QQ, ['x','y'], 'lex')
        sage: Z == X
        True
        sage: Z is X
        False

    The last parameter ('lex' in the above example) can also be 'deglex' or
    'degrevlex'; this would result in an Infinite Polynomial Ring in degree
    lexicographic or degree reverse lexicographic order.

    See :mod:`~sage.rings.polynomial.infinite_polynomial_ring` for more details.

    """
    def __init__(self, R, names, order):
        """
        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)

        Infinite Polynomial Rings are unique parent structures::

            sage: X is loads(dumps(X))
            True
            sage: p=x[10]*y[2]^3+2*x[1]*y[3]
            sage: p
            2*y3*x1 + y2^3*x10

        We define another Infinite Polynomial Ring with same generator
        names but a different order. These rings are different, but
        allow for coercion::

            sage: Y.<x,y> = InfinitePolynomialRing(QQ, order='deglex', implementation='dense')
            sage: Y is X
            False
            sage: q=y[2]^3*x[10]+x[1]*y[3]*2
            sage: q
            y2^3*x10 + 2*y3*x1
            sage: p==q
            True
            sage: X.gen(1)[2]*Y.gen(0)[1]
            y2*x1

        As usual, if there are coercion maps in both direction, the parent
        of the left operand is taken::

            sage: (X.gen(1)[2]*Y.gen(0)[1]).parent() is X
            True

        """
        for n in names:
            if (len(n) > 1) or not n.isalpha():
                raise ValueError, "variable names must be of length 1"
        if len(names)!=len(set(names)):
            raise ValueError, "variable names must be pairwise different"
        self._names = list(names)
        if not isinstance(order, basestring):
            raise TypeError, "The monomial order must be given as a string"
        if not (hasattr(R,'is_field') and R.is_field()):
            raise TypeError, "The base ring (= %s) must be a field"%(R)

        # now, the input is accepted
        self._order = order
        self._name_dict = dict([(names[i],i) for i in xrange(len(names))])
        CommutativeRing.__init__(self, base=R)
        self._populate_coercion_lists_()
        self._varpattern = re.compile('[%s]\\d+'%(''.join(names)))
        self._monpattern = re.compile('[%s]\\d+\\^?\\d*'%(''.join(names)))
        self._exppattern = re.compile('\\^')

    def __repr__(self):
        """
        EXAMPLES::

            sage: InfinitePolynomialRing(QQ)
            Infinite polynomial ring in x over Rational Field

            sage: X.<x,y> = InfinitePolynomialRing(QQ, order='deglex'); X
            Infinite polynomial ring in x, y over Rational Field

        """
        return "Infinite polynomial ring in %s over %s"%(", ".join(self._names), self._base)

    def _latex_(self):
        """
        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: latex(X)
            \Bold{Q}[x_{\ast}, y_{\ast}]
        """
        from sage.misc.latex import latex
        vars = ', '.join([latex(X) for X in self.gens()])
        return "%s[%s]"%(latex(self.base_ring()), vars)

    @cached_method
    def _coerce_map_from_(self, S):
        """
        Coerce things into self

        EXAMPLES::

        Here, we check to see that elements of a *finitely* generated polynomial ring
        with appropriate variable names coerce correctly into the Infinite Polynomial
        Ring::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: px0 = PolynomialRing(QQ,'x0').gen(0)
            sage: px0 + x[0]
            2*x0
            sage: px0==x[0]
            True

        Note that any coercion will preserve variable names.

        """
        # Case 1: another infinite polynomial ring
        if isinstance(S, InfinitePolynomialRing_sparse):
            return self._base.has_coerce_map_from(S.base_ring()) and set(S._names).issubset(set(self._names))
        # Case 2: a field (e.g. fraction field) compatible with the base field
        try:
            if S.fraction_field() is S: # Note that the base ring of self is a field
                return self._base.has_coerce_map_from(S)
        except (NotImplementedError, AttributeError):
            pass
        # Case 3: a plain ring that coerces into our base field (e.g. ZZ into QQ)
        try:
            if S.base_ring() is S:
                return self._base.has_coerce_map_from(S)
        except AttributeError:
            pass
        # Case 4: a polynomial ring
        try:
            for n in S.variable_names():
                if (not self._name_dict.has_key(n[0])) or (not n[1:].isdigit()):
                    return False
            return self._base.has_coerce_map_from(S.base_ring())
        except AttributeError:
            return False

    def _element_constructor_(self, x):
        """
        EXAMPLES::

            sage: X = InfinitePolynomialRing(QQ)
            sage: a = X(2); a
            2
            sage: a.parent()
            Infinite polynomial ring in x over Rational Field
            sage: R=PolynomialRing(ZZ,['x3'])
            sage: b = X(R.gen()); b
            x3
            sage: b.parent()
            Infinite polynomial ring in x over Rational Field
            sage: X('(x1^2+2/3*x4)*(x2+x5)')
            2/3*x5*x4 + x5*x1^2 + 2/3*x4*x2 + x2*x1^2

        """
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        try:
            if hasattr(x, 'parent'):
                R = x.parent()
                if R is self:
                    return x
                if self.has_coerce_map_from(R):
                    return InfinitePolynomial(self, x)
            return InfinitePolynomial(self,str(x))
        except ValueError:
            raise ValueError, "Can't convert %s into an element of %s"%(x,self)

    ## Auxiliary function for variable comparison
    def varname_cmp(self,x,y):
        """
        Comparison of two variable names

        INPUT:
            ``x,y`` should both be strings of the form ``a+str(n)``, where
            a is a single character appearing in the list of generator names,
            and n is an integer

        RETURN:
            -1,0,1 if x<y, x==y, x>y, respectively, where the order
            is defined as follows:
              x<y `\\iff` the letter ``x[0]`` is earlier in the list of
              generator names of self than ``y[0]``,
              or (``x[0]==y[0]`` and ``int(x[1:])<int(y[1:])``)

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: X.varname_cmp('x1','y10')
            -1
            sage: X.varname_cmp('y1','x10')
            1
            sage: X.varname_cmp('y1','y10')
            -1

        """
        try:
            return cmp([self._name_dict[x[0]],int(x[1:])],[self._name_dict[y[0]],int(y[1:])])
        except (KeyError, ValueError):
            raise ValueError, "%s or %s is not a valid variable name"%(x,y)

    def __cmp__(self, x):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: X2.<x> = InfinitePolynomialRing(QQ)
            sage: X3.<x> = InfinitePolynomialRing(QQ, order='deglex')
            sage: Y.<y> = InfinitePolynomialRing(QQ)
            sage: Z.<z> = InfinitePolynomialRing(GF(5))
            sage: X == X
            True
            sage: X == X2
            True
            sage: X == X3
            False
            sage: X == Y
            False
            sage: X == Z
            False

        """
        if not isinstance(x, InfinitePolynomialRing_sparse):
            return -1
        return cmp( (self._base, self._names, self._order), (x.base_ring(), x._names, x._order) )

    def ngens(self):
        """
        Returns the number of generators for this ring.  Since there
        are countably infinitely many variables in this polynomial
        ring, by 'generators' we mean the number of infinite families
        of variables. See :mod:`~sage.rings.polynomial.infinite_polynomial_ring`
        for more details.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: X.ngens()
            1

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: X.ngens()
            2

        """
        return len(self._names)

    @cached_method
    def gen(self, i=None):
        """
        Returns the `i^{th}` 'generator' (see the description in :meth:`.ngens`)
        of this infinite polynomial ring.

        EXAMPLES::

            sage: X = InfinitePolynomialRing(QQ)
            sage: x = X.gen()
            sage: x[1]
            x1
            sage: X.gen() is X.gen(0)
            True
            sage: XX = InfinitePolynomialRing(GF(5))
            sage: XX.gen(0) is XX.gen()
            True
        """
        if i > len(self._names):
            raise ValueError
        j = i if i is not None else 0
        res = InfinitePolynomialGen(self, self._names[j])
        if i is None:
            key = ((0,), ())
            if key in self._cache__gen:
                return self._cache__gen[key]
            else:
                self._cache__gen[key] = res
        return res

    def _ideal_class_(self):
        """
        Return :class:`SymmetricIdeals` (see there for further details).

        """
        import sage.rings.polynomial.symmetric_ideal
        return sage.rings.polynomial.symmetric_ideal.SymmetricIdeal

    def characteristic(self):
        """
        Return the characteristic of the base field.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(GF(25,'a'))
            sage: X
            Infinite polynomial ring in x, y over Finite Field in a of size 5^2
            sage: X.characteristic()
            5

        """
        return self._base.characteristic()

class InfinitePolynomialGen(SageObject):
    """
    This class provides the object which is responsible for returning
    variables in an infinite polynomial ring (implemented in
    :meth:`.__getitem__`).

    EXAMPLES::

        sage: X.<x,y> = InfinitePolynomialRing(RR)
        sage: x
        Generator for the x's in Infinite polynomial ring in x, y over Real Field with 53 bits of precision
        sage: x[5]
        x5
        sage: x == loads(dumps(x))
        True

    """

    def __init__(self, parent, name):
        """
        EXAPMLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: loads(dumps(x))
            Generator for the x's in Infinite polynomial ring in x over Rational Field

        """
        self._name = name
        self._parent = parent

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialGen
            sage: x2 = InfinitePolynomialGen(X, 'x')
            sage: x2 == x
            True

        """
        if not isinstance(other, InfinitePolynomialGen):
            return -1
        return cmp((self._name,self._parent),(other._name,other._parent))

    def _latex_(self):
        """
        EXAMPLES::

            sage: from sage.misc.latex import latex
            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: latex(x)
            x_{\ast}
            sage: latex(x[3])
            x_{3}
        """
        return self._name+'_{\\ast}'

    @cached_method
    def __getitem__(self, i):
        """
        Returns the the variable ``x[i]`` where ``x`` is this
        :class:`sage.rings.polynomial.infinite_polynomial_ring.InfinitePolynomialGen`,
        and i is a non-negative integer.

        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ)
            sage: x[1]
            x1

        """
        if i < 0:
            raise ValueError, "i (= %s) must be non-negative"%i
        P = self._parent
        from sage.rings.polynomial.infinite_polynomial_element import InfinitePolynomial
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        if hasattr(P,'_P'):
            if i <= P._max:
                return InfinitePolynomial(P, P._P(self._name+str(i)), is_good_poly=True)
            else:
                P._max = i
            #Calculate all of the new names needed
            names = [ [name+str(j) for name in P._names] for j in range(P._max+1)]
            names = reduce(operator.add, names)
            names.sort(cmp=P.varname_cmp,reverse=True)
            #Create the new polynomial ring
            P._P = PolynomialRing(P.base_ring(), names, order = P._order)
            ##Get the generators
            #P._pgens = P._P.gens()
            return InfinitePolynomial(P, P._P(self._name+str(i)), is_good_poly=True)
        return InfinitePolynomial(P, PolynomialRing(P.base_ring(), self._name+str(i), order=P._order).gen(), is_good_poly=True)

    def __repr__(self):
        """
        EXAPMLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ)
            sage: x
            Generator for the x's in Infinite polynomial ring in x, y over Rational Field

        """
        return "Generator for the %s's in %s"%(self._name, self._parent)

##############################################################
##  The dense implementation

class InfinitePolynomialRing_dense(InfinitePolynomialRing_sparse):
    """
    Dense implementation of Infinite Polynomial Rings

    Compared with :class:`~sage.rings.polynomial.infinite_polynomial_ring.InfinitePolynomialRing_sparse`,
    from which this class inherits, it keeps a polynomial ring that comprises all elements that have
    been created so far.
    """
    def __init__(self, R, names, order):
        """
        EXAMPLES::

            sage: X.<x> = InfinitePolynomialRing(QQ, implementation='sparse')
            sage: X == loads(dumps(X))
            True
        """
        #Generate the initial polynomial ring
        self._max = 0
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        InfinitePolynomialRing_sparse.__init__(self, R, names, order)
        #self._populate_coercion_lists_()
        VarList = [X+'0' for X in names]
        VarList.sort(cmp=self.varname_cmp, reverse=True)
        self._P = PolynomialRing(R, len(names), VarList)
        self._pgens = self._P.gens()

    def polynomial_ring(self):
        """
        Returns the underlying *finite* polynomial ring.

        .. note::

           This ring returned can change over time as more variables
           are used.

        EXAMPLES::

            sage: X.<x, y> = InfinitePolynomialRing(QQ)
            sage: X.polynomial_ring()
            Multivariate Polynomial Ring in y0, x0 over Rational Field
            sage: a = y[3]
            sage: X.polynomial_ring()
            Multivariate Polynomial Ring in y3, y2, y1, y0, x3, x2, x1, x0 over Rational Field

        """
        return self._P
