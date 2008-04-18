"""
Ring of Laurent Polynomials.

If R is a commutative ring, then the ring of laurent polynomials in n variables over R is $R[x_1^{\pm 1}, x_2^{\pm 1}, \ldots, x_n^{\pm 1}].$
We implement it as a quotient ring $R[x_1, xx_1, x_2, xx_2, \ldots, x_n, xx_n] / (x_1*xx_1 - 1, x_2*xx_2 - 1, \ldots, x_n, xx_n - 1)$.

AUTHORS:
  -- David Roe (2008-2-23)
"""

#################################################################################
#       Copyright (C) 2008 David Roe <roed@math.harvard.edu>, William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import weakref
from sage.structure.parent_gens import normalize_names
from sage.rings.quotient_ring import QuotientRing_generic
from sage.structure.element import is_Element
from sage.rings.ring import is_Ring
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import _single_variate as _single_variate_poly
from sage.rings.polynomial.polynomial_ring_constructor import _multi_variate as _multi_variate_poly
from sage.misc.latex import latex
from sage.rings.polynomial.laurent_polynomial import LaurentPolynomial_quotient, LaurentPolynomial_mpair
from sage.rings.polynomial.polydict import ETuple
from sage.rings.ring import CommutativeRing
from sage.structure.parent_gens import ParentWithGens

def is_LaurentPolynomialRing(R):
    return isinstance(R, LaurentPolynomialRing_generic)

def LaurentPolynomialRing(base_ring, arg1=None, arg2=None, sparse = False, order='degrevlex', names = None, name=None):
    r"""
    Return the globally unique univariate or multivariate laurent polynomial
    ring with given properties and variable name or names.

    There are four ways to call the polynomial ring constructor:
          1. LaurentPolynomialRing(base_ring, name,    sparse=False)
          2. LaurentPolynomialRing(base_ring, names,   order='degrevlex')
          3. LaurentPolynomialRing(base_ring, name, n, order='degrevlex')
          4. LaurentPolynomialRing(base_ring, n, name, order='degrevlex')

    The optional arguments sparse and order *must* be explicitly
    named, and the other arguments must be given positionally.

    INPUT:
         base_ring -- a commutative ring
         name -- a string
         names -- a list or tuple of names, or a comma separated string
         n -- a positive integer
         sparse -- bool (default: False), whether or not elements are sparse
         order -- string or TermOrder, e.g.,
                 'degrevlex' (default) -- degree reverse lexicographic
                 'lex'  -- lexicographic
                 'deglex' -- degree lexicographic
                 TermOrder('deglex',3) + TermOrder('deglex',3) -- block ordering

    OUTPUT:
        LaurentPolynomialRing(base_ring, name, sparse=False) returns a univariate
        laurent polynomial ring; all other input formats return a multivariate
        laurent polynomial ring.

    UNIQUENESS and IMMUTABILITY: In SAGE there is exactly one
    single-variate laurent polynomial ring over each base ring in each choice
    of variable and sparsenes.  There is also exactly one multivariate
    laurent polynomial ring over each base ring for each choice of names of
    variables and term order.  The names of the generators can only
    be temporarily changed after the ring has been created.  Do this
    using the localvars context:

        EXAMPLES of VARIABLE NAME CONTEXT:
            sage: R.<x,y> = LaurentPolynomialRing(QQ,2); R
            Multivariate Laurent Polynomial Ring in x, y over Rational Field
            sage: f = x^2 - 2*y^-2

        You can't just globally change the names of those variables.
        This is because objects all over SAGE could have pointers to
        that polynomial ring.
            sage: R._assign_names(['z','w'])
            Traceback (most recent call last):
            ...
            ValueError: variable names cannot be changed after object creation.

        However, you can very easily change the names within a "with" block:
            sage: with localvars(R, ['z','w']):
            ...     print f
            ...
            z^2 - 2*w^-2

        After the with block the names revert to what they were before.
            sage: print f
            x^2 - 2*y^-2

    EXAMPLES:
    1. LaurentPolynomialRing(base_ring, name,    sparse=False):
        sage: PolynomialRing(QQ, 'w')
        Univariate Laurent Polynomial Ring in w over Rational Field

    Use the diamond brackets notation to make the variable
    ready for use after you define the ring:
        sage: R.<w> = LaurentPolynomialRing(QQ)
        sage: (1 + w)^3
        w^3 + 3*w^2 + 3*w + 1

    You must specify a name:
        sage: LaurentPolynomialRing(QQ)
        Traceback (most recent call last):
        ...
        TypeError: You must specify the names of the variables.

        sage: R.<abc> = LaurentPolynomialRing(QQ, sparse=True); R
        Sparse Univariate Laurent Polynomial Ring in abc over Rational Field

        sage: R.<w> = LaurentPolynomialRing(PolynomialRing(GF(7),'k')); R
        Univariate Laurent Polynomial Ring in w over Univariate Polynomial Ring in k over Finite Field of size 7

    Rings with different variables are different:
        sage: LaurentPolynomialRing(QQ, 'x') == LaurentPolynomialRing(QQ, 'y')
        False

    2. LaurentPolynomialRing(base_ring, names,   order='degrevlex'):
        sage: R = LaurentPolynomialRing(QQ, 'a,b,c'); R
        Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

        sage: S = LaurentPolynomialRing(QQ, ['a','b','c']); S
        Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

        sage: T = LaurentPolynomialRing(QQ, ('a','b','c')); T
        Multivariate Laurent Polynomial Ring in a, b, c over Rational Field

    All three rings are identical.
        sage: (R is S) and  (S is T)
        True

    There is a unique laurent polynomial ring with each term order:
        sage: R = LaurentPolynomialRing(QQ, 'x,y,z', order='degrevlex'); R
        Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
        sage: S = LaurentPolynomialRing(QQ, 'x,y,z', order='invlex'); S
        Multivariate Laurent Polynomial Ring in x, y, z over Rational Field
        sage: S is LaurentPolynomialRing(QQ, 'x,y,z', order='invlex')
        True
        sage: R == S
        False


    3. LaurentPolynomialRing(base_ring, name, n, order='degrevlex'):

    If you specify a single name as a string and a number of
    variables, then variables labeled with numbers are created.
        sage: LaurentPolynomialRing(QQ, 'x', 10)
        Multivariate Laurent Polynomial Ring in x0, x1, x2, x3, x4, x5, x6, x7, x8, x9 over Rational Field

        sage: LaurentPolynomialRing(GF(7), 'y', 5)
        Multivariate Laurent Polynomial Ring in y0, y1, y2, y3, y4 over Finite Field of size 7

        sage: LaurentPolynomialRing(QQ, 'y', 3, sparse=True)
        Multivariate Laurent Polynomial Ring in y0, y1, y2 over Rational Field

    You can call \code{injvar}, which is a convenient shortcut for \code{inject_variables()}.
        sage: R = LaurentPolynomialRing(GF(7),15,'w'); R
        Multivariate Laurent Polynomial Ring in w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14 over Finite Field of size 7
        sage: R.injvar()
        Defining w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14
        sage: (w0 + 2*w8 + w13)^2
        w0^2 - 3*w0*w8 - 3*w8^2 + 2*w0*w13 - 3*w8*w13 + w13^2
    """
    if is_Element(arg1) and not isinstance(arg1, (int, long, Integer)):
        arg1 = repr(arg1)
    if is_Element(arg2) and not isinstance(arg2, (int, long, Integer)):
        arg2 = repr(arg2)

    if isinstance(arg1, (int, long, Integer)):
        arg1, arg2 = arg2, arg1

    if not names is None:
        arg1 = names
    elif not name is None:
        arg1 = name

    if not is_Ring(base_ring):
        raise TypeError, 'base_ring must be a ring'

    if arg1 is None:
        raise TypeError, "You must specify the names of the variables."

    R = None
    if isinstance(arg1, (list, tuple)):
        arg1 = [str(x) for x in arg1]
    if isinstance(arg2, (list, tuple)):
        arg2 = [str(x) for x in arg2]
    if isinstance(arg2, (int, long, Integer)):
        # 3. LaurentPolynomialRing(base_ring, names, n, order='degrevlex'):
        if not isinstance(arg1, (list, tuple, str)):
            raise TypeError, "You *must* specify the names of the variables."
        n = int(arg2)
        names = arg1
        R = _multi_variate(base_ring, names, n, sparse, order)

    elif isinstance(arg1, str) or (isinstance(arg1, (list,tuple)) and len(arg1) == 1) and isinstance(arg1[0], str):
        if isinstance(arg1, (list,tuple)):
            arg1 = arg1[0]
        if not ',' in arg1:
            # 1. LaurentPolynomialRing(base_ring, name, sparse=False):
            if not arg2 is None:
                raise TypeError, "if second arguments is a string with no commas, then there must be no other non-optional arguments"
            name = arg1
            R = _single_variate(base_ring, name, sparse)
        else:
            # 2-4. LaurentPolynomialRing(base_ring, names, order='degrevlex'):
            if not arg2 is None:
                raise TypeError, "invalid input to LaurentPolynomialRing function; please see the docstring for that function"
            names = arg1.split(',')
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)
    elif isinstance(arg1, (list, tuple)):
            # LaurentPolynomialRing(base_ring, names (list or tuple), order='degrevlex'):
            names = arg1
            n = len(names)
            R = _multi_variate(base_ring, names, n, sparse, order)

    if arg1 is None and arg2 is None:
        raise TypeError, "you *must* specify the indeterminates (as not None)."
    if R is None:
        raise TypeError, "invalid input (%s, %s, %s) to PolynomialRing function; please see the docstring for that function"%(base_ring, arg1, arg2)

    return R

_cache = {}
def _get_from_cache(key):
    try:
        if _cache.has_key(key):
            return _cache[key]()
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)
    return None

def _save_in_cache(key, R):
    try:
        _cache[key] = weakref.ref(R)
    except TypeError, msg:
        raise TypeError, 'key = %s\n%s'%(key,msg)

def _single_variate(base_ring, name, sparse):
    ####################################################################################################################################
    ##### This should later get moved to an actual single variate implementation with valuation tracking, but I don't want to right now.
    ####################################################################################################################################
    # We need to come up with a name for the inverse that is easy to search for in a string *and* doesn't overlap with the name that we already have.
    # For now, I'm going to use a name mangling with checking method.
    name = normalize_names(1, name)
    key = (base_ring, name, sparse)
    P = _get_from_cache(key)
    if P is not None:
        return P
    prepend_string = "qk"
    while True:
        if prepend_string in name:
            prepend_string += 'k'
        else:
            break
    R = _single_variate_poly(base_ring, name, sparse)
    RR = _multi_variate_poly(base_ring, (name, prepend_string + name), 2, False, 'degrevlex')
    P = LaurentPolynomialRing_quotient(base_ring, 1, R, RR, prepend_string, (name,))
    _save_in_cache(key, P)
    return P

def _multi_variate(base_ring, names, n, sparse, order):
    # We need to come up with names for the inverses that are easy to search for in a string *and* don't overlap with names that we already have.
    # For now, I'm going to use a name mangling with checking method.
    names = normalize_names(n, names)

    from term_order import TermOrder
    order = TermOrder(order, n)

    if isinstance(names, list):
        names = tuple(names)
    elif isinstance(names, str):
        if ',' in names:
            names = tuple(names.split(','))

    key = (base_ring, names, n, sparse, order)
    P = _get_from_cache(key)
    if P is not None:
        return P
    prepend_string = "qk"
    while True:
        for a in names:
            if prepend_string in a:
                prepend_string += 'k'
                break
        else:
            break
    R = _multi_variate_poly(base_ring, names, n, sparse, order)
    #RR = _multi_variate_poly(base_ring, list(names) + [prepend_string + a for a in names], 2*n, sparse, order)
    #P =  LaurentPolynomialRing_quotient(base_ring, n, R, RR, prepend_string, names)
    P = LaurentPolynomialRing_mpair(R, prepend_string, names)
    _save_in_cache(key, P)
    return P

class LaurentPolynomialRing_generic(CommutativeRing, ParentWithGens):
    def __init__(self, R, prepend_string, names):
        self._n = R.ngens()
        self._R = R
        self._prepend_string = prepend_string
        ParentWithGens.__init__(self, R.base_ring(), names)

    def _repr_(self):
        if self._n == 1:
            return "Univariate Laurent Polynomial Ring in %s over %s"%(self._R.variable_name(), self._R.base_ring())
        else:
            return "Multivariate Laurent Polynomial Ring in %s over %s"%(", ".join(self._R.variable_names()), self._R.base_ring())

    def ngens(self):
        return self._n

    def gen(self, i=0):
        if i < 0 or i >= self._n:
            raise ValueError, "Generator not defined"
        return self(self._R.gen(i))

    def is_integral_domain(self):
        return self.base_ring().is_integral_domain()

    def is_noetherian(self):
        raise NotImplementedError

    def construction(self):
        from sage.categories.pushout import LaurentPolynomialFunctor
        vars = self.variable_names()
        if len(vars) == 1:
            return LaurentPolynomialFunctor(vars[0], False), self.base_ring()
        else:
            return LaurentPolynomialFunctor(vars[-1], True), LaurentPolynomialRing(self.base_ring(), vars[:-1])

    def completion(self, p, prec=20, extras=None):
        raise NotImplementedError

    def remove_var(self, var):
        vars = list(self.variable_names())
        vars.remove(str(var))
        return LaurentPolynomialRing(self.base_ring(), vars)

    def coerce_map_from_impl(self, R):
        if R is self._R:
            from sage.categories.morphism import FormalCoercionMorphism
            return FormalCoercionMorphism(R, self)
        else:
            f = self._R.coerce_map_from(R)
            if f is not None:
                from sage.categories.homset import Hom
                from sage.categories.morphism import CallMorphism
                return CallMorphism(Hom(self._R, self)) * f
        return None

    def __cmp__(left, right):
        c = cmp(type(left), type(right))
        if c == 0:
            c = cmp(left._R, right._R)
        return c

    def _latex_(self):
        vars = [a + '^{\pm 1}' for a in self.latex_variable_names()]
        vars = vars.join(', ')
        return "%s[%s]"%(latex(self.base_ring()), vars)

    def _ideal_class_(self):
        # One may eventually want ideals in these guys.
        raise NotImplementedError

    def ideal(self):
        raise NotImplementedError

    def _is_valid_homomorphism_(self, codomain, im_gens):
        if not codomain.has_coerce_map_from(self.base_ring()):
            # we need that elements of the base ring
            # canonically coerce into codomain.
            return False
        for a in im_gens:
            # in addition, the image of each generator must be invertible.
            if not a.is_unit():
                return False
        return True

    def term_order(self):
        return self._R.term_order()

    def is_finite(self):
        return False

    def is_field(self):
        return False

    def polynomial_ring(self):
        return self._R

    def characteristic(self):
        return self._base_ring().characteristic()

    def krull_dimension(self):
        raise NotImplementedError

    def random_element(self, low_degree = -2, high_degree = 2, terms = 5, choose_degree=False,*args, **kwds):
        raise NotImplementedError

    def is_exact(self):
        return self.base_ring().is_exact()

    def change_ring(self, base_ring=None, names=None, sparse=False, order=None):
        if base_ring is None:
            base_ring = self.base_ring()
        if names is None:
            names = self.variable_names()
        if order is None:
            order = self._RR.term_order()
        if self._n == 1:
            return LaurentPolynomialRing(base_ring, names[0], sparse = sparse)
        else:
            return LaurentPolynomialRing(base_ring, self._n, names, order = order)

class LaurentPolynomialRing_quotient(LaurentPolynomialRing_generic, QuotientRing_generic):
    def __init__(self, base_ring, n, R, RR, prepend_string, names):
        if n <= 0:
            raise ValueError, "n must be positive"
        if not base_ring.is_integral_domain():
            raise ValueError, "base ring must be an integral domain"
        I = RR.ideal([RR.gen(i) * RR.gen(i + n) - 1 for i in range(n)])
        LaurentPolynomialRing_generic.__init__(self, R, prepend_string, names)
        QuotientRing_generic.__init__(self, RR, I, names)
        self._RR = RR

    def _representing_quotient_ring(self):
        return self._RR

    def __call__(self, x):
        return LaurentPolynomial_quotient(self, x)

    def _to_long_exp(self, E):
        if self._n == 1:
            L = [E]
        else:
            L = list(E)
        for m in range(self._n):
            if L[m] < 0:
                L.append(-L[m])
                L[m] = 0
            else:
                L.append(0)
        return ETuple(L)

    def _to_short_exp(self, E):
        if self._n == 1:
            return E[0] - E[1]
        else:
            return ETuple([E[i] - E[i+self._n] for i in range(self._n)])

class LaurentPolynomialRing_mpair(LaurentPolynomialRing_generic):
    def __init__(self, R, prepend_string, names):
        if R.ngens() <= 0:
            raise ValueError, "n must be positive"
        if not R.base_ring().is_integral_domain():
            raise ValueError, "base ring must be an integral domain"
        LaurentPolynomialRing_generic.__init__(self, R, prepend_string, names)

    def __call__(self, x):
        return LaurentPolynomial_mpair(self, x)

