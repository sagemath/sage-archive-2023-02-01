import sage.misc.misc as misc

include "sage/ext/stdsage.pxi"
from sage.rings.integer cimport Integer

from sage.misc.derivative import multi_derivative

def is_MPolynomial(x):
    return isinstance(x, MPolynomial)

cdef class MPolynomial(CommutativeRingElement):

    ####################
    # Some standard conversions
    ####################
    def __int__(self):
        if self.degree() == 0:
            return int(self.constant_coefficient())
        else:
            raise TypeError

    def __long__(self):
        if self.degree() == 0:
            return long(self.constant_coefficient())
        else:
            raise TypeError

    def __float__(self):
        if self.degree() == 0:
            return float(self.constant_coefficient())
        else:
            raise TypeError

    def _mpfr_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_mpfr_field_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _complex_double_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _real_double_(self, R):
        if self.degree() == 0:
            return R(self.constant_coefficient())
        else:
            raise TypeError

    def _rational_(self):
        if self.degree() == 0:
            from sage.rings.rational import Rational
            return Rational(repr(self))
        else:
            raise TypeError

    def _integer_(self):
        if self.degree() == 0:
            from sage.rings.integer import Integer
            return Integer(repr(self))
        else:
            raise TypeError

    def _polynomial_(self, R):
        var = R.variable_name()
        if var in self._parent.variable_names():
            return R(self.polynomial(self._parent(var)))
        else:
            return R([self])

    def coefficients(self):
        """
        Return the nonzero coefficients of this polynomial in a list.
        The returned list is decreasingly ordered by the term ordering
        of self.parent(), i.e. the list of coefficients matches the list
        of monomials returned by self.monomials().

        EXAMPLES:
            sage: R.<x,y,z> = MPolynomialRing(QQ,3,order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = MPolynomialRing(QQ,3,order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

            # Test the same stuff with ZZ -- different implementation
            sage: R.<x,y,z> = MPolynomialRing(ZZ,3,order='degrevlex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [23, 6, 1]
            sage: R.<x,y,z> = MPolynomialRing(ZZ,3,order='lex')
            sage: f=23*x^6*y^7 + x^3*y+6*x^7*z
            sage: f.coefficients()
            [6, 23, 1]

        AUTHOR:
            -- didier deshommes
        """
        degs = self.exponents()
        d = self.dict()
        return  [ d[i] for i in degs ]

    def truncate(self, var, n):
        """
        Returns a new multivariate polynomial obtained from self by
        deleting all terms that involve the given variable to a power
        at least n.
        """
        cdef int ind
        R = self.parent()
        G = R.gens()
        Z = list(G)
        try:
            ind = Z.index(var)
        except ValueError:
            raise ValueError, "var must be one of the generators of the parent polynomial ring."
        d = self.dict()
        return R(dict([(k, c) for k, c in d.iteritems() if k[ind] < n]))

    def _fast_float_(self, *vars):
        """
        Returns a quickly-evaluating function on floats.

        EXAMPLE:
            sage: K.<x,y,z> = QQ[]
            sage: f = (x+2*y+3*z^2)^2 + 42
            sage: f(1, 10, 100)
            901260483
            sage: ff = f._fast_float_()
            sage: ff(0, 0, 1)
            51.0
            sage: ff(0, 1, 0)
            46.0
            sage: ff(1, 10, 100)
            901260483.0
            sage: ff_swapped = f._fast_float_('z', 'y', 'x')
            sage: ff_swapped(100, 10, 1)
            901260483.0
            sage: ff_extra = f._fast_float_('x', 'A', 'y', 'B', 'z', 'C')
            sage: ff_extra(1, 7, 10, 13, 100, 19)
            901260483.0

        Currently, we use a fairly unoptimized method that evaluates one
        monomial at a time, with no sharing of repeated computations and
        with useless additions of 0 and multiplications by 1:
            sage: list(ff)
            ['push 0.0', 'push 12.0', 'load 1', 'load 2', 'dup', 'mul', 'mul', 'mul', 'add', 'push 4.0', 'load 0', 'load 1', 'mul', 'mul', 'add', 'push 42.0', 'add', 'push 1.0', 'load 0', 'dup', 'mul', 'mul', 'add', 'push 9.0', 'load 2', 'dup', 'mul', 'dup', 'mul', 'mul', 'add', 'push 6.0', 'load 0', 'load 2', 'dup', 'mul', 'mul', 'mul', 'add', 'push 4.0', 'load 1', 'dup', 'mul', 'mul', 'add']

        TESTS:
            sage: from sage.ext.fast_eval import fast_float
            sage: list(fast_float(K(0)))
            ['push 0.0']
            sage: list(fast_float(K(17)))
            ['push 0.0', 'push 17.0', 'add']
            sage: list(fast_float(y))
            ['push 0.0', 'push 1.0', 'load 1', 'mul', 'add']
        """
        from sage.ext.fast_eval import fast_float_arg, fast_float_constant
        my_vars = self.parent().variable_names()
        vars = list(vars)
        if len(vars) == 0:
            indices = range(len(my_vars))
        else:
            indices = [vars.index(v) for v in my_vars]
        x = [fast_float_arg(i) for i in indices]

        n = len(x)
        expr = fast_float_constant(0)
        for (m,c) in self.dict().iteritems():
            monom = misc.mul([ x[i]**m[i] for i in range(n) if m[i] != 0], fast_float_constant(c))
            expr = expr + monom
        return expr


    def derivative(self, *args):
        r"""
        The formal derivative of this polynomial, with respect to
        variables supplied in args.

        Multiple variables and iteration counts may be supplied; see
        documentation for the global derivative() function for more details.

        SEE ALSO:
            self._derivative()

        EXAMPLES:

        Polynomials implemented via Singular:
            sage: R.<x, y> = PolynomialRing(FiniteField(5))
            sage: f = x^3*y^5 + x^7*y
            sage: type(f)
            <type 'sage.rings.polynomial.multi_polynomial_libsingular.MPolynomial_libsingular'>
            sage: f.derivative(x)
            2*x^6*y - 2*x^2*y^5
            sage: f.derivative(y)
            x^7

        Generic multivariate polynomials:
            sage: R.<t> = PowerSeriesRing(QQ)
            sage: S.<x, y> = PolynomialRing(R)
            sage: f = (t^2 + O(t^3))*x^2*y^3 + (37*t^4 + O(t^5))*x^3
            sage: type(f)
            <class 'sage.rings.polynomial.multi_polynomial_element.MPolynomial_polydict'>
            sage: f.derivative(x)   # with respect to x
            (2*t^2 + O(t^3))*x*y^3 + (111*t^4 + O(t^5))*x^2
            sage: f.derivative(y)   # with respect to y
            (3*t^2 + O(t^3))*x^2*y^2
            sage: f.derivative(t)   # with respect to t (recurses into base ring)
            (2*t + O(t^2))*x^2*y^3 + (148*t^3 + O(t^4))*x^3
            sage: f.derivative(x, y) # with respect to x and then y
            (6*t^2 + O(t^3))*x*y^2
            sage: f.derivative(y, 3) # with respect to y three times
            (6*t^2 + O(t^3))*x^2
            sage: f.derivative()    # can't figure out the variable
            Traceback (most recent call last):
            ...
            ValueError: must specify which variable to differentiate with respect to

        Polynomials over the symbolic ring (just for fun....):
            sage: x = var("x")
            sage: S.<u, v> = PolynomialRing(SR)
            sage: f = u*v*x
            sage: f.derivative(x) == u*v
            True
            sage: f.derivative(u) == v*x
            True
        """
        return multi_derivative(self, args)


    def polynomial(self, var):
        """
        Let var be one of the variables of the parent of self.  This
        returns self viewed as a unvariate polynomial in var over the
        polynomial ring generated by all the other variables of the parent.

        EXAMPLES:
            sage: R.<x,w,z> = QQ[]
            sage: f = x^3 + 3*w*x + w^5 + (17*w^3)*x + z^5
            sage: f.polynomial(x)
            x^3 + (17*w^3 + 3*w)*x + w^5 + z^5
            sage: parent(f.polynomial(x))
            Univariate Polynomial Ring in x over Multivariate Polynomial Ring in w, z over Rational Field

            sage: f.polynomial(w)
            w^5 + 17*x*w^3 + 3*x*w + z^5 + x^3
            sage: f.polynomial(z)
            z^5 + w^5 + 17*x*w^3 + x^3 + 3*x*w
            sage: R.<x,w,z,k> = ZZ[]
            sage: f = x^3 + 3*w*x + w^5 + (17*w^3)*x + z^5 +x*w*z*k + 5
            sage: f.polynomial(x)
            x^3 + (17*w^3 + w*z*k + 3*w)*x + w^5 + z^5 + 5
            sage: f.polynomial(w)
            w^5 + 17*x*w^3 + (x*z*k + 3*x)*w + z^5 + x^3 + 5
            sage: f.polynomial(z)
            z^5 + x*w*k*z + w^5 + 17*x*w^3 + x^3 + 3*x*w + 5
            sage: f.polynomial(k)
            x*w*z*k + w^5 + z^5 + 17*x*w^3 + x^3 + 3*x*w + 5
        """
        cdef int ind
        R = self.parent()
        G = R.gens()
        Z = list(G)
        try:
            ind = Z.index(var)
        except ValueError:
            raise ValueError, "var must be one of the generators of the parent polynomial ring."

        if R.ngens() <= 1:
            return self.univariate_polynomial()

        other_vars = Z
        del other_vars[ind]

        # Make polynomial ring over all variables except var.
        S = R.base_ring()[tuple(other_vars)]
        ring = S[var]
        if not self:
            return ring(0)

        d = self.degree(var)
        B = ring.base_ring()
        w = dict([(remove_from_tuple(e, ind), val) for e, val in self.dict().iteritems() if not e[ind]])
        v = [B(w)]  # coefficients that don't involve var
        z = var
        for i in range(1,d+1):
            c = self.coefficient(z).dict()
            w = dict([(remove_from_tuple(e, ind), val) for e, val in c.iteritems()])
            v.append(B(w))
            z *= var
        return ring(v)

    def _mpoly_dict_recursive(self, vars=None, base_ring=None):
        """
        Return a dict of coefficent entries suitable for construction of a MPolynomial_polydict
        with the given variables.

        EXAMPLES:
            sage: R = Integers(10)['x,y,z']['t,s']
            sage: t,s = R.gens()
            sage: x,y,z = R.base_ring().gens()
            sage: (x+y+2*z*s+3*t)._mpoly_dict_recursive(['z','t','s'])
            {(1, 0, 1): 2, (0, 1, 0): 3, (0, 0, 0): x + y}

        TESTS:
            sage: R = Qp(7)['x,y,z,t,p']; S = ZZ['x,z,t']['p']
            sage: R(S.0)
            p
            sage: R = QQ['x,y,z,t,p']; S = ZZ['x']['y,z,t']['p']
            sage: z = S.base_ring().gen(1)
            sage: R(z)
            z
            sage: R = QQ['x,y,z,t,p']; S = ZZ['x']['y,z,t']['p']
            sage: z = S.base_ring().gen(1); p = S.0; x = S.base_ring().base_ring().gen()
            sage: R(z+p)
            z + p
            sage: R = Qp(7)['x,y,z,p']; S = ZZ['x']['y,z,t']['p'] # shouldn't work, but should throw a better error
            sage: R(S.0)
            p
        """
        from polydict import ETuple
        if not self:
            return {}

        if vars is None:
            vars = self.parent().variable_names_recursive()
        vars = list(vars)
        my_vars = self.parent().variable_names()
        if vars == list(my_vars):
            return self.dict()
        elif not my_vars[-1] in vars:
            x = base_ring(self) if base_ring is not None else self
#            print "vars", vars, type(vars)
            const_ix = ETuple((0,)*len(vars))
            return { const_ix: x }
        elif not set(my_vars).issubset(set(vars)):
            # we need to split it up
            return self.polynomial(self.parent().gen(len(my_vars)-1))._mpoly_dict_recursive(vars, base_ring)
        else:
            D = {}
            prev_vars = vars[:vars.index(my_vars[0])]
            var_range = range(len(my_vars))
            if len(prev_vars) > 0:
                mapping = [vars.index(v) - len(prev_vars) for v in my_vars]
                tmp = [0] * (len(vars) - len(prev_vars))
                try:
                    for ix,a in self.dict().iteritems():
                        for k in var_range:
                            tmp[mapping[k]] = ix[k]
                        postfix = ETuple(tmp)
                        mpoly = a._mpoly_dict_recursive(prev_vars, base_ring)
                        for prefix,b in mpoly.iteritems():
                            D[prefix+postfix] = b
                    return D

                except AttributeError:
                    pass

            if base_ring is self.base_ring():
                base_ring = None

            mapping = [vars.index(v) for v in my_vars]
            tmp = [0] * len(vars)
            for ix,a in self.dict().iteritems():
                for k in var_range:
                    tmp[mapping[k]] = ix[k]
                if base_ring is not None:
                    a = base_ring(a)
                D[ETuple(tmp)] = a
            return D

    cdef long _hash_c(self):
        """
        This hash incorporates the variable name in an effort to respect the obvious inclusions
        into multi-variable polynomial rings.

        The tuple algorithm is borrowed from http://effbot.org/zone/python-hash.htm.

        EXAMPLES:
            sage: T.<y>=QQ[]
            sage: R.<x>=ZZ[]
            sage: S.<x,y>=ZZ[]
            sage: hash(S.0)==hash(R.0)  # respect inclusions into mpoly rings (with matching base rings)
            True
            sage: hash(S.1)==hash(T.0)  # respect inclusions into mpoly rings (with unmatched base rings)
            True
            sage: hash(S(12))==hash(12)  # respect inclusions of the integers into an mpoly ring
            True
            sage: # the point is to make for more flexible dictionary look ups
            sage: d={S.0:12}
            sage: d[R.0]
            12
            sage: # or, more to the point, make subs in fraction field elements work
            sage: f=x/y
            sage: f.subs({x:1})
            1/y
        """
        cdef long result = 0 # store it in a c-int and just let the overflowing additions wrap
        cdef long result_mon
        var_name_hash = [hash(v) for v in self._parent.variable_names()]
        cdef long c_hash
        for m,c in self.dict().iteritems():
            #  I'm assuming (incorrectly) that hashes of zero indicate that the element is 0.
            # This assumption is not true, but I think it is true enough for the purposes and it
            # it allows us to write fast code that omits terms with 0 coefficients.  This is
            # important if we want to maintain the '==' relationship with sparse polys.
            c_hash = hash(c)
            if c_hash != 0: # this is always going to be true, because we are sparse (correct?)
                # Hash (self[i], gen_a, exp_a, gen_b, exp_b, gen_c, exp_c, ...) as a tuple according to the algorithm.
                # I omit gen,exp pairs where the exponent is zero.
                result_mon = c_hash
                for p in m.nonzero_positions():
                    result_mon = (1000003 * result_mon) ^ var_name_hash[p]
                    result_mon = (1000003 * result_mon) ^ m[p]
                result += result_mon
        if result == -1:
            return -2
        return result

    # you may have to replicate this boilerplate code in derived classes if you override
    # __richcmp__.  The python documentation at  http://docs.python.org/api/type-structs.html
    # explains how __richcmp__, __hash__, and __cmp__ are tied together.
    def __hash__(self):
        return self._hash_c()

    def args(self):
        r"""
        Returns the named of the arguments of \code{self}, in the
        order they are accepted from call.

        EXAMPLES:
            sage: R.<x,y> = ZZ[]
            sage: x.args()
            (x, y)
        """
        return self._parent.gens()

    def homogenize(self, var='h'):
        r"""
        Return \code{self} if \code{self} is homogeneous.  Otherwise
        return a homogenized polynomial for \code{self}. If a string
        is given, return a polynomial in one more variable named after
        the strig such that setting that variable equal to 1 yields
        self. This variable is added to the end of the variables. If a
        variable in \code{self.parent()} is given, this variable is
        used to homogenize the polynomial. If an integer is given, the
        variable with this index is used for homogenization.

        INPUT:
            var -- either a variable name, variable index or a
                   variable (default: 'h').

        OUTPUT:
            a multivariate polynomial

        EXAMPLES:
            sage: P.<x,y> = PolynomialRing(QQ,2)
            sage: f = x^2 + y + 1 + 5*x*y^10
            sage: g = f.homogenize('z'); g
            5*x*y^10 + x^2*z^9 + y*z^10 + z^11
            sage: g.parent()
            Multivariate Polynomial Ring in x, y, z over Rational Field

            sage: f.homogenize(x)
            2*x^11 + x^10*y + 5*x*y^10

            sage: f.homogenize(0)
            2*x^11 + x^10*y + 5*x*y^10

            sage: x, y = Zmod(3)['x', 'y'].gens()
            sage: (x + x^2).homogenize(y)
            x^2 + x*y

            sage: x, y = Zmod(3)['x', 'y'].gens()
            sage: (x + x^2).homogenize(y).parent()
            Multivariate Polynomial Ring in x, y over Ring of integers modulo 3

            sage: x, y = GF(3)['x', 'y'].gens()
            sage: (x + x^2).homogenize(y)
            x^2 + x*y

            sage: x, y = GF(3)['x', 'y'].gens()
            sage: (x + x^2).homogenize(y).parent()
            Multivariate Polynomial Ring in x, y over Finite Field of size 3

        """
        P = self.parent()

        if self.is_homogeneous():
            return self

        if PY_TYPE_CHECK(var, basestring):
            V = list(P.variable_names())
            try:
                i = V.index(var)
                return self._homogenize(i)
            except ValueError:
                P = P.__class__(P.base_ring(), len(V)+1, V + [var], order=P.term_order())
                return P(self)._homogenize(len(V))

        elif PY_TYPE_CHECK(var, MPolynomial) and \
             ((<MPolynomial>var)._parent is P or (<MPolynomial>var)._parent == P):
            V = list(P.gens())
            try:
                i = V.index(var)
                return self._homogenize(i)
            except ValueError:
                P = P.change_ring(names=P.variable_names() + [str(var)])
                return P(self)._homogenize(len(V))

        elif PY_TYPE_CHECK(var, int) or PY_TYPE_CHECK(var, Integer):
            if 0 <= var < P.ngens():
                return self._homogenize(var)
            else:
                raise TypeError, "Variable index %d must be < parent(self).ngens()."%var
        else:
            raise TypeError, "Parameter var must be either a variable, a string or an integer."

    def is_homogeneous(self):
        """
        Return True if self is a homogeneous polynomial.

        EXAMPLES:
            sage: x, y = MPolynomialRing(RationalField(), 2, names=['x', 'y']).gens()
            sage: (x+y).is_homogeneous()
            True
            sage: (x.parent()(0)).is_homogeneous()
            True
            sage: (x+y^2).is_homogeneous()
            False
            sage: (x^2 + y^2).is_homogeneous()
            True
            sage: (x^2 + y^2*x).is_homogeneous()
            False
            sage: (x^2*y + y^2*x).is_homogeneous()
            True
        """
        M = self.monomials()
        d = M.pop().degree()
        for m in M:
            if m.degree() != d:
                return False
        else:
            return True

    def __mod__(self, other):
        """
        EXAMPLE:
            sage: R.<x,y> = PolynomialRing(QQ)
            sage: f = (x^2*y + 2*x - 3)
            sage: g = (x + 1)*f
            sage: g % f
            0

            sage: (g+1) % f
            1

            sage: M = x*y
            sage: N = x^2*y^3
            sage: M.divides(N)
            True
        """
        q,r = self.quo_rem(other)
        return r

    def change_ring(self, R):
        """
        Return a copy of this polynomial but with coefficients in R,
        if at all possible.

        INPUT:
            R -- a ring

        EXAMPLE:
            sage: R.<x,y> = QQ[]
            sage: f = x^3 + 3/5*y + 1
            sage: f.change_ring(GF(7))
            x^3 + 2*y + 1
        """
        P = self._parent
        P = P.change_ring(R)
        return P(self)

cdef remove_from_tuple(e, int ind):
    w = list(e)
    del w[ind]
    return tuple(w)

