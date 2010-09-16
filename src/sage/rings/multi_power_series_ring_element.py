r"""
Multivariate Power Series

Construct and manipulate multivariate power series over a given commutative
ring.  Multivariate power series are implemented with total-degree
precision.

EXAMPLES:

Power series arithmetic, tracking precision::

    sage: R.<s,t> = PowerSeriesRing(ZZ); R
    Multivariate Power Series Ring in s, t over Integer Ring

    sage: f = 1 + s + 3*s^2; f
    1 + s + 3*s^2
    sage: g = t^2*s + 3*t^2*s^2 + R.O(5); g
    s*t^2 + 3*s^2*t^2 + O(s, t)^5
    sage: f = f.O(7); f
    1 + s + 3*s^2 + O(s, t)^7
    sage: f += s; f
    1 + 2*s + 3*s^2 + O(s, t)^7
    sage: f*g
    s*t^2 + 5*s^2*t^2 + O(s, t)^5
    sage: (f-1)*g
    2*s^2*t^2 + 9*s^3*t^2 + O(s, t)^6
    sage: f*g - g
    2*s^2*t^2 + O(s, t)^5

    sage: f*=s; f
    s + 2*s^2 + 3*s^3 + O(s, t)^8
    sage: f%2
    s + s^3 + O(s, t)^8
    sage: (f%2).parent()
    Multivariate Power Series Ring in s, t over Ring of integers modulo 2

Comparisons: As with univariate power series, comparison of f and g is done
up to the minimum precision of f and g.

::

    sage: f = 1 + t + s + s*t + R.O(3); f
    1 + s + t + s*t + O(s, t)^3
    sage: g = s^2 + 2*s^4 - s^5 + s^2*t^3 + R.O(6); g
    s^2 + 2*s^4 - s^5 + s^2*t^3 + O(s, t)^6
    sage: f == g
    False
    sage: g == g.add_bigoh(3)
    True
    sage: f < g
    False
    sage: f > g
    True

Calling::

    sage: f = s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + R.O(5); f
    s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
    sage: f(t,s)
    s*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5
    sage: f(t^2,s^2)
    s^2*t^2 + t^4 + s^2*t^4 + t^6 + 3*s^2*t^6 + 3*t^8 + O(s, t)^10

- Substitution is defined only for elements of positive valuation, unless f
  has infinite precision::

    sage: f(t^2,s^2+1)
    Traceback (most recent call last):
    ...
    TypeError: Substitution defined only for elements of positive valuation,
    unless self has infinite precision.

    sage: g = f.truncate()
    sage: g(t^2,s^2+1)
    t^2 + s^2*t^2 + 2*t^4 + s^2*t^4 + 4*t^6 + 3*s^2*t^6 + 3*t^8
    sage: g(t^2,(s^2+1).O(3))
    t^2 + s^2*t^2 + 2*t^4 + O(s, t)^5

- 0 has valuation +Infinity::

    sage: f(t^2,0)
    t^4 + t^6 + 3*t^8 + O(s, t)^10
    sage: f(t^2,s^2+s)
    s*t^2 + s^2*t^2 + t^4 + O(s, t)^5

- Substitution of power series with finite precision works too::

    sage: f(s.O(2),t)
    s^2 + s*t + O(s, t)^3
    sage: f(f,f)
    2*s^4 + 4*s^3*t + 2*s^2*t^2 + 4*s^5 + 8*s^4*t + 4*s^3*t^2 + 16*s^6 +
    34*s^5*t + 20*s^4*t^2 + 2*s^3*t^3 + O(s, t)^7
    sage: t(f,f)
    s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
    sage: t(0,f) == s(f,0)
    True

- subs works as expected::

    sage: r0 = -t^2 - s*t^3 - 2*t^6 + s^7 + s^5*t^2 + R.O(10)
    sage: r1 = s^4 - s*t^4 + s^6*t - 4*s^2*t^5 - 6*s^3*t^5 + R.O(10)
    sage: r2 = 2*s^3*t^2 - 2*s*t^4 - 2*s^3*t^4 + s*t^7 + R.O(10)
    sage: r0.subs({t:r2,s:r1})
    -4*s^6*t^4 + 8*s^4*t^6 - 4*s^2*t^8 + 8*s^6*t^6 - 8*s^4*t^8 - 4*s^4*t^9
    + 4*s^2*t^11 - 4*s^6*t^8 + O(s, t)^15
    sage: r0.subs({t:r2,s:r1}) == r0(r1,r2)
    True

- construct ring homomorphisms from one power series ring to another::

    sage: A.<a,b> = PowerSeriesRing(QQ)
    sage: X.<x,y> = PowerSeriesRing(QQ)

    sage: phi = Hom(A,X)([x,2*y]); phi
    Ring morphism:
      From: Multivariate Power Series Ring in a, b over Rational Field
      To:   Multivariate Power Series Ring in x, y over Rational Field
      Defn: a |--> x
            b |--> 2*y

    sage: phi(a+b+3*a*b^2 + A.O(5))
    x + 2*y + 12*x*y^2 + O(x, y)^5


Inversion::

    sage: h = 1 + s + t + s*t + s^2*t^2 + 3*s^4 + 3*s^3*t + R.O(5);
    sage: k = h^-1; k
    1 - s - t + s^2 + s*t + t^2 - s^3 - s^2*t - s*t^2 - t^3 - 2*s^4 -
    2*s^3*t + s*t^3 + t^4 + O(s, t)^5
    sage: h*k
    1 + O(s, t)^5

    sage: f = 1 - 5*s^29 - 5*s^28*t + 4*s^18*t^35 + \
    4*s^17*t^36 - s^45*t^25 - s^44*t^26 + s^7*t^83 + \
    s^6*t^84 + R.O(101)
    sage: h = 1/f; h
    1 + 5*s^29 + 5*s^28*t - 4*s^18*t^35 - 4*s^17*t^36 + 25*s^58 + 50*s^57*t
    + 25*s^56*t^2 + s^45*t^25 + s^44*t^26 - 40*s^47*t^35 - 80*s^46*t^36
    - 40*s^45*t^37 + 125*s^87 + 375*s^86*t + 375*s^85*t^2 + 125*s^84*t^3
    - s^7*t^83 - s^6*t^84 + 10*s^74*t^25 + 20*s^73*t^26 + 10*s^72*t^27
    + O(s, t)^101
    sage: h*f
    1 + O(s, t)^101



AUTHORS:

- Niles Johnson (07/2010): initial code

"""


#*****************************************************************************
#       Copyright (C) 2010 Niles Johnson <nilesj@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.power_series_ring_element import is_PowerSeries, PowerSeries

from sage.rings.polynomial.all import is_MPolynomialRing, is_PolynomialRing
from sage.rings.power_series_ring import is_PowerSeriesRing

from sage.rings.integer import Integer
from sage.rings.finite_rings.integer_mod_ring import Zmod

from sage.rings.infinity import infinity
from sage.structure.nonexact import Nonexact
from sage.structure.parent_gens import ParentWithGens
from sage.structure.category_object import CategoryObject


def is_MPowerSeries(f):
    """
        Return true if input is a multivariate power series.

    TESTS::

        sage: from sage.rings.power_series_ring_element import is_PowerSeries
        sage: from sage.rings.multi_power_series_ring_element import is_MPowerSeries
        sage: M = PowerSeriesRing(ZZ,4,'v');
        sage: is_PowerSeries(M.random_element(10))
        True
        sage: is_MPowerSeries(M.random_element(10))
        True
        sage: T.<v> = PowerSeriesRing(RR)
        sage: is_MPowerSeries(1 - v + v^2 +O(v^3))
        False
        sage: is_PowerSeries(1 - v + v^2 +O(v^3))
        True

        """

    return type(f) == MPowerSeries



class MPowerSeries(PowerSeries):
    ### methods from PowerSeries that we *don't* override:
    #
    # __hash__ : works just fine
    #
    # __reduce__ : don't really understand this
    #
    # is_sparse : works just fine
    #
    # is_dense : works just fine
    #
    # is_gen : works just fine
    #
    # base_extend : works just fine
    #
    # change_ring : works just fine
    #
    # _cmp_c_impl : don't understand this
    #
    # __copy__ : works just fine
    #
    # base_ring : works just fine
    #
    # common_prec : works just fine
    #
    # common_prec_c : seems fine
    #
    # _mul_prec : works just fine
    #
    # __nonzero__ : works just fine
    #
    """
    Multivariate power series; these are the elements of Multivariate Power
    Series Rings.

    INPUT

        - ``parent`` - A multivariate power series.

        - ``x`` - The element (default: 0).  This can be another
          MPowerSeries object, or an element of one of the following:

            - the background univariate power series ring

            - the foreground polynomial ring

            - a ring that coerces to one of the above two

        - ``prec`` - the precision (default: infinity)

        - ``is_gen`` - is this element one of the generators? (default: False)

        - ``check`` - needed by univariate power series class (default: False)


    EXAMPLES:

    Construct multivariate power series from generators::

        sage: S.<s,t> = PowerSeriesRing(ZZ)
        sage: f = s + 4*t + 3*s*t
        sage: f in S
        True
        sage: f = f.add_bigoh(4); f
        s + 4*t + 3*s*t + O(s, t)^4
        sage: g = 1 + s + t - s*t + S.O(5); g
        1 + s + t - s*t + O(s, t)^5


        sage: T = PowerSeriesRing(GF(3),5,'t'); T
        Multivariate Power Series Ring in t0, t1, t2, t3, t4 over Finite
        Field of size 3
        sage: t = T.gens()
        sage: w = t[0] - 2*t[1]*t[3] + 5*t[4]^3 - t[0]^3*t[2]^2; w
        t0 + t1*t3 - t4^3 - t0^3*t2^2
        sage: w = w.add_bigoh(5); w
        t0 + t1*t3 - t4^3 + O(t0, t1, t2, t3, t4)^5
        sage: w in T
        True

        sage: w = t[0] - 2*t[0]*t[2] + 5*t[4]^3 - t[0]^3*t[2]^2 + T.O(6)
        sage: w
        t0 + t0*t2 - t4^3 - t0^3*t2^2 + O(t0, t1, t2, t3, t4)^6


    Get random elements::

        sage: S.random_element(4) # random
        -2*t + t^2 - 12*s^3 + O(s, t)^4

        sage: T.random_element(10) # random
        -t1^2*t3^2*t4^2 + t1^5*t3^3*t4 + O(t0, t1, t2, t3, t4)^10


    Convert elements from polynomial rings::

        sage: R = PolynomialRing(ZZ,5,T.gens())
        sage: t = R.gens()
        sage: r = -t[2]*t[3] + t[3]^2 + t[4]^2
        sage: T(r)
        -t2*t3 + t3^2 + t4^2
        sage: r.parent()
        Multivariate Polynomial Ring in t0, t1, t2, t3, t4 over Integer Ring
        sage: r in T
        True

    """


    def __init__(self, parent, x=0, prec=infinity, is_gen=False, check=False):
        """
        input x can be an MPowerSeries, or an element of

          - the background univariate power series ring

          - the foreground polynomial ring

          - a ring that coerces to one of the above two


        TESTS::

            sage: S.<s,t> = PowerSeriesRing(ZZ)
            sage: f = s + 4*t + 3*s*t
            sage: f in S
            True
            sage: f = f.add_bigoh(4); f
            s + 4*t + 3*s*t + O(s, t)^4
            sage: g = 1 + s + t - s*t + S.O(5); g
            1 + s + t - s*t + O(s, t)^5

            sage: B.<s, t> = PowerSeriesRing(QQ)
            sage: C.<z> = PowerSeriesRing(QQ)
            sage: B(z)
            Traceback (most recent call last):
            ...
            TypeError: Cannot coerce input to polynomial ring.

            sage: D.<s> = PowerSeriesRing(QQ)
            sage: s.parent() is D
            True
            sage: B(s) in B
            True
            sage: d = D.random_element(20)
            sage: b = B(d) # test coercion from univariate power series ring
            sage: b in B
            True

        """
        PowerSeries.__init__(self, parent, prec, is_gen=is_gen)
        self._PowerSeries__is_gen = is_gen

        try:
            prec = min(prec, x.prec()) # use precision of input, if defined
        except AttributeError:
            pass


        # set the correct background value, depending on what type of input x is
        try:
            xparent = x.parent() # 'int' types have no parent
        except AttributeError:
            xparent = None

        # test whether x coerces to background univariate
        # power series ring of parent
        from sage.rings.multi_power_series_ring import is_MPowerSeriesRing
        if is_PowerSeriesRing(xparent) or is_MPowerSeriesRing(xparent):
            # x is either a multivariate or univariate power series
            #
            # test whether x coerces directly to designated parent
            if is_MPowerSeries(x):
                try:
                    self._bg_value = parent._bg_ps_ring(x._bg_value)
                except TypeError:
                    raise TypeError("Unable to coerce into background ring.")

            # test whether x coerces to background univariate
            # power series ring of parent
            elif xparent == parent._bg_ps_ring():
                self._bg_value = x
            elif parent._bg_ps_ring().has_coerce_map_from(xparent):
                # previous test may fail if precision or term orderings of
                # base rings do not match
                self._bg_value = parent._bg_ps_ring(x)
            else:
                # x is a univariate power series, but not from the
                # background power series ring
                #
                # convert x to a polynomial and send to background
                # ring of parent
                x = x.polynomial()
                self._bg_value = parent._send_to_bg(x).add_bigoh(prec)

        # test whether x coerces to underlying polynomial ring of parent
        elif is_PolynomialRing(xparent):
            self._bg_value = parent._send_to_bg(x).add_bigoh(prec)

        else:
            try:
                x = parent._poly_ring(x)
                #self._value = x
                self._bg_value = parent._send_to_bg(x).add_bigoh(prec)
            except (TypeError, AttributeError):
                raise TypeError("Input does not coerce to any of the expected rings.")

        self._go_to_fg = parent._send_to_fg
        self._prec = self._bg_value.prec()

        # self._parent is used a lot by the class PowerSeries
        self._parent = self.parent()

    def __reduce__(self):
        """
        EXAMPLES::

            sage: K.<s,t> = PowerSeriesRing(QQ)
            sage: f = 1 + t - s + s*t - s*t^3 + K.O(12)
            sage: loads(dumps(f)) == f
            True
        """

        from sage.rings.power_series_ring_element import make_element_from_parent_v0
        return make_element_from_parent_v0, (self._parent,self._bg_value,self._prec)

    def __call__(self, *x, **kwds):
        """
        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + R.O(5); f
            s^2 + s*t + s^3 + s^2*t + 3*s^4 + 3*s^3*t + O(s, t)^5
            sage: f(t,s)
            s*t + t^2 + s*t^2 + t^3 + 3*s*t^3 + 3*t^4 + O(s, t)^5

            sage: f(t,0)
            t^2 + t^3 + 3*t^4 + O(s, t)^5
            sage: f(t,2)
            Traceback (most recent call last):
            ...
            TypeError: Substitution defined only for elements of positive
            valuation, unless self has infinite precision.

            sage: f.truncate()(t,2)
            2*t + 3*t^2 + 7*t^3 + 3*t^4
        """
        if len(x) != self.parent().ngens():
            raise ValueError("Number of arguments does not match number of variables in parent.")

        sub_dict = {}
        valn_list =[]
        for i in range(len(x)):
            try:
                 xi = self.parent(x[i])
            except (AttributeError, TypeError):
                # Input does not coerce to parent ring of self
                # attempt formal substitution
                return self._subs_formal(*x,**kwds)
            if xi.valuation() == 0 and self.prec() is not infinity:
                raise TypeError("Substitution defined only for elements of positive valuation, unless self has infinite precision.")
            elif xi.valuation() > 0:
                sub_dict[self.parent()._poly_ring().gens()[i]] = xi.add_bigoh(xi.valuation()*self.prec())
                valn_list.append(xi.valuation())
            else:
                sub_dict[self.parent()._poly_ring().gens()[i]] = xi
        if self.prec() is infinity:
            newprec = infinity
        else:
            newprec = self.prec()*min(valn_list)
        return self._value().subs(sub_dict).add_bigoh(newprec)

    def _subs_formal(self,*x,**kwds):
        """
        Substitution of inputs as variables of self.  This is formal
        in the sense that the inputs do not need to be elements of
        same multivariate power series ring as self.  They can be any
        objects which support addition and multiplication with
        eachother and with the coefficients of self.  If self has
        finite precision, the inputs must also support an `add_bigoh`
        method.

        TESTS::

            sage: B.<s, t> = PowerSeriesRing(QQ)
            sage: C.<z> = PowerSeriesRing(QQ)
            sage: s(z,z)
            z

            sage: f = -2/33*s*t^2 - 1/5*t^5 - s^5*t + s^2*t^4
            sage: f(z,z) #indirect doctest
            -2/33*z^3 - 1/5*z^5
            sage: f(z,1) #indirect doctest
            -1/5 - 2/33*z + z^2 - z^5
            sage: RF = RealField(10)
            sage: f(z,RF(1)) #indirect doctest
            -0.20 - 0.061*z + 1.0*z^2 - 0.00*z^3 - 0.00*z^4 - 1.0*z^5

            sage: m = matrix(QQ,[[1,0,1],[0,2,1],[-1,0,0]])
            sage: m
            [ 1  0  1]
            [ 0  2  1]
            [-1  0  0]
            sage: f(m,m) #indirect doctest
            [     2/33         0       1/5]
            [   131/55 -1136/165    -24/11]
            [     -1/5         0   -23/165]
            sage: f(m,m) == -2/33*m^3 - 1/5*m^5 #indirect doctest
            True

            sage: f = f.add_bigoh(10)
            sage: f(z,z)
            -2/33*z^3 - 1/5*z^5 + O(z^10)
            sage: f(m,m)
            Traceback (most recent call last):
            ...
            AttributeError: 'sage.matrix.matrix_rational_dense.Matrix_rational_dense' object has no attribute 'add_bigoh'
        """
        from sage.misc.misc_c import prod

        if len(x) == 1 and isinstance(x[0], (list, tuple)):
            x = x[0]
        n = self.parent().ngens()
        if len(x) != n:
            raise ValueError("Input must be of correct length.")
        if n == 0:
            return self

        y = 0
        for (m,c) in self.dict().iteritems():
                y += c*prod([x[i]**m[i] for i in range(n) if m[i] != 0])
        if self.prec() == infinity:
            return y
        else:
            return y.add_bigoh(self.prec())

    def _value(self):
        """
        Return the value of self in the foreground polynomial ring.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(GF(5)); R
            Multivariate Power Series Ring in a, b, c over Finite Field of
            size 5
            sage: f = 1 + a + b - a*b + R.O(3); f
            1 + a + b - a*b + O(a, b, c)^3
            sage: f._value()
            1 + a + b - a*b
            sage: f._value().parent()
            Multivariate Polynomial Ring in a, b, c over Finite Field of size 5

        """
        return self._go_to_fg(self._bg_value)

    def _repr_(self):
        """
        Return string representation of self.

        EXAMPLES::

            sage: B.<s,t,v> = PowerSeriesRing(QQ)
            sage: e = 1 + s - s*t + t*v/2 - 2*s*t*v/8 + B.O(4)
            sage: e._repr_()
            '1 + s - s*t + 1/2*t*v - 1/4*s*t*v + O(s, t, v)^4'

        """
        if self._prec == infinity:
            return "%s" % self._value()
        return "%(val)s + O(%(gens)s)^%(prec)s" \
               %{'val':self._value(),
                 'gens':', '.join(str(g) for g in self.parent().gens()),
                 'prec':self._prec}

    def _latex_(self):
        """
        Return latex representation of this multivariate power series.

        EXAMPLES::
            sage: M = PowerSeriesRing(GF(5),3,'t'); M
            Multivariate Power Series Ring in t0, t1, t2 over Finite Field
            of size 5
            sage: t = M.gens()
            sage: f = -t[0]^4*t[1]^3*t[2]^4 - 2*t[0]*t[1]^4*t[2]^7 \
            + 2*t[1]*t[2]^12 + 2*t[0]^7*t[1]^5*t[2]^2 + M.O(15)
            sage: f
            -t0^4*t1^3*t2^4 - 2*t0*t1^4*t2^7 + 2*t1*t2^12 + 2*t0^7*t1^5*t2^2
            + O(t0, t1, t2)^15
            sage: f._latex_()
            '- t_{0}^{4} t_{1}^{3} t_{2}^{4} + 3 t_{0} t_{1}^{4} t_{2}^{7} +
            2 t_{1} t_{2}^{12} + 2 t_{0}^{7} t_{1}^{5} t_{2}^{2}
            + O(t0, t1, t2)^15'
        """
        if self._prec == infinity:
            return "%s" % self._value()
        return "%(val)s + O(%(gens)s)^%(prec)s" \
               %{'val':self._value()._latex_(),
                 'gens':', '.join(g._latex_() for g in self.parent().gens()),
                 'prec':self._prec}


    def _im_gens_(self, codomain, im_gens):
        """
        Returns the image of this series under the map that sends the
        generators to im_gens. This is used internally for computing
        homomorphisms.

        EXAMPLES::

            sage: A.<a,b> = PowerSeriesRing(QQ)
            sage: X.<x,y> = PowerSeriesRing(QQ)
            sage: phi = Hom(A,X)([x,2*y])
            sage: phi = Hom(A,X)([x,2*y]); phi
            Ring morphism:
              From: Multivariate Power Series Ring in a, b over Rational Field
              To:   Multivariate Power Series Ring in x, y over Rational Field
              Defn: a |--> x
                    b |--> 2*y
            sage: phi(a+b+3*a*b^2 + A.O(5))  # indirect doctest
            x + 2*y + 12*x*y^2 + O(x, y)^5
        """
        return codomain(self(*im_gens))

    def __getitem__(self,n):
        """
        Return summand of total degree n.

        TESTS::

            sage: R.<a,b> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(4)
            sage: f[0]
            1
            sage: f[2]
            -a*b
            sage: f[3]
            0
            sage: f[4]
            Traceback (most recent call last):
            ...
            IndexError: Cannot return terms of total degree greater than or
            equal to precision of self.
        """
        if n >= self.prec():
            raise IndexError("Cannot return terms of total degree greater than or equal to precision of self.")
        return self.parent(self._bg_value[n])

    def __invert__(self):
        """
        Return multiplicative inverse of this multivariate power series.

        Currently implemented only if constant coefficient is a unit in the
        base ring.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b - b*c - a*c + R.O(4)
            sage: 1/f
            1 - a - b + a^2 + 3*a*b + a*c + b^2 + b*c - a^3 - 5*a^2*b
            - 2*a^2*c - 5*a*b^2 - 4*a*b*c - b^3 - 2*b^2*c + O(a, b, c)^4
        """
        if self.valuation() == 0:
            return self.parent(self._bg_value.__invert__())
        else:
            raise NotImplementedError("Multiplicative inverse of multivariate power series currently implemented only if constant coefficient is a unit.")

    ## comparisons
    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(GF(5)); R
            Multivariate Power Series Ring in a, b, c over Finite Field of size 5
            sage: f = a + b + c + a^2*c
            sage: f == f^2
            False
            sage: f = f.truncate()
            sage: f == f.O(4)
            True

        Ordering is determined by underlying polynomial ring::

            sage: a > b
            True
            sage: a > a^2
            True
            sage: b > a^2
            True
            sage: (f^2).O(3)
            a^2 + 2*a*b + 2*a*c + b^2 + 2*b*c + c^2 + O(a, b, c)^3
            sage: f < f^2
            False
            sage: f > f^2
            True
            sage: f < 2*f
            True
        """
        return cmp(self._bg_value,other._bg_value)



    ## arithmetic
    def _add_(left, right):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0 + f1; g #indirect doctest
            -6*b*c^3 - 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 + a^6*b^2*c
            - 14*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == f0.polynomial() + f1.polynomial()
            True
        """
        f = left._bg_value + right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

    def _sub_(left, right):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0 - f1; g #indirect doctest
            6*b*c^3 + 4*a^2*b*c^2 - a^3*b*c^2 + a^2*b^2*c^4 - a^6*b^2*c
            - 10*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == f0.polynomial() - f1.polynomial()
            True
        """
        f = left._bg_value - right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

    def _mul_(left, right):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f0 = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: f1 = -6*b*c^3 - 4*a^2*b*c^2 + a^6*b^2*c - 2*a^3*b^3*c^3 + R.O(10)
            sage: g = f0*f1; g #indirect doctest
            6*a^3*b^2*c^5 + 4*a^5*b^2*c^4 - 6*a^2*b^3*c^7 - 4*a^4*b^3*c^6
            + 72*a^3*b^4*c^6 + O(a, b, c)^14
            sage: g in R
            True

        The power series product and polynomial product agree up to
        total degree < precision of g::

            sage: diff = g.polynomial() - f0.polynomial() * f1.polynomial()
            sage: all(S >= g.prec() for S in [sum(e) for e in diff.exponents()])
            True
        """
        f = left._bg_value * right._bg_value
        return MPowerSeries(left.parent(), f, prec=f.prec())

#    def _rmul_(self, c):
#        # multivariate power series rings are assumed to be commutative
#        return self._lmul_(c)

    def _lmul_(self, c):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: g = 3*f; g #indirect doctest
            -3*a^3*b*c^2 + 3*a^2*b^2*c^4 - 36*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == 3 * (f.polynomial())
            True
            sage: g = f*5; g #indirect doctest
            -5*a^3*b*c^2 + 5*a^2*b^2*c^4 - 60*a^3*b^3*c^3 + O(a, b, c)^10
            sage: g in R
            True
            sage: g.polynomial() == (f.polynomial()) * 5
            True
        """
        f = c * self._bg_value
        return MPowerSeries(self.parent(), f, prec=f.prec())

    def _div_(self, denom_r):
        """
        division by a unit works, but cancellation doesn't

        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = 1 + a + b - a*b + R.O(3)
            sage: g = 1/f; g #indirect doctest
            1 - a - b + a^2 + 3*a*b + b^2 + O(a, b, c)^3
            sage: g in R
            True
            sage: g = a/(a*f)
            Traceback (most recent call last):
            ...
            TypeError: denominator must be a unit

        """
        f = self._bg_value / denom_r._bg_value
        return MPowerSeries(self.parent(), f, prec=f.prec())

#    def _r_action_(self, c):
#        # multivariate power series rings are assumed to be commutative
#        return self._l_action_(c)

    def _l_action_(self, c):
        """
        Multivariate power series support multiplication by any ring for
        which there is a supported action on the base ring.

        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = 1 + t + s + s*t + R.O(3)
            sage: g = f._l_action_(1/2); g
            1/2 + 1/2*s + 1/2*t + 1/2*s*t + O(s, t)^3
            sage: g.parent()
            Multivariate Power Series Ring in s, t over Rational Field
            sage: g = (1/2)*f; g
            1/2 + 1/2*s + 1/2*t + 1/2*s*t + O(s, t)^3
            sage: g.parent()
            Multivariate Power Series Ring in s, t over Rational Field

            sage: K = NumberField(x-3,'a')
            sage: g = K.random_element()*f
            sage: g.parent()
            Multivariate Power Series Ring in s, t over Number Field in a with defining polynomial x - 3

        """
        try:
            f = c * self._bg_value
            if f.parent() == self.parent()._bg_ps_ring():
                return MPowerSeries(self.parent(), f, prec=f.prec())
            else:
                from sage.rings.all import PowerSeriesRing
                new_parent = PowerSeriesRing(f.base_ring().base_ring(), num_gens = f.base_ring().ngens(), names = f.base_ring().gens())
                return MPowerSeries(new_parent, f, prec=f.prec())
        except (TypeError, AttributeError):
            raise TypeError("Action not defined.")

    def __mod__(self, other):
        """
        TESTS::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ)
            sage: f = -a^3*b*c^2 + a^2*b^2*c^4 - 12*a^3*b^3*c^3 + R.O(10)
            sage: g = f % 2; g
            a^3*b*c^2 + a^2*b^2*c^4 + O(a, b, c)^10
            sage: g in R
            False
            sage: g in R.base_extend(Zmod(2))
            True
            sage: g.polynomial() == f.polynomial() % 2
            True
        """
        if isinstance(other,(int,Integer,long)):
            return self.change_ring(Zmod(other))
        raise NotImplementedError("Mod on multivariate power series ring elements not defined except modulo an integer.")



    def dict(self):
        """
        Return underlying dictionary with keys the exponents and values the
        coefficients of this power series.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t',sparse=True); M
            Sparse Multivariate Power Series Ring in t0, t1, t2, t3 over
            Rational Field

            sage: M.inject_variables()
            Defining t0, t1, t2, t3

            sage: m = 2/3*t0*t1^15*t3^48 \
            - t0^15*t1^21*t2^28*t3^5 \
            + 1/2*t0^12*t1^29*t2^46*t3^6 \
            - 1/4*t0^39*t1^5*t2^23*t3^30 \
            + M.O(100)
            sage: m.dict()
            {(1, 15, 0, 48): 2/3, (15, 21, 28, 5): -1, (12, 29, 46, 6): 1/2, (39, 5, 23, 30): -1/4}
        """
        out_dict = {}
        for j in self._bg_value.coefficients():
            out_dict.update(j.dict())
        return out_dict

    def polynomial(self):
        """
        Return underlying polynomial of self as an element of underlying
        multivariate polynomial ring.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t'); M
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational
            Field
            sage: t = M.gens()
            sage: f = 1/2*t[0]^3*t[1]^3*t[2]^2 + 2/3*t[0]*t[2]^6*t[3] \
            - t[0]^3*t[1]^3*t[3]^3 - 1/4*t[0]*t[1]*t[2]^7 + M.O(10)
            sage: f
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7 + O(t0, t1, t2, t3)^10

            sage: f.polynomial()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7

            sage: f.polynomial().parent()
            Multivariate Polynomial Ring in t0, t1, t2, t3 over Rational Field

        Contrast with truncate::

            sage: f.truncate()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3 - 1/4*t0*t1*t2^7
            sage: f.truncate().parent()
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field
        """
        return self._value()

    def variables(self):
        """
        Return tuple of variables occuring in self.

        EXAMPLES::

            sage: T = PowerSeriesRing(GF(3),5,'t'); T
            Multivariate Power Series Ring in t0, t1, t2, t3, t4 over
            Finite Field of size 3
            sage: t = T.gens()
            sage: w = t[0] - 2*t[0]*t[2] + 5*t[4]^3 - t[0]^3*t[2]^2 + T.O(6)
            sage: w
            t0 + t0*t2 - t4^3 - t0^3*t2^2 + O(t0, t1, t2, t3, t4)^6
            sage: w.variables()
            (t0, t2, t4)
        """
        return tuple(self.parent(v) for v in self._value().variables())

    def monomials(self):
        """
        Return a list of monomials of self.  These are the keys of the dict
        returned by 'coefficients'.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 1 + a + b - a*b - b*c - a*c + R.O(4)
            sage: f.monomials()
            [1, b*c, b, a, a*c, a*b]
            sage: f = 1 + 2*a + 7*b - 2*a*b - 4*b*c - 13*a*c + R.O(4)
            sage: f.monomials()
            [1, b*c, b, a, a*c, a*b]


        """
        return self.coefficients().keys()

    def coefficients(self):
        """
        Return a dict of monomials and coefficients.

        EXAMPLES::

            sage: R.<s,t> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in s, t over Integer Ring
            sage: f = 1 + t + s + s*t + R.O(3)
            sage: f.coefficients()
            {s*t: 1, 1: 1, s: 1, t: 1}
            sage: (f^2).coefficients()
            {t^2: 1, 1: 1, t: 2, s*t: 4, s^2: 1, s: 2}

            sage: g = f^2 + f - 2; g
            3*s + 3*t + s^2 + 5*s*t + t^2 + O(s, t)^3
            sage: cd = g.coefficients()
            sage: g2 = sum(k*v for (k,v) in cd.iteritems()); g2
            3*s + 3*t + s^2 + 5*s*t + t^2
            sage: g2 == g.truncate()
            True
        """
        if self.is_sparse():
            return self.dict()
        tmp = {}
        for j in self._bg_value.coefficients():
            for m in j.monomials():
                tmp[self.parent(m)]=j.monomial_coefficient(self.parent()._poly_ring(m))
        return tmp

    def constant_coefficient(self):
        """
        Return constant coefficient of self.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 3 + a + b - a*b - b*c - a*c + R.O(4)
            sage: f.constant_coefficient()
            3
            sage: f.constant_coefficient().parent()
            Integer Ring
        """
        return self.base_ring()(self._bg_value[0])

    def exponents(self):
        """
        Return a list of tuples which hold the exponents of each monomial of self.

        EXAMPLES::

            sage: H = QQ[['x,y']]
            sage: (x,y) = H.gens()
            sage: h = -y^2 - x*y^3 - 6/5*y^6 - x^7 + 2*x^5*y^2 + H.O(10)
            sage: h
            -y^2 - x*y^3 - 6/5*y^6 - x^7 + 2*x^5*y^2 + O(x, y)^10
            sage: h.exponents()
            [(0, 2), (1, 3), (0, 6), (7, 0), (5, 2)]
        """
        exp_list = []
        for m in self._bg_value.coefficients():
            exp_list += m.exponents()
        return exp_list

    def V(self, n):
        r"""
        If `f = \sum a_{m_0, \ldots, m_k} x_0^{m_0} \cdots x_k^{m_k}`,
        then this function returns
        `\sum a_{m_0, \ldots, m_k} x_0^{n m_0} \cdots x_n^{n m_k}`.

        The total-degree precision of the output is ``n`` times the precision of self.

        EXAMPLES::

            sage: H = QQ[['x,y,z']]
            sage: (x,y,z) = H.gens()
            sage: h = -x*y^4*z^7 - 1/4*y*z^12 + 1/2*x^7*y^5*z^2 \
            + 2/3*y^6*z^8 + H.O(15)
            sage: h.V(3)
            -x^3*y^12*z^21 - 1/4*y^3*z^36 + 1/2*x^21*y^15*z^6 + 2/3*y^18*z^24 + O(x, y, z)^45
        """
        cd = self.coefficients()
        Vs = sum(v*k**n for (k,v) in cd.iteritems())
        return Vs.add_bigoh(self.prec()*n)

    def prec(self):
        """
        Return precision of self.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b, c over Integer Ring
            sage: f = 3 + a + b - a*b - b*c - a*c + R.O(4)
            sage: f.prec()
            4
            sage: f.truncate().prec()
            +Infinity
        """
        return self._prec

    def add_bigoh(self,prec):
        """
        Return a multivariate power series of total precision prec obtained
        by truncating self at precision ``prec``.  This is the same as 'O'.

        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ); B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r.add_bigoh(4)
            1 + x^2 - x*y + O(x, y)^4
            sage: r.add_bigoh(2)
            1 + O(x, y)^2

        Note that this does not change self::

            sage: r
            1 + x^2 - x*y
        """
        return self.parent(self._bg_value.add_bigoh(prec))


    def O(self, prec):
        """
        Return a multivariate power series of total precision prec obtained
        by truncating self at precision ``prec``.  This is the same as
        'add_bigoh'.


        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ); B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r.O(4)
            1 + x^2 - x*y + O(x, y)^4
            sage: r.O(2)
            1 + O(x, y)^2

        Note that this does not change self::

            sage: r
            1 + x^2 - x*y
        """
        return self.parent(self._bg_value.O(prec))

    def truncate(self, prec=infinity):
        """
        Return infinite precision multivariate power series formed by
        truncating self at precision ``prec``.

        EXAMPLES::

            sage: M = PowerSeriesRing(QQ,4,'t'); M
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field
            sage: t = M.gens()
            sage: f = 1/2*t[0]^3*t[1]^3*t[2]^2 + 2/3*t[0]*t[2]^6*t[3] \
            - t[0]^3*t[1]^3*t[3]^3 - 1/4*t[0]*t[1]*t[2]^7 + M.O(10)
            sage: f
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7 + O(t0, t1, t2, t3)^10

            sage: f.truncate()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3
            - 1/4*t0*t1*t2^7
            sage: f.truncate().parent()
            Multivariate Power Series Ring in t0, t1, t2, t3 over Rational Field

        Contrast with polynomial::

            sage: f.polynomial()
            1/2*t0^3*t1^3*t2^2 + 2/3*t0*t2^6*t3 - t0^3*t1^3*t3^3 - 1/4*t0*t1*t2^7
            sage: f.polynomial().parent()
            Multivariate Polynomial Ring in t0, t1, t2, t3 over Rational Field
        """
        return self.parent((self.O(prec))._value())

    def valuation(self):
        """
        Return valuation of self.

        EXAMPLES::

            sage: R.<a,b> = PowerSeriesRing(GF(4949717)); R
            Multivariate Power Series Ring in a, b over Finite Field of
            size 4949717
            sage: f = a^2 + a*b + a^3 + R.O(9)
            sage: f.valuation()
            2
            sage: g = 1 + a + a^3
            sage: g.valuation()
            0
        """
        try:
            return self._bg_value.valuation()
        except (TypeError, AttributeError):
            if self._bg_value == 0:
                return infinity

            # at this stage, self is probably a non-zero
            # element of the base ring
            for a in range(len(self._bg_value.list())):
                if self._bg_value.list()[a] is not 0:
                    return a

    def is_nilpotent(self):
        """
        Return True if self is nilpotent.  This occurs if

            - self has finite precision and positive valuation

            - self is constant and nilpotent in base ring

        otherwise, return False.

        EXAMPLES::

            sage: R.<a,b,c> = PowerSeriesRing(Zmod(8)); R
            Multivariate Power Series Ring in a, b, c over Ring of integers
            modulo 8
            sage: f = a + b + c + a^2*c
            sage: f.is_nilpotent()
            False
            sage: f = f.O(4); f
            a + b + c + a^2*c + O(a, b, c)^4
            sage: f.is_nilpotent()
            True

            sage: g = R(2)
            sage: g.is_nilpotent()
            True
            sage: (g.O(4)).is_nilpotent()
            True

            sage: S = R.change_ring(QQ)
            sage: S(g).is_nilpotent()
            False
            sage: S(g.O(4)).is_nilpotent()
            False
        """
        if self.prec() < infinity and self.valuation() > 0:
            return True
        elif self == self.constant_coefficient() and \
           self.base_ring()(self.constant_coefficient()).is_nilpotent():
            return True
        else:
            return False


    def degree(self):
        """
        Return degree of underlying polynomial of self.

        EXAMPLES::

            sage: B.<x,y> = PowerSeriesRing(QQ)
            sage: B
            Multivariate Power Series Ring in x, y over Rational Field
            sage: r = 1 - x*y + x^2
            sage: r = r.add_bigoh(4); r
            1 + x^2 - x*y + O(x, y)^4
            sage: r.degree()
            2

        """
        return self._value().degree()

    def is_unit(self):
        """
        A multivariate power series is a unit if and only if its constant
        coefficient is a unit.

        EXAMPLES::

            sage: R.<a,b> = PowerSeriesRing(ZZ); R
            Multivariate Power Series Ring in a, b over Integer Ring
            sage: f = 2 + a^2 + a*b + a^3 + R.O(9)
            sage: f.is_unit()
            False
            sage: f.base_extend(QQ).is_unit()
            True
        """
        return self._bg_value[0].is_unit()

    ###
    ### the following could be implemented, but aren't
    ###

    def padded_list(self):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.padded_list()
            Traceback (most recent call last):
            ...
            NotImplementedError: padded_list
        """
        raise NotImplementedError("padded_list")

    def is_square(self):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.is_square()
            Traceback (most recent call last):
            ...
            NotImplementedError: is_square
        """
        raise NotImplementedError("is_square")

    def sqrt(self):
        """
        Method from univariate power series not yet implemented
        (An alias for ``square_root``.)

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.sqrt()
            Traceback (most recent call last):
            ...
            NotImplementedError: square_root
        """
        return self.square_root()

    def square_root(self):
        """
        Method from univariate power series not yet implemented.
        Depends on square root method for multivariate polynomials.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.square_root()
            Traceback (most recent call last):
            ...
            NotImplementedError: square_root
        """
        raise NotImplementedError("square_root")

    def derivative(self, v):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.derivative(a)
            Traceback (most recent call last):
            ...
            NotImplementedError: derivative
        """
        raise NotImplementedError("derivative")

    def ogf(self):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.ogf()
            Traceback (most recent call last):
            ...
            NotImplementedError: ogf
        """
        raise NotImplementedError("ogf")

    def egf(self):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.egf()
            Traceback (most recent call last):
            ...
            NotImplementedError: egf
        """
        raise NotImplementedError("egf")

    def _pari_(self):
        """
        Method from univariate power series not yet implemented

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f._pari_()
            Traceback (most recent call last):
            ...
            NotImplementedError: _pari_
        """
        raise NotImplementedError("_pari_")



    ###
    ### the following don't make sense for multivariable power series
    ###
    def list(self):
        """
        Doesn't make sense for multivariate power series.
        Multivariate polynomials don't have list of coefficients either.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.list()
            Traceback (most recent call last):
            ...
            NotImplementedError: Multivariate power series do not have list
            of coefficients; use 'coefficients' to get a dict of coefficients.
        """
        #return [self.parent(c) for c in self._bg_value.list()]
        raise NotImplementedError("Multivariate power series do not have list of coefficients; use 'coefficients' to get a dict of coefficients.")


    def variable(self):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.variable()
            Traceback (most recent call last):
            ...
            NotImplementedError: variable not defined for multivariate power
            series; use 'variables' instead.
        """
        raise NotImplementedError("variable not defined for multivariate power series; use 'variables' instead.")

    def shift(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.shift(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: shift not defined for multivariate power series.
        """
        raise NotImplementedError("shift not defined for multivariate power series.")

    def __lshift__(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.__lshift__(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: __lshift__ not defined for multivariate power series.
        """
        raise NotImplementedError("__lshift__ not defined for multivariate power series.")

    def __rshift__(self, n):
        """
        Doesn't make sense for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.__rshift__(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: __rshift__ not defined for multivariate power series.
        """
        raise NotImplementedError("__rshift__ not defined for multivariate power series.")

    def valuation_zero_part(self):
        """
        Doesn't make sense for multivariate power series;
        valuation zero with respect to which variable?

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.valuation_zero_part()
            Traceback (most recent call last):
            ...
            NotImplementedError: valuation_zero_part not defined for multivariate
            power series; perhaps 'constant_coefficient' is what you want.
        """
        raise NotImplementedError("valuation_zero_part not defined for multivariate power series; perhaps 'constant_coefficient' is what you want.")

    def solve_linear_de(self, prec = infinity, b=None, f0=None):
        """
        Not implemented for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.solve_linear_de()
            Traceback (most recent call last):
            ...
            NotImplementedError: solve_linear_de not defined for multivariate power series.
        """
        raise NotImplementedError("solve_linear_de not defined for multivariate power series.")

    def exp(self, prec=None):
        """
        Not implemented for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.exp()
            Traceback (most recent call last):
            ...
            NotImplementedError: exp not defined for multivariate power series.
        """
        raise NotImplementedError("exp not defined for multivariate power series.")

    def laurent_series(self):
        """
        Not implemented for multivariate power series.

        TESTS::

            sage: T.<a,b> = PowerSeriesRing(ZZ,2)
            sage: f = a + b + a*b + T.O(5)
            sage: f.laurent_series()
            Traceback (most recent call last):
            ...
            NotImplementedError: laurent_series not defined for multivariate power series.
        """
        raise NotImplementedError("laurent_series not defined for multivariate power series.")




