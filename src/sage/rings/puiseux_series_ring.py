r"""
Puiseux Series Ring
===================

The ring of Puiseux series.

This defines how to construct an element of the ring and
how other rings can coerce into this ring.

Classes
-------

.. autosummary::

    PuiseuxSeriesRing_generic
    PuiseuxSeriesRing_domain
    PuiseuxSeriesRing_field

Functions
---------

.. autosummary::

    PuiseuxSeriesRing
    is_PuiseuxSeriesRing

Examples
--------

See :mod:`sage.rings.puiseux_series_ring_element` for examples.

Contents
--------
"""
import weakref

from sage.all import parent
from sage.rings.ring import IntegralDomain, CommutativeRing, Field
from sage.rings.puiseux_series_ring_element import PuiseuxSeries
from sage.categories.complete_discrete_valuation import CompleteDiscreteValuationFields
from sage.categories.fields import Fields
from sage.rings.laurent_series_ring import (is_LaurentSeriesRing,
                                            LaurentSeriesRing)
from sage.rings.laurent_series_ring_element import LaurentSeries
from sage.rings.power_series_ring import is_PowerSeriesRing, PowerSeriesRing
from sage.rings.power_series_ring_element import PowerSeries

puiseux_series = {}


def PuiseuxSeriesRing(base_ring, name=None, names=None, default_prec=None,
                      sparse=False):
    """
    Rings of Puiseux series.

    EXAMPLES::

        sage: P = PuiseuxSeriesRing(QQ,'y')
        sage: y = P.gen()
        sage: f = y**(4/3)+y**(-5/6); f
        y^(-5/6) + y^(4/3)
        sage: f.add_bigoh(2)
        y^(-5/6) + y^(4/3) + O(y^2)
        sage: f.add_bigoh(1)
        y^(-5/6) + O(y)
    """
    if not names is None:
        name = names
    if name is None:
        raise TypeError('You must specify the name of the '
                        'indeterminate of the Puiseux series ring.')

    if default_prec is None:
        from sage.misc.defaults import series_precision
        default_prec = series_precision()

    global puiseux_series
    key = (base_ring, name, default_prec, sparse)
    if key in puiseux_series:
        x = puiseux_series[key]()
        if x is not None:
            return x

    if isinstance(base_ring, Field):
        R = PuiseuxSeriesRing_field(base_ring, name, default_prec, sparse)
    elif isinstance(base_ring, IntegralDomain):
        R = PuiseuxSeriesRing_domain(base_ring, name, default_prec, sparse)
    elif isinstance(base_ring, CommutativeRing):
        R = PuiseuxSeriesRing_generic(base_ring, name, default_prec, sparse)
    else:
        raise TypeError("base_ring must be a commutative ring")
    puiseux_series[key] = weakref.ref(R)
    return R


def is_PuiseuxSeriesRing(x):
    """
    Return whether this ring is a Puiseux series ring.

    EXAMPLES::

        sage: from sage.rings.puiseux_series_ring import is_PuiseuxSeriesRing
        sage: is_PuiseuxSeriesRing(QQ)
        False
        sage: P = PuiseuxSeriesRing(QQ,'y')
        sage: is_PuiseuxSeriesRing(P)
        True
    """
    return isinstance(x, PuiseuxSeriesRing_generic)


class PuiseuxSeriesRing_generic(CommutativeRing):
    def __init__(self, base_ring, name=None, default_prec=None, sparse=False,
                 category=None):
        """
        Generic class for Puiseux series rings.
        """
        CommutativeRing.__init__(self, base_ring, names=name,
                                 category=getattr(self, '_default_category',
                                                  Fields()))

        # If self is R(( x^(1/e) )) then the corresponding Laurent series
        # ring will be R(( x ))
        self._laurent_series_ring = LaurentSeriesRing(base_ring, name=name,
                                                      default_prec=default_prec, sparse=sparse)

    def base_extend(self, R):
        """
        Extend the coefficients.

        INPUT:

        R -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.base_extend(QQ)
            Puiseux Series Ring in y over Rational Field
        """
        if R.has_coerce_map_from(self.base_ring()):
            return self.change_ring(R)
        else:
            raise TypeError("no valid base extension defined")

    def change_ring(self, R):
        """
        Return a Puiseux series ring over another ring.

        INPUT:

        R -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.change_ring(QQ)
            Puiseux Series Ring in y over Rational Field
        """
        return PuiseuxSeriesRing(R, self.variable_name(),
                                 default_prec=self.default_prec(),
                                 sparse=self.is_sparse())

    def is_sparse(self):
        """
        Return whether self is sparse.

        INPUT:

        R -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_sparse()
            False
        """
        return self.laurent_series_ring().is_sparse()

    def is_field(self, proof=True):
        """
        Return whether self is a field.

        A Puiseux series ring is a field if and only
        its base ring is a field.

        INPUT:

        R -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_field()
            False
            sage: A.change_ring(QQ).is_field()
            True
        """
        return self.base_ring().is_field()

    def is_dense(self):
        """
        Return whether self is dense.

        INPUT:

        R -- a ring

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(ZZ, 'y')
            sage: A.is_dense()
            True
        """
        return self.laurent_series_ring().is_dense()

    def __reduce__(self):
        """
        Used for pickling.

        EXAMPLES::

            sage: P = PuiseuxSeriesRing(QQ,'y')
            sage: P.__reduce__()
            (<class 'sage.rings.puiseux_series_ring.PuiseuxSeriesRing_field_with_category'>,
             (Rational Field, 'y', 20, False))
        """
        return self.__class__, (self.base_ring(), self.variable_name(),
                                self.default_prec(), self.is_sparse())

    def _repr_(self):
        """
        String representation.

        EXAMPLES::

            sage: PuiseuxSeriesRing(AA, 'y')  # indirect doctest
            Puiseux Series Ring in y over Algebraic Real Field
        """
        s = "Puiseux Series Ring in %s over %s" % (self.variable_name(),
                                                   self.base_ring())
        if self.is_sparse():
            s = 'Sparse ' + s
        return s

    Element = PuiseuxSeries

    def _element_constructor_(self, x, e=1):
        r"""
        Construct a Puiseux series from `x`.

        INPUT:

        - ``x`` -- an object that can be converted into a Puiseux series in (x-a)

        - ``a`` -- (default: 0) the series is in powers of (var - a)

        - ``e`` -- (default: 1) the ramification index of the series

        EXAMPLES::

            sage: P = PuiseuxSeriesRing(QQ,'y')
            sage: y = P.gen()
            sage: P([1,3,5,7])
            1 + 3*y + 5*y^2 + 7*y^3
            sage: P(33/14)
            33/14

            sage: Q = PowerSeriesRing(QQ,'y')
            sage: z = Q([1,2,4,5]).O(6); z
            1 + 2*y + 4*y^2 + 5*y^3 + O(y^6)
            sage: P(z) + y**(1/2)
            1 + y^(1/2) + 2*y + 4*y^2 + 5*y^3 + O(y^6)

            sage: Q = LaurentSeriesRing(QQ,'y')
            sage: z = Q([3,2,1,2]).add_bigoh(5); z
            3 + 2*y + y^2 + 2*y^3 + O(y^5)
            sage: P(z) + y**(1/2)
            3 + y^(1/2) + 2*y + y^2 + 2*y^3 + O(y^5)
        """
        P = parent(x)

        # 1. x is a Puiseux series belonging to this ring
        if isinstance(x, self.element_class) and P is self:
            return x
        # 2. x is a Puiseux series but not an element of this ring. the laurent
        #    part should be coercible to the laurent series ring of self
        elif isinstance(x, self.element_class):
            l = self.laurent_series_ring()(x.laurent_part)
            e = x.ramification_index
        # 3. x is a member of the base ring then convert x to a laurent series
        #    and set the ramification index of the Puiseux series to 1.
        elif P is self.base_ring():
            l = self.laurent_series_ring()(x)
            e = 1
        # 4. x is a Laurent or power series with the same base ring
        elif (isinstance(x, (LaurentSeries, PowerSeries))
              and P is self.base_ring()):
            l = self.laurent_series_ring()(x)
        # 5. everything else: try to coerce to laurent series ring
        else:
            l = self.laurent_series_ring()(x)
            e = 1

        return self.element_class(self, l, e=e)

    def _coerce_map_from_(self, P):
        r"""
        Return a coercion map from `P` to `self`, or `True`, or `None`.

        The following rings admit a coercion map to the Puiseux series ring
        `A((x-a)^(1/e))`:

        - any ring that admits a coercion map to `A`

        - any Laurent series ring, power series ring, or polynomial ring in the
          variable `(x-a)` over a ring admitting a coercion map to `A`

        - any Puiseux series ring with the same center `a` and ramification
          index equal to a multiple of `self`'s ramification index. For
          example, Puiseux series in (x-a)^(1/2) can be interpreted as Puiseux
          series in (x-a)^(1/4).

        """
        # any ring that has a coercion map to A
        A = self.base_ring()
        if A is P:
            return True
        f = A.coerce_map_from(P)
        if f is not None:
            return self.coerce_map_from(A) * f

        # Laurent series rings, power series rings, and polynomial rings with
        # the same variable name and the base rings are coercible
        if ((is_PuiseuxSeriesRing(P) or is_LaurentSeriesRing(P) or
             is_PowerSeriesRing(P)) and
                P.variable_name() == self.variable_name() and
                A.has_coerce_map_from(P.base_ring())):
            return True

        # # other Puiseux series rings with the same variable name and
        # # center. Puiseux series rings with difference ramification indices are
        # # coercible to each other.
        # if (is_PuiseuxSeriesRing(P) and
        #     P.variable_name() == self.variable_name()):
        #     return True

    def gen(self, n=0):
        """
        Return the generator of self.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.gen()
            z
        """
        if n != 0:
            raise IndexError("Generator n not defined.")
        try:
            return self.__generator
        except AttributeError:
            #l = self.laurent_series_ring()([0,1])
            self.__generator = PuiseuxSeries(self, [0, 1], e=1)
            return self.__generator

    def ngens(self):
        """
        Return the number of generators of self, namely 1

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.ngens()
            1
        """
        return 1

    def laurent_series_ring(self):
        """
        Return the underlying Laurent series ring.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.laurent_series_ring()
            Laurent Series Ring in z over Algebraic Real Field
        """
        return self._laurent_series_ring

    def default_prec(self):
        """
        Return the default precision of self.

        EXAMPLES::

            sage: A = PuiseuxSeriesRing(AA, 'z')
            sage: A.default_prec()
            20
        """
        return self.laurent_series_ring().default_prec()


class PuiseuxSeriesRing_domain(PuiseuxSeriesRing_generic,
                               IntegralDomain):
    def __init__(self, base_ring, name=None, default_prec=None, sparse=False):
        PuiseuxSeriesRing_generic.__init__(self, base_ring, name,
                                           default_prec, sparse)


class PuiseuxSeriesRing_field(PuiseuxSeriesRing_generic, Field):
    _default_category = CompleteDiscreteValuationFields()

    def __init__(self, base_ring, name=None, default_prec=None, sparse=False):
        PuiseuxSeriesRing_generic.__init__(self, base_ring, name,
                                           default_prec, sparse)
