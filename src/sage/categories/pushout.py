"""
Coercion via Construction Functors
"""
import six
from sage.misc.lazy_import import lazy_import
from functor import Functor, IdentityFunctor_generic

lazy_import('sage.categories.commutative_additive_groups', 'CommutativeAdditiveGroups')
lazy_import('sage.categories.commutative_rings', 'CommutativeRings')
lazy_import('sage.categories.groups', 'Groups')
lazy_import('sage.categories.objects', 'Objects')
lazy_import('sage.categories.rings', 'Rings')

lazy_import('sage.structure.parent', 'CoercionException')

# TODO, think through the rankings, and override pushout where necessary.

class ConstructionFunctor(Functor):
    """
    Base class for construction functors.

    A construction functor is a functorial algebraic construction,
    such as the construction of a matrix ring over a given ring
    or the fraction field of a given ring.

    In addition to the class :class:`~sage.categories.functor.Functor`,
    construction functors provide rules for combining and merging
    constructions. This is an important part of Sage's coercion model,
    namely the pushout of two constructions: When a polynomial ``p`` in
    a variable ``x`` with integer coefficients is added to a rational
    number ``q``, then Sage finds that the parents ``ZZ['x']`` and
    ``QQ`` are obtained from ``ZZ`` by applying a polynomial ring
    construction respectively the fraction field construction. Each
    construction functor has an attribute ``rank``, and the rank of
    the polynomial ring construction is higher than the rank of the
    fraction field construction. This means that the pushout of ``QQ``
    and ``ZZ['x']``, and thus a common parent in which ``p`` and ``q``
    can be added, is ``QQ['x']``, since the construction functor with
    a lower rank is applied first.

    ::

        sage: F1, R = QQ.construction()
        sage: F1
        FractionField
        sage: R
        Integer Ring
        sage: F2, R = (ZZ['x']).construction()
        sage: F2
        Poly[x]
        sage: R
        Integer Ring
        sage: F3 = F2.pushout(F1)
        sage: F3
        Poly[x](FractionField(...))
        sage: F3(R)
        Univariate Polynomial Ring in x over Rational Field
        sage: from sage.categories.pushout import pushout
        sage: P.<x> = ZZ[]
        sage: pushout(QQ,P)
        Univariate Polynomial Ring in x over Rational Field
        sage: ((x+1) + 1/2).parent()
        Univariate Polynomial Ring in x over Rational Field

    When composing two construction functors, they are sometimes
    merged into one, as is the case in the Quotient construction::

        sage: Q15, R = (ZZ.quo(15*ZZ)).construction()
        sage: Q15
        QuotientFunctor
        sage: Q35, R = (ZZ.quo(35*ZZ)).construction()
        sage: Q35
        QuotientFunctor
        sage: Q15.merge(Q35)
        QuotientFunctor
        sage: Q15.merge(Q35)(ZZ)
        Ring of integers modulo 5

    Functors can not only be applied to objects, but also to morphisms in the
    respective categories. For example::

        sage: P.<x,y> = ZZ[]
        sage: F = P.construction()[0]; F
        MPoly[x,y]
        sage: A.<a,b> = GF(5)[]
        sage: f = A.hom([a+b,a-b],A)
        sage: F(A)
        Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Finite Field of size 5
        sage: F(f)
        Ring endomorphism of Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Finite Field of size 5
          Defn: Induced from base ring by
                Ring endomorphism of Multivariate Polynomial Ring in a, b over Finite Field of size 5
                  Defn: a |--> a + b
                        b |--> a - b
        sage: F(f)(F(A)(x)*a)
        (a + b)*x

    """
    def __mul__(self, other):
        """
        Compose ``self`` and ``other`` to a composite construction
        functor, unless one of them is the identity.

        NOTE:

        The product is in functorial notation, i.e., when applying the
        product to an object, the second factor is applied first.

        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: F*P
            FractionField(Poly[t](...))
            sage: P*F
            Poly[t](FractionField(...))
            sage: (F*P)(ZZ)
            Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            sage: I*P is P
            True
            sage: F*I is F
            True

        """
        if not isinstance(self, ConstructionFunctor) and not isinstance(other, ConstructionFunctor):
            raise CoercionException("Non-constructive product")
        if isinstance(other,IdentityConstructionFunctor):
            return self
        if isinstance(self,IdentityConstructionFunctor):
            return other
        return CompositeConstructionFunctor(other, self)

    def pushout(self, other):
        """
        Composition of two construction functors, ordered by their ranks.

        NOTE:

        - This method seems not to be used in the coercion model.

        - By default, the functor with smaller rank is applied first.

        TESTS::

            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: F.pushout(P)
            Poly[t](FractionField(...))
            sage: P.pushout(F)
            Poly[t](FractionField(...))

        """
        if self.rank > other.rank:
            return self * other
        else:
            return other * self

    def __cmp__(self, other):
        """
        Equality here means that they are mathematically equivalent, though they may have
        specific implementation data. This method will usually be overloaded in subclasses.
        by default, only the types of the functors are compared. Also see the \code{merge}
        function.

        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: I == F        # indirect doctest
            False
            sage: I == I        # indirect doctest
            True

        """
        return cmp(type(self), type(other))

    def _repr_(self):
        """
        NOTE:

        By default, it returns the name of the construction functor's class.
        Usually, this method will be overloaded.

        TEST::

            sage: F = QQ.construction()[0]
            sage: F                  # indirect doctest
            FractionField
            sage: Q = ZZ.quo(2).construction()[0]
            sage: Q                  # indirect doctest
            QuotientFunctor

        """
        s = str(type(self))
        import re
        return re.sub("<.*'.*\.([^.]*)'>", "\\1", s)

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return None.

        NOTE:

        The default is to merge only if the two functors coincide. But this
        may be overloaded for subclasses, such as the quotient functor.

        EXAMPLES::

            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: F.merge(F)
            FractionField
            sage: F.merge(P)
            sage: P.merge(F)
            sage: P.merge(P)
            Poly[t]

        """
        if self == other:
            return self
        else:
            return None

    def commutes(self, other):
        """
        Determine whether ``self`` commutes with another construction functor.

        NOTE:

        By default, ``False`` is returned in all cases (even if the two
        functors are the same, since in this case :meth:`merge` will apply
        anyway). So far there is no construction functor that overloads
        this method. Anyway, this method only becomes relevant if two
        construction functors have the same rank.

        EXAMPLES::

            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: F.commutes(P)
            False
            sage: P.commutes(F)
            False
            sage: F.commutes(F)
            False

        """
        return False

    def expand(self):
        """
        Decompose ``self`` into a list of construction functors.

        NOTE:

        The default is to return the list only containing ``self``.

        EXAMPLE::

            sage: F = QQ.construction()[0]
            sage: F.expand()
            [FractionField]
            sage: Q = ZZ.quo(2).construction()[0]
            sage: Q.expand()
            [QuotientFunctor]
            sage: P = ZZ['t'].construction()[0]
            sage: FP = F*P
            sage: FP.expand()
            [FractionField, Poly[t]]

        """
        return [self]

    # See the pushout() function below for explanation.
    coercion_reversed = False

    def common_base(self, other_functor, self_bases, other_bases):
        r"""
        This function is called by :func:`pushout` when no common parent
        is found in the construction tower.

        .. NOTE::

            The main use is for multivariate construction functors,
            which use this function to implement recursion for
            :func:`pushout`.

        INPUT:

        - ``other_functor`` -- a construction functor.

        - ``self_bases`` -- the arguments passed to this functor.

        - ``other_bases`` -- the arguments passed to the functor
          ``other_functor``.

        OUTPUT:

        Nothing, since a
        :class:`~sage.structure.coerce_exceptions.CoercionException`
        is raised.

        .. NOTE::

            Overload this function in derived class, see
            e.e. :class:`MultivariateConstructionFunctor`.

        TESTS::

            sage: from sage.categories.pushout import pushout
            sage: pushout(QQ, cartesian_product([ZZ]))  # indirect doctest
            Traceback (most recent call last):
            ...
            CoercionException: No common base ("join") found for
            FractionField(Integer Ring) and The cartesian_product functorial construction(Integer Ring).
        """
        self._raise_common_base_exception_(
            other_functor, self_bases, other_bases)

    def _raise_common_base_exception_(self, other_functor,
                                      self_bases, other_bases,
                                      reason=None):
        r"""
        Raise a coercion exception.

        INPUT:

        - ``other_functor`` -- a functor.

        - ``self_bases`` -- the arguments passed to this functor.

        - ``other_bases`` -- the arguments passed to the functor
          ``other_functor``.

        - ``reason`` -- a string or ``None`` (default).

        TESTS::

            sage: from sage.categories.pushout import pushout
            sage: pushout(QQ, cartesian_product([QQ]))  # indirect doctest
            Traceback (most recent call last):
            ...
            CoercionException: No common base ("join") found for
            FractionField(Integer Ring) and The cartesian_product functorial construction(Rational Field).
        """
        if not isinstance(self_bases, (tuple, list)):
            self_bases = (self_bases,)
        if not isinstance(other_bases, (tuple, list)):
            other_bases = (other_bases,)
        if reason is None:
            reason = '.'
        else:
            reason = ': ' + reason + '.'
        raise CoercionException(
            'No common base ("join") found for %s(%s) and %s(%s)%s' %
            (self, ', '.join(str(b) for b in self_bases),
             other_functor, ', '.join(str(b) for b in other_bases),
             reason))


class CompositeConstructionFunctor(ConstructionFunctor):
    """
    A Construction Functor composed by other Construction Functors.

    INPUT:

    ``F1, F2,...``: A list of Construction Functors. The result is the
    composition ``F1`` followed by ``F2`` followed by ...

    EXAMPLES::

        sage: from sage.categories.pushout import CompositeConstructionFunctor
        sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
        sage: F
        Poly[y](FractionField(Poly[x](FractionField(...))))
        sage: F == loads(dumps(F))
        True
        sage: F == CompositeConstructionFunctor(*F.all)
        True
        sage: F(GF(2)['t'])
        Univariate Polynomial Ring in y over Fraction Field of Univariate Polynomial Ring in x over Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)

    """

    def __init__(self, *args):
        """
        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F
            Poly[y](FractionField(Poly[x](FractionField(...))))
            sage: F == CompositeConstructionFunctor(*F.all)
            True

        """
        self.all = []
        for c in args:
            if isinstance(c, list):
                self.all += c
            elif isinstance(c, CompositeConstructionFunctor):
                self.all += c.all
            else:
                self.all.append(c)
        Functor.__init__(self, self.all[0].domain(), self.all[-1].codomain())

    def _apply_functor_to_morphism(self, f):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: R.<a,b> = QQ[]
            sage: f = R.hom([a+b, a-b])
            sage: F(f)           # indirect doctest
            Ring endomorphism of Univariate Polynomial Ring in y over Fraction Field of Univariate Polynomial Ring in x over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
              Defn: Induced from base ring by
                    Ring endomorphism of Fraction Field of Univariate Polynomial Ring in x over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
                      Defn: Induced from base ring by
                            Ring endomorphism of Univariate Polynomial Ring in x over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
                              Defn: Induced from base ring by
                                    Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field
                                      Defn: a |--> a + b
                                            b |--> a - b

        """
        for c in self.all:
            f = c(f)
        return f

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: R.<a,b> = QQ[]
            sage: F(R)       # indirect doctest
            Univariate Polynomial Ring in y over Fraction Field of Univariate Polynomial Ring in x over Fraction Field of Multivariate Polynomial Ring in a, b over Rational Field

        """
        for c in self.all:
            R = c(R)
        return R

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F == loads(dumps(F)) # indirect doctest
            True

        """
        if isinstance(other, CompositeConstructionFunctor):
            return cmp(self.all, other.all)
        else:
            return cmp(type(self), type(other))

    def __mul__(self, other):
        """
        Compose construction functors to a composit construction functor, unless one of them is the identity.

        NOTE:

        The product is in functorial notation, i.e., when applying the product to an object
        then the second factor is applied first.

        EXAMPLES::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F1 = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0])
            sage: F2 = CompositeConstructionFunctor(QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F1*F2
            Poly[x](FractionField(Poly[y](FractionField(...))))

        """
        if isinstance(self, CompositeConstructionFunctor):
            all = [other] + self.all
        elif isinstance(other,IdentityConstructionFunctor):
            return self
        else:
            all = other.all + [self]
        return CompositeConstructionFunctor(*all)

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F     # indirect doctest
            Poly[y](FractionField(Poly[x](FractionField(...))))

        """
        s = "..."
        for c in self.all:
            s = "%s(%s)" % (c,s)
        return s

    def expand(self):
        """
        Return expansion of a CompositeConstructionFunctor.

        NOTE:

        The product over the list of components, as returned by
        the ``expand()`` method, is equal to ``self``.

        EXAMPLES::

            sage: from sage.categories.pushout import CompositeConstructionFunctor
            sage: F = CompositeConstructionFunctor(QQ.construction()[0],ZZ['x'].construction()[0],QQ.construction()[0],ZZ['y'].construction()[0])
            sage: F
            Poly[y](FractionField(Poly[x](FractionField(...))))
            sage: prod(F.expand()) == F
            True

        """
        return list(reversed(self.all))


class IdentityConstructionFunctor(ConstructionFunctor):
    """
    A construction functor that is the identity functor.

    TESTS::

        sage: from sage.categories.pushout import IdentityConstructionFunctor
        sage: I = IdentityConstructionFunctor()
        sage: I(RR) is RR
        True
        sage: I == loads(dumps(I))
        True

    """
    rank = -100

    def __init__(self):
        """
        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: IdentityFunctor(Sets()) == I
            True
            sage: I(RR) is RR
            True

        """
        from sage.categories.sets_cat import Sets
        ConstructionFunctor.__init__(self, Sets(), Sets())

    def _apply_functor(self, x):
        """
        Return the argument unaltered.

        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: I(RR) is RR      # indirect doctest
            True
        """
        return x

    def _apply_functor_to_morphism(self, f):
        """
        Return the argument unaltered.

        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: f = ZZ['t'].hom(['x'],QQ['x'])
            sage: I(f) is f      # indirect doctest
            True
        """
        return f

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: I == IdentityFunctor(Sets())     # indirect doctest
            True
            sage: I == QQ.construction()[0]
            False

        """
        c = cmp(type(self),type(other))
        if c:
            from sage.categories.functor import IdentityFunctor_generic
            if isinstance(other,IdentityFunctor_generic):
               return 0
        return c

    def __mul__(self, other):
        """
        Compose construction functors to a composit construction functor, unless one of them is the identity.

        NOTE:

        The product is in functorial notation, i.e., when applying the product to an object
        then the second factor is applied first.

        TESTS::

            sage: from sage.categories.pushout import IdentityConstructionFunctor
            sage: I = IdentityConstructionFunctor()
            sage: F = QQ.construction()[0]
            sage: P = ZZ['t'].construction()[0]
            sage: I*F is F     # indirect doctest
            True
            sage: F*I is F
            True
            sage: I*P is P
            True
            sage: P*I is P
            True

        """
        if isinstance(self, IdentityConstructionFunctor):
            return other
        else:
            return self


class MultivariateConstructionFunctor(ConstructionFunctor):
    """
    An abstract base class for functors that take
    multiple inputs (e.g. Cartesian products).

    TESTS::

        sage: from sage.categories.pushout import pushout
        sage: A = cartesian_product((QQ['z'], QQ))
        sage: B = cartesian_product((ZZ['t']['z'], QQ))
        sage: pushout(A, B)
        The Cartesian product of (Univariate Polynomial Ring in z over
        Univariate Polynomial Ring in t over Rational Field,
        Rational Field)
        sage: A.construction()
        (The cartesian_product functorial construction,
         (Univariate Polynomial Ring in z over Rational Field, Rational Field))
        sage: pushout(A, B)
        The Cartesian product of (Univariate Polynomial Ring in z over Univariate Polynomial Ring in t over Rational Field, Rational Field)
    """
    def common_base(self, other_functor, self_bases, other_bases):
        r"""
        This function is called by :func:`pushout` when no common parent
        is found in the construction tower.

        INPUT:

        - ``other_functor`` -- a construction functor.

        - ``self_bases`` -- the arguments passed to this functor.

        - ``other_bases`` -- the arguments passed to the functor
          ``other_functor``.

        OUTPUT:

        A parent.

        If no common base is found a :class:`sage.structure.coerce_exceptions.CoercionException`
        is raised.

        .. NOTE::

            Overload this function in derived class, see
            e.g. :class:`MultivariateConstructionFunctor`.

        TESTS::

            sage: from sage.categories.pushout import pushout
            sage: pushout(cartesian_product([ZZ]), QQ)  # indirect doctest
            Traceback (most recent call last):
            ...
            CoercionException: No common base ("join") found for
            The cartesian_product functorial construction(Integer Ring) and FractionField(Integer Ring):
            (Multivariate) functors are incompatible.
            sage: pushout(cartesian_product([ZZ]), cartesian_product([ZZ, QQ]))  # indirect doctest
            Traceback (most recent call last):
            ...
            CoercionException: No common base ("join") found for
            The cartesian_product functorial construction(Integer Ring) and
            The cartesian_product functorial construction(Integer Ring, Rational Field):
            Functors need the same number of arguments.
        """
        if self != other_functor:
            self._raise_common_base_exception_(
                other_functor, self_bases, other_bases,
                '(Multivariate) functors are incompatible')
        if len(self_bases) != len(other_bases):
            self._raise_common_base_exception_(
                other_functor, self_bases, other_bases,
                'Functors need the same number of arguments')
        from sage.structure.element import get_coercion_model
        Z_bases = tuple(get_coercion_model().common_parent(S, O)
                        for S, O in zip(self_bases, other_bases))
        return self(Z_bases)


class PolynomialFunctor(ConstructionFunctor):
    """
    Construction functor for univariate polynomial rings.

    EXAMPLE::

        sage: P = ZZ['t'].construction()[0]
        sage: P(GF(3))
        Univariate Polynomial Ring in t over Finite Field of size 3
        sage: P == loads(dumps(P))
        True
        sage: R.<x,y> = GF(5)[]
        sage: f = R.hom([x+2*y,3*x-y],R)
        sage: P(f)((x+y)*P(R).0)
        (-x + y)*t

    By :trac:`9944`, the construction functor distinguishes sparse and
    dense polynomial rings. Before, the following example failed::

        sage: R.<x> = PolynomialRing(GF(5), sparse=True)
        sage: F,B = R.construction()
        sage: F(B) is R
        True
        sage: S.<x> = PolynomialRing(ZZ)
        sage: R.has_coerce_map_from(S)
        False
        sage: S.has_coerce_map_from(R)
        False
        sage: S.0 + R.0
        2*x
        sage: (S.0 + R.0).parent()
        Univariate Polynomial Ring in x over Finite Field of size 5
        sage: (S.0 + R.0).parent().is_sparse()
        False

    """
    rank = 9

    def __init__(self, var, multi_variate=False, sparse=False):
        """
        TESTS::

            sage: from sage.categories.pushout import PolynomialFunctor
            sage: P = PolynomialFunctor('x')
            sage: P(GF(3))
            Univariate Polynomial Ring in x over Finite Field of size 3

        There is an optional parameter ``multi_variate``, but
        apparently it is not used::

            sage: Q = PolynomialFunctor('x',multi_variate=True)
            sage: Q(ZZ)
            Univariate Polynomial Ring in x over Integer Ring
            sage: Q == P
            True

        """
        from rings import Rings
        Functor.__init__(self, Rings(), Rings())
        self.var = var
        self.multi_variate = multi_variate
        self.sparse = sparse

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST::

            sage: P = ZZ['x'].construction()[0]
            sage: P(GF(3))      # indirect doctest
            Univariate Polynomial Ring in x over Finite Field of size 3

        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(R, self.var, sparse=self.sparse)

    def _apply_functor_to_morphism(self, f):
        """
        Apply the functor ``self`` to the morphism `f`.

        TEST::

            sage: P = ZZ['x'].construction()[0]
            sage: P(ZZ.hom(GF(3)))  # indirect doctest
            Ring morphism:
              From: Univariate Polynomial Ring in x over Integer Ring
              To:   Univariate Polynomial Ring in x over Finite Field of size 3
              Defn: Induced from base ring by
                    Ring Coercion morphism:
                      From: Integer Ring
                      To:   Finite Field of size 3
        """
        from sage.rings.polynomial.polynomial_ring_homomorphism import PolynomialRingHomomorphism_from_base
        R = self._apply_functor(f.domain())
        S = self._apply_functor(f.codomain())
        return PolynomialRingHomomorphism_from_base(R.Hom(S), f)

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.categories.pushout import MultiPolynomialFunctor
            sage: Q = MultiPolynomialFunctor(('x',),'lex')
            sage: P = ZZ['x'].construction()[0]
            sage: P
            Poly[x]
            sage: Q
            MPoly[x]
            sage: P == Q
            True
            sage: P == loads(dumps(P))
            True
            sage: P == QQ.construction()[0]
            False

        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.var, other.var)
        elif isinstance(other, MultiPolynomialFunctor):
            return -cmp(other, self)
        return c

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return None.

        NOTE:

        Internally, the merging is delegated to the merging of
        multipolynomial construction functors. But in effect,
        this does the same as the default implementation, that
        returns ``None`` unless the to-be-merged functors coincide.

        EXAMPLE::

            sage: P = ZZ['x'].construction()[0]
            sage: Q = ZZ['y','x'].construction()[0]
            sage: P.merge(Q)
            sage: P.merge(P) is P
            True

        """
        if isinstance(other, MultiPolynomialFunctor):
            return other.merge(self)
        elif self == other:
            # i.e., they only differ in sparsity
            if not self.sparse:
                return self
            return other
        else:
            return None

    def _repr_(self):
        """
        TEST::

            sage: P = ZZ['x'].construction()[0]
            sage: P       # indirect doctest
            Poly[x]

        """
        return "Poly[%s]" % self.var

class MultiPolynomialFunctor(ConstructionFunctor):
    """
    A constructor for multivariate polynomial rings.

    EXAMPLES::

        sage: P.<x,y> = ZZ[]
        sage: F = P.construction()[0]; F
        MPoly[x,y]
        sage: A.<a,b> = GF(5)[]
        sage: F(A)
        Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Finite Field of size 5
        sage: f = A.hom([a+b,a-b],A)
        sage: F(f)
        Ring endomorphism of Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring in a, b over Finite Field of size 5
          Defn: Induced from base ring by
                Ring endomorphism of Multivariate Polynomial Ring in a, b over Finite Field of size 5
                  Defn: a |--> a + b
                        b |--> a - b
        sage: F(f)(F(A)(x)*a)
        (a + b)*x

    """

    rank = 9

    def __init__(self, vars, term_order):
        """
        EXAMPLES::

            sage: F = sage.categories.pushout.MultiPolynomialFunctor(['x','y'], None)
            sage: F
            MPoly[x,y]
            sage: F(ZZ)
            Multivariate Polynomial Ring in x, y over Integer Ring
            sage: F(CC)
            Multivariate Polynomial Ring in x, y over Complex Field with 53 bits of precision
        """
        Functor.__init__(self, Rings(), Rings())
        self.vars = vars
        self.term_order = term_order

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        EXAMPLES::

            sage: R.<x,y,z> = QQ[]
            sage: F = R.construction()[0]; F
            MPoly[x,y,z]
            sage: type(F)
            <class 'sage.categories.pushout.MultiPolynomialFunctor'>
            sage: F(ZZ)          # indirect doctest
            Multivariate Polynomial Ring in x, y, z over Integer Ring
            sage: F(RR)          # indirect doctest
            Multivariate Polynomial Ring in x, y, z over Real Field with 53 bits of precision
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        return PolynomialRing(R, self.vars)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: F = ZZ['x,y,z'].construction()[0]
            sage: G = QQ['x,y,z'].construction()[0]
            sage: F == G
            True
            sage: G == loads(dumps(G))
            True
            sage: G = ZZ['x,y'].construction()[0]
            sage: F == G
            False
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.vars, other.vars) or cmp(self.term_order, other.term_order)
        elif isinstance(other, PolynomialFunctor):
            c = cmp(self.vars, (other.var,))
        return c

    def __mul__(self, other):
        """
        If two MPoly functors are given in a row, form a single MPoly functor
        with all of the variables.

        EXAMPLES::

            sage: F = sage.categories.pushout.MultiPolynomialFunctor(['x','y'], None)
            sage: G = sage.categories.pushout.MultiPolynomialFunctor(['t'], None)
            sage: G*F
            MPoly[x,y,t]
        """
        if isinstance(other,IdentityConstructionFunctor):
            return self
        if isinstance(other, MultiPolynomialFunctor):
            if self.term_order != other.term_order:
                raise CoercionException("Incompatible term orders (%s,%s)." % (self.term_order, other.term_order))
            if set(self.vars).intersection(other.vars):
                raise CoercionException("Overlapping variables (%s,%s)" % (self.vars, other.vars))
            return MultiPolynomialFunctor(other.vars + self.vars, self.term_order)
        elif isinstance(other, CompositeConstructionFunctor) \
              and isinstance(other.all[-1], MultiPolynomialFunctor):
            return CompositeConstructionFunctor(other.all[:-1], self * other.all[-1])
        else:
            return CompositeConstructionFunctor(other, self)

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return None.

        EXAMPLES::

            sage: F = sage.categories.pushout.MultiPolynomialFunctor(['x','y'], None)
            sage: G = sage.categories.pushout.MultiPolynomialFunctor(['t'], None)
            sage: F.merge(G) is None
            True
            sage: F.merge(F)
            MPoly[x,y]
        """
        if self == other:
            return self
        else:
            return None

    def expand(self):
        """
        Decompose ``self`` into a list of construction functors.

        EXAMPLES::

            sage: F = QQ['x,y,z,t'].construction()[0]; F
            MPoly[x,y,z,t]
            sage: F.expand()
            [MPoly[t], MPoly[z], MPoly[y], MPoly[x]]

        Now an actual use case::

            sage: R.<x,y,z> = ZZ[]
            sage: S.<z,t> = QQ[]
            sage: x+t
            x + t
            sage: parent(x+t)
            Multivariate Polynomial Ring in x, y, z, t over Rational Field
            sage: T.<y,s> = QQ[]
            sage: x + s
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Multivariate Polynomial Ring in x, y, z over Integer Ring' and 'Multivariate Polynomial Ring in y, s over Rational Field'
            sage: R = PolynomialRing(ZZ, 'x', 500)
            sage: S = PolynomialRing(GF(5), 'x', 200)
            sage: R.gen(0) + S.gen(0)
            2*x0
        """
        if len(self.vars) <= 1:
            return [self]
        else:
            return [MultiPolynomialFunctor((x,), self.term_order) for x in reversed(self.vars)]

    def _repr_(self):
        """
        TEST::

            sage: QQ['x,y,z,t'].construction()[0]
            MPoly[x,y,z,t]
        """
        return "MPoly[%s]" % ','.join(self.vars)



class InfinitePolynomialFunctor(ConstructionFunctor):
    """
    A Construction Functor for Infinite Polynomial Rings (see :mod:`~sage.rings.polynomial.infinite_polynomial_ring`).

    AUTHOR:

    -- Simon King

    This construction functor is used to provide uniqueness of infinite polynomial rings as parent structures.
    As usual, the construction functor allows for constructing pushouts.

    Another purpose is to avoid name conflicts of variables of the to-be-constructed infinite polynomial ring with
    variables of the base ring, and moreover to keep the internal structure of an Infinite Polynomial Ring as simple
    as possible: If variables `v_1,...,v_n` of the given base ring generate an *ordered* sub-monoid of the monomials
    of the ambient Infinite Polynomial Ring, then they are removed from the base ring and merged with the generators
    of the ambient ring. However, if the orders don't match, an error is raised, since there was a name conflict
    without merging.

    EXAMPLES::

        sage: A.<a,b> = InfinitePolynomialRing(ZZ['t'])
        sage: A.construction()
        [InfPoly{[a,b], "lex", "dense"},
         Univariate Polynomial Ring in t over Integer Ring]
        sage: type(_[0])
        <class 'sage.categories.pushout.InfinitePolynomialFunctor'>
        sage: B.<x,y,a_3,a_1> = PolynomialRing(QQ, order='lex')
        sage: B.construction()
        (MPoly[x,y,a_3,a_1], Rational Field)
        sage: A.construction()[0]*B.construction()[0]
        InfPoly{[a,b], "lex", "dense"}(MPoly[x,y](...))

    Apparently the variables `a_1,a_3` of the polynomial ring are merged with the variables
    `a_0, a_1, a_2, ...` of the infinite polynomial ring; indeed, they form an ordered sub-structure.
    However, if the polynomial ring was given a different ordering, merging would not be allowed,
    resulting in a name conflict::

        sage: A.construction()[0]*PolynomialRing(QQ,names=['x','y','a_3','a_1']).construction()[0]
        Traceback (most recent call last):
        ...
        CoercionException: Incompatible term orders lex, degrevlex

    In an infinite polynomial ring with generator `a_\\ast`, the variable `a_3` will always be greater
    than the variable `a_1`. Hence, the orders are incompatible in the next example as well::

        sage: A.construction()[0]*PolynomialRing(QQ,names=['x','y','a_1','a_3'], order='lex').construction()[0]
        Traceback (most recent call last):
        ...
        CoercionException: Overlapping variables (('a', 'b'),['a_1', 'a_3']) are incompatible

    Another requirement is that after merging the order of the remaining variables must be unique.
    This is not the case in the following example, since it is not clear whether the variables `x,y`
    should be greater or smaller than the variables `b_\\ast`::

        sage: A.construction()[0]*PolynomialRing(QQ,names=['a_3','a_1','x','y'], order='lex').construction()[0]
        Traceback (most recent call last):
        ...
        CoercionException: Overlapping variables (('a', 'b'),['a_3', 'a_1']) are incompatible

    Since the construction functors are actually used to construct infinite polynomial rings, the following
    result is no surprise::

        sage: C.<a,b> = InfinitePolynomialRing(B); C
        Infinite polynomial ring in a, b over Multivariate Polynomial Ring in x, y over Rational Field

    There is also an overlap in the next example::

        sage: X.<w,x,y> = InfinitePolynomialRing(ZZ)
        sage: Y.<x,y,z> = InfinitePolynomialRing(QQ)

    `X` and `Y` have an overlapping generators `x_\\ast, y_\\ast`. Since the default lexicographic order is
    used in both rings, it gives rise to isomorphic sub-monoids in both `X` and `Y`. They are merged in the
    pushout, which also yields a common parent for doing arithmetic::

        sage: P = sage.categories.pushout.pushout(Y,X); P
        Infinite polynomial ring in w, x, y, z over Rational Field
        sage: w[2]+z[3]
        w_2 + z_3
        sage: _.parent() is P
        True

    """

    # We do provide merging with polynomial rings. However, it seems that it is better
    # to have a greater rank, since we want to apply InfinitePolynomialFunctor *after*
    # [Multi]PolynomialFunctor, which have rank 9. But there is the MatrixFunctor, which
    # has rank 10. So, do fine tuning...
    rank = 9.5

    def __init__(self, gens, order, implementation):
        """
        TEST::

            sage: F = sage.categories.pushout.InfinitePolynomialFunctor(['a','b','x'],'degrevlex','sparse'); F # indirect doctest
            InfPoly{[a,b,x], "degrevlex", "sparse"}
            sage: F == loads(dumps(F))
            True

        """
        if len(gens)<1:
            raise ValueError("Infinite Polynomial Rings have at least one generator")
        ConstructionFunctor.__init__(self, Rings(), Rings())
        self._gens = tuple(gens)
        self._order = order
        self._imple = implementation

    def _apply_functor_to_morphism(self, f):
        """
        Morphisms for inifinite polynomial rings are not implemented yet.

        TEST::

            sage: P.<x,y> = QQ[]
            sage: R.<alpha> = InfinitePolynomialRing(P)
            sage: f = P.hom([x+y,x-y],P)
            sage: R.construction()[0](f)     # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Morphisms for inifinite polynomial rings are not implemented yet.

        """
        raise NotImplementedError("Morphisms for inifinite polynomial rings are not implemented yet.")

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST::

            sage: F = sage.categories.pushout.InfinitePolynomialFunctor(['a','b','x'],'degrevlex','sparse'); F
            InfPoly{[a,b,x], "degrevlex", "sparse"}
            sage: F(QQ['t']) # indirect doctest
            Infinite polynomial ring in a, b, x over Univariate Polynomial Ring in t over Rational Field

        """
        from sage.rings.polynomial.infinite_polynomial_ring import InfinitePolynomialRing
        return InfinitePolynomialRing(R, self._gens, order=self._order, implementation=self._imple)

    def _repr_(self):
        """
        TEST::

            sage: F = sage.categories.pushout.InfinitePolynomialFunctor(['a','b','x'],'degrevlex','sparse'); F # indirect doctest
            InfPoly{[a,b,x], "degrevlex", "sparse"}

        """
        return 'InfPoly{[%s], "%s", "%s"}'%(','.join(self._gens), self._order, self._imple)

    def __cmp__(self, other):
        """
        TEST::

            sage: F = sage.categories.pushout.InfinitePolynomialFunctor(['a','b','x'],'degrevlex','sparse'); F # indirect doctest
            InfPoly{[a,b,x], "degrevlex", "sparse"}
            sage: F == loads(dumps(F)) # indirect doctest
            True
            sage: F == sage.categories.pushout.InfinitePolynomialFunctor(['a','b','x'],'deglex','sparse')
            False

        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self._gens, other._gens) or cmp(self._order, other._order) or cmp(self._imple, other._imple)
        return c

    def __mul__(self, other):
        """
        Compose construction functors to a composit construction functor, unless one of them is the identity.

        NOTE:

        The product is in functorial notation, i.e., when applying the product to an object
        then the second factor is applied first.

        TESTS::

            sage: F1 = QQ['a','x_2','x_1','y_3','y_2'].construction()[0]; F1
            MPoly[a,x_2,x_1,y_3,y_2]
            sage: F2 = InfinitePolynomialRing(QQ, ['x','y'],order='degrevlex').construction()[0]; F2
            InfPoly{[x,y], "degrevlex", "dense"}
            sage: F3 = InfinitePolynomialRing(QQ, ['x','y'],order='degrevlex',implementation='sparse').construction()[0]; F3
            InfPoly{[x,y], "degrevlex", "sparse"}
            sage: F2*F1
            InfPoly{[x,y], "degrevlex", "dense"}(Poly[a](...))
            sage: F3*F1
            InfPoly{[x,y], "degrevlex", "sparse"}(Poly[a](...))
            sage: F4 = sage.categories.pushout.FractionField()
            sage: F2*F4
            InfPoly{[x,y], "degrevlex", "dense"}(FractionField(...))

        """
        if isinstance(other,IdentityConstructionFunctor):
            return self
        if isinstance(other, self.__class__): #
            INT = set(self._gens).intersection(other._gens)
            if INT:
                # if there is overlap of generators, it must only be at the ends, so that
                # the resulting order after the merging is unique
                if other._gens[-len(INT):] != self._gens[:len(INT)]:
                    raise CoercionException("Overlapping variables (%s,%s) are incompatible" % (self._gens, other._gens))
                OUTGENS = list(other._gens) + list(self._gens[len(INT):])
            else:
                OUTGENS = list(other._gens) + list(self._gens)
            # the orders must coincide
            if self._order != other._order:
                return CompositeConstructionFunctor(other, self)
            # the implementations must coincide
            if self._imple != other._imple:
                return CompositeConstructionFunctor(other, self)
            return InfinitePolynomialFunctor(OUTGENS, self._order, self._imple)

        # Polynomial Constructor
        # Idea: We merge into self, if the polynomial functor really provides a substructure,
        # even respecting the order. Note that, if the pushout is computed, only *one* variable
        # will occur in the polynomial constructor. Hence, any order is fine, which is exactly
        # what we need in order to have coercion maps for different orderings.
        if isinstance(other, MultiPolynomialFunctor) or isinstance(other, PolynomialFunctor):
            if isinstance(other, MultiPolynomialFunctor):
                othervars = other.vars
            else:
                othervars = [other.var]
            OverlappingGens = [] ## Generator names of variable names of the MultiPolynomialFunctor
                              ## that can be interpreted as variables in self
            OverlappingVars = [] ## The variable names of the MultiPolynomialFunctor
                                 ## that can be interpreted as variables in self
            RemainingVars = [x for x in othervars]
            IsOverlap = False
            BadOverlap = False
            for x in othervars:
                if x.count('_') == 1:
                    g,n = x.split('_')
                    if n.isdigit():
                        if g.isalnum(): # we can interprete x in any InfinitePolynomialRing
                            if g in self._gens: # we can interprete x in self, hence, we will not use it as a variable anymore.
                                RemainingVars.pop(RemainingVars.index(x))
                                IsOverlap = True # some variables of other can be interpreted in self.
                                if OverlappingVars:
                                    # Is OverlappingVars in the right order?
                                    g0,n0 = OverlappingVars[-1].split('_')
                                    i = self._gens.index(g)
                                    i0 = self._gens.index(g0)
                                    if i<i0: # wrong order
                                        BadOverlap = True
                                    if i==i0 and int(n)>int(n0): # wrong order
                                        BadOverlap = True
                                OverlappingVars.append(x)
                            else:
                                if IsOverlap: # The overlap must be on the right end of the variable list
                                    BadOverlap = True
                        else:
                            if IsOverlap: # The overlap must be on the right end of the variable list
                                BadOverlap = True
                    else:
                        if IsOverlap: # The overlap must be on the right end of the variable list
                            BadOverlap = True
                else:
                    if IsOverlap: # The overlap must be on the right end of the variable list
                        BadOverlap = True

            if BadOverlap: # the overlapping variables appear in the wrong order
                raise CoercionException("Overlapping variables (%s,%s) are incompatible" % (self._gens, OverlappingVars))
            if len(OverlappingVars)>1: # multivariate, hence, the term order matters
                if other.term_order.name()!=self._order:
                    raise CoercionException("Incompatible term orders %s, %s" % (self._order, other.term_order.name()))
            # ok, the overlap is fine, we will return something.
            if RemainingVars: # we can only partially merge other into self
                if len(RemainingVars)>1:
                    return CompositeConstructionFunctor(MultiPolynomialFunctor(RemainingVars,term_order=other.term_order), self)
                return CompositeConstructionFunctor(PolynomialFunctor(RemainingVars[0]), self)
            return self
        return CompositeConstructionFunctor(other, self)

    def merge(self,other):
        """
        Merge two construction functors of infinite polynomial rings, regardless of monomial order and implementation.

        The purpose is to have a pushout (and thus, arithmetic) even in cases when the parents are isomorphic as
        rings, but not as ordered rings.

        EXAMPLES::

            sage: X.<x,y> = InfinitePolynomialRing(QQ,implementation='sparse')
            sage: Y.<x,y> = InfinitePolynomialRing(QQ,order='degrevlex')
            sage: X.construction()
            [InfPoly{[x,y], "lex", "sparse"}, Rational Field]
            sage: Y.construction()
            [InfPoly{[x,y], "degrevlex", "dense"}, Rational Field]
            sage: Y.construction()[0].merge(Y.construction()[0])
            InfPoly{[x,y], "degrevlex", "dense"}
            sage: y[3] + X(x[2])
            x_2 + y_3
            sage: _.parent().construction()
            [InfPoly{[x,y], "degrevlex", "dense"}, Rational Field]

        """
        # Merging is only done if the ranks of self and other are the same.
        # It may happen that other is a substructure of self up to the monomial order
        # and the implementation. And this is when we want to merge, in order to
        # provide multiplication for rings with different term orderings.
        if not isinstance(other, InfinitePolynomialFunctor):
            return None
        if set(other._gens).issubset(self._gens):
            return self
        return None
        try:
            OUT = self*other
            # The following happens if "other" has the same order type etc.
            if not isinstance(OUT, CompositeConstructionFunctor):
                return OUT
        except CoercionException:
            pass
        if isinstance(other,InfinitePolynomialFunctor):
            # We don't require that the orders coincide. This is a difference to self*other
            # We only merge if other's generators are an ordered subset of self's generators
            for g in other._gens:
                if g not in self._gens:
                    return None
            # The sequence of variables is part of the ordering. It must coincide in both rings
            Ind = [self._gens.index(g) for g in other._gens]
            if sorted(Ind)!=Ind:
                return None
            # OK, other merges into self. Now, chose the default dense implementation,
            # unless both functors refer to the sparse implementation
            if self._imple != other._imple:
                return InfinitePolynomialFunctor(self._gens, self._order, 'dense')
            return self
        return None

    def expand(self):
        """
        Decompose the functor `F` into sub-functors, whose product returns `F`.

        EXAMPLES::

            sage: F = InfinitePolynomialRing(QQ, ['x','y'],order='degrevlex').construction()[0]; F
            InfPoly{[x,y], "degrevlex", "dense"}
            sage: F.expand()
            [InfPoly{[y], "degrevlex", "dense"}, InfPoly{[x], "degrevlex", "dense"}]
            sage: F = InfinitePolynomialRing(QQ, ['x','y','z'],order='degrevlex').construction()[0]; F
            InfPoly{[x,y,z], "degrevlex", "dense"}
            sage: F.expand()
            [InfPoly{[z], "degrevlex", "dense"},
             InfPoly{[y], "degrevlex", "dense"},
             InfPoly{[x], "degrevlex", "dense"}]
            sage: prod(F.expand())==F
            True

        """
        if len(self._gens)==1:
            return [self]
        return [InfinitePolynomialFunctor((x,), self._order, self._imple) for x in reversed(self._gens)]



class MatrixFunctor(ConstructionFunctor):
    """
    A construction functor for matrices over rings.

    EXAMPLES::

        sage: MS = MatrixSpace(ZZ,2, 3)
        sage: F = MS.construction()[0]; F
        MatrixFunctor
        sage: MS = MatrixSpace(ZZ,2)
        sage: F = MS.construction()[0]; F
        MatrixFunctor
        sage: P.<x,y> = QQ[]
        sage: R = F(P); R
        Full MatrixSpace of 2 by 2 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
        sage: f = P.hom([x+y,x-y],P); F(f)
        Ring endomorphism of Full MatrixSpace of 2 by 2 dense matrices over Multivariate Polynomial Ring in x, y over Rational Field
          Defn: Induced from base ring by
                Ring endomorphism of Multivariate Polynomial Ring in x, y over Rational Field
                  Defn: x |--> x + y
                        y |--> x - y
        sage: M = R([x,y,x*y,x+y])
        sage: F(f)(M)
        [    x + y     x - y]
        [x^2 - y^2       2*x]

    """
    rank = 10

    def __init__(self, nrows, ncols, is_sparse=False):
        """
        TEST::

            sage: from sage.categories.pushout import MatrixFunctor
            sage: F = MatrixFunctor(2,3)
            sage: F == MatrixSpace(ZZ,2,3).construction()[0]
            True
            sage: F.codomain()
            Category of commutative additive groups
            sage: R = MatrixSpace(ZZ,2,2).construction()[0]
            sage: R.codomain()
            Category of rings
            sage: F(ZZ)
            Full MatrixSpace of 2 by 3 dense matrices over Integer Ring
            sage: F(ZZ) in F.codomain()
            True
            sage: R(GF(2))
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 2
            sage: R(GF(2)) in R.codomain()
            True
        """
        if nrows == ncols:
            Functor.__init__(self, Rings(), Rings()) # Algebras() takes a base ring
        else:
            # Functor.__init__(self, Rings(), MatrixAlgebras()) # takes a base ring
            Functor.__init__(self, Rings(), CommutativeAdditiveGroups()) # not a nice solution, but the best we can do.
        self.nrows = nrows
        self.ncols = ncols
        self.is_sparse = is_sparse

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST:

        The following is a test against a bug discussed at :trac:`8800`::

            sage: F = MatrixSpace(ZZ,2,3).construction()[0]
            sage: F(RR)         # indirect doctest
            Full MatrixSpace of 2 by 3 dense matrices over Real Field with 53 bits of precision
            sage: F(RR) in F.codomain()
            True

        """
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(R, self.nrows, self.ncols, sparse=self.is_sparse)

    def __cmp__(self, other):
        """
        TEST::

            sage: F = MatrixSpace(ZZ,2,3).construction()[0]
            sage: F == loads(dumps(F))
            True
            sage: F == MatrixSpace(ZZ,2,2).construction()[0]
            False

        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp((self.nrows, self.ncols), (other.nrows, other.ncols))
        return c

    def merge(self, other):
        """
        Merging is only happening if both functors are matrix functors of the same dimension.
        The result is sparse if and only if both given functors are sparse.

        EXAMPLE::

            sage: F1 = MatrixSpace(ZZ,2,2).construction()[0]
            sage: F2 = MatrixSpace(ZZ,2,3).construction()[0]
            sage: F3 = MatrixSpace(ZZ,2,2,sparse=True).construction()[0]
            sage: F1.merge(F2)
            sage: F1.merge(F3)
            MatrixFunctor
            sage: F13 = F1.merge(F3)
            sage: F13.is_sparse
            False
            sage: F1.is_sparse
            False
            sage: F3.is_sparse
            True
            sage: F3.merge(F3).is_sparse
            True

        """
        if self != other:
            return None
        else:
            return MatrixFunctor(self.nrows, self.ncols, self.is_sparse and other.is_sparse)

class LaurentPolynomialFunctor(ConstructionFunctor):
    """
    Construction functor for Laurent polynomial rings.

    EXAMPLES::

        sage: L.<t> = LaurentPolynomialRing(ZZ)
        sage: F = L.construction()[0]
        sage: F
        LaurentPolynomialFunctor
        sage: F(QQ)
        Univariate Laurent Polynomial Ring in t over Rational Field
        sage: K.<x> = LaurentPolynomialRing(ZZ)
        sage: F(K)
        Univariate Laurent Polynomial Ring in t over Univariate Laurent Polynomial Ring in x over Integer Ring
        sage: P.<x,y> = ZZ[]
        sage: f = P.hom([x+2*y,3*x-y],P)
        sage: F(f)
        Ring endomorphism of Univariate Laurent Polynomial Ring in t over Multivariate Polynomial Ring in x, y over Integer Ring
          Defn: Induced from base ring by
                Ring endomorphism of Multivariate Polynomial Ring in x, y over Integer Ring
                  Defn: x |--> x + 2*y
                        y |--> 3*x - y
        sage: F(f)(x*F(P).gen()^-2+y*F(P).gen()^3)
        (x + 2*y)*t^-2 + (3*x - y)*t^3

    """
    rank = 9

    def __init__(self, var, multi_variate=False):
        """
        INPUT:

        - ``var``, a string or a list of strings
        - ``multi_variate``, optional bool, default ``False`` if ``var`` is a string
          and ``True`` otherwise: If ``True``, application to a Laurent polynomial
          ring yields a multivariate Laurent polynomial ring.

        TESTS::

            sage: from sage.categories.pushout import LaurentPolynomialFunctor
            sage: F1 = LaurentPolynomialFunctor('t')
            sage: F2 = LaurentPolynomialFunctor('s', multi_variate=True)
            sage: F3 = LaurentPolynomialFunctor(['s','t'])
            sage: F1(F2(QQ))
            Univariate Laurent Polynomial Ring in t over Univariate Laurent Polynomial Ring in s over Rational Field
            sage: F2(F1(QQ))
            Multivariate Laurent Polynomial Ring in t, s over Rational Field
            sage: F3(QQ)
            Multivariate Laurent Polynomial Ring in s, t over Rational Field

        """
        Functor.__init__(self, Rings(), Rings())
        if not isinstance(var, (six.string_types,tuple,list)):
            raise TypeError("variable name or list of variable names expected")
        self.var = var
        self.multi_variate = multi_variate or not isinstance(var, six.string_types)

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: from sage.categories.pushout import LaurentPolynomialFunctor
            sage: F1 = LaurentPolynomialFunctor('t')
            sage: F2 = LaurentPolynomialFunctor('s', multi_variate=True)
            sage: F3 = LaurentPolynomialFunctor(['s','t'])
            sage: F1(F2(QQ))          # indirect doctest
            Univariate Laurent Polynomial Ring in t over Univariate Laurent Polynomial Ring in s over Rational Field
            sage: F2(F1(QQ))
            Multivariate Laurent Polynomial Ring in t, s over Rational Field
            sage: F3(QQ)
            Multivariate Laurent Polynomial Ring in s, t over Rational Field

        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing, is_LaurentPolynomialRing
        if self.multi_variate and is_LaurentPolynomialRing(R):
            return LaurentPolynomialRing(R.base_ring(), (list(R.variable_names()) + [self.var]))
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            return LaurentPolynomialRing(R, self.var)

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.categories.pushout import LaurentPolynomialFunctor
            sage: F1 = LaurentPolynomialFunctor('t')
            sage: F2 = LaurentPolynomialFunctor('t', multi_variate=True)
            sage: F3 = LaurentPolynomialFunctor(['s','t'])
            sage: F1 == F2
            True
            sage: F1 == loads(dumps(F1))
            True
            sage: F1 == F3
            False
            sage: F1 == QQ.construction()[0]
            False

        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.var, other.var)
        return c

    def merge(self, other):
        """
        Two Laurent polynomial construction functors merge if the variable names coincide.
        The result is multivariate if one of the arguments is multivariate.

        EXAMPLE::

            sage: from sage.categories.pushout import LaurentPolynomialFunctor
            sage: F1 = LaurentPolynomialFunctor('t')
            sage: F2 = LaurentPolynomialFunctor('t', multi_variate=True)
            sage: F1.merge(F2)
            LaurentPolynomialFunctor
            sage: F1.merge(F2)(LaurentPolynomialRing(GF(2),'a'))
            Multivariate Laurent Polynomial Ring in a, t over Finite Field of size 2
            sage: F1.merge(F1)(LaurentPolynomialRing(GF(2),'a'))
            Univariate Laurent Polynomial Ring in t over Univariate Laurent Polynomial Ring in a over Finite Field of size 2

        """
        if self == other or isinstance(other, PolynomialFunctor) and self.var == other.var:
            return LaurentPolynomialFunctor(self.var, (self.multi_variate or other.multi_variate))
        else:
            return None


class VectorFunctor(ConstructionFunctor):
    """
    A construction functor for free modules over commutative rings.

    EXAMPLE::

        sage: F = (ZZ^3).construction()[0]
        sage: F
        VectorFunctor
        sage: F(GF(2)['t'])
        Ambient free module of rank 3 over the principal ideal domain Univariate Polynomial Ring in t over Finite Field of size 2 (using NTL)


    """
    rank = 10 # ranking of functor, not rank of module.
    # This coincides with the rank of the matrix construction functor, but this is OK since they can not both be applied in any order

    def __init__(self, n, is_sparse=False, inner_product_matrix=None):
        """
        INPUT:

        - ``n``, the rank of the to-be-created modules (non-negative integer)
        - ``is_sparse`` (optional bool, default ``False``), create sparse implementation of modules
        - ``inner_product_matrix``: ``n`` by ``n`` matrix, used to compute inner products in the
          to-be-created modules

        TEST::

            sage: from sage.categories.pushout import VectorFunctor
            sage: F1 = VectorFunctor(3, inner_product_matrix = Matrix(3,3,range(9)))
            sage: F1.domain()
            Category of commutative rings
            sage: F1.codomain()
            Category of commutative additive groups
            sage: M1 = F1(ZZ)
            sage: M1.is_sparse()
            False
            sage: v = M1([3, 2, 1])
            sage: v*Matrix(3,3,range(9))*v.column()
            (96)
            sage: v.inner_product(v)
            96
            sage: F2 = VectorFunctor(3, is_sparse=True)
            sage: M2 = F2(QQ); M2; M2.is_sparse()
            Sparse vector space of dimension 3 over Rational Field
            True

        """
#        Functor.__init__(self, Rings(), FreeModules()) # FreeModules() takes a base ring
#        Functor.__init__(self, Objects(), Objects())   # Object() makes no sence, since FreeModule raises an error, e.g., on Set(['a',1]).
        ## FreeModule requires a commutative ring. Thus, we have
        Functor.__init__(self, CommutativeRings(), CommutativeAdditiveGroups())
        self.n = n
        self.is_sparse = is_sparse
        self.inner_product_matrix = inner_product_matrix

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: from sage.categories.pushout import VectorFunctor
            sage: F1 = VectorFunctor(3, inner_product_matrix = Matrix(3,3,range(9)))
            sage: M1 = F1(ZZ)   # indirect doctest
            sage: M1.is_sparse()
            False
            sage: v = M1([3, 2, 1])
            sage: v*Matrix(3,3,range(9))*v.column()
            (96)
            sage: v.inner_product(v)
            96
            sage: F2 = VectorFunctor(3, is_sparse=True)
            sage: M2 = F2(QQ); M2; M2.is_sparse()
            Sparse vector space of dimension 3 over Rational Field
            True
            sage: v = M2([3, 2, 1])
            sage: v.inner_product(v)
            14

        """
        from sage.modules.free_module import FreeModule
        return FreeModule(R, self.n, sparse=self.is_sparse, inner_product_matrix=self.inner_product_matrix)

    def _apply_functor_to_morphism(self, f):
        """
        This is not implemented yet.

        TEST::

            sage: F = (ZZ^3).construction()[0]
            sage: P.<x,y> = ZZ[]
            sage: f = P.hom([x+2*y,3*x-y],P)
            sage: F(f)       # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Can not create induced morphisms of free modules yet
        """
        ## TODO: Implement this!
        raise NotImplementedError("Can not create induced morphisms of free modules yet")

    def __cmp__(self, other):
        """
        Only the rank of the to-be-created modules is compared, *not* the inner product matrix.

        TESTS::

            sage: from sage.categories.pushout import VectorFunctor
            sage: F1 = VectorFunctor(3, inner_product_matrix = Matrix(3,3,range(9)))
            sage: F2 = (ZZ^3).construction()[0]
            sage: F1 == F2
            True
            sage: F1(QQ) == F2(QQ)
            True
            sage: F1(QQ).inner_product_matrix() == F2(QQ).inner_product_matrix()
            False
            sage: F1 == loads(dumps(F1))
            True
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.n, other.n)
        return c

    def merge(self, other):
        """
        Two constructors of free modules merge, if the module ranks coincide. If both
        have explicitly given inner product matrices, they must coincide as well.

        EXAMPLE:

        Two modules without explicitly given inner product allow coercion::

            sage: M1 = QQ^3
            sage: P.<t> = ZZ[]
            sage: M2 = FreeModule(P,3)
            sage: M1([1,1/2,1/3]) + M2([t,t^2+t,3])     # indirect doctest
            (t + 1, t^2 + t + 1/2, 10/3)

        If only one summand has an explicit inner product, the result will be provided
        with it::

            sage: M3 = FreeModule(P,3, inner_product_matrix = Matrix(3,3,range(9)))
            sage: M1([1,1/2,1/3]) + M3([t,t^2+t,3])
            (t + 1, t^2 + t + 1/2, 10/3)
            sage: (M1([1,1/2,1/3]) + M3([t,t^2+t,3])).parent().inner_product_matrix()
            [0 1 2]
            [3 4 5]
            [6 7 8]

        If both summands have an explicit inner product (even if it is the standard
        inner product), then the products must coincide. The only difference between
        ``M1`` and ``M4`` in the following example is the fact that the default
        inner product was *explicitly* requested for ``M4``. It is therefore not
        possible to coerce with a different inner product::

            sage: M4 = FreeModule(QQ,3, inner_product_matrix = Matrix(3,3,1))
            sage: M4 == M1
            True
            sage: M4.inner_product_matrix() == M1.inner_product_matrix()
            True
            sage: M4([1,1/2,1/3]) + M3([t,t^2+t,3])      # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Ambient quadratic space of dimension 3 over Rational Field
            Inner product matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]' and 'Ambient free quadratic module of rank 3 over the integral domain Univariate Polynomial Ring in t over Integer Ring
            Inner product matrix:
            [0 1 2]
            [3 4 5]
            [6 7 8]'

        """
        if self != other:
            return None
        if self.inner_product_matrix is None:
            return VectorFunctor(self.n, self.is_sparse and other.is_sparse, other.inner_product_matrix)
        if other.inner_product_matrix is None:
            return VectorFunctor(self.n, self.is_sparse and other.is_sparse, self.inner_product_matrix)
        # At this point, we know that the user wants to take care of the inner product.
        # So, we only merge if both coincide:
        if self.inner_product_matrix != other.inner_product_matrix:
            return None
        else:
            return VectorFunctor(self.n, self.is_sparse and other.is_sparse, self.inner_product_matrix)

class SubspaceFunctor(ConstructionFunctor):
    """
    Constructing a subspace of an ambient free module, given by a basis.

    NOTE:

    This construction functor keeps track of the basis. It can only be applied
    to free modules into which this basis coerces.

    EXAMPLES::

        sage: M = ZZ^3
        sage: S = M.submodule([(1,2,3),(4,5,6)]); S
        Free module of degree 3 and rank 2 over Integer Ring
        Echelon basis matrix:
        [1 2 3]
        [0 3 6]
        sage: F = S.construction()[0]
        sage: F(GF(2)^3)
        Vector space of degree 3 and dimension 2 over Finite Field of size 2
        User basis matrix:
        [1 0 1]
        [0 1 0]

    """
    rank = 11 # ranking of functor, not rank of module

    # The subspace construction returns an object admitting a coercion
    # map into the original, not vice versa.
    coercion_reversed = True

    def __init__(self, basis):
        """
        INPUT:

        ``basis``: a list of elements of a free module.

        TEST::

            sage: from sage.categories.pushout import SubspaceFunctor
            sage: M = ZZ^3
            sage: F = SubspaceFunctor([M([1,2,3]),M([4,5,6])])
            sage: F(GF(5)^3)
            Vector space of degree 3 and dimension 2 over Finite Field of size 5
            User basis matrix:
            [1 2 3]
            [4 0 1]
        """
##        Functor.__init__(self, FreeModules(), FreeModules()) # takes a base ring
##        Functor.__init__(self, Objects(), Objects())   # is too general
        ## It seems that the category of commutative additive groups
        ## currently is the smallest base ring free category that
        ## contains in- and output
        Functor.__init__(self, CommutativeAdditiveGroups(), CommutativeAdditiveGroups())
        self.basis = basis

    def _apply_functor(self, ambient):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: M = ZZ^3
            sage: S = M.submodule([(1,2,3),(4,5,6)]); S
            Free module of degree 3 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 2 3]
            [0 3 6]
            sage: F = S.construction()[0]
            sage: F(GF(2)^3)    # indirect doctest
            Vector space of degree 3 and dimension 2 over Finite Field of size 2
            User basis matrix:
            [1 0 1]
            [0 1 0]
        """
        return ambient.span_of_basis(self.basis)

    def _apply_functor_to_morphism(self, f):
        """
        This is not implemented yet.

        TEST::

            sage: F = (ZZ^3).span([(1,2,3),(4,5,6)]).construction()[0]
            sage: P.<x,y> = ZZ[]
            sage: f = P.hom([x+2*y,3*x-y],P)
            sage: F(f)      # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError: Can not create morphisms of free sub-modules yet
        """
        raise NotImplementedError("Can not create morphisms of free sub-modules yet")

    def __cmp__(self, other):
        """
        TEST::

            sage: F1 = (GF(5)^3).span([(1,2,3),(4,5,6)]).construction()[0]
            sage: F2 = (ZZ^3).span([(1,2,3),(4,5,6)]).construction()[0]
            sage: F3 = (QQ^3).span([(1,2,3),(4,5,6)]).construction()[0]
            sage: F4 = (ZZ^3).span([(1,0,-1),(0,1,2)]).construction()[0]
            sage: F1 == loads(dumps(F1))
            True

        The ``span`` method automatically transforms the given basis into
        echelon form. The bases look like that::

            sage: F1.basis
            [
            (1, 0, 4),
            (0, 1, 2)
            ]
            sage: F2.basis
            [
            (1, 2, 3),
            (0, 3, 6)
            ]
            sage: F3.basis
            [
            (1, 0, -1),
            (0, 1, 2)
            ]
            sage: F4.basis
            [
            (1, 0, -1),
            (0, 1, 2)
            ]


        The basis of ``F2`` is modulo 5 different from the other bases.
        So, we have::

            sage: F1 != F2 != F3
            True

        The bases of ``F1``, ``F3`` and ``F4`` are the same modulo 5; however,
        there is no coercion from ``QQ^3`` to ``GF(5)^3``. Therefore, we have::

            sage: F1 == F3
            False

        But there are coercions from ``ZZ^3`` to ``QQ^3`` and ``GF(5)^3``, thus::

            sage: F1 == F4 == F3
            True

        """
        c = cmp(type(self), type(other))
        if c == 0:
            # since comparing the basis involves constructing the pushout
            # of the ambient module, we can not do:
            #c = cmp(self.basis, other.basis)
            # Instead, we only test whether there are coercions.
            L = self.basis.universe()
            R = other.basis.universe()
            c = cmp(L,R)
            if L.has_coerce_map_from(R):
                c = cmp(tuple(self.basis),tuple(L(x) for x in other.basis))
            elif R.has_coerce_map_from(L):
                c = cmp(tuple(other.basis),tuple(R(x) for x in self.basis))
        return c

    def merge(self, other):
        """
        Two Subspace Functors are merged into a construction functor of the sum of two subspaces.

        EXAMPLE::

            sage: M = GF(5)^3
            sage: S1 = M.submodule([(1,2,3),(4,5,6)])
            sage: S2 = M.submodule([(2,2,3)])
            sage: F1 = S1.construction()[0]
            sage: F2 = S2.construction()[0]
            sage: F1.merge(F2)
            SubspaceFunctor
            sage: F1.merge(F2)(GF(5)^3) == S1+S2
            True
            sage: F1.merge(F2)(GF(5)['t']^3)
            Free module of degree 3 and rank 3 over Univariate Polynomial Ring in t over Finite Field of size 5
            User basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        TEST::

            sage: P.<t> = ZZ[]
            sage: S1 = (ZZ^3).submodule([(1,2,3),(4,5,6)])
            sage: S2 = (Frac(P)^3).submodule([(t,t^2,t^3+1),(4*t,0,1)])
            sage: v = S1([0,3,6]) + S2([2,0,1/(2*t)]); v   # indirect doctest
            (2, 3, (-12*t - 1)/(-2*t))
            sage: v.parent()
            Vector space of degree 3 and dimension 3 over Fraction Field of Univariate Polynomial Ring in t over Integer Ring
            User basis matrix:
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        if isinstance(other, SubspaceFunctor):
            # in order to remove linear dependencies, and in
            # order to test compatibility of the base rings,
            # we try to construct a sample submodule
            if not other.basis:
                return self
            if not self.basis:
                return other
            try:
                P = pushout(self.basis[0].parent().ambient_module(),other.basis[0].parent().ambient_module())
            except CoercionException:
                return None
            try:
                # Use span instead of submodule because we want to
                # allow denominators.
                submodule = P.span
            except AttributeError:
                return None
            S = submodule(self.basis+other.basis).echelonized_basis()
            return SubspaceFunctor(S)
        else:
            return None

class FractionField(ConstructionFunctor):
    """
    Construction functor for fraction fields.

    EXAMPLE::

        sage: F = QQ.construction()[0]
        sage: F
        FractionField
        sage: F.domain()
        Category of integral domains
        sage: F.codomain()
        Category of fields
        sage: F(GF(5)) is GF(5)
        True
        sage: F(ZZ['t'])
        Fraction Field of Univariate Polynomial Ring in t over Integer Ring
        sage: P.<x,y> = QQ[]
        sage: f = P.hom([x+2*y,3*x-y],P)
        sage: F(f)
        Ring endomorphism of Fraction Field of Multivariate Polynomial Ring in x, y over Rational Field
          Defn: x |--> x + 2*y
                y |--> 3*x - y
        sage: F(f)(1/x)
        1/(x + 2*y)
        sage: F == loads(dumps(F))
        True

    """
    rank = 5

    def __init__(self):
        """
        TEST::

            sage: from sage.categories.pushout import FractionField
            sage: F = FractionField()
            sage: F
            FractionField
            sage: F(ZZ['t'])
            Fraction Field of Univariate Polynomial Ring in t over Integer Ring
        """
        from sage.categories.integral_domains import IntegralDomains
        from sage.categories.fields import Fields
        Functor.__init__(self, IntegralDomains(), Fields())

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST::

            sage: F = QQ.construction()[0]
            sage: F(GF(5)['t'])      # indirect doctest
            Fraction Field of Univariate Polynomial Ring in t over Finite Field of size 5
        """
        return R.fraction_field()


# This isn't used anywhere in Sage, and so I remove it (Simon King, 2010-05)
#
#class LocalizationFunctor(ConstructionFunctor):
#
#    rank = 6
#
#    def __init__(self, t):
#        Functor.__init__(self, Rings(), Rings())
#        self.t = t
#    def _apply_functor(self, R):
#        return R.localize(t)
#    def __cmp__(self, other):
#        c = cmp(type(self), type(other))
#        if c == 0:
#            c = cmp(self.t, other.t)
#        return c

class CompletionFunctor(ConstructionFunctor):
    """
    Completion of a ring with respect to a given prime (including infinity).

    EXAMPLES::

        sage: R = Zp(5)
        sage: R
        5-adic Ring with capped relative precision 20
        sage: F1 = R.construction()[0]
        sage: F1
        Completion[5]
        sage: F1(ZZ) is R
        True
        sage: F1(QQ)
        5-adic Field with capped relative precision 20
        sage: F2 = RR.construction()[0]
        sage: F2
        Completion[+Infinity]
        sage: F2(QQ) is RR
        True
        sage: P.<x> = ZZ[]
        sage: Px = P.completion(x) # currently the only implemented completion of P
        sage: Px
        Power Series Ring in x over Integer Ring
        sage: F3 = Px.construction()[0]
        sage: F3(GF(3)['x'])
        Power Series Ring in x over Finite Field of size 3

    TEST::

        sage: R1.<a> = Zp(5,prec=20)[]
        sage: R2 = Qp(5,prec=40)
        sage: R2(1) + a
        (1 + O(5^20))*a + (1 + O(5^40))
        sage: 1/2 + a
        (1 + O(5^20))*a + (3 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + 2*5^10 + 2*5^11 + 2*5^12 + 2*5^13 + 2*5^14 + 2*5^15 + 2*5^16 + 2*5^17 + 2*5^18 + 2*5^19 + O(5^20))

    """
    rank = 4

    def __init__(self, p, prec, extras=None):
        """
        INPUT:

        - ``p``: A prime number, the generator of a univariate polynomial ring, or ``+Infinity``

        - ``prec``: an integer, yielding the precision in bits. Note that
          if ``p`` is prime then the ``prec`` is the *capped* precision,
          while it is the *set* precision if ``p`` is ``+Infinity``.

        - ``extras`` (optional dictionary): Information on how to print elements, etc.
          If 'type' is given as a key, the corresponding value should be a string among the following:

          - 'RDF', 'Interval', 'RLF', or 'RR' for completions at infinity

          - 'capped-rel', 'capped-abs', 'fixed-mod' or 'lazy' for completions at a finite place
            or ideal of a DVR.

        TESTS::

            sage: from sage.categories.pushout import CompletionFunctor
            sage: F1 = CompletionFunctor(5,100)
            sage: F1(QQ)
            5-adic Field with capped relative precision 100
            sage: F1(ZZ)
            5-adic Ring with capped relative precision 100
            sage: F2 = RR.construction()[0]
            sage: F2
            Completion[+Infinity]
            sage: F2.extras
            {'rnd': 0, 'sci_not': False, 'type': 'MPFR'}
        """
        Functor.__init__(self, Rings(), Rings())
        self.p = p
        self.prec = prec
        if extras is None:
            self.extras = {}
            self.type = None
        else:
            self.extras = dict(extras)
            self.type = extras.get('type', None)
            from sage.rings.infinity import Infinity
            if self.p == Infinity:
                if self.type not in self._real_types:
                    raise ValueError("completion type must be one of %s"%(", ".join(self._real_types)))
            else:
                if self.type not in self._dvr_types:
                    raise ValueError("completion type must be one of %s"%(", ".join(self._dvr_types)))

    def _repr_(self):
        """
        TEST::

            sage: Zp(7).construction()  # indirect doctest
            (Completion[7], Integer Ring)
        """
        return 'Completion[%s]'%repr(self.p)

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST::

            sage: R = Zp(5)
            sage: F1 = R.construction()[0]
            sage: F1(ZZ) is R  # indirect doctest
            True
            sage: F1(QQ)
            5-adic Field with capped relative precision 20

        """
        try:
            if len(self.extras) == 0:
                if self.type is None:
                    try:
                        return R.completion(self.p, self.prec)
                    except TypeError:
                        return R.completion(self.p, self.prec, {})
                else:
                    return R.completion(self.p, self.prec, {'type':self.type})
            else:
                extras = self.extras.copy()
                extras['type'] = self.type
                return R.completion(self.p, self.prec, extras)
        except (NotImplementedError,AttributeError):
            if R.construction() is None:
                raise NotImplementedError("Completion is not implemented for %s"%R.__class__)
            F, BR = R.construction()
            M = self.merge(F) or F.merge(self)
            if M is not None:
                return M(BR)
            if self.commutes(F) or F.commutes(self):
                return F(self(BR))
            raise NotImplementedError("Don't know how to apply %s to %s"%(repr(self),repr(R)))

    def __cmp__(self, other):
        """
        NOTE:

        Only the prime used in the completion is relevant to comparison
        of Completion functors, although the resulting rings also take
        the precision into account.

        TEST::

            sage: R1 = Zp(5,prec=30)
            sage: R2 = Zp(5,prec=40)
            sage: F1 = R1.construction()[0]
            sage: F2 = R2.construction()[0]
            sage: F1 == loads(dumps(F1))    # indirect doctest
            True
            sage: F1==F2
            True
            sage: F1(QQ)==F2(QQ)
            False
            sage: R3 = Zp(7)
            sage: F3 = R3.construction()[0]
            sage: F1==F3
            False
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.p, other.p)
        return c

    _real_types = ['Interval','Ball','MPFR','RDF','RLF']
    _dvr_types = [None, 'fixed-mod','capped-abs','capped-rel','lazy']

    def merge(self, other):
        """
        Two Completion functors are merged, if they are equal. If the precisions of
        both functors coincide, then a Completion functor is returned that results
        from updating the ``extras`` dictionary of ``self`` by ``other.extras``.
        Otherwise, if the completion is at infinity then merging does not increase
        the set precision, and if the completion is at a finite prime, merging
        does not decrease the capped precision.

        EXAMPLE::

            sage: R1.<a> = Zp(5,prec=20)[]
            sage: R2 = Qp(5,prec=40)
            sage: R2(1)+a         # indirect doctest
            (1 + O(5^20))*a + (1 + O(5^40))
            sage: R3 = RealField(30)
            sage: R4 = RealField(50)
            sage: R3(1) + R4(1)   # indirect doctest
            2.0000000
            sage: (R3(1) + R4(1)).parent()
            Real Field with 30 bits of precision

        TESTS:

        We check that #12353 has been resolved::

            sage: RealIntervalField(53)(-1) > RR(1)
            False
            sage: RealIntervalField(54)(-1) > RR(1)
            False
            sage: RealIntervalField(54)(1) > RR(-1)
            True
            sage: RealIntervalField(53)(1) > RR(-1)
            True

        We check that various pushouts work::

            sage: R0 = RealIntervalField(30)
            sage: R1 = RealIntervalField(30, sci_not=True)
            sage: R2 = RealIntervalField(53)
            sage: R3 = RealIntervalField(53, sci_not = True)
            sage: R4 = RealIntervalField(90)
            sage: R5 = RealIntervalField(90, sci_not = True)
            sage: R6 = RealField(30)
            sage: R7 = RealField(30, sci_not=True)
            sage: R8 = RealField(53, rnd = 'RNDD')
            sage: R9 = RealField(53, sci_not = True, rnd = 'RNDZ')
            sage: R10 = RealField(53, sci_not = True)
            sage: R11 = RealField(90, sci_not = True, rnd = 'RNDZ')
            sage: Rlist = [R0,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11]
            sage: from sage.categories.pushout import pushout
            sage: pushouts = [R0,R0,R0,R1,R0,R1,R0,R1,R0,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R1,R0,R1,R2,R2,R2,R3,R0,R1,R2,R3,R3,R3,R1,R1,R3,R3,R3,R3,R1,R1,R3,R3,R3,R3,R0,R1,R2,R3,R4,R4,R0,R1,R2,R3,R3,R5,R1,R1,R3,R3,R5,R5,R1,R1,R3,R3,R3,R5,R0,R1,R0,R1,R0,R1,R6,R6,R6,R7,R7,R7,R1,R1,R1,R1,R1,R1,R7,R7,R7,R7,R7,R7,R0,R1,R2,R3,R2,R3,R6,R7,R8,R9,R10,R9,R1,R1,R3,R3,R3,R3,R7,R7,R9,R9,R10,R9,R1,R1,R3,R3,R3,R3,R7,R7,R10,R10,R10,R10,R1,R1,R3,R3,R5,R5,R7,R7,R9,R9,R10,R11]
            sage: all([R is S for R, S in zip(pushouts, [pushout(a, b) for a in Rlist for b in Rlist])])
            True

        ::

            sage: P0 = ZpFM(5, 10)
            sage: P1 = ZpFM(5, 20)
            sage: P2 = ZpCR(5, 10)
            sage: P3 = ZpCR(5, 20)
            sage: P4 = ZpCA(5, 10)
            sage: P5 = ZpCA(5, 20)
            sage: P6 = Qp(5, 10)
            sage: P7 = Qp(5, 20)
            sage: Plist = [P2,P3,P4,P5,P6,P7]
            sage: from sage.categories.pushout import pushout
            sage: pushouts = [P2,P3,P4,P5,P6,P7,P3,P3,P5,P5,P7,P7,P4,P5,P4,P5,P6,P7,P5,P5,P5,P5,P7,P7,P6,P7,P6,P7,P6,P7,P7,P7,P7,P7,P7,P7]
            sage: all([P is Q for P, Q in zip(pushouts, [pushout(a, b) for a in Plist for b in Plist])])
            True
        """
        if self == other: # both are Completion functors with the same p
            from sage.all import Infinity
            if self.p == Infinity:
                new_prec = min(self.prec, other.prec)
                new_type = self._real_types[min(self._real_types.index(self.type), \
                                                self._real_types.index(other.type))]
                new_scinot = max(self.extras.get('sci_not',0), other.extras.get('sci_not',0))
                new_rnd = min(self.extras.get('rnd', 0), other.extras.get('rnd', 0))
                return CompletionFunctor(self.p, new_prec, {'type':new_type, 'sci_not':new_scinot, 'rnd':new_rnd})
            else:
                new_type = self._dvr_types[min(self._dvr_types.index(self.type), self._dvr_types.index(other.type))]
                if new_type == 'fixed-mod':
                    if self.type != 'fixed-mod' or other.type != 'fixed-mod':
                        return None # no coercion into fixed-mod
                    new_prec = min(self.prec, other.prec)
                else:
                    new_prec = max(self.prec, other.prec) # since elements track their own precision, we don't want to truncate them
                extras = self.extras.copy()
                extras.update(other.extras)
                extras['type'] = new_type
                return CompletionFunctor(self.p, new_prec, extras)

##   Completion has a lower rank than FractionField
##   and is thus applied first. However, fact is that
##   both commute. This is used in the call method,
##   since some fraction fields have no completion method
##   implemented.

    def commutes(self,other):
        """
        Completion commutes with fraction fields.

        EXAMPLE::

            sage: F1 = Qp(5).construction()[0]
            sage: F2 = QQ.construction()[0]
            sage: F1.commutes(F2)
            True

        TEST:

        The fraction field ``R`` in the example below has no completion
        method. But completion commutes with the fraction field functor,
        and so it is tried internally whether applying the construction
        functors in opposite order works. It does::

            sage: P.<x> = ZZ[]
            sage: C = P.completion(x).construction()[0]
            sage: R = FractionField(P)
            sage: hasattr(R,'completion')
            False
            sage: C(R) is Frac(C(P))
            True
            sage: F = R.construction()[0]
            sage: (C*F)(ZZ['x']) is (F*C)(ZZ['x'])
            True

        The following was fixed in :trac:`15329` (it used to result
        in an infinite recursion)::

            sage: from sage.categories.pushout import pushout
            sage: pushout(Qp(7),RLF)
            Traceback (most recent call last):
            ...
            CoercionException: ('Ambiguous Base Extension', 7-adic Field with capped relative precision 20, Real Lazy Field)

        """
        return isinstance(other,FractionField)

class QuotientFunctor(ConstructionFunctor):
    """
    Construction functor for quotient rings.

    NOTE:

    The functor keeps track of variable names.

    EXAMPLE::

        sage: P.<x,y> = ZZ[]
        sage: Q = P.quo([x^2+y^2]*P)
        sage: F = Q.construction()[0]
        sage: F(QQ['x','y'])
        Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + y^2)
        sage: F(QQ['x','y']) == QQ['x','y'].quo([x^2+y^2]*QQ['x','y'])
        True
        sage: F(QQ['x','y','z'])
        Traceback (most recent call last):
        ...
        CoercionException: Can not apply this quotient functor to Multivariate Polynomial Ring in x, y, z over Rational Field
        sage: F(QQ['y','z'])
        Traceback (most recent call last):
        ...
        TypeError: Could not find a mapping of the passed element to this ring.
    """
    rank = 4.5

    def __init__(self, I, names=None, as_field=False):
        """
        INPUT:

        - ``I``, an ideal (the modulus)
        - ``names`` (optional string or list of strings), the names for the quotient ring generators
        - ``as_field`` (optional bool, default false), return the quotient ring as field (if available).

        TESTS::

            sage: from sage.categories.pushout import QuotientFunctor
            sage: P.<t> = ZZ[]
            sage: F = QuotientFunctor([5+t^2]*P)
            sage: F(P)
            Univariate Quotient Polynomial Ring in tbar over Integer Ring with modulus t^2 + 5
            sage: F(QQ['t'])
            Univariate Quotient Polynomial Ring in tbar over Rational Field with modulus t^2 + 5
            sage: F = QuotientFunctor([5+t^2]*P,names='s')
            sage: F(P)
            Univariate Quotient Polynomial Ring in s over Integer Ring with modulus t^2 + 5
            sage: F(QQ['t'])
            Univariate Quotient Polynomial Ring in s over Rational Field with modulus t^2 + 5
            sage: F = QuotientFunctor([5]*ZZ,as_field=True)
            sage: F(ZZ)
            Finite Field of size 5
            sage: F = QuotientFunctor([5]*ZZ)
            sage: F(ZZ)
            Ring of integers modulo 5

        """
        Functor.__init__(self, Rings(), Rings()) # much more general...
        self.I = I
        if names is None:
            self.names = None
        elif isinstance(names, six.string_types):
            self.names = (names,)
        else:
            self.names = tuple(names)
        self.as_field = as_field

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: P.<x,y> = ZZ[]
            sage: Q = P.quo([2+x^2,3*x+y^2])
            sage: F = Q.construction()[0]; F
            QuotientFunctor
            sage: F(QQ['x','y'])     # indirect doctest
            Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x^2 + 2, y^2 + 3*x)

        Note that the ``quo()`` method of a field used to return the
        integer zero. That strange behaviour was removed in trac
        ticket :trac:`9138`. It now returns a trivial quotient ring
        when applied to a field::

            sage: F = ZZ.quo([5]*ZZ).construction()[0]
            sage: F(QQ)
            Ring of integers modulo 1
            sage: QQ.quo(5)
            Quotient of Rational Field by the ideal (1)
        """
        I = self.I
        from sage.all import QQ
        if not I.is_zero():
            from sage.categories.fields import Fields
            if R in Fields():
                from sage.all import Integers
                return Integers(1)
        if I.ring() != R:
            if I.ring().has_coerce_map_from(R):
                R = I.ring()
            else:
                R = pushout(R,I.ring().base_ring())
                I = [R(1)*t for t in I.gens()]*R
        try:
            Q = R.quo(I,names=self.names)
        except IndexError: # That may happen!
            raise CoercionException("Can not apply this quotient functor to %s"%R)
        if self.as_field:# and hasattr(Q, 'field'):
            try:
                Q = Q.field()
            except AttributeError:
                pass
        return Q

    def __cmp__(self, other):
        """
        The types, the names and the moduli are compared.

        TESTS::

            sage: P.<x> = QQ[]
            sage: F = P.quo([(x^2+1)^2*(x^2-3),(x^2+1)^2*(x^5+3)]).construction()[0]
            sage: F == loads(dumps(F))
            True
            sage: P2.<x,y> = QQ[]
            sage: F == P2.quo([(x^2+1)^2*(x^2-3),(x^2+1)^2*(x^5+3)]).construction()[0]
            False
            sage: P3.<x> = ZZ[]
            sage: F == P3.quo([(x^2+1)^2*(x^2-3),(x^2+1)^2*(x^5+3)]).construction()[0]
            True
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.names, other.names)
        if c == 0:
            c = cmp(self.I, other.I)
        return c

    def merge(self, other):
        """
        Two quotient functors with coinciding names are merged by taking the gcd of their moduli.

        EXAMPLE::

            sage: P.<x> = QQ[]
            sage: Q1 = P.quo([(x^2+1)^2*(x^2-3)])
            sage: Q2 = P.quo([(x^2+1)^2*(x^5+3)])
            sage: from sage.categories.pushout import pushout
            sage: pushout(Q1,Q2)    # indirect doctest
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^4 + 2*x^2 + 1

        The following was fixed in :trac:`8800`::

            sage: pushout(GF(5), Integers(5))
            Finite Field of size 5

        """
        if type(self) is not type(other):
            return None
        if self.names != other.names:
            return None
        if self == other:
            if self.as_field == other.as_field:
                return self
            return QuotientFunctor(self.I, names=self.names, as_field=True) # one of them yields a field!
        try:
            gcd = self.I + other.I
        except (TypeError, NotImplementedError):
            try:
                gcd = self.I.gcd(other.I)
            except (TypeError, NotImplementedError):
                return None
        if gcd.is_trivial() and not gcd.is_zero():
            # quotient by gcd would result in the trivial ring/group/...
            # Rather than create the zero ring, we claim they can't be merged
            # TODO: Perhaps this should be detected at a higher level...
            raise TypeError("Trivial quotient intersection.")
        # GF(p) has a coercion from Integers(p). Hence, merging should
        # yield a field if either self or other yields a field.
        return QuotientFunctor(gcd, names=self.names, as_field=self.as_field or other.as_field)

class AlgebraicExtensionFunctor(ConstructionFunctor):
    """
    Algebraic extension (univariate polynomial ring modulo principal ideal).

    EXAMPLE::

        sage: K.<a> = NumberField(x^3+x^2+1)
        sage: F = K.construction()[0]
        sage: F(ZZ['t'])
        Univariate Quotient Polynomial Ring in a over Univariate Polynomial Ring in t over Integer Ring with modulus a^3 + a^2 + 1

    Note that, even if a field is algebraically closed, the algebraic
    extension will be constructed as the quotient of a univariate
    polynomial ring::

        sage: F(CC)
        Univariate Quotient Polynomial Ring in a over Complex Field with 53 bits of precision with modulus a^3 + a^2 + 1.00000000000000
        sage: F(RR)
        Univariate Quotient Polynomial Ring in a over Real Field with 53 bits of precision with modulus a^3 + a^2 + 1.00000000000000

    Note that the construction functor of a number field applied to
    the integers returns an order (not necessarily maximal) of that
    field, similar to the behaviour of ``ZZ.extension(...)``::

        sage: F(ZZ)
        Order in Number Field in a with defining polynomial x^3 + x^2 + 1

    This also holds for non-absolute number fields::

        sage: K.<a,b> = NumberField([x^3+x^2+1,x^2+x+1])
        sage: F = K.construction()[0]
        sage: O = F(ZZ); O
        Relative Order in Number Field in a with defining polynomial x^3 + x^2 + 1 over its base field

    Unfortunately, the relative number field is not a unique parent::

        sage: O.ambient() is K
        False
        sage: O.ambient() == K
        True

    """
    rank = 3

    def __init__(self, polys, names, embeddings, cyclotomic=None, **kwds):
        """
        INPUT:

        - ``polys``: a list of polynomials (or of integers, for
          finite fields and unramified local extensions)
        - ``names``: a list of strings of the same length as the
          list ``polys``
        - ``embeddings``: a list of approximate complex values,
          determining an embedding of the generators into the
          complex field, or ``None`` for each generator whose
          embedding is not prescribed.
        - ``cyclotomic``: optional integer. If it is provided,
          application of the functor to the rational field yields
          a cyclotomic field, rather than just a number field.
        - ``**kwds``: further keywords; when the functor is applied to
          a ring `R`, these are passed to the ``extension()`` method
          of `R`.

        REMARK:

        Currently, an embedding can only be provided for the last
        generator, and only when the construction functor is applied
        to the rational field. There will be no error when constructing
        the functor, but when applying it.

        TESTS::

            sage: from sage.categories.pushout import AlgebraicExtensionFunctor
            sage: P.<x> = ZZ[]
            sage: F1 = AlgebraicExtensionFunctor([x^3 - x^2 + 1], ['a'], [None])
            sage: F2 = AlgebraicExtensionFunctor([x^3 - x^2 + 1], ['a'], [0])
            sage: F1==F2
            False
            sage: F1(QQ)
            Number Field in a with defining polynomial x^3 - x^2 + 1
            sage: F1(QQ).coerce_embedding()
            sage: phi = F2(QQ).coerce_embedding().__copy__(); phi
            Generic morphism:
              From: Number Field in a with defining polynomial x^3 - x^2 + 1
              To:   Real Lazy Field
              Defn: a -> -0.7548776662466928?
            sage: F1(QQ)==F2(QQ)
            False
            sage: F1(GF(5))
            Univariate Quotient Polynomial Ring in a over Finite Field of size 5 with modulus a^3 + 4*a^2 + 1
            sage: F2(GF(5))
            Traceback (most recent call last):
            ...
            NotImplementedError: ring extension with prescripted embedding is not implemented

        When applying a number field constructor to the ring of
        integers, an order (not necessarily maximal) of that field is
        returned, similar to the behaviour of ``ZZ.extension``::

            sage: F1(ZZ)
            Order in Number Field in a with defining polynomial x^3 - x^2 + 1

        The cyclotomic fields form a special case of number fields
        with prescribed embeddings::

            sage: C = CyclotomicField(8)
            sage: F,R = C.construction()
            sage: F
            AlgebraicExtensionFunctor
            sage: R
            Rational Field
            sage: F(R)
            Cyclotomic Field of order 8 and degree 4
            sage: F(ZZ)
            Maximal Order in Cyclotomic Field of order 8 and degree 4

        """
        Functor.__init__(self, Rings(), Rings())
        if not (isinstance(polys,(list,tuple)) and isinstance(names,(list,tuple)) and isinstance(embeddings,(list,tuple))):
            raise ValueError("Arguments must be lists or tuples")
        if not (len(names)==len(polys)==len(embeddings)):
            raise ValueError("The three arguments must be of the same length")
        self.polys = list(polys)
        self.names = list(names)
        self.embeddings = list(embeddings)
        self.cyclotomic = int(cyclotomic) if cyclotomic is not None else None
        self.kwds = kwds

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: K.<a>=NumberField(x^3+x^2+1)
            sage: F = K.construction()[0]
            sage: F(ZZ)       # indirect doctest
            Order in Number Field in a with defining polynomial x^3 + x^2 + 1
            sage: F(ZZ['t'])  # indirect doctest
            Univariate Quotient Polynomial Ring in a over Univariate Polynomial Ring in t over Integer Ring with modulus a^3 + a^2 + 1
            sage: F(RR)       # indirect doctest
            Univariate Quotient Polynomial Ring in a over Real Field with 53 bits of precision with modulus a^3 + a^2 + 1.00000000000000

        Check that :trac:`13538` is fixed::

            sage: K = Qp(3,3)
            sage: R.<a> = K[]
            sage: AEF = sage.categories.pushout.AlgebraicExtensionFunctor([a^2-3], ['a'], [None])
            sage: AEF(K)
            Eisenstein Extension of 3-adic Field with capped relative precision 3 in a defined by (1 + O(3^3))*a^2 + (O(3^4))*a + (2*3 + 2*3^2 + 2*3^3 + O(3^4))

        """
        from sage.all import QQ, ZZ, CyclotomicField
        if self.cyclotomic:
            if R==QQ:
                return CyclotomicField(self.cyclotomic)
            if R==ZZ:
                return CyclotomicField(self.cyclotomic).maximal_order()
        if len(self.polys) == 1:
            return R.extension(self.polys[0], names=self.names[0], embedding=self.embeddings[0], **self.kwds)
        return R.extension(self.polys, names=self.names, embedding=self.embeddings)

    def __cmp__(self, other):
        """
        TEST::

            sage: K.<a>=NumberField(x^3+x^2+1)
            sage: F = K.construction()[0]
            sage: F == loads(dumps(F))
            True
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.polys, other.polys)
        if c == 0:
            c = cmp(self.embeddings, other.embeddings)
        return c

    def merge(self,other):
        """
        Merging with another :class:`AlgebraicExtensionFunctor`.

        INPUT:

        ``other`` -- Construction Functor.

        OUTPUT:

        - If ``self==other``, ``self`` is returned.
        - If ``self`` and ``other`` are simple extensions
          and both provide an embedding, then it is tested
          whether one of the number fields provided by
          the functors coerces into the other; the functor
          associated with the target of the coercion is
          returned. Otherwise, the construction functor
          associated with the pushout of the codomains
          of the two embeddings is returned, provided that
          it is a number field.
        - If these two extensions are defined by Conway polynomials
          over finite fields, merges them into a single extension of
          degree the lcm of the two degrees.
        - Otherwise, None is returned.

        REMARK:

        Algebraic extension with embeddings currently only
        works when applied to the rational field. This is
        why we use the admittedly strange rule above for
        merging.

        EXAMPLES:

        The following demonstrate coercions for finite fields using Conway or
        pseudo-Conway polynomials::

            sage: k = GF(3^2, prefix='z'); a = k.gen()
            sage: l = GF(3^3, prefix='z'); b = l.gen()
            sage: a + b # indirect doctest
            z6^5 + 2*z6^4 + 2*z6^3 + z6^2 + 2*z6 + 1

        Note that embeddings are compatible in lattices of such finite fields::

            sage: m = GF(3^5, prefix='z'); c = m.gen()
            sage: (a+b)+c == a+(b+c) # indirect doctest
            True
            sage: from sage.categories.pushout import pushout
            sage: n = pushout(k, l)
            sage: o = pushout(l, m)
            sage: q = pushout(n, o)
            sage: q(o(b)) == q(n(b)) # indirect doctest
            True

        Coercion is also available for number fields::

            sage: P.<x> = QQ[]
            sage: L.<b> = NumberField(x^8-x^4+1, embedding=CDF.0)
            sage: M1.<c1> = NumberField(x^2+x+1, embedding=b^4-1)
            sage: M2.<c2> = NumberField(x^2+1, embedding=-b^6)
            sage: M1.coerce_map_from(M2)
            sage: M2.coerce_map_from(M1)
            sage: c1+c2; parent(c1+c2)    #indirect doctest
            -b^6 + b^4 - 1
            Number Field in b with defining polynomial x^8 - x^4 + 1
            sage: pushout(M1['x'],M2['x'])
            Univariate Polynomial Ring in x over Number Field in b with defining polynomial x^8 - x^4 + 1

        In the previous example, the number field ``L`` becomes the pushout
        of ``M1`` and ``M2`` since both are provided with an embedding into
        ``L``, *and* since ``L`` is a number field. If two number fields
        are embedded into a field that is not a numberfield, no merging
        occurs::

            sage: K.<a> = NumberField(x^3-2, embedding=CDF(1/2*I*2^(1/3)*sqrt(3) - 1/2*2^(1/3)))
            sage: L.<b> = NumberField(x^6-2, embedding=1.1)
            sage: L.coerce_map_from(K)
            sage: K.coerce_map_from(L)
            sage: pushout(K,L)
            Traceback (most recent call last):
            ...
            CoercionException: ('Ambiguous Base Extension', Number Field in a with defining polynomial x^3 - 2, Number Field in b with defining polynomial x^6 - 2)

        """
        if isinstance(other, AlgebraicClosureFunctor):
            return other
        elif not isinstance(other, AlgebraicExtensionFunctor):
            return None
        if self == other:
            return self
        # This method is supposed to be used in pushout(),
        # *after* expanding the functors. Hence, we can
        # assume that both functors have a single variable.
        # But for being on the safe side...:
        if len(self.names)!=1 or len(other.names)!=1:
            return None
##       We don't accept a forgetful coercion, since, together
##       with bidirectional coercions between two embedded
##       number fields, it would yield to contradictions in
##       the coercion system.
#        if self.polys==other.polys and self.names==other.names:
#            # We have a forgetful functor:
#            if self.embeddings==[None]:
#                return self
#            if  other.embeddings==[None]:
#                return other
        # ... or we may use the given embeddings:
        if self.embeddings!=[None] and other.embeddings!=[None]:
            from sage.all import QQ
            KS = self(QQ)
            KO = other(QQ)
            if KS.has_coerce_map_from(KO):
                return self
            if KO.has_coerce_map_from(KS):
                return other
            # nothing else helps, hence, we move to the pushout of the codomains of the embeddings
            try:
                P = pushout(self.embeddings[0].parent(), other.embeddings[0].parent())
                from sage.rings.number_field.number_field import is_NumberField
                if is_NumberField(P):
                    return P.construction()[0]
            except CoercionException:
                return None
        # Finite fields and unramified local extensions may use
        # integers to encode degrees of extensions.
        from sage.rings.integer import Integer
        if (isinstance(self.polys[0], Integer) and isinstance(other.polys[0], Integer)
            and self.embeddings == [None] and other.embeddings == [None] and self.kwds == other.kwds):
            return AlgebraicExtensionFunctor([self.polys[0].lcm(other.polys[0])], [None], [None], **self.kwds)

    def __mul__(self, other):
        """
        Compose construction functors to a composit construction functor, unless one of them is the identity.

        NOTE:

        The product is in functorial notation, i.e., when applying the product to an object
        then the second factor is applied first.

        TESTS::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-5,embedding=0)
            sage: L.<b> = K.extension(x^2+a)
            sage: F,R = L.construction()
            sage: prod(F.expand())(R) == L #indirect doctest
            True

        """
        if isinstance(other,IdentityConstructionFunctor):
            return self
        if isinstance(other, AlgebraicExtensionFunctor):
            if set(self.names).intersection(other.names):
                raise CoercionException("Overlapping names (%s,%s)" % (self.names, other.names))
            return AlgebraicExtensionFunctor(self.polys + other.polys, self.names + other.names,
                                             self.embeddings + other.embeddings, **self.kwds)
        elif isinstance(other, CompositeConstructionFunctor) \
              and isinstance(other.all[-1], AlgebraicExtensionFunctor):
            return CompositeConstructionFunctor(other.all[:-1], self * other.all[-1])
        else:
            return CompositeConstructionFunctor(other, self)

    def expand(self):
        """
        Decompose the functor `F` into sub-functors, whose product returns `F`.

        EXAMPLES::

            sage: P.<x> = QQ[]
            sage: K.<a> = NumberField(x^3-5,embedding=0)
            sage: L.<b> = K.extension(x^2+a)
            sage: F,R = L.construction()
            sage: prod(F.expand())(R) == L
            True
            sage: K = NumberField([x^2-2, x^2-3],'a')
            sage: F, R = K.construction()
            sage: F
            AlgebraicExtensionFunctor
            sage: L = F.expand(); L
            [AlgebraicExtensionFunctor, AlgebraicExtensionFunctor]
            sage: L[-1](QQ)
            Number Field in a1 with defining polynomial x^2 - 3
        """
        if len(self.polys)==1:
            return [self]
        return [AlgebraicExtensionFunctor([self.polys[i]], [self.names[i]], [self.embeddings[i]], **self.kwds)
                for i in xrange(len(self.polys))]

class AlgebraicClosureFunctor(ConstructionFunctor):
    """
    Algebraic Closure.

    EXAMPLE::

        sage: F = CDF.construction()[0]
        sage: F(QQ)
        Algebraic Field
        sage: F(RR)
        Complex Field with 53 bits of precision
        sage: F(F(QQ)) is F(QQ)
        True

    """
    rank = 3

    def __init__(self):
        """
        TEST::

            sage: from sage.categories.pushout import AlgebraicClosureFunctor
            sage: F = AlgebraicClosureFunctor()
            sage: F(QQ)
            Algebraic Field
            sage: F(RR)
            Complex Field with 53 bits of precision
            sage: F == loads(dumps(F))
            True

        """
        Functor.__init__(self, Rings(), Rings())

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TEST::

            sage: F = CDF.construction()[0]
            sage: F(QQ)       # indirect doctest
            Algebraic Field
        """
        try:
            c = R.construction()
            if c is not None and c[0]==self:
                return R
        except AttributeError:
            pass
        return R.algebraic_closure()

    def merge(self, other):
        """
        Mathematically, Algebraic Closure subsumes Algebraic Extension.
        However, it seems that people do want to work with algebraic
        extensions of ``RR``. Therefore, we do not merge with algebraic extension.

        TEST::

            sage: K.<a>=NumberField(x^3+x^2+1)
            sage: CDF.construction()[0].merge(K.construction()[0]) is None
            True
            sage: CDF.construction()[0].merge(CDF.construction()[0])
            AlgebraicClosureFunctor

        """
        if self==other:
            return self
        return None
        # Mathematically, Algebraic Closure subsumes Algebraic Extension.
        # However, it seems that people do want to work with
        # algebraic extensions of RR (namely RR/poly*RR). So, we don't do:
        # if isinstance(other,AlgebraicExtensionFunctor):
        #     return self

class PermutationGroupFunctor(ConstructionFunctor):

    rank = 10

    def __init__(self, gens, domain):
        """
        EXAMPLES::

            sage: from sage.categories.pushout import PermutationGroupFunctor
            sage: PF = PermutationGroupFunctor([PermutationGroupElement([(1,2)])], [1,2]); PF
            PermutationGroupFunctor[(1,2)]
        """
        Functor.__init__(self, Groups(), Groups())
        self._gens = gens
        self._domain = domain

    def _repr_(self):
        """
        EXAMPLES::

            sage: P1 = PermutationGroup([[(1,2)]])
            sage: PF, P = P1.construction()
            sage: PF
            PermutationGroupFunctor[(1,2)]
        """
        return "PermutationGroupFunctor%s"%self.gens()

    def __call__(self, R):
        """
        EXAMPLES::

            sage: P1 = PermutationGroup([[(1,2)]])
            sage: PF, P = P1.construction()
            sage: PF(P)
            Permutation Group with generators [(1,2)]
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        return PermutationGroup([g for g in (R.gens() + self.gens()) if not g.is_one()],
                                domain=self._domain)

    def gens(self):
        """
        EXAMPLES::

            sage: P1 = PermutationGroup([[(1,2)]])
            sage: PF, P = P1.construction()
            sage: PF.gens()
            [(1,2)]
        """
        return self._gens

    def merge(self, other):
        """
        Merge ``self`` with another construction functor, or return None.

        EXAMPLES::

            sage: P1 = PermutationGroup([[(1,2)]])
            sage: PF1, P = P1.construction()
            sage: P2 = PermutationGroup([[(1,3)]])
            sage: PF2, P = P2.construction()
            sage: PF1.merge(PF2)
            PermutationGroupFunctor[(1,2), (1,3)]
        """
        if self.__class__ != other.__class__:
            return None
        from sage.sets.all import FiniteEnumeratedSet

        new_domain = set(self._domain).union(set(other._domain))
        new_domain = FiniteEnumeratedSet(sorted(new_domain))
        return PermutationGroupFunctor(self.gens() + other.gens(),
                                       new_domain)

class BlackBoxConstructionFunctor(ConstructionFunctor):
    """
    Construction functor obtained from any callable object.

    EXAMPLES::

        sage: from sage.categories.pushout import BlackBoxConstructionFunctor
        sage: FG = BlackBoxConstructionFunctor(gap)
        sage: FS = BlackBoxConstructionFunctor(singular)
        sage: FG
        BlackBoxConstructionFunctor
        sage: FG(ZZ)
        Integers
        sage: FG(ZZ).parent()
        Gap
        sage: FS(QQ['t'])
        //   characteristic : 0
        //   number of vars : 1
        //        block   1 : ordering lp
        //                  : names    t
        //        block   2 : ordering C
        sage: FG == FS
        False
        sage: FG == loads(dumps(FG))
        True
    """
    rank = 100

    def __init__(self, box):
        """
        TESTS::

            sage: from sage.categories.pushout import BlackBoxConstructionFunctor
            sage: FG = BlackBoxConstructionFunctor(gap)
            sage: FM = BlackBoxConstructionFunctor(maxima)
            sage: FM == FG
            False
            sage: FM == loads(dumps(FM))
            True
        """
        ConstructionFunctor.__init__(self,Objects(),Objects())
        if not callable(box):
            raise TypeError("input must be callable")
        self.box = box

    def _apply_functor(self, R):
        """
        Apply the functor to an object of ``self``'s domain.

        TESTS::

            sage: from sage.categories.pushout import BlackBoxConstructionFunctor
            sage: f = lambda x: x^2
            sage: F = BlackBoxConstructionFunctor(f)
            sage: F(ZZ)           # indirect doctest
            Ambient free module of rank 2 over the principal ideal domain Integer Ring

        """
        return self.box(R)

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.categories.pushout import BlackBoxConstructionFunctor
            sage: FG = BlackBoxConstructionFunctor(gap)
            sage: FM = BlackBoxConstructionFunctor(maxima)
            sage: FM == FG       # indirect doctest
            False
            sage: FM == loads(dumps(FM))
            True
        """
        c = cmp(type(self), type(other))
        if c == 0:
            c = cmp(self.box, other.box)
        #return self.box == other.box
        return c

def pushout(R, S):
    r"""
    Given a pair of objects `R` and `S`, try to construct a
    reasonable object `Y` and return maps such that
    canonically `R \leftarrow Y \rightarrow S`.

    ALGORITHM:

    This incorporates the idea of functors discussed at Sage Days 4.
    Every object `R` can be viewed as an initial object and a series
    of functors (e.g. polynomial, quotient, extension, completion,
    vector/matrix, etc.). Call the series of increasingly simple
    objects (with the associated functors) the "tower" of `R`. The
    construction method is used to create the tower.

    Given two objects `R` and `S`, try to find a common initial object
    `Z`. If the towers of `R` and `S` meet, let `Z` be their join.
    Otherwise, see if the top of one coerces naturally into the other.

    Now we have an initial object and two ordered lists of functors to
    apply. We wish to merge these in an unambiguous order, popping
    elements off the top of one or the other tower as we apply them to
    `Z`.

    - If the functors are of distinct types, there is an absolute
      ordering given by the rank attribute. Use this.

    - Otherwise:

      - If the tops are equal, we (try to) merge them.

      - If exactly one occurs lower in the other tower, we may
        unambiguously apply the other (hoping for a later merge).

      - If the tops commute, we can apply either first.

      - Otherwise fail due to ambiguity.

    The algorithm assumes by default that when a construction `F` is
    applied to an object `X`, the object `F(X)` admits a coercion map
    from `X`.  However, the algorithm can also handle the case where
    `F(X)` has a coercion map *to* `X` instead.  In this case, the
    attribute ``coercion_reversed`` of the class implementing `F`
    should be set to ``True``.

    EXAMPLES:

    Here our "towers" are `R = Complete_7(Frac(\ZZ))` and `Frac(Poly_x(\ZZ))`,
    which give us `Frac(Poly_x(Complete_7(Frac(\ZZ))))`::

        sage: from sage.categories.pushout import pushout
        sage: pushout(Qp(7), Frac(ZZ['x']))
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20

    Note we get the same thing with
    ::

        sage: pushout(Zp(7), Frac(QQ['x']))
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20
        sage: pushout(Zp(7)['x'], Frac(QQ['x']))
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20

    Note that polynomial variable ordering must be unambiguously determined.
    ::

        sage: pushout(ZZ['x,y,z'], QQ['w,z,t'])
        Traceback (most recent call last):
        ...
        CoercionException: ('Ambiguous Base Extension', Multivariate Polynomial Ring in x, y, z over Integer Ring, Multivariate Polynomial Ring in w, z, t over Rational Field)
        sage: pushout(ZZ['x,y,z'], QQ['w,x,z,t'])
        Multivariate Polynomial Ring in w, x, y, z, t over Rational Field

    Some other examples::

        sage: pushout(Zp(7)['y'], Frac(QQ['t'])['x,y,z'])
        Multivariate Polynomial Ring in x, y, z over Fraction Field of Univariate Polynomial Ring in t over 7-adic Field with capped relative precision 20
        sage: pushout(ZZ['x,y,z'], Frac(ZZ['x'])['y'])
        Multivariate Polynomial Ring in y, z over Fraction Field of Univariate Polynomial Ring in x over Integer Ring
        sage: pushout(MatrixSpace(RDF, 2, 2), Frac(ZZ['x']))
        Full MatrixSpace of 2 by 2 dense matrices over Fraction Field of Univariate Polynomial Ring in x over Real Double Field
        sage: pushout(ZZ, MatrixSpace(ZZ[['x']], 3, 3))
        Full MatrixSpace of 3 by 3 dense matrices over Power Series Ring in x over Integer Ring
        sage: pushout(QQ['x,y'], ZZ[['x']])
        Univariate Polynomial Ring in y over Power Series Ring in x over Rational Field
        sage: pushout(Frac(ZZ['x']), QQ[['x']])
        Laurent Series Ring in x over Rational Field

    A construction with ``coercion_reversed = True`` (currently only
    the :class:`SubspaceFunctor` construction) is only applied if it
    leads to a valid coercion::

        sage: A = ZZ^2
        sage: V = span([[1, 2]], QQ)
        sage: P = sage.categories.pushout.pushout(A, V)
        sage: P
        Vector space of dimension 2 over Rational Field
        sage: P.has_coerce_map_from(A)
        True

        sage: V = (QQ^3).span([[1, 2, 3/4]])
        sage: A = ZZ^3
        sage: pushout(A, V)
        Vector space of dimension 3 over Rational Field
        sage: B = A.span([[0, 0, 2/3]])
        sage: pushout(B, V)
        Vector space of degree 3 and dimension 2 over Rational Field
        User basis matrix:
        [1 2 0]
        [0 0 1]

    Some more tests with ``coercion_reversed = True``::

        sage: from sage.categories.pushout import ConstructionFunctor
        sage: class EvenPolynomialRing(type(QQ['x'])):
        ....:     def __init__(self, base, var):
        ....:         super(EvenPolynomialRing, self).__init__(base, var)
        ....:         self.register_embedding(base[var])
        ....:     def __repr__(self):
        ....:         return "Even Power " + super(EvenPolynomialRing, self).__repr__()
        ....:     def construction(self):
        ....:         return EvenPolynomialFunctor(), self.base()[self.variable_name()]
        ....:     def _coerce_map_from_(self, R):
        ....:         return self.base().has_coerce_map_from(R)
        ....:
        sage: class EvenPolynomialFunctor(ConstructionFunctor):
        ....:     rank = 10
        ....:     coercion_reversed = True
        ....:     def __init__(self):
        ....:         ConstructionFunctor.__init__(self, Rings(), Rings())
        ....:     def _apply_functor(self, R):
        ....:         return EvenPolynomialRing(R.base(), R.variable_name())
        ....:
        sage: pushout(EvenPolynomialRing(QQ, 'x'), ZZ)
        Even Power Univariate Polynomial Ring in x over Rational Field
        sage: pushout(EvenPolynomialRing(QQ, 'x'), QQ)
        Even Power Univariate Polynomial Ring in x over Rational Field
        sage: pushout(EvenPolynomialRing(QQ, 'x'), RR)
        Even Power Univariate Polynomial Ring in x over Real Field with 53 bits of precision

        sage: pushout(EvenPolynomialRing(QQ, 'x'), ZZ['x'])
        Univariate Polynomial Ring in x over Rational Field
        sage: pushout(EvenPolynomialRing(QQ, 'x'), QQ['x'])
        Univariate Polynomial Ring in x over Rational Field
        sage: pushout(EvenPolynomialRing(QQ, 'x'), RR['x'])
        Univariate Polynomial Ring in x over Real Field with 53 bits of precision

        sage: pushout(EvenPolynomialRing(QQ, 'x'), EvenPolynomialRing(QQ, 'x'))
        Even Power Univariate Polynomial Ring in x over Rational Field
        sage: pushout(EvenPolynomialRing(QQ, 'x'), EvenPolynomialRing(RR, 'x'))
        Even Power Univariate Polynomial Ring in x over Real Field with 53 bits of precision

        sage: pushout(EvenPolynomialRing(QQ, 'x')^2, RR^2)
        Ambient free module of rank 2 over the principal ideal domain Even Power Univariate Polynomial Ring in x over Real Field with 53 bits of precision
        sage: pushout(EvenPolynomialRing(QQ, 'x')^2, RR['x']^2)
        Ambient free module of rank 2 over the principal ideal domain Univariate Polynomial Ring in x over Real Field with 53 bits of precision

    Some more tests related to univariate/multivariate
    constructions. We consider a generalization of polynomial rings,
    where in addition to the coefficient ring `C` we also specify
    an additive monoid `E` for the exponents of the indeterminate.
    In particular, the elements of such a parent are given by

    .. MATH::

        \sum_{i=0}^I c_i X^{e_i}

    with `c_i \in C` and `e_i \in E`. We define
    ::

        sage: class GPolynomialRing(Parent):
        ....:     def __init__(self, coefficients, var, exponents):
        ....:         self.coefficients = coefficients
        ....:         self.var = var
        ....:         self.exponents = exponents
        ....:         super(GPolynomialRing, self).__init__(category=Rings())
        ....:     def _repr_(self):
        ....:         return 'Generalized Polynomial Ring in %s^(%s) over %s' % (
        ....:                self.var, self.exponents, self.coefficients)
        ....:     def construction(self):
        ....:         return GPolynomialFunctor(self.var, self.exponents), self.coefficients
        ....:     def _coerce_map_from_(self, R):
        ....:         return self.coefficients.has_coerce_map_from(R)

    and
    ::

        sage: class GPolynomialFunctor(ConstructionFunctor):
        ....:     rank = 10
        ....:     def __init__(self, var, exponents):
        ....:         self.var = var
        ....:         self.exponents = exponents
        ....:         ConstructionFunctor.__init__(self, Rings(), Rings())
        ....:     def _repr_(self):
        ....:         return 'GPoly[%s^(%s)]' % (self.var, self.exponents)
        ....:     def _apply_functor(self, coefficients):
        ....:         return GPolynomialRing(coefficients, self.var, self.exponents)
        ....:     def merge(self, other):
        ....:         if isinstance(other, GPolynomialFunctor) and self.var == other.var:
        ....:             exponents = pushout(self.exponents, other.exponents)
        ....:             return GPolynomialFunctor(self.var, exponents)

    We can construct a parent now in two different ways::

        sage: GPolynomialRing(QQ, 'X', ZZ)
        Generalized Polynomial Ring in X^(Integer Ring) over Rational Field
        sage: GP_ZZ = GPolynomialFunctor('X', ZZ); GP_ZZ
        GPoly[X^(Integer Ring)]
        sage: GP_ZZ(QQ)
        Generalized Polynomial Ring in X^(Integer Ring) over Rational Field

    Since the construction
    ::

        sage: GP_ZZ(QQ).construction()
        (GPoly[X^(Integer Ring)], Rational Field)

    uses the coefficient ring, we have the usual coercion with respect
    to this parameter::

        sage: pushout(GP_ZZ(ZZ), GP_ZZ(QQ))
        Generalized Polynomial Ring in X^(Integer Ring) over Rational Field
        sage: pushout(GP_ZZ(ZZ['t']), GP_ZZ(QQ))
        Generalized Polynomial Ring in X^(Integer Ring) over Univariate Polynomial Ring in t over Rational Field
        sage: pushout(GP_ZZ(ZZ['a,b']), GP_ZZ(ZZ['b,c']))
        Generalized Polynomial Ring in X^(Integer Ring)
          over Multivariate Polynomial Ring in a, b, c over Integer Ring
        sage: pushout(GP_ZZ(ZZ['a,b']), GP_ZZ(QQ['b,c']))
        Generalized Polynomial Ring in X^(Integer Ring)
          over Multivariate Polynomial Ring in a, b, c over Rational Field
        sage: pushout(GP_ZZ(ZZ['a,b']), GP_ZZ(ZZ['c,d']))
        Traceback (most recent call last):
        ...
        CoercionException: ('Ambiguous Base Extension', ...)

    ::

        sage: GP_QQ = GPolynomialFunctor('X', QQ)
        sage: pushout(GP_ZZ(ZZ), GP_QQ(ZZ))
        Generalized Polynomial Ring in X^(Rational Field) over Integer Ring
        sage: pushout(GP_QQ(ZZ), GP_ZZ(ZZ))
        Generalized Polynomial Ring in X^(Rational Field) over Integer Ring

    ::

        sage: GP_ZZt = GPolynomialFunctor('X', ZZ['t'])
        sage: pushout(GP_ZZt(ZZ), GP_QQ(ZZ))
        Generalized Polynomial Ring in X^(Univariate Polynomial Ring in t
          over Rational Field) over Integer Ring

    ::

        sage: pushout(GP_ZZ(ZZ), GP_QQ(QQ))
        Generalized Polynomial Ring in X^(Rational Field) over Rational Field
        sage: pushout(GP_ZZ(QQ), GP_QQ(ZZ))
        Generalized Polynomial Ring in X^(Rational Field) over Rational Field
        sage: pushout(GP_ZZt(QQ), GP_QQ(ZZ))
        Generalized Polynomial Ring in X^(Univariate Polynomial Ring in t
          over Rational Field) over Rational Field
        sage: pushout(GP_ZZt(ZZ), GP_QQ(QQ))
        Generalized Polynomial Ring in X^(Univariate Polynomial Ring in t
          over Rational Field) over Rational Field
        sage: pushout(GP_ZZt(ZZ['a,b']), GP_QQ(ZZ['c,d']))
        Traceback (most recent call last):
        ...
        CoercionException: ('Ambiguous Base Extension', ...)
        sage: pushout(GP_ZZt(ZZ['a,b']), GP_QQ(ZZ['b,c']))
        Generalized Polynomial Ring in X^(Univariate Polynomial Ring in t over Rational Field)
          over Multivariate Polynomial Ring in a, b, c over Integer Ring

    Some tests with Cartesian products::

        sage: from sage.sets.cartesian_product import CartesianProduct
        sage: A = CartesianProduct((ZZ['x'], QQ['y'], QQ['z']), Sets().CartesianProducts())
        sage: B = CartesianProduct((ZZ['x'], ZZ['y'], ZZ['t']['z']), Sets().CartesianProducts())
        sage: A.construction()
        (The cartesian_product functorial construction,
         (Univariate Polynomial Ring in x over Integer Ring,
          Univariate Polynomial Ring in y over Rational Field,
          Univariate Polynomial Ring in z over Rational Field))
        sage: pushout(A, B)
        The Cartesian product of
         (Univariate Polynomial Ring in x over Integer Ring,
          Univariate Polynomial Ring in y over Rational Field,
          Univariate Polynomial Ring in z over Univariate Polynomial Ring in t over Rational Field)
        sage: pushout(ZZ, cartesian_product([ZZ, QQ]))
        Traceback (most recent call last):
        ...
        CoercionException: 'NoneType' object is not iterable

    ::

        sage: from sage.categories.pushout import PolynomialFunctor
        sage: from sage.sets.cartesian_product import CartesianProduct
        sage: class CartesianProductPoly(CartesianProduct):
        ....:     def __init__(self, polynomial_rings):
        ....:         sort = sorted(polynomial_rings, key=lambda P: P.variable_name())
        ....:         super(CartesianProductPoly, self).__init__(sort, Sets().CartesianProducts())
        ....:     def vars(self):
        ....:         return tuple(P.variable_name() for P in self.cartesian_factors())
        ....:     def _pushout_(self, other):
        ....:         if isinstance(other, CartesianProductPoly):
        ....:             s_vars = self.vars()
        ....:             o_vars = other.vars()
        ....:             if s_vars == o_vars:
        ....:                 return
        ....:             return pushout(CartesianProductPoly(
        ....:                     self.cartesian_factors() +
        ....:                     tuple(f for f in other.cartesian_factors()
        ....:                           if f.variable_name() not in s_vars)),
        ....:                 CartesianProductPoly(
        ....:                     other.cartesian_factors() +
        ....:                     tuple(f for f in self.cartesian_factors()
        ....:                           if f.variable_name() not in o_vars)))
        ....:         C = other.construction()
        ....:         if C is None:
        ....:             return
        ....:         elif isinstance(C[0], PolynomialFunctor):
        ....:             return pushout(self, CartesianProductPoly((other,)))

    ::

        sage: pushout(CartesianProductPoly((ZZ['x'],)),
        ....:         CartesianProductPoly((ZZ['y'],)))
        The Cartesian product of
         (Univariate Polynomial Ring in x over Integer Ring,
          Univariate Polynomial Ring in y over Integer Ring)
        sage: pushout(CartesianProductPoly((ZZ['x'], ZZ['y'])),
        ....:         CartesianProductPoly((ZZ['x'], ZZ['z'])))
        The Cartesian product of
         (Univariate Polynomial Ring in x over Integer Ring,
          Univariate Polynomial Ring in y over Integer Ring,
          Univariate Polynomial Ring in z over Integer Ring)
        sage: pushout(CartesianProductPoly((QQ['a,b']['x'], QQ['y'])),
        ....:         CartesianProductPoly((ZZ['b,c']['x'], SR['z'])))
        The Cartesian product of
         (Univariate Polynomial Ring in x over
            Multivariate Polynomial Ring in a, b, c over Rational Field,
          Univariate Polynomial Ring in y over Rational Field,
          Univariate Polynomial Ring in z over Symbolic Ring)

    ::

        sage: pushout(CartesianProductPoly((ZZ['x'],)), ZZ['y'])
        The Cartesian product of
         (Univariate Polynomial Ring in x over Integer Ring,
          Univariate Polynomial Ring in y over Integer Ring)
        sage: pushout(QQ['b,c']['y'], CartesianProductPoly((ZZ['a,b']['x'],)))
        The Cartesian product of
         (Univariate Polynomial Ring in x over
            Multivariate Polynomial Ring in a, b over Integer Ring,
          Univariate Polynomial Ring in y over
            Multivariate Polynomial Ring in b, c over Rational Field)

    ::

        sage: pushout(CartesianProductPoly((ZZ['x'],)), ZZ)
        Traceback (most recent call last):
        ...
        CoercionException: No common base ("join") found for
        The cartesian_product functorial construction(...) and None(Integer Ring):
        (Multivariate) functors are incompatible.

    AUTHORS:

    - Robert Bradshaw
    - Peter Bruin
    - Simon King
    - Daniel Krenn
    - David Roe
    """
    if R is S or R == S:
        return R

    if hasattr(R, '_pushout_'):
        P = R._pushout_(S)
        if P is not None:
            return P

    if hasattr(S, '_pushout_'):
        P = S._pushout_(R)
        if P is not None:
            return P

    if isinstance(R, type):
        R = type_to_parent(R)

    if isinstance(S, type):
        S = type_to_parent(S)

    R_tower = construction_tower(R)
    S_tower = construction_tower(S)
    Rs = [c[1] for c in R_tower]
    Ss = [c[1] for c in S_tower]

    # If there is a multivariate construction functor in the tower, we must chop off the end
    # because tuples don't have has_coerce_map_from functions and to align with the
    # modification of Rs and Ss below
    from sage.structure.parent import Parent
    if not isinstance(Rs[-1], Parent):
        Rs = Rs[:-1]
    if not isinstance(Ss[-1], Parent):
        Ss = Ss[:-1]

    if R in Ss:
        if not any(c[0].coercion_reversed for c in S_tower[1:]):
            return S
    elif S in Rs:
        if not any(c[0].coercion_reversed for c in R_tower[1:]):
            return R

    if Rs[-1] in Ss:
        Rs, Ss = Ss, Rs
        R_tower, S_tower = S_tower, R_tower

    # look for join
    Z = None
    if Ss[-1] in Rs:
        if Rs[-1] == Ss[-1]:
            while Rs and Ss and Rs[-1] == Ss[-1]:
                Rs.pop()
                Z = Ss.pop()
        else:
            Rs = Rs[:Rs.index(Ss[-1])]
            Z = Ss.pop()

    # look for topmost coercion
    elif S.has_coerce_map_from(Rs[-1]):
        while not Ss[-1].has_coerce_map_from(Rs[-1]):
            Ss.pop()
        while len(Rs) > 0 and Ss[-1].has_coerce_map_from(Rs[-1]):
            Rs.pop()
        Z = Ss.pop()

    elif R.has_coerce_map_from(Ss[-1]):
        while not Rs[-1].has_coerce_map_from(Ss[-1]):
            Rs.pop()
        while len(Ss) > 0 and Rs[-1].has_coerce_map_from(Ss[-1]):
            Ss.pop()
        Z = Rs.pop()

    if Z is None and R_tower[-1][0] is not None:
        Z = R_tower[-1][0].common_base(S_tower[-1][0], R_tower[-1][1], S_tower[-1][1])
        R_tower = expand_tower(R_tower[:len(Rs)])
        S_tower = expand_tower(S_tower[:len(Ss)])
    else:
        # Rc is a list of functors from Z to R and Sc is a list of functors from Z to S
        R_tower = expand_tower(R_tower[:len(Rs)+1])
        S_tower = expand_tower(S_tower[:len(Ss)+1])
    Rc = [c[0] for c in R_tower[1:]]
    Sc = [c[0] for c in S_tower[1:]]

    all = IdentityConstructionFunctor()

    def apply_from(Xc):
        c = Xc.pop()
        if c.coercion_reversed:
            Yc = Sc if Xc is Rc else Rc
            Y_tower = S_tower if Xc is Rc else R_tower
            Y_partial = Y_tower[len(Yc)][1]
            if not (c * all)(Z).has_coerce_map_from(Y_partial):
                return all
        return c * all

    try:

        while len(Rc) > 0 or len(Sc) > 0:
            # print Z
            # if we are out of functors in either tower, there is no ambiguity
            if len(Sc) == 0:
                all = apply_from(Rc)
            elif len(Rc) == 0:
                all = apply_from(Sc)
            # if one of the functors has lower rank, do it first
            elif Rc[-1].rank < Sc[-1].rank:
                all = apply_from(Rc)
            elif Sc[-1].rank < Rc[-1].rank:
                all = apply_from(Sc)
            else:
                # the ranks are the same, so things are a bit subtler
                if Rc[-1] == Sc[-1]:
                    # If they are indeed the same operation, we only do it once.
                    # The \code{merge} function here takes into account non-mathematical
                    # distinctions (e.g. single vs. multivariate polynomials).
                    cR = Rc.pop()
                    cS = Sc.pop()
                    c = cR.merge(cS) or cS.merge(cR)
                    if c:
                        all = c * all
                    else:
                        raise CoercionException("Incompatible Base Extension %r, %r (on %r, %r)" % (R, S, cR, cS))
                else:
                    # Now we look ahead to see if either top functor is
                    # applied later on in the other tower.
                    # If this is the case for exactly one of them, we unambiguously
                    # postpone that operation, but if both then we abort.
                    if Rc[-1] in Sc:
                        if Sc[-1] in Rc:
                            raise CoercionException("Ambiguous Base Extension", R, S)
                        else:
                            all = apply_from(Sc)
                    elif Sc[-1] in Rc:
                        all = apply_from(Rc)
                    # If, perchance, the two functors commute, then we may do them in any order.
                    elif Rc[-1].commutes(Sc[-1]) or Sc[-1].commutes(Rc[-1]):
                        all = Sc.pop() * Rc.pop() * all
                    else:
                        # try and merge (default merge is failure for unequal functors)
                        cR = Rc.pop()
                        cS = Sc.pop()
                        c = cR.merge(cS) or cS.merge(cR)
                        if c is not None:
                            all = c * all
                        else:
                            # Otherwise, we cannot proceed.
                            raise CoercionException("Ambiguous Base Extension", R, S)

        return all(Z)

    except CoercionException:
        raise
    except (TypeError, ValueError, AttributeError, NotImplementedError) as ex:
        # We do this because we may be trying all kinds of things that don't
        # make sense, and in this case simply want to return that a pushout
        # couldn't be found.
        raise CoercionException(ex)



def pushout_lattice(R, S):
    r"""
    Given a pair of objects `R` and `S`, try to construct a
    reasonable object `Y` and return maps such that
    canonically `R \leftarrow Y \rightarrow S`.

    ALGORITHM:

    This is based on the model that arose from much discussion at
    Sage Days 4.  Going up the tower of constructions of `R` and `S`
    (e.g. the reals come from the rationals come from the integers),
    try to find a common parent, and then try to fill in a lattice
    with these two towers as sides with the top as the common ancestor
    and the bottom will be the desired ring.

    See the code for a specific worked-out example.

    EXAMPLES::

        sage: from sage.categories.pushout import pushout_lattice
        sage: A, B = pushout_lattice(Qp(7), Frac(ZZ['x']))
        sage: A.codomain()
        Fraction Field of Univariate Polynomial Ring in x over 7-adic Field with capped relative precision 20
        sage: A.codomain() is B.codomain()
        True
        sage: A, B = pushout_lattice(ZZ, MatrixSpace(ZZ[['x']], 3, 3))
        sage: B
        Identity endomorphism of Full MatrixSpace of 3 by 3 dense matrices over Power Series Ring in x over Integer Ring

    AUTHOR:

    - Robert Bradshaw

    """
    R_tower = construction_tower(R)
    S_tower = construction_tower(S)
    Rs = [c[1] for c in R_tower]
    Ss = [c[1] for c in S_tower]

    # look for common ancestor
    start = None
    for Z in Rs:
        if Z in Ss:
            start = Z
    if start is None:
        # Should I test for a map between the tops of the towers?
        # Or, if they're both not ZZ, is it hopeless?
        return None

    # truncate at common ancestor
    R_tower = list(reversed(R_tower[:Rs.index(start)+1]))
    S_tower = list(reversed(S_tower[:Ss.index(start)+1]))
    Rs = [c[1] for c in R_tower] # the list of objects
    Ss = [c[1] for c in S_tower]
    Rc = [c[0] for c in R_tower] # the list of functors
    Sc = [c[0] for c in S_tower]

    # Here we try and construct a 2-dimensional lattice as follows.
    # Suppose our towers are Z -> Q -> Qp = R and Z -> Z[t] -> Frac(Z[t]) = S
    lattice = {}
    # First we fill in the sides
    #
    #         Z
    #       /   \
    #      Q    Z[t]
    #    /         \
    #   Qp       Frac(Z[t])
    #
    for i in range(len(Rs)):
        lattice[i,0] = Rs[i]
    for j in range(len(Ss)):
        lattice[0,j] = Ss[j]

    # Now we attempt to fill in the center, one (diagonal) row at a time,
    # one commuting square at a time.
    #
    #          Z
    #       /    \
    #      Q     Z[t]
    #    /   \  /    \
    #   Qp   Q[t]   Frac(Z[t])
    #    \   /
    #    Qp[t]
    #
    # There is always exactly one "correct" path/order in which to apply operations
    # from the top to the bottom. In our example, this is down the far left side.
    # We keep track of which that is by clearing out Rc and Sc as we go along.
    #
    # Note that when applying the functors in the correct order, base extension
    # is not needed (though it may occur in the resulting morphisms).
    #
    for i in range(len(Rc)-1):
        for j in range(len(Sc)-1):
            try:
                if lattice[i,j+1] == lattice[i+1,j]:
                    # In this case we have R <- S -> R
                    # We don't want to perform the operation twice
                    # and all subsequent squares will come from objects
                    # where the operation was already performed (either
                    # to the left or right)
                    Rc[i] = Sc[j] = None # IdentityConstructionFunctor()
                    lattice[i+1,j+1] = lattice[i,j+1]
                elif Rc[i] is None and Sc[j] is None:
                    lattice[i+1,j+1] = lattice[i,j+1]
                elif Rc[i] is None:
                    lattice[i+1,j+1] = Sc[j](lattice[i+1,j])
                elif Sc[j] is None:
                    lattice[i+1,j+1] = Rc[i](lattice[i,j+1])
                else:
                    # For now, we just look at the rank.
                    # TODO: be more sophisticated and query the functors themselves
                    if Rc[i].rank < Sc[j].rank:
                        lattice[i+1,j+1] = Sc[j](lattice[i+1,j])
                        Rc[i] = None # force us to use pre-applied Rc[i]
                    else:
                        lattice[i+1,j+1] = Rc[i](lattice[i,j+1])
                        Sc[j] = None # force us to use pre-applied Sc[i]
            except (AttributeError, NameError):
                # print i, j
                # pp(lattice)
                for i in range(100):
                    for j in range(100):
                        try:
                            R = lattice[i,j]
                            print i, j, R
                        except KeyError:
                            break
                raise CoercionException("%s does not support %s" % (lattice[i,j], 'F'))

    # If we are successful, we should have something that looks like this.
    #
    #          Z
    #       /    \
    #      Q     Z[t]
    #    /   \  /    \
    #   Qp   Q[t]   Frac(Z[t])
    #    \   /  \    /
    #    Qp[t]  Frac(Q[t])
    #      \      /
    #     Frac(Qp[t])
    #
    R_loc = len(Rs)-1
    S_loc = len(Ss)-1

    # Find the composition coercion morphisms along the bottom left...
    if S_loc > 0:
        R_map = lattice[R_loc,1].coerce_map_from(R)
        for i in range(1, S_loc):
            map = lattice[R_loc, i+1].coerce_map_from(lattice[R_loc, i]) # The functor used is implicit here, should it be?
            R_map = map * R_map
    else:
        R_map = R.coerce_map_from(R) # id

    # ... and bottom right
    if R_loc > 0:
        S_map = lattice[1, S_loc].coerce_map_from(S)
        for i in range(1, R_loc):
            map = lattice[i+1, S_loc].coerce_map_from(lattice[i, S_loc])
            S_map = map * S_map
    else:
        S_map = S.coerce_map_from(S) # id

    return R_map, S_map


## def pp(lattice):
##     """
##     Used in debugging to print the current lattice.
##     """
##     for i in range(100):
##         for j in range(100):
##             try:
##                 R = lattice[i,j]
##                 print i, j, R
##             except KeyError:
##                 break

def construction_tower(R):
    """
    An auxiliary function that is used in :func:`pushout` and :func:`pushout_lattice`.

    INPUT:

    An object

    OUTPUT:

    A constructive description of the object from scratch, by a list of pairs
    of a construction functor and an object to which the construction functor
    is to be applied. The first pair is formed by ``None`` and the given object.

    EXAMPLE::

        sage: from sage.categories.pushout import construction_tower
        sage: construction_tower(MatrixSpace(FractionField(QQ['t']),2))
        [(None, Full MatrixSpace of 2 by 2 dense matrices over Fraction Field of Univariate Polynomial Ring in t over Rational Field), (MatrixFunctor, Fraction Field of Univariate Polynomial Ring in t over Rational Field), (FractionField, Univariate Polynomial Ring in t over Rational Field), (Poly[t], Rational Field), (FractionField, Integer Ring)]

    """
    tower = [(None, R)]
    c = R.construction()
    from sage.structure.parent import Parent
    while c is not None:
        f, R = c
        if not isinstance(f, ConstructionFunctor):
            f = BlackBoxConstructionFunctor(f)
        tower.append((f,R))
        if not isinstance(R, Parent):
            break
        c = R.construction()
    return tower

def expand_tower(tower):
    """
    An auxiliary function that is used in :func:`pushout`.

    INPUT:

    A construction tower as returned by :func:`construction_tower`.

    OUTPUT:

    A new construction tower with all the construction functors expanded.

    EXAMPLE::

        sage: from sage.categories.pushout import construction_tower, expand_tower
        sage: construction_tower(QQ['x,y,z'])
        [(None, Multivariate Polynomial Ring in x, y, z over Rational Field),
         (MPoly[x,y,z], Rational Field),
         (FractionField, Integer Ring)]
        sage: expand_tower(construction_tower(QQ['x,y,z']))
        [(None, Multivariate Polynomial Ring in x, y, z over Rational Field),
         (MPoly[z], Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field),
         (MPoly[y], Univariate Polynomial Ring in x over Rational Field),
         (MPoly[x], Rational Field),
         (FractionField, Integer Ring)]
    """
    new_tower = []
    for f, R in reversed(tower):
        if f is None:
            new_tower.append((f, R))
        else:
            fs = f.expand()
            for ff in reversed(fs[1:]):
                new_tower.append((ff, R))
                R = ff(R)
            new_tower.append((fs[0], R))
    return list(reversed(new_tower))

def type_to_parent(P):
    """
    An auxiliary function that is used in :func:`pushout`.

    INPUT:

    A type

    OUTPUT:

    A Sage parent structure corresponding to the given type

    TEST::

        sage: from sage.categories.pushout import type_to_parent
        sage: type_to_parent(int)
        Integer Ring
        sage: type_to_parent(float)
        Real Double Field
        sage: type_to_parent(complex)
        Complex Double Field
        sage: type_to_parent(list)
        Traceback (most recent call last):
        ...
        TypeError: Not a scalar type.
    """
    import sage.rings.all
    if P in [int, long]:
        return sage.rings.all.ZZ
    elif P is float:
        return sage.rings.all.RDF
    elif P is complex:
        return sage.rings.all.CDF
    else:
        raise TypeError("Not a scalar type.")
