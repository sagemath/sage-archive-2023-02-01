r"""
Asymptotic Growth Group
"""

import re

from sage.structure.element import MultiplicativeGroupElement
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.real_mpfr import RR



class GenericGrowthElement(MultiplicativeGroupElement):
    r"""
    Class for a generic asymptotic growth group element. These
    elements hold exactly one asymptotic term and can be compared to
    each other, multiplied and divided, but possess no explicit
    coefficient. At this stage, just the order of magnitude shall be
    managed. In this class, only base structure is handled.
    For a concrete realization see :class:`GrowthElementPower`.

    INPUT:

    - ``parent`` -- a GenericGrowthGroup, the parent of the
      element.

    OUTPUT:

    A generic element of a generic growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GenericGrowthGroup()
        sage: e = agg.GenericGrowthElement(P); e
        Generic element of a Generic Asymptotic Growth Group
        sage: e.parent()
        Generic Asymptotic Growth Group
    """

    def __init__(self, parent):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: e = agg.GenericGrowthElement(P)
            sage: e.category()
            Category of elements of Generic Asymptotic Growth Group

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e = agg.GenericGrowthElement(None)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError("The parent must be provided")
        super(GenericGrowthElement, self).__init__(parent=parent)


    def _mul_(self, other):
        r"""
        Abstract multiplication method for generic asymptotic growth
        group elements. For a concrete realization see
        :meth:`GrowthElementPower._mul_`.

        INPUT:

        - ``other`` -- a GenericGrowthElement from the same parent as
          ``self``.

        OUTPUT:

        A GenericGrowthElement representing the product of ``self``
        and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e1._mul_(e2)
            x^5
            sage: e1 * e2 * e1
            x^7
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def _repr_(self):
        r"""
        Represent the abstract element of an abstract asymptotic growth
        group as "Generic element of a Generic Asymptotic Growth
        Group".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: e = agg.GenericGrowthElement(P); e._repr_()
            'Generic element of a Generic Asymptotic Growth Group'
        """
        return "Generic element of a Generic Asymptotic Growth Group"


    def is_le_one(self):
        r"""
        Abstract method for comparison with one. See
        :meth:`GrowthElementPower.is_le_one` for a
        concrete implementation.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: (~P.gen()).is_le_one()
            True
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def __le__(self, other):
        r"""
        Abstract method for comparison of ``self`` and ``other``. See
        :meth:`GrowthElementPower.__le__` and
        :meth:`GrowthElementPower._le_` for a concrete realization.

        INPUT:

        - ``other`` -- a generic asymptotic growth element.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: (~P.gen()).__le__(P.gen())
            True
        """
        raise NotImplementedError("Only implemented in concrete realizations")


    def __hash__(self):
        r"""
        Return the hash of the representation of the element produced
        by :meth:`GenericGrowthElement._repr_`

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup();
            sage: e = agg.GenericGrowthElement(P); e.__hash__()  # random
            -7923874249531374658
        """
        return hash(repr(self))


class GenericGrowthGroup(Parent, UniqueRepresentation):
    r"""
    Class for the generic asymptotic growth group. Only base
    structure is handled here, for a concrete realization see
    :class:`GrowthGroupPower`.

    INPUT:

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a
      subcategory of ``Join of Category of groups and Category
      of posets``. This is also the default category if ``None``
      is specified.

    - ``base`` -- The base of the parent associated to concrete
      implementations of this abstract base class.

    OUTPUT:

    A generic asymptotic growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GenericGrowthGroup(); P
        Generic Asymptotic Growth Group
    """

    # TODO: implement some sort of "assume", where basic assumptions
    # for the variables can be stored. --> within the cartesian product

    # TODO: implement a cartesian product class for the asymptotic
    # growth group. the "standard" cart. prod. produces an element
    # of the cart. prod. class, which does not have a order structure.

    # enable the category framework for elements
    Element = GenericGrowthElement


    def __init__(self, category=None, base=None):
        r"""
        See :class:`GenericGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup().category()
            Join of Category of groups and Category of posets

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup()
            sage: P.is_parent_of(agg.GenericGrowthElement(P))
            True
            sage: P2 = agg.GenericGrowthGroup(category=FiniteGroups()\
            ....:                               & Posets())
            sage: P2.category()
            Join of Category of finite groups and Category of finite posets
            sage: P3 = agg.GenericGrowthGroup(category=Rings())
            Traceback (most recent call last):
            ...
            ValueError: (Category of rings,) is not a subcategory of Join of Category of groups and Category of posets
        """
        from sage.categories.groups import Groups
        from sage.categories.posets import Posets

        if category is None:
            category = Groups() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category, )
            if not any(cat.is_subcategory(Groups() & Posets()) for cat in
                       category):
                raise ValueError("%s is not a subcategory of %s"
                                 % (category, Groups() & Posets()))
        super(GenericGrowthGroup, self).__init__(category=category,
                                                 base=base)


    def le(self, x, y):
        r"""
        Return whether the asymptotic order of magnitude of `x` is less
        than or equal to the asymptotic order of magnitude of `y`.

        INPUT:

        - ``x``, ``y`` -- elements of ``self``.


        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P.<x> = agg.GrowthGroupPower()
            sage: P.le(x, x^2)
            True
            sage: P.le(x^2, x)
            False
            sage: P.le(x^0,1)
            True
        """
        return (self(x) / self(y)).is_le_one()


    def _repr_(self):
        r"""
        Represent the generic asymptotic growth group as
        "Generic Asymptotic Growth Group".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GenericGrowthGroup()._repr_()
            'Generic Asymptotic Growth Group'
        """
        return "Generic Asymptotic Growth Group"


    def __hash__(self):
        r"""
        Return the hash of the representation of the group produced by
        :meth:`GenericGrowthGroup._repr_`

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GenericGrowthGroup(); P.__hash__()  # random
            8493512696244708699
        """
        return hash(repr(self))


class GrowthElementPower(GenericGrowthElement):
    r"""
    Class for the concrete realization of asymptotic growth group
    elements in the case of polynomial growth. These elements hold
    exactly one asymptotic term.

    A power growth element represents a polynomial term
    `\operatorname{variable}^{\operatorname{exponent}}`.
    More complex constructions including logarithmic or exponential
    terms can be constructed via a cartesian product of the related
    growth groups. Asymptotic growth elements can be multiplied,
    divided, inverted, and compared to each other. However, they
    possess no explicit coefficient.

    The elements can be specified by either an expression ``x`` being
    a string, an element from the symbolic or a polynomial ring or the
    integer `1`. On the other hand, elements can also be specified
    by their exponent.

    INPUT:

    - ``parent`` -- a :class:`GrowthGroupPower`, the
      parent of the element.
    - ``x`` -- an expression (string, polynomial ring element,
      symbolic ring element, or the integer `1`) representing
      the element to be initialized.
    - ``exponent`` -- the exponent of the power element.

    OUTPUT:

    An asymptotic power growth element with the specified
    parent and magnitude of growth, i.e. exponent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GrowthGroupPower("x")
        sage: e1 = P(x=1); e1
        1
        sage: e2 = P(x=None, exponent=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
        sage: P.le(e1, P.gen()) and P.le(P.gen(), e2)
        True
    """
    # TODO: implement comparison for the cartesian product structure.


    def __init__(self, parent, x=None, exponent=None):
        r"""
        See :class:`GrowthElementPower` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P.gen(); e1
            x
            sage: e1.is_idempotent()
            False
            sage: e1.is_one()
            False
            sage: e1.parent()
            Asymptotic Power Growth Group in x over Integer Ring
            sage: e2 = P.one(); e2
            1
            sage: e2.is_idempotent() and e2.is_one()
            True
        """
        if exponent not in RR:
            raise NotImplementedError("Non-real exponents are not supported.")
        else:
            self.exponent = parent.base()(exponent)
        super(GrowthElementPower, self).__init__(parent=parent)


    def _repr_(self):
        r"""
        Represent the asymptotic power growth element as
        ``variable^{exponent}``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x", base=QQ)
            sage: P(x=1)._repr_()
            '1'
            sage: P(x=None, exponent=5)._repr_()
            'x^5'
            sage: P(x=None, exponent=1/2)._repr_()
            'x^(1/2)'
        """
        if self.exponent == 0:
            return "1"
        elif self.exponent == 1:
            return self.parent().variable
        elif self.exponent in ZZ:
            return self.parent().variable + "^" + str(self.exponent)
        else:
            return self.parent().variable + "^(" + str(self.exponent) + ")"


    def _mul_(self, other):
        r"""
        Multiply two asymptotic power growth elements from the
        same parent by adding their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element to be
          multiplied with ``self``.

        OUTPUT:

        An asymptotic power growth element representing the product
        of ``self`` and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e3 = e1._mul_(e2); e3
            x^5
            sage: e3 == e1 * e2
            True
        """
        cls = self.__class__
        return cls(self.parent(), exponent=self.exponent + other.exponent)


    def __invert__(self):
        r"""
        Return the multiplicative inverse from a given asymptotic power
        growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse asymptotic power growth element
        of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = e1.__invert__(); e2
            x^-2
            sage: e2 == ~e1
            True
        """
        cls = self.__class__
        return cls(self.parent(), exponent=-self.exponent)


    def _div_(self, other):
        r"""
        Divide two asymptotic power growth elements from the same
        parent by subtracting their exponents.

        INPUT:

        - ``other`` -- the asymptotic growth element which ``self``
          is divided by.

        OUTPUT:

        The result of the division of ``self`` by ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e1._div_(P.gen())
            x
            sage: e1._div_(P.gen()) == e1 / P.gen()
            True
        """
        cls = self.__class__
        return cls(self.parent(), exponent=self.exponent - other.exponent)


    def __pow__(self, power):
        r"""
        Return a asymptotic power element to the power of
        ``power``.

        INPUT:

        - ``power`` -- a rational number.

        OUTPUT:

        The asymptotic power growth element ``self`` to the power of
        ``power``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: P.gen().__pow__(5)
            x^5
            sage: P.gen().__pow__(1/2)
            x^(1/2)
            sage: P.gen()^7
            x^7
        """
        new_exponent = self.exponent * power
        if new_exponent in self.parent().base():
            return self.parent()(None, exponent=self.exponent * power)

        if new_exponent in RR:
            pnt = GrowthGroupPower(self.parent().variable,
                                   base=new_exponent.parent())
            return pnt(None, exponent=new_exponent)
        else:
            raise NotImplementedError("Only real exponents are implemented.")


    def __eq__(self, other):
        r"""
        Tests for equality of the given elements (with taking care of
        different parents by using the coercion model).

        INPUT:

        - ``other`` -- power element to be compared with ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: P2 = agg.GrowthGroupPower("x", base=QQ)
            sage: P1.gen() == P2.gen()
            True
        """
        from sage.structure.sage_object import have_same_parent
        if have_same_parent(self, other):
            return self._eq_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.eq)
        except TypeError:
            return False


    def _eq_(self, other):
        r"""
        Return whether the asymptotic power growth elements
        ``self`` and ``other`` are equal and have the same parent.

        INPUT:

        - ``other`` -- an asymptotic power growth element to
          be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P(x=None, exponent=1)
            sage: e1._eq_(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self.exponent == other.exponent


    def __le__(self, other):
        r"""
        Return whether the growth of the asymptotic power growth
        element ``self`` is less than or equal to the growth of the
        asymptotic power growth element ``other``.

        INPUT:

        - ``other`` -- a growth power element to be compared
          to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: P2 = agg.GrowthGroupPower("x", base=QQ)
            sage: P1.gen() <= P2.gen()^2
            True
        """
        from sage.structure.sage_object import have_same_parent
        if have_same_parent(self, other):
            return self._le_(other)

        from sage.structure.element import get_coercion_model
        import operator
        try:
            return get_coercion_model().bin_op(self, other, operator.le)
        except TypeError:
            return False


    def _le_(self, other):
        r"""
        Return whether the exponent of ``self`` is less than or equal
        to the exponent of ``other`` by calling
        :meth:`GrowthGroupPower.le`.

        INPUT:

        - ``other`` -- a growth power element from the same parent to
          be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P.gen()
            sage: e2 = P(None, exponent=2)
            sage: e1._le_(e2)
            True
            sage: e2._le_(e1)
            False
        """
        return self.parent().le(self, other)


    def is_le_one(self):
        r"""
        Return whether or not the growth of the asymptotic power
        growth element ``self`` is less than or equal to the
        (constant) growth of `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P = agg.GrowthGroupPower("x")
            sage: e1 = P.gen()
            sage: e1.is_le_one()
            False
            sage: (P.one() / P.gen()).is_le_one()
            True
        """
        return self.exponent <= 0

    def __hash__(self):
        r"""
        Return the hash of the tuple containing the exponent and the
        parent.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P.<x> = agg.GrowthGroupPower()
            sage: hash(x)  # random
            -3234005094684624010
        """
        return hash((self.exponent, self.parent()))


class GrowthGroupPower(GenericGrowthGroup):
    r"""
    Class for the concrete realization of asymptotic power growth
    groups. These are the parents for the
    :class:`GrowthElementPower` elements. These parents
    contain single monomial terms in a specified variable with
    exponents from a specified base.

    More complex growth elements can be constructed by constructing the
    cartesian product of various asymptotic growth groups.

    INPUT:

    - ``variable`` -- either a symbol from the symbolic ring, a
      generator from a polynomial ring, or a alphanumeric string
      starting with a letter and optionally containing underscores.

    - ``base`` -- the base structure containing the exponents.

    - ``category`` -- The category of the parent can be specified
      in order to broaden the base structure. Has to be a
      subcategory of ``Join of Category of groups and Category of
      posets``. This is also the default category if ``None`` is
      specified.

    OUTPUT:

    An asymptotic power growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as agg
        sage: P = agg.GrowthGroupPower("x"); P
        Asymptotic Power Growth Group in x over Integer Ring
    """
    # TODO: implement the cartesian product structure


    @staticmethod
    def __classcall__(cls, variable=None, base=ZZ, category=None, names=None):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: x1, x2, x3 = "x", PolynomialRing(ZZ, "x").gen(), var("x")
            sage: P1 = agg.GrowthGroupPower(x1)
            sage: P2 = agg.GrowthGroupPower(x2)
            sage: P3 = agg.GrowthGroupPower(x3)
            sage: P1 is P2 and P2 is P3
            True
            sage: x4 = buffer("xylophone", 0, 1)
            sage: P4 = agg.GrowthGroupPower(x4)
            sage: P1 is P4
            True
        """
        if names is not None and len(names) == 1:
            variable = names[0]

        if hasattr(variable, "is_symbol") and variable.is_symbol():
            variable = repr(variable)
        elif hasattr(variable, "is_gen") and variable.is_gen():
            variable = repr(variable)
        elif isinstance(variable, buffer):
            variable = str(variable)
        return super(GrowthGroupPower, cls).\
            __classcall__(cls, variable, base, category)


    def __init__(self, variable, base=ZZ, category=None):
        r"""
        For more information see :class:`GrowthGroupPower`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x"); P1
            Asymptotic Power Growth Group in x over Integer Ring
            sage: var('n')
            n
            sage: P2 = agg.GrowthGroupPower(n, base=QQ); P2
            Asymptotic Power Growth Group in n over Rational Field
            sage: y = PolynomialRing(ZZ, "y").gen()
            sage: P3 = agg.GrowthGroupPower(y); P3
            Asymptotic Power Growth Group in y over Integer Ring
            sage: P4 = agg.GrowthGroupPower("y"); P4
            Asymptotic Power Growth Group in y over Integer Ring
            sage: P3 is P4
            True
        """
        if variable is None:
            raise ValueError("Variable for initialization required.")
        else:
            if re.match("^[A-Za-z][A-Za-z0-9_]*$", variable):
                self.variable = str(variable)
            else:
                raise ValueError("Only alphanumeric strings starting with a "
                                 "letter may be variables")
        super(GrowthGroupPower, self).\
            __init__(category=category, base=base)
        self._populate_coercion_lists_()

    # enable the category framework for elements
    Element = GrowthElementPower


    def _coerce_map_from_(self, S):
        r"""
        Another GrowthGroupPower ``S`` coerces into
        ``self`` if and only if the base of ``S`` coerces into the
        base of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: P2 = agg.GrowthGroupPower("x", base=QQ)
            sage: bool(P1._coerce_map_from_(P2))
            False
            sage: bool(P2._coerce_map_from_(P1))
            True
        """
        if isinstance(S, GrowthGroupPower):
            if self.base().coerce_map_from(S.base()) is not None:
                return True


    def _element_constructor_(self, x, exponent=None):
        r"""
        Coerce object to this asymptotic power growth group.

        INPUT:

        - ``x`` -- an expression (string, polynomial ring element,
          symbolic ring element, or some integer) representing the
          element to be initialized.
        - ``exponent`` -- the exponent of the asymptotic element.

        OUTPUT:

        An element of an asymptotic power growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P1 = agg.GrowthGroupPower("x", base=ZZ)
            sage: P2 = agg.GrowthGroupPower("x", base=QQ)
            sage: e = P2(None, exponent=3/2) / P1.gen(); e
            x^(1/2)
            sage: e.parent() is P2
            True
        """
        if isinstance(x, GrowthElementPower):
            return self.element_class(self, None, exponent=x.exponent)

        if x is None and exponent is None:
            raise ValueError("Neither x nor exponent are specified.")
        elif x is not None and exponent is not None:
            raise ValueError("Both x and exponent are specified.")
        elif exponent is None:
            if x == 1:
                exponent = 0
            else:
                raise NotImplementedError("Parsing of %s is not yet "
                                          "implemented" % x)
        return self.element_class(self, None, exponent=exponent)


    def _repr_(self):
        r"""
        Represent the asymptotic power growth group as
        "Asymptotic Power Growth Group in ``variable``
        over ``base``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: agg.GrowthGroupPower("x")._repr_()
            'Asymptotic Power Growth Group in x over Integer Ring'
            sage: agg.GrowthGroupPower("v_107", base=QQ)._repr_()
            'Asymptotic Power Growth Group in v_107 over Rational Field'
        """
        return "Asymptotic Power Growth Group in %s over %s" \
               % (self.variable, self.base())


    def gens(self):
        r"""
        Return a list with all generators of ``self``. For power growth
        groups, this is a list with only one element: the variable
        to the power `1`.

        INPUT:

        Nothing.

        OUTPUT:

        A list of generators.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P.<x> = agg.GrowthGroupPower()
            sage: P.gens()
            (x,)
        """
        return (self.gen(), )


    def ngens(self):
        r"""
        Return the number of generators of ``self``. For power growth
        groups, this is exactly `1`.

        INPUT:

        Nothing.

        OUTPUT:

        The number of generators of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P.<x> = agg.GrowthGroupPower()
            sage: P.ngens()
            1
        """
        return 1


    def one(self):
        r"""
        Return the neutral element of the asymptotic power growth
        group.

        INPUT:

        Nothing.

        OUTPUT:

        The neutral element of the asymptotic power growth group
        ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.GrowthGroupPower("x").one(); e1
            1
            sage: e1.is_idempotent()
            True
        """
        return self(1)


    def gen(self):
        r"""
        Return the asymptotic power growth element with
        exponent 1.

        INPUT:

        Nothing.

        OUTPUT:

        The asymptotic power growth element with exponent 1.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: e1 = agg.GrowthGroupPower("x").gen(); e1
            x
            sage: e1.exponent == 1
            True
        """
        return self(x=None, exponent=1)

    def __hash__(self):
        r"""
        Return the hash of the tuple containing the variable and the
        base.

        INPUT:

        Nothing.

        OUTPUT:

        An integer.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as agg
            sage: P.<x> = agg.GrowthGroupPower()
            sage: hash(P)  # random
            -8144479309627091876
        """
        return hash((self.variable, self.base()))
