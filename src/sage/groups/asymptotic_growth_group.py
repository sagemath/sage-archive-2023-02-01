import re

from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.groups import Groups
from sage.categories.posets import Posets

from sage.misc.functional import parent
from sage.categories.partially_ordered_monoids import PartiallyOrderedMonoids
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.symbolic.constants import e
from sage.symbolic.ring import SR

from sage.structure.unique_representation import UniqueRepresentation


class AsymptoticGrowthElement(Element):
    r"""
    Class for an abstract asymptotic growth group element. These elements hold exactly one asymptotic term. These
    elements can be compared to each other, multiplied and divided, but possess no explicit coefficient. At this stage,
    just the order of magnitude shall be managed. In this class, only base structure is handled. For a concrete
    realization see :class:`AsymptoticGrowthElementUnivariate`.

    INPUT:

        - ``parent`` -- The parent of the element.

    OUTPUT:

    Returns an abstract element of an abstract growth group.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as ar
        sage: P = ar.AsymptoticGrowthGroupParent()
        sage: e = ar.AsymptoticGrowthElement(P); e
        Generic element of a structure
        sage: e.parent()
        Abstract Asymptotic Growth Group
    """

    def __init__(self, parent):
        r"""
        See :class:`AsymptoticGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParent()
            sage: e = ar.AsymptoticGrowthElement(P)
            sage: e.category()
            Category of elements of Abstract Asymptotic Growth Group

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: e = ar.AsymptoticGrowthElement(None)
            Traceback (most recent call last):
            ...
            ValueError: The parent must be provided
        """
        if parent is None:
            raise ValueError("The parent must be provided")
        super(AsymptoticGrowthElement, self).__init__(parent=parent)

    def _mul_(self, other):
        r"""
        Abstract multiplication method for abstract asymptotic growth group elements. For a concrete realization see
        :meth:`~sage.groups.asymptotic_growth_group.AsymptoticGrowthElementUnivariate._mul_`.

        INPUT:

            - ``self``, ``other`` -- two AsymptoticGrowthElement objects over the same parent.

        OUTPUT:

        An AsymptoticGrowthElement representing the product of ``self`` and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e1._mul_(e2)
            x^5
            sage: e1 * e2 * e1
            x^7
        """
        pass


class AsymptoticGrowthGroupParent(Parent, UniqueRepresentation):
    r"""
    Class for the abstract asymptotic growth group parent. Only base structure is handled here, for a concrete
    realization see :class:`AsymptoticGrowthGroupParentUnivariate`.

    INPUT:

        - ``category`` -- The category of the parent can be specified in order to broaden the base structure. Has to
          be a subcategory of ``Join of Category of groups and Category of posets``.

    OUTPUT:

    Returns an abstract asymptotic growth group parent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as ar
        sage: P = ar.AsymptoticGrowthGroupParent(); P
        Abstract Asymptotic Growth Group
    """

    # TODO: implement some sort of "assume", where basic assumptions for the variables can be stored. (probably in
    #       the multivariate implementation of the asymptotic growth group).

    # TODO: implement a cartesian product class for the asymptotic growth group. the "standard" cart. prod. produces
    #       an element of the cart. prod. class, which does not have a order structure.

    # enable the category framework for elements
    Element = AsymptoticGrowthElement

    def __init__(self, category=None):
        r"""
        See :class:`AsymptoticGrowthElement` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParent()
            sage: P.category()
            Join of Category of groups and Category of posets

        TESTS::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParent()
            sage: P.is_parent_of(ar.AsymptoticGrowthElement(P))
            True
            sage: P2 = ar.AsymptoticGrowthGroupParent(category=FiniteGroups() & Posets())
            sage: P2.category()
            Join of Category of finite groups and Category of finite posets
            sage: P3 = ar.AsymptoticGrowthGroupParent(category=Rings())
            Traceback (most recent call last):
            ...
            ValueError: (Category of rings,) is not a subcategory of Join of Category of groups and Category of posets
        """
        if category is None:
            category = Groups() & Posets()
        else:
            if not isinstance(category, tuple):
                category = (category, )
            if not any(cat.is_subcategory(Groups() & Posets()) for cat in category):
                raise ValueError("%s is not a subcategory of %s" % (category, Groups() & Posets()))
        super(AsymptoticGrowthGroupParent, self).__init__(category=category)

    def le(self, x, y):
        r"""
        Return whether the asymptotic order of magnitude of `x` is less than or equal to the asymptotic order of
        magnitude of `y`.

        INPUT:

            - ``x``, ``y`` -- elements of ``self``.


        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: x = P.gen()
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
        Represents the abstract asymptotic growth group as "Abstract Asymptotic Growth Group".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: ar.AsymptoticGrowthGroupParent()._repr_()
            'Abstract Asymptotic Growth Group'
        """
        return "Abstract Asymptotic Growth Group"


class AsymptoticGrowthElementUnivariate(AsymptoticGrowthElement):
    r"""
    Class for the concrete realization of asymptotic growth group elements in the univariate case. These elements
    hold exactly one asymptotic term. The elements can be compared to each other, multiplied and divided, but possess
    no explicit coefficient.

    Basically, a univariate asymptotic growth element represents a polynomial term ``var^{exponent}``. More complex
    constructions including logarithmic or exponential terms can be constructed via a cartesian product. Asymptotic
    growth elements can be multiplied, divided, inverted, and compared to each other. However, they possess no explicit
    coefficient.

    The elements can be specified by either an expression ``x`` being a string, an element from the symbolic or a
    polynomial ring or the integer ``1``. On the other hand, elements can also be specified by their exponent.

    INPUT:

        - ``parent`` -- The parent of the element.
        - ``x`` -- an expression (string, polynomial ring element, symbolic ring element, or the integer ``1``)
            representing the element to be initialized.
        - ``exponent`` -- the exponent of the asymptotic element, directly controlling the growth of the asymptotic
            growth element.

    OUTPUT:

        Returns a univariate asymptotic growth element with the specified parent and magnitude of growth.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as ar
        sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
        sage: e1 = P(x=1); e1
        1
        sage: e2 = P(x=None, exponent=2); e2
        x^2
        sage: e1 == e2
        False
        sage: P.le(e1, e2)
        True
    """
    # TODO: implement comparison for the cartesian product structure.

    def __init__(self, parent, x=None, exponent=None):
        r"""
        See :class:`AsymptoticGrowthElementUnivariate` for more information.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P.gen(); e1
            x
            sage: e1.is_idempotent()
            False
            sage: e1.is_one()
            False
            sage: e1.parent()
            Univariate Asymptotic Growth Group in x
            sage: e2 = P.one(); e2
            1
            sage: e2.is_idempotent() and e2.is_one()
            True
        """
        if x is None and exponent is None:
            raise ValueError("Neither x nor exponent are specified.")
        elif x is not None and exponent is not None:
            raise ValueError("Both x and exponent are specified.")
        elif exponent is None:
            if x == 1:
                self.exponent = 0
            else:
                print(x)
                raise NotImplementedError("Parsing of x is NYI")
        elif exponent is not None:
            if exponent not in QQ:
                raise TypeError("Non-rational exponents are not supported")
            else:
                self.exponent = exponent
        super(AsymptoticGrowthElementUnivariate, self).__init__(parent=parent)

    def _repr_(self):
        r"""
        Represents the univariate asymptotic growth element as ``variable^{exponent}``.

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
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
        Multiplies two univariate asymptotic growth elements from the same parent by adding their exponents.

        INPUT:

            - ``other`` -- the asymptotic growth element to be multiplied with ``self``.

        OUTPUT:

        A univariate asymptotic growth element representing the product of ``self`` and ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = P(x=None, exponent=3)
            sage: e3 = e1._mul_(e2); e3
            x^5
            sage: e3 == e1 * e2
            True
        """
        C = self.__class__
        return C(self.parent(), exponent=self.exponent + other.exponent)

    def __invert__(self):
        r"""
        Returns the multiplicative inverse from a given univariate asymptotic growth element.

        INPUT:

        Nothing.

        OUTPUT:

        The multiplicative inverse univariate asymptotic growth element of ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e2 = e1.__invert__(); e2
            x^-2
            sage: e2 == ~e1
            True
        """
        C = self.__class__
        return C(self.parent(), exponent=-self.exponent)

    def _div_(self, other):
        r"""
        Divides two univariate asymptotic growth elements from the same parent by subtracting their exponents.

        INPUT:

            - ``other`` -- the asymptotic growth element which ``self`` is divided by.

        OUTPUT:

        The result of the division of ``self`` by ``other``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P(x=None, exponent=2)
            sage: e1._div_(P.gen())
            x
            sage: e1._div_(P.gen()) == e1 / P.gen()
            True
        """
        C = self.__class__
        return C(self.parent(), exponent=self.exponent - other.exponent)
    __div__ = _div_

    def _pow_(self, power):
        r"""
        Returns a univariate asymptotic element to the power of ``power``.

        INPUT:

            - ``power`` -- a rational number.

        OUTPUT:

        The univariate asymptotic element ``self`` to the power of ``power``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: P.gen()._pow_(5)
            x^5
            sage: P.gen()._pow_(1/2)
            x^(1/2)
            sage: P.gen()^7
            x^7
        """
        C = self.__class__
        if power not in QQ:
            raise TypeError("Non-rational exponents are not supported")
        else:
            return C(self.parent(), exponent=self.exponent * power)
    __pow__ = _pow_

    def __eq__(self, other):
        r"""
        Returns whether the univariate asymptotic growth elements ``self`` and ``other`` are equal and have the same
        parent.

        INPUT:

            - ``other`` -- a univariate asymptotic growth element to be compared to ``self``.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P(x=None, exponent=1)
            sage: e1.__eq__(P.gen())
            True
            sage: e2 = e1^4
            sage: e2 == e1^2 * e1 * e1
            True
            sage: e2 == e1
            False
        """
        return self.parent() is other.parent() and self.exponent == other.exponent

    def is_le_one(self):
        r"""
        Returns whether or not the univariate asymptotic growth element ``self`` is in O(1).

        INPUT:

        Nothing.

        OUTPUT:

        A boolean.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x")
            sage: e1 = P.gen()
            sage: e1.is_le_one()
            False
            sage: (P.one() / P.gen()).is_le_one()
            True
        """
        return self.exponent <= 0


class AsymptoticGrowthGroupParentUnivariate(AsymptoticGrowthGroupParent):
    r"""
    Class for the concrete realization of asymptotic growth group parents in the univariate case. These are the
    parents for the :class:`AsymptoticGrowthElementUnivariate` elements. These parents contain single polynomial terms
    in a specified variable.

    More complex univariate growth elements can be constructed by constructing the cartesian product of various
    univariate asymptotic growth group parents in the same variable (which can be constructed with
    the :meth:`sage.groups.asymptotic_growth_group.AsymptoticGrowthGroupParentUnivariate.create_exp_parent` and
    :meth:`sage.groups.asymptotic_growth_group.AsymptoticGrowthGroupParentUnivariate.create_log_parent` method.

    INPUT:

        - ``variable`` -- either a symbol from the symbolic ring, a generator from a polynomial ring, or a alphanumeric
          string starting with a letter and optionally containing underscores.

        - ``category`` -- The category of the parent can be specified in order to broaden the base structure. Has to
          be a subcategory of ``Join of Category of groups and Category of posets``.

    OUTPUT:

    Returns a univariate asymptotic growth group parent.

    EXAMPLES::

        sage: import sage.groups.asymptotic_growth_group as ar
        sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x"); P
        Univariate Asymptotic Growth Group in x
    """
    # TODO: implement the cartesian product structure, as well as the exp_parent and log_parent methods.

    def __init__(self, variable, category=None):
        r"""
        For more information see :class:`AsymptoticGrowthGroupParentUnivariate`.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P1 = ar.AsymptoticGrowthGroupParentUnivariate("x"); P1
            Univariate Asymptotic Growth Group in x
            sage: var('n')
            n
            sage: P2 = ar.AsymptoticGrowthGroupParentUnivariate(n); P2
            Univariate Asymptotic Growth Group in n
            sage: y = PolynomialRing(ZZ, "y").gen()
            sage: P3 = ar.AsymptoticGrowthGroupParentUnivariate(y); P3
            Univariate Asymptotic Growth Group in y
            sage: P4 = ar.AsymptoticGrowthGroupParentUnivariate("y"); P4
            Univariate Asymptotic Growth Group in y
        """
        if parent(variable) is SR and variable.is_symbol():
            self.variable = repr(variable)
        elif hasattr(variable, "is_gen") and variable.is_gen():
            self.variable = repr(variable)
        elif variable is None:
            raise ValueError("Variable for initialization required.")
        else:
            if re.match("^[A-Za-z][A-Za-z0-9_]*$", variable):
                self.variable = str(variable)
            else:
                raise ValueError("Only alphanumeric strings starting with a letter may be variables.")
        super(AsymptoticGrowthGroupParentUnivariate, self).__init__(category=category)

    # enable the category framework for elements
    Element = AsymptoticGrowthElementUnivariate

    def _repr_(self):
        r"""
        Represents the univariate asymptotic growth group as "Univariate Asymptotic Growth Group in ``variable``".

        INPUT:

        Nothing.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: ar.AsymptoticGrowthGroupParentUnivariate("x")._repr_()
            'Univariate Asymptotic Growth Group in x'
            sage: ar.AsymptoticGrowthGroupParentUnivariate("v_107")._repr_()
            'Univariate Asymptotic Growth Group in v_107'
        """
        return "Univariate Asymptotic Growth Group in %s" % self.variable

    def one(self):
        r"""
        Returns the neutral element of the univariate asymptotic growth group.

        INPUT:

        Nothing.

        OUTPUT:

        The neutral element of the univariate asymptotic growth group ``self``.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: e1 = ar.AsymptoticGrowthGroupParentUnivariate("x").one(); e1
            1
            sage: e1.is_idempotent()
            True
        """
        return self(x=None, exponent=0)

    def gen(self):
        r"""
        Returns the univariate asymptotic growth element with exponent 1.

        INPUT:

        Nothing.

        OUTPUT:

        The univariate asymptotic growth element with exponent 1.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: e1 = ar.AsymptoticGrowthGroupParentUnivariate("x").gen(); e1
            x
            sage: e1.exponent == 1
            True
        """
        return self(x=None, exponent=1)

    def create_exp_parent(self, base=e):
        # TODO: Not entirely sure that this is the best approach. Has to be discussed.
        r"""
        Returns a univariate asymptotic growth group which can be used to contain elements of the form
        ``base^(variable^exponent)``. By calling this method on the univariate asymptotic growth group constructed
        in this manor, iterated exponentiation can also be modeled.

        .. TODO::

            Essentially, this should only modify the :meth:`_repr_` method of the univariate asymptotic growth group
            and its elements. Comparison between different parents is then realized by comparison over the respective
            cartesian product.

        INPUT:

            - ``base`` -- a real number, the base of the exponentiation. Defaults to ``e``.

        OUTPUT:

        A univariate asymptotic growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x").create_exp_parent(); P  # TODO: not implemented
            Univariate Asymptotic Growth Group in x with exponential degree 1
            sage: P(x=None, exponent=2)  # TODO: not implemented
            exp(x^2)
        """
        raise NotImplementedError("NYI")

    def create_log_parent(self):
        # TODO: see create_exp_parent above.
        r"""
        Returns a univariate asymptotic growth group which can be used to contain elements of the form
        ``log^{exponent}(variable)``. By calling this method on the univariate asymptotic growth group constructed
        in this manor, iterated logarithms can also be modeled.

        .. TODO::

            See :meth:`create_exp_parent`.

        INPUT:

        Nothing.

        OUTPUT:

        A univariate asymptotic growth group.

        EXAMPLES::

            sage: import sage.groups.asymptotic_growth_group as ar
            sage: P = ar.AsymptoticGrowthGroupParentUnivariate("x").create_log_parent(); P  # TODO: not implemented
            Univariate Asymptotic Growth Group in x with logarithmic degree 1
            sage: P(x=None, exponent=2)  # TODO: not implemented
            log^2(x)
        """
        raise NotImplementedError("NYI")

    def cartesian_product(self, other):
        # TODO: This has to be thoroughly discussed!
        r"""
        Returns the cartesian product of ``self`` and ``other`` (possibly also a list of univariate asymptotic growth
        groups).

        .. TODO::

            Thoroughly discuss the implementation of this. The standard framework for cartesian products probably cannot
            be used due to the fact that the :class:`sage.sets.cartesian_product.CartesianProduct` has no function for
            comparing elements of the resulting parent.

        INPUT:

            - ``other`` -- a univariate asymptotic growth group parent.

        OUTPUT:

        The cartesian product of ``self`` and ``other``.

        EXAMPLES::

            sage: ???  # TODO: not implemented
        """
        raise NotImplementedError("NYI")


