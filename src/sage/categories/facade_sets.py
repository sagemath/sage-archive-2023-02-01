r"""
Facade Sets
"""
#*****************************************************************************
#  Copyright (C) 2010 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category import Category

class FacadeSets(Category_singleton):
    r"""
    The category of facade sets

    A *facade set* is a parent ``P`` whose elements actually belong to
    some other parent::

        sage: P = Sets().example(); P
        Set of prime numbers (basic implementation)
        sage: p = Sets().example().an_element(); p
        47
        sage: p in P
        True
        sage: p.parent()
        Integer Ring

    Typical use cases include modeling a subset of an existing
    parent::

        sage: Sets().Facades().example()
        An example of facade set: the monoid of positive integers

    or the union of several parents::

        sage: Sets().Facades().example("union")
        An example of a facade set: the integers completed by +-infinity

    or endowing a parent with more (or less!) structure::

        sage: Posets().example("facade")
        An example of a facade poset: the positive integers ordered by divisibility

    Let us consider one of the examples above in detail: the partially ordered
    set `P` of positive integers w.r.t. divisibility order. There are two
    options for representing its elements:

     1. as plain integers
     2. as integers, modified to be aware that their parent is `P`

    The advantage of 1. is that one needs not to do conversions back and
    forth between `P` and `\ZZ`. The disadvantage is that this
    introduces an ambiguity when writing `2 < 3`::

        sage: 2 < 3
        True

    To raise this ambiguity, one needs to explicitely specify the order
    as in `2 <_P 3`::


        sage: P = Posets().example("facade")
        sage: P.lt(2,3)
        False

    In short `P` being a facade parent is one of the programmatic
    counterpart (with e.g. coercions) of the usual mathematical idiom:
    "for ease of notation, we identify an element of `P` with the
    corresponding integer". Too many identifications lead to
    confusion; the lack thereof leads to heavy, if not obfuscated,
    notations. Finding the right balance is an art, and even though
    there are common guidelines, it is ultimately up to the writer to
    choose which identifications to do. This is no different in code.

    .. seealso::

       ::

        sage: Sets().example("facade")
        Set of prime numbers (facade implementation)
        sage: Sets().example("inherits")
        Set of prime numbers
        sage: Sets().example("wrapper")
        Set of prime numbers (wrapper implementation)

    .. rubric:: Specifications

    A parent which is a facade must either:

    - call :meth:`Parent.__init__` using the ``facade`` parameter to
      specify a parent, or tuple thereof.
    - overload the method :meth:`~FacadeSets.ParentMethods.facade_for`.

    .. note:: the concept of facade parents was originally introduced
       in the computer algebra system MuPAD.

    TESTS:

    Check that multiple categories initialisation works (:trac:`13801`)::

        sage: class A(Parent):
        ...     def __init__(self):
        ...         Parent.__init__(self, category=(FiniteEnumeratedSets(),Monoids()), facade=True)
        sage: a = A()
    """

    @cached_method
    def super_categories(self):
        r"""
        Returns the super categories of ``self``, as per
        :meth:`Category.super_categories`.

        EXAMPLES::

            sage: Sets().Facades().super_categories()
            [Category of sets]
        """
        from sage.categories.sets_cat import Sets
        return [Sets()]

    def example(self, choice='subset'):
        r"""
        Returns an example of facade set, as per
        :meth:`Category.example()
        <sage.categories.category.Category.example>`.

        INPUT:

        - ``choice`` -- 'union' or 'subset' (default: 'subset').

        EXAMPLES::

            sage: Sets().Facades().example()
            An example of facade set: the monoid of positive integers
            sage: Sets().Facades().example(choice='union')
            An example of a facade set: the integers completed by +-infinity
            sage: Sets().Facades().example(choice='subset')
            An example of facade set: the monoid of positive integers
        """
        import sage.categories.examples.facade_sets as examples
        if choice == "union":
            return examples.IntegersCompletion()
        elif choice == 'subset':
            return examples.PositiveIntegerMonoid()
        else:
            raise TypeError, "choice should be 'union' or 'subset'"

    class ParentMethods:

        def _element_constructor_(self, element):
            """
            Coerce ``element`` into ``self``

            INPUT:

            - ``element`` -- any object

            This default implementation returns ``element`` if
            ``self`` is a facade for ``parent(element)`. Otherwise it
            attempts in turn to coerce ``element`` into each parent
            ``self`` is a facade for.

            This implementation is only valid for a facade parent
            which models the full union of the parents it is a facade
            for. Other facade parents should redefine
            :meth:`element_constructor` appropriately.

            EXAMPLES::

                sage: S = Sets().Facades().example("union"); S
                An example of a facade set: the integers completed by +-infinity
                sage: S(1)
                1
                sage: S(1/2)
                Traceback (most recent call last):
                ...
                ValueError: Can't coerce `1/2` in any parent `An example of a facade set: the integers completed by +-infinity` is a facade for
                sage: S(2/1)
                2
                sage: S(2/1).parent()
                Integer Ring
                sage: S(int(1))
                1
                sage: S(int(1)).parent()
                Integer Ring

            Facade parents that model strict subsets should redefine
            :meth:`element_constructor`::

                sage: S = Sets().Facades().example(); S
                An example of facade set: the monoid of positive integers
                sage: S(-1)
                Traceback (most recent call last):
                ...
                ValueError: %s should be positive
            """
            if self.is_parent_of(element):
                return element
            else:
                parents = self.facade_for()
                if parents is True:
                    return NotImplementedError
                for parent in self.facade_for():
                    try:
                        return parent(element)
                    except Exception:
                        pass
            raise ValueError, "Can't coerce `%s` in any parent `%s` is a facade for"%(element, self)

        def facade_for(self):
            """
            Returns the parents this set is a facade for

            This default implementation assumes that ``self`` has
            an attribute ``_facade_for``, typically initialized by
            :meth:`Parent.__init__`. If the attribute is not present, the method
            raises a NotImplementedError.

            EXAMPLES::

                sage: S = Sets().Facades().example(); S
                An example of facade set: the monoid of positive integers
                sage: S.facade_for()
                (Integer Ring,)

            Check that :trac:`13801` is corrected::

                sage: class A(Parent):
                ...     def __init__(self):
                ...         Parent.__init__(self, category=Sets(), facade=True)
                sage: a = A()
                sage: a.facade_for()
                Traceback (most recent call last):
                ...
                NotImplementedError: this parent did not specify which parents it is a facade for
            """
            try:
                return self._facade_for
            except AttributeError:
                raise NotImplementedError("this parent did not specify which parents it is a facade for")

        def is_parent_of(self, element):
            """
            Returns whether ``self`` is the parent of ``element``

            INPUT:

            - ``element`` -- any object

            Since ``self`` is a facade domain, this actually tests
            whether the parent of ``element`` is any of the parent
            ``self`` is a facade for.

            EXAMPLES::

                sage: S = Sets().Facades().example(); S
                An example of facade set: the monoid of positive integers
                sage: S.is_parent_of(1)
                True
                sage: S.is_parent_of(1/2)
                False

            This method differs from :meth:`__contains__` in two
            ways.  First, this does not take into account the fact
            that ``self`` may be a strict subset of the parent(s)
            it is a facade for::

                sage: -1 in S, S.is_parent_of(-1)
                (False, True)

            Furthermore, there is no coercion attempted::

                sage: int(1) in S, S.is_parent_of(int(1))
                (True, False)

            .. warning::

               this implementation does not handle facade parents of facade
               parents. Is this a feature we want generically?
            """
            parents = self.facade_for()
            if parents is True:
                return True
            from sage.structure.element import parent
            return parent(element) in parents

        def __contains__(self, element):
            """
            Membership testing

            Returns whether ``element`` is in one of the parents
            ``self`` is a facade for.

            .. warning:: this default implementation is currently
            overriden by :meth:`Parent.__contains__`.

            EXAMPLES::

                sage: S = Sets().Facades().example("union"); S
                An example of a facade set: the integers completed by +-infinity
                sage: 1 in S, -5 in S, oo in S, -oo in S, int(1) in S, 2/1 in S
                (True, True, True, True, True, True)
                sage: 1/2 in S, "bla" in S
                (False, False)
            """
            return any(element in parent for parent in self.facade_for())

        def _an_element_(self):
            """
            Try to return an element of ``self``, as per
            :meth:`Sets.ParentMethods.an_element`.

            For each parent ``self`` is a facade for, this default
            implementation tries the method ``an_element`` until it finds an
            element in ``self``. If none is found raise a
            ``NotImplementedError``.

            EXAMPLES::

                sage: S = Sets().Facades().example(); S
                An example of facade set: the monoid of positive integers
                sage: S.an_element()
                1
            """
            for parent in self.facade_for():
                x = parent.an_element()
                if x in self:
                    return x
            raise NotImplementedError


