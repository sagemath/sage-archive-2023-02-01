from sage.misc.all import cached_method
from sage.categories.all import Monoids
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.sets.family import Family
from sage.rings.integer import Integer
import operator


class AutomaticMonoid(UniqueRepresentation, Parent):
    r"""

    Construct (lazily) a monoid from a set of concrete generators
    living in an ambient monoid.

    EXAMPLES::

        sage: R = IntegerModRing(12)
        sage: M = AutomaticMonoid(Family({1: R(3), 2: R(5)}), one = R.one())
        sage: M.one()
        []
        sage: M.one() in M
        True

    Elements are represented by default by their reduced word, or by
    the corresponding element in the ambient monoid if the reduced
    word is not yet known.

    ::

        sage: g = M.generators; g
        Finite family {1: 3, 2: 5}
        sage: g[1]*g[2]
        3

    Calling cardinality, or list, or iterating through the monoid will
    trigger its full construction and, as a side effect, compute all
    the reduced words. The order of the elements, and the induced
    choice of reduced word is currently length-lexicographic
    (i.e. the chosen reduced word is of minimal length, and then
    minimal lexicographically w.r.t. the order of the indices of the
    generators).

        sage: M.cardinality()
        4
        sage: M.list()
        [[], [1], [2], [1, 1]]
        sage: g = M.generators

        sage: g[1]*g[2]
        [1]

        sage: g[1].transition(1)
        [1, 1]
        sage: g[1] * g[1]
        [1, 1]
        sage: g[1] * g[1] * g[1]
        [1]
        sage: g[1].transition(2)
        [1]
        sage: g[1] * g[2]
        [1]
        sage: M.list()
        [[], [1], [2], [1, 1]]
        sage: [ x.lift() for x in M.list() ]
        [1, 3, 5, 9]
        sage: M.idempotents()
        [[], [1, 1]]
        sage: G = M.cayley_graph(side = "twosided"); G
        Looped multi-digraph on 4 vertices
        sage: G.edges()
        [([], [1], (1, 'left')), ([], [1], (1, 'right')), ([], [2], (2, 'left')), ([], [2], (2, 'right')), ([1], [1], (2, 'left')), ([1], [1], (2, 'right')), ([1], [1, 1], (1, 'left')), ([1], [1, 1], (1, 'right')), ([2], [], (2, 'left')), ([2], [], (2, 'right')), ([2], [1], (1, 'left')), ([2], [1], (1, 'right')), ([1, 1], [1], (1, 'left')), ([1, 1], [1], (1, 'right')), ([1, 1], [1, 1], (2, 'left')), ([1, 1], [1, 1], (2, 'right'))]

        sage: map(sorted, M.j_classes())
        [[[], [2]], [[1, 1], [1]]]
        sage: M.j_classes_of_idempotents()
        [[[]], [[1, 1]]]
        sage: M.j_transversal_of_idempotents()
        [[], [1, 1]]

        sage: map(attrcall('pseudo_order'), M.list())
        [[1, 0], [3, 1], [2, 0], [2, 1]]

    TESTS::

        sage: (g[1]).__hash__() == (g[1]*g[1]*g[1]).__hash__()
        True
        sage: g[1] == g[1]*g[1]*g[1]
        True
        sage: copy(g[1]) is g[1]
        True
        sage: from copy import deepcopy
        sage: deepcopy(g[1]) is g[1]
        True
        sage: M.__class__.mro() # todo: not implemented (wait until CombinatorialClass is a category)
        [<class 'sage.combinat.automatic_monoid.AutomaticMonoid_with_categories'>, <class 'sage.categories.category.uniq'>, <type 'sage.structure.parent.Parent'>, <type 'sage.structure.category_object.CategoryObject'>, <type 'sage.structure.sage_object.SageObject'>, <class 'sage.categories.finite_monoids.FiniteMonoids.parent_class'>, <class 'sage.categories.category.Monoids.parent_class'>, <class 'sage.categories.category.Semigroups.parent_class'>, <class 'sage.categories.category_types.Sets.parent_class'>, <type 'object'>]

        sage: TestSuite(M).run()

    CAVEATS:

     - AutomaticMonoid is designed primarily for finite monoids. This
       property is not checked automatically (this would be too costly
       if not impossible). Some of the features should still work with
       infinite monoids. In that case, the category Monoids() should
       be passed as extra argument, instead of the default
       Monoids().Finite().

        sage: AutomaticMonoid(Family({1:2}), category = Monoids().Finite())
        The (automatic) monoid with generators Finite family {1: 2}

    BUGS:

    - Running a command like [ f(p) for p in M ] where f involves
      computing products of elements and M.list() has not yet been
      computed caused an infinite loop at least in one occasion looking like::

         sage: all(operator.sub(*p.pseudo_order())==1 for p in M)
         False

    """

    @staticmethod
    def __classcall_private__(cls, generators, one=None, mul=operator.mul, category=Monoids().Finite()):
        generators = Family(generators)
        return super(AutomaticMonoid, cls).__classcall__(cls, generators, one=one, mul=mul, category=category)

    def __init__(self, generators, one=None, mul=operator.mul, category=Monoids().Finite()):
        """
        TESTS::

            sage: R = IntegerModRing(12)
            sage: M = AutomaticMonoid(Family(()), one = R.one())
            sage: M.ambient() == R
            True
            sage: M = AutomaticMonoid(Family(()))
            Traceback (most recent call last):
            ...
            ValueError: AutomaticMonoid requires at least one generator or `one` to determine the ambient space

        """
        Parent.__init__(self, category=category)
        self.generators = generators                    # todo: rename to self._generators?
        if self.generators.cardinality() > 0:
            self._ambient = self.generators.first().parent()
        elif one is not None:
            self._ambient = one.parent()
        else:
            raise ValueError("AutomaticMonoid requires at least one generator or `one` to determine the ambient space")

        if one is None:
            one = self.generators.first().parent().one()
        self._one = self.retract(one)
        self.mul = mul
        self.generators_in_ambient = generators
        self.generators = generators.map(self.retract)

        self.elements     = [self.one()]
        self.elements_set = set(self.elements)
        self.one()._reduced_word = []
        self.done         = 0

    def one(self):
        return self._one

    def _repr_(self):
        return "The (automatic) monoid with generators %s"%self.generators

    def an_element(self):
        return self.generators.first()

    def some_elements(self):
        return self.generators

    def ambient(self):
        return self._ambient

    @cached_method
    def retract(self, ambient_element):
        return self.element_class(self, ambient_element)

    def lift(self, x):
        assert(x in self)
        return x.lift()

    def semigroup_generators(self):
        return self.generators

    def __iter__(self):
        for x in self.elements:
            yield x
        while self.done < len(self.elements):
            x = self.elements[self.done]
            for i in self.generators.keys():
                y = x.transition(i)
                if y in self.elements_set:
                    continue
                self.elements.append(y)
                self.elements_set.add(y)
                y._reduced_word = x.reduced_word()+[i]
                yield y
            self.done+=1
        return

    def cardinality(self):
        """
            sage: R = IntegerModRing(12)
            sage: M = AutomaticMonoid(Family({1: R(3), 2: R(5)}), one = R.one())
            sage: M.cardinality()
            4

        TESTS::

            sage: assert isinstance(M.cardinality(), Integer)  # This did fail at some point
        """
        if self.done == len(self.elements):
            return Integer(self.done)
        else:
            return self._cardinality_from_iterator()

    def __contains__(self, x):
        return x.parent() is self

    def product(self, x, y):
        assert(x in self)
        assert(y in self)
        red = y.reduced_word()
        if red is None:
            return self.retract(self.mul(x.lift(), y.lift()))
        else:
            for i in red:
                x = x.transition(i)
        return x

    def from_reduced_word(self, l):
        result = self.one()
        for i in l:
            result = result.transition(i)
        return result

    class Element(ElementWrapper):
        def __init__(self, ambient_element, parent):
            ElementWrapper.__init__(self, ambient_element, parent)
            self._reduced_word = None

        # should we inherit from UniqueRepresentation instead?
        def __eq__      (self, other): return self is other
        def __hash__    (self): return id(self)
        def __copy__    (self): return self
        def __deepcopy__(self, memo): return self

        def __reduce__(self):
            """
            sage: R = IntegerModRing(12)
            sage: M = AutomaticMonoid(Family({1: R(3), 2: R(5)}), one = R.one())
            sage: bla = M.list()[2]
            sage: loads(dumps(bla)) is bla   # indirect doctest
            True
            """
            return (unreduce_retract, (self.parent(), self.value,))

        def lift(self):
            return self.value

        def reduced_word(self):
            """
            Returns the length-lexicographic shortest reduced word for
            self, or None if it has not yet been computed

            TODO: should it instead trigger the computation?
            """
            return self._reduced_word

        @cached_method
        def transition(self, i):
            """
            The multiplication on the right by a generator.
            Namely, this returns x * self.generators[i]
            """
            parent = self.parent()
            assert(i in parent.generators.keys())
            # Do it in to steps, as a temporary workaround for #5843
            res = parent.mul(self.lift(), parent.generators_in_ambient[i])
            return parent.retract(res)

        def _repr_(self):
            rep = self.reduced_word()
            if rep is None:
                rep = self.lift()
            return "%s"%rep


def unreduce_retract(parent, ambient_element):
    return parent.retract(ambient_element)
