r"""
Subsets of a Universe Defined by a Predicate
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.symbolic.callable import is_CallableSymbolicExpression
from sage.symbolic.ring import SymbolicRing

from .set import Set

class ConditionSet(Parent, UniqueRepresentation):
    r"""
    Set of elements of a universe that satisfy a predicate

    EXAMPLES::

        sage: Evens = ConditionSet(ZZ, is_even); Evens
        { x ∈ Integer Ring : <function is_even at 0x...>(x) }
        sage: 2 in Evens
        True
        sage: 3 in Evens
        False
        sage: 2.0 in Evens
        True

        sage: P = polytopes.cube(); P
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
        sage: P.rename("P")
        sage: P_inter_B = ConditionSet(P, lambda x: x.norm() < 1.2); P_inter_B
        { x ∈ P : <function <lambda> at 0x...>(x) }
        sage: vector([1, 0, 0]) in P_inter_B
        True
        sage: vector([1, 1, 1]) in P_inter_B
        False

        sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 1.2; predicate
        (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 1.20000000000000
        sage: P_inter_B_again = ConditionSet(P, predicate); P_inter_B_again
        { (x, y, z) ∈ P : sqrt(x^2 + y^2 + z^2) < 1.20000000000000 }
        sage: vector([1, 0, 0]) in P_inter_B
        True
        sage: vector([1, 1, 1]) in P_inter_B
        False
    """
    @staticmethod
    def __classcall_private__(cls, universe, predicate, category=None):
        if category is None:
            category = Sets()
        if isinstance(universe, Parent):
            if universe in Sets().Finite():
                category = category & Sets().Finite()
        if is_CallableSymbolicExpression(predicate):
            return ConditionSet_callable_symbolic_expression(universe, predicate, category)
        return super().__classcall__(cls, universe, predicate, category)

    def __init__(self, universe, predicate, category):
        self._universe = universe
        self._predicate = predicate
        facade = None
        if isinstance(universe, Parent):
            facade = universe
        super().__init__(facade=facade, category=category)

    def _repr_(self):
        universe = self._universe
        predicate = self._predicate
        var = 'x'
        return '{ ' + f'{var} ∈ {universe} : {predicate}({var})' + ' }'

    def _element_constructor_(self, *args, **kwds):
        try:
            universe_element_constructor = self._universe._element_constructor_
        except AttributeError:
            if len(args) != 1 or kwds:
                raise ValueError('element constructor only takes 1 argument')
            element = args[0]
            if element not in self._universe:
                raise ValueError(f'{element} is not an element of the universe')
        else:
            element = universe_element_constructor(*args, **kwds)
        if not self._predicate(element):
            raise ValueError(f'{element} does not satisfy the condition')
        return element

    def ambient(self):
        return self._universe

class ConditionSet_callable_symbolic_expression(ConditionSet):

    def _repr_(self):
        universe = self._universe
        predicate = self._predicate
        args = predicate.arguments()
        predicate_expr = SymbolicRing._repr_element_(predicate.parent(), predicate)
        return '{ ' + f'{args} ∈ {universe} : {predicate_expr}' + ' }'

    def _sympy_(self):
        r"""
        EXAMPLES::

            sage: predicate(x, y, z) = sqrt(x^2 + y^2 + z^2) < 12; predicate
            (x, y, z) |--> sqrt(x^2 + y^2 + z^2) < 12
            sage: (ZZ^3).rename('ZZ^3')
            sage: SmallTriples = ConditionSet(ZZ^3, predicate); SmallTriples
            { (x, y, z) ∈ ZZ^3 : sqrt(x^2 + y^2 + z^2) < 12 }
            sage: ST = SmallTriples._sympy_(); ST
            ConditionSet((x, y, z), sqrt(x**2 + y**2 + z**2) < 12,
                         ProductSet(Integers, Integers, Integers))
            sage: (1, 3, 5) in ST
            True
            sage: (5, 7, 9) in ST
            False
        """
        import sympy
        predicate = self._predicate
        args = predicate.arguments()
        if len(args) == 1:
            sym = args._sympy_()
        else:
            sym = tuple(x._sympy_() for x in args)
        return sympy.ConditionSet(sym,
                                  predicate._sympy_(),
                                  base_set=self._universe._sympy_())
