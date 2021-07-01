r"""
Subsets of a Universe Defined by Predicates
"""

from sage.structure.category_object import normalize_names
from sage.structure.parent import Parent, Set_generic
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.symbolic.expression import is_Expression
from sage.symbolic.callable import is_CallableSymbolicExpression
from sage.symbolic.ring import SymbolicRing, SR, is_SymbolicVariable

from .set import Set, Set_base, Set_boolean_operators, Set_add_sub_operators

class ConditionSet(Set_generic, Set_base, Set_boolean_operators, Set_add_sub_operators,
                   UniqueRepresentation):
    r"""
    Set of elements of a universe that satisfy given predicates

    EXAMPLES::

        sage: Evens = ConditionSet(ZZ, is_even); Evens
        { x ∈ Integer Ring : <function is_even at 0x...>(x) }
        sage: 2 in Evens
        True
        sage: 3 in Evens
        False
        sage: 2.0 in Evens
        True

        sage: Odds = ConditionSet(ZZ, is_odd); Odds
        { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
        sage: EvensAndOdds = Evens | Odds; EvensAndOdds
        Set-theoretic union of
         { x ∈ Integer Ring : <function is_even at 0x...>(x) } and
         { x ∈ Integer Ring : <function is_odd at 0x...>(x) }
        sage: 5 in EvensAndOdds
        True
        sage: 7/2 in EvensAndOdds
        False

        sage: var('y')
        y
        sage: SmallOdds = ConditionSet(ZZ, is_odd, abs(y) <= 11, vars=[y]); SmallOdds
        { y ∈ Integer Ring : abs(y) <= 11, <function is_odd at 0x3c6cb8310>(y) }

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
        sage: vector([1, 0, 0]) in P_inter_B_again
        True
        sage: vector([1, 1, 1]) in P_inter_B_again
        False

    TESTS::

        sage: TestSuite(Evens).run()
        sage: TestSuite(P_inter_B).run(skip='_test_pickling')  # cannot pickle lambdas
        sage: TestSuite(P_inter_B_again).run()
    """
    @staticmethod
    def __classcall_private__(cls, universe, *predicates, vars=None, names=None, category=None):
        if category is None:
            category = Sets()
        if isinstance(universe, Parent):
            if universe in Sets().Finite():
                category = category & Sets().Finite()

        if vars is not None:
            if names is not None:
                raise ValueError('cannot use names and vars at the same time; they are aliases')
            names, vars = vars, None

        if names is not None:
            names = normalize_names(-1, names)

        callable_symbolic_predicates = []
        other_predicates = []

        for predicate in predicates:
            if is_CallableSymbolicExpression(predicate):
                if names is None:
                    names = tuple(str(var) for var in predicate.args())
                elif len(names) != predicates.args():
                    raise TypeError('mismatch in number of arguments')
                callable_symbolic_predicates.append(predicate)
            elif is_Expression(predicate):
                if names is None:
                    raise TypeError('use callable symbolic expressions or provide variable names')
                if vars is None:
                    vars = tuple(SR.var(name) for name in names)
                callable_symbolic_predicates.append(predicate.function(*vars))
            else:
                other_predicates.append(predicate)

        predicates = callable_symbolic_predicates + other_predicates

        if not other_predicates:
            if not callable_symbolic_predicates:
                if names is None and category is None:
                    # No conditions, no variable names, no category, just use Set.
                    return Set(universe)
            # Use ConditionSet_callable_symbolic_expression even if no conditions
            # are present; this will make the _sympy_ method available.
            return ConditionSet_callable_symbolic_expression(universe, *predicates,
                                                             names=names, category=category)

        if names is None:
            names = ("x",)
        return super().__classcall__(cls, universe, *predicates,
                                     names=names, category=category)

    def __init__(self, universe, *predicates, names=None, category=None):
        self._universe = universe
        self._predicates = predicates
        facade = None
        if isinstance(universe, Parent):
            facade = universe
        super().__init__(facade=facade, category=category,
                         names=names, normalize=False) # names already normalized by classcall

    def _repr_(self):
        s = "{ "
        names = self.variable_names()
        comma_sep_names = ", ".join(str(name) for name in names)
        if len(names) == 1:
            s += f"{comma_sep_names}"
        else:
            s += f"({comma_sep_names})"
        universe = self._universe
        s += f" ∈ {universe}"
        sep = " : "
        for predicate in self._predicates:
            s += sep + self._repr_condition(predicate)
            sep = ", "
        s += " }"
        return s

    @cached_method
    def _repr_condition(self, predicate):
        if is_CallableSymbolicExpression(predicate):
            args = self.arguments()
            if len(args) == 1:
                args = args[0]
            condition = self._call_predicate(predicate, args)
            return str(condition)
        comma_sep_names = ", ".join(str(name)
                                    for name in self.variable_names())
        return f"{predicate}({comma_sep_names})"

    @cached_method
    def arguments(self):
        return SR.var(self.variable_names())

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
        if not all(self._call_predicate(predicate, element)
                   for predicate in self._predicates):
            raise ValueError(f'{element} does not satisfy the condition')
        return element

    def _call_predicate(self, predicate, element):
        if is_CallableSymbolicExpression(predicate):
            if len(predicate.arguments()) != 1:
                return predicate(*element)
        return predicate(element)

    def _an_element_(self):
        for element in self._universe.some_elements():
            if element in self:
                return element
        raise NotImplementedError

    def ambient(self):
        return self._universe


class ConditionSet_callable_symbolic_expression(ConditionSet):

    @cached_method
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

            sage: Interval = ConditionSet(RR, x >= -7, x <= 4, vars=[x]); Interval
            { x ∈ Real Field with 53 bits of precision : x >= -7, x <= 4 }
            sage: Interval._sympy_()
            ConditionSet(x, (x >= -7) & (x <= 4), SageSet(Real Field with 53 bits of precision))
        """
        import sympy
        args = self.arguments()
        if len(args) == 1:
            args = args[0]
            sym = args._sympy_()
        else:
            sym = tuple(x._sympy_() for x in args)
        conditions = [self._call_predicate(predicate, args)
                      for predicate in self._predicates]
        return sympy.ConditionSet(sym,
                                  sympy.And(*[condition._sympy_()
                                              for condition in conditions]),
                                  base_set=self._universe._sympy_())
