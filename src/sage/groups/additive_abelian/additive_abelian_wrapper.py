r"""
Wrapper class for abelian groups

This class is intended as a template for anything in Sage that needs the
functionality of abelian groups. One can create an AdditiveAbelianGroupWrapper
object from any given set of elements in some given parent, as long as an
``_add_`` method has been defined.


EXAMPLES:

We create a toy example based on the Mordell-Weil group of an elliptic curve over `\QQ`::

    sage: E = EllipticCurve('30a2')
    sage: pts = [E(4,-7,1), E(7/4, -11/8, 1), E(3, -2, 1)]
    sage: M = AdditiveAbelianGroupWrapper(pts[0].parent(), pts, [3, 2, 2])
    sage: M
    Additive abelian group isomorphic to Z/2 + Z/6 embedded in Abelian group of
    points on Elliptic Curve defined by y^2 + x*y + y = x^3 - 19*x + 26 over
    Rational Field
    sage: M.gens()
    ((4 : -7 : 1), (7/4 : -11/8 : 1), (3 : -2 : 1))
    sage: 3*M.0
    (0 : 1 : 0)
    sage: 3000000000000001 * M.0
    (4 : -7 : 1)

We check that ridiculous operations are being avoided::

    sage: set_verbose(2, 'additive_abelian_wrapper.py')
    sage: 300001 * M.0
    verbose 1 (...: additive_abelian_wrapper.py, _discrete_exp) Calling discrete exp on (1, 0, 0)
    (4 : -7 : 1)
    sage: set_verbose(0, 'additive_abelian_wrapper.py')


TODO:

- Implement proper black-box discrete logarithm (using baby-step giant-step).
  The discrete_exp function can also potentially be speeded up substantially
  via caching.

- Think about subgroups and quotients, which probably won't work in the current
  implementation -- some fiddly adjustments will be needed in order to be able
  to pass extra arguments to the subquotient's init method.
"""

import additive_abelian_group as addgp
from sage.rings.all import ZZ
from sage.misc.misc import verbose
from sage.categories.morphism import Morphism

class UnwrappingMorphism(Morphism):
    r"""
    The embedding into the ambient group. Used by the coercion framework.
    """
    def __init__(self, domain):
        Morphism.__init__(self, domain.Hom(domain.universe()))

    def _call_(self, x):
        r"""
        TEST::

            sage: E = EllipticCurve("65a1")
            sage: G = E.torsion_subgroup()
            sage: isinstance(G, sage.groups.additive_abelian.additive_abelian_wrapper.AdditiveAbelianGroupWrapper)
            True
            sage: P1 = E([1,-1,1])
            sage: P2 = E([0,1,0])
            sage: P1 in G
            False
            sage: P2 in G
            True
            sage: (G(P2) + P1) in G
            False
            sage: (G(P2) + P1).parent()
            Abelian group of points on Elliptic Curve defined by y^2 + x*y = x^3 - x over Rational Field
        """
        return self.codomain()(x.element())

class AdditiveAbelianGroupWrapper(addgp.AdditiveAbelianGroup_fixed_gens):

    def __init__(self, universe, gens, invariants):
        self._universe = universe
        self._gen_elements = tuple(universe(x) for x in gens)
        cover,rels = addgp.cover_and_relations_from_invariants(invariants)
        addgp.AdditiveAbelianGroup_fixed_gens.__init__(self, cover, rels, cover.gens())
        self._gen_orders = invariants
        try:
            self.register_embedding(UnwrappingMorphism(self))
        except AssertionError:
            # coercion already exists -- this can happen, and we can't rely on
            # has_coerce_map_from at this stage in the game (try it, it doesn't
            # work!). So we just catch the AssertionError and ignore it.
            pass

    def universe(self):
        return self._universe

    def generator_orders(self):
        return tuple(self._gen_orders)

    def _repr_(self):
        return addgp.AdditiveAbelianGroup_fixed_gens._repr_(self) + " embedded in " + self.universe()._repr_()

    def _element_class(self):
        return AdditiveAbelianGroupWrapperElement

    def _discrete_exp(self, v):
        v = self.V()(v)
        verbose("Calling discrete exp on %s" % v)
        # DUMB IMPLEMENTATION!
        return sum([self._gen_elements[i] * ZZ(v[i]) for i in xrange(len(v))], self.universe()(0))

    def _discrete_log(self,x):
        # EVEN DUMBER IMPLEMENTATION!
        u = [y for y in self.list() if y.element() == x]
        if len(u) == 0: raise TypeError, "Not in group"
        if len(u) > 1: raise NotImplementedError
        return u[0]

    def __call__(self, x, check=False):
        if hasattr(x,"parent"):
            if x.parent() is self:
                return x
            elif x.parent() is self.universe():
                return AdditiveAbelianGroupWrapperElement(self, self._discrete_log(x), element = x)
        return addgp.AdditiveAbelianGroup_fixed_gens.__call__(self, x, check)

class AdditiveAbelianGroupWrapperElement(addgp.AdditiveAbelianGroupElement):

    def __init__(self, parent, vector, element=None, check=False):
        addgp.AdditiveAbelianGroupElement.__init__(self, parent, vector, check)
        if element is not None:
            element = self.parent().universe()(element)
        self._element = element

    def element(self):
        if self._element is None:
            self._element = self.parent()._discrete_exp(self._hermite_lift())
        return self._element

    def _repr_(self):
        return repr(self.element())

