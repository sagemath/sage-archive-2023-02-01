r"""
Crystal of Bernstien-Zelevinsky Multisegments
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

class InfinityCrystalOfMultisegments(Parent, UniqueRepresentation):
    r"""
    The type `A_n^{(1)}` crystal `B(\infty)` realized using
    Bernstein-Zelevinsky (BZ) multisegments.

    INPUT:

    - ``n`` -- for type `A_n^{(1)}`

    EXAMPLES::

    REFERENCES:

    - [Vazirani2002]_
    - [TingleyLN]_

    .. [LTV1999] Bernard Leclerc, Jean-Yves Thibon, and Eric Vasserot.
       *Zelevinsky's involution at roots of unity*. J. Reine Angew. Math.
       513:33-51 (1999).
    """
    def __init__(self, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = InfinityCrystalOfMultisegments(2)
        """
        self._cartan_type = CartanType(['A', n, 1])
        Parent.__init__(self, category=(HighestWeightCrystals(), InfiniteEnumeratedSets()))
        self.module_generators = (self.module_generator(),)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "The infinity crystal of BZ-multisegments of type {}".format(self._cartan_type)

    def module_generator(self):
        """
        Return the module generator of ``self``.
        """
        return self.element_class(self, ())

    class Element(ElementWrapper):
        """
        An element in a BZ multisegments crystal.
        """
        def __init__(self, parent, value):
            """
            Initialize ``self``.
            """
            ElementWrapper.__init__(self, parent, tuple(sorted(value, reverse=True)))

        def e(self, i):
            r"""
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::
            """
            M = self.value
            pos = [] # The positions of the uncancelled minuses
            for j,(a,b) in enumerate(M):
                if a == i:
                    pos.append(j)
                elif a + 1 == i and len(pos) > 0:
                    pos.pop()

            if len(pos) == 0:
                return None
            j = pos[0]
            a,b = M[j]
            if b == 1:
                return self.__class__(self.parent(), M[:j] + M[j+1:])
            return self.__class__(self.parent(), M[:j] + ((a-1,b-1),) + M[j+1:])

        def f(self, i):
            r"""
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::
            """
            M = self.value
            pos = [] # The positions of the uncancelled minuses
            for j,(a,b) in reversed(list(enumerate(M))):
                if a + 1 == i:
                    pos.append(j)
                elif a == i and len(pos) > 0:
                    pos.pop()

            if len(pos) == 0:
                Z = IntegerModRing(self.cartan_type().rank())
                return self.__class__(self.parent(), ((Z(1), i),) + M)
            j = pos[0]
            a,b = M[j]
            return self.__class__(self.parent(), M[:j] + ((a+1,b+1),) + M[j+1:])

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set

            EXAMPLES::
            """
            M = self.value
            epsilon = 0
            for j,(a,b) in enumerate(M):
                if b == i:
                    epsilon += 1
                elif b == i + 1 and epsilon > 0:
                    epsilon -= 1

            return epsilon

        def phi(self,i):
            r"""
            Return `\varphi_i` of ``self``.

            Let `T \in \mathcal{B}(\infty)` Define `\varphi_i(T) :=
            \varepsilon_i(T) + \langle h_i, \mathrm{wt}(T) \rangle`, where `h_i`
            is the `i`-th simple coroot and `\mathrm{wt}(T)` is the :meth:`weight`
            of `T`.

            INPUT:

            - ``i`` -- an element of the index set
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return self.epsilon(i) + P(self.weight()).scalar(h[i])

        def weight(self):
            """
            Return the weight of ``self``.
            """
            WLR = self.parent().weight_lattice_reallization()
            alpha = WLR.simple_roots()
            return WLR.sum(alpha[j] for a,b in self.value for j in range(a,b+1))

