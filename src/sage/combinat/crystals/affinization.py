r"""
Affinization Crystals
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
from sage.structure.element import Element
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.all import ZZ
from sage.combinat.partition_tuple import PartitionTuple
from sage.combinat.root_system.cartan_type import CartanType

class AffinizationCrystal(Parent, UniqueRepresentation):
    r"""
    An affiniziation crystal.

    Let `\mathfrak{g}` be of affine type. An affinization crystal is a
    `U_q(\mathfrak{g})`-crystal formed from a
    `U_q^{\prime}(\mathfrak{g})`-crystal `B` by taking the set:

    .. MATH::

        B^{aff} = \{ b(m) \mid b \in B, m \in \ZZ \}

    and defining the crystal structure by:

    .. MATH::

        \begin{aligned}
        e_i(b(m)) & = \begin{cases} (e_0 b)(m+1) & i = 0
        \\ (e_i b)(m) & i \neq 0 \end{cases} \\

        \\ f_i(b(m)) & = \begin{cases} (f_0 b)(m-1) & i = 0
        \\ (f_i b)(m) & i \neq 0 \end{cases}

        \\ \mathrm{wt}(b(m)) & = \wt(b) + m \delta.
        \end{aligned}
    """
    def __init__(self, B):
        """
        Initialize ``self``.
        """
        self._B = B
        Parent.__init__(self, category=(RegularCrystals(), InfiniteEnumeratedSets()))
        self.moodule_generators = (self.element_class(self, self._B.module_generator(), 0),)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Affinization of {}".format(self._B)

    def weight_lattice_realization(self):
        """
        Return the weight lattice realization of ``self``.
        """
        return self.cartan_type().ambient_space()

    class Element(Element):
        """
        An element in an affinization crystal.
        """
        def __init__(self, parent b, m):
            """
            Initialize ``self``.
            """
            self._b = b
            self._m = m
            Element.__init__(self, parent)

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.
            """
            if i == 0:
                return self.__class__(self.parent(), self._b.e(0), self._m+1)
            return self.__class__(self.parent(), self._b.e(i), self._m)

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.
            """
            if i == 0:
                return self.__class__(self.parent(), self._b.f(0), self._m+1)
            return self.__class__(self.parent(), self._b.f(i), self._m)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.
            """
            return self._b.epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.
            """
            return self._b.phi(i)

        def weight(self):
            r"""
            Return the weight of ``self``.
            """
            WLR = self.parent().weight_lattice_realization()
            La = WLR.fundamental_weights()
            return WLR.sum(c*La[i] for i,c in self._b.weight()) + self._m * WLR.null_root()

