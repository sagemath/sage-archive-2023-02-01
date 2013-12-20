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

    EXAMPLES::

        sage: K = KirillovReshetikhinCrystal(['A',2,1], 1, 1)
        sage: A = K.affinization()

    REFERENCES:

    - [HK02]_
    """
    def __init__(self, B):
        """
        Initialize ``self``.
        """
        self._B = B
        self._cartan_type = B.cartan_type()
        Parent.__init__(self, category=(RegularCrystals(), InfiniteEnumeratedSets()))
        self.module_generators = (self.element_class(self, self._B.module_generator(), 0),)

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
        def __init__(self, parent, b, m):
            """
            Initialize ``self``.
            """
            self._b = b
            self._m = m
            Element.__init__(self, parent)

        def _repr_(self):
            """
            Return a string representation of ``self``.
            """
            return "{}({})".format(self._b, self._m)

        def _latex_(self):
            r"""
            Return a LaTeX representation of ``self``.
            """
            from sage.misc.latex import latex
            return latex(self._b) + "({})".format(self._m)

        def e(self, i):
            """
            Return the action of `e_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index se
            """
            bp = self._b.e(i)
            if bp is None:
                return None
            if i == 0:
                return self.__class__(self.parent(), bp, self._m+1)
            return self.__class__(self.parent(), bp, self._m)

        def f(self, i):
            """
            Return the action of `f_i` on ``self``.

            INPUT:

            - ``i`` -- an element of the index set
            """
            bp = self._b.f(i)
            if bp is None:
                return None
            if i == 0:
                return self.__class__(self.parent(), bp, self._m+1)
            return self.__class__(self.parent(), bp, self._m)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set
            """
            return self._b.epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            INPUT:

            - ``i`` -- an element of the index set
            """
            return self._b.phi(i)

        def weight(self):
            r"""
            Return the weight of ``self``.
            """
            WLR = self.parent().weight_lattice_realization()
            La = WLR.fundamental_weights()
            return WLR.sum(c*La[i] for i,c in self._b.weight()) + self._m * WLR.null_root()

