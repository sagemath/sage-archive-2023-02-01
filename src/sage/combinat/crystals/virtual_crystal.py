r"""
Virtual Crystals

These are the crystals that are subsets of a larger ambient crystal with
virtual crystal operators.

AUTHORS:

- Travis Scrimshaw (2013-10-16): Initial implementation
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

from sage.misc.cachefunc import cached_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.categories.crystals import Crystals
from sage.categories.sets_cat import Sets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.subcrystal import Subcrystal

class VirtualCrystal(Subcrystal):
    r"""
    A virtual crystal `B` of an ambient crystal `\widehat{B}` is a crystal
    formed by taking a subset of `\widehat{B}` and whose crystal structure is
    given by

    .. MATH::

        \widehat{e}_i = \prod_{j \in \sigma_i} e_j, \quad
        \widehat{f}_i = \prod_{j \in \sigma_i} f_j,

    The below is wrong!!!

    .. MATH::

        \widehat{\varepsilon}_i = \sum_{j \in \sigma_i} \varepsilon_j, \quad
        \widehat{\varphi}_i = \sum_{j \in \sigma_i} \varphi_j,

    where `\sigma_i` is a sequence of indices in `\widehat{B}`.

    INPUT:

    - ``ambient`` -- the ambient crystal
    - ``virtualization`` -- a dictionary whose key `i` corresponds
      to `\sigma_i`
    - ``generators`` -- (optional) the generators for the virtual crystal; the
      default is the generators for the ambient crystal
    - ``cartan_type`` -- (optional) the Cartan type for the virtual crystal;
      the default is the Cartan type for the ambient crystal
    - ``index_set`` -- (optional) the index set for the virtual crystal; the
      default is the index set for the Cartan type
    - ``category`` -- (optional) the category for the virtual crystal; the
      default is the :class:`~sage.categories.crystals.Crystals` category
    """
    def __init__(self, ambient, virtualization, cartan_type=None,
                 generators=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

        """
        self._virtualization = virtualization
        Subcrystal.__init__(ambient, None, generators, cartan_type, index_set, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

        """
        return "Virtual crystal of {}".format(self._ambient)

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
        """
        if not Subcrystal.__contains__(self, x):
            return False
        # TODO: check that it is in the vitual crystal
        return False

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: S.index_set()
            (1, 2)
        """
        return self._index_set

    class Element(ElementWrapper):
        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::
            """
            ret = self.value.e_string(virtualization[i])
            if ret is None:
                return None
            return self.__class__(self.parent(), ret)

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::
            """
            ret = self.value.f_string(virtualization[i])
            if ret is None:
                return None
            return self.__class__(self.parent(), ret)

        # Phi and Epsilon are wrong
        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::
            """
            return sum(self.value.epsilon(j) for j in virtualization[i])

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::
            """
            return sum(self.value.phi(j) for j in virtualization[i])

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.weight()
                (0, 1, 0, 1, 1)
            """
            return self.value.weight()

