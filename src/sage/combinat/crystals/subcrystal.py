r"""
Subcrystals

These are the crystals that are subsets of a larger ambient crystal.

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

class Subcrystal(Parent):
    """
    A subcrystal `X` of an ambient crystal `Y` is a crystal formed by taking a
    subset of `Y` and whose crystal structure is induced by `Y`.

    INPUT:

    - ``ambient`` -- the ambient crystal
    - ``contained`` -- (optional) a function or set which specifies when an
      element is contained in the subcrystal; the default is everything is
      the function ``lambda x: True``
    - ``generators`` -- (optional) the generators for the subcrystal; the
      default is the generators for the ambient crystal
    - ``cartan_type`` -- (optional) the Cartan type for the subcrystal; the
      default is the Cartan type for the ambient crystal
    - ``index_set`` -- (optional) the index set for the subcrystal; the
      default is the index set for the Cartan type
    - ``category`` -- (optional) the category for the subcrystal; the
      default is the :class:`~sage.categories.crystals.Crystals` category
    """
    def __init__(self, ambient, contained=None, generators=None,
                 cartan_type=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            sage: TestSuite(S).run()
        """
        if contained is None:
            contained = lambda x: True
        elif isinstance(contained, (list, tuple, set, frozenset)) \
                or contained in Sets():
            contained = lambda x: x in contained

        if cartan_type is None:
            cartan_type = ambient.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
        if index_set is None:
            index_set = cartan_type.index_set
        if generators is None:
            generators = ambient.module_generators

        category = Crystals().or_subcategory(category)

        self._ambient = ambient
        self._contained = contained
        self._cartan_type = cartan_type
        self._index_set = tuple(index_set)
        Parent.__init__(self, category=category)
        self.module_generators = tuple(map(self, generators))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
            sage: B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
            Subcrystal of The crystal of tableaux of type ['A', 4] and shape(s) [[2, 1]]
        """
        return "Subcrystal of {}".format(self._ambient)

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
            sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
        """
        if isinstance(x, self.element_class) and x.parent() is self:
            return True
        if x in self._ambient:
            return self._contained(x)
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

                sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.e(2)
                sage: mg.e(1)
                [[1, 4], [5]]
            """
            ret = self.value.e(i)
            if ret is None or not self.parent()._contained(ret):
                return None
            return self.__class__(self.parent(), ret)

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.f(1)
                sage: mg.f(2)
                [[3, 4], [5]]
            """
            ret = self.value.f(i)
            if ret is None or not self.parent()._contained(ret):
                return None
            return self.__class__(self.parent(), ret)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.epsilon(1)
                1
                sage: mg.epsilon(2)
                0
            """
            return self.value.epsilon(i)

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['A',4], shape=[2,1])
                sage: S = B.subcrystal(generators=(B(2,1,1), B(5,2,4)), index_set=[1,2])
                sage: mg = S.module_generators[1]
                sage: mg.phi(1)
                0
                sage: mg.phi(2)
                1
            """
            return self.value.phi(i)

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

