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
from sage.categories.finite_crystals import FiniteCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.crystals.subcrystal import Subcrystal
from sage.sets.family import Family

class VirtualCrystal(Subcrystal):
    r"""
    A virtual crystal `\widehat{B}` of an ambient crystal `B` is a crystal
    formed by taking a subset of `B` and whose crystal structure is given by

    .. MATH::

        \widehat{e}_i = \prod_{j \in \sigma_i} e_j^{\gamma_i}, \quad
        \widehat{f}_i = \prod_{j \in \sigma_i} f_j^{\gamma_i},

    .. MATH::

        \widehat{\varepsilon}_i = \frac{\varepsilon_j}{\gamma_j}, \quad
        \widehat{\varphi}_i = \frac{\varphi_j}{\gamma_j}, \quad

    where `\sigma_i` is a set of indices in `B` and the *scaling factors*
    `\gamma_i \in \ZZ`. We note that for the crystal to be well-defined,
    we must have

    .. MATH::

        \varepsilon_j = \varepsilon_{j^{\prime}},
        \quad \varphi_j = \varphi_{j^{\prime}}

    for all `j, j^{\prime} \in \sigma_i` and that the order that the Kashiwara
    operators in the ambient space are applied does not affect the result.
    We define the weight `\widehat{\mathrm{wt}}` to be such that

    .. MATH::

        \widehat{\mathrm{wt}}(b) = c_i \widehat{\Lambda}_i =
        c_i \sum_{j \in \sigma_i} \gamma_j \Lambda_j = \mathrm{wt}(b).

    INPUT:

    - ``ambient`` -- the ambient crystal
    - ``virtualization`` -- a dictionary whose key `i` corresponds
      to the set `\sigma_i`
    - ``scaling_factors`` -- a dictionary whose key `i` corresponds to
      the scaling factor `\gamma_i`
    - ``contained`` -- (optional) a set (or function) which specifies when an
      element is contained in the subcrystal; the default is everything
      possible is included
    - ``generators`` -- (optional) the generators for the virtual crystal; the
      default is the generators for the ambient crystal
    - ``cartan_type`` -- (optional) the Cartan type for the virtual crystal;
      the default is the Cartan type for the ambient crystal
    - ``index_set`` -- (optional) the index set for the virtual crystal; the
      default is the index set for the Cartan type
    - ``category`` -- (optional) the category for the virtual crystal; the
      default is the :class:`~sage.categories.crystals.Crystals` category

    REFERENCES:

    - [FOS09]_
    - [OSS03]_
    - [OSS2003]_
    """
    @staticmethod
    def __classcall_private__(cls, ambient, virtualization, scaling_factors,
                              contained=None, generators=None,
                              cartan_type=None, index_set=None, category=None):
        """
        Normalize arguments to ensure a unique representation.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi1 = B.crystal_morphism(C.module_generators)
            sage: V1 = psi1.image()
            sage: psi2 = B.crystal_morphism(C.module_generators, index_set=[1,2,3])
            sage: V2 = psi2.image()
            sage: V1 is V2
            True
        """
        if cartan_type is None:
            cartan_type = ambient.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
        if index_set is None:
            index_set = cartan_type.index_set
        if generators is None:
            generators = ambient.module_generators
        virtualization = Family(virtualization)
        scaling_factors = Family(scaling_factors)

        category = Crystals().or_subcategory(category)

        return super(Subcrystal, cls).__classcall__(cls, ambient, virtualization, scaling_factors,
                                                    contained, tuple(generators), cartan_type,
                                                    tuple(index_set), category)

    def __init__(self, ambient, virtualization, scaling_factors,
                 contained, generators, cartan_type, index_set, category):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: V = psi.image()
            sage: TestSuite(V).run()
        """
        self._virtualization = virtualization
        self._scaling_factors = scaling_factors
        Subcrystal.__init__(self, ambient, contained, generators,
                            cartan_type, index_set, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: psi.image()
            Virtual crystal of The crystal of tableaux of type ['D', 4] and shape(s) [[2]] of type ['B', 3]
        """
        return "Virtual crystal of {} of type {}".format(self._ambient, self._cartan_type)

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: V = psi.image()
            sage: mg = C.module_generators[0]
            sage: mg in V
            True
            sage: mg.f(1) in V
            False
            sage: mg.f(1).f(1) in V
            True
        """
        if not Subcrystal.__contains__(self, x):
            return False
        if self in FiniteCrystals():
            if isinstance(x, self._ambient.element_class):
                if x.parent() == self:
                    x = self.element_class(self, self._ambient(x))
                elif x.parent() == self._ambient:
                    x = self.element_class(self, self._ambient(x))
            elif isinstance(x, self.element_class) and x.parent() != self:
                x = self.element_class(self, x.value)
            return x in self.list()
        return True

    def virtualization(self):
        """
        Return the virtualization sets `\sigma_i`.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: V = psi.image()
            sage: V.virtualization()
            Finite family {1: (1,), 2: (2,), 3: (3, 4)}
        """
        return self._virtualization

    def scaling_factors(self):
        """
        Return the scaling factors `\gamma_i`.

        EXAMPLES::

            sage: B = CrystalOfTableaux(['B',3], shape=[1])
            sage: C = CrystalOfTableaux(['D',4], shape=[2])
            sage: psi = B.crystal_morphism(C.module_generators)
            sage: V = psi.image()
            sage: V.scaling_factors()
            Finite family {1: 2, 2: 2, 3: 1}
        """
        return self._scaling_factors

    class Element(Subcrystal.Element):
        """
        An element of a virtual (sub)crystal. Wraps an element in the
        ambient crystal.
        """
        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['B',3], shape=[1])
                sage: C = CrystalOfTableaux(['D',4], shape=[2])
                sage: psi = B.crystal_morphism(C.module_generators)
                sage: V = psi.image()
                sage: mg = V.module_generators[0]
                sage: mg.e(1)
                sage: b = psi(B.module_generators[0].f(1))
                sage: V(b).e(1)
                [[1, 1]]
            """
            s = []
            P = self.parent()
            sf = P._scaling_factors[i]
            for j in P._virtualization[i]:
                s += [j]*sf
            ret = self.value.e_string(s)
            if ret is None:
                return None
            return self.__class__(P, ret)

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['B',3], shape=[1])
                sage: C = CrystalOfTableaux(['D',4], shape=[2])
                sage: psi = B.crystal_morphism(C.module_generators)
                sage: V = psi.image()
                sage: mg = V.module_generators[0]
                sage: mg.f(1)
                [[2, 2]]
                sage: mg.f(2)
            """
            s = []
            P = self.parent()
            sf = P._scaling_factors[i]
            for j in P._virtualization[i]:
                s += [j]*sf
            ret = self.value.f_string(s)
            if ret is None:
                return None
            return self.__class__(P, ret)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['B',3], shape=[1])
                sage: C = CrystalOfTableaux(['D',4], shape=[2])
                sage: psi = B.crystal_morphism(C.module_generators)
                sage: V = psi.image()
                sage: mg = V.module_generators[0]
                sage: mg.epsilon(2)
                0
                sage: mg.f(1).epsilon(1)
                1
            """
            P = self.parent()
            return self.value.epsilon(P._virtualization[i][0]) / P._scaling_factors[i]

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['B',3], shape=[1])
                sage: C = CrystalOfTableaux(['D',4], shape=[2])
                sage: psi = B.crystal_morphism(C.module_generators)
                sage: V = psi.image()
                sage: mg = V.module_generators[0]
                sage: mg.phi(1)
                1
                sage: mg.phi(2)
                0
            """
            P = self.parent()
            return self.value.phi(P._virtualization[i][0]) / P._scaling_factors[i]

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::

                sage: B = CrystalOfTableaux(['B',3], shape=[1])
                sage: C = CrystalOfTableaux(['D',4], shape=[2])
                sage: psi = B.crystal_morphism(C.module_generators)
                sage: V = psi.image()
                sage: mg = V.module_generators[0]
                sage: mg.weight()
                (1, 0, 0)
                sage: mg.f(1).weight()
                (0, 1, 0)
            """
            P = self.parent()
            WLR = P.weight_lattice_realization()
            wt = self.value.weight()
            ac = P._ambient.weight_lattice_realization().simple_coroots()
            La = WLR.fundamental_weights()
            v = P._virtualization
            sf = P._scaling_factors
            return WLR.sum(wt.scalar(ac[v[i][0]]) / sf[i] * La[i] for i in self.index_set())

# TODO: implement a devirtualization map

