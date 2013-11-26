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
    A virtual crystal `\widehat{B}` of an ambient crystal `B` is a crystal
    formed by taking a subset of `B` and whose crystal structure is given by

    .. MATH::

        \widehat{e}_i = \prod_{j \in \sigma_i} e_j^{\gamma_i}, \quad
        \widehat{f}_i = \prod_{j \in \sigma_i} f_j^{\gamma_i},

    .. MATH::

        \widehat{\varepsilon}_i = \frac{\varepsilon_j}{\gamma_j}, \quad
        \widehat{\varphi}_i = \frac{\varphi_j}{\gamma_j}, \quad

    where `\sigma_i` is a set of indices in `B` and the *scaling factors*
    `\gamma_i`. We note that for the crystal to be well-defined, we must have

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
    - ``cartan_type`` -- (optional) the Cartan type for the virtual crystal;
      the default is the Cartan type for the ambient crystal
    - ``generators`` -- (optional) the generators for the virtual crystal; the
      default is the generators for the ambient crystal
    - ``index_set`` -- (optional) the index set for the virtual crystal; the
      default is the index set for the Cartan type
    - ``category`` -- (optional) the category for the virtual crystal; the
      default is the :class:`~sage.categories.crystals.Crystals` category
    """
    def __init__(self, ambient, virtualization, scaling_factors,
                 cartan_type=None, generators=None, index_set=None, category=None):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._virtualization = Family(virtualization)
        self._scaling_factors = Family(scaling_factors)
        Subcrystal.__init__(ambient, None, generators, cartan_type, index_set, category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

        """
        return "Virtual crystal of {} of type {}".format(self._ambient, self._cartan_type)

    def __contains__(self, x):
        """
        Check if ``x`` is in ``self``.

        EXAMPLES::
        """
        if not Subcrystal.__contains__(self, x):
            return False
        # TODO: check that it is in the vitual crystal
        return False

    def virtualization(self):
        """
        Return the virtualization sets `\sigma_i`.

        EXAMPLES::
        """
        return self._virtualization

    def scaling_factors(self):
        """
        Return the scaling factors `\gamma_i`.

        EXAMPLES::
        """
        return self._scaling_factors

    class Element(ElementWrapper):
        def e(self, i):
            """
            Return `e_i` of ``self``.

            EXAMPLES::
            """
            s = []
            sf = self._scaling_factors[i]
            for j in self._virtualization[i]:
                s += [j]*sf
            ret = self.value.e_string(s)
            if ret is None:
                return None
            return self.__class__(self.parent(), ret)

        def f(self, i):
            """
            Return `f_i` of ``self``.

            EXAMPLES::
            """
            s = []
            sf = self._scaling_factors[i]
            for j in self._virtualization[i]:
                s += [j]*sf
            ret = self.value.f_string(s)
            if ret is None:
                return None
            return self.__class__(self.parent(), ret)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            EXAMPLES::
            """            
            return self.value.epsilon(self._virtualization[i][0]) / self._scaling_factors[i]

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            EXAMPLES::
            """
            return self.value.phi(self._virtualization[i][0]) / self._scaling_factors[i]

        def weight(self):
            """
            Return the weight of ``self``.

            EXAMPLES::
            """
            WLR = self.weight_lattice_realization()
            wt = self.value.weight()
            La = WLR.fundamental_weights()
            sf = self.parent().scaling_factors()
            return WLR.sum_of_terms((i, wt / sf[i]) for i in self.index_set())

