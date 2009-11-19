"""
Root lattices and root spaces
"""
#*****************************************************************************
#       Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from root_lattice_realization import RootLatticeRealizationElement
from weight_lattice_realization import WeightLatticeRealization


class WeightSpace(CombinatorialFreeModule, WeightLatticeRealization):
    r"""
    INPUT:
     - root_system: a root system
     - base_ring: a ring `R`

    The weight space (or lattice if base_ring is ZZ) of a root system,
    that is the formal free module `\bigoplus_i R \Lambda_i` generated
    by the simple fundamental weights `\Lambda_i`.

    This class is also used for coweight spaces (or lattices).

    TESTS::
        sage: for ct in CartanType.samples():
        ...       if ct.is_implemented():
        ...           P = ct.root_system().weight_space()
        ...           TestSuite(P).run()
        ...
    """

    def __init__(self, root_system, base_ring):
        """
        EXAMPLES::

            sage: Q = RootSystem(['A',4]).weight_lattice()
            sage: s = Q.simple_reflections()

        """
        self.root_system = root_system
        basis_name = "Lambdacheck" if root_system.dual_side else "Lambda"
        CombinatorialFreeModule.__init__(self, base_ring,
                                         root_system.index_set(),
                                         element_class = WeightSpaceElement,
                                         prefix=basis_name)

    def _repr_(self):
        """
        TESTS::

            sage: RootSystem(['A',4]).weight_lattice()
            Weight lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).weight_space()
            Weight space over the Rational Field of the Root system of type ['B', 4]
            sage: RootSystem(['A',4]).coweight_lattice()
            Coweight lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).coweight_space()
            Coweight space over the Rational Field of the Root system of type ['B', 4]

        """
        return self._name_string()

    def _name_string(self, capitalize=True, base_ring=True, type=True):
        """
        EXAMPLES::

            sage: RootSystem(['A',4]).weight_lattice()._name_string()
            "Weight lattice of the Root system of type ['A', 4]"
        """
        return self._name_string_helper("weight", capitalize=capitalize, base_ring=base_ring, type=type)

    fundamental_weight = CombinatorialFreeModule.term

    def simple_root(self, j):
        """
        Returns the `j^{th}` simple root

        TESTS::

            sage: R = RootSystem(["C",4])
            sage: R.weight_lattice().simple_root(3)
            -Lambda[2] + 2*Lambda[3] - Lambda[4]

            sage: R.weight_lattice().simple_roots()
            Finite family {1: 2*Lambda[1] - Lambda[2], 2: -Lambda[1] + 2*Lambda[2] - Lambda[3], 3: -Lambda[2] + 2*Lambda[3] - Lambda[4], 4: -2*Lambda[3] + 2*Lambda[4]}

        """
        assert(j in self.index_set())
        return self._from_dict(dict([(i,c) for (i,c) in self.root_system.dynkin_diagram().column(j)]))


class WeightSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):
    def scalar(self, lambdacheck):
        """
        The canonical scalar product between the weight lattice and
        the coroot lattice.

        TODO: merge with_apply_multi_module_morphism

        EXAMPLES::

            sage: R = RootSystem(["C",4])
            sage: w = R.weight_lattice().simple_root(1)
            sage: cw = R.coweight_lattice().simple_root(1)
            sage: w.scalar(cw)
            5
        """
        zero = self.parent().base_ring().zero()
        if len(self) < len(lambdacheck):
            return sum( (lambdacheck[i]*c for (i,c) in self), zero)
        else:
            return sum( (self[i]*c for (i,c) in lambdacheck), zero)

    def is_dominant(self):
	"""
	Checks whether an element in the weight space lies in the positive cone spanned
	by the basis elements (fundamental weights).

	EXAMPLES::

	    sage: W=RootSystem(['A',3]).weight_space()
	    sage: Lambda=W.basis()
	    sage: w=Lambda[1]+Lambda[3]
	    sage: w.is_dominant()
	    True
	    sage: w=Lambda[1]-Lambda[2]
	    sage: w.is_dominant()
	    False
	"""
	for c in self.coefficients():
	    if c < 0:
		return False
	return True
