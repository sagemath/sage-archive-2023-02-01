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
from root_lattice_realization import RootLatticeRealization, RootLatticeRealizationElement
from sage.misc.cachefunc import ClearCacheOnPickle

# TODO: inheriting from ClearCacheOnPickle is a technical detail unrelated to root spaces
# could we abstract this somewhere higher?

class RootSpace(ClearCacheOnPickle, CombinatorialFreeModule, RootLatticeRealization):
    r"""
    INPUT:
     - root_system: a root system
     - base_ring: a ring `R`

    The root space (or lattice if base_ring is ZZ) of a root system,
    that is the formal free module `\bigoplus_i R \alpha_i` generated
    by the simple roots `\alpha_i`.

    This class is also used for coroot spaces (or lattices).

    Todo: standardize the variable used for the root space in the examples (P?)

    TESTS::

        sage: for ct in CartanType.samples():
        ...       if ct.is_implemented():
        ...           P = ct.root_system().root_space()
        ...           TestSuite(P).run()
        ...
        sage: r = RootSystem(['A',4]).root_lattice()
        sage: r.simple_root(1)
        alpha[1]
        sage: latex(r.simple_root(1))
        \alpha_{1}

    """

    def __init__(self, root_system, base_ring):
        """
        EXAMPLES::

            sage: P = RootSystem(['A',4]).root_space()
            sage: s = P.simple_reflections()

        """
        self.root_system = root_system
        basis_name = "alphacheck" if root_system.dual_side else "alpha"
        basis_name_latex  = "\\alpha^\\vee" if root_system.dual_side else "\\alpha"
        CombinatorialFreeModule.__init__(self, base_ring,
                                         root_system.index_set(),
                                         element_class = RootSpaceElement,
                                         prefix=basis_name,
                                         latex_prefix=basis_name_latex)

    def _repr_(self):
        """
        TESTS::
            sage: RootSystem(['A',4]).root_lattice()
            Root lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).root_space()
            Root space over the Rational Field of the Root system of type ['B', 4]
            sage: RootSystem(['A',4]).coroot_lattice()
            Coroot lattice of the Root system of type ['A', 4]
            sage: RootSystem(['B',4]).coroot_space()
            Coroot space over the Rational Field of the Root system of type ['B', 4]

        """
        return self._name_string()

    def _name_string(self, capitalize=True, base_ring=True, type=True):
        """
        EXAMPLES::

            sage: RootSystem(['A',4]).root_space()._name_string()
            "Root space over the Rational Field of the Root system of type ['A', 4]"
        """
        return self._name_string_helper("root", capitalize=capitalize, base_ring=base_ring, type=type)

    simple_root = CombinatorialFreeModule.monomial


class RootSpaceElement(CombinatorialFreeModuleElement, RootLatticeRealizationElement):
    def scalar(self, lambdacheck):
        """
        The scalar product between the root lattice and
        the coroot lattice.

        EXAMPLES::

            sage: r = RootSystem(['A',4]).root_lattice()
            sage: cr = RootSystem(['A',4]).coroot_lattice()
            sage: a1 = r.simple_root(1)
            sage: ac1 = cr.simple_root(1)
            sage: a1.scalar(ac1)
            2

        """
        # Find some better test
        assert lambdacheck in self.parent().coroot_lattice() or lambdacheck in self.parent().coroot_space()
        zero = self.parent().base_ring().zero()
        cartan_matrix = self.parent().dynkin_diagram()
        return sum( (sum( (lambdacheck[i]*s for i,s in cartan_matrix.column(j)), zero) * c for j,c in self), zero)

    def is_positive_root(self):
        """
        Checks whether an element in the root space lies in the
        nonnegative cone spanned by the simple roots.

        EXAMPLES::

            sage: R=RootSystem(['A',3,1]).root_space()
            sage: B=R.basis()
            sage: w=B[0]+B[3]
            sage: w.is_positive_root()
            True
            sage: w=B[1]-B[2]
            sage: w.is_positive_root()
            False
        """
        for c in self.coefficients():
            if c < 0:
                return False
        return True

