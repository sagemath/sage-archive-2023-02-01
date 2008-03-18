"""
Spaces of homomorphisms between modular abelian varieties.

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.homset import HomsetWithBase

import abvar as abelian_variety
import morphism

import sage.rings.integer_ring

from sage.rings.ring import Ring
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.element import is_Matrix

ZZ = sage.rings.integer_ring.ZZ

class Homspace(HomsetWithBase):
    """
    A space of homomorphisms between two modular abelian varieties.
    """
    def __init__(self, domain, codomain, cat):
        """
        Create a homspace.

        INPUT:
            domain, codomain -- modular abelian varieties
            cat -- category

        EXAMPLES:
            sage: H = Hom(J0(11), J0(22)); H
            Space of homomorphisms from Jacobian of the modular curve associated to the congruence subgroup Gamma0(11) to Jacobian of the modular curve associated to the congruence subgroup Gamma0(22)
            sage: type(H)
            <class 'sage.modular.abvar.homspace.Homspace'>
            sage: H.homset_category()
            Category of modular abelian varieties over Rational Field
        """
        if not abelian_variety.is_ModularAbelianVariety(domain):
            raise TypeError, "domain must be a modular abelian variety"
        if not abelian_variety.is_ModularAbelianVariety(codomain):
            raise TypeError, "codomain must be a modular abelian variety"
        HomsetWithBase.__init__(self, domain, codomain, cat)

    def _repr_(self):
        """
        String representation of a modular abelian variety homspace.

        EXAMPLES:
            sage: J = J0(11)
            sage: End(J)._repr_()
            'Space of homomorphisms from Jacobian of the modular curve associated to the congruence subgroup Gamma0(11) to Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)'
        """
        return "Space of homomorphisms from %s to %s"%\
               (self.domain(), self.codomain())


class EndomorphismSubring(Homspace, Ring):

    def __init__(self, A):
        """
        """
        self._J = A.ambient_variety()
        self._A = A
        self._E = A.Hom(A)
        self._gens_set = False
        self._matrix_space = MatrixSpace(ZZ,2*A.dimension())

    def _repr_(self):
        return "Subring of endomorphism ring of %s"%self._A

    def domain(self):
        return self._A

    def codomain(self):
        return self._A

    def abelian_variety(self):
        return self._A

    def _End(self):
        return self._E

    def matrix_space(self):
        return self._matrix_space

    def _set_generators(self, gens):
        self._gens_set = True
        self._gens = tuple([ self.matrix_space()(g.list()) for g in gens ])

    def index_in(self, other, check=True):
        if check:
            if not isinstance(other, EndomorphismSubring):
                raise ValueError, "other must be a subring of an endomorphism ring of an abelian variety."
            if not (self.abelian_variety() == other.abelian_variety()):
                raise ValueError, "self and other must be endomorphisms of the same abelian variety"

        M = self.free_module()
        N = other.free_module()
        return M.index_in(N)


    def free_module(self):
        if not self._gens_set:
            raise ValueError, "generators of self unknown"
        V = ZZ**(4*self.abelian_variety().dimension())
        return V.submodule([ V(m.list()) for m in self.gens() ])

    def gen(self, i=0):
        if not self._gens_set:
            raise ValueError, "generators of self unknown"
        if i > self.ngens():
            raise ValueError, "self only has %s generators"%self.ngens()
        return self._gens[i]

    def ngens(self):
        if not self._gens_set:
            raise ValueError, "number of generators unknown"
        return len(self._gens)

    def __call__(self, M, check=True):
        if check:
            if not is_Matrix(M):
                raise ValueError, "can only coerce in matrices"
            if (M.nrows() != M.ncols()) or (M.nrows() != 2*self.abelian_variety().dimension()):
                raise ValueError, "incorrect matrix dimensions"
            if M.base_ring() != ZZ:
                M = M.change_ring(ZZ)

        return morphism.Morphism(self._End(), M)

    def image_of_hecke_algebra(self):
        A = self.abelian_variety()
        #if A.is_hecke_stable():
        return

class EndomorphismSubAlgebra(EndomorphismSubring):
    def __init__(self, A):
        EndomorphismSubring.__init__(self, A)


