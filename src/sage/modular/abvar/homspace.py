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
from sage.misc.functional import parent

import abvar as abelian_variety
import morphism

import sage.rings.integer_ring

from sage.rings.ring import Ring
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix
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
            Space of homomorphisms from Abelian variety J0(11) of dimension 1 to Abelian variety J0(22) of dimension 2
            sage: Hom(J0(11), J0(11))
            Endomorphism ring of Abelian variety J0(11) of dimension 1
            sage: type(H)
            <class 'sage.modular.abvar.homspace.Homspace'>
            sage: H.homset_category()
            Category of modular abelian varieties over Rational Field
        """
        if not abelian_variety.is_ModularAbelianVariety(domain):
            raise TypeError, "domain must be a modular abelian variety"
        if not abelian_variety.is_ModularAbelianVariety(codomain):
            raise TypeError, "codomain must be a modular abelian variety"
        self._matrix_space = MatrixSpace(ZZ,2*domain.dimension(), 2*codomain.dimension())
        self._gens = None
        HomsetWithBase.__init__(self, domain, codomain, cat)

    def __call__(self, M):
        r"""
        Create a homomorphism in this space from M. M can
        be any of the following:

          * a Morphism of abelian varieties
          * a matrix of the appropriate size (i.e.
            2*self.domain().dimension() x
            2*self.codomain().dimension()) whose entries
            are coercible into self.base_ring()
          * anything that can be coerced into self.matrix_space()

        EXAMPLES:
            sage: H = Hom(J0(11), J0(22))
            sage: phi = H(matrix(ZZ,2,4,[5..12])) ; phi
            Abelian variety morphism:
              From: Abelian variety J0(11) of dimension 1
              To:   Abelian variety J0(22) of dimension 2
            sage: phi.matrix()
            [ 5  6  7  8]
            [ 9 10 11 12]
            sage: phi.matrix().parent()
            Full MatrixSpace of 2 by 4 dense matrices over Integer Ring

            sage: H = J0(22).Hom(J0(11)*J0(11))
            sage: m1 = J0(22).degeneracy_map(11,1).matrix() ; m1
            [ 0  1]
            [-1  1]
            [-1  0]
            [ 0 -1]
            sage: m2 = J0(22).degeneracy_map(11,2).matrix() ; m2
            [ 1 -2]
            [ 0 -2]
            [ 1 -1]
            [ 0 -1]
            sage: m = m1.transpose().stack(m2.transpose()).transpose() ; m
            [ 0  1  1 -2]
            [-1  1  0 -2]
            [-1  0  1 -1]
            [ 0 -1  0 -1]
            sage: phi = H(m) ; phi
            Abelian variety morphism:
              From: Abelian variety J0(22) of dimension 2
              To:   Abelian variety J0(11) x J0(11) of dimension 2
            sage: phi.matrix()
            [ 0  1  1 -2]
            [-1  1  0 -2]
            [-1  0  1 -1]
            [ 0 -1  0 -1]
        """
        if isinstance(M, morphism.Morphism):
            if M.parent() is self:
                return M
            elif M.domain() == self.domain() and M.codomain() == self.codomain():
                M = M.matrix()
            else:
                raise ValueError, "cannot convert %s into %s" % (M, self)
        elif is_Matrix(M):
            if M.base_ring() != ZZ:
                M = M.change_ring(ZZ)
            if M.nrows() != 2*self.domain().dimension() or M.ncols() != 2*self.codomain().dimension():
                raise TypeError, "matrix has wrong dimension"
        elif self.matrix_space().has_coerce_map_from(parent(M)):
            M = self.matrix_space()(M)
        else:
            raise TypeError, "can only coerce in matrices or morphisms"
        return morphism.Morphism(self, M)

    def _coerce_impl(self, x):
        """
        Coerce x into self, if possible.

        EXAMPLES:
            sage: J = J0(37) ; J.Hom(J)._coerce_impl(matrix(ZZ,4,[5..20]))
            Abelian variety endomorphism of Abelian variety J0(37) of dimension 2
            sage: K = J0(11) * J0(11) ; J.Hom(K)._coerce_impl(matrix(ZZ,4,[5..20]))
            Abelian variety morphism:
              From: Abelian variety J0(37) of dimension 2
              To:   Abelian variety J0(11) x J0(11) of dimension 2
        """
        if self.matrix_space().has_coerce_map_from(parent(x)):
            return self(x)
        else:
            return HomsetWithBase._coerce_impl(self, x)

    def _repr_(self):
        """
        String representation of a modular abelian variety homspace.

        EXAMPLES:
            sage: J = J0(11)
            sage: End(J)._repr_()
            'Space of homomorphisms from Abelian variety J0(11) of dimension 1 to Abelian variety J0(11) of dimension 1'
        """
        return "Space of homomorphisms from %s to %s"%\
               (self.domain(), self.codomain())

    def _get_matrix(self, g):
        """
        Given an object g, try to return a matrix corresponding to g
        with dimensions the same as those of self.matrix_space().

        INPUT:
            g -- a matrix or morphism or object with a list method

        OUTPUT:
            a matrix

        EXAMPLES:
            sage: E = End(J0(11))
            sage: E._get_matrix(matrix(QQ,2,[1,2,3,4]))
            [1 2]
            [3 4]
            sage: E._get_matrix(J0(11).hecke_operator(2))
            [-2  0]
            [ 0 -2]

            sage: H = Hom(J0(11) * J0(17), J0(22))

            sage: H._get_matrix(tuple([8..23]))

            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]

            sage: H._get_matrix(tuple([8..23]))

            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]

            sage: H._get_matrix([8..23])

            [ 8  9 10 11]
            [12 13 14 15]
            [16 17 18 19]
            [20 21 22 23]

        """
        try:
            if g.parent() is self.matrix_space():
                return g
        except AttributeError:
            pass

        if isinstance(g, morphism.Morphism):
            return g.matrix()
        elif hasattr(g, 'list'):
            return self.matrix_space()(g.list())
        else:
            return self.matrix_space()(g)

    def free_module(self):
        r"""
        Return this endomorphism ring as a free submodule of a big $\ZZ^{4nm}$,
        where $n$ is the dimension of the domain abelian variety and $m$ the
        dimension of the codomain.

        OUTPUT:
            free module

        EXAMPLES:
            sage: E = Hom(J0(11), J0(22))
            sage: E.free_module()
            Free module of degree 8 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0 0 0 1 1 0 0]
            [0 0 0 1 0 0 1 1]
        """
        self.calculate_generators()
        V = ZZ**(4*self.domain().dimension() * self.codomain().dimension())
        return V.submodule([ V(m.matrix().list()) for m in self.gens() ])

    def gen(self, i=0):
        """
        Return i-th generator of self.

        INPUT:
            i -- an integer

        OUTPUT:
            a morphism

        EXAMPLES:
            sage: E = End(J0(22))
            sage: E.gen(0).matrix()
            [3 0 0 0]
            [0 3 0 0]
            [0 0 0 0]
            [0 0 0 0]
        """
        self.calculate_generators()
        if i > self.ngens():
            raise ValueError, "self only has %s generators"%self.ngens()
        return morphism.Morphism(self, self._gens[i])

    def ngens(self):
        """
        Return number of generators of self.

        OUTPUT:
            integer

        EXAMPLES:
            sage: E = End(J0(22))
            sage: E.ngens()
            4
        """
        self.calculate_generators()
        return len(self._gens)

    def gens(self):
        """
        Return tuple of generators for this endomorphism ring.

        EXAMPLES:
            sage: E = End(J0(22))
            sage: E.gens()
            (Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2,
             Abelian variety endomorphism of Abelian variety J0(22) of dimension 2)
        """
        try:
            return self._gen_morphisms
        except AttributeError:
            self.calculate_generators()
            self._gen_morphisms = tuple([self.gen(i) for i in range(self.ngens())])
            return self._gen_morphisms

    def matrix_space(self):
        """
        Return the underlying matrix space that we view this endomorphism ring as
        being embedded into.

        EXAMPLES:
            sage: E = End(J0(22))
            sage: E.matrix_space()
            Full MatrixSpace of 4 by 4 dense matrices over Integer Ring
        """
        return self._matrix_space

    def calculate_generators(self):
        """
        If generators haven't already been computed, calculate generators
        for this homspace.   If they have been computed, do nothing.

        EXAMPLES:
            sage: E = End(J0(11))
            sage: E.calculate_generators()
        """

        if self._gens is not None:
            return

        phi = self.domain()._isogeny_to_product_of_simples()
        psi = self.codomain()._isogeny_to_product_of_simples()

        H_simple = phi.codomain().Hom(psi.codomain())
        im_gens = H_simple._calculate_product_gens()

        M = phi.matrix()
        Mt = psi.complementary_isogeny().matrix()

        R = ZZ**(4*self.domain().dimension()*self.codomain().dimension())
        gens = R.submodule([ (M*self._get_matrix(g)*Mt).list() for g in im_gens ]).saturation().basis()
        self._gens = tuple([ self._get_matrix(g) for g in gens ])

    def _calculate_product_gens(self):
        """
        If generators haven't already been computed, calculate generators
        for this homspace.   If they have been computed, do nothing.

        EXAMPLES:
            sage: E = End(J0(11))
            sage: E.calculate_generators()
        """

        Afactors = self.domain().decomposition(simple=False)
        Bfactors = self.codomain().decomposition(simple=False)
        matrix_space = self.matrix_space()
        if len(Afactors) == 1 and len(Bfactors) == 1:
            Asimples = Afactors[0].decomposition()
            Bsimples = Bfactors[0].decomposition()
            if len(Asimples) == 1 and len(Bsimples) == 1:
                # Handle the base case of A, B simple
                gens = self._calculate_simple_gens()

            else:
                # Handle the case of A, B simple powers
                gens = []
                for i in range(len(Asimples)):
                    for j in range(len(Bsimples)):
                        hom_gens = Asimples[i].Hom(Bsimples[j]).gens()
                        for sub_gen in hom_gens:
                            sub_mat = sub_gen.matrix()
                            M = self.matrix_space()(0)
                            M.set_block(sub_mat.nrows()*i, sub_mat.ncols()*j, sub_mat)
                            gens.append(M)

        else:
            # Handle the case of A, B generic
            gens = []
            cur_row = 0
            for Afactor in Afactors:
                cur_row += Afactor.dimension() * 2
                cur_col = 0
                for Bfactor in Bfactors:
                    cur_col += Bfactor.dimension() * 2
                    Asimple = Afactor[0]
                    Bsimple = Bfactor[0]
                    if Asimple.newform_label() == Bsimple.newform_label():
                        for sub_gen in Afactor.Hom(Bfactor).gens():
                            sub_mat = sub_gen.matrix()
                            M = self.matrix_space()(0)
                            M.set_block(cur_row - sub_mat.nrows(),
                                        cur_col - sub_mat.ncols(),
                                        sub_mat)
                            gens.append(M)

        return gens

        # set the gens
        #R = ZZ**(4*self.domain().dimension()*self.codomain().dimension())
        #gens = R.submodule([ self._get_matrix(g).list() for g in gens ]).saturation().basis()
        #self._gens = tuple([ self._get_matrix(g) for g in gens ])

    def _calculate_simple_gens(self):
        """
        Calculate generators for self, where both the domain and
        codomain for self are assumed to be simple abelian varieties.
        The saturation of the span of these generators in self will be
        the full space of homomorphisms from the domain of self to its
        codomain.

        EXAMPLES:
            sage: H = Hom(J0(11), J0(22)[0])
            sage: H._calculate_simple_gens()
            [[1 0]
            [1 1]]
            sage: J = J0(11) * J0(33) ; J.decomposition()
            [
            Simple abelian subvariety 11a(1,11) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 11a(1,33) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 11a(3,33) of dimension 1 of J0(11) x J0(33),
            Simple abelian subvariety 33a(1,33) of dimension 1 of J0(11) x J0(33)
            ]
            sage: J[0].Hom(J[1])._calculate_simple_gens()
            [[ 0 -1]
            [ 1 -1]]
            sage: J[0].Hom(J[2])._calculate_simple_gens()
            [[-1  0]
            [-1 -1]]
            sage: J[0].Hom(J[0])._calculate_simple_gens()
            [[1 0]
            [0 1]]
            sage: J[1].Hom(J[2])._calculate_simple_gens()
            [[ 0 -4]
            [ 4  0]]

            sage: J = J0(23) ; J.decomposition()
            [
            Simple abelian variety J0(23) of dimension 2
            ]
            sage: J[0].Hom(J[0])._calculate_simple_gens()
            [[1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1],
             [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]]
            sage: J.hecke_operator(2).matrix()
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]

            sage: H = Hom(J0(11), J0(22)[0])
            sage: H._calculate_simple_gens()
            [[1 0]
            [1 1]]
        """
        A = self.domain()
        B = self.codomain()

        # TODO: POSSIBLE BUG ALERT!!! If the groups() are different this can be wrong, e.g., J0(11) and J1(11)
        if A.newform_label() != B.newform_label():
            return []

        f = A._isogeny_to_newform_abelian_variety()
        g = B._isogeny_to_newform_abelian_variety().complementary_isogeny()

        Af = f.codomain()
        ls = Af._calculate_endomorphism_generators()

        Mf = f.matrix()
        Mg = g.matrix()

        return [ Mf * self._get_matrix(e) * Mg for e in ls ]

# NOTE/WARNING/TODO:  Below in the __init__, etc. we do *not* check
# that the input gens are give something that spans a sub*ring*, as apposed
# to just a subgroup.
class EndomorphismSubring(Homspace, Ring):

    def __init__(self, A, gens=None):
        """
        A subring of the endomorphism ring.

        INPUT:
            A -- an abelian variety
            gens -- (default: None); optional; if given should be a
                 tuple of the generators as matrices

        EXAMPLES:
            sage: J0(23).endomorphism_ring()
            Endomorphism ring of Abelian variety J0(23) of dimension 2
            sage: sage.modular.abvar.homspace.EndomorphismSubring(J0(25))
            Endomorphism ring of Abelian variety J0(25) of dimension 0
            sage: type(J0(11).endomorphism_ring())
            <class 'sage.modular.abvar.homspace.EndomorphismSubring'>
        """
        self._J = A.ambient_variety()
        self._A = A
        Homspace.__init__(self, A, A, A.category())
        Ring.__init__(self, A.base_ring())

        if gens is None:
            self._gens = None
        else:
            self._gens = tuple([ self._get_matrix(g) for g in gens ])
        self._is_full_ring = gens is None

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: J0(31).endomorphism_ring()._repr_()
            'Endomorphism ring of Abelian variety J0(31) of dimension 2'
            sage: J0(31).endomorphism_ring().image_of_hecke_algebra()._repr_()
            'Subring of endomorphism ring of Abelian variety J0(31) of dimension 2'
        """
        if self._is_full_ring:
            return "Endomorphism ring of %s" % self._A
        else:
            return "Subring of endomorphism ring of %s" % self._A

    def abelian_variety(self):
        """
        Return the abelian variety that this endomorphism ring
        is attached to.

        EXAMPLES:
            sage: J0(11).endomorphism_ring().abelian_variety()
            Abelian variety J0(11) of dimension 1
        """
        return self._A


#     def calculate_generators(self):
#         """
#         Calculate a set of generators for self.
#
#         EXAMPLES:
#             sage: J0(11).endomorphism_ring().calculate_generators()
#             BOOM! This is broken .  TODO -- fix me.
#         """
#         if self._gens is None:
#             gens = self._A._calculate_endomorphism_generators()
#             M = ZZ**(4*self._A.dimension()**2)
#             gens = M.submodule([ x.matrix().list() for x in gens ]).saturation().basis()
#             self._gens = tuple([ self._get_matrix(g) for g in gens ])


    def index_in(self, other, check=True):
        """
        Return the index of self in other.

        INPUT:
            other -- another endomorphism subgring of the same abelian variety
            check -- bool (default: True); whether to do some type and other consistency checks

        EXAMPLES:
            sage: R = J0(33).endomorphism_ring()
            sage: R.index_in(R)
            boom.

            AND I would like to do an example with the Hecke algebra.
        """
        if check:
            if not isinstance(other, EndomorphismSubring):
                raise ValueError, "other must be a subring of an endomorphism ring of an abelian variety."
            if not (self.abelian_variety() == other.abelian_variety()):
                raise ValueError, "self and other must be endomorphisms of the same abelian variety"

        M = self.free_module()
        N = other.free_module()
        return M.index_in(N)

    def discriminant(self):
        """
        Return the discriminant of this ring, which is the
        discriminant of the trace pairing.

        EXAMPLES:
             sage: J0(33).endomorphism_ring().discriminant()
             boom!
        """
        g = self.gens()
        M = Matrix(ZZ,len(g), [ (g[i]*g[j]).trace()
                                for i in range(len(g)) for j in range(len(g)) ])
        return M.determinant()

    def image_of_hecke_algebra(self):
        """
        Compute the image of the Hecke algebra inside this endomorphism
        subring.

        EXAMPLES:
            sage: E = J0(33).endomorphism_ring()
            sage: E.image_of_hecke_algebra()
            boom!! it's broken.
        """
        try:
            return self.__hecke_algebra_image
        except AttributeError:
            pass

        A = self.abelian_variety()
        if not A.is_hecke_stable():
            raise ValueError, "ambient variety is not Hecke stable"

        M = A.modular_symbols()

        d = A.dimension()
        EndVecZ = ZZ**(4*d**2)
        T_matrices = [ A.hecke_operator(n).matrix().list() for n in
                       range(1,M.sturm_bound()+1) ]
        W = EndVecZ.submodule(T_matrices)

        T = EndomorphismSubring(A, W.basis())
        self.__hecke_algebra_image = T
        return self.__hecke_algebra_image



