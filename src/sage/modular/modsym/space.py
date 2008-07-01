"""
Space of modular symbols (base class)

All the spaces of modular symbols derive from this class.  This class
is an abstract base class.
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
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
#*****************************************************************************


import math
import random
import weakref

import sage.modules.free_module as free_module
import sage.matrix.matrix_space as matrix_space
import sage.modules.free_module_morphism as free_module_morphism
from   sage.modules.free_module_element  import is_FreeModuleElement
import sage.misc.misc as misc
import sage.modular.dims as dims
import sage.modular.hecke.all as hecke
import sage.modular.modsym.element
import sage.structure.parent_gens as gens
import sage.rings.arith as arith
from   sage.rings.all import PowerSeriesRing, Integer, O, QQ, ZZ, is_NumberField
from   sage.structure.all import Sequence, SageObject
import sage.modular.modsym.ambient

def is_ModularSymbolsSpace(x):
    return isinstance(x, ModularSymbolsSpace)

class ModularSymbolsSpace(hecke.HeckeModule_free_module):
    def __init__(self, group, weight, character, sign, base_ring):
        """
        Create a space of modular symbols.

        EXAMPLES:
            sage: M = ModularSymbols(22,6) ; M
            Modular Symbols space of dimension 30 for Gamma_0(22) of weight 6 with sign 0 over Rational Field
            sage: M == loads(dumps(M))
            True
        """
        self.__group = group
        self.__character = character
        self.__sign = sign
        hecke.HeckeModule_free_module.__init__(self, base_ring, group.level(), weight)

    def __cmp__(self, other):
        """
        Compare self and other.

        When spaces are in a common ambient space, we order
        lexicographically by the sequence of traces of Hecke operators
        $T_p$, for all primes $p$.  In general we order first by the
        group, then the weight, then the character, then the sign then
        the base ring, then the dimension.

        EXAMPLES:
            sage: M = ModularSymbols(21,4) ; N = ModularSymbols(Gamma1(5),6)
            sage: M.cuspidal_submodule().__cmp__(N)
            1
            sage: M.cuspidal_submodule() == N
            False
        """
        if not isinstance(other, ModularSymbolsSpace):
            return cmp(type(self), type(other))
        c = cmp(self.__group,other.__group)
        if c: return c
        c = cmp(self.weight(), other.weight())
        if c: return c
        c = cmp(self.__character, other.__character)
        if c: return c
        c = cmp(self.__sign, other.__sign)
        if c: return c
        c = cmp(self.base_ring(), other.base_ring())
        if c: return c
        if self.is_ambient() or other.is_ambient():
            # if one is ambient they are equal iff they have same
            # dimension at this point, since all defining properties
            # are the same, so they are in the same ambient space.
            return cmp(self.dimension(), other.dimension())

        c = cmp(self.ambient_hecke_module(), other.ambient_hecke_module())
        if c: return c
        c = cmp(self.dimension(), other.dimension())
        if c: return c
        d = cmp(self.free_module(), other.free_module())
        if d == 0:
            return d
        # distinguish using Hecke operators, if possible.
        try:
            for p in arith.prime_range(self.hecke_bound()):
                ap = self.hecke_matrix(p).trace()
                bp = other.hecke_matrix(p).trace()
                c = cmp(ap, bp)
                if c: return c
        except ArithmeticError:
            pass
        # fallback on subspace comparison
        return d


    def character(self):
        """
        Return the character associated to self.

        EXAMPLES:
            sage: ModularSymbols(12,8).character()
            [1, 1]
            sage: ModularSymbols(DirichletGroup(25).0, 4).character()
            [zeta20]
        """
        return self.__character

    def cuspidal_submodule(self):
        """
        Return the cuspidal submodule of self.

        NOTE: This should be overridden by all derived classes.

        EXAMPLES:
            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).cuspidal_submodule()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of cuspidal submodule not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).cuspidal_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
        """
        raise NotImplementedError, "computation of cuspidal submodule not yet implemented for this class"

    def cuspidal_subspace(self):
        """
        Synonym for cuspidal_submodule.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.cuspidal_subspace().new_subspace().dimension()
            2
        """
        return self.cuspidal_submodule()

    def new_subspace(self, p=None):
        """
        Synonym for new_submodule.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma0(5),12); m.dimension()
            12
            sage: m.new_subspace().dimension()
            6
            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.new_subspace().dimension()
            2
        """
        return self.new_submodule(p)

    def old_subspace(self, p=None):
        """
        Synonym for old_submodule.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.old_subspace().dimension()
            6
        """
        return self.old_submodule(p)

    def eisenstein_subspace(self):
        """
        Synonym for eisenstein_submodule.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma1(3),12); m.dimension()
            8
            sage: m.eisenstein_subspace().dimension()
            2
            sage: m.cuspidal_subspace().dimension()
            6
        """
        return self.eisenstein_submodule()

    def dimension_of_associated_cuspform_space(self):
        """
        Return the dimension of the corresponding space of cusp forms.

        The input space must be cuspidal, otherwise there is no
        corresponding space of cusp forms.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma0(389),2).cuspidal_subspace(); m.dimension()
            64
            sage: m.dimension_of_associated_cuspform_space()
            32
            sage: m = ModularSymbols(Gamma0(389),2,sign=1).cuspidal_subspace(); m.dimension()
            32
            sage: m.dimension_of_associated_cuspform_space()
            32
        """
        if not self.is_cuspidal():
            raise ArithmeticError, "space must be cuspidal"
        if self.sign() == 0:
            return self.dimension() // 2
        return self.dimension()

    def dual_star_involution_matrix(self):
        """
        Return the matrix of the dual star involution, which is
        induced by complex conjugation on the linear dual of modular
        symbols.

        NOTE: This should be overridden in all derived classes.

        EXAMPLES:
            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).dual_star_involution_matrix()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of dual star involution matrix not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).dual_star_involution_matrix()
            [ 1  0  0]
            [ 0 -1  0]
            [ 0  1  1]
        """
        raise NotImplementedError, "computation of dual star involution matrix not yet implemented for this class"

    def group(self):
        """
        Returns the group of this modular symbols space.

        INPUT:
           ModularSymbols self -- an arbitrary space of modular symbols

        OUTPUT:
           CongruenceSubgroup -- the congruence subgroup that this is a space
                              of modular symbols for.

        ALGORITHM:
           The group is recorded when this space is created.

        EXAMPLES:
            sage: m = ModularSymbols(20)
            sage: m.group()
            Congruence Subgroup Gamma0(20)
        """
        return self.__group

    def is_ambient(self):
        """
        Return True if self is an ambient space of modular symbols.

        EXAMPLES:
            sage: ModularSymbols(21,4).is_ambient()
            True
            sage: ModularSymbols(21,4).cuspidal_submodule().is_ambient()
            False
        """
        return isinstance(self, sage.modular.modsym.ambient.ModularSymbolsAmbient)

    def is_cuspidal(self):
        """
        Return True if self is a cuspidal space of modular symbols.

        NOTE: This should be overridden in all derived classes.

        EXAMPLES:
            sage: sage.modular.modsym.space.ModularSymbolsSpace(Gamma0(11),2,DirichletGroup(11).gens()[0]**10,0,QQ).is_cuspidal()
            Traceback (most recent call last):
            ...
            NotImplementedError: computation of cuspidal subspace not yet implemented for this class
            sage: ModularSymbols(Gamma0(11),2).is_cuspidal()
            False
        """
        raise NotImplementedError, "computation of cuspidal subspace not yet implemented for this class"

    def is_simple(self):
        """
        Return whether not this modular symbols space is simple as a module
        over the anemic Hecke algebra adjoin *.

        EXAMPLES:
            sage: m = ModularSymbols(Gamma0(33),2,sign=1)
            sage: m.is_simple()
            False
            sage: o = m.old_subspace()
            sage: o.decomposition()
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field,
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field
            ]
            sage: C=ModularSymbols(1,14,0,GF(5)).cuspidal_submodule()
            sage: C
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(1) of weight 14 with sign 0 over Finite Field of size 5
            sage: C.is_simple()
            True
        """
        try:
            return self._is_simple
        except AttributeError:
            D = self.factorization()
            if len(D) == 0 or len(D) == 1 and D[0][1] == 1:
                self._is_simple = True
            else:
                self._is_simple = False
            return self._is_simple

    def multiplicity(self, S, check_simple=True):
        """
        Return the multiplicity of the simple modular symbols space S
        in self.  S must be a simple anemic Hecke module.

        ASSUMPTION: self is an anemic Hecke module with the same weight and
        group as S, and S is simple.

        EXAMPLES:
            sage: M = ModularSymbols(11,2,sign=1)
            sage: N1, N2 = M.decomposition()
            sage: N1.multiplicity(N2)
            0
            sage: M.multiplicity(N1)
            1
            sage: M.multiplicity(ModularSymbols(14,2))
            0
        """
        if self.level() % S.level() != 0 or S.weight() != self.weight():
            return 0
        if check_simple and not S.is_simple():
            raise ArithmeticError, "S must be simple"
        A = self.ambient_hecke_module()
        B = A.submodule_generated_by_images(S)
        C = self.intersection(B)
        d = C.rank()
        n = S.rank()
        assert d % n == 0, "the dimension of intersection must be a multiple of dimension of simple space.  bug!"
        return d//n

    def ngens(self):
        """
        The number of generators of self.

        INPUT:
           ModularSymbols self -- arbitrary space of modular symbols.
        OUTPUT:
           int -- the number of generators, which is the same as the
                  dimension of self.

        ALGORITHM:
           Call the dimension function.

        EXAMPLES:
            sage: m = ModularSymbols(33)
            sage: m.ngens()
            9
            sage: m.rank()
            9
            sage: ModularSymbols(100, weight=2, sign=1).ngens()
            18
        """
        return self.rank()


    #########################################################################
    #
    #  Computation of q-expansions
    #
    #########################################################################

    def default_prec(self):
        r"""
        Get the default precision for computation of $q$-expansion
        associated to the ambient space of this space of modular
        symbols (and all subspaces).  Use \code{set_default_prec} to
        change the default precision.

        EXAMPLES:
            sage: M = ModularSymbols(15)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [
            q - q^2 - q^3 - q^4 + q^5 + q^6 + O(q^8)
            ]
            sage: M.set_default_prec(20)

        Notice that setting the default precision of the ambient space
        affects the subspaces.

            sage: M.cuspidal_submodule().q_expansion_basis()
            [
            q - q^2 - q^3 - q^4 + q^5 + q^6 + 3*q^8 + q^9 - q^10 - 4*q^11 + q^12 - 2*q^13 - q^15 - q^16 + 2*q^17 - q^18 + 4*q^19 + O(q^20)
            ]
            sage: M.cuspidal_submodule().default_prec()
            20
        """
        if not self.is_ambient():
            return self.ambient_hecke_module().default_prec()
        try:
            return self.__default_prec
        except AttributeError:
            self.__default_prec = Integer(8)
            return self.__default_prec

    def set_default_prec(self, prec):
        """
        Set the default precision for computation of $q$-expansion
        associated to the ambient space of this space of modular
        symbols (and all subspaces).

        EXAMPLES:
            sage: M = ModularSymbols(Gamma1(13),2)
            sage: M.set_default_prec(5)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [
            q - 4*q^3 - q^4 + O(q^5),
            q^2 - 2*q^3 - q^4 + O(q^5)
            ]
        """
        if not self.is_ambient():
            return self.ambient_hecke_module().set_default_prec(prec)
        else:
            self.__default_prec = Integer(prec)

    def set_precision(self, prec):
        """
        Same as self.set_default_prec(prec).

        EXAMPLES:
            sage: M = ModularSymbols(17,2)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [
            q - q^2 - q^4 - 2*q^5 + 4*q^7 + O(q^8)
            ]
            sage: M.set_precision(10)
            sage: M.cuspidal_submodule().q_expansion_basis()
            [
            q - q^2 - q^4 - 2*q^5 + 4*q^7 + 3*q^8 - 3*q^9 + O(q^10)
            ]
        """
        self.set_default_prec(prec)

    def q_expansion_basis(self, prec=None, algorithm='default'):
        r"""
        Returns a basis of q-expansions (as power series) to precision
        prec of the space of modular forms associated to self.  The
        q-expansions are defined over the same base ring as self, and
        a put in echelon form.

        INPUT:
            self -- a space of CUSPIDAL modular symbols
            prec -- an integer
            algorithm -- string:
                    'default' (default) -- decide which algorithm to use based on heuristics
                    'hecke' -- compute basis by computing homomorphisms
                               T --> K, where T is the Hecke algebra
                    'eigen' -- compute basis using eigenvectors for the Hecke action
                               and Atkin-Lehner-Li theory to patch them together
                    'all'   -- compute using hecke_dual and eigen algorithms and verify
                               that the results are the same.

        The computed basis is \emph{not} cached, though of course Hecke
        operators used in computing the basis are cached.

        EXAMPLES:
            sage: M = ModularSymbols(1, 12).cuspidal_submodule()
            sage: M.q_expansion_basis(8)
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)
            ]

            sage: M.q_expansion_basis(8, algorithm='eigen')
            [
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)
            ]


            sage: M = ModularSymbols(1, 24).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [
            q + 195660*q^3 + 12080128*q^4 + 44656110*q^5 - 982499328*q^6 - 147247240*q^7 + O(q^8),
            q^2 - 48*q^3 + 1080*q^4 - 15040*q^5 + 143820*q^6 - 985824*q^7 + O(q^8)
            ]

            sage: M = ModularSymbols(11, 2, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 + O(q^8)
            ]

            sage: M = ModularSymbols(Gamma1(13), 2, sign=1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')
            [
            q - 4*q^3 - q^4 + 3*q^5 + 6*q^6 + O(q^8),
            q^2 - 2*q^3 - q^4 + 2*q^5 + 2*q^6 + O(q^8)
            ]


            sage: M = ModularSymbols(Gamma1(5), 3, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8, algorithm='eigen')   # dimension is 0
            []

            sage: M = ModularSymbols(Gamma1(7), 3, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_basis(8)
            [
            q - 3*q^2 + 5*q^4 - 7*q^7 + O(q^8)
            ]

            sage: M = ModularSymbols(43, 2, sign=0).cuspidal_submodule()
            sage: M[0]
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field
            sage: M[0].q_expansion_basis()
            [
            q - 2*q^2 - 2*q^3 + 2*q^4 - 4*q^5 + 4*q^6 + O(q^8)
            ]
            sage: M[1]
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field
            sage: M[1].q_expansion_basis()
            [
            q + 2*q^5 - 2*q^6 - 2*q^7 + O(q^8),
            q^2 - q^3 - q^5 + q^7 + O(q^8)
            ]
        """
        if prec is None:
            prec = self.default_prec()
        else:
            prec = Integer(prec)

        if prec < 1:
            raise ValueError, "prec (=%s) must be >= 1"%prec

        if not self.is_cuspidal():
            raise ArithmeticError, "space must be cuspidal"

        if self.sign() == 0:
            P = self.plus_submodule(compute_dual=True)
            return Sequence(P.q_expansion_basis(prec=prec, algorithm=algorithm), cr=True)

        if self.dimension() == 0:
            return Sequence([])

        if algorithm == 'default':
            algorithm = 'hecke'
        if algorithm == 'hecke':
            return Sequence(self._q_expansion_basis_hecke_dual(prec), cr=True)
        elif algorithm == 'eigen':
            return Sequence(self._q_expansion_basis_eigen(prec, 'alpha'), cr=True)
        elif algorithm == 'all':
            B1 = self._q_expansion_basis_hecke_dual(prec)
            B2 = self._q_expansion_basis_eigen(prec, 'alpha')
            if B1 != B2:
                raise RuntimeError, "There is a bug in q_expansion_basis -- basis computed differently with two algorithms:\n%s\n%s\n"%(B1, B2,)
            return Sequence(B1, cr=True)
        else:
            raise ValueError, "no algorithm '%s'"%algorithm

    def q_expansion_module(self, prec = None, R=None):
        r"""
        Return a basis over R for the space spanned by the coefficient
        vectors of the $q$-expansions corresponding to self.  If R is
        not the base ring of self, returns the restriction of scalars
        down to R (for this, self must have base ring $\QQ$ or a
        number field).

        INPUT:
             self -- must be cuspidal
             prec -- an integer (default: self.default_prec())
             R -- either ZZ, QQ, or the base_ring of self (which is the default)

        OUTPUT:
             A free module over R.

        TODO -- extend to more general R (though that is fairly easy
        for the user to get by just doing base_extend or change_ring
        on the output of this function).

        Note that the prec needed to distinguish elements of the
        restricted-down-to-R basis may be bigger than
        \code{self.hecke_bound()}, since one must use the Sturm bound
        for modular forms on $\Gamma_H(N)$.

        INPUT:
            prec -- integer

        OUTPUT:
            A QQ-vector space


        EXAMPLES WITH SIGN 1 and R=QQ:

        Basic example with sign 1:
            sage: M = ModularSymbols(11, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        Same example with sign -1:
            sage: M = ModularSymbols(11, sign=-1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        An example involving old forms:
            sage: M = ModularSymbols(22, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -2]
            [ 0  0  1  0 -2]

        An example that (somewhat spuriously) is over a number field:
            sage: x = polygen(QQ)
            sage: k = NumberField(x^2+1, 'a')
            sage: M = ModularSymbols(11, base_ring=k, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(5, QQ)
            Vector space of degree 5 and dimension 1 over Rational Field
            Basis matrix:
            [ 0  1 -2 -1  2]

        An example that involves an eigenform with coefficients in a number field:
            sage: M = ModularSymbols(23, sign=1).cuspidal_submodule()
            sage: M.q_eigenform(4, 'gamma')
            q + gamma*q^2 + (-2*gamma - 1)*q^3 + O(q^4)
            sage: M.q_expansion_module(11, QQ)
            Vector space of degree 11 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -1  0 -2  2 -1  2  2]
            [ 0  0  1 -2 -1  2  1  2 -2  0 -2]

        An example that is genuinely over a base field besides QQ.
            sage: eps = DirichletGroup(11).0
            sage: M = ModularSymbols(eps,3,sign=1).cuspidal_submodule(); M
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 and level 11, weight 3, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4
            sage: M.q_eigenform(4, 'beta')
            q + (-zeta10^3 + 2*zeta10^2 - 2*zeta10)*q^2 + (2*zeta10^3 - 3*zeta10^2 + 3*zeta10 - 2)*q^3 + O(q^4)
            sage: M.q_expansion_module(7, QQ)
            Vector space of degree 7 and dimension 4 over Rational Field
            Basis matrix:
            [  0   1   0   0   0 -40  64]
            [  0   0   1   0   0 -24  41]
            [  0   0   0   1   0 -12  21]
            [  0   0   0   0   1  -4   4]

        An example involving an eigenform rational over the base, but the base is not QQ.
            sage: k.<a> = NumberField(x^2-5)
            sage: M = ModularSymbols(23, base_ring=k, sign=1).cuspidal_submodule()
            sage: D = M.decomposition(); D
            [
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Number Field in a with defining polynomial x^2 - 5,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Number Field in a with defining polynomial x^2 - 5
            ]
            sage: M.q_expansion_module(8, QQ)
            Vector space of degree 8 and dimension 2 over Rational Field
            Basis matrix:
            [ 0  1  0 -1 -1  0 -2  2]
            [ 0  0  1 -2 -1  2  1  2]

        An example involving an eigenform not rational over the base and for which the base is not QQ.
            sage: eps = DirichletGroup(25).0^2
            sage: M = ModularSymbols(eps,2,sign=1).cuspidal_submodule(); M
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 25, weight 2, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4
            sage: D = M.decomposition(); D
            [
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 and level 25, weight 2, character [zeta10], sign 1, over Cyclotomic Field of order 10 and degree 4
            ]
            sage: D[0].q_eigenform(4, 'mu')
            q + mu*q^2 + ((zeta10^3 + zeta10 - 1)*mu + zeta10^2 - 1)*q^3 + O(q^4)
            sage: D[0].q_expansion_module(11, QQ)
            Vector space of degree 11 and dimension 8 over Rational Field
            Basis matrix:
            [  0   1   0   0   0   0   0   0 -20  -3   0]
            [  0   0   1   0   0   0   0   0 -16  -1   0]
            [  0   0   0   1   0   0   0   0 -11  -2   0]
            [  0   0   0   0   1   0   0   0  -8  -1   0]
            [  0   0   0   0   0   1   0   0  -5  -1   0]
            [  0   0   0   0   0   0   1   0  -3  -1   0]
            [  0   0   0   0   0   0   0   1  -2   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1]
            sage: D[0].q_expansion_module(11)
            Vector space of degree 11 and dimension 2 over Cyclotomic Field of order 10 and degree 4
            Basis matrix:
            [                                  0                                   1                                   0                        zeta10^2 - 1                       -zeta10^2 - 1                -zeta10^3 - zeta10^2                   zeta10^2 - zeta10           2*zeta10^3 + 2*zeta10 - 1    zeta10^3 - zeta10^2 - zeta10 + 1        zeta10^3 - zeta10^2 + zeta10   -2*zeta10^3 + 2*zeta10^2 - zeta10]
            [                                  0                                   0                                   1               zeta10^3 + zeta10 - 1                         -zeta10 - 1                -zeta10^3 - zeta10^2 -2*zeta10^3 + zeta10^2 - zeta10 + 1                            zeta10^2                                   0                        zeta10^3 + 1  2*zeta10^3 - zeta10^2 + zeta10 - 1]

        EXAMPLES WITH SIGN 0 and R=QQ:
        ** TODO -- this doesn't work yet -- not implemented!! **
             M = ModularSymbols(11,2).cuspidal_submodule()
             M.q_expansion_module()
            ... boom ...


        EXAMPLES WITH SIGN 1 and R=ZZ (computes saturation):

            sage: M = ModularSymbols(43,2, sign=1).cuspidal_submodule()
            sage: M.q_expansion_module(8, QQ)
            Vector space of degree 8 and dimension 3 over Rational Field
            Basis matrix:
            [   0    1    0    0    0    2   -2   -2]
            [   0    0    1    0 -1/2    1 -3/2    0]
            [   0    0    0    1 -1/2    2 -3/2   -1]
            sage: M.q_expansion_module(8, ZZ)
            Free module of degree 8 and rank 3 over Integer Ring
            Echelon basis matrix:
            [ 0  1  0  0  0  2 -2 -2]
            [ 0  0  1  1 -1  3 -3 -1]
            [ 0  0  0  2 -1  4 -3 -2]
        """
        if prec is None:
            prec = self.default_prec()
        if R == ZZ:
            return self._q_expansion_module_integral(prec)
        elif R == QQ:
            return self._q_expansion_module_rational(prec)
        elif R is None or R == self.base_ring():
            ## names is never used in this case
            return self._q_expansion_module(prec)
        else:
            raise NotImplementedError, "R must be ZZ, QQ, or the base ring of the modular symbols space."

    def _q_eigenform_images(self, A, prec, names):
        """
        Return list of images in space corresponding to self of
        eigenform corresponding to A under the degeneracy maps. This
        is mainly a helper function for other internal functions.

        INPUT:
             self -- space of modular symbols
             A    -- cuspidal simple space of level dividing the level of self
                     and the same weight
             prec -- a positive integer

        EXAMPLES:
            sage: M = ModularSymbols(33,2,sign=1)
            sage: A = M.modular_symbols_of_level(11).cuspidal_submodule()
            sage: M._q_eigenform_images(A, 10, names='a')
            [q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 - 2*q^9 + O(q^10),
             q^3 - 2*q^6 - q^9 + 2*q^12 + q^15 + 2*q^18 - 2*q^21 - 2*q^27 + O(q^30)]
        """
        f = A.q_eigenform(prec, names)
        if A.level() == self.level():
            return [f]
        D = arith.divisors(self.level() // A.level())
        q = f.parent().gen()
        return [f] + [f(q**d) for d in D if d > 1]

    def _q_expansion_module(self, prec, algorithm='hecke'):
        """
        Return module spanned by the $q$-expansions corresponding to
        self.
        """
        if not self.is_cuspidal():
            raise ValueError, "self must be cuspidal"
        R = self.base_ring()
        if not R.is_field():
            if R == ZZ:
                return self._q_expansion_module_integral(prec)
            raise NotImplementedError, "base ring must be a field (or ZZ)."

        if algorithm == 'hecke' or algorithm == 'default':
            A = R ** prec
            return A.span([f.padded_list(prec) for f in self.q_expansion_basis(prec, algorithm)])

        if algorithm != 'eigen':
            raise ValueError, "unknown algorithm '%s'"%algorithm

        V = R ** prec
        def q_eigen_gens(f):
            X = f.padded_list(prec)
            d = A.dimension()
            if d == 1:
                # X is just a list of elements of R
                return [X]
            else:
                # X is a list of elements of a poly quotient ring
                return [[X[i][j] for i in xrange(prec)] for j in xrange(d)]

        if self.sign() == 0:
            X = self.plus_submodule(compute_dual=True)
        else:
            X = self

        B = [sum([q_eigen_gens(f) for f in self._q_eigenform_images(A, prec, 'zeta')], []) for A, _ in X.factorization()]

        A = R ** prec
        return A.span(sum(B, []))

    def _q_expansion_module_rational(self, prec):
        """
        Return a free module over $\ZZ$ for the space spanned by the
        $q$-expansions corresponding to self.  The base ring of self
        must be $\QQ$ or a number field, and self must be cuspidal.
        The returned space is a $\ZZ$-module, where the coordinates
        are the coefficients of $q$-expansions.
        """
        if not self.is_cuspidal():
            raise ValueError, "self must be cuspidal"
        K = self.base_ring()
        if not (K == QQ or is_NumberField(K)):
            raise TypeError, "self must be over QQ or a number field."
        n = K.degree()
        if n == 1:
            return self._q_expansion_module(prec)

        # Construct the vector space over QQ of dimension equal to
        # the degree of the base field times the dimension over C
        # of the space of cusp forms corresponding to self.
        V = QQ**prec ## is this needed?
        def q_eigen_gens(f):
            # Return restricted down to QQ gens for cusp space corresponding
            # to the simple factor A.
            X = f.padded_list(prec)
            d = A.dimension()
            if d == 1:
                return [[X[i][j] for i in xrange(prec)] for j in xrange(n)]
            else:
                # This looks like it might be really slow -- though
                # perhaps it's nothing compared to the time taken by
                # whatever computed this in the first place.
                return [[(X[i].list())[j][k] for i in xrange(prec)] for j in xrange(d) for k in range(n)]
        if self.sign() == 0:
            X = self.plus_submodule(compute_dual=True)
        else:
            X = self

        B = [sum([q_eigen_gens(f) for f in self._q_eigenform_images(A, prec, 'alpha')], []) for A, _ in X.factorization()]
        A = QQ**prec
        W = A.span(sum(B, []))
        return W


    def _q_expansion_module_integral(self, prec):
        r"""
        Return module over $\ZZ$ for the space spanned by the
        $q$-expansions corresponding to self.  The base ring of self
        must be $\QQ$ or a number field, and self must be cuspidal.
        The returned space is a $\ZZ$-module, where the coordinates
        are the coefficients of $q$-expansions.
        """
        V = self.q_expansion_module(prec, QQ)
        return free_module.span(V.basis(), ZZ).saturation()


    def congruence_number(self, other, prec=None):
        r"""
        Given two cuspidal spaces of modular symbols, compute the
        congruence number, using prec terms of the $q$-expansions.

        The congruence number is defined as follows.  If $V$ is the
        submodule of integral cusp forms corresponding to self
        (satured in $\Z[[q]]$, by definition) and $W$ is the submodule
        corresponding to other, each computed to precision prec, the
        congruence number is the index of $V+W$ in its saturation in
        $\Z[[q]]$.

        If prec is not given it is set equal to the max of the
        \code{hecke_bound} function called on each space.
        """
        if not self.is_cuspidal():
            raise ValueError, "self must be cuspidal"
        if not other.is_cuspidal():
            raise ValueError, "right must be cuspidal"
        if prec is None:
            prec = max(self.hecke_bound(), other.hecke_bound())
        prec = int(prec)

        V =  self.q_expansion_module(prec, ZZ)
        W = other.q_expansion_module(prec, ZZ)
        K = V+W
        return K.index_in_saturation()

    #########################################################################
    #
    #  Computation of a basis using eigenforms
    #
    #########################################################################

    def q_eigenform(self, prec, names=None):
        """
        Returns the q-expansion to precision prec of a new eigenform
        associated to self, where self must be new, cuspidal, and
        simple.
        """
        if self.dimension() > 1 and names is None:
            raise ValueError, "please specify a name to use for the field of eigenvalues"

        if prec is None:
            prec = self.default_prec()
        try:
            f = self._q_expansion_dict[names]
        except (AttributeError, KeyError):
            self._q_expansion_dict = {}
            if not self.is_cuspidal():
                raise ArithmeticError, "self must be cuspidal."

            if not self.is_simple():
                if self.sign() == 0:
                    return self.plus_submodule(compute_dual=True).q_eigenform(prec, names)
                raise ArithmeticError, "self must be simple."
            a2 = self.eigenvalue(2, names)
            R = PowerSeriesRing(a2.parent(), "q")
            q = R.gen(0)
            f = q + a2*q**2 + O(q**3)

        if f.prec() < prec:
            R = f.parent()
            ext = [self.eigenvalue(n, names) for n in range(f.prec(),prec)]
            f = R(f.padded_list(f.prec()) + ext)
            self._q_expansion_dict[names] = f.add_bigoh(prec)
            return self._q_expansion_dict[names]
        else:
            return f.O(prec)

    def _q_expansion_basis_eigen(self, prec, names):
        if self.is_simple():
            f = self.q_eigenform(prec, names)
            R = PowerSeriesRing(self.base_ring(), 'q')
            B = [R([f[i][j] for i in xrange(prec)],prec) for j in range(self.rank())]
            return B
        else:
            raise NotImplementedError


    #########################################################################
    #
    #  Computation of a basis using linear functionals on the Hecke algebra.
    #
    #########################################################################

    def q_expansion_cuspforms(self, prec=None):
        """
        Returns a function f(i,j) such that each value f(i,j) is the
        q-expansion, to the given precision, of an element of the
        corresponding space~$S$ of cusp forms.  Together these
        functions span~$S$.  Here $i,j$ are integers with
        $0\leq i,j < d$, where $d$ is the dimension of self.

        For a reduced echelon basis, use the function
        \code{q_expansion_basis} instead.

        More precisely, this function returns the $q$-expansions
        obtained by taking the $ij$ entry of the matrices of the Hecke
        operators $T_n$ acting on the subspace of the linear dual of
        modular symbols corresponding to self.

        EXAMPLES:
            sage: S = ModularSymbols(11,2, sign=1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 2*q^2 - q^3 + 2*q^4 + q^5 + 2*q^6 - 2*q^7 + O(q^8)

            sage: S = ModularSymbols(37,2).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q + q^3 - 2*q^4 - q^7 + O(q^8)
            sage: f(3,3)
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + 6*q^6 - q^7 + O(q^8)
            sage: f(1,2)
            q^2 + 2*q^3 - 2*q^4 + q^5 - 3*q^6 + O(q^8)

            sage: S = ModularSymbols(Gamma1(13),2,sign=-1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 2*q^2 + q^4 - q^5 + 2*q^6 + O(q^8)
            sage: f(0,1)
            q^2 - 2*q^3 - q^4 + 2*q^5 + 2*q^6 + O(q^8)

            sage: S = ModularSymbols(1,12,sign=-1).cuspidal_submodule()
            sage: f = S.q_expansion_cuspforms(8)
            sage: f(0,0)
            q - 24*q^2 + 252*q^3 - 1472*q^4 + 4830*q^5 - 6048*q^6 - 16744*q^7 + O(q^8)
        """
        if prec is None:
            prec = self.default_prec()
        if not self.is_cuspidal():
            raise ArithmeticError, "self must be cuspidal"
        K = self.base_ring()
        M = matrix_space.MatrixSpace(K, prec-1, self.dimension())
        T = [self.dual_hecke_matrix(n) for n in range(1,prec)]
        R = PowerSeriesRing(self.base_ring(), 'q')
        def form(i, j):
            return R([0] + [t[i,j] for t in T], prec)
        return form

    def _q_expansion_basis_hecke_dual(self, prec):
        d = self.dimension_of_associated_cuspform_space()
        prec = Integer(prec)
        if prec < 1:
            raise ValueError, "prec (=%s) must be >= 1"%prec
        if d >= prec-1:
            d = prec-1
        K = self.base_ring()

        A = free_module.VectorSpace(K, prec-1)
        M = matrix_space.MatrixSpace(K, prec-1, self.dimension())

        V = A.zero_submodule()
        i = self.dimension()-1
        j = 0

        t = misc.verbose('computing basis to precision %s'%prec)
        while V.dimension() < d and i >= 0:
            v = [self.dual_hecke_matrix(n).column(i) for n in range(1,prec)]
            t = misc.verbose('iteration: %s'%j,t)
            X = M(v).transpose()
            V += X.row_space()
            t = misc.verbose('addition of row space: %s'%j,t)
            i -= 1
            j += 1

        R = PowerSeriesRing(K, 'q')
        B = V.basis()
        if len(B) < d:
            B += [V(0)] * (d-len(B))
        return [R([0] + b.list(), prec) for b in B]


    #########################################################################
    #
    #  Decomposition of spaces
    #
    ##########################################################################

    def factorization(self):
        """
        Returns a list of pairs $(S,e)$ where $S$ is simple spaces of
        modular symbols and self is isomorphic to the direct sum of
        the $S^e$ as a module over the \emph{anemic} Hecke algebra
        adjoin the star involution.

        ASSUMPTION: self is a module over the anemic Hecke algebra.
        """
        try:
            return self._factorization
        except AttributeError:
            raise NotImplementedError

    def hecke_module_of_level(self, level):
        r"""
        See the documentation for \code{self.modular_symbols_of_level(level)}.
        """
        return self.modular_symbols_of_level(Integer(level))

    def sign(self):
        """
        Returns the sign of self.

        For efficiency reasons, it is often useful to compute in the
        (largest) quotient of modular symbols where the * involution
        acts as +1, or where it acts as -1.


        INPUT:
           ModularSymbols self -- arbitrary space of modular symbols.

        OUTPUT:
           int -- the sign of self, either -1, 0, or 1.
                  -1 -- if this is factor of quotient where * acts as -1,
                  +1 -- if this is factor of quotient where * acts as +1,
                   0 -- if this is full space of modular symbols (no quotient).

        EXAMPLES:
            sage: m = ModularSymbols(33)
            sage: m.rank()
            9
            sage: m.sign()
            0
            sage: m = ModularSymbols(33, sign=0)
            sage: m.sign()
            0
            sage: m.rank()
            9
            sage: m = ModularSymbols(33, sign=-1)
            sage: m.sign()
            -1
            sage: m.rank()
            3
        """
        return self.__sign

    def simple_factors(self):
        """
        Returns a list modular symbols spaces $S$ where $S$ is simple
        spaces of modular symbols (for the anemic Hecke algebra) and
        self is isomorphic to the direct sum of the $S$ with some
        multiplicities, as a module over the \emph{anemic} Hecke
        algebra.  For the multiplicities use factorization() instead.

        ASSUMPTION: self is a module over the anemic Hecke algebra.

        EXAMPLES:
            sage: ModularSymbols(1,100,sign=-1).simple_factors()
            [Modular Symbols subspace of dimension 8 of Modular Symbols space of dimension 8 for Gamma_0(1) of weight 100 with sign -1 over Rational Field]
            sage: ModularSymbols(1,16,0,GF(5)).simple_factors()
            [Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5,
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 3 for Gamma_0(1) of weight 16 with sign 0 over Finite Field of size 5]
        """
        return [S for S,_ in self.factorization()]

    def star_eigenvalues(self):
        """
        Returns the eigenvalues of the star involution acting on self.

        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: D = M.decomposition()
            sage: M.star_eigenvalues()
            [1, -1]
            sage: D[0].star_eigenvalues()
            [1]
            sage: D[1].star_eigenvalues()
            [1, -1]
            sage: D[1].plus_submodule().star_eigenvalues()
            [1]
            sage: D[1].minus_submodule().star_eigenvalues()
            [-1]
        """
        try:
            return self.__star_eigenvalues
        except AttributeError:
            pass
        if self.sign() != 0:
            return [self.sign()]
        M = self.star_involution().matrix()
        R = self.base_ring()
        if M == R(1):
            self.__star_eigenvalues = [R(1)]
        elif M == R(-1):
            self.__star_eigenvalues = [R(-1)]
        else:
            self.__star_eigenvalues = [R(1), R(-1)]
        return self.__star_eigenvalues

    def star_decomposition(self):
        """
        """
        S = self.star_involution()
        return S.decomposition()

    def integral_structure(self):
        r"""
        Return the $\Z$-structure of this modular symbols spaces
        generated by all integral modular symbols.

        EXAMPLES:
            sage: M = ModularSymbols(11,4)
            sage: M.integral_structure()
            Free module of degree 6 and rank 6 over Integer Ring
            Echelon basis matrix:
            [    1     0     0     0     0     0]
            [    0  1/14   1/7  5/14   1/2 13/14]
            [    0     0   1/2     0     0   1/2]
            [    0     0     0     1     0     0]
            [    0     0     0     0     1     0]
            [    0     0     0     0     0     1]
            sage: M.cuspidal_submodule().integral_structure()
            Free module of degree 6 and rank 4 over Integer Ring
            Echelon basis matrix:
            [     0   1/14    1/7   5/14    1/2 -15/14]
            [     0      0    1/2      0      0   -1/2]
            [     0      0      0      1      0     -1]
            [     0      0      0      0      1     -1]
        """
        try:
            return self.__integral_structure
        except AttributeError:
            pass
        A = self.ambient_hecke_module()
        I = A.integral_structure()
        J = self.free_module().intersection(I)
        self.__integral_structure = J
        return J

    def intersection_number(self, M):
        """
        Given modular symbols spaces self and M in some common ambient
        space, returns the intersection number of these two spaces.
        This is the index in their saturation of the sum of their
        underlying integral structures.

        If self and M are of weight two and defined over QQ, and
        correspond to newforms f and g, then this number equals the
        order of the intersection of the modular abelian varieties
        attached to f and g.

        EXAMPLES:
            sage: m = ModularSymbols(389,2)
            sage: d = m.decomposition(2)
            sage: eis = d[0]
            sage: ell = d[1]
            sage: af = d[-1]
            sage: af.intersection_number(eis)
            97
            sage: af.intersection_number(ell)
            400
        """
        if not isinstance(M, ModularSymbolsSpace):
            raise TypeError, "M must be a modular symbols space"
        if M.ambient() != self.ambient():
            raise ValueError, "self and M must be in the same ambient space."
        A = self.integral_structure()
        B = M.integral_structure()
        return (A+B).index_in_saturation()

    def integral_basis(self):
        r"""
        Return a basis for the $\Z$-submodule of this modular symbols
        space spanned by the generators.

        Modular symbols spaces for congruence subgroups have a
        $\Z$-structure.  Computing this $\Z$-structure is expensive,
        so by default modular symbols spaces for congruence subgroups
        in \sage are defined over $\Q$.  This function returns a tuple
        of independent elements in this modular symbols space whose
        $\Z$-span is the corresponding space of modular symbols over
        $\Z$.

        EXAMPLES:
            sage: M = ModularSymbols(11)
            sage: M.basis()
            ((1,0), (1,8), (1,9))
            sage: M.integral_basis()
            ((1,0), (1,8), (1,9))
            sage: S = M.cuspidal_submodule()
            sage: S.basis()
            ((1,8), (1,9))
            sage: S.integral_basis()
            ((1,8), (1,9))

            sage: M = ModularSymbols(13,4)
            sage: M.basis()
            ([X^2,(0,1)], [X^2,(1,4)], [X^2,(1,5)], [X^2,(1,7)], [X^2,(1,9)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
            sage: M.integral_basis()
            ([X^2,(0,1)], 1/28*[X^2,(1,4)] + 2/7*[X^2,(1,5)] + 3/28*[X^2,(1,7)] + 11/14*[X^2,(1,9)] + 2/7*[X^2,(1,10)] + 11/28*[X^2,(1,11)] + 3/28*[X^2,(1,12)], [X^2,(1,5)], 1/2*[X^2,(1,7)] + 1/2*[X^2,(1,9)], [X^2,(1,9)], [X^2,(1,10)], [X^2,(1,11)], [X^2,(1,12)])
            sage: S = M.cuspidal_submodule()
            sage: S.basis()
            ([X^2,(1,4)] - [X^2,(1,12)], [X^2,(1,5)] - [X^2,(1,12)], [X^2,(1,7)] - [X^2,(1,12)], [X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,10)] - [X^2,(1,12)], [X^2,(1,11)] - [X^2,(1,12)])
            sage: S.integral_basis()
            (1/28*[X^2,(1,4)] + 2/7*[X^2,(1,5)] + 3/28*[X^2,(1,7)] + 11/14*[X^2,(1,9)] + 2/7*[X^2,(1,10)] + 11/28*[X^2,(1,11)] - 53/28*[X^2,(1,12)], [X^2,(1,5)] - [X^2,(1,12)], 1/2*[X^2,(1,7)] + 1/2*[X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,9)] - [X^2,(1,12)], [X^2,(1,10)] - [X^2,(1,12)], [X^2,(1,11)] - [X^2,(1,12)])

        This function currently raises a NotImplementedError on
        modular symbols spaces with character of order bigger than $2$:

        EXAMPLES:
            sage: M = ModularSymbols(DirichletGroup(13).0^2, 2); M
            Modular Symbols space of dimension 4 and level 13, weight 2, character [zeta6], sign 0, over Cyclotomic Field of order 6 and degree 2
            sage: M.basis()
            ((1,0), (1,5), (1,10), (1,11))
            sage: M.integral_basis()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        try:
            return self.__integral_basis
        except AttributeError:
            pass
        B = self.integral_structure().basis()
        self.__integral_basis = tuple([self(b) for b in B])
        return self.__integral_basis


    def integral_hecke_matrix(self, n):
        r"""
        Return the matrix of the $n$th Hecke operator acting on the integral
        structure on self (as returned by \code{self.integral_structure()}.
        """
        n = int(n)
        try:
            return self.__integral_hecke_matrix[n]
        except AttributeError:
            self.__integral_hecke_matrix = {}
        except KeyError:
            pass
        #raise NotImplementedError, "code past this point is broken / not done"  # todo
        A = self.ambient_hecke_module()
        T = A.hecke_matrix(n)
        S = T.restrict(self.integral_structure()).change_ring(ZZ)
        self.__integral_hecke_matrix[n] = S
        return S

    def sturm_bound(self):
        r"""
        Returns the Sturm bound for this space of modular symbols.

        Type \code{sturm\_bound?} for more details.

        EXAMPLES:
            sage: ModularSymbols(11,2).sturm_bound()
            2
            sage: ModularSymbols(389,2).sturm_bound()
            65
            sage: ModularSymbols(1,12).sturm_bound()
            1
            sage: ModularSymbols(1,36).sturm_bound()
            3
        """
        # For Gamma_0(N), n = \frac{k}{12}[\SL_2(\Z):\Gamma_0(N)]
        try:
            return self.__sturm_bound
        except:
            self.__sturm_bound = dims.sturm_bound(self.level(), self.weight())
        return self.__sturm_bound


    def plus_submodule(self, compute_dual=True):
        """
        Return the subspace of self on which the star involution acts as +1.

        INPUT:
            compute_dual -- bool (default: True) also compute dual subspace.
                            This are useful for many algorithms.

        OUTPUT:
            subspace of modular symbols

        EXAMPLES:
            sage: ModularSymbols(17,2)
            Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
            sage: ModularSymbols(17,2).plus_submodule()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(17) of weight 2 with sign 0 over Rational Field
        """
        return self.sign_submodule(+1, compute_dual)

    def minus_submodule(self, compute_dual=True):
        """
        Return the subspace of self on which the star involution acts as -1.

        INPUT:
            compute_dual -- bool (default: True) also compute dual subspace.
                            This are useful for many algorithms.
        OUTPUT:
            subspace of modular symbols

        EXAMPLES:
            sage: ModularSymbols(14,4)
            Modular Symbols space of dimension 12 for Gamma_0(14) of weight 4 with sign 0 over Rational Field
            sage: ModularSymbols(14,4).minus_submodule()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 12 for Gamma_0(14) of weight 4 with sign 0 over Rational Field
        """
        return self.sign_submodule(-1, compute_dual)

    def _compute_sign_submodule(self, sign, compute_dual=True):
        A = self.ambient()
        S = A.sign_submodule(sign, compute_dual=compute_dual)
        V = S.free_module().intersection(self.free_module())
        if compute_dual:
            W = S.dual_free_module()
            Y = self.dual_free_module()
            D = W.intersection(Y)
            M = A.submodule(V, D, check=False)
        else:
            M = A.submodule(V, check=False)
        M._set_sign(sign)
        return M

    def _set_sign(self, sign):
        sign = int(sign)
        if not (sign in [-1,0,1]):
            raise ValueError, "sign (=%s) must be -1, 0, or 1"%sign
        self.__sign = sign

    def sign_submodule(self, sign, compute_dual=True):
        """
        Return the subspace of self that is fixed under the star involution.

        INPUT:
            sign -- int (either -1, 0 or +1)
            compute_dual -- bool (default: True) also compute dual subspace.
                            This are useful for many algorithms.
        OUTPUT:
            subspace of modular symbols

        EXAMPLES:
            sage: M = ModularSymbols(29,2)
            sage: M.sign_submodule(1)
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 5 for Gamma_0(29) of weight 2 with sign 0 over Rational Field
            sage: M.sign_submodule(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(29) of weight 2 with sign 0 over Rational Field
            sage: M.sign_submodule(-1).sign()
            -1
        """
        sign = int(sign)
        if not sign in [-1, 0, 1]:
            raise ValueError, "sign must be -1, 0 or 1"
        if self.sign() == sign:  # an easy case
            return self
        if self.sign() == -sign:  # another easy case
            return self.zero_submodule()
        if sign == 0:
            # if sign is zero then self.sign() isn't 0 because
            # of the above checks.
            raise ArithmeticError, "There is no sign 0 subspace of a space of modular symbols with nonzero sign."
        try:
            return self.__plus_submodule[(sign, compute_dual)]
        except AttributeError:
            self.__plus_submodule = {}
        except KeyError:
            pass
        P = self._compute_sign_submodule(sign, compute_dual)
        P.__star_eigenvalue = sign
        self.__plus_submodule[(sign,compute_dual)] = P
        return P

    def star_involution(self):
        """
        Return the star involution on self, which is induced by complex
        conjugation on modular symbols.
        """
        raise NotImplementedError

    def abelian_variety(self):
        """
        Return the corresponding abelian variety.

        INPUT:
            self -- modular symbols space of weight 2 for a congruence subgroup
                    such as Gamma0, Gamma1 or GammaH.

        EXAMPLES:
            sage: ModularSymbols(Gamma0(11)).cuspidal_submodule().abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: ModularSymbols(Gamma1(11)).cuspidal_submodule().abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: ModularSymbols(GammaH(11,[3])).cuspidal_submodule().abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1

        The abelian variety command only works on cuspidal modular symbols spaces:
            sage: M = ModularSymbols(37)
            sage: M[0].abelian_variety()
            Traceback (most recent call last):
            ...
            ValueError: self must be cuspidal
            sage: M[1].abelian_variety()
            Abelian subvariety of dimension 1 of J0(37)
            sage: M[2].abelian_variety()
            Abelian subvariety of dimension 1 of J0(37)
        """
        try:
            return self.__modular_abelian_variety
        except AttributeError:
            if not self.is_cuspidal():
                raise ValueError, "self must be cuspidal"
            from sage.modular.abvar.abvar import ModularAbelianVariety_modsym
            A = ModularAbelianVariety_modsym(self, check=False)
            self.__modular_abelian_variety = A
            return A

    def rational_period_mapping(self):
        r"""
        Return the rational period mapping associated to self.

        This is a homomorphism to a vector space whose kernel is the
        same as the kernel of the period mapping associated to self.
        For this to exist, self must be Hecke equivariant.

        Use \code{integral_period_mapping} to obtain a homomorphism to
        a $\ZZ$-module, normalized so the image of integral modular
        symbols is exactly $\ZZ^n$.

        EXAMPLES:
            sage: M = ModularSymbols(37)
            sage: A = M[1]; A
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r = A.rational_period_mapping(); r
            Rational period mapping associated to Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r(M.0)
            (0, 0)
            sage: r(M.1)
            (1, 0)
            sage: r.matrix()
            [ 0  0]
            [ 1  0]
            [ 0  1]
            [-1 -1]
            [ 0  0]
            sage: r.domain()
            Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: r.codomain()
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        try:
            return self.__rational_period_mapping
        except AttributeError:
            pass
        V = self.dual_free_module()
        # the rational period mapping is just dotting with
        # each element of V.
        A = V.basis_matrix().transpose()
        A.set_immutable()
        R = RationalPeriodMapping(self, A)
        self.__rational_period_mapping = R
        return R


    def integral_period_mapping(self):
        r"""
        Return the integral period mapping associated to self.

        This is a homomorphism to a vector space whose kernel is the
        same as the kernel of the period mapping associated to self,
        normalized so the image of integral modular symbols is
        exactly $\ZZ^n$.

        EXAMPLES:
            sage: m = ModularSymbols(23).cuspidal_submodule()
            sage: i = m.integral_period_mapping()
            sage: i
            Integral period mapping associated to Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: i.matrix()
            [-1/11  1/11     0  3/11]
            [    1     0     0     0]
            [    0     1     0     0]
            [    0     0     1     0]
            [    0     0     0     1]
            sage: [i(b) for b in m.integral_structure().basis()]
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
            sage: [i(b) for b in m.ambient_module().basis()]
            [(-1/11, 1/11, 0, 3/11),
             (1, 0, 0, 0),
             (0, 1, 0, 0),
             (0, 0, 1, 0),
             (0, 0, 0, 1)]

        We compute the image of the winding element:
            sage: m = ModularSymbols(37,sign=1)
            sage: a = m[1]
            sage: f = a.integral_period_mapping()
            sage: e = m([0,oo])
            sage: f(e)
            (-2/3)

        The input space must be cuspidal:
            sage: m = ModularSymbols(37,2,sign=1)
            sage: m.integral_period_mapping()
            Traceback (most recent call last):
            ...
            ValueError: integral mapping only defined for cuspidal spaces
        """
        try:
            return self.__integral_period_mapping
        except AttributeError:
            pass
        if self.base_ring() != QQ:
            raise ValueError, "integral mapping only defined for spaces over QQ"
        if not self.is_cuspidal():
            raise ValueError, "integral mapping only defined for cuspidal spaces"
        D = self.dual_free_module().basis_matrix().transpose()
        I = self.ambient_module().cuspidal_submodule().integral_structure().basis_matrix()
        # image of cuspidal integral submodule
        C = I * D
        if not C.is_one():
            if not C.is_square():
                C = (ZZ**C.ncols()).span(C.rows()).basis_matrix()
            D = D * C**(-1)
        D.set_immutable()
        R = IntegralPeriodMapping(self, D)
        self.__integral_period_mapping = R
        return R

    def modular_symbols_of_sign(self, sign, bound=None):
        """
        Returns a space of modular symbols with the same defining
        properties (weight, level, etc.) and Hecke eigenvalues as this
        space except with given sign.

        INPUT:
            self -- a cuspidal space of modular symbols
            sign -- an integer, one of -1, 0, or 1
            bound -- integer (default: None); if specified only use
                 Hecke operators up to the given bound.

        EXAMPLES:
            sage: S = ModularSymbols(Gamma0(11),2,sign=0).cuspidal_subspace()
            sage: S
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

            sage: S = ModularSymbols(43,2,sign=1)[2]; S
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 4 for Gamma_0(43) of weight 2 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(43) of weight 2 with sign -1 over Rational Field

            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(43) of weight 2 with sign 0 over Rational Field


            sage: S = ModularSymbols(389,sign=1)[3]; S
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 33 for Gamma_0(389) of weight 2 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 3 of Modular Symbols space of dimension 32 for Gamma_0(389) of weight 2 with sign -1 over Rational Field
            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 65 for Gamma_0(389) of weight 2 with sign 0 over Rational Field

            sage: S = ModularSymbols(23,sign=1,weight=4)[2]; S
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 7 for Gamma_0(23) of weight 4 with sign 1 over Rational Field
            sage: S.modular_symbols_of_sign(1) is S
            True
            sage: S.modular_symbols_of_sign(0)
            Modular Symbols subspace of dimension 8 of Modular Symbols space of dimension 12 for Gamma_0(23) of weight 4 with sign 0 over Rational Field
            sage: S.modular_symbols_of_sign(-1)
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(23) of weight 4 with sign -1 over Rational Field
        """
        if sign == self.sign():
            return self
        if not self.is_cuspidal():
            raise ValueError, "self must be cuspidal for modular symbols space with given sign to be defined."
        d = self.dimension()
        if d == 0:
            return self
        if sign != 0:
            if self.sign() == 0:
                d = d//2
        elif sign == 0:      # self has nonzero sign
            d = 2*d
        B = self.ambient_module().modular_symbols_of_sign(sign)
        p = 2
        if bound is None:
            bound = self.hecke_bound()
        while B.dimension() > d and p <= bound:
            while self.level() % p == 0:
                p = arith.next_prime(p)
            f = self.hecke_polynomial(p)
            g = misc.prod(g for g,_ in f.factor())   # square free part
            t = B.hecke_operator(p)
            s = g(t)
            B = s.kernel()
            p = arith.next_prime(p)
        return B


class PeriodMapping(SageObject):
    def __init__(self, modsym, A):
        self.__modsym = modsym
        self.__domain = modsym.ambient_module()
        self.__A = A
        A.set_immutable()

    def modular_symbols_space(self):
        return self.__modsym

    def __call__(self, x):
        if is_FreeModuleElement(x):
            v = x
        else:
            v = self.__domain(x).element()
        return v*self.__A

    def matrix(self):
        return self.__A

    def domain(self):
        return self.__domain

    def codomain(self):
        return self.__A.row_module()

class RationalPeriodMapping(PeriodMapping):
    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: ModularSymbols(40,2).rational_period_mapping()._repr_()
            'Rational period mapping associated to Modular Symbols space of dimension 13 for Gamma_0(40) of weight 2 with sign 0 over Rational Field'
        """
        return "Rational period mapping associated to %s"%self.modular_symbols_space()


class IntegralPeriodMapping(PeriodMapping):
    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: ModularSymbols(40,2).cuspidal_submodule().integral_period_mapping()._repr_()
            'Integral period mapping associated to Modular Symbols subspace of dimension 6 of Modular Symbols space of dimension 13 for Gamma_0(40) of weight 2 with sign 0 over Rational Field'
        """
        return "Integral period mapping associated to %s"%self.modular_symbols_space()
