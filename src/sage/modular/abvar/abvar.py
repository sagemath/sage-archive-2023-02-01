"""
Base class for modular abelian varieties

AUTHOR:
    -- William Stein (2007-03)

TESTS:
    sage: A = J0(33)
    sage: D = A.decomposition(); D
    [
    Abelian variety factor of dimension 1 of J0(33),
    Abelian variety factor of dimension 2 of J0(33)
    ]
    sage: loads(dumps(D)) == D
    True
    sage: loads(dumps(A)) == A
    True
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.all        import ModularAbelianVarieties
from sage.structure.sequence    import Sequence
from sage.structure.parent_base import ParentWithBase
from hecke_operator             import HeckeOperator
from torsion_subgroup           import TorsionSubgroup
from finite_subgroup            import FiniteSubgroup_gens, FiniteSubgroup, FiniteSubgroupElement
from cuspidal_subgroup          import CuspidalSubgroup, RationalCuspidalSubgroup
from sage.rings.all             import ZZ, QQ, QQbar, is_Ring, LCM, divisors
from sage.modules.all           import is_FreeModule
from sage.modular.congroup      import is_CongruenceSubgroup, is_Gamma0, is_Gamma1, is_GammaH
from sage.modular.modsym.all    import ModularSymbols
from sage.modular.modsym.space  import ModularSymbolsSpace
from sage.matrix.all            import matrix
from sage.groups.all            import AbelianGroup

import homology
import homspace
import lseries

def is_ModularAbelianVariety(x):
    """
    Return True if x is a modular abelian variety.

    INPUT:
        x -- object

    EXAMPLES:
        sage: is_ModularAbelianVariety(5)
        False
        sage: is_ModularAbelianVariety(J0(37))
        True

    Returning True is a statement about the data type not
    whether or not some abelian variety is modular:
        sage: is_ModularAbelianVariety(EllipticCurve('37a'))
        False
    """
    return isinstance(x, ModularAbelianVariety_abstract)


class ModularAbelianVariety_abstract(ParentWithBase):
    def __init__(self, base_field, check=True):
        if check and not is_Ring(base_field) and base_field.is_field():
            raise TypeError, "base_field must be a field"
        ParentWithBase.__init__(self, base_field)

    # groups() and lattice() *must* be defined by every derived class!!!!
    def groups(self):
        raise NotImplementedError

    def lattice(self):
        raise NotImplementedError

    def base_field(self):
        r"""
        Synonym for \code{self.base_ring()}.

        EXAMPLES:
            sage: J0(11).base_field()
            Rational Field
        """
        return self.base_ring()

    def __contains__(self, x):
        """
            sage: J = J0(67); G = (J[0] + J[1]).intersection(J[1] + J[2])
            sage: G[0]
            Finite subgroup with invariants [5, 10] over QQbar of Abelian variety factor of dimension 3 of J0(67)
            sage: a = G[0].0; a
            [(1/10, 1/10, 3/10, 1/2, 1/5, 4/5)]
            sage: a in J[0]
            False
            sage: a in (J[0]+J[1])
            True
            sage: a in (J[1]+J[2])
            True
            sage: C = G[1]   # abelian variety in kernel
            sage: G[0].0
            [(1/10, 1/10, 3/10, 1/2, 1/5, 4/5)]
            sage: 5*G[0].0
            [(1/2, 1/2, 3/2, 5/2, 1, 4)]
            sage: 5*G[0].0 in C
            True
        """
        if not isinstance(x, FiniteSubgroupElement):
            return False
        if x.parent().abelian_variety().groups() != self.groups():
            return False
        v = x.ambient_element()
        n = v.denominator()
        nLambda = self.ambient_variety().lattice().scale(n)
        return n*v in self.lattice() + nLambda

    def __cmp__(self, other):
        if not isinstance(other, ModularAbelianVariety_abstract):
            return cmp(type(self), type(other))
        if self is other:
            return 0
        c = cmp(self.groups(), other.groups())
        if c: return c
        return cmp(self.lattice(), other.lattice())

    def _repr_(self):
        """
        Return string representation of this modular abelian variety.

        This is just the generic base class, so it's unlikely to be called in practice.

        EXAMPLES:
            sage: A = J0(23)
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety_abstract._repr_(A)
            'Abelian variety J0(23)'
        """
        if self.is_ambient():
            return 'Abelian variety %s'%self._ambient_repr()
        return "Abelian variety factor of dimension %s of %s%s"%(
            self.dimension(),
            self._ambient_repr(),
            '' if self.base_field() == QQ else ' over %s'%self.base_field())

    def _Hom_(self, B, cat=None):
        """
        INPUT:
            B -- modular abelian varieties
            cat -- category

        EXAMPLES:
            sage: J0(37)._Hom_(J1(37))
            Space of homomorphisms from Jacobian of the modular curve associated to the congruence subgroup Gamma0(37) to Jacobian of the modular curve associated to the congruence subgroup Gamma1(37)
        """
        if cat is None:
            K = self.base_field(); L = B.base_field()
            if K == L:
                F = K
            elif K == QQbar or L == QQbar:
                F = QQbar
            else:
                # TODO -- improve this
                raise ValueError, "please specify a category"
            cat = ModularAbelianVarieties(F)
        return homspace.Homspace(self, B, cat)

    def intersection(self, other):
        """
        Returns the intersection of self and other inside the ambient
        Jacobian product.  self and other must be abelian subvarieties
        of the ambient Jacobian product.

        This function returns a finite group $H$ and an abelian variety $A$
        such that the intersection is $H + A$, which need not be connected.

        INPUT:
            other -- a modular abelian variety

        OUTPUT:
            G -- finite subgroup of self
            A -- abelian variety

        EXAMPLES:
        We intersect some abelian varieties with finite intersection.
            sage: J = J0(37)
            sage: J[0].intersection(J[1])
            (Finite subgroup with invariants [2, 2] over QQ of Abelian variety factor of dimension 1 of J0(37), Abelian variety factor of dimension 0 of J0(37))

            sage: J = J0(65)
            sage: D = J.decomposition(); D
            [
            Abelian variety factor of dimension 1 of J0(65),
            Abelian variety factor of dimension 2 of J0(65),
            Abelian variety factor of dimension 2 of J0(65)
            ]
            sage: A = D[0] + D[1]; B = D[1] + D[2]
            sage: A.intersection(B)
            (Finite subgroup with invariants [2] over QQbar of Abelian variety factor of dimension 3 of J0(65),
             Abelian variety factor of dimension 2 of J0(65))
            sage: D[0].intersection(D[2])
            (Finite subgroup with invariants [2] over QQ of Abelian variety factor of dimension 1 of J0(65), Abelian variety factor of dimension 0 of J0(65))

            sage: J = J0(33)
            sage: J[0].intersection(J[1])
            (Finite subgroup with invariants [3, 3] over QQ of Abelian variety factor of dimension 1 of J0(33), Abelian variety factor of dimension 0 of J0(33))

        Next we intersect two abelian varieties with non-finite intersection:
            sage: J = J0(67); D = J.decomposition(); D
            [
            Abelian variety factor of dimension 1 of J0(67),
            Abelian variety factor of dimension 2 of J0(67),
            Abelian variety factor of dimension 2 of J0(67)
            ]
            sage: (D[0] + D[1]).intersection(D[1] + D[2])
            (Finite subgroup with invariants [5, 10] over QQbar of Abelian variety factor of dimension 3 of J0(67), Abelian variety factor of dimension 2 of J0(67))

        """
        if not is_ModularAbelianVariety(other):
            raise TypeError, "other must be a modular abelian variety"
        if self.groups() != other.groups():
            raise ValueError, "incompatible ambient Jacobians"
        if not self.is_subvariety_of_ambient_jacobian() or not other.is_subvariety_of_ambient_jacobian():
            raise ValueError, "self and other must be subvarieties of the ambient product Jacobian"

        # 1. find the abvar part
        L = self.lattice().intersection(other.lattice())
        if L.dimension() > 0:
            L = L.intersection(self._ambient_lattice())
            finitegroup_base_field = QQbar
        else:
            if self.base_field() == other.base_field():
                finitegroup_base_field = self.base_field()
            else:
                # todo -- compositum
                finitegroup_base_field = QQbar

        A = ModularAbelianVariety(self.groups(), L, self.base_field(), check=False)

        # 2. find the finite component group, as a finite subgroup of self.
        L = self.lattice().basis_matrix()
        M = other.lattice().basis_matrix()

        if False:
            LM = L.stack(M)
            P = LM.transpose().pivots()
            V = (QQ**L.ncols()).span_of_basis([LM.row(p) for p in P])
            S = (self.lattice() + other.lattice()).saturation().scale(5)
            n = self.lattice().rank()
            gens = [V.coordinates(w)[:n] for w in S.basis()]

        if True:
            B = L.stack(M)
            piv0 = B.pivots()
            S = B.change_ring(ZZ).saturation().matrix_from_columns(piv0).transpose()

            B = B.transpose()
            piv1 = B.pivots()
            B = B.matrix_from_rows_and_columns(piv0, piv1)

            # Write each column of S in terms of the columns of B.
            X = B.solve_right(S.change_ring(QQ))

            # Finally, project to the L factor
            gens = X.matrix_from_rows(range(L.nrows())).columns()


        G = self.finite_subgroup(gens, base_field=finitegroup_base_field)

        return G, A


    def __add__(self, other):
        """
        Returns the sum of the images of self and other inside the
        ambient Jacobian product.   self and other must be abelian subvarieties
        of the ambient Jacobian product.

        EXAMPLES:

        """
        if not is_ModularAbelianVariety(other):
            raise TypeError, "other must be a modular abelian variety"
        if self.groups() != other.groups():
            raise ValueError, "incompatible ambient Jacobians"
        if not self.is_subvariety_of_ambient_jacobian() or not other.is_subvariety_of_ambient_jacobian():
            raise ValueError, "self and other must be subvarieties of the ambient product Jacobian"
        L = self.lattice() + other.lattice()
        M = L.intersection(self._ambient_lattice())
        return ModularAbelianVariety(self.groups(), M, self.base_field(), check=False)

    def direct_product(self, other):
        """
        Compute the direct product of self and other.
        """
        return self * other

    def __mul__(self, other):
        """
        Compute the direct product of self and other.

        EXAMPLES:
        Some modular Jacobians:
            sage: J0(11) * J0(33)
            Abelian variety J0(11) x J0(33)
            sage: J0(11) * J0(33) * J0(11)
            Abelian variety J0(11) x J0(33) x J0(11)

        We multiply some factors of $J_0(65)$:
            sage: d = J0(65).decomposition()
            sage: d[0] * d[1] * J0(11)
            Abelian variety factor of dimension 4 of J0(65) x J0(65) x J0(11)
        """
        if not is_ModularAbelianVariety(other):
            raise TypeError, "other must be a modular abelian variety"
        if other.base_ring() != self.base_ring():
            raise TypeError, "self and other must have the same base ring"
        groups = tuple(list(self.groups()) + list(other.groups()))
        lattice = self.lattice().direct_sum(other.lattice())
        base_field = self.base_ring()
        return ModularAbelianVariety(groups, lattice, base_field, check=False)

    def __div__(self, other):
        """
        Compute the quotient of self and other, where other is either
        an abelian subvariety of self or a finite subgroup of self.

        INPUT:
            other -- a finite subgroup or subvariety

        EXAMPLES:
            sage: J = J0(67); G = (J[0] + J[1]).intersection(J[1] + J[2])
            sage: Q = J/G[0]; Q
            Abelian variety factor of dimension 5 of J0(67) over Algebraic Field
            sage: Q.base_field()
            Algebraic Field
            sage: Q.lattice()
            Free module of degree 10 and rank 10 over Integer Ring
            Echelon basis matrix:
            [1/10 1/10 3/10  1/2    0    0    0 3/10    0  1/2]
            [   0  1/5  4/5  4/5    0    0    0    0    0  3/5]
            ...
        """
        if isinstance(other, FiniteSubgroup):
            if other.abelian_variety() != self:
                other = self.finite_subgroup(other)
            return self._quotient_by_finite_subgroup(other)
        elif isinstance(other, ModularAbelianVariety_abstract):
            if other.is_subvariety(self):
                return self._quotient_by_abelian_subvariety(other)
        else:
            raise TypeError, "other must be a subgroup or abelian subvariety"

    def _quotient_by_finite_subgroup(self, G):
        if G.order() == 1:
            return self
        return ModularAbelianVariety(self.groups(), self.lattice() + G.lattice(), G.base_field())

    def is_subvariety_of_ambient_jacobian(self):
        try:
            return self.__is_sub_ambient
        except AttributeError:
            self.__is_sub_ambient = (self.lattice().denominator() == 1)
            return self.__is_sub_ambient

    def ambient_variety(self):
        """
        Return the ambient modular abelian variety that contains this
        abelian variety.  The ambient variety is always a product of
        Jacobians of modular curves.
        """
        try:
            return self.__ambient_variety
        except AttributeError:
            A = ModularAbelianVariety(self.groups(), ZZ**(2*self._ambient_dimension()),
                                     self.base_field(), check=False)
            self.__ambient_variety = A
            return A

    def is_ambient(self):
        try:
            return self.__is_ambient
        except AttributeError:
            pass
        L = self.lattice()
        self.__is_ambient = (self.lattice() == ZZ**L.degree())
        return self.__is_ambient

    def dimension(self):
        """
        Return the dimension of this abelian variety.

        EXAMPLES:
            sage: A = J0(23)
            sage: A.dimension()
            2
        """
        return self.lattice().rank() // 2

    def degree(self):
        """
        Return the degree of this abelian variety, which
        is the dimension of the ambient Jacobian product.

        EXAMPLES:
            sage: A = J0(23)
            sage: A.dimension()
            2
        """
        return self._ambient_dimension()

    def endomorphism_ring(self):
        try:
            return self.__endomorphism_ring
        except AttributeError:
            pass

        self.__endomorphism_ring = homspace.EndomorphismSubring(self)
        return self.__endomorphism_ring

    def is_hecke_stable(self):
        """
        Return True if self is stable under the Hecke operators of
        its ambient Jacobian.
        """
        try:
            return self._is_hecke_stable
        except AttributeError:
            pass

        b = self.modular_symbols().sturm_bound()
        J = self.ambient_variety()
        L = self.lattice()
        B = self.lattice().basis()

        for n in range(1,b+1):
            Tn_matrix = J.hecke_operator(n).matrix()
            for v in B:
                if not (v*Tn_matrix in L):
                    self._is_hecke_stable = False
                    return False

        self._is_hecke_stable = True
        return True

    def is_subvariety(self, other):
        """
        Return True if self is a subvariety of other as they sit in a
        common ambient modular Jacobian.  In particular, this function
        will only return True if self and other have exactly the same
        ambient Jacobians.

        EXAMPLES:
            sage: J = J0(37); J
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: A = J[0]; A
            Abelian variety factor of dimension 1 of J0(37)
            sage: A.is_subvariety(A)
            True
            sage: A.is_subvariety(J)
            True
        """
        if not is_ModularAbelianVariety(other):
            return False
        if self is other:
            return True
        if self.groups() != other.groups():
            return False
        L = self.lattice()
        M = other.lattice()
        # self is an abelian subvariety of other if and only if
        #   1. L is a subset of M (so the abelian subvarieties of the ambient J are equal), and
        #   2. L is relatively saturated in M, i.e., M/L is torsion free.
        if not L.is_submodule(M):
            return False
        # To determine if L is relatively staturated we compute the intersection
        # of M with (L tensor Q) and see if that equals L.
        return L.change_ring(QQ).intersection(M) == L

    def change_ring(self, R):
        """
        Change the base ring of this modular abelian variety.

        EXAMPLES:
            sage: A = J0(23)
            sage: A.change_ring(QQ)
            Abelian variety J0(23)
        """
        return ModularAbelianVariety(self.groups(), self.lattice(), R, check=False)

    def category(self):
        """
        Return the category of modular abelian varieties that contains
        this modular abelian variety.

        EXAMPLES:
            sage: J0(23).category()
            Category of modular abelian varieties over Rational Field
        """
        try:
            return self.__category
        except AttributeError:
            C = ModularAbelianVarieties(self.base_ring())
            self.__category = C
            return C

    def level(self):
        """
        TODO: Rewrite
        Return the level of this modular abelian variety, which is an integer
        N (usually minimal) such that this modular abelian variety is a quotient
        of $J_1(N)$.

        EXAMPLES:
            sage: J1(5077).level()
            5077
            sage: JH(389,[4]).level()
            389
        """
        try:
            return self.__level
        except AttributeError:
            self.__level = LCM([G.level() for G in self.groups()])
            return self.__level

    ###############################################################################
    # Properties of the ambient product of Jacobians
    ###############################################################################
    def _ambient_repr(self):
        v = []
        for G in self.groups():
            if is_Gamma0(G):
                v.append('J0(%s)'%G.level())
            elif is_Gamma1(G):
                v.append('J1(%s)'%G.level())
            elif is_GammaH(G):
                v.append('JH(%s,%s)'%(G.level(), G._generators_for_H()))
        return ' x '.join(v)

    def _ambient_lattice(self):
        try:
            return self.__ambient_lattice
        except AttributeError:
            self.__ambient_lattice = ZZ**(2*self.degree())
            return self.__ambient_lattice

    def _ambient_modular_symbols_spaces(self):
        try:
            return self.__ambient_modular_symbols_spaces
        except AttributeError:
            self.__ambient_modular_symbols_spaces = tuple([ModularSymbols(G).cuspidal_subspace() for G in self.groups()])
            return self.__ambient_modular_symbols_spaces

    def _ambient_dimension(self):
        try:
            return self.__ambient_dimension
        except AttributeError:
            self.__ambient_dimension = sum([S.dimension() for S in self._ambient_modular_symbols_spaces()]) // 2
            return self.__ambient_dimension

    def _ambient_hecke_matrix_on_modular_symbols(self, n):
        r"""
        Return block direct sum of the matrix of the Hecke operator $T_n$ acting
        on each of the ambient modular symbols spaces.

        INPUT:
            n -- an integer $\geq 1$.

        OUTPUT:
            a matrix
        """
        try:
            return self.__ambient_hecke_matrix_on_modular_symbols[n]
        except AttributeError:
            self.__ambient_hecke_matrix_on_modular_symbols = {}
        except KeyError:
            pass
        M = self._ambient_modular_symbols_spaces()
        if len(M) == 0:
            return matrix(QQ,0)
        T = M[0].hecke_matrix(n)
        for i in range(1,len(M)):
            T = T.block_sum(M[i].hecke_matrix(n))
        self.__ambient_hecke_matrix_on_modular_symbols[n] = T
        return T

    ###############################################################################
    # Rational and Integral Homology
    ###############################################################################
    def _rational_homology_space(self):
        """
        Return the rational homology of this modular abelian variety.

        EXAMPLES:
            sage: J = J0(11)
            sage: J._rational_homology_space()
            Vector space of dimension 2 over Rational Field

        The result is cached:
            sage: J._rational_homology_space() is J._rational_homology_space()
            True
        """
        try:
            return self.__rational_homology_space
        except AttributeError:
            HQ = self.rational_homology().free_module()
            self.__rational_homology_space = HQ
            return HQ

    def homology(self, base_ring=ZZ):
        """
        Return the homology of this modular abelian variety.

        WARNING: For efficiency reasons the basis of the integral
        homology need not be the same as the basis for the rational
        homology.

        EXAMPLES:
            sage: J0(389).homology(GF(7))
            Homology with coefficients in Finite Field of size 7 of Jacobian of the modular curve associated to the congruence subgroup Gamma0(389)
            sage: J0(389).homology(QQ)
            Rational Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(389)
            sage: J0(389).homology(ZZ)
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(389)
        """
        try:
            return self._homology[base_ring]
        except AttributeError:
            self._homology = {}
        except KeyError:
            pass
        if base_ring == ZZ:
            H = homology.IntegralHomology(self)
        elif base_ring == QQ:
            H = homology.RationalHomology(self)
        else:
            H = homology.Homology_over_base(self, base_ring)
        self._homology[base_ring] = H
        return H

    def integral_homology(self):
        """
        Return the integral homology of this modular abelian variety.

        EXAMPLES:
            sage: H = J0(43).integral_homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(43)
            sage: H.rank()
            6
            sage: H = J1(17).integral_homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma1(17)
            sage: H.rank()
            10

        If you just ask for the rank of the homology, no serious
        calculations are done, so the following is fast:
            sage: H = J0(50000).integral_homology(); H
            Integral Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(50000)
            sage: H.rank()
            14702
        """
        return self.homology(ZZ)

    def rational_homology(self):
        """
        Return the rational homology of this modular abelian variety.

        EXAMPLES:
            sage: H = J0(37).rational_homology(); H
            Rational Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: H.rank()
            4
            sage: H.base_ring()
            Rational Field
            sage: H = J1(17).rational_homology(); H
            Rational Homology of Jacobian of the modular curve associated to the congruence subgroup Gamma1(17)
            sage: H.rank()
            10
            sage: H.base_ring()
            Rational Field
        """
        return self.homology(QQ)

    ###############################################################################
    # L-series
    ###############################################################################
    def lseries(self):
        """
        Return the complex $L$-series of this modular abelian variety.

        EXAMPLES:
            sage: A = J0(37)
            sage: A.lseries()
            Complex L-series attached to Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
        """
        try:
            return self.__lseries
        except AttributeError:
            pass
        self.__lseries = lseries.Lseries_complex(self)
        return self.__lseries

    def padic_lseries(self, p):
        """
        Return the $p$-adic $L$-series of this modular abelian variety.

        EXAMPLES:
            sage: A = J0(37)
            sage: A.padic_lseries(7)
            7-adic L-series attached to Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
        """
        p = int(p)
        try:
            return self.__lseries_padic[p]
        except AttributeError:
            self.__lseries_padic = {}
        except KeyError:
            pass
        self.__lseries_padic[p] = lseries.Lseries_padic(self, p)
        return self.__lseries_padic[p]

    ###############################################################################
    # Hecke Operators
    ###############################################################################
    def hecke_operator(self, n):
        """
        Return the $n$-th Hecke operator on the modular abelian
        variety, if this makes sense [[ellaborate]].  Otherwise raise
        a ValueError.

        EXAMPLES:
        We compute $T_2$ on $J_0(37)$.
            sage: t2 = J0(37).hecke_operator(2); t2
            Hecke operator T_2 on Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: t2.charpoly().factor()
            x^2 * (x + 2)^2
            sage: t2.index()
            2

        Note that there is no matrix associated to Hecke operators on
        modular abelian varieties.  For a matrix, instead consider, e.g.,
        the Hecke operator on integral or rational homology.
            sage: t2.action_on_homology().matrix()
            [-1  1  1 -1]
            [ 1 -1  1  0]
            [ 0  0 -2  1]
            [ 0  0  0  0]
        """
        try:
            return self._hecke_operator[n]
        except AttributeError:
            self._hecke_operator = {}
        except KeyError:
            pass
        Tn = HeckeOperator(self, n)
        self._hecke_operator[n] = Tn
        return Tn

    def hecke_polynomial(self, n, var='x'):
        return self.hecke_operator(n).charpoly(var='x')

    def _integral_hecke_matrix(self, n):
        """
        Return the matrix of the Hecke operator $T_n$ acting on the
        integral homology of this modular abelian variety, if the
        modular abelian variety is stable under $T_n$.  Otherwise,
        raise an ArithmeticError.

        EXAMPLES:
            sage: A = J0(23)
            sage: t = A._integral_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Integer Ring
        """
        A = self._ambient_hecke_matrix_on_modular_symbols(n)
        return A.restrict(self.lattice())

    def _rational_hecke_matrix(self, n):
        r"""
        Return the matrix of the Hecke operator $T_n$ acting on the
        rational homology $H_1(A,\Q)$ of this modular abelian variety,
        if this action is defined.  Otherwise, raise an
        ArithmeticError.

        EXAMPLES:
            sage: A = J0(23)
            sage: t = A._rational_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field
        """
        return self._integral_hecke_matrix(n)

    ###############################################################################
    # Finite Subgroups
    ###############################################################################
    def torsion_subgroup(self):
        """
        EXAMPLES:
            sage: J = J0(33)
            sage: A = J.new_quotient()
            sage: A
            Abelian variety factor of dimension 1 of J0(33)
            sage: t = A.torsion_subgroup()
            sage: t.multiple_of_order()
            4
            sage: t.divisor_of_order()
            4
            sage: t.order()
            4
            sage: t.gens()
            [[(1/2, 0)], [(0, 1/2)]]
            sage: t
            Torsion subgroup of Abelian variety factor of dimension 1 of J0(33)
        """
        try:
            return self._torsion_subgroup
        except AttributeError:
            T = TorsionSubgroup(self)
            self._torsion_subgroup = T
            return T

    def cuspidal_subgroup(self):
        """
        Return the cuspidal subgroup of this modular abelian variety.
        This is the subgroup generated by rational cusps.

        EXAMPLES:
            sage: J = J0(54)
            sage: C = J.cuspidal_subgroup()
            sage: C.gens()
            [[(1/3, 0, 0, 0, 0, 1/3, 0, 2/3)], [(0, 1/3, 0, 0, 0, 2/3, 0, 1/3)], [(0, 0, 1/9, 1/9, 1/9, 1/9, 1/9, 2/9)], [(0, 0, 0, 1/3, 0, 1/3, 0, 0)], [(0, 0, 0, 0, 1/3, 1/3, 0, 1/3)], [(0, 0, 0, 0, 0, 0, 1/3, 2/3)]]
            sage: C.invariants()
            [3, 3, 3, 3, 3, 9]
        """
        try:
            return self._cuspidal_subgroup
        except AttributeError:
            T = CuspidalSubgroup(self)
            self._cuspidal_subgroup = T
            return T

    def rational_cuspidal_subgroup(self):
        """
        Return the subgroup of this modular abelian variety generated
        by rational cusps.

        EXAMPLES:
            sage: J = J0(54)
            sage: CQ = J.rational_cuspidal_subgroup(); CQ
            Rational cuspidal subgroup with invariants [3, 3, 9] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(54)
            sage: CQ.gens()
            [[(1/3, 0, 0, 1/3, -1/3, -2/3, 1/3, 0)], [(0, 0, 1/9, 1/9, -2/9, -2/9, 1/9, -1/9)], [(0, 0, 0, 1, -1, -1, 2/3, -2/3)]]

            sage: factor(CQ.order())
            3^4
            sage: CQ.invariants()
            [3, 3, 9]
        """
        try:
            return self._rational_cuspidal_subgroup
        except AttributeError:
            T = RationalCuspidalSubgroup(self)
            self._rational_cuspidal_subgroup = T
            return T

    def zero_subgroup(self):
        """
        Return the zero subgroup of this modular abelian variety, as a
        finite group.

        EXAMPLES:
            sage: A =J0(54); G = A.zero_subgroup(); G
            Finite subgroup with invariants [] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(54)
            sage: G.is_subgroup(A)
            True
        """
        try:
            return self._zero_subgroup
        except AttributeError:
            G = FiniteSubgroup_gens(self, [], base_field=QQ)
            self._zero_subgroup = G
            return G

    def finite_subgroup(self, X, base_field=None):
        """
        Return a finite subgroup of this modular abelian variety.

        INPUT:
            X -- list of elements of other finite subgroups of
                 this modular abelian variety or elements that
                 coerce into the rational homology (viewed as
                 a rational vector space); also X could be
                 a finite subgroup itself that is contained
                 in this abelian variety.
            base_field -- (default: None) field over which this
                 group is defined.  If None try to figure out
                 the best base field.

        OUTPUT:
            a finite subgroup of a modular abelian variety

        EXAMPLES:
            sage: J = J0(11)
            sage: J.finite_subgroup([[1/5,0], [0,1/3]])
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)

            sage: J = J0(33); C = J[0].cuspidal_subgroup(); C
            Cuspidal subgroup with invariants [2, 2] over QQ of Abelian variety factor of dimension 1 of J0(33)
            sage: J.finite_subgroup([[0,0,0,0,0,1/6]])
            Finite subgroup with invariants [6] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(33)
            sage: J.finite_subgroup(C)
            Finite subgroup with invariants [2, 2] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(33)

        """
        if isinstance(X, FiniteSubgroup):
            if base_field is None:
                base_field = X.base_field()
            A = X.abelian_variety()
            if A.groups() != self.groups():
                raise ValueError, "ambient product Jacobians must be equal"
            if A == self:
                X = [v.element() for v in X.gens()]
            else:
                L = self.lattice()
                B = A.lattice().matrix()
                try:
                    # BROKEN
                    X = [L.coordinates(v.element()*B) for v in X.gens()]
                except ValueError:
                    raise TypeError, "unable to coerce subgroup into abelian variety."

        if base_field is None:
            base_field = QQbar

        return FiniteSubgroup_gens(self, X, base_field=base_field, check=True)


    def n_torsion_subgroup(self, n):
        """
        Return the $n$-torsion subgroup of elements of order dividing $n$
        of this modular abelian variety $A$, i.e., the group $A[n]$.

        EXAMPLES:
            sage: A = J0(23)
            sage: G = A.n_torsion_subgroup(5); G
            Finite subgroup with invariants [5, 5, 5, 5] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)
            sage: G.order()
            625
            sage: G.gens()
            [[(1/5, 0, 0, 0)], [(0, 1/5, 0, 0)], [(0, 0, 1/5, 0)], [(0, 0, 0, 1/5)]]
            sage: A = J0(23)
            sage: A.n_torsion_subgroup(2).order()
            16
        """
        n = ZZ(n)
        try:
            return self.__n_torsion_subgroup[n]
        except KeyError:
            pass
        except AttributeError:
            self.__n_torsion_subgroup = {}
        G = self.zero_subgroup()
        H = G.multiply(1/n)
        self.__n_torsion_subgroup[n] = H
        return H


    ###############################################################################
    # Decomposition
    ###############################################################################
    def decomposition(self):
        """
        Return a sequence of abelian subvarieties of self that are all simple,
        have finite intersection and sum to self.
        """
        try:
            return self.__decomposition
        except AttributeError:
            pass


        intersect = (self.dimension() < self._ambient_dimension())

        L = self.lattice()

        lattices = []
        S = self._ambient_modular_symbols_spaces()

        for i in range(len(S)):
            before = sum(S[j].dimension() for j in range(i))
            after  = sum(S[j].dimension()  for j in range(i+1,len(S)))
            M = S[i]
            for N in divisors(M.level()):
                P = M.ambient_module().modular_symbols_of_level(N)
                PS = P.cuspidal_subspace()
                zero_module = (QQ**M.ambient_module().dimension()).zero_submodule()
                D = PS.new_subspace().decomposition()
                for A in D:
                    # Now let B be the sum in the big ambient space
                    if N == M.level():
                        B = A
                    else:
                        # take all images of A at higher level
                        B = zero_module
                        for t in divisors(M.level()//N):
                            delta = A.degeneracy_map(M.level(), t).matrix()
                            B += delta.image()
                    # Figure out coordinates of this sum of images of A
                    # in terms of coordinates for the cuspidal subspace
                    # of modular symbols.
                    V = M.free_module().coordinate_module(B.free_module())
                    # Embed V in the space with 0's everywhere except at
                    # M factor.
                    AV = V.basis_matrix()
                    big = matrix(QQ,AV.nrows(), before).augment(AV).augment(matrix(QQ,AV.nrows(),after))
                    V_embed = big.row_module(QQ)
                    Z = V_embed.intersection(L)
                    if Z.dimension() > 0:
                        lattices.append((Z, Z.dimension() // A.dimension()))

        groups = self.groups()
        X = [ModularAbelianVariety(groups, L, QQ, check=False) for L, i in lattices]
        self.__decomposition = Sequence(X, immutable=True, cr=True, universe=self.category())
        return self.__decomposition

    def __getitem__(self, i):
        """
        Return the i-th decomposition factor of self.

        EXAMPLES:
            sage: J = J0(389)
            sage: J.decomposition()
            [
            Abelian variety factor of dimension 1 of J0(389),
            Abelian variety factor of dimension 2 of J0(389),
            Abelian variety factor of dimension 3 of J0(389),
            Abelian variety factor of dimension 6 of J0(389),
            Abelian variety factor of dimension 20 of J0(389)
            ]
            sage: J[2]
            Abelian variety factor of dimension 3 of J0(389)
            sage: J[-1]
            Abelian variety factor of dimension 20 of J0(389)
        """
        return self.decomposition()[i]

    def __getslice__(self, i, j):
        """
        The slice i:j of decompositions of self.

        EXAMPLES:
            sage: J = J0(125); J.decomposition()
            [
            Abelian variety factor of dimension 2 of J0(125),
            Abelian variety factor of dimension 2 of J0(125),
            Abelian variety factor of dimension 4 of J0(125)
            ]
            sage: J[:2]
            [
            Abelian variety factor of dimension 2 of J0(125),
            Abelian variety factor of dimension 2 of J0(125)
            ]
        """
        return self.decomposition()[i:j]



class ModularAbelianVariety(ModularAbelianVariety_abstract):
    def __init__(self, groups, lattice, base_field, check=True):
        r"""
        Create a modular abelian variety with given level and base field.

        INPUT:
            groups -- a tuple of congruence subgroups
            lattice -- a full lattice in $\ZZ^n$, where $n$ is the sum of
                       the dimensions of the spaces of cuspidal modular
                       symbols corresponding to each $\Gamma \in$ groups
            base_field -- a field

        EXAMPLES:
            sage: J0(23)
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)
        """
        if check:
            if not isinstance(groups, tuple):
                raise TypeError, "groups must be a tuple"
            for G in groups:
                if not is_CongruenceSubgroup(G):
                    raise TypeError, "each element of groups must be a congruence subgroup"
        self.__groups = groups

        if check:
            n = self._ambient_dimension()
            if not is_FreeModule(lattice):
                raise TypeError, "lattice must be a free module"
            if lattice.base_ring() != ZZ:
                raise TypeError, "lattice must be over ZZ"
            if lattice.degree() != 2*n:
                raise ValueError, "lattice must have degree n (=%s)"%n
            if not lattice.saturation().is_submodule(lattice):  # potentially expensive
                raise ValueError, "lattice must be full"
        self.__lattice = lattice

        ModularAbelianVariety_abstract.__init__(self, base_field, check=check)

    def groups(self):
        return self.__groups

    def lattice(self):
        return self.__lattice


class ModularAbelianVariety_modsym_abstract(ModularAbelianVariety_abstract):
    # Anything that derives from this class must define the
    # modular_symbols method, which returns a cuspidal modular
    # symbols space over QQ.  It can have any sign.
    def _modular_symbols(self):
        raise NotImplementedError, "bug -- must define this"


    def __add__(self, other):
        """
        Add two modular abelian variety factors.

        EXAMPLES:
            sage: A = J0(42); D = A.decomposition(); D
            [
            Abelian variety factor of dimension 1 of J0(42),
            Abelian variety factor of dimension 2 of J0(42),
            Abelian variety factor of dimension 2 of J0(42)
            ]
            sage: D[0] + D[1]
            Abelian variety factor of dimension 3 of J0(42)
            sage: D[1].is_subvariety(D[0] + D[1])
            True
            sage: D[0] + D[1] + D[2]
            Abelian variety J0(42)
            sage: D[0] + D[0]
            Abelian variety factor of dimension 1 of J0(42)
            sage: D[0] + D[0] == D[0]
            True
            sage: sum(D, D[0]) == A
            True
        """
        if not is_ModularAbelianVariety(other):
            raise TypeError, "sum not defined"
        if not isinstance(other, ModularAbelianVariety_modsym_abstract):
            return ModularAbelianVariety_abstract.__add__(self, other)
        if self.groups() != other.groups():
            raise TypeError, "sum not defined since ambient spaces different"
        M = self.modular_symbols(1) + other.modular_symbols(1)
        return ModularAbelianVariety_modsym(M)

    def groups(self):
        return (self._modular_symbols().group(), )

    def lattice(self):
        try:
            return self.__lattice
        except AttributeError:
            M = self.modular_symbols(0)
            S = M.ambient_module().cuspidal_submodule()
            if M.dimension() == S.dimension():
                s = 1 if M.sign() == 0 else 2
                L = ZZ**(s*M.dimension())
            else:
                K0 = M.integral_structure()
                K1 = S.integral_structure()
                L = K1.coordinate_module(K0)
            self.__lattice = L
            return self.__lattice

    def modular_symbols(self, sign=0):
        """
        Return space of modular symbols (with given sign) associated
        to this modular abelian variety.

        EXAMPLES:
            sage: A = J0(37)
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: A.modular_symbols(1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1 over Rational Field


        More examples:
            sage: J0(11).modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: J0(11).modular_symbols(sign=1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 2 for Gamma_0(11) of weight 2 with sign 1 over Rational Field
            sage: J0(11).modular_symbols(sign=0)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(11) of weight 2 with sign 0 over Rational Field
            sage: J0(11).modular_symbols(sign=-1)
            Modular Symbols space of dimension 1 for Gamma_0(11) of weight 2 with sign -1 over Rational Field

        Even more examples:
            sage: A = J0(33)[1]; A
            Abelian variety factor of dimension 2 of J0(33)
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 9 for Gamma_0(33) of weight 2 with sign 0 over Rational Field
            sage: A.modular_symbols(1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 6 for Gamma_0(33) of weight 2 with sign 1 over Rational Field
            sage: A.modular_symbols(-1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(33) of weight 2 with sign -1 over Rational Field
        """
        return self._modular_symbols().modular_symbols_of_sign(sign)

    def hecke_polynomial(self, n, var='x'):
        """
        Return the characteristic polynomial of the $n$-th Hecke
        operator on self.

        NOTE: If self has dimension d, then this is a polynomial of
        degree d.  It is not of degree 2*d, so it is the square root
        of the characteristic polynomial of the Hecke operator on
        integral or rational homology (which has degree 2*d).

        EXAMPLES:
            sage: factor(J0(11).hecke_polynomial(2))
            x + 2
            sage: factor(J0(23).hecke_polynomial(2))
            x^2 + x - 1
            sage: factor(J1(13).hecke_polynomial(2))
            x^2 + 3*x + 3
            sage: factor(J0(43).hecke_polynomial(2))
            (x + 2) * (x^2 - 2)

        The Hecke polynomial is the square root of the characteristic
        polynomial:
            sage: factor(J0(43).hecke_operator(2).charpoly())
            (x + 2)^2 * (x^2 - 2)^2
        """
        return self.modular_symbols(sign=1).hecke_polynomial(n, var)

    def __cmp__(self, other):
        """
        Compare two modular abelian varieties associated to spaces of
        cuspidal modular symbols if possible; otherwise, fallback to
        generic comparison.

        If other is a modular abelian variety attached to modular
        symbols, then this function compares the underlying +1 modular
        symbols spaces.  Otherwise it just compares the underlying
        types.

        EXAMPLES:
            sage: A = J0(37)
            sage: cmp(A,A)
            0
            sage: cmp(A,J0(43))
            -1
            sage: cmp(J0(43),A)
            1

        cmp also works when other is not a modular abelian variety.
            sage: cmp(A,17) #random (meaningless since it depends on memory layout)
            1
            sage: cmp(17,A) #random (meaningless since it depends on memory layout)
            -1
        """
        if isinstance(other, ModularAbelianVariety_modsym):
            return cmp(self.modular_symbols(1), other.modular_symbols(1))
        else:
            return ModularAbelianVariety_abstract.__cmp__(self, other)

    def _integral_hecke_matrix(self, n, sign=0):
        """
        Return the action of the Hecke operator $T_n$ on the
        integral homology of self.

        INPUT:
            n -- a positive integer
            sign -- 0, +1, or -1; if 1 or -1 act on the +1 or
                   -1 quotient of the integral homology.

        EXAMPLES:
            sage: J1(13)._integral_hecke_matrix(2)     # slightly random choice of basis
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J1(13)._integral_hecke_matrix(2,sign=1)  # slightly random choice of basis
            [-1  1]
            [-1 -2]
            sage: J1(13)._integral_hecke_matrix(2,sign=-1)  # slightly random choice of basis
            [-2 -1]
            [ 1 -1]
        """
        return self.modular_symbols(sign).integral_hecke_matrix(n)

    def _rational_hecke_matrix(self, n, sign=0):
        """
        Return the action of the Hecke operator $T_n$ on the
        rational homology of self.

        INPUT:
            n -- a positive integer
            sign -- 0, +1, or -1; if 1 or -1 act on the +1 or
                   -1 quotient of the rational homology.

        EXAMPLES:
            sage: J1(13)._rational_hecke_matrix(2)    # slightly random choice of basis
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J0(43)._rational_hecke_matrix(2,sign=1)  # slightly random choice of basis
            [-2  0  1]
            [-1 -2  2]
            [-2  0  2]
        """
        return self.modular_symbols(sign).hecke_matrix(n)

    def group(self):
        """
        Return the congruence subgroup associated that this modular abelian
        variety is associated to.

        EXAMPLES:
            sage: J0(13).group()
            Congruence Subgroup Gamma0(13)
            sage: J1(997).group()
            Congruence Subgroup Gamma1(997)
            sage: JH(37,[3]).group()
            Congruence Subgroup Gamma_H(37) with H generated by [3]
            sage: J0(37)[1].groups()
            (Congruence Subgroup Gamma0(37),)
        """
        return self.modular_symbols(1).group()

    def is_subvariety(self, other):
        """
        Return True if self is a subvariety of other.

        EXAMPLES:
            sage: J = J0(37); J
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: A = J[0]; A
            Abelian variety factor of dimension 1 of J0(37)
            sage: A.is_subvariety(J)
            True
            sage: A.is_subvariety(J0(11))
            False

        There may be a way to map $A$ into $J_0(74)$, but $A$ is
        not equipped with any special structure of an embedding.
            sage: A.is_subvariety(J0(74))
            False

        Some ambient examples:
            sage: J = J0(37)
            sage: J.is_subvariety(J)
            True
            sage: J.is_subvariety(25)
            False

        More examples:
            sage: A = J0(42); D = A.decomposition(); D
            [
            Abelian variety factor of dimension 1 of J0(42),
            Abelian variety factor of dimension 2 of J0(42),
            Abelian variety factor of dimension 2 of J0(42)
            ]
            sage: D[0].is_subvariety(A)
            True
            sage: D[1].is_subvariety(D[0] + D[1])
            True
            sage: D[2].is_subvariety(D[0] + D[1])
            False
        """
        if not is_ModularAbelianVariety(other):
            return False
        if not isinstance(other, ModularAbelianVariety_modsym_abstract):
            return ModularAbelianVariety_abstract.is_subvariety(self, other)
        return self.modular_symbols(1).is_submodule(other.modular_symbols(1))

    def is_ambient(self):
        return self.degree() == self.dimension()

    def dimension(self):
        """
        Return the dimension of this modular abelian variety.

        EXAMPLES:
            sage: J0(37)[0].dimension()
            1
            sage: J0(43)[1].dimension()
            2
            sage: J1(17)[1].dimension()
            4
        """
        try:
            return self._dimension
        except AttributeError:
            M = self._modular_symbols()
            if M.sign() == 0:
                d = M.dimension() // 2
            else:
                d = M.dimension()
            self._dimension = d
            return d

    def new_quotient(self, p=None):
        """
        Return the new or $p$-new quotient variety of self.

        INPUT:
            self -- a modular abelian variety
            p -- prime number or None (default); if p is a prime,
                 return the p-new quotient.  Otherwise return the
                 full new quotient.

        EXAMPLES:
            sage: J0(33).new_quotient()
            Abelian variety factor of dimension 1 of J0(33)
            sage: J0(100).new_quotient()
            Abelian variety factor of dimension 1 of J0(100)
            sage: J1(13).new_quotient()
            Abelian variety J1(13)
        """
        try:
            return self.__new_quotient[p]
        except AttributeError:
            self.__new_quotient = {}
        except KeyError:
            pass
        A = self.modular_symbols(sign=1)
        N = A.new_submodule(p=p)
        B = ModularAbelianVariety_modsym(N)
        self.__new_quotient[p] = B
        return B

    def old_quotient(self, p=None):
        """
        Return the old or $p$-old quotient variety of self.

        INPUT:
            self -- a modular abelian variety
            p -- prime number or None (default); if p is a prime,
                 return the p-old quotient.  Otherwise return the
                 full old quotient.

        EXAMPLES:
            sage: J0(33).old_quotient()
            Abelian variety factor of dimension 2 of J0(33)
            sage: J0(100).old_quotient()
            Abelian variety factor of dimension 6 of J0(100)
            sage: J1(13).old_quotient()
            Abelian variety factor of dimension 0 of J1(13)
        """
        try:
            return self.__old_quotient[p]
        except AttributeError:
            self.__old_quotient = {}
        except KeyError:
            pass
        A = self.modular_symbols(sign=1)
        N = A.old_submodule(p=p)
        B = ModularAbelianVariety_modsym(N)
        self.__old_quotient[p] = B
        return B

    def decomposition(self, bound=None):
        r"""
        Decompose this modular abelian variety as a product of Hecke
        equivariant modular abelian quotient varieties, up to isogeny.
        Each factor is an \emph{abelian subvariety} of self that
        corresponds to a newform of level dividing the level of self;
        in particular, each space of modular symbols can be cut out
        using Hecke operators of index coprime to the level.

        The old factors are \emph{not} simple!

        EXAMPLES:
            sage: J = J0(33)
            sage: J.decomposition()
            [
            Abelian variety factor of dimension 1 of J0(33),
            Abelian variety factor of dimension 2 of J0(33)
            ]
            sage: J1(17).decomposition()
            [
            Abelian variety factor of dimension 1 of J1(17),
            Abelian variety factor of dimension 4 of J1(17)
            ]
        """
        try:
            return self.__decomposition
        except AttributeError:
            pass
        A = self.modular_symbols(sign=1)

        D = Sequence([ModularAbelianVariety_modsym(B) for B in A.decomposition(bound = bound)],
                     immutable=True, cr=True, universe=self.category())

        self.__decomposition = D
        return D

class ModularAbelianVariety_modsym(ModularAbelianVariety_modsym_abstract):

    def __init__(self, modsym, check=True):
        """
        Modular abelian variety that corresponds to a Hecke stable
        space of cuspidal modular symbols.
        """
        if check:
            if not isinstance(modsym, ModularSymbolsSpace):
                raise TypeError, "modsym must be a modular symbols space"
            if not modsym.is_cuspidal():
                raise ValueError, "modsym must be cuspidal"

        ModularAbelianVariety_abstract.__init__(self, modsym.base_ring())
        self.__modsym = modsym

    def _modular_symbols(self):
        return self.__modsym

