"""
Base class for modular abelian varieties

AUTHOR:
    -- William Stein (2007-03)

TESTS:
    sage: A = J0(33)
    sage: D = A.decomposition(); D
    [
    Modular abelian variety quotient of dimension 1 and level 33,
    Modular abelian variety quotient of dimension 2 and level 33
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
from finite_subgroup            import FiniteSubgroup_gens
from cuspidal_subgroup          import CuspidalSubgroup, RationalCuspidalSubgroup
from sage.rings.all             import ZZ, QQ, QQbar, is_Field

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
    return isinstance(x, ModularAbelianVariety)

class ModularAbelianVariety(ParentWithBase):
    """
    A modular abelian variety.
    """
    def __init__(self, level, base_field):
        """
        Create a modular abelian variety with given level and base field.

        EXAMPLES:
            sage: J0(23)
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)
        """
        self.__level = level
        if not is_Field(base_field):
            raise TypeError, "base_field must be a field"
        ParentWithBase.__init__(self, base_field)

    def _repr_(self):
        """
        Return string representation of this modular abelian variety.

        This is just the generic base class, so it's unlikely to be called in practice.

        EXAMPLES:
            sage: A = J0(23)
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety._repr_(A)
            'Modular abelian variety of level 23 over Rational Field'
        """
        return "Modular abelian variety of level %s over %s"%(self.__level, self.base_ring())

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
        Return the level of this modular abelian variety, which is an integer
        N (usually minimal) such that this modular abelian variety is a quotient
        of $J_1(N)$.

        EXAMPLES:
            sage: J1(5077).level()
            5077
            sage: JH(389,[4]).level()
            389
        """
        return self.__level

    def base_field(self):
        r"""
        Synonym for \code{self.base_ring()}.

        EXAMPLES:
            sage: J0(11).base_field()
            Rational Field
        """
        return self.base_ring()

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

    def hecke_operator(self, n):
        """
        Return the n-th Hecke operator on the modular abelian variety.

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
        """
        Return the characteristic polynomial of the n-th Hecke
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

    def torsion_subgroup(self):
        """
        EXAMPLES:
            sage: J = J0(33)
            sage: A = J.new_quotient()
            sage: A
            Modular abelian variety quotient of dimension 1 and level 33
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
            Torsion subgroup of Modular abelian variety quotient of dimension 1 and level 33
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

    def finite_subgroup(self, gens):
        """
        Return a finite subgroup of this modular abelian variety.

        INPUT:
            gens -- either elements of other finite subgroups of
                    this modular abelian variety or elements that
                    coerce into the rational homology (viewed as
                    a rational vector space).

        OUTPUT:
            a finite subgroup of a modular abelian variety

        EXAMPLES:
            sage: J = J0(11)
            sage: J.finite_subgroup([[1/5,0], [0,1/3]])
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
        """
        return FiniteSubgroup_gens(self, gens, base_field=QQbar, check=True)


    def n_torsion_subgroup(self, n):
        """
        Return the n-torsion subgroup of elements of order dividing n
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


    def dimension(self):   # Derived classes *must* define this:
        """
        Return the dimension of this abelian variety.

        This function must be overloaded in the derived class.
        It just raises a NotImplementedError if called directly or
        not overloaded (which is a bug).

        EXAMPLES:
            sage: A = J0(23)
            sage: import sage.modular.abvar.abvar as abvar
            sage: A.dimension()
            2
            sage: abvar.ModularAbelianVariety.dimension(A)
            Traceback (most recent call last):
            ...
            NotImplementedError: bug in Sage; dimension method must be defined in the derived class
        """
        raise NotImplementedError, "bug in Sage; dimension method must be defined in the derived class"

    def is_subvariety(self, other):             # Derived classes *must* define this:
        """
        Return True if self is known to be a subvariety of other.

        EXAMPLES:
            sage: J = J0(37); J
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: A = J[0]; A
            Modular abelian variety quotient of dimension 1 and level 37
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety.is_subvariety(A, A)
            True

        For anything nontrivial the derived class must implement the
        functionality:
            sage: abvar.ModularAbelianVariety.is_subvariety(A, J)
            Traceback (most recent call last):
            ...
            NotImplementedError
            sage: A.is_subvariety(J)
            True
        """
        if self is other:
            return True
        raise NotImplementedError


    def change_ring(self, R):                   # Derived classes *must* define this:
        """
        Change the base ring of this modular abelian variety.

        This must be defined in the derived class.

        EXAMPLES:
            sage: A = J0(23)
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety.change_ring(A, QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError: bug in Sage; change_ring must be defined in the derived class
        """
        raise NotImplementedError, "bug in Sage; change_ring must be defined in the derived class"

    def _integral_hecke_matrix(self, n):        # derived classes may define
        """
        Return the matrix of the Hecke operator $T_n$ acting on the
        integral homology of this modular abelian variety, if this
        action is defined.  Otherwise, raise a ValueError.

        EXAMPLES:
            sage: A = J0(23)
            sage: t = A._integral_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Integer Ring

        This is usually defined in the derived class.   The base
        class method raises a ValueError, as illustrated below:
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety._integral_hecke_matrix(A, 2)
            Traceback (most recent call last):
            ...
            ValueError: no action of Hecke operators over ZZ
        """
        # this is allowed to raise an error if associated modular
        # symbols space has no Hecke action, e.g., isn't hecke
        # invariant.
        raise ValueError, "no action of Hecke operators over ZZ"

    def _rational_hecke_matrix(self, n):        # derived classes may define
        r"""
        Return the matrix of the Hecke operator $T_n$ acting on the
        rational homology $H_1(A,\Q)$ of this modular abelian variety,
        if this action is defined.  Otherwise, raise a ValueError.

        EXAMPLES:
            sage: A = J0(23)
            sage: t = A._rational_hecke_matrix(2); t
            [ 0  1 -1  0]
            [ 0  1 -1  1]
            [-1  2 -2  1]
            [-1  1  0 -1]
            sage: t.parent()
            Full MatrixSpace of 4 by 4 dense matrices over Rational Field

        This is usually defined in the derived class.   The base
        class method raises a ValueError, as illustrated below:
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety._rational_hecke_matrix(A, 2)
            Traceback (most recent call last):
            ...
            ValueError: no action of Hecke operators over QQ
        """
        raise ValueError, "no action of Hecke operators over QQ"


class ModularAbelianVariety_modsym(ModularAbelianVariety):
    """
    Abstract base class for modular abelian variety that corresponds
    to a space of cuspidal modular symbols.
    """
    def modular_symbols(self, sign=0):   # derived classes must overload this.
        """
        Return the modular symbols space associated to self.

        This must be defined in the derived class.

        EXAMPLES:
            sage: A = J0(37)

        The base class function must be redefined in the derived
        class.  Since A is in a derived class this is the case here:
            sage: type(A)
            <class 'sage.modular.abvar.abvar_ambient_jacobian.ModAbVar_ambient_jacobian_class'>
            sage: import sage.modular.abvar.abvar as abvar
            sage: abvar.ModularAbelianVariety_modsym.modular_symbols(A)
            Traceback (most recent call last):
            ...
            NotImplementedError: bug -- no associated modular symbols space

        The function is defined in the derived class.
            sage: A.modular_symbols()
            Modular Symbols subspace of dimension 4 of Modular Symbols space of dimension 5 for Gamma_0(37) of weight 2 with sign 0 over Rational Field
            sage: A.modular_symbols(1)
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 3 for Gamma_0(37) of weight 2 with sign 1 over Rational Field
        """
        raise NotImplementedError, "bug -- no associated modular symbols space"

    def __cmp__(self, other):
        """
        Compare two modular abelian varieties associated to spaces of
        cuspidal modular symbols.

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
        if not isinstance(other, ModularAbelianVariety_modsym):
            return cmp(type(self), type(other))
        return cmp(self.modular_symbols(1), other.modular_symbols(1))

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
            sage: J0(37)[1].group()
            Congruence Subgroup Gamma0(37)
        """
        return self.modular_symbols(1).group()

    def is_subvariety(self, other):
        """
        Return True if self is known to be a subvariety of other.

        EXAMPLES:
            sage: J = J0(37); J
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: A = J[0]; A
            Modular abelian variety quotient of dimension 1 and level 37
            sage: A.is_subvariety(J)
            True
            sage: A.is_subvariety(J0(11))
            False

        There may be a way to map $A$ into $J_0(74)$, but $A$ is
        not equipped with any special structure of an embedding.
            sage: A.is_subvariety(J0(74))
            False
        """
        if not isinstance(other, ModularAbelianVariety):
            return False
        if not isinstance(other, ModularAbelianVariety_modsym):
            return NotImplementedError, "general inclusion checking not implemented"
        return self.modular_symbols(1).is_submodule(other.modular_symbols(1))


    def dimension(self):
        """
        Return the dimension of this modular abelian variety.

        EXAMPLES:
            sage: J0(37)[1].dimension()
            1
            sage: J0(43)[1].dimension()
            2
            sage: J1(17)[1].dimension()
            4
        """
        try:
            return self._dimension
        except AttributeError:
            d = self.modular_symbols(sign=1).dimension()
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
            Modular abelian variety quotient of dimension 1 and level 33
            sage: J0(100).new_quotient()
            Modular abelian variety quotient of dimension 1 and level 100
            sage: J1(13).new_quotient()
            Modular abelian variety quotient of dimension 2 and level 13
        """
        try:
            return self.__new_quotient[p]
        except AttributeError:
            self.__new_quotient = {}
        except KeyError:
            pass
        A = self.modular_symbols(sign=1)
        N = A.new_submodule(p=p)

        from abvar_modsym_factor import ModAbVar_modsym_factor
        B = ModAbVar_modsym_factor(self, N)
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
            Modular abelian variety quotient of dimension 2 and level 33
            sage: J0(100).old_quotient()
            Modular abelian variety quotient of dimension 6 and level 100
            sage: J1(13).old_quotient()
            Modular abelian variety quotient of dimension 0 and level 13
        """
        try:
            return self.__old_quotient[p]
        except AttributeError:
            self.__old_quotient = {}
        except KeyError:
            pass
        A = self.modular_symbols(sign=1)
        N = A.old_submodule(p=p)

        from abvar_modsym_factor import ModAbVar_modsym_factor
        B = ModAbVar_modsym_factor(self, N)
        self.__old_quotient[p] = B
        return B

    def decomposition(self, bound=None):
        r"""
        Decompose this modular abelian variety as a product of Hecke
        equivariant modular abelian quotient varieties, up to isogeny.
        Each factor is an \emph{abelian subvariety} of self that
        corresponds to a newform of level dividing the level of self;
        in particular, each space of modular symbols that can be cut
        out using Hecke operators of index coprime to the level.

        EXAMPLES:
            sage: J = J0(33)
            sage: J.decomposition()
            [
            Modular abelian variety quotient of dimension 1 and level 33,
            Modular abelian variety quotient of dimension 2 and level 33
            ]
            sage: J1(17).decomposition()
            [
            Modular abelian variety quotient of dimension 1 and level 17,
            Modular abelian variety quotient of dimension 4 and level 17
            ]
        """
        try:
            return self.__decomposition
        except AttributeError:
            pass
        A = self.modular_symbols(sign=1)
        from abvar_modsym_factor import ModAbVar_modsym_factor
        D = Sequence([ModAbVar_modsym_factor(self, B) for B in A.decomposition(bound = bound)],
                     immutable=True, cr=True, universe=self.category())
        self.__decomposition = D
        return D

    def __getitem__(self, i):
        """
        Return the i-th decomposition factor of self.

        EXAMPLES:
            sage: J = J0(389)
            sage: J.decomposition()
            [
            Modular abelian variety quotient of dimension 1 and level 389,
            Modular abelian variety quotient of dimension 2 and level 389,
            Modular abelian variety quotient of dimension 3 and level 389,
            Modular abelian variety quotient of dimension 6 and level 389,
            Modular abelian variety quotient of dimension 20 and level 389
            ]
            sage: J[2]
            Modular abelian variety quotient of dimension 3 and level 389
            sage: J[-1]
            Modular abelian variety quotient of dimension 20 and level 389
        """
        return self.decomposition()[i]

    def __getslice__(self, i, j):
        """
        The slice i:j of decompositions of self.

        EXAMPLES:
            sage: J = J0(125); J.decomposition()
            [
            Modular abelian variety quotient of dimension 2 and level 125,
            Modular abelian variety quotient of dimension 2 and level 125,
            Modular abelian variety quotient of dimension 4 and level 125
            ]
            sage: J[:2]
            [
            Modular abelian variety quotient of dimension 2 and level 125,
            Modular abelian variety quotient of dimension 2 and level 125
            ]
            sage: J[1:]
            [
            Modular abelian variety quotient of dimension 2 and level 125,
            Modular abelian variety quotient of dimension 4 and level 125
            ]
        """
        return self.decomposition()[i:j]
