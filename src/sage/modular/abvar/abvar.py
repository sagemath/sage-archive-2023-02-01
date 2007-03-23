"""
Base class for modular abelian varieties

AUTHOR:
    -- William Stein (2007-03)
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.categories.all        import ModularAbelianVarieties
from sage.structure.sequence    import Sequence
from sage.structure.sage_object import SageObject
from hecke_operator             import HeckeOperator
from torsion_subgroup           import TorsionSubgroup
from finite_subgroup            import FiniteSubgroup_gens
from cuspidal_subgroup          import CuspidalSubgroup, RationalCuspidalSubgroup
from sage.rings.all             import ZZ, QQ

import homology

def is_ModularAbelianVariety(x):
    return isinstance(x, ModularAbelianVariety)

class ModularAbelianVariety(SageObject):
    def __init__(self, level, base_ring):
        self._level = level
        self._base_ring = base_ring


    def _repr_(self):
        return "Modular abelian variety of level %s over %s"%(self._level, self._base_ring)

    def category(self):
        try:
            return self.__category
        except AttributeError:
            C = ModularAbelianVarieties(self._base_ring)
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
        return self._level

    def base_ring(self):
        """
        Return the ring that this modular abelian varety is defined over.

        EXAMPLES:
            sage: J0(11).base_ring()
            Rational Field
        """
        return self._base_ring

    def base_field(self):
        r"""
        Synonym for \code{self.base_ring()}.

        EXAMPLES:
            sage: J0(11).base_field()
            Rational Field
        """
        return self.base_ring()

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

        If you just ask for the rank of the homology, no serious calculations are done, so the
        following is fast:
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
        integral or rational homology (which has degree degree 2*d).

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
            Rational cuspidal subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(54)
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
        try:
            return self._zero_subgroup
        except AttributeError:
            G = FiniteSubgroup_gens(self, [], base_field=QQ)
            self._zero_subgroup = G
            return G

    def n_torsion_subgroup(self, n):
        try:
            return self._n_torsion_subgroup
        except AttributeError:
            G = self.zero_subgroup()
            H = G.multiply(1/ZZ(n))
            self._n_torsion_subgroup = H
            return H


    def dimension(self):   # Derived classes *must* overload this:
        raise NotImplementedError


    def change_ring(self, R):                   # Derived classes *must* overload this:
        raise NotImplementedError

    def _integral_hecke_matrix(self, n):        # derived classes should overload
        # this is allowed to raise an error if
        # associated modular symbols space.
        raise ValueError, "no action of Hecke operators over ZZ"

    def _rational_hecke_matrix(self, n):        # derived classes should overload
        raise ValueError, "no action of Hecke operators over QQ"


class ModularAbelianVariety_modsym(ModularAbelianVariety):
    """
    Modular abelian variety that corresponds to a space of
    cuspidal modular symbols.
    """
    # derived classes must overload this.
    def modular_symbols(self, sign=0):
        """
        Return the modular symbols space associated to self.

        This raises a ValueError if there is no associated modular
        symbols space.
        """
        raise ValueError, "no associated modular symbols space"

    def __cmp__(self, other):
        if not isinstance(other, ModularAbelianVariety_modsym):
            return cmp(type(self), type(other))
        return cmp(self.modular_symbols(1), other.modular_symbols(1))

    def _integral_hecke_matrix(self, n, sign=0):
        """
        EXAMPLES:
            sage: J1(13)._integral_hecke_matrix(2)
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J1(13)._integral_hecke_matrix(2,sign=1)
            [-1  1]
            [-1 -2]
            sage: J1(13)._integral_hecke_matrix(2,sign=-1)
            [-2 -1]
            [ 1 -1]
        """
        return self.modular_symbols(sign).integral_hecke_matrix(n)

    def _rational_hecke_matrix(self, n, sign=0):
        """
        EXAMPLES:
            sage: J1(13)._rational_hecke_matrix(2)
            [-2  0 -1  1]
            [ 1 -1  0 -1]
            [ 1  1 -2  0]
            [ 0  1 -1 -1]
            sage: J0(43)._rational_hecke_matrix(2,sign=1)
            [-2  0  1]
            [-1 -2  2]
            [-2  0  2]
        """
        return self.modular_symbols(sign).hecke_matrix(n)

    def group(self):
        """
        Return the congruence subgroup associated this this modular abelian
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
        Return the new or p-new quotient variety of self.

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
        Return the old or p-old quotient variety of self.

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
        """
        Decompose this modular abelian variety as a product of Hecke
        equivariant modular abelian quotient varieties, up to isogeny.
        Each factor is an optimal quotient of self that corresponds to
        a new form of level dividing the level of self; in particular,
        each space of modular symbols that can be cut out using Hecke
        operators of index coprime to the level.

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
