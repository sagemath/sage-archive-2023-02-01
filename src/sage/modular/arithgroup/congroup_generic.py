r"""
Congruence arithmetic subgroups of `{\rm SL}_2(\mathbb{Z})`

Sage can compute extensively with the standard congruence subgroups
`\Gamma_0(N)`, `\Gamma_1(N)`, and `\Gamma_H(N)`.

AUTHOR:
    -- William Stein
"""

################################################################################
#
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from sage.rings.all import QQ, ZZ, divisors, euler_phi
from sage.matrix.matrix_space import MatrixSpace
#import sage.modular.cusps # circular!
from sage.misc.misc import ellipsis_range
from sage.misc.cachefunc import cached_method

from arithgroup_element import ArithmeticSubgroupElement
from arithgroup_generic import ArithmeticSubgroup

def is_CongruenceSubgroup(x):
    """
    Return True if x is of type CongruenceSubgroup.

    EXAMPLES::

        sage: from sage.modular.arithgroup.congroup_generic import is_CongruenceSubgroup
        sage: is_CongruenceSubgroup(SL2Z)
        True
        sage: is_CongruenceSubgroup(Gamma0(13))
        True
        sage: is_CongruenceSubgroup(Gamma1(6))
        True
        sage: is_CongruenceSubgroup(GammaH(11, [3]))
        True
        sage: is_CongruenceSubgroup(SymmetricGroup(3))
        False
    """
    return isinstance(x, CongruenceSubgroup)

class CongruenceSubgroup(ArithmeticSubgroup):
    def __init__(self, level):
        """
        Create a congruence subgroup with given level.

        EXAMPLES::

            sage: Gamma0(500)
            Congruence Subgroup Gamma0(500)
        """
        level = ZZ(level)
        if level <= 0:
            raise ArithmeticError, "Congruence groups only defined for positive levels."
        self.__level = level

    def _repr_(self):
        """
        Return the string representation of self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5)._repr_()
            'Generic congruence subgroup'
        """
        return "Generic congruence subgroup"

    def is_congruence(self):
        r"""
        Return True, since this is a congruence subgroup.

        EXAMPLE::

            sage: Gamma0(7).is_congruence()
            True
        """

        return True

    def modular_symbols(self, sign=0, weight=2, base_ring=QQ):
        """
        Return the space of modular symbols of the specified weight and sign
        on the congruence subgroup self.

        EXAMPLES::

            sage: G = Gamma0(23)
            sage: G.modular_symbols()
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: G.modular_symbols(weight=4)
            Modular Symbols space of dimension 12 for Gamma_0(23) of weight 4 with sign 0 over Rational Field
            sage: G.modular_symbols(base_ring=GF(7))
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Finite Field of size 7
            sage: G.modular_symbols(sign=1)
            Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Rational Field
        """
        from sage.modular.modsym.modsym import ModularSymbols
        return ModularSymbols(self, sign=sign, weight=weight, base_ring=base_ring)

    def modular_abelian_variety(self):
        """
        Return the modular abelian variety corresponding to the congruence
        subgroup self.

        EXAMPLES::

            sage: Gamma0(11).modular_abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: Gamma1(11).modular_abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: GammaH(11,[3]).modular_abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1
        """
        from sage.modular.abvar.abvar_ambient_jacobian import ModAbVar_ambient_jacobian
        return ModAbVar_ambient_jacobian(self)

    def level(self):
        """
        Return the level of this congruence subgroup.

        EXAMPLES::

            sage: SL2Z.level()
            1
            sage: Gamma0(20).level()
            20
            sage: Gamma1(11).level()
            11
            sage: GammaH(14, [2]).level()
            14
        """
        return self.__level

    def _new_group_from_level(self, level):
        r"""
        Return a new group of the same type (Gamma0, Gamma1, or
        GammaH) as self of the given level. In the case that self is
        of type GammaH, we take the largest H inside
        $(Z/\text{level}Z)^\times$ which maps to H, namely its inverse
        image under the natural reduction map.

        EXAMPLES::

            sage: G = Gamma0(20)
            sage: G._new_group_from_level(4)
            Congruence Subgroup Gamma0(4)
            sage: G._new_group_from_level(40)
            Congruence Subgroup Gamma0(40)

            sage: G = Gamma1(10)
            sage: G._new_group_from_level(6)
            Traceback (most recent call last):
            ...
            ValueError: one level must divide the other

            sage: G = GammaH(50,[3,37])
            sage: G
            Congruence Subgroup Gamma_H(50) with H generated by [3, 37]
            sage: G._new_group_from_level(25)
            Congruence Subgroup Gamma_H(25) with H generated by [3, 12]
            sage: G._new_group_from_level(100)
            Congruence Subgroup Gamma_H(100) with H generated by [3, 37, 53, 87]
        """

        from congroup_gamma0 import is_Gamma0
        from congroup_gamma1 import is_Gamma1
        from congroup_gammaH import is_GammaH
        from all import Gamma0, Gamma1, GammaH
        N = self.level()
        if (level%N) and (N%level):
            raise ValueError, "one level must divide the other"
        if is_Gamma0(self):
            return Gamma0(level)
        elif is_Gamma1(self):
            return Gamma1(level)
        elif is_GammaH(self):
            H = self._generators_for_H()
            if level > N:
                d = level // N
                diffs = [ N*i for i in range(d) ]
                return GammaH(level, [ h + diff for h in H for diff in diffs ])
            else:
                return GammaH(level, [ h%level for h in H ])
        else:
            raise NotImplementedError

    def __call__(self, x, check=True):
        """
        Coerce x into self.

        NOTE: This function should be overridden by any subclass the
        user will interact with directly.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).__call__(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if isinstance(x, ArithmeticSubgroupElement) and x.parent() == self:
            return x
        raise NotImplementedError


    def sturm_bound(self, weight=2):
        r"""
        Returns the Sturm bound for modular forms of the given weight and level
        this congruence subgroup.

        INPUT:

        -  ``weight`` - an integer `\geq 2` (default: 2)

        EXAMPLES::
            sage: Gamma0(11).sturm_bound(2)
            2
            sage: Gamma0(389).sturm_bound(2)
            65
            sage: Gamma0(1).sturm_bound(12)
            1
            sage: Gamma0(100).sturm_bound(2)
            30
            sage: Gamma0(1).sturm_bound(36)
            3
            sage: Gamma0(11).sturm_bound()
            2
            sage: Gamma0(13).sturm_bound()
            3
            sage: Gamma0(16).sturm_bound()
            4
            sage: GammaH(16,[13]).sturm_bound()
            8
            sage: GammaH(16,[15]).sturm_bound()
            16
            sage: Gamma1(16).sturm_bound()
            32
            sage: Gamma1(13).sturm_bound()
            28
            sage: Gamma1(13).sturm_bound(5)
            70

        FURTHER DETAILS: This function returns a positive integer
        `n` such that the Hecke operators
        `T_1,\ldots, T_n` acting on *cusp  forms* generate the
        Hecke algebra as a `\mathbb{Z}`-module when the character
        is trivial or quadratic. Otherwise, `T_1,\ldots,T_n`
        generate the Hecke algebra at least as a
        `\mathbb{Z}[\varepsilon]`-module, where
        `\mathbb{Z}[\varepsilon]` is the ring generated by the
        values of the Dirichlet character `\varepsilon`.
        Alternatively, this is a bound such that if two cusp forms
        associated to this space of modular symbols are congruent modulo
        `(\lambda, q^n)`, then they are congruent modulo
        `\lambda`.

        REFERENCES:

        - See the Agashe-Stein appendix to Lario and Schoof,
          *Some computations with Hecke rings and deformation rings*,
          Experimental Math., 11 (2002), no. 2, 303-311.

        - This result originated in the paper Sturm,
          *On the congruence of modular forms*,
          Springer LNM 1240, 275-280, 1987.

        REMARK: Kevin Buzzard pointed out to me (William Stein) in Fall
        2002 that the above bound is fine for `\Gamma_1(N)` with
        character, as one sees by taking a power of `f`. More
        precisely, if `f \cong 0 \pmod{p}` for first
        `s` coefficients, then `f^r \cong 0 \pmod{p}` for
        first `sr` coefficents. Since the weight of `f^r`
        is `r\cdot k(f)`, it follows that if
        `s \geq b`, where `b` is the Sturm bound for
        `\Gamma_0(N)` at weight `k(f)`, then `f^r`
        has valuation large enough to be forced to be `0` at
        `r*k(f)` by Sturm bound (which is valid if we choose
        `r` correctly). Thus `f \cong 0 \pmod{p}`.
        Conclusion: For `\Gamma_1(N)` with fixed character, the
        Sturm bound is *exactly* the same as for `\Gamma_0(N)`.

        A key point is that we are finding
        `\mathbb{Z}[\varepsilon]` generators for the Hecke algebra
        here, not `\mathbb{Z}`-generators. So if one wants
        generators for the Hecke algebra over `\mathbb{Z}`, this
        bound must be suitably modified (and I'm not sure what the
        modification is).

        AUTHORS:

        - William Stein
        """
        return ZZ((self.index() * weight / ZZ(12)).ceil())

#import congroup_pyx
#degeneracy_coset_representatives_gamma0 = congroup_pyx.degeneracy_coset_representatives_gamma0
#degeneracy_coset_representatives_gamma1 = congroup_pyx.degeneracy_coset_representatives_gamma1

