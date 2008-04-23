r"""
Congruence subgroups of SL2(Z)

\sage can compute with the congruence subgroups $\Gamma_0(N)$,
$\Gamma_1(N)$, and $\Gamma_H(N)$.

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

import random

import sage.rings.arith as arith

from sage.groups.group import Group

from sage.rings.integer_ring import IntegerRing, ZZ
from sage.rings.integer_mod_ring import IntegerModRing
from sage.matrix.matrix_space import MatrixSpace

from sage.rings.infinity import infinity


from congroup_element import CongruenceSubgroupElement

import sage.modular.modsym.p1list

from math import sqrt
import sage.ext.arith as fast_arith

from sage.rings.all import QQ, divisors

# Just for now until we make an SL_2 group type.
from sage.matrix.matrix_space import MatrixSpace
Mat2Z = MatrixSpace(IntegerRing(),2)

def get_gcd(order):
    if order <= 46340:   # todo: don't hard code
        return fast_arith.arith_int().gcd_int
    elif order <= 2147483647:   # todo: don't hard code
        return fast_arith.arith_llong().gcd_longlong
    raise NotImplementedError

def get_inverse_mod(order):
    if order <= 46340:   # todo: don't hard code
        return fast_arith.arith_int().inverse_mod_int
    elif order <= 2147483647:   # todo: don't hard code
        return fast_arith.arith_llong().inverse_mod_longlong
    raise NotImplementedError


def is_CongruenceSubgroup(x):
    """
    Return True if x is of type CongruenceSubgroup.

    EXAMPLES:
        sage: is_CongruenceSubgroup(SL2Z())
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

class CongruenceSubgroup(Group):
    def __init__(self, level):
        level = ZZ(level)
        if level <= 0:
            raise ArithmeticError, "Congruence groups only defined for positive levels."
        self.__level = level

    def _repr_(self):
        return "Generic congruence subgroup"

    def __hash__(self):
        return hash(str(self))

    def modular_symbols(self, sign=0, weight=2, base_ring=QQ):
        """
        Return the space of modular symbols of the specified weight and sign
        on the congruence subgroup self.

        EXAMPLES:
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

        EXAMPLES:
            sage: Gamma0(11).modular_abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: Gamma1(11).modular_abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: GammaH(11,[3]).modular_abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1
        """
        from sage.modular.abvar.abvar_ambient_jacobian import ModAbVar_ambient_jacobian
        return ModAbVar_ambient_jacobian(self)

    def are_equivalent(self, x, y):
        raise NotImplementedError

    def coset_reps(self):
        raise NotImplementedError

    def generators(self):
        raise NotImplementedError

    def gens(self):
        """
        Return tuple of generators for this congruence subgroup.

        The generators need not be minimal.

        EXAMPLES:
            sage: SL2Z().gens()
            ([ 0 -1]
            [ 1  0], [1 1]
            [0 1])
            sage: Gamma0(3).gens()
            ([1 1]
            [0 1],
            [-1  0]
            [ 0 -1],
            [ 1 -1]
            [ 0  1],
            [ 1 -1]
            [ 3 -2],
            [ 2 -1]
            [ 3 -1],
            [-2  1]
            [-3  1])
            sage: Gamma1(2).gens()
            ([1 1]
            [0 1],
            [-1  0]
            [ 0 -1],
            [ 1 -1]
            [ 0  1],
            [ 1 -1]
            [ 2 -1],
            [-1  1]
            [-2  1])
            sage: GammaH(3, [2]).gens()
            ([1 1]
            [0 1],
            [-1  0]
            [ 0 -1],
            [ 1 -1]
            [ 0  1],
            [ 1 -1]
            [ 3 -2],
            [ 2 -1]
            [ 3 -1],
            [-2  1]
            [-3  1])
        """
        return tuple(self.generators())

    def gen(self, i):
        """
        Return the i-th generator of self, i.e. the i-th element of the
        tuple self.gens().

        EXAMPLES:
            sage: SL2Z().gen(1)
            [1 1]
            [0 1]
            sage: Gamma0(20).gen(7)
            [17 11]
            [20 13]
            sage: Gamma1(16).gen(5)
            [ 1329  -758]
            [ 3056 -1743]
            sage: GammaH(13, [8]).gen(3)
            [ 885 -218]
            [1157 -285]
        """
        return self.generators()[i]

    def ngens(self):
        r"""
        Return the number of generators for this congruence subgroup.

        This need not be the minimal number of generators of self.

        EXAMPLES:
            sage: Gamma0(22).ngens()
            42
            sage: Gamma1(14).ngens()
            250
            sage: GammaH(11, [3]).ngens()
            31
            sage: SL2Z().ngens()
            2
        """
        return len(self.generators())

    def level(self):
        """
        Return the level of this congruence subgroup.

        EXAMPLES:
            sage: SL2Z().level()
            1
            sage: Gamma0(20).level()
            20
            sage: Gamma1(11).level()
            11
            sage: GammaH(14, [2]).level()
            14
        """
        return self.__level

    def __cmp__(self, right):
        raise NotImplementedError

    def is_abelian(self):
        """
        Return True if this congruence subgroup is abelian.

        Since congruence subgroups are always nonabelian, this always
        returns False.

        EXAMPLES:
            sage: SL2Z().is_abelian()
            False
            sage: Gamma0(3).is_abelian()
            False
            sage: Gamma1(12).is_abelian()
            False
            sage: GammaH(4, [2]).is_abelian()
            False
        """
        return False

    def is_finite(self):
        """
        Return True if this congruence subgroup is finite.

        Since congruence subgroups are always infinite, this always
        returns False.

        EXAMPLES:
            sage: SL2Z().is_finite()
            False
            sage: Gamma0(3).is_finite()
            False
            sage: Gamma1(12).is_finite()
            False
            sage: GammaH(4, [2]).is_finite()
            False
        """
        return False

    def is_subgroup(self, right):
        raise NotImplementedError

    def is_odd(self):
        """
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES:
            sage: SL2Z().is_odd()
            False
            sage: Gamma0(20).is_odd()
            False
            sage: Gamma1(5).is_odd()
            True
            sage: GammaH(11, [3]).is_odd()
            True
        """
        return not self.is_even()

    def is_even(self):
        """
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES:
            sage: SL2Z().is_even()
            True
            sage: Gamma0(20).is_even()
            True
            sage: Gamma1(5).is_even()
            False
            sage: GammaH(11, [3]).is_even()
            False
        """
        return [-1, 0, 0, -1] in self

    def _new_group_from_level(self, level):
        """
        Return a new group of the same type (Gamma0, Gamma1, or GammaH)
        as self of the given level. In the case that self is of type
        GammaH, we take the largest H inside $(Z/\text{level}Z)^\times$
        which maps to H, namely its inverse image under the natural map.

        EXAMPLES:
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
        N = self.level()
        if (level%N) and (N%level):
            raise ValueError, "one level must divide the other"
        if is_Gamma0(self):
            return Gamma0(level)
        elif is_Gamma1(self):
            return Gamma1(level)
        else:
            H = self._generators_for_H()
            if level > N:
                d = level // N
                diffs = [ N*i for i in range(d) ]
                return GammaH(level, [ h + diff for h in H for diff in diffs ])
            else:
                return GammaH(level, [ h%level for h in H ])

    def order(self):
        """
        Return the number of elements in this congruence subgroup.

        Since congruence subgroups are always infinite, this always returns
        infinity.

        EXAMPLES:
            sage: SL2Z().order()
            +Infinity
            sage: Gamma0(5).order()
            +Infinity
            sage: Gamma1(2).order()
            +Infinity
            sage: GammaH(12, [5]).order()
            +Infinity
        """
        return infinity

    def __call__(self, x, check=True):
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        raise NotImplementedError

def lift_to_sl2z(c, d, N):
    """
    If this Manin symbol is (c,d) viewed modulo N, this function
    computes and returns a list [a, b, c', d'] that defines a 2x2
    matrix with determinant 1 and integer entries, such that
    c=c'(mod N) and d=d'(mod N).

    EXAMPLES:
        sage: lift_to_sl2z(1, 0, 12)
        [0, -1, 1, 0]
        sage: lift_to_sl2z(29, 100, 7)
        [9, 31, 29, 100]
        sage: lift_to_sl2z(5, 10, 3)
        [2, 5, 5, 13]
        sage: lift_to_sl2z(1000, 10000, 10000)
        [0, -1, 1000, 20000]
    """
    if N == 1:
        return [1,0,0,1]
    g, z1, z2 = arith.XGCD(c,d)

    # We're lucky: z1*c + z2*d = 1.
    if g==1:
        return [z2, -z1, c, d]

    # Have to try harder.
    if c == 0:
        c += N;
    if d == 0:
        d += N;
    m = c;

    # compute prime-to-d part of m.
    while True:
        g = arith.GCD(m,d)
        if g == 1:
            break
        m //= g

    # compute prime-to-N part of m.
    while True:
        g = arith.GCD(m,N);
        if g == 1:
            break
        m //= g
    d += N*m
    g, z1, z2 = arith.XGCD(c,d)

    assert g==1

    return [z2, -z1, c, d]


def is_Gamma0(x):
    """
    Return True if x is a congruence subgroup of type Gamma0.

    EXAMPLES:
        sage: is_Gamma0(SL2Z())
        True
        sage: is_Gamma0(Gamma0(13))
        True
        sage: is_Gamma0(Gamma1(6))
        False
    """
    return isinstance(x, Gamma0)

class Gamma0(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_0(N)$.

    EXAMPLES:
        sage: G = Gamma0(11); G
        Congruence Subgroup Gamma0(11)
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self, level):
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        return "Congruence Subgroup Gamma0(%s)"%self.level()

    def _latex_(self):
        """
        EXAMPLES:
            sage: Gamma0(20)._latex_()
            '\\Gamma_0(20)'
            sage: latex(Gamma0(20))
            \Gamma_0(20)
        """
        return "\\Gamma_0(%s)"%self.level()

    def _generators_for_H(self):
        """
        EXAMPLES:
            sage: Gamma0(15)._generators_for_H()
            [11, 7]
        """
        try:
            return self.__generators_for_H
        except AttributeError:
            self.__generators_for_H = [int(x) for x in IntegerModRing(self.level()).unit_gens()]
            return self.__generators_for_H

    def _list_of_elements_in_H(self):
        """
        Returns a sorted list of Python ints that are representatives
        between 0 and N-1 of the elements of H.

        EXAMPLES:
            sage: G = Gamma0(11)
            sage: G._list_of_elements_in_H()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        return range(1, self.level())

    def __cmp__(self, right):
        if not isinstance(right, Gamma0):
            if isinstance(right, CongruenceSubgroup):
                c = cmp(self.level(), right.level())
                if c: return c
            return cmp(type(self), type(right))
        return cmp(self.level(), right.level())

    def is_even(self):
        """
        Return True precisely if this subgroup contains the matrix -1.

        Since Gamma0(N) always, contains the matrix -1, this always returns
        True.

        EXAMPLES:
            sage: Gamma0(12).is_even()
            True
            sage: SL2Z().is_even()
            True
        """
        return True

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES:
            sage: G = Gamma0(20)
            sage: G.is_subgroup(SL2Z())
            True
            sage: G.is_subgroup(Gamma0(4))
            True
            sage: G.is_subgroup(Gamma0(20))
            True
            sage: G.is_subgroup(Gamma0(7))
            False
            sage: Gamma0(2).is_subgroup(Gamma1(2))
            True
        """
        if right.level() == 1:
            return True
        if isinstance(right, Gamma0):
            return self.level() % right.level() == 0
        if isinstance(right, Gamma1):
            if right.level() >= 3:
                return False
            elif right.level() == 2:
                return self.level() == 2
            # case level 1 dealt with above
        raise NotImplementedError

    def coset_reps(self):
        r"""
        Return representatives for the right cosets of this
        congruence subgroup in ${\rm SL}_2(\Z)$ as a generator.

        Use \code{list(self.coset_reps())} to obtain coset reps as a
        list.

        EXAMPLES:
            sage: list(Gamma0(5).coset_reps())
            [[1, 0, 0, 1],
            [0, -1, 1, 0],
            [1, 0, 1, 1],
            [1, 1, 1, 2],
            [1, 2, 1, 3],
            [1, 3, 1, 4]]
            sage: list(Gamma0(4).coset_reps())
            [[1, 0, 0, 1],
            [0, -1, 1, 0],
            [1, 0, 1, 1],
            [1, 1, 1, 2],
            [1, 2, 1, 3],
            [-1, -1, 2, 1]]
            sage: list(Gamma0(1).coset_reps())
            [[1, 0, 0, 1]]
        """
        N = self.level()
        for z in sage.modular.modsym.p1list.P1List(N):
            yield lift_to_sl2z(z[0], z[1], N)

    def generators(self):
        r"""
        Return generators for this congruence subgroup.

        The result is cached.

        EXAMPLE:
            sage: for g in Gamma0(3).generators():
            ...     print g
            ...     print '---'
            [1 1]
            [0 1]
            ---
            [-1  0]
            [ 0 -1]
            ---
            ...
            ---
            [-2  1]
            [-3  1]
            ---

        """
        try:
            return self.__gens
        except AttributeError:
            from sage.modular.modsym.p1list import P1List
            from congroup_pyx import generators_helper
            level = self.level()
            gen_list = generators_helper(P1List(level), level, Mat2Z)
            self.__gens = [self(g, check=False) for g in gen_list]
            return self.__gens

    def gamma_h_subgroups(self):
        r"""
        Return the subgroups of the form $\Gamma_H(N)$ contained
        in self, where $N$ is the level of self.

        EXAMPLES:
            sage: G = Gamma0(11)
            sage: G.gamma_h_subgroups()
            [Congruence Subgroup Gamma_H(11) with H generated by [2], Congruence Subgroup Gamma_H(11) with H generated by [4], Congruence Subgroup Gamma_H(11) with H generated by [10], Congruence Subgroup Gamma_H(11) with H generated by []]
            sage: G = Gamma0(12)
            sage: G.gamma_h_subgroups()
            [Congruence Subgroup Gamma_H(12) with H generated by [5, 7], Congruence Subgroup Gamma_H(12) with H generated by [7], Congruence Subgroup Gamma_H(12) with H generated by [5], Congruence Subgroup Gamma_H(12) with H generated by []]
        """
        N = self.level()
        R = IntegerModRing(N)
        return [GammaH(N, H) for H in R.multiplicative_subgroups()]

    def __call__(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES:
            sage: G = Gamma0(12)
            sage: G([1, 0, 24, 1])
            [ 1  0]
            [24  1]
            sage: G(matrix(ZZ, 2, [1, 1, -12, -11]))
            [  1   1]
            [-12 -11]
            sage: G([1, 0, 23, 1])
            Traceback (most recent call last):
            ...
            TypeError: matrix must have lower left entry (=23) divisible by 12
        """
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        x = CongruenceSubgroupElement(self, x, check=check)
        if not check:
            return x

        c = x.c()
        N = self.level()
        if c%N == 0:
            return x
        else:
            raise TypeError, "matrix must have lower left entry (=%s) divisible by %s" %(c, N)

def is_SL2Z(x):
    """
    Return True if x is the modular group ${\rm SL}_2(\Z)$.

    EXAMPLES:
        sage: is_SL2Z(SL2Z())
        True
        sage: is_SL2Z(Gamma0(6))
        False
    """
    return isinstance(x, SL2Z)

class SL2Z(Gamma0):
    r"""
    The modular group ${\rm SL}_2(\Z)$.

    EXAMPLES:
        sage: G = SL2Z(); G
        Modular Group SL(2,Z)
        sage: G.gens()
        ([ 0 -1]
        [ 1  0], [1 1]
        [0 1])
        sage: G.0
        [ 0 -1]
        [ 1  0]
        sage: G.1
        [1 1]
        [0 1]
        sage: latex(G)
        \mbox{\rm SL}_2(\mathbf{Z})
        sage: G([1,-1,0,1])
        [ 1 -1]
        [ 0  1]
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self):
        Gamma0.__init__(self, 1)

    def _repr_(self):
        return "Modular Group SL(2,Z)"

    def _latex_(self):
        """
        EXAMPLES:
            sage: SL2Z()._latex_()
            '\\mbox{\\rm SL}_2(\\mathbf{Z})'
            sage: latex(SL2Z())
            \mbox{\rm SL}_2(\mathbf{Z})
        """
        return "\\mbox{\\rm SL}_2(%s)"%(IntegerRing()._latex_())

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES:
            sage: SL2Z().is_subgroup(SL2Z())
            True
            sage: SL2Z().is_subgroup(Gamma1(1))
            True
            sage: SL2Z().is_subgroup(Gamma0(6))
            False
        """
        return right.level() == 1

def is_Gamma1(x):
    """
    Return True if x is a congruence subgroup of type Gamma1.

    EXAMPLES:
        sage: is_Gamma1(SL2Z())
        True
        sage: is_Gamma1(Gamma1(13))
        True
        sage: is_Gamma1(Gamma0(6))
        False
    """
    return (isinstance(x, Gamma1) or isinstance(x, SL2Z))

class Gamma1(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_1(N)$.

    EXAMPLES:
        sage: G = Gamma1(11); G
        Congruence Subgroup Gamma1(11)
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self, level):
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        return "Congruence Subgroup Gamma1(%s)"%self.level()

    def _latex_(self):
        """
        EXAMPLES:
            sage: Gamma1(3)._latex_()
            '\\Gamma_1(3)'
            sage: latex(Gamma1(3))
            \Gamma_1(3)
        """
        return "\\Gamma_1(%s)"%self.level()

    def __cmp__(self, right):
        if not isinstance(right, Gamma1):
            if isinstance(right, CongruenceSubgroup):
                c = cmp(self.level(), right.level())
                if c: return c
            return cmp(type(self), type(right))
        return cmp(self.level(), right.level())

    def is_even(self):
        """
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES:
            sage: Gamma1(1).is_even()
            True
            sage: Gamma1(2).is_even()
            True
            sage: Gamma1(15).is_even()
            False
        """
        return self.level() in [1,2]

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES:
            sage: Gamma1(3).is_subgroup(SL2Z())
            True
            sage: Gamma1(3).is_subgroup(Gamma1(5))
            False
            sage: Gamma1(3).is_subgroup(Gamma1(6))
            False
            sage: Gamma1(6).is_subgroup(Gamma1(3))
            True
            sage: Gamma1(6).is_subgroup(Gamma0(2))
            True
        """
        if right.level() == 1:
            return True
        if isinstance(right, (Gamma1, Gamma0)):
            return self.level() % right.level() == 0
        raise NotImplementedError

    def generators(self):
        r"""
        Return generators for this congruence subgroup.

        The result is cached.

        EXAMPLE:
            sage: for g in Gamma1(3).generators():
            ...     print g
            ...     print '---'
            [1 1]
            [0 1]
            ---
            [ 31 -14]
            [ 51 -23]
            ---
            [-5  4]
            [-9  7]
            ---
            ...
            ---
            [4 3]
            [9 7]
            ---
            [ -5  -2]
            [-12  -5]
            ---

        """
        try:
            return self.__gens
        except AttributeError:
            from sage.modular.modsym.g1list import G1list
            from congroup_pyx import generators_helper
            level = self.level()
            gen_list = generators_helper(G1list(level), level, Mat2Z)
            self.__gens = [self(g, check=False) for g in gen_list]
            return self.__gens

    def __call__(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES:
            sage: G = Gamma1(5)
            sage: G([1, 0, -10, 1])
            [ 1   0]
            [-10  1]
            sage: G(matrix(ZZ, 2, [6, 1, 5, 1]))
            [6  1]
            [5  1]
            sage: G([1, 1, 6, 7])
            Traceback (most recent call last):
            ...
            TypeError: matrix must have diagonal entries (=1, 7) congruent to 1 modulo 5, and lower left entry (=6) divisible by 5
        """
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        x = CongruenceSubgroupElement(self, x, check=check)
        if not check:
            return x

        a = x.a()
        c = x.c()
        d = x.d()
        N = self.level()
        if (a%N == 1) and (c%N == 0) and (d%N == 1):
            return x
        else:
            raise TypeError, "matrix must have diagonal entries (=%s, %s) congruent to 1 modulo %s, and lower left entry (=%s) divisible by %s" %(a, d, N, c, N)

def GammaH(level, H):
    r"""
    Return the congruence subgroup $\Gamma_H(N)$.

    INPUT:
        level -- an integer
        H -- either 0, 1, or a list
             * If H is a list, return $\Gamma_H(N)$, where $H$
               is the subgroup of $(\Z/N\Z)^*$ *generated* by the
               elements of the list.
             * If H = 0, returns $\Gamma_0(N)$.
             * If H = 1, returns $\Gamma_1(N)$.

    EXAMPLES:
        sage: GammaH(11,0)
        Congruence Subgroup Gamma0(11)
        sage: GammaH(11,1)
        Congruence Subgroup Gamma1(11)
        sage: GammaH(11,[2])
        Congruence Subgroup Gamma_H(11) with H generated by [2]
        sage: GammaH(11,[2,1])
        Congruence Subgroup Gamma_H(11) with H generated by [2]
    """
    if H == 0:
        return Gamma0(level)
    elif H == 1:
        return Gamma1(level)
    return GammaH_class(level, H)

def is_GammaH(x):
    """
    Return True if x is a congruence subgroup of type GammaH.

    EXAMPLES:
        sage: is_GammaH(GammaH(13, [2]))
        True
        sage: is_GammaH(Gamma0(6))
        False
    """
    return isinstance(x, GammaH_class)

class GammaH_class(CongruenceSubgroup):
    r"""
    The congruence subgroup $\Gamma_H(N)$.
    """
    def __init__(self, level, H):
        CongruenceSubgroup.__init__(self, level)
        if not isinstance(H, list):
            raise TypeError, "H must be a list."
        H = list(set(H))
        H.sort()
        if 1 in H:
            H.remove(1)
        self.__H = H

    def restrict(self, M):
        r"""
        Return the subgroup of $\Gamma_0(M)$ obtained by taking $H$ to
        be the image of the $H$ at level $N$ modulo $M$.

        EXAMPLES:
            sage: G = GammaH(33,[2])
            sage: G.restrict(11)
            Congruence Subgroup Gamma_H(11) with H generated by [2]
            sage: G.restrict(1)
            Congruence Subgroup Gamma_H(1) with H generated by []
            sage: G.restrict(15)
            Traceback (most recent call last):
            ...
            ValueError: M (=15) must be a divisor of the level (33) of self
        """
        N = self.level()
        M = ZZ(M)
        if N % M:
            raise ValueError, "M (=%s) must be a divisor of the level (%s) of self"%(M, N)
        if N == M:
            return self
        v = self.__H
        w = [a % M for a in v if a%M]
        return GammaH(M, w)

    def divisor_subgroups(self):
        r"""
        Given this congruence subgroup $\Gamma_H(N)$, return all
        subgroups $\Gamma_G(M)$ for $M$ a divisor of $N$ and such that
        $G$ is equal to the image of $H$ modulo $M$.

        EXAMPLES:
            sage: G = GammaH(33,[2]); G
            Congruence Subgroup Gamma_H(33) with H generated by [2]
            sage: G._list_of_elements_in_H()
            [1, 2, 4, 8, 16, 17, 25, 29, 31, 32]
            sage: G.divisor_subgroups()
            [Congruence Subgroup Gamma_H(1) with H generated by [],
             Congruence Subgroup Gamma_H(3) with H generated by [2],
             Congruence Subgroup Gamma_H(11) with H generated by [2],
             Congruence Subgroup Gamma_H(33) with H generated by [2]]
        """
        v = self.__H
        ans = []
        for M in divisors(self.level()):
            w = [a % M for a in v if a%M]
            ans.append(GammaH(M, w))
        return ans

    def __cmp__(self, other):
        if not isinstance(other, CongruenceSubgroup):
            return cmp(type(self), type(other))
        if isinstance(other, GammaH_class):
            c = cmp(self.level(), other.level())
            if c: return c
            return cmp(self._list_of_elements_in_H(), other._list_of_elements_in_H())
        return cmp(type(self), type(other))

    def _generators_for_H(self):
        return self.__H

    def _repr_(self):
        return "Congruence Subgroup Gamma_H(%s) with H generated by %s"%(self.level(), self.__H)

    def _latex_(self):
        return "\\Gamma_H(%s)"%self.level()

    def _list_of_elements_in_H(self):
        """
        Returns a sorted list of Python ints that are representatives
        between 1 and N-1 of the elements of H.

        WARNING: Do not change this returned list.

        EXAMPLES:
            sage: G = GammaH(11,[3]); G
            Congruence Subgroup Gamma_H(11) with H generated by [3]
            sage: G._list_of_elements_in_H()
            [1, 3, 4, 5, 9]
        """
        try:
            return self.__list_of_elements_in_H
        except AttributeError:
            pass
        N = self.level()
        if N == 1:
            self.__list_of_elements_in_H = [1]
            return [1]
        gens = self.__H

        H = set([1])
        N = int(N)
        for g in gens:
            if sage.rings.arith.gcd(g, N) != 1:
                raise ValueError, "gen (=%s) is not in (Z/%sZ)^*"%(g,N)
            gk = int(g) % N
            sbgrp = [gk]
            while not (gk in H):
                gk = (gk * g)%N
                sbgrp.append(gk)
            H = set([(x*h)%N for x in sbgrp for h in H])
        H = list(H)
        H.sort()
        self.__list_of_elements_in_H = H
        return H

    def is_even(self):
        """
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES:
            sage: GammaH(10, [3]).is_even()
            True
            sage: GammaH(14, [1]).is_even()
            False
        """
        if self.level() == 1:
            return True
        v = self._list_of_elements_in_H()
        return int(self.level() - 1) in v

    def generators(self):
        r"""
        Return generators for this congruence subgroup.

        The result is cached.

        EXAMPLE:
            sage: for g in GammaH(3, [2]).generators():
            ...     print g
            ...     print '---'
            [1 1]
            [0 1]
            ---
            [-1  0]
            [ 0 -1]
            ---
            [ 1 -1]
            [ 0  1]
            ---
            [ 1 -1]
            [ 3 -2]
            ---
            [ 2 -1]
            [ 3 -1]
            ---
            [-2  1]
            [-3  1]
            ---

        """
        try:
            return self.__gens
        except AttributeError:
            from sage.modular.modsym.ghlist import GHlist
            from congroup_pyx import generators_helper
            level = self.level()
            gen_list = generators_helper(GHlist(self), level, Mat2Z)
            self.__gens = [self(g, check=False) for g in gen_list]
            return self.__gens

    def _coset_reduction_data_first_coord(G):
        """
        Compute data used for determining the canonical coset
        representative of an element of SL_2(Z) modulo G.  This function
        specfically returns data needed for the first part of the
        reduction step (the first coordinate).

        INPUT:
            G -- a congruence subgroup Gamma_0(N), Gamma_1(N), or Gamma_H(N).

        OUTPUT:
            A list v such that
                v[u] = (min(u*h: h in H),
                        gcd(u,N) ,
                        an h such that h*u = min(u*h: h in H)).

        EXAMPLES:
            sage: G = GammaH(12,[-1,5]); G
            Congruence Subgroup Gamma_H(12) with H generated by [-1, 5]
            sage: G._coset_reduction_data_first_coord()
            [(0, 12, 0), (1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1), (1, 1, 5), (6, 6, 1),
            (1, 1, 7), (4, 4, 5), (3, 3, 7), (2, 2, 5), (1, 1, 11)]
        """
        H = G._list_of_elements_in_H()
        N = int(G.level())
        SQN = int(sqrt(N))

        # Get some useful fast functions for inverse and gcd
        inverse_mod = get_inverse_mod(N)   # optimal gcd function
        gcd = get_gcd(N)   # optimal gcd function

        # We will be adding to this list below.
        reduct_data = [(0,0,N,0)] + [ (u, 1, 1, inverse_mod(u, N)) for u in H ]
        not_yet_done = list(set(range(1,N)).difference(set(H)))
        not_yet_done.sort()

        # Make a table of the reduction of H (mod N/d),
        # one for each divisor d <= sqrt(N).
        repr_H_mod_N_over_d = {}
        for d in divisors(N):
            N_over_d = N//d
            w = [1]
            z = [1]
            for x in H:
                if not (x % N_over_d  in  w):
                    w.append(x % N_over_d)
                    z.append(x)
            repr_H_mod_N_over_d[d] = z

        # Compute the rest of the tuples
        while len(not_yet_done) > 0:
            u = not_yet_done[0]
            d = gcd(u, N)
            z = repr_H_mod_N_over_d[d]
            v = [((u*x) % N,  u, d, inverse_mod(x, N)) for x in z]
            reduct_data += v
            not_yet_done = list(set(not_yet_done).difference(set([a[0] for a in v])))
            not_yet_done.sort()

        # delete first entry.
        reduct_data.sort()
        reduct_data = [(a,b,c) for _,a,b,c in reduct_data]
        return reduct_data

    def _coset_reduction_data_second_coord(G):
        """
        Compute data used for determining the canonical coset
        representative of an element of SL_2(Z) modulo G.  This function
        specfically returns data needed for the first part of the
        reduction step (the first coordinate).

        INPUT:
            G -- a congruence subgroup Gamma_0(N), Gamma_1(N), or Gamma_H(N).

        OUTPUT:
            dictionary v with keys the divisors of N that are < sqrt(N)
            such that v[d] is the subgroup {h in H : h = 1 (mod N/d)}.

        EXAMPLES:
            sage: G = GammaH(1200,[-1,7]); G
            Congruence Subgroup Gamma_H(1200) with H generated by [-1, 7]
            sage: G._coset_reduction_data_second_coord()
            {1: [1], 2: [1], 3: [1], 4: [1], 5: [1], 6: [1], 8: [1], 10: [1], 12: [1],
            15: [1], 16: [1], 20: [1], 24: [1, 1151], 25: [1, 49], 30: [1]}
        """
        H = G._list_of_elements_in_H()
        N = int(G.level())
        SQN = int(sqrt(N))
        v = { 1: [1] }
        for d in [x for x in divisors(N) if x > 1 and x <= SQN]:
            N_over_d = N // d
            v[d] = [x for x in H if x % N_over_d == 1]
        return v

    def _coset_reduction_data(self):
        """
        Compute data used for determining the canonical coset
        representative of an element of SL_2(Z) modulo G.

        EXAMPLES:
            sage: G = GammaH(12,[-1,7]); G
            Congruence Subgroup Gamma_H(12) with H generated by [-1, 7]
            sage: G._coset_reduction_data()
            ([(0, 12, 0), (1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1), (1, 1, 5), (6, 6, 1),
              (1, 1, 7), (4, 4, 5), (3, 3, 7), (2, 2, 5), (1, 1, 11)],
              {1: [1], 2: [1, 7], 3: [1, 5]})
        """
        try:
            return self.__coset_reduction_data
        except AttributeError:
            pass
        self.__coset_reduction_data = (self._coset_reduction_data_first_coord(),
                                       self._coset_reduction_data_second_coord())
        return self.__coset_reduction_data


    def _reduce_coset(self, uu, vv):
        """
        Compute a canonical form for a given Manin symbol.

        INPUT:
        Two integers (uu,vv) that define an element of $(Z/NZ)^2$.
            uu -- an integer
            vv -- an integer

        OUTPUT:
           pair of integers that are equivalent to (uu,vv).

        NOTE: We do *not* require that gcd(uu,vv,N) = 1.  If the gcd is
        not 1, we return (0,0).

        EXAMPLE:
        An example at level 9.
            sage: G = GammaH(9,[7]); G
            Congruence Subgroup Gamma_H(9) with H generated by [7]
            sage: a = []
            sage: for i in range(G.level()):
            ...     for j in range(G.level()):
            ...       a.append(G._reduce_coset(i,j))
            sage: v = list(set(a))
            sage: v.sort()
            sage: v
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (3, 1), (3, 2), (6, 1), (6, 2)]

        An example at level 100.
            sage: G = GammaH(100,[3,7]); G
            Congruence Subgroup Gamma_H(100) with H generated by [3, 7]
            sage: a = []
            sage: for i in range(G.level()):
            ...   for j in range(G.level()):
            ...       a.append(G._reduce_coset(i,j))
            sage: v = list(set(a))
            sage: v.sort()
            sage: len(v)
            361
        """
        N = int(self.level())
        u = uu % N
        v = vv % N
        first, second = self._coset_reduction_data()
        gcd_u_N = first[u][1]
        gcd_v_N = first[v][1]
        if arith.gcd(gcd_u_N, gcd_v_N) != 1:
            return (0, 0)

        if u == 0:
            return (0, first[v][0])
        if v == 0:
            return (first[u][0], 0)

        def f(u,v):
            v = (v*first[u][2]) % N
            d = first[u][1]   # gcd(u,N)
            u = first[u][0]
            H_cong1_mod_N_over_d = second[d]
            if len(H_cong1_mod_N_over_d) > 1:
                vmin = v
                for x in H_cong1_mod_N_over_d:
                    v1 = (v*x) % N
                    if v1 < vmin:
                        vmin = v1
                    v = vmin
            return u,v

        if gcd_u_N <= gcd_v_N:
            (u,v) = f(u,v)
        else:
            (v,u) = f(v,u)

        return (u,v)


    def _reduce_cusp(self, c):
        """
        Compute a canonical form for the cusp $uu/vv$.

        INPUT:
            cusp
        OUTPUT:
            cusp


        OUTPUT:
            bool -- True if self and other are equivalent
            int -- $1$, $0$, or $-1$,
                   If the two cusps are $u1/v1$ and $u2/v2$, then they are
                   equivalent modulo Gamma_H(N) if and only if
                        $v1 =  h*v2 (mod N)$ and $u1 =  h^(-1)*u2 (mod gcd(v1,N))$
                   or
                        $v1 = -h*v2 (mod N)$ and $u1 = -h^(-1)*u2 (mod gcd(v1,N))$
                   where $h \in H$.   In the first case we return $1$, and in
                   the second we return $-1$.   We return $0$ if the two
                   cusps are not equivalent.

        EXAMPLE:
        """
        from sage.rings.integer import Integer
        N = int(self.level())
        Cusps = c.parent()
        u = int(c.numerator() % N)
        v = int(c.denominator() % N)
        first, second = self._coset_reduction_data()
        gcd_u_N = first[u][1]
        gcd_v_N = first[v][1]

        if u == 0:
            return Cusps(0), 1
        if v == 0:
            return Cusps((1,0)), 1

        def f(u,v):
            h = first[v][2]
            v = first[v][0]
            hinv = int(Integer(h).inverse_mod(gcd_v_N))     # optimize
            d = first[v][1]   # d = gcd(v,N)
            if d == 1:
                return 0, 1, h
            u = (hinv * u) % d
            H_cong1_mod_N_over_d = second[d]
            if len(H_cong1_mod_N_over_d) > 1:
                umin = u
                for x in H_cong1_mod_N_over_d:
                    u1 = (u*x) % d
                    if u1 < umin:
                        umin = u1
                    u = umin
            return u,v, h

        u, v, h = f(u,v)
        N_over_2 = N//2
        if u > N_over_2:
            u = N_over_2 - u
            v = N_over_2 - v
            t = -1
        else:
            t = 1

        return Cusps((u,v)), t



    def __call__(self, x, check=True):
        r"""
        Create an element of this congruence subgroup from x.

        If the optional flag check is True (default), check whether
        x actually gives an element of self.

        EXAMPLES:
            sage: G = GammaH(10, [3])
            sage: G([1, 0, -10, 1])
            [ 1   0]
            [-10  1]
            sage: G(matrix(ZZ, 2, [7, 1, 20, 3]))
            [ 7  1]
            [20  3]
            sage: GammaH(10, [9])([7, 1, 20, 3])
            Traceback (most recent call last):
            ...
            TypeError: matrix must have lower right entry (=3) congruent modulo 10 to some element of H
        """
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        x = CongruenceSubgroupElement(self, x, check=check)
        if not check:
            return x

        c = x.c()
        d = x.d()
        N = self.level()
        if c%N != 0:
            raise TypeError, "matrix must have lower left entry (=%s) divisible by %s" %(c, N)
        elif d%N in self._list_of_elements_in_H():
            return x
        else:
            raise TypeError, "matrix must have lower right entry (=%s) congruent modulo %s to some element of H" %(d, N)


import congroup_pyx
degeneracy_coset_representatives_gamma0 = congroup_pyx.degeneracy_coset_representatives_gamma0
degeneracy_coset_representatives_gamma1 = congroup_pyx.degeneracy_coset_representatives_gamma1



## def xxx_degeneracy_coset_representatives_gamma0(N, M, t):
##     r"""
##     Let $N$ be a positive integer and $M$ a multiple of $N$.  Let $t$ be a
##     divisor of $N/M$, and let $T$ be the $2x2$ matrix $T=[0,0; 0,t]$.
##     This function returns representatives for the orbit set
##     $\Gamma_0(N) \backslash T \Gamma_0(M)$, where $\Gamma_0(N)$
##     acts on the left on $T \Gamma_0(M)$.

##     INPUT:
##         N -- int
##         M -- int (divisor of N)
##         t -- int (divisors of N/M)

##     OUTPUT:
##         list -- list of lists [a,b,c,d], where [a,b,c,d] should
##                 be viewed as a 2x2 matrix.

##     This function is used for computation of degeneracy maps between
##     spaces of modular symbols, hence its name.

##     We use that $T^{-1}*(a,b;c,d)*T = (a,bt,c/t,d),$
##     that the group $T^{-1}Gamma_0(N) T$ is contained in $\Gamma_0(M)$,
##     and that $\Gamma_0(N) T$ is contained in $T \Gamma_0(M)$.

##     ALGORITHM:
##     \begin{enumerate}
##     \item Compute representatives for $\Gamma_0(N/t,t)$ inside of $\Gamma_0(M)$:
##           COSET EQUIVALENCE:
##            Two right cosets represented by $[a,b;c,d]$ and
##            $[a',b';c',d']$ of $\Gamma_0(N/t,t)$ in ${\rm SL}_2(\Z)$
##            are equivalent if and only if
##            $(a,b)=(a',b')$ as points of $\P^1(\Z/t\Z)$,
##            i.e., $ab' \con ba' \pmod{t}$,
##            and $(c,d) = (c',d')$ as points of $\P^1(\Z/(N/t)\Z)$.

##         ALGORITHM to list all cosets:
##         \begin{enumerate}
##            \item Compute the number of cosets.
##            \item Compute a random element x of Gamma_0(M).
##            \item Check if x is equivalent to anything generated so
##                far; if not, add x to the list.
##            \item Continue until the list is as long as the bound
##                computed in A.
##         \end{enumerate}

##     \item There is a bijection between $\Gamma_0(N)\backslash T \Gamma_0(M)$
##           and $\Gamma_0(N/t,t) \Gamma_0(M)$ given by $T r$ corresponds to $r$.
##           Consequently we obtain coset representatives for
##           $\Gamma_0(N)\backslash T \Gamma_0(M)$ by left multiplying by $T$ each
##           coset representative of $\Gamma_0(N/t,t) \Gamma_0(M)$ found
##           in step 1.
##     \end{enumerate}
##     """
##     import sage.modular.dims as dims
##     N = int(N); M = int(M); t = int(t)
##     if N % M != 0:
##         raise ArithmeticError, "M (=%s) must be a divisor of N (=%s)"%(M,N)
##     n     = dims.idxG0(N) // dims.idxG0(M)          # number of coset representatives.
##     Ndivd = N // t
##     R     = []                        # coset reps found so far.
##     halfmax = 2*(n+10)
##     while len(R) < n:
##         # try to find another coset representative.
##         cc = M*random.randrange(-halfmax, halfmax+1)
##         dd =   random.randrange(-halfmax, halfmax+1)
##         g, bb, aa = arith.xgcd(-cc,dd)
##         if g == 0: continue
##         cc //= g
##         dd //= g
##         if cc % M != 0: continue
##         # Test if we've found a new coset representative.
##         is_new = True
##         for r in R:
##             if ((r[1]*aa - r[0]*bb) % t == 0) and \
##                      ((r[3]*cc - r[2]*dd) % Ndivd == 0):
##                 is_new = False
##                 break
##         # If our matrix is new add it to the list.
##         if is_new: R.append([aa,bb,cc,dd])
##     # Return the list left multiplied by T.
##     return [[r[0], r[1], r[2]*t, r[3]*t] for r in R]
