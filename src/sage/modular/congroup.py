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

import sage.groups.group as group

import sage.rings.arith as arith
from sage.rings.integer_mod_ring import IntegerModRing
from sage.rings.all import QQ, ZZ, divisors

from sage.matrix.matrix_space import MatrixSpace

from congroup_element import CongruenceSubgroupElement

import sage.modular.modsym.p1list

# Just for now until we make an SL_2 group type.
from sage.matrix.matrix_space import MatrixSpace
Mat2Z = MatrixSpace(ZZ,2)

def is_CongruenceSubgroup(x):
    """
    Return True if x is of type CongruenceSubgroup.

    EXAMPLES:
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

class CongruenceSubgroup(group.Group):
    def __init__(self, level):
        """
        Create a congruence subgroup with given level.

        EXAMPLES:
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

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5)._repr_()
            'Generic congruence subgroup'
        """
        return "Generic congruence subgroup"

    def __reduce__(self):
        """
        Used for pickling self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: all subclasses must define a __reduce__ method
        """
        raise NotImplementedError, "all subclasses must define a __reduce__ method"

    def __hash__(self):
        """
        Return a hash of self.

        EXAMPLES:
            sage: Gamma0(11).__hash__()
            -545929996 # 32-bit
            466678398374495476 # 64-bit
            sage: Gamma1(11).__hash__()
            -830809815 # 32-bit
            4909638266971150633 # 64-bit
        """
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
        r"""
        Determine whether $x$ and $y$ are equivalent by an element of
        self, i.e. whether or not there exists an element $g$ of self
        such that $g\cdot x = y$.

        NOTE: This function must be overridden by all subclasses.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).are_equivalent(0, 0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def coset_reps(self):
        """
        Return coset representatives for this congruence subgroup.

        NOTE: This function must be overridden by all subclasses.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).coset_reps()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def generators(self):
        """
        Return generators for this congruence subgroup.

        NOTE: This function must be overridden by all subclasses.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).generators()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def gens(self):
        """
        Return a tuple of generators for this congruence subgroup.

        The generators need not be minimal.

        EXAMPLES:
            sage: SL2Z.gens()
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
            sage: SL2Z.gen(1)
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
            sage: SL2Z.ngens()
            2
        """
        return len(self.generators())

    def level(self):
        """
        Return the level of this congruence subgroup.

        EXAMPLES:
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

    def __cmp__(self, right):
        """
        Compare self to right.

        NOTE: This function must be overridden by all subclasses.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).__cmp__(ZZ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_abelian(self):
        """
        Return True if this congruence subgroup is abelian.

        Since congruence subgroups are always nonabelian, this always
        returns False.

        EXAMPLES:
            sage: SL2Z.is_abelian()
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
            sage: SL2Z.is_finite()
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
        """
        Return True if self is a subgroup of right, and False
        otherwise.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).is_subgroup(SL2Z)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def is_odd(self):
        """
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES:
            sage: SL2Z.is_odd()
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
            sage: SL2Z.is_even()
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
        r"""
        Return a new group of the same type (Gamma0, Gamma1, or
        GammaH) as self of the given level. In the case that self is
        of type GammaH, we take the largest H inside
        $(Z/\text{level}Z)^\times$ which maps to H, namely its inverse
        image under the natural reduction map.

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
            sage: SL2Z.order()
            +Infinity
            sage: Gamma0(5).order()
            +Infinity
            sage: Gamma1(2).order()
            +Infinity
            sage: GammaH(12, [5]).order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def __call__(self, x, check=True):
        """
        Coerce x into self.

        NOTE: This function should be overridden by any subclass the
        user will interact with directly.

        EXAMPLES:
            sage: sage.modular.congroup.CongruenceSubgroup(5).__call__(0)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if isinstance(x, CongruenceSubgroupElement) and x.parent() == self:
            return x
        raise NotImplementedError

def lift_to_sl2z(c, d, N):
    """
    Given a vector (c, d) in (Z/NZ)^2, this function computes and
    returns a list [a, b, c', d'] that defines a 2x2 matrix with
    determinant 1 and integer entries, such that c=c'(mod N) and
    d=d'(mod N).

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
        sage: is_Gamma0(SL2Z)
        True
        sage: is_Gamma0(Gamma0(13))
        True
        sage: is_Gamma0(Gamma1(6))
        False
    """
    return isinstance(x, Gamma0_class)

_gamma0_cache = {}
def Gamma0(N):
    """
    Return the congruence subgroup Gamma0(N).

    EXAMPLES:
        sage: G = Gamma0(51) ; G
        Congruence Subgroup Gamma0(51)
        sage: G == Gamma0(51)
        True
        sage: G is Gamma0(51)
        True
    """
    try:
        return _gamma0_cache[N]
    except KeyError:
        _gamma0_cache[N] = Gamma0_class(N)
        return _gamma0_cache[N]

class Gamma0_class(CongruenceSubgroup):
    def __init__(self, level):
        r"""
        The congruence subgroup $\Gamma_0(N)$.

        EXAMPLES:
            sage: G = Gamma0(11); G
            Congruence Subgroup Gamma0(11)
            sage: loads(G.dumps()) == G
            True
        """
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: Gamma0(98)._repr_()
            'Congruence Subgroup Gamma0(98)'
        """
        return "Congruence Subgroup Gamma0(%s)"%self.level()

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES:
            sage: Gamma0(22).__reduce__()
            (<function Gamma0 at ...>, (22,))
        """
        return Gamma0, (self.level(),)

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES:
            sage: Gamma0(20)._latex_()
            '\\Gamma_0(20)'
            sage: latex(Gamma0(20))
            \Gamma_0(20)
        """
        return "\\Gamma_0(%s)"%self.level()

    def _generators_for_H(self):
        """
        Return generators for the subgroup H of the units mod
        self.level() that defines self.

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

            sage: G = Gamma0(6)
            sage: G._list_of_elements_in_H()
            [1, 5]

            sage: G = Gamma0(1)
            sage: G._list_of_elements_in_H()
            [1]
        """
        N = self.level()
        if N != 1:
            gcd = arith.gcd
            return [ x for x in range(1, N) if gcd(x, N) == 1 ]
        else:
            return [1]

    def __cmp__(self, right):
        """
        Compare self to right.

        EXAMPLES:
            sage: Gamma0(21).__cmp__(Gamma0(21))
            0
            sage: Gamma0(21) < Gamma0(32)
            True
        """
        if not is_Gamma0(right):
            if is_CongruenceSubgroup(right):
                c = cmp(self.level(), right.level())
                if c: return c
            return cmp(type(self), type(right))
        return cmp(self.level(), right.level())

    def is_even(self):
        """
        Return True precisely if this subgroup contains the matrix -1.

        Since Gamma0(N) always, contains the matrix -1, this always
        returns True.

        EXAMPLES:
            sage: Gamma0(12).is_even()
            True
            sage: SL2Z.is_even()
            True
        """
        return True

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES:
            sage: G = Gamma0(20)
            sage: G.is_subgroup(SL2Z)
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
        if is_Gamma0(right):
            return self.level() % right.level() == 0
        if is_Gamma1(right):
            if right.level() >= 3:
                return False
            elif right.level() == 2:
                return self.level() == 2
            # case level 1 dealt with above
        raise NotImplementedError

    def coset_reps(self):
        r"""
        Return representatives for the right cosets of this congruence
        subgroup in ${\rm SL}_2(\Z)$ as a generator object.

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
        sage: is_SL2Z(SL2Z)
        True
        sage: is_SL2Z(Gamma0(6))
        False
    """
    return isinstance(x, SL2Z_class)

class SL2Z_class(Gamma0_class):
    def __init__(self):
        r"""
        The modular group ${\rm SL}_2(\Z)$.

        EXAMPLES:
            sage: G = SL2Z; G
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
            sage: SL2Z.0 * SL2Z.1
            [ 0 -1]
            [ 1  1]

            sage: SL2Z == loads(dumps(SL2Z))
            True
            sage: SL2Z is loads(dumps(SL2Z))
            True
        """
        Gamma0_class.__init__(self, 1)

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES:
            sage: SL2Z.__reduce__()
            (<function _SL2Z_ref at ...>, ())
        """
        return _SL2Z_ref, ()

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: SL2Z._repr_()
            'Modular Group SL(2,Z)'
        """
        return "Modular Group SL(2,Z)"

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES:
            sage: SL2Z._latex_()
            '\\mbox{\\rm SL}_2(\\mathbf{Z})'
            sage: latex(SL2Z)
            \mbox{\rm SL}_2(\mathbf{Z})
        """
        return "\\mbox{\\rm SL}_2(%s)"%(ZZ._latex_())

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES:
            sage: SL2Z.is_subgroup(SL2Z)
            True
            sage: SL2Z.is_subgroup(Gamma1(1))
            True
            sage: SL2Z.is_subgroup(Gamma0(6))
            False
        """
        return right.level() == 1

SL2Z = SL2Z_class()

def _SL2Z_ref():
    """
    Return SL2Z. (Used for pickling SL2Z.)

    EXAMPLES:
        sage: sage.modular.congroup._SL2Z_ref()
        Modular Group SL(2,Z)
        sage: sage.modular.congroup._SL2Z_ref() is SL2Z
        True
    """
    return SL2Z

def is_Gamma1(x):
    """
    Return True if x is a congruence subgroup of type Gamma1.

    EXAMPLES:
        sage: is_Gamma1(SL2Z)
        True
        sage: is_Gamma1(Gamma1(13))
        True
        sage: is_Gamma1(Gamma0(6))
        False
    """
    return (isinstance(x, Gamma1_class) or is_SL2Z(x))

_gamma1_cache = {}
def Gamma1(N):
    r"""
    Return the congruence subgroup $\Gamma_1(N)$.

    EXAMPLES:
        sage: Gamma1(5)
        Congruence Subgroup Gamma1(5)
        sage: G = Gamma1(23)
        sage: G is Gamma1(23)
        True
        sage: G == loads(dumps(G))
        True
        sage: G is loads(dumps(G))
        True
    """
    try:
        return _gamma1_cache[N]
    except KeyError:
        _gamma1_cache[N] = Gamma1_class(N)
        return _gamma1_cache[N]

class Gamma1_class(CongruenceSubgroup):
    def __init__(self, level):
        r"""
        The congruence subgroup $\Gamma_1(N)$.

        EXAMPLES:
            sage: G = Gamma1(11); G
            Congruence Subgroup Gamma1(11)
            sage: loads(G.dumps()) == G
            True
        """
        CongruenceSubgroup.__init__(self, level)

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: Gamma1(133)._repr_()
            'Congruence Subgroup Gamma1(133)'
        """
        return "Congruence Subgroup Gamma1(%s)"%self.level()

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES:
            sage: Gamma1(82).__reduce__()
            (<function Gamma1 at ...>, (82,))
        """
        return Gamma1, (self.level(),)

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES:
            sage: Gamma1(3)._latex_()
            '\\Gamma_1(3)'
            sage: latex(Gamma1(3))
            \Gamma_1(3)
        """
        return "\\Gamma_1(%s)"%self.level()

    def __cmp__(self, right):
        """
        Compare self to right.

        EXAMPLES:
            sage: G = Gamma1(111)
            sage: G.__cmp__(Gamma1(111))
            0
            sage: G.__cmp__(135) is not 0
            True
        """
        if not is_Gamma1(right):
            if is_CongruenceSubgroup(right):
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
            sage: Gamma1(3).is_subgroup(SL2Z)
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
        if is_Gamma0(right) or is_Gamma1(right):
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

_gammaH_cache = {}
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

    H = _normalize_H(H, level)
    key = (level, tuple(H))
    try:
        return _gammaH_cache[key]
    except KeyError:
        _gammaH_cache[key] = GammaH_class(level, H)
        return _gammaH_cache[key]

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

def _normalize_H(H, level):
    """
    Normalize representatives for a given subgroup H of the units
    modulo level.

    NOTE: This function does *not* make any attempt to find a minimal
    set of generators for H. It simply normalizes the inputs for use
    in hashing.

    EXAMPLES:
        sage: sage.modular.congroup._normalize_H([23], 10)
        [3]
        sage: sage.modular.congroup._normalize_H([1,5], 7)
        [5]
        sage: sage.modular.congroup._normalize_H([4,18], 14)
        [4]
        sage: sage.modular.congroup._normalize_H([-1,7,9], 10)
        [7, 9]
    """
    if not isinstance(H, list):
        raise TypeError, "H must be a list."
    H = list(set([h%level for h in H]))
    H.sort()
    if 1 in H:
        H.remove(1)
    return H

class GammaH_class(CongruenceSubgroup):
    def __init__(self, level, H):
        r"""
        The congruence subgroup $\Gamma_H(N)$. The subgroup H
        must be input as a list.

        EXAMPLES:
            sage: GammaH(117, [4])
            Congruence Subgroup Gamma_H(117) with H generated by [4]
            sage: G = GammaH(16, [7])
            sage: G == loads(dumps(G))
            True
            sage: G is loads(dumps(G))
            True
        """
        CongruenceSubgroup.__init__(self, level)
        self.__H = _normalize_H(H, level)

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

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES:
            sage: GammaH(92,[5,11]).__reduce__()
            (<function GammaH at ...>, (92, [5, 11]))
        """
        return GammaH, (self.level(), self.__H)

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
        """
        Compare self to right.

        EXAMPLES:
            sage: G = GammaH(86, [9])
            sage: G.__cmp__(G)
            0
            sage: G.__cmp__(GammaH(86, [11])) is not 0
            True
        """
        if not is_CongruenceSubgroup(other):
            return cmp(type(self), type(other))
        if is_GammaH(other):
            c = cmp(self.level(), other.level())
            if c: return c
            return cmp(self._list_of_elements_in_H(), other._list_of_elements_in_H())
        return cmp(type(self), type(other))

    def _generators_for_H(self):
        """
        Return generators for the subgroup H of the units mod
        self.level() that defines self.

        EXAMPLES:
            sage: GammaH(17,[4])._generators_for_H()
            [4]
            sage: GammaH(12,[-1])._generators_for_H()
            [11]
        """
        return self.__H

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES:
            sage: GammaH(123, [55])._repr_()
            'Congruence Subgroup Gamma_H(123) with H generated by [55]'
        """
        return "Congruence Subgroup Gamma_H(%s) with H generated by %s"%(self.level(), self.__H)

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES:
            sage: GammaH(3,[2])._latex_()
            '\\Gamma_H(3)'
        """
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
        representative of an element of SL_2(Z) modulo G. This
        function specifically returns data needed for the first part
        of the reduction step (the first coordinate).

        INPUT:
            G -- a congruence subgroup Gamma_0(N), Gamma_1(N), or Gamma_H(N).

        OUTPUT:
            A list v such that
                v[u] = (min(u*h: h in H),
                        gcd(u,N) ,
                        an h such that h*u = min(u*h: h in H)).

        EXAMPLES:
            sage: G = GammaH(12,[-1,5]); G
            Congruence Subgroup Gamma_H(12) with H generated by [5, 11]
            sage: G._coset_reduction_data_first_coord()
            [(0, 12, 0), (1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1), (1, 1, 5), (6, 6, 1),
            (1, 1, 7), (4, 4, 5), (3, 3, 7), (2, 2, 5), (1, 1, 11)]
        """
        H = [ int(x) for x in G._list_of_elements_in_H() ]
        N = int(G.level())

        # Get some useful fast functions for inverse and gcd
        inverse_mod = arith.get_inverse_mod(N)   # optimal inverse function
        gcd = arith.get_gcd(N)   # optimal gcd function

        # We will be filling this list in below.
        reduct_data = [0] * N

        # We can fill in 0 and all elements of H immediately
        reduct_data[0] = (0,N,0)
        for u in H:
            reduct_data[u] = (1, 1, inverse_mod(u, N))

        # Make a table of the reduction of H (mod N/d), one for each
        # divisor d.
        repr_H_mod_N_over_d = {}
        for d in divisors(N):
            # We special-case N == d because in this case,
            # 1 % N_over_d is 0
            if N == d:
                repr_H_mod_N_over_d[d] = [1]
                break
            N_over_d = N//d
            # For each element of H, we look at its image mod
            # N_over_d. If we haven't yet seen it, add it on to
            # the end of z.
            w = [0] * N_over_d
            z = [1]
            for x in H:
                val = x%N_over_d
                if not w[val]:
                    w[val] = 1
                    z.append(x)
            repr_H_mod_N_over_d[d] = z

        # Compute the rest of the tuples. The values left to process
        # are those where reduct_data has a 0. Note that several of
        # these values are processed on each loop below, so re-index
        # each time.
        while True:
            try:
                u = reduct_data.index(0)
            except ValueError:
                break
            d = gcd(u, N)
            for x in repr_H_mod_N_over_d[d]:
                reduct_data[(u*x)%N] = (u, d, inverse_mod(x,N))

        return reduct_data

    def _coset_reduction_data_second_coord(G):
        """
        Compute data used for determining the canonical coset
        representative of an element of SL_2(Z) modulo G. This
        function specifically returns data needed for the second part
        of the reduction step (the second coordinate).

        INPUT:
            self

        OUTPUT:
            a dictionary v with keys the divisors of N such that v[d]
            is the subgroup {h in H : h = 1 (mod N/d)}.

        EXAMPLES:
            sage: G = GammaH(240,[7,239])
            sage: G._coset_reduction_data_second_coord()
            {1: [1], 2: [1], 3: [1], 4: [1], 5: [1, 49], 6: [1], 48: [1, 191], 8: [1], 80: [1, 7, 49, 103], 10: [1, 49], 12: [1], 15: [1, 49], 240: [1, 7, 49, 103, 137, 191, 233, 239], 40: [1, 7, 49, 103], 20: [1, 49], 24: [1, 191], 120: [1, 7, 49, 103, 137, 191, 233, 239], 60: [1, 49, 137, 233], 30: [1, 49, 137, 233], 16: [1]}
            sage: G = GammaH(1200,[-1,7]); G
            Congruence Subgroup Gamma_H(1200) with H generated by [7, 1199]
            sage: K = G._coset_reduction_data_second_coord().keys() ; K.sort()
            sage: K == divisors(1200)
            True
        """
        H = G._list_of_elements_in_H()
        N = G.level()
        v = { 1: [1] , N: H }
        for d in [x for x in divisors(N) if x > 1 and x < N ]:
            N_over_d = N // d
            v[d] = [x for x in H if x % N_over_d == 1]
        return v

    def _coset_reduction_data(self):
        """
        Compute data used for determining the canonical coset
        representative of an element of SL_2(Z) modulo G.

        EXAMPLES:
            sage: G = GammaH(12,[-1,7]); G
            Congruence Subgroup Gamma_H(12) with H generated by [7, 11]
            sage: G._coset_reduction_data()
            ([(0, 12, 0), (1, 1, 1), (2, 2, 1), (3, 3, 1), (4, 4, 1), (1, 1, 5), (6, 6, 1), (1, 1, 7), (4, 4, 5), (3, 3, 7), (2, 2, 5), (1, 1, 11)],
            {1: [1], 2: [1, 7], 3: [1, 5],  4: [1, 7], 6: [1, 5, 7, 11], 12: [1, 5, 7, 11]})
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

        This demonstrates the problem underlying trac #1220:
            sage: G = GammaH(99, [67])
            sage: G._reduce_coset(11,-3)
            (11, 96)
            sage: G._reduce_coset(77, -3)
            (11, 96)
        """
        N = int(self.level())
        u = uu % N
        v = vv % N
        first, second = self._coset_reduction_data()

        if arith.gcd(first[u][1], first[v][1]) != 1:
            return (0,0)
        if not u:
            return (0, first[v][0])
        if not v:
            return (first[u][0], 0)

        new_u = first[u][0]
        d = first[u][1]
        new_v = (first[u][2] * v) % N
        H_ls = second[d]
        if len(H_ls) > 1:
            new_v = min([ (new_v * h)%N for h in H_ls ])

        return (new_u, new_v)

    def _reduce_cusp(self, c):
        r"""
        Compute a minimal representative for the given cusp c.
        Returns a pair (c', t), where c' is the minimal representative
        for the given cusp, and t is either 1 or -1, as explained
        below.

        The minimal representative for a cusp is the element in $P^1(Q)$
        in lowest terms with minimal positive denominator, and minimal
        positive numerator for that denominator.

        Two cusps $u1/v1$ and $u2/v2$ are equivalent modulo $\Gamma_H(N)$
        if and only if
            $v1 =  h*v2 (mod N)$ and $u1 =  h^(-1)*u2 (mod gcd(v1,N))$
        or
            $v1 = -h*v2 (mod N)$ and $u1 = -h^(-1)*u2 (mod gcd(v1,N))$
        for some $h \in H$. Then t is 1 or -1 as c and c' fall into
        the first or second case, respectively.

        EXAMPLES:
            sage: GammaH(6,[5])._reduce_cusp(Cusp(5,3))
            (1/3, -1)
            sage: GammaH(12,[5])._reduce_cusp(Cusp(8,9))
            (1/3, -1)
            sage: GammaH(12,[5])._reduce_cusp(Cusp(5,12))
            (Infinity, 1)
            sage: GammaH(12,[])._reduce_cusp(Cusp(5,12))
            (5/12, 1)
            sage: GammaH(21,[5])._reduce_cusp(Cusp(-9/14))
            (1/7, 1)
        """
        N = int(self.level())
        Cusps = c.parent()
        v = int(c.denominator() % N)
        H = self._list_of_elements_in_H()

        # First, if N | v, take care of this case. If u is in \pm H,
        # then we return Infinity. If not, let u_0 be the minimum
        # of \{ h*u | h \in \pm H \}. Then return u_0/N.
        if not v:
            u = c.numerator() % N
            if u in H:
                return Cusps((1,0)), 1
            if (N-u) in H:
                return Cusps((1,0)), -1
            ls = [ (u*h)%N for h in H ]
            m1 = min(ls)
            m2 = N-max(ls)
            if m1 < m2:
                return Cusps((m1,N)), 1
            else:
                return Cusps((m2,N)), -1

        u = int(c.numerator() % v)
        gcd = arith.get_gcd(N)
        d = gcd(v,N)

        # If (N,v) == 1, let v_0 be the minimal element
        # in \{ v * h | h \in \pm H \}. Then we either return
        # Infinity or 1/v_0, as v is or is not in \pm H,
        # respectively.
        if d == 1:
            if v in H:
                return Cusps((0,1)), 1
            if (N-v) in H:
                return Cusps((0,1)), -1
            ls = [ (v*h)%N for h in H ]
            m1 = min(ls)
            m2 = N-max(ls)
            if m1 < m2:
                return Cusps((1,m1)), 1
            else:
                return Cusps((1,m2)), -1

        val_min = v
        inv_mod = arith.get_inverse_mod(N)

        # Now we're in the case (N,v) > 1. So we have to do several
        # steps: first, compute v_0 as above. While computing this
        # minimum, keep track of *all* pairs of (h,s) which give this
        # value of v_0.
        hs_ls = [(1,1)]
        for h in H:
            tmp = (v*h)%N

            if tmp < val_min:
                val_min = tmp
                hs_ls = [(inv_mod(h,N), 1)]
            elif tmp == val_min:
                hs_ls.append((inv_mod(h,N), 1))

            if (N-tmp) < val_min:
                val_min = N - tmp
                hs_ls = [(inv_mod(h,N), -1)]
            elif (N-tmp) == val_min:
                hs_ls.append((inv_mod(h,N), -1))

        # Finally, we find our minimal numerator. Let u_1 be the
        # minimum of s*h^-1*u mod d as (h,s) ranges over the elements
        # of hs_ls. We must find the smallest integer u_0 which is
        # smaller than v_0, congruent to u_1 mod d, and coprime to
        # v_0. Then u_0/v_0 is our minimal representative.
        u_min = val_min
        sign = None
        for h_inv,s in hs_ls:
            tmp = (h_inv * s * u)%d
            while gcd(tmp, val_min) > 1 and tmp < u_min:
                tmp += d
            if tmp < u_min:
                u_min = tmp
                sign = s

        return Cusps((u_min, val_min)), sign

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
