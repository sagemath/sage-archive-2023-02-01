r"""
Finite subgroups of modular abelian varieties

\sage can compute with fairly general finite subgroups of modular
abelian varieties.  Elements of finite order are represented by
equivalence classes of elements in $H_1(A,\QQ)$ modulo $H_1(A,\ZZ)$.
A finite subgroup can be defined by giving generators and via various
other constructions.  Given a finite subgroup, one can compute
generators, the structure as an abstract group.  Arithmetic on
subgroups is also supported, including adding two subgroups together,
checking inclusion, etc.

TODO: Intersection, action of Hecke operators.

AUTHOR:
    -- William Stein (2007-03)

EXAMPLES:
    sage: J = J0(33)
    sage: C = J.cuspidal_subgroup()
    sage: C
    Cuspidal subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(33)
    sage: C.order()
    100
    sage: C.gens()
    [[(1/10, 0, 1/10, 1/10, 1/10, -7/10)], [(0, 1/5, 1/10, 0, 1/10, -1/10)]]
    sage: C.0 + C.1
    [(1/10, 1/5, 1/5, 1/10, 1/5, -4/5)]
    sage: 10*(C.0 + C.1)
    [(0, 0, 0, 0, 0, 0)]
    sage: G = C.subgroup([C.0 + C.1]); G
    Finite subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(33) defined over Complex Field with 53 bits of precision
    sage: G.gens()
    [[(1/10, 1/5, 1/5, 1/10, 1/5, -4/5)]]
    sage: G.order()
    10
    sage: G <= C
    True
    sage: G >= C
    False


We make a table of the order of the cuspidal subgroup for the first
few levels:

    sage: for N in range(11,40): print N, J0(N).cuspidal_subgroup().order()
    ...
    11 5
    12 1
    13 1
    14 6
    15 8
    16 1
    17 4
    18 1
    19 3
    20 6
    21 8
    22 25
    23 11
    24 8
    25 1
    26 21
    27 9
    28 36
    29 7
    30 192
    31 5
    32 8
    33 100
    34 48
    35 48
    36 12
    37 3
    38 135
    39 56
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.modules.module      import Module
from sage.structure.element   import ModuleElement
from sage.structure.sequence  import Sequence
from sage.rings.all           import gcd, lcm, QQ, ZZ, CC
from sage.misc.misc           import prod

class FiniteSubgroup(Module):
    def __init__(self, abvar, base_field=QQ):
        Module.__init__(self, base_field)
        self._abvar = abvar
        self._base_field = base_field

    def __cmp__(self, other):
        if not isinstance(other, FiniteSubgroup):
            return cmp(type(self), type(other))
        c = cmp(self.abelian_variety(), other.abelian_variety())
        if c: return c
        # Minus sign because order gets reversed in passing to lattices.
        return -cmp(self._full_module(), other._full_module())

    def is_subgroup(self, other):
        """
        Return True exactly if self is a subgroup of other, and both
        are defined as subgroups of the same ambient abelian variety.

        EXAMPLES:
            sage: C = J0(22).cuspidal_subgroup()
            sage: H = C.subgroup([C.0])
            sage: K = C.subgroup([C.1])
            sage: H.is_subgroup(K)
            False
            sage: K.is_subgroup(H)
            False
            sage: K.is_subgroup(C)
            True
            sage: H.is_subgroup(C)
            True
        """
        if not isinstance(other, FiniteSubgroup):
            return False
        if self.abelian_variety() != other.abelian_variety():
            return False
        if self._full_module().is_submodule(other._full_module()):
            return True
        return False

    def __add__(self, other):
        """
        Return the sum of two subgroups.

        EXAMPLES:
            sage: C = J0(22).cuspidal_subgroup()
            sage: C.gens()
            [[(1/5, 1/5, -1/5, 0)], [(0, 0, 0, 1/5)]]
            sage: A = C.subgroup([C.0]); B = C.subgroup([C.1])
            sage: A + B == C
            True
        """
        if not isinstance(other, FiniteSubgroup):
            raise TypeError, "only addition of two finite subgroups is defined"
        if self.abelian_variety() != other.abelian_variety():
            raise TypeError, "finite subgroups must be in the same ambient abelian variety"
        K = Sequence([self.base_field()(0), other.base_field()(0)]).universe()
        return FiniteSubgroup_gens(self.abelian_variety(),
                        self._generators() + other._generators(), base_field=K)

    def __mul__(self, right):
        r = QQ(right)
        G = [r*v for v in self._generators()]
        if r.denominator() != 1:
            L = self._lattice()
            d = 1/r.denominator()
            G += [d * v for v in L.basis()]
        return FiniteSubgroup_gens(self.abelian_variety(), G, base_field = self.base_field())

    def __rmul__(self, left):
        return self * left

    def multiply(self, right):
        return self * right

    def _generators(self):
        """
        Return a list of vectors that define elements of the rational
        homology that generate this finite subgroup.

        Raises a ValueError if no explicit presentation of this finite
        subgroup is known.
        """
        raise ValueError, "no explicit presentation of this finite subgroup is known"

    def abelian_variety(self):
        return self._abvar

    def base_field(self):
        return self._base_field

    def base_ring(self):
        return self.base_field()

    def _repr_(self):
        return "Finite subgroup of %s defined over %s"%(self._abvar, self._base_field)

    def order(self):
        try:
            return self.__order
        except AttributeError:
            if self._abvar.dimension() == 0:
                self.__order = ZZ(1)
                return self.__order
            o = prod(self.invariants())
            self.__order = o
            return o

    def _lattice(self):
        try:
            return self.__lattice
        except AttributeError:
            A = self._abvar
            L = ZZ**(2*A.dimension())
            self.__lattice = L
            return L


    def _full_module(self):
        """
        Return the ZZ-module that corresponds to this finite subgroup.
        This is a ZZ-module that contains the integral homology of the
        ambient variety with finite index.

        EXAMPLES:
            sage: C = J0(11).cuspidal_subgroup()
            sage: C._full_module()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [  1   0]
            [  0 1/5]
        """
        try:
            return self.__full_module
        except AttributeError:
            pass
        G = self._generators()
        V = self._lattice()
        W = V + V.span(G)
        self.__full_module = W
        return W


    def _rescaled_module(self):
        r"""
        Return d * gens as a module, where gens is a list of generators
        of self modulo the $\ZZ^n$.
        """
        try:
            return self.__rescaled_module, self.__denom
        except AttributeError:
            pass
        G = self._generators()
        V = self._lattice()
        if len(G) == 0:
            self.__rescaled_module = V
            self.__denom = ZZ(1)
            return self.__rescaled_module, self.__denom

        d = lcm([v.denominator() for v in G])
        self.__denom = d

        if d == 1:
            B = G
        else:
            B = [d*v for v in G]

        W = V.span(B)
        self.__rescaled_module = W
        return W, d

    def gens(self):
        """
        Return generators for this finite subgroup.

        EXAMPLES:
        We list generators for several cuspidal subgroups:
            sage: J0(11).cuspidal_subgroup().gens()
            [[(0, 1/5)]]
            sage: J0(37).cuspidal_subgroup().gens()
            [[(0, 0, 0, 1/3)]]
            sage: J0(43).cuspidal_subgroup().gens()
            [[(0, 1/7, 0, -1/7, 0, -2/7)]]
            sage: J1(13).cuspidal_subgroup().gens()
            [[(1/19, 0, 0, 9/19)], [(0, 1/19, 1/19, 18/19)]]
        """

        try:
            return self.__gens
        except AttributeError:
            pass

        W, d = self._rescaled_module()

        e = 1/d
        B = [FiniteSubgroupElement(self, e*v) for
                 v in W.basis() if gcd(v.list()) % d != 0]

        #endif
        self.__gens = Sequence(B, immutable=True)
        return self.__gens

    def gen(self, n):
        return self.gens()[n]

    def __call__(self, x):
        if isinstance(x, FiniteSubgroupElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return FiniteSubgroupElement(self, x._element)
            elif x.parent()._abvar == self._abvar:
                return self(x._element)
            else:
                raise TypeError, "no known way to coerce x yet."
        else:
            W, d = self._rescaled_module()
            x = W.ambient_vector_space(x)
            if d * x in W:
                return FiniteSubgroupElement(self, x)
            else:
                raise TypeError, "x does not define an element of self."

    def subgroup(self, gens):
        """
        Return the subgroup of self spanned by the given generators,
        which all must be elements of self.
        """
        if not isinstance(gens, (tuple, list)):
            raise TypeError, "gens must be a list or tuple"
        G = [self(g).element() for g in gens]
        return FiniteSubgroup_gens(self.abelian_variety(), G)

    def invariants(self):
        """
        Return the elementary invariants of this abelian group.
        """
        try:
            return self._invariants
        except AttributeError:
            pass
        W, d = self._rescaled_module()
        B = W.basis_matrix().change_ring(ZZ)
        E = B.elementary_divisors()
        # That the formula below is right is a somewhat tricky diagram change
        # involving the snake lemma and commutative algebra over ZZ.
        # The key conclusion is that the structure of the group is
        # realized as the kernel of a natural map
        #   (Z^n)/(d*Z^n) ---> (Z^n)/(d*Lambda + d*Z^n)
        #
        I = Sequence([d // gcd(e,d) for e in E if e != 0 and e != d])
        I.sort()
        I.set_immutable()
        self._invariants = I
        return I


class FiniteSubgroup_gens(FiniteSubgroup):
    def __init__(self, abvar, gens, base_field=CC):
        FiniteSubgroup.__init__(self, abvar, base_field)
        self._gens = gens

    def _generators(self):
        return self._gens


class FiniteSubgroupElement(ModuleElement):
    def __init__(self, parent, element):
        ModuleElement.__init__(self, parent)
        if element.denominator() == 1:
            element = element.parent().zero_vector()
        self._element = element

    def element(self):
        return self._element

    def _repr_(self):
        return '[%s]'%self._element

    def _add_(self, other):
        return FiniteSubgroupElement(self.parent(), self._element + other._element)

    def _sub_(self, other):
        return FiniteSubgroupElement(self.parent(), self._element - other._element)

    def _neg_(self):
        return FiniteSubgroupElement(self.parent(), -self._element)

    def _rmul_(self, left):
        return FiniteSubgroupElement(self.parent(), left * self._element)

    def _lmul_(self, right):
        return FiniteSubgroupElement(self.parent(), self._element * right)

    def __cmp__(self, right):
        v = self._element - right._element
        if v.denominator() == 1:
            # two elements are equal modulo the lattice
            return 0
        return cmp(self._element, right._element)

    def additive_order(self):
        return self._element.denominator()
