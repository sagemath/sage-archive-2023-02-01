r"""
Finite subgroups of modular abelian varieties

\sage can compute with fairly general finite subgroups of modular
abelian varieties.  Elements of finite order are represented by
equivalence classes of elements in $H_1(A,\QQ)$ modulo $H_1(A,\ZZ)$.
A finite subgroup can be defined by giving generators and via various
other constructions.  Given a finite subgroup, one can compute
generators, as well as the structure as an abstract group.  Arithmetic
on subgroups is also supported, including adding two subgroups
together, checking inclusion, etc.

TODO: Intersection, action of Hecke operators.

AUTHOR:
    -- William Stein (2007-03)

EXAMPLES:
    sage: J = J0(33)
    sage: C = J.cuspidal_subgroup()
    sage: C
    Finite subgroup with invariants [10, 10] over QQbar of Abelian variety J0(33) of dimension 3
    sage: C.order()
    100
    sage: C.gens()
    [[(1/10, 0, 1/10, 1/10, 1/10, 3/10)], [(0, 1/5, 1/10, 0, 1/10, 9/10)], [(0, 0, 1/2, 0, 1/2, 1/2)]]
    sage: C.0 + C.1
    [(1/10, 1/5, 1/5, 1/10, 1/5, 6/5)]
    sage: 10*(C.0 + C.1)
    [(0, 0, 0, 0, 0, 0)]
    sage: G = C.subgroup([C.0 + C.1]); G
    Finite subgroup with invariants [10] over QQbar of Abelian variety J0(33) of dimension 3
    sage: G.gens()
    [[(1/10, 1/5, 1/5, 1/10, 1/5, 6/5)]]
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

TESTS:
    sage: G = J0(11).finite_subgroup([[1/3,0], [0,1/5]]); G
    Finite subgroup with invariants [15] over QQbar of Abelian variety J0(11) of dimension 1
    sage: loads(dumps(G)) == G
    True
    sage: loads(dumps(G.0)) == G.0
    True
"""

###########################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.modules.module      import Module
from sage.structure.element   import ModuleElement
from sage.structure.sequence  import Sequence
from sage.rings.all           import gcd, lcm, QQ, ZZ, QQbar, is_Field, Integer
from sage.misc.misc           import prod

import abvar as abelian_variety

# TODO: obviously this goes somewhere else in the rings module.
def composite_field(K,L):
    """
    Return a field that contains both $K$ and $L$.

    INPUT:
        K -- field
        L -- field
    OUTPUT:
        field

    EXAMPLES:
    """
    if K == L:
        return K
    if K == QQbar:
        return QQbar
    raise NotImplementedError, "need to implement this"

class QQbarTorsionSubgroup(Module):
    def __init__(self, abvar):
        self.__abvar = abvar
        Module.__init__(self, ZZ)

    def _repr_(self):
        return 'Group of all torsion points in QQbar on %s'%self.__abvar

    def field_of_definition(self):
        return self.__abvar.base_field()

    def __call__(self, x):
        v = self.__abvar.vector_space()(x)
        return FiniteSubgroupElement(self, v)

    def abelian_variety(self):
        return self.__abvar


class FiniteSubgroup(Module):
    """
    A finite subgroup of a modular abelian variety.
    """
    def __init__(self, abvar, field_of_definition=QQ):
        """
        Create a finite subgroup of a modular abelian variety.

        INPUT:
            abvar -- a modular abelian variety
            field_of_definition -- a field over which this group is defined.

        EXAMPLES:
        This is an abstract base class, so there are no instances
        of this class itself.
            sage: A = J0(37)
            sage: G = A.n_torsion_subgroup(3); G
            Finite subgroup with invariants [3, 3, 3, 3] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(37)
            sage: type(G)
            <class 'sage.modular.abvar.finite_subgroup.FiniteSubgroup_gens'>
            sage: from sage.modular.abvar.finite_subgroup import FiniteSubgroup
            sage: isinstance(G, FiniteSubgroup)
            True
        """
        if not is_Field(field_of_definition):
            raise TypeError, "field_of_definition must be a field"
        if not abelian_variety.is_ModularAbelianVariety(abvar):
            raise TypeError, "abvar must be a modular abelian variety"
        Module.__init__(self, ZZ)
        self.__abvar = abvar
        self.__field_of_definition = field_of_definition

    def __cmp__(self, other):
        """
        Compare this finite subgroup to other.

        If other is not a modular abelian variety finite subgroup,
        then the types of self and other are compared.  If other is a
        finite subgroup, and the ambient abelian varieties are equal,
        then the subgroups themselves are compared, by comparing their
        full modules.  If the containing abelian varieties are not
        equal and their ambient varieties are different they are
        compared; if they are the same, then a NotImplemnetedError
        is raised (this is temporary).


        EXAMPLES:
        We first compare to subgroups of $J_0(37)$:
            sage: A = J0(37)
            sage: G = A.n_torsion_subgroup(3); G.order()
            81
            sage: H = A.cuspidal_subgroup(); H.order()
            3
            sage: H < G
            True
            sage: H.is_subgroup(G)
            True
            sage: H < 5 #random (meaningless since it depends on memory layout)
            False
            sage: 5 < H #random (meaningless since it depends on memory layout)
            True

        The ambient varieties are compared:
            sage: cmp(A[0].cuspidal_subgroup(), J0(11).cuspidal_subgroup())
            1

        The following case is not implemented yet, but should be:
            sage: cmp(A[0].cuspidal_subgroup(), A[1].cuspidal_subgroup())
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if not isinstance(other, FiniteSubgroup):
            return cmp(type(self), type(other))
        A = self.abelian_variety()
        B = other.abelian_variety()
        c = cmp(A,B)
        if c == 0:
            # Minus sign because order gets reversed in passing to lattices.
            return -cmp(self.full_module(), other.full_module())
        c = cmp(A.ambient_variety(), B.ambient_variety())
        if c:
            return c
        # Ambient varieties are the same, but the varieties that
        # contain these groups aren't -- we should embed each group
        # into ambient variety, but this isn't implemented yet.
        raise NotImplementedError

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
            if abelian_variety.is_ModularAbelianVariety(other):
                if self.abelian_variety().is_subvariety(other):
                    return True
                raise NotImplementedError, "determining general inclusions not completely implemented yet."
            return False
        if self.abelian_variety() != other.abelian_variety():
            return False
        if self.full_module().is_submodule(other.full_module()):
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
        K = Sequence([self.field_of_definition()(0), other.field_of_definition()(0)]).universe()
        return FiniteSubgroup_gens(self.abelian_variety(),
                        self._generators() + other._generators(), field_of_definition=K)

    def exponent(self):
        """
        Return the exponent of this finite abelian group.

        OUTPUT:
            Integer

        EXAMPLES:
            sage: t = J0(33).hecke_operator(7)
            sage: G = t.kernel()[0]; G
            Finite subgroup with invariants [2, 2, 2, 2, 4, 4] over QQ of Abelian variety J0(33) of dimension 3
            sage: G.exponent()
            4
        """
        try:
            return self.__exponent
        except AttributeError:
            e = lcm(self.invariants())
            self.__exponent = e
            return e

    def intersection(self, other):
        """
        Return the intersection of the finite subgroups self and other.

        INPUT:
            other -- a finite group
        OUTPUT:
            a finite group

        EXAMPLES:
            sage: E11a0, E11a1, B = J0(33)
            sage: G = E11a0.n_torsion_subgroup(6); H = E11a0.n_torsion_subgroup(9)
            sage: G.intersection(H)
            Finite subgroup with invariants [3, 3] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: W = E11a1.n_torsion_subgroup(15)
            sage: G.intersection(W)
            Finite subgroup with invariants [3] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: E11a0.intersection(E11a1)[0]
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)

        We intersect subgroups of different abelian varieties.
            sage: E11a0, E11a1, B = J0(33)
            sage: G = E11a0.n_torsion_subgroup(5); H = E11a1.n_torsion_subgroup(5)
            sage: G.intersection(H)
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)
            sage: E11a0.intersection(E11a1)[0]
            Finite subgroup with invariants [5] over QQ of Simple abelian subvariety 11a(1,33) of dimension 1 of J0(33)

        We intersect abelian varieties with subgroups:
            sage: t = J0(33).hecke_operator(7)
            sage: G = t.kernel()[0]; G
            Finite subgroup with invariants [2, 2, 2, 2, 4, 4] over QQ of Abelian variety J0(33) of dimension 3
            sage: A = J0(33).old_subvariety()
            sage: A.intersection(G)
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian variety J0(33) of dimension 3
            sage: A.hecke_operator(7).kernel()[0]
            Finite subgroup with invariants [2, 2, 2, 2] over QQ of Abelian subvariety of dimension 2 of J0(33)
            sage: B = J0(33).new_subvariety()
            sage: B.intersection(G)
            Finite subgroup with invariants [4, 4] over QQ of Abelian variety J0(33) of dimension 3
            sage: B.hecke_operator(7).kernel()[0]
            Finite subgroup with invariants [4, 4] over QQ of Abelian subvariety of dimension 1 of J0(33)
            sage: A.intersection(B)[0]
            Finite subgroup with invariants [3, 3] over QQ of Abelian subvariety of dimension 2 of J0(33)
        """
        A = self.abelian_variety()
        if abelian_variety.is_ModularAbelianVariety(other):
            B = other
            M = B.lattice().scale(Integer(1)/self.exponent())
            K = composite_field(self.field_of_definition(), other.base_field())
        else:
            if not isinstance(other, FiniteSubgroup):
                raise TypeError, "only addition of two finite subgroups is defined"
            B = other.abelian_variety()
            if A.ambient_variety() != B.ambient_variety():
                raise TypeError, "finite subgroups must be in the same ambient product Jacobian"
            M = other.lattice()
            K = composite_field(self.field_of_definition(), other.field_of_definition())

        L = self.lattice()
        if A != B:
            # TODO: This might be way slower than what we could do if
            # we think more carefully.
            C = A + B
            L = L + C.lattice()
            M = M + C.lattice()
        W = L.intersection(M).intersection(A.vector_space())
        # Now write each of a basis for W in terms of the basis for L.
        gens = A.lattice().coordinate_module(W).basis()
        return FiniteSubgroup_gens(self.abelian_variety(), gens, field_of_definition=K)

    def __mul__(self, right):
        """
        Multiply this subgroup by the rational number right.

        If right is an integer the result is a subgroup of self.  If
        right is a rational number $n/m$, then this group is first
        divided by $m$ then multiplied by $n$.

        INPUT:
            right -- a rational number

        OUTPUT:
            a subgroup

        EXAMPLES:
            sage: J = J0(37)
            sage: H = J.cuspidal_subgroup(); H.order()
            3
            sage: G = H * 3; G.order()
            1
            sage: G = H * (1/2); G.order()
            48
            sage: J.n_torsion_subgroup(2) + H == G
            True
            sage: G = H*(3/2); G.order()
            16
        """
        r = QQ(right)
        G = [r*v for v in self._generators()]
        if r.denominator() != 1:
            L = self._ambient_lattice()
            d = 1/r.denominator()
            G += [d * v for v in L.basis()]
        return FiniteSubgroup_gens(self.abelian_variety(), G, field_of_definition = self.field_of_definition())

    def __rmul__(self, left):
        """
        Multiply this finite subgroup on the left by an integer.

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.cuspidal_subgroup(); factor(G.order())
            2^8 * 3^2
            sage: H = G.__rmul__(2)
            sage: H.order().factor()
            2^4 * 3^2

        Automatic calling is currently broken because of trac \#2283.
            sage: 2*G
            Finite subgroup with invariants [6, 24] over QQbar of Abelian variety J0(42) of dimension 5
        """
        return self * left

    def multiply(self, right):
        """
        Multiply this finite group by the rational number right.

        INPUT:
            right -- rational number

        OUTPUT:
            a new finite group

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.cuspidal_subgroup(); factor(G.order())
            2^8 * 3^2
            sage: G.multiply(3).order()
            256
            sage: G.multiply(0).order()
            1
            sage: G.multiply(1/5).order()
            22500000000
        """

        return self * right

    def _generators(self):
        """
        Return a tuple of vectors that define elements of the rational
        homology that generate this finite subgroup.

        Raises a ValueError if no explicit presentation of this finite
        subgroup is known, i.e., if this function hasn't been
        overridden in a derived class.

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.torsion_subgroup(); G
            Torsion subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(42)
            sage: G._generators()
            Traceback (most recent call last):
            ...
            ValueError: no explicit presentation of this finite subgroup is known (unable to compute explicitly)
        """
        raise ValueError, "no explicit presentation of this finite subgroup is known"

    def abelian_variety(self):
        """
        Return the abelian variety that this is a finite subgroup of.

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.torsion_subgroup(); G
            Torsion subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(42)
            sage: G.abelian_variety()
            Jacobian of the modular curve associated to the congruence subgroup Gamma0(42)
        """
        return self.__abvar

    def field_of_definition(self):
        """
        Return the field over which this finite modular abelian
        variety subgroup is defined.  This is a field over which
        this subgroup is defined.

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.torsion_subgroup(); G
            Torsion subgroup of Jacobian of the modular curve associated to the congruence subgroup Gamma0(42)
            sage: G.field_of_definition()
            Rational Field
        """
        return self.__field_of_definition

    def _repr_(self):
        """
        Return string representation of this finite subgroup.

        EXAMPLES:
            sage: J = J0(42)
            sage: G = J.n_torsion_subgroup(3); G._repr_()
            'Finite subgroup with invariants [3, 3, 3, 3, 3, 3, 3, 3, 3, 3] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(42)'
        """
        K = self.__field_of_definition
        if K == QQbar:
            field = "QQbar"
        elif K == QQ:
            field = "QQ"
        else:
            field = str(K)
        return "Finite subgroup %sover %s of %s"%(self._invariants_repr(), field, self.__abvar)

    def _invariants_repr(self):
        """
        The string representation of the 'invariants' part of this group.

        We make this a separate function so it is possible to create
        finite subgroups that don't print their invariants, since
        printing them could be expensive.

        EXAMPLES:
            sage: J0(42).cuspidal_subgroup()._invariants_repr()
            'with invariants [2, 2, 12, 48] '
        """
        return 'with invariants %s '%(self.invariants(), )

    def order(self):
        """
        Return the order (number of elements) of this finite subgroup.

        EXAMPLES:
            sage: J = J0(42)
            sage: C = J.cuspidal_subgroup()
            sage: C.order()
            2304
        """
        try:
            return self.__order
        except AttributeError:
            if self.__abvar.dimension() == 0:
                self.__order = ZZ(1)
                return self.__order
            o = prod(self.invariants())
            self.__order = o
            return o

    def _ambient_lattice(self):
        r"""
        Return the lattice corresponding to this finite subgroup --
        this is just a free module of rank equal to the rank of the
        ambient abelian variety, i.e., it's $H_1(A, \Z)$.  This is
        used internally by other algorithms.

        EXAMPLES:
            sage: J = J0(23)
            sage: C = J.cuspidal_subgroup()
            sage: C._ambient_lattice().basis()
            [
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, 1)
            ]

        We emphasize that lattice has essentially nothing to do with
        this finite subgroup -- it's really a function of the ambient
        variety.  You really want \code{lattice} for something
        about $C$ itself.
            sage: C.full_module().basis()
            [
            (1/11, 10/11, 0, 8/11),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, 1)
            ]
        """
        try:
            return self.__ambient_lattice
        except AttributeError:
            A = self.__abvar
            L = ZZ**(2*A.dimension())
            self.__ambient_lattice = L
            return L


    def full_module(self):
        r"""
        Return the $\ZZ$-module that corresponds to this finite subgroup.
        This is a $\ZZ$-module that contains the integral homology of the
        ambient variety with finite index.  It has degree equal to twice
        the dimension of the ambient abelian variety.

        EXAMPLES:
            sage: C = J0(11).cuspidal_subgroup()
            sage: C.full_module()
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
        V = self._ambient_lattice()
        W = V + V.span(G)
        self.__full_module = W
        return W

    def lattice(self):
        """
        Return the lattice in the homology the modular Jacobian
        product corresponding to this subgroup.  The elements of the
        subgroup are represented by vecotrs in the ambient vector
        space (the rational homology), and this returns the lattice
        they span.

        EXAMPLES:
            sage: J = J0(33); C = J[0].cuspidal_subgroup(); C
            Cuspidal subgroup with invariants [2, 2] over QQ of Abelian variety factor of dimension 1 of J0(33)
            sage: C.lattice()
            Free module of degree 6 and rank 2 over Integer Ring
            Echelon basis matrix:
            [ 1/2    0    0 -1/2    0    0]
            [   0    0  1/2    0  1/2 -1/2]
        """
        try:
            return self.__lattice
        except AttributeError:
            L = self.full_module()
            B = L.basis_matrix() * self.abelian_variety().lattice().basis_matrix()
            self.__lattice = B.row_module(ZZ)
            return self.__lattice



    def rescaled_module(self):
        r"""
        Return d * M and the integer d as a module, where M is the
        ZZ-module generated by gens modulo $\ZZ^n$.  This is the
        analogue of \code{full_module}, except it returns a
        not-necessarily full rank integral full_module and the denominator
        instead of a (possibly) non-integral full rank full_module.

        EXAMPLES:
            sage: J = J0(23)
            sage: C = J.cuspidal_subgroup()
            sage: C.full_module()
            Free module of degree 4 and rank 4 over Integer Ring
            Echelon basis matrix:
            [ 1/11 10/11     0  8/11]
            [    0     1     0     0]
            [    0     0     1     0]
            [    0     0     0     1]
            sage: C.rescaled_module()
            (Free module of degree 4 and rank 1 over Integer Ring
            Echelon basis matrix:
            [ 1 -1  0 -3], 11)
        """
        try:
            return self.__rescaled_module, self.__denom
        except AttributeError:
            pass
        G = self._generators()
        V = self._ambient_lattice()
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

        W, d = self.rescaled_module()

        e = 1/d
        B = [FiniteSubgroupElement(self, e*v, check=False) for
                 v in W.basis() if gcd(v.list()) % d != 0]

        #endif
        self.__gens = Sequence(B, immutable=True)
        return self.__gens

    def gen(self, n):
        r"""
        Return $n$th generator of self.

        EXAMPLES:
            sage: J = J0(23)
            sage: C = J.n_torsion_subgroup(3)
            sage: C.gens()
            [[(1/3, 0, 0, 0)], [(0, 1/3, 0, 0)], [(0, 0, 1/3, 0)], [(0, 0, 0, 1/3)]]
            sage: C.gen(0)
            [(1/3, 0, 0, 0)]
            sage: C.gen(3)
            [(0, 0, 0, 1/3)]
            sage: C.gen(4)
            Traceback (most recent call last):
            ...
            IndexError: list index out of range

        Negative indices wrap around:
            sage: C.gen(-1)
            [(0, 0, 0, 1/3)]
        """
        return self.gens()[n]

    def __call__(self, x):
        r"""
        Coerce $x$ into this finite subgroup.

        This works when the abelian varieties that contains x and
        self are the same, or if $x$ is coercible into the rational
        homology (viewed as an abstract $\QQ$-vector space).

        EXAMPLES:
        We first construct the $11$-torsion subgroup of $J_0(23)$:
            sage: J = J0(23)
            sage: G = J.n_torsion_subgroup(11)
            sage: G.invariants()
            [11, 11, 11, 11]

        We also construct the cuspidal subgroup.
            sage: C = J.cuspidal_subgroup()
            sage: C.invariants()
            [11]

        Coercing something into its parent returns it:
            sage: G(G.0) is G.0
            True

        We coerce an element from the cuspidal subgroup into the
        $11$-torsion subgroup:
            sage: z = G(C.0); z
            [(1/11, -1/11, 0, -3/11)]
            sage: z.parent() == G
            True

        We coerce a list, which defines an element of the underlying
        \code{full_module} into $G$, and verify an equality:
            sage: x = G([1/11, 1/11, 0, -1/11])
            sage: x == G([1/11, 1/11, 0, 10/11])
            True

        Finally we coerce in two elements that shouldn't work, since
        they do not define elements of $G$:
            sage: G(J.n_torsion_subgroup(3).0)
            Traceback (most recent call last):
            ...
            TypeError: x does not define an element of self
            sage: G(J0(11).n_torsion_subgroup(11).0)
            Traceback (most recent call last):
            ...
            TypeError: x does not define an element of self
        """
        if isinstance(x, FiniteSubgroupElement):
            if x.parent() is self:
                return x
            elif x.parent() == self:
                return FiniteSubgroupElement(self, x.element(), check=False)
            elif x.parent().abelian_variety() == self.abelian_variety():
                return self(x.element())
            else:
                raise TypeError, "x does not define an element of self"
        else:
            W, d = self.rescaled_module()
            x = W.ambient_vector_space()(x)
            if d * x in W:
                return FiniteSubgroupElement(self, x, check=False)
            else:
                raise TypeError, "x does not define an element of self"

    def __contains__(self, x):
        try:
            self(x)
        except TypeError:
            return False
        return True

    def subgroup(self, gens):
        """
        Return the subgroup of self spanned by the given generators,
        which all must be elements of self.

        EXAMPLES:
            sage: J = J0(23)
            sage: G = J.n_torsion_subgroup(11); G
            Finite subgroup with invariants [11, 11, 11, 11] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)

        We create the subgroup of the 11-torsion subgroup of $J_0(23)$ generated
        by the first $11$-torsion point:
            sage: H = G.subgroup([G.0]); H
            Finite subgroup with invariants [11] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(23)
            sage: H.invariants()
            [11]

        We can also create a subgroup from a list of objects that
        coerce into the ambient rational homology.
            sage: H == G.subgroup([[1/11,0,0,0]])
            True
        """
        if not isinstance(gens, (tuple, list)):
            raise TypeError, "gens must be a list or tuple"
        G = [self(g).element() for g in gens]
        return FiniteSubgroup_gens(self.abelian_variety(), G)

    def invariants(self):
        """
        Return elementary invariants of this abelian group, by which
        we mean a nondecreasing (immutable) sequence of integers
        $n_i$, $1 \leq i \leq k$, with $n_i$ dividing $n_{i+1}$, and
        such that this group is abstractly isomorphic to
        $\ZZ/n_1\ZZ \times\cdots\times \ZZ/n_k\ZZ.$

        EXAMPLES:
            sage: J = J0(38)
            sage: C = J.cuspidal_subgroup(); C
            Cuspidal subgroup with invariants [3, 45] over QQ of Jacobian of the modular curve associated to the congruence subgroup Gamma0(38)
            sage: v = C.invariants(); v
            [3, 45]
            sage: v[0] = 5
            Traceback (most recent call last):
            ...
            ValueError: object is immutable; please change a copy instead.
            sage: type(v[0])
            <type 'sage.rings.integer.Integer'>
        """
        try:
            #print "getting", self.__invariants
            return self.__invariants
        except AttributeError:
            pass
        W, d = self.rescaled_module()
        B = W.basis_matrix().change_ring(ZZ)
        E = B.elementary_divisors()
        # That the formula below is right is a somewhat tricky diagram
        # chase involving the snake lemma and commutative algebra over
        # ZZ.  The key conclusion is that the structure of the group
        # is realized as the kernel of a natural map
        #   (Z^n)/(d*Z^n) ---> (Z^n)/(d*Lambda + d*Z^n)
        #
        v = [d // gcd(e,d) for e in E if e != 0 and e != d]
        v = [x for x in v if v != 1]
        I = Sequence(v)
        I.sort()
        I.set_immutable()
        print "setting", I
        self.__invariants = I
        return I


class FiniteSubgroup_gens(FiniteSubgroup):
    """
    A finite subgroup of a modular abelian variety that is generated
    by given generators.
    """
    def __init__(self, abvar, gens, field_of_definition=QQbar, check=True):
        """
        Create a finite subgroup with given generators.

        INPUT:
            abvar -- a modular abelian variety
            gens -- a list or tuple of generators of a subgroup
            check -- bool (default: True) whether or not to check that each
                     generator is a FiniteSubgroupElement with
                     the abvar as abelian variety.

        EXAMPLES:
            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
        """
        if check:
            HQ = abvar._rational_homology_space()
            v = []
            for g in gens:
                if not isinstance(g, FiniteSubgroupElement):
                    g = HQ(g)
                else:
                    # it is a FiniteSubgroupElement
                    if g.parent().abelian_variety() != abvar:
                        # TODO
                        raise NotImplementedError, "coercion of elements into abelian varieties not implemented in general"
                    # Now g is a FiniteSubgroupElement with correct ambient variety
                    # but not in this subgroup.
                    g = g.element()
                # done
                v.append(g)
        else:
            v = gens

        FiniteSubgroup.__init__(self, abvar, field_of_definition)
        self.__v = tuple(v)

    def _generators(self):
        r"""
        Return a tuple of elements of the vector space underlying the
        rational homology that together generator this abelian
        variety. They need not be independent.  This is mainly for
        internal use; use the \code{gens()} method for independent
        \code{FiniteSubgroupElement}s.

        EXAMPLES:
            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
            sage: G._generators()
            ((1/3, 0), (0, 1/5))
        """
        return self.__v


class FiniteSubgroupElement(ModuleElement):
    """
    An element of a finite subgroup of a modular abelian variety.
    """
    def __init__(self, parent, element, check=True):
        """
        Create a finite subgroup element.

        INPUT:
            parent  -- a finite subgroup of a modular abelian variety
            element -- a QQ vector space element that represents this
                       element in rational homology modulo integral
                       homology.
            check   -- bool (default: True) whether to check that element
                       is in the appropriate vector space

        EXAMPLES:
        The following calls the FiniteSubgroupElement constructor implicitly:
            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
            sage: type(G.0)
            <class 'sage.modular.abvar.finite_subgroup.FiniteSubgroupElement'>
        """
        ModuleElement.__init__(self, parent)
        # TODO: check
        if element.denominator() == 1:
            element = element.parent().zero_vector()
        self.__element = element

    def ambient_element(self):
        try:
            return self.__ambient_element
        except AttributeError:
            B = self.parent().abelian_variety().lattice().basis_matrix()
            self.__ambient_element = self.__element * B
            return self.__ambient_element

    def element(self):
        """
        Return an underlying QQ-vector space element that defines this
        element of a modular abelian variety.

        EXAMPLES:
        We create some elements of $J_0(11)$:
            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
            sage: G.0.element()
            (1/3, 0)

        The underlying element is a vector over the rational numbers:
            sage: v = (G.0-G.1).element(); v
            (1/3, -1/5)
            sage: type(v)
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        return self.__element

    def _repr_(self):
        r"""
        Return string representation of this finite subgroup element.
        Since they are represented as equivalences classes of rational
        homology modulo integral homology, we represent an element
        corresponding to $v$ in the rational homology by \code{[v]}.

        EXAMPLES:
            sage: J = J0(11)
            sage: G = J.finite_subgroup([[1/3,0], [0,1/5]]); G
            Finite subgroup with invariants [15] over QQbar of Jacobian of the modular curve associated to the congruence subgroup Gamma0(11)
            sage: G.0._repr_()
            '[(1/3, 0)]'
        """
        return '[%s]'%self.__element

    def _add_(self, other):
        """
        Add two finite subgroup elements with the same parent.  This
        is called implicitly by +.

        INPUT:
            other -- a FiniteSubgroupElement with the same parent as self

        OUTPUT:
            a FiniteSubgroupElement

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._add_(G.1)
            [(1/3, 1/5)]
            sage: G.0 + G.1
            [(1/3, 1/5)]
        """
        return FiniteSubgroupElement(self.parent(), self.__element + other.__element, check=False)

    def _sub_(self, other):
        """
        Subtract two finite subgroup elements with the same parent.  This
        is called implicitly by +.

        INPUT:
            other -- a FiniteSubgroupElement with the same parent as self

        OUTPUT:
            a FiniteSubgroupElement

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._sub_(G.1)
            [(1/3, -1/5)]
            sage: G.0 - G.1
            [(1/3, -1/5)]
        """
        return FiniteSubgroupElement(self.parent(), self.__element - other.__element, check=False)

    def _neg_(self):
        """
        Negate a finite subgroup element.

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._neg_()
            [(-1/3, 0)]
        """
        return FiniteSubgroupElement(self.parent(), -self.__element, check=False)

    def _rmul_(self, left):
        """
        Left multiply a finite subgroup element by an integer.

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._rmul_(2)
            [(2/3, 0)]
            sage: 2*G.0
            [(2/3, 0)]
        """
        return FiniteSubgroupElement(self.parent(), ZZ(left) * self.__element, check=False)

    def _lmul_(self, right):
        """
        Right multiply a finite subgroup element by an integer.

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0._lmul_(2)
            [(2/3, 0)]
            sage: G.0 * 2
            [(2/3, 0)]
        """
        return FiniteSubgroupElement(self.parent(), self.__element * ZZ(right), check=False)

    def __cmp__(self, right):
        """
        Compare self and right.

        INPUT:
            self, right -- elements of the same finite abelian variety
            subgroup.

        OUTPUT:
            -1, 0, or 1

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: cmp(G.0, G.1)
            1
            sage: cmp(G.0, G.0)
            0
            sage: 3*G.0 == 0
            True
            sage: 3*G.0 == 5*G.1
            True

        We make sure things that shouldn't be equal aren't:
            sage: H = J0(14).finite_subgroup([[1/3,0]])
            sage: G.0 == H.0
            False
            sage: cmp(G.0, H.0)
            -1
            sage: G.0
            [(1/3, 0)]
            sage: H.0
            [(1/3, 0)]
        """
        if self.parent() != right.parent():
            return cmp(self.parent(), right.parent())
        v = self.__element - right.__element
        if v.denominator() == 1:
            # two elements are equal if they are equal modulo the lattice
            return 0
        return cmp(self.__element, right.__element)

    def additive_order(self):
        """
        Return the additive order of this element.

        EXAMPLES:
            sage: J = J0(11); G = J.finite_subgroup([[1/3,0], [0,1/5]])
            sage: G.0.additive_order()
            3
            sage: G.1.additive_order()
            5
            sage: (G.0 + G.1).additive_order()
            15
            sage: (3*G.0).additive_order()
            1
        """
        return self.__element.denominator()
