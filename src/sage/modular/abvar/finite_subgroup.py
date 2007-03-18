"""
Finite subgroups of modular abelian varieties

EXAMPLES:

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


from sage.modules.module      import Module
from sage.structure.element   import ModuleElement
from sage.structure.sequence  import Sequence
from sage.rings.all           import gcd, lcm, ZZ

class FiniteSubgroup(Module):
    def __init__(self, abvar):
        self._abvar = abvar

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

    def _repr_(self):
        return "Finite subgroup of %s"%self._abvar

    def order(self):
        try:
            return self.__order
        except AttributeError:
            if self._abvar.dimension() == 0:
                self.__order = ZZ(1)
                return self.__order
            W, d = self._rescaled_module()
            # compute the order
            s = ZZ(W.index_in_saturation())
            o = (d**W.rank())//s
            self.__order = o
            return o

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
        A = self._abvar
        V = ZZ**(2*A.dimension())
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
        gens = [self(g) for g in gens]




class FiniteSubgroupElement(ModuleElement):
    def __init__(self, parent, element):
        ModuleElement.__init__(self, parent)
        self._element = element

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
