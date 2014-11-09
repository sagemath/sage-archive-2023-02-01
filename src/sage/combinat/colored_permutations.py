r"""
Colored Permutations
"""
from sage.categories.groups import Groups
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from sage.combinat.permutation import Permutations
from sage.combinat.cartesian_product import CartesianProduct
from sage.matrix.constructor import diagonal_matrix
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.rings.number_field.number_field import CyclotomicField

class ColoredPermutation(Element):
    """
    A colored permutation.
    """
    def __init__(self, parent, colors, perm):
        """
        Initialize ``self``.
        """
        self._colors = tuple(colors)
        self._perm = perm
        Element.__init__(self, parent=parent)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return repr([list(self._colors), self._perm])

    def _mul_(self, other):
        """
        Multiply ``self`` and ``other``.
        """
        colors = tuple(self._colors[i] + other._colors[val-1]
                       for i,val in enumerate(~self._perm))
        return self.__class__(self.parent(), colors, self._perm * other._perm)

    def __eq__(self, other):
        """
        Check equality.
        """
        if not isinstance(other, ColoredPermutation):
            return False
        return (self.parent() is other.parent()
                and self._colors == other._colors
                and self._perm == other._perm)

    def __ne__(self, other):
        """
        Check inequality.
        """
        return not self.__eq__(other)

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        for i,p in enumerate(self._perm):
            yield (self._colors[i], p)

    def one_line_form(self):
        """
        Return the one line form of ``self``.
        """
        return list(self)

    def colors(self):
        """
        Return the colors of ``self``.
        """
        return list(self._colors)

    def permutation(self):
        """
        Return the permutation of ``self``.
        """
        return self._perm

    def to_matrix(self):
        """
        Return a matrix of ``self``.
        """
        Cp = CyclotomicField(self.parent()._m)
        g = Cp.gen()
        return diagonal_matrix(Cp, [g**i for i in self._colors]) * self._perm.to_matrix()

class ColoredPermutations(Parent, UniqueRepresentation):
    r"""
    The group of `m`-colored permutations on `\{1, 2, \ldots, n\}`.

    Let `S_n` be the symmetric group on `n` letters and `C_m` be the cyclic
    group of order `m`. The `m`-colored permutation group on `n` letters
    is given by `P_n^m = C_m \wr S_n`. This is also the complex reflection
    group `G(m, 1, n)`.

    We define our multiplication by

    .. MATH::

        ((s_1, \ldots s_n), \sigma) \cdot ((t_1, \ldots, t_n), \tau)
        = ((s_1 t_{\sigma^{-1}(1)}, \ldots, s_n t_{\sigma^{-1}(n)}),
           \sigma \tau).

    REFERENCES:

    - :wikipedia:`Generalized_symmetric_group`
    - :wikipedia:`Complex_reflection_group`
    """
    def __init__(self, m, n):
        """
        Initialize ``self``.

        EXAMPLES::
        """
        self._m = m
        self._n = n
        self._C = IntegerModRing(self._m)
        self._P = Permutations(self._n)
        Parent.__init__(self, category=(Groups(), FiniteEnumeratedSets()))

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "{}-colored permutations of size {}".format(self._m, self._n)

    def one(self):
        """
        Return the identity element of ``self``.
        """
        return self.element_class(self, [self._C.zero()]*self._n, self._P.identity())

    @cached_method
    def gens(self):
        """
        Return the generators of ``self``.
        """
        zero = [self._C.zero()]*self._n
        g = []
        for i in range(self._n-1):
            p = range(1, self._n+1)
            p[i] = i+2
            p[i+1] = i+1
            g.append( self.element_class(self, zero, self._P(p)) )
        zero[-1] = self._C.one()
        g.append( self.element_class(self, zero, self._P.identity()) )
        return tuple(g)

    def matrix_group(self):
        """
        Return the matrix group corresponding to ``self``.
        """
        from sage.groups.matrix_gps.finitely_generated import MatrixGroup
        return MatrixGroup([g.to_matrix() for g in self.gens()])

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.
        """
        if isinstance(x, list):
            if isinstance(x[0], tuple):
                c = []
                p = []
                for k in x:
                    if len(k) != 2:
                        raise ValueError("input must be pairs (color, element)")
                    c.append( self._C(k[0]) )
                    p.append(k[1])
                return self.element_class(self, c, self._P(p))

            if len(x) != 2:
                raise ValueError("input must be a pair of a list of colors and a permutation")
            return self.element_class(self, map(self._C, x[0]), self._P(x[1]))

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        C = CartesianProduct(*[self._C]*self._n)
        for p in self._P:
            for c in C:
                yield self.element_class(self, c, p)

    def cardinality(self):
        """
        Return the cardinality of ``self``.
        """
        return self._m**self._n * self._P.cardinality()

    def rank(self):
        """
        Return the rank of ``self``.

        The rank of a complex reflection group is equal to the dimension
        of the complex vector space the group acts on.
        """
        if self._m == 1:
            return self._n - 1
        return self._n

    def degrees(self):
        """
        Return the degrees of ``self``.

        The degrees of a complex reflection group are the degrees of
        the fundamental invariants of the ring of polynomial invariants.

        If `m = 1`, then we are in the special case of the symmetric group
        and the degrees are `(2, 3, \ldots, n, n+1)`. Otherwise the degrees
        are `(m, 2m, \ldots, nm)`.

        EXAMPLES::

            sage: CP = ColoredPermutations(4, 3)
            sage: CP.degrees()
            [4, 8, 12]
            sage: S = ColoredPermutations(1, 3)
            sage: S.degrees()
            [2, 3]

        We now check that the product of the degrees is equal to the
        cardinality of ``self``::

            sage: prod(CP.degrees()) == CP.cardinality()
            True
            sage: prod(S.degrees()) == S.cardinality()
            True
        """
        if self._m == 1: # Special case for the usual symmetric group
            return range(2, self._n+1)
        return [self._m * i for i in range(1, self._n+1)]

    def codegrees(self):
        """
        Return the codegrees of ``self``.

        Let `G` be a complex reflection group. The codegrees
        `d_1^* \leq d_2^* \leq \cdots \leq d_{\ell}^*` of `G` can be
        defined in terms of the fixed point polynomial:

            f_G(q) = \prod_{i=1}^{\ell} (q - d_i^* - 1).

        If `m = 1`, then we are in the special case of the symmetric group
        and the codegrees are `(n-2, n-3, \ldots 1, 0)`. Otherwise the degrees
        are `((n-1)m, (n-2)m, \ldots, m, 0)`.
        """
        if self._m == 1: # Special case for the usual symmetric group
            return list(reversed(range(self._n-1)))
        return [self._m * i for i in reversed(range(self._n))]

    def number_reflection_hyperplanes(self):
        """
        Return the number of reflection hyperplanes of ``self``.

        The number of reflection hyperplanes of a complex reflection
        group is equal to the sume of the codegrees plus the rank.
        """
        return sum(self.codegrees()) + self.rank()

    def fixed_point_polynomial(self, q):
        r"""
        The fixed point polynomial of ``self``.

        The fixed point polynomial `f_G` of a complex reflection group `G` is
        counting the dimesions fixed points subspaces:

        .. MATH::

            f_G(q) = \sum_{w \in W} q^{\dim V^w}.

        Furthermore, let `d_1, d_2, \ldots, d_{\ell}` be the degrees of `G`,
        then the fixed point polynomial is given by

        .. MATH::

            f_G(q) = \prod_{i=1}^{\ell} (q + d_i - 1).
        """
        return prod(q + d - 1 for d in self.degrees())

    def is_well_generated(self):
        """
        Return if ``self`` is a well-generated complex reflection group.

        A complex reflection group `G` is well-generated if it is
        generated by `\ell` reflections. Equivalently, `G` is well-generated
        if `d_i + d_i^* = d_{\ell}` for all `1 \leq i \leq \ell`
        """
        deg = self.degrees()
        dstar = self.codegrees()
        return all(deg[-1] == d + dstar[i] for i,d in enumerate(deg))

    Element = ColoredPermutation

#####################################################################
## Signed permutations

class SignedPermutation(ColoredPermutation):
    """
    A signed permutation.
    """
    def _mul_(self, other):
        """
        Multiply ``self`` and ``other``.
        """
        colors = tuple(self._colors[i] * other._colors[val-1]
                       for i,val in enumerate(~self._perm))
        return self.__class__(self.parent(), colors, self._perm * other._perm)

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        for i,p in enumerate(self._perm):
            yield self._colors[i] * p

    def to_matrix(self):
        """
        Return a matrix of ``self``.
        """
        return identity_matrix(self._colors) * self._perm.to_matrix()

class SignedPermutations(ColoredPermutations):
    r"""
    Group of signed permutations.

    The group of signed permutations is also known as the hyperoctahedral
    group and the 2-colored permutation group. Thus it can be constructed
    as the wreath product `S_2 \wr S_n`.

    REFERENCES:

    - :wikipedia:`Hyperoctahedral_group`
    """
    def __init__(self, n):
        """
        Initialize ``self``.
        """
        ColoredPermutations.__init__(self, 2, n)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Signed permutations of {}".format(self._n)

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        C = CartesianProduct(*[[-1,1]]*self._n)
        for p in self._P:
            for c in C:
                yield self.element_class(self, c, p)

    Element = SignedPermutation

class EvenSignedPermutations(SignedPermutations):
    """
    Group of even signed permutations.
    """
    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Even signed permtuations of {}".format(self._n)

    def __iter__(self):
        """
        Iterate over ``self``.
        """
        for s in SignedPermutations.__iter__(self):
            total = 0
            for pm in s._colors:
                if pm == -1:
                    total += 1

            if total % 2 == 0:
                yield s

