#*****************************************************************************
#       Copyright (C) 2008 Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.sage_object import SageObject
from sage.interfaces.gap import gap, GapElement, gfq_gap_to_sage
from sage.rings.all import Integer
from sage.rings.all import CyclotomicField

class ClassFunction(SageObject):
    """
    A wrapper of GAP's ClassFunction function.

    NOTES:
        It is *not* checked whether the given values describes a character,
        since GAP does not do this.

    EXAMPLES:
        sage: G = CyclicPermutationGroup(4)
        sage: values  = [1, -1, 1, -1]
        sage: chi = ClassFunction(G, values); chi
        Character of Cyclic group of order 4 as a permutation group

    AUTHOR: Franco Saliola (November 2008)
    """
    def __init__(self, G, values):
        r"""
        Returns the character of the group G with values given by the list
        values. The order of the values must correspond to the output of
        G.conjugacy_classes_representatives().

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: chi = ClassFunction(G, values); chi
            Character of Cyclic group of order 4 as a permutation group
        """
        self._group = G
        if isinstance(values, GapElement) and gap.IsClassFunction(values):
            self._gap_classfunction = values
        else:
            self._gap_classfunction = gap.ClassFunction(G, list(values))
        e = self._gap_classfunction.Conductor()
        self._base_ring = CyclotomicField(e)

    def _gap_init_(self):
        r"""
        Returns a string showing how to declare / initialize self in Gap.
        Stored in the \code{self._gap_string} attribute.

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: ClassFunction(G, values)._gap_init_()
            'ClassFunction( CharacterTable( Group( [ (1,2,3,4) ] ) ), [ 1, -1, 1, -1 ] )'
        """
        return str(self._gap_classfunction)

    def _gap_(self, *args):
        r"""
        Coerce self into a GAP element.

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
            sage: values  = [1, -1, 1, -1]
            sage: gap(ClassFunction(G, values))
            ClassFunction( CharacterTable( Group( [ (1,2,3,4) ] ) ), [ 1, -1, 1, -1 ] )
        """
        return self._gap_classfunction

    def __repr__(self):
        r"""
        EXAMPLES:
            sage: G = SymmetricGroup(4)
            sage: values  = [1, -1, 1, 1, -1]
            sage: ClassFunction(G, values)
            Character of Symmetric group of order 4! as a permutation group
        """
        return "Character of %s" % repr(self._group)

    def __iter__(self):
        r"""
        Iterate through the values of self evaluated on the conjugacy
        classes.

        EXAMPLES:
            sage: xi = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: list(xi)
            [1, -1, 1, 1, -1]
        """
        for v in self._gap_classfunction:
            yield self._base_ring(v)

    def domain(self):
        r"""
        Returns the domain of the self.

        EXAMPLES:
            sage: ClassFunction(SymmetricGroup(4), [1,-1,1,1,-1]).domain()
            Symmetric group of order 4! as a permutation group
        """
        return self._group

    def __call__(self, g):
        """
        Evaluate the character on the group element g. Returns an error if
        g is not in G.

        EXAMPLES:
            sage: G = GL(2,7)
            sage: values = list(gap(G).CharacterTable().Irr()[2])
            sage: chi = ClassFunction(G, values)
            sage: z = G([[3,0],[0,3]]); z
            [3 0]
            [0 3]
            sage: chi(z)
            1
            sage: G = GL(2,3)
            sage: chi = G.irreducible_characters()[3]
            sage: g = G.conjugacy_class_representatives()[6]
            sage: chi(g)
            zeta8^3 + zeta8
        """
        rep = self._group._gap_().ConjugacyClass(g).Representative()
        reps = [x.Representative() for x in gap(self._group).ConjugacyClasses()]
        return self._base_ring(self._gap_classfunction[reps.index(rep)+1])

    def __add__(self, other):
        r"""
        Returns the sum of the characters self and other.

        EXAMPLES:
            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: s = chi+chi
            sage: s
            Character of Symmetric group of order 4! as a permutation group
            sage: s.values()
            [6, 2, -2, 0, -2]
        """
        if not isinstance(other, ClassFunction):
            NotImplementedError
        else:
            s = self._gap_classfunction + other._gap_classfunction
            return ClassFunction(self._group, s)

    def __mul__(self, other):
        r"""
        Returns the sum of the characters self and other.

        EXAMPLES:
            sage: G = SymmetricGroup(4)
            sage: chi1 = ClassFunction(G, [3, 1, -1, 0, -1])
            sage: chi2 = ClassFunction(G, [1, -1, 1, 1, -1])
            sage: p = chi1*chi2
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [3, -1, -1, 0, 1]
        """
        if not isinstance(other, ClassFunction):
            NotImplementedError
        else:
            p = self._gap_classfunction * other._gap_classfunction
            return ClassFunction(self._group, p)

    def __pow__(self, other):
        r"""
        Returns the product of self with itself other times.

        EXAMPLES:
            sage: chi = ClassFunction(SymmetricGroup(4), [3, 1, -1, 0, -1])
            sage: p = chi**3
            sage: p
            Character of Symmetric group of order 4! as a permutation group
            sage: p.values()
            [27, 1, -1, 0, -1]
        """
        if isinstance(other, (int,Integer)):
            return ClassFunction(self._group, self._gap_classfunction ** other)
        else:
            NotImplementedError

    def scalar_product(self, other):
        r"""
        Returns the scalar product of self with other.

        EXAMPLES:
            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [[x.scalar_product(y) for x in irr] for y in irr]
            [[1, 0, 0, 0, 0],
             [0, 1, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 0, 0, 1, 0],
             [0, 0, 0, 0, 1]]
        """
        return self._gap_classfunction.ScalarProduct(other)

    def is_irreducible(self):
        r"""
        Returns True if self cannot be written as the sum of two nonzero
        characters of self.

        EXAMPLES:
            sage: S4 = SymmetricGroup(4)
            sage: irr = S4.irreducible_characters()
            sage: [x.is_irreducible() for x in irr]
            [True, True, True, True, True]
        """
        return bool(self._gap_classfunction.IsIrreducible())

    def degree(self):
        r"""
        Returns the degree of the character self.

        EXAMPLES:
            sage: S5 = SymmetricGroup(5)
            sage: irr = S5.irreducible_characters()
            sage: [x.degree() for x in irr]
            [1, 4, 5, 6, 5, 4, 1]
        """
        return Integer(self._gap_classfunction.DegreeOfCharacter())

    def irreducible_constituents(self):
        r"""
        Returns a list of the characters that appear in the decomposition
        of chi.

        EXAMPLES:
            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: irr = chi.irreducible_constituents(); irr
            [Character of Symmetric group of order 5! as a permutation group, Character of Symmetric group of order 5! as a permutation group]
            sage: map(list, irr)
            [[4, -2, 0, 1, 1, 0, -1], [5, -1, 1, -1, -1, 1, 0]]

            sage: G = GL(2,3)
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, -1, -1, -1, -1])
            sage: chi.irreducible_constituents()
            [Character of General Linear Group of degree 2 over Finite Field of size 3]
            sage: chi = ClassFunction(G, [1, 1, 1, 1, 1, 1, 1, 1])
            sage: chi.irreducible_constituents()
            [Character of General Linear Group of degree 2 over Finite Field of size 3]
            sage: chi = ClassFunction(G, [2, 2, 2, 2, 2, 2, 2, 2])
            sage: chi.irreducible_constituents()
            [Character of General Linear Group of degree 2 over Finite Field of size 3]
            sage: chi = ClassFunction(G, [-1, -1, -1, -1, 3, -1, -1, 1])
            sage: ic = chi.irreducible_constituents(); ic
            [Character of General Linear Group of degree 2 over Finite Field of size 3, Character of General Linear Group of degree 2 over Finite Field of size 3]
            sage: map(list, ic)
            [[2, -1, 2, -1, 2, 0, 0, 0], [3, 0, 3, 0, -1, 1, 1, -1]]
        """
        L = self._gap_classfunction.ConstituentsOfCharacter()
        return [ClassFunction(self._group, list(l)) for l in L]

    def decompose(self):
        r"""
        Returns a list of the characters that appear in the decomposition
        of chi.

        EXAMPLES:
            sage: S5 = SymmetricGroup(5)
            sage: chi = ClassFunction(S5, [22, -8, 2, 1, 1, 2, -3])
            sage: chi.decompose()
            [(3, Character of Symmetric group of order 5! as a permutation group), (2, Character of Symmetric group of order 5! as a permutation group)]
        """
        L = []
        for irr in self.irreducible_constituents():
            L.append((self.scalar_product(irr), irr))
        return L

    def norm(self):
        r"""
        Returns the norm of self.

        EXAMPLES:
            sage: A5 = AlternatingGroup(5)
            sage: [x.norm() for x in A5.irreducible_characters()]
            [1, 1, 1, 1, 1]
        """
        return self._gap_classfunction.Norm()

    def values(self):
        r"""
        Return the list of values of self on the conjugacy classes.

        EXAMPLES:
            sage: G = GL(2,3)
            sage: [x.values() for x in G.irreducible_characters()]
            [[1, 1, 1, 1, 1, 1, 1, 1], [1, 1, 1, 1, 1, -1, -1, -1], [2, -1, 2, -1, 2, 0, 0, 0], [2, 1, -2, -1, 0, -zeta8^3 - zeta8, zeta8^3 + zeta8, 0], [2, 1, -2, -1, 0, zeta8^3 + zeta8, -zeta8^3 - zeta8, 0], [3, 0, 3, 0, -1, -1, -1, 1], [3, 0, 3, 0, -1, 1, 1, -1], [4, -1, -4, 1, 0, 0, 0, 0]]
        """
        return list(self)

    def central_character(self):
        r"""
        Returns the central character of self.

        EXAMPLES:
            sage: t = SymmetricGroup(4).trivial_character()
            sage: t.central_character().values()
            [1, 6, 3, 8, 6]
        """
        return ClassFunction(self._group, self._gap_classfunction.CentralCharacter())

    def determinant_character(self):
        r"""
        Returns the determinant character of self.

        EXAMPLES:
            sage: t = ClassFunction(SymmetricGroup(4), [1, -1, 1, 1, -1])
            sage: t.determinant_character().values()
            [1, -1, 1, 1, -1]
        """
        return ClassFunction(self._group, self._gap_classfunction.DeterminantOfCharacter())

    def tensor_product(self, other):
        r"""
        EXAMPLES:
            sage: S3 = SymmetricGroup(3)
            sage: chi1, chi2, chi3 = S3.irreducible_characters()
            sage: chi1.tensor_product(chi3).values()
            [1, -1, 1]
        """
        return ClassFunction(self._group, gap.Tensored([self],[other])[1])

    def restrict(self, H):
        r"""
        EXAMPLES:
            sage: G = SymmetricGroup(5)
            sage: chi = ClassFunction(G, [3, -3, -1, 0, 0, -1, 3]); chi
            Character of Symmetric group of order 5! as a permutation group
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: chi.restrict(H)
            Character of Subgroup of SymmetricGroup(5) generated by [(1,2,3), (1,2), (4,5)]
            sage: chi.restrict(H).values()
            [3, -3, -3, -1, 0, 0]
        """
        rest = self._gap_classfunction.RestrictedClassFunction(H._gap_())
        return ClassFunction(H, rest)

    def induct(self, G):
        r"""
        EXAMPLES:
            sage: G = SymmetricGroup(5)
            sage: H = G.subgroup([(1,2,3), (1,2), (4,5)])
            sage: xi = H.trivial_character(); xi
            Character of Subgroup of SymmetricGroup(5) generated by [(1,2,3), (1,2), (4,5)]
            sage: xi.induct(G)
            Character of Symmetric group of order 5! as a permutation group
            sage: xi.induct(G).values()
            [10, 4, 2, 1, 1, 0, 0]
        """
        rest = self._gap_classfunction.InducedClassFunction(G._gap_())
        return ClassFunction(G, rest)
