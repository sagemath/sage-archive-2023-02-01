"""
Matrix Groups

AUTHORS:

- William Stein: initial version

- David Joyner (2006-03-15): degree, base_ring, _contains_, list,
  random, order methods; examples

- William Stein (2006-12): rewrite

- David Joyner (2007-12): Added invariant_generators (with Martin
  Albrecht and Simon King)

- David Joyner (2008-08): Added module_composition_factors (interface
  to GAP's MeatAxe implementation) and as_permutation_group (returns
  isomorphic PermutationGroup).

This class is designed for computing with matrix groups defined by
a (relatively small) finite set of generating matrices.

EXAMPLES::

    sage: F = GF(3)
    sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
    sage: G = MatrixGroup(gens)
    sage: G.conjugacy_class_representatives()
    [
    [1 0]
    [0 1],
    [0 1]
    [2 1],
    [0 1]
    [2 2],
    [0 2]
    [1 1],
    [0 2]
    [1 2],
    [0 1]
    [2 0],
    [2 0]
    [0 2]
    ]

Loading and saving work::

    sage: G = GL(2,5); G
    General Linear Group of degree 2 over Finite Field of size 5
    sage: loads(dumps(G)) == G
    True
    sage: g = G.1; g
    [4 1]
    [4 0]
    sage: loads(dumps(g)) == g
    True
"""

##############################################################################
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.misc.randstate import current_randstate
from matrix_group_element import MatrixGroupElement
from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring, infinity
from sage.misc.functional import is_field
from sage.rings.finite_field import is_FiniteField
from sage.interfaces.gap import gap, GapElement
from sage.matrix.all import MatrixSpace, is_MatrixSpace, is_Matrix
import sage.rings.integer as integer
from sage.misc.latex import latex
from sage.structure.sequence import Sequence
from sage.structure.sage_object import SageObject
from sage.groups.class_function import ClassFunction
from sage.misc.randstate import current_randstate

#################################################################

class MatrixGroup_generic(Group):
    pass

def is_MatrixGroup(x):
    """
    EXAMPLES::

        sage: from sage.groups.matrix_gps.matrix_group import is_MatrixGroup
        sage: is_MatrixGroup(MatrixSpace(QQ,3))
        False
        sage: is_MatrixGroup(Mat(QQ,3))
        False
        sage: is_MatrixGroup(GL(2,ZZ))
        True
        sage: is_MatrixGroup(MatrixGroup([matrix(2,[1,1,0,1])]))
        True
    """
    return isinstance(x, MatrixGroup_generic)

def MatrixGroup(gens):
    r"""
    Return the matrix group with given generators.

    INPUT:


    -  ``gens`` - list of matrices in a matrix space or
       matrix group


    EXAMPLES::

        sage: F = GF(5)
        sage: gens = [matrix(F,2,[1,2, -1, 1]), matrix(F,2, [1,1, 0,1])]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 5 with 2 generators:
        [[[1, 2], [4, 1]], [[1, 1], [0, 1]]]

    In the second example, the generators are a matrix over
    `\mathbb{Z}`, a matrix over a finite field, and the integer
    `2`. Sage determines that they both canonically map to
    matrices over the finite field, so creates that matrix group
    there.

    ::

        sage: gens = [matrix(2,[1,2, -1, 1]), matrix(GF(7), 2, [1,1, 0,1]), 2]
        sage: G = MatrixGroup(gens); G
        Matrix group over Finite Field of size 7 with 3 generators:
        [[[1, 2], [6, 1]], [[1, 1], [0, 1]], [[2, 0], [0, 2]]]

    Each generator must be invertible::

        sage: G = MatrixGroup([matrix(ZZ,2,[1,2,3,4])])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix but one is not:
        [1 2]
        [3 4]

    Some groups aren't supported::

        sage: SL(2, CC).gens()
        Traceback (most recent call last):
        ...
        NotImplementedError: Matrix group over Complex Field with 53 bits of precision not implemented.
        sage: G = SL(0, QQ)
        Traceback (most recent call last):
        ...
        ValueError: The degree must be at least 1
    """
    if len(gens) == 0:
        raise ValueError, "gens must have positive length"
    try:
        R = gens[0].base_ring()
    except AttributeError:
        raise TypeError, "gens must be a list of matrices"
    for i in range(len(gens)):
        if not is_Matrix(gens[i]):
            try:
                gens[i] = gens[i].matrix()
            except AttributeError:
                pass
    if is_FiniteField(R):
        return MatrixGroup_gens_finite_field(gens)
    else:
        return MatrixGroup_gens(gens)

class MatrixGroup_gap(MatrixGroup_generic):
    def __init__(self, n, R, var='a'):
        """
        INPUT:


        -  ``n`` - the degree

        -  ``R`` - the base ring

        -  ``var`` - variable used to define field of
           definition of actual matrices in this group.
        """
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R


        self._var = var
        self.__n = integer.Integer(n)
        if self.__n <= 0:
            raise ValueError, "The degree must be at least 1"
        self.__R = R

    def _gap_(self, G=None):
        try:
            return SageObject._gap_(self, G)
        except TypeError:
            raise NotImplementedError, "Matrix group over %s not implemented."%self.__R

    def __cmp__(self, H):
        if not isinstance(H, MatrixGroup_gap):
            return cmp(type(self), type(H))
        return cmp(gap(self), gap(H))

    def __call__(self, x):
        """
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
            sage: G(1)
            [1 0]
            [0 1]
        """
        if isinstance(x, MatrixGroupElement) and x.parent() is self:
            return x
        M = self.matrix_space()(x)
        g = MatrixGroupElement(M, self)
        if not gap(g) in gap(self):
            raise TypeError, "no way to coerce element to self."
        return g

    def _Hom_(self, G, cat=None):
        if not (cat is None or (cat is G.category() and cat is self.category())):
            raise NotImplementedError
        if not is_MatrixGroup(G):
            raise TypeError, "G (=%s) must be a matrix group."%G
        import homset
        return homset.MatrixGroupHomset(self, G)

    def hom(self, x):
        v = Sequence(x)
        U = v.universe()
        if not is_MatrixGroup(U):
            raise TypeError, "u (=%s) must have universe a matrix group."%U
        return self.Hom(U)(x)

    def matrix_space(self):
        """
        Return the matrix space corresponding to this matrix group.

        This is a matrix space over the field of definition of this matrix
        group.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G.matrix_space()
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 5
        """
        try:
            return self.__matrix_space
        except AttributeError:
            pass
        self.__matrix_space = MatrixSpace(self.field_of_definition(), self.__n)
        return self.__matrix_space

    def degree(self):
        """
        Return the degree of this matrix group.

        EXAMPLES::

            sage: SU(5,5).degree()
            5
        """
        return self.__n

    def field_of_definition(self, var='a'):
        """
        Return a field that contains all the matrices in this matrix
        group.

        EXAMPLES::

            sage: G = SU(3,GF(5))
            sage: G.base_ring()
            Finite Field of size 5
            sage: G.field_of_definition()
            Finite Field in a of size 5^2
            sage: G = GO(4,GF(7),1)
            sage: G.field_of_definition()
            Finite Field of size 7
            sage: G.base_ring()
            Finite Field of size 7
        """
        return self.__R

    def base_ring(self):
        """
        Return the base ring of this matrix group.

        EXAMPLES::

            sage: GL(2,GF(3)).base_ring()
            Finite Field of size 3
            sage: G = SU(3,GF(5))
            sage: G.base_ring()
            Finite Field of size 5
            sage: G.field_of_definition()
            Finite Field in a of size 5^2
        """
        return self.__R

    base_field = base_ring

    def is_finite(self):
        """
        Return True if this matrix group is finite.

        EXAMPLES::

            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
            sage: SL(2,ZZ).is_finite()
            False
        """
        if self.base_ring().is_finite():
            return True
        return self._gap_().IsFinite().bool()

    def order(self):
        """
        EXAMPLES::

            sage: G = Sp(4,GF(3))
            sage: G.order()
            51840
            sage: G = SL(4,GF(3))
            sage: G.order()
            12130560
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            480
            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.order()
            +Infinity
        """
        g = self._gap_()
        if g.IsFinite().bool():
            return integer.Integer(gap(self).Size())
        return infinity

    def gens(self):
        """
        Return generators for this matrix group.

        EXAMPLES::

            sage: G = GO(3,GF(5))
            sage: G.gens()
            [
            [2 0 0]
            [0 3 0]
            [0 0 1],
            [0 1 0]
            [1 4 4]
            [0 2 1]
            ]
        """
        try:
            return self.__gens
        except AttributeError:
            pass
        F = self.field_of_definition()
        gap_gens = list(gap(self).GeneratorsOfGroup())
        gens = Sequence([MatrixGroupElement(g._matrix_(F), self, check=False) for g in gap_gens],
                        cr=True, universe=self, check=False)
        self.__gens = gens
        return gens

    def ngens(self):
        """
        Return the number of generators of this linear group.

        EXAMPLES::

            sage: G = GO(3,GF(5))
            sage: G.ngens()
            2
        """
        return len(self.gens())


    def gen(self, n):
        """
        Return the n-th generator.

        EXAMPLES::

            sage: G = GU(4,GF(5), var='beta')
            sage: G.gen(0)
            [  beta      0      0      0]
            [     0      1      0      0]
            [     0      0      1      0]
            [     0      0      0 3*beta]
        """
        return self.gens()[n]

    def as_matrix_group(self):
        """
        Return this group, but as a general matrix group, i.e., throw away
        the extra structure of general unitary group.

        EXAMPLES::

            sage: G = SU(4,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field in a of size 5^2 with 2 generators:
            [[[a, 0, 0, 0], [0, 2*a + 3, 0, 0], [0, 0, 4*a + 1, 0], [0, 0, 0, 3*a]], [[1, 0, 4*a + 3, 0], [1, 0, 0, 0], [0, 2*a + 4, 0, 1], [0, 3*a + 1, 0, 0]]]

        ::

            sage: G = GO(3,GF(5))
            sage: G.as_matrix_group()
            Matrix group over Finite Field of size 5 with 2 generators:
            [[[2, 0, 0], [0, 3, 0], [0, 0, 1]], [[0, 1, 0], [1, 4, 4], [0, 2, 1]]]
        """
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        return MatrixGroup([g.matrix() for g in self.gens()])

    def list(self):
        """
        Return list of all elements of this group.

        Always returns a new list, so it is safe to change the returned
        list.

        EXAMPLES::

            sage: F = GF(3)
            sage: gens = [matrix(F,2, [1,0, -1,1]), matrix(F, 2, [1,1,0,1])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            24
            sage: v = G.list()
            sage: len(v)
            24
            sage: v[:2]
            [[0 1]
            [2 0], [0 1]
            [2 1]]
            sage: G.list()[0] in G
            True

        ::

            sage: GL(2,ZZ).list()
            Traceback (most recent call last):
            ...
            ValueError: group must be finite
        """
        # We check the cache for the result
        try:
            return list(self.__list)
        except AttributeError:
            pass
        if not self.is_finite():
            raise ValueError, "group must be finite"

        # Get basic properties of the field over which we are working
        F = self.field_of_definition()
        n = F.degree()
        p = F.characteristic()
        a = F.prime_subfield().multiplicative_generator()
        b = F.multiplicative_generator()

        # Get string representation of the list of elements of self.
        # Since the output is usually big, we use a file, which can
        # easily give us a hundred-times speedup for at all large output.
        s = self._gap_().Elements().str(use_file=True)
        s = ''.join(s.split())

        # Replace the two types of gap-style 'power of generator' notation
        s = s.replace('Z(%s^%s)'%(p,n),'b')
        s = s.replace('Z(%s)'%p,'a')
        s = s.replace('^','**')
        # Then eval the string with a and b set to the corresponding
        # multiplicative generators.
        v = eval(s, {'a':a, 'b':b})

        # Finally, create the matrix space in which all these matrices live,
        # and make each element as a MatrixGroupElement.
        MS = self.matrix_space()
        v = [MatrixGroupElement(MS(x), self, check=False) for x in v]
        self.__list = v
        return list(v)

    def irreducible_characters(self):
        """
        Returns the list of irreducible characters of the group.

        EXAMPLES::

            sage: G = GL(2,2)
            sage: G.irreducible_characters()
            [Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2,
             Character of General Linear Group of degree 2 over Finite Field of size 2]

        """
        current_randstate().set_seed_gap()
        Irr = self._gap_().Irr()
        L = []
        for irr in Irr:
            L.append(ClassFunction(self,irr))
        return L

class MatrixGroup_gap_finite_field(MatrixGroup_gap):
    """
    Python class for matrix groups over a finite field.
    """
    def order(self):
        """
        EXAMPLES::

            sage: G = Sp(4,GF(3))
            sage: G.order()
            51840
            sage: G = SL(4,GF(3))
            sage: G.order()
            12130560
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.order()
            480
            sage: G = MatrixGroup([matrix(ZZ,2,[1,1,0,1])])
            sage: G.order()
            +Infinity
        """
        return integer.Integer(gap(self).Size())

    def random_element(self):
        """
        Return a random element of this group.

        EXAMPLES::

            sage: G = Sp(4,GF(3))
            sage: G.random_element()
            [2 1 1 1]
            [1 0 2 1]
            [0 1 1 0]
            [1 0 0 1]

        ::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.random_element()
            [1 3]
            [0 3]
            sage: G.random_element()
            [2 2]
            [1 0]
            sage: G.random_element()
            [4 0]
            [1 4]
        """
        current_randstate().set_seed_gap()
        from matrix_group_element import MatrixGroupElement
        F = self.field_of_definition()
        return MatrixGroupElement(gap(self).Random()._matrix_(F), self, check=False)

    def random(self):
        """
        Deprecated. Use self.random_element() instead.
        """
        raise NotImplementedError, "Deprecated: use random_element() instead"


    def __contains__(self, x):
        """
        Return True if `x` is an element of this abelian group.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 0], [0, 1]], [[1, 2], [3, 4]]]
            sage: G.order()
            8
            sage: G(1)
            [1 0]
            [0 1]
            sage: G.1 in G
            True
            sage: 1 in G
            True
            sage: [1,2,3,4] in G
            True
            sage: matrix(GF(5),2,[1,2,3,5]) in G
            False
            sage: G(matrix(GF(5),2,[1,2,3,5]))
            Traceback (most recent call last):
            ...
            TypeError: no way to coerce element to self.
        """
        from matrix_group_element import MatrixGroupElement
        if isinstance(x, MatrixGroupElement):
            if x.parent() == self:
                return True
            return gap(x) in gap(self)
        try:
            self(x)
            return True
        except TypeError:
            return False

    def conjugacy_class_representatives(self):
        """
        Return a set of representatives for each of the conjugacy classes
        of the group.

        EXAMPLES::

            sage: G = SU(3,GF(2))
            sage: len(G.conjugacy_class_representatives())
            16
            sage: len(GL(2,GF(3)).conjugacy_class_representatives())
            8
            sage: len(GU(2,GF(5)).conjugacy_class_representatives())
            36
        """
        current_randstate().set_seed_gap()
        try:
            return self.__reps
        except AttributeError:
            pass
        from matrix_group_element import MatrixGroupElement
        G    = self._gap_().ConjugacyClasses()
        reps = list(gap.List(G, 'x -> Representative(x)'))
        F    = self.field_of_definition()
        self.__reps = Sequence([self(g._matrix_(F)) for g in reps], cr=True, universe=self, check=False)
        return self.__reps

    def center(self):
        """
        Return the center of this linear group as a matrix group.

        EXAMPLES::

            sage: G = SU(3,GF(2))
            sage: G.center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators:
             [[[a, 0, 0], [0, a, 0], [0, 0, a]]]
            sage: GL(2,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators:
             [[[2, 0], [0, 2]]]
            sage: GL(3,GF(3)).center()
            Matrix group over Finite Field of size 3 with 1 generators:
             [[[2, 0, 0], [0, 2, 0], [0, 0, 2]]]
            sage: GU(3,GF(2)).center()
            Matrix group over Finite Field in a of size 2^2 with 1 generators:
             [[[a + 1, 0, 0], [0, a + 1, 0], [0, 0, a + 1]]]
        """
        try:
            return self.__center
        except AttributeError:
            pass
        G = list(self._gap_().Center().GeneratorsOfGroup())
        F = self.field_of_definition()
        from sage.groups.matrix_gps.matrix_group import MatrixGroup
        self.__center = MatrixGroup([g._matrix_(F) for g in G])
        return self.__center


class MatrixGroup_gens(MatrixGroup_gap):
    """
    EXAMPLES:

    A ValueError is raised if one of the generators is not invertible.

    ::

        sage: F = GF(5); MS = MatrixSpace(F,2,2)
        sage: G = MatrixGroup([MS.0])
        Traceback (most recent call last):
        ...
        ValueError: each generator must be an invertible matrix but one is not:
        [1 0]
        [0 0]
    """
    def __init__(self, gensG):
        v = Sequence(gensG, immutable=True)
        M = v.universe()
        if not is_MatrixSpace(M):
            raise TypeError, "universe of sequence (=%s) of generators must be a matrix space"%M
        if M.nrows() != M.ncols():
            raise ValueError, "matrices must be square."
        for x in v:
            if not x.is_invertible():
                raise ValueError, "each generator must be an invertible matrix but one is not:\n%s"%x
        self._gensG = v
        MatrixGroup_gap.__init__(self, M.nrows(), M.base_ring())

    def as_permutation_group(self, method =None):
        r"""
        This returns a permutation group representation for the group. In
        most cases occurring in practice, this is a permutation group of
        minimal degree (the degree begin determined from orbits under the
        group action). When these orbits are hard to compute, the procedure
        can be time-consuming and the degree may not be minimal. The
        "method=smaller" option tries return an isomorphic group of lower
        degree.

        EXAMPLES::

            sage: MS = MatrixSpace( GF(2), 5, 5)
            sage: A = MS([[0,0,0,0,1],[0,0,0,1,0],[0,0,1,0,0],[0,1,0,0,0],[1,0,0,0,0]])
            sage: G = MatrixGroup([A])
            sage: G.as_permutation_group()
            Permutation Group with generators [(1,2)]
            sage: MS = MatrixSpace( GF(7), 12, 12)
            sage: GG = gap("ImfMatrixGroup( 12, 3 )")
            sage: GG.GeneratorsOfGroup().Length()
            3
            sage: g1 = MS(eval(str(GG.GeneratorsOfGroup()[1]).replace("\n","")))
            sage: g2 = MS(eval(str(GG.GeneratorsOfGroup()[2]).replace("\n","")))
            sage: g3 = MS(eval(str(GG.GeneratorsOfGroup()[3]).replace("\n","")))
            sage: G = MatrixGroup([g1, g2, g3])
            sage: G.order()
            21499084800
            sage: set_random_seed(0); current_randstate().set_seed_gap()
            sage: G.as_permutation_group()
            Permutation Group with generators [(2,3,5,11,20,38,14,26,43,66)(4,8,16,29,47,63,48,74,56,86)(6,13,24,41,22,27,45,69,98,67)(9,19,35,58,88,94,70,100,84,104)(17,32,52,36,60,75,108,133,101,128)(30,50,77,110,117,106,61,92,123,89)(33,55,83,81,72,103,115,137,112,119)(53,80), (2,4,9,20)(3,6,14,27)(8,17,33,56)(11,22,16,30)(13,19,36,61)(29,48,75,103)(32,53,81,92)(38,63,94,66)(43,67,74,106)(45,70,101,50)(47,72,104,128)(55,84,117,137)(58,89,83,115)(60,86,119,88)(77,108,80,112), (1,2)(3,7)(4,10)(5,12)(6,15)(8,18)(9,21)(11,23)(13,25)(14,28)(16,31)(17,34)(19,37)(20,39)(22,40)(24,42)(26,44)(27,46)(29,49)(30,51)(32,54)(33,57)(35,59)(36,62)(38,64)(41,65)(43,68)(45,71)(47,73)(48,76)(50,78)(52,79)(53,82)(55,85)(56,87)(58,90)(60,91)(61,93)(63,95)(66,96)(67,97)(69,99)(70,102)(72,105)(74,107)(75,109)(77,111)(80,113)(81,114)(83,116)(84,118)(86,120)(88,121)(89,122)(92,124)(94,125)(98,126)(100,127)(101,129)(103,130)(104,131)(106,132)(108,134)(110,135)(112,136)(115,138)(117,139)(119,140)(123,141)(128,142)(133,143)(137,144)]
            sage: set_random_seed(3); current_randstate().set_seed_gap()
            sage: G.as_permutation_group(method="smaller")
            Permutation Group with generators [(1,2)(3,7,13,25,45,5,10,19,35,60)(8,16,30,41,69,11,22,40,31,54)(14,28,49,80,91,20,38,64,93,78)(17,33,57,26,48,23,43,72,36,63)(46,61)(52,73,55,76,94,67,58,70,89,81), (1,3,8,17)(2,5,11,23)(7,14)(10,20)(13,26,49,81)(16,31,55,72)(19,36,64,94)(22,41,70,57)(25,46,76,93)(28,38)(30,52,63,91)(33,58)(35,61,89,80)(40,67,48,78)(43,73)(45,69)(54,60), (1,4)(2,6)(3,9)(5,12)(7,15)(8,18)(10,21)(11,24)(13,27)(14,29)(16,32)(17,34)(19,37)(20,39)(22,42)(23,44)(25,47)(26,50)(28,51)(30,53)(31,56)(33,59)(35,62)(36,65)(38,66)(40,68)(41,71)(43,74)(45,75)(46,77)(48,79)(49,82)(52,83)(54,84)(55,85)(57,86)(58,87)(60,88)(61,90)(63,92)(64,95)(67,96)(69,97)(70,98)(72,99)(73,100)(76,101)(78,102)(80,103)(81,104)(89,105)(91,106)(93,107)(94,108)]

        In this case, the "smaller" option returned an isomorphic group of
        lower degree. The above example used GAP's library of irreducible
        maximal finite ("imf") integer matrix groups to construct the
        MatrixGroup G over GF(7). The section "Irreducible Maximal Finite
        Integral Matrix Groups" in the GAP reference manual has more
        details.
        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        F = self.base_ring()
        if not(F.is_finite()):
            raise NotImplementedError, "Base ring must be finite."
        q = F.order()
        gens = self.gens()
        n = self.degree()
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        W = self.matrix_space().row_space()
        for g in gens:
            p = MS(g.matrix())
            m = p.rows()
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        gap.eval("iso:=IsomorphismPermGroup(Group("+mats_str+"))")
        C = gap("Image( iso )")
        if method == "smaller":
            gap.eval("small:= SmallerDegreePermutationRepresentation( Image( iso ) );")
            C = gap("Image( small )")
        return PermutationGroup(gap_group=C)

    def module_composition_factors(self, method=None):
        r"""
        Returns a list of triples consisting of [base field, dimension,
        irreducibility], for each of the Meataxe composition factors
        modules. The method="verbose" option returns more information, but
        in Meataxe notation.

        EXAMPLES::

            sage: F=GF(3);MS=MatrixSpace(F,4,4)
            sage: M=MS(0)
            sage: M[0,1]=1;M[1,2]=1;M[2,3]=1;M[3,0]=1
            sage: G = MatrixGroup([M])
            sage: G.module_composition_factors()
            [[Finite Field of size 3, 1, True],
             [Finite Field of size 3, 1, True],
             [Finite Field of size 3, 2, True]]
            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.module_composition_factors()
            [[Finite Field of size 7, 2, True]]

        Type "G.module_composition_factors(method='verbose')" to get a
        more verbose version.

        For more on MeatAxe notation, see
        http://www.gap-system.org/Manuals/doc/htm/ref/CHAP067.htm
        """
        from sage.misc.sage_eval import sage_eval
        F = self.base_ring()
        if not(F.is_finite()):
            raise NotImplementedError, "Base ring must be finite."
        q = F.order()
        gens = self.gens()
        n = self.degree()
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        W = self.matrix_space().row_space()
        for g in gens:
            p = MS(g.matrix())
            m = p.rows()
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        gap.eval("MCFs := MTX.CompositionFactors( M )")
        N = eval(gap.eval("Length(MCFs)"))
        if method == "verbose":
            print gap.eval('MCFs')+"\n"
        L = []
        for i in range(1,N+1):
            gap.eval("MCF := MCFs[%s]"%i)
            L = L + [[sage_eval(gap.eval("MCF.field")), eval(gap.eval("MCF.dimension")), sage_eval(gap.eval("MCF.IsIrreducible"))]]
        return L

    def gens(self):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: gens[0] in G
            True
            sage: gens = G.gens()
            sage: gens[0] in G
            True
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]

        ::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: G = MatrixGroup([MS(1), MS([1,2,3,4])])
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 0], [0, 1]], [[1, 2], [3, 4]]]
            sage: G.gens()
            [[1 0]
            [0 1], [1 2]
            [3 4]]
        """
        try:
            return self.__gens
        except AttributeError:
            t = Sequence([MatrixGroupElement(x, self) for x in self._gensG],
                         immutable=True, universe=self)
            self.__gens = t
            return t

    def _gap_init_(self):
        """
        Returns a string representation of the corresponding GAP object.

        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G._gap_init_() # The variable $sage11 belongs to gap(F) and is somehow random
            'Group([[Z(5)^0,Z(5)^1],[Z(5)^2,Z(5)^0]]*One($sage11),[[Z(5)^0,Z(5)^0],[0*Z(5),Z(5)^0]]*One($sage11))'
            sage: gap(G._gap_init_())
            Group([ [ [ Z(5)^0, Z(5) ], [ Z(5)^2, Z(5)^0 ] ],
              [ [ Z(5)^0, Z(5)^0 ], [ 0*Z(5), Z(5)^0 ] ] ])
        """
        gens_gap = ','.join([x._gap_init_() for x in self._gensG])
        return 'Group(%s)'%gens_gap

    def _repr_(self):
        """
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: G
            Matrix group over Finite Field of size 5 with 2 generators:
             [[[1, 2], [4, 1]], [[1, 1], [0, 1]]]
        """
        gns = [x.list() for x in self.gens()]
        return "Matrix group over %s with %s generators: \n %s"%(self.base_ring(), self.ngens(), gns)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: latex(G)
            \left\langle \left(\begin{array}{rr}
            1 & 2 \\
            4 & 1
            \end{array}\right), \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right) \right\rangle
        """
        gens = ', '.join([latex(x) for x in self.gens()])
        return '\\left\\langle %s \\right\\rangle'%gens

    def invariant_generators(self):
        """
        Wraps Singular's invariant_algebra_reynolds and invariant_ring
        in finvar.lib, with help from Simon King and Martin Albrecht.
        Computes generators for the polynomial ring
        `F[x_1,\ldots,x_n]^G`, where G in GL(n,F) is a finite
        matrix group.

        In the "good characteristic" case the polynomials returned form a
        minimal generating set for the algebra of G-invariant polynomials.
        In the "bad" case, the polynomials returned are primary and
        secondary invariants, forming a not necessarily minimal generating
        set for the algebra of G-invariant polynomials.

        EXAMPLES::

            sage: F = GF(7); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[0,1],[-1,0]]),MS([[1,1],[2,3]])]
            sage: G = MatrixGroup(gens)
            sage: G.invariant_generators()
            [x1^7*x2 - x1*x2^7, x1^12 - 2*x1^9*x2^3 - x1^6*x2^6 + 2*x1^3*x2^9 + x2^12, x1^18 + 2*x1^15*x2^3 + 3*x1^12*x2^6 + 3*x1^6*x2^12 - 2*x1^3*x2^15 + x2^18]
            sage: q = 4; a = 2
            sage: MS = MatrixSpace(QQ, 2, 2)
            sage: gen1 = [[1/a,(q-1)/a],[1/a, -1/a]]; gen2 = [[1,0],[0,-1]]; gen3 = [[-1,0],[0,1]]
            sage: G = MatrixGroup([MS(gen1),MS(gen2),MS(gen3)])
            sage: G.order()
            12
            sage: G.invariant_generators()
            [x1^2 + 3*x2^2, x1^6 + 15*x1^4*x2^2 + 15*x1^2*x2^4 + 33*x2^6]
            sage: F = GF(5); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,2],[-1,1]]),MS([[1,1],[-1,1]])]
            sage: G = MatrixGroup(gens)
            sage: G.invariant_generators()  ## takes a long time (several mins)
            [x1^20 + x1^16*x2^4 + x1^12*x2^8 + x1^8*x2^12 + x1^4*x2^16 + x2^20, x1^20*x2^4 + x1^16*x2^8 + x1^12*x2^12 + x1^8*x2^16 + x1^4*x2^20]
            sage: F=CyclotomicField(8)
            sage: z=F.gen()
            sage: a=z+1/z
            sage: b=z^2
            sage: MS=MatrixSpace(F,2,2)
            sage: g1=MS([[1/a,1/a],[1/a,-1/a]])
            sage: g2=MS([[1,0],[0,b]])
            sage: g3=MS([[b,0],[0,1]])
            sage: G=MatrixGroup([g1,g2,g3])
            sage: G.invariant_generators()
            [x1^8 + 14*x1^4*x2^4 + x2^8,
             x1^24 + 10626/1025*x1^20*x2^4 + 735471/1025*x1^16*x2^8 + 2704156/1025*x1^12*x2^12 + 735471/1025*x1^8*x2^16 + 10626/1025*x1^4*x2^20 + x2^24]

        AUTHORS:

        - David Joyner, Simon King and Martin Albrecht.

        REFERENCES:

        - Singular reference manual

        - B. Sturmfels, "Algorithms in invariant theory", Springer-Verlag, 1993.

        - S. King, "Minimal Generating Sets of non-modular invariant rings of
          finite groups", arXiv:math.AC/0703035
        """
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.interfaces.singular import singular
        gens = self.gens()
        singular.LIB("finvar.lib")
        n = len((gens[0].matrix()).rows())
        F = self.base_ring()
        q = F.characteristic()
        ## test if the field is admissible
        if F.gen()==1: # we got the rationals or GF(prime)
            FieldStr = str(F.characteristic())
        elif hasattr(F,'polynomial'): # we got an algebraic extension
            if len(F.gens())>1:
                raise NotImplementedError, "can only deal with finite fields and (simple algebraic extensions of) the rationals"
            FieldStr = '(%d,%s)'%(F.characteristic(),str(F.gen()))
        else: # we have a transcendental extension
            FieldStr = '(%d,%s)'%(F.characteristic(),','.join([str(p) for p in F.gens()]))

        ## Setting Singular's variable names
        ## We need to make sure that field generator and variables get different names.
        if str(F.gen())[0]=='x':
            VarStr = 'y'
        else:
            VarStr = 'x'
        VarNames='('+','.join((VarStr+str(i+1) for i in range(n)))+')'
        R=singular.ring(FieldStr,VarNames,'dp')
        if hasattr(F,'polynomial') and F.gen()!=1: # we have to define minpoly
            singular.eval('minpoly = '+str(F.polynomial()).replace('x',str(F.gen())))
        A = [singular.matrix(n,n,str((x.matrix()).list())) for x in gens]
        Lgens = ','.join((x.name() for x in A))
        PR = PolynomialRing(F,n,[VarStr+str(i) for i in range(1,n+1)])
        if q == 0 or (q > 0 and self.order()%q != 0):
            ReyName = 't'+singular._next_var_name()
            singular.eval('list %s=group_reynolds((%s))'%(ReyName,Lgens))
            IRName = 't'+singular._next_var_name()
            singular.eval('matrix %s = invariant_algebra_reynolds(%s[1])'%(IRName,ReyName))
            OUT = [singular.eval(IRName+'[1,%d]'%(j)) for j in range(1,1+singular('ncols('+IRName+')'))]
            return [PR(gen) for gen in OUT]
        if self.order()%q == 0:
            PName = 't'+singular._next_var_name()
            SName = 't'+singular._next_var_name()
            singular.eval('matrix %s,%s=invariant_ring(%s)'%(PName,SName,Lgens))
            OUT = [singular.eval(PName+'[1,%d]'%(j)) for j in range(1,1+singular('ncols('+PName+')'))] + [singular.eval(SName+'[1,%d]'%(j)) for j in range(2,1+singular('ncols('+SName+')'))]
            return [PR(gen) for gen in OUT]


class MatrixGroup_gens_finite_field(MatrixGroup_gens, MatrixGroup_gap_finite_field):
    pass

##     def conjugacy_class_representatives_gap(self):
##         """
##         Wraps GAP Representative+ConjugactClasses but returns a list of
##         strings representing the GAP matrices which form a complete
##         set of representatives of the conjugacy classes of the group.

##         EXAMPLES:
##             sage: F = GF(3); MS = MatrixSpace(F,2,2)
##             sage: gens = [MS([[1,0],[-1,1]]),MS([[1,1],[0,1]])]
##             sage: G = MatrixGroup(gens)
##             sage: G.conjugacy_class_representatives_gap()
## 	    ['[ [ Z(3)^0, 0*Z(3) ], [ 0*Z(3), Z(3)^0 ] ]',
##  	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3)^0 ] ]',
##             '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), Z(3) ] ]',
##             '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3)^0 ] ]',
##      	    '[ [ 0*Z(3), Z(3) ], [ Z(3)^0, Z(3) ] ]',
##      	    '[ [ 0*Z(3), Z(3)^0 ], [ Z(3), 0*Z(3) ] ]',
##  	    '[ [ Z(3), 0*Z(3) ], [ 0*Z(3), Z(3) ] ]']

##         AUTHOR: David Joyner (1-2006)
##         """
##         F = self.__R
##         deg = self.__n
##         gp_gens = self.gens()
##         L = [gap(A) for A in gp_gens]
##         sL = ','.join(str(x) for x in L)
##         if is_FiniteField(F):
##             q = F.order()
##             gap.eval("cl:=ConjugacyClasses(Group(["+sL+"]))")
##             m = eval(gap.eval("Length(cl)"))
##             gap.eval("reps:=List(cl,x->Representative(x))")
##             sreps = [gap.eval("reps["+str(i+1)+"]") for i in range(m)]
##             return sreps
##         raise TypeError, "R (=%s) must be a finite field"%R



