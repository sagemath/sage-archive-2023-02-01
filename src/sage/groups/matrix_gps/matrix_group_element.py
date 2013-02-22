"""
Matrix Group Elements

AUTHORS:

- David Joyner (2006-05): initial version David Joyner

- David Joyner (2006-05): various modifications to address William
  Stein's TODO's.

- William Stein (2006-12-09): many revisions.

EXAMPLES::

    sage: F = GF(3); MS = MatrixSpace(F,2,2)
    sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
    sage: G = MatrixGroup(gens); G
    Matrix group over Finite Field of size 3 with 2 generators:
     [[[1, 0], [0, 1]], [[1, 1], [0, 1]]]
    sage: g = G([[1,1],[0,1]])
    sage: h = G([[1,2],[0,1]])
    sage: g*h
    [1 0]
    [0 1]

You cannot add two matrices, since this is not a group operation.
You can coerce matrices back to the matrix space and add them
there::

    sage: g + h
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand type(s) for +: 'MatrixGroup_gens_finite_field_with_category.element_class' and 'MatrixGroup_gens_finite_field_with_category.element_class'

::

    sage: g.matrix() + h.matrix()
    [2 0]
    [0 2]

::

    sage: 2*g
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '*': 'Integer Ring' and 'Matrix group over Finite Field of size 3 with 2 generators:
     [[[1, 0], [0, 1]], [[1, 1], [0, 1]]]'
    sage: 2*g.matrix()
    [2 2]
    [0 2]
"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import Integer, Infinity
from sage.interfaces.gap import gap
import sage.structure.element as element
from sage.matrix.matrix import Matrix
from sage.structure.factorization import Factorization
from sage.structure.sage_object import have_same_parent

def is_MatrixGroupElement(x):
    return isinstance(x, MatrixGroupElement)

class MatrixGroupElement(element.MultiplicativeGroupElement):
    """
    An element of a matrix group.

    EXAMPLES::

        sage: F = GF(3); MS = MatrixSpace(F,2,2)
               sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
        sage: G = MatrixGroup(gens)
        sage: g = G.random_element()
        sage: type(g)
        <class 'sage.groups.matrix_gps.matrix_group_element.MatrixGroup_gens_finite_field_with_category.element_class'>
    """
    def __init__(self, g, parent, check = True):
        r"""
        Create element of a matrix group.

        INPUT:


        -  ``g`` - a Matrix

        -  ``parent`` - defines parent group (g must be in
           parent or a TypeError is raised).

        -  ``check`` - bool (default: True), if true does some
           type checking.


        TESTS::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.random_element()
            sage: TestSuite(g).run()
        """
        if check:
            if not isinstance(g, Matrix):
                raise TypeError, "g must be a matrix"
            if not (g.parent() is parent):
                g = parent.matrix_space()(g)
        element.Element.__init__(self, parent)
        self.__mat = g

    def matrix(self):
        """
        Obtain the usual matrix (as an element of a matrix space)
        associated to this matrix group element.

        One reason to compute the associated matrix is that matrices
        support a huge range of functionality.

        EXAMPLES::

            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: g = G.0
            sage: g.matrix()
            [1 1]
            [0 1]
            sage: parent(g.matrix())
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 7

        Matrices have extra functionality that matrix group elements do not
        have.

        ::

            sage: g.matrix().charpoly('t')
            t^2 + 5*t + 1
        """
        return self.__mat.__copy__()

    def _gap_init_(self):
        """
        Return a string representation of a Gap object corresponding to
        this matrix group element.

        EXAMPLES::

            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])]); g = G.1
            sage: g._gap_init_() # The variable $sage27 belongs to gap(k) and is somehow random
            '[[Z(7)^0,0*Z(7)],[0*Z(7),Z(7)^2]]*One($sage27)'
            sage: gap(g._gap_init_())
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]

        It may be better to use gap(the matrix), since the result is
        cached.

        ::

            sage: gap(G.1)
            [ [ Z(7)^0, 0*Z(7) ], [ 0*Z(7), Z(7)^2 ] ]
            sage: gap(G.1).IsMatrix()
            true
        """
        return self.__mat._gap_init_()

    def _gap_latex_(self):
        r"""
        Return the GAP latex version of this matrix.

        AUTHORS:

        - S. Kohl: Wrote the GAP function.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g._gap_latex_()
            \left(\begin{array}{rr}%
            Z(3)^{0}&Z(3)^{0}\\%
            0*Z(3)&Z(3)^{0}\\%
            \end{array}\right)%

        Type view(g._latex_()) to see the object in an xdvi window
        (assuming you have latex and xdvi installed).
        """
        s1 = self.__mat._gap_init_()
        s2 = gap.eval("LaTeX("+s1+")")
        return eval(s2)

    def _repr_(self):
        """
        Return string representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: g                        # indirect doctest
            [1 1]
            [0 1]
        """
        return str(self.__mat)

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G([[1, 1], [0, 1]])
            sage: print g._latex_()
            \left(\begin{array}{rr}
            1 & 1 \\
            0 & 1
            \end{array}\right)

        Type ``view(g._latex_())`` to see the object in an
        xdvi window (assuming you have latex and xdvi installed).
        """
        return self.__mat._latex_()

    def _mul_(self,other):
        """
        Return the product of self and other, which must have identical
        parents.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g*h         # indirect doctest
            [1 2]
            [0 1]
        """
        parent = self.parent()
        return parent.element_class(self.__mat * other.__mat, parent, check=False)

    def _act_on_(self, x, self_on_left):
        """
        EXAMPLES::

            sage: G = GL(4,7)
            sage: G.0 * vector([1,2,3,4])
            (3, 2, 3, 4)
            sage: v = vector(GF(7), [3,2,1,-1])
            sage: g = G.1
            sage: v * g == v * g.matrix()   # indirect doctest
            True
        """
        if not is_MatrixGroupElement(x) and x not in self.parent().base_ring():
            try:
                if self_on_left:
                    return self.matrix() * x
                else:
                    return x * self.matrix()
            except TypeError:
                return None

    def __invert__(self):
        parent = self.parent()
        return parent.element_class(~self.__mat, parent, check=False)

    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = G([1,1, 0,1])
            sage: g == h
            True
        """
        return cmp(self.__mat, other.__mat)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: H = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = H([1,1, 0,1])
            sage: g == h
            True
        """
        return have_same_parent(self, other) and self.__mat == other.__mat

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2)
            sage: gens = [MS([1,0, 0,1]), MS([1,1, 0,1])]
            sage: G = MatrixGroup(gens)
            sage: H = MatrixGroup(gens)
            sage: g = G([1,1, 0,1])
            sage: h = H([1,1, 0,1])
            sage: g != h
            False
        """
        return not self.__eq__(other)

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer `n` such that `g^n = 1`, or
        +Infinity if no such integer exists.

        EXAMPLES::

            sage: k = GF(7); G = MatrixGroup([matrix(k,2,[1,1,0,1]), matrix(k,2,[1,0,0,2])])
            sage: G
            Matrix group over Finite Field of size 7 with 2 generators:
             [[[1, 1], [0, 1]], [[1, 0], [0, 2]]]
            sage: G.order()
            21

        See trac #1170

        ::

            sage: gl=GL(2,ZZ); gl
            General Linear Group of degree 2 over Integer Ring
            sage: g=gl.gens()[2]; g
            [1 1]
            [0 1]
            sage: g.order()
            +Infinity
        """
        try:
            return self.__order
        except AttributeError:
            try:
                self.__order = Integer(self._gap_().Order())
            except TypeError:
                self.__order = Infinity
            return self.__order

    def word_problem(self, gens=None):
        r"""
        Write this group element in terms of the elements of the list
        ``gens``, or the standard generators of the parent of self if
        ``gens`` is None.

        INPUT:

        - ``gens`` - a list of elements that can be coerced into the parent of
          self, or None.

        ALGORITHM: Use GAP, which has optimized algorithms for solving the
        word problem (the GAP functions EpimorphismFromFreeGroup and
        PreImagesRepresentative).

        EXAMPLE::

            sage: G = GL(2,5); G
            General Linear Group of degree 2 over Finite Field of size 5
            sage: G.gens()
            [
            [2 0]
            [0 1],
            [4 1]
            [4 0]
            ]
            sage: G(1).word_problem([G.1, G.0])
            1

        Next we construct a more complicated element of the group from the
        generators::

            sage: s,t = G.0, G.1
            sage: a = (s * t * s); b = a.word_problem(); b
            ([2 0]
            [0 1]) *
            ([4 1]
            [4 0]) *
            ([2 0]
            [0 1])
            sage: b.prod() == a
            True

        We solve the word problem using some different generators::

            sage: s = G([2,0,0,1]); t = G([1,1,0,1]); u = G([0,-1,1,0])
            sage: a.word_problem([s,t,u])
            ([2 0]
            [0 1])^-1 *
            ([1 1]
            [0 1])^-1 *
            ([0 4]
            [1 0]) *
            ([2 0]
            [0 1])^-1

        We try some elements that don't actually generate the group::

            sage: a.word_problem([t,u])
            Traceback (most recent call last):
            ...
            ValueError: Could not solve word problem

        AUTHORS:

        - David Joyner and William Stein
        - David Loeffler (2010): fixed some bugs
        """
        G = self.parent()
        gg = gap(G)
        if not gens:
            gens = gg.GeneratorsOfGroup()
        else:
            gens = gap([G(x) for x in gens])
        H = gens.Group()
        hom = H.EpimorphismFromFreeGroup()
        in_terms_of_gens = str(hom.PreImagesRepresentative(self))
        if 'identity' in in_terms_of_gens:
            return Factorization([])
        elif in_terms_of_gens == 'fail':
            raise ValueError, "Could not solve word problem"
        v = list(H.GeneratorsOfGroup())
        F = G.field_of_definition()
        w = [G(x._matrix_(F)) for x in v]
        factors = [Factorization([(x,1)]) for x in w]
        d = dict([('x%s'%(i+1),factors[i]) for i in range(len(factors))])

        from sage.misc.sage_eval import sage_eval
        F = sage_eval(str(in_terms_of_gens), d)
        F._set_cr(True)
        return F

    def list(self):
        """
        Return list representation of this matrix.

        EXAMPLES::

            sage: F = GF(3); MS = MatrixSpace(F,2,2)
            sage: gens = [MS([[1,0],[0,1]]),MS([[1,1],[0,1]])]
            sage: G = MatrixGroup(gens)
            sage: g = G.0
            sage: g.list()
            [[1, 0], [0, 1]]
        """
        rws = (self.__mat).rows()
        lrws = [list(x) for x in rws]
        return lrws

    def conjugacy_class(self):
        r"""
        Returns the conjugacy class of ``self``

        EXAMPLES::

            sage: G = SL(2, GF(2))
            sage: g = G.gens()[0]
            sage: g.conjugacy_class()
            Conjugacy class of [1 1]
            [0 1] in Special Linear Group of degree 2 over Finite Field of size 2
        """
        from sage.groups.conjugacy_classes import ConjugacyClassGAP
        return ConjugacyClassGAP(self.parent(), self)
