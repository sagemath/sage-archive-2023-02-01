"""
Contains general base classes for the classical groups
GL, SL, Sp, SO, SU

TODO:
   Implement "twisted" groups.

AUTHORS:
   William Stein -- initial version
   David Joyner  -- degree, base_ring, random, order methods; examples
   David Joyner (2006-05) -- added center, more examples,
                             renamed random attributes, bug fixes.

REFERENCES:
    [KL] Peter Kleidman and Martin Liebeck. The subgroup structure of the finite
    classical groups. Cambridge University Press, 1990.
    [C] R. W. Carter. Simple groups of Lie type, volume 28 of Pure and
    Applied Mathematics. John Wiley and Sons, 1972.

"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.group import Group
from sage.rings.all import IntegerRing, is_Ring, Integer
from sage.interfaces.gap import gap
from sage.rings.finite_field import FiniteField

class LinearGroup_generic(Group):
    def __init__(self, n, R, e=0):
        """
        n -- the degree
        R -- the base ring
        e -- a parameter for orthogonal groups only depending
             on the invariant form



        """
        self.__n = Integer(n)
        self.__R = R
        self.__form = e
        if not is_Ring(R):
            raise TypeError, "R (=%s) must be a ring"%R

    def degree(self):
        return self.__n

    def field_of_definition(self):
        """
        This is only used for unitary groups at the moment.
        It can eventually be used for inner forms.

        EXAMPLES:
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
        from sage.groups.matrix_gps.unitary import SpecialUnitaryGroup_finite_field,GeneralUnitaryGroup_finite_field
        if type(self)==SpecialUnitaryGroup_finite_field:
            q = (self.__R).order()
            return FiniteField(q**2)
        if type(self)==GeneralUnitaryGroup_finite_field:
            q = (self.__R).order()
            return FiniteField(q**2)
        return self.__R

    def base_ring(self):
        return self.__R

    def invariant_form(self):
        return self.__form

    def is_finite(self):
        """
        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: G.is_finite()
            True
        """
        return self.base_ring().is_finite()

class LinearGroup_finite_field(LinearGroup_generic):
    def order(self):
        """
        EXAMPLES:
            sage: G = Sp(4,GF(3))
	    sage: G.order()
            51840
            sage: G = SL(4,GF(3))
            sage: G.order()
            12130560

        """
        cmd = self._gap_init_()
        return eval(gap.eval("Size("+cmd+")"))

    def random(self):
        """
        Wraps GAP's Random function.

        EXAMPLES:
            sage: G = Sp(4,GF(3))
            sage: G.random()        ## random output
            [0 2 1 1]
            [1 1 1 2]
            [2 2 2 2]
            [0 2 0 2]

        """
        from matrix_group_element import MatrixGroupElement
        F = self.field_of_definition()
        cmd = self._gap_init_()
        s = gap("Random("+cmd+")")
        return MatrixGroupElement(s._matrix_(F),self, check=False)

    def random_gap(self):
        """
        Return a random element of this group, using GAP with
        the output in GAP notation.

        EXAMPLES:
            sage: G = Sp(4,GF(3))
            sage: G.random_gap()        ## random output
            [ [ Z(3), 0*Z(3), Z(3), 0*Z(3) ], [ Z(3), Z(3), Z(3)^0, Z(3) ],
              [ Z(3), 0*Z(3), Z(3)^0, 0*Z(3) ], [ Z(3)^0, Z(3), Z(3)^0, Z(3)^0 ] ]

        """
        from matrix_group_element import MatrixGroupElement
        #print self._gap_().Random()
        return MatrixGroupElement(self._gap_().Random(),
                                       self, check=False)

    def __contains__(self, x):
        """
        Return True if $x$ is an element of this group.

        EXAMPLES:
            sage: G = GL(3,GF(4))
            sage: g = G.random()
            sage: g in G
            True
        """
        from matrix_group_element import MatrixGroupElement
        return isinstance(x, MatrixGroupElement) and x.parent() == self

    def conjugacy_class_representatives(self):
        """
        Return a complete set of representatives of the conjugacy
        classes of the group.

        WARNING: Does not work for unitary groups.

        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: C = G.conjugacy_class_representatives()
            sage: len(C)
            8
            sage: C[0]
            [1 0]
            [0 1]
            sage: [g.list() for g in C]     # prints more nicely
            [[[1, 0], [0, 1]],
             [[0, 2], [1, 1]],
             [[2, 0], [0, 2]],
             [[0, 2], [1, 2]],
             [[0, 2], [1, 0]],
             [[0, 1], [1, 2]],
             [[0, 1], [1, 1]],
             [[2, 0], [0, 1]]]
            sage: G = GL(2,GF(4))
            sage: C = G.conjugacy_class_representatives()
            sage: [g.list() for g in C]      # prints more nicely
            [[[1, 0], [0, 1]],
             [[0, 1], [1, 0]],
             [[a, 0], [0, a]],
             [[0, a + 1], [1, 0]],
             [[a + 1, 0], [0, a + 1]],
             [[0, a], [1, 0]],
             [[0, 1], [1, a]],
             [[0, 1], [1, a + 1]],
             [[0, a], [1, 1]],
             [[0, a], [1, a]],
             [[0, a + 1], [1, 1]],
             [[0, a + 1], [1, a + 1]],
             [[1, 0], [0, a]],
             [[1, 0], [0, a + 1]],
             [[a, 0], [0, a + 1]]]

        """
        from matrix_group_element import MatrixGroupElement
        G = self._gap_()
        C = G.ConjugacyClasses()
        gap = G.parent()
        reps = gap.List(C, 'x -> Representative(x)')
        K = self.base_ring()
        self.__reps = [MatrixGroupElement(x._matrix_(K),self, check=False) for x in reps]
        return self.__reps

    def center(self):
        """
        Return the center the group (wraps GAP's Center function).

        WARNING: Does not work for unitary groups.

        EXAMPLES:
            sage: G = GL(2,GF(3))
            sage: C = G.center()
            [[1 0]
             [0 1], [2 0]
                    [0 2]]
            sage: C
            Matrix group over Finite Field of size 3 with 2 generators:
             [[[1, 0], [0, 1]], [[2, 0], [0, 2]]]
            sage: print C
            MatrixGroup( [[[1, 0], [0, 1]], [[2, 0], [0, 2]]] )

        """
        from sage.misc.sage_eval import sage_eval as seval
        from matrix_group import MatrixGroup
        F = self.base_ring()
        G = self._gap_()
        C = G.Center()
        reps = C.Elements()
        repns = [x._matrix_(F) for x in reps]
        print repns
        return MatrixGroup(repns)
