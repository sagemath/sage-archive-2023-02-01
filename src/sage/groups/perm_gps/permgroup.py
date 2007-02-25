r"""
Permutation groups

A {\it permutation group} is a finite group G whose elements are permutations
of a given finite set X (i.e., bijections X --> X) and whose group operation is
the composition of permutations. The number of elements of $X$ is called the
{\it degree} of G.

In \sage a permutation is represented as either a string that defines a
permutation using disjoint cycle notation, or a list of tuples, which represent
disjoint cycles.

\begin{verbatim}
(a,...,b)(c,...,d)...(e,...,f)  <--> [(a,...,b), (c,...,d),..., (e,...,f)]
                  () = identity <--> []
\end{verbatim}

You can construct the following permutation groups:

-- SymmetricGroup, $S_n$ of order $n!$ (n can also be a list $X$ of distinct
                   positive integers, in which case it returns $S_X$)

-- AlternatingGroup, $A_n$ or order $n!/2$ (n can also be a list $X$
                   of distinct positive integers, in which case it returns
                   $A_X$)

-- DihedralGroup, $D_n$ of order $2n$

-- CyclicPermutationGroup, $C_n$ of order $n$

-- TransitiveGroup, $i^{th}$ transitive group of degree $n$
                      from the GAP tables of transitive groups (requires
                      the "optional" package database_gap)

-- Mathieu groups, of degrees 9, 10, 11, 12, 21, 22, 23, or 24.

-- KleinFourGroup, subgroup of $S_4$ of order $4$ which is not $C_2 \times C_2$

-- PGL(n,q), projective general linear group of $n\times n$ matrices over
             the finite field GF(q)

-- PSL(n,q), projective special linear group of $n\times n$ matrices over
             the finite field GF(q)

-- PSp(2n,q), projective symplectic linear group of $2n\times 2n$ matrices
              over the finite field GF(q)

-- PSU(n,q), projective special unitary group of $n\times n$ matrices having
             coefficients in the finite field $GF(q^2)$ that respect a
             fixed nondegenerate sesquilinear form, of determinant 1.

-- PGU(n,q), projective general unitary group of $n\times n$ matrices having
             coefficients in the finite field $GF(q^2)$ that respect a
             fixed nondegenerate sesquilinear form, modulo the centre.

-- Suzuki(q), Suzuki group over GF(q), $^2 B_2(2^{2k+1}) = Sz(2^{2k+1})$.

-- direct_product_permgroups, which takes a list of permutation
             groups and returns their direct product.

JOKE:
    Q: What's hot, chunky, and acts on a polygon? A: Dihedral soup.
    Renteln, P. and Dundes, A. "Foolproof: A Sampling of Mathematical
    Folk Humor." Notices Amer. Math. Soc. 52, 24--34, 2005.

AUTHOR:
    - David Joyner (2005-10-14): first version
    - David Joyner (2005-11-17)
    - William Stein (2005-11-26): rewrite to better wrap Gap
    - David Joyner (2005-12-21)
    - Stein and Joyner (2006-01-04): added conjugacy_class_representatives
    - David Joyner (2006-03): reorganization into subdirectory perm_gps;
                              added __contains__, has_element; fixed _cmp_;
                              added subgroup class+methods, PGL,PSL,PSp, PSU classes,
    - David Joyner (2006-06): added PGU, functionality to SymmetricGroup, AlternatingGroup,
                                    direct_product_permgroups
    - David Joyner (2006-08): added degree, ramification_module_decomposition_modular_curve
                              and ramification_module_decomposition_hurwitz_curve
                              methods to PSL(2,q), MathieuGroup, is_isomorphic
    - Bobby Moretti (2006)-10): Added KleinFourGroup, fixed bug in DihedralGroup
    - David Joyner (2006-10): added is_subgroup (fixing a bug found by Kiran Kedlaya),
                              is_solvable, normalizer, is_normal_subgroup, Suzuki
    - David Kohel (2007-02): fixed __contains__ to not enumerate group elements,
                              following the convention for __call__

REFERENCES:
    Cameron, P., Permutation Groups. New York: Cambridge University Press, 1999.
    Wielandt, H., Finite Permutation Groups. New York: Academic Press, 1964.
    Dixon, J. and Mortimer, B., Permutation Groups, Springer-Verlag, Berlin/New York, 1996.

TODO:
  Wrap Ree groups (using IsPermGroup filter)?

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                          David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import random
from types import ListType

import sage.structure.element as element
import sage.groups.group as group

from sage.rings.all      import RationalField, Integer
#from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
import sage.structure.coerce as coerce
from sage.rings.finite_field import GF
from sage.rings.arith import factor
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.misc.functional import is_even, log

def gap_format(x):
    """
    Put a permutation in Gap format, as a string.
    """
    if isinstance(x,str): return x
    x = str(tuple(x)).replace(' ','')
    return x.replace('),(',')(').replace(',)',')')

def PermutationGroup(x, from_group=False, check=True):
    """
    Return the permutation group associated to $x$ (typically a list
    of generators).

    EXAMPLES:
        sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5), (3,4)]

    We can also make permutation groups from PARI groups:
        sage: H = pari('x^4 - 2*x^3 - 2*x + 1').polgalois()
        sage: G = PariGroup(H, 4); G
        PARI group [8, -1, 3, "D(4)"] of degree 4
        sage: H = PermutationGroup(G); H          # requires optional database_gap
        Transitive group number 3 of degree 4
        sage: H.gens()                            # requires optional database_gap
        ((1,2,3,4), (1,3))

    EXAMPLES:
    There is an underlying gap object that implements each permutation group.

        sage: G = PermutationGroup([[(1,2,3,4)]])
        sage: G._gap_()
        Group([ (1,2,3,4) ])
        sage: gap(G)
        Group([ (1,2,3,4) ])
        sage: gap(G) is G._gap_()
        True
        sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
        sage: G._gap_().DerivedSeries()   # output somewhat random
        [ Group([ (1,2,3)(4,5), (3,4) ]), Group([ (1,5)(3,4), (1,5)(2,4), (1,4,5) ]) ]
    """
    if not is_ExpectElement(x) and hasattr(x, '_permgroup_'):
        return x._permgroup_()
    return PermutationGroup_generic(x, from_group, check)


class PermutationGroup_generic(group.FiniteGroup):
    """
    EXAMPLES:
        sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5), (3,4)]
        sage: G.center()
        Permutation Group with generators [()]
        sage: G.group_id()          # requires optional database_gap
        [120, 34]
        sage: n = G.order(); n
        120
        sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
        sage: loads(G.dumps()) == G
        True
    """
    def __init__(self, gens, from_group = False, check=True):
        if from_group and isinstance(gens, str):
            self.__gap = gens
            self.gens()  # so will check that group can be defined in GAP (e.g., no missing packages, etc.)
            return
        if is_GapElement(gens):
            if from_group:
                # the variable gens is actually a Gap group or string
                # representation of one.
                self.__gap = str(gens)
            else:
                # gens is a Gap object that represents a list of generators for a group
                self.__gap = 'Group(%s)'%gens
            self.gens() # so will check that group can be defined in GAP (e.g., no missing packages, etc.)
            return

        if check:
            if isinstance(gens, tuple):
                gens = list(gens)
            elif not isinstance(gens, list):
                raise TypeError, "gens must be a tuple or list"
            gens = [gap_format(x) for x in gens]

        cmd = 'Group(%s)'%gens
        cmd = cmd.replace("'","")  # get rid of quotes
        self.__gap = cmd
        self.gens()

    def _gap_init_(self):
        return self.__gap

    def _magma_init_(self):
        g = str(self.gens())[1:-1]
        return 'PermutationGroup<%s | %s>'%(self.degree(), g)

    def __cmp__(self, right):
        """
        Compare self and right.

        The ordering is whatever it is in Gap.

        EXAMPLES:
            sage: G1 = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G2 = PermutationGroup([[(1,2,3),(4,5)]])
            sage: G1 > G2
            True
            sage: G1 < G2
            False
        """
        if not isinstance(right, PermutationGroup_generic):
            return -1
        return right._gap_().__cmp__(self._gap_())


    def __call__(self, x):
        """
        Coerce x into this permutation group.

        The input can be either a string that defines a permutation in
        cycle notation, a permutation group element, a list of
        integers that gives the permutation as a mapping, a list of
        tuples, or the integer 1.

        EXAMPLES:
        We illustrate each way to make a permutation in S4:
            sage: G = SymmetricGroup(4)
            sage: G((1,2,3,4))
            (1,2,3,4)
            sage: G([(1,2),(3,4)])
            (1,2)(3,4)
            sage: G('(1,2)(3,4)')
            (1,2)(3,4)
            sage: G([1,2,4,3])
            (3,4)
            sage: G([2,3,4,1])
            (1,2,3,4)
            sage: G(G((1,2,3,4)))
            (1,2,3,4)
            sage: G(1)
            ()

        Some more examples:
            sage: G = PermutationGroup([(1,2,3,4)])
            sage: G([(1,3), (2,4)])
            (1,3)(2,4)
            sage: G(G.0^3)
            (1,4,3,2)
            sage: G(1)
            ()
            sage: G((1,4,3,2))
            (1,4,3,2)
            sage: G([(1,2)])
            Traceback (most recent call last):
            ...
            TypeError: permutation (1,2) not in Permutation Group with generators [(1,2,3,4)]
        """
        if isinstance(x, PermutationGroupElement):
            if x.parent() is self:
                return x
            else:
                return PermutationGroupElement(x._gap_(), self, check = True)
        elif isinstance(x, (list, str)):
            if isinstance(x, list) and len(x) > 0 and not isinstance(x[0], tuple):
                # todo: This is ok, but is certainly not "industry strength" fast
                # compared to what is possible.
                x = gap.eval('PermList(%s)'%x)
            return PermutationGroupElement(x, self, check = True)
        elif isinstance(x, tuple):
            return PermutationGroupElement([x], self, check = True)
        elif isinstance(x, (int, long, Integer)) and x == 1:
            return self.identity()
        else:
            raise TypeError, "unable to coerce %s to permutation in %s"%(x, self)

    def _coerce_impl(self, x):
        r"""
        EXAMPLES:
        The following test verifies that Trac #270 has been fixed:

        sage: g1 = PermutationGroupElement([(1,2),(3,4,5)])
        sage: g1.parent()
        Symmetric group of order 5! as a permutation group
        sage: g2 = PermutationGroupElement([(1,2)])
        sage: g2.parent()
        Symmetric group of order 2! as a permutation group
        sage: g1*g2
        (3,4,5)
        sage: g2*g2
        ()
        sage: g2*g1
        (3,4,5)
        """
        if isinstance(x, PermutationGroupElement):
            if x.parent() is self:
                return x
            elif x.parent().degree() <= self.degree() and x._gap_() in self._gap_():
                return PermutationGroupElement(x._gap_(), self, check = False)
        raise TypeError

    def list(self):
        """
        Return list of all elements of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3,4)], [(1,2)]])
            sage: G.list()
            [(), (3,4), (2,3), (2,3,4), (2,4,3), (2,4), (1,2), (1,2)(3,4), (1,2,3), (1,2,3,4), (1,2,4,3), (1,2,4), (1,3,2), (1,3,4,2), (1,3), (1,3,4), (1,3)(2,4), (1,3,2,4), (1,4,3,2), (1,4,2), (1,4,3), (1,4), (1,4,2,3), (1,4)(2,3)]
        """
        X = self._gap_().Elements()
        n = X.Length()
        return [PermutationGroupElement(X[i], self, check = False)
                            for i in range(1,n+1)]

    def __contains__(self, item):
        """
        Returns boolean value of "item in self"

        EXAMPLES:
            sage: G = SymmetricGroup(16)
	    sage: g = G.gen(0)
	    sage: h = G.gen(1)
	    sage: g^7*h*g*h in G
	    True
	    sage: G = SymmetricGroup(4)
	    sage: g = G((1,2,3,4))
	    sage: h = G((1,2))
            sage: H = PermutationGroup([[(1,2,3,4)], [(1,2),(3,4)]])
	    sage: g in H
	    True
	    sage: h in H
	    False
        """
        if isinstance(item, PermutationGroupElement):
            if not item.parent() is self:
                try:
		    item = PermutationGroupElement(item._gap_(), self, check = True)
		except:
		    return False
 	    return True
        elif isinstance(item, (list, str)):
            if isinstance(item, list) and len(item) > 0 and not isinstance(item[0], tuple):
                item = gap.eval('PermList(%s)'%item)
            try:
	        item = PermutationGroupElement(item, self, check = True)
	    except:
	        return False
            return True
        elif isinstance(item, tuple):
            if isinstance(item, list) and len(item) > 0 and not isinstance(item[0], tuple):
                item = gap.eval('PermList(%s)'%item)
            try:
	        item = PermutationGroupElement(item, self, check = True)
	    except:
	        return False
            return True
        elif isinstance(item, (int, long, Integer)) and item == 1:
            return True
        return False

    def has_element(self,item):
        """
        Returns boolean value of "item in self" -- however *ignores* parentage.

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
	    sage: gens = G.gens()
	    sage: H = DihedralGroup(4)
	    sage: g = G([(1,2,3,4)]); g
            (1,2,3,4)
	    sage: G.has_element(g)
            True
	    sage: h = H([(1,2),(3,4)]); h
            (1,2)(3,4)
	    sage: G.has_element(h)
            False

        """
        L = [str(x) for x in self.list()]
        return (str(item) in L)

    def __iter__(self):
        """
        Return an iterator over the elements of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3)], [(1,2)]])
            sage: [a for a in G]
            [(), (2,3), (1,2), (1,2,3), (1,3,2), (1,3)]
        """
        for g in self.list():
            yield g

    def gens(self):
        """
        Return tuple of generators of this group.  These need not be
        minimal, as they are the generators used in defining this
        group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3)], [(1,2)]])
            sage: G.gens()
            ((1,2,3), (1,2))

        Note that the generators need not be minimal.
            sage: G = PermutationGroup([[(1,2)], [(1,2)]])
            sage: G.gens()
            ((1,2), (1,2))

            sage: G = PermutationGroup([[(1,2,3,4), (5,6)], [(1,2)]])
            sage: g = G.gens()
            sage: g[0]
            (1,2,3,4)(5,6)
            sage: g[1]
            (1,2)
        """
        try:
            return self.__gens
        except AttributeError:
            try:
                gens = self._gap_().GeneratorsOfGroup()
            except TypeError, s:
                raise RuntimeError, "(It might be necessary to install the database_gap optional SAGE package, if you haven't already.)\n%s"%s
            self.__gens = tuple([PermutationGroupElement(gens[n],
                                    self, check = False) for n in \
                                 range(1, int(gens.Length())+1)])
            return self.__gens

    def gen(self, i):
        return self.gens()[int(i)]

    def identity(self):
        """
        Return the identity element of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)]])
            sage: e = G.identity()
            sage: e
            ()
            sage: g = G.gen(0)
            sage: g*e
            (1,2,3)(4,5)
            sage: e*g
            (1,2,3)(4,5)
        """
        return PermutationGroupElement('()', self, check=True)

    def degree(self):
        """
        Synonym for largest_moved_point().

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G.degree()
            5
        """
        return self.largest_moved_point()

    def largest_moved_point(self):
        """
        Return the largest point moved by a permutation in this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
            sage: G.largest_moved_point()
            4
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4,10)]])
            sage: G.largest_moved_point()
            10
        """
        try:
            return self.__largest_moved_point
        except AttributeError:
            n = int(self._gap_().LargestMovedPoint())
        self.__largest_moved_point = n
        return n

    def smallest_moved_point(self):
        """
        Return the smallest point moved by a permutation in this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(3,4)], [(2,3,4)]])
            sage: G.smallest_moved_point()
            2
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4,10)]])
            sage: G.smallest_moved_point()
            1
        """
        return self._gap_().SmallestMovedPoint()

    def _repr_(self):
        G = self.gens()
        return "Permutation Group with generators %s"%list(self.gens())

    def _latex_(self):
        return '\\langle ' + \
               ', '.join([x._latex_() for x in self.gens()]) + ' \\rangle'

    def order(self):
        """
        Return the number of elements of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)], [(1,2)]])
            sage: G.order()
            12
        """
        return Integer(self._gap_().Size())

    def random_element(self):
        """
        Return a random element of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)], [(1,2)]])
            sage: G.random_element()         ## random output
            (1, 3)(4, 5)
        """
        return PermutationGroupElement(self._gap_().Random(),
                                       self, check=False)

    def group_id(self):
        """
        Return the ID code of this group, which is a list of two
        integers. Requires "optional" database_gap-4.4.x package.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)], [(1,2)]])
            sage: G.group_id()    # requires optional database_gap-4.4.6 package
            [12, 4]
        """
        return [Integer(n) for n in eval(str(self._gap_().IdGroup()))]

    def id(self):
        """
        Same as self.group_id()
        """
        return self.group_id()

    def direct_product(self,other,maps=True):
        """
        Wraps GAP's DirectProduct, Embedding, and Projection.

        SAGE calls GAP's DirectProduct, which chooses an efficient representation for the direct product.
        The direct product of permutation groups will be a permutation group again.
        For a direct product D, the GAP operation Embedding(D,i) returns the homomorphism embedding the
        i-th factor into D. The GAP operation Projection(D,i) gives the projection of D onto the
        i-th factor.

        INPUT:
            self, other -- permutation groups

        This method returns a 5-tuple - a permutation groups and 4 morphisms.

        OUTPUT:
            D     -- a direct product of the inputs, returned as a permutation group as well
            iota1 -- an embedding of self into D
            iota2 -- an embedding of other into D
            pr1   -- the projection of D onto self  (giving a splitting 1 -> other -> D ->> self -> 1)
            pr2   -- the projection of D onto other (giving a splitting 1 -> self -> D ->> other -> 1)

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
            sage: D = G.direct_product(G,False)
            sage: D
            Permutation Group with generators [(1,2,3,4), (5,6,7,8)]
            sage: D,iota1,iota2,pr1,pr2 = G.direct_product(G)
            sage: D; iota1; iota2; pr1; pr2
            Permutation Group with generators [(1,2,3,4), (5,6,7,8)]
            Homomorphism : Cyclic group of order 4 as a permutation group --> Permutation Group with generators [(1,2,3,4), (5,6,7,8)]
            Homomorphism : Cyclic group of order 4 as a permutation group --> Permutation Group with generators [(1,2,3,4), (5,6,7,8)]
            Homomorphism : Permutation Group with generators [(1,2,3,4), (5,6,7,8)] --> Cyclic group of order 4 as a permutation group
            Homomorphism : Permutation Group with generators [(1,2,3,4), (5,6,7,8)] --> Cyclic group of order 4 as a permutation group

            sage: g=D([(1,3),(2,4)]); g
            (1,3)(2,4)
            sage: d=D([(1,4,3,2),(5,7),(6,8)]); d
            (1,4,3,2)(5,7)(6,8)
            sage: iota1(g); iota2(g); pr1(d); pr2(d)
            (1,3)(2,4)
            (5,7)(6,8)
            (1,4,3,2)
            (1,3)(2,4)

        """
        from sage.groups.perm_gps.permgroup_morphism import PermutationGroupMorphism_from_gap
        G1 = self._gap_init_()
        #print G1
        G2 = other._gap_init_()
        cmd1 = "G:=DirectProduct("+G1+","+G2+")"
        cmd2 = "iota1:=Embedding(G,1)"
        cmd3 = "iota2:=Embedding(G,2)"
        cmd4 = "pr1:=Projection(G,1)"
        cmd5 = "pr2:=Projection(G,2)"
        if not(maps):
            return PermutationGroup(gap.eval(cmd1), from_group = True)
        else:
            D = PermutationGroup_generic(gap.eval(cmd1), from_group = True)
            iota1 = PermutationGroupMorphism_from_gap(self,D, cmd2, "iota1")
            iota2 = PermutationGroupMorphism_from_gap(other,D, cmd3, "iota2")
            pr1 = PermutationGroupMorphism_from_gap(D,self, cmd4, "pr1")
            pr2 = PermutationGroupMorphism_from_gap(D,other, cmd5, "pr2")
            return D,iota1,iota2,pr1,pr2

    def center(self):
        """
        Return the subgroup of elements of that commute with
        every element of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3,4)]])
            sage: G.center()
            Permutation Group with generators [(1,2,3,4)]
            sage: G = PermutationGroup([[(1,2,3,4)], [(1,2)]])
            sage: G.center()
            Permutation Group with generators [()]
        """
        C = self._gap_().Center()
        return PermutationGroup(C, from_group = True)

    def derived_series(self):
        """
        Return the derived series of this group as a list of
        permutation groups.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G.derived_series()        # somewhat random output
            [Permutation Group with generators [(1,2,3)(4,5), (3,4)],
             Permutation Group with generators [(1,5)(3,4), (1,5)(2,4), (2,4)(3,5)]]
        """
        ans = []
        DS = self._gap_().DerivedSeries()
        n = DS.Length()
        for i in range(1,n+1):
            ans.append(PermutationGroup(DS[i], from_group = True))
        return ans

    def character_table(self):
        r"""
        Returns the matrix of values of the irreducible characters of
        a permutation group $G$ at the conjugacy classes of $G$. The
        columns represent the the conjugacy classes of $G$ and the
        rows represent the different irreducible characters in the
        ordering given by GAP.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3)]])
            sage: G.order()
            12
            sage: G.character_table()
            [         1          1          1          1]
            [         1          1 -zeta3 - 1      zeta3]
            [         1          1      zeta3 -zeta3 - 1]
            [         3         -1          0          0]
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3)]])
            sage: CT = gap(G).CharacterTable()

        Type print gap.eval("Display(%s)"%CT.name()) to display this nicely.

            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
            sage: G.order()
            8
            sage: G.character_table()
            [ 1  1  1  1  1]
            [ 1 -1 -1  1  1]
            [ 1 -1  1 -1  1]
            [ 1  1 -1 -1  1]
            [ 2  0  0  0 -2]
            sage: CT = gap(G).CharacterTable()

        Again, type print gap.eval("Display(%s)"%CT.name()) to display this.

            sage: SymmetricGroup(2).character_table()
            [ 1 -1]
            [ 1  1]
            sage: SymmetricGroup(3).character_table()
            [ 1 -1  1]
            [ 2  0 -1]
            [ 1  1  1]
            sage: SymmetricGroup(5).character_table()
            [ 1 -1  1  1 -1 -1  1]
            [ 4 -2  0  1  1  0 -1]
            [ 5 -1  1 -1 -1  1  0]
            [ 6  0 -2  0  0  0  1]
            [ 5  1  1 -1  1 -1  0]
            [ 4  2  0  1 -1  0 -1]
            [ 1  1  1  1  1  1  1]
            sage: list(AlternatingGroup(6).character_table())
            [(1, 1, 1, 1, 1, 1, 1), (5, 1, 2, -1, -1, 0, 0), (5, 1, -1, 2, -1, 0, 0), (8, 0, -1, -1, 0, zeta5^3 + zeta5^2 + 1, -zeta5^3 - zeta5^2), (8, 0, -1, -1, 0, -zeta5^3 - zeta5^2, zeta5^3 + zeta5^2 + 1), (9, 1, 0, 0, 1, -1, -1), (10, -2, 1, 1, 0, 0, 0)]

        Suppose that you have a class function $f(g)$ on $G$ and you
        know the values $v_1, ..., v_n$ on the conjugacy class
        elements in \code{conjugacy_classes_representatives(G)} =
        $[g_1, \ldots, g_n]$.  Since the irreducible characters
        $\rho_1, \ldots, \rho_n$ of $G$ form an $E$-basis of the space
        of all class functions ($E$ a ``sufficiently large''
        cyclotomic field), such a class function is a linear
        combination of these basis elements, $f = c_1\rho_1 + \cdots +
        c_n\rho_n$. To find the coefficients $c_i$, you simply solve
        the linear system \code{character_table_values(G)}*$[v_1, ...,
        v_n] = [c_1, ..., c_n]$,
        where $[v_1, ...,v_n]$ = \code{character_table_values(G)}$^{-1}[c_1, ...,c_n]$.

        AUTHORS:
            - David Joyner and William Stein (2006-01-04)
        """
        G    = self._gap_()
        cl   = G.ConjugacyClasses()
        n    = int(cl.Length())
        irrG = G.Irr()
        ct   = [[irrG[i+1,j+1] for j in range(n)] for i in range(n)]

        # Now we have the tricky task of figuring out what common
        # cyclotomic field to put all these cyclotomic elements in.
        # Our trick is to compute a list of strings that begin with
        # the order of the root of unit for each element followed by ")junk",
        # then get rid of the junk.
        s = str(ct).split('E(')[1:]   # with junk
        s = [eval(x.split(')')[0]) for x in s]  # get rid of trailing junk

        from sage.rings.all import lcm, CyclotomicField
        e = lcm(s)

        # Now make field and coerce all elements into it.
        K = CyclotomicField(e)
        ct = [[K(x) for x in v] for v in ct]

        # Finally return the result as a matrix.
        from sage.matrix.all import MatrixSpace
        MS = MatrixSpace(K,n)
        return MS(ct)


    def is_abelian(self):
        """
        Return True if this group is abelian.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: G.is_abelian()
            False
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G.is_abelian()
            True
        """
        t = self._gap_().IsAbelian()
        return t.bool()

    def is_commutative(self):
        """
        Return True if this group is commutative.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: G.is_commutative()
            False
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G.is_commutative()
            True
        """
        return self.is_abelian()

    def is_isomorphic(self,right,mode=None):
        """
        Return True if the groups are isomorphic. If mode="verbose" then
        an isomorphism is printed.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: H = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G.is_isomorphic(H)
            False

        Now, if you type G.is_isomorphic(G,mode="verbose"), of course you will get
        "True" but in addition an isomorphism (determined by a map of group generators)
        is printed. The isomorphism is random, but in this case you will probably
        get something like [ (2,3), (1,2)(3,4,5) ] -> [ (2,3), (1,2)(3,5,4) ],
        which indicates that g1 = (2,3) and g2 = (1,2)(3,4,5) are generators of G and
        the map found is the identity map.
        """
        G = self._gap_init_()
        H = right._gap_init_()
        gap.eval("x:=IsomorphismGroups( %s, %s )"%(G,H))
        if mode=="verbose":
            print gap.eval("x")
            return not(gap.eval("x")=="fail")
        return not(gap.eval("x")=="fail")

    def is_normal(self,other):
        """
        Return True if the group self is a normal subgroup of other.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: H = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: G.is_normal(H)
            True

        """
        t = self._gap_().IsNormal(other._gap_())
        return t.bool()

    def is_solvable(self):
        """
        Returns True if the group is solvable.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G.is_solvable()
            true

        """
        ans = self._gap_().IsSolvableGroup()
        return ans

    def is_subgroup(self,other):
        """
        Returns true if self is a subgroup of other.

        EXAMPLES:
            sage: G = AlternatingGroup(5)
            sage: H = SymmetricGroup(5)
            sage: G.is_subgroup(H)
            True
        """
        G = other
        gens = self.gens()
        for i in range(len(gens)):
            x = gens[i]
            if not (G.has_element(x)):
                return False
        return True

    def conjugacy_classes_representatives(self):
        """
        Returns a complete list of representatives of conjugacy classes
        in a permutation group G. The ordering is that given by GAP.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
            sage: cl = G.conjugacy_classes_representatives(); cl
            [(), (2,4), (1,2)(3,4), (1,2,3,4), (1,3)(2,4)]
            sage: cl[3] in G
            True

            sage: G = SymmetricGroup(5)
            sage: G.conjugacy_classes_representatives ()
            [(), (1,2), (1,2)(3,4), (1,2,3), (1,2,3)(4,5), (1,2,3,4), (1,2,3,4,5)]


        AUTHOR: David Joyner and William Stein (2006-01-04)
        """
        cl = self._gap_().ConjugacyClasses()
        n = int(cl.Length())
        L = gap("List([1..Length(%s)], i->Representative(%s[i]))"%(
            cl.name(),  cl.name()))
        return [PermutationGroupElement(L[i], self, check=False) \
                for i in range(1,n+1)]

    def conjugacy_classes_subgroups(self):
        """
        Returns a complete list of representatives of conjugacy classes
        of subgroups in a permutation group G. The ordering is that given by GAP.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
            sage: cl = G.conjugacy_classes_subgroups()
            sage: cl[3]
            Group([ (2,4) ])

            sage: G = SymmetricGroup(3)
            sage: G.conjugacy_classes_subgroups()
            [Group(()), Group([ (2,3) ]), Group([ (1,2,3) ]), Group([ (1,3,2), (1,2) ])]

        AUTHOR: David Joyner (2006-10)
        """
        G = self._gap_()
        cl = G.ConjugacyClassesSubgroups()
        n = int(cl.Length())
        L = gap("List([1..Length(%s)], i->Representative(%s[i]))"%(
            cl.name(),  cl.name()))
        return [PermutationGroupElement(L[i], self, check=False) \
                for i in range(1,n+1)]

    def normalizer(self,g):
        """
        Returns the normalizer of g in self.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3,4)]])
            sage: g = G.random_element()
            sage: G.normalizer(g)  ## random output
            Group([ (2,4), (1,3) ])
            sage: g = G([(1,2,3,4)])
            sage: G.normalizer(g)
            Group([ (1,2,3,4), (1,3)(2,4), (2,4) ])

        """
        N = self._gap_().Normalizer(str(g))
        return N

def direct_product_permgroups(P):
    """
    Takes the direct product of the permutation groups listed in P.

    EXAMPLES:
        sage: G1 = AlternatingGroup([1,2,4,5])
        sage: G2 = AlternatingGroup([3,4,6,7])
        sage: D = direct_product_permgroups([G1,G2,G1])
        sage: D.order()
        1728
        sage: D = direct_product_permgroups([G1])
        sage: D==G1
        True
        sage: direct_product_permgroups([])
        Symmetric group of order 0! as a permutation group
    """
    n = len(P)
    if n == 0:
        return SymmetricGroup(1)
    if n == 1:
        return P[0]
    from sage.groups.perm_gps.permgroup_morphism import PermutationGroupMorphism_from_gap
    G = [H._gap_init_() for H in P]
    Glist = ""
    for H in G:
        Glist = Glist + H + ","
    cmd = "G:=DirectProduct([" + Glist[:-1] + "])"
    gap.eval(cmd)
    return PermutationGroup(gap.eval("G"), from_group = True)

class SymmetricGroup(PermutationGroup_generic):
    """
    The full symmetric group of order $n!$, as a permutation group.
    (If n is a list of positive integers then it returns the
    symmetric group of the associated set.)

    INPUT:
       n -- a positive integer

    EXAMPLE:
        sage: G = SymmetricGroup(8)
        sage: G.order()
        40320
        sage: G
        Symmetric group of order 8! as a permutation group
        sage: G.degree()
        8
        sage: S8 = SymmetricGroup(8)
        sage: loads(dumps(S8)) == S8
        True
        sage: G = SymmetricGroup([1,2,4,5])
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.set()
        [1, 2, 4, 5]
        sage: G = SymmetricGroup(4)
        sage: G
        Symmetric group of order 4! as a permutation group
        sage: G.set()
        [1, 2, 3, 4]
    """
    def __init__(self, n):
        self._deg = n
        if isinstance(n, ListType):
            PermutationGroup_generic.__init__(self, 'SymmetricGroup(%s)'%n, from_group = True)
        else:
            try:
                n = Integer(n)
                if n < 1:
                    raise ValueError, "n (=%s) must be >= 1"%n
                PermutationGroup_generic.__init__(self, 'SymmetricGroup(%s)'%n, from_group = True)
            except TypeError, msg:
                raise ValueError, "%s\nn (=%s) must be an integer >= 1 or a list (but n has type %s)"%(msg, n,type(n))

    def _repr_(self):
        if isinstance(self._deg,ListType):
            deg = len(self._deg)
        else:
            deg = self.degree()
        return "Symmetric group of order %s! as a permutation group"%deg

    def __str__(self):
        if isinstance(self._deg,ListType):
            deg = len(self._deg)
        else:
            deg = self.degree()
        return "SymmetricGroup(%s)"%deg

    def set(self):
        if isinstance(self._deg,ListType):
            X = self._deg
        else:
            X = range(1,self._deg + 1)
        return X

class AlternatingGroup(PermutationGroup_generic):
    """
    The alternating group of order $n!/2$, as a permutation group.

    INPUT:
        n -- integer $n \geq 1$

    EXAMPLE:
        sage: G = AlternatingGroup(8)
        sage: G.order()
        20160
        sage: G
        Alternating group of order 8!/2 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: G = AlternatingGroup([1,2,4,5])
        sage: G
        Alternating group of order 4!/2 as a permutation group
        sage: G.set()
        [1, 2, 4, 5]
    """
    def __init__(self, n):
        self._deg = n
        if isinstance(n,ListType):
            PermutationGroup_generic.__init__(self, 'AlternatingGroup(%s)'%n, from_group = True)
        elif isinstance(n,Integer):
            n = Integer(n)
            if n < 1:
                raise ValueError, "n (=%s) must be >= 1"%n
            PermutationGroup_generic.__init__(self, 'AlternatingGroup(%s)'%n, from_group = True)
        else:
            raise ValueError, "n (=%s) must be an integer >= 1 or a list"%n

    def _repr_(self):
        if isinstance(self._deg,ListType):
            deg = len(self._deg)
        else:
            deg = self.degree()
        return "Alternating group of order %s!/2 as a permutation group"%deg

    def __str__(self):
        if isinstance(self._deg,ListType):
            deg = len(self._deg)
        else:
            deg = self.degree()
        return "AlternatingGroup(%s)"%deg

    def set(self):
        if isinstance(self._deg,ListType):
            X = self._deg
        else:
            X = range(1,self._deg + 1)
        return X

class CyclicPermutationGroup(PermutationGroup_generic):
    """
    A cyclic group of order n, as a permutation group.

    INPUT:
        n -- a positive integer

    EXAMPLE:
        sage: G = CyclicPermutationGroup(8)
        sage: G.order()
        8
        sage: G
        Cyclic group of order 8 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: C = CyclicPermutationGroup(10)
        sage: C.is_abelian()
        True
        sage: C = CyclicPermutationGroup(10)
        sage: C.as_AbelianGroup()
        Multiplicative Abelian Group isomorphic to C2 x C5

    """
    def __init__(self, n):
        """
        """
        n = Integer(n)
        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n
        gens = tuple(range(1, n+1))
        PermutationGroup_generic.__init__(self, [gens], n)

    def _repr_(self):
        return "Cyclic group of order %s as a permutation group"%self.order()

    def is_commutative(self):
        return True

    def is_abelian(self):
        return True

    def as_AbelianGroup(self):
	"""
	Returns the corresponding Abelian Group instance.

	EXAMPLES:

	"""
	n = self.order()
        a = list(factor(n))
        invs = [x[0]**x[1] for x in a]
        G = AbelianGroup(len(a),invs)
	return G

class KleinFourGroup(PermutationGroup_generic):
    r"""
    The Klein 4 Group, which has order $4$ and exponent $2$, viewed
    as a subgroup of $S_4$.

    OUTPUT:
        -- the Klein 4 group of order 4, as a permutation group of degree 4.

    EXAMPLES:
        sage: G = KleinFourGroup(); G
        The Klein 4 group of order 4, as a permutation group
        sage: list(G)
        [(), (3,4), (1,2), (1,2)(3,4)]

    AUTHOR:
        -- Bobby Moretti (2006-10)
    """
    def __init__(self):
        gens = ((1,2),(3,4))
        PermutationGroup_generic.__init__(self, gens, from_group=True)

    def _repr_(self):
        return 'The Klein 4 group of order 4, as a permutation group'


class DihedralGroup(PermutationGroup_generic):
    """
    The Dihedral group of order $2n$ for any integer $n\geq 1$.

    INPUT:
        n -- a positive integer

    OUTPUT:
        -- the dihedral group of order 2*n, as a permutation group

    EXAMPLE:
        sage: DihedralGroup(1)
        Dihedral group of order 2 as a permutation group

        sage: DihedralGroup(2)
        Dihedral group of order 4 as a permutation group
        sage: DihedralGroup(2).gens()
        ((1,2), (3,4))

        sage: DihedralGroup(5).gens()
        ((1,2,3,4,5), (1,5)(2,4))
        sage: list(DihedralGroup(5))
        [(), (2,5)(3,4), (1,2)(3,5), (1,2,3,4,5), (1,3)(4,5), (1,3,5,2,4), (1,4)(2,3), (1,4,2,5,3), (1,5,4,3,2), (1,5)(2,4)]

        sage: G = DihedralGroup(6)
        sage: G.order()
        12
        sage: G = DihedralGroup(5)
        sage: G.order()
        10
        sage: G
        Dihedral group of order 10 as a permutation group
        sage: loads(G.dumps()) == G
        True
        sage: G.gens()
        ((1,2,3,4,5), (1,5)(2,4))

        sage: DihedralGroup(0)
        Traceback (most recent call last):
        ...
        ValueError: n must be positive
    """
    def __init__(self, n):
        n = Integer(n)
        if n <= 0:
            raise ValueError, "n must be positive"

        # the first generator generates the cyclic subgroup of D_n, <(1...n)> in
        # cycle notation
        gen0 = range(1,n+1)

        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n

        # D_1 is a subgroup of S_2, we need the cyclic group of order 2
        if n == 1:
            gens = CyclicPermutationGroup(2).gens()
        elif n == 2:
            gens = ((1,2),(3,4))
        else:
            #gens = tuple(gen0.append([(i, n+1-i) for i in range(1, n//2 + 1)]))
            gen1 = [(i, n-i+1) for i in range(1, n//2 +1)]
            gens = tuple([tuple(gen0),tuple(gen1)])
        # send this off to the parent's class __init__()
        PermutationGroup_generic.__init__(self, gens, from_group = True)

    def _repr_(self):
        return "Dihedral group of order %s as a permutation group"%self.order()

class MathieuGroup(PermutationGroup_generic):
    """
    The Mathieu group of degree $n$.

    INPUT:
        n -- a positive integer in  {9, 10, 11, 12, 21, 22, 23, 24}.

    OUTPUT:
        -- the Mathieu group of degree n, as a permutation group

    EXAMPLE:
        sage: G = MathieuGroup(12)
        sage: G
        Mathieu group of degree 12 and order 95040 as a permutation group
    """
    def __init__(self, n):
        n = Integer(n)
        self._n = n
        if not(n in [9, 10, 11, 12, 21, 22, 23, 24]):
            raise ValueError,"argument must belong to {9, 10, 11, 12, 21, 22, 23, 24}."
        id = 'MathieuGroup(%s)'%n
        PermutationGroup_generic.__init__(self, id, from_group=True, check=False)

    def _repr_(self):
        return "Mathieu group of degree %s and order %s as a permutation group"%(self._n,self.order())

class TransitiveGroup(PermutationGroup_generic):
    """
    The transitive group from the GAP tables of transitive groups.
    """
    def __init__(self, d, n):
        """
        INPUT:
            d -- positive integer; the degree
            n -- positive integer; the number

        OUTPUT:
            the n-th transitive group of degree d

        EXAMPLE:
            sage: G = TransitiveGroup(1,1); G
            Transitive group number 1 of degree 1
            sage: G = TransitiveGroup(5, 2); G         # requires optional database_gap
            Transitive group number 2 of degree 5
            sage: G.gens()                             # requires optional database_gap
            ((1,2,3,4,5), (1,4)(2,3))

            sage: loads(G.dumps()) == G                # requires optional database_gap
            True
        """
        if d == 1:
            id = 'Group([()])'
        else:
            id = 'TransitiveGroup(%s,%s)'%(d,n)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._d = d
        self._n = n

    def _repr_(self):
        return "Transitive group number %s of degree %s"%(self._n, self._d)

class PGL(PermutationGroup_generic):
    """
    The projective general linear groups over GF(q).
    """
    def __init__(self, n, q):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field

        OUTPUT:
            PGL(n,q)

        EXAMPLE:
            sage: G = PGL(2,3); G
            Permutation Group with generators [(3,4), (1,2,4)]
            sage: print G
            The projective general linear group of degree 2 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3
            sage: G.order()
            24

        """
        if n == 1:
            id = 'Group([()])'
        else:
            id = 'PGL(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)
        self._n = n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective general linear group of degree %s over %s"%(self._n, self.base_ring())

class PSL(PermutationGroup_generic):
    """
    The projective special linear groups over GF(q).
    """
    def __init__(self, n, q):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field

        OUTPUT:
            PSL(n,q)

        EXAMPLE:
            sage: G = PSL(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G.base_ring()
            Finite Field of size 3
            sage: print G
            The projective special linear group of degree 2 over Finite Field of size 3

        """
        if n == 1:
            id = 'Group([()])'
        else:
            id = 'PSL(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)
        self._n = n

    def matrix_degree(self):
        return self._n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective special linear group of degree %s over %s"%(self._n, self.base_ring())

    def ramification_module_decomposition_hurwitz_curve(self):
        """
        Helps compute the decomposition of the ramification module
        for the Hurwitz curves X (over CC say) with automorphism group
        G = PSL(2,q), q a "Hurwitz prime" (ie, p is $\pm 1 \pmod 7$).
        Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        Hurwitz curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE: David Joyner, Amy Ksir, Roger Vogeler,
                   "Group representations on Riemann-Roch spaces of some
                   Hurwitz curves," preprint, 2006.

        EXAMPLES:
            sage: G = PSL(2,13)
            sage: G.ramification_module_decomposition_hurwitz_curve() #random
            [0, 7, 7, 12, 12, 12, 13, 15, 14]

        This means, for example, that the trivial representation does not
        occur in the ramification module of a Hurwitz curve with automorphism
        group PSL(2,13), since the trivial representation is listed first
        and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.

        If you try to use this function on a group PSL(2,q) where q is
        not a (smallish) "Hurwitz prime", an error mesage will be printed.
        """
        if self.matrix_degree()!=2:
            return ValueError, "Degree must be 2."
        F = self.base_ring()
        q = F.order()
        from sage.misc.misc import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/joyner/hurwitz_crv_rr_sp.gap'
        gap.eval('Read("'+gapcode+'")')
        mults = gap.eval("ram_module_hurwitz("+str(q)+")")
        return eval(mults)

    def ramification_module_decomposition_modular_curve(self):
        """
        Helps compute the decomposition of the ramification module
        for the modular curve X(p) (over CC say) with automorphism group G = PSL(2,q),
        q a prime > 5. Using this computation and Borne's formula helps determine the
        G-module structure of the RR spaces of equivariant
        divisors can be determined explicitly.

        The output is a list of integer multiplicities: [m1,...,mn],
        where n is the number of conj classes of G=PSL(2,p) and mi is the
        multiplicity of pi_i in the ramification module of a
        modular curve with automorphism group G.
        Here IrrRepns(G) = [pi_1,...,pi_n] (in the order listed in the
        output of self.character_table()).

        REFERENCE: D. Joyner and A. Ksir, 'Modular representations
                   on some Riemann-Roch spaces of modular curves
                   $X(N)$', Computational Aspects of Algebraic Curves,
                   (Editor: T. Shaska) Lecture Notes in Computing, WorldScientific,
                   2005.)

        EXAMPLES:
            sage: G = PSL(2,7)
            sage: G.ramification_module_decomposition_modular_curve() ## random
            [0, 4, 3, 6, 7, 8]

        This means, for example, that the trivial representation does not
        occur in the ramification module of X(7), since the trivial representation
        is listed first and that entry has multiplicity 0. The "randomness" is due to the
        fact that GAP randomly orders the conjugacy classes of the same order
        in the list of all conjugacy classes. Similarly, there is some
        randomness to the ordering of the characters.
        """
        if self.matrix_degree()!=2:
            return ValueError, "Degree must be 2."
        F = self.base_ring()
        q = F.order()
        from sage.misc.misc import SAGE_EXTCODE
        gapcode = SAGE_EXTCODE + '/gap/joyner/modular_crv_rr_sp.gap'
        gap.eval('Read("'+gapcode+'")')
        mults = gap.eval("ram_module_X("+str(q)+")")
        return eval(mults)

class PSp(PermutationGroup_generic):
    """
    The projective symplectic linear groups over GF(q).
    """
    def __init__(self, n, q):
        """
        INPUT:
            n -- positive integer; the degree
            q -- prime power; the size of the ground field

        OUTPUT:
            PSp(n,q)

        EXAMPLE:
            sage: G = PSp(2,3); G
            Permutation Group with generators [(2,3,4), (1,2)(3,4)]
            sage: G.order()
            12
            sage: G = PSp(4,3); G
            Permutation Group with generators [(3,4)(6,7)(9,10)(12,13)(17,20)(18,21)(19,22)(23,32)(24,33)(25,34)(26,38)(27,
            39)(28,40)(29,35)(30,36)(31,37), (1,5,14,17,27,22,19,36,3)(2,6,32)(4,7,23,20,37,13,16,26,40)(8,24,29,30,39,10,
            33,11,34)(9,15,35)(12,25,38)(21,28,31)]
            sage: G.order()
            25920
            sage: print G
            The projective symplectic linear group of degree 4 over Finite Field of size 3
            sage: G.base_ring()
            Finite Field of size 3

        """
        if n%2 == 1:
            raise TypeError, "The degree n must be even"
        else:
            id = 'PSp(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)
        self._n = n

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The projective symplectic linear group of degree %s over %s"%(self._n, self.base_ring())

PSP = PSp

class PSU(PermutationGroup_generic):
    """
    The projective special unitary groups over GF(q).

    INPUT:
        n -- positive integer; the degree
        q -- prime power; the size of the ground field

    OUTPUT:
        PSU(n,q)

    EXAMPLE:
        sage: PSU(2,3)
        The projective special unitary group of degree 2 over Finite Field of size 3
    """
    def __init__(self, n, q, var='a'):
        id = 'PSU(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)
        self._field_of_definition = GF(q**2, var)
        self._n = n

    def field_of_definition(self):
        return self._field_of_definition

    def base_ring(self):
        return self._base_ring

    def _repr_(self):
        return "The projective special unitary group of degree %s over %s"%(self._n, self.base_ring())

class PGU(PermutationGroup_generic):
    """
    The projective general unitary groups over GF(q).

    INPUT:
        n -- positive integer; the degree
        q -- prime power; the size of the ground field

    OUTPUT:
        PGU(n,q)

    EXAMPLE:
        sage: PGU(2,3)
        The projective general unitary group of degree 2 over Finite Field of size 3
    """
    def __init__(self, n, q, var='a'):
        id = 'PGU(%s,%s)'%(n,q)
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)
        self._field_of_definition = GF(q**2, var)
        self._n = n

    def field_of_definition(self):
        return self._field_of_definition

    def base_ring(self):
        return self._base_ring

    def _repr_(self):
        return "The projective general unitary group of degree %s over %s"%(self._n, self.base_ring())


class Suzuki(PermutationGroup_generic):
    r"""
    The Suzuki group over GF(q), $^2 B_2(2^{2k+1}) = Sz(2^{2k+1})$. A wrapper for the GAP function SuzukiGroup.

    INPUT:
        $q = 2^n$ -- an odd power of 2; the size of the ground field. (Strictly speaking, n should be
	             greater than 1, or else this group os not simple.)
    OUTPUT:
        Sz(q)

    EXAMPLE:
        sage: Suzuki(8)
        Permutation Group with generators [(1,28,10,44)(3,50,11,42)(4,43,53,64)(5,9,39,52)(6,36,63,13)(7,51,60,57)(8,33,
        37,16)(12,24,55,29)(14,30,48,47)(15,19,61,54)(17,59,22,62)(18,23,34,31)(20,38,
        49,25)(21,26,45,58)(27,32,41,65)(35,46,40,56), (1,2)(3,10)(4,42)(5,18)(6,50)(7,26)(8,58)(9,34)(12,28)(13,45)(14,44)(15,
        23)(16,31)(17,21)(19,39)(20,38)(22,25)(24,61)(27,60)(29,65)(30,55)(32,33)(35,
        52)(36,49)(37,59)(40,54)(41,62)(43,53)(46,48)(47,56)(51,63)(57,64)]
        sage: print Suzuki(8)
        The Suzuki group over Finite Field in a of size 2^3


    REFERENCES:
        http://en.wikipedia.org/wiki/Group_of_Lie_type\#Suzuki-Ree_groups
    """
    def __init__(self, q, var='a'):
	if is_even(round(log(q,2))):
	    raise ValueError,"The ground field size %s must be an odd power of 2."%q
        id = 'SuzukiGroup(IsPermGroup,%s)'%q
        PermutationGroup_generic.__init__(self, id,
                                          from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q, var)

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The Suzuki group over %s"%self.base_ring()

"""
class Ree(PermutationGroup_generic):
    #
    #The Ree group over GF(q), $^2 G_2(3^{2k+1}) = Ree(3^{2k+1})$. A wrapper for the GAP function ReeGroup.
    #
    #INPUT:
    #    q = 3^n -- an odd power of 3; the size of the ground field. (Strictly speaking, n should be
    #	             greater than 1, or else this group os not simple.)
    #
    #OUTPUT:
    #    Ree(q)
    #
    #EXAMPLE:
    #    sage: Ree(27)
    #    ???
    #    sage: print Ree(27)
    #    ???
    #
    #REFERENCES:
    #    http://en.wikipedia.org/wiki/Group_of_Lie_type#Suzuki-Ree_groups
    #
    def __init__(self, q):
	if is_even(int(log(q,3))):
	    raise ValueError,"The ground field size %s must be an odd power of 3."%q
	id = 'ReeGroup(IsPermGroup,%s)'%q   ######## It appears IsPermGroup
	                                    ######## has not been implemented in GAP yet!
        PermutationGroup_generic.__init__(self, id, from_group=True, check=False)
        self._q = q
        self._base_ring = GF(q)

    def base_ring(self):
        return self._base_ring

    def __str__(self):
        return "The Ree group over %s"%self.base_ring()
"""

class PermutationGroup_subgroup(PermutationGroup_generic):
    """
    Subgroup subclass of PermutationGroup_generic, so instance methods are
    inherited.

    EXAMPLES:
        sage: G = CyclicPermutationGroup(4)
        sage: gens = G.gens()
        sage: H = DihedralGroup(4)
        sage: PermutationGroup_subgroup(H,list(gens))
        Subgroup of Dihedral group of order 8 as a permutation group generated by [(1,2,3,4)]
        sage: K=PermutationGroup_subgroup(H,list(gens))
        sage: K.list()
        [(), (1,2,3,4), (1,3)(2,4), (1,4,3,2)]
        sage: K.ambient_group()
        Dihedral group of order 8 as a permutation group
        sage: K.gens()
        [(1,2,3,4)]
    """
    def __init__(self, ambient, gens, from_group = False,
                 check=True):
        if not isinstance(ambient, PermutationGroup_generic):
            raise TypeError, "ambient (=%s) must be perm group."%ambient
        if not isinstance(gens, list):
            raise TypeError, "gens (=%s) must be a list"%gens
        if check:
            pass

        self.__ambient_group = ambient
        self.__gens = gens
        cmd = 'Group(%s)'%gens
        cmd = cmd.replace("'","")  # get rid of quotes
        self.__gap = cmd

        G = ambient
        if check:
            for i in range(len(gens)):
                x = gens[i]
                if not (G.has_element(x)):
                    raise TypeError, "each generator must be in the ambient group"
        self.__ambient_group = G

        PermutationGroup_generic.__init__(self, gens, from_group, check)

    def __cmp__(self, other):
        r"""
        Compare self and other.  If self and other are in a common ambient group,
        then self <= other precisely if self is contained in other.

        EXAMPLES:
            sage: G = CyclicPermutationGroup(4)
	    sage: gens = G.gens()
	    sage: H = DihedralGroup(4)
 	    sage: PermutationGroup_subgroup(H,list(gens))
            Subgroup of Dihedral group of order 8 as a permutation group generated by [(1,2,3,4)]
	    sage: K=PermutationGroup_subgroup(H,list(gens))
            sage: G<K
            False
            sage: G>K
            False
        """
        if self is other:
            return 0
        if not isinstance(other, PermutationGroup_generic):
            return -1
        c = cmp(self.ambient_group(), other.ambient_group())
        if c: return c
        if self.is_subgroup(other):
            return -1
        else:
            return 1

    def _repr_(self):
        s = "Subgroup of %s generated by %s"%(self.ambient_group(), self.gens())
        return s

    def _latex_(self):
        r"""
        Return latex representation of this group.
        """
        return self._repr_()


    def ambient_group(self):
        """
        Return the ambient group related to self.
        """
        return self.__ambient_group

    def gens(self):
        """
        Return the generators for this subgroup.
        """
        return self.__gens

    def is_normal(self):
	"""
	Returns True if self is normal in ambient_group.

       EXAMPLES:

       """
        amb = self.ambient_group()
        U = self._gap_init_()
        G = amb._gap_init_()
        return gap.eval("IsNormal("+G+","+U)

    def is_subgroup(self,ambient):
        """
        Returns true if self is a subgroup of ambient.

        EXAMPLES:

        """
        G = ambient
        gens = self.gens()
        for i in range(len(gens)):
            x = gens[i]
            if not (G.has_element(x)):
                return 0
        return 1



