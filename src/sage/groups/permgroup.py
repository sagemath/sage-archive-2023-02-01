r"""
Permutation groups

In \sage a permutation is represented as either a string
that defines a permutation using disjoint cycle notation,
or a list of tuples, which represent disjoint cycles.

\begin{verbatim}
(a,...,b)(c,...,d)...(e,...,f)  <--> [(a,...,b), (c,...,d),..., (e,...,f)]
                  () = identity <--> []
\end{verbatim}

JOKE:
    Q: What's hot, chunky, and acts on a polygon? A: Dihedral soup.
    Renteln, P. and Dundes, A. "Foolproof: A Sampling of Mathematical
    Folk Humor." Notices Amer. Math. Soc. 52, 24--34, 2005.

AUTHOR:
    - David Joyner (2005-10-14): first version
    - David Joyner (2005-11-17)
    - William Stein (2005-11-26): rewrite to better wrap Gap
    - David Joyner (2005-12-21)
    - Stein and Joyner (2005-01-04): added conjugacy_class_representatives

"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import random

import sage.structure.element as element
import group

from sage.rings.all      import RationalField, Integer, MPolynomial, MPolynomialRing, Polynomial
from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement

import sage.ext.coerce as coerce

def gap_format(x):
    """
    Put a permutation in Gap format, as a string.
    """
    x = str(x).replace(' ','')
    return x.replace('),(',')(').replace('[','').replace(']','')

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
        sage: G.random()
        (2,5)

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

    def __cmp__(self, right):
        """
        Compare self and right.

        The ordering is whatever it is in Gap.

        EXAMPLES:
            sage: G1 = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G2 = PermutationGroup([[(1,2,3),(4,5)]])
            sage: G1 < G2
            True
            sage: G1 > G2
            False
        """
        if not isinstance(right, PermutationGroup_generic):
            return -1
        return self._gap_().__cmp__(right._gap_())


    def __call__(self, x):
        """
        Coerce x into this permutation group.

        EXAMPLES:
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
        if isinstance(x, (int, long, Integer)) and x == 1:
            return self.identity()
        if isinstance(x, PermutationGroupElement):
            if x.parent() is self:
                return x
            else:
                return PermutationGroupElement(x._gap_(), self, check = True)
        elif isinstance(x, (list, str)):
            return PermutationGroupElement(x, self, check = True)
        elif isinstance(x, tuple):
            return PermutationGroupElement([x], self, check = True)
        else:
            raise TypeError, "unable to coerce %s to permutation in %s"%(x, self)


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
        return self._gap_().LargestMovedPoint()

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

    def random(self):
        """
        Return a random element of this group.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)], [(1,2)]])
            sage: G.random()
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
            [Permutation Group with generators [(1,2,3)(4,5), (3,4)], Permutation Group with generators [(1,5)(3,4), (1,5)(2,4), (2,4)(3,5)]]
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
            [          1           1           1           1]
            [          1           1 -zeta_3 - 1      zeta_3]
            [          1           1      zeta_3 -zeta_3 - 1]
            [          3          -1           0           0]
            sage: G = PermutationGroup([[(1,2),(3,4)], [(1,2,3)]])
            sage: CT = gap(G).CharacterTable()
            sage: print gap.eval("Display(%s)"%CT.name())
            CT1
            <BLANKLINE>
                2  2  2  .  .
                3  1  .  1  1
            <BLANKLINE>
                   1a 2a 3a 3b
                2P 1a 1a 3b 3a
                3P 1a 2a 1a 1a
            <BLANKLINE>
            X.1     1  1  1  1
            X.2     1  1  A /A
            X.3     1  1 /A  A
            X.4     3 -1  .  .
            <BLANKLINE>
            A = E(3)^2
              = (-1-ER(-3))/2 = -1-b3

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
            sage: print gap.eval("Display(%s)"%CT.name())
            CT2
            <BLANKLINE>
             2  3  2  2  2  3
            <BLANKLINE>
               1a 2a 2b 4a 2c
            2P 1a 1a 1a 2c 1a
            3P 1a 2a 2b 4a 2c
            <BLANKLINE>
            X.1     1  1  1  1  1
            X.2     1 -1 -1  1  1
            X.3     1 -1  1 -1  1
            X.4     1  1 -1 -1  1
            X.5     2  .  .  . -2

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
            [(1, 1, 1, 1, 1, 1, 1),
             (5, 1, -1, 2, -1, 0, 0),
             (5, 1, 2, -1, -1, 0, 0),
             (8, 0, -1, -1, 0, zeta_5^3 + zeta_5^2 + 1, -zeta_5^3 - zeta_5^2),
             (8, 0, -1, -1, 0, -zeta_5^3 - zeta_5^2, zeta_5^3 + zeta_5^2 + 1),
             (9, 1, 0, 0, 1, -1, -1),
             (10, -2, 1, 1, 0, 0, 0)]

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

class PermutationGroupElement(element.Element_cmp_,
                              element.MultiplicativeGroupElement):
    """
    An element of a permutation group.

    EXAMPLES:
        sage: G = PermutationGroup(['(1,2,3)(4,5)'])
        sage: G
        Permutation Group with generators [(1,2,3)(4,5)]
        sage: g = G.gen(0); g
        (1,2,3)(4,5)
        sage: print g
        (1,2,3)(4,5)
        sage: g*g
        (1,3,2)
        sage: g**(-1)
        (1,3,2)(4,5)
        sage: g**2
        (1,3,2)
        sage: G = PermutationGroup([(1,2,3)])
        sage: g = G.gen(0); g
        (1,2,3)
        sage: g.order()
        3

    This example illustrates how permutations act on multivariate
    polynomials.

        sage: R = MPolynomialRing(RationalField(), 5, ["x","y","z","u","v"])
        sage: x, y, z, u, v = R.gens()
        sage: f = x**2 - y**2 + 3*z**2
        sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
        sage: sigma = G.gen(0)
        sage: f * sigma
        -1*z^2 + y^2 + 3*x^2

    """
    def __init__(self, g, parent = None, check = True):
        r"""
        Create element of a permutation group.

        There are several ways to define a permutation group element:
        \begin{itemize}
        \item  Define a permutation group $G$, then use
               \code{G.gens()} and multiplication * to construct
               elements.
        \item Define a permutation group $G$, then use e.g.,
               \code{G([(1,2),(3,4,5)])} to construct an element of
               the group.  You could also use \code{G('(1,2)(3,4,5)')}
        \item Use e.g., \code{PermutationGroupElement([(1,2),(3,4,5)])}
        or \code{PermutationGroupElement('(1,2)(3,4,5)')}
              to make a permutation group element with parent $S_5$.
        \end{itemize}

        INPUT:
            g -- defines element
            parent (optional) -- defines parent group (g must be in parent if
                                 specified, or a TypeError is raised).
            check -- bool (default: True), if False assumes g is a
                     gap element in parent (if specified).

        EXAMPLES:
        We illustrate construction of permutation using several
        different methods.

        First we construct elements by multiplying together generators
        for a group.

            sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'])
            sage: s = G.gens()
            sage: s[0]
            (1,2)(3,4)
            sage: s[1]
            (3,4,5,6)
            sage: s[0]*s[1]
            (1,2)(3,5,6)
            sage: (s[0]*s[1]).parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]

        Next we illustrate creation of a permutation using
        coercion into an already-created group.

            sage: g = G([(1,2),(3,5,6)])
            sage: g
            (1,2)(3,5,6)
            sage: g.parent()
            Permutation Group with generators [(1,2)(3,4), (3,4,5,6)]
            sage: g == s[0]*s[1]
            True

        We can also use a string instead of a list to specify
        the permutation.

            sage: h = G('(1,2)(3,5,6)')
            sage: g == h
            True

        We can also make a permutation group element directly
        using the \code{PermutationGroupElement} command.  Note
        that the parent is then the full symmetric group $S_n$,
        where $n$ is the largest integer that is moved by the
        permutation.

            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k
            (1,2)(3,5,6)
            sage: k.parent()
            Symmetric group of order 6! as a permutation group

        Note the comparison of permutations doesn't require that the
        parent groups are the same.

            sage: k == g
            True

        Arithmetic with permutations having different parents is also defined:

            sage: k*g
            (3,6,5)
            sage: (k*g).parent()
            Symmetric group of order 6! as a permutation group

            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: loads(dumps(G.0)) == G.0
            True

        EXAMPLES:
            sage: k = PermutationGroupElement('(1,2)(3,5,6)')
            sage: k._gap_()
            (1,2)(3,5,6)
            sage: k._gap_().parent()
            Gap
        """
        if check:
            if not (parent is None or isinstance(parent, PermutationGroup_generic)):
                raise TypeError, 'parent must be a permutation group'
            self.__gap = gap_format(g)
            if not parent is None:
                P = parent._gap_()
                if not self._gap_(P.parent()) in P:
                    raise TypeError, 'permutation %s not in %s'%(self.__gap, parent)
        else:
            self.__gap = str(g)
        if parent is None:
            parent = SymmetricGroup(self._gap_().LargestMovedPoint())
        element.Element.__init__(self, parent)

    def _gap_init_(self):
        return self.__gap


    def _repr_(self):
        """
        Return string representation of this permutation.

        EXAMPLES:
        We create the permutation $(1,2,3)(4,5)$ and print it.

            sage: g = PermutationGroupElement([(1,2,3),(4,5)])
            sage: g._repr_()
            '(1,2,3)(4,5)'

        Permutation group elements support renaming them so
        they print however you want, as illustred below:

            sage: g.rename('sigma')
            sage: g
            sigma
            sage: g.rename()
            sage: g
            (1,2,3)(4,5)
        """
        return self.__gap

    def _latex_(self):
        return str(self)

    def __getitem__(self, i):
        """
        Return the ith permutation cycle in the disjoint cycle
        representation of self.

        INPUT:
            i -- integer

        OUTPUT:
            a permutation group element

        EXAMPLE:
            sage: G = PermutationGroup([[(1,2,3),(4,5)]],5)
            sage: g = G.gen(0)
            sage: g[0]
            (1,2,3)
            sage: g[1]
            (4,5)
        """
        S = str(self).split(')(')
        if i >= 0 and i < len(S):
            T = S[i]
            if i > 0:
                T = '(' + T
            if i < len(S)-1:
                T += ')'
        else:
            raise IndexError, "i (=%s) must be between 0 and %s, inclusive"%(i, len(S)-1)
        return PermutationGroupElement(gap(T), check = False)

    def _cmp_(self, right):
        """
        Compare group elements self and right.

        EXAMPLES:
            sage: G = PermutationGroup([[(1,2,3),(4,5)],[(3,4)]])
            sage: G.gen(0) < G.gen(1)
            False
            sage: G.gen(0) > G.gen(1)
            True
        """
        r = right._gap_()
        G = r.parent()
        l = self._gap_(G)
        return cmp(l,r)

    def __call__(self, i):
        """
        Returns the image of the integer i under this permutation.

        EXAMPLE:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: G
            Permutation Group with generators [(1,2,3)(4,5)]
            sage: g = G.gen(0)
            sage: g(5)
            4
        """
        return int(gap.eval('%s^%s'%(i, self._gap_().name())))

    #def __mul__(self, other):
    #    """
    #    This overloaded operator implements multiplication *and*
    #    permutation action on polynomials.
    #    EXAMPLES:
    #        sage: G = PermutationGroup(['(1,2)(3,4)', '(3,4,5,6)'])
    #        sage: g = G.gens()
    #        sage: g[0] * g[1]
    #        (1,2)(3,5,6)
    #    """
    #    if isinstance(other, MPolynomial):
    #        return self.right_action_on_polynomial(other)
    #    else:
    #        return element.MultiplicativeGroupElement.__mul__(self, other)

    def _r_action(self, left):
        """
        Return the right action of self on left.

        For example, if f=left is a polynomial, then this function
        returns f(sigma*x), which is image of f under the right action
        of sigma on the indeterminates.  This is a right action since
        the image of f(sigma*x) under tau is f(sigma*tau*x).

        INPUT:
            left -- element of space on which permutations act from the right

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)', '(1,2,3,4,5)'])
            sage: R = MPolynomialRing(RationalField(), 5, ["x","y","z","u","v"])
            sage: x,y,z,u,v = R.gens()
            sage: f = x**2 + y**2 - z**2 + 2*u**2
            sage: sigma, tau = G.gens()
            sage: f*sigma
            2*v^2 + z^2 + y^2 - x^2
            sage: f*tau
            2*v^2 - u^2 + z^2 + y^2
            sage: f*(sigma*tau)
            u^2 + z^2 - y^2 + 2*x^2
            sage: (f*sigma)*tau
            u^2 + z^2 - y^2 + 2*x^2
        """
        if isinstance(left, Polynomial):
            if not (self == 1):
                raise ValueError, "%s does not act on %s"%(self, left.parent())
            return left
        elif isinstance(left, PermutationGroupElement):
            return PermutationGroupElement(self._gap_()*left._gap_(),
                                           parent = None, check = True)
        elif isinstance(left, MPolynomial):
            F = left.base_ring()
            R = left.parent()
            x = R.gens()
            vars = list(x)
            try:
                sigma_x  = [vars[int(self(i+1)-1)] for i in range(len(x))]
            except IndexError:
                raise ValueError, "%s does not act on %s"%(self, left.parent())
            return left(tuple(sigma_x))
        else:
            raise TypeError, "left (=%s) must be a polynomial."%left


    def _mul_(self, other):
        return PermutationGroupElement(self._gap_()*other._gap_(),
                                       self.parent(), check = True)

    def _div_(self, other):
        """
        Returns self divided by other, i.e., self times the inverse
        of other.

        EXAMPLES:
            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: h = PermutationGroupElement('(1,2,3)')
            sage: g/h
            (4,5)
        """
        return PermutationGroupElement(self._gap_()/other._gap_(),
                                           self.parent(),
                                           check = True)

    def __invert__(self):
        """
        Return the inverse of this permutation.

        EXAMPLES:
            sage: g = PermutationGroupElement('(1,2,3)(4,5)')
            sage: ~g
            (1,3,2)(4,5)
            sage: (~g) * g
            ()
        """
        return PermutationGroupElement(self._gap_().Inverse(),
                          self.parent(), check=False)

    def order(self):
        """
        Return the order of this group element, which is the smallest
        positive integer $n$ for which $g^n = 1$.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.order()
            6
        """
        return int(self._gap_().Order())

    def sign(self):
        """
        Returns the sign of self, which is $(-1)^{s}$, where $s$ is
        the number of swaps.

        EXAMPLES:
            sage: s = PermutationGroupElement('(1,2)(3,5,6)')
            sage: s.sign()
            -1
        """
        return int(self._gap_().SignPerm())


    def orbit(self, n):
        """
        Returns the orbit of the integer $n$ under this group element,
        as a sorted list of integers.

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.orbit(4)
            [4, 5]
            sage: g.orbit(3)
            [1, 2, 3]
            sage: g.orbit(10)
            [10]
        """
        n = Integer(n)
        # We use eval to avoid creating intermediate gap objects (so
        # this is slightly faster.
        s = gap.eval('Orbit(Group(%s), %s)'%(self._gap_().name(), Integer(n)))
        v = [Integer(k) for k in eval(s)]
        v.sort()
        return v

    def matrix(self):
        """
        Returns deg x deg permutation matrix associated
        to the permutation self

        EXAMPLES:
            sage: G = PermutationGroup(['(1,2,3)(4,5)'])
            sage: g = G.gen(0)
            sage: g.matrix()
            [0 1 0 0 0]
            [0 0 1 0 0]
            [1 0 0 0 0]
            [0 0 0 0 1]
            [0 0 0 1 0]
        """
        deg = self.parent().degree()
        M = MatrixSpace(RationalField(),deg,deg)
        A = M(0)
        for i in range(deg):
            A[i, self(i+1) - 1] = 1
        return A

class SymmetricGroup(PermutationGroup_generic):
    """
    The full symmetric group of order $n!$, as a permutation group.
    """
    def __init__(self, n):
        """
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
        """
        n = Integer(n)
        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n
        PermutationGroup_generic.__init__(self, 'SymmetricGroup(%s)'%n, from_group = True)

    def _repr_(self):
        return "Symmetric group of order %s! as a permutation group"%self.degree()


class AlternatingGroup(PermutationGroup_generic):
    """
    The alternating group of order $n!/2$, as a permutation group.
    """
    def __init__(self, n):
        """
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
        """
        n = Integer(n)
        if n < 1:
            raise ValueError, "n (=%s) must be >= 1"%n
        PermutationGroup_generic.__init__(self, 'AlternatingGroup(%s)'%n, from_group = True)

    def _repr_(self):
        return "Alternating group of order %s!/2 as a permutation group"%self.degree()


class CyclicPermutationGroup(PermutationGroup_generic):
    """
    A cyclic group of order n, as a permutation group.

    EXAMPLES:
        sage: C = CyclicPermutationGroup(10)
        sage: C.is_abelian()
        True
    """
    def __init__(self, n):
        """
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

class DihedralGroup(PermutationGroup_generic):
    """
    The Dihedral group of degree $n$ and order $2n$.
    """
    def __init__(self, n):
        """
        INPUT:
            n -- a positive integer

        OUTPUT:
            -- the dihedral group of order 2*n, as a permutation group

        EXAMPLE:
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
        """
        n = Integer(n)
        n0 = n//2
        m = 2*n0 + 1
        if n % 2 != 0:
            m += 1
        gen0 = tuple(range(1,m))
        gen1 = tuple([(i,m-i) for i in range(1,n0+1)])
        PermutationGroup_generic.__init__(self, [gen0, gen1], from_group = True)

    def _repr_(self):
        return "Dihedral group of order %s as a permutation group"%self.order()


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

