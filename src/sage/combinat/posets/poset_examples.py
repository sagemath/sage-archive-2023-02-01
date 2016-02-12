r"""
A catalog of posets and lattices.

Some common posets can be accessed through the ``posets.<tab>`` object::

    sage: posets.PentagonPoset()
    Finite lattice containing 5 elements

Moreover, the set of all posets of order `n` is represented by ``Posets(n)``::

    sage: Posets(5)
    Posets containing 5 vertices

**Catalog of common posets:**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Posets.AntichainPoset` | Return an antichain on `n` elements.
    :meth:`~Posets.BooleanLattice` | Return the Boolean lattice on `2^n` elements.
    :meth:`~Posets.ChainPoset` | Return a chain on `n` elements.
    :meth:`~Posets.DiamondPoset` | Return the lattice of rank two on `n` elements.
    :meth:`~Posets.IntegerCompositions` | Return the poset of integer compositions of `n`.
    :meth:`~Posets.IntegerPartitions` | Return the poset of integer partitions of ``n``.
    :meth:`~Posets.PentagonPoset` | Return the Pentagon poset.
    :meth:`~Posets.RandomPoset` | Return a random poset on `n` vertices according to a probability `p`.
    :meth:`~Posets.RestrictedIntegerPartitions` | Return the poset of integer partitions of `n`, ordered by restricted refinement.
    :meth:`~Posets.ShardPoset` | Return the shard intersection order.
    :meth:`~Posets.SSTPoset` | Return the poset on semistandard tableaux of shape `s` and largest entry `f` that is ordered by componentwise comparison.
    :meth:`~Posets.StandardExample` | Return the standard example of a poset with dimension `n`.
    :meth:`~Posets.SymmetricGroupBruhatIntervalPoset` | The poset of permutations with respect to Bruhat order.
    :meth:`~Posets.SymmetricGroupBruhatOrderPoset` | The poset of permutations with respect to Bruhat order.
    :meth:`~Posets.SymmetricGroupWeakOrderPoset` | The poset of permutations of `\{ 1, 2, \ldots, n \}` with respect to the weak order.
    :meth:`~Posets.TamariLattice` | Return the Tamari lattice.
    :meth:`~Posets.TetrahedralPoset` | Return the Tetrahedral poset with `n-1` layers based on the input colors.

Constructions
-------------
"""
#*****************************************************************************
#       Copyright (C) 2008 Peter Jipsen <jipsen@chapman.edu>,
#                          Franco Saliola <saliola@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.classcall_metaclass import ClasscallMetaclass
import sage.categories.posets
from sage.combinat.permutation import Permutations, Permutation
from sage.combinat.posets.posets import Poset, FinitePosets_n
from sage.combinat.posets.lattices import LatticePoset
from sage.graphs.digraph import DiGraph
from sage.rings.integer import Integer

class Posets(object):
    r"""
    A collection of posets and lattices.

    EXAMPLES::

        sage: Posets.BooleanLattice(3)
        Finite lattice containing 8 elements
        sage: Posets.ChainPoset(3)
        Finite lattice containing 3 elements
        sage: Posets.RandomPoset(17,.15)
        Finite poset containing 17 elements

    The category of all posets::

        sage: Posets()
        Category of posets

    The enumerated set of all posets on `3` vertices, up to an
    isomorphism::

        sage: Posets(3)
        Posets containing 3 vertices

    .. seealso:: :class:`~sage.categories.posets.Posets`, :class:`FinitePosets`, :func:`Poset`

    TESTS::

        sage: P = Posets
        sage: TestSuite(P).run()
    """

    __metaclass__ = ClasscallMetaclass
    @staticmethod
    def __classcall__(cls, n = None):
        r"""
        Return either the category of all posets, or the finite
        enumerated set of all finite posets on ``n`` elements up to an
        isomorphism.

        EXAMPLES::

            sage: Posets()
            Category of posets
            sage: Posets(4)
            Posets containing 4 vertices
        """
        if n is None:
            return sage.categories.posets.Posets()
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        return FinitePosets_n(n)

    @staticmethod
    def BooleanLattice(n):
        """
        Returns the Boolean lattice containing `2^n` elements.

        EXAMPLES::

            sage: Posets.BooleanLattice(5)
            Finite lattice containing 32 elements
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        if n==0:
            return LatticePoset( ([0], []) )
        if n==1:
            return LatticePoset( ([0,1], [[0,1]]) )
        return LatticePoset([[Integer(x|(1<<y)) for y in range(0,n) if x&(1<<y)==0] for
            x in range(0,2**n)])

    @staticmethod
    def ChainPoset(n):
        """
        Returns a chain (a totally ordered poset) containing ``n`` elements.

        EXAMPLES::

            sage: C = Posets.ChainPoset(6); C
            Finite lattice containing 6 elements
            sage: C.linear_extension()
            [0, 1, 2, 3, 4, 5]
            sage: for i in range(5):
            ...       for j in range(5):
            ...           if C.covers(C(i),C(j)) and j != i+1:
            ...              print "TEST FAILED"

        TESTS:

        Check that :trac:`8422` is solved::

            sage: Posets.ChainPoset(0)
            Finite lattice containing 0 elements
            sage: C = Posets.ChainPoset(1); C
            Finite lattice containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = Posets.ChainPoset(2); C
            Finite lattice containing 2 elements
            sage: C.cover_relations()
            [[0, 1]]
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        return LatticePoset((range(n), [[x,x+1] for x in range(n-1)]))

    @staticmethod
    def AntichainPoset(n):
        """
        Returns an antichain (a poset with no comparable elements)
        containing `n` elements.

        EXAMPLES::

            sage: A = Posets.AntichainPoset(6); A
            Finite poset containing 6 elements
            sage: for i in range(5):
            ...       for j in range(5):
            ...           if A.covers(A(i),A(j)):
            ...              print "TEST FAILED"

        TESTS:

        Check that :trac:`8422` is solved::

            sage: Posets.AntichainPoset(0)
            Finite poset containing 0 elements
            sage: C = Posets.AntichainPoset(1); C
            Finite poset containing 1 elements
            sage: C.cover_relations()
            []
            sage: C = Posets.AntichainPoset(2); C
            Finite poset containing 2 elements
            sage: C.cover_relations()
            []
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        return Poset((range(n), []))

    @staticmethod
    def PentagonPoset(facade = None):
        """
        Returns the Pentagon poset.

        INPUT:

        - ``facade`` (boolean) -- whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`). The
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor).

        EXAMPLES::

            sage: P = Posets.PentagonPoset(); P
            Finite lattice containing 5 elements
            sage: P.cover_relations()
            [[0, 1], [0, 2], [1, 4], [2, 3], [3, 4]]

        This is smallest lattice that is not modular::

            sage: P.is_modular()
            False

        This poset and the :meth:`DiamondPoset` are the two smallest
        lattices which are not distributive::

            sage: P.is_distributive()
            False
            sage: Posets.DiamondPoset(5).is_distributive()
            False
        """
        p = LatticePoset([[1,2],[4],[3],[4],[]], facade = facade)
        p.hasse_diagram()._pos = {0:[2,0],1:[0,2],2:[3,1],3:[3,3],4:[2,4]}
        return p

    @staticmethod
    def DiamondPoset(n, facade = None):
        """
        Return the lattice of rank two containing ``n`` elements.

        INPUT:

        - ``n`` - number of vertices, an integer at least 3.

        - ``facade`` (boolean) -- whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`). The
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor).

        EXAMPLES::

            sage: Posets.DiamondPoset(7)
            Finite lattice containing 7 elements
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n <= 2:
            raise ValueError("n must be an integer at least 3")
        c = [[n-1] for x in range(n)]
        c[0] = [x for x in range(1,n-1)]
        c[n-1] = []
        return LatticePoset(c, facade = facade)

    @staticmethod
    def IntegerCompositions(n):
        """
        Returns the poset of integer compositions of the integer ``n``.

        A composition of a positive integer `n` is a list of positive
        integers that sum to `n`. The order is reverse refinement:
        `[p_1,p_2,...,p_l] < [q_1,q_2,...,q_m]` if `q` consists
        of an integer composition of `p_1`, followed by an integer
        composition of `p_2`, and so on.

        EXAMPLES::

            sage: P = Posets.IntegerCompositions(7); P
            Finite poset containing 64 elements
            sage: len(P.cover_relations())
            192
        """
        from sage.combinat.composition import Compositions
        C = Compositions(n)
        return Poset((C, [[c,d] for c in C for d in C if d.is_finer(c)]), cover_relations=False)

    @staticmethod
    def IntegerPartitions(n):
        """
        Returns the poset of integer partitions on the integer ``n``.

        A partition of a positive integer `n` is a non-increasing list
        of positive integers that sum to `n`. If `p` and `q` are
        integer partitions of `n`, then `p` covers `q` if and only
        if `q` is obtained from `p` by joining two parts of `p`
        (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.IntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            28
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers
            of elements in the poset of integer partitions.
            """
            lc = []
            for i in range(0,len(partition)-1):
                for j in range(i+1,len(partition)):
                    new_partition = partition[:]
                    del new_partition[j]
                    del new_partition[i]
                    new_partition.append(partition[i]+partition[j])
                    new_partition.sort(reverse=True)
                    tup = tuple(new_partition)
                    if tup not in lc:
                        lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def RestrictedIntegerPartitions(n):
        """
        Returns the poset of integer partitions on the integer `n`
        ordered by restricted refinement. That is, if `p` and `q`
        are integer partitions of `n`, then `p` covers `q` if and
        only if `q` is obtained from `p` by joining two distinct
        parts of `p` (and sorting, if necessary).

        EXAMPLES::

            sage: P = Posets.RestrictedIntegerPartitions(7); P
            Finite poset containing 15 elements
            sage: len(P.cover_relations())
            17
        """
        def lower_covers(partition):
            r"""
            Nested function for computing the lower covers of elements in the
            restricted poset of integer partitions.
            """
            lc = []
            for i in range(0,len(partition)-1):
                for j in range(i+1,len(partition)):
                    if partition[i] != partition[j]:
                        new_partition = partition[:]
                        del new_partition[j]
                        del new_partition[i]
                        new_partition.append(partition[i]+partition[j])
                        new_partition.sort(reverse=True)
                        tup = tuple(new_partition)
                        if tup not in lc:
                            lc.append(tup)
            return lc
        from sage.combinat.partition import Partitions
        H = DiGraph(dict([[tuple(p),lower_covers(p)] for p in Partitions(n)]))
        return Poset(H.reverse())

    @staticmethod
    def RandomPoset(n,p):
        r"""
        Generate a random poset on ``n`` vertices according to a
        probability ``p``.

        INPUT:

        - ``n`` - number of vertices, a non-negative integer

        - ``p`` - a probability, a real number between 0 and 1 (inclusive)

        OUTPUT:

        A poset on ``n`` vertices.  The construction decides to make an
        ordered pair of vertices comparable in the poset with probability
        ``p``, however a pair is not made comparable if it would violate
        the defining properties of a poset, such as transitivity.

        So in practice, once the probability exceeds a small number the
        generated posets may be very similar to a chain.  So to create
        interesting examples, keep the probability small, perhaps on the
        order of `1/n`.

        EXAMPLES::

            sage: Posets.RandomPoset(17,.15)
            Finite poset containing 17 elements

        TESTS::

            sage: Posets.RandomPoset('junk', 0.5)
            Traceback (most recent call last):
            ...
            TypeError: number of elements must be an integer, not junk

            sage: Posets.RandomPoset(-6, 0.5)
            Traceback (most recent call last):
            ...
            ValueError: number of elements must be non-negative, not -6

            sage: Posets.RandomPoset(6, 'garbage')
            Traceback (most recent call last):
            ...
            TypeError: probability must be a real number, not garbage

            sage: Posets.RandomPoset(6, -0.5)
            Traceback (most recent call last):
            ...
            ValueError: probability must be between 0 and 1, not -0.5
        """
        from sage.misc.prandom import random
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("number of elements must be an integer, not {0}".format(n))
        if n < 0:
            raise ValueError("number of elements must be non-negative, not {0}".format(n))
        try:
            p = float(p)
        except Exception:
            raise TypeError("probability must be a real number, not {0}".format(p))
        if p < 0 or p> 1:
            raise ValueError("probability must be between 0 and 1, not {0}".format(p))

        D = DiGraph(loops=False,multiedges=False)
        D.add_vertices(range(n))
        for i in range(n):
            for j in range(n):
                if random() < p:
                    D.add_edge(i,j)
                    if not D.is_directed_acyclic():
                        D.delete_edge(i,j)
        return Poset(D,cover_relations=False)

    @staticmethod
    def SSTPoset(s,f=None):
        """
        The poset on semistandard tableaux of shape ``s`` and largest
        entry ``f`` that is ordered by componentwise comparison of the
        entries.

        INPUT:

        - ``s`` - shape of the tableaux

        - ``f`` - maximum fill number.  This is an optional
          argument.  If no maximal number is given, it will use
          the number of cells in the shape.

        NOTE: This is a basic implementation and most certainly
        not the most efficient.

        EXAMPLES::

            sage: Posets.SSTPoset([2,1])
            Finite poset containing 8 elements

            sage: Posets.SSTPoset([2,1],4)
            Finite poset containing 20 elements

            sage: Posets.SSTPoset([2,1],2).cover_relations()
            [[[[1, 1], [2]], [[1, 2], [2]]]]

            sage: Posets.SSTPoset([3,2]).bottom()  # long time (6s on sage.math, 2012)
            [[1, 1, 1], [2, 2]]

            sage: Posets.SSTPoset([3,2],4).maximal_elements()
            [[[3, 3, 4], [4, 4]]]
        """
        from sage.combinat.tableau import SemistandardTableaux
        def tableaux_is_less_than(a,b):
            atstring = []
            btstring = []
            for i in a:
                atstring += i
            for i in b:
                btstring += i
            for i in range(len(atstring)):
                if atstring[i] > btstring[i]:
                    return False
            return True
        if f is None:
            f=0
            for i in s:
                f += i
        E = SemistandardTableaux(s, max_entry=f)
        return Poset((E, tableaux_is_less_than))

    @staticmethod
    def StandardExample(n, facade=None):
        r"""
        Return the partially ordered set on ``2n`` elements with
        dimension ``n``.

        Let `P` be the poset on `\{0, 1, 2, \ldots, 2n-1\}` whose defining
        relations are that `i < j` for every `0 \leq i < n \leq j < 2n`
        except when `i + n = j`. The poset `P` is the so-called
        *standard example* of a poset with dimension `n`.

        INPUT:

        - ``n`` -- an integer `\ge 2`, dimension of the constructed poset
        - ``facade`` (boolean) -- whether to make the returned poset a
          facade poset (see :mod:`sage.categories.facade_sets`); the
          default behaviour is the same as the default behaviour of
          the :func:`~sage.combinat.posets.posets.Poset` constructor

        OUTPUT:

        The standard example of a poset of dimension `n`.

        EXAMPLES::

            sage: A = Posets.StandardExample(3); A
            Finite poset containing 6 elements
            sage: A.dimension()
            3

        REFERENCES:

        .. [Rosen] K. Rosen *Handbook of Discrete and Combinatorial
           Mathematics* (1999), Chapman and Hall.

        .. [Garg] V. Garg *Introduction to Lattice Theory with Computer
           Science Applications* (2015), Wiley.

        TESTS::

            sage: A = Posets.StandardExample(10); A
            Finite poset containing 20 elements
            sage: len(A.cover_relations())
            90

            sage: P = Posets.StandardExample(5, facade=False)
            sage: P(4) < P(3), P(4) > P(3)
            (False, False)
        """
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("dimension must be an integer, not {0}".format(n))
        if n < 2:
            raise ValueError("dimension must be at least 2, not {0}".format(n))
        return Poset((range(2*n), [[i, j+n] for i in range(n)
                                   for j in range(n) if i != j]),
                     facade=facade)

    @staticmethod
    def SymmetricGroupBruhatOrderPoset(n):
        """
        The poset of permutations with respect to Bruhat order.

        EXAMPLES::

            sage: Posets.SymmetricGroupBruhatOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10:
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        return Poset(dict([[s,s.bruhat_succ()]
                for s in Permutations(n)]),element_labels)

    @staticmethod
    def SymmetricGroupBruhatIntervalPoset(start, end):
        """
        The poset of permutations with respect to Bruhat order.

        INPUT:

        - ``start`` - list permutation

        - ``end`` - list permutation (same n, of course)

        .. note::

           Must have ``start`` <= ``end``.

        EXAMPLES:

        Any interval is rank symmetric if and only if it avoids these
        permutations::

            sage: P1 = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [3,4,1,2])
            sage: P2 = Posets.SymmetricGroupBruhatIntervalPoset([1,2,3,4], [4,2,3,1])
            sage: ranks1 = [P1.rank(v) for v in P1]
            sage: ranks2 = [P2.rank(v) for v in P2]
            sage: [ranks1.count(i) for i in uniq(ranks1)]
            [1, 3, 5, 4, 1]
            sage: [ranks2.count(i) for i in uniq(ranks2)]
            [1, 3, 5, 6, 4, 1]

        """
        start = Permutation(start)
        end = Permutation(end)
        if len(start) != len(end):
            raise TypeError("Start (%s) and end (%s) must have same length."%(start, end))
        if not start.bruhat_lequal(end):
            raise TypeError("Must have start (%s) <= end (%s) in Bruhat order."%(start, end))
        unseen = [start]
        nodes = {}
        while len(unseen) > 0:
            perm = unseen.pop(0)
            nodes[perm] = [succ_perm for succ_perm in perm.bruhat_succ()
                                if succ_perm.bruhat_lequal(end)]
            for succ_perm in nodes[perm]:
                if succ_perm not in nodes:
                    unseen.append(succ_perm)
        return Poset(nodes)

    @staticmethod
    def SymmetricGroupWeakOrderPoset(n, labels="permutations", side="right"):
        r"""
        The poset of permutations of `\{ 1, 2, \ldots, n \}` with respect
        to the weak order (also known as the permutohedron order, cf.
        :meth:`~sage.combinat.permutation.Permutation.permutohedron_lequal`).

        The optional variable ``labels`` (default: ``"permutations"``)
        determines the labelling of the elements if `n < 10`. The optional
        variable ``side`` (default: ``"right"``) determines whether the
        right or the left permutohedron order is to be used.

        EXAMPLES::

            sage: Posets.SymmetricGroupWeakOrderPoset(4)
            Finite poset containing 24 elements
        """
        if n < 10 and labels == "permutations":
            element_labels = dict([[s,"".join(map(str,s))] for s in Permutations(n)])
        if n < 10 and labels == "reduced_words":
            element_labels = dict([[s,"".join(map(str,s.reduced_word_lexmin()))] for s in Permutations(n)])
        if side == "left":
            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of left weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                    s.length() + (s.inverse().right_action_product(v)).length() == v.length()]
        else:
            def weak_covers(s):
                r"""
                Nested function for computing the covers of elements in the
                poset of right weak order for the symmetric group.
                """
                return [v for v in s.bruhat_succ() if
                    s.length() + (s.inverse().left_action_product(v)).length() == v.length()]
        return Poset(dict([[s,weak_covers(s)] for s in Permutations(n)]),element_labels)
        
    @staticmethod
    def TetrahedralPoset(n, *colors, **labels):
        r"""
        Return the tetrahedral poset based on the input colors. 
        
        This method will return the tetrahedral poset with n-1 layers and 
        covering relations based on the input colors of 'green', 'red', 
        'orange', 'silver', 'yellow' and 'blue' as defined in [Striker2011]_.  
        For particular color choices, the order ideals of the resulting 
        tetrahedral poset will be isomorphic to known combinatorial objects.
        
        For example, for the colors 'blue', 'yellow', 'orange', and 'green', 
        the order ideals will be in bijection with alternating sign matrices.
        For the colors 'yellow', 'orange', and 'green', the order ideals will 
        be in bijection with semistandard Young tableaux of staircase shape.
        For the colors 'red', 'orange', 'green', and optionally 'yellow', the 
        order ideals will be in bijection with totally symmetric 
        self-complementary plane partitions in a `2n \times 2n \times 2n` box.

        INPUT:

        - ``n`` - Defines the number (n-1) of layers in the poset.

        - ``colors`` - The colors that define the covering relations of the 
          poset. Colors used are 'green', 'red', 'yellow', 'orange', 'silver', 
          and 'blue'.
          
        - ``labels`` - Keyword variable used to determine whether the poset
          is labeled with integers or tuples.  To label with integers, the
          method should be called with ``labels='integers'``.  Otherwise, the 
          labeling will default to tuples.

        EXAMPLES::

            sage: Posets.TetrahedralPoset(4,'green','red','yellow','silver','blue','orange')
            Finite poset containing 10 elements
            
            sage: Posets.TetrahedralPoset(4,'green','red','yellow','silver','blue','orange', labels='integers')
            Finite poset containing 10 elements
            
            sage: A = AlternatingSignMatrices(3)
            sage: p = A.lattice()
            sage: ji = p.join_irreducibles_poset()
            sage: tet = Posets.TetrahedralPoset(3, 'green','yellow','blue','orange')
            sage: ji.is_isomorphic(tet)
            True
        
        REFERENCES:

        .. [Striker2011] J. Striker. *A unifying poset perpective on 
           alternating sign matrices, plane partitions, Catalan objects, 
           tournaments, and tableaux*, Advances in Applied Mathematics 46 
           (2011), no. 4, 583-609. :arXiv:`1408.5391`
        """
        n=n-1
        try:
            n = Integer(n)
        except TypeError:
            raise TypeError("n must be an integer.")
        if n < 2:
            raise ValueError("n must be greater than 2.")
        for c in colors:
            if(c not in ('green', 'red', 'yellow', 'orange', 'silver', 'blue')):
                raise ValueError("Color input must be from the following: 'green', 'red', 'yellow', 'orange', 'silver', and 'blue'.")
        elem=[(i,j,k) for i in range (n) for j in range (n-i) for k in range (n-i-j)]
        rels = []
        elem_labels = {}
        if 'labels' in labels:
            if labels['labels'] == 'integers':
                labelcount = 0
                for (i,j,k) in elem:
                    elem_labels[(i,j,k)] = labelcount
                    labelcount += 1
        for c in colors:
            for (i,j,k) in elem:
                if(i+j+k < n-1):
                    if(c=='green'):
                        rels.append([(i,j,k),(i+1,j,k)])
                    if(c=='red'):
                        rels.append([(i,j,k),(i,j,k+1)])
                    if(c=='yellow'):
                        rels.append([(i,j,k),(i,j+1,k)])
                if(j<n-1 and k>0):
                    if(c=='orange'):
                        rels.append([(i,j,k),(i,j+1,k-1)])
                if(i<n-1 and j>0):
                    if(c=='silver'):
                        rels.append([(i,j,k),(i+1,j-1,k)])
                if(i<n-1 and k>0):
                    if(c=='blue'):
                        rels.append([(i,j,k),(i+1,j,k-1)])
        return Poset([elem,rels], elem_labels)

    # shard intersection order
    import sage.combinat.shard_order
    ShardPoset = staticmethod(sage.combinat.shard_order.shard_poset)

    # Tamari lattices
    import sage.combinat.tamari_lattices
    TamariLattice = staticmethod(sage.combinat.tamari_lattices.TamariLattice)


posets = Posets
