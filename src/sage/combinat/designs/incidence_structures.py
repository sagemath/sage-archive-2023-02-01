r"""
Incidence structures (i.e. hypergraphs, i.e. set systems)

An incidence structure is specified by a list of points, blocks, or an incidence
matrix ([1]_, [2]_). :class:`IncidenceStructure` instances have the following methods:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~IncidenceStructure.ground_set` | Return the ground set (i.e the list of points).
    :meth:`~IncidenceStructure.num_points` | Return the size of the ground set.
    :meth:`~IncidenceStructure.num_blocks` | Return the number of blocks.
    :meth:`~IncidenceStructure.blocks` | Return the list of blocks.
    :meth:`~IncidenceStructure.block_sizes` | Return the set of block sizes.
    :meth:`~IncidenceStructure.degree` | Return the degree of a point ``p``
    :meth:`~IncidenceStructure.is_connected` | Test whether the design is connected.
    :meth:`~IncidenceStructure.is_simple` | Test whether this design is simple (i.e. no repeated block).
    :meth:`~IncidenceStructure.incidence_matrix` | Return the incidence matrix `A` of the design
    :meth:`~IncidenceStructure.incidence_graph` | Return the incidence graph of the design
    :meth:`~IncidenceStructure.packing` | Return a maximum packing
    :meth:`~IncidenceStructure.relabel` | Relabel the ground set
    :meth:`~IncidenceStructure.is_resolvable` | Test whether the hypergraph is resolvable
    :meth:`~IncidenceStructure.is_t_design` | Test whether ``self`` is a `t-(v,k,l)` design.
    :meth:`~IncidenceStructure.dual` | Return the dual design.
    :meth:`~IncidenceStructure.automorphism_group` | Return the automorphism group
    :meth:`~IncidenceStructure.canonical_label` | Return a canonical label for the incidence structure.
    :meth:`~IncidenceStructure.is_isomorphic` | Return whether the two incidence structures are isomorphic.
    :meth:`~IncidenceStructure.edge_coloring` | Return an optimal edge coloring`
    :meth:`~IncidenceStructure.copy` | Return a copy of the incidence structure.


REFERENCES:

.. [1] Block designs and incidence structures from wikipedia,
  :wikipedia:`Block_design`
  :wikipedia:`Incidence_structure`

.. [2] E. Assmus, J. Key, Designs and their codes, CUP, 1992.

AUTHORS:

- Peter Dobcsanyi and David Joyner (2007-2008)

  This is a significantly modified form of part of the module block_design.py
  (version 0.6) written by Peter Dobcsanyi peter@designtheory.org.

- Vincent Delecroix (2014): major rewrite

Methods
-------
"""
#***************************************************************************
#                              Copyright (C) 2007                          #
#                                                                          #
#                Peter Dobcsanyi       and         David Joyner            #
#           <peter@designtheory.org>          <wdjoyner@gmail.com>         #
#                                                                          #
#                                                                          #
#    Distributed under the terms of the GNU General Public License (GPL)   #
#    as published by the Free Software Foundation; either version 2 of     #
#    the License, or (at your option) any later version.                   #
#                    http://www.gnu.org/licenses/                          #
#***************************************************************************


from sage.misc.superseded import deprecated_function_alias
from sage.misc.cachefunc import cached_method

from sage.rings.all import ZZ
from sage.rings.integer import Integer
from sage.misc.latex import latex
from sage.sets.set import Set

def IncidenceStructureFromMatrix(M, name=None):
    """
    Deprecated function that builds an incidence structure from a matrix.

    You should now use ``designs.IncidenceStructure(incidence_matrix=M)``.

    INPUT:

    - ``M`` -- a binary matrix. Creates a set of "points" from the rows and a
      set of "blocks" from the columns.

    EXAMPLES::

        sage: BD1 = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        sage: M = BD1.incidence_matrix()
        sage: BD2 = IncidenceStructureFromMatrix(M)
        doctest:...: DeprecationWarning: IncidenceStructureFromMatrix is deprecated.
        Please use designs.IncidenceStructure(incidence_matrix=M) instead.
        See http://trac.sagemath.org/16553 for details.
        sage: BD1 == BD2
        True
    """
    from sage.misc.superseded import deprecation
    deprecation(16553, 'IncidenceStructureFromMatrix is deprecated. Please use designs.IncidenceStructure(incidence_matrix=M) instead.')
    return IncidenceStructure(incidence_matrix=M, name=name)

class IncidenceStructure(object):
    r"""
    A base class for incidence structures (i.e. hypergraphs, i.e. set systems)

    An incidence structure (i.e. hypergraph, i.e. set system) can be defined
    from a collection of blocks (i.e. sets, i.e. edges), optionally with an
    explicit ground set (i.e. point set, i.e. vertex set). Alternatively they
    can be defined from a binary incidence matrix.

    INPUT:

    - ``points`` -- (i.e. ground set, i.e. vertex set) the underlying set. If
      ``points`` is an integer `v`, then the set is considered to be `\{0, ...,
      v-1\}`.

      .. NOTE::

          The following syntax, where ``points`` is ommitted, automatically
          defines the ground set as the union of the blocks::

              sage: H = IncidenceStructure([['a','b','c'],['c','d','e']])
              sage: H.ground_set()
              ['a', 'b', 'c', 'd', 'e']

    - ``blocks`` -- (i.e. edges, i.e. sets) the blocks defining the incidence
      structure. Can be any iterable.

    - ``incidence_matrix`` -- a binary incidence matrix. Each column represents
      a set.

    - ``name`` (a string, such as "Fano plane").

    - ``check`` -- whether to check the input

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the
      :class:`IncidenceStructure` instance's internal data.

    EXAMPLES:

    An incidence structure can be constructed by giving the number of points and
    the list of blocks::

        sage: designs.IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
        Incidence structure with 7 points and 7 blocks

    Only providing the set of blocks is sufficient. In this case, the ground set
    is defined as the union of the blocks::

        sage: IncidenceStructure([[1,2,3],[2,3,4]])
        Incidence structure with 4 points and 2 blocks

    Or by its adjacency matrix (a `\{0,1\}`-matrix in which rows are indexed by
    points and columns by blocks)::

        sage: m = matrix([[0,1,0],[0,0,1],[1,0,1],[1,1,1]])
        sage: designs.IncidenceStructure(m)
        Incidence structure with 4 points and 3 blocks

    The points can be any (hashable) object::

        sage: V = [(0,'a'),(0,'b'),(1,'a'),(1,'b')]
        sage: B = [(V[0],V[1],V[2]), (V[1],V[2]), (V[0],V[2])]
        sage: I = designs.IncidenceStructure(V, B)
        sage: I.ground_set()
        [(0, 'a'), (0, 'b'), (1, 'a'), (1, 'b')]
        sage: I.blocks()
        [[(0, 'a'), (0, 'b'), (1, 'a')], [(0, 'a'), (1, 'a')], [(0, 'b'), (1, 'a')]]

    The order of the points and blocks does not matter as they are sorted on
    input (see :trac:`11333`)::

        sage: A = designs.IncidenceStructure([0,1,2], [[0],[0,2]])
        sage: B = designs.IncidenceStructure([1,0,2], [[0],[2,0]])
        sage: B == A
        True

        sage: C = designs.BlockDesign(2, [[0], [1,0]])
        sage: D = designs.BlockDesign(2, [[0,1], [0]])
        sage: C == D
        True

    If you care for speed, you can set ``copy`` to ``False``, but in that
    case, your input must be a list of lists and the ground set must be `{0,
    ..., v-1}`::

        sage: blocks = [[0,1],[2,0],[1,2]]  # a list of lists of integers
        sage: I = designs.IncidenceStructure(3, blocks, copy=False)
        sage: I.blocks(copy=False) is blocks
        True
    """
    def __init__(self, points=None, blocks=None, incidence_matrix=None,
            name=None, check=True, test=None, copy=True):
        r"""
        TESTS::

            sage: designs.IncidenceStructure(3, [[4]])
            Traceback (most recent call last):
            ...
            ValueError: Block [4] is not contained in the point set

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]], test=True)
            doctest:...: DeprecationWarning: the keyword test is deprecated,
            use check instead
            See http://trac.sagemath.org/16553 for details.
            Incidence structure with 3 points and 2 blocks

            sage: designs.IncidenceStructure(2, [[0,1,2,3,4,5]], test=False)
            Incidence structure with 2 points and 1 blocks

        We avoid to convert to integers when the points are not (but compare
        equal to integers because of coercion)::

            sage: V = GF(5)
            sage: e0,e1,e2,e3,e4 = V
            sage: [e0,e1,e2,e3,e4] == range(5)   # coercion makes them equal
            True
            sage: blocks = [[e0,e1,e2],[e0,e1],[e2,e4]]
            sage: I = designs.IncidenceStructure(V, blocks)
            sage: type(I.ground_set()[0])
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
            sage: type(I.blocks()[0][0])
            <type 'sage.rings.finite_rings.integer_mod.IntegerMod_int'>
        """
        if test is not None:
            from sage.misc.superseded import deprecation
            deprecation(16553, "the keyword test is deprecated, use check instead")
            check = test

        from sage.matrix.constructor import matrix
        from sage.structure.element import Matrix

        # Reformatting input
        if isinstance(points, Matrix):
            assert incidence_matrix is None, "'incidence_matrix' cannot be defined when 'points' is a matrix"
            assert blocks is None, "'blocks' cannot be defined when 'points' is a matrix"
            incidence_matrix = points
            points = blocks = None
        elif points and blocks is None:
            blocks = points
            points = set().union(*blocks)
        if points:
            assert incidence_matrix is None, "'incidence_matrix' cannot be defined when 'points' is defined"

        if incidence_matrix:
            M = matrix(incidence_matrix)
            v = M.nrows()
            self._points = range(v)
            self._point_to_index = None
            self._blocks = sorted(M.nonzero_positions_in_column(i) for i in range(M.ncols()))

        else:
            if isinstance(points, (int,Integer)):
                self._points = range(points)
                self._point_to_index = None
            else:
                self._points = sorted(points)
                if self._points == range(len(points)) and all(isinstance(x,(int,Integer)) for x in self._points):
                    self._point_to_index = None
                else:
                    self._point_to_index = {e:i for i,e in enumerate(self._points)}

            if check:
                for block in blocks:
                    if any(x not in self._points for x in block):
                        raise ValueError("Block {} is not contained in the point set".format(block))
                    if len(block) != len(set(block)):
                        raise ValueError("Repeated element in block {}".format(block))

            if self._point_to_index:
                # translate everything to integers between 0 and v-1
                blocks = [sorted(self._point_to_index[e] for e in block) for block in blocks]
            elif copy:
                # create a new list made of sorted blocks
                blocks = [sorted(block) for block in blocks]
            else:
                # sort the data but avoid copying it
                for b in blocks:
                    b.sort()

            blocks.sort()
            self._blocks = blocks

        self._name = str(name) if name is not None else 'IncidenceStructure'
        self._classes = None
        self._canonical_label = None

    def __iter__(self):
        """
        Iterator over the blocks.

        Note that it is faster to call for the method ``.blocks(copy=True)``
        (but in that case the output should not be modified).

        EXAMPLES::

            sage: sts = designs.steiner_triple_system(9)
            sage: list(sts)
            [[0, 1, 5], [0, 2, 4], [0, 3, 6], [0, 7, 8], [1, 2, 3], [1, 4, 7],
            [1, 6, 8], [2, 5, 8], [2, 6, 7], [3, 4, 8], [3, 5, 7], [4, 5, 6]]

            sage: b = designs.IncidenceStructure('ab', ['a','ab'])
            sage: it = iter(b)
            sage: it.next()
            ['a']
            sage: it.next()
            ['a', 'b']
        """
        if self._point_to_index is None:
            for b in self._blocks: yield b[:]
        else:
            for b in self._blocks:
                yield [self._points[i] for i in b]

    def __repr__(self):
        """
        A print method.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD
            Incidence structure with 7 points and 7 blocks
        """
        return 'Incidence structure with {} points and {} blocks'.format(
                self.num_points(), self.num_blocks())

    def __str__(self):
        """
        A print method.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD = designs.IncidenceStructure(range(7),[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: print BD
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
        """
        return '{}<points={}, blocks={}>'.format(
                self._name, self.ground_set(), self.blocks())

    def __eq__(self, other):
        """
        Tests is the two incidence structures are equal

        TESTS::

            sage: blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD1 = designs.IncidenceStructure(7, blocks)
            sage: M = BD1.incidence_matrix()
            sage: BD2 = designs.IncidenceStructure(incidence_matrix=M)
            sage: BD1 == BD2
            True

            sage: e1 = frozenset([0,1])
            sage: e2 = frozenset([2])
            sage: sorted([e1,e2]) == [e1,e2]
            True
            sage: sorted([e2,e1]) == [e2,e1]
            True
            sage: I1 = designs.IncidenceStructure([e1,e2], [[e1],[e1,e2]])
            sage: I2 = designs.IncidenceStructure([e1,e2], [[e2,e1],[e1]])
            sage: I3 = designs.IncidenceStructure([e2,e1], [[e1,e2],[e1]])
            sage: I1 == I2 and I2 == I1 and I1 == I3 and I3 == I1 and I2 == I3 and I3 == I2
            True
        """
        # We are extra careful in this method since we cannot assume that a
        # total order is defined on the point set.
        if not isinstance(other, IncidenceStructure):
            return False

        if self._points == other._points:
            return self._blocks == other._blocks

        if (self.num_points() != other.num_points() or
            self.num_blocks() != other.num_blocks()):
            return False

        p_to_i = self._point_to_index if self._point_to_index else range(self.num_points())

        if any(p not in p_to_i for p in other.ground_set()):
            return False

        other_blocks = sorted(sorted(p_to_i[p] for p in b) for b in other.blocks())
        return self._blocks == other_blocks

    def __ne__(self, other):
        r"""
        Difference test.

        EXAMPLES::

            sage: BD1 = designs.IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: M = BD1.incidence_matrix()
            sage: BD2 = designs.IncidenceStructure(incidence_matrix=M)
            sage: BD1 != BD2
            False
        """
        return not self.__eq__(other)

    def __contains__(self, block):
        r"""
        Tests if a block belongs to the incidence structure

        INPUT:

        - ``block`` -- a block.

        EXAMPLES::

            sage: [1,2,3,4] in IncidenceStructure([[1,2,3,4]])
            True
            sage: [1,2,4,3] in IncidenceStructure([[1,2,3,4]])
            True
            sage: [1,2,"3",4] in IncidenceStructure([[1,2,3,4]])
            False
            sage: [1,2,"3",4] in IncidenceStructure([[1,2,"3",4]])
            True

        More complicated examples::

            sage: str="I had a dream of a time when a 3-lines patch does not kill one hour"
            sage: sets = Subsets(str.split(), 4)
            sage: IS = IncidenceStructure(sets) # a complete 4-uniform hypergraph
            sage: ["I", "dream", "of", "one"] in IS
            True
            sage: ["does", "patch", "kill", "dream"] in IS
            True
            sage: ["Am", "I", "finally", "done ?"] in IS
            False
            sage: IS = designs.ProjectiveGeometryDesign(3, 1, GF(2))
            sage: [3,8,7] in IS
            True
            sage: [3,8,9] in IS
            False
        """
        try:
            iter(block)
        except TypeError:
            return False

        # Relabel to 0,...,n-1 if necessary
        if self._point_to_index is not None:
            try:
                block = [self._point_to_index[x] for x in block]
            except KeyError:
                return False

        return sorted(block) in self._blocks

    def canonical_label(self):
        r"""
        Return a canonical label for the incidence structure.

        A canonical label is relabeling of the points into integers
        `\{0,...,n-1\}` such that isomorphic incidence structures are
        relabelled to equal objects.

        EXAMPLE::

            sage: fano1 = designs.balanced_incomplete_block_design(7,3)
            sage: fano2 = designs.projective_plane(2)
            sage: fano1 == fano2
            False
            sage: fano1.relabel(fano1.canonical_label())
            sage: fano2.relabel(fano2.canonical_label())
            sage: fano1 == fano2
            True
        """
        if self._canonical_label is None:
            from sage.graphs.graph import Graph
            g = Graph()
            n = self.num_points()
            g.add_edges((i+n,x) for i,b in enumerate(self._blocks) for x in b)
            canonical_label = g.canonical_label([range(n),range(n,n+self.num_blocks())],certify=True)[1]
            canonical_label = [canonical_label[x] for x in range(n)]
            self._canonical_label = canonical_label

        return dict(zip(self._points,self._canonical_label))

    def is_isomorphic(self, other, certificate=False):
        r"""
        Return whether the two incidence structures are isomorphic.

        .. NOTE::

            If you need to test isomorphisms between one incidence
            structure and many others, you should consider using
            :meth:`canonical_label` instead of this function.

        INPUT:

        - ``other`` -- an incidence structure.

        - ``certificate`` (boolean) -- whether to return an
          insomorphism from ``self`` to ``other`` instead of a boolean
          answer.

        EXAMPLE::

            sage: fano1 = designs.balanced_incomplete_block_design(7,3)
            sage: fano2 = designs.projective_plane(2)
            sage: fano1.is_isomorphic(fano2)
            True
            sage: fano1.is_isomorphic(fano2,certificate=True)
            {0: 0, 1: 1, 2: 2, 3: 6, 4: 4, 5: 3, 6: 5}

        TESTS::

            sage: IS  = IncidenceStructure([["A",5,pi],["A",5,"Wouhou"],["A","Wouhou",(9,9)],[pi,12]])
            sage: IS2 = IS.copy()
            sage: IS2.relabel(IS2.canonical_label())
            sage: IS.is_isomorphic(IS2)
            True
            sage: canon = IS.is_isomorphic(IS2,certificate=True)
            sage: IS.relabel(canon)
            sage: IS==IS2
            True

            sage: IS2 = IncidenceStructure([[1,2]])
            sage: IS2.is_isomorphic(IS)
            False
            sage: IS2.is_isomorphic(IS,certificate=True)
            {}
        """
        if (self.num_points() != other.num_points() or
            self.num_blocks() != other.num_blocks() or
            sorted(self.block_sizes()) != sorted(other.block_sizes())):
            return {} if certificate else False

        A = self.copy()
        B = other.copy()

        A_canon = A.canonical_label()
        B_canon = B.canonical_label()
        A.relabel(A_canon)
        B.relabel(B_canon)

        if A == B:
            if certificate:
                B_canon_rev = {y:x for x,y in B_canon.iteritems()}
                return {x:B_canon_rev[xint] for x,xint in A_canon.iteritems()}
            else:
                return True
        else:
            return {} if certificate else False

    def copy(self):
        r"""
        Return a copy of the incidence structure.

        EXAMPLE::

            sage: IS = IncidenceStructure([[1,2,3,"e"]],name="Test")
            sage: IS
            Incidence structure with 4 points and 1 blocks
            sage: copy(IS)
            Incidence structure with 4 points and 1 blocks
            sage: [1, 2, 3, 'e'] in copy(IS)
            True
            sage: copy(IS)._name
            'Test'
        """
        IS = IncidenceStructure(self._blocks,
                                name=self._name,
                                check=False)
        IS.relabel(dict(zip(range(self.num_points()),self._points)))
        IS._canonical_label = None if self._canonical_label is None else self._canonical_label[:]

        return IS

    __copy__ = copy


    def ground_set(self, copy=True):
        r"""
        Return the ground set (i.e the list of points).

        INPUT:

        - ``copy`` (boolean) -- ``True`` by default. When set to ``False``, a
          pointer toward the object's internal data is given. Set it to
          ``False`` only if you know what you are doing.

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]]).ground_set()
            [0, 1, 2]
        """
        if copy:
            return self._points[:]
        return self._points

    def num_points(self):
        r"""
        Return the size of the ground set.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_points()
            7
            sage: B = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_points()
            4
        """
        return len(self._points)

    def num_blocks(self):
        r"""
        Return the number of blocks.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_blocks()
            7
            sage: B = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_blocks()
            5
        """
        return len(self._blocks)

    def blocks(self, copy=True):
        """
        Return the list of blocks.

        INPUT:

        - ``copy`` (boolean) -- ``True`` by default. When set to ``False``, a
          pointer toward the object's internal data is given. Set it to
          ``False`` only if you know what you are doing.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.blocks()
            [[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]

        What you should pay attention to::

            sage: blocks = BD.blocks(copy=False)
            sage: del blocks[0:6]
            sage: BD
            Incidence structure with 7 points and 1 blocks

        """
        if copy:
            if self._point_to_index is None:
                return [b[:] for b in self._blocks]
            else:
                return [[self._points[i] for i in b] for b in self._blocks]
        else:
            return self._blocks

    def block_sizes(self):
        r"""
        Return the set of block sizes.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(8, [[0,1,3],[1,4,5,6],[1,2],[5,6,7]])
            sage: BD.block_sizes()
            [3, 2, 4, 3]
            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
        """
        return map(len, self._blocks)

    def degree(self, p=None):
        r"""
        Return the degree of a point ``p``

        The degree of a point `p` is the number of blocks that contain it.

        INPUT:

        - ``p`` -- a point. If set to ``None`` (default), a dictionary
          associating the points with their degrees is returned.

        EXAMPLES::

            sage: designs.steiner_triple_system(9).degree(3)
            4
            sage: designs.steiner_triple_system(9).degree()
            {0: 4, 1: 4, 2: 4, 3: 4, 4: 4, 5: 4, 6: 4, 7: 4, 8: 4}
        """
        if p is None:
            d = [0]*self.num_points()
            for b in self._blocks:
                for x in b:
                    d[x] += 1
            return {p: d[i] for i, p in enumerate(self._points)}
        else:
            p = self._point_to_index[p] if self._point_to_index else p
            return sum(1 for b in self._blocks if p in b)

    def is_connected(self):
        r"""
        Test whether the design is connected.

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[0,2]]).is_connected()
            True
            sage: designs.IncidenceStructure(4, [[0,1],[2,3]]).is_connected()
            False
        """
        return self.incidence_graph().is_connected()

    def is_simple(self):
        r"""
        Test whether this design is simple (i.e. no repeated block).

        EXAMPLES::

            sage: designs.IncidenceStructure(3, [[0,1],[1,2],[0,2]]).is_simple()
            True
            sage: designs.IncidenceStructure(3, [[0],[0]]).is_simple()
            False

            sage: V = [(0,'a'),(0,'b'),(1,'a'),(1,'b')]
            sage: B = [[V[0],V[1]], [V[1],V[2]]]
            sage: I = designs.IncidenceStructure(V, B)
            sage: I.is_simple()
            True
            sage: I2 = designs.IncidenceStructure(V, B*2)
            sage: I2.is_simple()
            False
        """
        B = self._blocks
        return all(B[i] != B[i+1] for i in xrange(len(B)-1))

    def _gap_(self):
        """
        Return the GAP string describing the design.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD._gap_()
            'BlockDesign(7,[[1, 2, 3], [1, 4, 5], [1, 6, 7], [2, 4, 6], [2, 5, 7], [3, 4, 7], [3, 5, 6]])'
        """
        B = self.blocks()
        v = self.num_points()
        gB = [[x+1 for x in b] for b in self._blocks]
        return "BlockDesign("+str(v)+","+str(gB)+")"

    def incidence_matrix(self):
        r"""
        Return the incidence matrix `A` of the design. A is a `(v \times b)`
        matrix defined by: ``A[i,j] = 1`` if ``i`` is in block ``B_j`` and 0
        otherwise.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.block_sizes()
            [3, 3, 3, 3, 3, 3, 3]
            sage: BD.incidence_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [1 0 0 0 0 1 1]
            [0 1 0 1 0 1 0]
            [0 1 0 0 1 0 1]
            [0 0 1 1 0 0 1]
            [0 0 1 0 1 1 0]

            sage: I = designs.IncidenceStructure('abc', ('ab','abc','ac','c'))
            sage: I.incidence_matrix()
            [1 1 1 0]
            [1 1 0 0]
            [0 1 1 1]
        """
        from sage.matrix.constructor import Matrix
        from sage.rings.all import ZZ
        A = Matrix(ZZ, self.num_points(), self.num_blocks(), sparse=True)
        for j, b in enumerate(self._blocks):
            for i in b:
                A[i, j] = 1
        return A

    def incidence_graph(self):
        """
        Return the incidence graph of the design, where the incidence
        matrix of the design is the adjacency matrix of the graph.

        EXAMPLE::

            sage: BD = designs.IncidenceStructure(7, [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.incidence_graph()
            Bipartite graph on 14 vertices
            sage: A = BD.incidence_matrix()
            sage: Graph(block_matrix([[A*0,A],[A.transpose(),A*0]])) == BD.incidence_graph()
            True

        REFERENCE:

        - Sage Reference Manual on Graphs
        """
        from sage.graphs.bipartite_graph import BipartiteGraph
        A = self.incidence_matrix()
        return BipartiteGraph(A)

    def relabel(self, perm=None, inplace=True):
        r"""
        Relabel the ground set

        INPUT:

        - ``perm`` -- can be one of

            - a dictionary -- then each point ``p`` (which should be a key of
              ``d``) is relabeled to ``d[p]``

            - a list or a tuple of length ``n`` -- the first point returned by
              :meth:`ground_set` is relabeled to ``l[0]``, the second to
              ``l[1]``, ...

            - ``None`` -- the incidence structure is relabeled to be on
              `\{0,1,...,n-1\}` in the ordering given by :meth:`ground_set`.

        - ``inplace`` -- If ``True`` then return a relabeled graph and does not
          touch ``self`` (default is ``False``).


        EXAMPLES::

            sage: TD=designs.transversal_design(5,5)
            sage: TD.relabel({i:chr(97+i) for i in range(25)})
            sage: print TD.ground_set()
            ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y']
            sage: print TD.blocks()[:3]
            [['a', 'f', 'k', 'p', 'u'], ['a', 'g', 'm', 's', 'y'], ['a', 'h', 'o', 'q', 'x']]

        Relabel to integer points::

            sage: TD.relabel()
            sage: print TD.blocks()[:3]
            [[0, 5, 10, 15, 20], [0, 6, 12, 18, 24], [0, 7, 14, 16, 23]]

        TESTS:

        Check that the relabel is consistent on a fixed incidence structure::

            sage: I = designs.IncidenceStructure([0,1,2,3,4],
            ....:               [[0,1,3],[0,2,4],[2,3,4],[0,1]])
            sage: I.relabel()
            sage: from itertools import permutations
            sage: for p in permutations([0,1,2,3,4]):
            ....:     J = I.relabel(p,inplace=False)
            ....:     if I == J: print p
            (0, 1, 2, 3, 4)
            (0, 1, 4, 3, 2)

        And one can also verify that we have exactly two automorphisms::

            sage: I.automorphism_group()
            Permutation Group with generators [(2,4)]
        """
        if not inplace:
            from copy import copy
            G = copy(self)
            G.relabel(perm=perm, inplace=True)
            return G

        if perm is None:
            self._points = range(self.num_points())
            self._point_to_index = None
            return

        if isinstance(perm, (list,tuple)):
            perm = dict(zip(self._points, perm))

        if not isinstance(perm, dict):
            raise ValueError("perm argument must be None, a list or a dictionary")

        if len(set(perm.values())) != len(perm):
            raise ValueError("Two points are getting relabelled with the same name !")

        self._points = [perm[x] for x in self._points]
        if self._points == range(self.num_points()):
            self._point_to_index  = None
        else:
            self._point_to_index = {v:i for i,v in enumerate(self._points)}

    def __hash__(self):
        r"""
        Not Implemented

        This object is mutable because of .relabel()

        EXAMPLE::

            sage: TD=designs.transversal_design(5,5)
            sage: hash(TD)
            Traceback (most recent call last):
            ...
            RuntimeError: This object is mutable !
        """
        raise RuntimeError("This object is mutable !")

    #####################
    # real computations #
    #####################

    def packing(self, solver=None, verbose=0):
        r"""
        Return a maximum packing

        A maximum packing in a hypergraph is collection of disjoint sets/blocks
        of maximal cardinality. This problem is NP-complete in general, and in
        particular on 3-uniform hypergraphs. It is solved here with an Integer
        Linear Program.

        For more information, see the :wikipedia:`Packing_in_a_hypergraph`.

        INPUT:

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP)
          solver to be used. If set to ``None``, the default one is used. For
          more information on LP solvers and which default solver is used, see
          the method
          :meth:`solve <sage.numerical.mip.MixedIntegerLinearProgram.solve>`
          of the class
          :class:`MixedIntegerLinearProgram <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.
          Only useful when ``algorithm == "LP"``.

        EXAMPLE::

            sage; IncidenceStructure([[1,2],[3,"A"],[2,3]]).packing()
            [[1, 2], [3, 'A']]
            sage: len(designs.steiner_triple_system(9).packing())
            3
        """
        from sage.numerical.mip import MixedIntegerLinearProgram

        # List of blocks containing a given point x
        d = [[] for x in self._points]
        for i,B in enumerate(self._blocks):
            for x in B:
                d[x].append(i)

        p = MixedIntegerLinearProgram(solver=solver)
        b = p.new_variable(binary=True)
        for x,L in enumerate(d): # Set of disjoint blocks
            p.add_constraint(p.sum([b[i] for i in L]) <= 1)

        # Maximum number of blocks
        p.set_objective(p.sum([b[i] for i in range(self.num_blocks())]))

        p.solve(log=verbose)

        return [[self._points[x] for x in self._blocks[i]]
                for i,v in p.get_values(b).iteritems() if v]

    def is_t_design(self, t=None, v=None, k=None, l=None, return_parameters=False):
        r"""
        Test whether ``self`` is a `t-(v,k,l)` design.

        A `t-(v,k,\lambda)` (sometimes called `t`-design for short) is a block
        design in which:

        - the underlying set has cardinality `v`
        - the blocks have size `k`
        - each `t`-subset of points is covered by `\lambda` blocks

        INPUT:

        - ``t,v,k,l`` (integers) -- their value is set to ``None`` by
          default. The function tests whether the design is a ``t-(v,k,l)``
          design using the provided values and guesses the others. Note that
          `l`` cannot be specified if ``t`` is not.

        - ``return_parameters`` (boolean)-- whether to return the parameters of
          the `t`-design. If set to ``True``, the function returns a pair
          ``(boolean_answer,(t,v,k,l))``.

        EXAMPLES::

            sage: fano_blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = designs.IncidenceStructure(7, fano_blocks)
            sage: BD.is_t_design()
            True
            sage: BD.is_t_design(return_parameters=True)
            (True, (2, 7, 3, 1))
            sage: BD.is_t_design(2, 7, 3, 1)
            True
            sage: BD.is_t_design(1, 7, 3, 3)
            True
            sage: BD.is_t_design(0, 7, 3, 7)
            True

            sage: BD.is_t_design(0,6,3,7) or BD.is_t_design(0,7,4,7) or BD.is_t_design(0,7,3,8)
            False

            sage: BD = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: BD.is_t_design(1)
            True
            sage: BD.is_t_design(2)
            True

        Steiner triple and quadruple systems are other names for `2-(v,3,1)` and
        `3-(v,4,1)` designs::

            sage: S3_9 = designs.steiner_triple_system(9)
            sage: S3_9.is_t_design(2,9,3,1)
            True

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = designs.IncidenceStructure(8, blocks)
            sage: S4_8.is_t_design(3,8,4,1)
            True

            sage: blocks = designs.steiner_quadruple_system(14)
            sage: S4_14 = designs.IncidenceStructure(14, blocks)
            sage: S4_14.is_t_design(3,14,4,1)
            True

        Some examples of Witt designs that need the gap database::

            sage: BD = designs.WittDesign(9)         # optional - gap_packages
            sage: BD.is_t_design(2,9,3,1)            # optional - gap_packages
            True
            sage: W12 = designs.WittDesign(12)       # optional - gap_packages
            sage: W12.is_t_design(5,12,6,1)          # optional - gap_packages
            True
            sage: W12.is_t_design(4)                 # optional - gap_packages
            True

        Further examples::

            sage: D = designs.IncidenceStructure(4,[[],[]])
            sage: D.is_t_design(return_parameters=True)
            (True,  (0, 4, 0, 2))

            sage: D = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (0, 4, 2, 3))

            sage: D = designs.IncidenceStructure(4, [[0],[1],[2],[3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (1, 4, 1, 1))

            sage: D = designs.IncidenceStructure(4,[[0,1],[2,3]])
            sage: D.is_t_design(return_parameters=True)
            (True, (1, 4, 2, 1))

            sage: D = designs.IncidenceStructure(4, [range(4)])
            sage: D.is_t_design(return_parameters=True)
            (True, (4, 4, 4, 1))

        TESTS::

            sage: blocks = designs.steiner_quadruple_system(8)
            sage: S4_8 = designs.IncidenceStructure(8, blocks)
            sage: R = range(15)
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(3,v,k,l)]
            [(8, 4, 1)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(2,v,k,l)]
            [(8, 4, 3)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(1,v,k,l)]
            [(8, 4, 7)]
            sage: [(v,k,l) for v in R for k in R for l in R if S4_8.is_t_design(0,v,k,l)]
            [(8, 4, 14)]
            sage: A = designs.AffineGeometryDesign(3, 1, GF(2))
            sage: A.is_t_design(return_parameters=True)
            (True, (2, 8, 2, 1))
            sage: A = designs.AffineGeometryDesign(4, 2, GF(2))
            sage: A.is_t_design(return_parameters=True)
            (True, (3, 16, 4, 1))
            sage: I = designs.IncidenceStructure(2, [])
            sage: I.is_t_design(return_parameters=True)
            (True, (0, 2, 0, 0))
            sage: I = designs.IncidenceStructure(2, [[0],[0,1]])
            sage: I.is_t_design(return_parameters=True)
            (False, (0, 0, 0, 0))
        """
        from sage.rings.arith import binomial

        # Missing parameters ?
        if v is None:
            v = self.num_points()

        if k is None:
            k = len(self._blocks[0]) if self._blocks else 0

        if l is not None and t is None:
            raise ValueError("t must be set when l=None")

        b = self.num_blocks()

        # Trivial wrong answers
        if (any(len(block) != k for block in self._blocks) or # non k-uniform
            v != self.num_points()):
            return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case t>k
        if (t is not None and t>k):
            if (l is None or l == 0):
                return (True, (t,v,k,0)) if return_parameters else True
            else:
                return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case k=0
        if k==0:
            if (l is None or l == 0):
                return (True, (0,v,k,b)) if return_parameters else True
            else:
                return (False, (0,0,0,0)) if return_parameters else False

        # Trivial case k=v (includes v=0)
        if k == v:
            if t is None:
                t = v
            if l is None or b == l:
                return (True, (t,v,k,b)) if return_parameters else True
            else:
                return (True, (0,0,0,0)) if return_parameters else False

        # Handbook of combinatorial design theorem II.4.8:
        #
        # a t-(v,k,l) is also a t'-(v,k,l')
        # for t' < t and l' = l* binomial(v-t',t-t') / binomial(k-t',t-t')
        #
        # We look for the largest t such that self is a t-design
        from itertools import combinations
        for tt in (range(1,k+1) if t is None else [t]):
            # is lambda an integer?
            if (b*binomial(k,tt)) % binomial(v,tt) != 0:
                tt -= 1
                break

            s = {}
            for block in self._blocks:
                for i in combinations(block,tt):
                    s[i] = s.get(i,0) + 1

            if len(set(s.values())) != 1:
                tt -= 1
                break

            ll = b*binomial(k,tt) // binomial(v,tt)

        if ((t is not None and t!=tt) or
            (l is not None and l!=ll)):
            return (False, (0,0,0,0)) if return_parameters else False
        else:
            if tt == 0:
                ll = b
            return (True, (tt,v,k,ll)) if return_parameters else True

    def dual(self, algorithm=None):
        """
        Return the dual of the incidence structure.

        INPUT:

        - ``algorithm`` -- whether to use Sage's implementation
          (``algorithm=None``, default) or use GAP's (``algorithm="gap"``).

          .. NOTE::

              The ``algorithm="gap"`` option requires GAP's Design package
              (included in the gap_packages Sage spkg).

        EXAMPLES:

        The dual of a projective plane is a projective plane::

            sage: PP = designs.DesarguesianProjectivePlaneDesign(4)
            sage: PP.dual().is_t_design(return_parameters=True)
            (True, (2, 21, 5, 1))

        TESTS::

            sage: D = designs.IncidenceStructure(4, [[0,2],[1,2,3],[2,3]])
            sage: D
            Incidence structure with 4 points and 3 blocks
            sage: D.dual()
            Incidence structure with 3 points and 4 blocks
            sage: print D.dual(algorithm="gap")       # optional - gap_packages
            IncidenceStructure<points=[0, 1, 2], blocks=[[0], [0, 1, 2], [1], [1, 2]]>
            sage: blocks = [[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]]
            sage: BD = designs.IncidenceStructure(7, blocks, name="FanoPlane");
            sage: BD
            Incidence structure with 7 points and 7 blocks
            sage: print BD.dual(algorithm="gap")         # optional - gap_packages
            IncidenceStructure<points=[0, 1, 2, 3, 4, 5, 6], blocks=[[0, 1, 2], [0, 3, 4], [0, 5, 6], [1, 3, 5], [1, 4, 6], [2, 3, 6], [2, 4, 5]]>
            sage: BD.dual()
            Incidence structure with 7 points and 7 blocks

        REFERENCE:

        - Soicher, Leonard, Design package manual, available at
          http://www.gap-system.org/Manuals/pkg/design/htm/CHAP003.htm
        """
        if algorithm == "gap":
            from sage.interfaces.gap import gap
            gap.load_package("design")
            gD = self._gap_()
            gap.eval("DD:=DualBlockDesign("+gD+")")
            v = eval(gap.eval("DD.v"))
            gblcks = eval(gap.eval("DD.blocks"))
            gB = []
            for b in gblcks:
                gB.append([x-1 for x in b])
            return IncidenceStructure(range(v), gB, name=None, check=False)
        else:
            return IncidenceStructure(
                          incidence_matrix=self.incidence_matrix().transpose(),
                          check=False)

    def automorphism_group(self):
        r"""
        Return the subgroup of the automorphism group of the incidence graph
        which respects the P B partition. It is (isomorphic to) the automorphism
        group of the block design, although the degrees differ.

        EXAMPLES::

            sage: P = designs.DesarguesianProjectivePlaneDesign(2); P
            Incidence structure with 7 points and 7 blocks
            sage: G = P.automorphism_group()
            sage: G.is_isomorphic(PGL(3,2))
            True
            sage: G
            Permutation Group with generators [(2,3)(4,5), (2,4)(3,5), (1,2)(4,6), (0,1)(4,5)]

        A non self-dual example::

            sage: IS = designs.IncidenceStructure(range(4), [[0,1,2,3],[1,2,3]])
            sage: IS.automorphism_group().cardinality()
            6
            sage: IS.dual().automorphism_group().cardinality()
            1

        Examples with non-integer points::

            sage: I = designs.IncidenceStructure('abc', ('ab','ac','bc'))
            sage: I.automorphism_group()
            Permutation Group with generators [('b','c'), ('a','b')]
            sage: designs.IncidenceStructure([[(1,2),(3,4)]]).automorphism_group()
            Permutation Group with generators [((1,2),(3,4))]
        """
        from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_element import standardize_generator
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        M1 = self.incidence_matrix().transpose()
        M2 = MatrixStruct(M1)
        M2.run()
        gens = M2.automorphism_group()[0]
        gens = [standardize_generator([x+1 for x in g]) for g in gens]
        if self._point_to_index:
            gens = [[tuple([self._points[i-1] for i in cycle]) for cycle in g] for g in gens]
        else:
            gens = [[tuple([i-1 for i in cycle]) for cycle in g] for g in gens]
        return PermutationGroup(gens, domain=self._points)

    def is_resolvable(self, certificate=False, solver=None, verbose=0, copy=True, check=True):
        r"""
        Test whether the hypergraph is resolvable

        A hypergraph is said to be resolvable if its sets can be partitionned
        into classes, each of which is a partition of the ground set.

        .. NOTE::

            This problem is solved using an Integer Linear Program, and GLPK
            (the default LP solver) has been reported to be very slow on some
            instances. If you hit this wall, consider installing a more powerful
            LP solver (CPLEX, Gurobi, ...).

        INPUT:

        - ``certificate`` (boolean) -- whether to return the classes along with
          the binary answer (see examples below).

        - ``solver`` -- (default: ``None``) Specify a Linear Program (LP) solver
          to be used. If set to ``None``, the default one is used. For more
          information on LP solvers and which default solver is used, see the
          method :meth:`solve
          <sage.numerical.mip.MixedIntegerLinearProgram.solve>` of the class
          :class:`MixedIntegerLinearProgram
          <sage.numerical.mip.MixedIntegerLinearProgram>`.

        - ``verbose`` -- integer (default: ``0``). Sets the level of
          verbosity. Set to 0 by default, which means quiet.

        - ``copy`` (boolean) -- ``True`` by default. When set to ``False``, a
          pointer toward the object's internal data is given. Set it to
          ``False`` only if you know what you are doing.

        - ``check`` (boolean) -- whether to check that output is correct before
          returning it. As this is expected to be useless (but we are cautious
          guys), you may want to disable it whenever you want speed. Set to ``True``
          by default.

        EXAMPLES:

        Some resolvable designs::

            sage: TD = designs.transversal_design(2,2,resolvable=True)
            sage: TD.is_resolvable()
            True

            sage: AG = designs.AffineGeometryDesign(3,1,GF(2))
            sage: AG.is_resolvable()
            True

        Their classes::

            sage: b,cls = TD.is_resolvable(True)
            sage: b
            True
            sage: cls # random
            [[[0, 3], [1, 2]], [[1, 3], [0, 2]]]

            sage: b,cls = AG.is_resolvable(True)
            sage: b
            True
            sage: cls # random
            [[[6, 7], [4, 5], [0, 1], [2, 3]],
             [[5, 7], [0, 4], [3, 6], [1, 2]],
             [[0, 2], [4, 7], [1, 3], [5, 6]],
             [[3, 4], [0, 7], [1, 5], [2, 6]],
             [[3, 7], [1, 6], [0, 5], [2, 4]],
             [[0, 6], [2, 7], [1, 4], [3, 5]],
             [[4, 6], [0, 3], [2, 5], [1, 7]]]

        A non-resolvable design::

            sage: Fano = designs.balanced_incomplete_block_design(7,3)
            sage: Fano.is_resolvable()
            False
            sage: Fano.is_resolvable(True)
            (False, [])

        TESTS::

            sage: _,cls1 = AG.is_resolvable(certificate=True, copy=True)
            sage: _,cls2 = AG.is_resolvable(certificate=True, copy=True)
            sage: cls1 is cls2
            False

            sage: _,cls1 = AG.is_resolvable(certificate=True, copy=False)
            sage: _,cls2 = AG.is_resolvable(certificate=True, copy=False)
            sage: cls1 is cls2
            True
        """
        if self._classes is None:
            degrees = set(self.degree().itervalues())
            if len(degrees) != 1:
                self._classes = False
            else:
                from sage.numerical.mip import MixedIntegerLinearProgram
                from sage.numerical.mip import MIPSolverException
                n_classes = degrees.pop()
                p = MixedIntegerLinearProgram(solver=solver)
                b = p.new_variable(binary=True)
                domain = range(self.num_points())

                # Lists of blocks containing i for every i
                dual = [[] for i in domain]
                for i,B in enumerate(self._blocks):
                    for x in B:
                        dual[x].append(i)

                # Each class is a partition
                for t in range(n_classes):
                    for x in domain:
                        p.add_constraint(p.sum(b[t,i] for i in dual[x]) == 1)

                # Each set appears exactly once
                for i in range(len(self._blocks)):
                    p.add_constraint(p.sum(b[t,i] for t in range(n_classes)) == 1)

                try:
                    p.solve(log=verbose)
                except MIPSolverException:
                    self._classes = False
                else:
                    # each class is stored as the list of indices of its blocks
                    self._classes = [[] for _ in range(n_classes)]
                    for (t,i),v in p.get_values(b).iteritems():
                        if v:
                            self._classes[t].append(self._blocks[i])

        if check and self._classes is not False:
            assert sorted(id(c) for cls in self._classes for c in cls) == sorted(id(b) for b in self._blocks), "some set does not appear exactly once"
            domain = range(self.num_points())
            for i,c in enumerate(self._classes):
                assert sorted(sum(c,[])) == domain, "class {} is not a partition".format(i)

        if self._classes is False:
            return (False, []) if certificate else False

        if certificate:
            if copy:
                if self._point_to_index is None:
                    classes = [[block[:] for block in classs] for classs in self._classes]
                else:
                    classes = [[[self._points[i] for i in block] for block in classs] for classs in self._classes]
            else:
                classes = self._classes

            return (True, classes)

        else:
            return True

    ###############
    # Deprecation #
    ###############

    def parameters(self):
        r"""
        Deprecated function. You should use :meth:`is_t_design` instead.

        EXAMPLES::

            sage: I = designs.IncidenceStructure('abc', ['ab','ac','bc'])
            sage: I.parameters()
            doctest:...: DeprecationWarning: .parameters() is deprecated. Use
            `is_t_design` instead
            See http://trac.sagemath.org/16553 for details.
            (2, 3, 2, 1)
        """
        from sage.misc.superseded import deprecation
        deprecation(16553, ".parameters() is deprecated. Use `is_t_design` instead")
        return self.is_t_design(return_parameters=True)[1]

    dual_design = deprecated_function_alias(16553, dual)
    dual_incidence_structure = deprecated_function_alias(16553, dual)
    is_block_design = deprecated_function_alias(16553, is_t_design)
    points = deprecated_function_alias(16553, ground_set)

    def block_design_checker(self, t, v, k, lmbda, type=None):
        """
        This method is deprecated and will soon be removed (see :trac:`16553`).
        You could use :meth:`is_t_design` instead.

        This is *not* a wrapper for GAP Design's IsBlockDesign. The GAP
        Design function IsBlockDesign
        http://www.gap-system.org/Manuals/pkg/design/htm/CHAP004.htm
        apparently simply checks the record structure and no mathematical
        properties. Instead, the function below checks some necessary (but
        not sufficient) "easy" identities arising from the identity.

        INPUT:

        - ``t`` - the t as in "t-design"

        - ``v`` - the number of points

        - ``k`` - the number of blocks incident to a point

        - ``lmbda`` - each t-tuple of points should be incident with
          lmbda blocks

        - ``type`` - can be 'simple' or 'binary' or 'connected'
          Depending on the option, this wraps IsBinaryBlockDesign,
          IsSimpleBlockDesign, or IsConnectedBlockDesign.

          - Binary: no block has a repeated element.

          - Simple: no block is repeated.

          - Connected: its incidence graph is a connected graph.

        WARNING: This is very fast but can return false positives.

        EXAMPLES::

            sage: BD = designs.IncidenceStructure(7,[[0,1,2],[0,3,4],[0,5,6],[1,3,5],[1,4,6],[2,3,6],[2,4,5]])
            sage: BD.is_t_design(return_parameters=True)
            (True, (2, 7, 3, 1))
            sage: BD.block_design_checker(2, 7, 3, 1)
            doctest:...: DeprecationWarning: .block_design_checker(v,t,k,lmbda) is deprecated; please use
            .is_t_design(v,t,k,lmbda) instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"binary")
            doctest:...: DeprecationWarning: .block_design_checker(type='binary') is
            deprecated; use .is_binary() instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"connected")
            doctest:...: DeprecationWarning: block_design_checker(type='connected') is
            deprecated, please use .is_connected() instead
            See http://trac.sagemath.org/16553 for details.
            True

            sage: BD.block_design_checker(2, 7, 3, 1,"simple")
            doctest:...: DeprecationWarning: .block_design_checker(type='simple')
            is deprecated; all designs here are simple!
            See http://trac.sagemath.org/16553 for details.
            True
        """
        from sage.misc.superseded import deprecation

        ans = self.is_t_design(t,v,k,lmbda)

        if type is None:
            deprecation(16553, ".block_design_checker(v,t,k,lmbda) is deprecated; please use .is_t_design(v,t,k,lmbda) instead")
            return ans

        if type == "binary":
            deprecation(16553, ".block_design_checker(type='binary') is deprecated; use .is_binary() instead")
            return True
        if type == "simple":
            deprecation(16553, ".block_design_checker(type='simple') is deprecated; all designs here are simple!")
            return True
        if type == "connected":
            deprecation(16553, "block_design_checker(type='connected') is deprecated, please use .is_connected() instead")
            return self.incidence_graph().is_connected()

    def edge_coloring(self):
        r"""
        Compute a proper edge-coloring.

        A proper edge-coloring is an assignment of colors to the sets of the
        incidence structure such that two sets with non-empty intersection
        receive different colors. The coloring returned minimizes the number of
        colors.

        OUTPUT:

        A partition of the sets into color classes.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: C = H.edge_coloring()
            sage: C # random
            [[[3, 4, 5]], [[2, 3, 4]], [[4, 5, 6], [1, 2, 3]]]
            sage: Set(map(Set,sum(C,[]))) == Set(map(Set,H.blocks()))
            True
        """
        from sage.graphs.graph import Graph
        blocks = self.blocks()
        blocks_sets = map(frozenset,blocks)
        g = Graph([range(self.num_blocks()),lambda x,y : len(blocks_sets[x]&blocks_sets[y])],loops = False)
        return [[blocks[i] for i in C] for C in g.coloring(algorithm="MILP")]

    def _spring_layout(self):
        r"""
        Return a spring layout for the points.

        The layout is computed by creating a graph `G` on the points *and* sets
        of the incidence structure. Each set is then made adjacent in `G` with
        all points it contains before a spring layout is computed for this
        graph. The position of the points in the graph gives the position of the
        points in the final drawing.

        .. NOTE::

            This method also returns the position of the "fake" points,
            i.e. those representing the sets.

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: L = H._spring_layout()
            sage: L # random
            {1: (0.238, -0.926),
             2: (0.672, -0.518),
             3: (0.449, -0.225),
             4: (0.782, 0.225),
             5: (0.558, 0.518),
             6: (0.992, 0.926),
             {3, 4, 5}: (0.504, 0.173),
             {2, 3, 4}: (0.727, -0.173),
             {4, 5, 6}: (0.838, 0.617),
             {1, 2, 3}: (0.393, -0.617)}
            sage: all(v in L for v in H.ground_set())
            True
            sage: all(v in L for v in map(Set,H.blocks()))
            True
        """
        from sage.graphs.graph import Graph

        g = Graph()
        for s in map(Set,self.blocks()):
            for x in s:
                g.add_edge(s,x)

        _ = g.plot(iterations = 50000,save_pos=True)

        # The values are rounded as TikZ does not like accuracy.
        return {k:(round(x,3),round(y,3)) for k,(x,y) in g.get_pos().items()}

    def _latex_(self):
        r"""
        Return a TikZ representation of the incidence structure

        EXAMPLES::

            sage: H = Hypergraph([{1,2,3},{2,3,4},{3,4,5},{4,5,6}]); H
            Incidence structure with 6 points and 4 blocks
            sage: view(H) # not tested

        With sets of size 4::

            sage: g = graphs.Grid2dGraph(5,5)
            sage: C4 = graphs.CycleGraph(4)
            sage: sets = Set(map(Set,list(g.subgraph_search_iterator(C4))))
            sage: H = Hypergraph(sets)
            sage: view(H) # not tested
        """
        from sage.rings.integer import Integer
        from sage.functions.trig import arctan2

        from sage.misc.misc import warn
        warn("\nThe hypergraph is drawn as a set of closed curves. The curve "
             "representing a set S go **THROUGH** the points contained "
             "in S.\n A point which is encircled by a curve but is not located "
             "on its boundary is **NOT** included in the corresponding set.\n"
             "\n"
             "The colors are picked for readability and have no other meaning.")

        latex.add_package_to_preamble_if_available("tikz")
        latex.add_to_mathjax_avoid_list("tikz")

        if not latex.has_file("tikz.sty"):
            raise RuntimeError("You must have TikZ installed in order "
                               "to draw a hypergraph.")

        domain = self.ground_set()
        pos = self._spring_layout()
        tex = "\\begin{tikzpicture}[scale=3]\n"

        colors = ["black", "red", "green", "blue", "cyan", "magenta", "yellow","pink","brown"]
        colored_sets = [(s,i) for i,S in enumerate(self.edge_coloring()) for s in S]

        # Prints each set with its color
        for s,i in colored_sets:
            current_color = colors[i%len(colors)]

            if len(s) == 2:
                s = list(s)
                tex += ("\\draw[color="+str(current_color)+","+
                        "line width=.1cm,opacity = .6] "+
                        str(pos[s[0]])+" -- "+str(pos[s[1]])+";\n")
                continue

            tex += ("\\draw[color="+str(current_color)+","
                    "line width=.1cm,opacity = .6,"
                    "line cap=round,"
                    "line join=round]"
                    "plot [smooth cycle,tension=1] coordinates {")

            # Reorders the vertices of s according to their angle with the
            # "center", i.e. the vertex representing the set s
            cx, cy = pos[Set(s)]
            s = map(lambda x: pos[x], s)
            s = sorted(s, key = lambda x_y: arctan2(x_y[0] - cx, x_y[1] - cy))

            for x in s:
                tex += str(x)+" "
            tex += "};\n"

        # Prints each vertex
        for v in domain:
            tex += "\\draw node[fill,circle,scale=.5,label={90:$"+latex(v)+"$}] at "+str(pos[v])+" {};\n"

        tex += "\\end{tikzpicture}"
        return tex

class GroupDivisibleDesign(IncidenceStructure):
    r"""
    Group Divisible Design (GDD)

    Let `K` and `G` be sets of positive integers and let `\lambda` be a positive
    integer. A Group Divisible Design of index `\lambda` and order `v` is a
    triple `(V,\mathcal G,\mathcal B)` where:

    - `V` is a set of cardinality `v`

    - `\mathcal G` is a partition of `V` into groups whose size belongs to `G`

    - `\mathcal B` is a family of subsets of `V` whose size belongs to `K` such
      that any two points `p_1,p_2\in V` from different groups appear
      simultaneously in exactly `\lambda` elements of `mathcal B`. Besides, a
      group and a block intersect on at most one point.

    If `K=\{k_1,...,k_k\}` and `G` has exactly `m_i` groups of cardinality `k_i`
    then `G` is said to have type `k_1^{m_1}...k_k^{m_k}`.

    INPUT:

    - ``points`` -- the underlying set. If ``points`` is an integer `v`, then
      the set is considered to be `\{0, ..., v-1\}`.

    - ``groups`` -- the groups of the design

    - ``blocks`` -- collection of blocks

    - ``G`` -- list of integers of which the sizes of the groups must be
      elements. Set to ``None`` (automatic guess) by default.

    - ``K`` -- list of integers of which the sizes of the blocks must be
      elements. Set to ``None`` (automatic guess) by default.

    - ``lambd`` (integer) -- value of `\lambda`, set to `1` by default.

    - ``check`` (boolean) -- whether to check that the design is indeed a `GDD`
      with the right parameters. Set to ``True`` by default.

    - ``copy`` -- (use with caution) if set to ``False`` then ``blocks`` must be
      a list of lists of integers. The list will not be copied but will be
      modified in place (each block is sorted, and the whole list is
      sorted). Your ``blocks`` object will become the instance's internal data.

    EXAMPLE::

        sage: from sage.combinat.designs.incidence_structures import GroupDivisibleDesign
        sage: TD = designs.transversal_design(4,10)
        sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
        sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
        Group Divisible Design on 40 points of type 10^4
    """
    def __init__(self, points, groups, blocks, G=None, K=None, lambd=1, check=True, copy=True,**kwds):
        r"""
        Constructor function

        EXAMPLE::

            sage: from sage.combinat.designs.incidence_structures import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
        """
        from designs_pyx import is_group_divisible_design

        self._lambd = lambd

        IncidenceStructure.__init__(self,
                                    points,
                                    blocks,
                                    copy=copy,
                                    check=False,
                                    **kwds)

        if copy is False and self._point_to_index is None:
            self._groups = groups
        elif self._point_to_index is None:
            self._groups = [g[:] for g in groups]
        else:
            self._groups = [[self._point_to_index[x] for x in g] for g in groups]

        if check:
            assert is_group_divisible_design(self._groups,self._blocks,self.num_points(),G,K,lambd)


    def groups(self, copy=True):
        r"""
        Return the groups of the Group-Divisible Design.

        INPUT:

        - ``copy`` (boolean) -- ``True`` by default. When set to ``False``, a
          pointer toward the object's internal data is given. Set it to
          ``False`` only if you know what you are doing.

        EXAMPLE::

            sage: from sage.combinat.designs.incidence_structures import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
            sage: GDD.groups()
            [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
             [10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
             [20, 21, 22, 23, 24, 25, 26, 27, 28, 29],
             [30, 31, 32, 33, 34, 35, 36, 37, 38, 39]]

        TESTS:

        Non-integer ground set::

            sage: TD=designs.transversal_design(5,5)
            sage: TD.relabel({i:chr(97+i) for i in range(25)})
            sage: TD.groups()
            [['a', 'b', 'c', 'd', 'e'],
             ['f', 'g', 'h', 'i', 'j'],
             ['k', 'l', 'm', 'n', 'o'],
             ['p', 'q', 'r', 's', 't'],
             ['u', 'v', 'w', 'x', 'y']]
        """
        if copy is False:
            return self._groups
        elif self._point_to_index is None:
            return [list(g) for g in self._groups]
        else:
            return [[self._points[i] for i in g] for g in self._groups]

    def __repr__(self):
        r"""
        Returns a string that describes self

        EXAMPLE::

            sage: from sage.combinat.designs.incidence_structures import GroupDivisibleDesign
            sage: TD = designs.transversal_design(4,10)
            sage: groups = [range(i*10,(i+1)*10) for i in range(4)]
            sage: GDD = GroupDivisibleDesign(40,groups,TD); GDD
            Group Divisible Design on 40 points of type 10^4
        """
        from string import join
        group_sizes = map(len, self.groups(copy=False))

        gdd_type = ["{}^{}".format(s,group_sizes.count(s))
                    for s in sorted(set(group_sizes))]
        gdd_type = join(gdd_type,".")

        if not gdd_type:
            gdd_type = "1^0"

        v = self.num_points()

        return "Group Divisible Design on {} points of type {}".format(v,gdd_type)
