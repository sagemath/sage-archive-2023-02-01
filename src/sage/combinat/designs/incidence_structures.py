"""
Incidence structures (i.e. hypergraphs, i.e. set systems)

An incidence structure is specified by a list of points, blocks, and
an incidence matrix ([1]_, [2]_).

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
        The number of points in that design.

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
        The number of blocks.

        EXAMPLES::

            sage: designs.DesarguesianProjectivePlaneDesign(2).num_blocks()
            7
            sage: B = designs.IncidenceStructure(4, [[0,1],[0,2],[0,3],[1,2], [1,2,3]])
            sage: B.num_blocks()
            5
        """
        return len(self._blocks)

    def blocks(self, copy=True):
        """Return the list of blocks.

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
                from copy import deepcopy
                return deepcopy(self._blocks)
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
        Returns the degree of a point ``p``

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
        Returns the incidence graph of the design, where the incidence
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

    #####################
    # real computations #
    #####################

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
        Returns the dual of the incidence structure.

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
        Returns the subgroup of the automorphism group of the incidence graph
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

        An example with points other than integers::

            sage: I = designs.IncidenceStructure('abc', ('ab','ac','bc'))
            sage: I.automorphism_group()
            Permutation Group with generators [('b','c'), ('a','b')]
        """
        from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        M1 = self.incidence_matrix().transpose()
        M2 = MatrixStruct(M1)
        M2.run()
        gens = M2.automorphism_group()[0]
        if self._point_to_index:
            gens = [[self._points[i] for i in p] for p in gens]
        return PermutationGroup(gens, domain=self._points)

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
