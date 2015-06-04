"""
Coxeter matrices
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#                     2015 Travis Scrimshaw <tscrim at ucdavis.edu>
#                     2015 Jean-Philippe Labbe <labbe at math.huji.ac.il>
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

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.matrix.constructor import matrix
from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.classcall_metaclass import ClasscallMetaclass, typecall
from sage.matrix.matrix_integer_dense import Matrix_integer_dense
from sage.matrix.matrix_generic_dense import Matrix_generic_dense
from sage.graphs.graph import Graph
from sage.graphs.generators.basic import CycleGraph
from sage.rings.all import ZZ, QQ, RR
from sage.rings.infinity import infinity
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.coxeter_type import CoxeterType
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.family import Family


class CoxeterMatrix(CoxeterType):
    """
    A Coxeter matrix.

    .. TODO::

        Because there is no object `\ZZ \cup \{ \infty \}`, we define `-1`
        to represent `\infty`.

    EXAMPLES::

        sage: CoxeterMatrix(['A', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
        sage: CoxeterMatrix(['B', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: CoxeterMatrix(['C', 4])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 4]
        [2 2 4 1]
        sage: CoxeterMatrix(['D', 4])
        [1 3 2 2]
        [3 1 3 3]
        [2 3 1 2]
        [2 3 2 1]

        sage: CoxeterMatrix(['E', 6])
        [1 2 3 2 2 2]
        [2 1 2 3 2 2]
        [3 2 1 3 2 2]
        [2 3 3 1 3 2]
        [2 2 2 3 1 3]
        [2 2 2 2 3 1]

        sage: CoxeterMatrix(['F', 4])
        [1 3 2 2]
        [3 1 4 2]
        [2 4 1 3]
        [2 2 3 1]

        sage: CoxeterMatrix(['G', 2])
        [1 6]
        [6 1]

    Because there currently is no class for `\ZZ \cup \{ \infty \}`, labels
    of `\infty` are given by `-1` in the Coxeter matrix::

        sage: G = Graph([(0,1,None), (1,2,4), (0,2,oo)])
        sage: CoxeterMatrix(G)
        [ 1  3 -1]
        [ 3  1  4]
        [-1  4  1]

    It is possible to give a number `\leq -1` to represent an infinite label

        sage: CoxeterMatrix([[1,-1],[-1,1]])
        [ 1 -1]
        [-1  1]
        sage: CoxeterMatrix([[1,-3/2],[-3/2,1]])
        [   1 -3/2]
        [-3/2    1]

    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **kwds):
        """
        A Coxeter matrix can we created via a graph, a Coxeter type, or
        a matrix.

        .. NOTE::

            To disable the Coxeter type check, use the optional argument
            ``coxeter_type_check = False``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1],['a','b'])
            sage: C2 = CoxeterMatrix([[1, -1], [-1, 1]])
            sage: C3 = CoxeterMatrix(matrix([[1, -1], [-1, 1]]), [0, 1])
            sage: C == C2 and C == C3
            True

        Check with `\infty` because of the hack of using `-1` to represent
        `\infty` in the Coxeter matrix::

            sage: G = Graph([(0, 1, 3), (1, 2, oo)])
            sage: W1 = CoxeterMatrix([[1, 3, 2], [3, 1, -1], [2, -1, 1]])
            sage: W2 = CoxeterMatrix(G)
            sage: W1 == W2
            True
            sage: CoxeterMatrix(W1.coxeter_graph()) == W1
            True

        The base ring of the matrix depends on the entries given::

            sage: CoxeterMatrix([[1,-1],[-1,1]])._matrix.base_ring()
            Integer Ring
            sage: CoxeterMatrix([[1,-3/2],[-3/2,1]])._matrix.base_ring()
            Rational Field
            sage: CoxeterMatrix([[1,-1.5],[-1.5,1]])._matrix.base_ring()
            Real Field with 53 bits of precision

        """

        # Special cases with 0 args
        if not args:
            if "coxeter_type" in kwds:  # kwds has Coxeter type
                args = ( CoxeterType(kwds["coxeter_type"]), )
            elif "cartan_type" in kwds:  # kwds has Cartan type
                args = ( CoxeterType(CartanType(kwds["cartan_type"])), )

            data = []
            n = 0
            index_set = tuple()
            coxeter_type = None
            base_ring = ZZ
            mat = typecall(cls, MatrixSpace(base_ring, n, sparse=False), data, coxeter_type, index_set)
            mat._subdivisions = None

            return mat

        elif len(args) == 4 and isinstance(args[0], MatrixSpace):  # For pickling
            return typecall(cls, args[0], args[1], args[2], args[3])
        elif isinstance(args[0], CoxeterMatrix):  # Initiate from itself
            return args[0]
        else:
            # Get the type check
            if kwds.get("coxeter_type_check", True):
                coxeter_type_check = True
            else:
                coxeter_type_check = False

            # Initiate from a graph:
            if isinstance(args[0], Graph):
                return cls._from_graph(args[0], coxeter_type_check)

            # Get the Coxeter type
            coxeter_type = None
            from sage.combinat.root_system.cartan_type import CartanType_abstract
            if isinstance(args[0], CartanType_abstract):
                coxeter_type = args[0].coxeter_type()
            else:
                try:
                    coxeter_type = CoxeterType(args[0])
                except (TypeError, ValueError, NotImplementedError):
                    pass

            # Initiate from a Coxeter type
            if coxeter_type:
                return cls._from_coxetertype(coxeter_type)

            # Get the index set
            n = len(list(args[0]))
            index_set = None
            if kwds.get("index_set", None):
                index_set = tuple(kwds["index_set"])
            if len(args) == 2:
                index_set = tuple(args[1])
            elif len(args) > 2:
                raise ValueError("too many arguments")
            if index_set and len(set(index_set)) != n:
                    raise ValueError("the given index set is not valid")

            return cls._from_matrix(args[0], coxeter_type, index_set, coxeter_type_check)

    def __init__(self, parent, data, coxeter_type, index_set):
        """
        Initialize ``self``.

        TESTS::

            sage: C = CoxeterMatrix(['A', 2, 1])
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])
        """

        self._matrix = Matrix_generic_dense(parent, data, False, True)
        self._matrix.set_immutable()

        if self._matrix.base_ring() not in [ZZ, QQ]:
            self._is_cyclotomic = False
        else:
            self._is_cyclotomic = True
        self._coxeter_type = coxeter_type

        if self._coxeter_type is not None:
            if self._coxeter_type.is_finite():
                self._is_finite = True
                self._is_affine = False
            elif self._coxeter_type.is_affine():
                self._is_finite = False
                self._is_affine = True
            else:
                self._is_finite = False
                self._is_affine = False
        else:
            self._is_finite = False
            self._is_affine = False

        self._index_set = index_set
        self._rank = self._matrix.nrows()

        self._dict = {(self._index_set[i], self._index_set[j]): self._matrix[i, j]
                for i in range(self._rank) for j in range(self._rank)}

        for key in self._index_set:
            index_key = self._index_set.index(key)
            self._dict[key] = {i: self._matrix[index_key, self._index_set.index(i)] for i in self._index_set}

    @classmethod
    def _from_matrix(cls, data, coxeter_type, index_set, coxeter_type_check):
        """
        Initiate the Coxeter matrix from a matrix.

        TESTS::

            sage: CM = CoxeterMatrix([[1,2],[2,1]]);CM
            [1 2]
            [2 1]
            sage: CM = CoxeterMatrix([[1,-1],[-1,1]]);CM
            [ 1 -1]
            [-1  1]
            sage: CM = CoxeterMatrix([[1,-1.5],[-1.5,1]]);CM
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
            sage: CM = CoxeterMatrix([[1,-3/2],[-3/2,1]]);CM
            [   1 -3/2]
            [-3/2    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,-1],[5,-1,1]]);CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
            sage: CM = CoxeterMatrix([[1,-3/2,5],[-3/2,1,oo],[5,oo,1]]);CM
            [   1 -3/2    5]
            [-3/2    1   -1]
            [   5   -1    1]
        """

        # Check that the data is valid
        check_coxeter_matrix(data)

        M = matrix(data)
        n = M.ncols()

        # TODO:: remove when oo is possible in matrices.
        entries = []
        for r in data:
            entries += list(r)
        raw_data = map(lambda x: x if x != infinity else -1, entries)
        M = matrix(n, n, raw_data)
        # until here

        base_ring = M.base_ring()

        if not coxeter_type:
            if n == 1:
                coxeter_type = CoxeterType(['A', 1])
            elif coxeter_type_check:
                coxeter_type = recognize_coxeter_type_from_matrix(M)
            else:
                coxeter_type = None
        if not index_set:
            index_set = tuple(range(1,n+1))

        raw_data = M.list()

        mat = typecall(cls, MatrixSpace(base_ring, n, sparse=False), raw_data,
                       coxeter_type, index_set)
        mat._subdivisions = M._subdivisions

        return mat

    @classmethod
    def _from_graph(cls, graph, coxeter_type_check):
        """
        Initiate the Coxeter matrix from a graph.

        TESTS::

            sage: CoxeterMatrix(CoxeterMatrix(['A',4,1]).coxeter_graph())
            [1 3 2 2 3]
            [3 1 3 2 2]
            [2 3 1 3 2]
            [2 2 3 1 3]
            [3 2 2 3 1]
            sage: CoxeterMatrix(CoxeterMatrix(['B',4,1]).coxeter_graph())
            [1 2 3 2 2]
            [2 1 3 2 2]
            [3 3 1 3 2]
            [2 2 3 1 4]
            [2 2 2 4 1]
            sage: CoxeterMatrix(CoxeterMatrix(['F',4]).coxeter_graph())
            [1 3 2 2]
            [3 1 4 2]
            [2 4 1 3]
            [2 2 3 1]

            sage: G=Graph()
            sage: G.add_edge([0,1,oo])
            sage: CoxeterMatrix(G)
            [ 1 -1]
            [-1  1]
            sage: H = Graph()
            sage: H.add_edge([0,1,-1.5])
            sage: CoxeterMatrix(H)
            [ 1.00000000000000 -1.50000000000000]
            [-1.50000000000000  1.00000000000000]
        """

        verts = sorted(graph.vertices())
        index_set = tuple(verts)
        n = len(index_set)

        # Setup the basis matrix as all 2 except 1 on the diagonal
        data = []
        for i in range(n):
            data += [[]]
            for j in range(n):
                if i == j:
                    data[-1] += [ZZ.one()]
                else:
                    data[-1] += [2]

        for e in graph.edges():
            label = e[2]
            if label is None:
                label = 3
            elif label == infinity:
                label = -1
            elif label not in ZZ and label > -1:
                raise ValueError("invalid Coxeter graph label")
            elif label == 0 or label == 1:
                raise ValueError("invalid Coxeter graph label")
            i = verts.index(e[0])
            j = verts.index(e[1])
            data[j][i] = data[i][j] = label

        return cls._from_matrix(data, None, index_set, coxeter_type_check)

    @classmethod
    def _from_coxetertype(cls, coxeter_type):
        """
        Initiate the Coxeter matrix from a Coxeter type.

        TESTS::

            sage: CoxeterMatrix(['A',4]).coxeter_type()
            Coxeter type of ['A', 4]
            sage: CoxeterMatrix(['A',4,1]).coxeter_type()
            Coxeter type of ['A', 4, 1]
            sage: CoxeterMatrix(['D',4,1]).coxeter_type()
            Coxeter type of ['D', 4, 1]
        """
        index_set = coxeter_type.index_set()
        n = len(index_set)
        reverse = {index_set[i]: i for i in range(n)}
        data = [[1 if i == j else 2 for j in range(n)] for i in range(n)]
        for (i, j, l) in coxeter_type.coxeter_graph().edge_iterator():
            if l == infinity:
                l = -1
            data[reverse[i]][reverse[j]] = l
            data[reverse[j]][reverse[i]] = l

        return cls._from_matrix(data, coxeter_type, index_set, False)

    @classmethod
    def samples(self, finite=None, affine=None, crystallographic=None, higher_rank=None):
        """
        Return a sample of the available Coxeter types.

        INPUT:

        - ``finite`` -- a boolean or ``None`` (default: ``None``)

        - ``affine`` -- a boolean or ``None`` (default: ``None``)

        - ``crystallographic`` -- a boolean or ``None`` (default: ``None``)

        - ``higher_rank`` -- a boolean or ``None`` (default: ``None``)

        The sample contains all the exceptional finite and affine
        Coxeter types, as well as typical representatives of the
        infinite families.

        Here the ``higher_rank`` term denotes non-finite, non-affine, 
        Coxeter groups (including hyperbolic types).

        .. TODO:: Implement the hyperbolic and compact hyperbolic in the
        samples.

        EXAMPLES::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples()]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            Coxeter type of ['F', 4], Coxeter type of ['H', 3],
            Coxeter type of ['H', 4], Coxeter type of ['I', 10],
            Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
            Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
            Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
            Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
            Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1],
            [ 1 -2]
            [-2  1],
            [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1],
            [1 2 3]
            [2 1 7]
            [3 7 1],
            [ 1 -2  3  2]
            [-2  1  2  3]
            [ 3  2  1 -8]
            [ 2  3 -8  1]]

        The finite, affine and crystallographic options allow
        respectively for restricting to (non) finite, (non) affine,
        and (non) crystallographic Cartan types::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(finite=True)]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            Coxeter type of ['F', 4], Coxeter type of ['H', 3],
            Coxeter type of ['H', 4], Coxeter type of ['I', 10]]

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(affine=True)]
            [Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
            Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
            Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
            Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
            Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1]]

            sage: [CM.coxeter_type() for CM in CoxeterMatrix.samples(crystallographic=True)]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            Coxeter type of ['F', 4], Coxeter type of ['A', 2, 1],
            Coxeter type of ['B', 5, 1], Coxeter type of ['C', 5, 1],
            Coxeter type of ['D', 5, 1], Coxeter type of ['E', 6, 1],
            Coxeter type of ['E', 7, 1], Coxeter type of ['E', 8, 1],
            Coxeter type of ['F', 4, 1], Coxeter type of ['G', 2, 1]]

            sage: CoxeterMatrix.samples(crystallographic=False)
            [[1 3 2]
            [3 1 5]
            [2 5 1], [1 3 2 2]
            [3 1 3 2]
            [2 3 1 5]
            [2 2 5 1], [ 1 10]
            [10  1], [ 1 -1]
            [-1  1], [ 1 -2]
            [-2  1], [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1], [1 2 3]
            [2 1 7]
            [3 7 1], [ 1 -2  3  2]
            [-2  1  2  3]
            [ 3  2  1 -8]
            [ 2  3 -8  1]]

        .. TODO:: add some reducible Coxeter types (suggestions?)

        TESTS::

            sage: for ct in CoxeterMatrix.samples(): TestSuite(ct).run()
        """
        result = self._samples()
        if crystallographic is not None:
            result = [t for t in result if t.is_crystallographic() == crystallographic]
        if finite is not None:
            result = [t for t in result if t.is_finite() == finite]
        if affine is not None:
            result = [t for t in result if t.is_affine() == affine]
        if higher_rank is not None:
            result = [t for t in result if not t.is_affine() and not t.is_finite()]
        return result

    @cached_method
    def _samples(self):
        """
        Return a sample of all implemented Coxeter types.

        .. NOTE:: This is intended to be used through :meth:`samples`.

        EXAMPLES::

            sage: [CM.coxeter_type() for CM in CoxeterMatrix._samples()]
            [Coxeter type of ['A', 1], Coxeter type of ['A', 5],
            Coxeter type of ['B', 5], Coxeter type of ['D', 4],
            Coxeter type of ['D', 5], Coxeter type of ['E', 6],
            Coxeter type of ['E', 7], Coxeter type of ['E', 8],
            Coxeter type of ['F', 4], Coxeter type of ['H', 3],
            Coxeter type of ['H', 4], Coxeter type of ['I', 10],
            Coxeter type of ['A', 2, 1], Coxeter type of ['B', 5, 1],
            Coxeter type of ['C', 5, 1], Coxeter type of ['D', 5, 1],
            Coxeter type of ['E', 6, 1], Coxeter type of ['E', 7, 1],
            Coxeter type of ['E', 8, 1], Coxeter type of ['F', 4, 1],
            Coxeter type of ['G', 2, 1], Coxeter type of ['A', 1, 1],
            [ 1 -2]
            [-2  1],
            [ 1 -1 -1]
            [-1  1 -1]
            [-1 -1  1],
            [1 2 3]
            [2 1 7]
            [3 7 1],
            [ 1 -2  3  2]
            [-2  1  2  3]
            [ 3  2  1 -8]
            [ 2  3 -8  1]]
        """
        finite = [CoxeterMatrix(t)       for t in [['A', 1], ['A', 5], ['B', 5],
                                            ['D', 4], ['D', 5], ['E', 6], ['E', 7],
                                            ['E', 8], ['F', 4], ['H', 3], ['H', 4],
                                            ['I', 10]]]

        affine = [CoxeterMatrix(t)      for t in ['A', 2, 1], ['B', 5, 1],
                                            ['C', 5, 1], ['D', 5, 1], ['E', 6, 1],
                                            ['E', 7, 1], ['E', 8, 1], ['F', 4, 1],
                                            ['G', 2, 1], ['A', 1, 1]]

        higher_matrices = [[[1, -2], [-2, 1]],
                [[1, -1, -1], [-1, 1, -1], [-1, -1, 1]],
                [[1, 2, 3], [2, 1, 7], [3, 7, 1]],
                [[1, -2, 3, 2], [-2, 1, 2, 3], [3, 2, 1, -8], [2, 3, -8, 1]]]

        higher = [CoxeterMatrix(m) for m in higher_matrices]

        return finite + affine + higher

    def relabel(self, relabelling):
        """
        Return a relabelled copy of this Coxeter matrix.

        INPUT:

        - ``relabelling`` -- a function (or dictionary)

        OUTPUT:

        an isomorphic Coxeter type obtained by relabelling the nodes of
        the Coxeter graph. Namely, the node with label ``i`` is
        relabelled ``f(i)`` (or, by ``f[i]`` if ``f`` is a dictionary).

        EXAMPLES::

            sage: CoxeterMatrix(['F',4]).relabel({ 1:2, 2:3, 3:4, 4:1})
            [1 4 2 3]
            [4 1 3 2]
            [2 3 1 2]
            [3 2 2 1]
            sage: CoxeterMatrix(['F',4]).relabel(lambda x: x+1 if x<4 else 1)
            [1 4 2 3]
            [4 1 3 2]
            [2 3 1 2]
            [3 2 2 1]
        """

        if isinstance(relabelling,type({})):
            data = [[self[relabelling[i]][relabelling[j]] for j in self.index_set()] for i in self.index_set()]
        else:
            data = [[self[relabelling(i)][relabelling(j)] for j in self.index_set()] for i in self.index_set()]

        return CoxeterMatrix(data)

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: C = CoxeterMatrix(['A',4])
            sage: M = loads(dumps(C))
            sage: M._index_set
            (1, 2, 3, 4)
        """
        return (CoxeterMatrix, (self._matrix.parent(), self._matrix.list(),
                                self._coxeter_type, self._index_set))

    def __repr__(self):
        """
        String representation of the Coxeter matrix.
        """

        return self._matrix.__repr__()

    def __iter__(self):
        """
        Return an iterator for the rows of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,8],[8,1]])
            sage: CM.__iter__().next()
            (1, 8)
        """

        return self._matrix.__iter__()

    def __getitem__(self, key):
        """
        Return a dictionary of labels adjacent to a node or
        the label of an edge in the Coxeter graph.

        EXAMPLES::
            
            sage: CM = CoxeterMatrix([[1,-2],[-2,1]])
            sage: CM = CoxeterMatrix([[1,-2],[-2,1]],['a','b'])
            sage: CM['a']
            {'a': 1, 'b': -2}
            sage: CM['b']
            {'a': -2, 'b': 1}
            sage: CM['a','b']
            -2
            sage: CM['a','a']
            1
        """

        return self._dict[key]

    def __hash__(self):
        r"""
        Return hash of the Coxeter matrix.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-2],[-2,1]],['a','b'])
            sage: CM.__hash__()
            1
            sage: CM = CoxeterMatrix([[1,-3],[-3,1]],['1','2'])
            sage: CM.__hash__()
            4
                                                                                                                                                                """
        
        return self._matrix.__hash__()

    def __eq__(self, other):
        r"""
        Return if ``self`` and ``other`` are equal, ``False`` otherwise.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-2],[-2,1]],['a','b'])
            sage: CM.__hash__()
            1
            sage: CM = CoxeterMatrix([[1,-3],[-3,1]],['1','2'])
            sage: CM.__hash__()
            4
                                                                                                                                                                """

        return self._matrix.__eq__(other._matrix)

    def _matrix_(self, R = None):
        """
        Return ``self`` as a matrix over the ring ``R``.

        EXAMPLES::

            sage: CM = CoxeterMatrix([[1,-3],[-3,1]])
            sage: matrix(CM)
            [ 1 -3]
            [-3  1]
            sage: matrix(CM,RR)
            [ 1.00000000000000 -3.00000000000000]
            [-3.00000000000000  1.00000000000000]
        """

        if R is not None:
            return self._matrix.change_ring(R)
        else:
            return self._matrix

    ##########################################################################
    # Coxeter type methods

    def index_set(self):
        """
        Return the index set of ``self``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1])
            sage: C.index_set()
            (0, 1)
            sage: C = CoxeterMatrix(['E',6])
            sage: C.index_set()
            (1, 2, 3, 4, 5, 6)
        """
        return self._index_set

    def coxeter_type(self):
        """
        Return the Coxeter type of ``self`` or ``self`` if unknown.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',4,1])
            sage: C.coxeter_type()
            Coxeter type of ['A', 4, 1]

        If the Coxeter type is unknown::

            sage: C = CoxeterMatrix([[1,3,4], [3,1,-1], [4,-1,1]])
            sage: C.coxeter_type()
            [ 1  3  4]
            [ 3  1 -1]
            [ 4 -1  1]
        """
        if self._coxeter_type is None:
            return self
        return self._coxeter_type

    def rank(self):
        r"""
        Return the rank of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).rank()
            3
            sage: CoxeterMatrix(["A2","B2","F4"]).rank()
            8
        """
        return self._rank

    def coxeter_matrix(self):
        r"""
        Return the Coxeter matrix of ``self``.

        EXAMPLES::

            sage: CoxeterMatrix(['C',3]).coxeter_matrix()
            [1 3 2]
            [3 1 4]
            [2 4 1]
        """
        return self

    def bilinear_form(self):
        r"""
        Return the bilinear form of ``self``.

        EXAMPLES::
        """

        return CoxeterType.bilinear_form(self)

    @cached_method
    def coxeter_graph(self):
        """
        Return the Coxeter graph of ``self``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',3])
            sage: C.coxeter_graph()
            Graph on 3 vertices

            sage: C = CoxeterMatrix([['A',3],['A',1]])
            sage: C.coxeter_graph()
            Graph on 4 vertices
        """
        n = self.rank()
        I = self.index_set()
        val = lambda x: infinity if x == -1 else x
        G = Graph([(I[i], I[j], val((self._matrix)[i, j]))
                   for i in range(n) for j in range(i)
                   if (self._matrix)[i, j] not in [1, 2]])
        G.add_vertices(I)
        return G.copy(immutable = True)

    def is_simply_laced(self):
        """
        Return if ``self`` is simply-laced.

        A Coxeter matrix is simply-laced if all non-diagonal entries are
        either 2 or 3.

        EXAMPLES::

            sage: cm = CoxeterMatrix([[1,3,3,3], [3,1,3,3], [3,3,1,3], [3,3,3,1]])
            sage: cm.is_simply_laced()
            True
        """
        # We include 1 in this list to account for the diagonal
        L = [1, 2, 3]
        return all(x in L for row in self for x in row)

    def is_crystallographic(self):
        """
        Return if ``self`` is crystallographic.

        A Coxeter matrix is crystallographic if all non-diagonal entries
        are either 2, 4, or 6.

        EXAMPLES::

            sage: CoxeterMatrix(['F',4]).is_crystallographic()
            True
            sage: CoxeterMatrix(['H',3]).is_crystallographic()
            False
        """
        # We include 1 in this list to account for the diagonal
        L = [1, 2, 3, 4, 6]
        return all(x in L for row in self for x in row)

    def is_finite(self):
        """
        Return if ``self`` is a finite type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_finite()
            True
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_finite()
            False
            sage: M = CoxeterMatrix([[1, -1], [-1, 1]])
            sage: M.is_finite()
            False
        """
        return self._is_finite

    def is_affine(self):
        """
        Return if ``self`` is an affine type or ``False`` if unknown.

        EXAMPLES::

            sage: M = CoxeterMatrix(['C',4])
            sage: M.is_affine()
            False
            sage: M = CoxeterMatrix(['D',4,1])
            sage: M.is_affine()
            True
            sage: M = CoxeterMatrix([[1, 3],[3,1]])
            sage: M.is_affine()
            False
            sage: M = CoxeterMatrix([[1, -1, 7], [-1, 1, 3], [7, 3, 1]])
            sage: M.is_affine()
            False
        """
        return self._is_affine


#####################################################################
## Type check functions

def recognize_coxeter_type_from_matrix(coxeter_matrix):
    """
    Return the Coxeter type of ``coxeter_matrix`` if known,
    otherwise return ``None``.

    EXAMPLES:

    Some infinite ones::

        sage: C = CoxeterMatrix([[1,3,2],[3,1,-1],[2,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False
        sage: C = CoxeterMatrix([[1,-1,-1],[-1,1,-1],[-1,-1,1]])
        sage: C.is_finite()  # indirect doctest
        False

    Some finite ones::

        sage: m = matrix(CoxeterMatrix(['D', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True
        sage: m = matrix(CoxeterMatrix(['H', 4]))
        sage: CoxeterMatrix(m).is_finite()  # indirect doctest
        True

        sage: CoxeterMatrix(CoxeterType(['A',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10]
        sage: CoxeterMatrix(CoxeterType(['B',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['C',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10]
        sage: CoxeterMatrix(CoxeterType(['D',10]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10]
        sage: CoxeterMatrix(CoxeterType(['E',6]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6]
        sage: CoxeterMatrix(CoxeterType(['E',7]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7]
        sage: CoxeterMatrix(CoxeterType(['E',8]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8]
        sage: CoxeterMatrix(CoxeterType(['F',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4]
        sage: CoxeterMatrix(CoxeterType(['G',2]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2]
        sage: CoxeterMatrix(CoxeterType(['H',3]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 3]
        sage: CoxeterMatrix(CoxeterType(['H',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 4]
        sage: CoxeterMatrix(CoxeterType(['I',100]).coxeter_graph()).coxeter_type()
        Coxeter type of ['I', 100]

    Some affine graphs::

        sage: CoxeterMatrix(CoxeterType(['A',1,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 1, 1]
        sage: CoxeterMatrix(CoxeterType(['A',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['B',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['C',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['C', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['D',10,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['D', 10, 1]
        sage: CoxeterMatrix(CoxeterType(['E',6,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 6, 1]
        sage: CoxeterMatrix(CoxeterType(['E',7,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 7, 1]
        sage: CoxeterMatrix(CoxeterType(['E',8,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['E', 8, 1]
        sage: CoxeterMatrix(CoxeterType(['F',4,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4, 1]
        sage: CoxeterMatrix(CoxeterType(['G',2,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['G', 2, 1]

    TESTS::

        sage: from sage.combinat.root_system.coxeter_matrix import recognize_coxeter_type_from_matrix
        sage: for C in CoxeterMatrix.samples():
        ....:     relabelling_perm = Permutations(C.index_set()).random_element()
        ....:     relabelling_dict = {C.index_set()[i]: relabelling_perm[i] for i in range(C.rank())}
        ....:     relabeled_matrix = C.relabel(relabelling_dict)._matrix
        ....:     recognized_type = recognize_coxeter_type_from_matrix(relabeled_matrix)
        ....:     if C.is_finite() or C.is_affine():
        ....:         assert recognized_type == C.coxeter_type()

    """
    # First, we build the Coxeter graph of the group without the edge labels
    n = ZZ(coxeter_matrix.nrows())
    G = Graph([range(n), lambda i, j: coxeter_matrix[i, j] not in [1, 2]])
    # Coxeter graphs of finite Coxeter groups are forests
    # Coxeter graphs of affine Coxeter groups are forests possibly with cycle
    # graphs

    comps = G.connected_components()
    # The group is finite if and only if for every connected
    # component ``comp`` of its Coxeter graph, the submatrix of
    # the Coxeter matrix corresponding to ``comp`` is one of the
    # type-A,B,D,E,F,H,I matrices (up to permutation). So we
    # shall check this condition on every ``comp``.
    types = []
    for comp in comps:
        l = len(comp)
        G0 = G.subgraph(comp)
        if l == 1:
            # Any `1 \times 1` Coxeter matrix gives a finite group.
            types.append(['A', 1])
            continue  # A1
        elif l == 2:
            # A finite dihedral group iff there is no `\infty` in its
            # Coxeter matrix. Otherwise it is affine iff the
            # off-diagonal value is -1.
            c0, c1 = comp
            if coxeter_matrix[c0, c1] > 0:
                types.append(['I', coxeter_matrix[c0, c1]])
                continue
            elif coxeter_matrix[c0, c1] == -1:
                types.append(['A', 1, 1])
                continue
            else:
                return None  # TODO: return hyperbolic type once implemented
        elif l == 3:
            # The `3`-node case. The finite groups to check for
            # here are `A_3`, `B_3` and `H_3`.
            # The affine groups to check for here are `A_{2,1}`, `G_{2,1}`,
            # `B_{2,1}`.
            c0, c1, c2 = comp
            s = sorted([coxeter_matrix[c0, c1],
                        coxeter_matrix[c0, c2],
                        coxeter_matrix[c1, c2]])
            if s[0] == 2:  # Have a tree
                if s[1] == 3:
                    if s[2] == 3:
                        types.append(['A', 3])
                        continue
                    elif s[2] == 4:
                        types.append(['B', 3])
                        continue
                    elif s[2] == 5:
                        types.append(['H', 3])
                        continue
                    elif s[2] == 6:
                        types.append(['G', 2, 1])
                        continue
                    else:
                        return None
                elif s[1] == 4 and s[2] == 4:
                    types.append(['B', 2, 1])
                    continue
                else:
                    return None
            elif s[0] == 3 and s[1] == 3 and s[2] == 3:
                types.append(['A', 2, 1])
                continue
            else:  # Have a hyperbolic type of rank 3
                return None
        elif l == 4:
            # The `4`-node case. The finite groups to check for
            # here are `A_4`, `B_4`, `D_4`, `F_4` and `H_4`.
            # The affine groups to check for here are `A_{3,1}`, `B_{3,1}`,
            # and `C_{3,1}`.

            c0, c1, c2, c3 = comp
            u = [coxeter_matrix[c0, c1],
                 coxeter_matrix[c0, c2],
                 coxeter_matrix[c0, c3],
                 coxeter_matrix[c1, c2],
                 coxeter_matrix[c1, c3],
                 coxeter_matrix[c2, c3]]
            s = sorted(u)
            # ``s`` is the list of all off-diagonal entries of
            # the ``comp``-submatrix of the Coxeter matrix,
            # sorted in increasing order.

            if s[:3] == [2, 2, 2]:
                if s[3:] == [3, 3, 3]:
                    if max(G0.degree()) == 2:
                        types.append(['A', 4])
                        continue
                    else:
                        types.append(['D', 4])
                        continue
                elif s[3:] == [3, 3, 4] or s[3:] == [3, 3, 5] or s[3:] == [3,
                        4, 4]:
                    if max(G0.degree()) == 3:
                        if s[4:] == [3, 4]:
                            types.append(['B', 3, 1])
                            continue
                        else:
                            return None
                    else:  # The graph is a path
                           # Differenciate using sum of edge labels
                        u0 = u[0] + u[1] + u[2]
                        u1 = u[0] + u[3] + u[4]
                        u2 = u[1] + u[3] + u[5]
                        u3 = u[2] + u[4] + u[5]
                        ss = sorted([u0, u1, u2, u3])
                        if s[5] == 4:
                            if ss == [7, 7, 9, 9]:
                                types.append(['F', 4])
                                continue
                            elif ss == [7, 8, 8, 9]:
                                types.append(['B', 4])
                                continue
                            elif ss == [8, 8, 9, 9]:
                                types.append(['C', 3, 1])
                                continue
                            else:
                                return None
                        elif ss == [7, 8, 9, 10]:
                            types.append(['H', 4])
                            continue
                        else:
                            return None

                else:
                    return None

            elif s == [2, 2, 3, 3, 3, 3] and max(G0.degree()) == 2:
                types.append(['A', 3, 1])
                continue
            else:
                return None
            return None
        else:
            # The case of `l \geq 5` nodes. The finite
            # groups to check for here are `A_l`, `B_l`, `D_l`,
            # and `E_l` (for `l = 6, 7, 8`).

            # The affine groups to check for here are `A_{l-1,1}`, `B_{l-1,1}`,
            # `C_{l-1,1}`, `D_{l-1,1}`, `E_{l_1,1}` (for `l = 6, 7, 8`),
            # or `F_{4,1}`.

            degrees = G0.degree()
            vertices = G0.vertices()
            sub_cox_matrix = coxeter_matrix.matrix_from_rows_and_columns(vertices, vertices)
            vertices_labels = [set(filter(lambda x: x != 1, r)) for r in
                    sub_cox_matrix.rows()]
            label_list = filter(lambda x: x != 1, sub_cox_matrix.list())
            labels = sorted(set(label_list))

            vertices_4 = [index for index in range(l) if 4 in
                    vertices_labels[index]]
            occur_4 = label_list.count(4)/2  # Each label appear twice

            if not G0.is_tree():
                if G0.is_isomorphic(CycleGraph(l)):  # Type `A_{l-1,1}`
                    types.append(['A', l-1, 1])
                    continue
                else:
                    return None

            elif max(degrees) == 2:  # The component is a path

                if labels[-1] == 3:  # Highest label is 3
                    types.append(['A', l])
                    continue
                elif labels[-1] == 4:  # Highest label is 4

                    if occur_4 == 1:  # There is 1 edge with label 4
                        if not (3 in vertices_labels[vertices_4[0]] and 3 in
                        vertices_labels[vertices_4[1]]):  # The edge is at the end of the path
                            types.append(['B', l])
                            continue
                        elif l == 5:
                            types.append(['F', 4, 1])
                            continue
                        else:
                            return None
                    elif occur_4 == 2:  # There are 2 edges labeled 4
                        if len(filter(lambda x: 2 in x and 3 not in x and 4 in
                                x, [vertices_labels[i] for i in vertices_4])) == 2:
                            # The edges with 4 are at the ends of the path
                            types.append(['C', l-1, 1])
                            continue
                    else:
                        return None
                else:
                    return None

            else:  # The graph contains branching vertices

                #Finite D_n,E_6,E_7,E_8
                #Affine B_n,D_n, E_6,E_7,E_8

                if max(degrees) == 3:  # Branching degree is 3

                    highest_label = labels[-1]

                    nb_branching = degrees.count(3)
                    ecc = sorted(G0.eccentricity())

                    if nb_branching == 1:

                        if highest_label == 3:

                            if ecc[-3] == l - 2:
                                types.append(['D', l])
                                continue  # Dl
                            elif l <= 9 and ecc[-2] == l - 2 and ecc[-5] == l - 3:
                                if l <= 8:  # E_{6,7,8}
                                    types.append(['E', l])
                                    continue
                                else:  # E_{8,1}
                                    types.append(['E', l-1, 1])
                                    continue
                            elif l <= 8 and ecc[0] == l-5 and ecc[-3] == l-3:
                                # TO TEST
                                types.append(['E', l-1, 1])
                                continue
                            else:
                                return None

                        elif highest_label == 4:
                            if ecc[-3] == l - 2 and occur_4 == 1:  # D_n graph
                                # and 1 occurence of label 4
                                if not (3 in vertices_labels[vertices_4[0]] and 3 in vertices_labels[vertices_4[1]]) and 3 not in G0.degree(vertices=vertices_4):
                                    # Edge with label 4 is at the end and not
                                    # on a degree 3 vertex.
                                    types.append(['B', l-1, 1])
                                    continue
                                else:
                                    return None
                            else:
                                return None
                        else:
                            return None

                    elif nb_branching == 2:
                        if highest_label == 3:
                            if ecc[-4] == l - 3:
                                types.append(['D', l-1, 1])
                                continue
                            else:
                                return None

                        else:  # Highest label too high
                            return None

                    else:  # Nb of branching too high
                        return None

                else:  # Branching degree too high
                    return None

    if len(types) == 1:
        types = types[0]
    return CoxeterType(types)

#####################################################################
## Other functions

def check_coxeter_matrix(m):
    """
    Check if ``m`` represents a generalized Coxeter matrix and raise
    and error if not.

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_matrix import check_coxeter_matrix
        sage: m = matrix([[1,3,2],[3,1,-1],[2,-1,1]])
        sage: check_coxeter_matrix(m)

        sage: m = matrix([[1,3],[3,1],[2,-1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: not a square matrix

        sage: m = matrix([[1,3,2],[3,1,-1],[2,-1,2]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: the matrix diagonal is not all 1

        sage: m = matrix([[1,3,3],[3,1,-1],[2,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: the matrix is not symmetric

        sage: m = matrix([[1,3,1/2],[3,1,-1],[1/2,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: invalid Coxeter label 1/2

        sage: m = matrix([[1,3,1],[3,1,-1],[1,-1,1]])
        sage: check_coxeter_matrix(m)
        Traceback (most recent call last):
        ...
        ValueError: invalid Coxeter label 1

    """

    mat = matrix(m)
    if not mat.is_square():
        raise ValueError("not a square matrix")
    for i, row in enumerate(m):
        if mat[i, i] != 1:
            raise ValueError("the matrix diagonal is not all 1")
        for j, val in enumerate(row[i+1:]):
            if val != m[j+i+1][i]:
                raise ValueError("the matrix is not symmetric")
            if val not in ZZ:
                if val > -1 and val in RR and val != infinity:
                    raise ValueError("invalid Coxeter label {}".format(val))
            else:
                if val == 1 or val == 0:
                    raise ValueError("invalid Coxeter label {}".format(val))

def coxeter_matrix_as_function(t):
    """
    Return the Coxeter matrix, as a function

    INPUT:

    - ``t`` -- a Cartan type

    EXAMPLES::

        sage: from sage.combinat.root_system.coxeter_matrix import coxeter_matrix_as_function
        sage: f = coxeter_matrix_as_function(['A',4])
        sage: matrix([[f(i,j) for j in range(1,5)] for i in range(1,5)])
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
    """
    t = CartanType(t)
    m = t.coxeter_matrix()
    index_set = t.index_set()
    reverse = dict((index_set[i], i) for i in range(len(index_set)))
    return lambda i, j: m[i, j]

def coxeter_matrix(t):
    """
    This was deprecated in :trac:`17798` for :class:`CartanMatrix`.

    EXAMPLES::

        sage: coxeter_matrix(['A', 4])
        doctest:...: DeprecationWarning: coxeter_matrix() is deprecated. Use CoxeterMatrix() instead
        See http://trac.sagemath.org/17798 for details.
        [1 3 2 2]
        [3 1 3 2]
        [2 3 1 3]
        [2 2 3 1]
    """
    from sage.misc.superseded import deprecation
    deprecation(17798, 'coxeter_matrix() is deprecated. Use CoxeterMatrix() instead')
    return CoxeterMatrix(t)
