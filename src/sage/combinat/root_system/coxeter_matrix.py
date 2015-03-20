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
from sage.rings.all import ZZ, RR
from sage.rings.infinity import infinity
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.coxeter_type import CoxeterType
from sage.combinat.root_system.root_system import RootSystem
from sage.sets.family import Family

class CoxeterMatrix(Matrix_generic_dense, CoxeterType):
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
        Normalize input so we can inherit from a dense matrix with the 
        appropriate ring.

        .. NOTE::

            To disable the Coxeter type check, use the optional argument
            ``coxeter_type_check = False``.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',1,1])
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

            sage: CoxeterMatrix([[1,-1],[-1,1]]).base_ring()
            Integer Ring
            sage: CoxeterMatrix([[1,-3/2],[-3/2,1]]).base_ring()
            Rational Field
            sage: CoxeterMatrix([[1,-1.5],[-1.5,1]]).base_ring()
            Real Field with 53 bits of precision

        """
        # Special cases with 0 args
        if not args:
            if "coxeter_type" in kwds: # kwds has Coxeter type
                args = ( CoxeterType(kwds["coxeter_type"]), )
            elif "cartan_type" in kwds: # kwds has Cartan type
                args = ( CoxeterType(CartanType(kwds["cartan_type"])), )
            base_ring = ZZ

        if not args:
            data = []
            n = 0
            index_set = tuple()
            coxeter_type = None
            subdivisions = None
            base_ring = ZZ

        elif len(args) == 4 and isinstance(args[0], MatrixSpace): # For pickling
            return typecall(cls, args[0], args[1], args[2], args[3])
        elif isinstance(args[0], CoxeterMatrix):
            return args[0]

        else:
            coxeter_type = None
            subdivisions = None
            index_set = None
            base_ring = ZZ

            from sage.combinat.root_system.cartan_type import CartanType_abstract
            if isinstance(args[0], CartanType_abstract):
                coxeter_type = args[0].coxeter_type()
            elif isinstance(args[0], Graph):
                G = args[0]
                n = G.num_verts()

                # Setup the basis matrix as all 2 except 1 on the diagonal
                data = []
                for i in range(n):
                    data += [[]]
                    for j in range(n):
                        if i == j:
                            data[-1] += [ZZ.one()]
                        else:
                            data[-1] += [2]

                verts = sorted(G.vertices())
                for e in G.edges():
                    m = e[2]
                    if m is None:
                        m = 3
                    elif m == infinity: 
                        m = -1
                    elif m not in ZZ and m > -1:
                        raise ValueError("invalid Coxeter graph label")
                    elif m == 0 or m == 1:
                        raise ValueError("invalid Coxeter graph label")
                    i = verts.index(e[0])
                    j = verts.index(e[1])
                    data[j][i] = data[i][j] = m
                    base_ring = (base_ring.one()*m).parent()
                data = MatrixSpace(base_ring,n)(data)
                index_set = tuple(verts)
                args = [data]
            else:
                try:
                    coxeter_type = CoxeterType(args[0])
                except (TypeError, ValueError, NotImplementedError):
                    pass

            if coxeter_type:
                index_set = coxeter_type.index_set()
                n = len(index_set)
                reverse = {index_set[i]: i for i in range(n)}
                data = [[1 if i == j else 2 for j in range(n)] for i in range(n)]
                for (i,j,l) in coxeter_type.coxeter_graph().edge_iterator():
                    if l == infinity:
                        l = -1
                    data[reverse[i]][reverse[j]] = l
                    data[reverse[j]][reverse[i]] = l
                data = [val for row in data for val in row]
            else:
                M = matrix(args[0])
                base_ring = M.base_ring()
                check_coxeter_matrix(M)
                n = M.ncols()
                if "coxeter_type" in kwds:
                    coxeter_type = CoxeterType(kwds["coxeter_type"])
                elif n == 1:
                    coxeter_type = CoxeterType(['A', 1])
                elif kwds.get("coxeter_type_check", True):
                    coxeter_type = find_coxeter_type_from_matrix(M)
                data = M.list()
                subdivisions = M._subdivisions

            if kwds.get("index_set", None):
                index_set = tuple(kwds["index_set"])

            if len(args) == 1:
                if coxeter_type is not None and index_set is None:
                    index_set = tuple(coxeter_type.index_set())
            elif len(args) == 2:
                index_set = tuple(args[1])
            else:
                raise ValueError("too many arguments")

            if index_set is None:
                index_set = tuple(range(n))

            if len(set(index_set)) != n:
                raise ValueError("the given index set is not valid")

        mat = typecall(cls, MatrixSpace(base_ring, n, sparse=False), data,
                       coxeter_type, index_set)
        mat._subdivisions = subdivisions
        return mat

    def __init__(self, parent, data, coxeter_type, index_set):
        """
        Initialize ``self``.

        TESTS::

            sage: C = CoxeterMatrix(['A', 2, 1])
            sage: TestSuite(C).run(skip=["_test_category", "_test_change_ring"])
        """
        Matrix_generic_dense.__init__(self, parent, data, False, True)

        if min(data)<-1:
            self._is_cyclotomic = False
        else:
            self._is_cyclotomic = True
        self._coxeter_type = coxeter_type
        self._index_set = index_set
        self.set_immutable()

    def __reduce__(self):
        """
        Used for pickling.

        TESTS::

            sage: C = CoxeterMatrix(['A',4])
            sage: M = loads(dumps(C))
            sage: M._index_set
            (1, 2, 3, 4)
        """
        return (CoxeterMatrix, (self.parent(), self.list(),
                                self._coxeter_type, self._index_set))

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
        Return the Cartan type of ``self`` or ``self`` if unknown.

        EXAMPLES::

            sage: C = CoxeterMatrix(['A',4,1])
            sage: C.coxeter_type()
            Coxeter type of ['A', 4, 1]

        If the Cartan type is unknown::

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
        return len(self._index_set)

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

    @cached_method # Transferred the caching here from CoxeterType
    def bilinear_form(self):
        r"""
        Return the bilinear form of ``self``.

        EXAMPLES::
        """

        if self._is_cyclotomic:
            return super(CoxeterMatrix, self).bilinear_form()
        else:
            return super(CoxeterMatrix, self).bilinear_form(R=RR) #FIXME put the right ring

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
        n = self.nrows()
        I = self.index_set()
        val = lambda x: infinity if x == -1 else x
        G = Graph([(I[i], I[j], val(self[i, j]))
                   for i in range(n) for j in range(i)
                   if self[i, j] not in [1, 2]])
        G.add_vertices(I)
        return G.copy(immutable=True)

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
        if self._coxeter_type is not None:
            return self._coxeter_type.is_finite()
        # The type checker should catch all finite types, and checking for
        #   positive definiteness is difficult, so we just assume it is
        #   not finite if the type checker cannot determine the type.
        #return self.bilinear_form().is_positive_definite()
        return False

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
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self._coxeter_type is not None: 
            return self._coxeter_type.is_affine()
        else:
            return False
            # if self.bilinear_form().rank() == self.rank() - 1 and self.bilinear_form().determinant() == 0:
            #     return True
            # else:
            #     return False

#####################################################################
## Type check functions

def find_coxeter_type_from_matrix(coxeter_matrix):
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
        sage:CoxeterMatrix(CoxeterType(['I',100]).coxeter_graph()).coxeter_type()
        Coxeter type of ['I', 100]
        sage: CoxeterMatrix(CoxeterType(['H',3]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 3]
        sage: CoxeterMatrix(CoxeterType(['H',4]).coxeter_graph()).coxeter_type()
        Coxeter type of ['H', 4]

        
    Some affine graphs::
        
        sage: CoxeterMatrix(CoxeterType(['A',3,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['A', 3, 1]
        sage: CoxeterMatrix(CoxeterType(['B',3,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['B', 3, 1]
        sage: CoxeterMatrix(CoxeterType(['C',3,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['C', 3, 1]
        sage: CoxeterMatrix(CoxeterType(['F',4,1]).coxeter_graph()).coxeter_type()
        Coxeter type of ['F', 4, 1]


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

    """
    # First, we build the Coxeter graph of the group without the edge labels
    n = ZZ(coxeter_matrix.nrows())
    G = Graph([range(n), lambda i,j: coxeter_matrix[i, j] not in [1,2]])
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
            types.append(['A',1])
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
                return None #TODO: return hyperbolic type once implemented
        elif l == 3:
            # The `3`-node case. The finite groups to check for
            # here are `A_3`, `B_3` and `H_3`.
            # The affine groups to check for here are `A_{2,1}`, `G_{2,1}`,
            # `B_{2,1}`.
            c0, c1, c2 = comp
            s = sorted([coxeter_matrix[c0, c1],
                        coxeter_matrix[c0, c2],
                        coxeter_matrix[c1, c2]])
            if s[0] == 2: #Have a tree
                if s[1] == 3:
                    if s[2] == 3:
                        types.append(['A',3])
                        continue
                    elif s[2] == 4:
                        types.append(['B',3])
                        continue
                    elif s[2] == 5:
                        types.append(['H',3])
                        continue
                    elif s[2] == 6:
                        types.append(['G',2,1])
                        continue
                    else:
                        return None
                elif s[1] == 4 and s[2] == 4:
                    types.append(['B',2,1])
                    continue
                else:
                    return None
            elif s[0] == 3 and s[1] == 3 and s[2] == 3:
                types.append(['A',2,1])
                continue
            else:
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
                        types.append(['A',4])
                        continue
                    else:
                        types.append(['D',4])
                        continue
                elif s[3:] == [3, 3, 4] or s[3:] == [3, 3, 5] or s[3:] == [3,
                        4, 4]: 
                    if max(G0.degree()) == 3:
                        if s[4:] == [3, 4]:
                            types.append(['B', 3, 1])
                            continue
                        else:
                            return None
                    else: # The graph is a path
                          # Differenciate using sum of edge labels
                        u0 = u[0] + u[1] + u[2]
                        u1 = u[0] + u[3] + u[4]
                        u2 = u[1] + u[3] + u[5]
                        u3 = u[2] + u[4] + u[5]
                        ss = sorted([u0, u1, u2, u3])
                        if s[5] == 4:
                            if ss == [7, 7, 9, 9]:
                                types.append(['F',4])
                                continue
                            elif ss == [7, 8, 8, 9]:
                                types.append(['B',4])
                                continue
                            elif ss == [8, 8, 9, 9]:
                                types.append(['C',3,1])
                                continue
                            else:
                                return None
                        elif ss == [7, 8, 9, 10]:
                            types.append(['H',4])
                            continue
                        else:
                            return None

                else:
                    return None

            elif s == [2, 2, 3, 3, 3, 3] and max(G0.degree()) == 2:
                types.append(['A',3,1])
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
            vertices= G0.vertices()
            sub_cox_matrix = coxeter_matrix.matrix_from_rows_and_columns(vertices,vertices)
            vertices_labels = [set(filter(lambda x: x != 1, r)) for r in
                    sub_cox_matrix.rows()]
            label_list = filter(lambda x: x != 1, sub_cox_matrix.list())
            labels = sorted(set(label_list))
            
            vertices_4 = [index for index in range(l) if 4 in
                    vertices_labels[index]]
            occur_4 = label_list.count(4)/2 # Each label appear twice

            if not G0.is_tree():
                if G0.is_isomorphic(CycleGraph(l)): # Type `A_{l-1,1}`
                    types.append(['A',l-1,1])
                    continue
                else:
                    return None
            
            elif max(degrees) == 2: # The component is a path
                
                
                if labels[-1] == 3: # Highest label is 3
                    types.append(['A',l])
                    continue
                elif labels[-1] == 4: # Highest label is 4
                    
                    if occur_4 == 1: # There is 1 edge with label 4
                        if not (3 in vertices_labels[vertices_4[0]] and 3 in
                        vertices_labels[vertices_4[1]]): # The edge is at the end of the path
                            types.append(['B',l])
                            continue
                        elif l == 5:
                            types.append(['F',4,1])
                            continue
                        else:
                            return None
                    elif occur_4 == 2: # There are 2 edges labeled 4
                        if len(filter(lambda x: 2 in x and 3 not in x and 4 in
                            x, [vertices_labels[i] for i in vertices_4])) == 2:
                            # The edges with 4 are at the ends of the path
                            types.append(['C',l-1,1])
                            continue
                    else:
                        return None
                else:
                    return None


            else: # The graph contains branching vertices

                #Finite D_n,E_6,E_7,E_8
                #Affine B_n,D_n, E_6,E_7,E_8

                if max(degrees) == 3: # Branching degree is 3

                    highest_label = labels[-1]

                    nb_branching = degrees.count(3)
                    ecc = sorted(G0.eccentricity())
                    
                    if nb_branching == 1:

                        if highest_label == 3:

                            if ecc[-3] == l - 2:
                                types.append(['D', l])
                                continue  # Dl
                            elif l <= 9 and ecc[-2] == l - 2 and ecc[-5] == l - 3:
                                if l <= 8: #E_{6,7,8}
                                    types.append(['E', l])
                                    continue
                                else: #E_{8,1}
                                    types.append(['E',l-1,1])
                                    continue
                            elif l <= 8 and ecc[0] == l-5 and ecc[-3] == l-3:
                                #TO TEST
                                types.append(['E',l-1,1])
                                continue
                            else:
                                return None

                        elif highest_label == 4:
                            if ecc[-3] == l - 2 and occur_4 == 1: # D_n graph
                                # and 1 occurence of label 4
                                if not (3 in vertices_labels[vertices_4[0]] and 3 in vertices_labels[vertices_4[1]]) and 3 not in G0.degree(vertices=vertices_4):
                                    # Edge with label 4 is at the end and not
                                    # on a degree 3 vertex.
                                    types.append(['B',l-1,1])
                                    continue
                        else:
                            return None

                    elif nb_branching == 2:
                        if highest_label == 3:
                            if ecc[-4] == l - 3:
                                types.append(['D',l-1,1])
                                continue
                            else:
                                return None

                        else: # Highest label too high
                            return None

                    else: # Nb of branching too high
                        return None

                else: # Branching degree too high
                    return None


            # # Checking that the Coxeter matrix of the subgroup
            # # corresponding to the vertices ``comp`` has all its
            # # off-diagonal entries equal to 2, 3 or at most once 4
            # found_a_4 = False
            # for j in range(l):
            #     for i in range(j):
            #         coxeter_entry = coxeter_matrix[comp[i], comp[j]]
            #         if coxeter_entry in [2, 3]:
            #             continue
            #         if coxeter_entry == 4 and not found_a_4:
            #             found_a_4 = True
            #             continue
            #         return None
            #
            # G0 = G.subgraph(comp)
            # if found_a_4:
            #     # The case when a `4` has been found in the
            #     # Coxeter matrix. This needs only to be checked
            #     # against `B_l`. We use the observation that
            #     # the group is `B_l` if and only if the Coxeter
            #     # graph is an `l`-path (i.e., has diameter
            #     # `l - 1`) and the `4` corresponds to one of
            #     # its two outermost edges.
            #     diameter = G0.diameter()
            #     if diameter != l - 1:
            #         return None
            #
            #     ecc = sorted(((u, v) for (v, u) in G0.eccentricity(with_labels=True).items()))
            #     left_end = ecc[-1][1]
            #     right_end = ecc[-2][1]
            #     left_almost_end = G0.neigbors(left_end)[0]
            #     right_almost_end = G0.neigbors(right_end)[0]
            #     if (coxeter_matrix[left_end, left_almost_end] == 4
            #         or coxeter_matrix[right_end, right_almost_end] == 4):
            #         types.append(['B', l])
            #         continue  # Bl
            #     return None
            #
            # # Now, all off-diagonal entries of the Coxeter matrix
            # # are 2's and 3's. We need to check our group against
            # # `A_l`, `D_l` and `E_l`. Knowing that the Coxeter
            # # graph is a tree, we can use its vertex
            # # eccentricities to check this.
            # ecc = sorted(G0.eccentricity())
            # if ecc[-1] == l - 1:
            #     types.append(['A', l])
            #     continue  # Al
            # if ecc[-3] == l - 2:
            #     types.append(['D', l])
            #     continue  # Dl
            # if l <= 8 and ecc[-2] == l - 2 and ecc[-5] == l - 3:
            #     types.append(['E', l])
            #     continue  # El
            # return None

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

   """
    if not m.is_square():
        raise ValueError("not a square matrix")
    for i, row in enumerate(m.rows()):
        if m[i,i] != 1:
            raise ValueError("the matrix diagonal is not all 1")
        for j, val in enumerate(row[i+1:]):
            if val != m[j+i+1,i]:
                raise ValueError("the matrix is not symmetric")
            if val not in ZZ:
                if val > -1:
                    raise ValueError("invalid Coxeter label {}".format(val))
            else:
                if val == 1 or val == 0:
                    raise ValueError("invalid Coxeter label {}".format(val))

def coxeter_matrix_as_function(t):
    """
    Returns the Coxeter matrix, as a function

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
    return lambda i,j: m[reverse[i], reverse[j]]

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

