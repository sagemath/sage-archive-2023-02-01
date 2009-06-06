r"""
Finite simplicial complexes

AUTHORS:

- John H. Palmieri (2009-04)

This module implements the basic structure of finite simplicial
complexes. Given a set `V` of "vertices", a simplicial complex on `V`
is a collection `K` of subsets of `V` satisfying the condition that if
`S` is one of the subsets in `K`, then so is every subset of `S`.  The
subsets `S` are called the 'simplices' of `K`.

A simplicial complex `K` can be viewed as a purely combinatorial
object, as described above, but it also gives rise to a topological
space `|K|` (its *geometric realization*) as follows: first, the
points of `V` should be in general position in euclidean space.  Next,
if `\{v\}` is in `K`, then the vertex `v` is in `|K|`.  If `\{v, w\}`
is in `K`, then the line segment from `v` to `w` is in `|K|` If `\{u,
v, w\}` is in `K`, then the triangle with vertices `u`, `v`, and `w`
is in `|K|`.  In general, `|K|` is the union of the convex hulls of
simplices of `K`.  Frequently, one abuses notation and uses `K` to
denote both the simplicial complex and the associated topological
space.

For any simplicial complex `K` and any commutative ring `R` there is
an associated chain complex, with differential of degree `-1`.  The
`n^{th}` term is the free `R`-module with basis given by the
`n`-simplices of `K`.  The differential is determined by its value on
any simplex: on the `n`-simplex with vertices `(v_0, v_1, ..., v_n)`,
the differential is the alternating sum with `i^{th}` summand `(-1)^i`
multiplied by the `(n-1)`-simplex obtained by omitting vertex `v_i`.

In the implementation here, the vertex set must be finite. To define a
simplicial complex, specify its vertex set: this should be a list,
tuple, or set, or it can be a non-negative integer `n`, in which case
the vertex set is `(0, ..., n)`.  Also specify the facets: the maximal
faces.

.. note::

   The elements of the vertex set are not automatically contained in
   the simplicial complex: each one is only included if and only if it
   is a vertex of at least one of the specified facets.

EXAMPLES::

    sage: SimplicialComplex([1, 3, 7], [[1], [3, 7]])
    Simplicial complex with vertex set (1, 3, 7) and facets {(3, 7), (1,)}
    sage: SimplicialComplex(2)   # the empty simplicial complex
    Simplicial complex with vertex set (0, 1, 2) and facets {()}
    sage: X = SimplicialComplex(3, [[0,1], [1,2], [2,3], [3,0]])
    sage: X
    Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2), (2, 3), (0, 3), (0, 1)}
    sage: X.stanley_reisner_ring()
    Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)
    sage: X.is_pure()
    True

Sage can perform a number of operations on simplicial complexes, such
as the join and the product, and it can also compute homology::

    sage: S = SimplicialComplex(3, [[0,1], [1,2], [0,2]]) # circle
    sage: T = S.product(S)  # torus
    sage: T
    Simplicial complex with 16 vertices and 18 facets
    sage: T.homology()   # this computes reduced homology
    {0: 0, 1: Z x Z, 2: Z}
    sage: T.euler_characteristic()
    0

Sage knows about some basic combinatorial data associated to a
simplicial complex::

    sage: X = SimplicialComplex(3, [[0,1], [1,2], [2,3], [0,3]])
    sage: X.f_vector()
    [1, 4, 4]
    sage: X.face_poset()
    Finite poset containing 8 elements
    sage: X.stanley_reisner_ring()
    Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)
"""

# possible future directions for SimplicialComplex:
#
#  make compatible with GAP (see http://linalg.org/gap.html)
#  compare to and make compatible with polymake
#  see Macaulay: http://www.math.uiuc.edu/Macaulay2/doc/Macaulay2-1.1/share/doc/Macaulay2/SimplicialComplexes/html/___Simplicial__Complex.html; compare performance and make compatible
#  define morphisms: they are determined by where each vertex goes
#  should + have any meaning?
#  cohomology: compute cup products (and Massey products?)

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import normalize_names
from sage.misc.latex import latex
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex, HomologyGroup
from sage.graphs.graph import Graph
from sage.misc.flatten import flatten
from sage.rings.all import GF

def lattice_paths(t1,t2):
    """
    Given lists (or tuples or ...) ``t1`` and ``t2``, think of them as
    labelings for vertices: ``t1`` labeling points on the x-axis,
    ``t2`` labeling points on the y-axis, both increasing.  Return the
    list of rectilinear paths along the grid defined by these points
    in the plane, starting from ``(t1[0], t2[0])``, ending at
    ``(t1[last], t2[last])``, and at each grid point, going either
    right or up.  See the examples.

    INPUT:

    -  ``t1``, ``t2`` - tuples, lists, other iterables.

    OUTPUT: list of lists of vertices making up the paths as described
    above

    This is used when triangulating the product of simplices.

    EXAMPLES::

        sage: from sage.homology.simplicial_complex import lattice_paths
        sage: lattice_paths([0,1,2], [0,1,2])
        [[(0, 0), (0, 1), (0, 2), (1, 2), (2, 2)],
         [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)],
         [(0, 0), (1, 0), (1, 1), (1, 2), (2, 2)],
         [(0, 0), (0, 1), (1, 1), (2, 1), (2, 2)],
         [(0, 0), (1, 0), (1, 1), (2, 1), (2, 2)],
         [(0, 0), (1, 0), (2, 0), (2, 1), (2, 2)]]
        sage: lattice_paths(('a', 'b', 'c'), (0, 3, 5))
        [[('a', 0), ('a', 3), ('a', 5), ('b', 5), ('c', 5)],
         [('a', 0), ('a', 3), ('b', 3), ('b', 5), ('c', 5)],
         [('a', 0), ('b', 0), ('b', 3), ('b', 5), ('c', 5)],
         [('a', 0), ('a', 3), ('b', 3), ('c', 3), ('c', 5)],
         [('a', 0), ('b', 0), ('b', 3), ('c', 3), ('c', 5)],
         [('a', 0), ('b', 0), ('c', 0), ('c', 3), ('c', 5)]]
    """
    if len(t1) == 0 or len(t2) == 0:
        return [[]]
    elif len(t1) == 1:
        return [[(t1[0], w) for w in t2]]
    elif len(t2) == 1:
        return [[(v, t2[0]) for v in t1]]
    else:
        return ([path + [(t1[-1], t2[-1])] for path
                in lattice_paths(t1[:-1], t2)] +
                [path + [(t1[-1], t2[-1])] for path
                in lattice_paths(t1, t2[:-1])])

class Simplex(SageObject):
    """
    Define a simplex.

    Topologically, a simplex is the convex hull of a collection of
    vertices in general position.  Combinatorially, it is defined just
    by specifying a set of vertices.  It is represented in Sage by the
    tuple of the vertices.

    INPUT:

    -  ``X`` - set of vertices

    OUTPUT: simplex with those vertices

    ``X`` may be a non-negative integer `n`, in which case the
    simplicial complex will have `n+1` vertices `(0, 1, ..., n)`, or
    it may be anything which may be converted to a tuple, in which
    case the vertices will be that tuple.

    .. warning::

       The vertices should be distinct, and no error checking is done
       to make sure this is the case.

    EXAMPLES::

        sage: from sage.homology.simplicial_complex import Simplex
        sage: Simplex(4)
        (0, 1, 2, 3, 4)
        sage: Simplex([3, 4, 1])
        (3, 4, 1)
        sage: X = Simplex((3, 'a', 'vertex')); X
        (3, 'a', 'vertex')
        sage: X == loads(dumps(X))
        True
    """

    def __init__(self, X):
        """
        Define a simplex.  See ``Simplex`` for full documentation.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(2)
            (0, 1, 2)
            sage: Simplex(('a', 'b', 'c'))
            ('a', 'b', 'c')
        """
        try:
            N = int(X) + 1
            self.__tuple = tuple(range(N))
        except TypeError:
            self.__tuple = tuple(X)
        self.__set = frozenset(self.__tuple)

    def tuple(self):
        """
        The tuple attached to this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(3).tuple()
            (0, 1, 2, 3)

        Although simplices are printed as if they were tuples, they
        are not the same type::

            sage: type(Simplex(3).tuple())
            <type 'tuple'>
            sage: type(Simplex(3))
            <class 'sage.homology.simplicial_complex.Simplex'>
        """
        return self.__tuple

    def set(self):
        """
        The frozenset attached to this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(3).set()
            frozenset([0, 1, 2, 3])
        """
        return self.__set

    def is_face(self, other):
        """
        Return True iff this simplex is a face of other.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(3).is_face(Simplex(5))
            True
            sage: Simplex(5).is_face(Simplex(2))
            False
            sage: Simplex(['a', 'b', 'c']).is_face(Simplex(8))
            False
        """
        return self.__set.issubset(other.__set)

    def __contains__(self, x):
        """
        Return True iff ``x`` is a vertex of this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: 3 in Simplex(5)
            True
            sage: 3 in Simplex(2)
            False
        """
        return self.__set.__contains__(x)

    def __getitem__(self, n):
        """
        Return the nth vertex in this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(5)[2]
            2
            sage: Simplex(['a', 'b', 'c'])[1]
            'b'
        """
        return self.__tuple.__getitem__(n)

    def __iter__(self):
        """
        Iterator for the vertices of this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: [v**2 for v in Simplex(3)]
            [0, 1, 4, 9]
        """
        return self.__tuple.__iter__()

    def __add__(self, other):
        """
        Simplex obtained by concatenating the underlying tuples of the
        two arguments.

        INPUT:

        -  ``other`` - another simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex((1,2,3)) + Simplex((5,6))
            (1, 2, 3, 5, 6)
        """
        return Simplex(self.__tuple.__add__(other.__tuple))

    def face(self, n):
        """
        The nth face of this simplex.

        INPUT:

        -  ``n`` - an integer between 0 and the dimension of this simplex

        OUTPUT: the simplex obtained by removing the nth vertex from
        this simplex

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: S = Simplex(4)
            sage: S.face(0)
            (1, 2, 3, 4)
            sage: S.face(3)
            (0, 1, 2, 4)
        """
        if n >= 0 and n <= self.dimension():
            return Simplex(self.__tuple[:n] + self.__tuple[n+1:])
        else:
            raise IndexError, "%s does not have an nth face for n=%s." % (self, n)

    def faces(self):
        """
        The list of faces (of codimension 1) of this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: S = Simplex(4)
            sage: S.faces()
            [(1, 2, 3, 4), (0, 2, 3, 4), (0, 1, 3, 4), (0, 1, 2, 4), (0, 1, 2, 3)]
            sage: len(Simplex(10).faces())
            11
        """
        return [self.face(i) for i in range(self.dimension()+1)]

    def dimension(self):
        """
        The dimension of this simplex: the number of vertices minus 1.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(5).dimension() == 5
            True
            sage: Simplex(5).face(1).dimension()
            4
        """
        return len(self.__tuple) - 1

    def is_empty(self):
        """
        Return True iff this simplex is the empty simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: [Simplex(n).is_empty() for n in range(-1,4)]
            [True, False, False, False, False]
        """
        return self.dimension() < 0

    def join(self, right, rename_vertices=True):
        """
        The join of this simplex with another one.

        The join of two simplices `[v_0, ..., v_k]` and `[w_0, ...,
        w_n]` is the simplex `[v_0, ..., v_k, w_0, ..., w_n]`.

        INPUT:

        -  ``right`` - the other simplex (the right-hand factor)

        -  ``rename_vertices`` -- boolean (optional, default True).  If
           this is True, the vertices in the join will be renamed by
           this formula: vertex "v" in the left-hand factor --> vertex
           "Lv" in the join, vertex "w" in the right-hand factor -->
           vertex "Rw" in the join.  If this is false, this tries to
           construct the join without renaming the vertices; this may
           cause problems if the two factors have any vertices with
           names in common.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(2).join(Simplex(3))
            ('L0', 'L1', 'L2', 'R0', 'R1', 'R2', 'R3')
            sage: Simplex(['a', 'b']).join(Simplex(['x', 'y', 'z']))
            ('La', 'Lb', 'Rx', 'Ry', 'Rz')
            sage: Simplex(['a', 'b']).join(Simplex(['x', 'y', 'z']), rename_vertices=False)
            ('a', 'b', 'x', 'y', 'z')
        """
        if rename_vertices:
            vertex_set = (["L" + str(v) for v in self]
                          + ["R" + str(w) for w in right])
        else:
            vertex_set = self.__tuple + right.__tuple
        return Simplex(vertex_set)

    def product(self, other, rename_vertices=False):
        """
        The product of this simplex with another one, as a list of simplices.

        INPUT:

        -  ``other`` - the other simplex

        - ``rename_vertices`` -- boolean (optional, default True). If
           this is False, then the vertices in the product are the set
           of ordered pairs `(v,w)` where `v` is a vertex in the
           left-hand factor (``self``) and `w` is a vertex in the
           right-hand factor (``other``). If this is True, then the
           vertices are renamed as "LvRw" (e.g., the vertex (1,2)
           would become "L1R2").  This is useful if you want to define
           the Stanley-Reisner ring of the complex: vertex names like
           (0,1) are not suitable for that, while vertex names like
           "L0R1" are.

        Algorithm: see Hatcher, p. 277-278 (who in turn refers to
        Eilenberg-Steenrod, p. 68): given ``Simplex(m)`` and
        ``Simplex(n)``, then ``Simplex(m)`` x ``Simplex(n)`` can be
        triangulated as follows: for each path `f` from `(0,0)` to
        `(m,n)` along the integer grid in the plane, going up or right
        at each lattice point, associate an `(m+n)`-simplex with
        vertices `v_0`, `v_1`, ..., where `v_k` is the `k^{th}` vertex
        in the path `f`.

        Note that there are `m+n` choose `n` such paths.  Note also
        that each vertex in the product is a pair of vertices `(v,w)`
        where `v` is a vertex in the left-hand factor and `w`
        is a vertex in the right-hand factor.

        .. note::

           This produces a list of simplices -- not a ``Simplex``, not
           a ``SimplicialComplex``.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: len(Simplex(2).product(Simplex(2)))
            6
            sage: Simplex(1).product(Simplex(1))
            [((0, 0), (0, 1), (1, 1)), ((0, 0), (1, 0), (1, 1))]

        REFERENCES:

        - A. Hatcher, "Algebraic Topology", Cambridge University Press
          (2002).
        """
        if not rename_vertices:
            return [Simplex(x) for x in lattice_paths(self.tuple(), other.tuple())]

        answer = []
        for x in lattice_paths(self.tuple(), other.tuple()):
            new = tuple(["L" + str(v) + "R" + str(w) for (v,w) in x])
            answer.append(Simplex(new))
        return answer


    def __cmp__(self, other):
        """
        Return True iff this simplex is the same as ``other``: that
        is, if the vertices of the two are the same, even with a
        different ordering

        INPUT:

        -  ``other`` - the other simplex

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex([0,1,2]) == Simplex([0,2,1])
            True
            sage: Simplex([0,1,2]) == Simplex(['a','b','c'])
            False
        """
        if self.__set == other.__set:
            return 0
        else:
            return -1

    def __hash__(self):
        """
        Hash value for this simplex.  This computes the hash value of
        the Python frozenset of the underlying tuple, since this is
        what's important when testing equality.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex([1,2,0]).__hash__() == Simplex(2).__hash__()
            True
            sage: Simplex([1,2,0,1,1,2]).__hash__() == Simplex(2).__hash__()
            True
        """
        return hash(self.__set)

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: S = Simplex(5)
            sage: S._repr_()
            '(0, 1, 2, 3, 4, 5)'
        """
        return self.__tuple.__repr__()

    def _latex_(self):
        """
        LaTeX representation.

        EXAMPLES::

            sage: from sage.homology.simplicial_complex import Simplex
            sage: Simplex(3)._latex_()
            \left(0,
            1,
            2,
            3\right)
        """
        return latex(self.__tuple)


class SimplicialComplex(SageObject):
    """
    Define a simplicial complex.

    INPUT:

    -  ``vertex_set`` - set of vertices

    -  ``maximal_faces`` - set of maximal faces

    -  ``vertex_check`` - boolean (optional, default True)

    -  ``maximality_check`` - boolean (optional, default True)

    -  ``sort_facets`` - boolean (optional, default True)

    -  ``name_check`` - boolean (optional, default False)

    OUTPUT: a simplicial complex

    ``vertex_set`` may be a non-negative integer `n` (in which case
    the simplicial complex will have `n+1` vertices `\{0, 1, ...,
    n\}`), or it may be anything which may be converted to a tuple.
    Call the elements of this 'vertices'.

    ``maximal_faces`` should be a list or tuple or set (indeed,
    anything which may be converted to a set) whose elements are lists
    (or tuples, etc.) of vertices.

    If ``vertex_check`` is True, check to see that each given maximal
    face is a subset of the vertex set. Raise an error for any bad
    face.

    If ``maximality_check`` is True, check that each maximal face is,
    in fact, maximal. In this case, when producing the internal
    representation of the simplicial complex, omit those that are not.
    It is highly recommended that this be True; various methods for
    this class may fail if faces which are claimed to be maximal are
    in fact not.

    If ``sort_facets`` is True, sort the vertices in each facet.  If
    the vertices in different facets are not ordered compatibly (e.g.,
    if you have facets (1, 3, 5) and (5, 3, 8)), then homology
    calculations may have unpredictable results.

    If ``name_check`` is True, check the names of the vertices to see
    if they can be easily converted to generators of a polynomial ring
    -- use this if you plan to use the Stanley-Reisner ring for the
    simplicial complex.

    .. note::

       The elements of ``vertex_set`` are not automatically in the
       simplicial complex: each one is only included if it is a vertex
       of at least one of the specified facets.

    EXAMPLES::

        sage: SimplicialComplex(4, [[1,2], [1,4]])
        Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(1, 2), (1, 4)}
        sage: SimplicialComplex(3, [[0,2], [0,3], [0]])
        Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2), (0, 3)}
        sage: SimplicialComplex(3, [[0,2], [0,3], [0]], maximality_check=False)
        Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2), (0, 3), (0,)}
        sage: S = SimplicialComplex(['a', 'b', 'c'], (('a', 'b'), ('a', 'c'), ('b', 'c')))
        sage: S
        Simplicial complex with vertex set ('a', 'b', 'c') and facets {('b', 'c'), ('a', 'c'), ('a', 'b')}
        sage: S == loads(dumps(S))
        True
        """
    def __init__(self, vertex_set=[], maximal_faces=[], **kwds):
        """
        Define a simplicial complex.  See ``SimplicialComplex`` for more
        documentation.

        EXAMPLES::

            sage: SimplicialComplex(3, [[0,2], [0,3], [0]])
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2), (0, 3)}
            sage: SimplicialComplex(['a', 'b', 'c'], (('a', 'b'), ('a', 'c'), ('b', 'c')))
            Simplicial complex with vertex set ('a', 'b', 'c') and facets {('b', 'c'), ('a', 'c'), ('a', 'b')}
        """
        # process kwds
        if 'sort_facets' in kwds:
            sort_facets = kwds['sort_facets']
        else:
            sort_facets = True
        if 'vertex_check' in kwds:
            vertex_check = kwds['vertex_check']
        else:
            vertex_check = True
        if 'maximality_check' in kwds:
            maximality_check = kwds['maximality_check']
        else:
            maximality_check = True
        if 'name_check' in kwds:
            name_check = kwds['name_check']
        else:
            name_check = False
        # done with kwds
        if sort_facets:
            try:  # vertex_set is an iterable
                vertices = Simplex(sorted(vertex_set))
            except TypeError:  # vertex_set is an integer
                vertices = Simplex(vertex_set)
        else:
            vertices = Simplex(vertex_set)
        gen_dict = {}
        for v in vertices:
            if name_check:
                try:
                    if int(v) < 0:
                        raise ValueError, "The vertex %s does not have an appropriate name."%v
                except ValueError:  # v is not an integer
                    try:
                        normalize_names(1, v)
                    except ValueError:
                        raise ValueError, "The vertex %s does not have an appropriate name."%v
            # build dictionary of generator names
            try:
                gen_dict[v] = 'x%s'%int(v)
            except:
                gen_dict[v] = v
        # build set of facets
        good_faces = []
        maximal_simplices = [Simplex(f) for f in maximal_faces]
        for face in maximal_simplices:
            # check whether vertices of each face are contained in vertex set
            if vertex_check:
                if not face.is_face(vertices):
                    raise ValueError, "The face %s is not a subset of the vertex set." % face
            # check whether each given face is actually maximal
            face_is_maximal = True
            if maximality_check:
                faces_to_be_removed = []
                for other in good_faces:
                    if other.is_face(face):
                        faces_to_be_removed.append(other)
                    elif face_is_maximal:
                        face_is_maximal = not face.is_face(other)
                for x in faces_to_be_removed:
                    good_faces.remove(x)
            if sort_facets:
                face = Simplex(sorted(face.tuple()))
            if face_is_maximal:
                good_faces += [face]
        # if no maximal faces, add the empty face as a facet
        if len(maximal_simplices) == 0:
            good_faces.append(Simplex(-1))
        # now record the attributes for self
        # self._vertex_set: the Simplex formed by the vertices
        self._vertex_set = vertices
        # self._facets: list of facets
        self._facets = good_faces
        # self._faces: dictionary of dictionaries of faces.  The main
        # dictionary is keyed by subcomplexes, and each value is a
        # dictionary keyed by dimension.  This should be empty until
        # needed -- that is, until the faces method is called
        self._faces = {}
        # self._gen_dict: dictionary of names for the polynomial
        # generators of the Stanley-Reisner ring
        self._gen_dict = gen_dict
        # self._complex: dictionary indexed by dimension d, base_ring,
        # etc.: differential from dim d to dim d-1 in the associated
        # chain complex.  thus to get the differential in the cochain
        # complex from dim d-1 to dim d, take the transpose of this
        # one.
        self._complex = {}
        # self.__contractible: if not None, a contractible subcomplex
        # of self, as found by the _contractible_subcomplex method.
        self.__contractible = None
        # self.__enlarged: dictionary of enlarged subcomplexes,
        # indexed by subcomplexes.  For use in the _enlarge_subcomplex
        # method.
        self.__enlarged = {}

    def __cmp__(self,right):
        """
        Two simplicial complexes are equal iff their vertex sets are
        equal and their sets of facets are equal.

        EXAMPLES::

            sage: SimplicialComplex(4, [[1,2], [2,3], [4]]) == SimplicialComplex(4, [[4], [2,3], [3], [2,1]])
            True
            sage: X = SimplicialComplex(4)
            sage: X.add_face([1,3])
            sage: X == SimplicialComplex(4, [[1,3]])
            True
        """
        if (self.vertices() == right.vertices() and
            set(self._facets) == set(right._facets)):
            return 0
        else:
            return -1

    def vertices(self):
        """
        The vertex set of this simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex(15, [[0,1], [1,2]])
            sage: S
            Simplicial complex with 16 vertices and facets {(1, 2), (0, 1)}
            sage: S.vertices()
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)

        Note that this actually returns a simplex::

            sage: type(S.vertices())
            <class 'sage.homology.simplicial_complex.Simplex'>
        """
        return self._vertex_set

    def maximal_faces(self):
        """
        The maximal faces (a.k.a. facets) of this simplicial complex.

        This just returns the set of facets used in defining the
        simplicial complex, so if the simplicial complex was defined
        with no maximality checking, none is done here, either.

        EXAMPLES::

            sage: Y = SimplicialComplex(5, [[0,2], [1,4]])
            sage: Y.maximal_faces()
            {(1, 4), (0, 2)}

        ``facets`` is a synonym for ``maximal_faces``::

            sage: S = SimplicialComplex(2, [[0,1], [0,1,2]])
            sage: S.facets()
            {(0, 1, 2)}
        """
        return Set(self._facets)

    facets = maximal_faces

    def faces(self, subcomplex=None):
        """
        The faces of this simplicial complex, in the form of a
        dictionary of sets keyed by dimension.  If the optional
        argument ``subcomplex`` is present, then return only the
        faces which are *not* in the subcomplex.

        INPUT:

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default None).  Return faces which are not in
           this subcomplex.

        EXAMPLES::

            sage: Y = SimplicialComplex(5, [[1,2], [1,4]])
            sage: Y.faces()
            {0: set([(4,), (2,), (1,)]), 1: set([(1, 2), (1, 4)]), -1: set([()])}
            sage: L = SimplicialComplex(5, [[1,2]])
            sage: Y.faces(subcomplex=L)
            {0: set([(4,)]), 1: set([(1, 4)]), -1: set([])}
        """
        if subcomplex not in self._faces:
            # Faces is the dictionary of faces in self but not in
            # subcomplex, indexed by dimension
            Faces = {}
            # sub_facets is the dictionary of facets in the subcomplex
            sub_facets = {}
            for i in range(-1,self.dimension()+1):
                Faces[i] = set([])
                sub_facets[i] = set([])
            for f in self._facets:
                dim = f.dimension()
                Faces[dim].add(f)
            if subcomplex is not None:
                for g in subcomplex._facets:
                    dim = g.dimension()
                    Faces[dim].discard(g)
                    sub_facets[dim].add(g)
            # bad_faces is the set of faces in the subcomplex in the
            # current dimension
            bad_faces = sub_facets[self.dimension()]
            for dim in range(self.dimension(), -1, -1):
                # bad_bdries = boundaries of bad_faces: things to be
                # discarded in dim-1
                bad_bdries = sub_facets[dim-1]
                for f in bad_faces:
                    bad_bdries.update(f.faces())
                for f in Faces[dim]:
                    Faces[dim-1].update(set(f.faces()).difference(bad_bdries))
                bad_faces = bad_bdries
            self._faces[subcomplex] = Faces
        return self._faces[subcomplex]

    def _flattened_faces(self):
        """
        The faces of this simplicial complex, as a list.  (Flattened
        from the dictionary returned by the ``faces`` method.)

        EXAMPLES::

            sage: Y = SimplicialComplex(5, [[1,2], [1,4]])
            sage: Y.faces()
            {0: set([(4,), (2,), (1,)]), 1: set([(1, 2), (1, 4)]), -1: set([()])}
            sage: Y._flattened_faces()
            [(4,), (2,), (1,), (1, 2), (1, 4), ()]
        """
        return flatten([list(f) for f in self.faces().values()])

    def n_faces(self, n, subcomplex=None):
        """
        The set of faces of dimension n of this simplicial complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional faces which are *not* in the
        subcomplex.

        INPUT:

        -  ``n`` - non-negative integer

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default None).  Return ``n``-dimensional faces
           which are not in this subcomplex.

        EXAMPLES::

            sage: S = Set(range(1,5))
            sage: Z = SimplicialComplex(S, S.subsets())
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.n_faces(2)
            set([(1, 3, 4), (1, 2, 3), (2, 3, 4), (1, 2, 4)])
            sage: K = SimplicialComplex(S, [[1,2,3], [2,3,4]])
            sage: Z.n_faces(2, subcomplex=K)
            set([(1, 3, 4), (1, 2, 4)])
            """
        if n in self.faces(subcomplex):
            return self.faces(subcomplex)[n]
        else:
            return set([])

    def f_vector(self):
        """
        The `f`-vector of this simplicial complex: a list whose
        `n^{th}` item is the number of `(n-1)`-faces.  Note that, like
        all lists in Sage, this is indexed starting at 0: the 0th
        element in this list is the number of -1 faces.

        EXAMPLES::

            sage: S = Set(range(1,5))
            sage: Z = SimplicialComplex(S, S.subsets()); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.f_vector()
            [1, 4, 6, 4, 1]
            sage: Y = SimplicialComplex(5, [[1,2], [1,4]])
            sage: Y.f_vector()[2]
            2
        """
        return [self._f_dict()[n] for n in range(-1, self.dimension()+1)]

    def _f_dict(self):
        """
        The `f`-vector of this simplicial complex as a dictionary: the
        item associated to an integer `n` is the number of the
        `n`-faces.

        EXAMPLES::

            sage: S = Set(range(1,5))
            sage: Z = SimplicialComplex(S, S.subsets()); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: [Z._f_dict()[n] for n in range(-1, 4)]
            [1, 4, 6, 4, 1]
            sage: Y = SimplicialComplex(5, [[1,2], [1,4]])
            sage: Y._f_dict()[1]
            2
        """
        answer = {}
        answer[-1] = 1
        for n in range(self.dimension() + 1):
            answer[n] = len(self.n_faces(n))
        return answer

    def euler_characteristic(self):
        r"""
        The Euler characteristic of this simplicial complex: the
        alternating sum over `n \geq 0` of the number of
        `n`-simplices.

        EXAMPLES::

            sage: Y = SimplicialComplex(5, [[1,2], [1,4]])
            sage: Y.euler_characteristic()
            1
            sage: X = SimplicialComplex(3, [[0,1], [0,2], [1,2]])
            sage: X.euler_characteristic()
            0
        """
        return sum([(-1)**n * len(self.n_faces(n)) for n in range(self.dimension() + 1)])

    def dimension(self):
        """
        The dimension of this simplicial complex: the maximum
        dimension of its faces.

        EXAMPLES::

            sage: U = SimplicialComplex(5, [[1,2], [1, 3, 4]])
            sage: U.dimension()
            2
            sage: X = SimplicialComplex(3, [[0,1], [0,2], [1,2]])
            sage: X.dimension()
            1
        """
        return max([face.dimension() for face in self._facets])

    def is_pure(self):
        """
        True iff this simplicial complex is pure: a simplicial complex
        is pure iff all of its maximal faces have the same dimension.

        .. warning::

           This may give the wrong answer if the simplicial complex
           was constructed with ``maximality_check`` set to False.

        EXAMPLES::

            sage: U = SimplicialComplex(5, [[1,2], [1, 3, 4]])
            sage: U.is_pure()
            False
            sage: X = SimplicialComplex(3, [[0,1], [0,2], [1,2]])
            sage: X.is_pure()
            True
        """
        dims = [face.dimension() for face in self._facets]
        return max(dims) == min(dims)

    def product(self, right, rename_vertices=True):
        """
        The product of this simplicial complex with another one.

        INPUT:

        -  ``right`` - the other simplicial complex (the right-hand
           factor)

        - ``rename_vertices`` -- boolean (optional, default True). If
           this is False, then the vertices in the product are the set
           of ordered pairs `(v,w)` where `v` is a vertex in ``self``
           and `w` is a vertex in ``right``. If this is True, then the
           vertices are renamed as "LvRw" (e.g., the vertex (1,2)
           would become "L1R2").  This is useful if you want to define
           the Stanley-Reisner ring of the complex: vertex names like
           (0,1) are not suitable for that, while vertex names like
           "L0R1" are.

        The vertices in the product will be the set of ordered pairs
        `(v,w)` where `v` is a vertex in self and `w` is a vertex in
        right.

        .. warning::

           If ``X`` and ``Y`` are simplicial complexes, then ``X*Y``
           returns their join, not their product.

        EXAMPLES::

            sage: S = SimplicialComplex(3, [[0,1], [1,2], [0,2]]) # circle
            sage: K = SimplicialComplex(1, [[0,1]])   # edge
            sage: S.product(K).vertices()  # cylinder
            ('L0R0', 'L0R1', 'L1R0', 'L1R1', 'L2R0', 'L2R1', 'L3R0', 'L3R1')
            sage: S.product(K, rename_vertices=False).vertices()
            ((0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1))
            sage: T = S.product(S)  # torus
            sage: T
            Simplicial complex with 16 vertices and 18 facets
            sage: T.homology()
            {0: 0, 1: Z x Z, 2: Z}

        These can get large pretty quickly::

            sage: T = simplicial_complexes.Torus(); T
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 14 facets
            sage: K = simplicial_complexes.KleinBottle(); K
            Simplicial complex with 9 vertices and 18 facets
            sage: T.product(K)      # long time: 5 or 6 seconds
            Simplicial complex with 63 vertices and 1512 facets
        """
        vertices = []
        for v in self.vertices():
            for w in right.vertices():
                if rename_vertices:
                    vertices.append("L" + str(v) + "R" + str(w))
                else:
                    vertices.append((v,w))
        facets = []
        for f in self._facets:
            for g in right._facets:
                facets.extend(f.product(g, rename_vertices))
        return SimplicialComplex(vertices, facets)

    def join(self, right, rename_vertices=True):
        """
        The join of this simplicial complex with another one.

        The join of two simplicial complexes `S` and `T` is the
        simplicial complex `S*T` with simplices of the form `[v_0,
        ..., v_k, w_0, ..., w_n]` for all simplices `[v_0, ..., v_k]` in
        `S` and `[w_0, ..., w_n]` in `T`.

        INPUT:

        -  ``right`` - the other simplicial complex (the right-hand
           factor)

        -  ``rename_vertices`` -- boolean (optional, default True).  If
           this is True, the vertices in the join will be renamed by
           the formula: vertex "v" in the left-hand factor --> vertex
           "Lv" in the join, vertex "w" in the right-hand factor -->
           vertex "Rw" in the join.  If this is false, this tries to
           construct the join without renaming the vertices; this will
           cause problems if the two factors have any vertices with
           names in common.

        EXAMPLES::

            sage: S = SimplicialComplex(1, [[0], [1]])
            sage: T = SimplicialComplex([2, 3], [[2], [3]])
            sage: S.join(T)
            Simplicial complex with vertex set ('L0', 'L1', 'R2', 'R3') and 4 facets
            sage: S.join(T, rename_vertices=False)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 3), (1, 2), (0, 2), (0, 3)}

        The notation '*' may be used, as well::

            sage: S * S
            Simplicial complex with vertex set ('L0', 'L1', 'R0', 'R1') and 4 facets
            sage: S * S * S * S * S * S * S * S
            Simplicial complex with 16 vertices and 256 facets
        """
        if rename_vertices:
            vertex_set = (["L" + str(v) for v in self.vertices()]
                          + ["R" + str(w) for w in right.vertices()])
        else:
            vertex_set = tuple(self._vertex_set) + tuple(right.vertices())
        facets = []
        for f in self._facets:
            for g in right._facets:
                facets.append(f.join(g, rename_vertices))
        return SimplicialComplex(vertex_set, facets)

    __mul__ = join

    def cone(self):
        """
        The cone on this simplicial complex.

        The cone is the simplicial complex formed by adding a new
        vertex `C` and simplices of the form `[C, v_0, ..., v_k]` for
        every simplex `[v_0, ..., v_k]` in the original simplicial
        complex.  That is, the cone is the join of the original
        complex with a one-point simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex(1, [[0], [1]])
            sage: S.cone()
            Simplicial complex with vertex set ('L0', 'L1', 'R0') and facets {('L0', 'R0'), ('L1', 'R0')}
        """
        return self.join(SimplicialComplex(["0"], [["0"]]),
                         rename_vertices = True)

    def suspension(self, n=1):
        """
        The suspension of this simplicial complex.

        INPUT:

        -  ``n`` - positive integer (optional, default 1): suspend this
           many times.

        The suspension is the simplicial complex formed by adding two
        new vertices `S_0` and `S_1` and simplices of the form `[S_0,
        v_0, ..., v_k]` and `[S_1, v_0, ..., v_k]` for every simplex
        `[v_0, ..., v_k]` in the original simplicial complex.  That
        is, the cone is the join of the original complex with a
        two-point simplicial complex.

        EXAMPLES::

            sage: S = SimplicialComplex(1, [[0], [1]])
            sage: S.suspension()
            Simplicial complex with vertex set ('L0', 'L1', 'R0', 'R1') and 4 facets
            sage: S3 = S.suspension(3)  # the 3-sphere
            sage: S3.homology()
            {0: 0, 1: 0, 2: 0, 3: Z}
        """
        if n<0:
            raise ValueError, "n must be non-negative."
        if n==0:
            return self
        if n==1:
            return self.join(SimplicialComplex(["0", "1"], [["0"], ["1"]]),
                             rename_vertices = True)
        return self.suspension().suspension(int(n-1))

    def chain_complex(self, dimensions=None, base_ring=ZZ, subcomplex=None,
                      augmented=False, cochain=False, verbose=False,
                      check_diffs=False):
        """
        The chain complex associated to this simplicial complex.

        INPUT:

        -  ``dimensions`` - if None, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.

        -  ``base_ring`` - commutative ring (optional, default ZZ)

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default empty).  Compute the chain complex
           relative to this subcomplex.

        -  ``augmented`` - boolean (optional, default False).  If True,
           return the augmented chain complex (that is, include a class
           in dimension `-1` corresponding to the empty cell).  This is
           ignored if ``dimensions`` is specified.

        -  ``cochain`` - boolean (optional, default False).  If True,
           return the cochain complex (that is, the dual of the chain
           complex).

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the chain complex is computed.

        -  ``check_diffs`` - boolean (optional, default False).  If True,
           make sure that the chain complex is actually a chain complex:
           the differentials are composable and their product is zero.

        .. note::

           If subcomplex is nonempty, then the argument ``augmented``
           has no effect: the chain complex relative to a nonempty
           subcomplex is zero in dimension `-1`.

        EXAMPLES::

            sage: circle = SimplicialComplex(2, [[0,1], [1,2], [0, 2]])
            sage: circle.chain_complex()
            Chain complex with at most 2 nonzero terms over Integer Ring.
            sage: circle.chain_complex()._latex_()
            '\\Bold{Z}^{3} \\xrightarrow{d_{1}} \\Bold{Z}^{3}'
            sage: circle.chain_complex(base_ring=QQ, augmented=True)
            Chain complex with at most 3 nonzero terms over Rational Field.
        """
        # initialize subcomplex
        if subcomplex is None:
            subcomplex = SimplicialComplex(self.vertices())
        else:
            # subcomplex is not empty, so don't augment the chain complex
            augmented = False
        # now construct the range of dimensions in which to compute
        if dimensions is None:
            dimensions = range(0, self.dimension()+1)
            first = 0
        else:
            augmented = False
            first = dimensions[0]
        differentials = {}
        # in the chain complex, compute the first dimension by hand,
        # and don't cache it: it may be differ from situation to
        # situation because of boundary effects.
        current = None
        current_dim = None
        if augmented:  # then first == 0
            current = list(self.n_faces(0, subcomplex=subcomplex))
            current_dim = 0
            if cochain:
                differentials[-1] = matrix(base_ring, len(current), 1,
                                           [1]*len(current))
            else:
                differentials[0] = matrix(base_ring, 1, len(current),
                                          [1]*len(current))
        elif first == 0 and not augmented:
            current = list(self.n_faces(0, subcomplex=subcomplex))
            current_dim = 0
            if not cochain:
                differentials[0] = matrix(base_ring, 0, len(current))
        else:  # first > 0
            current = list(self.n_faces(first, subcomplex=subcomplex))
            current_dim = first
            if not cochain:
                differentials[first] = matrix(base_ring, 0, len(current))
        for n in dimensions[1:]:
            if verbose:
                print "  starting dimension %s" % n
            if (n, subcomplex) in self._complex:
                if cochain:
                    differentials[n-1] = self._complex[(n, subcomplex)].transpose().change_ring(base_ring)
                    mat = differentials[n-1]
                else:
                    differentials[n] = self._complex[(n, subcomplex)].change_ring(base_ring)
                    mat = differentials[n]
                if verbose:
                    print "    boundary matrix (cached): it's %s by %s." % (mat.nrows(), mat.ncols())
            else:
                # 'current' is the list of faces in dimension n
                #
                # 'old' is a dictionary, with keys the faces in the
                # previous dimension (dim n-1 for the chain complex,
                # n+1 for the cochain complex), values the integers 0,
                # 1, 2, ... (the index of the face).  finding an entry
                # in a dictionary seems to be faster than finding the
                # index of an entry in a list.
                if current_dim == n-1:
                    old = dict(zip(current, range(len(current))))
                else:
                    set_of_faces = list(self.n_faces(n-1, subcomplex=subcomplex))
                    old = dict(zip(set_of_faces, range(len(set_of_faces))))
                current = list(self.n_faces(n, subcomplex=subcomplex))
                current_dim = n
                # construct matrix.  it is easiest to construct it as
                # a sparse matrix, specifying which entries are
                # nonzero via a dictionary.
                matrix_data = {}
                col = 0
                if len(old) and len(current):
                    for simplex in current:
                        for i in range(n+1):
                            face_i = simplex.face(i)
                            try:
                                matrix_data[(old[face_i], col)] = (-1)**i
                            except KeyError:
                                pass
                        col += 1
                mat = matrix(ZZ, len(old), len(current), matrix_data)
                if cochain:
                    self._complex[(n, subcomplex)] = mat
                    differentials[n-1] = mat.transpose().change_ring(base_ring)
                else:
                    self._complex[(n, subcomplex)] = mat
                    differentials[n] = mat.change_ring(base_ring)
                if verbose:
                    print "    boundary matrix computed: it's %s by %s." % (mat.nrows(), mat.ncols())
        # now for the cochain complex, compute the last dimension by
        # hand, and don't cache it.
        if cochain:
            n = dimensions[-1] + 1
            if current_dim != n-1:
                current = list(self.n_faces(n-1, subcomplex=subcomplex))
            differentials[n-1] = matrix(base_ring, 0, len(current))
        # finally, return the chain complex
        if cochain:
            return ChainComplex(data=differentials, base_ring=base_ring,
                                degree=1, check_products=check_diffs)
        else:
            return ChainComplex(data=differentials, base_ring=base_ring,
                                degree=-1, check_products=check_diffs)

    def homology(self, dim=None, base_ring=ZZ, subcomplex=None,
                 cohomology=False, enlarge=True, algorithm='auto',
                 verbose=False):
        """
        The reduced homology of this simplicial complex.

        INPUT:

        -  ``dim`` - integer or list of integers or None (optional,
           default None).  If None, then return the homology in every
           dimension.  If ``dim`` is an integer or list, return the
           homology in the given dimensions.  (Actually, if ``dim`` is
           a list, return the homology in the range from ``min(dim)``
           to ``max(dim)``.)

        -  ``base_ring`` - commutative ring (optional, default ZZ).
           Must be ZZ or a field.

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default None).  Compute homology relative to
           this subcomplex.

        -  ``cohomology`` - boolean (optional, default False).  If True,
           compute cohomology rather than homology.

        -  ``enlarge`` - boolean (optional, default True).  If True,
           find a new subcomplex homotopy equivalent to, and probably
           larger than, the given one.

        -  ``algorithm`` - string (optional, default 'auto').  This
           only has an effect if working over the integers.  If 'dhsw',
           then preprocess each boundary matrix using the Dumas,
           Heckenbach, Saunders, and Welker elimination algorithm.  If
           'pari', then compute elementary divisors using Pari.  If
           'linbox', then use LinBox.  If 'auto', then use 'dhsw' for
           large matrices and Pari for small ones.

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the homology is computed.

        Algorithm: if ``subcomplex`` is None, replace it with a facet
        -- a contractible subcomplex of the original complex.  Then no
        matter what ``subcomplex`` is, replace it with a subcomplex
        `L` which is homotopy equivalent and as large as possible.
        Compute the homology of the original complex relative to `L`:
        if `L` is large, then the relative chain complex will be small
        enough to speed up computations considerably.

        EXAMPLES::

            sage: circle = SimplicialComplex(2, [[0,1], [1,2], [0, 2]])
            sage: circle.homology()
            {0: 0, 1: Z}
            sage: sphere = SimplicialComplex(3, [[0,1,2,3]])
            sage: sphere.remove_face([0,1,2,3])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: sphere.homology()
            {0: 0, 1: 0, 2: Z}

        Another way to get a two-sphere: take a two-point space and take its
        three-fold join with itself::

            sage: S = SimplicialComplex(1, [[0], [1]])
            sage: (S*S*S).homology(dim=2, cohomology=True)
            Z

        Relative homology::

            sage: T = SimplicialComplex(2, [[0,1,2]])
            sage: U = SimplicialComplex(2, [[0,1], [1,2], [0,2]])
            sage: T.homology(subcomplex=U)
            {0: 0, 1: 0, 2: Z}
        """
        if not (base_ring.is_field() or base_ring == ZZ):
            raise NotImplementedError, "Can't compute homology if the base ring is not the integers or a field."
        if dim is not None:
            if isinstance(dim, (list, tuple)):
                low = min(dim) - 1
                high = max(dim) + 2
            else:
                low = dim - 1
                high = dim + 2
            dims = range(low, high)
        else:
            dims = None
        if verbose:
            print "starting calculation of the homology of this",
            print "%s-dimensional simplicial complex" % self.dimension()
        if subcomplex is None:
            if enlarge:
                if verbose:
                    print "Constructing contractible subcomplex..."
                L = self._contractible_subcomplex(verbose=verbose)
                if verbose:
                    print "Done finding contractible subcomplex."
                    vec = [len(self.n_faces(n-1, subcomplex=L)) for n in range(self.dimension()+2)]
                    print "The difference between the f-vectors is:"
                    print "  %s" % vec
            else:
                L = SimplicialComplex(self.vertices(), [[self.vertices().tuple()[0]]])
        else:
            if enlarge:
                if verbose:
                    print "Enlarging subcomplex..."
                L = self._enlarge_subcomplex(subcomplex, verbose=verbose)
                if verbose:
                    print "Done enlarging subcomplex:"
            else:
                L = subcomplex
        if verbose:
            print "Computing the chain complex..."
        C = self.chain_complex(base_ring=base_ring, subcomplex=L,
                               dimensions=dims, augmented=True,
                               cochain=cohomology, verbose=verbose)
        if verbose:
            print " Done computing the chain complex. "
            print "Now computing homology..."
        answer = C.homology(algorithm=algorithm, verbose=verbose)
        if isinstance(answer, dict):
            if cohomology:
                too_big = self.dimension() + 1
                if (not ((isinstance(dim, (list, tuple)) and too_big in dim)
                        or too_big == dim)
                    and too_big in answer):
                    del answer[too_big]
            if -2 in answer:
                del answer[-2]
            if -1 in answer:
                del answer[-1]
            if dim is not None:
                if isinstance(dim, (list, tuple)):
                    temp = {}
                    for n in dim:
                        temp[n] = answer[n]
                    answer = temp
                else:  # just a single dimension
                    answer = answer[dim]
        return answer

    def cohomology(self, dim=None, base_ring=ZZ, subcomplex=None,
                   enlarge=True, algorithm='auto', verbose=False):
        """
        The reduced cohomology of this simplicial complex.

        INPUT:

        -  ``dim`` - integer or list of integers or None (optional,
           default None).  If None, then return the cohomology in
           every dimension.  If ``dim`` is an integer or list, return
           the cohomology in the given dimensions.  (Actually, if
           ``dim`` is a list, return the cohomology in the range from
           ``min(dim)`` to ``max(dim)``.)

        -  ``base_ring`` - commutative ring (optional, default ZZ).
           Must be ZZ or a field.

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default empty).  Compute cohomology relative to
           this subcomplex.

        -  ``enlarge`` - boolean (optional, default True).  If True,
           find a new subcomplex homotopy equivalent to, and probably
           larger than, the given one.

        -  ``algorithm`` - string (optional, default 'auto').  This
           only has an effect if working over the integers.  If
           'dhsw', then preprocess each boundary matrix using the
           Dumas, Heckenbach, Saunders, and Welker elimination
           algorithm.  If 'pari', then compute elementary divisors
           using Pari.  If 'linbox', then use LinBox.  If 'auto', then
           use 'dhsw' for large matrices and Pari for small ones.

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the cohomology is computed.

        EXAMPLES::

            sage: circle = SimplicialComplex(2, [[0,1], [1,2], [0, 2]])
            sage: circle.cohomology()
            {0: 0, 1: Z}
            sage: P2 = SimplicialComplex(5, [[0,1,2], [0,2,3], [0,1,5], [0,4,5], [0,3,4], [1,2,4], [1,3,4], [1,3,5], [2,3,5], [2,4,5]])   # projective plane
            sage: P2.cohomology()
            {0: 0, 1: 0, 2: C2}
            sage: P2.cohomology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
            1: Vector space of dimension 1 over Finite Field of size 2,
            2: Vector space of dimension 1 over Finite Field of size 2}
            sage: P2.cohomology(base_ring=GF(3))
            {0: Vector space of dimension 0 over Finite Field of size 3,
            1: Vector space of dimension 0 over Finite Field of size 3,
            2: Vector space of dimension 0 over Finite Field of size 3}

        Relative cohomology::

            sage: T = SimplicialComplex(1, [[0,1]])
            sage: U = SimplicialComplex(1, [[0], [1]])
            sage: T.cohomology(subcomplex=U)
            {0: 0, 1: Z}
        """
        return self.homology(dim=dim, base_ring=base_ring,
                             subcomplex=subcomplex,
                             cohomology=True, enlarge=enlarge,
                             algorithm=algorithm, verbose=verbose)

    def betti(self, dim=None, subcomplex=None):
        """
        The Betti numbers of this simplicial complex as a dictionary
        (or a single Betti number, if only one dimension is given).

        INPUT:

        -  ``dim`` - integer or list of integers or None (optional,
           default None).  If None, then return every Betti number, as
           a dictionary with keys the non-negative integers.  If
           ``dim`` is an integer or list, return the Betti number for
           each given dimension.  (Actually, if ``dim`` is a list,
           return the Betti numbers, as a dictionary, in the range
           from ``min(dim)`` to ``max(dim)``.  If ``dim`` is a number,
           return the Betti number in that dimension.)

        -  ``subcomplex`` - a subcomplex of this simplicial complex
           (optional, default None).  Compute the Betti numbers of the
           homology relative to this subcomplex.

        EXAMPLES: Build the two-sphere as a three-fold join of a
        two-point space with itself::

            sage: S = SimplicialComplex(1, [[0], [1]])
            sage: (S*S*S).betti()
            {0: 0, 1: 0, 2: 1}
            sage: (S*S*S).betti([1,2])
            {1: 0, 2: 1}
            sage: (S*S*S).betti(2)
            1
        """
        dict = {}
        H = self.homology(dim, base_ring=QQ, subcomplex=subcomplex)
        try:
            for n in H.keys():
                dict[n] = H[n].dimension()
            return dict
        except AttributeError:
            return H.dimension()

    def add_face(self, face):
        """
        Add a face to this simplicial complex

        INPUT:

        -  ``face`` - a subset of the vertex set

        This changes the simplicial complex, adding a new face and all
        of its subfaces.

        EXAMPLES::

            sage: X = SimplicialComplex(2, [[0,1], [0,2]])
            sage: X.add_face([0,1,2,]); X
            Simplicial complex with vertex set (0, 1, 2) and facets {(0, 1, 2)}
            sage: Y = SimplicialComplex(3); Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {()}
            sage: Y.add_face([0,1])
            sage: Y.add_face([1,2,3])
            sage: Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2, 3), (0, 1)}

        If you add a face which is already present, there is no effect::

            sage: Y.add_face([1,3]); Y
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2, 3), (0, 1)}
        """
        new_face = Simplex(face)
        if not new_face.is_face(self.vertices()):
            raise ValueError, "The face to be added is not a subset of the vertex set."
        else:
            face_is_maximal = True
            for other in self._facets:
                if face_is_maximal:
                    face_is_maximal = not new_face.is_face(other)
            if face_is_maximal:
                # remove any old facets which are no longer maximal
                Facets = list(self._facets)
                for old_face in self._facets:
                    if old_face.is_face(new_face):
                        Facets.remove(old_face)
                # add new_face to facet list
                Facets.append(new_face)
                self._facets = Facets
                # update self._faces if necessary
                if None in self._faces:
                    all_new_faces = SimplicialComplex(self.vertices(),
                                                      [new_face]).faces()
                    for dim in range(0, new_face.dimension()+1):
                        if dim in self._faces[None]:
                            self._faces[None][dim] = self._faces[None][dim].union(all_new_faces[dim])
                        else:
                            self._faces[None][dim] = all_new_faces[dim]
            return None

    def remove_face(self, face):
        """
        Remove a face from this simplicial complex

        INPUT:

        -  ``face`` - a face of the simplicial complex

        This changes the simplicial complex, removing the given face
        any face which contains it.

        Algorithm: take the Alexander dual, add the complement of
        ``face``, and then take the Alexander dual again.

        EXAMPLES::

            sage: S = range(1,5)
            sage: Z = SimplicialComplex(S, [S]); Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 2, 3, 4)}
            sage: Z.remove_face([1,2])
            sage: Z
            Simplicial complex with vertex set (1, 2, 3, 4) and facets {(1, 3, 4), (2, 3, 4)}
        """
        if not Simplex(face).is_face(self.vertices()):
            raise ValueError, "The face to be removed is not a subset of the vertex set."
        else:
            X = self.alexander_dual()
            X.add_face(self._complement(face))
            self._facets = X.alexander_dual()._facets
            if None in self._faces:
                s = Simplex(face)
                bad_faces = SimplicialComplex(self.vertices(), [s]).faces()
                for dim in range(0, s.dimension()+1):
                    self._faces[None][dim] = self._faces[None][dim].difference(bad_faces[dim])

    def link(self, simplex):
        """
        The link of a simplex in this simplicial complex.

        The link of a simplex `F` is the simplicial complex formed by
        all simplices `G` which are disjoint from `F` but for which `F
        \cup G` is a simplex.

        INPUT:

        -  ``simplex`` - a simplex in this simplicial complex.

        EXAMPLES::

            sage: X = SimplicialComplex(4, [[0,1,2], [1,2,3]])
            sage: from sage.homology.simplicial_complex import Simplex
            sage: X.link(Simplex([0]))
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(1, 2)}
            sage: X.link([1,2])
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {(3,), (0,)}
            sage: Y = SimplicialComplex(3, [[0,1,2,3]])
            sage: Y.link([1])
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3)}
        """
        faces = []
        s = Simplex(simplex)
        for f in self._facets:
            if s.is_face(f):
                faces.append(Simplex(list(f.set().difference(s.set()))))
        return SimplicialComplex(self.vertices(), faces)

    def _complement(self, simplex):
        """
        Return the complement of a simplex in the vertex set of this
        simplicial complex.

        INPUT:

        -  ``simplex`` - a simplex (need not be in the simplicial complex)

        OUTPUT: its complement: the simplex formed by the vertices not
        contained in ``simplex``.

        Note that this only depends on the vertex set of the
        simplicial complex, not on its simplices.

        EXAMPLES::

            sage: X = SimplicialComplex(5)
            sage: X._complement([1,2,3])
            (0, 4, 5)
            sage: X._complement([0,1,3,4])
            (2, 5)
            sage: X._complement([0,4,1,3])
            (2, 5)
        """
        return Simplex(set(self.vertices()).difference(simplex))

    def _transpose_simplices(self, *simplices):
        """
        Given tuple ``L`` of simplices, returns new list, where each
        simplex is formed by taking a vertex from each simplex from
        ``L``.

        INPUT:

        -  ``simplices`` - a bunch of simplices

        If ``simplices`` consists of `(f_0, f_1, f_2, ...)`, then the
        output consists of all possible simplices of the form `(v_0,
        v_1, v_2, ...)`, where `v_i` is a vertex of `f_i`.  If a
        vertex appears more than once in such a simplex, remove all
        but one of its appearances.  If such a simplex contains others
        already produced, then ignore that larger simplex -- the
        output should be a list of minimal simplices constructed in
        this way.

        This is used in computing the minimal nonfaces and hence the
        Stanley-Reisner ring.

        Note that this only depends on the vertex set of the
        simplicial complex, not on its simplices.

        I don't know if there is a standard name for this, but it
        looked sort of like the transpose of a matrix; hence the name
        for this method.

        EXAMPLES::

            sage: X = SimplicialComplex(5)
            sage: X._transpose_simplices([1,2])
            [(1,), (2,)]
            sage: X._transpose_simplices([1,2], [3,4])
            [(1, 3), (1, 4), (2, 3), (2, 4)]

        In the following example, one can construct the simplices
        (1,2) and (1,3), but you can also construct (1,1) = (1,),
        which is a face of both of the others.  So the answer omits
        (1,2) and (1,3)::

            sage: X._transpose_simplices([1,2], [1,3])
            [(1,), (2, 3)]
        """
        answer = []
        if len(simplices) == 1:
            answer = [Simplex((v,)) for v in simplices[0]]
        elif len(simplices) > 1:
            face = simplices[0]
            rest = simplices[1:]
            new_simplices = []
            for v in face:
                for partial in self._transpose_simplices(*rest):
                    if v not in partial:
                        L = [v] + list(partial)
                        L.sort()
                        simplex = Simplex(L)
                    else:
                        simplex = partial
                    add_simplex = True
                    simplices_to_delete = []
                    for already in answer:
                        if add_simplex:
                            if already.is_face(simplex):
                                add_simplex = False
                            if add_simplex and simplex.is_face(already):
                                simplices_to_delete.append(already)
                    if add_simplex:
                        answer.append(simplex)
                    for x in simplices_to_delete:
                        answer.remove(x)
        return answer

    def minimal_nonfaces(self):
        """
        Set consisting of the minimal subsets of the vertex set of
        this simplicial complex which do not form faces.

        Algorithm: first take the complement (within the vertex set)
        of each facet, obtaining a set `(f_1, f_2, ...)` of simplices.
        Now form the set of all simplices of the form `(v_1, v_2,
        ...)` where vertex `v_i` is in face `f_i`.  This set will
        contain the minimal nonfaces and may contain some non-minimal
        nonfaces also, so loop through the set to find the minimal
        ones.  (The last two steps are taken care of by the
        ``_transpose_simplices`` routine.)

        This is used in computing the Stanley-Reisner ring and the
        Alexander dual.

        EXAMPLES::

            sage: X = SimplicialComplex(4)
            sage: X.minimal_nonfaces()
            {(4,), (2,), (3,), (0,), (1,)}
            sage: X.add_face([1,2])
            sage: X.add_face([1,3])
            sage: X.minimal_nonfaces()
            {(4,), (2, 3), (0,)}
            sage: Y = SimplicialComplex(3, [[0,1], [1,2], [2,3], [3,0]])
            sage: Y.minimal_nonfaces()
            {(1, 3), (0, 2)}
        """
        complements = [self._complement(facet) for facet in self._facets]
        return Set(self._transpose_simplices(*complements))

    def _stanley_reisner_base_ring(self, base_ring=ZZ):
        """
        The polynomial algebra of which the Stanley-Reisner ring is a
        quotient.

        INPUT:

        -  ``base_ring`` - a commutative ring (optional, default ZZ)

        OUTPUT: a polynomial algebra with coefficients in base_ring,
        with one generator for each vertex in the simplicial complex.

        See the documentation for ``stanley_reisner_ring`` for a
        warning about the names of the vertices.

        EXAMPLES::

            sage: X = SimplicialComplex(3, [[1,2], [0]])
            sage: X._stanley_reisner_base_ring()
            Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring
            sage: Y = SimplicialComplex(['a', 'b', 'c'])
            sage: Y._stanley_reisner_base_ring(base_ring=QQ)
            Multivariate Polynomial Ring in a, c, b over Rational Field
        """
        return PolynomialRing(base_ring, self._gen_dict.values())

    def stanley_reisner_ring(self, base_ring=ZZ):
        """
        The Stanley-Reisner ring of this simplicial complex.

        INPUT:

        -  ``base_ring`` - a commutative ring (optional, default ZZ)

        OUTPUT:

        -  a quotient of a polynomial algebra with coefficients in
           ``base_ring``, with one generator for each vertex in the
           simplicial complex, by the ideal generated by the products
           of those vertices which do not form faces in it.

        Thus the ideal is generated by the products corresponding to
        the minimal nonfaces of the simplicial complex.

        .. warning::

           This may be quite slow!

           Also, this may behave badly if the vertices have the
           'wrong' names. To avoid this, define the simplicial complex
           at the start with the flag ``name_check`` set to True.

           More precisely, this is a quotient of a polynomial ring
           with one generator for each vertex.  If the name of a
           vertex is a non-negative integer, then the corresponding
           polynomial generator is named 'x' followed by that integer
           (e.g., 'x2', 'x3', 'x5', ...).  Otherwise, the polynomial
           generators are given the same names as the vertices.  Thus
           if the vertex set is (2, 'x2'), there will be problems.

        EXAMPLES::

            sage: X = SimplicialComplex(3, [[0,1], [1,2], [2,3], [0,3]])
            sage: X.stanley_reisner_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3 over Integer Ring by the ideal (x1*x3, x0*x2)
            sage: Y = SimplicialComplex(4); Y
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {()}
            sage: Y.stanley_reisner_ring()
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Integer Ring by the ideal (x4, x2, x3, x0, x1)
            sage: Y.add_face([0,1,2,3,4])
            sage: Y.stanley_reisner_ring(base_ring=QQ)
            Quotient of Multivariate Polynomial Ring in x0, x1, x2, x3, x4 over Rational Field by the ideal (0)
        """
        R = self._stanley_reisner_base_ring(base_ring)
        products = []
        for f in self.minimal_nonfaces():
            prod = 1
            for v in f:
                prod *= R(self._gen_dict[v])
            products.append(prod)
        return R.quotient(products)

    def alexander_dual(self):
        """
        The Alexander dual of this simplicial complex: according to
        the Macaulay2 documentation, this is the simplicial complex
        whose faces are the complements of its nonfaces.

        Thus find the minimal nonfaces and take their complements to
        find the facets in the Alexander dual.

        EXAMPLES::

            sage: Y = SimplicialComplex(4); Y
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and facets {()}
            sage: Y.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3, 4) and 5 facets
            sage: X = SimplicialComplex(3, [[0,1], [1,2], [2,3], [3,0]])
            sage: X.alexander_dual()
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 3), (0, 2)}
        """
        nonfaces = self.minimal_nonfaces()
        return SimplicialComplex(self.vertices(), [self._complement(f) for f in nonfaces])

    def face_poset(self):
        """
        The face poset of this simplicial complex, the poset of
        nonempty faces, ordered by inclusion.

        EXAMPLES::

            sage: P = SimplicialComplex(3, [[0, 1], [1,2], [2,3]]).face_poset(); P
            Finite poset containing 7 elements
            sage: P.list()
            [(3,), (2,), (2, 3), (1,), (0,), (0, 1), (1, 2)]
        """
        from sage.combinat.posets.posets import Poset
        covers = {}
        # The code for posets seems to work better if the elements are
        # from the class tuple instead of Simplex; thus each face is
        # converted to a tuple here.
        for simplex in self._flattened_faces():
            if not simplex.is_empty():
                covers[tuple(simplex)] = []
        for simplex in self._flattened_faces():
            for face in simplex.faces():
                if not face.is_empty():
                    covers[tuple(face)].append(tuple(simplex))
        return Poset(covers)

    def barycentric_subdivision(self):
        """
        The barycentric subdivision of this simplicial complex.

        See http://en.wikipedia.org/wiki/Barycentric_subdivision for a
        definition.

        EXAMPLES::

            sage: triangle = SimplicialComplex(2, [[0,1], [1,2], [0, 2]])
            sage: hexagon = triangle.barycentric_subdivision()
            sage: hexagon
            Simplicial complex with 6 vertices and 6 facets
            sage: hexagon.homology(1) == triangle.homology(1)
            True

        Barycentric subdivisions can get quite large, since each
        `n`-dimensional facet in the original complex produces
        `(n+1)!` facets in the subdivision::

            sage: S4 = simplicial_complexes.Sphere(4)
            sage: S4
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5) and 6 facets
            sage: S4.barycentric_subdivision()
            Simplicial complex with 62 vertices and 720 facets
        """
        return self.face_poset().order_complex()

    def graph(self):
        """
        The 1-skeleton of this simplicial complex, as a graph.

        EXAMPLES::

            sage: S = SimplicialComplex(3, [[0,1,2,3]])
            sage: G = S.graph(); G
            Graph on 4 vertices
            sage: G.edges()
            [(0, 1, None), (0, 2, None), (0, 3, None), (1, 2, None), (1, 3, None), (2, 3, None)]
        """
        edges = self.n_faces(1)
        d = {}
        for e in edges:
            v = min(e)
            if v in d:
                d[v].append(max(e))
            else:
                d[v] = [max(e)]
        return Graph(d)

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this simplicial complex: the simplicial
        complex obtained by discarding all of the simplices in
        dimensions larger than `n`.

        INPUT:

        -  ``n`` - non-negative integer

        EXAMPLES::

            sage: X = SimplicialComplex(3, [[0,1], [1,2,3], [0,2,3]])
            sage: X.n_skeleton(1)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(2, 3), (0, 2), (1, 3), (1, 2), (0, 3), (0, 1)}
            sage: X.n_skeleton(2)
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (1, 2, 3), (0, 1)}
        """
        facets = filter(lambda f: f.dimension()<n, self._facets)
        facets.extend(self.n_faces(n))
        return SimplicialComplex(self.vertices(), facets)

    def _contractible_subcomplex(self, verbose=False):
        """
        Find a contractible subcomplex `L` of this simplicial complex,
        preferably one which is as large as possible.

        INPUT:

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the simplicial complex is computed.

        Motivation: if `K` is the original complex and if `L` is
        contractible, then the relative homology `H_*(K,L)` is
        isomorphic to the reduced homology of `K`.  If `L` is large,
        then the relative chain complex will be a good deal smaller
        than the augmented chain complex for `K`, and this leads to a
        speed improvement for computing the homology of `K`.

        This just passes a subcomplex consisting of a facet to the
        method ``_enlarge_subcomplex``.

        .. note::

           Thus when the simplicial complex is empty, so is the
           resulting 'contractible subcomplex', which is therefore not
           technically contractible.  In this case, that doesn't
           matter because the homology is computed correctly anyway.

        EXAMPLES::

            sage: sphere = SimplicialComplex(3, [[0,1,2,3]])
            sage: sphere.remove_face([0,1,2,3])
            sage: sphere
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 2, 3), (0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: L = sphere._contractible_subcomplex(); L
            Simplicial complex with vertex set (0, 1, 2, 3) and facets {(0, 1, 2), (1, 2, 3), (0, 1, 3)}
            sage: L.homology()
            {0: 0, 1: 0, 2: 0}
        """
        vertices = self.vertices()
        facets = [self._facets[0]]
        return self._enlarge_subcomplex(SimplicialComplex(vertices, facets), verbose=verbose)

    def _enlarge_subcomplex(self, subcomplex, verbose=False):
        """
        Given a subcomplex `S` of this simplicial complex `K`, find a
        subcomplex `L`, as large as possible, containing `S` which is
        homotopy equivalent to `S` (so that `H_{*}(K,S)` is isomorphic
        to `H_{*}(K,L)`).  This way, the chain complex for computing
        `H_{*}(K,L)` will be smaller than that for computing
        `H_{*}(K,S)`, so the computations should be faster.

        INPUT:

        -  ``subcomplex`` - a subcomplex of this simplicial complex

        -  ``verbose`` - boolean (optional, default False).  If True,
           print some messages as the simplicial complex is computed.

        OUTPUT: a complex `L` containing ``subcomplex`` and contained
        in ``self``, homotopy equivalent to ``subcomplex``.

        Algorithm: start with the subcomplex `S` and loop through the
        facets of `K` which are not in `S`.  For each one, see whether
        its intersection with `S` is contractible, and if so, add it.
        This is recursive: testing for contractibility calls this
        routine again, via ``_contractible_subcomplex``.

        EXAMPLES::

            sage: T = simplicial_complexes.Torus(); T
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 14 facets

        Inside the torus, define a subcomplex consisting of a loop:

            sage: S = SimplicialComplex(T.vertices(), [[0,1], [1,2], [0,2]])
            sage: S.homology()
            {0: 0, 1: Z}
            sage: L = T._enlarge_subcomplex(S)
            sage: L
            Simplicial complex with vertex set (0, 1, 2, 3, 4, 5, 6) and 8 facets
            sage: L.facets()
            {(0, 1, 5), (1, 3, 6), (1, 2), (1, 2, 4), (1, 3, 4), (0, 2), (1, 5, 6), (0, 1)}
            sage: L.homology()
            {0: 0, 1: Z, 2: 0}
        """
        if subcomplex in self.__enlarged:
            return self.__enlarged[subcomplex]
        faces = filter(lambda x: x not in subcomplex._facets, self._facets)
        done = False
        new_facets = list(subcomplex._facets)
        while not done:
            done = True
            remove_these = []
            if verbose:
                print "  looping through %s facets" % len(faces)
            for f in faces:
                f_set = f.set()
                int_facets = []
                for a in new_facets:
                    int_facets.append(a.set().intersection(f_set))
                intersection = SimplicialComplex(self.vertices(), int_facets)
                if not intersection._facets[0].is_empty():
                    if (len(intersection._facets) == 1 or
                        intersection == intersection._contractible_subcomplex()):
                        new_facets.append(f)
                        remove_these.append(f)
                        done = False
            if verbose and not done:
                print "    added %s facets" % len(remove_these)
            for f in remove_these:
                faces.remove(f)
        if verbose:
            print "  now constructing a simplicial complex with %s vertices and %s facets" % (self.vertices().dimension()+1, len(new_facets))
        L = SimplicialComplex(self.vertices(), new_facets, maximality_check=False,
                              vertex_check=False, sort_facets=False)
        self.__enlarged[subcomplex] = L
        return L

    def category(self):
        """
        Return the category to which this chain complex belongs: the
        category of all simplicial complexes.

        EXAMPLES::

            sage: SimplicialComplex(5, [[0,1], [1,2,3,4,5]]).category()
            Category of simplicial complexes
        """
        import sage.categories.all
        return sage.categories.all.SimplicialComplexes()

    def _repr_(self):
        """
        Print representation. If there are only a few vertices or
        faces, they are listed.  If there are lots, the number is
        given.

        EXAMPLES::

            sage: X = SimplicialComplex(3, [[0,1], [1,2]])
            sage: X._repr_()
            'Simplicial complex with vertex set (0, 1, 2, 3) and facets {(1, 2), (0, 1)}'
            sage: SimplicialComplex(15)
            Simplicial complex with 16 vertices and facets {()}
        """
        vertex_limit = 45
        facet_limit = 55
        vertices = self.vertices()
        facets = Set(self._facets)
        vertex_string = "with vertex set %s" % vertices
        if len(vertex_string) > vertex_limit:
            vertex_string = "with %s vertices" % str(1+vertices.dimension())
        facet_string = "facets %s" % facets
        if len(facet_string) > facet_limit:
            facet_string = "%s facets" % len(facets)
        return "Simplicial complex " + vertex_string + " and " + facet_string
