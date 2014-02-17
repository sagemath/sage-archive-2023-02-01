#*****************************************************************************
#  Copyright (C) 2012 Jim Stark <jstarx@gmail.com>
#                2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.factory import UniqueFactory
from sage.modules.module import Module
from sage.modules.module_element import ModuleElement

class QuiverRepFactory(UniqueFactory):
    """
    A quiver representation is a diagram in the category of vector spaces whose
    underlying graph is the quiver.  Giving a finite dimensional representation
    is equivalent to giving a finite dimensional right module for the path
    algebra of the quiver.

    INPUT:

    The first two arguments specify the base ring and the Quiver, and they are
    always required:

    - ``k`` - ring, the base ring of the representation

    - ``P`` - the partial semigroup formed by the paths of the quiver of the
      representation

    Then to specify the spaces and maps associated to the Quiver there are three
    possible options.  The first is the ``values`` option, where the next two
    arguments give the data to be assigned.  The following can either be the next
    two entries in the argument list or they can be passed by keyword.  If the
    argument list is long enough the keywords are ignored; the keywords are only
    checked in the event that the argument list does not have enough entries after
    ``P``.

    - ``spaces`` - dict (default: empty), a dictionary associating to each vertex a
      free module over the base ring `k`.  Not all vertices must be specified;
      unspecified vertices are automatically set to `k^0`.  Keys of the dictionary
      that don't correspond to vertices are ignored.

    - ``maps`` - dict (default: empty), a dictionary associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges must be specified;
      unspecified edges are automatically set to the zero map.  Keys of the
      dictionary that don't correspond to edges are ignored.

    The second option is the ``paths`` option which creates a module by generating a
    right ideal from a list of paths.  Thus the basis elements of this module
    correspond to paths of the quiver and the maps are given by right
    multiplication by the corresponding edge.  As above this can be passed either
    as the next entry in the argument list or as a keyword.  The keyword is only
    checked if there is no entry in the argument list after ``Q``.

    - ``basis`` - list, a nonempty list of paths in the quiver Q.  Entries that do not
      represent valid paths are ignored and duplicate paths are deleted.  There
      must be at least one valid path in the list or a ValueError is raised.  The
      closure of this list under right multiplication forms the basis of the
      resulting representation.

    The third option is the ``dual paths`` option which creates the dual of a left
    ideal in the quiver algebra.  Thus the basis elements of this module correspond
    to paths of the quiver and the maps are given by deleting the corresponding
    edge from the start of the path (the edge map is zero on a path if that edge is
    not the initial edge of the path).  As above this can be passed either as the
    next entry in the argument list or as a keyword.

    - ``basis`` - list, a nonempty list of paths in the quiver Q.  Entries that do not
      represent valid paths are ignored and duplicate paths are deleted.  There
      must be at least one valid path in the list or a ValueError is raised.  The
      closure of this list under left multiplication of edges forms the basis of
      the resulting representation.

    Using the second and third options requires that the following keyword be
    passed to the constructor.  This must be passed as a keyword.

    - ``option`` - string (default: ``None``), either 'values' or 'paths' or
      'dual paths'.  ``None`` is equivalent to 'values'.

    OUTPUT:

    - QuiverRep

    EXAMPLES::

        sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()

    When the ``option`` keyword is not supplied the constructor uses the 'values'
    option and expects the spaces and maps to be specified.  If no maps or spaces
    are given the zero module is created::

        sage: M = Q1.representation(GF(5))
        sage: M.is_zero()
        True

    The simple modules, indecomposable projectives, and indecomposable injectives are
    examples of quiver representations::

        sage: S = Q1.S(GF(3), 1)
        sage: I = Q1.I(QQ, 2)
        sage: P = Q1.P(GF(3), 1)

    Various standard submodules can be computed, such as radicals and socles.  We can
    also form quotients and test for certain attributes such as semisimplicity::

        sage: R = P.radical()
        sage: R.is_zero()
        False
        sage: (P/R).is_simple()
        True
        sage: P == R
        False

    With the option 'paths' the input data should be a list of QuiverPaths or
    things that QuiverPaths can be constructed from.  The resulting module is the
    submodule generated by these paths in the quiver algebra, when considered as a
    right module over itself::

        sage: P1 = Q1.representation(QQ, [[(1, 1)]], option='paths')
        sage: P1.dimension()
        2

    Notice that the 3rd and 4th paths are actually the same, the duplicate is
    removed::

        sage: N = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='paths')
        sage: N.dimension()
        3

    The dimension at each vertex equals the number of paths in the closed basis whose
    terminal point is that vertex::

        sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
        sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
        sage: M.dimension_vector()
        (0, 2, 2)
        sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
        sage: N.dimension_vector()
        (0, 1, 2)
    """

    def create_key(self, k, P, *args, **kwds):
        """
        Returns a key for the specified module.

        The key is a tuple.  The first and second entries are the base ring
        ``k`` and the partial semigroup ``P`` formed by the paths of a quiver.
        The third entry is the 'option' and the remaining entries depend on
        that option.  If the option is 'values' and the Quiver has n vertices
        then the next n entries are the vector spaces to be assigned to those
        vertices.  After that are the matrices of the maps assigned to edges,
        listed in the same order that ``Q.edges()`` uses.  If the option is
        'paths' or 'dual paths' then the next entry is a tuple containing a
        sorted list of the paths that form a basis of the Quiver.

        INPUT:

        - See the class documentation

        OUTPUT:

        - tuple

        EXAMPLES::

            sage: P = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: from sage.quivers.representation import QuiverRep
            sage: QuiverRep.create_key(GF(5), P)
            (Finite Field of size 5, Partial semigroup formed by the directed paths of Multi-digraph on 2 vertices, 'values', Vector space of dimension 0 over Finite Field of size 5, Vector space of dimension 0 over Finite Field of size 5, [])
        """
        key = [k, P]
        Q = P.quiver()
        if 'option' in kwds and (kwds['option'] == 'paths' or kwds['option'] == 'dual paths'):
            # Follow the 'paths' specification for the input
            key.append(kwds['option'])
            if args:
                basis = args[0]
            else:
                basis = kwds['basis']

            # Add as QuiverPaths to a set
            paths = set()
            for p in basis:
                paths.add(P(p))

            if kwds['option'] == 'paths':
                # Close the set under right mult by edges
                edges = [P(e) for e in Q.edges()]
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            pe = p*e
                            if pe is not None and pe not in paths:
                                to_be_added.append(pe)

                    paths.update(to_be_added)
                    just_added = to_be_added

            if kwds['option'] == 'dual paths':
                # Close the set under left mult by edges
                edges = [P(e) for e in Q.edges()]
                just_added = paths
                while just_added:
                    to_be_added = []
                    for e in edges:
                        for p in just_added:
                            ep = e*p
                            if ep is not None and ep not in paths:
                                to_be_added.append(ep)

                    paths.update(to_be_added)
                    just_added = to_be_added

            # Add to the key
            key.append(tuple(sorted(paths)))

        else:
            # Assume the input type the 'values' option
            key.append('values')
            if args:
                spaces = args[0]
            elif 'spaces' in kwds:
                spaces = kwds['spaces']
            else:
                spaces = {}
            if len(args) > 1:
                maps = args[1]
            elif 'maps' in kwds:
                maps = kwds['maps']
            else:
                maps = {}

            # If the vertex is not specified set it as a free module of rank 0, if
            # an integer is given set it as a free module of that rank, otherwise
            # assume the object is a module and assign it to the vertex.
            from sage.rings.finite_rings.integer_mod_ring import Integers
            verts = Q.vertices()
            for x in verts:
                if x not in spaces:
                    key.append(k**0)
                elif spaces[x] in Integers():
                    key.append(k**spaces[x])
                else:
                    key.append(spaces[x])

            # The preferred method of specifying an edge is as a tuple
            # (i, t, l) where i is the initial vertex, t is the terminal
            # vertex, and l is the label.  This is the form in which
            # quiver.edges() and other such functions give the edge.  But here
            # edges can be specified by giving only the two vertices or giving
            # only the edge label.
            #
            # Note that the first space is assigned to key[3] and the first
            # vertex is 1 so the space assigned to vertex v is key[2 + v]
            from sage.matrix.constructor import Matrix
            from sage.categories.morphism import is_Morphism
            for x in Q.edges():
                if x in maps:
                    e = maps[x]
                elif (x[0], x[1]) in maps:
                    e = maps[(x[0], x[1])]
                elif x[2] in maps:
                    e = maps[x[2]]
                else:
                    e = Matrix(k, key[3 + verts.index(x[0])].dimension(), key[3 + verts.index(x[1])].dimension())

                # If a morphism is specified take it's matrix.  Create one if
                # needed.  Otherwise assume the Matrix function can convert the
                # object to a Matrix.
                if is_Morphism(e):
                    if hasattr(e, 'matrix'):
                        key.append(e.matrix())
                    else:
                        gens_images = [key[3 + verts.index(x[1])].coordinate_vector(e(x))
                                        for x in key[3 + verts.index(x[0])].gens()]
                        key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                          key[3 + verts.index(x[1])].dimension(), gens_images))
                else:
                    key.append(Matrix(k, key[3 + verts.index(x[0])].dimension(),
                                      key[3 + verts.index(x[1])].dimension(), e))

                # Make sure the matrix is immutable so it hashes
                key[-1].set_immutable()

        # Wrap as a tuple and return
        return tuple(key)

    def create_object(self, version, key, **extra_args):
        """
        Creates a QuiverRep_generic or QuiverRep_with_path_basis object from the key.

        The key is a tuple.  The first and second entries are the base ring
        ``k`` and the Quiver ``Q``.  The third entry is the 'option' and the
        remaining entries depend on that option.  If the option is 'values'
        and the Quiver has n vertices then the next n entries are the vector
        spaces to be assigned to those vertices.  After that are the matrices
        of the maps assigned to edges, listed in the same order that
        ``Q.edges()`` uses.  If the option is 'paths' or 'dual paths' then the
        next entry is a tuple containing a sorted list of the paths that form
        a basis of the Quiver.

        INPUT:

        - ``version`` - The version of sage, this is currently ignored.

        - ``key`` - tuple

        OUTPUT:

        - QuiverRep_generic or QuiverRep_with_path_basis object

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: from sage.quivers.representation import QuiverRep
            sage: key = QuiverRep.create_key(GF(5), Q)
            sage: QuiverRep.create_object(Q.version(), key)
            Representation with dimension vector (0, 0)
        """

        if len(key) < 4:
            raise ValueError("Invalid key used in QuiverRepFactory!")

        # Get the quiver
        P = key[1]
        Q = P.quiver()

        if key[2] == 'values':
            # Get the spaces
            spaces = {}
            i = 3
            for v in Q:
                spaces[v] = key[i]
                i += 1

            # Get the maps
            maps = {}
            for e in Q.edges():
                maps[e] = key[i]
                i += 1

            # Create and return the module
            return QuiverRep_generic(key[0], P, spaces, maps)

        elif key[2] == 'paths':
            # Create and return the module
            return QuiverRep_with_path_basis(key[0], P, key[3])

        elif key[2] == 'dual paths':
            # Create and return the module
            return QuiverRep_with_dual_path_basis(key[0], P, key[3])

        else:
            raise ValueError("Invalid key used in QuiverRepFactory!")

QuiverRep = QuiverRepFactory("sage.quivers.representation.QuiverRep")

#########################################################################
##  Elements

class QuiverRepElement(ModuleElement):
    """
    An element of a quiver representation is a choice of element from each of the
    spaces assigned to the vertices of the quiver.  Addition, subtraction, and
    scalar multiplication of these elements is done pointwise within these spaces.

    INPUT:

    - ``module`` - QuiverRep (default: None), the module to which the element belongs

    - ``elements`` - dict (default: empty), a dictionary associating to each vertex a
      vector or an object from which sage can create a vector.  Not all vertices must
      be specified, unspecified vertices will be assigned the zero vector of the space
      associated to that vertex in the given module.  Keys that do not correspond to a
      vertex are ignored.

    - ``name`` - string (default: None), the name of the element

    OUTPUT:

    - QuiverRepElement

    .. NOTES::

        The constructor needs to know the quiver in order to create an element of a
        representation over that quiver.  The default is to read this information
        from ``module`` as well as to fill in unspecified vectors with the zeros of the
        spaces in ``module``.  If ``module``=None then ``quiver`` MUST be a quiver and each
        vertex MUST be specified or an error will result.  If both ``module`` and
        ``quiver`` are given then ``quiver`` is ignored.

    EXAMPLES::

        sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
        sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
        sage: M = Q.representation(GF(3), spaces)
        sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
        sage: M(elems)
        Element of quiver representation
        sage: v = M(elems, 'v')
        sage: v
        v
        sage: (v + v + v).is_zero()
        True
    """

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, parent, elements=None, name=None):
        """
        Type QuiverRepElement? for more information.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: M()
            Element of quiver representation
            sage: v = M(elems, 'v')
            sage: v
            v
            sage: (v + v + v).is_zero()
            True
        """
        # In the default call method, the default value of the first argument is zero
        if not elements:
            elements = {}
        # The data describing an element is held in the following private
        # variables:
        #
        # * _elems
        #      A dictionary that assigns to each vertex of the quiver a choice
        #      of element from the space assigned to that vertex in the parent
        #      representation.
        # * _quiver
        #      The quiver of the representation.

        super(QuiverRepElement, self).__init__(parent)

        self._elems = {}
        self._quiver = parent._quiver
        for v in self._quiver:
            if v in elements:
                self._elems[v] = parent._spaces[v](elements[v])
            else:
                self._elems[v] = parent._spaces[v].zero()

        # Assign a name if supplied
        if name is not None:
            self.rename(name)

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: Q.P(QQ, 3).an_element() # indirect doctest
            Element of quiver representation
        """
        return "Element of quiver representation"

    def _add_(left, right):
        """
        This overrides the + operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: (v + v + v).is_zero() # indirect doctest
            True
        """

        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] + right._elems[v]

        return left.parent()(elements)

    def _iadd_(self, right):
        """
        This overrides the += operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = M(elems)
            sage: v += w # indirect doctest
            sage: v == w
            False
            sage: v.is_zero()
            False
        """

        for v in self._quiver:
            self._elems[v] += right._elems[v]

        return self

    def _sub_(left, right):
        """
        This overrides the - operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: (v - v).is_zero() # indirect doctest
            True
        """

        elements = {}
        for v in left._quiver:
            elements[v] = left._elems[v] - right._elems[v]

        return left.parent()(elements)

    def _isub_(self, right):
        """
        This overrides the -= operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v -= v # indirect doctest
            sage: v.is_zero()
            True
        """

        for v in self._quiver:
            self._elems[v] -= right._elems[v]

        return self

    def _neg_(self):
        """
        This overrides the unary - operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v == -v # indirect doctest
            False
        """

        elements = {}
        for v in self._quiver:
            elements[v] = -self._elems[v]

        return self.parent()(elements)

    def __mul__(self, other):
        """
        Implements * for right multiplication by Quiver Algebra elements.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: A = Q.algebra(QQ)
            sage: m = P.an_element()
            sage: a = A((1, 2, 'a'))
            sage: e1 = A((1, 1))
            sage: m.support()
            [1, 2]
            sage: (m*a).support()
            [2]
            sage: (m*e1).support()
            [1]
            sage: (m*(a + e1)).support()
            [1, 2]
        """

        # Make sure the input is an element of the quiver algebra and get the
        # coefficients of the monomials in it
        parent = self.parent()
        mons = parent._base(other).monomial_coefficients()
        result = parent()   # this must not be the cached parent.zero(),
                             # since otherwise it gets changed in place!!

        for path in mons:
            # Multiply by the scalar
            x = mons[path]*self._elems[path.initial_vertex()]

            # If the edge isn't trivial apply the corresponding maps
            if len(path) > 0:
                for edge in path:
                    x = parent.get_map(edge)(x)

            # Sum the results of acting by each monomial
            result._elems[path.terminal_vertex()] += x

        return result

    def __eq__(self, other):
        """
        This overrides the == operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = M(elems)
            sage: v == w
            True
            sage: v += w
            sage: v == w
            False
        """

        # Return False if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return False

        # Return False if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return False

        return True

    def __ne__(self, other):
        """
        This overrides the != operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = M(elems)
            sage: v != w
            False
            sage: v += w
            sage: v != w
            True
        """

        # Return True if being compared to something other than a
        # QuiverRepElement or if comparing two elements from representations
        # with different quivers
        if not isinstance(other, QuiverRepElement) or self._quiver != other._quiver:
            return True

        # Return True if the elements differ at any vertex
        for v in self._quiver:
            if self._elems[v] != other._elems[v]:
                return True

        return False

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view and modify the representation data. #
    #                                                                         #
    ###########################################################################

    def quiver(self):
        """
        Returns the quiver of the representation.

        OUTPUT:

        - Quiver, the quiver of the representation

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: P = Q.P(QQ, 1)
            sage: v = P.an_element()
            sage: v.quiver() is Q.quiver()
            True
        """

        return self._quiver

    def get_element(self, vertex):
        """
        Returns the element at the given vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver

        OUTPUT:

        - vector, the vaector assigned to the given vertex

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v.get_element(3)
            (2, 1)
        """

        return self._elems[vertex]

    def set_element(self, vector, vertex):
        """
        Sets the element at the given vertex.

        INPUT:

        - ``vector`` - a vector or an object from which the space associated to the given
          vertex in the parent can create a vector

        - ``vertex`` - integer, a vertex of the quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.get_element(1)
            (1, 0)
            sage: v.set_element((1, 1), 1)
            sage: v.get_element(1)
            (1, 1)
          """

        self._elems[vertex] = self.parent()._spaces[vertex](vector)

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the homomorphism.         #
    #                                                                         #
    ###########################################################################

    def is_zero(self):
        """
        Tests whether self is zero.

        OUTPUT:

        - bool, True is the element is the zero element, False otherwise

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.is_zero()
            False
            sage: w = M()
            sage: w.is_zero()
            True

        TESTS::

            sage: M = Q.P(QQ, 1)
            sage: M.zero().is_zero()
            True

        """

        for v in self._quiver:
            if not self._elems[v].is_zero():
                return False

        return True

    def support(self):
        """
        Returns the support of self as a list.

        The support is the set of vertices to which a nonzero vector is associated.

        OUTPUT

        - list, the support

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 0), 3: (2, 1)}
            sage: v = M(elems)
            sage: v.support()
            [1, 3]
        """

        sup = []
        for v in self._quiver:
            if not self._elems[v].is_zero():
                sup.append(v)

        return sup

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def copy(self):
        """
        Returns a copy of self.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{3:['c']}}).path_semigroup()
            sage: spaces = dict((v, GF(3)^2) for v in Q.quiver())
            sage: M = Q.representation(GF(3), spaces)
            sage: elems = {1: (1, 0), 2: (0, 1), 3: (2, 1)}
            sage: v = M(elems)
            sage: w = v.copy()
            sage: w.set_element((0, 0), 1)
            sage: w.get_element(1)
            (0, 0)
            sage: v.get_element(1)
            (1, 0)
        """

        if hasattr(self, '__custom_name'):
            name = self.__custom_name
        else:
            name = None
        return self.parent()(self._elems.copy(), name)

####################################################################
# The representations

class QuiverRep_generic(Module):
    """
    This function should not be called by the user.

    Call QuiverRep with option='values' (which is the default) instead.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``spaces`` - dict (default: empty), a dictionary associating to each vertex a
      free module over the base ring k.  Not all vertices need to be specified,
      unspecified vertices are automatically set to k^0.  Keys of the dictionary
      that don't correspond to vertices are ignored.

    - ``maps`` - dict (default: empty), a dictionary associating to each edge a map
      whose domain and codomain are the spaces associated to the initial and
      terminal vertex of the edge respectively.  Not all edges need to be specified,
      unspecified edges are automatically set to the zero map.  Keys of the
      dictionary that don't correspond to edges are ignored.

    OUTPUT:

    - QuiverRep

    EXAMPLES::

        sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
        sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
        sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
        sage: M = Q.representation(QQ, spaces, maps)

    ::

        sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
        sage: P = Q.P(GF(3), 1)
        sage: I = Q.I(QQ, 1)
        sage: P.an_element() in P
        True
        sage: I.an_element() in P
        False

    TESTS::

        sage: TestSuite(M).run()
        sage: TestSuite(P).run()
        sage: TestSuite(I).run()
        sage: TestSuite(Q.S(ZZ,2)).run()

    """
    Element = QuiverRepElement

    ###########################################################################
    #                                                                         #
    # PRIVATE FUNCTIONS                                                       #
    #    These functions are not meant to be seen by the end user.            #
    #                                                                         #
    ###########################################################################

    def __init__(self, k, P, spaces, maps):
        """
        Type QuiverRep? for more information.

        NOTE:

        Use :meth:`~sage.quivers.quiver.Quiver.representation` to create a
        representation. The following is just a test, but it is not the
        recommended way of creating a representation.

        TESTS::

            sage: from sage.quivers.representation import QuiverRep_generic
            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: QuiverRep_generic(GF(5), Q, {1: GF(5)^2, 2: GF(5)^3}, {})
            Representation with dimension vector (2, 3)

        """
        # This class can handle representations over
        # an arbitrary base ring, not necessarily a field, so long as sage can
        # construct free modules over that ring.  The data of a representation is held
        # in the following private variables:
        #
        # * _quiver
        #      The quiver of the representation.
        # * _base
        #       The quiver algebra of _quiver over _base_ring
        # * _base_ring
        #      The base ring of the representation.
        # * _spaces
        #      A dictionary which associates to each vertex of the quiver a free
        #      module over the base ring.
        # * _maps
        #      A dictionary which associates to each edge of the quiver a homomorphism
        #      whose domain and codomain equal the initial and terminal vertex of the
        #      edge.
        Q = P.quiver()
        self._semigroup = P
        self._quiver = Q
        self._base_ring = k
        self._spaces = {}
        self._maps = {}

        # If the vertex is not specified set it as a free module of rank 0, if
        # an integer is given set it as a free module of that rank, otherwise
        # assume the object is a module and assign it to the vertex.
        from sage.rings.finite_rings.integer_mod_ring import Integers
        for x in Q:
            if x not in spaces:
                self._spaces[x] = k**0
            elif spaces[x] in Integers():
                self._spaces[x] = k**spaces[x]
            else:
                self._spaces[x] = spaces[x]

        # The preferred method of specifying an edge is as a tuple (i, t, l)
        # where i is the initial vertex, t is the terminal vertex, and l is the
        # label.  This is the form in which quiver.edges() and other such
        # functions give the edge.  But here edges can be specified by giving
        # only the two vertices or giving only the edge label.
        for x in Q.edges():
            if x in maps:
                e = maps[x]
            elif (x[0], x[1]) in maps:
                e = maps[(x[0], x[1])]
            elif x[2] in maps:
                e = maps[x[2]]
            else:
                e = self._spaces[x[0]].Hom(self._spaces[x[1]]).zero_element()

            #If a morphism is specified use it, otherwise assume the hom
            # function can convert the object to a morphism.  Matrices and the
            # zero and one of the base ring are valid inputs (one is valid only
            # when the domain and codomain are equal.
            from sage.categories.morphism import Morphism
            if isinstance(e, Morphism):
                self._maps[x] = e
            else:
                self._maps[x] = self._spaces[x[0]].hom(e, self._spaces[x[1]])

        self._assert_valid_quiverrep()

        super(QuiverRep_generic, self).__init__(P.algebra(k))

    def _assert_valid_quiverrep(self):
        """
        Raises an error if the representation is not well defined.

        Specifically it checks the map assigned to each edge.  The domain and codomain
        must equal the modules assigned to the initial and terminal vertices of the
        edge.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.P(GF(3), 2) # indirect doctest

        Due to unique representation, we will cause bugs in later code if we modify M
        to be an invalid representation.  So we make sure to store the original values
        and replace them after testing::

            sage: sv = M._spaces[1]
            sage: M._spaces[1] = 0
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: Domain of map at edge 'a' does not match.
            sage: M._spaces[1] = sv
            sage: M = Q.P(GF(3), 2)
            sage: sv = M._maps[(1, 2, 'a')]
            sage: M._maps[(1, 2, 'a')] = (QQ^2).Hom(QQ^1).zero_element()
            sage: M._assert_valid_quiverrep()
            Traceback (most recent call last):
            ...
            ValueError: Domain of map at edge 'a' does not match.
            sage: M._maps[(1, 2, 'a')] = sv
        """

        for x in self._quiver.edges():
            if self._maps[x].domain() != self._spaces[x[0]]:
                raise ValueError("Domain of map at edge '" + str(x[2]) + "' does not match.")
            if self._maps[x].codomain() != self._spaces[x[1]]:
                raise ValueError("Codomain of map at edge '" + str(x[2]) + "' does not match.")

    def _repr_(self):
        """
        Default string representation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: Q.P(GF(3), 2) # indirect doctest
            Representation with dimension vector (0, 1)
        """

        return "Representation with dimension vector " + str(self.dimension_vector())

    def __div__(self, sub):
        """
        This and :meth:`__Truediv__` below together overload the /
        operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(GF(3), 1)
            sage: R = P.radical()
            sage: (P/R).is_simple()
            True
        """
        return self.quotient(sub)

    def __truediv__(self, sub):
        """
        This and :meth:`__div__` above together overload the / operator.

        TESTS::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P = Q.P(GF(3), 1)
            sage: R = P.radical()
            sage: (P/R).is_simple()
            True
        """
        return self.quotient(sub)

    def _submodule(self, spaces={}):
        """
        Returns the submodule specified by the data.

        This differs from self.submodule in that it assumes the data correctly
        specifies a submodule whereas self.submodule returns the smallest submodule
        containing the data.

        TESTS::

            sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v]) # indirect doctest
            Representation with dimension vector (1, 1, 1)
            sage: M.submodule(spaces={1: QQ^2}) # indirect doctest
            Representation with dimension vector (2, 0, 2)
            sage: M.submodule(spaces={2: QQ^3}) # indirect doctest
            Representation with dimension vector (0, 3, 1)
            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3}) # indirect doctest
            Representation with dimension vector (1, 3, 1)
            sage: M.submodule().is_zero() # indirect doctest
            True
            sage: M.submodule(M.gens()) is M # indirect doctest
            True
        """

        # Add zero submodules
        for v in self._quiver:
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

        # Create edge homomorphisms restricted to the new domains and codomains
        maps = {}
        for e in self._quiver.edges():
            maps[e] = self._maps[e].restrict_domain(spaces[e[0]]).restrict_codomain(spaces[e[1]])

        return self._semigroup.representation(self._base_ring, spaces, maps)

    def _coerce_map_from_(self, domain):
        """
        Return either a QuiverRepHom from ``domain`` to ``self``, or
        ``False``.

        .. NOTE::

            This function simply tries to coerce a map at each vertex and then check if
            the result is a valid homomorphism.  If it is, then that homomorphism is
            returned.  If it is not or if no coercion was possible then it returns
            ``False``.

        INPUT:

        - ``domain`` - a Sage object

        OUTPUT:

        - QuiverRepHom or bool

        TESTS::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: S = M.radical()
            sage: M.coerce_map_from(S) # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices
            sage: (M/S).coerce_map_from(M) # indirect doctest
            Homomorphism of representations of Multi-digraph on 3 vertices

        Here sage coerces a map but the result is not a homomorphism::

            sage: S.coerce_map_from(M) # indirect doctest

        Here sage cannot coerce a map::

            sage: N = Q.P(QQ, 3)
            sage: N.coerce_map_from(M) # indirect doctest
        """

        from sage.quivers.morphism import QuiverRepHom
        # Domain must be a QuiverRep
        if not isinstance(domain, QuiverRep_generic):
            return False

        # Coerce a map at each vertex, return false if it fails
        maps = {}
        for v in self._quiver:
            maps[v] = self._spaces[v].coerce_map_from(domain._spaces[v])
            if maps[v] is None or maps[v] is False:
                return False

        # Create and return the hom, return False if it wasn't valid
        try:
            return domain.hom(maps, self)
        except ValueError:
            return False

    def _Hom_(self, codomain, category):
        """
        This function is used by the coercion model.

        INPUT:

        - ``codomain`` - QuiverRepHom

        - ``category`` - This input is ignored

        EXAMPLES::

            sage: from sage.categories.homset import Hom
            sage: Q = DiGraph({1:{2:['a', 'b']}, 2:{3:['c', 'd']}}).path_semigroup()
            sage: P = Q.P(GF(3), 2)
            sage: S = P/P.radical()
            sage: Hom(P, S) # indirect doctest
            Dimension 1 QuiverHomSpace
        """
        from sage.quivers.homspace import QuiverHomSpace

        if not isinstance(codomain, QuiverRep_generic):
            raise TypeError("codomain must be of type QuiverRep_generic")

        return QuiverHomSpace(self, codomain)

    ###########################################################################
    #                                                                         #
    # ACCESS FUNCTIONS                                                        #
    #    These functions are used to view the representation data.            #
    #                                                                         #
    ###########################################################################

    def get_space(self, vertex):
        """
        Returns the module associated to the given vertex.

        INPUT:

        - ``vertex`` - integer, a vertex of the quiver of the module

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}}).path_semigroup()
            sage: Q.P(QQ, 1).get_space(1)
            Vector space of dimension 1 over Rational Field
        """

        return self._spaces[vertex]

    def get_map(self, edge):
        """
        Returns the map associated to the given edge.

        INPUT:

        - ``edge`` - tuple of the form (initial vertex, terminal vertex, label) specifying
          the edge whose map is returned

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.P(ZZ, 1).get_map((1, 2, 'a'))
            Free module morphism defined by the matrix
            [1 0]
            Domain: Ambient free module of rank 1 over the principal ideal domain Integer Ring
            Codomain: Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """

        return self._maps[edge]

    def quiver(self):
        """
        Returns the quiver of the representation.

        OUTPUT:

        - Quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.representation(GF(5))
            sage: M.quiver() is Q.quiver()
            True
        """

        return self._quiver

    def base_ring(self):
        """
        Returns the base ring of the representation.

        OUTPUT:

        - ring

        .. TODO::

            Get rid of this method. The convention for modules in Sage is
            that the base ring is the same as the base of the module.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: M = Q.representation(GF(5))
            sage: M.base_ring() is GF(5)
            True
        """

        return self._base_ring

    ###########################################################################
    #                                                                         #
    # DATA FUNCTIONS                                                          #
    #    These functions return data collected from the representation.       #
    #                                                                         #
    ###########################################################################

    def dimension(self, vertex=None):
        """
        Returns the dimension of the space associated to the given vertex
        ``vertex``.

        INPUT:

        - ``vertex`` - integer or ``None`` (default: ``None``), the given
          vertex

        OUTPUT:

        - integer, the dimension over the base ring of the space
          associated to the given vertex.  If ``vertex=None`` then the
          dimension over the base ring of the module is returned

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension(1)
            1
            sage: P.dimension(2)
            2

        The total dimension of the module is the sum of the dimensions
        at each vertex::

            sage: P.dimension()
            3
        """

        if vertex is None:
            # Sum the dimensions of each vertex
            dim = 0
            for x in self._quiver:
                dim += self._spaces[x].dimension()
            return dim
        else:
            # Return the dimension of just the one vertex
            return self._spaces[vertex].dimension()

    def dimension_vector(self):
        """
        Return the dimension vector of the representation.

        OUTPUT:

        - tuple

        .. NOTE::

            The order of the entries in the tuple matches the order given by calling
            the vertices() method on the quiver.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: P = Q.P(GF(2), 1)
            sage: P.dimension_vector()
            (1, 2)

        Each coordinate of the dimension vector is the dimension of the space
        associated to that coordinate::

            sage: P.get_space(2).dimension()
            2
        """

        return tuple(self._spaces[x].dimension() for x in self._quiver)

    def is_zero(self):
        """
        Tests whether the representation is zero.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(ZZ)
            sage: N = Q.representation(ZZ, {1: 1})
            sage: M
            Representation with dimension vector (0, 0)
            sage: N
            Representation with dimension vector (1, 0)
            sage: M.is_zero()
            True
            sage: N.is_zero()
            False
        """

        return self.dimension() == 0

    def is_simple(self):
        """
        Tests whether the representation is simple.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.P(RR, 1).is_simple()
            False
            sage: Q.S(RR, 1).is_simple()
            True
        """

        # A module for an acyclic quiver is simple if and only if it has total
        # dimension 1.
        return self.dimension() == 1

    def is_semisimple(self):
        """
        Tests whether the representation is semisimple.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: (M/M.radical()).is_semisimple()
            True
        """

        # A quiver representation is semisimple if and only if the zero map is
        # assigned to each edge.
        for x in self._quiver.edges():
            if not self._maps[x].is_zero():
                return False
        return True

    def an_element(self):
        """
        Return an element of ``self``.

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.an_element()
            Element of quiver representation
        """

        # Here we just use the an_element function from each space.
        elements = dict((v, self._spaces[v].an_element()) for v in self._quiver)
        return self(elements)

    def support(self):
        """
        Returns the support of self as a list.

        OUTPUT:

        - list, the vertices of the representation that have nonzero spaces associated
          to them

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a']}, 3:{2:['b'], 4:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 3)
            sage: M
            Representation with dimension vector (0, 1, 1, 1)
            sage: M.support()
            [2, 3, 4]
        """

        sup = []
        for v in self._quiver:
            if self._spaces[v].dimension() != 0:
                sup.append(v)

        return sup

    def gens(self, names='v'):
        """
        Returns a list of generators.

        INPUT:

        - ``names`` - an iterable variable of length equal to the number
          of generators, or a string (default: 'v'); gives the names of
          the generators either by giving a name to each generator or by
          giving a name to which an index will be appended

        OUTPUT:

        - list of QuiverRepElement objects, the linear generators of the module (over
          the base ring)

        .. NOTES::

            The generators are ordered first by vertex and then by the order given by
            the gens() method of the space associated to that vertex.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.gens()
            [v_0, v_1, v_2]

        If a string is given then it is used as the name of each generator, with the
        index of the generator appended in order to differentiate them::

            sage: M.gens('generator')
            [generator_0, generator_1, generator_2]

        If a list or other iterable variable is given then each generator is named
        using the appropriate entry.  The length of the variable must equal the number
        of generators (the dimension of the module)::

            sage: M.gens(['w', 'x', 'y', 'z'])
            Traceback (most recent call last):
            ...
            TypeError: can only concatenate list (not "str") to list
            sage: M.gens(['x', 'y', 'z'])
            [x, y, z]

        Strings are iterable, so if the length of the string is equal to the number of
        generators then the characters of the string will be used as the names::

            sage: M.gens('xyz')
            [x, y, z]
        """

        # Use names as a list if and only if it is the correct length
        uselist = (len(names) == self.dimension())
        i = 0

        # Create bases for each space and construct QuiverRepElements from
        # them.  Other functions depend on the ordering of generators produced
        # by this code.  Make sure you know which they are before you change
        # anything
        basis = []
        for v in self._quiver:
            for m in self._spaces[v].gens():
                if uselist:
                    basis.append(self({v: m}, names[i]))
                else:
                    basis.append(self({v: m}, names + "_" + str(i)))
                i += 1

        return basis

    def coordinates(self, vector):
        """
        Returns the coordinates when ``vector`` is expressed in terms of
        the gens.

        INPUT:

        - ``vector`` - QuiverRepElement

        OUTPUT:

        - list, the coefficients when the vector is expressed as a linear combination
          of the generators of the module

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: M.coordinates(x - y + z)
            [1, -1, 1]
            sage: M.coordinates(M.an_element())
            [1, 1, 0]
            sage: M.an_element() == x + y
            True
        """

        # Just use the coordinates functions on each space
        coords = []
        for v in self._quiver:
            coords += self._spaces[v].coordinates(vector._elems[v])

        return coords

    def linear_combination_of_basis(self, coordinates):
        """
        Return the linear combination of the basis for ``self`` given
        by ``coordinates``.

        INPUT:

        - ``coordinates`` - list, a list whose length is the dimension of
          ``self``.  The `i`-th element of this list defines the
          coefficient of the `i`-th basis vector in the linear
          combination.

        OUTPUT:

        - QuiverRepElement

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: x, y, z = M.gens('xyz')
            sage: e = x - y + z
            sage: M.coordinates(e)
            [1, -1, 1]
            sage: M.linear_combination_of_basis([1, -1, 1]) == e
            True
        """

        # Make sure the input is valid
        gens = self.gens()
        if len(gens) != len(coordinates):
            raise ValueError("The coordinates do not match the dimension of the module.")

        result = self()  # this must not be self.zero(), which is cached
        for i in range(0, len(gens)):
            result += coordinates[i]*gens[i]

        return result

    ###########################################################################
    #                                                                         #
    # CONSTRUCTION FUNCTIONS                                                  #
    #    These functions create and return submodules and homomorphisms.      #
    #                                                                         #
    ###########################################################################

    def submodule(self, elements=[], spaces=None):
        """
        Returns the submodule generated by the data.

        INPUT:

        - ``elements`` -- a collection of QuiverRepElements (default:
          empty list), each should be an element of ``self``

        - ``spaces`` -- dictionary (default: empty), this dictionary
          should contain entries of the form ``{v: S}`` where `v` is a
          vertex of the quiver and `S` is a subspace of the vector space
          associated to `v`

        OUTPUT:

        - QuiverRep, the smallest subrepresentation of ``self``
          containing the given elements and the given subspaces

        .. NOTE::

            This function returns only a QuiverRep object ``sub``.  The inclusion map
            of ``sub`` into ``M``=self can be obtained by calling ``M.coerce_map_from(sub)``.

        EXAMPLES::

            sage: Q = DiGraph({1:{3:['a']}, 2:{3:['b']}}).path_semigroup()
            sage: spaces = {1: QQ^2, 2: QQ^3, 3: QQ^2}
            sage: maps = {(1, 3, 'a'): (QQ^2).Hom(QQ^2).identity(), (2, 3, 'b'): [[1, 0], [0, 0], [0, 0]]}
            sage: M = Q.representation(QQ, spaces, maps)
            sage: v = M.an_element()
            sage: M.submodule([v])
            Representation with dimension vector (1, 1, 1)

        The smallest submodule containing the vector space at vertex 1 also contains
        the image of the rank 1 homomorphism associated to the edge (1, 3, 'a')::

            sage: M.submodule(spaces={1: QQ^2})
            Representation with dimension vector (2, 0, 2)

        The smallest submodule containing the vector space at vertex 2 also contains
        the entire vector space associated to vertex 3 because there is an isomorphism
        associated to the edge (2, 3, 'b')::

            sage: M.submodule(spaces={2: QQ^3})
            Representation with dimension vector (0, 3, 1)

        As v is not already contained in this submodule adding it as a generator yields
        a larger submodule::

            sage: v.support()
            [1, 2, 3]
            sage: M.submodule([v], {2: QQ^3})
            Representation with dimension vector (1, 3, 1)

        Giving no generating data yields the zero submodule::

            sage: M.submodule().is_zero()
            True

        If the given data generates all of M then the result is M::

            sage: M.submodule(M.gens()) is M
            True
        """

        if spaces is None:
            spaces = {}

        # For each vertex generate a submodule from the given data
        dim = old_dim = 0
        for v in self._quiver:
            #Start with the zero submodule if no space is specified
            if v not in spaces:
                spaces[v] = self._spaces[v].zero_submodule()

            # Sum this with the submodule generated by the given elements.
            # Note that we are only taking the part of the element at the
            # vertex v.  We can always multiply an element of a quiver
            # representation by a constant path so we don't need to worry about
            # subspaces being embedded diagonally across multiple vertices.
            spaces[v] += self._spaces[v].submodule([m._elems[v] for m in elements])
            dim += spaces[v].dimension()

        #Now to enlarge the subspace to a submodule we sum a subspace at a
        # vertex with the images of the subspaces at adjacent vertices.  The
        # dimension of the subspace will strictly increase until we generate a
        # submodule.  At that point the dimension stabilizes and we can exit
        # the loop.
        while old_dim != dim:
            old_dim, dim = dim, 0

            # First sum the subspaces
            for e in self._quiver.edges():
                spaces[e[1]] += self._maps[e](spaces[e[0]])

            # Then get the resulting dimensions
            for v in self._quiver:
                dim += spaces[v].dimension()

        # Return self if the entire module was generated, otherwise return a
        # submodule
        if dim == self.dimension():
            return self
        else:
            return self._submodule(spaces)

    def quotient(self, sub, check=True):
        """
        Return the quotient of ``self`` by the submodule ``sub``.

        INPUT:

        - ``sub`` -- QuiverRep; this must be a submodule of ``self``,
          meaning the space associated to each vertex `v` of ``sub`` is a
          subspace of the space associated to `v` in ``self`` and the map
          associated to each edge `e` of ``sub`` is the restriction of
          the map associated to `e` in ``self``

        - ``check`` -- bool; if ``True`` then ``sub`` is checked to verify
          that it is indeed a submodule of ``self`` and an error is raised
          if it is not

        OUTPUT:

        - QuiverRep, the quotient module ``self``/``sub``

        .. NOTE::

            This function returns only a QuiverRep object ``quot``.  The
            projection map from ``self`` to ``quot`` can be obtained by
            calling ``quot.coerce_map_from(self)``.

        EXAMPLES:

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.I(GF(3), 3)
            sage: N = Q.S(GF(3), 3)
            sage: M.quotient(N)
            Representation with dimension vector (2, 1, 0)
            sage: M.quotient(M.radical())
            Representation with dimension vector (2, 0, 0)
            sage: M.quotient(M)
            Representation with dimension vector (0, 0, 0)
        """

        # First form the quotient space at each vertex
        spaces = {}
        for v in self._quiver:
            spaces[v] = self._spaces[v].quotient(sub._spaces[v], check)

        # Check the maps of sub if desired
        if check:
            for e in self._quiver.edges():
                for x in sub._spaces[e[0]].gens():
                    if sub._maps[e](x) != self._maps[e](x):
                        raise ValueError("The quotient method was not passed a submodule.")

        # Then pass the edge maps to the quotient
        maps = {}
        for e in self._quiver.edges():
            # Sage can automatically coerce an element of a module to an
            # element of a quotient of that module but not the other way
            # around.  So in order to pass a map to the quotient we need to
            # construct the quotient map for the domain so that we can take
            # inverse images to lift elments.  As sage can coerce to a quotient
            # this is easy, we just send generators to themselves and set the
            # domain to be the quotient.

            # TODO: This 'if' shouldn't need to be here, but sage crashes when
            # coercing elements into a quotient that is zero.  Once Trac ticket
            # 12413 gets fixed only the else should need to execute.
            # NOTE: This is no longer necessary. Keep around for speed or
            # remove? -- darij, 16 Feb 2014
            if spaces[e[1]].dimension() == 0:
                maps[e] = spaces[e[0]].Hom(spaces[e[1]]).zero_element()
            else:
                factor = self._spaces[e[0]].hom(self._spaces[e[0]].gens(), spaces[e[0]])
                # Now we create a homomorphism by specifying the images of
                # generators.  Each generator is lifted to the original domain and
                # mapped over using the original map.  The codomain is set as the
                # quotient so sage will take care of pushing the result to the
                # quotient in the codomain.
                maps[e] = spaces[e[0]].hom([self._maps[e](factor.lift(x)) for x in spaces[e[0]].gens()], spaces[e[1]])

        return self._semigroup.representation(self._base_ring, spaces, maps)

    def socle(self):
        """
        The socle of ``self``.

        OUTPUT:

        - QuiverRep, the socle

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.socle()
            Representation with dimension vector (0, 0, 2)
        """

        # The socle of a representation is the intersection of the kernels of
        # all the edge maps.  The empty intersection is defined to be the
        # entire space so this is what we start with.
        spaces = self._spaces.copy()
        for e in self._quiver.edges():
            spaces[e[0]] = spaces[e[0]].intersection(self._maps[e].kernel())

        return self._submodule(spaces)

    def radical(self):
        """
        Returns the Jacobson radical of self.

        OUTPUT:

        - QuiverRep, the socle

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.radical()
            Representation with dimension vector (0, 2, 2)
        """

        #The Jacobson radical of a representation is the sum of the images of
        # all of the edge maps.  The empty sum is defined to be zero so this is
        # what we start with.
        spaces = dict((v, self._spaces[v].zero_submodule()) for v in self._quiver)
        for e in self._quiver.edges():
            spaces[e[1]] += self._maps[e].image()

        return self._submodule(spaces)

    def top(self):
        """
        Returns the top of self.

        OUTPUT:

        - QuiverRep, the quotient of self by its radical

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.top()
            Representation with dimension vector (1, 0, 0)
            sage: M.top() == M/M.radical()
            True
        """

        return self.quotient(self.radical())

    def zero_submodule(self):
        """
        Returns the zero submodule.

        OUTPUT:

        - QuiverRep, the quotient of self by its radical

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.zero_submodule()
            Representation with dimension vector (0, 0, 0)
            sage: M.zero_submodule().is_zero()
            True
        """

        # When no data is specified this constructor automatically returns the
        # zero submodule
        return self._submodule()

    def linear_dual(self):
        """
        Computes the linear dual Hom_k(M, k) of the module M=self over the base ring k.

        OUTPUT:

        - QuiverRep, the dual representation

        .. NOTES::

            If e is an edge of the quiver Q then we let (fe)(m) = f(me).  This gives
            Hom_k(M, k) a module structure over the opposite quiver Q.reverse().

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: M.linear_dual()
            Representation with dimension vector (1, 2, 2)
            sage: M.linear_dual().quiver() is Q.reverse().quiver()
            True
        """

        # This module is formed by taking the transpose of the edge maps.
        spaces = self._spaces.copy()
        maps = dict(((e[1], e[0], e[2]), self._spaces[e[1]].hom(self._maps[e].matrix().transpose(), self._spaces[e[0]]))
                    for e in self._quiver.edges())

        # Reverse the bases if present
        if hasattr(self, '_bases'):
            bases = {}
            basis = []
            for v in self._bases:
                bases[v] = [p.reverse() for p in self._bases[v]]
                basis.extend(bases[v])

        if isinstance(self, QuiverRep_with_path_basis):
            result = self._semigroup.reverse().representation(self._base_ring, basis, option='dual paths')
            result._maps = maps
            result._bases = bases
        elif isinstance(self, QuiverRep_with_dual_path_basis):
            result = self._semigroup.reverse().representation(self._base_ring, basis, option='paths')
            result._maps = maps
            result._bases = bases
        else:
            result = self._semigroup.reverse().representation(self._base_ring, spaces, maps)

        return result

    def algebraic_dual(self, basis=False):
        """
        Computes the algebraic dual Hom_Q(M, kQ) of the module M=self.

        INPUT:

        - ``basis`` - bool, if false then only the module is returned.  If true then a
          tuple is returned.  The first element is the QuiverRep and the second element
          is a dictionary which associates to each vertex a list.  The elements of this
          list are the homomorphisms which correspond to the basis elements of that
          vertex in the module.

        OUTPUT:

        - QuiverRep or tuple

        .. NOTES::

            Here kQ is the path algebra considered as a right module over itself.  If e
            is an edge of the quiver Q then we let (fe)(m) = ef(m).  This gives
            Hom_Q(M, kQ) a module structure over the opposite quiver Q.reverse().

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b'], 3: ['c', 'd']}, 2:{3:['e']}}).path_semigroup()
            sage: Q.free_module(GF(7)).algebraic_dual().dimension_vector()
            (7, 2, 1)
        """
        from sage.quivers.homspace import QuiverHomSpace

        return QuiverHomSpace(self, self._semigroup.free_module(self._base_ring)).left_module(basis)

    def Hom(self, codomain):
        """
        Returns the hom space from self to codomain.

        For more information see the QuiverHomSpace documentation.

        TESTS::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: Q.S(QQ, 2).Hom(Q.P(QQ, 1))
            Dimension 2 QuiverHomSpace
        """
        from sage.quivers.homspace import QuiverHomSpace

        return QuiverHomSpace(self, codomain)

    def direct_sum(self, modules, return_maps=False):
        """
        Returns the direct sum of self with the given modules.

        The modules must be modules over the same quiver and base ring.

        INPUT:

        - ``modules`` - QuiverRep or list of QuiverReps

        - ``return_maps`` - Boolean (default: False), if False then the output is a single
          QuiverRep object which is the direct sum of self with the given module or
          modules.  If True then the output is a list ``[sum, iota, pi]``.  The
          first entry ``sum`` is the direct sum of self with the given module or
          modules.  Both ``iota`` and ``pi`` are lists of QuiverRepHoms with one entry for
          each summand; ``iota[i]`` is the inclusion map and ``pi[i]`` is the projection
          map of the ith summand.  The summands are ordered as given with self
          being the zeroth summand.

        OUTPUT:

        - QuiverRep or tuple

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a'], 3:['b']}, 2:{4:['c']}, 3:{4:['d']}}).path_semigroup()
            sage: P1 = Q.P(QQ, 1)
            sage: P2 = Q.P(QQ, 2)
            sage: S = P1.direct_sum(P2)
            sage: P1.dimension_vector()
            (1, 1, 1, 2)
            sage: P2.dimension_vector()
            (0, 1, 0, 1)
            sage: S.dimension_vector()
            (1, 2, 1, 3)
            sage: S, iota, pi = P1.direct_sum(P2, return_maps=True)
            sage: iota[0].domain() is P1
            True
            sage: iota[1].domain() is P2
            True
            sage: pi[0].codomain() is P1
            True
            sage: pi[1].codomain() is P2
            True
            sage: iota[0].codomain() is S
            True
            sage: iota[1].codomain() is S
            True
            sage: pi[0].domain() is S
            True
            sage: pi[1].domain() is S
            True
            sage: iota[0].get_matrix(4)
            [1 0 0]
            [0 1 0]
            sage: pi[0].get_matrix(4)
            [1 0]
            [0 1]
            [0 0]
            sage: P1prime = S/iota[1].image()
            sage: f = P1prime.coerce_map_from(S)
            sage: g = f*iota[0]
            sage: g.is_isomorphism()
            True
        """

        from sage.quivers.morphism import QuiverRepHom
        # Convert a single module into a length 1 list and check validity
        if isinstance(modules, QuiverRep_generic):
            mods = [self, modules]
        else:
            mods = [self] + list(modules)
        for i in range(1, len(mods)):
            if not isinstance(mods[i], QuiverRep_generic):
                raise ValueError("modules must be a QuiverRep or a list of QuiverReps")
            if self._quiver is not mods[i]._quiver:
                raise ValueError("Cannot direct sum modules of different quivers")
            if self._base_ring is not mods[i]._base_ring:
                raise ValueError("Cannot direct sum modules with different base rings")

        # Get the dimensions of all spaces at each vertex
        dims = dict((v, [x._spaces[v].dimension() for x in mods]) for v in self._quiver)

        # Create spaces of the correct dimensions
        spaces = dict((v, self._base_ring**sum(dims[v])) for v in self._quiver)

        # Take block sums of matrices to form the maps
        from sage.matrix.constructor import block_diagonal_matrix
        maps = {}
        for e in self._quiver.edges():
            maps[e] = block_diagonal_matrix([x._maps[e].matrix() for x in mods], subdivide=False)

        # Create the QuiverRep, return if the maps aren't wanted
        result = self._semigroup.representation(self._base_ring, spaces, maps)
        if not return_maps:
            return result

        # Create the inclusions and projections
        from sage.matrix.constructor import block_matrix, Matrix
        from sage.rings.integer import Integer
        iota = []
        pi = []
        for i in range(0, len(mods)):
            incl_maps = {}
            proj_maps = {}
            for v in self._quiver:
                # Create the maps using block matrices
                pre_dims = sum(dims[v][:i])
                post_dims = sum(dims[v][i + 1:])
                incl_maps[v] = block_matrix(1, 3, [Matrix(dims[v][i], pre_dims),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(dims[v][i], post_dims)])
                proj_maps[v] = block_matrix(3, 1, [Matrix(pre_dims, dims[v][i]),
                                                   Matrix(dims[v][i], dims[v][i], Integer(1)),
                                                   Matrix(post_dims, dims[v][i])])

            # Create and save the QuiverRepHom

            iota.append(mods[i].hom(incl_maps, result))
            pi.append(result.hom(proj_maps, mods[i]))

        # Return all the data
        return [result, iota, pi]

    def projective_cover(self, return_maps=False):
        """
        Return the projective cover of self.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c','d']}}).path_semigroup()
            sage: S1 = Q.S(GF(3), 1)
            sage: f1 = S1.projective_cover()
            sage: f1.is_surjective()
            True
            sage: f1._domain
            Representation with dimension vector (1, 2, 4)
            sage: Q.P(GF(3), 1)
            Representation with dimension vector (1, 2, 4)
        """

        # Deal with the zero module
        if self.dimension() == 0:
            return self.coerce_map_from(self)

        # Get the projection onto the top and generators
        top = self.top()
        proj_to_top = top.coerce_map_from(self)
        gens = top.gens()

        # Lift each of the generators
        lifts = [proj_to_top.lift(x) for x in gens]

        # Get projective summands of the cover
        Ps = [self._semigroup.P(self._base_ring, x.support()[0]) for x in gens]

        # Factor the maps through self
        maps = [Ps[i].hom(lifts[i], self) for i in range(0, len(gens))]

        # Sum them and return
        return maps[0].direct_sum(maps[1:], return_maps, 'codomain')

    def transpose(self):
        """
        Returns the transpose of self.

        The transpose, Tr M, of a module M is defined as follows.  Let p: P1 -> P2 be
        the second map in a minimal projective presentation P1 -> P2 -> M -> 0 of M.
        If p^t is the algebraic dual of p then define Tr M = coker p^t.

        OUTPUT:

        - QuiverRep

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(GF(3), {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: M.transpose()
            Representation with dimension vector (1, 1)
        """

        # Get the second map in the minimal projective resolution of self.
        # Note we create this map as the projective cover of a kernel so we
        # need to alter the codomain to be the projective and not the kernel.
        p1 = self.projective_cover()
        k = p1.kernel()
        p = p1.domain().coerce_map_from(k)*k.projective_cover()

        # Return the cokernel
        return p.algebraic_dual().cokernel()

    def AR_translate(self):
        """
        Return the Auslander-Reiten translate of self.

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a', 'b']}}).path_semigroup()
            sage: M = Q.representation(GF(3), {1: 1, 2: 1}, {(1, 2, 'a'): 1})
            sage: tauM = M.AR_translate()
            sage: tauM
            Representation with dimension vector (1, 1)
            sage: tauM.get_map((1, 2, 'a')).matrix()
            [1]
            sage: tauM.get_map((1, 2, 'b')).matrix()
            [0]

        The module M above is its own AR translate.  This is not always true::

            sage: Q2 = DiGraph({3:{1:['b']}, 5:{3:['a']}}).path_semigroup()
            sage: Q2.S(QQ, 5).AR_translate()
            Representation with dimension vector (0, 1, 0)
        """

        return self.transpose().linear_dual()

    ###########################################################################
    #                                                                         #
    # ADDITIONAL OPERATIONS                                                   #
    #    These functions operations that are not implemented via binary       #
    #    operators.                                                           #
    #                                                                         #
    ###########################################################################

    def right_edge_action(self, element, path):
        """
        Returns the result of element*edge.

        INPUT:

        - ``element`` - QuiverRepElement, an element of self

        - ``path`` - QuiverPath or list of tuples

        OUTPUT:

        - QuiverRepElement, the result of element*path when path is considered an
          element of the path algebra of the quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.P(QQ, 1)
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]
            sage: M.right_edge_action(v, (1, 1)).support()
            [1]
            sage: M.right_edge_action(v, [(1, 1)]).support()
            [1]
            sage: M.right_edge_action(v, [(1, 1), (1, 2, 'a')]).support()
            [2]
            sage: M.right_edge_action(v, (1, 2, 'a')) == M.right_edge_action(v, [(1, 1), (1, 2, 'a'), (2, 2)])
            True
        """

        # Convert to a QuiverPath
        qpath = self._semigroup(path)

        # Invalid paths are zero in the quiver algebra
        result = self()  # this must not be self.zero(), which is cached
        if not qpath:
            return result

        # Start with the element at the initial vertex
        x = element._elems[qpath.initial_vertex()]

        # Act by each edge
        for e in qpath:
            x = self.get_map(e)(x)

        # Assign the result to the terminal vertex and return
        result._elems[qpath.terminal_vertex()] = x
        return result

class QuiverRep_with_path_basis(QuiverRep_generic):
    """
    The basis of the module must be closed under right multiplication by an
    edge; that is, appending any edge to the end of any path in the basis must
    result in either an invalid path or a valid path also contained in the
    basis of the module.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``basis`` - list (default: empty), should be a list of paths (also lists) in the
      quiver Q.  Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  The closure of this list under right multiplication forms the
      basis of the resulting representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are path multiplication.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the codomain.  This class derives
    # from the QuiverRep class and the following private methods and variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths (also lists)
    #       which correspond to the basis elements of the space assigned to that
    #       vertex.
    #
    # If the right closure of the basis is also closed under left mult by an
    # edge then the object will have the following private variable:
    #
    # * _left_action_mats
    #       A dictionary of dictionaries.  If e is an edge or trivial path of
    #       the quiver then _left_action_mats[e] is a dictionary that
    #       associates a matrix to each vertex.  The left action of e on a
    #       QuiverRepElement is given by multiplying, on the left, the vector
    #       associated to vertex v in the QuiverRepElement by the matrix
    #       _left_action_mats[e][v]

    def __init__(self, k, P, basis):
        """
        Type QuiverRep_with_path_basis? for more information.

        TESTS::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P1 = Q1.representation(QQ, [[(1, 1)]], option='paths')
            sage: P1.dimension()
            2
            sage: kQ = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='paths')
            sage: kQ.dimension()
            3
            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.dimension_vector()
            (0, 2, 2)
            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.dimension_vector()
            (0, 1, 2)
        """

        self._quiver = Q = P.quiver()
        self._base_ring = k

        # Add the paths to the basis dictionary.  The terminal vertex is the
        # key
        self._bases = dict((v, []) for v in Q)
        for path in basis:
            if path:
                self._bases[path.terminal_vertex()].append(path)

        # Create the matrices of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in Q.edges():
            arrow = P(e)
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(self._base_ring, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(0, len(self._bases[e[0]])):
                # Add an entry to the matrix coresponding to where the new path is found
                j = self._bases[e[1]].index(self._bases[e[0]][i]*arrow)
                maps[e][i, j] = self._base_ring(1)

        # Create the spaces and then the representation
        spaces = dict((v, len(self._bases[v])) for v in Q)
        super(QuiverRep_with_path_basis, self).__init__(k, P, spaces, maps)

        # Try and create the matrices for the left edge action of edges.  If it
        # fails just return, there's no edge action and the construction is
        # done
        action_mats = {}
        for e in self._quiver.edges():
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(self._base_ring, l, l)

                for j in range(0, l):
                    if e[1] == self._bases[v][j].initial_vertex():
                        try:
                            action_mats[e][v][self._bases[v].index(P(e)*self._bases[v][j]), j] = self._base_ring(1)
                        except ValueError:
                            # There is no left action
                            return

        # If haven't returned yet then there is a left action.  Create the
        # matrices for acting by trivial paths
        for vert in self._quiver:
            e = (vert, vert)
            action_mats[e] = {}
            for v in self._quiver:
                # Start with the zero matrix and fill in
                l = len(self._bases[v])
                action_mats[e][v] = Matrix(self._base_ring, l, l)

                # Paths not beginning at vert are sent to zero, paths beginning
                # at vert are fixed
                for i in range(0, l):
                    if self._bases[v][i].initial_vertex() == vert:
                        action_mats[e][v][i, i] = self._base_ring(1)

        # Define the method and save the matrices
        self.left_edge_action = self._left_edge_action
        self._left_action_mats = action_mats

    def _left_edge_action(self, edge, element):
        """
        Returns the result of edge*element.

        INPUT:

        - ``element`` - QuiverRepElement, an element of self

        - ``edge`` - an edge of the quiver (a tuple) or a list of edges in the quiver.
          Such a list can be empty (in which case no action is performed) and can
          contain trivial paths (tuples of the form (v, v) where v is a vertex of the
          quiver)

        OUTPUT:

        - QuiverRepElement, the result of edge*element when edge is considered an
          element of the path algebra of the quiver

        EXAMPLES::

            sage: Q = DiGraph({1:{2:['a','b']}, 2:{3:['c']}}).path_semigroup()
            sage: M = Q.representation(QQ, [[(1, 1)], [(2, 2)], [(3, 3)]], option='paths')
            sage: v = M.an_element()
            sage: v.support()
            [1, 2, 3]

        The sum of all trivial paths is the identity element under this action::

            sage: x = M.left_edge_action((1, 1), v) # indirect doctest
            sage: y = M.left_edge_action((2, 2), v) # indirect doctest
            sage: z = M.left_edge_action((3, 3), v) # indirect doctest
            sage: x + y + z == v
            True

        Note that the action only depends on the element of the path algebra
        that the input specifies::

            sage: a = M.left_edge_action([(1, 1), (1, 2, 'a'), (2, 2)], v) # indirect doctest
            sage: b = M.left_edge_action((1, 2, 'a'), v) # indirect doctest
            sage: a == b
            True

        These two edges multiply to zero in the path algebra::

            sage: M.left_edge_action([(1, 2, 'a'), (1, 2, 'b')], v).is_zero() # indirect doctest
            True
        """

        # Deal with lists by calling this function recursively
        if isinstance(edge, list):
            if len(edge) == 0:
                return element;
            else:
                return self.left_edge_action(edge[:-1], self.left_edge_action(edge[-1], element))

        # Now we are just acting by a single edge
        elems = dict((v, self._left_action_mats[edge][v]*element._elems[v]) for v in self._quiver)
        return self(elems)

    def is_left_module(self):
        """
        Tests whether the basis is closed under left multiplication.

        EXAMPLES::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: P2 = Q1.representation(QQ, [[(2, 2)]], option='paths')
            sage: P2.is_left_module()
            False

        The supplied basis is not closed under left multiplication, but it's not closed
        under right multiplication either.  When the closure under right multiplication
        is taken the result is also closed under left multiplication and therefore
        produces a left module structure::

            sage: kQ = Q1.representation(QQ, [[(1, 1)], [(2, 2)]], option='paths')
            sage: kQ.is_left_module()
            True

        Taking the right closure of a left closed set produces another left closed set::

            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a')]], option='paths')
            sage: M.is_left_module()
            True

        Note that the second path is length 2, so even though the edge (1, 2, 'a')
        appears in the input the path [(1, 2, 'a')] is not in the right closure::

            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='paths')
            sage: N.is_left_module()
            False
        """

        return hasattr(self, 'left_edge_action')

class QuiverRep_with_dual_path_basis(QuiverRep_generic):
    """
    The basis of the module must be closed under left deletion of an edge; that
    is, deleting any edge from the beginning of any path in the basis must
    result in a path also contained in the basis of the module.

    INPUT:

    - ``k`` - ring, the base ring of the representation

    - ``Q`` - Quiver, the quiver of the representation

    - ``basis`` - list (default: empty), should be a list of paths (also lists) in the
      quiver Q.  Entries that do not represent valid paths are ignored and duplicate
      paths are deleted.  The closure of this list under left deletion forms the
      basis of the resulting representation.
    """
    # This class implements quiver representations whose bases correspond to
    # paths in the path algebra and whose maps are edge deletion.  The
    # main advantage to having such a basis is that a homomorphism can be
    # defined by giving a single element in the domain.  This class derives
    # from the QuiverRep class and the following methods and private variables
    # have been added:
    #
    # * _bases
    #       A dictionary associating to each vertex a list of paths which
    #       correspond to the basis elements of the space assigned to that
    #       vertex.

    def __init__(self, k, P, basis):
        """
        Type QuiverRep_with_dual_path_basis? for more information.

        TESTS::

            sage: Q1 = DiGraph({1:{2:['a']}}).path_semigroup()
            sage: I2 = Q1.representation(QQ, [(2, 2)], option='dual paths')
            sage: I2.dimension()
            2
            sage: kQdual = Q1.representation(QQ, [[(1, 1)], [(2, 2)], [(1, 1), (1, 2, 'a'), (2, 2)], [(1, 2, 'a')]], option='dual paths')
            sage: kQdual.dimension()
            3
            sage: Q2 = DiGraph({1:{2:['a'], 3:['b', 'c']}, 2:{3:['d']}}).path_semigroup()
            sage: M = Q2.representation(QQ, [[(1, 2, 'a'), (2, 3, 'd')], [(1, 3, 'b')]], option='dual paths')
            sage: M.dimension_vector()
            (2, 0, 0)
            sage: N = Q2.representation(QQ, [[(2, 2)], [(1, 2, 'a'), (2, 3, 'd')]], option='dual paths')
            sage: N.dimension_vector()
            (2, 1, 0)
        """

        self._quiver = Q = P.quiver()
        self._base_ring = k

        # Add the paths to the basis dictionary.  The initial vertex is the
        # key
        self._bases = dict((v, []) for v in Q)
        for path in basis:
            if path:
                self._bases[path.initial_vertex()].append(path)

        # Create the matrices of the maps
        from sage.matrix.constructor import Matrix
        maps = {}
        for e in Q.edges():
            arrow = P(e)
            # Start with the zero matrix and fill in from there
            maps[e] = Matrix(self._base_ring, len(self._bases[e[0]]), len(self._bases[e[1]]))
            for i in range(0, len(self._bases[e[0]])):
                # Add an entry to the matrix coresponding to where the new path is found
                if self._bases[e[0]][i] % arrow in self._bases[e[1]]:
                    j = self._bases[e[1]].index(self._bases[e[0]][i] % e)
                    maps[e][i, j] = self._base_ring(1)

        # Create the spaces and then the representation
        spaces = dict((v, len(self._bases[v])) for v in Q)
        super(QuiverRep_with_dual_path_basis, self).__init__(k, P, spaces, maps)

