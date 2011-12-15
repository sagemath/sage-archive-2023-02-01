r"""
Polyhedra

This module contains functions for computing with exact (rational)
or floating-point convex polyhedra in arbitrary dimensions. The bulk
of this functionality is provided through the cddlib library of Komei
Fukuda.

There seems to be some inconsistency in the use of the word polyhedra.
In this module, a polyhedron is a convex (possibly unbounded) set in
Euclidean space cut out by a finite set of linear inequalities and
linear equations. Note that the dimension of the polyhedron can be
less than the dimension of the ambient space. There are two
complementary representations of the same data:

**H(alf-space/Hyperplane)-representation**
    This describes a polyhedron as the common solution set of a
    finite number of

        * linear inequalities `A \vec{x} + b \geq 0`, and
        * linear equations  `C \vec{x} + d \geq 0`.


**V(ertex)-representation**
    The other representation is as the convex hull of vertices (and
    rays and lines to all for unbounded polyhedra) as generators. The
    polyhedron is then the Minkowski sum

    .. MATH::

        P = \text{conv}\{v_1,\dots,v_k\} +
        \sum_{i=1}^m \RR_+ r_i +
        \sum_{j=1}^n \RR \ell_j

    where

        * `v_1`, `\dots`, `v_k` are a finite number of vertices,
        * `r_1`, `\dots`, `r_m` are generators of rays,
        * and `\ell_1`, `\dots`, `\ell_n` are generators of full lines.


A polytope is defined as a bounded polyhedron.

EXAMPLES::

    sage: trunc_quadr = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
    sage: trunc_quadr
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays.
    sage: v = trunc_quadr.vertex_generator().next()  # the first vertex in the internal enumeration
    sage: v
    A vertex at (1, 0)
    sage: v()
    (1, 0)
    sage: list(v)
    [1, 0]
    sage: len(v)
    2
    sage: v[0] + v[1]
    1
    sage: v.is_vertex()
    True
    sage: type(v)
    <class 'sage.geometry.polyhedra.Vertex'>
    sage: type( v() )
    <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
    sage: v.polyhedron()
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices and 2 rays.
    sage: r = trunc_quadr.ray_generator().next()
    sage: r
    A ray in the direction (1, 0)
    sage: r()
    (1, 0)
    sage: [x for x in v.neighbors()]
    [A ray in the direction (1, 0), A vertex at (0, 1)]

Inequalities `A \vec{x} + b \geq 0` (and, similarly, equations) are specified by a list ``[b, A]``::

    sage: Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,-1]]).Hrepresentation()
    [An inequality (1, 0) x + 0 >= 0, An inequality (0, 1) x + 0 >= 0, An inequality (-1, -1) x + 1 >= 0]

In addition to polyhedra, this module provides the function
:func:`Hasse_diagram_from_incidences` for computing Hasse diagrams of
finite atomic and coatomic lattices in the sense of partially ordered
sets where any two elements have meet and joint. For example, the face
lattice of a polyhedron.

REFERENCES:

    Komei Fukuda's `FAQ in Polyhedral Computation
    <http://www.ifor.math.ethz.ch/~fukuda/polyfaq/polyfaq.html>`_

AUTHORS:

    - Marshall Hampton: first version, bug fixes, and various improvements, 2008 and 2009
    - Arnaud Bergeron: improvements to triangulation and rendering, 2008
    - Sebastien Barthelemy: documentation improvements, 2008
    - Volker Braun: refactoring, handle non-compact case, 2009 and 2010
    - Andrey Novoseltsev: added Hasse_diagram_from_incidences, 2010

TESTS::

    TestSuite(s).run()
"""


########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence
from sage.categories.objects import Objects

from subprocess import Popen, PIPE
from sage.misc.all import tmp_filename, cached_method, prod
from sage.misc.functional import norm
from sage.misc.package import is_package_installed

from sage.rings.all import Integer, QQ, ZZ, primes_first_n
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.functions.other import sqrt, floor, ceil
from sage.functions.trig import sin, cos

from sage.plot.all import point2d, line2d, arrow, polygon2d
from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
from sage.graphs.graph import Graph

from sage.combinat.combinat import permutations
from sage.combinat.cartesian_product import CartesianProduct
from sage.groups.perm_gps.permgroup_named import AlternatingGroup




#########################################################################
def _common_length_of(l1, l2=None, l3=None):
    """
    The arguments are containers or ``None``. The function applies
    ``len()`` to each element, and returns the common length. If the
    length differs, ``ValueError`` is raised. Used to check arguments.

    OUTPUT:

    A tuple (number of entries, common length of the entries)

    EXAMPLES::

        sage: import sage.geometry.polyhedra as P
        sage: P._common_length_of([[1,2,3],[1,3,34]])
        (2, 3)
    """
    args = [];
    if l1 != None: args.append(l1)
    if l2 != None: args.append(l2)
    if l3 != None: args.append(l3)

    length = None
    num = 0
    for l in args:
        for i in l:
            num += 1
            length_i = len(i)
            if length!=None and length_i!=length:
                raise ValueError, "Argument lengths differ!"
            length = length_i

    return num, length


#########################################################################
def _set_to_None_if_empty(x):
    """
    Helper function to clean up arguments: Returns None if x==None or
    x is an empty container.

    EXAMPLES::

        sage: import sage.geometry.polyhedra as P
        sage: None == P._set_to_None_if_empty([])
        True
        sage: P._set_to_None_if_empty([1])
        [1]
    """
    if x==None: return x
    if len(x)==0: return None
    return x

#########################################################################
def _to_space_separated_string(l):
    """
    Converts a container to a space-separated string.

    EXAMPLES::

        sage: import sage.geometry.polyhedra as P
        sage: P._to_space_separated_string([2,3])
        '2 3'
    """
    s = '';
    for x in l:
        if len(s)>0: s += ' '
        s += repr(x)
    return s


#########################################################################
def cdd_Vrepresentation(cdd_type, vertices, rays, lines):
    r"""
    Return a string containing the V-representation in cddlib's ext format.

    NOTE:

    If there is no vertex given, then the origin will be implicitly
    added. You cannot write the empty V-representation (which cdd
    would refuse to process).

    EXAMPLES::

        sage: from sage.geometry.polyhedra import cdd_Vrepresentation
        sage: print cdd_Vrepresentation('rational', [[0,0]], [[1,0]], [[0,1]])
        V-representation
        linearity 1 1
        begin
          3 3 rational
          0 0 1
          0 1 0
          1 0 0
        end
    """
    vertices = _set_to_None_if_empty(vertices)
    rays     = _set_to_None_if_empty(rays)
    lines    = _set_to_None_if_empty(lines)

    num, ambient_dim = _common_length_of(vertices, rays, lines)

    # cdd implicitly assumes that the origin is a vertex if none is given
    if vertices==None:
        vertices = [[0]*ambient_dim]
        num += 1

    s = 'V-representation\n'
    if lines!=None:
        n = len(lines)
        s += "linearity " + repr(n) + ' '
        s += _to_space_separated_string(range(1,n+1)) + '\n'
    s += 'begin\n'
    s += ' ' + repr(num) + ' ' + repr(ambient_dim+1) + ' ' + cdd_type + '\n'
    if lines!=None:
        for l in lines:
            s += ' 0 ' + _to_space_separated_string(l) + '\n'
    if rays!=None:
        for r in rays:
            s += ' 0 ' + _to_space_separated_string(r) + '\n'
    if vertices!=None:
        for v in vertices:
            s += ' 1 ' + _to_space_separated_string(v) + '\n'
    s += 'end\n'
    return s


#########################################################################
def cdd_Hrepresentation(cdd_type, ieqs, eqns):
    r"""
    Return a string containing the H-representation in cddlib's ine format.

    EXAMPLES::

        sage: from sage.geometry.polyhedra import cdd_Hrepresentation
        sage: cdd_Hrepresentation('rational', None, [[0,1]])
        'H-representation\nlinearity 1 1\nbegin\n 1 2 rational\n 0 1\nend\n'
    """
    ieqs = _set_to_None_if_empty(ieqs)
    eqns  = _set_to_None_if_empty(eqns)

    num, ambient_dim = _common_length_of(ieqs, eqns)
    ambient_dim -= 1

    s = 'H-representation\n'
    if eqns!=None:
        assert len(eqns)>0
        n = len(eqns)
        s += "linearity " + repr(n) + ' '
        s += _to_space_separated_string(range(1,n+1)) + '\n'
    s += 'begin\n'
    s += ' ' + repr(num) + ' ' + repr(ambient_dim+1) + ' ' + cdd_type + '\n'
    if eqns!=None:
        for e in eqns:
            s += ' ' + _to_space_separated_string(e) + '\n'
    if ieqs!=None:
        for i in ieqs:
            s += ' ' + _to_space_separated_string(i) + '\n'
    s += 'end\n'
    return s


#########################################################################
#                      PolyhedronRepresentation
#                       /                     \
#                      /                       \
#              Hrepresentation            Vrepresentation
#                   /     \                 /     |   \
#                  /       \               /      |    \
#           Inequality  Equation        Vertex   Ray   Line

class PolyhedronRepresentation(SageObject):
    """
    The internal base class for all representation objects of
    ``Polyhedron`` (vertices/rays/lines and inequalites/equations)

    NOTES:

    You should not (and cannot) instantiate it yourself. You can only
    obtain them from a Polyhedron() class.

    """

    # Numeric values for the output of the type() method
    INEQUALITY = 0
    EQUATION = 1
    VERTEX = 2
    RAY = 3
    LINE = 4

    def __init__(self, polyhedron, data):
        """
        Initializes the PolyhedronRepresentation object.

        TESTS::

            sage: import sage.geometry.polyhedra as P
            sage: pr = P.PolyhedronRepresentation(Polyhedron(vertices = [[1,2,3]]), [1,2,3])
            sage: pr.vector()
            (1, 2, 3)
        """
        self._representation_vector = vector(polyhedron.field(), data);
        self._polyhedron = polyhedron;

    def __len__(self):
        """
        Returns the length of the representation data.

        TESTS::

            sage: import sage.geometry.polyhedra as P
            sage: pr = P.PolyhedronRepresentation(Polyhedron(vertices = [[1,2,3]]), [1,2,3])
            sage: pr.__len__()
            3
        """
        return len(self._representation_vector)

    def __getitem__(self, i):
        """
        Supports indexing.

        TESTS::

            sage: import sage.geometry.polyhedra as P
            sage: pr = P.PolyhedronRepresentation(Polyhedron(vertices = [[1,2,3]]), [1,2,3])
            sage: pr.__getitem__(1)
            2
        """
        return self._representation_vector[i]

    def polyhedron(self):
        """
        Returns the underlying polyhedron.

        TESTS::

            sage: import sage.geometry.polyhedra as P
            sage: pr = P.PolyhedronRepresentation(Polyhedron(vertices = [[1,2,3]]), [1,2,3])
            sage: pr.polyhedron()
            A 0-dimensional polyhedron in QQ^3 defined as the convex hull of 1 vertex.
        """
        return self._polyhedron

    def __call__(self):
        """
        Returns the vector representation of the representation
        object. Shorthand for the vector() method.

        TESTS::

            sage: import sage.geometry.polyhedra as P
            sage: pr = P.PolyhedronRepresentation(Polyhedron(vertices = [[1,2,3]]), [1,2,3])
            sage: pr.__call__()
            (1, 2, 3)
        """
        return self._representation_vector

    def vector(self):
        """
        Returns the vector representation of the representation object.

        NOTES:

          * For vertices/lines/rays this is a vector of length
            ``ambient_dim()``. For inequalities and equations this is
            a vector of length ``ambient_dim()+1``

        EXAMPLES::

            sage: s = polytopes.cuboctahedron()
            sage: v = s.vertex_generator().next()
            sage: v
            A vertex at (0, -1/2, -1/2)
            sage: v.vector()
            (0, -1/2, -1/2)
            sage: v()
            (0, -1/2, -1/2)
            sage: type(v())
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        """
        return self._representation_vector

    def index(self):
        """
        Returns an arbitrary but fixed number according to the internal
        storage order.

        NOTES:

        H-representation and V-representation objects are enumerated
        independently. That is, amongst all vertices/rays/lines there
        will be one with ``index()==0``, and amongs all
        inequalities/equations there will be one with ``index()==0``,
        unless the polyhedron is empty or spans the whole space.

        EXAMPLES::

            sage: s = Polyhedron(vertices=[[1],[-1]])
            sage: first_vertex = s.vertex_generator().next()
            sage: first_vertex.index()
            0
            sage: first_vertex == s.Vrepresentation(0)
            True
        """
        return self._index



class Hrepresentation(PolyhedronRepresentation):
    """
    The internal base class for H-representation objects of
    a polyhedron. Inherits from ``PolyhedronRepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initialization function.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
        """
        if len(data) != 1+polyhedron.ambient_dim():
            raise ValueError, \
                'H-representation data requires a list of length ambient_dim+1'
        super(Hrepresentation, self).__init__(polyhedron, data)
        self._index = len(polyhedron._Hrepresentation)
        polyhedron._Hrepresentation.append(self)

    def is_H(self):
        """
        Returns True if the object is part of a H-representation
        (inequality or equation).

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_H()
            True
        """
        return True

    def is_inequality(self):
        """
        Returns True if the object is an inequality of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_inequality()
            True
        """
        return False

    def is_equation(self):
        """
        Returns True if the object is an equation of the H-representation.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]], eqns = [[1,1,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_equation()
            True
        """
        return False

    def A(self):
        r"""
        Returns the coefficient vector `A` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.A()
            (1, 0)
        """
        return self._representation_vector[1:]

    def b(self):
        r"""
        Returns the constant `b` in `A\vec{x}+b`.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,0,],[1,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.b()
            0
        """
        return self._representation_vector[0]

    def neighbors(self):
        """
        Iterate over the adjacent facets (i.e. inequalities/equations)

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],[1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: a[0]
            An inequality (0, 1, 0) x + 0 >= 0
            sage: list(a[0])
            [0, 0, 1, 0]
        """
        adjacency_matrix = self.polyhedron().facet_adjacency_matrix()
        for x in self.polyhedron().Hrep_generator():
            if adjacency_matrix[self.index(), x.index()] == 1:
                yield x

    def adjacent(self):
        """
        Alias for neighbors().

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,0,2],[0,0,1,0,],[0,10,0,0],[1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: a = list(pH.neighbors())
            sage: b = list(pH.adjacent())
            sage: a==b
            True
        """
        return self.neighbors()

    def is_incident(self, Vobj):
        """
        Returns whether the incidence matrix element (Vobj,self) == 1

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],[1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH.is_incident(p.Vrepresentation(1))
            True
            sage: pH.is_incident(p.Vrepresentation(5))
            False
        """
        return self.polyhedron().incidence_matrix()[Vobj.index(), self.index()] == 1

    def __mul__(self, Vobj):
        """
        Shorthand for ``self.eval(x)``

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0,],[0,1,0,0],[1,-1,0,0],[1,0,-1,0,],[1,0,0,-1]])
            sage: pH = p.Hrepresentation(0)
            sage: pH*p.Vrepresentation(5)
            1
        """
        return self.eval(Vobj)

    def eval(self, Vobj):
        r"""
        Evaluates the left hand side `A\vec{x}+b` on the given
        vertex/ray/line.

        NOTES:

          * Evaluating on a vertex returns `A\vec{x}+b`
          * Evaluating on a ray returns `A\vec{r}`. Only the sign or
            whether it is zero is meaningful.
          * Evaluating on a line returns `A\vec{l}`. Only whether it
            is zero or not is meaningful.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = triangle.inequality_generator().next()
            sage: ineq
            An inequality (-1, -1) x + 1 >= 0
            sage: [ ineq.eval(v) for v in triangle.vertex_generator() ]
            [0, 0, 3]
            sage: [ ineq * v for v in triangle.vertex_generator() ]
            [0, 0, 3]

        If you pass a vector, it is assumed to be the coordinate vector of a point::

            sage: ineq.eval( vector(ZZ, [3,2]) )
            -4
        """
        try:
            return Vobj.evaluated_on(self)
        except AttributeError:
            return self.A() * Vobj + self.b()

    def incident(self):
        """
        Returns a generator for the incident H-representation objects,
        that is, the vertices/rays/lines satisfying the (in)equality.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = triangle.inequality_generator().next()
            sage: ineq
            An inequality (-1, -1) x + 1 >= 0
            sage: [ v for v in ineq.incident()]
            [A vertex at (1, 0), A vertex at (0, 1)]
            sage: p = Polyhedron(vertices=[[0,0,0],[0,1,0],[0,0,1]], rays=[[1,-1,-1]])
            sage: ineq = p.Hrepresentation(2)
            sage: ineq
            An inequality (1, 1, 0) x + 0 >= 0
            sage: [ x for x in ineq.incident() ]
            [A ray in the direction (1, -1, -1),
             A vertex at (0, 0, 0),
             A vertex at (0, 0, 1)]
        """
        incidence_matrix = self.polyhedron().incidence_matrix()
        for V in self.polyhedron().Vrep_generator():
            if incidence_matrix[V.index(), self.index()] == 1:
                yield V


class Inequality(Hrepresentation):
    """
    A linear inequality (supporting hyperplane) of the
    polyhedron. Inherits from ``Hrepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        A linear inequality (supporting hyperplane) of the
        polyhedron. Inherits from ``Hrepresentation``.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.inequality_generator().next()
            sage: a
            An inequality (-1, 0, 0) x + 1 >= 0
        """
        super(Inequality, self).__init__(polyhedron, data)


    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = p.inequality_generator().next()
            sage: repr_obj.type()
            0
            sage: repr_obj.type() == repr_obj.INEQUALITY
            True
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.INEQUALITY


    def is_inequality(self):
        """
        Returns True since this is, by construction, an inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.inequality_generator().next()
            sage: a.is_inequality()
            True
        """
        return True

    def _repr_(self):
        """
        The string representation of the inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.inequality_generator().next()
            sage: a._repr_()
            'An inequality (-1, 0, 0) x + 1 >= 0'
            sage: Polyhedron(ieqs=[(1,-1),(-1,2)]).Hrepresentation()
            [An inequality (-1) x + 1 >= 0, An inequality (2) x - 1 >= 0]
            sage: Polyhedron(eqns=[(1,0)]).Hrepresentation()
            [An equation 1 == 0]
            sage: Polyhedron(eqns=[(-1,0)]).Hrepresentation()
            [An equation -1 == 0]
        """
        s = 'An inequality '
        have_A = not self.A().is_zero()
        if have_A:
            s += repr(self.A()) + ' x '
        if self.b()>=0:
            if have_A:
                s += '+'
        else:
            s += '-'
        if have_A:
            s += ' '
        s += repr(abs(self.b())) + ' >= 0'
        return s

    def contains(self, Vobj):
        """
        Tests whether the halfspace (including its boundary) defined
        by the inequality contains the given vertex/ray/line.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: i1 = p.inequality_generator().next()
            sage: [i1.contains(q) for q in p.vertex_generator()]
            [True, True, True, True, True, True]
            sage: p2 = 3*polytopes.n_cube(3)
            sage: [i1.contains(q) for q in p2.vertex_generator()]
            [True, False, True, True, False, False, True, False]
        """
        try:
            if Vobj.is_vector(): # assume we were passed a point
                return self.polyhedron()._is_nonneg( self.eval(Vobj) )
        except AttributeError:
            pass

        if Vobj.is_line():
            return self.polyhedron()._is_zero( self.eval(Vobj) )
        else:
            return self.polyhedron()._is_nonneg( self.eval(Vobj) )

    def interior_contains(self, Vobj):
        """
        Tests whether the interior of the halfspace (excluding its
        boundary) defined by the inequality contains the given
        vertex/ray/line.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: i1 = p.inequality_generator().next()
            sage: [i1.interior_contains(q) for q in p.vertex_generator()]
            [True, False, True, False, True, False]
            sage: p2 = 3*polytopes.n_cube(3)
            sage: [i1.interior_contains(q) for q in p2.vertex_generator()]
            [True, False, True, True, False, False, True, False]

        If you pass a vector, it is assumed to be the coordinate vector of a point::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[-1,1],[-1,-1]])
            sage: p = vector(ZZ, [1,0] )
            sage: [ ieq.interior_contains(p) for ieq in P.inequality_generator() ]
            [True, True, True, False]
        """
        try:
            if Vobj.is_vector(): # assume we were passed a point
                return self.polyhedron()._is_positive( self.eval(Vobj) )
        except AttributeError:
            pass

        if Vobj.is_line():
            return self.polyhedron()._is_zero( self.eval(Vobj) )
        elif Vobj.is_vertex():
            return self.polyhedron()._is_positive( self.eval(Vobj) )
        else: # Vobj.is_ray()
            return self.polyhedron()._is_nonneg( self.eval(Vobj) )


class Equation(Hrepresentation):
    """
    A linear equation of the polyhedron. That is, the polyhedron is
    strictly smaller-dimensional than the ambient space, and contained
    in this hyperplane. Inherits from ``Hrepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initializes an equation object.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.equation_generator().next()
            sage: a
            An equation (0, 0, 1) x + 0 == 0
        """
        super(Equation, self).__init__(polyhedron, data)


    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = p.equation_generator().next()
            sage: repr_obj.type()
            1
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            True
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.EQUATION


    def is_equation(self):
        """
        Tests if this object is an equation.  By construction, it must be.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.equation_generator().next()
            sage: a.is_equation()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of this object.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.equation_generator().next()
            sage: a._repr_()
            'An equation (0, 0, 1) x + 0 == 0'
            sage: Polyhedron().Hrepresentation(0)
            An equation -1 == 0
        """
        s = 'An equation '
        have_A = not self.A().is_zero()
        if have_A:
            s += repr(self.A()) + ' x '
        if self.b()>=0:
            if have_A:
                s += '+'
        else:
            s += '-'
        if have_A:
            s += ' '
        s += repr(abs(self.b())) + ' == 0'
        return s

    def contains(self, Vobj):
        """
        Tests whether the hyperplane defined by the equation contains
        the given vertex/ray/line.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: v = p.vertex_generator().next()
            sage: v
            A vertex at (0, 0, 0)
            sage: a = p.equation_generator().next()
            sage: a
            An equation (0, 0, 1) x + 0 == 0
            sage: a.contains(v)
            True
        """
        return self.polyhedron()._is_zero( self.eval(Vobj) )

    def interior_contains(self, Vobj):
        """
        Tests whether the interior of the halfspace (excluding its
        boundary) defined by the inequality contains the given
        vertex/ray/line.

        NOTE:

        Returns False for any equation.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: v = p.vertex_generator().next()
            sage: v
            A vertex at (0, 0, 0)
            sage: a = p.equation_generator().next()
            sage: a
            An equation (0, 0, 1) x + 0 == 0
            sage: a.interior_contains(v)
            False
        """
        return False


class Vrepresentation(PolyhedronRepresentation):
    """
    The base class for V-representation objects of a
    polyhedron. Inherits from ``PolyhedronRepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initialization function for the vertex representation.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: from sage.geometry.polyhedra import Vrepresentation
            sage: p._Vrepresentation = Sequence([])
            sage: vr = Vrepresentation(p,[1,2,3,4])
            sage: vr.__init__(p,[1,2,3,5])
            sage: vr.vector()
            (1, 2, 3, 5)
        """
        if len(data) != polyhedron.ambient_dim():
            raise ValueError, \
                'V-representation data requires a list of length ambient_dim'
        super(Vrepresentation, self).__init__(polyhedron, data)
        self._index = len(polyhedron._Vrepresentation)
        polyhedron._Vrepresentation.append(self)

    def is_V(self):
        """
        Returns True if the object is part of a V-representation
        (a vertex, ray, or line).

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,3]])
            sage: v = p.vertex_generator().next()
            sage: v.is_V()
            True
        """
        return True

    def is_vertex(self):
        """
        Returns True if the object is a vertex of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Vertex.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,3]])
            sage: v = p.vertex_generator().next()
            sage: v.is_vertex()
            True
            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: r1 = p.ray_generator().next()
            sage: r1.is_vertex()
            False
        """
        return False

    def is_ray(self):
        """
        Returns True if the object is a ray of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Ray.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: r1 = p.ray_generator().next()
            sage: r1.is_ray()
            True
            sage: v1 = p.vertex_generator().next()
            sage: v1
            A vertex at (-1, -1, 0, -1)
            sage: v1.is_ray()
            False
        """
        return False

    def is_line(self):
        """
        Returns True if the object is a line of the V-representation.
        This method is over-ridden by the corresponding method in the
        derived class Line.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 0, 1], [1, 1, 0, 0, 0], [1, 0, 1, 0, 0]])
            sage: line1 = p.line_generator().next()
            sage: line1.is_line()
            True
            sage: v1 = p.vertex_generator().next()
            sage: v1.is_line()
            False
        """
        return False

    def neighbors(self):
        """
        Returns a generator for the adjacent vertices/rays/lines.

        EXAMPLES::

             sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,4]])
             sage: v = p.vertex_generator().next()
             sage: v.neighbors().next()
             A vertex at (1, 0)
        """
        adjacency_matrix = self.polyhedron().vertex_adjacency_matrix()
        for x in self.polyhedron().Vrep_generator():
            if adjacency_matrix[self.index(), x.index()] == 1:
                yield x

    def adjacent(self):
        """
        Alias for neighbors().

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,3],[1,4]])
            sage: v = p.vertex_generator().next()
            sage: a = v.neighbors().next()
            sage: b = v.adjacent().next()
            sage: a==b
            True
        """
        return self.neighbors()

    def is_incident(self, Hobj):
        """
        Returns whether the incidence matrix element (self,Hobj) == 1

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: h1 = p.inequality_generator().next()
            sage: h1
            An inequality (0, 0, 1) x + 1 >= 0
            sage: v1 = p.vertex_generator().next()
            sage: v1
            A vertex at (1, 1, 1)
            sage: v1.is_incident(h1)
            False
        """
        return self.polyhedron().incidence_matrix()[self.index(), Hobj.index()] == 1

    def __mul__(self, Hobj):
        """
        Shorthand for self.evaluated_on(Hobj)

        TESTS::

            sage: p = polytopes.n_cube(3)
            sage: h1 = p.inequality_generator().next()
            sage: v1 = p.vertex_generator().next()
            sage: v1.__mul__(h1)
            2
        """
        return self.evaluated_on(Hobj)

    def incident(self):
        """
        Returns a generator for the equations/inequalities that are satisfied on the given
        vertex/ray/line.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[-1,-1]])
            sage: ineq = triangle.inequality_generator().next()
            sage: ineq
            An inequality (-1, -1) x + 1 >= 0
            sage: [ v for v in ineq.incident()]
            [A vertex at (1, 0), A vertex at (0, 1)]
            sage: p = Polyhedron(vertices=[[0,0,0],[0,1,0],[0,0,1]], rays=[[1,-1,-1]])
            sage: ineq = p.Hrepresentation(2)
            sage: ineq
            An inequality (1, 1, 0) x + 0 >= 0
            sage: [ x for x in ineq.incident() ]
            [A ray in the direction (1, -1, -1),
             A vertex at (0, 0, 0),
             A vertex at (0, 0, 1)]
        """
        incidence_matrix = self.polyhedron().incidence_matrix()
        for H in self.polyhedron().Hrep_generator():
            if incidence_matrix[self.index(), H.index()] == 1:
                yield H


class Vertex(Vrepresentation):
    """
    A vertex of the polyhedron. Inherits from ``Vrepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initializes the vertex object.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: v1 = p.vertex_generator().next()
            sage: v1
            A vertex at (1, 0)
            sage: p._Vrepresentation = Sequence([])
            sage: v1.__init__(p,[1,2])
            sage: v1
            A vertex at (1, 2)
        """
        super(Vertex, self).__init__(polyhedron, data)


    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: repr_obj = p.vertex_generator().next()
            sage: repr_obj.type()
            2
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            True
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.VERTEX


    def is_vertex(self):
        """
        Tests if this object is a vertex.  By construction it always is.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = p.vertex_generator().next()
            sage: a.is_vertex()
            True
        """
        return True

    def _repr_(self):
        """
        Returns a string representation of the vertex.

        OUTPUT:

        String.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: v = p.vertex_generator().next()
            sage: v.__repr__()
            'A vertex at (1, 0)'
        """
        return 'A vertex at ' + repr(self._representation_vector);

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{x}+b`

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: v = p.vertex_generator().next()
            sage: h = p.inequality_generator().next()
            sage: v
            A vertex at (1, 1, 1)
            sage: h
            An inequality (0, 0, 1) x + 1 >= 0
            sage: v.evaluated_on(h)
            2
        """
        return Hobj.A() * self._representation_vector + Hobj.b()

    def is_integral(self):
        r"""
        Return whether the coordinates of the vertex are all integral.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = Polyhedron([(1/2,3,5), (0,0,0), (2,3,7)])
            sage: [ v.is_integral() for v in p.vertex_generator() ]
            [False, True, True]
        """
        return all(x in ZZ for x in self._representation_vector)


class Ray(Vrepresentation):
    """
    A ray of the polyhedron. Inherits from ``Vrepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initializes a ray object.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: r1 = p.ray_generator().next()
            sage: r1
            A ray in the direction (0, 1)
            sage: p._Vrepresentation = Sequence([])
            sage: r1.__init__(p,[1,2])
            sage: p.ray_generator().next()
            A ray in the direction (1, 2)
        """
        super(Ray, self).__init__(polyhedron, data)


    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: repr_obj = p.ray_generator().next()
            sage: repr_obj.type()
            3
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            True
            sage: repr_obj.type() == repr_obj.LINE
            False
        """
        return self.RAY


    def is_ray(self):
        """
        Tests if this object is a ray.  Always True by construction.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = p.ray_generator().next()
            sage: a.is_ray()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of the ray.

        TESTS::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = p.ray_generator().next()
            sage: a._repr_()
            'A ray in the direction (0, 1)'
        """
        return 'A ray in the direction ' + repr(self._representation_vector);

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{r}`

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,1],[0,1,0],[1,-1,0]])
            sage: a = p.ray_generator().next()
            sage: h = p.inequality_generator().next()
            sage: a.evaluated_on(h)
            1
        """
        return Hobj.A() * self._representation_vector


class Line(Vrepresentation):
    r"""
    A line (Minkowski summand `\simeq\RR`) of the
    polyhedron. Inherits from ``Vrepresentation``.
    """
    def __init__(self, polyhedron, data):
        """
        Initializes the line object.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = p.line_generator().next()
            sage: p._Vrepresentation = Sequence([])
            sage: a.__init__(p,[1,2,3])
            sage: a
            A line in the direction (1, 2, 3)
        """
        super(Line, self).__init__(polyhedron, data)


    def type(self):
        r"""
        Returns the type (equation/inequality/vertex/ray/line) as an
        integer.

        OUTPUT:

        Integer. One of ``PolyhedronRepresentation.INEQUALITY``,
        ``.EQUATION``, ``.VERTEX``, ``.RAY``, or ``.LINE``.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: repr_obj = p.line_generator().next()
            sage: repr_obj.type()
            4
            sage: repr_obj.type() == repr_obj.INEQUALITY
            False
            sage: repr_obj.type() == repr_obj.EQUATION
            False
            sage: repr_obj.type() == repr_obj.VERTEX
            False
            sage: repr_obj.type() == repr_obj.RAY
            False
            sage: repr_obj.type() == repr_obj.LINE
            True
        """
        return self.LINE


    def is_line(self):
        """
        Tests if the object is a line.  By construction it must be.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = p.line_generator().next()
            sage: a.is_line()
            True
        """
        return True

    def _repr_(self):
        """
        A string representation of the line.

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = p.line_generator().next()
            sage: a.__repr__()
            'A line in the direction (0, 1, 0)'
        """
        return 'A line in the direction ' + repr(self._representation_vector);

    def evaluated_on(self, Hobj):
        r"""
        Returns `A\vec{\ell}`

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: a = p.line_generator().next()
            sage: h = p.inequality_generator().next()
            sage: a.evaluated_on(h)
            0
        """
        return Hobj.A() * self._representation_vector


#############################################################
def render_2d(projection, **kwds):
    """
    Return 2d rendering of the projection of a polyhedron into
    2-dimensional ambient space.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[[1,1]], rays=[[1,1]])
        sage: q1 = p1.projection()
        sage: p2 = Polyhedron(vertices=[[1,0], [0,1], [0,0]])
        sage: q2 = p2.projection()
        sage: p3 = Polyhedron(vertices=[[1,2]])
        sage: q3 = p3.projection()
        sage: p4 = Polyhedron(vertices=[[2,0]], rays=[[1,-1]], lines=[[1,1]])
        sage: q4 = p4.projection()
        sage: q1.show() + q2.show() + q3.show() + q4.show()
        sage: from sage.geometry.polyhedra import render_2d
        sage: q = render_2d(p1.projection())
        sage: q._Graphics__objects
        [Point set defined by 1 point(s), Arrow from (1.0,1.0) to (2.0,2.0), Polygon defined by 3 points]
    """
    if isinstance(projection, Polyhedron): projection = Projection(projection)
    return \
        projection.render_points_2d(zorder=2, pointsize=10, **kwds) + \
        projection.render_outline_2d(zorder=1, **kwds) + \
        projection.render_fill_2d(zorder=0, rgbcolor=(0,1,0), **kwds)


def render_3d(projection, **kwds):
    """
    Return 3d rendering of a polyhedron projected into
    3-dimensional ambient space.

    NOTE:

    This method, ``render_3d``, is used in the ``show()``
    method of a polyhedron if it is in 3 dimensions.

    EXAMPLES::

        sage: p1 = Polyhedron(vertices=[[1,1,1]], rays=[[1,1,1]])
        sage: p2 = Polyhedron(vertices=[[2,0,0], [0,2,0], [0,0,2]])
        sage: p3 = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]], rays=[[-1,-1,-1]])
        sage: p1.projection().show() + p2.projection().show() + p3.projection().show()

    It correctly handles various degenerate cases::

        sage: Polyhedron(lines=[[1,0,0],[0,1,0],[0,0,1]]).show()                              # whole space
        sage: Polyhedron(vertices=[[1,1,1]], rays=[[1,0,0]], lines=[[0,1,0],[0,0,1]]).show()  # half space
        sage: Polyhedron(vertices=[[1,1,1]], lines=[[0,1,0],[0,0,1]]).show()                  # R^2 in R^3
        sage: Polyhedron(rays=[[0,1,0],[0,0,1]], lines=[[1,0,0]]).show()                      # quadrant wedge in R^2
        sage: Polyhedron(rays=[[0,1,0]], lines=[[1,0,0]]).show()                              # upper half plane in R^3
        sage: Polyhedron(lines=[[1,0,0]]).show()                                              # R^1 in R^2
        sage: Polyhedron(rays=[[0,1,0]]).show()                                               # Half-line in R^3
        sage: Polyhedron(vertices=[[1,1,1]]).show()                                           # point in R^3
    """
    if isinstance(projection, Polyhedron): projection = Projection(projection)
    return \
        projection.render_vertices_3d(width=3, color='green', **kwds) +\
        projection.render_wireframe_3d(width=3, color='green', **kwds) + \
        projection.render_solid_3d(**kwds)


def render_4d(polyhedron, **kwds):
    """
    Return a 3d rendering of the Schlegel projection of a 4d
    polyhedron projected into 3-dimensional space.

    NOTES:

    The ``show()`` method of ``Polyhedron()`` uses this to draw itself
    if the ambient dimension is 4.

    INPUT:

    - ``polyhedron`` - A ``Polyhedron()`` object.
    - ``kwds`` - plot keywords. Passing
      ``projection_direction=<list>`` sets the projetion direction of
      the Schlegel projection. If it is not given, the center of a
      facet is used.

    EXAMPLES::

        sage: poly = polytopes.twenty_four_cell()
        sage: poly
        A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 24 vertices.
        sage: poly.show()
        sage: poly.show(projection_direction=[2,5,11,17])
        sage: type( poly.show() )
        <class 'sage.plot.plot3d.base.Graphics3dGroup'>

    TESTS::

        sage: from sage.geometry.polyhedra import render_4d
        sage: p = polytopes.n_cube(4)
        sage: q = render_4d(p)
        sage: tach_str = q.tachyon()
        sage: tach_str.count('FCylinder')
        32
    """
    projection_direction = None
    try:
        projection_direction = kwds.pop('projection_directior')
    except KeyError:
        for ineq in polyhedron.inequality_generator():
            center = [v() for v in ineq.incident() if v.is_vertex()]
            center = sum(center) / len(center)
            if not center.is_zero():
                projection_direction = center
                break
    projection_3d = Projection(polyhedron).schlegel(projection_direction)
    return render_3d(projection_3d, **kwds)


#########################################################################
class Polyhedron(SageObject):
    """
    A class for polyhedral objects. You may either instantiate with
    vertex/ray/line or inequalities/equations data, but not
    both. Redundant data will automatically be removed, and the
    complementary representation will be computed.

    EXAMPLES:

    Construct some polyhedra::

        sage: square_from_vertices = Polyhedron(vertices = [[1, 1], [1, -1], [-1, 1], [-1, -1]])
        sage: square_from_ieqs = Polyhedron(ieqs = [[1, 0, 1], [1, 1, 0], [1, 0, -1], [1, -1, 0]])
        sage: list(square_from_ieqs.vertex_generator())
        [A vertex at (1, -1),
         A vertex at (1, 1),
         A vertex at (-1, 1),
         A vertex at (-1, -1)]
        sage: list(square_from_vertices.inequality_generator())
        [An inequality (0, 1) x + 1 >= 0,
         An inequality (1, 0) x + 1 >= 0,
         An inequality (0, -1) x + 1 >= 0,
         An inequality (-1, 0) x + 1 >= 0]
        sage: p = Polyhedron(vertices = [[1.1, 2.2], [3.3, 4.4]], field=RDF)
        sage: p.n_inequalities()
        2

    The same polyhedron given in two ways::

        sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0]])
        sage: p.Vrepresentation()
        [A line in the direction (0, 0, 1),
         A ray in the direction (1, 0, 0),
         A ray in the direction (0, 1, 0),
         A vertex at (0, 0, 0)]
        sage: q = Polyhedron(vertices=[[0,0,0]], rays=[[1,0,0],[0,1,0]], lines=[[0,0,1]])
        sage: q.Hrepresentation()
        [An inequality (1, 0, 0) x + 0 >= 0,
         An inequality (0, 1, 0) x + 0 >= 0]

    Finally, a more complicated example. Take `\mathbb{R}_{\geq 0}^6` with
    coordinates `a, b, \dots, f` and

      * The inequality `e+b \geq c+d`
      * The inequality `e+c \geq b+d`
      * The equation `a+b+c+d+e+f = 31`

    ::

        sage: positive_coords = Polyhedron(ieqs=[[0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
        sage: P = Polyhedron(ieqs=positive_coords.inequalities() + [[0,0,1,-1,-1,1,0], [0,0,-1,1,-1,1,0]], eqns=[[-31,1,1,1,1,1,1]])
        sage: P
        A 5-dimensional polyhedron in QQ^6 defined as the convex hull of 7 vertices.
        sage: P.dim()
        5
        sage: P.Vrepresentation()
        [A vertex at (0, 31/2, 31/2, 0, 0, 0), A vertex at (0, 31/2, 0, 0, 31/2, 0), A vertex at (0, 0, 0, 0, 31, 0), A vertex at (0, 0, 31/2, 0, 31/2, 0), A vertex at (0, 0, 0, 31/2, 31/2, 0), A vertex at (31, 0, 0, 0, 0, 0), A vertex at (0, 0, 0, 0, 0, 31)]

    Finally, a more complicated example. Take `\mathbb{R}_{\geq 0}^6` with
    coordinates `a, b, \dots, f` and

      * The inequality `e+b \geq c+d`
      * The inequality `e+c \geq b+d`
      * The equation `a+b+c+d+e+f = 31`

    ::

        sage: positive_coords = Polyhedron(ieqs=[[0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
        sage: P = Polyhedron(ieqs=positive_coords.inequalities() + [[0,0,1,-1,-1,1,0], [0,0,-1,1,-1,1,0]], eqns=[[-31,1,1,1,1,1,1]])
        sage: P
        A 5-dimensional polyhedron in QQ^6 defined as the convex hull of 7 vertices.
        sage: P.dim()
        5
        sage: P.Vrepresentation()
        [A vertex at (0, 31/2, 31/2, 0, 0, 0), A vertex at (0, 31/2, 0, 0, 31/2, 0), A vertex at (0, 0, 0, 0, 31, 0), A vertex at (0, 0, 31/2, 0, 31/2, 0), A vertex at (0, 0, 0, 31/2, 31/2, 0), A vertex at (31, 0, 0, 0, 0, 0), A vertex at (0, 0, 0, 0, 0, 31)]

    NOTES:

      * Once constructed, a ``Polyhedron`` object is immutable.
      * Although the option ``field = RDF`` allows numerical data to
        be used, it might not give the right answer for degenerate
        input data - the results can depend upon the tolerance
        setting of cddlib.
    """
    _render_method = [ None, None, render_2d, render_3d, render_4d ]


    def __init__(self,
                 vertices = None, rays = None, lines = None,
                 ieqs = None, eqns = None,
                 field = QQ, compute_adjacency=False, verbose = False):
        """
        Initializes the polyhedron.

        INPUT:

        - ``vertices`` -- list of points. Each point can be specified
          as any iterable container of ``field`` elements.

        - ``rays`` -- list of rays. Each ray can be specified as any
          iterable container of ``field`` elements.

        - ``lines`` -- list of lines. Each line can be specified as
          any iterable container of ``field`` elements.

        - ``ieqs`` -- list of inequalities. Each line can be specified as
          any iterable container of ``field`` elements.

        - ``eqns`` -- list of equalities. Each line can be specified
          as any iterable container of ``field`` elements.

        - ``field`` -- either ``QQ`` or ``RDF``. The field over which
          the polyhedron will be defined. For ``QQ``, exact arithmetic
          will be used. For ``RDF``, floating point numbers will be
          used. Floating point arithmetic is faster but might give the
          wrong result for degenerate input.

        - ``compute_adjacency`` -- boolean. Whether to compute the
          adjacency data of vertices and facets during
          initialization. This might take considerably longer than
          just computing a minimal H/V-representation pair. If
          adjacency data is not computed during initialization it will
          be computed as needed, but this then requires another call
          do cdd.

        .. NOTE::

            You can either specify a subset of (``vertices``,
            ``rays``, ``lines``) or a subset of (``ieqs``,
            ``eqns``). But you must not combine V-representation and
            H-representation objects.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[1,1]])
            sage: p
            A 1-dimensional polyhedron in QQ^2 defined as the convex hull of 2 vertices.
            sage: del p._n_vertices
            sage: del p._equations
            sage: del p._n_equations
            sage: p.__init__(vertices = [[1,1]])
            sage: p
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex.

        Degenerate input::

            sage: p = Polyhedron(ieqs=[(-1,0),(1,0)]); p  # incompatible inequalities
            The empty polyhedron in QQ^1.
            sage: p.Hrepresentation()
            [An equation -1 == 0]
            sage: p.Vrepresentation()
            []
            sage: p = Polyhedron(eqns=[(0,1,0),(0,0,1)]); p  # equations cutting out a point
            A 0-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex.
            sage: p.Hrepresentation()
            [An equation (1, 0) x + 0 == 0, An equation (0, 1) x + 0 == 0]
            sage: p.Vrepresentation()
            [A vertex at (0, 0)]
            sage: p = Polyhedron(ieqs=[(0,1,0)]); p  # half plane
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 vertex, 1 ray, 1 line.
            sage: p.Hrepresentation()
            [An inequality (1, 0) x + 0 >= 0]
            sage: p.Vrepresentation()
            [A line in the direction (0, 1), A ray in the direction (1, 0), A vertex at (0, 0)]
        """
        self._init_field(field)

        # Clean up the arguments
        vertices = _set_to_None_if_empty(vertices)
        rays     = _set_to_None_if_empty(rays)
        lines    = _set_to_None_if_empty(lines)
        ieqs     = _set_to_None_if_empty(ieqs)
        eqns     = _set_to_None_if_empty(eqns)

        got_Vrep = (vertices!=None or rays!=None or lines!=None)
        got_Hrep = (ieqs!=None or eqns!=None)

        if got_Vrep and got_Hrep:
            raise ValueError, "Must not be initialized with H-rep and V-rep."
        elif got_Vrep:
            s = cdd_Vrepresentation(self._cdd_type, vertices, rays, lines)
        elif got_Hrep:
            s = cdd_Hrepresentation(self._cdd_type, ieqs, eqns)
        else:
            self._init_empty_polyhedron()
            return

        # Compute complementary representation via cdd
        if compute_adjacency:
            cmdline_arg = '--all'
        else:
            cmdline_arg = '--reps'
        self._init_from_cdd_input(s, cmdline_arg, verbose)


    def __getstate__(self):
        r"""
        Return the dictionary that should be pickled.

        OUTPUT:

        A :class:`dict`.

        EXAMPLES::

            sage: P = polytopes.dodecahedron(QQ)
            sage: P == loads(dumps(P))
            True
        """
        exclude_keys = ['_is_zero', '_is_positive', '_is_nonneg']
        state = {}
        for key,value in self.__dict__.items():
            if key in exclude_keys:
                continue
            state[key] = value
        return state


    def __setstate__(self, state):
        r"""
        Return the dictionary that should be pickled.

        OUTPUT:

        A :class:`dict`.

        EXAMPLES::

            sage: P = polytopes.dodecahedron(QQ)
            sage: Q = loads(dumps(P))
            sage: P.Hrepresentation(0).polyhedron() is P
            True
            sage: Q.Hrepresentation(0).polyhedron() is Q
            True
            sage: P == Q
            True
            sage: P is Q
            False
        """
        exclude_keys = ['_is_zero', '_is_positive', '_is_nonneg']
        for key,value in state.items():
            if key in exclude_keys:
                continue
            if key == '_field':
                self._init_field(value)
                continue
            self.__dict__[key] = value
        return state


    def __lt__(self, other):
        """
        Test whether ``self`` is a strict sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P < Q   # indirect doctest
            False
            sage: P < P   # indirect doctest
            False
            sage: Q < P   # indirect doctest
            True
        """
        return self._is_subpolyhedron(other) and not other._is_subpolyhedron(self)


    def __le__(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P <= Q   # indirect doctest
            False
            sage: P <= P   # indirect doctest
            True
            sage: Q <= P   # indirect doctest
            True
        """
        return self._is_subpolyhedron(other)


    def __eq__(self, other):
        """
        Test whether ``self`` is a strict sub-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P == Q   # indirect doctest
            False
            sage: P == P   # indirect doctest
            True
            sage: Q == P   # indirect doctest
            False
        """
        return self._is_subpolyhedron(other) and other._is_subpolyhedron(self)


    def __ne__(self, other):
        """
        Test whether ``self`` is not equal to ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P != Q   # indirect doctest
            True
            sage: P != P   # indirect doctest
            False
            sage: Q != P   # indirect doctest
            True
        """
        return not self.__eq__(other)


    def __gt__(self, other):
        """
        Test whether ``self`` is a strict super-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P > Q   # indirect doctest
            True
            sage: P > P   # indirect doctest
            False
            sage: Q > P   # indirect doctest
            False
        """
        return other._is_subpolyhedron(self) and not self._is_subpolyhedron(other)


    def __ge__(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        super-polyhedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P >= Q   # indirect doctest
            True
            sage: P >= P   # indirect doctest
            True
            sage: Q >= P   # indirect doctest
            False
        """
        return other._is_subpolyhedron(self)


    def _is_subpolyhedron(self, other):
        """
        Test whether ``self`` is a (not necessarily strict)
        sub-polyhdedron of ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (0,1)], rays=[(1,1)])
            sage: Q = Polyhedron(vertices=[(1,0), (0,1)])
            sage: P._is_subpolyhedron(Q)
            False
            sage: Q._is_subpolyhedron(P)
            True
        """
        if not isinstance(other,Polyhedron):
            raise ValueError('Can only compare Polyhedron objects.')
        return all( other_H.contains(self_V)
                    for other_H, self_V in
                    CartesianProduct(other.Hrep_generator(), self.Vrep_generator()) )


    def plot(self, **kwds):
        """
        Return a graphical representation.

        INPUT:

        - ``**kwds`` -- optional keyword parameters.

        See :func:`render_2d`, :func:`render_3d`, :func:`render_4d`
        for a description of available options for different ambient
        space dimensions.

        OUTPUT:

        A graphics object.

        TESTS::

            sage: polytopes.n_cube(2).plot()
        """
        if self.ambient_dim() < len(Polyhedron._render_method):
            render = Polyhedron._render_method[self.ambient_dim()]
            if render != None:
                return render(self,**kwds)
        raise NotImplementedError('Plotting of '+str(self.ambient_dim())+
                                  '-dimensional polyhedra not implemented')


    show = plot


    def _init_field(self, field):
        """
        Internal method: initialize the cdd numeric type and
        executable depending on the desired field. Supported are
        rational numbers (using the gmp library) or floating-point
        numbers (which can yield wrong results for degenerate
        inequalities/equations due to rounding errors)

        TESTS::

            sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]])
            sage: p._field==QQ
            True
            sage: p._init_field(RDF)
            sage: p._field
            Real Double Field
        """
        if field == QQ:
            self._cdd_type = 'rational'
            self._cdd_executable = 'cdd_both_reps_gmp'
            self._is_zero = lambda x: x==0
            self._is_nonneg = lambda x: x>=0
            self._is_positive = lambda x: x>0
        elif field == RDF:
            self._cdd_type = 'real'
            self._cdd_executable = 'cdd_both_reps'
            # 1e-6 is the cddf+ default fuzzy zero cutoff
            self._is_zero = lambda x: abs(x)<=1e-6
            self._is_nonneg = lambda x: x>=-1e-6
            self._is_positive = lambda x: x>=-1e-6
        else:
            raise ValueError, 'Field must be QQ or RDF.'
        self._field = field


    def _init_empty_polyhedron(self):
        """
        Initializes an empty polyhedron in a zero-dimensional ambient space.

        TESTS::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0.
            sage: empty.Vrepresentation()
            []
            sage: empty.Hrepresentation()
            [An equation -1 == 0]
            sage: Polyhedron(vertices = [])
            The empty polyhedron in QQ^0.
            sage: Polyhedron()._init_empty_polyhedron()
        """
        self._ambient_dim = 0

        self._Vrepresentation = Sequence([])
        self._Vrepresentation.set_immutable()

        self._Hrepresentation = Sequence([])
        Equation(self, [-1]);
        self._Hrepresentation.set_immutable()

        self._V_adjacency_matrix = matrix(ZZ, 0, 0, 0)
        self._V_adjacency_matrix.set_immutable()

        self._H_adjacency_matrix = matrix(ZZ, 1, 1, 0)
        self._H_adjacency_matrix.set_immutable()


    def _init_from_cdd_input(self, cdd_input_string,
                             cmdline_argument='--all', verbose=False):
        """
        Internal method: run cdd on a cdd H- or V-representation
        and initialize ourselves with the output.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: from sage.geometry.polyhedra import cdd_Vrepresentation
            sage: s = cdd_Vrepresentation('rational', [[0,0,1],[0,1,0],[1,0,0]], [], [])
            sage: p._init_from_cdd_input(s)
            sage: p.dim()
            2
        """
        if verbose: print cdd_input_string

        cdd_proc = Popen([self._cdd_executable, cmdline_argument],
                         stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = cdd_proc.communicate(input=cdd_input_string)

        # FIXME: check err
        if verbose: print ans

        self._init_from_cdd_output(ans)


    def _init_from_cdd_output(self, cdd_output_string):
        """
        Initialize ourselves with the output from cdd.

        TESTS::

            sage: from sage.geometry.polyhedra import cdd_Vrepresentation
            sage: s = cdd_Vrepresentation('rational',[[0,0],[1,0],[0,1],[1,1]], [], [])
            sage: from subprocess import Popen, PIPE
            sage: cdd_proc = Popen(['cdd_both_reps_gmp', '--all'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            sage: ans, err = cdd_proc.communicate(input=s)
            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,1],[1,1]])
            sage: p._init_from_cdd_output(ans)
            sage: p.vertices()
            [[0, 0], [1, 0], [0, 1], [1, 1]]
        """
        cddout=cdd_output_string.splitlines()
        suppressed_vertex = False   # whether cdd suppressed the vertex in output

        # nested function
        def expect_in_cddout(expected_string):
            l = cddout.pop(0).strip()
            if l!=expected_string:
                raise ValueError, ('Error while parsing cdd output: expected "'
                                   +expected_string+'" but got "'+l+'".\n' )
        # nested function
        def cdd_linearities():
            l = cddout[0].split()
            if l[0] != "linearity":
                return []
            cddout.pop(0)
            assert len(l) == int(l[1])+2, "Not enough linearities given"
            return [int(i)-1 for i in l[2:]]  # make indices pythonic

        # nested function
        def cdd_convert(string, field=self.field()):
            """
            Converts the cdd output string to a numerical value.
            """
            return [field(x) for x in string.split()]

        # nested function
        def find_in_cddout(expected_string):
            """
            Find the expected string in a list of strings, and
            truncates ``cddout`` to start at that point. Returns
            ``False`` if search fails.
            """
            for pos in range(0,len(cddout)):
                l = cddout[pos].strip();
                if l==expected_string:
                    # must not assign to cddout in nested function
                    for i in range(0,pos+1):
                        cddout.pop(0)
                    return True
            return False

        if find_in_cddout('V-representation'):
            self._Vrepresentation = Sequence([])
            lines = cdd_linearities()
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            self._ambient_dim = int(l[1])-1
            suppressed_vertex = True
            for i in range(int(l[0])):
                l = cddout.pop(0).strip()
                l_type = l[0]
                l = l[1:]
                if i in lines:
                    Line(self, cdd_convert(l));
                elif l_type == '0':
                    Ray(self, cdd_convert(l));
                else:
                    Vertex(self, cdd_convert(l));
                    suppressed_vertex = False
            if suppressed_vertex and self.n_Vrepresentation()>0:
                # cdd does not output the vertex if it is only the origin
                Vertex(self, [0] * self._ambient_dim)
            self._Vrepresentation.set_immutable()
            expect_in_cddout('end')

        if find_in_cddout('H-representation'):
            self._Hrepresentation = Sequence([])
            equations = cdd_linearities()
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert self._ambient_dim == int(l[1])-1, "Different ambient dimension?"
            for i in range(int(l[0])):
                l = cddout.pop(0)
                if i in equations:
                    Equation(self, cdd_convert(l));
                else:
                    Inequality(self, cdd_convert(l));
            self._Hrepresentation.set_immutable()
            expect_in_cddout('end')

        # nested function
        def cdd_adjacencies():
            l = cddout.pop(0).split()
            assert l[2] == ':', "Not a line of the adjacency data?"
            return [int(i)-1 for i in l[3:]]

        if find_in_cddout('Vertex graph'):
            n = len(self._Vrepresentation);
            if suppressed_vertex:
                n_cdd=n-1;
            else:
                n_cdd=n;
            self._V_adjacency_matrix = matrix(ZZ, n, n, 0)
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert int(l[0]) == n_cdd, "Not enough V-adjacencies in cdd output?"
            for i in range(n_cdd):
                for a in cdd_adjacencies():
                    self._V_adjacency_matrix[i,a] = 1
                # cdd reports that lines are never adjacent to anything.
                # I disagree, they are adjacent to everything!
                if self._Vrepresentation[i].is_line():
                    for j in range(n):
                        self._V_adjacency_matrix[i,j] = 1
                        self._V_adjacency_matrix[j,i] = 1
                    self._V_adjacency_matrix[i,i] = 0
            if suppressed_vertex: # cdd implied that there is only one vertex
                for i in range(n-1):
                    self._V_adjacency_matrix[i,n-1] = 1
                    self._V_adjacency_matrix[n-1,i] = 1
            self._V_adjacency_matrix.set_immutable()
            expect_in_cddout('end')

        if find_in_cddout('Facet graph'):
            n = len(self._Hrepresentation);
            self._H_adjacency_matrix = matrix(ZZ, n, n, 0)
            expect_in_cddout('begin')
            l = cddout.pop(0).split()
            assert int(l[0]) == n, "Not enough H-adjacencies in cdd output?"
            for i in range(n):
                for a in cdd_adjacencies():
                    self._H_adjacency_matrix[i,a] = 1
            self._H_adjacency_matrix.set_immutable()
            expect_in_cddout('end')


    def _init_facet_adjacency_matrix(self, verbose=False):
        """
        Compute the facet adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], compute_adjacency=False)
            sage: '_H_adjacency_matrix' in p.__dict__
            False
            sage: p._init_facet_adjacency_matrix()
            sage: p._H_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        self._init_from_cdd_input(self.cdd_Hrepresentation(),
                                  '--adjacency', verbose)


    def _init_vertex_adjacency_matrix(self, verbose=False):
        """
        Compute the vertex adjacency matrix in case it has not been
        computed during initialization.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[(0,0),(1,0),(0,1)], compute_adjacency=False)
            sage: '_V_adjacency_matrix' in p.__dict__
            False
            sage: p._init_vertex_adjacency_matrix()
            sage: p._V_adjacency_matrix
            [0 1 1]
            [1 0 1]
            [1 1 0]
        """
        self._init_from_cdd_input(self.cdd_Vrepresentation(),
                                  '--adjacency', verbose)


    def _repr_(self):
        """
        Return a description of the polyhedron.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1]])
            sage: poly_test._repr_()
            'A 2-dimensional polyhedron in QQ^4 defined as the convex hull of 3 vertices.\n'
            sage: grammar_test = Polyhedron(vertices = [[1,1,1,1,1,1]])
            sage: grammar_test._repr_()
            'A 0-dimensional polyhedron in QQ^6 defined as the convex hull of 1 vertex.\n'
        """
        desc = ''
        if self.n_vertices()==0:
            desc += 'The empty polyhedron'
        else:
            desc += 'A ' + repr(self.dim()) + '-dimensional polyhedron'
        desc += ' in '
        if self.field()==QQ: desc += 'QQ'
        else:                desc += 'RDF'
        desc += '^' + repr(self.ambient_dim())

        if self.n_vertices()>0:
            desc += ' defined as the convex hull of '
            desc += repr(self.n_vertices())
            if self.n_vertices()==1: desc += ' vertex'
            else:                    desc += ' vertices'

            if self.n_rays()>0:
                if self.n_lines()>0: desc += ", "
                else:                desc += " and "
                desc += repr(self.n_rays())
                if self.n_rays()==1: desc += ' ray'
                else:                desc += ' rays'

            if self.n_lines()>0:
                if self.n_rays()>0: desc += ", "
                else:               desc += " and "
                desc += repr(self.n_lines())
                if self.n_lines()==1: desc +=' line'
                else:                 desc +=' lines'

        return desc + ".\n";


    def cdd_Hrepresentation(self):
        """
        Write the inequalities/equations data of the polyhedron in
        cdd's H-representation format.

        OUTPUT:

        A string. If you save the output to filename.ine then you can
        run the stand-alone cdd via ``cddr+ filename.ine``

        EXAMPLES::

            sage: p = polytopes.n_cube(2)
            sage: print p.cdd_Hrepresentation()
            H-representation
            begin
             4 3 rational
             1 0 1
             1 1 0
             1 0 -1
             1 -1 0
            end
        """
        return cdd_Hrepresentation(self._cdd_type,
                                   [i for i in self.inequalities()],
                                   [e for e in self.equation_generator()] )


    def cdd_Vrepresentation(self):
        """
        Write the vertices/rays/lines data of the polyhedron in cdd's
        V-representation format.

        OUTPUT:

        A string. If you save the output to filename.ext then you can
        run the stand-alone cdd via ``cddr+ filename.ext``

        EXAMPLES::

            sage: q = Polyhedron(vertices = [[1,1],[0,0],[1,0],[0,1]])
            sage: print q.cdd_Vrepresentation()
            V-representation
            begin
             4 3 rational
             1 1 1
             1 0 0
             1 1 0
             1 0 1
            end
        """
        return cdd_Vrepresentation(self._cdd_type,
                            self.vertices(),
                            [r for r in self.ray_generator()],
                            [l for l in self.line_generator()] )


    def n_equations(self):
        """
        Return the number of equations. The representation will
        always be minimal, so the number of equations is the
        codimension of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_equations()
            1
        """
        try:
            return self._n_equations
        except AttributeError:
            self._n_equations = len(self.equations())
            return self._n_equations


    def n_inequalities(self):
        """
        Return the number of inequalities. The representation will
        always be minimal, so the number of inequalities is the
        number of facets of the polyhedron in the ambient space.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0,0],[0,1,0],[0,0,1]])
            sage: p.n_inequalities()
            3
        """
        try:
            return self._n_inequalities
        except AttributeError:
            self._n_inequalities = 0
            for i in self.inequalities(): self._n_inequalities += 1
            return self._n_inequalities


    def n_facets(self):
        """
        Return the number of facets in the polyhedron.  This is the
        same as the n_inequalities function.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        return self.n_inequalities()


    def n_vertices(self):
        """
        Return the number of vertices. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1],[1,1]], rays=[[1,1]])
            sage: p.n_vertices()
            2
        """
        try:
            return self._n_vertices
        except AttributeError:
            self._n_vertices = 0
            for v in self.vertex_generator(): self._n_vertices += 1
            return self._n_vertices


    def n_rays(self):
        """
        Return the number of rays. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1,0],[0,1]], rays=[[1,1]])
            sage: p.n_rays()
            1
        """
        try:
            return self._n_rays
        except AttributeError:
            self._n_rays = 0
            for r in self.rays(): self._n_rays += 1
            return self._n_rays


    def n_lines(self):
        """
        Return the number of lines. The representation will
        always be minimal.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0]], rays=[[0,1],[0,-1]])
            sage: p.n_lines()
            1
        """
        try:
            return self._n_lines
        except AttributeError:
            self._n_lines = len(self.lines())
            return self._n_lines


    def Hrepresentation(self, index=None):
        """
        Return the objects of the H-representaton. Each entry is
        either an inequality or a equation.

        INPUT:

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index in
        `0...n_Hrepresentations()`. If present, the H-representation
        object at the given index will be returned. Without an
        argument, returns the list of all H-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p.Hrepresentation(0)
            An inequality (0, 0, 1) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation() [0]
            True
        """
        if index==None:
            return self._Hrepresentation
        else:
            return self._Hrepresentation[index]


    def Hrep_generator(self):
        """
        Return an iterator over the objects of the H-representation
        (inequalities or equations).

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p.Hrep_generator().next()
            An inequality (0, 0, 1) x + 1 >= 0
        """
        for H in self.Hrepresentation():
            yield H


    def n_Hrepresentation(self):
        """
        Return the number of objects that make up the
        H-representation of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: p.n_Hrepresentation()
            16
            sage: p.n_Hrepresentation() == p.n_inequalities() + p.n_equations()
            True
        """
        return len(self.Hrepresentation())


    def Vrepresentation(self, index=None):
        """
        Return the objects of the V-representation. Each entry is
        either a vertex, a ray, or a line.

        INPUT:

        - ``index`` -- either an integer or ``None``.

        OUTPUT:

        The optional argument is an index in
        `0...n_Vrepresentation()`. If present, the V-representation
        object at the given index will be returned. Without an
        argument, returns the list of all V-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4)
            sage: p.Vrepresentation(0)
            A vertex at (0, 0, 0, -44721/50000)
            sage: p.Vrepresentation(0) == p.Vrepresentation() [0]
            True
        """
        if index==None:
            return self._Vrepresentation
        else:
            return self._Vrepresentation[index]


    def n_Vrepresentation(self):
        """
        Return the number of objects that make up the
        V-representation of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4)
            sage: p.n_Vrepresentation()
            5
            sage: p.n_Vrepresentation() == p.n_vertices() + p.n_rays() + p.n_lines()
            True
        """
        return len(self.Vrepresentation())


    def Vrep_generator(self):
        """
        Returns an iterator over the objects of the V-representation
        (vertices, rays, and lines).

        EXAMPLES::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: vg = p.Vrep_generator()
            sage: vg.next()
            A vertex at (0, 0, 0)
            sage: vg.next()
            A vertex at (1, 1, 1)
        """
        for V in self.Vrepresentation():
            yield V


    def facial_adjacencies(self):
        """
        Return the list of face indices (i.e. indices of
        H-representation objects) and the indices of faces adjacent to
        them.

        .. NOTE::

            Instead of working with face indices, you can use the
            H-representation objects directly (see example).

        EXAMPLES::

            sage: p = polytopes.permutahedron(4)
            sage: p.facial_adjacencies()[0:3]
            [[0, [1, 2, 12, 13]], [1, [0, 2, 3, 6, 8, 13]], [2, [0, 1, 3, 4, 7, 12]]]
            sage: f0 = p.Hrepresentation(0)
            sage: f0.index() == 0
            True
            sage: f0_adjacencies = [f0.index(), [n.index() for n in f0.neighbors()]]
            sage: p.facial_adjacencies()[0] == f0_adjacencies
            True
        """
        try:
            return self._facial_adjacencies
        except AttributeError:
            self._facial_adjacencies = \
                [ [ h.index(),
                    [n.index() for n in h.neighbors()]
                  ] for h in self.Hrepresentation() ]
            return self._facial_adjacencies


    def facial_incidences(self):
        """
        Return the face-vertex incidences in the form `[f_i, [v_{i_0}, v_{i_1},\dots ,v_{i_2}]]`.

        .. NOTE::

            Instead of working with face/vertex indices, you can use
            the H-representation/V-representation objects directly
            (see examples).

        OUTPUT:

        The face indices are the indices of the H-representation
        objects, and the vertex indices are the indices of the
        V-representation objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[0,0,0],[2,2,5]])
            sage: p.facial_incidences()
            [[0, [1, 3, 4]],
            [1, [0, 3, 4]],
            [2, [0, 2, 4]],
            [3, [1, 2, 4]],
            [4, [0, 1, 2, 3]]]
            sage: f0 = p.Hrepresentation(0)
            sage: f0.index() == 0
            True
            sage: f0_incidences = [f0.index(), [v.index() for v in f0.incident()]]
            sage: p.facial_incidences()[0] == f0_incidences
            True
        """
        try:
            return self._facial_incidences
        except AttributeError:
            self._facial_incidences = \
                [ [ h.index(),
                    [v.index() for v in h.incident()]
                  ] for h in self.Hrepresentation() ]
            return self._facial_incidences


    def vertex_adjacencies(self):
        """
        Return a list of vertex indices and their adjacent vertices.

        .. NOTE::

            Instead of working with vertex indices, you can use the
            V-representation objects directly (see examples).

        OUTPUT:

        The vertex indices are the indices of the V-representation
        objects.

        EXAMPLES::

            sage: permuta3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: permuta3.vertex_adjacencies()[0:3]
            [[0, [1, 2, 6]], [1, [0, 3, 7]], [2, [0, 4, 8]]]
            sage: v0 = permuta3.Vrepresentation(0)
            sage: v0.index() == 0
            True
            sage: list( v0.neighbors() )
            [A vertex at (1, 2, 4, 3), A vertex at (1, 3, 2, 4), A vertex at (2, 1, 3, 4)]
            sage: v0_adjacencies = [v0.index(), [v.index() for v in v0.neighbors()]]
            sage: permuta3.vertex_adjacencies()[0] == v0_adjacencies
            True
        """
        try:
            return self._vertex_adjacencies
        except AttributeError:
            self._vertex_adjacencies = \
                [ [ v.index(),
                    [n.index() for n in v.neighbors()]
                  ] for v in self.Vrepresentation() ]
            return self._vertex_adjacencies


    def vertex_incidences(self):
        """
        Return the vertex-face incidences in the form `[v_i, [f_{i_0}, f_{i_1},\dots ,f_{i_2}]]`.

        .. NOTE::

            Instead of working with face/vertex indices, you can use
            the H-representation/V-representation objects directly
            (see examples).

        EXAMPLES::

            sage: from sage.geometry.polyhedra import polytopes
            sage: p = polytopes.n_simplex(3)
            sage: p.vertex_incidences()
            [[0, [0, 1, 3]], [1, [0, 2, 3]], [2, [1, 2, 3]], [3, [0, 1, 2]]]
            sage: v0 = p.Vrepresentation(0)
            sage: v0.index() == 0
            True
            sage: p.vertex_incidences()[0] == [ v0.index(), [h.index() for h in v0.incident()] ]
            True
        """
        try:
            return self._vertex_incidences
        except AttributeError:
            self._vertex_incidences = \
                [ [ v.index(),
                    [h.index() for h in v.incident()]
                  ] for v in self.Vrepresentation() ]
            return self._vertex_incidences


    def inequality_generator(self):
        """
        Return  a generator for the defining inequalities of the
        polyhedron.

        OUTPUT:

        A generator of the inequality Hrepresentation objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.inequality_generator(): print(v)
            An inequality (-1, 0) x + 1 >= 0
            An inequality (0, -1) x + 1 >= 0
            An inequality (1, 1) x - 1 >= 0
            sage: [ v for v in triangle.inequality_generator() ]
            [An inequality (-1, 0) x + 1 >= 0,
             An inequality (0, -1) x + 1 >= 0,
             An inequality (1, 1) x - 1 >= 0]
            sage: [ [v.A(), v.b()] for v in triangle.inequality_generator() ]
            [[(-1, 0), 1], [(0, -1), 1], [(1, 1), -1]]
        """
        for H in self.Hrepresentation():
            if H.is_inequality():
                yield H


    def inequalities(self):
        """
        Return a list of inequalities as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`inequality_generator`
            instead to iterate over the list of :class:`Inequality`
            objects.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities()[0:3]
            [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]
            sage: p3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: ieqs = p3.inequalities()
            sage: ieqs[0]
            [-1, 0, 0, 1, 0]
            sage: ieqs[-1]
            [4, -1, 0, 0, 0]
            sage: ieqs == [list(x) for x in p3.inequality_generator()]
            True
        """
        try:
            return self._inequalities
        except AttributeError:
            self._ieqs = [list(x) for x in self.inequality_generator()]
            return self._ieqs


    def ieqs(self):
        """
        Deprecated. Alias for inequalities()

        EXAMPLES::

            sage: p3 = Polyhedron(vertices = permutations([1,2,3,4]))
            sage: p3.ieqs() == p3.inequalities()
            True
        """
        return self.inequalities()


    def equation_generator(self):
        """
        Return a generator for the linear equations satisfied by the
        polyhedron.

        EXAMPLES::

            sage: p = polytopes.regular_polygon(8,field=RDF)
            sage: p3 = Polyhedron(vertices = [x+[0] for x in p.vertices()], field=RDF)
            sage: p3.equation_generator().next()
            An equation (0.0, 0.0, 1.0) x + 0.0 == 0
        """
        for H in self.Hrepresentation():
            if H.is_equation():
                yield H


    def equations(self):
        """
        Return the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        .. NOTE::

            It is recommended to use :meth:`equation_generator()` instead
            to iterate over the list of :class:`Equation` objects.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.equations()
            [[-10, 1, 1, 1, 1]]
        """
        try:
            return self._equations
        except:
            self._equations = [list(eq) for eq in self.equation_generator()]
            return self._equations


    def linearities(self):
        """
        Deprecated.  Use equations() instead.
        Returns the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.linearities()
            [[-10, 1, 1, 1, 1]]
            sage: test_p.linearities() == test_p.equations()
            True
        """
        return self.equations()


    def vertices(self):
        """
        Return a list of vertices of the polyhedron.

        .. NOTE::

            It is recommended to use :meth:`vertex_generator` instead to
            iterate over the list of :class:`Vertex` objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: triangle.vertices()
            [[1, 0], [0, 1], [1, 1]]
            sage: a_simplex = Polyhedron(ieqs = [[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0],[0,0,0,0,1]], eqns = [[1,-1,-1,-1,-1]])
            sage: a_simplex.vertices()
            [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
            sage: a_simplex.vertices() == [list(v) for v in a_simplex.vertex_generator()]
            True
        """
        try:
            return self._vertices
        except:
            self._vertices = [list(x) for x in self.vertex_generator()]
            return self._vertices


    def vertex_generator(self):
        """
        Return a generator for the vertices of the polyhedron.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.vertex_generator(): print(v)
            A vertex at (1, 0)
            A vertex at (0, 1)
            A vertex at (1, 1)
            sage: v_gen = triangle.vertex_generator()
            sage: v_gen.next()   # the first vertex
            A vertex at (1, 0)
            sage: v_gen.next()   # the second vertex
            A vertex at (0, 1)
            sage: v_gen.next()   # the third vertex
            A vertex at (1, 1)
            sage: try: v_gen.next()   # there are only three vertices
            ... except StopIteration: print "STOP"
            STOP
            sage: type(v_gen)
            <type 'generator'>
            sage: [ v for v in triangle.vertex_generator() ]
            [A vertex at (1, 0), A vertex at (0, 1), A vertex at (1, 1)]

        """
        for V in self.Vrepresentation():
            if V.is_vertex():
                yield V


    def ray_generator(self):
        """
        Return a generator for the rays of the polyhedron.

        EXAMPLES::

            sage: pi = Polyhedron(ieqs = [[1,1,0],[1,0,1]])
            sage: pir = pi.ray_generator()
            sage: [x.vector() for x in pir]
            [(1, 0), (0, 1)]
        """
        for V in self.Vrepresentation():
            if V.is_ray():
                yield V


    def rays(self):
        """
        Return a list of rays as coefficient lists.

        .. NOTE::

            It is recommended to use :meth:`ray_generator` instead to
            iterate over the list of :class:`Ray` objects.

        OUTPUT:

        A list of rays as lists of coordinates.

        EXAMPLES::

            sage: p = Polyhedron(ieqs = [[0,0,0,1],[0,0,1,0],[1,1,0,0]])
            sage: p.rays()
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            sage: p.rays() == [list(r) for r in p.ray_generator()]
            True
        """
        try:
            return self._rays
        except:
            self._rays = [list(x) for x in self.ray_generator()]
            return self._rays


    def line_generator(self):
        """
        Return a generator for the lines of the polyhedron.

        EXAMPLES::

            sage: pr = Polyhedron(rays = [[1,0],[-1,0],[0,1]], vertices = [[-1,-1]])
            sage: pr.line_generator().next().vector()
            (1, 0)
        """
        for V in self.Vrepresentation():
            if V.is_line():
                yield V


    def lines(self):
        """
        Return a list of lines of the polyhedron.  The line data is given
        as a list of coordinates rather than as a Hrepresentation object.

        .. NOTE::

            It is recommended to use :meth:`line_generator` instead to
            iterate over the list of :class:`Line` objects.

        EXAMPLES::

            sage: p = Polyhedron(rays = [[1,0],[-1,0],[0,1],[1,1]], vertices = [[-2,-2],[2,3]])
            sage: p.lines()
            [[1, 0]]
            sage: p.lines() == [list(x) for x in p.line_generator()]
            True
        """
        try:
            return self._lines
        except:
            self._lines = [list(x) for x in self.line_generator()]
            return self._lines


    def bounded_edges(self):
        """
        Return the bounded edges (excluding rays and lines).

        OUTPUT:

        A generator for pairs of vertices, one pair per edge.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
            sage: [ e for e in p.bounded_edges() ]
            [(A vertex at (1, 0), A vertex at (0, 1))]
            sage: for e in p.bounded_edges(): print e
            (A vertex at (1, 0), A vertex at (0, 1))
        """
        obj = self.Vrepresentation()
        edges = []
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(i+1,len(obj)):
                if not obj[j].is_vertex(): continue
                if self.vertex_adjacency_matrix()[i,j] == 0: continue
                yield (obj[i], obj[j])


    def ambient_dim(self):
        r"""
        Return the dimension of the ambient space.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_dim()
            4
        """
        return self._ambient_dim


    def dim(self):
        """
        Return the dimension of the polyhedron.

        EXAMPLES::

            sage: simplex = Polyhedron(vertices = [[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
            sage: simplex.dim()
            3
            sage: simplex.ambient_dim()
            4
       """
        return self.ambient_dim() - self.n_equations()


    def adjacency_matrix(self):
        """
        This is an alias for :meth:`vertex_adjacency_matrix`

        EXAMPLES::

            sage: polytopes.n_cube(3).adjacency_matrix()
            [0 1 1 0 1 0 0 0]
            [1 0 0 1 0 1 0 0]
            [1 0 0 1 0 0 1 0]
            [0 1 1 0 0 0 0 1]
            [1 0 0 0 0 1 1 0]
            [0 1 0 0 1 0 0 1]
            [0 0 1 0 1 0 0 1]
            [0 0 0 1 0 1 1 0]
        """
        return self.vertex_adjacency_matrix()


    def vertex_adjacency_matrix(self):
        """
        Return the binary matrix of vertex adjacencies.

        EXAMPLES::

            sage: polytopes.n_simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        if '_V_adjacency_matrix' not in self.__dict__:
            self._init_vertex_adjacency_matrix()
        return self._V_adjacency_matrix;


    def facet_adjacency_matrix(self):
        """
        Return the adjacency matrix for the facets and hyperplanes.

        EXAMPLES::

            sage: polytopes.n_simplex(4).facet_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        if '_H_adjacency_matrix' not in self.__dict__:
            self._init_facet_adjacency_matrix()
        return self._H_adjacency_matrix;


    def incidence_matrix(self):
        """
        Return the incidence matrix.

        .. NOTE::

            The columns correspond to inequalities/equations in the
            order :meth:`Hrepresentation`, the rows correspond to
            vertices/rays/lines in the order
            :meth:`Vrepresentation`

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: p.incidence_matrix()
            [0 0 1 0 0 1 0 0 1 1 0 0 0 0]
            [0 0 0 1 0 0 0 0 0 1 1 0 1 0]
            [0 1 1 0 0 0 0 0 1 0 0 0 0 1]
            [0 0 0 0 0 0 0 0 0 0 1 1 1 1]
            [1 1 0 0 0 0 0 0 0 0 0 1 0 1]
            [0 0 0 0 0 0 0 0 1 1 1 0 0 1]
            [1 0 0 0 0 0 0 1 0 0 0 1 1 0]
            [1 1 1 0 0 0 1 0 0 0 0 0 0 0]
            [1 0 0 0 1 0 1 1 0 0 0 0 0 0]
            [0 0 0 1 1 0 0 1 0 0 0 0 1 0]
            [0 0 0 1 1 1 0 0 0 1 0 0 0 0]
            [0 0 1 0 1 1 1 0 0 0 0 0 0 0]
            sage: v = p.Vrepresentation(0)
            sage: v
            A vertex at (0, -1/2, -1/2)
            sage: h = p.Hrepresentation(2)
            sage: h
            An inequality (0, 2, 0) x + 1 >= 0
            sage: h.eval(v)        # evaluation (0, 2, 0) * (0, -1/2, -1/2) + 1
            0
            sage: h*v              # same as h.eval(v)
            0
            sage: p.incidence_matrix() [0,2]   # this entry is (v,h)
            1
            sage: h.contains(v)
            True
            sage: p.incidence_matrix() [2,0]   # note: not symmetric
            0
        """
        try:
            return self._incidence_matrix
        except AttributeError:
            self._incidence_matrix = matrix(ZZ, len(self.Vrepresentation()),
                                                len(self.Hrepresentation()), 0)
            for V in self.Vrep_generator():
                for H in self.Hrep_generator():
                    if self._is_zero(H*V):
                        self._incidence_matrix[V.index(),H.index()] = 1

            return self._incidence_matrix


    def field(self):
        """
        Return the number type that we are working with, either
        ``QQ`` (exact arithmetic using gmp, default) or ``RDF``
        (double precision floating-point arithmetic)

        EXAMPLES::

            sage: triangle = Polyhedron(vertices = [[1,0],[0,1],[1,1]])
            sage: triangle.field() == QQ
            True
        """
        return self._field


    def coerce_field(self, other):
        """
        Return the common field for both ``self`` and ``other``.

        INPUT:

        - ``other`` -- must be either:

            * another ``Polyhedron`` object

            * `\QQ` or `RDF`

            * a constant that can be coerced to `\QQ` or `RDF`

        OUTPUT:

        Either `\QQ` or `RDF`. Raises ``TypeError`` if ``other`` is not a
        suitable input.

        .. NOTE::

            "Real" numbers in sage are not necessarily elements of
            `RDF`. For example, the literal `1.0` is not.

        EXAMPLES::

            sage: triangle_QQ  = Polyhedron(vertices = [[1,0],[0,1],[1,1]], field=QQ)
            sage: triangle_RDF = Polyhedron(vertices = [[1,0],[0,1],[1,1]], field=RDF)
            sage: triangle_QQ.coerce_field(QQ)
            Rational Field
            sage: triangle_QQ.coerce_field(triangle_RDF)
            Real Double Field
            sage: triangle_RDF.coerce_field(triangle_QQ)
            Real Double Field
            sage: triangle_QQ.coerce_field(RDF)
            Real Double Field
            sage: triangle_QQ.coerce_field(ZZ)
            Rational Field
            sage: triangle_QQ.coerce_field(1/2)
            Rational Field
            sage: triangle_QQ.coerce_field(0.5)
            Real Double Field
        """
        try:
            # other is a Polyhedron object?
            other_field = other.field()
        except AttributeError:

            try:
                # other is a constant?
                other_parent = other.parent()
            except AttributeError:
                other_parent = other

            # other is a field?
            if QQ.coerce_map_from(other_parent) != None:
                other_field = QQ
            elif RDF.coerce_map_from(other_parent) != None:
                other_field = RDF
            else:
                raise TypeError("cannot determine field from %s!" % other)

        assert other_field==QQ or other_field==RDF

        if self.field()==RDF or other_field==RDF:
            return RDF
        else:
            return QQ


    def center(self):
        """
        Return the average of the vertices.

        OUTPUT:

        The center of the polyhedron. All rays and lines are
        ignored. Raises a ``ZeroDivisionError`` for the empty
        polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p = p + vector([1,0,0])
            sage: p.center()
            (1, 0, 0)
        """
        try:
            return self._center
        except AttributeError:
            self._center = vector(ZZ, [0]*self.ambient_dim())
            for v in self.vertex_generator(): self._center += v.vector()
            self._center /= self.n_vertices()
            return self._center


    def radius_square(self):
        """
        Return the square of the maximal distance from the center to
        a vertex. All rays and lines are ignored.

        OUTPUT:

        The square of the radius, which is in :meth:`field`.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            5
        """
        try:
            return self._radius_2
        except AttributeError:
            center = self.center()
            self._radius_2 = max( (v.vector() - center).dot_product(
                                   v.vector() - center) for v in
                                   self.vertex_generator() )
            return self._radius_2


    def radius(self):
        """
        Return the maximal distance from the center to a vertex. All
        rays and lines are ignored.

        OUTPUT:

        The radius for a rational polyhedron is, in general, not
        rational.  use :meth:`radius_square` if you need a rational
        distance measure.

        EXAMPLES::

            sage: p = polytopes.n_cube(4)
            sage: p.radius()
            2
        """
        return sqrt(self.radius_square())


    def is_compact(self):
        """
        Test for boundedness of the polytope.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: p.is_compact()
            True
            sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,-1,0,0]])
            sage: p.is_compact()
            False
        """
        return self.n_rays()==0 and self.n_lines()==0


    def is_simple(self):
        """
        Test for simplicity of a polytope.

        EXAMPLES::

            sage: p = Polyhedron([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
            sage: p.is_simple()
            True
            sage: p = Polyhedron([[0,0,0],[4,4,0],[4,0,0],[0,4,0],[2,2,2]])
            sage: p.is_simple()
            False

        REFERENCES:

            http://en.wikipedia.org/wiki/Simple_polytope
        """
        if not self.is_compact(): return False

        for v in self.vertex_generator():
            adj = [a for a in v.neighbors()]
            if len(adj) != self.dim():
                return False

        return True

    def gale_transform(self):
        """
        Return the Gale transform of a polytope as described in the
        reference below.

        OUTPUT:

        A list of vectors, the Gale transform.  The dimension is the
        dimension of the affine dependencies of the vertices of the
        polytope.

        EXAMPLES:

        This is from the reference, for a triangular prism::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0]])
            sage: p2 = p.prism()
            sage: p2.gale_transform()
            [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]

        REFERENCES:

            Lectures in Geometric Combinatorics, R.R.Thomas, 2006, AMS Press.
        """
        if not self.is_compact(): raise ValueError('Not a polytope.')

        A = matrix(self.n_vertices(),
                   [ [1]+list(x) for x in self.vertex_generator()])
        A = A.transpose()
        A_ker = A.right_kernel()
        return A_ker.basis_matrix().transpose().rows()

    def triangulate(self, engine='auto', connected=True, fine=False, regular=None, star=None):
        r"""
        Returns a triangulation of the polytope.

        INPUT:

        - ``engine`` -- either 'auto' (default), 'internal', or
          'TOPCOM'.  The latter two instruct this package to always
          use its own triangulation algorithms or TOPCOM's algorithms,
          respectively. By default ('auto'), TOPCOM is used if it is
          available and internal routines otherwise.

        The remaining keyword parameters are passed through to the
        :class:`~sage.geometry.triangulation.point_configuration.PointConfiguration`
        constructor:

        - ``connected`` -- boolean (default: ``True``). Whether the
          triangulations should be connected to the regular
          triangulations via bistellar flips. These are much easier to
          compute than all triangulations.

        - ``fine`` -- boolean (default: ``False``). Whether the
          triangulations must be fine, that is, make use of all points
          of the configuration.

        - ``regular`` -- boolean or ``None`` (default:
          ``None``). Whether the triangulations must be regular. A
          regular triangulation is one that is induced by a
          piecewise-linear convex support function. In other words,
          the shadows of the faces of a polyhedron in one higher
          dimension.

          * ``True``: Only regular triangulations.

          * ``False``: Only non-regular triangulations.

          * ``None`` (default): Both kinds of triangulation.

        - ``star`` -- either ``None`` (default) or a point. Whether
          the triangulations must be star. A triangulation is star if
          all maximal simplices contain a common point. The central
          point can be specified by its index (an integer) in the
          given points or by its coordinates (anything iterable.)

        OUTPUT:

        A triangulation of the convex hull of the vertices as a
        :class:`~sage.geometry.triangulation.point_configuration.Triangulation`. The
        indices in the triangulation correspond to the
        :meth:`Vrepresentation` objects.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: triangulation = cube.triangulate(engine='internal') # to make doctest independent of TOPCOM
            sage: triangulation
            (<0,1,2,7>, <0,1,4,7>, <0,2,4,7>, <1,2,3,7>, <1,4,5,7>, <2,4,6,7>)
            sage: simplex_indices = triangulation[0]; simplex_indices
            (0, 1, 2, 7)
            sage: simplex_vertices = [ cube.Vrepresentation(i) for i in simplex_indices ]
            sage: simplex_vertices
            [A vertex at (1, 1, 1), A vertex at (-1, 1, 1),
             A vertex at (1, -1, 1), A vertex at (-1, -1, -1)]
            sage: Polyhedron(simplex_vertices)
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 4 vertices.
        """
        if not self.is_compact():
            raise NotImplementedError('I can only triangulate compact polytopes.')
        from sage.geometry.triangulation.point_configuration import PointConfiguration
        pc = PointConfiguration((v.vector() for v in self.vertex_generator()),
                                connected=connected, fine=fine, regular=regular, star=star)
        pc.set_engine(engine)
        return pc.triangulate()

    def triangulated_facial_incidences(self):
        """
        Return a list of the form [face_index, [v_i_0,
        v_i_1,...,v_i_{n-1}]] where the face_index refers to the
        original defining inequality.  For a given face, the
        collection of triangles formed by each list of v_i should
        triangulate that face.

        In dimensions greater than 3, this is computed by randomly
        lifting each face up a dimension; this does not always work!
        This should eventually be fixed by using lrs or another
        program that computes triangulations.

        EXAMPLES:

        If the figure is already composed of triangles, then all is well::

            sage: Polyhedron(vertices = [[5,0,0],[0,5,0],[5,5,0],[2,2,5]]).triangulated_facial_incidences()
            doctest:...: DeprecationWarning: (Since Sage Version 4.7.1)
            This method is deprecated. Use triangulate() instead.
            [[0, [0, 2, 3]], [1, [0, 1, 2]], [2, [1, 2, 3]], [3, [0, 1, 3]]]

        Otherwise some faces get split up to triangles::

            p = polytopes.regular_polygon(5)
            sage: Polyhedron(vertices = [[2,0,0],[4,1,0],[0,5,0],[5,5,0],[1,1,0],[0,0,1]]).triangulated_facial_incidences()
            [[0, [0, 1, 5]], [1, [0, 4, 5]], [2, [2, 4, 5]], [3, [2, 3, 5]], [4, [1, 3, 5]], [5, [0, 1, 4]], [5, [1, 4, 3]], [5, [4, 3, 2]]]
        """
        from sage.misc.misc import deprecation
        deprecation('This method is deprecated. Use triangulate() instead.', 'Sage Version 4.7.1')
        try:
            return self._triangulated_facial_incidences
        except AttributeError:
            t_fac_incs = []
            for a_face in self.facial_incidences():
                vert_number = len(a_face[1])
                if vert_number == self.dim():
                    t_fac_incs.append(a_face)
                elif self.dim() >= 4:
                    lifted_verts = []
                    for vert_index in a_face[1]:
                        lifted_verts.append(self.vertices()[vert_index] +
                                            [randint(-vert_index,5000+vert_index + vert_number**2)])
                    temp_poly = Polyhedron(vertices = lifted_verts)
                    for t_face in temp_poly.facial_incidences():
                        if len(t_face[1]) != self.dim():
                            print 'Failed for face: ' + str(a_face)
                            print 'Attempted simplicial face: ' + str(t_face)
                            print 'Attempted lifted vertices: ' + str(lifted_verts)
                            raise RuntimeError, "triangulation failed"
                        normal_fdir = temp_poly.ieqs()[t_face[0]][-1]
                        if normal_fdir >= 0:
                            t_fac_verts = [temp_poly.vertices()[i] for i in t_face[1]]
                            proj_verts = [q[0:self.dim()] for q in t_fac_verts]
                            t_fac_incs.append([a_face[0],
                                               [self.vertices().index(q) for q in proj_verts]])
                else:
                    vs = a_face[1][:]
                    adj = dict([a[0], filter(lambda p: p in a_face[1], a[1])]
                               for a in filter(lambda va: va[0] in a_face[1],
                                               self.vertex_adjacencies()))
                    t = vs[0]
                    vs.remove(t)
                    ts = adj[t]
                    for v in ts:
                        vs.remove(v)
                    t_fac_incs.append([a_face[0], [t] + ts])
                    while vs:
                        t = ts[0]
                        ts = ts[1:]
                        for v in adj[t]:
                            if v in vs:
                                vs.remove(v)
                                ts.append(v)
                                t_fac_incs.append([a_face[0], [t] + ts])
                                break
        self._triangulated_facial_incidences = t_fac_incs
        return t_fac_incs


    def simplicial_complex(self):
        """
        Return a simplicial complex from a triangulation of the polytope.

        Warning: This first triangulates the polytope using
        ``triangulated_facial_incidences``, and this function may fail
        in dimensions greater than 3, although it usually doesn't.

        OUTPUT:

        A simplicial complex.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: sc = p.simplicial_complex()
            doctest:...: DeprecationWarning: (Since Sage Version 4.7.1)
            This method is deprecated. Use triangulate().simplicial_complex() instead.
            doctest:...: DeprecationWarning: (Since Sage Version 4.7.1)
            This method is deprecated. Use triangulate() instead.
            sage: sc
            Simplicial complex with 13 vertices and 20 facets
        """
        from sage.misc.misc import deprecation
        deprecation('This method is deprecated. Use triangulate().simplicial_complex() instead.',
                    'Sage Version 4.7.1')
        from sage.homology.simplicial_complex import SimplicialComplex
        return SimplicialComplex(vertex_set = self.n_vertices(),
                                 maximal_faces = [x[1] for x in self.triangulated_facial_incidences()])

    def __add__(self, other):
        """
        The Minkowski sum of ``self`` and ``other``.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        EXAMPLES::

            sage: four_cube = polytopes.n_cube(4)
            sage: four_simplex = Polyhedron(vertices = [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]])
            sage: unholy_union = four_cube + four_simplex
            sage: unholy_union.dim()
            4
            sage: poly_spam = Polyhedron([[3,4,5,2],[1,0,0,1],[0,0,0,0],[0,4,3,2],[-3,-3,-3,-3]])
            sage: poly_eggs = Polyhedron([[5,4,5,4],[-4,5,-4,5],[4,-5,4,-5],[0,0,0,0]])
            sage: poly_spam_and_eggs = poly_spam + poly_spam + poly_eggs
            sage: poly_spam_and_eggs.n_vertices()
            12
        """
        if isinstance(other,Polyhedron):
            new_vertices = []
            for v1 in self.vertex_generator():
                for v2 in other.vertex_generator():
                    new_vertices.append(list(v1() + v2()))
            new_rays = self.rays() + other.rays()
            new_lines = self.lines() + other.lines()
            other_field = other.field()

        else:  # assume other is a vector and try to add vertices
            displacement = vector(other)
            new_vertices = [list(x() + displacement) for x in self.vertex_generator()]
            new_rays = self.rays()
            new_lines = self.lines()
            other_field = displacement.base_ring()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          field=self.coerce_field(other_field))


    def __mul__(self, other):
        """
        Multiplication by ``other``.

        INPUT:

        - ``other`` -- A scalar, not necessarily in :meth:`field`, or
          a :class:`Polyhedron`.

        OUTPUT:

        Multiplication by another polyhedron returns the product
        polytope. Multiplication by a scalar returns the polytope
        dilated by that scalar, possibly coerced to the bigger field.

        EXAMPLES::

             sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,6)])
             sage: p.vertex_generator().next()
             A vertex at (2, 4, 8)
             sage: p2 = p*2
             sage: p2.vertex_generator().next()
             A vertex at (4, 8, 16)
        """
        if isinstance(other,Polyhedron):
            new_vertices = [ list(x)+list(y)
                             for x in self.vertex_generator() for y in other.vertex_generator()]
            new_rays = []
            new_rays.extend( [ list(r)+[0]*other.ambient_dim()
                               for r in self.ray_generator() ] )
            new_rays.extend( [ [0]*self.ambient_dim()+list(r)
                               for r in other.ray_generator() ] )
            new_lines = []
            new_lines.extend( [ list(l)+[0]*other.ambient_dim()
                                for l in self.line_generator() ] )
            new_lines.extend( [ [0]*self.ambient_dim()+list(l)
                               for l in other.line_generator() ] )
        else:
            new_vertices = [ list(other*v()) for v in self.vertex_generator()]
            new_rays =  self.rays()
            new_lines = self.lines()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines,
                          field=self.coerce_field(other))


    def __rmul__(self,other):
        """
        Right multiplication.

        See :meth:`__mul__` for details.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: p2 = 3*p + p
            sage: p2.vertex_generator().next()
            A vertex at (8, 16, 32)
        """
        return self.__mul__(other)


    def union(self, other):
        """
        Deprecated.  Use ``self.convex_hull(other)`` instead.

        EXAMPLES::

            sage: Polyhedron(vertices=[[0]]).union( Polyhedron(vertices=[[1]]) )
            doctest:...: DeprecationWarning: (Since Sage Version 4.4.4) The function union is replaced by convex_hull.
            A 1-dimensional polyhedron in QQ^1 defined as the convex hull of 2 vertices.
        """
        from sage.misc.misc import deprecation
        deprecation('The function union is replaced by convex_hull.', 'Sage Version 4.4.4')
        return self.convex_hull(other)


    def convex_hull(self, other):
        """
        Return the convex hull of the set-theoretic union of the two
        polyhedra.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        The convex hull.

        EXAMPLES::

            sage: a_simplex = polytopes.n_simplex(3)
            sage: verts = a_simplex.vertices()
            sage: verts = [[x[0]*3/5+x[1]*4/5, -x[0]*4/5+x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices = verts)
            sage: simplex_union = a_simplex.convex_hull(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        hull_vertices = self.vertices() + other.vertices()
        hull_rays = self.rays() + other.rays()
        hull_lines = self.lines() + other.lines()
        hull_field = self.coerce_field(other)
        return Polyhedron(vertices=hull_vertices,
                          rays=hull_rays, lines=hull_lines,
                          field=hull_field)


    def intersection(self, other):
        """
        Return the intersection of one polyhedron with another.

        INPUT:

        - ``other`` -- a :class:`Polyhedron`.

        OUTPUT:

        The intersection.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: oct = polytopes.cross_polytope(3)
            sage: cube_oct = cube.intersection(oct*2)
            sage: len(list( cube_oct.vertex_generator() ))
            12
            sage: cube_oct
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 12 vertices.
        """
        new_ieqs = []
        new_ieqs.extend(self.inequalities())
        new_ieqs.extend(other.inequalities())

        new_eqns = []
        new_eqns.extend(self.equations())
        new_eqns.extend(other.equations())

        return Polyhedron(ieqs = new_ieqs, eqns = new_eqns,
                          field=self.coerce_field(other))


    def edge_truncation(self, cut_frac = Integer(1)/3):
        r"""
        Return a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

        - ``cut_frac`` -- integer. how deeply to cut into the edge.
            Default is `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: trunc_cube = cube.edge_truncation()
            sage: trunc_cube.n_vertices()
            24
            sage: trunc_cube.n_inequalities()
            14
        """
        new_vertices = []
        for e in self.bounded_edges():
            new_vertices.append((1-cut_frac)*e[0]() + cut_frac *e[1]())
            new_vertices.append(cut_frac *e[0]() + (1-cut_frac)*e[1]())

        new_vertices = [list(v) for v in new_vertices]
        new_rays =  self.rays()
        new_lines = self.lines()

        return Polyhedron(vertices=new_vertices, rays=new_rays,
                          lines=new_lines,
                          field=self.coerce_field(cut_frac))


    def face_lattice(self):
        """
        Return the face-lattice poset.

        OUTPUT:

        A :class:`~sage.combinat.posets.posets.FinitePoset`. Elements
        are given as :class:`~sage.geometry.polyhedra.PolyhedronFace`.

        In the case of a full-dimensional polytope, the faces are
        pairs (vertices, inequalities) of the spanning vertices and
        corresponding saturated inequalities. In general, a face is
        defined by a pair (V-rep. objects, H-rep. objects). The
        V-representation objects span the face, and the corresponding
        H-representation objects are those inequalities and equations
        that are saturated on the face.

        The bottom-most element of the face lattice is the "empty
        face". It contains no V-representation object. All
        H-representation objects are incident.

        The top-most element is the "full face". It is spanned by all
        V-representation objects. The incident H-representation
        objects are all equations and no inequalities.

        In the case of a full-dimensional polytope, the "empty face"
        and the "full face" are the empty set (no vertices, all
        inequalities) and the full polytope (all vertices, no
        inequalities), respectively.

        ALGORITHM:

        For a full-dimensional polytope, the basic algorithm is
        described in :func:`Hasse_diagram_from_incidences`. There are
        three generalizations of [KP2002]_ necessary to deal with more
        general polytopes, corresponding to the extra
        H/V-representation objects:

        * Lines are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face V-representation except for the "empty face".

        * Equations are removed before calling
          :func:`Hasse_diagram_from_incidences`, and then added back
          to each face H-representation.

        * Rays: Consider the half line as an example. The
          V-representation objects are a point and a ray, which we can
          think of as a point at infinity. However, the point at
          infinity has no inequality associated to it, so there is
          only one H-representation object alltogether. The face
          lattice does not contain the "face at infinity". This means
          that in :func:`Hasse_diagram_from_incidences`, one needs to
          drop faces with V-representations that have no matching
          H-representation. In addition, one needs to ensure that
          every non-empty face contains at least one vertex.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: square.face_lattice()
            Finite poset containing 10 elements
            sage: list(_)
            [<>, <0>, <1>, <2>, <3>, <2,3>, <1,3>, <0,1>, <0,2>, <0,1,2,3>]
            sage: poset_element = _[6]
            sage: a_face = poset_element.element
            sage: a_face
            <1,3>
            sage: a_face.dim()
            1
            sage: set(a_face.Vrepresentation()) == set([square.Vrepresentation(1), square.Vrepresentation(3)])
            True
            sage: a_face.Vrepresentation()
            [A vertex at (-1, 1), A vertex at (-1, -1)]
            sage: a_face.Hrepresentation()
            [An inequality (1, 0) x + 1 >= 0]

        A more complicated example::

            sage: c5_10 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,11)])
            sage: c5_10_fl = c5_10.face_lattice()
            sage: [len(x) for x in c5_10_fl.level_sets()]
            [1, 10, 45, 100, 105, 42, 1]

        Note that if the polyhedron contains lines then there is a
        dimension gap between the empty face and the first non-empty
        face in the face lattice::

            sage: line = Polyhedron(vertices=[(0,)], lines=[(1,)])
            sage: [ fl.element.dim() for fl in line.face_lattice() ]
            [-1, 1]

        TESTS::

            sage: c5_20 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,21)]) # not tested - very long time
            sage: c5_20_fl = c5_20.face_lattice() # not tested - very long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # not tested - very long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.n_cube(2).face_lattice().plot()
            sage: level_sets = polytopes.cross_polytope(2).face_lattice().level_sets()
            sage: print level_sets[0], level_sets[-1]
            [<>] [<0,1,2,3>]

        Various degenerate polyhedra::

            sage: Polyhedron(vertices=[[0,0,0],[1,0,0],[0,1,0]]).face_lattice().level_sets()
            [[<>], [<0>, <1>, <2>], [<0,1>, <1,2>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], rays=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<1>, <2>], [<0,1>, <0,2>, <1,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0,0),(0,1,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<2>], [<1,2>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(rays=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<2>], [<1,2>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(1,),(0,)]).face_lattice().level_sets()
            [[<>], [<0>, <1>], [<0,1>]]
            sage: Polyhedron(vertices=[(1,0,0),(0,1,0)], lines=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>, <0,2>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0,0)], vertices=[(0,0,1)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0),(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1,2>]]
            sage: Polyhedron(lines=[(1,0)], rays=[(0,1)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,2>], [<0,1,2>]]
            sage: Polyhedron(vertices=[(0,)], lines=[(1,)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
            sage: Polyhedron(lines=[(1,0)], vertices=[(0,0)]).face_lattice().level_sets()
            [[<>], [<0,1>]]
        """
        try:
            return self._face_lattice
        except AttributeError:
            pass

        coatom_to_Hindex = [ h.index() for h in self.inequality_generator() ]
        Hindex_to_coatom = [None] * self.n_Hrepresentation()
        for i in range(0,len(coatom_to_Hindex)):
            Hindex_to_coatom[ coatom_to_Hindex[i] ] = i

        atom_to_Vindex = [ v.index() for v in self.Vrep_generator() if not v.is_line() ]
        Vindex_to_atom = [None] * self.n_Vrepresentation()
        for i in range(0,len(atom_to_Vindex)):
                        Vindex_to_atom[ atom_to_Vindex[i] ] = i

        atoms_incidences   = [ tuple([ Hindex_to_coatom[h.index()]
                                       for h in v.incident() if h.is_inequality() ])
                               for v in self.Vrepresentation() if not v.is_line() ]

        coatoms_incidences = [ tuple([ Vindex_to_atom[v.index()]
                                       for v in h.incident() if not v.is_line() ])
                               for h in self.Hrepresentation() if h.is_inequality() ]

        atoms_vertices = [ Vindex_to_atom[v.index()] for v in self.vertex_generator() ]
        equations = [ e.index() for e in self.equation_generator() ]
        lines     = [ l.index() for l in self.line_generator() ]

        def face_constructor(atoms,coatoms):
            if len(atoms)==0:
                Vindices = ()
            else:
                Vindices = tuple(sorted([   atom_to_Vindex[i] for i in   atoms ]+lines))
            Hindices = tuple(sorted([ coatom_to_Hindex[i] for i in coatoms ]+equations))
            return PolyhedronFace(self, Vindices, Hindices)

        self._face_lattice = Hasse_diagram_from_incidences\
            (atoms_incidences, coatoms_incidences, face_constructor=face_constructor, required_atoms=atoms_vertices)
        return self._face_lattice


    def f_vector(self):
        r"""
        Return the f-vector.

        OUTPUT:

        Returns a vector whose ``i``-th entry is the number of
        ``i``-dimensional faces of the polytope.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            (1, 7, 12, 7, 1)
        """
        try:
            return self._f_vector
        except AttributeError:
            self._f_vector = vector(ZZ,[len(x) for x in self.face_lattice().level_sets()])
            return self._f_vector


    def vertex_graph(self):
        """
        Return a graph in which the vertices correspond to vertices
        of the polyhedron, and edges to edges.

        EXAMPLES::

            sage: g3 = polytopes.n_cube(3).vertex_graph()
            sage: len(g3.automorphism_group())
            48
            sage: s4 = polytopes.n_simplex(4).vertex_graph()
            sage: s4.is_eulerian()
            True
        """
        try:
            return self._graph
        except AttributeError:
            self._graph = Graph(self.vertex_adjacency_matrix(), loops=True)
            return self._graph


    graph = vertex_graph


    def polar(self):
        """
        Return the polar (dual) polytope.  The original vertices are
        translated so that their barycenter is at the origin, and then
        the vertices are used as the coefficients in the polar inequalities.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,0],[1,0,0],[0,0,0],[1,1,1]])
            sage: p
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 5 vertices.
            sage: p.polar()
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 6 vertices.
        """
        assert self.is_compact(), "Not a polytope."

        verts = [list(v() - self.center()) for v in self.vertex_generator()]
        return Polyhedron(ieqs = [[1] + list(v) for v in verts],
                          field = self.field())


    def pyramid(self):
        """
        Returns a polyhedron that is a pyramid over the original.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: egyptian_pyramid = square.pyramid()
            sage: egyptian_pyramid.n_vertices()
            5
            sage: for v in egyptian_pyramid.vertex_generator(): print v
            A vertex at (0, 1, 1)
            A vertex at (0, -1, 1)
            A vertex at (0, 1, -1)
            A vertex at (0, -1, -1)
            A vertex at (1, 0, 0)
        """
        new_verts = \
            [[0] + list(x) for x in self.Vrep_generator()] + \
            [[1] + list(self.center())]

        return Polyhedron(vertices = new_verts, field=self.field())


    def bipyramid(self):
        """
        Return a polyhedron that is a bipyramid over the original.

        EXAMPLES::

            sage: octahedron = polytopes.cross_polytope(3)
            sage: cross_poly_4d = octahedron.bipyramid()
            sage: cross_poly_4d.n_vertices()
            8
            sage: q = [list(v) for v in cross_poly_4d.vertex_generator()]
            sage: q
            [[0, 0, 0, 1],
             [0, 0, 1, 0],
             [0, 1, 0, 0],
             [0, 0, 0, -1],
             [0, 0, -1, 0],
             [0, -1, 0, 0],
             [1, 0, 0, 0],
             [-1, 0, 0, 0]]

        Now check that bipyramids of cross-polytopes are cross-polytopes::

            sage: q2 = [list(v) for v in polytopes.cross_polytope(4).vertex_generator()]
            sage: [v in q2 for v in q]
            [True, True, True, True, True, True, True, True]
        """
        new_verts = \
            [[ 0] + list(x) for x in self.vertex_generator()] + \
            [[ 1] + list(self.center())] + \
            [[-1] + list(self.center())]
        new_rays = [[0] + r for r in self.rays()]
        new_lines = [[0] + list(l) for l in self.lines()]
        return Polyhedron(vertices=new_verts,
                          rays=new_rays, lines=new_lines, field=self.field())


    def prism(self):
        """
        Return a prism of the original polyhedron.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: cube = square.prism()
            sage: cube
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 8 vertices.
            sage: hypercube = cube.prism()
            sage: hypercube.n_vertices()
            16
        """
        new_verts = []
        new_verts.extend( [ [0] + v for v in self.vertices()] )
        new_verts.extend( [ [1] + v for v in self.vertices()] )
        new_rays =        [ [0] + r for r in self.rays()]
        new_lines =       [ [0] + list(l) for l in self.lines()]
        return Polyhedron(vertices=new_verts,
                          rays=new_rays, lines=new_lines, field=self.field())


    def projection(self):
        """
        Return a projection object.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj
            The projection of a polyhedron into 3 dimensions.
        """
        self.projection = Projection(self)
        return self.projection


    def render_solid(self, **kwds):
        """
        Return a solid rendering of a 2- or 3-d polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p_solid = p.render_solid(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        proj = self.projection()
        if self.ambient_dim()==3:
            return proj.render_solid_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_fill_2d(**kwds)
        raise ValueError, "render_solid is only defined for 2 and 3 dimensional polyhedra."


    def render_wireframe(self, **kwds):
        """
        For polytopes in 2 or 3 dimensions, return the edges
        as a list of lines.

        EXAMPLES::

            sage: p = Polyhedron([[1,2,],[1,1],[0,0]])
            sage: p_wireframe = p.render_wireframe()
            sage: p_wireframe._Graphics__objects
            [Line defined by 2 points, Line defined by 2 points, Line defined by 2 points]
        """
        proj = self.projection()
        if self.ambient_dim()==3:
            return proj.render_wireframe_3d(**kwds)
        if self.ambient_dim()==2:
            return proj.render_outline_2d(**kwds)
        raise ValueError, "render_wireframe is only defined for 2 and 3 dimensional polyhedra."


    def schlegel_projection(self, projection_dir = None, height = 1.1):
        """
        Returns a projection object whose transformed coordinates are
        a Schlegel projection of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: sch_proj = p.schlegel_projection()
            sage: schlegel_edge_indices = sch_proj.lines
            sage: schlegel_edges = [sch_proj.coordinates_of(x) for x in schlegel_edge_indices]
            sage: len([x for x in schlegel_edges if x[0][0] > 0])
            8
        """
        proj = self.projection()
        if projection_dir == None:
            v = self.vertices()
            f0 = (self.facial_incidences()[0])[1]
            projection_dir = [sum([v[f0[i]][j]/len(f0) for i in range(len(f0))])
                              for j in range(len(v[0]))]
        return proj.schlegel(projection_direction = projection_dir, height = height)


    def lrs_volume(self, verbose = False):
        """
        Computes the volume of a polytope.

        OUTPUT:

        The volume, cast to RDF (although lrs seems to output a
        rational value this must be an approximation in some cases).

        EXAMPLES::

            sage: polytopes.n_cube(3).lrs_volume() #optional, needs lrs package installed
            8.0
            sage: (polytopes.n_cube(3)*2).lrs_volume() #optional, needs lrs package installed
            64.0
            sage: polytopes.twenty_four_cell().lrs_volume() #optional, needs lrs package installed
            2.0

        REFERENCES:

             David Avis's lrs program.
        """
        if is_package_installed('lrs') != True:
            print 'You must install the optional lrs package ' \
                  'for this function to work'
            raise NotImplementedError

        in_str = self.cdd_Vrepresentation()
        in_str += 'volume'
        in_filename = tmp_filename()
        in_file = file(in_filename,'w')
        in_file.write(in_str)
        in_file.close()
        if verbose: print in_str

        lrs_procs = Popen(['lrs',in_filename],
                          stdin = PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = lrs_procs.communicate()
        if verbose:
            print ans
        # FIXME: check err

        for a_line in ans.splitlines():
            if 'Volume=' in a_line:
                volume = a_line.split('Volume=')[1]
                volume = RDF(QQ(volume))
                return volume

        raise ValueError, "lrs did not return a volume"


    def contains(self, point):
        """
        Test whether the polyhedron contains the given ``point``.

        See also :meth:`interior_contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point (an iterable).

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[1,1],[1,-1],[0,0]])
            sage: P.contains( [1,0] )
            True
            sage: P.contains( P.center() )  # true for any convex set
            True

        The point need not have coordinates in the same field as the
        polyhedron::

            sage: ray = Polyhedron(vertices=[(0,0)], rays=[(1,0)], field=QQ)
            sage: ray.contains([sqrt(2)/3,0])        # irrational coordinates are ok
            True
            sage: a = var('a')
            sage: ray.contains([a,0])                # a might be negative!
            False
            sage: assume(a>0)
            sage: ray.contains([a,0])
            True
            sage: ray.contains(['hello', 'kitty'])   # no common ring for coordinates
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0.
            sage: empty.contains([])
            False
            sage: empty.contains([0])               # not a point in QQ^0
            False
            sage: full = Polyhedron(vertices=[()]); full
            A 0-dimensional polyhedron in QQ^0 defined as the convex hull of 1 vertex.
            sage: full.contains([])
            True
            sage: full.contains([0])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.contains(p):
                return False
        return True


    def interior_contains(self, point):
        """
        Test whether the interior of the polyhedron contains the
        given ``point``.

        See also :meth:`contains` and
        :meth:`relative_interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[[0,0],[1,1],[1,-1]])
            sage: P.contains( [1,0] )
            True
            sage: P.interior_contains( [1,0] )
            False

        If the polyhedron is of strictly smaller dimension than the
        ambient space, its interior is empty::

            sage: P = Polyhedron(vertices=[[0,1],[0,-1]])
            sage: P.contains( [0,0] )
            True
            sage: P.interior_contains( [0,0] )
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0.
            sage: empty.interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for H in self.Hrep_generator():
            if not H.interior_contains(p):
                return False
        return True


    def relative_interior_contains(self, point):
        """
        Test whether the relative interior of the polyhedron
        contains the given ``point``.

        See also :meth:`contains` and :meth:`interior_contains`.

        INPUT:

        - ``point`` -- coordinates of a point.

        OUTPUT:

        ``True`` or ``False``.

        EXAMPLES::

            sage: P = Polyhedron(vertices=[(1,0), (-1,0)])
            sage: P.contains( (0,0) )
            True
            sage: P.interior_contains( (0,0) )
            False
            sage: P.relative_interior_contains( (0,0) )
            True
            sage: P.relative_interior_contains( (1,0) )
            False

        The empty polyhedron needs extra care, see trac #10238::

            sage: empty = Polyhedron(); empty
            The empty polyhedron in QQ^0.
            sage: empty.relative_interior_contains([])
            False
        """
        try:
            p = vector(point)
        except TypeError: # point not iterable or no common ring for elements
            if len(point)>0:
                return False
            else:
                p = vector(self.field(), [])

        if len(p)!=self.ambient_dim():
            return False

        for eq in self.equation_generator():
            if not eq.contains(p):
                return False

        for ine in self.inequality_generator():
            if not ine.interior_contains(p):
                return False

        return True

    def is_simplex(self):
        r"""
        Return whether the polyhedron is a simplex.

        EXAMPLES::

            sage: Polyhedron([(0,0,0), (1,0,0), (0,1,0)]).is_simplex()
            True
            sage: polytopes.n_simplex(3).is_simplex()
            True
            sage: polytopes.n_cube(3).is_simplex()
            False
        """
        return self.is_compact() and (self.dim()+1==self.n_vertices())

    def is_lattice_polytope(self):
        r"""
        Return whether the polyhedron is a lattice polytope.

        OUTPUT:

        ``True`` if the polyhedron is compact and has only integral
        vertices, ``False`` otherwise.

        EXAMPLES::

            sage: polytopes.cross_polytope(3).is_lattice_polytope()
            True
            sage: polytopes.regular_polygon(5).is_lattice_polytope()
            False
        """
        try:
            return self._is_lattice_polytope
        except AttributeError:
            pass
        self._is_lattice_polytope = self.is_compact() and \
            all(v.is_integral() for v in self.vertex_generator())
        return self._is_lattice_polytope

    def lattice_polytope(self, envelope=False):
        r"""
        Return an encompassing lattice polytope.

        INPUT:

        - ``envelope`` -- boolean (default: ``False``). If the
          polyhedron has non-integral vertices, this option decides
          whether to return a strictly larger lattice polytope or
          raise a ``ValueError``. This option has no effect if the
          polyhedron has already integral vertices.

        OUTPUT:

        A :class:`LatticePolytope
        <sage.geometry.lattice_polytope.LatticePolytopeClass>`. If the
        polyhedron is compact and has integral vertices, the lattice
        polytope equals the polyhedron. If the polyhedron is compact
        but has at least one non-integral vertex, a strictly larger
        lattice polytope is returned.

        If the polyhedron is not compact, a ``NotImplementedError`` is
        raised.

        If the polyhedron is not integral and ``envelope=False``, a
        ``ValueError`` is raised.

        ALGORITHM:

        For each non-integral vertex, a bounding box of integral
        points is added and the convex hull of these integral points
        is returned.

        EXAMPLES:

        First, a polyhedron with integral vertices::

            sage: P = Polyhedron( vertices = [(1, 0), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope(); lp
            A lattice polytope: 2-dimensional, 4 vertices.
            sage: lp.vertices()
            [ 1  0 -1  0]
            [ 0  1  0 -1]

        Here is a polyhedron with non-integral vertices::

            sage: P = Polyhedron( vertices = [(1/2, 1/2), (0, 1), (-1, 0), (0, -1)])
            sage: lp = P.lattice_polytope()
            Traceback (most recent call last):
            ...
            ValueError: Some vertices are not integral. You probably want
            to add the argument "envelope=True" to compute an enveloping
            lattice polytope.
            sage: lp = P.lattice_polytope(True); lp
            A lattice polytope: 2-dimensional, 5 vertices.
            sage: lp.vertices()
            [ 0  1  1 -1  0]
            [ 1  0  1  0 -1]
        """
        if not self.is_compact():
            raise NotImplementedError, 'Only compact lattice polytopes are allowed.'

        def nonintegral_error():
            raise ValueError, 'Some vertices are not integral. '+\
                'You probably want to add the argument '+\
                '"envelope=True" to compute an enveloping lattice polytope.'

        # try to make use of cached values, if possible
        if envelope:
            try:
                return self._lattice_polytope
            except AttributeError:
                pass
        else:
            try:
                assert self._is_lattice_polytope
                return self._lattice_polytope
            except AttributeError:
                pass
            except AssertionError:
                nonintegral_error()

        # find the integral vertices
        try:
            vertices = matrix(ZZ, self.vertices()).transpose()
            self._is_lattice_polytope = True
        except TypeError:
            self._is_lattice_polytope = False
            if envelope==False: nonintegral_error()
            vertices = []
            for v in self.vertex_generator():
                vbox = [ set([floor(x),ceil(x)]) for x in v ]
                vertices.extend( CartesianProduct(*vbox) )
            vertices = matrix(ZZ, vertices).transpose()

        # construct the (enveloping) lattice polytope
        from sage.geometry.lattice_polytope import LatticePolytope
        self._lattice_polytope = LatticePolytope(vertices)
        return self._lattice_polytope

    def _integral_points_PALP(self):
        r"""
        Return the integral points in the polyhedron using PALP.

        This method is for testing purposes and will eventually be removed.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [(-1, -1), (1, 0), (1, 1), (0, 1), (0, 0)]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)]).lattice_polytope(True).points()
            [ 0 -1 -1  1  1  0  0]
            [-1  0 -1  0  1  1  0]
            sage: Polyhedron(vertices=[(-1/2,-1/2),(1,0),(1,1),(0,1)])._integral_points_PALP()
            [(1, 0), (1, 1), (0, 1), (0, 0)]
        """
        if not self.is_compact():
            raise ValueError, 'Can only enumerate points in a compact polyhedron.'
        lp = self.lattice_polytope(True)
        # remove cached values to get accurate timings
        try:
            del lp._points
            del lp._npoints
        except AttributeError:
            pass
        if self.is_lattice_polytope():
            return lp.points().columns()
        points = filter(lambda p: self.contains(p),
                        lp.points().columns())
        return points

    @cached_method
    def bounding_box(self, integral=False):
        r"""
        Return the coordinates of a rectangular box containing the non-empty polytope.

        INPUT:

        - ``integral`` -- Boolean (default: ``False``). Whether to
          only allow integral coordinates in the bounding box.

        OUTPUT:

        A pair of tuples ``(box_min, box_max)`` where ``box_min`` are
        the coordinates of a point bounding the coordinates of the
        polytope from below and ``box_max`` bounds the coordinates
        from above.

        EXAMPLES::

            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box()
            ((1/3, 1/3), (2/3, 2/3))
            sage: Polyhedron([ (1/3,2/3), (2/3, 1/3) ]).bounding_box(integral=True)
            ((0, 0), (1, 1))
            sage: polytopes.buckyball().bounding_box()
            ((-1059/1309, -1059/1309, -1059/1309), (1059/1309, 1059/1309, 1059/1309))
        """
        box_min = []
        box_max = []
        if self.n_vertices==0:
            raise ValueError('Empty polytope is not allowed')
        for i in range(0,self.ambient_dim()):
            coords = [ v[i] for v in self.Vrep_generator() ]
            max_coord = max(coords)
            min_coord = min(coords)
            if integral:
                box_max.append(ceil(max_coord))
                box_min.append(floor(min_coord))
            else:
                box_max.append(max_coord)
                box_min.append(min_coord)
        return (tuple(box_min), tuple(box_max))

    def integral_points(self, threshold=100000):
        r"""
        Return the integral points in the polyhedron.

        Uses either the naive algorithm (iterate over a rectangular
        bounding box) or triangulation + Smith form.

        INPUT:

        - ``threshold`` -- integer (default: 100000). Use the naive
          algorith as long as the bounding box is smaller than this.

        OUTPUT:

        The list of integral points in the polyhedron. If the
        polyhedron is not compact, a ``ValueError`` is raised.

        EXAMPLES::

            sage: Polyhedron(vertices=[(-1,-1),(1,0),(1,1),(0,1)]).integral_points()
            ((-1, -1), (0, 0), (0, 1), (1, 0), (1, 1))

            sage: simplex = Polyhedron([(1,2,3), (2,3,7), (-2,-3,-11)])
            sage: simplex.integral_points()
            ((-2, -3, -11), (0, 0, -2), (1, 2, 3), (2, 3, 7))

        The polyhedron need not be full-dimensional::

            sage: simplex = Polyhedron([(1,2,3,5), (2,3,7,5), (-2,-3,-11,5)])
            sage: simplex.integral_points()
            ((-2, -3, -11, 5), (0, 0, -2, 5), (1, 2, 3, 5), (2, 3, 7, 5))

            sage: point = Polyhedron([(2,3,7)])
            sage: point.integral_points()
            ((2, 3, 7),)

            sage: empty = Polyhedron()
            sage: empty.integral_points()
            ()

        Here is a simplex where the naive algorithm of running over
        all points in a rectangular bounding box no longer works fast
        enough::

            sage: v = [(1,0,7,-1), (-2,-2,4,-3), (-1,-1,-1,4), (2,9,0,-5), (-2,-1,5,1)]
            sage: simplex = Polyhedron(v); simplex
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 5 vertices.
            sage: len(simplex.integral_points())
            49

        Finally, the 3-d reflexive polytope number 4078::

            sage: v = [(1,0,0), (0,1,0), (0,0,1), (0,0,-1), (0,-2,1),
            ...        (-1,2,-1), (-1,2,-2), (-1,1,-2), (-1,-1,2), (-1,-3,2)]
            sage: P = Polyhedron(v)
            sage: pts1 = P.integral_points()                     # Sage's own code
            sage: all(P.contains(p) for p in pts1)
            True
            sage: pts2 = LatticePolytope(v).points().columns()   # PALP
            sage: for p in pts1: p.set_immutable()
            sage: for p in pts2: p.set_immutable()
            sage: set(pts1) == set(pts2)
            True

            sage: timeit('Polyhedron(v).integral_points()')   # random output
            sage: timeit('LatticePolytope(v).points()')       # random output
        """
        if not self.is_compact():
            raise ValueError('Can only enumerate points in a compact polyhedron.')
        if self.n_vertices() == 0:
            return tuple()

        # for small bounding boxes, it is faster to naively iterate over the points of the box
        box_min, box_max = self.bounding_box(integral=True)
        box_points = prod(max_coord-min_coord+1 for min_coord, max_coord in zip(box_min, box_max))
        if  not self.is_lattice_polytope() or \
                (self.is_simplex() and box_points<1000) or \
                box_points<threshold:
            from sage.geometry.integral_points import rectangular_box_points
            return rectangular_box_points(box_min, box_max, self)

        # for more complicate polytopes, triangulate & use smith normal form
        from sage.geometry.integral_points import simplex_points
        if self.is_simplex():
            return simplex_points(self.Vrepresentation())
        triangulation = self.triangulate()
        points = set()
        for simplex in triangulation:
            triang_vertices = [ self.Vrepresentation(i) for i in simplex ]
            new_points = simplex_points(triang_vertices)
            for p in new_points:
                p.set_immutable()
            points.update(new_points)
        # assert all(self.contains(p) for p in points)   # slow
        return tuple(points)

    def combinatorial_automorphism_group(self):
        """
        Computes the combinatorial automorphism group of the vertex
        graph of the polyhedron.

        OUTPUT:

        A
        :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the combinatorial automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        EXAMPLES::

            sage: quadrangle = Polyhedron(vertices=[(0,0),(1,0),(0,1),(2,3)])
            sage: quadrangle.combinatorial_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)(3,4)]
            sage: quadrangle.restricted_automorphism_group()
            Permutation Group with generators [()]

        Permutations can only exchange vertices with vertices, rays
        with rays, and lines with lines::

            sage: P = Polyhedron(vertices=[(1,0,0), (1,1,0)], rays=[(1,0,0)], lines=[(0,0,1)])
            sage: P.combinatorial_automorphism_group()
            Permutation Group with generators [(3,4)]
        """
        if '_combinatorial_automorphism_group' in self.__dict__:
            return self._combinatorial_automorphism_group

        from sage.groups.perm_gps.permgroup import PermutationGroup

        G = Graph(sparse=True)
        for edge in self.vertex_graph().edges():
            i = edge[0]
            j = edge[1]
            G.add_edge(i, j, (self.Vrepresentation(i).type(), self.Vrepresentation(j).type()) )

        group, node_dict = G.automorphism_group(edge_labels=True, translation=True)

        # Relabel the permutation group
        perm_to_vertex = dict( (i,v+1) for v,i in node_dict.items() )
        group = PermutationGroup([ [ tuple([ perm_to_vertex[i] for i in cycle ])
                                     for cycle in generator.cycle_tuples() ]
                                   for generator in group.gens() ])

        self._combinatorial_automorphism_group = group
        return group


    def _affine_coordinates(self, Vrep_object):
        r"""
        Return affine coordinates for a V-representation object.

        INPUT:

        - ``v`` -- a V-representation object or any iterable
          containing ``self.ambient_dim()`` coordinates. The
          coordinates must specify a point in the affine plane
          containing the polyhedron, or the output will be invalid. No
          checks on the input are performed.

        OUTPUT:

        A ``self.dim()``-dimensional coordinate vector. It contains
        the coordinates of ``v`` in an arbitrary but fixed basis for
        the affine span of the polyhedron.

        EXAMPLES::

            sage: P = Polyhedron(rays=[(1,0,0),(0,1,0)])
            sage: P._affine_coordinates( (-1,-2,-3) )
            (-1, -2)
            sage: [ P._affine_coordinates(v) for v in P.Vrep_generator() ]
            [(1, 0), (0, 1), (0, 0)]
        """
        if '_affine_coordinates_pivots' not in self.__dict__:
            v_list = [ vector(v) for v in self.Vrepresentation() ]
            if len(v_list)>0:
                origin = v_list[0]
                v_list = [ v - origin for v in v_list ]
            coordinates = matrix(v_list)
            self._affine_coordinates_pivots = coordinates.pivots()

        v = list(Vrep_object)
        if len(v) != self.ambient_dim():
            raise ValueError('Incorrect dimension: '+str(v))

        return vector(self.field(), [ v[i] for i in self._affine_coordinates_pivots ])


    def restricted_automorphism_group(self):
        r"""
        Return the restricted automorphism group.

        First, let the linear automorphism group be the subgroup of
        the Euclidean group `E(d) = GL(d,\RR) \ltimes \RR^d`
        preserving the `d`-dimensional polyhedron. The Euclidean group
        acts in the usual way `\vec{x}\mapsto A\vec{x}+b` on the
        ambient space.

        The restricted automorphism group is the subgroup of the linear
        automorphism group generated by permutations of the generators
        of the same type. That is, vertices can only be permuted with
        vertices, ray generators with ray generators, and line
        generators with line generators.

        For example, take the first quadrant

        .. MATH::

            Q = \Big\{ (x,y) \Big| x\geq 0,\; y\geq0 \Big\}
            \subset \QQ^2

        Then the linear automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            a & 0 \\ 0 & b
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & c \\ d & 0
            \end{pmatrix}
            :~
            a, b, c, d \in \QQ_{>0}
            \right\}
            \subset
            GL(2,\QQ)
            \subset
            E(d)

        Note that there are no translations that map the quadrant `Q`
        to itself, so the linear automorphism group is contained in
        the subgroup of rotations of the whole Euclidean group. The
        restricted automorphism group is

        .. MATH::

            \mathrm{Aut}(Q) =
            \left\{
            \begin{pmatrix}
            1 & 0 \\ 0 & 1
            \end{pmatrix}
            ,~
            \begin{pmatrix}
            0 & 1 \\ 1 & 0
            \end{pmatrix}
            \right\}
            \simeq \ZZ_2

        OUTPUT:

        A :class:`PermutationGroup<sage.groups.perm_gps.permgroup.PermutationGroup_generic>`
        that is isomorphic to the restricted automorphism group is
        returned.

        Note that in Sage, permutation groups always act on positive
        integers while ``self.Vrepresentation()`` is indexed by
        nonnegative integers. The indexing of the permutation group is
        chosen to be shifted by ``+1``. That is, ``i`` in the
        permutation group corresponds to the V-representation object
        ``self.Vrepresentation(i-1)``.

        REFERENCES:

        ..  [BSS]
            David Bremner, Mathieu Dutour Sikiric, Achill Schuermann:
            Polyhedral representation conversion up to symmetries.
            http://arxiv.org/abs/math/0702239

        EXAMPLES::

            sage: P = polytopes.cross_polytope(3)
            sage: AutP = P.restricted_automorphism_group();  AutP
            Permutation Group with generators [(3,6), (2,3)(5,6), (2,5), (1,2)(4,5), (1,4)]
            sage: P24 = polytopes.twenty_four_cell()
            sage: AutP24 = P24.restricted_automorphism_group()
            sage: PermutationGroup([   # computed with sympow
            ...    [(1,2),(3,4),(5,6),(7,8),(9,10),(11,12),(13,14),(15,16),(17,21)],
            ...    [(1,5),(2,6),(3,7),(4,8),(9,13),(10,14),(11,15),(12,16),(19,23)],
            ...    [(2,5),(4,7),(10,13),(12,15),(17,19),(21,23)],
            ...    [(1,3),(2,4),(5,7),(6,8),(9,11),(10,12),(13,15),(14,16),(18,22)],
            ...    [(3,5),(4,6),(11,13),(12,14),(18,19),(22,23)],
            ...    [(1,9),(2,10),(3,11),(4,12),(5,13),(6,14),(7,15),(8,16),(20,24)],
            ...    [(3,9),(4,10),(7,13),(8,14),(18,20),(22,24)],
            ...    [(1,3,11,9),(2,4,12,10),(5,7,15,13),(6,8,16,14),(18,20,22,24)],
            ...    [(3,9,5),(4,10,6),(7,11,13),(8,12,14),(18,20,19),(22,24,23)],
            ...    [(1,5,7,15,11,9),(2,6,8,16,12,10),(3,13),(4,14),(18,20,23,22,24,19)],
            ...    [(2,5,3,9),(4,13),(6,7,11,10),(8,15,12,14),(17,19,18,20),(21,23,22,24)],
            ...    [(1,2,6,8,16,15,11,9),(3,10,5,4,14,7,12,13),(17,19,18,20,21,23,22,24)],
            ...    [(1,12,20,16,5,24),(2,21,6,14,18,10),(3,22,7,15,17,11),(4,8,23,13,9,19)]
            ...   ]) == AutP24
            True

        Here is the quadrant example mentioned in the beginning::

            sage: P = Polyhedron(rays=[(1,0),(0,1)])
            sage: P.Vrepresentation()
            [A ray in the direction (1, 0), A ray in the direction (0, 1), A vertex at (0, 0)]
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]

        Also, the polyhedron need not be full-dimensional::

            sage: P = Polyhedron(vertices=[(1,2,3,4,5),(7,8,9,10,11)])
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(1,2)]

        Translations do not change the restricted automorphism
        group. For example, any non-degenerate triangle has the
        dihedral group with 6 elements, `D_6`, as its automorphism
        group::

            sage: initial_points = [vector([1,0]), vector([0,1]), vector([-2,-1])]
            sage: points = initial_points
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[0] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]
            sage: points = [pt - 2*initial_points[1] for pt in initial_points]
            sage: Polyhedron(vertices=points).restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        Floating-point computations are supported with a simple fuzzy
        zero implementation::

            sage: P = Polyhedron(vertices=[(1.0/3.0,0,0),(0,1.0/3.0,0),(0,0,1.0/3.0)], field=RDF)
            sage: P.restricted_automorphism_group()
            Permutation Group with generators [(2,3), (1,2)]

        TESTS::

            sage: p = Polyhedron(vertices=[(1,0), (1,1)], rays=[(1,0)])
            sage: p.restricted_automorphism_group()
            Permutation Group with generators [(2,3)]
        """
        if '_restricted_automorphism_group' in self.__dict__:
            return self._restricted_automorphism_group

        from sage.groups.perm_gps.permgroup import PermutationGroup

        if self.field() is QQ:
            def rational_approximation(c):
                return c

        else:  # self.field() is RDF
            c_list = []
            def rational_approximation(c):
                # Implementation detail: Return unique integer if two
                # c-values are the same up to machine precision. But
                # you can think of it as a uniquely-chosen rational
                # approximation.
                for i,x in enumerate(c_list):
                    if self._is_zero(x-c):
                        return i
                c_list.append(c)
                return len(c_list)-1

        # The algorithm identifies the restricted automorphism group
        # with the automorphism group of a edge-colored graph. The
        # nodes of the graph are the V-representation objects. If all
        # V-representation objects are vertices, the edges are
        # labelled by numbers (to be computed below). Roughly
        # speaking, the edge label is the inner product of the
        # coordinate vectors with some orthogonalization thrown in
        # [BSS].
        def edge_label_compact(i,j,c_ij):
            return c_ij

        # In the non-compact case we also label the edges by the type
        # of the V-representation object. This ensures that vertices,
        # rays, and lines are only permuted amongst themselves.
        def edge_label_noncompact(i,j,c_ij):
            return (self.Vrepresentation(i).type(), c_ij, self.Vrepresentation(j).type())

        if self.is_compact():
            edge_label = edge_label_compact
        else:
            edge_label = edge_label_noncompact

        # good coordinates for the V-representation objects
        v_list = []
        for v in self.Vrepresentation():
            v_coords = list(self._affine_coordinates(v))
            if v.is_vertex():
                v_coords = [1]+v_coords
            else:
                v_coords = [0]+v_coords
            v_list.append(vector(v_coords))

        # Finally, construct the graph
        Qinv = sum( v.column() * v.row() for v in v_list ).inverse()

        # Was set to sparse = False, but there is a problem with Graph
        # backends. It should probably be set back to sparse = False as soon as
        # the backends are fixed.
        G = Graph(sparse=True)
        for i in range(0,len(v_list)):
            for j in range(i+1,len(v_list)):
                v_i = v_list[i]
                v_j = v_list[j]
                c_ij = rational_approximation( v_i * Qinv * v_j )
                G.add_edge(i,j, edge_label(i,j,c_ij))

        group, node_dict = G.automorphism_group(edge_labels=True, translation=True)

        # Relabel the permutation group
        perm_to_vertex = dict( (i,v+1) for v,i in node_dict.items() )
        group = PermutationGroup([ [ tuple([ perm_to_vertex[i] for i in cycle ])
                                     for cycle in generator.cycle_tuples() ]
                                   for generator in group.gens() ])

        self._restricted_automorphism_group = group
        return group


#############################################################
def cyclic_sort_vertices_2d(Vlist):
    """
    Return the vertices/rays in cyclic order if possible.

    NOTES:

    This works if and only if each vertex/ray is adjacent to exactly
    two others. For example, any 2-dimensional polyhedron satisfies
    this.

    EXAMPLES::

        sage: from sage.geometry.polyhedra import cyclic_sort_vertices_2d
        sage: square = Polyhedron([[1,0],[-1,0],[0,1],[0,-1]])
        sage: vertices = [v for v in square.vertex_generator()]
        sage: vertices
        [A vertex at (1, 0),
         A vertex at (-1, 0),
         A vertex at (0, 1),
         A vertex at (0, -1)]
        sage: cyclic_sort_vertices_2d(vertices)
        [A vertex at (0, -1),
         A vertex at (1, 0),
         A vertex at (0, 1),
         A vertex at (-1, 0)]
    """
    if len(Vlist)==0: return Vlist

    adjacency_matrix = Vlist[0].polyhedron().vertex_adjacency_matrix()
    result = [ Vlist.pop() ]
    while len(Vlist)>0:
        for i in range(len(Vlist)):
            if adjacency_matrix[Vlist[i].index(), result[-1].index()] == 1:
                result.append( Vlist.pop(i) )
                break;
        else:
            raise ValueError
    return result




#########################################################################
def projection_func_identity(x):
    """
    The identity projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedra import projection_func_identity
        sage: projection_func_identity((1,2,3))
        [1, 2, 3]
    """
    return list(x)



class ProjectionFuncStereographic():
    """
    The stereographic (or perspective) projection.

    EXAMPLES::

        sage: from sage.geometry.polyhedra import ProjectionFuncStereographic
        sage: cube = polytopes.n_cube(3).vertices()
        sage: proj = ProjectionFuncStereographic([1.1,1.1,1.1])
        sage: ppoints = [proj(vector(x)) for x in cube]
        sage: ppoints[1]
        (-0.3182829598..., 1.18784817...)
    """
    def __init__(self, projection_point):
        """
        Create a stereographic projection function.

        INPUT:

        - ``projection_point`` -- a list of coordinates in the
          appropriate dimension, which is the point projected from.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__init__([1.0,1.0])
            sage: proj.house
            [-0.7071067811...  0.7071067811...]
            [ 0.7071067811...  0.7071067811...]
        """
        self.projection_point = vector(projection_point)
        self.dim = self.projection_point.degree()

        pproj = vector(RDF,self.projection_point)
        self.psize = norm(pproj)
        if (self.psize).is_zero():
            raise ValueError, "projection direction must be a non-zero vector."
        v = vector(RDF, [0.0]*(self.dim-1) + [self.psize]) - pproj
        polediff = matrix(RDF,v).transpose()
        denom = RDF((polediff.transpose()*polediff)[0][0])
        if denom.is_zero():
            self.house = identity_matrix(RDF,self.dim)
        else:
            self.house = identity_matrix(RDF,self.dim) \
            - 2*polediff*polediff.transpose()/denom   # Householder reflector

    def __call__(self, x):
        """
        Action of the stereographic projection.

        INPUT:

        - ``x`` -- a vector or anything convertible to a vector.

        OUTPUT:

        First reflects ``x`` with a Householder reflection which takes
        the projection point to ``(0,...,0,self.psize)`` where
        ``psize`` is the length of the projection point, and then
        dilates by ``1/(zdiff)`` where ``zdiff`` is the difference
        between the last coordinate of ``x`` and ``psize``.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import ProjectionFuncStereographic
            sage: proj = ProjectionFuncStereographic([1.0,1.0])
            sage: proj.__call__(vector([1,2]))
            (-1.0)
            sage: proj = ProjectionFuncStereographic([2.0,1.0])
            sage: proj.__call__(vector([1,2]))
            (3.0)
            sage: proj = ProjectionFuncStereographic([0,0,2])
            sage: proj.__call__(vector([0,0,1]))
            (0.0, 0.0)
            sage: proj.__call__(vector([1,0,0]))
            (0.5, 0.0)
        """
        img = self.house * x
        denom = self.psize-img[self.dim-1]
        if denom.is_zero():
            raise ValueError, 'Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x)
        return vector(RDF, [img[i]/denom for i in range(self.dim-1)])


class ProjectionFuncSchlegel():
    """
    The Schlegel projection from the given input point.

    EXAMPLES::

        sage: from sage.geometry.polyhedra import ProjectionFuncSchlegel
        sage: proj = ProjectionFuncSchlegel([2,2,2])
        sage: proj(vector([1.1,1.1,1.11]))[0]
        0.0302...
    """
    def __init__(self, projection_direction, height = 1.1):
        """
        Initializes the projection.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__init__([2,2,2])
            sage: proj(vector([1.1,1.1,1.11]))[0]
            0.0302...
        """
        self.projection_dir = vector(RDF, projection_direction)
        if norm(self.projection_dir).is_zero():
            raise ValueError, "projection direction must be a non-zero vector."
        self.dim = self.projection_dir.degree()
        spcenter = height * self.projection_dir/norm(self.projection_dir)
        self.height = height
        v = vector(RDF, [0.0]*(self.dim-1) + [self.height]) - spcenter
        polediff = matrix(RDF,v).transpose()
        denom = (polediff.transpose()*polediff)[0][0]
        if denom.is_zero():
            self.house = identity_matrix(RDF,self.dim)
        else:
            self.house = identity_matrix(RDF,self.dim) \
            - 2*polediff*polediff.transpose()/denom #Householder reflector

    def __call__(self, x):
        """
        Apply the projection to a vector.

        - ``x`` -- a vector or anything convertible to a vector.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import ProjectionFuncSchlegel
            sage: proj = ProjectionFuncSchlegel([2,2,2])
            sage: proj.__call__([1,2,3])
            (0.56162854..., 2.09602626...)
        """
        v = vector(RDF,x)
        if v.is_zero():
            raise ValueError, "The origin must not be a vertex."
        v = v/norm(v)         # normalize vertices to unit sphere
        v = self.house*v      # reflect so self.projection_dir is at "north pole"
        denom = self.height-v[self.dim-1]
        if denom.is_zero():
            raise ValueError, 'Point cannot coincide with ' \
                'coordinate singularity at ' + repr(x)
        return vector(RDF, [ v[i]/denom for i in range(self.dim-1) ])



#########################################################################
class Projection(SageObject):
    """
    The projection of a :class:`Polyhedron`.

    This class keeps track of the necessary data to plot the input
    polyhedron.
    """

    def __init__(self, polyhedron, proj=projection_func_identity):
        """
        Initialize the projection of a Polyhedron() object.

        INPUT:

          - ``polyhedron`` - a ``Polyhedron()`` object
          - ``proj`` - a projection function for the points

        NOTES:

        Once initialized, the polyhedral data is fixed. However, the
        projection can be changed later on.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedra import Projection
            sage: Projection(p)
            The projection of a polyhedron into 3 dimensions.
            sage: def pr_12(x): return [x[1],x[2]]
            sage: Projection(p, pr_12)
            The projection of a polyhedron into 2 dimensions.
            sage: Projection(p,  lambda x: [x[1],x[2]] )   # another way of doing the same projection
            The projection of a polyhedron into 2 dimensions.
            sage: _.show()   # plot of the projected icosahedron in 2d
            sage: proj = Projection(p)
            sage: proj.stereographic([1,2,3])
            The projection of a polyhedron into 2 dimensions.
            sage: proj.show()
        """
        self.coords = Sequence([])
        self.points = Sequence([])
        self.lines  = Sequence([])
        self.arrows = Sequence([])
        self.polygons = Sequence([])
        self.polyhedron_ambient_dim = polyhedron.ambient_dim()

        if polyhedron.ambient_dim() == 2:
            self._init_from_2d(polyhedron)
        elif polyhedron.ambient_dim() == 3:
            self._init_from_3d(polyhedron)
        else:
            self._init_points(polyhedron)
            self._init_lines_arrows(polyhedron)

        self.coords.set_immutable()
        self.points.set_immutable()
        self.lines.set_immutable()
        self.arrows.set_immutable()
        self.polygons.set_immutable()

        self.__call__(proj)


    def _repr_(self):
        """
        Return a string describing the projection.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: from sage.geometry.polyhedra import Projection
            sage: proj = Projection(p)
            sage: print proj._repr_()
            The projection of a polyhedron into 3 dimensions.
        """
        s = 'The projection of a polyhedron into ' + repr(self.dimension)
        s += ' dimensions.'
        return s + "\n"


    def __call__(self, proj=projection_func_identity):
        """
        Apply a projection.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedra import Projection
            sage: pproj = Projection(p)
            sage: from sage.geometry.polyhedra import ProjectionFuncStereographic
            sage: pproj_stereo = pproj.__call__(proj = ProjectionFuncStereographic([1,2,3]))
            sage: pproj_stereo.polygons[0]
            [8, 1, 3]
        """
        self.transformed_coords = \
            Sequence([proj(p) for p in self.coords])
        self._init_dimension()
        return self


    def identity(self):
        """
        Return the identity projection of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedra import Projection
            sage: pproj = Projection(p)
            sage: ppid = pproj.identity()
            sage: ppid.dimension
            3
        """
        return self.__call__(projection_func_identity)


    def stereographic(self, projection_point=None):
        r"""
        Rteurn the stereographic projection.

        INPUT:

        - ``projection_point`` - The projection point. This must be
          distinct from the polyhedron's vertices. Default is `(1,0,\dots,0)`

        EXAMPLES::

            sage: from sage.geometry.polyhedra import Projection
            sage: proj = Projection(polytopes.buckyball())  #long time
            sage: proj                                      #long time
            The projection of a polyhedron into 3 dimensions.
            sage: proj.stereographic([5,2,3]).show()        #long time
            sage: Projection( polytopes.twenty_four_cell() ).stereographic([2,0,0,0])
            The projection of a polyhedron into 3 dimensions.

        """
        if projection_point == None:
            projection_point = [1] + [0]*(self.polyhedron_ambient_dim-1)
        return self.__call__(ProjectionFuncStereographic(projection_point))


    def schlegel(self, projection_direction=None, height = 1.1):
        """
        Return the Schlegel projection.

        The vertices are normalized to the unit sphere, and
        stereographically projected from a point slightly outside of
        the sphere.

        INPUT:

          - ``projection_direction`` - The direction of the Schlegel
            projection. By default, the vector consisting of the first
            n primes is chosen.

        EXAMPLES::

            sage: cube4 = polytopes.n_cube(4)
            sage: from sage.geometry.polyhedra import Projection
            sage: Projection(cube4).schlegel([1,0,0,0])
            The projection of a polyhedron into 3 dimensions.
            sage: _.show()
        """
        if projection_direction == None:
            for poly in self.polygons:
                center = sum([self.coords[i] for i in poly]) / len(poly)
                print center, "\n"
                if not center.is_zero():
                    projection_direction = center
                    break
        if projection_direction == None:
            projection_direction = primes_first_n(self.polyhedron_ambient_dim)
        return self.__call__(ProjectionFuncSchlegel(projection_direction, height = height))


    def coord_index_of(self, v):
        """
        Convert a coordinate vector to its internal index.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj.coord_index_of(vector((1,1,1)))
            0
        """
        try:
            return self.coords.index(v)
        except ValueError:
            self.coords.append(v)
            return len(self.coords)-1


    def coord_indices_of(self, v_list):
        """
        Convert list of coordinate vectors to the corresponding list
        of internal indices.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: proj = p.projection()
            sage: proj.coord_indices_of([vector((1,1,1)),vector((1,-1,1))])
            [0, 2]
        """
        return [self.coord_index_of(v) for v in v_list]


    def coordinates_of(self, coord_index_list):
        """
        Given a list of indices, return the projected coordinates.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4).projection()
            sage: p.coordinates_of([1])
            [[0, 0, -43301/50000, 22361/100000]]
        """
        return [self.transformed_coords[i] for i in coord_index_list]


    def _init_dimension(self):
        """
        Internal function: Initialize from 2d polyhedron. Must always
        be called after a coordinate projection.

        TESTS::

            sage: from sage.geometry.polyhedra import Projection, render_2d
            sage: p = polytopes.n_simplex(2).projection()
            sage: test = p._init_dimension()
            sage: p.show.__doc__ == render_2d.__doc__
            True
        """
        self.dimension = len(self.transformed_coords[0])

        if self.dimension == 2:
            self.show = lambda **kwds: render_2d(self,**kwds)
            self.show.__doc__ = render_2d.__doc__
        elif self.dimension == 3:
            self.show = lambda **kwds: render_3d(self,**kwds)
            self.show.__doc__ = render_3d.__doc__
        else:
            try:
                del self.show
            except AttributeError:
                pass


    def _init_from_2d(self, polyhedron):
        """
        Internal function: Initialize from 2d polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0],[0,1],[1,0],[1,1]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0]], [[0, 1]], [[1, 0]], [[1, 1]]]
            sage: proj._init_from_2d
            <bound method Projection._init_from_2d of The projection of a polyhedron into 2 dimensions.
            >
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        self.dimension = 2
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_area_2d(polyhedron)


    def _init_from_3d(self, polyhedron):
        """
        Internal function: Initialize from 3d polyhedron.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,1],[0,1,2],[1,0,3],[1,1,5]])
            sage: proj = p.projection()
            sage: [proj.coordinates_of([i]) for i in proj.points]
            [[[0, 0, 1]], [[0, 1, 2]], [[1, 0, 3]], [[1, 1, 5]]]
            sage: proj._init_from_3d
            <bound method Projection._init_from_3d of The projection of a polyhedron into 3 dimensions.
            >
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"
        self.dimension = 3
        self._init_points(polyhedron)
        self._init_lines_arrows(polyhedron)
        self._init_solid_3d(polyhedron)


    def _init_points(self, polyhedron):
        """
        Internal function: Initialize points (works in arbitrary
        dimensions).

        TESTS::

            sage: p = polytopes.n_cube(2)
            sage: pp = p.projection()
            sage: del pp.points
            sage: pp.points = Sequence([])
            sage: pp._init_points(p)
            sage: pp.points
            [0, 1, 2, 3]
        """
        for v in polyhedron.vertex_generator():
            self.points.append( self.coord_index_of(v.vector()) )


    def _init_lines_arrows(self, polyhedron):
        """
        Internal function: Initialize compact and non-compact edges
        (works in arbitrary dimensions).

        TESTS::

            sage: p = Polyhedron(ieqs = [[1, 0, 0, 1],[1,1,0,0]])
            sage: pp = p.projection()
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]
            sage: del pp.arrows
            sage: pp.arrows = Sequence([])
            sage: pp._init_lines_arrows(p)
            sage: pp.arrows
            [[0, 1], [0, 2], [0, 3], [0, 4]]
        """
        obj = polyhedron.Vrepresentation()
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(len(obj)):
                if polyhedron.vertex_adjacency_matrix()[i,j] == 0: continue
                if i<j and obj[j].is_vertex():
                    l = [obj[i].vector(), obj[j].vector()]
                    self.lines.append( [ self.coord_index_of(l[0]),
                                         self.coord_index_of(l[1]) ] )
                if obj[j].is_ray():
                    l = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l[0]),
                                          self.coord_index_of(l[1]) ] )
                if obj[j].is_line():
                    l1 = [obj[i].vector(), obj[i].vector() + obj[j].vector()]
                    l2 = [obj[i].vector(), obj[i].vector() - obj[j].vector()]
                    self.arrows.append( [ self.coord_index_of(l1[0]),
                                          self.coord_index_of(l1[1]) ] )
                    self.arrows.append( [ self.coord_index_of(l2[0]),
                                          self.coord_index_of(l2[1]) ] )


    def _init_area_2d(self, polyhedron):
        """
        Internal function: Initialize polygon area for 2d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(2,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_area_2d(p)
            sage: proj.polygons
            [[3, 0, 1, 2]]
        """
        assert polyhedron.ambient_dim() == 2, "Requires polyhedron in 2d"
        vertices = [v for v in polyhedron.Vrep_generator()]
        vertices = cyclic_sort_vertices_2d(vertices)
        coords = []

        def adjacent_vertices(i):
            n = len(vertices)
            if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
            if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

        for i in range(len(vertices)):
            v = vertices[i]
            if v.is_vertex():
                coords.append(v())
            if v.is_ray():
                for a in adjacent_vertices(i):
                    coords.append(a() + v())

        if polyhedron.n_lines() == 0:
            self.polygons.append( self.coord_indices_of(coords) )
            return

        polygons = []

        if polyhedron.n_lines() == 1:
            aline = polyhedron.line_generator().next()
            for shift in [aline(), -aline()]:
                for i in range(len(coords)):
                    polygons.append( [ coords[i-1],coords[i],
                                       coords[i]+shift, coords[i-1]+shift ] )

        if polyhedron.n_lines() == 2:
            [line1, line2] = [l for l in polyhedron.lines()]
            assert len(coords)==1, "Can have only a single vertex!"
            v = coords[0]
            l1 = line1()
            l2 = line2()
            polygons = [ [v-l1-l2, v+l1-l2, v+l1+l2, v-l1+l2] ]

        polygons = [ self.coord_indices_of(p) for p in polygons ]
        self.polygons.extend(polygons)



    def _init_solid_3d(self, polyhedron):
        """
        Internal function: Initialize facet polygons for 3d polyhedron.

        TESTS::

            sage: p = polytopes.cyclic_polytope(3,4)
            sage: proj = p.projection()
            sage: proj.polygons = Sequence([])
            sage: proj._init_solid_3d(p)
            sage: proj.polygons
            [[3, 1, 2], [3, 0, 2], [3, 0, 1], [2, 0, 1]]
        """
        assert polyhedron.ambient_dim() == 3, "Requires polyhedron in 3d"

        if polyhedron.dim() <= 1: # empty or 0d or 1d polyhedron => no polygon
            return None

        def defining_equation():  # corresponding to a polygon
            if polyhedron.dim() < 3:
                yield polyhedron.equation_generator().next()
            else:
                for ineq in polyhedron.inequality_generator():
                    yield ineq

        faces = []
        for facet_equation in defining_equation():
            vertices = [v for v in facet_equation.incident()]
            vertices = cyclic_sort_vertices_2d(vertices)
            coords = []

            def adjacent_vertices(i):
                n = len(vertices)
                if vertices[(i-1) % n].is_vertex(): yield vertices[(i-1) % n]
                if vertices[(i+1) % n].is_vertex(): yield vertices[(i+1) % n]

            for i in range(len(vertices)):
                v = vertices[i]
                if v.is_vertex():
                    coords.append(v())
                if v.is_ray():
                    for a in adjacent_vertices(i):
                        coords.append(a() + v())

            faces.append(coords)

        if polyhedron.n_lines() == 0:
            assert len(faces)>0, "no vertices?"
            self.polygons.extend( [self.coord_indices_of(f) for f in faces] )
            return

        # now some special cases if there are lines (dim < ambient_dim)
        polygons = []

        if polyhedron.n_lines()==1:
            assert len(faces)>0, "no vertices?"
            aline = polyhedron.line_generator().next()
            for shift in [aline(), -aline()]:
                for coords in faces:
                    assert len(coords)==2, "There must be two points."
                    polygons.append( [ coords[0],coords[1],
                                       coords[1]+shift, coords[0]+shift ] )

        if polyhedron.n_lines()==2:
            [line1, line2] = [l for l in polyhedron.line_generator()]
            l1 = line1()
            l2 = line2()
            for v in polyhedron.vertex_generator():
                polygons.append( [v()-l1-l2, v()+l1-l2, v()+l1+l2, v()-l1+l2] )

        self.polygons.extend( [self.coord_indices_of(p) for p in polygons] )


    def render_points_2d(self, **kwds):
        """
        Return the points of a polyhedron in 2d.

        EXAMPLES::

            sage: hex = polytopes.regular_polygon(6)
            sage: proj = hex.projection()
            sage: hex_points = proj.render_points_2d()
            sage: hex_points._Graphics__objects
            [Point set defined by 6 point(s)]
        """
        return point2d(self.coordinates_of(self.points), **kwds)


    def render_outline_2d(self, **kwds):
        """
        Return the outline (edges) of a polyhedron in 2d.

        EXAMPLES::

            sage: penta = polytopes.regular_polygon(5)
            sage: outline = penta.projection().render_outline_2d()
            sage: outline._Graphics__objects[0]
            Line defined by 2 points
        """
        wireframe = [];
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line2d(l_coords, **kwds) )
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append( arrow(a_coords[0], a_coords[1], **kwds) )
        return sum(wireframe)


    def render_fill_2d(self, **kwds):
        """
        Return the filled interior (a polygon) of a polyhedron in 2d.

        EXAMPLES::

            sage: cps = [i^3 for i in srange(-2,2,1/5)]
            sage: p = Polyhedron(vertices = [[(t^2-1)/(t^2+1),2*t/(t^2+1)] for t in cps])
            sage: proj = p.projection()
            sage: filled_poly = proj.render_fill_2d()
            sage: filled_poly.axes_width()
            0.8
        """
        poly = [polygon2d(self.coordinates_of(p), **kwds)
                 for p in self.polygons]
        return sum(poly)


    def render_vertices_3d(self, **kwds):
        """
        Return the 3d rendering of the vertices.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(3)
            sage: proj = p.projection()
            sage: verts = proj.render_vertices_3d()
            sage: verts.bounding_box()
            ((-1.0, -1.0, -1.0), (1.0, 1.0, 1.0))
        """
        return point3d(self.coordinates_of(self.points), **kwds)


    def render_wireframe_3d(self, **kwds):
        r"""
        Return the 3d wireframe rendering.

        EXAMPLES::

            sage: cube = polytopes.n_cube(3)
            sage: cube_proj = cube.projection()
            sage: wire = cube_proj.render_wireframe_3d()
            sage: print wire.tachyon().split('\n')[77]  # for testing
            FCylinder base 1.0 -1.0 1.0 apex 1.0 1.0 1.0 rad 0.005 texture...
        """
        wireframe = [];
        for l in self.lines:
            l_coords = self.coordinates_of(l)
            wireframe.append( line3d(l_coords, **kwds))
        for a in self.arrows:
            a_coords = self.coordinates_of(a)
            wireframe.append(arrow3d(a_coords[0], a_coords[1], **kwds))
        return sum(wireframe)


    def render_solid_3d(self, **kwds):
        """
        Return solid 3d rendering of a 3d polytope.

        EXAMPLES::

            sage: p = polytopes.n_cube(3).projection()
            sage: p_solid = p.render_solid_3d(opacity = .7)
            sage: type(p_solid)
            <class 'sage.plot.plot3d.base.Graphics3dGroup'>
        """
        return sum([ polygon3d(self.coordinates_of(f), **kwds)
                     for f in self.polygons ])


#########################################################################
class Polytopes():
    """
    A class of constructors for commonly used, famous, or interesting
    polytopes.
    """

    @staticmethod
    def orthonormal_1(dim_n=5):
        """
        A matrix of rational approximations to orthonormal vectors to
        ``(1,...,1)``.

        INPUT:

        - ``dim_n`` - the dimension of the vectors

        OUTPUT:

        A matrix over ``QQ`` whose rows are close to an orthonormal
        basis to the subspace normal to ``(1,...,1)``.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import Polytopes
            sage: m = Polytopes.orthonormal_1(5)
            sage: m
            [ 70711/100000   -7071/10000             0             0             0]
            [    1633/4000     1633/4000 -81649/100000             0             0]
            [   7217/25000    7217/25000    7217/25000  -43301/50000             0]
            [ 22361/100000  22361/100000  22361/100000  22361/100000  -44721/50000]
        """
        pb = []
        for i in range(0,dim_n-1):
            pb.append([1.0/(i+1)]*(i+1) + [-1] + [0]*(dim_n-i-2))
        m = matrix(RDF,pb)
        new_m = []
        for i in range(0,dim_n-1):
            new_m.append([RDF(100000*q/norm(m[i])).ceil()/100000 for q in m[i]])
        return matrix(QQ,new_m)

    @staticmethod
    def project_1(fpoint):
        """
        Take a ndim-dimensional point and projects it onto the plane
        perpendicular to (1,1,...,1).

        INPUT:

          - ``fpoint`` - a list of ndim numbers

        EXAMPLES::

            sage: from sage.geometry.polyhedra import Polytopes
            sage: Polytopes.project_1([1,1,1,1,2])
            [1/100000, 1/100000, 1/50000, -559/625]
        """
        dim_n = len(fpoint)
        p_basis = [list(q) for q in Polytopes.orthonormal_1(dim_n)]
        out_v = []
        for v in p_basis:
            out_v.append(sum([fpoint[ind]*v[ind] for ind in range(dim_n)]))
        return out_v

    @staticmethod
    def _pfunc(i,j,perm):
        """
        An internal utility function for constructing the Birkhoff polytopes.

        EXAMPLES::

            sage: from sage.geometry.polyhedra import Polytopes
            sage: Polytopes._pfunc(1,2,permutations(3)[0])
            0
        """
        if perm[i-1] == j:
            return 1
        else:
            return 0


    def regular_polygon(self, n, field = QQ):
        """
        Return a regular polygon with n vertices.  Over the rational
        field the vertices may not be exact.

        INPUT:

        - ``n`` -- a positive integer, the number of vertices.

        - ``field`` -- either ``QQ`` or ``RDF``.

        EXAMPLES::

            sage: octagon = polytopes.regular_polygon(8)
            sage: len(octagon.vertices())
            8
        """
        npi = 3.14159265359
        verts = []
        for i in range(n):
            t = 2*npi*i/n
            verts.append([sin(t),cos(t)])
        verts = [[field(RDF(x)) for x in y] for y in verts]
        return Polyhedron(vertices = verts, field = field)


    def Birkhoff_polytope(self, n):
        """
        Return the Birkhoff polytope with n! vertices.  Each vertex
        is a (flattened) n by n permutation matrix.

        INPUT:

        - ``n`` -- a positive integer giving the size of the permutation matrices.

        EXAMPLES::

            sage: b3 = polytopes.Birkhoff_polytope(3)
            sage: b3.n_vertices()
            6
        """
        perms = permutations(range(1,n+1))
        verts = []
        for p in perms:
            verts += [ [Polytopes._pfunc(i,j,p) for j in range(1,n+1)
                        for i in range(1,n+1) ] ]
        return Polyhedron(vertices = verts)


    def n_simplex(self, dim_n=3, project = True):
        """
        Return a rational approximation to a regular simplex in
        dimension ``dim_n``.

        INPUT:

        - ``dim_n`` -- The dimension of the cross-polytope, a positive
          integer.

        - ``project`` -- Optional argument, whether to project
          orthogonally.  Default is True.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional simplex.

        EXAMPLES::

            sage: s5 = polytopes.n_simplex(5)
            sage: s5.dim()
            5
        """
        verts = permutations([0 for i in range(dim_n)] + [1])
        if project: verts = [Polytopes.project_1(x) for x in verts]
        return Polyhedron(vertices = verts)


    def icosahedron(self, field = QQ):
        """
        Return an icosahedron with edge length 1.

        INPUT:

        - ``field`` -- Either ``QQ`` or ``RDF``.

        OUTPUT:

        A Polyhedron object of a floating point or rational
        approximation to the regular 3d icosahedron.

        If ``field=QQ``, a rational approximation is used and the
        points are not exactly the vertices of the icosahedron. The
        icosahedron's coordinates contain the golden ratio, so there
        is no exact representation possible.

        EXAMPLES::

            sage: ico = polytopes.icosahedron()
            sage: sum(sum( ico.vertex_adjacency_matrix() ))/2
            30
        """
        if field == QQ:
            g = QQ(1618033)/1000000 # Golden ratio approximation
            r12 = QQ(1)/2
        elif field == RDF:
            g = RDF( (1 + sqrt(5))/2 )
            r12 = RDF( QQ(1)/2 )
        else:
            raise ValueError, "field must be QQ or RDF."
        verts = [i([0,r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,r12,-g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,g/2]) for i in AlternatingGroup(3)]
        verts = verts + [i([0,-r12,-g/2]) for i in AlternatingGroup(3)]
        return Polyhedron(vertices = verts, field = field)


    def dodecahedron(self, field=QQ):
        """
        Return a dodecahedron.

        INPUT:

        - ``field`` -- Either ``QQ`` (in which case a rational
          approximation to the golden ratio is used) or ``RDF``.

        EXAMPLES::

            sage: d12 = polytopes.dodecahedron()
            sage: d12.n_inequalities()
            12
        """
        return self.icosahedron(field=field).polar()


    def small_rhombicuboctahedron(self):
        """
        Return an Archimedean solid with 24 vertices and 26 faces.

        EXAMPLES::

            sage: sr = polytopes.small_rhombicuboctahedron()
            sage: sr.n_vertices()
            24
            sage: sr.n_inequalities()
            26
        """
        verts = [ [-3/2, -1/2, -1/2], [-3/2, -1/2, 1/2], [-3/2, 1/2, -1/2],
                  [-3/2, 1/2, 1/2], [-1/2, -3/2, -1/2], [-1/2, -3/2, 1/2],
                  [-1/2, -1/2, -3/2], [-1/2,-1/2, 3/2], [-1/2, 1/2, -3/2],
                  [-1/2, 1/2, 3/2], [-1/2, 3/2, -1/2], [-1/2, 3/2, 1/2],
                  [1/2, -3/2, -1/2], [1/2, -3/2, 1/2], [1/2, -1/2,-3/2],
                  [1/2, -1/2, 3/2], [1/2, 1/2, -3/2], [1/2, 1/2, 3/2],
                  [1/2, 3/2,-1/2], [1/2, 3/2, 1/2], [3/2, -1/2, -1/2],
                  [3/2, -1/2, 1/2], [3/2, 1/2,-1/2], [3/2, 1/2, 1/2] ]
        return Polyhedron(vertices = verts)


    def great_rhombicuboctahedron(self, field=QQ):
        """
        Return an Archimedean solid with 48 vertices and 26 faces.

        EXAMPLES::

            sage: gr = polytopes.great_rhombicuboctahedron()
            sage: gr.n_vertices()
            48
            sage: gr.n_inequalities()
            26
        """
        v1 = QQ(131739771357/54568400000)
        v2 = QQ(104455571357/27284200000)
        verts = [ [1, v1, v2],
                  [1, v2, v1],
                  [v1, 1, v2],
                  [v1, v2, 1],
                  [v2, 1, v1],
                  [v2, v1, 1] ]
        verts = verts + [[x[0],x[1],-x[2]] for x in verts]
        verts = verts + [[x[0],-x[1],x[2]] for x in verts]
        verts = verts + [[-x[0],x[1],x[2]] for x in verts]
        if field!=QQ:
            verts = [field(v) for v in verts]
        return Polyhedron(vertices=verts, field=field)


    def rhombic_dodecahedron(self):
        """
        This face-regular, vertex-uniform polytope is dual to the
        cuboctahedron. It has 14 vertices and 12 faces.

        EXAMPLES::

            sage: rd = polytopes.rhombic_dodecahedron()
            sage: rd.n_vertices()
            14
            sage: rd.n_inequalities()
            12
        """
        v = [ [1, 1, 1], [1, 1, -1], [1, -1, 1], [1, -1, -1], [-1, 1, 1],
              [-1, 1, -1], [-1, -1, 1], [-1, -1, -1], [0, 0, 2], [0, 2, 0],
              [2, 0, 0], [0, 0, -2], [0, -2, 0], [-2, 0, 0] ]
        return Polyhedron(vertices = v)


    def cuboctahedron(self):
        """
        An Archimedean solid with 12 vertices and 14 faces.  Dual to
        the rhombic dodecahedron.

        EXAMPLES::

            sage: co = polytopes.cuboctahedron()
            sage: co.n_vertices()
            12
            sage: co.n_inequalities()
            14
        """
        one = Integer(1)
        v = [ [0, -one/2, -one/2], [0, one/2, -one/2], [one/2, -one/2, 0],
              [one/2, one/2, 0], [one/2, 0, one/2], [one/2, 0, -one/2],
              [0, one/2, one/2], [0, -one/2, one/2], [-one/2, 0, one/2],
              [-one/2, one/2, 0], [-one/2, 0, -one/2], [-one/2, -one/2, 0] ]
        return Polyhedron(vertices = v)


    def buckyball(self, field = QQ):
        """
        Also known as the truncated icosahedron, an Archimedean solid.
        It has 32 faces and 60 vertices.  Rational coordinates are not
        exact.

        EXAMPLES::

            sage: bb = polytopes.buckyball()
            sage: bb.n_vertices()
            60
            sage: bb.n_inequalities()   # number of facets
            32
            sage: bb.field()
            Rational Field
        """
        # Note: field=QQ would give some incorrecty subdivided facets
        p = self.icosahedron(field=RDF).edge_truncation()
        if field==RDF:
            return p
        # Converting with low precision to save time.
        new_ieqs = [[int(1000*x)/QQ(1000) for x in y] for y in p.inequalities()]
        return Polyhedron(ieqs = new_ieqs)


    def pentakis_dodecahedron(self):
        """
        This face-regular, vertex-uniform polytope is dual to the
        truncated icosahedron.  It has 60 faces and 32 vertices.

        EXAMPLES::

            sage: pd = polytopes.pentakis_dodecahedron()
            sage: pd.n_vertices()
            32
            sage: pd.n_inequalities()   # number of facets
            60
        """
        return self.buckyball().polar()


    def twenty_four_cell(self):
        """
        Return the standard 24-cell polytope.

        OUTPUT:

        A Polyhedron object of the 4-dimensional 24-cell, a regular
        polytope. The coordinates of this polytope are exact.

        EXAMPLES::

            sage: p24 = polytopes.twenty_four_cell()
            sage: v = p24.vertex_generator().next()
            sage: for adj in v.neighbors(): print adj
            A vertex at (1/2, 1/2, 1/2, -1/2)
            A vertex at (1/2, 1/2, -1/2, 1/2)
            A vertex at (1/2, -1/2, 1/2, 1/2)
            A vertex at (-1/2, 1/2, 1/2, 1/2)
            A vertex at (0, 0, 0, 1)
            A vertex at (0, 0, 1, 0)
            A vertex at (0, 1, 0, 0)
            A vertex at (1, 0, 0, 0)
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        verts = verts + permutations([0,0,0,1])
        verts = verts + permutations([0,0,0,-1])
        return Polyhedron(vertices = verts)


    def six_hundred_cell(self):
        """
        Return the standard 600-cell polytope.

        OUTPUT:

        A Polyhedron object of the 4-dimensional 600-cell, a regular
        polytope.  In many ways this is an analogue of the
        icosahedron.  The coordinates of this polytope are rational
        approximations of the true coordinates of the 600-cell, some
        of which involve the (irrational) golden ratio.

        EXAMPLES::

            sage: p600 = polytopes.six_hundred_cell() # not tested - very long time
            sage: len(list(p600.bounded_edges())) # not tested - very long time
            120
        """
        verts = []
        q12 = QQ(1)/2
        base = [q12,q12,q12,q12]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        verts.append([x for x in base])
                        base[3] = base[3]*(-1)
                    base[2] = base[2]*(-1)
                base[1] = base[1]*(-1)
            base[0] = base[0]*(-1)
        for x in permutations([0,0,0,1]):
            verts.append(x)
        for x in permutations([0,0,0,-1]):
            verts.append(x)
        g = QQ(1618033)/1000000 # Golden ratio approximation
        verts = verts + [i([q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,1/(g*2),0]) for i in AlternatingGroup(4)]
        verts = verts + [i([-q12,-g/2,-1/(g*2),0]) for i in AlternatingGroup(4)]
        return Polyhedron(vertices = verts)


    def cyclic_polytope(self, dim_n, points_n, field=QQ):
        """
        Return a cyclic polytope.

        INPUT:

        - ``dim_n`` -- positive integer. the dimension of the polytope.

        - ``points_n`` -- positive integer. the number of vertices.

        - ``field`` -- either ``QQ`` (default) or ``RDF``.

        OUTPUT:

        A cyclic polytope of dim_n with points_n vertices on the
        moment curve ``(t,t^2,...,t^n)``, as Polyhedron object.

        EXAMPLES::

            sage: c = polytopes.cyclic_polytope(4,10)
            sage: c.n_inequalities()
            35
        """
        verts = [[t**i for i in range(1,dim_n+1)] for t in range(points_n)]
        return Polyhedron(vertices=verts, field=field)


    def hypersimplex(self, dim_n, k, project = True):
        """
        The hypersimplex in dimension dim_n with d choose k vertices,
        projected into (dim_n - 1) dimensions.

        INPUT:

        - ``n`` -- the numbers ``(1,...,n)`` are permuted

        - ``project`` -- If ``False``, the polyhedron is left in
          dimension ``n``.

        OUTPUT:

        A Polyhedron object representing the hypersimplex.

        EXAMPLES::

            sage: h_4_2 = polytopes.hypersimplex(4,2) # combinatorially equivalent to octahedron
            sage: h_4_2.n_vertices()
            6
            sage: h_4_2.n_inequalities()
            8
        """
        vert0 = [0]*(dim_n-k) + [1]*k
        verts = permutations(vert0)
        if project:
            verts = [Polytopes.project_1(x) for x in verts]
        return Polyhedron(vertices = verts)


    def permutahedron(self, n, project = True):
        """
        The standard permutahedron of (1,...,n) projected into n-1
        dimensions.

        INPUT:

        - ``n`` -- the numbers ``(1,...,n)`` are permuted

        - ``project`` -- If ``False`` the polyhedron is left in dimension ``n``.

        OUTPUT:

        A Polyhedron object representing the permutahedron.

        EXAMPLES::

            sage: perm4 = polytopes.permutahedron(4)
            sage: perm4
            A 3-dimensional polyhedron in QQ^3 defined as the convex hull of 24 vertices.
            sage: polytopes.permutahedron(5).show()    # long time
        """
        verts = range(1,n+1)
        verts = permutations(verts)
        if project:
            verts = [Polytopes.project_1(x) for x in verts]
        p = Polyhedron(vertices = verts)
        return p


    def n_cube(self, dim_n):
        """
        Return a cube in the given dimension

        INPUT:

        - ``dim_n`` -- integer. The dimension of the cube.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional cube, with
        exact coordinates.

        EXAMPLES::

            sage: four_cube = polytopes.n_cube(4)
            sage: four_cube.is_simple()
            True
        """
        if dim_n == 1:
            return Polyhedron(vertices = [[1],[-1]])

        pre_cube = polytopes.n_cube(dim_n-1)
        vertices = [];
        for pre_v in pre_cube.vertex_generator():
            vertices.append( [ 1] + [v for v in pre_v] );
            vertices.append( [-1] + [v for v in pre_v] );
        return Polyhedron(vertices = vertices)


    def cross_polytope(self, dim_n):
        """
        Return a cross-polytope in dimension ``dim_n``. These are
        the generalization of the octahedron.

        INPUT:

        - ``dim_n`` -- integer. The dimension of the cross-polytope.

        OUTPUT:

        A Polyhedron object of the ``dim_n``-dimensional cross-polytope,
        with exact coordinates.

        EXAMPLES::

            sage: four_cross = polytopes.cross_polytope(4)
            sage: four_cross.is_simple()
            False
            sage: four_cross.n_vertices()
            8
        """
        verts = permutations([0 for i in range(dim_n-1)] + [1])
        verts = verts + permutations([0 for i in range(dim_n-1)] + [-1])
        return Polyhedron(vertices = verts)


    def parallelotope(self, generators):
        r"""
        Return the parallelotope spanned by the generators.

        INPUT:

        - ``generators`` -- an iterable of anything convertible to vector
          (for example, a list of vectors) such that the vectors all
          have the same dimension.

        OUTPUT:

        The parallelotope. This is the multi-dimensional
        generalization of a parallelogram (2 generators) and a
        parallelepiped (3 generators).

        EXAMPLES::

            sage: polytopes.parallelotope([ (1,0), (0,1) ])
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 vertices.
            sage: polytopes.parallelotope([[1,2,3,4],[0,1,0,7],[3,1,0,2],[0,0,1,0]])
            A 4-dimensional polyhedron in QQ^4 defined as the convex hull of 16 vertices.
        """
        try:
            generators = [ vector(QQ,v) for v in generators ]
            field = QQ
        except TypeError:
            generators = [ vector(RDF,v) for v in generators ]
            field = RDF

        from sage.combinat.combination import Combinations
        par =  [ 0*generators[0] ]
        par += [ sum(c) for c in Combinations(generators) if c!=[] ]
        return Polyhedron(vertices=par, field=field)



polytopes = Polytopes()



#########################################################################
class PolyhedronFace(SageObject):
    r"""
    A face of a polyhedron.

    This class is for use in
    :meth:`~sage.geometry.polyhedra.Polyhedron.face_lattice`.

    INPUT:

    No checking is performed whether the H/V-representation indices
    actually determine a face of the polyhedron. You should not
    manually create :class:`PolyhedronFace` objects unless you know
    what you are doing.

    OUTPUT:

    A :class:`PolyhedronFace`.

    EXAMPLES::

        sage: octahedron = polytopes.cross_polytope(3)
        sage: inequality = octahedron.Hrepresentation(2)
        sage: face_h = tuple([ inequality ])
        sage: face_v = tuple( inequality.incident() )
        sage: face_h_indices = [ h.index() for h in face_h ]
        sage: face_v_indices = [ v.index() for v in face_v ]
        sage: from sage.geometry.polyhedra import PolyhedronFace
        sage: face = PolyhedronFace(octahedron, face_v_indices, face_h_indices)
        sage: face
        <0,4,5>
        sage: face.dim()
        2
        sage: face.Hrepresentation()
        [An inequality (1, 1, -1) x + 1 >= 0]
        sage: face.Vrepresentation()
        [A vertex at (0, 0, 1), A vertex at (0, -1, 0), A vertex at (-1, 0, 0)]
    """

    def __init__(self, polyhedron, V_indices, H_indices):
        r"""
        The constructor.

        See :class:`PolyhedronFace` for more information.

        INPUT:

        - ``polyhedron`` -- a :class:`Polyhedron`. The ambient
          polyhedron.

        - ``V_indices`` -- list of integers. The indices of the
          face-spanning V-representation objects in the ambient
          polyhedron.

        - ``H_indices`` -- list of integers. The indices of the
          H-representation objects of the ambient polyhedron that are
          saturated on the face.

        TESTS::

            sage: from sage.geometry.polyhedra import PolyhedronFace
            sage: PolyhedronFace(Polyhedron(), [], [])   # indirect doctest
            <>
        """
        self._polyhedron = polyhedron
        self._Vrepresentation = Sequence([ polyhedron.Vrepresentation(i) for i in V_indices ],
                                         universe=Objects())
        self._Hrepresentation = Sequence([ polyhedron.Hrepresentation(i) for i in H_indices ],
                                         universe=Objects())
        self._Vrepresentation.set_immutable()
        self._Hrepresentation.set_immutable()


    def Hrepresentation(self, index=None):
        r"""
        Return the H-representation objects defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a sequence of
        H-representation objects. Each entry is either an inequality
        or an equation.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the sequence is returned.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: for fl in square.face_lattice():
            ...       print fl.element.Hrepresentation()
            ...
            [An inequality (0, 1) x + 1 >= 0, An inequality (1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0]
            [An inequality (0, -1) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0]
            [An inequality (1, 0) x + 1 >= 0, An inequality (0, -1) x + 1 >= 0]
            [An inequality (0, 1) x + 1 >= 0, An inequality (-1, 0) x + 1 >= 0]
            [An inequality (0, 1) x + 1 >= 0, An inequality (1, 0) x + 1 >= 0]
            [An inequality (0, 1) x + 1 >= 0]
            [An inequality (1, 0) x + 1 >= 0]
            [An inequality (0, -1) x + 1 >= 0]
            [An inequality (-1, 0) x + 1 >= 0]
            []
        """
        if index==None:
            return self._Hrepresentation
        else:
            return self._Hrepresentation[index]


    def Vrepresentation(self, index=None):
        r"""
        Return the V-representation objects defining the face.

        INPUT:

        - ``index`` -- optional. Either an integer or ``None``
          (default).

        OUTPUT:

        If the optional argument is not present, a sequence of
        V-representation objects. Each entry is either a vertex, a
        ray, or a line.

        If the optional integer ``index`` is specified, the
        ``index``-th element of the sequence is returned.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: for fl in square.face_lattice():
            ...       print fl.element.Vrepresentation()
            ...
            []
            [A vertex at (1, 1)]
            [A vertex at (-1, 1)]
            [A vertex at (1, -1)]
            [A vertex at (-1, -1)]
            [A vertex at (1, -1), A vertex at (-1, -1)]
            [A vertex at (-1, 1), A vertex at (-1, -1)]
            [A vertex at (1, 1), A vertex at (-1, 1)]
            [A vertex at (1, 1), A vertex at (1, -1)]
            [A vertex at (1, 1), A vertex at (-1, 1), A vertex at (1, -1), A vertex at (-1, -1)]
        """
        if index==None:
            return self._Vrepresentation
        else:
            return self._Vrepresentation[index]


    def dim(self):
        """
        Return the dimension of the face.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: fl = polytopes.dodecahedron().face_lattice()
            sage: [ x.element.dim() for x in fl ]
            [-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3]
        """
        if '_dim' in self.__dict__:
            return self._dim

        if len(self._Vrepresentation)==0:
            self._dim = -1
            return self._dim

        origin = vector(self._Vrepresentation[0])
        v_list = [ vector(v)-origin for v in self._Vrepresentation ]
        self._dim = matrix(v_list).rank()
        return self._dim


    def _repr_(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string listing the V-representation indices of the face.

        EXAMPLES::

            sage: square = polytopes.n_cube(2)
            sage: a_face = list( square.face_lattice() )[8].element
            sage: a_face.__repr__()
            '<0,2>'
        """
        s = '<'
        s += ','.join([ str(v.index()) for v in self.Vrepresentation() ])
        s += '>'
        return s


def Hasse_diagram_from_incidences(atom_to_coatoms, coatom_to_atoms,
                                  face_constructor=None,
                                  required_atoms=None,
                                  key = None,
                                  **kwds):
    r"""
    Compute the Hasse diagram of an atomic and coatomic lattice.

    INPUT:

    - ``atom_to_coatoms`` -- list, ``atom_to_coatom[i]`` should list all
      coatoms over the ``i``-th atom;

    - ``coatom_to_atoms`` -- list, ``coatom_to_atom[i]`` should list all
      atoms under the ``i``-th coatom;

    - ``face_constructor`` -- function or class taking as the first two
      arguments sorted :class:`tuple` of integers and any keyword arguments.
      It will be called to construct a face over atoms passed as the first
      argument and under coatoms passed as the second argument. Default
      implementation will just return these two tuples as a tuple;

    - ``required_atoms`` -- list of atoms (default:None). Each
      non-empty "face" requires at least on of the specified atoms
      present. Used to ensure that each face has a vertex.

    - ``key`` -- any hashable value (default: None). It is passed down
      to :class:`sage.combinat.posets.posets.FinitePoset`.

    - all other keyword arguments will be passed to ``face_constructor`` on
      each call.

    OUTPUT:

    - :class:`finite poset <sage.combinat.posets.posets.FinitePoset>` with
      elements constructed by ``face_constructor``.

    .. NOTE::

        In addition to the specified partial order, finite posets in Sage have
        internal total linear order of elements which extends the partial one.
        This function will try to make this internal order to start with the
        bottom and atoms in the order corresponding to ``atom_to_coatoms`` and
        to finish with coatoms in the order corresponding to
        ``coatom_to_atoms`` and the top. This may not be possible if atoms and
        coatoms are the same, in which case the preference is given to the
        first list.

    ALGORITHM:

    The detailed description of the used algorithm is given in [KP2002]_.

    The code of this function follows the pseudo-code description in the
    section 2.5 of the paper, although it is mostly based on frozen sets
    instead of sorted lists - this makes the implementation easier and should
    not cost a big performance penalty. (If one wants to make this function
    faster, it should be probably written in Cython.)

    While the title of the paper mentions only polytopes, the algorithm (and
    the implementation provided here) is applicable to any atomic and coatomic
    lattice if both incidences are given, see Section 3.4.

    In particular, this function can be used for strictly convex cones and
    complete fans.

    REFERENCES:

    ..  [KP2002]
        Volker Kaibel and Marc E. Pfetsch,
        "Computing the Face Lattice of a Polytope from its Vertex-Facet
        Incidences", Computational Geometry: Theory and Applications,
        Volume 23, Issue 3 (November 2002), 281-290.
        Available at http://portal.acm.org/citation.cfm?id=763203
        and free of charge at http://arxiv.org/abs/math/0106043

    AUTHORS:

    - Andrey Novoseltsev (2010-05-13) with thanks to Marshall Hampton for the
      reference.

    EXAMPLES:

    Let's construct the Hasse diagram of a lattice of subsets of {0, 1, 2}.
    Our atoms are {0}, {1}, and {2}, while our coatoms are {0,1}, {0,2}, and
    {1,2}. Then incidences are ::

        sage: atom_to_coatoms = [(0,1), (0,2), (1,2)]
        sage: coatom_to_atoms = [(0,1), (0,2), (1,2)]

    and we can compute the Hasse diagram as ::

        sage: L = sage.geometry.cone.Hasse_diagram_from_incidences(
        ...                       atom_to_coatoms, coatom_to_atoms)
        sage: L
        Finite poset containing 8 elements
        sage: for level in L.level_sets(): print level
        [((), (0, 1, 2))]
        [((0,), (0, 1)), ((1,), (0, 2)), ((2,), (1, 2))]
        [((0, 1), (0,)), ((0, 2), (1,)), ((1, 2), (2,))]
        [((0, 1, 2), ())]

    For more involved examples see the *source code* of
    :meth:`sage.geometry.cone.ConvexRationalPolyhedralCone.face_lattice` and
    :meth:`sage.geometry.fan.RationalPolyhedralFan._compute_cone_lattice`.
    """
    from sage.graphs.all import DiGraph
    from sage.combinat.posets.posets import FinitePoset
    def default_face_constructor(atoms, coatoms, **kwds):
        return (atoms, coatoms)
    if face_constructor is None:
        face_constructor = default_face_constructor
    atom_to_coatoms = [frozenset(atc) for atc in atom_to_coatoms]
    A = frozenset(range(len(atom_to_coatoms)))  # All atoms
    coatom_to_atoms = [frozenset(cta) for cta in coatom_to_atoms]
    C = frozenset(range(len(coatom_to_atoms)))  # All coatoms
    # Comments with numbers correspond to steps in Section 2.5 of the article
    L = DiGraph(1)       # 3: initialize L
    faces = dict()
    atoms = frozenset()
    coatoms = C
    faces[atoms, coatoms] = 0
    next_index = 1
    Q = [(atoms, coatoms)]              # 4: initialize Q with the empty face
    while Q:                            # 5
        q_atoms, q_coatoms = Q.pop()    # 6: remove some q from Q
        q = faces[q_atoms, q_coatoms]
        # 7: compute H = {closure(q+atom) : atom not in atoms of q}
        H = dict()
        candidates = set(A.difference(q_atoms))
        for atom in candidates:
            coatoms = q_coatoms.intersection(atom_to_coatoms[atom])
            atoms = A
            for coatom in coatoms:
                atoms = atoms.intersection(coatom_to_atoms[coatom])
            H[atom] = (atoms, coatoms)
        # 8: compute the set G of minimal sets in H
        minimals = set([])
        while candidates:
            candidate = candidates.pop()
            atoms = H[candidate][0]
            if atoms.isdisjoint(candidates) and atoms.isdisjoint(minimals):
                minimals.add(candidate)
        # Now G == {H[atom] : atom in minimals}
        for atom in minimals:   # 9: for g in G:
            g_atoms, g_coatoms = H[atom]
            if not required_atoms is None:
                if g_atoms.isdisjoint(required_atoms):
                    continue
            if (g_atoms, g_coatoms) in faces:
                g = faces[g_atoms, g_coatoms]
            else:               # 11: if g was newly created
                g = next_index
                faces[g_atoms, g_coatoms] = g
                next_index += 1
                Q.append((g_atoms, g_coatoms))  # 12
            L.add_edge(q, g)                    # 14
    # End of algorithm, now construct a FinitePoset.
    # In principle, it is recommended to use Poset or in this case perhaps
    # even LatticePoset, but it seems to take several times more time
    # than the above computation, makes unnecessary copies, and crashes.
    # So for now we will mimic the relevant code from Poset.

    # Enumeration of graph vertices must be a linear extension of the poset
    new_order = L.topological_sort()
    # Make sure that coatoms are in the end in proper order
    tail = [faces[atoms, frozenset([coatom])]
            for coatom, atoms in enumerate(coatom_to_atoms)]
    tail.append(faces[A, frozenset()])
    new_order = [n for n in new_order if n not in tail] + tail
    # Make sure that atoms are in the beginning in proper order
    head = [0] # We know that the empty face has index 0
    head.extend(faces[frozenset([atom]), coatoms]
                for atom, coatoms in enumerate(atom_to_coatoms)
                if required_atoms is None or atom in required_atoms)
    new_order = head + [n for n in new_order if n not in head]
    # "Invert" this list to a dictionary
    labels = dict()
    for new, old in enumerate(new_order):
        labels[old] = new
    L.relabel(labels)
    # Construct the actual poset elements
    elements = [None] * next_index
    for face, index in faces.items():
        atoms, coatoms = face
        elements[labels[index]] = face_constructor(
                        tuple(sorted(atoms)), tuple(sorted(coatoms)), **kwds)
    return FinitePoset(L, elements, key = key)
