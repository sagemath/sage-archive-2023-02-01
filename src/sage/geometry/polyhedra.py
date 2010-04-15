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


REFERENCES:

    Komei Fukuda's `FAQ in Polyhedral Computation
    <http://www.ifor.math.ethz.ch/~fukuda/polyfaq/polyfaq.html>`_

AUTHORS:

    - Marshall Hampton: first version, bug fixes, and various improvements, 2008 and 2009
    - Arnaud Bergeron: improvements to triangulation and rendering, 2008
    - Sebastien Barthelemy: documentation improvements, 2008
    - Volker Braun: refactoring, handle non-compact case, 2009

TESTS::

    TestSuite(s).run()
"""


########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence

from subprocess import Popen, PIPE
from sage.misc.all import SAGE_TMP, tmp_filename, union
from sage.misc.functional import norm
from sage.misc.package import is_package_installed

from sage.rings.all import Integer, QQ, ZZ, primes_first_n
from sage.rings.rational import Rational
from sage.rings.real_double import RDF
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix, identity_matrix
from sage.functions.other import sqrt
from sage.functions.trig import sin, cos

from sage.plot.all import point2d, line2d, arrow, polygon2d
from sage.plot.plot3d.all import point3d, line3d, arrow3d, polygon3d
from sage.graphs.graph import Graph

from sage.combinat.posets.posets import Poset
from sage.combinat.combinat import permutations
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
        """
        return Vobj.evaluated_on(self)

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
    is_inequality.__doc__ = Hrepresentation.is_inequality.__doc__

    def _repr_(self):
        """
        The string representation of the inequality.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.inequality_generator().next()
            sage: a._repr_()
            'An inequality (-1, 0, 0) x + 1 >= 0'
        """
        s = 'An inequality '
        s += repr(self.A()) + ' x + ' + repr(self.b()) + ' >= 0'
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
        """
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

    is_equation.__doc__ = Hrepresentation.is_equation.__doc__

    def _repr_(self):
        """
        A string representation of this object.

        TESTS::

            sage: p = Polyhedron(vertices = [[0,0,0],[1,1,0],[1,2,0]])
            sage: a = p.equation_generator().next()
            sage: a._repr_()
            'An equation (0, 0, 1) x + 0 == 0'
        """
        s = 'An equation '
        s += repr(self.A()) + ' x + ' + repr(self.b()) + ' == 0'
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
    is_vertex.__doc__ = Vrepresentation.is_vertex.__doc__

    def _repr_(self):
        """
        Returns a string representation of the vertex.

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
    is_ray.__doc__ = Vrepresentation.is_ray.__doc__

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

    is_line.__doc__ = Vrepresentation.is_line.__doc__

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
    Returns 2d rendering of the projection of a polyhedron into
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
    Returns 3d rendering of a polyhedron projected into
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
        projection.render_vertices_3d(thickness=3, color='green', **kwds) +\
        projection.render_wireframe_3d(thickness=3, color='green', **kwds) + \
        projection.render_solid_3d(**kwds)


def render_4d(polyhedron, **kwds):
    """
    Returns a 3d rendering of the Schlegel projection of a 4d
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
                 field = QQ, verbose = False):
        """
        Initializes the polyhedron.

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
        self._init_from_cdd_input(s, verbose=verbose)

        # Only add a show() method if we have one
        if self.ambient_dim() < len(Polyhedron._render_method):
            render = Polyhedron._render_method[self.ambient_dim()]
            if render != None:
                self.show = lambda **kwds: render(self,**kwds)
                self.show.__doc__ = render.__doc__

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
        Initializes an empty polyhedron.  This is not implemented yet.

        TESTS::

            sage: Polyhedron()
            The empty polyhedron in QQ^0.
            sage: Polyhedron(vertices = [])
            The empty polyhedron in QQ^0.
            sage: Polyhedron()._init_empty_polyhedron()
        """
        self._ambient_dim = 0

        self._Vrepresentation = Sequence([])
        self._Vrepresentation.set_immutable()

        self._Hrepresentation = Sequence([])
        self._Hrepresentation.set_immutable()

        self._V_adjacency_matrix = matrix(ZZ, 0, 0, 0)
        self._V_adjacency_matrix.set_immutable()

        self._H_adjacency_matrix = matrix(ZZ, 0, 0, 0)
        self._H_adjacency_matrix.set_immutable()


    def _init_from_cdd_input(self, cdd_input_string, verbose=False):
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

        cdd_proc = Popen(self._cdd_executable,
                         stdin=PIPE, stdout=PIPE, stderr=PIPE)
        ans, err = cdd_proc.communicate(input=cdd_input_string)

        # FIXME: check err
        if verbose: print ans

        self._init_from_cdd_output(ans)


    def _init_from_cdd_output(self, cdd_output_string):
        """
        Internal method: initialize ourselves with the output from cdd.

        TESTS::

            sage: from sage.geometry.polyhedra import cdd_Vrepresentation
            sage: s = cdd_Vrepresentation('rational',[[0,0],[1,0],[0,1],[1,1]], [], [])
            sage: from subprocess import Popen, PIPE
            sage: cdd_proc = Popen('cdd_both_reps_gmp', stdin=PIPE, stdout=PIPE, stderr=PIPE)
            sage: ans, err = cdd_proc.communicate(input=s)
            sage: p = Polyhedron(vertices = [[0,0],[1,0],[0,1],[1,1]])
            sage: p._init_from_cdd_output(ans)
            sage: p.vertices()
            [[0, 0], [1, 0], [0, 1], [1, 1]]
        """
        cddout = cdd_output_string.splitlines()

        def skip_to(expected_string):
            """
            Finds the expected string in a list of strings.
            """
            while True:
                if len(cddout)==0:
                    raise ValueError, "Did not find expected string."
                l = cddout.pop(0)
                l.strip();
                if l==expected_string: return

        def cdd_linearities():
            l = cddout[0].split()
            if l[0] != "linearity":
                return []
            cddout.pop(0)
            assert len(l) == int(l[1])+2, "Not enough linearities given"
            return [int(i)-1 for i in l[2:]]  # make indices pythonic

        def cdd_convert(string, field=self.field()):
            """
            Converts the cdd output string to a numerical value.
            """
            return [field(x) for x in string.split()]

        skip_to('V-representation')
        self._Vrepresentation = Sequence([])
        lines = cdd_linearities()
        skip_to('begin')
        l = cddout.pop(0).split()
        self._ambient_dim = int(l[1])-1
        no_vertex = True
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
                no_vertex = False
        if no_vertex:
            Vertex(self, [0] * self._ambient_dim)
        self._Vrepresentation.set_immutable()
        skip_to('end')

        skip_to('H-representation')
        self._Hrepresentation = Sequence([])
        equations = cdd_linearities()
        skip_to('begin')
        l = cddout.pop(0).split()
        assert self._ambient_dim == int(l[1])-1, "Different ambient dimension?"
        for i in range(int(l[0])):
            l = cddout.pop(0)
            if i in equations:
                Equation(self, cdd_convert(l));
            else:
                Inequality(self, cdd_convert(l));
        self._Hrepresentation.set_immutable()
        skip_to('end')

        def cdd_adjacencies():
            l = cddout.pop(0).split()
            assert l[2] == ':', "Not a line of the adjacency data?"
            return [int(i)-1 for i in l[3:]]

        skip_to('Vertex graph')
        n = len(self._Vrepresentation);
        if no_vertex:
            n_cdd=n-1;
        else:
            n_cdd=n;
        self._V_adjacency_matrix = matrix(ZZ, n, n, 0)
        skip_to('begin')
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
        if no_vertex: # cdd implied that there is only one vertex
            for i in range(n-1):
                self._V_adjacency_matrix[i,n-1] = 1
                self._V_adjacency_matrix[n-1,i] = 1
        self._V_adjacency_matrix.set_immutable()
        skip_to('end')

        skip_to('Facet graph')
        n = len(self._Hrepresentation);
        self._H_adjacency_matrix = matrix(ZZ, n, n, 0)
        skip_to('begin')
        l = cddout.pop(0).split()
        assert int(l[0]) == n, "Not enough H-adjacencies in cdd output?"
        for i in range(n):
            for a in cdd_adjacencies():
                self._H_adjacency_matrix[i,a] = 1
        self._H_adjacency_matrix.set_immutable()
        skip_to('end')


    def _repr_(self):
        """
        Returns a description of the polyhedron.

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
        Writes the inequalities/equations data of the polyhedron in
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
        Writes the vertices/rays/lines data of the polyhedron in cdd's
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
        Returns the number of equations. The representation will
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
        Returns the number of inequalities. The representation will
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
        Returns the number of facets in the polyhedron.  This is the
        same as the n_inequalities function.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        return self.n_inequalities()


    def n_vertices(self):
        """
        Returns the number of vertices. The representation will
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
        Returns the number of rays. The representation will
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
        Returns the number of lines. The representation will
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


    def Hrepresentation(self, *args):
        """
        Returns the objects of the H-representaton. Each entry is
        either an inequality or a equation.

        INPUT:

        The optional argument is an index in
        0...n_Hrepresentations(). Without an argument, returns the
        list of all H-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p.Hrepresentation(0)
            An inequality (0, 0, 1) x + 1 >= 0
            sage: p.Hrepresentation(0) == p.Hrepresentation() [0]
            True
        """
        nargs=len(args)
        if nargs==0:
            return self._Hrepresentation
        elif nargs==1:
            return self._Hrepresentation[args[0]]
        raise ValueError, "Hrepresentation() requires 0 or 1 argument."


    def Hrep_generator(self):
        """
        Returns an iterator over the objects of the H-representation
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
        Returns the number of objects that make up the
        H-representation of the polyhedron.

        EXAMPLES::

            sage: p = polytopes.cross_polytope(4)
            sage: p.n_Hrepresentation()
            16
            sage: p.n_Hrepresentation() == p.n_inequalities() + p.n_equations()
            True
        """
        return len(self.Hrepresentation())


    def Vrepresentation(self, *args):
        """
        Returns the objects of the V-representation. Each entry is
        either a vertex, a ray, or a line.

        INPUT:

        The optional argument is an index in
        0...n_Hrepresentations(). Without an argument, returns the
        list of all H-representation objects.

        EXAMPLES::

            sage: p = polytopes.n_simplex(4)
            sage: p.Vrepresentation(0)
            A vertex at (0, 0, 0, -44721/50000)
            sage: p.Vrepresentation(0) == p.Vrepresentation() [0]
            True
        """
        nargs=len(args)
        if nargs==0:
            return self._Vrepresentation
        elif nargs==1:
            return self._Vrepresentation[args[0]]
        raise ValueError, "Vrepresentation() requires 0 or 1 argument."


    def n_Vrepresentation(self):
        """
        Returns the number of objects that make up the
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
        Returns a list of face indices (i.e. indices of
        H-representation objects) and the indices of faces adjacent to
        them.

        NOTES:

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
        Returns the face-vertex incidences in the form `[f_i, [v_{i_0}, v_{i_1},\dots ,v_{i_2}]]`.

        NOTES:

        Instead of working with face/vertex indices, you can use the
        H-representation/V-representation objects directly (see
        examples).

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
        Returns a list of vertex indices and their adjacent vertices.

        NOTES:

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
        Returns the vertex-face incidences in the form `[v_i, [f_{i_0}, f_{i_1},\dots ,f_{i_2}]]`.

        NOTES:

        Instead of working with face/vertex indices, you can use the
        H-representation/V-representation objects directly (see
        examples).

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


    def graph(self):
        """
        Returns a graph in which the vertices correspond to vertices of the polyhedron,
        and edges to edges.

        OUTPUT:

        A graph

        EXAMPLES::

            sage: g3 = polytopes.n_cube(3).graph()
            sage: len(g3.automorphism_group())
            48
            sage: s4 = polytopes.n_simplex(4).graph()
            sage: s4.is_eulerian()
            True
        """
        return Graph(self.vertex_adjacency_matrix(), loops=True)


    def inequality_generator(self):
        """
        Returns  a generator for the defining inequalities of the
        polyhedron.

        OUTPUT:

        A generator of the inequality Hrepresentation objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.inequality_generator(): print(v)
            An inequality (-1, 0) x + 1 >= 0
            An inequality (0, -1) x + 1 >= 0
            An inequality (1, 1) x + -1 >= 0
            sage: [ v for v in triangle.inequality_generator() ]
            [An inequality (-1, 0) x + 1 >= 0,
             An inequality (0, -1) x + 1 >= 0,
             An inequality (1, 1) x + -1 >= 0]
            sage: [ [v.A(), v.b()] for v in triangle.inequality_generator() ]
            [[(-1, 0), 1], [(0, -1), 1], [(1, 1), -1]]
        """
        for H in self.Hrepresentation():
            if H.is_inequality():
                yield H


    def inequalities(self):
        """
        Returns a list of inequalities as coefficient lists.

        NOTES:

        It is recommended to use ``inequality_generator()`` instead to
        iterate over the list of ``Inequality`` objects.

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
        Returns a generator for the linear equations satisfied by the
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
        Returns the linear constraints of the polyhedron. As with
        inequalities, each constraint is given as [b -a1 -a2 ... an]
        where for variables x1, x2,..., xn, the polyhedron satisfies
        the equation b = a1*x1 + a2*x2 + ... + an*xn.

        NOTES:

        It is recommended to use ``equation_generator()`` instead to
        iterate over the list of ``Equation`` objects.

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
        Returns a list of vertices of the polyhedron.

        NOTES:

        It is recommended to use ``vertex_generator()`` instead to
        iterate over the list of ``Vertex`` objects.

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
        Returns a generator for the vertices of the polyhedron.

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
        Returns a generator for the rays of the polyhedron.

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
        Returns a list of rays as coefficient lists.

        NOTES:

        It is recommended to use ``ray_generator()`` instead to
        iterate over the list of ``Ray`` objects.

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
        Returns a generator for the lines of the polyhedron.

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
        Returns a list of lines of the polyhedron.  The line data is given
        as a list of coordinates rather than as a Hrepresentation object.

        NOTES:

        It is recommended to use ``line_generator()`` instead to
        iterate over the list of ``Line`` objects.

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
        Returns the bounded edges (excluding rays and lines).

        OUTPUT:

            A generator for pairs of vertices, one pair per edge.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0],[0,1]], rays=[[1,0],[0,1]])
            sage: [ e for e in p.bounded_edges() ]
            [[A vertex at (1, 0), A vertex at (0, 1)]]
            sage: for e in p.bounded_edges(): print e
            [A vertex at (1, 0), A vertex at (0, 1)]
        """
        obj = self.Vrepresentation()
        edges = []
        for i in range(len(obj)):
            if not obj[i].is_vertex(): continue
            for j in range(i+1,len(obj)):
                if not obj[j].is_vertex(): continue
                if self.vertex_adjacency_matrix()[i,j] == 0: continue
                yield [obj[i], obj[j]]


    def ambient_dim(self):
        r"""
        Returns the dimension of the ambient space.

        EXAMPLES::

            sage: poly_test = Polyhedron(vertices = [[1,0,0,0],[0,1,0,0]])
            sage: poly_test.ambient_dim()
            4
        """
        return self._ambient_dim


    def dim(self):
        """
        Returns the dimension of the polyhedron.

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
        This is an alias for ``vertex_adjacency_matrix()``

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
        Returns the binary matrix of vertex adjacencies.

        EXAMPLES::

            sage: polytopes.n_simplex(4).vertex_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        return self._V_adjacency_matrix;


    def facet_adjacency_matrix(self):
        """
        Returns the adjacency matrix for the facets and hyperplanes.

        EXAMPLES::

            sage: polytopes.n_simplex(4).facet_adjacency_matrix()
            [0 1 1 1 1]
            [1 0 1 1 1]
            [1 1 0 1 1]
            [1 1 1 0 1]
            [1 1 1 1 0]
        """
        return self._H_adjacency_matrix;


    def incidence_matrix(self):
        """
        Returns the incidence matrix.

        NOTES:

        The columns correspond to inequalities/equations in the order
        ``self.Hrepresentation()``, the rowns correspond to
        vertices/rays/lines in the order ``self.Vrepresentation()``

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


    def center(self):
        """
        Returns the average of the vertices. Returns the origin for
        the empty polytope. All rays and lines are ignored.

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
        Returns the square of the maximal distance from the center to
        a vertex. All rays and lines are ignored.

        EXAMPLES::

            sage: p = polytopes.permutahedron(4, project = False)
            sage: p.radius_square()
            720
        """
        try:
            return self._radius_2
        except AttributeError:
            self._radius_2 = 0
            for v in self.vertex_generator():
                self._radius_2 += v.vector()*v.vector()
            return self._radius_2


    def radius(self):
        """
        Returns the maximal distance from the center to a vertex. All
        rays and lines are ignored.

        NOTES:

        The radius for a rational vector is, in general, not rational.
        use ``radius_square()`` if you need a rational distance
        measure.

        EXAMPLES::

            sage: p = polytopes.n_cube(4)
            sage: p.radius()
            8
        """
        return sqrt(self.radius_square())


    def is_compact(self):
        """
        Tests for boundedness of the polytope.

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
        Tests for simplicity of a polytope.

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
        Returns the Gale transform of a polytope as described in the
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
        if not self.is_compact(): raise ValueError, "Not a polytope"

        A = matrix(self.n_vertices(),
                   [ [1]+list(x) for x in self.vertex_generator()])
        A = A.transpose()
        A_ker = A.right_kernel()
        return A_ker.basis_matrix().transpose().rows()


    def triangulated_facial_incidences(self):
        """
        Returns a list of the form [face_index, [v_i_0,
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
            [[0, [0, 2, 3]], [1, [0, 1, 2]], [2, [1, 2, 3]], [3, [0, 1, 3]]]

        Otherwise some faces get split up to triangles::

            p = polytopes.regular_polygon(5)
            sage: Polyhedron(vertices = [[2,0,0],[4,1,0],[0,5,0],[5,5,0],[1,1,0],[0,0,1]]).triangulated_facial_incidences()
            [[0, [0, 1, 5]], [1, [0, 4, 5]], [2, [2, 4, 5]], [3, [2, 3, 5]], [4, [1, 3, 5]], [5, [0, 1, 4]], [5, [1, 4, 3]], [5, [4, 3, 2]]]
        """
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
        Returns a simplicial complex from a triangulation of the polytope.

        Warning: This first triangulates the polytope using
        ``triangulated_facial_incidences``, and this function may fail
        in dimensions greater than 3, although it usually doesn't.

        OUTPUT:

            A simplicial complex.

        EXAMPLES::

            sage: p = polytopes.cuboctahedron()
            sage: sc = p.simplicial_complex()
            sage: sc
            Simplicial complex with 13 vertices and 20 facets
        """
        from sage.homology.simplicial_complex import SimplicialComplex
        return SimplicialComplex(vertex_set = self.n_vertices(),
                                 maximal_faces = [x[1] for x in self.triangulated_facial_incidences()])

    def __add__(self, other):
        """
        Addition of two polyhedra is defined as their Minkowski sum.

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

        else:  # assume other is a vector and try to add vertices
            displacement = vector(other)
            new_vertices = [list(x() + displacement) for x in self.vertex_generator()]
            new_rays = self.rays()
            new_lines = self.lines()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines, field=self.field())


    def __mul__(self, other):
        """
        Multiplication by another polyhedron returns the product
        polytope.  Multiplication by a scalar returns the polytope
        dilated by that scalar.

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
                          rays=new_rays, lines=new_lines, field=self.field())


    def __rmul__(self,other):
        """
        Multiplication by another polyhedron returns the product
        polytope.  Multiplication by a scalar returns the polytope
        dilated by that scalar.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[t,t^2,t^3] for t in srange(2,4)])
            sage: p2 = 3*p + p
            sage: p2.vertex_generator().next()
            A vertex at (8, 16, 32)
        """
        return self.__mul__(other)


    def union(self, other):
        """
        Returns the convex hull of the set-theoretic union of the two
        polyhedra.

        EXAMPLES::

            sage: a_simplex = polytopes.n_simplex(3)
            sage: verts = a_simplex.vertices()
            sage: verts = [[x[0]*3/5+x[1]*4/5, -x[0]*4/5+x[1]*3/5, x[2]] for x in verts]
            sage: another_simplex = Polyhedron(vertices = verts)
            sage: simplex_union = a_simplex.union(another_simplex)
            sage: simplex_union.n_vertices()
            7
        """
        new_vertices = self.vertices() + other.vertices()

        new_rays = self.rays() + other.rays()

        new_lines = self.lines() + other.lines()

        return Polyhedron(vertices=new_vertices,
                          rays=new_rays, lines=new_lines, field=self.field())


    def intersection(self, other):
        """
        Returns the intersection of one polyhedron with another.

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
                          field=self.field())


    def edge_truncation(self, cut_frac = Integer(1)/3):
        r"""
        Returns a new polyhedron formed from two points on each edge
        between two vertices.

        INPUT:

          - ``cut_frac`` - how deeply to cut into the edge.  Default is
            `\frac{1}{3}`.

        OUTPUT:

        A Polyhedron object, truncated as described above.

        EXAMPLES:

        Truncating a cube::

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
                          lines=new_lines, field=self.field())


    def _v_closure(self, vertex_indices):
        """
        Returns the list of vertex indices for the vertices in the smallest
        face containing the input vertices.

        EXAMPLES::

            sage: p = polytopes.n_cube(3)
            sage: p._v_closure([1,2])
            [0, 1, 2, 3]
        """
        return self._face_vertex_indexset(self._vertex_face_indexset(vertex_indices))


    def _vertex_face_indexset(self, Vrep_indices):
        """
        Returns the maximal list of facet indices whose intersection contains
        the given vertices, i.e. a list of all the facets which contain all of
        the given vertices.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0, 0, 0], [0, 0, 1], [0, 1, 1], [1, 0, 0], [1, 1, 1], [1, 0, 1], [1, 1, 0]])
            sage: list( p.Vrep_generator() )
            [A vertex at (0, 0, 0),
             A vertex at (0, 0, 1),
             A vertex at (0, 1, 1),
             A vertex at (1, 0, 0),
             A vertex at (1, 1, 1),
             A vertex at (1, 0, 1),
             A vertex at (1, 1, 0)]
            sage: list( p.Hrep_generator() )
            [An inequality (0, 0, 1) x + 0 >= 0,
             An inequality (1, -1, 1) x + 0 >= 0,
             An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (-1, 0, 0) x + 1 >= 0]
            sage: p._vertex_face_indexset([2])
            [1, 2, 4, 5]
        """
        if Vrep_indices == []:
            return range(self.n_Vrepresentation())
        Vrep_list = [self.Vrepresentation(i) for i in Vrep_indices]
        faces = []
        for aface in self.Hrep_generator():
            if all([aface.is_incident(v) for v in Vrep_list]):
                faces.append(aface.index())
        return faces


    def _face_vertex_indexset(self, Hrep_indices):
        """
        Returns the maximal list of vertex indices who lie on the given facets.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: p._face_vertex_indexset([0])
            [1, 3, 8]
        """
        if Hrep_indices == []:
            return range(self.n_Hrepresentation())
        Hrep_list = [self.Hrepresentation(i) for i in Hrep_indices]
        verts = []
        for avert in self.Vrep_generator():
            if all([avert.is_incident(h) for h in Hrep_list]):
                verts.append(avert.index())
        return verts


    def face_lattice(self):
        """
        Computes the face-lattice poset. Elements are tuples of
        (vertices, facets) - i.e. this keeps track of both the
        vertices in each face, and all the facets containing them.

        EXAMPLES::

            sage: c5_10 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,11)])
            sage: c5_10_fl = c5_10.face_lattice()
            sage: [len(x) for x in c5_10_fl.level_sets()]
            [1, 10, 45, 100, 105, 42, 1]

        TESTS::

            sage: c5_20 = Polyhedron(vertices = [[i,i^2,i^3,i^4,i^5] for i in range(1,21)]) # not tested - very long time
            sage: c5_20_fl = c5_20.face_lattice() # not tested - very long time
            sage: [len(x) for x in c5_20_fl.level_sets()] # not tested - very long time
            [1, 20, 190, 580, 680, 272, 1]
            sage: polytopes.n_cube(2).face_lattice().plot()

        REFERENCES:

            'Computing the Face Lattice of a Polytope from its
            Vertex-Facet Incidences', by V. Kaibel and M.E. Pfetsch.
        """
        try:
            return self._face_lattice
        except AttributeError:
            pass
        # dictionary of the form: (vertices,faces):(upper cover keys)
        f_l_dict = {(None,None):[(tuple([x]),tuple(self._vertex_face_indexset([x])))
                                 for x in range(self.n_Vrepresentation())]}
        todolist = f_l_dict.values()[0][:]
        while todolist != []:
            todo = todolist.pop()
            f_l_dict[todo] = []
            candidates = []
            for aface in todo[1]:
                H = self.Hrepresentation(aface)
                candidates = union(candidates + [v.index() for v in H.incident()])
            candidates = [x for x in candidates if not x in todo[0]]
            minimals = []
            if candidates == []:
                f_l_dict[todo].append((tuple(range(self.n_Vrepresentation())),
                                       tuple(range(self.n_Hrepresentation()))))
            while candidates != []:
                c = candidates[0]
                closure = self._v_closure(list(todo[0])+[c])
                clos_comp = [x for x in closure if not x in list(todo[0])+[c]]
                clos_comp = [x for x in clos_comp if x in candidates+minimals]
                if clos_comp == []:
                    minimals.append(c)
                    candidates.remove(c)
                    newelm = (tuple(closure), tuple(self._vertex_face_indexset(closure)))
                    f_l_dict[todo].append(newelm)
                    if not f_l_dict.has_key(newelm):
                        todolist.append(newelm)
                else:
                    candidates.remove(c)
        self._face_lattice = Poset(f_l_dict)
        return self._face_lattice


    def f_vector(self):
        r"""
        Returns the f-vector (number of faces in each dimension) of
        the polytope as a list.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2], [3, 2, 1], [0, 0, 0]])
            sage: p.f_vector()
            [1, 7, 12, 7, 1]
        """
        try:
            return self._f_vector
        except AttributeError:
            self._f_vector = [len(x) for x in self.face_lattice().level_sets()]
            return self._f_vector


    def vertex_graph(self):
        """
        Returns a graph in which the vertices correspond to vertices
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


    def polar(self):
        """
        Returns the polar (dual) polytope.  The original vertices are
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
        Returns a polyhedron that is a bipyramid over the original.

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
        Returns a prism of the original polyhedron.

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
        Returns a projection object.

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
        Returns solid rendering of a 2- or 3-d polytope.

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
        For polytopes in 2 or 3 dimensions, returns the edges
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
            raise NotImplemented

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


#############################################################
def cyclic_sort_vertices_2d(Vlist):
    """
    Returns the vertices/rays in cyclic order if possible.

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
        sage: ppoints[0]
        (0.0, 0.0)
    """
    def __init__(self, projection_point):
        """
        Creates a stereographic projection function.  Input is a list
        of coordinates in the appropriate dimension, which is the point
        projected from.

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
        First reflects x with a Householder reflection which takes
        the projection point to (0,...0,self.psize) where psize is
        the length of the projection point, and then dilates by
        1/(zdiff) where zdiff is the difference between the last
        coordinate of x and psize.

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
        Applies the projection to a vector (can be input as a list).

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
        Returns a string describing the projection.

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
        Returns the identity projection.

        EXAMPLES::

            sage: p = polytopes.icosahedron()
            sage: from sage.geometry.polyhedra import Projection
            sage: pproj = Projection(p)
            sage: ppid = pproj.identity()
            sage: ppid.dimension
            3
        """
        return self.__call__(projection_func_identity)
    identity.__doc__ = projection_func_identity.__doc__


    def stereographic(self, projection_point=None):
        r"""
        The stereographic projection.

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
        The Schlegel projection.

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
        Converts a coordinate vector to its internal index.

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
        Converts list of coordinate vectors to the corresponding list
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
        Given a list of indices, returns the projected coordinates.

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
        Returns the points of a polyhedron in 2d.

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
        Returns the outline (edges) of a polyhedron in 2d.

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
        Returns the filled interior (a polygon) of a polyhedron in 2d.

        EXAMPLES::

            sage: cps = [i^3 for i in srange(-2,2,1/5)]
            sage: p = Polyhedron(vertices = [[(t^2-1)/(t^2+1),2*t/(t^2+1)] for t in cps])
            sage: proj = p.projection()
            sage: filled_poly = proj.render_fill_2d()
            sage: filled_poly.axes_width()
            0.8000...
        """
        poly = [polygon2d(self.coordinates_of(p), **kwds)
                 for p in self.polygons]
        return sum(poly)


    def render_vertices_3d(self, **kwds):
        """
        Returns the 3d rendering of the vertices.

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
        Returns the 3d wireframe rendering.

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
        Returns solid 3d rendering of a 3d polytope.

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
        Takes a ndim-dimensional point and projects it onto the plane
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
        Returns a regular polygon with n vertices.  Over the rational
        field the vertices may not be exact.

        INPUT:

          - ``n`` - a positive integer, the number of vertices.

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
        Returns the Birkhoff polytope with n! vertices.  Each vertex
        is a (flattened) n by n permutation matrix.

        INPUT:

          - ``n`` - a positive integer giving the size of the permutation matrices.

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
        Returns a rational approximation to a regular simplex in
        dimension ``dim_n``.

        INPUT:

          - ``dim_n`` - The dimension of the cross-polytope, a positive
                        integer.
          - ``project`` - Optional argument, whether to project
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
        Returns an icosahedron with edge length 1.  The vertices are
        rational, so a rational approximation of the golden ratio
        is used.

        INPUT:

          - ``field`` - Either ``QQ`` or ``RDF``. The icosahedron's
            coordinates contain the golden ratio, so there is no exact
            representation possible.

        OUTPUT:

        A Polyhedron object of a floating point or rational
        approximation to the regular 3d icosahedron.

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
        Returns a dodecahedron.

        INPUT:

          - ``field`` - Either ``QQ`` (in which case a rational
            approximation to the golden ratio is used) or ``RDF``.

        EXAMPLES::

            sage: d12 = polytopes.dodecahedron()
            sage: d12.n_inequalities()
            12
        """
        return self.icosahedron(field=field).polar()


    def small_rhombicuboctahedron(self):
        """
        An Archimedean solid with 24 vertices and 26 faces.

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
        An Archimedean solid with 48 vertices and 26 faces.

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
        Returns the standard 24-cell polytope.

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
        Returns the standard 600-cell polytope.

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
        Returns a cyclic polytope in dimension=dim_n and with points_n
        number of points.

        INPUT:

          - ``dim_n`` - a positive integer, the dimension of the polytope.
          - ``points_n`` - the number of vertices.

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

          - ``n`` - the numbers (1,...,n) are permuted
          - ``project`` - If False the polyhedron is left in dimension ``n``.

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

          - ``n`` - the numbers (1,...,n) are permuted
          - ``project`` - If False the polyhedron is left in dimension ``n``.

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
        Returns a cube in the given dimension

        INPUT:

          - ``dim_n`` - The dimension of the cube, a positive integer.

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
        Returns a cross-polytope in dimension ``dim_n``. These are
        the generalization of the octahedron.

        INPUT:

          - ``dim_n`` - The dimension of the cross-polytope, a positive integer.

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


polytopes = Polytopes()
