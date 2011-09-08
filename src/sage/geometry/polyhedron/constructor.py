r"""
Polyhedra

In this module, a polyhedron is a convex (possibly unbounded) set in
Euclidean space cut out by a finite set of linear inequalities and
linear equations. Note that the dimension of the polyhedron can be
less than the dimension of the ambient space. There are two
complementary representations of the same data:

**H(alf-space/Hyperplane)-representation**
    This describes a polyhedron as the common solution set of a
    finite number of

        * linear inequalities `A \vec{x} + b \geq 0`, and
        * linear equations  `C \vec{x} + d = 0`.


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
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 2 rays
    sage: v = trunc_quadr.vertex_generator().next()  # the first vertex in the internal enumeration
    sage: v
    A vertex at (0, 1)
    sage: v.vector()
    (0, 1)
    sage: list(v)
    [0, 1]
    sage: len(v)
    2
    sage: v[0] + v[1]
    1
    sage: v.is_vertex()
    True
    sage: type(v)
    <class 'sage.geometry.polyhedron.representation.Vertex'>
    sage: type( v() )
    <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>
    sage: v.polyhedron()
    A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 2 vertices and 2 rays
    sage: r = trunc_quadr.ray_generator().next()
    sage: r
    A ray in the direction (0, 1)
    sage: r.vector()
    (0, 1)
    sage: [x for x in v.neighbors()]
    [A ray in the direction (0, 1), A ray in the direction (1, 0), A vertex at (1, 0)]

Inequalities `A \vec{x} + b \geq 0` (and, similarly, equations) are
specified by a list ``[b, A]``::

    sage: Polyhedron(ieqs = [[0,1,0],[0,0,1],[1,-1,-1]]).Hrepresentation()
    (An inequality (-1, -1) x + 1 >= 0,
     An inequality (1, 0) x + 0 >= 0,
     An inequality (0, 1) x + 0 >= 0)

See :func:`Polyhedron` for a detailed description of all possible ways
to construct a polyhedron.

REFERENCES:

    Komei Fukuda's `FAQ in Polyhedral Computation
    <http://www.ifor.math.ethz.ch/~fukuda/polyfaq/polyfaq.html>`_

AUTHORS:

    - Marshall Hampton: first version, bug fixes, and various improvements, 2008 and 2009
    - Arnaud Bergeron: improvements to triangulation and rendering, 2008
    - Sebastien Barthelemy: documentation improvements, 2008
    - Volker Braun: refactoring, handle non-compact case, 2009 and 2010
    - Andrey Novoseltsev: added Hasse_diagram_from_incidences, 2010
    - Volker Braun: rewrite to use PPL instead of cddlib, 2011
"""

########################################################################
#       Copyright (C) 2008 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.rings.all import QQ, ZZ, RDF
from sage.misc.decorators import rename_keyword

from misc import _make_listlist, _common_length_of






#########################################################################
@rename_keyword(deprecation=11634, field='base_ring')
def Polyhedron(vertices=None, rays=None, lines=None,
               ieqs=None, eqns=None,
               ambient_dim=None, base_ring=None, minimize=True, verbose=False,
               backend=None):
    """
    Construct a polyhedron object.

    You may either define it with vertex/ray/line or
    inequalities/equations data, but not both. Redundant data will
    automatically be removed (unless ``minimize=False``), and the
    complementary representation will be computed.

    INPUT:

    - ``vertices`` -- list of point. Each point can be specified as
      any iterable container of ``base_ring`` elements. If ``rays`` or
      ``lines`` are specified but no ``vertices``, the origin is
      taken to be the single vertex.

    - ``rays`` -- list of rays. Each ray can be specified as any
      iterable container of ``base_ring`` elements.

    - ``lines`` -- list of lines. Each line can be specified as any
      iterable container of ``base_ring`` elements.

    - ``ieqs`` -- list of inequalities. Each line can be specified as
      any iterable container of ``base_ring`` elements.

    - ``eqns`` -- list of equalities. Each line can be specified as
      any iterable container of ``base_ring`` elements.

    - ``base_ring`` -- either ``QQ`` or ``RDF``. The field over which
      the polyhedron will be defined. For ``QQ``, exact arithmetic
      will be used. For ``RDF``, floating point numbers will be
      used. Floating point arithmetic is faster but might give the
      wrong result for degenerate input.

    - ``ambient_dim`` -- integer. The ambient space dimension. Usually
      can be figured out automatically from the H/Vrepresentation
      dimensions.

    - ``backend`` -- string or ``None`` (default). The backend to use. Valid choices are

      * ``'cdd'``: use cdd
        (:mod:`~sage.geometry.polyhedron.backend_cdd`) with `\QQ` or
        `\RDF` coefficients depending on ``base_ring``.


      * ``'ppl'``: use ppl
        (:mod:`~sage.geometry.polyhedron.backend_ppl`) with `\ZZ` or
        `\QQ` coefficients depending on ``base_ring``.

    Some backends support further optional arguments:

    - ``minimize`` -- boolean (default: ``True``). Whether to
      immediately remove redundant H/V-representation data. Currently
      not used.

    - ``verbose`` -- boolean (default: ``False``). Whether to print
      verbose output for debugging purposes. Only supported by the cdd
      backends.

    OUTPUT:

    The polyhedron defined by the input data.

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
        [An inequality (1, 0) x + 1 >= 0,
         An inequality (0, 1) x + 1 >= 0,
         An inequality (-1, 0) x + 1 >= 0,
         An inequality (0, -1) x + 1 >= 0]
        sage: p = Polyhedron(vertices = [[1.1, 2.2], [3.3, 4.4]], base_ring=RDF)
        sage: p.n_inequalities()
        2

    The same polyhedron given in two ways::

        sage: p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0]])
        sage: p.Vrepresentation()
        (A line in the direction (0, 0, 1),
         A ray in the direction (1, 0, 0),
         A ray in the direction (0, 1, 0),
         A vertex at (0, 0, 0))
        sage: q = Polyhedron(vertices=[[0,0,0]], rays=[[1,0,0],[0,1,0]], lines=[[0,0,1]])
        sage: q.Hrepresentation()
        (An inequality (1, 0, 0) x + 0 >= 0,
         An inequality (0, 1, 0) x + 0 >= 0)

    Finally, a more complicated example. Take `\mathbb{R}_{\geq 0}^6` with
    coordinates `a, b, \dots, f` and

      * The inequality `e+b \geq c+d`
      * The inequality `e+c \geq b+d`
      * The equation `a+b+c+d+e+f = 31`

    ::

        sage: positive_coords = Polyhedron(ieqs=[
        ...       [0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0],
        ...       [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1]])
        sage: P = Polyhedron(ieqs=positive_coords.inequalities() + (
        ...       [0,0,1,-1,-1,1,0], [0,0,-1,1,-1,1,0]), eqns=[[-31,1,1,1,1,1,1]])
        sage: P
        A 5-dimensional polyhedron in QQ^6 defined as the convex hull of 7 vertices
        sage: P.dim()
        5
        sage: P.Vrepresentation()
        (A vertex at (31, 0, 0, 0, 0, 0), A vertex at (0, 0, 0, 0, 0, 31),
         A vertex at (0, 0, 0, 0, 31, 0), A vertex at (0, 0, 31/2, 0, 31/2, 0),
         A vertex at (0, 31/2, 31/2, 0, 0, 0), A vertex at (0, 31/2, 0, 0, 31/2, 0),
         A vertex at (0, 0, 0, 31/2, 31/2, 0))

    .. NOTE::

      * Once constructed, a ``Polyhedron`` object is immutable.
      * Although the option ``field=RDF`` allows numerical data to
        be used, it might not give the right answer for degenerate
        input data - the results can depend upon the tolerance
        setting of cdd.
    """
    # Clean up the arguments
    vertices = _make_listlist(vertices)
    rays     = _make_listlist(rays)
    lines    = _make_listlist(lines)
    ieqs     = _make_listlist(ieqs)
    eqns     = _make_listlist(eqns)

    got_Vrep = (len(vertices+rays+lines) > 0)
    got_Hrep = (len(ieqs+eqns) > 0)

    if got_Vrep and got_Hrep:
        raise ValueError('You cannot specify both H- and V-representation.')
    elif got_Vrep:
        deduced_ambient_dim = _common_length_of(vertices, rays, lines)[1]
    elif got_Hrep:
        deduced_ambient_dim = _common_length_of(ieqs, eqns)[1] - 1
    else:
        if ambient_dim is None:
            deduced_ambient_dim = 0
        else:
            deduced_ambient_dim = ambient_dim
        if base_ring is None:
            base_ring = ZZ

    # set ambient_dim
    if ambient_dim is not None and deduced_ambient_dim!=ambient_dim:
        raise ValueError('Ambient space dimension mismatch. Try removing the "ambient_dim" parameter.')
    ambient_dim = deduced_ambient_dim

    # figure out base_ring
    from sage.misc.flatten import flatten
    values = flatten(vertices+rays+lines+ieqs+eqns)
    if base_ring is not None:
        try:
            convert = not all(x.parent() is base_ring for x in values)
        except AttributeError:   # No x.parent() method?
            convert = True
    else:
        from sage.rings.integer import is_Integer
        from sage.rings.rational import is_Rational
        from sage.rings.real_double import is_RealDoubleElement
        if all(is_Integer(x) for x in values):
            if got_Vrep:
                base_ring = ZZ
            else:   # integral inequalities usually do not determine a latice polytope!
                base_ring = QQ
            convert=False
        elif all(is_Rational(x) for x in values):
            base_ring = QQ
            convert=False
        elif all(is_RealDoubleElement(x) for x in values):
            base_ring = RDF
            convert=False
        else:
            try:
                map(ZZ, values)
                if got_Vrep:
                    base_ring = ZZ
                else:
                    base_ring = QQ
                convert = True
            except TypeError:
                from sage.structure.sequence import Sequence
                values = Sequence(values)
                if QQ.has_coerce_map_from(values.universe()):
                    base_ring = QQ
                    convert = True
                else:
                    base_ring = RDF
                    convert = True

    # Add the origin if necesarry
    if got_Vrep and len(vertices)==0:
        vertices = [ [0]*ambient_dim ]

    # Specific backends can override the base_ring
    from sage.geometry.polyhedron.parent import Polyhedra
    parent = Polyhedra(base_ring, ambient_dim, backend=backend)
    base_ring = parent.base_ring()

    # Convert into base_ring if necessary
    def convert_base_ring(lstlst):
        return [ [base_ring(x) for x in lst] for lst in lstlst]
    Hrep = Vrep = None
    if got_Hrep:
        if convert:
            Hrep = [convert_base_ring(ieqs), convert_base_ring(eqns)]
        else:
            Hrep = [ieqs, eqns]
    if got_Vrep:
        if convert:
            Vrep = [convert_base_ring(vertices), convert_base_ring(rays), convert_base_ring(lines)]
        else:
            Vrep = [vertices, rays, lines]

    # finally, construct the Polyhedron
    return parent.element_class(parent, Vrep, Hrep)
