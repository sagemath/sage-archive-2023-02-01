r"""
Fano toric varieties

This module provides support for (Crepant Partial Resolutions of) Fano toric
varieties, corresponding to crepant subdivisions of face fans of reflexive
:class:`lattice polytopes
<sage.geometry.lattice_polytope.LatticePolytopeClass>`.
The interface is provided via :func:`CPRFanoToricVariety`.

A careful exposition of different flavours of Fano varieties can be found in
the paper by Benjamin Nill [Nill2005]_. The main goal of this module is to
support work with **Gorenstein weak Fano toric varieties**. Such a variety
corresponds to a **coherent crepant refinement of the normal fan of a
reflexive polytope** `\Delta`, where crepant means that primitive generators
of the refining rays lie on the facets of the polar polytope `\Delta^\circ`
and coherent (a.k.a. regular or projective) means that there exists a strictly
upper convex piecewise linear function whose domains of linearity are
precisely the maximal cones of the subdivision. These varieties are important
for string theory in physics, as they serve as ambient spaces for mirror pairs
of Calabi-Yau manifolds via constructions due to Victor V. Batyrev
[Batyrev1994]_ and Lev A. Borisov [Borisov1993]_.

From the combinatorial point of view "crepant" requirement is much more simple
and natural to work with than "coherent." For this reason, the code in this
module will allow work with arbitrary crepant subdivisions without checking
whether they are coherent or not. We refer to corresponding toric varieties as
**CPR-Fano toric varieties**.

REFERENCES:

..  [Batyrev1994]
    Victor V. Batyrev,
    "Dual polyhedra and mirror symmetry for Calabi-Yau hypersurfaces in toric
    varieties",
    J. Algebraic Geom. 3 (1994), no. 3, 493-535.
    arXiv:alg-geom/9310003v1

..  [Borisov1993]
    Lev A. Borisov,
    "Towards the mirror symmetry for Calabi-Yau complete intersections in
    Gorenstein Fano toric varieties", 1993.
    arXiv:alg-geom/9310001v1

..  [CD2007]
    Adrian Clingher and Charles F. Doran,
    "Modular invariants for lattice polarized K3 surfaces",
    Michigan Math. J. 55 (2007), no. 2, 355-393.
    arXiv:math/0602146v1 [math.AG]

..  [Nill2005]
    Benjamin Nill,
    "Gorenstein toric Fano varieties",
    Manuscripta Math. 116 (2005), no. 2, 183-210.
    arXiv:math/0405448v1 [math.AG]

AUTHORS:

- Andrey Novoseltsev (2010-05-18): initial version.

EXAMPLES:

Most of the functions available for Fano toric varieties are the same as
for general toric varieties, so here we will concentrate only on
Calabi-Yau subvarieties, which were the primary goal for creating this
module.

For our first example we realize the projective plane as a Fano toric
variety::

    sage: simplex = lattice_polytope.projective_space(2)
    sage: P2 = CPRFanoToricVariety(Delta_polar=simplex)

Its anticanonical "hypersurface" is a one-dimensional Calabi-Yau
manifold::

    sage: P2.anticanonical_hypersurface(
    ...         monomial_points="all")
    Closed subscheme of 2-d CPR-Fano toric variety
    covered by 3 affine patches defined by:
      a0*z0^3 + a9*z0^2*z1 + a7*z0*z1^2
    + a1*z1^3 + a8*z0^2*z2 + a6*z0*z1*z2
    + a4*z1^2*z2 + a5*z0*z2^2
    + a3*z1*z2^2 + a2*z2^3

In many cases it is sufficient to work with the "simplified polynomial
moduli space" of anticanonical hypersurfaces::

    sage: P2.anticanonical_hypersurface(
    ...         monomial_points="simplified")
    Closed subscheme of 2-d CPR-Fano toric variety
    covered by 3 affine patches defined by:
      a0*z0^3 + a1*z1^3 + a6*z0*z1*z2 + a2*z2^3

The mirror family to these hypersurfaces lives inside the Fano toric
variety obtained using ``simplex`` as ``Delta`` instead of ``Delta_polar``::

    sage: FTV = CPRFanoToricVariety(Delta=simplex,
    ...         coordinate_points="all")
    sage: FTV.anticanonical_hypersurface(
    ...         monomial_points="simplified")
    Closed subscheme of 2-d CPR-Fano toric variety
    covered by 9 affine patches defined by:
      a2*z2^3*z3^2*z4*z5^2*z8
    + a1*z1^3*z3*z4^2*z7^2*z9
    + a3*z0*z1*z2*z3*z4*z5*z7*z8*z9
    + a0*z0^3*z5*z7*z8^2*z9^2

Here we have taken the resolved version of the ambient space for the
mirror family, but in fact we don't have to resolve singularities
corresponding to the interior points of facets - they are singular
points which do not lie on a generic anticanonical hypersurface::

    sage: FTV = CPRFanoToricVariety(Delta=simplex,
    ...         coordinate_points="all but facets")
    sage: FTV.anticanonical_hypersurface(
    ...         monomial_points="simplified")
    Closed subscheme of 2-d CPR-Fano toric variety
    covered by 3 affine patches defined by:
      a0*z0^3 + a1*z1^3 + a3*z0*z1*z2 + a2*z2^3

This looks very similar to our second version of the anticanonical
hypersurface of the projective plane, as expected, since all
one-dimensional Calabi-Yau manifolds are elliptic curves!

Now let's take a look at a toric realization of `M`-polarized K3 surfaces
studied by Adrian Clingher and Charles F. Doran in [CD2007]_::

    sage: p4318 = ReflexivePolytope(3, 4318)  # long
    sage: FTV = CPRFanoToricVariety(Delta_polar=p4318)  # long
    sage: FTV.anticanonical_hypersurface()  # long
    Closed subscheme of 3-d CPR-Fano toric variety
    covered by 4 affine patches defined by:
      a3*z2^12 + a4*z2^6*z3^6 + a2*z3^12
    + a8*z0*z1*z2*z3 + a0*z1^3 + a1*z0^2

Below you will find detailed descriptions of available functions. Current
functionality of this module is very basic, but it is under active
development and hopefully will improve in future releases of Sage. If there
are some particular features that you would like to see implemented ASAP,
please consider reporting them to the Sage Development Team or even
implementing them on your own as a patch for inclusion!
"""
# The first example of the tutorial is taken from
# CPRFanoToricVariety_field.anticanonical_hypersurface


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.geometry.all import Cone, FaceFan, Fan
from sage.misc.all import latex, prod
from sage.rings.all import (PolynomialRing, QQ,
                            is_FractionField, is_Field,
                            is_MPolynomialRing, is_PolynomialRing)
from sage.schemes.generic.toric_variety import (
                                            AlgebraicScheme_subscheme_toric,
                                            ToricVariety_field,
                                            normalize_names)


# Default coefficient for anticanonical hypersurfaces
DEFAULT_COEFFICIENT = "a"


def is_CPRFanoToricVariety(x):
    r"""
    Check if ``x`` is a CPR-Fano toric variety.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a :class:`CPR-Fano toric variety
      <CPRFanoToricVariety_field>` and ``False`` otherwise.

    .. NOTE::

        While projective spaces are Fano toric varieties mathematically, they
        are not toric varieties in Sage due to efficiency considerations, so
        this function will return ``False``.

    EXAMPLES::

        sage: from sage.schemes.generic.fano_toric_variety import (
        ...     is_CPRFanoToricVariety)
        sage: is_CPRFanoToricVariety(1)
        False
        sage: FTV = CPRFanoToricVariety(lattice_polytope.octahedron(2))
        sage: FTV
        2-d CPR-Fano toric variety covered by 4 affine patches
        sage: is_CPRFanoToricVariety(FTV)
        True
        sage: is_CPRFanoToricVariety(ProjectiveSpace(2))
        False
    """
    return isinstance(x, CPRFanoToricVariety_field)


def CPRFanoToricVariety(Delta=None,
                        Delta_polar=None,
                        coordinate_points=None,
                        charts=None,
                        coordinate_names=None,
                        coordinate_name_indices=None,
                        make_simplicial=False,
                        base_field=None,
                        check=True):
    r"""
    Construct a CPR-Fano toric variety.

    .. NOTE::

        See documentation of the module
        :mod:`~sage.schemes.generic.fano_toric_variety` for the used
        definitions and supported varieties.

    Due to the large number of available options, it is recommended to always
    use keyword parameters.

    INPUT:

    - ``Delta`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *normal fan* of ``Delta``. Either ``Delta`` or ``Delta_polar`` must be
      given, but not both at the same time, since one is completely determined
      by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method;

    - ``Delta_polar`` -- reflexive :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The fan of the
      constructed CPR-Fano toric variety will be a crepant subdivision of the
      *face fan* of ``Delta_polar``. Either ``Delta`` or ``Delta_polar`` must
      be given, but not both at the same time, since one is completely
      determined by another via :meth:`polar
      <sage.geometry.lattice_polytope.LatticePolytopeClass.polar>` method;

    - ``coordinate_points`` -- list of integers or string. A list will be
      interpreted as indices of (boundary) points of ``Delta_polar`` which
      should be used as rays of the underlying fan. It must include all
      vertices of ``Delta_polar`` and no repetitions are allowed. A string
      must be one of the following descriptions of points of ``Delta_polar``:

      * "vertices" (default),
      * "all" (will not include the origin),
      * "all but facets" (will not include points in the relative interior of
        facets);

    - ``charts`` -- list of lists of elements from ``coordinate_points``. Each
      of these lists must define a generating cone of a fan subdividing the
      normal fan of ``Delta``. Default ``charts`` correspond to the normal fan
      of ``Delta`` without subdivision. The fan specified by ``charts`` will
      be subdivided to include all of the requested ``coordinate_points``;

    - ``coordinate_names`` -- names of variables for the coordinate ring, see
      :func:`~sage.schemes.generic.toric_variety.normalize_names`
      for acceptable formats. If not given, indexed variable names will be
      created automatically;

    - ``coordinate_name_indices`` -- list of integers, indices for indexed
      variables. If not given, the index of each variable will coincide with
      the index of the corresponding point of ``Delta_polar``;

    - ``make_simplicial`` -- if ``True``, the underlying fan will be made
      simplicial (default: ``False``);

    - ``base_field`` -- base field of the CPR-Fano toric variety
      (default: `\QQ`);

    - ``check`` -- by default the input data will be checked for correctness
      (e.g. that ``charts`` do form a subdivision of the normal fan of
      ``Delta``). If you know for sure that the input is valid, you may
      significantly decrease construction time using ``check=False`` option.

    OUTPUT:

    - :class:`CPR-Fano toric variety <CPRFanoToricVariety_field>`.

    EXAMPLES:

    We start with the product of two projective lines::

        sage: diamond = lattice_polytope.octahedron(2)
        sage: diamond.vertices()
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        sage: P1xP1 = CPRFanoToricVariety(Delta_polar=diamond)
        sage: P1xP1
        2-d CPR-Fano toric variety covered by 4 affine patches

    "Unfortunately," this variety is smooth to start with and we cannot
    perform any subdivisions of the underlying fan without leaving the
    category of CPR-Fano toric varieties. Our next example starts with a
    square::

        sage: square = diamond.polar()
        sage: square.vertices()
        [-1  1 -1  1]
        [ 1  1 -1 -1]
        sage: square.points()
        [-1  1 -1  1 -1  0  0  0  1]
        [ 1  1 -1 -1  0 -1  0  1  0]

    We will construct several varieties associated to it::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square)
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1]
        [ 1  1 -1 -1]
        sage: FTV.gens()
        (z0, z1, z2, z3)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,8])
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1  1]
        [ 1  1 -1 -1  0]
        sage: FTV.gens()
        (z0, z1, z2, z3, z8)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[8,0,2,1,3],
        ...         coordinate_names="x+")
        sage: FTV.fan().ray_matrix()
        [ 1 -1 -1  1  1]
        [ 0  1 -1  1 -1]
        sage: FTV.gens()
        (x8, x0, x2, x1, x3)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points="all",
        ...         coordinate_names="x y Z+")
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1 -1  0  0  1]
        [ 1  1 -1 -1  0 -1  1  0]
        sage: FTV.gens()
        (x, y, Z2, Z3, Z4, Z5, Z7, Z8)

    Note that ``Z6`` is "missing". This is due to the fact that the 6-th point
    of ``square`` is the origin, and all automatically created names have the
    same indices as corresponding points of
    :meth:`~CPRFanoToricVariety_field.Delta_polar`. This is usually very
    convenient, especially if you have to work with several partial
    resolutions of the same Fano toric variety. However, you can change it, if
    you want::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points="all",
        ...         coordinate_names="x y Z+",
        ...         coordinate_name_indices=range(8))
        sage: FTV.gens()
        (x, y, Z2, Z3, Z4, Z5, Z6, Z7)

    Note that you have to provide indices for *all* variables, including those
    that have "completely custom" names. Again, this is usually convenient,
    because you can add or remove "custom" variables without disturbing too
    much "automatic" ones::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points="all",
        ...         coordinate_names="x Z+",
        ...         coordinate_name_indices=range(8))
        sage: FTV.gens()
        (x, Z1, Z2, Z3, Z4, Z5, Z6, Z7)

    If you prefer to always start from zero, you will have to shift indices
    accordingly::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points="all",
        ...         coordinate_names="x Z+",
        ...         coordinate_name_indices=[0] + range(7))
        sage: FTV.gens()
        (x, Z0, Z1, Z2, Z3, Z4, Z5, Z6)

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points="all",
        ...         coordinate_names="x y Z+",
        ...         coordinate_name_indices=[0]*2 + range(6))
        sage: FTV.gens()
        (x, y, Z0, Z1, Z2, Z3, Z4, Z5)

    So you always can get any names you want, somewhat complicated default
    behaviour was designed with the hope that in most cases you will have no
    desire to provide different names.

    Now we will use the possibility to specify initial charts::

        sage: charts = [(0,1), (1,3), (3,2), (2,0)]

    (these charts actually form exactly the face fan of our square) ::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,4],
        ...         charts=charts)
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1 -1]
        [ 1  1 -1 -1  0]
        sage: [cone.ambient_ray_indices() for cone in FTV.fan()]
        [(0, 1), (1, 3), (2, 3), (2, 4), (0, 4)]

    If charts are wrong, it should be detected::

        sage: bad_charts = charts + [(2,0)]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,4],
        ...         charts=bad_charts)
        d...: UserWarning: you have provided a non-minimal set
        of generating cones, 1 of them were discarded!
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1 -1]
        [ 1  1 -1 -1  0]
        sage: [cone.ambient_ray_indices() for cone in FTV.fan()]
        [(0, 1), (1, 3), (2, 3), (2, 4), (0, 4)]

    OK, everything still works but that's because these charts still define a
    valid fan, they just happened to list one cone twice and it was discarded
    by the fan constructor. However, when you try to speed up your code such
    charts are indeed bad::

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,4],
        ...         charts=bad_charts,
        ...         check=False)
        sage: FTV.fan().ray_matrix()
        [-1  1 -1  1 -1]
        [ 1  1 -1 -1  0]
        sage: [cone.ambient_ray_indices() for cone in FTV.fan()]
        [(0, 1), (1, 3), (2, 3), (2, 4), (0, 4), (2, 4), (0, 4)]

    The last line shows two of the generating cones twice. While "everything
    still works" in the sense "it does not crash," any work with such a
    variety may lead to mathematically wrong results, so use ``check=False``
    carefully!

    Let's also look at some of the mistakes which are reported, rather than
    fixed::

        sage: bad_charts = charts + [(0,3)]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,4],
        ...         charts=bad_charts)
        Traceback (most recent call last):
        ...
        ValueError: (0, 3) does not form a chart of a
        subdivision of the face fan of A polytope polar
        to An octahedron: 2-dimensional, 4 vertices.!

        sage: bad_charts = charts[:-1]
        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,4],
        ...         charts=bad_charts)
        Traceback (most recent call last):
        ...
        ValueError: given charts do not form a complete fan!

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: all 4 vertices of Delta_polar
        must be used for coordinates!
        Got: [1, 2, 3, 4]

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,0,1,2,3,4])
        Traceback (most recent call last):
        ...
        ValueError: no repetitions are
        allowed for coordinate points!
        Got: [0, 0, 1, 2, 3, 4]

        sage: FTV = CPRFanoToricVariety(Delta_polar=square,
        ...         coordinate_points=[0,1,2,3,6])
        Traceback (most recent call last):
        ...
        ValueError: the origin (point #6)
        cannot be used for a coordinate!
        Got: [0, 1, 2, 3, 6]
    """
    # Check/normalize Delta_polar
    if Delta is None and Delta_polar is None:
        raise ValueError("either Delta or Delta_polar must be given!")
    elif Delta is not None and Delta_polar is not None:
        raise ValueError("Delta and Delta_polar cannot be given together!")
    elif Delta_polar is None:
        Delta_polar = Delta.polar()
    elif not Delta_polar.is_reflexive():
        raise ValueError("Delta_polar must be reflexive!")
    # Check/normalize coordinate_points and construct fan rays
    if coordinate_points is None:
        coordinate_points = range(Delta_polar.nvertices())
        if charts is not None:
            for chart in charts:
                for point in chart:
                    if point not in coordinate_points:
                        coordinate_points.append(point)
    elif coordinate_points == "vertices":
        coordinate_points = range(Delta_polar.nvertices())
    elif coordinate_points == "all":
        coordinate_points = range(Delta_polar.npoints())
        coordinate_points.remove(Delta_polar.origin())
    elif coordinate_points == "all but facets":
        coordinate_points = Delta_polar.skeleton_points(Delta_polar.dim() - 2)
    elif isinstance(coordinate_points, str):
        raise ValueError("unrecognized description of the coordinate points!"
                         "\nGot: %s" % coordinate_points)
    elif check:
        cp_set = set(coordinate_points)
        if len(cp_set) != len(coordinate_points):
            raise ValueError(
                "no repetitions are allowed for coordinate points!\nGot: %s"
                % coordinate_points)
        if not cp_set.issuperset(range(Delta_polar.nvertices())):
            raise ValueError("all %d vertices of Delta_polar must be used "
                "for coordinates!\nGot: %s"
                % (Delta_polar.nvertices(), coordinate_points))
        if Delta_polar.origin() in cp_set:
            raise ValueError("the origin (point #%d) cannot be used for a "
                "coordinate!\nGot: %s"
                % (Delta_polar.origin(), coordinate_points))
    point_to_ray = dict()
    for n, point in enumerate(coordinate_points):
        point_to_ray[point] = n
    # This can be simplified if LatticePolytopeClass is adjusted.
    rays = [Delta_polar.point(p) for p in coordinate_points]
    for ray in rays:
        ray.set_immutable()
    # Check/normalize charts and construct the fan based on them.
    if charts is None:
        # Start with the face fan
        fan = FaceFan(Delta_polar)
    else:
        # First of all, check that each chart is completely contained in a
        # single facet of Delta_polar, otherwise they do not form a
        # subdivision of the face fan of Delta_polar
        if check:
            facet_sets = [frozenset(facet.points())
                          for facet in Delta_polar.facets()]
            for chart in charts:
                is_bad = True
                for fset in facet_sets:
                    if fset.issuperset(chart):
                        is_bad = False
                        break
                if is_bad:
                    raise ValueError(
                        "%s does not form a chart of a subdivision of the "
                        "face fan of %s!" % (chart, Delta_polar))
        # We will construct the initial fan from Cone objects: since charts
        # may not use all of the necessary rays, alternative form is tedious
        # With check=False it should not be long anyway.
        cones = [Cone((rays[point_to_ray[point]] for point in chart),
                      check=check)
                 for chart in charts]
        fan = Fan(cones, check=check)
        if check and not fan.is_complete():
            raise ValueError("given charts do not form a complete fan!")
    # Subdivide this fan to use all required points
    fan = fan.subdivide(new_rays=(ray for ray in rays
                                      if ray not in fan.ray_set()),
                        make_simplicial=make_simplicial)
    # Now create yet another fan making sure that the order of the rays is
    # the same as requested (it is a bit difficult to get it from the start)
    trans = dict()
    for n, ray in enumerate(fan.rays()):
        trans[n] = rays.index(ray)
    cones = tuple(tuple(sorted(trans[r] for r in cone.ambient_ray_indices()))
                  for cone in fan)
    fan = Fan(cones, rays, check=False, normalize=False)
    # Check/normalize base_field
    if base_field is None:
        base_field = QQ
    elif not is_Field(base_field):
        raise TypeError("need a field to construct a Fano toric variety!"
                        "\n Got %s" % base_field)
    fan._is_complete = True     # At this point it must be for sure
    return CPRFanoToricVariety_field(Delta_polar, fan, coordinate_points,
        point_to_ray, coordinate_names, coordinate_name_indices, base_field)


class CPRFanoToricVariety_field(ToricVariety_field):
    r"""
    Construct a CPR-Fano toric variety associated to a reflexive polytope.

    .. WARNING::

        This class does not perform any checks of correctness of input and it
        does assume that the internal structure of the given parameters is
        coordinated in a certain way. Use
        :func:`CPRFanoToricVariety` to construct CPR-Fano toric varieties.

    .. NOTE::

        See documentation of the module
        :mod:`~sage.schemes.generic.fano_toric_variety` for the used
        definitions and supported varieties.

    INPUT:

    - ``Delta_polar`` -- reflexive polytope;

    - ``fan`` -- rational polyhedral fan subdividing the face fan of
      ``Delta_polar``;

    - ``coordinate_points`` -- list of indices of points of ``Delta_polar``
      used for rays of ``fan``;

    - ``point_to_ray`` -- dictionary mapping the index of a coordinate point
      to the index of the corresponding ray;

    - ``coordinate_names`` -- names of the variables of the coordinate ring in
      the format accepted by
      :func:`~sage.schemes.generic.toric_variety.normalize_names`;

    - ``coordinate_name_indices`` -- indices for indexed variables,
      if ``None``, will be equal to ``coordinate_points``;

    - ``base_field`` -- base field of the CPR-Fano toric variety.

    OUTPUT:

    - :class:`CPR-Fano toric variety <CPRFanoToricVariety_field>`.

    TESTS::

        sage: P1xP1 = CPRFanoToricVariety(
        ...       Delta_polar=lattice_polytope.octahedron(2))
        sage: P1xP1
        2-d CPR-Fano toric variety covered by 4 affine patches
    """

    def __init__(self, Delta_polar, fan, coordinate_points, point_to_ray,
                 coordinate_names, coordinate_name_indices, base_field):
        r"""
        See :class:`CPRFanoToricVariety_field` for documentation.

        Use ``CPRFanoToricVariety`` to construct CPR-Fano toric varieties.

        TESTS::

            sage: P1xP1 = CPRFanoToricVariety(
            ...       Delta_polar=lattice_polytope.octahedron(2))
            sage: P1xP1
            2-d CPR-Fano toric variety covered by 4 affine patches
        """
        self._Delta_polar = Delta_polar
        self._coordinate_points = tuple(coordinate_points)
        self._point_to_ray = point_to_ray
        # Check/normalize coordinate_indices
        if coordinate_name_indices is None:
            coordinate_name_indices = coordinate_points
        super(CPRFanoToricVariety_field, self).__init__(fan, coordinate_names,
                                        coordinate_name_indices, base_field)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: P1xP1 = CPRFanoToricVariety(
            ...       Delta_polar=lattice_polytope.octahedron(2))
            sage: P1xP1._latex_()
            '\\mathbb{P}_{\\Delta^{2}}'
        """
        return r"\mathbb{P}_{%s}" % latex(self.Delta())

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: P1xP1 = CPRFanoToricVariety(
            ...       Delta_polar=lattice_polytope.octahedron(2))
            sage: P1xP1._repr_()
            '2-d CPR-Fano toric variety covered by 4 affine patches'
        """
        return ("%d-d CPR-Fano toric variety covered by %d affine patches"
                % (self.dimension_relative(), self.fan().ngenerating_cones()))

    def anticanonical_hypersurface(self, **kwds):
        r"""
        Return the anticanonical hypersurface of ``self``.

        .. NOTE::

            The returned hypersurface may be actually a subscheme of
            **another** CPR-Fano toric variety: if the base field of ``self``
            does not include all of the required names for generic monomial
            coefficients, it will be automatically extended.

        Below `\Delta` is the reflexive polytope corresponding to ``self``,
        i.e. the fan of ``self`` is a refinement of the normal fan of
        `\Delta`. This function accepts only keyword parameters.

        INPUT:

        - ``monomial points`` -- list of integers or string. A list will be
          interpreted as indices of points of `\Delta` which should be used
          for monomials of this hypersurface. A string must be one of the
          following descriptions of points of `\Delta`:

          * "vertices",
          * "vertices+origin",
          * "all",
          * "simplified" (default) -- all points of `\Delta` except for
            the interior points of facets, this choice corresponds to working
            with the "simplified polynomial moduli space" of anticanonical
            hypersurfaces;

        - ``coefficient_names`` -- names for the monomial coefficients, see
          :func:`~sage.schemes.generic.toric_variety.normalize_names`
          for acceptable formats. If not given, indexed coefficient names will
          be created automatically;

        - ``coefficient_name_indices`` -- list of integers, indices for
          indexed coefficients. If not given, the index of each coefficient
          will coincide with the index of the corresponding point of `\Delta`.

        OUTPUT:

        - :class:`anticanonical hypersurface <AnticanonicalHypersurface>` of
          ``self`` (with the extended base field, if necessary).

        EXAMPLES:

        We realize the projective plane as a Fano toric variety::

            sage: simplex = lattice_polytope.projective_space(2)
            sage: P2 = CPRFanoToricVariety(Delta_polar=simplex)

        Its anticanonical "hypersurface" is a one-dimensional Calabi-Yau
        manifold::

            sage: P2.anticanonical_hypersurface(
            ...         monomial_points="all")
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 3 affine patches defined by:
              a0*z0^3 + a9*z0^2*z1 + a7*z0*z1^2
            + a1*z1^3 + a8*z0^2*z2 + a6*z0*z1*z2
            + a4*z1^2*z2 + a5*z0*z2^2
            + a3*z1*z2^2 + a2*z2^3

        In many cases it is sufficient to work with the "simplified polynomial
        moduli space" of anticanonical hypersurfaces::

            sage: P2.anticanonical_hypersurface(
            ...         monomial_points="simplified")
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 3 affine patches defined by:
              a0*z0^3 + a1*z1^3 + a6*z0*z1*z2 + a2*z2^3

        The mirror family to these hypersurfaces lives inside the Fano toric
        variety obtained using ``simplex`` as ``Delta`` instead of
        ``Delta_polar``::

            sage: FTV = CPRFanoToricVariety(Delta=simplex,
            ...         coordinate_points="all")
            sage: FTV.anticanonical_hypersurface(
            ...         monomial_points="simplified")
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 9 affine patches defined by:
              a2*z2^3*z3^2*z4*z5^2*z8
            + a1*z1^3*z3*z4^2*z7^2*z9
            + a3*z0*z1*z2*z3*z4*z5*z7*z8*z9
            + a0*z0^3*z5*z7*z8^2*z9^2

        Here we have taken the resolved version of the ambient space for the
        mirror family, but in fact we don't have to resolve singularities
        corresponding to the interior points of facets - they are singular
        points which do not lie on a generic anticanonical hypersurface::

            sage: FTV = CPRFanoToricVariety(Delta=simplex,
            ...         coordinate_points="all but facets")
            sage: FTV.anticanonical_hypersurface(
            ...         monomial_points="simplified")
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 3 affine patches defined by:
              a0*z0^3 + a1*z1^3 + a3*z0*z1*z2 + a2*z2^3

        This looks very similar to our second version of the anticanonical
        hypersurface of the projective plane, as expected, since all
        one-dimensional Calabi-Yau manifolds are elliptic curves!
        """
        # The example above is also copied to the tutorial section in the
        # main documentation of the module.
        return AnticanonicalHypersurface(self, **kwds)

    def change_ring(self, F):
        r"""
        Return a CPR-Fano toric variety over field ``F``, otherwise the same
        as ``self``.

        INPUT:

        - ``F`` -- field.

        OUTPUT:

        - :class:`CPR-Fano toric variety <CPRFanoToricVariety_field>` over
          ``F``.

        .. NOTE::

            There is no need to have any relation between ``F`` and the base
            field of ``self``. If you do want to have such a relation, use
            :meth:`base_extend` instead.

        EXAMPLES::

            sage: P1xP1 = CPRFanoToricVariety(
            ...       Delta_polar=lattice_polytope.octahedron(2))
            sage: P1xP1.base_ring()
            Rational Field
            sage: P1xP1_RR = P1xP1.change_ring(RR)
            sage: P1xP1_RR.base_ring()
            Real Field with 53 bits of precision
            sage: P1xP1_QQ = P1xP1_RR.change_ring(QQ)
            sage: P1xP1_QQ.base_ring()
            Rational Field
            sage: P1xP1_RR.base_extend(QQ)
            Traceback (most recent call last):
            ...
            ValueError: no natural map from the base ring
            (=Real Field with 53 bits of precision)
            to R (=Rational Field)!
        """
        if self.base_ring() == F:
            return self
        else:
            return CPRFanoToricVariety_field(self._Delta_polar, self._fan,
                self._coordinate_points, self._point_to_ray,
                self.variable_names(), None, F)
                # coordinate_name_indices do not matter, we give explicit
                # names for all variables

    def coordinate_point_to_coordinate(self, point):
        r"""
        Return the variable of the coordinate ring corresponding to ``point``.

        INPUT:

        - ``point`` -- integer from the list of :meth:`coordinate_points`.

        OUTPUT:

        - the corresponding generator of the coordinate ring of ``self``.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: FTV = CPRFanoToricVariety(diamond,
            ...         coordinate_points=[0,1,2,3,8])
            sage: FTV.coordinate_points()
            (0, 1, 2, 3, 8)
            sage: FTV.gens()
            (z0, z1, z2, z3, z8)
            sage: FTV.coordinate_point_to_coordinate(8)
            z8
        """
        return self.gen(self._point_to_ray[point])

    def coordinate_points(self):
        r"""
        Return indices of points of :meth:`Delta_polar` used for coordinates.

        OUTPUT:

        - :class:`tuple` of integers.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: square = diamond.polar()
            sage: FTV = CPRFanoToricVariety(Delta_polar=square,
            ...         coordinate_points=[0,1,2,3,8])
            sage: FTV.coordinate_points()
            (0, 1, 2, 3, 8)
            sage: FTV.gens()
            (z0, z1, z2, z3, z8)

            sage: FTV = CPRFanoToricVariety(Delta_polar=square,
            ...         coordinate_points="all")
            sage: FTV.coordinate_points()
            (0, 1, 2, 3, 4, 5, 7, 8)
            sage: FTV.gens()
            (z0, z1, z2, z3, z4, z5, z7, z8)

        Note that one point is missing, namely ::

            sage: square.origin()
            6
        """
        return self._coordinate_points

    def Delta(self):
        r"""
        Return the reflexive polytope associated to ``self``.

        OUTPUT:

        - reflexive :class:`lattice polytope
          <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The
          underlying fan of ``self`` is a coherent subdivision of the
          *normal fan* of this polytope.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: P1xP1 = CPRFanoToricVariety(Delta_polar=diamond)
            sage: P1xP1.Delta()
            A polytope polar to An octahedron: 2-dimensional, 4 vertices.
            sage: P1xP1.Delta() is diamond.polar()
            True
        """
        return self._Delta_polar.polar()

    def Delta_polar(self):
        r"""
        Return polar of :meth:`Delta`.

        OUTPUT:

        - reflexive :class:`lattice polytope
          <sage.geometry.lattice_polytope.LatticePolytopeClass>`. The
          underlying fan of ``self`` is a coherent subdivision of the
          *face fan* of this polytope.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: P1xP1 = CPRFanoToricVariety(Delta_polar=diamond)
            sage: P1xP1.Delta_polar()
            An octahedron: 2-dimensional, 4 vertices.
            sage: P1xP1.Delta_polar() is diamond
            True
            sage: P1xP1.Delta_polar() is P1xP1.Delta().polar()
            True
        """
        return self._Delta_polar

    def plot3d(self):
        r"""
        Return the 3D plot of the underlying fan of ``self``.

        Since this fan corresponds to a triangulation of the associated
        reflexive polytope, this function actually plots this polytope with
        the triangulation.

        OUTPUT:

        - plot.

        .. NOTE::

            This function is currently present to allow quick control over
            subdivisions. It will be shortly replaced with a superior one
            allowing great customization of the output.

        EXAMPLES::

            sage: P1xP1xP1 = CPRFanoToricVariety(
            ...          Delta_polar=lattice_polytope.octahedron(3))
            sage: P1xP1xP1.plot3d()

        A more subdivided one::

            sage: FTV = CPRFanoToricVariety(
            ...         Delta=lattice_polytope.octahedron(3),
            ...         coordinate_points="all but facets")
            sage: FTV.plot3d()
        """
        # This function should be greatly improved! Probably after
        # plotting capabilities of lattice polytopes are switched to
        # dictionaries to allow uniformity.
        #
        # Get cone polytopes and compute their faces in a single run of PALP
        dim = self.Delta_polar().dim()
        lps = [cone.lattice_polytope()
               for cone in self.fan()]
        for lp in lps:
            lp._dim = dim
        # Put import here in case we actually want to move this function to
        # another module
        from sage.geometry.lattice_polytope import all_faces
        all_faces(lps)
        # Now build the plot
        p = None
        for lp in lps:
            p += lp.plot3d(show_vindices=False, show_points=False)
        p += self.Delta_polar().plot3d()
        return p

    def resolve(self, **kwds):
        r"""
        Construct a toric variety whose fan subdivides the fan of ``self``.

        This function accepts only keyword arguments, none of which are
        mandatory.

        INPUT:

        - ``new_points`` -- list of integers, indices of boundary points of
          :meth:`Delta_polar`, which should be added as rays to the
          subdividing fan;

        - all other arguments will be passed to
          :meth:`~sage.schemes.generic.toric_variety.ToricVariety_field.resolve`
          method of (general) toric varieties, see its documentation for
          details.

        OUTPUT:

        - :class:`CPR-Fano toric variety <CPRFanoToricVariety_field>` if there
          was no ``new_rays`` argument and :class:`toric variety
          <sage.schemes.generic.toric_variety.ToricVariety_field>` otherwise.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: FTV = CPRFanoToricVariety(Delta=diamond)
            sage: FTV.coordinate_points()
            (0, 1, 2, 3)
            sage: FTV.gens()
            (z0, z1, z2, z3)
            sage: FTV_res = FTV.resolve(new_points=[6,8])
            Traceback (most recent call last):
            ...
            ValueError: the origin (point #6)
            cannot be used for subdivision!
            sage: FTV_res = FTV.resolve(new_points=[8,5])
            sage: FTV_res
            2-d CPR-Fano toric variety covered by 6 affine patches
            sage: FTV_res.coordinate_points()
            (0, 1, 2, 3, 8, 5)
            sage: FTV_res.gens()
            (z0, z1, z2, z3, z8, z5)

            sage: TV_res = FTV.resolve(new_rays=[(1,2)])
            sage: TV_res
            2-d toric variety covered by 5 affine patches
            sage: TV_res.gens()
            (z0, z1, z2, z3, z4)
        """
        # Reasons to override the base class:
        # - allow using polytope point indices for subdivision
        # - handle automatic name creation in a different fashion
        # - return CPR-Fano toric variety if the above feature was used and
        #   just toric variety if subdivision involves rays
        if "new_rays" in kwds:
            if "new_points" in kwds:
                raise ValueError("you cannot give new_points and new_rays at "
                                 "the same time!")
            return super(CPRFanoToricVariety_field, self).resolve(**kwds)
        # Now we need to construct another Fano variety
        new_points = kwds.pop("new_points", ())
        coordinate_points = self.coordinate_points()
        new_points = tuple(point for point in new_points
                                 if point not in coordinate_points)
        Delta_polar = self._Delta_polar
        if Delta_polar.origin() in new_points:
            raise ValueError("the origin (point #%d) cannot be used for "
                             "subdivision!" % Delta_polar.origin())
        if new_points:
            coordinate_points = coordinate_points + new_points
            point_to_ray = dict()
            for n, point in enumerate(coordinate_points):
                point_to_ray[point] = n
        else:
            point_to_ray = self._point_to_ray
        new_rays = [Delta_polar.point(point) for point in new_points]
        coordinate_name_indices = kwds.pop("coordinate_name_indices",
                                           coordinate_points)
        fan = self.fan()
        if "coordinate_names" in kwds:
            coordinate_names = kwds.pop("coordinate_names")
        else:
            coordinate_names = list(self.variable_names())
            coordinate_names.extend(normalize_names(ngens=len(new_rays),
                                indices=coordinate_name_indices[fan.nrays():],
                                prefix=self._coordinate_prefix))
            coordinate_names.append(self._coordinate_prefix + "+")
        rfan = fan.subdivide(new_rays=new_rays, **kwds)
        resolution = CPRFanoToricVariety_field(Delta_polar, rfan,
                            coordinate_points, point_to_ray, coordinate_names,
                            coordinate_name_indices, self.base_ring())
        R = self.coordinate_ring()
        R_res = resolution.coordinate_ring()
        resolution_map = resolution.hom(R.hom(R_res.gens()[:R.ngens()]), self)
        resolution._resolution_map = resolution_map
        return resolution


class AnticanonicalHypersurface(AlgebraicScheme_subscheme_toric):
    r"""
    Construct the anticanonical hypersurface of a CPR-Fano toric variety.

    INPUT:

    - ``P_Delta`` -- :class:`CPR-Fano toric variety
      <CPRFanoToricVariety_field>` associated to a reflexive polytope
      `\Delta`;

    -  see :meth:`CPRFanoToricVariety_field.anticanonical_hypersurface` for
       documentation on all other acceptable parameters.

    OUTPUT:

    - :class:`anticanonical hypersurface <AnticanonicalHypersurface>` of
      ``P_Delta`` (with the extended base field, if necessary).

    EXAMPLES::

        sage: P1xP1 = CPRFanoToricVariety(
        ...       Delta_polar=lattice_polytope.octahedron(2))
        sage: import sage.schemes.generic.fano_toric_variety as ftv
        sage: ftv.AnticanonicalHypersurface(P1xP1)
        Closed subscheme of 2-d CPR-Fano toric variety
        covered by 4 affine patches defined by:
          a1*z0^2*z1^2 + a0*z1^2*z2^2 + a6*z0*z1*z2*z3
        + a3*z0^2*z3^2 + a2*z2^2*z3^2

    See :meth:`~CPRFanoToricVariety_field.anticanonical_hypersurface()` for a
    more elaborate example.
    """
    def __init__(self, P_Delta, monomial_points=None, coefficient_names=None,
                 coefficient_name_indices=None):
        r"""
        See :meth:`CPRFanoToricVariety_field.anticanonical_hypersurface` for
        documentation.

        TESTS::

            sage: P1xP1 = CPRFanoToricVariety(
            ...       Delta_polar=lattice_polytope.octahedron(2))
            sage: import sage.schemes.generic.fano_toric_variety as ftv
            sage: ftv.AnticanonicalHypersurface(P1xP1)
            Closed subscheme of 2-d CPR-Fano toric variety
            covered by 4 affine patches defined by:
              a1*z0^2*z1^2 + a0*z1^2*z2^2 + a6*z0*z1*z2*z3
            + a3*z0^2*z3^2 + a2*z2^2*z3^2
        """
        if not is_CPRFanoToricVariety(P_Delta):
            raise TypeError("anticanonical hypersurfaces can only be "
                            "constructed for CPR-Fano toric varieties!"
                            "\nGot: %s" % P_Delta)
        Delta = P_Delta.Delta()
        Delta_polar = Delta.polar()
        # Monomial points normalization
        if monomial_points == "vertices":
            monomial_points = range(Delta.nvertices())
        elif monomial_points == "all":
            monomial_points = range(Delta.npoints())
        elif monomial_points == "vertices+origin":
            monomial_points = range(Delta.nvertices())
            monomial_points.append(Delta.origin())
        elif monomial_points == "simplified" or monomial_points is None:
            monomial_points = Delta.skeleton_points(Delta.dim() - 2)
            monomial_points.append(Delta.origin())
        elif isinstance(monomial_points, str):
            raise ValueError("%s is an unsupported description of monomial "
                             "points!" % monomial_points)
        monomial_points = tuple(monomial_points)
        self._monomial_points = monomial_points
        # Make the necessary ambient space
        if coefficient_name_indices is None:
            coefficient_name_indices = monomial_points
        coefficient_names = normalize_names(
                                coefficient_names, len(monomial_points),
                                DEFAULT_COEFFICIENT, coefficient_name_indices)
        self._coefficient_names = coefficient_names
        F = P_Delta.base_ring()
        if is_FractionField(F):
            # Q(a) ---> Q(a, b) rather than Q(a)(b)
            R = F.ring()
            if is_PolynomialRing(R) or is_MPolynomialRing(R):
                parameter_names = R.variable_names()
                new_names = [name for name in coefficient_names
                                  if name not in parameter_names]
                if new_names:
                    parameter_names.extend(new_names)
                    R = PolynomialRing(R.base_ring(), parameter_names)
                    F = R.fraction_field()
        else:
            F = PolynomialRing(F, coefficient_names).fraction_field()
        P_Delta = P_Delta.base_extend(F)
        # Defining polynomial
        h = sum(F(coef) * prod(P_Delta.coordinate_point_to_coordinate(n)
                                ** (Delta.point(m) * Delta_polar.point(n) + 1)
                          for n in P_Delta.coordinate_points())
            for m, coef in zip(monomial_points, coefficient_names))
        super(AnticanonicalHypersurface, self).__init__(P_Delta, h)
