r"""
Toric varieties

This module provides support for (normal) toric varieties, corresponding to
:class:`rational polyhedral fans <sage.geometry.fan.RationalPolyhedralFan>`.
See also :mod:`~sage.schemes.generic.fano_toric_variety` for a more
restrictive class of (weak) Fano toric varieties.

An **excellent reference on toric varieties** is the book "Toric Varieties" by
David A. Cox, John B. Little, and Hal Schenck [CLS]_. Its draft **is freely
available** at
http://www.cs.amherst.edu/~dac/toric.html
**but will be removed** from this site once it is published, so hurry up!

The interface to this module is provided through functions
:func:`AffineToricVariety` and :func:`ToricVariety`, although you may also be
interested in :func:`normalize_names`.

.. NOTE::

    We do NOT build "general toric varieties" from affine toric varieties.
    Instead, we are using the quotient representation of toric varieties with
    the homogeneous coordinate ring (a.k.a. Cox's ring or the total coordinate
    ring). This description works best for simplicial fans of the full
    dimension.

REFERENCES:

..  [CLS]
    David A. Cox, John B. Little,  Hal Schenck,
    "Toric Varieties", Graduate Studies in Mathematics,
    Amer. Math. Soc., Providence, RI, to appear.

AUTHORS:

- Andrey Novoseltsev (2010-05-17): initial version.

EXAMPLES:

We start with constructing the affine plane as an affine toric variety. First,
we need to have a corresponding cone::

    sage: quadrant = Cone([(1,0), (0,1)])

If you don't care about variable names and the base field, that's all we need
for now::

    sage: A2 = AffineToricVariety(quadrant)
    sage: A2
    2-d affine toric variety
    sage: origin = A2(0,0)
    sage: origin
    [0 : 0]

Only affine toric varieties have points whose (homogeneous) coordinates
are all zero. ::

    sage: parent(origin)
    Set of Rational Points of 2-d affine toric variety

As you can see, by default toric varieties live over the field of rational
numbers::

    sage: A2.base_ring()
    Rational Field

While usually toric varieties are considered over the field of complex
numbers, for computational purposes it is more convenient to work with fields
that have exact representation on computers. You can also always do ::

    sage: C2 = AffineToricVariety(quadrant, base_field=CC)
    sage: C2.base_ring()
    Complex Field with 53 bits of precision
    sage: C2(1,2+i)
    [1.00000000000000 : 2.00000000000000 + 1.00000000000000*I]

or even ::

    sage: F = CC["a, b"].fraction_field()
    sage: F.inject_variables()
    Defining a, b
    sage: A2 = AffineToricVariety(quadrant, base_field=F)
    sage: A2(a,b)
    [a : b]

OK, if you need to work only with affine spaces,
:func:`~sage.schemes.generic.affine_space.AffineSpace` may be a better way to
construct them. Our next example is the product of two projective lines
realized as the toric variety associated to the
:func:`face fan <sage.geometry.fan.FaceFan>` of the "diamond"::

    sage: diamond = lattice_polytope.octahedron(2)
    sage: diamond.vertices()
    [ 1  0 -1  0]
    [ 0  1  0 -1]
    sage: fan = FaceFan(diamond)
    sage: P1xP1 = ToricVariety(fan)
    sage: P1xP1
    2-d toric variety covered by 4 affine patches
    sage: P1xP1.fan().ray_matrix()
    [ 1  0 -1  0]
    [ 0  1  0 -1]
    sage: P1xP1.gens()
    (z0, z1, z2, z3)

We got four coordinates - two for each of the projective lines, but their
names are perhaps not very well chosen. Let's make `(x,y)` to be coordinates
on the first line and `(s,t)` on the second one::

    sage: P1xP1 = ToricVariety(fan, coordinate_names="x s y t")
    sage: P1xP1.gens()
    (x, s, y, t)

Now, if we want to define subschemes of this variety, the defining polynomials
must be homogeneous in each of these pairs::

    sage: P1xP1.inject_variables()
    Defining x, s, y, t
    sage: P1xP1.subscheme(x)
    Closed subscheme of 2-d toric variety
    covered by 4 affine patches defined by:
      x
    sage: P1xP1.subscheme(x^2 + y^2)
    Closed subscheme of 2-d toric variety
    covered by 4 affine patches defined by:
      x^2 + y^2
    sage: P1xP1.subscheme(x^2 + s^2)
    Traceback (most recent call last):
    ...
    ValueError: x^2 + s^2 is not homogeneous
    on 2-d toric variety covered by 4 affine patches!
    sage: P1xP1.subscheme([x^2*s^2 + x*y*t^2 +y^2*t^2, s^3 + t^3])
    Closed subscheme of 2-d toric variety
    covered by 4 affine patches defined by:
      x^2*s^2 + x*y*t^2 + y^2*t^2,
      s^3 + t^3

While we don't build toric varieties from affine toric varieties, we still can
access the "building pieces"::

    sage: patch = P1xP1.affine_patch(2)
    sage: patch
    2-d affine toric variety
    sage: patch.fan().ray_matrix()
    [1 0]
    [0 1]
    sage: patch.embedding_morphism()
    Scheme morphism:
      From: 2-d affine toric variety
      To:   2-d toric variety covered by 4 affine patches
      Defn: Defined on coordinates by sending [x : s] to
            [x : s : 1 : 1]

The patch above was specifically chosen to coincide with our representation of
the affine plane before, but you can get the other three patches as well.
(While any cone of a fan will correspond to an affine toric variety, the main
interest is usually in the generating fans as "the biggest" affine
subvarieties, and these are precisely the patches that you can get from
:meth:`~ToricVariety_field.affine_patch`.)

All two-dimensional toric varieties are "quite nice" because any
two-dimensional cone is generated by exactly two rays. From the point of view
of the corresponding toric varieties, this means that they have at worst
quotient singularities::

    sage: P1xP1.is_orbifold()
    True
    sage: P1xP1.is_smooth()
    True
    sage: TV = ToricVariety(NormalFan(diamond))
    sage: TV.fan().ray_matrix()
    [-1  1 -1  1]
    [ 1  1 -1 -1]
    sage: TV.is_orbifold()
    True
    sage: TV.is_smooth()
    False

In higher dimensions worse things can happen::

    sage: TV3 = ToricVariety(NormalFan(lattice_polytope.octahedron(3)))
    sage: TV3.fan().ray_matrix()
    [-1  1 -1  1 -1  1 -1  1]
    [-1 -1  1  1 -1 -1  1  1]
    [ 1  1  1  1 -1 -1 -1 -1]
    sage: TV3.is_orbifold()
    False

Fortunately, we can perform a (partial) resolution::

    sage: TV3_res = TV3.resolve_to_orbifold()
    sage: TV3_res.is_orbifold()
    True
    sage: TV3_res.fan().ngenerating_cones()
    12
    sage: TV3.fan().ngenerating_cones()
    6

In this example we had to double the number of affine patches. The result is
still singular::

    sage: TV3_res.is_smooth()
    False

You can resolve it further using :meth:`~ToricVariety_field.resolve` method,
but (at least for now) you will have to specify which rays should be inserted
into the fan. See also
:func:`~sage.schemes.generic.fano_toric_variety.FanoToricVariety`,
which can construct some other "nice partial resolutions."

Below you will find detailed descriptions of available functions. If you are
familiar with toric geometry, you will likely see that many important objects
and operations are unavailable. However, this module is under active
development and hopefully will improve in future releases of Sage. If there
are some particular features that you would like to see implemented ASAP,
please consider reporting them to the Sage Development Team or even
implementing them on your own as a patch for inclusion!
"""


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.geometry.all import Cone, Fan
from sage.matrix.all import matrix
from sage.misc.all import latex
from sage.modules.free_module_element import vector
from sage.rings.all import PolynomialRing, QQ, is_Field, is_FractionField
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme
from sage.schemes.generic.ambient_space import AmbientSpace
from sage.schemes.generic.homset import SchemeHomset_coordinates
from sage.schemes.generic.morphism import (SchemeMorphism_coordinates,
                                           SchemeMorphism_on_points,
                                           is_SchemeMorphism)
from sage.schemes.generic.scheme import is_Scheme
from sage.structure.sequence import Sequence


# Default prefix for indexed coordinates
DEFAULT_PREFIX = "z"


def is_ToricVariety(x):
    r"""
    Check if ``x`` is a toric variety.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a :class:`toric variety <ToricVariety_field>` and
      ``False`` otherwise.

    .. NOTE::

        While projective spaces are toric varieties mathematically, they are
        not toric varieties in Sage due to efficiency considerations, so this
        function will return ``False``.

    EXAMPLES::

        sage: from sage.schemes.generic.toric_variety import is_ToricVariety
        sage: is_ToricVariety(1)
        False
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P = ToricVariety(fan)
        sage: P
        2-d toric variety covered by 4 affine patches
        sage: is_ToricVariety(P)
        True
        sage: is_ToricVariety(ProjectiveSpace(2))
        False
    """
    return isinstance(x, ToricVariety_field)


def ToricVariety(fan,
                 coordinate_names=None,
                 coordinate_indices=None,
                 base_field=QQ):
    r"""
    Construct a toric variety.

    INPUT:

    - ``fan`` -- :class:`rational polyhedral fan
      <sage.geometry.fan.RationalPolyhedralFan>`;

    - ``coordinate_names`` -- names of variables for the coordinate ring, see
      :func:`normalize_names` for acceptable formats. If not given, indexed
      variable names will be created automatically;

    - ``coordinate_indices`` -- list of integers, indices for indexed
      variables. If not given, the index of each variable will coincide with
      the index of the corresponding ray of the fan;

    - ``base_field`` -- base field of the toric variety (default: `\QQ`).

    OUTPUT:

    - :class:`toric variety <ToricVariety_field>`.

    EXAMPLES:

    We will create the product of two projective lines::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: fan.ray_matrix()
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1.gens()
        (z0, z1, z2, z3)

    Let's create some points::

        sage: P1xP1(1,1,1,1)
        [1 : 1 : 1 : 1]
        sage: P1xP1(0,1,1,1)
        [0 : 1 : 1 : 1]
        sage: P1xP1(0,1,0,1)
        Traceback (most recent call last):
        ...
        TypeError: coordinates (0, 1, 0, 1)
        are in the exceptional set!

    We cannot set to zero both coordinates of the same projective line!

    Let's change the names of the variables. We have to re-create our toric
    variety::

        sage: P1xP1 = ToricVariety(fan, "x s y t")
        sage: P1xP1.gens()
        (x, s, y, t)

    Now `(x, y)` correspond to one line and `(s, t)` to the other one. ::

        sage: P1xP1.inject_variables()
        Defining x, s, y, t
        sage: P1xP1.subscheme(x*s-y*t)
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s - y*t
    """
    if not is_Field(base_field):
        raise TypeError("need a field to construct a toric variety!\n Got %s"
                        % base_field)
    return ToricVariety_field(fan, coordinate_names, coordinate_indices,
                             base_field)


def AffineToricVariety(cone, *args, **kwds):
    r"""
    Construct an affine toric variety.

    INPUT:

    - ``cone`` -- :class:`strictly convex rational polyhedral cone
      <sage.geometry.cone.ConvexRationalPolyhedralCone>`.

    This cone will be used to construct a :class:`rational polyhedral fan
    <sage.geometry.fan.RationalPolyhedralFan>`, which will be passed to
    :func:`ToricVariety` with the rest of positional and keyword arguments.

    OUTPUT:

    - :class:`toric variety <ToricVariety_field>`.

    .. NOTE::

        The generating rays of the fan of this variety are guaranteed to be
        listed in the same order as the rays of the original cone.

    EXAMPLES:

    We will create the affine plane as an affine toric variety::

        sage: quadrant = Cone([(1,0), (0,1)])
        sage: A2 = AffineToricVariety(quadrant)
        sage: origin = A2(0,0)
        sage: origin
        [0 : 0]
        sage: parent(origin)
        Set of Rational Points of 2-d affine toric variety

    Only affine toric varieties have points whose (homogeneous) coordinates
    are all zero.
    """
    if not cone.is_strictly_convex():
        raise ValueError("affine toric varieties are defined for strictly "
                         "convex cones only!")
    # We make sure that Fan constructor does not meddle with the order of
    # rays, this is very important for affine patches construction
    fan = Fan([tuple(range(cone.nrays()))], cone.rays(),
              check=False, normalize=False)
    return ToricVariety(fan, *args, **kwds)


class ToricVariety_field(AmbientSpace):
    r"""
    Construct a toric variety associated to a rational polyhedral fan.

    .. WARNING::

        This class does not perform any checks of correctness of input. Use
        :func:`ToricVariety` and :func:`AffineToricVariety` to construct toric
        varieties.

    INPUT:

    - ``fan`` -- :class:`rational polyhedral fan
      <sage.geometry.fan.RationalPolyhedralFan>`;

    - ``coordinate_names`` -- names of variables, see :func:`normalize_names`
      for acceptable formats. If ``None``, indexed variable names will be
      created automatically;

    - ``coordinate_indices`` -- list of integers, indices for indexed
      variables. If ``None``, the index of each variable will coincide with
      the index of the corresponding ray of the fan;

    - ``base_field`` -- base field of the toric variety.

    OUTPUT:

    - :class:`toric variety <ToricVariety_field>`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
    """

    def __init__(self, fan, coordinate_names, coordinate_indices, base_field):
        r"""
        See :class:`ToricVariety_field` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
        """
        self._fan = fan
        super(ToricVariety_field, self).__init__(fan.lattice_dim(),
                                                 base_field)
        self._torus_factor_dim = fan.lattice_dim() - fan.dim()
        coordinate_names = normalize_names(coordinate_names,
                        fan.nrays() + self._torus_factor_dim, DEFAULT_PREFIX,
                        coordinate_indices, return_prefix=True)
        # Save the prefix for use in resolutions
        self._coordinate_prefix = coordinate_names.pop()
        self._assign_names(names=coordinate_names, normalize=False)

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is of the same type as ``self``, their fans are the
          same, names of variables are the same and stored in the same order,
          and base fields are the same. 1 or -1 otherwise.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1a = ToricVariety(fan, "x s y t")
            sage: P1xP1b = ToricVariety(fan)
            sage: cmp(P1xP1, P1xP1a)
            1
            sage: cmp(P1xP1a, P1xP1)
            -1
            sage: cmp(P1xP1, P1xP1b)
            0
            sage: P1xP1 is P1xP1b
            False
            sage: cmp(P1xP1, 1)
            -1
        """
        c = cmp(type(self), type(right))
        if c:
            return c
        return cmp([self.fan(),
                    self.variable_names(),
                    self.base_ring()],
                   [right.fan(),
                    right.variable_names(),
                    right.base_ring()])

    def _an_element_(self):
        r"""
        Construct an element of ``self``.

        This function is needed (in particular) for the test framework.

        OUTPUT:

        - a point of ``self`` with coordinates [1 : 2: ... : n].

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._an_element_()
            [1 : 2 : 3 : 4]
        """
        return self(range(1, self.ngens() + 1))

    def _check_satisfies_equations(self, coordinates):
        r"""
        Check if ``coordinates`` define a valid point of ``self``.

        INPUT:

        - ``coordinates`` -- list of elements of the base field of ``self``.

        OUTPUT:

        - ``True`` if ``coordinates`` do define a valid point of ``self``,
          otherwise a ``TypeError`` or ``ValueError`` exception is raised.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._check_satisfies_equations([1,1,1,1])
            True
            sage: P1xP1._check_satisfies_equations([0,0,1,1])
            True
            sage: P1xP1._check_satisfies_equations([0,1,0,1])
            Traceback (most recent call last):
            ...
            TypeError: coordinates (0, 1, 0, 1)
            are in the exceptional set!
            sage: P1xP1._check_satisfies_equations([1,1,1])
            Traceback (most recent call last):
            ...
            TypeError: coordinates (1, 1, 1) must have 4 components!
            sage: P1xP1._check_satisfies_equations(1)
            Traceback (most recent call last):
            ...
            TypeError: 1 can not be used as coordinates!
            Use a list or a tuple.
            sage: P1xP1._check_satisfies_equations([1,1,1,fan])
            Traceback (most recent call last):
            ...
            TypeError: coordinate Rational polyhedral fan
            in 2-d lattice N is not an element of Rational Field!
        """
        try:
            coordinates = tuple(coordinates)
        except TypeError:
            raise TypeError("%s can not be used as coordinates! "
                            "Use a list or a tuple." % coordinates)
        n = self.ngens()
        if len(coordinates) != n:
            raise TypeError("coordinates %s must have %d components!"
                            % (coordinates, n))
        base_field = self.base_ring()
        for coordinate in coordinates:
            if coordinate not in base_field:
                raise TypeError("coordinate %s is not an element of %s!"
                                % (coordinate, base_field))
        zero_positions = set(position
                            for position, coordinate in enumerate(coordinates)
                            if coordinate == 0)
        if not zero_positions:
            return True
        for i in range(n - self._torus_factor_dim, n):
            if i in zero_positions:
                raise ValueError("coordinates on the torus factor cannot be "
                                 "zero! Got %s" % str(coordinates))
        if len(zero_positions) == 1:
            return True
        fan = self.fan()
        possible_charts = set(fan._ray_to_cones(zero_positions.pop()))
        for i in zero_positions:
            possible_charts.intersection_update(fan._ray_to_cones(i))
        if possible_charts:
            return True     # All zeros are inside one generating cone
        raise TypeError("coordinates %s are in the exceptional set!"
                        % str(coordinates)) # Need str, coordinates is a tuple

    def _homset_class(self, *args, **kwds):
        r"""
        Construct a `Hom`-space for ``self``.

        INPUT:

        - same as for :class:`SchemeHomset_toric_coordinates_field`.

        OUPUT:

        - :class:`SchemeHomset_toric_coordinates_field`.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._homset_class(P1xP1, P1xP1.base_ring())
            Set of Rational Points of 2-d toric variety
            covered by 4 affine patches
        """
        return SchemeHomset_toric_coordinates_field(*args, **kwds)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._latex_()
            '\\mathbb{X}_{\\Sigma^{2}}'
        """
        return r"\mathbb{X}_{%s}" % latex(self.fan())

    def _latex_generic_point(self, coordinates=None):
        """
        Return a LaTeX representation of a point of ``self``.

        INPUT:

        - ``coordinates`` -- list of coordinates of a point of ``self``.
          If not given, names of coordinates of ``self`` will be used.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._latex_generic_point()
            '\\left[z_{0} : z_{1} : z_{2} : z_{3}\\right]'
            sage: P1xP1._latex_generic_point([1,2,3,4])
            '\\left[1 : 2 : 3 : 4\\right]'
        """
        if coordinates is None:
            coordinates = self.gens()
        return r"\left[%s\right]" % (" : ".join(str(latex(coord))
                                                for coord in coordinates))

    def _point_class(self, *args, **kwds):
        r"""
        Construct a point of ``self``.

        INPUT:

        - same as for :class:`SchemeMorphism_toric_coordinates_field`.

        OUPUT:

        - :class:`SchemeMorphism_toric_coordinates_field`.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._point_class(P1xP1, [1,2,3,4])
            [1 : 2 : 3 : 4]
        """
        return SchemeMorphism_toric_coordinates_field(*args, **kwds)

    def _point_morphism_class(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        INPUT:

        - same as for :class:`SchemeMorphism_on_points_toric_variety`.

        OUPUT:

        - :class:`SchemeMorphism_on_points_toric_variety`.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: H = P1xP1.Hom(P1)
            sage: P1xP1._point_morphism_class(H, [z0,z1,z0,z3])
            Scheme morphism:
              From: 2-d toric variety covered by 4 affine patches
              To:   Closed subscheme of 2-d toric variety
              covered by 4 affine patches defined by:
              z0 - z2
              Defn: Defined on coordinates by sending
                    [z0 : z1 : z2 : z3] to [z0 : z1 : z0 : z3]
        """
        return SchemeMorphism_on_points_toric_variety(*args, **kwds)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._repr_()
            '2-d toric variety covered by 4 affine patches'
        """
        result = "%d-d" % self.dimension_relative()
        if self.fan().ngenerating_cones() == 1:
            result += " affine toric variety"
        else:
            result += (" toric variety covered by %d affine patches"
                       % self.fan().ngenerating_cones())
        return result

    def _repr_generic_point(self, coordinates=None):
        r"""
        Return a string representation of a point of ``self``.

        INPUT:

        - ``coordinates`` -- list of coordinates of a point of ``self``.
          If not given, names of coordinates of ``self`` will be used.

        OUTPUT:

        - string.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1._repr_generic_point()
            '[z0 : z1 : z2 : z3]'
            sage: P1xP1._repr_generic_point([1,2,3,4])
            '[1 : 2 : 3 : 4]'
        """
        if coordinates is None:
            coordinates = self.gens()
        return "[%s]" % (" : ".join(str(coord) for coord in coordinates))

    def _validate(self, polynomials):
        """
        Check if ``polynomials`` define valid functions on ``self``.

        Since this is a toric variety, polynomials must be homogeneous in the
        total coordinate ring of ``self``.

        INPUT:

        - ``polynomials`` -- list of polynomials in the coordinate ring of
          ``self`` (this function does not perform any conversions).

        OUTPUT:

        - ``polynomials`` (the input parameter without any modifications) if
          ``polynomials`` do define valid polynomial functions on ``self``,
          otherwise a ``ValueError`` exception is raised.

        TESTS:

        We will use the product of two projective lines with coordinates
        `(x, y)` for one and `(s, t)` for the other::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.ray_matrix()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: P1xP1._validate([x - y, x*s + y*t])
            [x - y, x*s + y*t]
            sage: P1xP1._validate([x + s])
            Traceback (most recent call last):
            ...
            ValueError: x + s is not homogeneous on
            2-d toric variety covered by 4 affine patches!
        """
        for p in polynomials:
            if not self.is_homogeneous(p):
                raise ValueError("%s is not homogeneous on %s!" % (p, self))
        return polynomials

    def affine_patch(self, i):
        r"""
        Return the ``i``-th affine patch of ``self``.

        INPUT:

        - ``i`` -- integer, index of a generating cone of the fan of ``self``.

        OUTPUT:

        - affine :class:`toric variety <ToricVariety_field>` corresponding to
          the ``i``-th generating cone of the fan of ``self``.

        The result is cached, so the ``i``-th patch is always the same object
        in memory.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: patch0 = P1xP1.affine_patch(0)
            sage: patch0
            2-d affine toric variety
            sage: patch0.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [x : t] to
                    [x : 1 : 1 : t]
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : t] to
                    [1 : 1 : y : t]
            sage: patch1 is P1xP1.affine_patch(1)
            True
        """
        i = int(i)   # implicit type checking
        try:
            return self._affine_patches[i]
        except AttributeError:
            self._affine_patches = dict()
        except KeyError:
            pass
        cone = self.fan().generating_cone(i)
        names = self.variable_names()
        # Number of "honest fan coordinates"
        n = self.fan().nrays()
        # Number of "torus factor coordinates"
        t = self._torus_factor_dim
        names = ([names[ray] for ray in cone.ambient_ray_indices()]
                 + list(names[n:]))
        patch = AffineToricVariety(cone, names, base_field=self.base_ring())
        embedding_coordinates = [1] * n
        for k, ray in enumerate(cone.ambient_ray_indices()):
            embedding_coordinates[ray] = patch.gen(k)
        if t > 0: # Passing "-0" gives unintended result
            embedding_coordinates.extend(patch.gens()[-t:])
        patch._embedding_morphism = patch.hom(embedding_coordinates, self)
        self._affine_patches[i] = patch
        return patch

    def change_ring(self, F):
        r"""
        Return a toric variety over ``F`` and otherwise the same as ``self``.

        INPUT:

        - ``F`` -- field.

        OUTPUT:

        - :class:`toric variety <ToricVariety_field>` over ``F``.

        .. NOTE::

            There is no need to have any relation between ``F`` and the base
            field of ``self``. If you do want to have such a relation, use
            :meth:`base_extend` instead.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
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
            return ToricVariety(self.fan(), self.variable_names(),
                                base_field=F)

    def coordinate_ring(self):
        r"""
        Return the coordinate ring of ``self``.

        For toric varieties this is the homogeneous coordinate ring (a.k.a.
        Cox's ring and total ring).

        OUTPUT:

        - polynomial ring.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.coordinate_ring()
            Multivariate Polynomial Ring in z0, z1, z2, z3
            over Rational Field
        """
        if "_coordinate_ring" not in self.__dict__:
            self._coordinate_ring = PolynomialRing(self.base_ring(),
                                                   self.variable_names())
        return self._coordinate_ring

    def embedding_morphism(self):
        r"""
        Return the default embedding morphism of ``self``.

        Such a morphism is always defined for an affine patch of a toric
        variety (which is also a toric varieties itself).

        OUTPUT:

        - :class:`scheme morphism <SchemeMorphism_on_points_toric_variety>`
          if the default embedding morphism was defined for ``self``,
          otherwise a ``ValueError`` exception is raised.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.embedding_morphism()
            Traceback (most recent call last):
            ...
            ValueError: no default embedding was
            defined for this toric variety!
            sage: patch = P1xP1.affine_patch(0)
            sage: patch
            2-d affine toric variety
            sage: patch.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [x : t] to
                    [x : 1 : 1 : t]
        """
        try:
            return self._embedding_morphism
        except AttributeError:
            raise ValueError("no default embedding was defined for this "
                             "toric variety!")

    def fan(self, dim=None, codim=None):
        r"""
        Return the underlying fan of ``self`` or its cones.

        INPUT:

        - ``dim`` -- dimension of the requested cones;

        - ``codim`` -- codimension of the requested cones.

        OUTPUT:

        - :class:`rational polyhedral fan
          <sage.geometry.fan.RationalPolyhedralFan>` if no parameters were
          given, :class:`tuple` of :class:`cones
          <sage.geometry.cone.ConvexRationalPolyhedralCone>` otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.fan()
            Rational polyhedral fan in 2-d lattice N
            sage: P1xP1.fan() is fan
            True
            sage: P1xP1.fan(1)[0]
            1-d cone of Rational polyhedral fan in 2-d lattice N
        """
        return self._fan(dim, codim)

    def inject_coefficients(self, scope=None, verbose=True):
        r"""
        Inject generators of the base field of ``self`` into ``scope``.

        This function is useful if the base field is the field of rational
        functions.

        INPUT:

        - ``scope`` -- namespace (default: global);

        - ``verbose`` -- if ``True`` (default), names of injected generators
          will be printed.

        OUTPUT:

        - none.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_coefficients()

        The last command does nothing, since ``P1xP1`` is defined over `\QQ`.
        Let's construct a toric variety over a more complicated field::

            sage: F = QQ["a, b"].fraction_field()
            sage: P1xP1 = ToricVariety(fan, base_field=F)
            sage: P1xP1.inject_coefficients()
            Defining a, b
        """
        if is_FractionField(self.base_ring()):
            self.base_ring().inject_variables(scope, verbose)

    def is_homogeneous(self, polynomial):
        r"""
        Check if ``polynomial`` is homogeneous.

        The coordinate ring of a toric variety is multigraded by relations
        between generating rays of the underlying fan.

        INPUT:

        - ``polynomial`` -- polynomial in the coordinate ring of ``self`` or
          its quotient.

        OUTPUT:

        - ``True`` if ``polynomial`` is homogeneous and ``False`` otherwise.

        EXAMPLES:

        We will use the product of two projective lines with coordinates
        `(x, y)` for one and `(s, t)` for the other::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.ray_matrix()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: P1xP1.is_homogeneous(x - y)
            True
            sage: P1xP1.is_homogeneous(x*s + y*t)
            True
            sage: P1xP1.is_homogeneous(x - t)
            False
        """
        if "_relation_matrix" not in self.__dict__:
            m = self.fan().ray_matrix().transpose().kernel().matrix()
            # We ignore degrees of torus factor coordinates
            m = m.augment(matrix(m.nrows(), self._torus_factor_dim))
            self._relation_matrix = m
        relation_matrix = self._relation_matrix
        if polynomial not in self.coordinate_ring():
            # Then it should be in the quotient corresponding to a subscheme
            polynomial = polynomial.lift()
        monomials = polynomial.monomials()
        if not monomials:
            return True
        degree = relation_matrix * vector(monomials[0].degrees())
        for monomial in monomials:
            if relation_matrix * vector(monomial.degrees()) != degree:
                return False
        return True

    def is_isomorphic(self, another):
        r"""
        Check if ``self`` is isomorphic to ``another``.

        INPUT:

        - ``another`` - :class:`toric variety <ToricVariety_field>`.

        OUTPUT:

        - ``True`` if ``self`` and ``another`` are isomorphic,
          ``False`` otherwise.

        EXAMPLES::

            sage: fan1 = FaceFan(lattice_polytope.octahedron(2))
            sage: TV1 = ToricVariety(fan1)
            sage: fan2 = NormalFan(lattice_polytope.octahedron(2))
            sage: TV2 = ToricVariety(fan2)

        Only the most trivial case is implemented so far::

            sage: TV1.is_isomorphic(TV1)
            True
            sage: TV1.is_isomorphic(TV2)
            Traceback (most recent call last):
            ...
            NotImplementedError:
            isomorphism check is not yet implemented!
        """
        if self is another:
            return True
        if not is_ToricVariety(another):
            raise TypeError(
                "only another toric variety can be checked for isomorphism! "
                "Got %s" % another)
        raise NotImplementedError("isomorphism check is not yet implemented!")

    def is_complete(self):
        r"""
        Check if ``self`` is complete.

        OUTPUT:

        - ``True`` if ``self`` is complete and ``False`` otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.is_complete()
            True
            sage: P1xP1.affine_patch(0).is_complete()
            False
        """
        return self.fan().is_complete()

    def is_orbifold(self):
        r"""
        Check if ``self`` has only quotient singularities.

        OUTPUT:

        - ``True`` if ``self`` has only quotient singularities, ``False``
          otherwise.

        EXAMPLES::

            sage: fan1 = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan1)
            sage: P1xP1.is_orbifold()
            True
            sage: fan2 = NormalFan(lattice_polytope.octahedron(3))
            sage: TV = ToricVariety(fan2)
            sage: TV.is_orbifold()
            False
        """
        return self.fan().is_simplicial()

    def is_smooth(self):
        r"""
        Check if ``self`` is smooth.

        OUTPUT:

        - ``True`` if ``self`` is smooth and ``False`` otherwise.

        EXAMPLES::

            sage: diamond = lattice_polytope.octahedron(2)
            sage: fan1 = FaceFan(diamond)
            sage: P1xP1 = ToricVariety(fan1)
            sage: P1xP1.is_smooth()
            True
            sage: fan2 = NormalFan(lattice_polytope.octahedron(2))
            sage: TV = ToricVariety(fan2)
            sage: TV.is_smooth()
            False
        """
        return self.fan().is_smooth()

    def resolve(self, **kwds):
        r"""
        Construct a toric variety whose fan subdivides the fan of ``self``.

        The name of this function reflects the fact that usually such
        subdivisions are done for resolving singularities of the original
        variety.

        INPUT:

        This function accepts only keyword arguments, none of which are
        mandatory.

        - ``coordinate_names`` -- names for coordinates of the new variety. If
          not given, will be constructed from the coordinate names of ``self``
          and necessary indexed ones. See :func:`normalize_names` for the
          description of acceptable formats;

        - ``coordinate_indices`` -- coordinate indices which should be used
          for indexed variables of the new variety;

        - all other arguments will be passed to
          :meth:`~sage.geometry.fan.RationalPolyhedralFan.subdivide` method of
          the underlying :class:`rational polyhedral fan
          <sage.geometry.fan.RationalPolyhedralFan>`, see its documentation
          for the available options.

        OUTPUT:

        - :class:`toric variety <ToricVariety_field>`.

        EXAMPLES:

        First we will "manually" resolve a simple orbifold singularity::

            sage: cone = Cone([(1,1), (-1,1)])
            sage: fan = Fan([cone])
            sage: TV = ToricVariety(fan)
            sage: TV.is_smooth()
            False
            sage: TV_res = TV.resolve(new_rays=[(0,1)])
            sage: TV_res.is_smooth()
            True
            sage: TV_res.fan().ray_matrix()
            [-1  1  0]
            [ 1  1  1]
            sage: [cone.ambient_ray_indices() for cone in TV_res.fan()]
            [(1, 2), (0, 2)]

        Now let's "automatically" partially resolve a more complicated fan::

            sage: fan = NormalFan(lattice_polytope.octahedron(3))
            sage: TV = ToricVariety(fan)
            sage: TV.is_smooth()
            False
            sage: TV.is_orbifold()
            False
            sage: TV.fan().nrays()
            8
            sage: TV.fan().ngenerating_cones()
            6
            sage: TV_res = TV.resolve(make_simplicial=True)
            sage: TV_res.is_smooth()
            False
            sage: TV_res.is_orbifold()
            True
            sage: TV_res.fan().nrays()
            8
            sage: TV_res.fan().ngenerating_cones()
            12
            sage: TV.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)
            sage: TV_res.gens()
            (z0, z1, z2, z3, z4, z5, z6, z7)
            sage: TV_res = TV.resolve(coordinate_names="x+",
            ...                       make_simplicial=True)
            sage: TV_res.gens()
            (x0, x1, x2, x3, x4, x5, x6, x7)
        """
        # If you are changing this function, check out resolve in Fano toric
        # varieties to see if it should be changed too
        #
        # Currently the resolution of fans works for full-dimensional ones
        # only, so there is no point to deal with the general case here, since
        # we will not be able to check that it works.
        coordinate_names = kwds.pop("coordinate_names", None)
        coordinate_indices = kwds.pop("coordinate_indices", None)
        fan = self.fan()
        if fan.dim() != fan.lattice_dim():
            raise NotImplementedError("resolution of toric varieties with "
                                      "torus factors is not yet implemented!")
            # When it is implemented, should be careful with the torus factor
        rfan = fan.subdivide(**kwds)
        if coordinate_names is None:
            coordinate_names = list(self.variable_names())
            if coordinate_indices is None:
                coordinate_indices = range(fan.nrays(), rfan.nrays())
            else:
                coordinate_indices = coordinate_indices[fan.nrays():]
            coordinate_names.extend(normalize_names(
                                    ngens=rfan.nrays() - fan.nrays(),
                                    indices=coordinate_indices,
                                    prefix=self._coordinate_prefix))
            coordinate_names.append(self._coordinate_prefix + "+")
        resolution = ToricVariety(rfan, coordinate_names=coordinate_names,
                                  coordinate_indices=coordinate_indices,
                                  base_field=self.base_ring())
        R = self.coordinate_ring()
        R_res = resolution.coordinate_ring()
        resolution_map = resolution.hom(R.hom(R_res.gens()[:R.ngens()]), self)
        resolution._resolution_map = resolution_map
        # The above map does not have (yet) public methods to access it.
        # While this map is defined correctly, base classes of schemes and
        # morphisms do not treat it as they should. The plan is to fix this
        # situation soon and to be able to use this map!
        return resolution

    def resolve_to_orbifold(self, **kwds):
        r"""
        Construct an orbifold whose fan subdivides the fan of ``self``.

        It is a synonym for :meth:`resolve` with ``make_simplicial=True``
        option.

        INPUT:

        - this function accepts only keyword arguments. See :meth:`resolve`
          for documentation.

        OUTPUT:

        - :class:`toric variety <ToricVariety_field>`.

        EXAMPLES::

            sage: fan = NormalFan(lattice_polytope.octahedron(3))
            sage: TV = ToricVariety(fan)
            sage: TV.is_orbifold()
            False
            sage: TV.fan().nrays()
            8
            sage: TV.fan().ngenerating_cones()
            6
            sage: TV_res = TV.resolve_to_orbifold()
            sage: TV_res.is_orbifold()
            True
            sage: TV_res.fan().nrays()
            8
            sage: TV_res.fan().ngenerating_cones()
            12
        """
        return self.resolve(make_simplicial=True, **kwds)

    def subscheme(self, polynomials):
        r"""
        Return the subscheme of ``self`` defined by ``polynomials``.

        INPUT:

        - ``polynomials`` -- list of polynomials in the coordinate ring of
          ``self``.

        OUTPUT:

        - :class:`subscheme of a toric variety
          <AlgebraicScheme_subscheme_toric>`.

        EXAMPLES:

        We will construct a subscheme of the product of two projective lines
        with coordinates `(x, y)` for one and `(s, t)` for the other::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.ray_matrix()
            [ 1  0 -1  0]
            [ 0  1  0 -1]
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: X = P1xP1.subscheme([x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d toric variety
            covered by 4 affine patches defined by:
              x*s + y*t,
              x^3 + y^3
            sage: X.defining_polynomials()
            (x*s + y*t, x^3 + y^3)
            sage: X.defining_ideal()
            Ideal (x*s + y*t, x^3 + y^3)
            of Multivariate Polynomial Ring in x, s, y, t
            over Rational Field
            sage: X.base_ring()
            Rational Field
            sage: X.base_scheme()
            Spectrum of Rational Field
            sage: X.structure_morphism()
            Scheme morphism:
              From: Closed subscheme of
              2-d toric variety covered by 4 affine patches defined by:
              x*s + y*t,
              x^3 + y^3
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        return AlgebraicScheme_subscheme_toric(self, polynomials)


class AlgebraicScheme_subscheme_toric(AlgebraicScheme_subscheme):
    r"""
    Construct an algebraic subscheme of a toric variety.

    .. WARNING::

        You should not create objects of this class directly. The preferred
        method to construct such subschemes is to use
        :meth:`~ToricVariety_field.subscheme` method of
        :class:`toric varieties <ToricVariety_field>`.

    INPUT:

    - ``toric_variety`` -- ambient :class:`toric variety
      <ToricVariety_field>`;

    - ``polynomials`` -- single polynomial, list, or ideal of defining
      polynomials in the coordinate ring of ``toric_variety``.

    OUTPUT:

    - :class:`algebraic subscheme of a toric variety
      <AlgebraicScheme_subscheme_toric>`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan, "x s y t")
        sage: P1xP1.inject_variables()
        Defining x, s, y, t
        sage: import sage.schemes.generic.toric_variety as tv
        sage: X = tv.AlgebraicScheme_subscheme_toric(
        ...         P1xP1, [x*s + y*t, x^3+y^3])
        sage: X
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3

    A better way to construct the same scheme as above::

        sage: P1xP1.subscheme([x*s + y*t, x^3+y^3])
        Closed subscheme of 2-d toric variety
        covered by 4 affine patches defined by:
          x*s + y*t,
          x^3 + y^3
    """

    def __init__(self, toric_variety, polynomials):
        r"""
        See :class:`AlgebraicScheme_subscheme_toric` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: import sage.schemes.generic.toric_variety as tv
            sage: X = tv.AlgebraicScheme_subscheme_toric(
            ...         P1xP1, [x*s + y*t, x^3+y^3])
            sage: X
            Closed subscheme of 2-d toric variety
            covered by 4 affine patches defined by:
              x*s + y*t,
              x^3 + y^3
        """
        # Just to make sure that keyword arguments will be passed correctly
        super(AlgebraicScheme_subscheme_toric, self).__init__(toric_variety,
                                                              polynomials)

    def _point_morphism_class(self, *args, **kwds):
        r"""
        Construct a morphism determined by action on points of ``self``.

        INPUT:

        - same as for :class:`SchemeMorphism_on_points_toric_variety`.

        OUPUT:

        - :class:`SchemeMorphism_on_points_toric_variety`.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: H = P1.Hom(P1xP1)
            sage: P1._point_morphism_class(H, [z0,z1,z0,z3])
            Scheme morphism:
              From: Closed subscheme of 2-d toric variety
              covered by 4 affine patches defined by:
              z0 - z2
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [z0 : z1 : z2 : z3] to
                    [z2bar : z1bar : z2bar : z3bar]
        """
        return SchemeMorphism_on_points_toric_variety(*args, **kwds)

    def affine_patch(self, i):
        r"""
        Return the ``i``-th affine patch of ``self``.

        INPUT:

        - ``i`` -- integer, index of a generating cone of the fan of the
          ambient space of ``self``.

        OUTPUT:

        - subscheme of an affine :class:`toric variety <ToricVariety_field>`
          corresponding to the pull-back of ``self`` by the embedding morphism
          of the ``i``-th :meth:`affine patch of the ambient space
          <ToricVariety_field.affine_patch>` of ``self``.

        The result is cached, so the ``i``-th patch is always the same object
        in memory.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : t] to
                    [1 : 1 : y : t]
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: P1 = P1xP1.subscheme(x-y)
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
        """
        i = int(i)   # implicit type checking
        try:
            return self._affine_patches[i]
        except AttributeError:
            self._affine_patches = dict()
        except KeyError:
            pass
        ambient_patch = self.ambient_space().affine_patch(i)
        phi_p = ambient_patch.embedding_morphism().defining_polynomials()
        patch = ambient_patch.subscheme(
                            [p(phi_p) for p in self.defining_polynomials()])
        patch._embedding_morphism = patch.hom(phi_p, self)
        self._affine_patches[i] = patch
        return patch

    def dimension(self):
        """
        Return the dimension of ``self``.

        .. NOTE::

            Currently the dimension of subschemes of toric varieties can be
            returned only if it was somehow set before.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: P1.dimension()
            Traceback (most recent call last):
            ...
            NotImplementedError:
            cannot compute dimension of this scheme!
        """
        if "_dimension" not in self.__dict__:
            raise NotImplementedError(
                                "cannot compute dimension of this scheme!")
        return self._dimension

    def embedding_morphism(self):
        r"""
        Return the default embedding morphism of ``self``.

        Such a morphism is always defined for an affine patch of a subscheme
        of a toric variety (which is a subscheme of a toric variety itself).

        OUTPUT:

        - :class:`scheme morphism <SchemeMorphism_on_points_toric_variety>`
          if the default embedding morphism was defined for ``self``,
          otherwise a ``ValueError`` exception is raised.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan, "x s y t")
            sage: patch1 = P1xP1.affine_patch(1)
            sage: patch1.embedding_morphism()
            Scheme morphism:
              From: 2-d affine toric variety
              To:   2-d toric variety covered by 4 affine patches
              Defn: Defined on coordinates by sending [y : t] to
                    [1 : 1 : y : t]
            sage: P1xP1.inject_variables()
            Defining x, s, y, t
            sage: P1 = P1xP1.subscheme(x-y)
            sage: P1.embedding_morphism()
            Traceback (most recent call last):
            ...
            ValueError: no default embedding was defined
            for this subscheme of a toric variety!
            sage: subpatch = P1.affine_patch(1)
            sage: subpatch
            Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
            sage: subpatch.embedding_morphism()
            Scheme morphism:
              From: Closed subscheme of 2-d affine toric variety defined by:
              -y + 1
              To:   Closed subscheme of 2-d toric variety
              covered by 4 affine patches defined by:
              x - y
              Defn: Defined on coordinates by sending [y : t] to
                    [1 : 1 : 1 : tbar]
        """
        try:
            return self._embedding_morphism
        except AttributeError:
            raise ValueError("no default embedding was defined for this "
                             "subscheme of a toric variety!")


class SchemeHomset_toric_coordinates_field(SchemeHomset_coordinates):
    """
    Construct the `Hom`-space of morphisms given on coordinates.

    .. WARNING::

        You should not create objects of this class directly.


    INPUT:

    - same as for :class:`sage.schemes.generic.SchemeHomset_coordinates`.

    OUPUT:

    - :class:`SchemeHomset_toric_coordinates_field`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: import sage.schemes.generic.toric_variety as tv
        sage: tv.SchemeHomset_toric_coordinates_field(P1xP1, QQ)
        Set of Rational Points of 2-d toric variety
        covered by 4 affine patches

    A better way to construct the same `Hom`-space as above::

        sage: P1xP1(QQ)
        Set of Rational Points of 2-d toric variety
        covered by 4 affine patches
    """
    # Mimicking SchemeHomset_projective_coordinates_field,
    # affine spaces implement only "except" case
    def __call__(self, *arg):
        r"""
        Construct a morphism from given parameters.

        INPUT:

        - data determining a morphism.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: import sage.schemes.generic.toric_variety as tv
            sage: H = tv.SchemeHomset_toric_coordinates_field(P1xP1, QQ)
            sage: H(1,2,3,4)
            [1 : 2 : 3 : 4]
        """
        # This may break for one-dimensional varieties.
        if len(arg) == 1:
            arg = arg[0]
        X = self.codomain()
        try:
            return X._point_class(X, arg)
        except AttributeError:  # should be very rare
            return SchemeMorphism_toric_coordinates_field(self, arg)


class SchemeMorphism_toric_coordinates_field(SchemeMorphism_coordinates):
    """
    Construct a morphism determined by giving coordinates in a field.

    .. WARNING::

        You should not create objects of this class directly.

    INPUT:

    - ``X`` -- subscheme of a toric variety.

    - ``coordinates`` -- list of coordinates in the base field of ``X``.

    - ``check`` -- if ``True`` (default), the input will be checked for
      correctness.

    OUTPUT:

    - :class:`SchemeMorphism_toric_coordinates_field`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1(1,2,3,4)
        [1 : 2 : 3 : 4]
    """
    # Mimicking affine/projective classes
    def __init__(self, X, coordinates, check=True):
        r"""
        See :class:`SchemeMorphism_toric_coordinates_field` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1(1,2,3,4)
            [1 : 2 : 3 : 4]
        """
        # Convert scheme to its set of points over the base ring
        if is_Scheme(X):
            X = X(X.base_ring())
        super(SchemeMorphism_toric_coordinates_field, self).__init__(X)
        if check:
            # Verify that there are the right number of coords
            # Why is it not done in the parent?
            if is_SchemeMorphism(coordinates):
                coordinates = list(coordinates)
            if not isinstance(coordinates, (list, tuple)):
                raise TypeError("coordinates must be a scheme point, list, "
                                "or tuple. Got %s" % coordinates)
            d = X.codomain().ambient_space().ngens()
            if len(coordinates) != d:
                raise ValueError("there must be %d coordinates! Got only %d: "
                                 "%s" % (d, len(coordinates), coordinates))
            # Make sure the coordinates all lie in the appropriate ring
            coordinates = Sequence(coordinates, X.value_ring())
            # Verify that the point satisfies the equations of X.
            X.codomain()._check_satisfies_equations(coordinates)
        self._coords = coordinates


class SchemeMorphism_on_points_toric_variety(SchemeMorphism_on_points):
    """
    Construct a morphism determined by action on points.

    .. WARNING::

        You should not create objects of this class directly.

    INPUT:

    - same as for :class:`sage.schemes.generic.SchemeMorphism_on_points`.

    OUPUT:

    - :class:`SchemeMorphism_on_points_toric_variety`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1.inject_variables()
        Defining z0, z1, z2, z3
        sage: P1 = P1xP1.subscheme(z0-z2)
        sage: H = P1xP1.Hom(P1)
        sage: import sage.schemes.generic.toric_variety as tv
        sage: tv.SchemeMorphism_on_points_toric_variety(H, [z0,z1,z0,z3])
        Scheme morphism:
          From: 2-d toric variety covered by 4 affine patches
          To:   Closed subscheme of 2-d toric variety
                covered by 4 affine patches defined by:
          z0 - z2
          Defn: Defined on coordinates by sending
                [z0 : z1 : z2 : z3] to [z0 : z1 : z0 : z3]
    """

    def __init__(self, parent, polynomials, check=True):
        r"""
        See :class:`SchemeMorphism_on_points_toric_variety` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: H = P1xP1.Hom(P1)
            sage: import sage.schemes.generic.toric_variety as tv
            sage: tv.SchemeMorphism_on_points_toric_variety(H, [z0,z1,z0,z3])
            Scheme morphism:
              From: 2-d toric variety covered by 4 affine patches
              To:   Closed subscheme of 2-d toric variety
                    covered by 4 affine patches defined by:
              z0 - z2
              Defn: Defined on coordinates by sending
                    [z0 : z1 : z2 : z3] to [z0 : z1 : z0 : z3]
        """
        SchemeMorphism_on_points.__init__(self, parent, polynomials, check)
        if check:
            # Check that defining polynomials are homogeneous (degrees can be
            # different if the target uses weighted coordinates)
            for p in self.defining_polynomials():
                if not self.domain().ambient_space().is_homogeneous(p):
                    raise ValueError("%s is not homogeneous!" % p)


def normalize_names(names=None, ngens=None, prefix=None, indices=None,
                    return_prefix=False):
    r"""
    Return a list of names in the standard form.

    INPUT:

    All input parameters are optional.

    - ``names`` -- names given either as a single string (with individual
      names separated by commas or spaces) or a list of strings with each
      string specifying a name. If the last name ends with the plus sign,
      "+", this name will be used as ``prefix`` (even if ``prefix`` was
      given explicitly);

    - ``ngens`` -- number of names to be returned;

    - ``prefix`` -- prefix for the indexed names given as a string;

    - ``indices`` -- list of integers (default: ``range(ngens)``) used as
      indices for names with ``prefix``. If given, must be of length
      ``ngens``;

    - ``return_prefix`` -- if ``True``, the last element of the returned list
      will contain the prefix determined from ``names`` or given as the
      parameter ``prefix``. This is useful if you may need more names in the
      future.

    OUTPUT:

    - list of names given as strings.

    These names are constructed in the following way:

    #. If necessary, split ``names`` into separate names.
    #. If the last name ends with "+", put it into ``prefix``.
    #. If ``ngens`` was given, add to the names obtained so far as many
       indexed names as necessary to get this number. If the ``k``-th name of
       the *total* list of names is indexed, it is
       ``prefix + str(indices[k])``. If there were already more names than
       ``ngens``, discard "extra" ones.
    #. Check if constructed names are valid. See :func:`certify_names` for
       details.
    #. If the option ``return_prefix=True`` was given, add ``prefix`` to the
       end of the list.

    EXAMPLES:

    As promised, all parameters are optional::

        sage: from sage.schemes.generic.toric_variety import normalize_names
        sage: normalize_names()
        []

    One of the most common uses is probably this one::

        sage: normalize_names("x+", 4)
        ['x0', 'x1', 'x2', 'x3']

    Now suppose that you want to enumerate your variables starting with one
    instead of zero::

        sage: normalize_names("x+", 4, indices=range(1,5))
        ['x1', 'x2', 'x3', 'x4']

    You may actually have an arbitrary enumeration scheme::

        sage: normalize_names("x+", 4, indices=[1, 10, 100, 1000])
        ['x1', 'x10', 'x100', 'x1000']

    Now let's add some "explicit" names::

        sage: normalize_names("x y z t+", 4)
        ['x', 'y', 'z', 't3']

    Note that the "automatic" name is ``t3`` instead of ``t0``. This may seem
    weird, but the reason for this behaviour is that the fourth name in this
    list will be the same no matter how many explicit names were given::

        sage: normalize_names("x y t+", 4)
        ['x', 'y', 't2', 't3']

    This is especially useful if you get ``names`` from a user but want to
    specify all default names::

        sage: normalize_names("x, y", 4, prefix="t")
        ['x', 'y', 't2', 't3']

    In this format, the user can easily override your choice for automatic
    names::

        sage: normalize_names("x y s+", 4, prefix="t")
        ['x', 'y', 's2', 's3']

    Let's now use all parameters at once::

        sage: normalize_names("x, y, s+", 4, prefix="t",
        ...       indices=range(1,5), return_prefix=True)
        ['x', 'y', 's3', 's4', 's']

    Note that you still need to give indices for all names, even if some of
    the first ones will be "wasted" because of the explicit names. The reason
    is the same as before - this ensures consistency of automatically
    generated names, no matter how many explicit names were given.

    The prefix is discarded if ``ngens`` was not given::

        sage: normalize_names("alpha, beta, gamma, zeta+")
        ['alpha', 'beta', 'gamma']

    Finally, let's take a look at some possible mistakes::

        sage: normalize_names("123")
        Traceback (most recent call last):
        ...
        ValueError: name must start with a letter! Got 123

    A more subtle one::

        sage: normalize_names("x1", 4, prefix="x")
        Traceback (most recent call last):
        ...
        ValueError: names must be distinct! Got: ['x1', 'x1', 'x2', 'x3']
    """
    if names is None:
        names = []
    elif isinstance(names, str):
        names = names.replace(",", " ").split()
    else:
        try:
            names = list(names)
        except TypeError:
            raise TypeError(
                    "names must be a string or a list or tuple of them!")
        for name in names:
            if not isinstance(name, str):
                raise TypeError(
                    "names must be a string or a list or tuple of them!")
    if names and names[-1].endswith("+"):
        prefix = names.pop()[:-1]
    if ngens is None:
        ngens = len(names)
    if len(names) < ngens:
        if prefix is None:
            raise IndexError("need %d names but only %d are given!"
                             % (ngens, len(names)))
        if indices is None:
            indices = range(ngens)
        elif len(indices) != ngens:
            raise ValueError("need exactly %d indices, but got %d!"
                             % (ngens, len(indices)))
        names += [prefix + str(i) for i in indices[len(names):]]
    if len(names) > ngens:
        names = names[:ngens]
    # Check that all given and constructed names are valid
    certify_names(names)
    if return_prefix:
        names.append(prefix)
    return names


def certify_names(names):
    r"""
    Make sure that ``names`` are valid in Python.

    INPUT:

    - ``names`` -- list of strings.

    OUTPUT:

    - none, but a ``ValueError`` exception is raised if ``names`` are invalid.

    Each name must satisfy the following requirements:

    * Be non-empty.
    * Contain only (Latin) letters, digits, and underscores ("_").
    * Start with a letter.

    In addition, all names must be distinct.

    EXAMPLES::

        sage: from sage.schemes.generic.toric_variety import certify_names
        sage: certify_names([])
        sage: certify_names(["a", "x0", "x_45"])
        sage: certify_names(["", "x0", "x_45"])
        Traceback (most recent call last):
        ...
        ValueError: name must be nonempty!
        sage: certify_names(["a", "0", "x_45"])
        Traceback (most recent call last):
        ...
        ValueError: name must start with a letter! Got 0
        sage: certify_names(["a", "x0", "Щ_45"])
        Traceback (most recent call last):
        ...
        ValueError: name must be alphanumeric! Got Щ_45
        sage: certify_names(["a", "x0", "x0"])
        Traceback (most recent call last):
        ...
        ValueError: names must be distinct! Got: ['a', 'x0', 'x0']
    """
    for name in names:
        if not name:
            raise ValueError("name must be nonempty!")
        if not name.isalnum() and not name.replace("_","").isalnum():
            # Must be alphanumeric except for non-leading '_'
            raise ValueError("name must be alphanumeric! Got %s" % name)
        if not name[0].isalpha():
            raise ValueError("name must start with a letter! Got %s" % name)
    if len(set(names)) != len(names):
        raise ValueError("names must be distinct! Got: %s" % names)
