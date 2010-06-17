r"""
Rational polyhedral fans

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.generic.toric_variety`,
:mod:`~sage.schemes.generic.fano_toric_variety`). While the emphasis is on
complete full-dimensional fans, arbitrary fans are supported. Work
with distinct lattices. The default lattice is :class:`ToricLattice
<sage.geometry.toric_lattice.ToricLatticeFactory>` `N` of the appropriate
dimension. The only case when you must specify lattice explicitly is creation
of a 0-dimensional fan, where dimension of the ambient space cannot be
guessed.

A **rational polyhedral fan** is a *finite* collection of *strictly* convex
rational polyhedral cones, such that the intersection of any two cones of the
fan is a face of each of them and each face of each cone is also a cone of the
fan.

AUTHORS:

- Andrey Novoseltsev (2010-05-15): initial version.

- Andrey Novoseltsev (2010-06-17): substantial improvement during review by
  Volker Braun.

EXAMPLES:

Use :func:`Fan` to construct fans "explicitly"::

    sage: fan = Fan(cones=[(0,1), (1,2)],
    ...             rays=[(1,0), (0,1), (-1,0)])
    sage: fan
    Rational polyhedral fan in 2-d lattice N

In addition to giving such lists of cones and rays you can also create cones
first using :func:`~sage.geometry.cone.Cone` and then combine them into a fan.
See the documentation of :func:`Fan` for details.

Instead of building a fan from scratch, for this tutorial we will use an easy
way to get two fans assosiated to :class:`lattice polytopes
<sage.geometry.lattice_polytope.LatticePolytopeClass>`: :func:`FaceFan` and
:func:`NormalFan`::

    sage: fan1 = FaceFan(lattice_polytope.octahedron(3))
    sage: fan2 = NormalFan(lattice_polytope.octahedron(3))

Given such "automatic" fans, you may wonder what are their rays and cones::

    sage: fan1.rays()
    (N(1, 0, 0), N(0, 1, 0), N(0, 0, 1),
     N(-1, 0, 0), N(0, -1, 0), N(0, 0, -1))
    sage: fan1.ray_matrix()
    [ 1  0  0 -1  0  0]
    [ 0  1  0  0 -1  0]
    [ 0  0  1  0  0 -1]
    sage: fan1.generating_cones()
    (3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N)

The last output is not very illuminating. Let's try to improve it::

    sage: for cone in fan1: print cone.rays()
    (N(1, 0, 0), N(0, 1, 0), N(0, 0, -1))
    (N(0, 1, 0), N(-1, 0, 0), N(0, 0, -1))
    (N(1, 0, 0), N(0, -1, 0), N(0, 0, -1))
    (N(-1, 0, 0), N(0, -1, 0), N(0, 0, -1))
    (N(1, 0, 0), N(0, 1, 0), N(0, 0, 1))
    (N(0, 1, 0), N(0, 0, 1), N(-1, 0, 0))
    (N(1, 0, 0), N(0, 0, 1), N(0, -1, 0))
    (N(0, 0, 1), N(-1, 0, 0), N(0, -1, 0))

You can also do ::

    sage: for cone in fan1: print cone.ambient_ray_indices()
    (0, 1, 5)
    (1, 3, 5)
    (0, 4, 5)
    (3, 4, 5)
    (0, 1, 2)
    (1, 2, 3)
    (0, 2, 4)
    (2, 3, 4)

to see indices of rays of the fan corresponding to each cone.

While the above cycles were over "cones in fan", it is obvious that we did not
get ALL the cones: every face of every cone in a fan must also be in the fan,
but all of the above cones were of dimension three. The reason for this
behaviour is that in many cases it is enough to work with generating cones of
the fan, i.e. cones which are not faces of bigger cones. When you do need to
work with lower dimensional cones, you can easily get access to them using
:meth:`~sage.geometry.fan.RationalPolyhedralFan.cones`::

    sage: [cone.ambient_ray_indices() for cone in fan1.cones(2)]
    [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4),
     (2, 4), (3, 4), (3, 5), (4, 5), (0, 5), (1, 5)]

In fact, you don't have to type ``.cones``::

    sage: [cone.ambient_ray_indices() for cone in fan1(2)]
    [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (0, 4),
     (2, 4), (3, 4), (3, 5), (4, 5), (0, 5), (1, 5)]

You may also need to know the inclusion relations between all of the cones of
the fan. In this case check out
:meth:`~sage.geometry.fan.RationalPolyhedralFan.cone_lattice`::

    sage: L = fan1.cone_lattice()
    sage: L
    Finite poset containing 28 elements
    sage: L.bottom()
    0-d cone of Rational polyhedral fan in 3-d lattice N
    sage: L.top()
    Rational polyhedral fan in 3-d lattice N
    sage: cone = L.level_sets()[2][0]
    sage: cone
    2-d cone of Rational polyhedral fan in 3-d lattice N
    sage: L.hasse_diagram().neighbors(cone)
    [3-d cone of Rational polyhedral fan in 3-d lattice N,
     3-d cone of Rational polyhedral fan in 3-d lattice N,
     1-d cone of Rational polyhedral fan in 3-d lattice N,
     1-d cone of Rational polyhedral fan in 3-d lattice N]

Note, that while ``cone`` above seems to be a "cone", it is not::

    sage: cone.rays()
    Traceback (most recent call last):
    ...
    AttributeError: 'PosetElement' object has no attribute 'rays'

To get your hands on the "real" cone, you need to do one more step::

    sage: cone = cone.element
    sage: cone.rays()
    (N(1, 0, 0), N(0, 1, 0))

You can check how "good" a fan is::

    sage: fan1.is_complete()
    True
    sage: fan1.is_simplicial()
    True
    sage: fan1.is_smooth()
    True

The face fan of the octahedron is really good! Time to remember that we have
also constructed its normal fan::

    sage: fan2.is_complete()
    True
    sage: fan2.is_simplicial()
    False
    sage: fan2.is_smooth()
    False

This one does have some "problems," but we can fix them::

    sage: fan3 = fan2.make_simplicial()
    sage: fan3.is_simplicial()
    True
    sage: fan3.is_smooth()
    False

Note that we had to save the result of
:meth:`~sage.geometry.fan.RationalPolyhedralFan.make_simplicial` in a new fan.
Fans in Sage are immutable, so any operation that does change them constructs
a new fan.

We can also make ``fan3`` smooth, but it will take a bit more work::

    sage: cube = lattice_polytope.octahedron(3).polar()
    sage: sk = cube.skeleton_points(2)
    sage: rays = [cube.point(p) for p in sk]
    sage: fan4 = fan3.subdivide(new_rays=rays)
    sage: fan4.is_smooth()
    True

Let's see how "different" are ``fan2`` and ``fan4``::

    sage: fan2.ngenerating_cones()
    6
    sage: fan2.nrays()
    8
    sage: fan4.ngenerating_cones()
    48
    sage: fan4.nrays()
    26

Smoothness does not come for free!

Please take a look at the rest of the available functions below and their
complete descriptions. If you need any features that are missing, feel free to
suggest them. (Or implement them on your own and submit a patch to Sage for
inclusion!)
"""


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import collections
import warnings

from sage.combinat.posets.posets import FinitePoset
from sage.geometry.cone import (Cone,
                                ConvexRationalPolyhedralCone,
                                IntegralRayCollection,
                                hasse_diagram_from_incidences,
                                is_Cone,
                                normalize_rays)
from sage.geometry.lattice_polytope import (LatticePolytope,
                                            all_faces,
                                            all_facet_equations)
from sage.geometry.toric_lattice import is_ToricLattice
from sage.graphs.all import DiGraph
from sage.matrix.all import matrix
from sage.misc.all import flatten, walltime
from sage.rings.all import QQ, RR, ZZ
from sage.structure.all import Sequence
from sage.structure.coerce import parent


def is_Fan(x):
    r"""
    Check if ``x`` is a Fan.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a fan and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.fan import is_Fan
        sage: is_Fan(1)
        False
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: fan
        Rational polyhedral fan in 2-d lattice N
        sage: is_Fan(fan)
        True
    """
    return isinstance(x, RationalPolyhedralFan)


def Fan(cones, rays=None, lattice=None, check=True, normalize=True,
        is_complete=None):
    r"""
    Construct a rational polyhedral fan.

    .. NOTE::

        Approximate time to construct a fan consisting of `n` cones is `n^2/5`
        seconds. That is half an hour for 100 cones. This time can be
        significantly reduced in the future, but it is still likely to be
        `\sim n^2` (with, say, `/500` instead of `/5`). If you know that your
        input does form a valid fan, use ``check=False`` option to skip
        consistency checks.

    INPUT:

    - ``cones`` -- list of either
      :class:`Cone<sage.geometry.cone.ConvexRationalPolyhedralCone>` objects
      or lists of integers interpreted as indices of generating rays in
      ``rays``;

    - ``rays`` -- list of rays given as list or vectors convertible to the
      rational extension of ``lattice``. If ``cones`` are given by
      :class:`Cone<sage.geometry.cone.ConvexRationalPolyhedralCone>` objects
      ``rays`` may be determined automatically. You still may give them
      explicitly to ensure a particular order of rays in the fan. In this case
      you must list all rays that appear in ``cones``. You can give "extra"
      ones if it is convenient (e.g. if you have a big list of rays for
      several fans), but all "extra" rays will be discarded;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically;

    - ``check`` -- by default the input data will be checked for correctness
      (e.g. that intersection of any two given cones is a face of each) and
      generating cones will be selected. If you know for sure that the input
      is the set of generating cones of a fan, you may significantly
      decrease construction time using ``check=False`` option;

    - ``normalize`` -- you can further speed up construction using
      ``normalize=False`` option. In this case ``cones`` must be a list of
      **sorted** :class:`tuples` and ``rays`` must be immutable primitive
      vectors in ``lattice``. In general, you should not use this option, it
      is designed for code optimization and does not give as drastic
      improvement in speed as the previous one;

    - ``is_complete`` -- every fan can determine on its own if it is complete
      or not, however it can take quite a bit of time for "big" fans with many
      generating cones. On the other hand, in some situations it is known in
      advance that a certain fan is complete. In this case you can pass
      ``is_complete=True`` option to speed up some computations. You may also
      pass ``is_complete=False`` option, although it is less likely to be
      beneficial. Of course, passing a wrong value can compromise the
      integrity of data structures of the fan and lead to wrong results, so
      you should be very careful if you decide to use this option.

    OUTPUT:

    - :class:`rational polyhedral fan <RationalPolyhedralFan>`.

    EXAMPLES:

    Let's construct a fan corresponding to the projective plane in several
    ways::

        sage: cone1 = Cone([(1,0), (0,1)])
        sage: cone2 = Cone([(0,1), (-1,-1)])
        sage: cone3 = Cone([(-1,-1), (1,0)])
        sage: P2 = Fan([cone1, cone2, cone2]) # random output (warning)
        sage: P2.ngenerating_cones()
        2

    Oops! There was a typo and the second ``cone2`` got discarded (you will
    actually see a warning the first time you do such a thing yourself, i.e.
    provide a non-minimal generating set of cones of a fan). Let's fix it::

        sage: P2 = Fan([cone1, cone2, cone3])
        sage: P2.ngenerating_cones()
        3

    Looks better. An alternative way is ::

        sage: rays = [(1,0), (0,1), (-1,-1)]
        sage: cones = [(0,1), (1,2), (2,0)]
        sage: P2a = Fan(cones, rays)
        sage: P2a.ngenerating_cones()
        3
        sage: P2 == P2a
        False

    That may seem wrong, but it is not::

        sage: P2.is_equivalent(P2a)
        True

    See :meth:`~RationalPolyhedralFan.is_equivalent` for details.

    Yet another way to construct this fan is ::

        sage: P2b = Fan(cones, rays, check=False)
        sage: P2b.ngenerating_cones()
        3
        sage: P2a == P2b
        True

    If you try the above examples, you are likely to notice the difference in
    speed, so when you are sure that everything is correct, it is a good idea
    to use ``check=False`` option. On the other hand, it is usually **NOT** a
    good idea to use ``normalize=False`` option::

        sage: P2c = Fan(cones, rays, check=False, normalize=False)
        Traceback (most recent call last):
        ...
        AttributeError: 'tuple' object has no attribute 'parent'

    Yet another way is to use functions :func:`FaceFan` and :func:`NormalFan`
    to construct fans from :class:`lattice polytopes
    <sage.geometry.lattice_polytope.LatticePolytopeClass>`.

    We have not yet used ``lattice`` argument, since if was determined
    automatically::

        sage: P2.lattice()
        2-d lattice N
        sage: P2b.lattice()
        2-d lattice N

    However, it is necessary to specify it explicitly if you want to construct
    a fan without rays or cones::

        sage: Fan([], [])
        Traceback (most recent call last):
        ...
        ValueError: you must specify the lattice
        when you construct a fan without rays and cones!
        sage: F = Fan([], [], lattice=ToricLattice(2, "L"))
        sage: F
        Rational polyhedral fan in 2-d lattice L
        sage: F.lattice_dim()
        2
        sage: F.dim()
        0
    """
    if not check and not normalize:
        return RationalPolyhedralFan(cones, rays, lattice, is_complete)
    rays = normalize_rays(rays, lattice)
    if not isinstance(cones, list):
        try:
            cones = list(cones)
        except TypeError:
            raise TypeError(
                "cones must be given as an iterable!"
                "\nGot: %s" % cones)
    if not cones:
        if not rays and lattice is None:
            raise ValueError("you must specify the lattice when you "
                             "construct a fan without rays and cones!")
        return RationalPolyhedralFan(
                                (), (), rays[0].parent() if rays else lattice)
    if is_Cone(cones[0]):
        # Construct the fan from Cone objects
        if check:
            for cone in cones:
                if not cone.is_strictly_convex():
                    raise ValueError(
                                    "cones of a fan must be strictly convex!")
        ray_set = set([])
        for cone in cones:
            ray_set.update(cone.rays())
        if rays:    # Preserve the initial order of rays, if they were given
            new_rays = []
            for ray in rays:
                if ray in ray_set and ray not in new_rays:
                    new_rays.append(ray)
            if len(new_rays) != len(ray_set):
                raise ValueError(
                  "if rays are given, they must include all rays of the fan!")
            rays = new_rays
        else:
            rays = tuple(ray_set)
        if check:
            # Maybe we should compute all faces of all cones and save them for
            # later if we are doing this check?
            # Should we try to keep the order of cones when we select
            # generating ones? Sorting by dimension is necessary here.
            cones.sort(key=lambda cone: cone.dim(), reverse=True)
            generating_cones = []
            for cone in cones:
                is_generating = True
                for g_cone in generating_cones:
                    i_cone = cone.intersection(g_cone)
                    if i_cone.is_face_of(cone) and i_cone.is_face_of(g_cone):
                        if i_cone.dim() == cone.dim():
                            is_generating = False # cone is a face of g_cone
                            break
                    else:
                        raise ValueError(
                                "these cones cannot belong to the same fan!"
                                "\nCone 1 rays: %s\nCone 2 rays: %s"
                                % (g_cone.rays(), cone.rays()))
                if is_generating:
                    generating_cones.append(cone)
            if len(cones) > len(generating_cones):
                warnings.warn("you have provided a non-minimal set of "
                    "generating cones, %d of them were discarded!"
                    % (len(cones) - len(generating_cones)),
                    stacklevel=2)
            cones = generating_cones
        return RationalPolyhedralFan(
            (tuple(sorted(rays.index(ray) for ray in cone.rays()))
            for cone in cones), rays, lattice, is_complete)
    # Construct the fan from rays and "tuple cones"
    for n, cone in enumerate(cones):
        try:
            cones[n] = sorted(cone)
        except TypeError:
            raise TypeError("since rays are given, cones must be iterables!"
                            "\nGot: %s" % cone)
    if not check:
        return RationalPolyhedralFan(cones, rays, lattice, is_complete)
    # If we do need to make all the check, build explicit cone objects first
    return Fan((Cone([rays[n] for n in cone], lattice) for cone in cones),
               rays, lattice, is_complete=is_complete)


def FaceFan(polytope, lattice=None):
    r"""
    Construct the face fan of the given lattice ``polytope``.

    INPUT:

    - ``polytope`` -- :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically.

    OUTPUT:

    - :class:`rational polyhedral fan <RationalPolyhedralFan>`.

    See also :func:`NormalFan`.

    EXAMPLES:

    Let's construct the fan corresponding to the product of two projective
    lines::

        sage: diamond = lattice_polytope.octahedron(2)
        sage: P1xP1 = FaceFan(diamond)
        sage: P1xP1.ray_matrix()
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        sage: for cone in P1xP1: print cone.rays()
        (N(1, 0), N(0, -1))
        (N(-1, 0), N(0, -1))
        (N(1, 0), N(0, 1))
        (N(0, 1), N(-1, 0))
     """
    if any(d <= 0 for d in polytope.distances([0]*polytope.dim())):
        raise ValueError("face fans are defined only for polytopes containing"
                         "the origin as an interior point!")
    cones = (facet.vertices() for facet in polytope.facets())
    rays = polytope.vertices().columns(copy=False)
    fan = Fan(cones, rays, lattice=lattice, check=False)
    fan._is_complete = polytope.dim() == polytope.ambient_dim()
    return fan


def NormalFan(polytope, lattice=None):
    r"""
    Construct the normal fan of the given lattice ``polytope``.

    INPUT:

    - ``polytope`` -- :class:`lattice polytope
      <sage.geometry.lattice_polytope.LatticePolytopeClass>`;

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If not specified, an attempt will
      be made to determine an appropriate toric lattice automatically.

    OUTPUT:

    - :class:`rational polyhedral fan <RationalPolyhedralFan>`.

    See also :func:`FaceFan`.

    EXAMPLES:

    Let's construct the fan corresponding to the product of two projective
    lines::

        sage: square = lattice_polytope.octahedron(2).polar()
        sage: P1xP1 = NormalFan(square)
        sage: P1xP1.ray_matrix()
        [ 1  0 -1  0]
        [ 0  1  0 -1]
        sage: for cone in P1xP1: print cone.rays()
        (N(1, 0), N(0, -1))
        (N(-1, 0), N(0, -1))
        (N(1, 0), N(0, 1))
        (N(0, 1), N(-1, 0))
    """
    rays = (polytope.facet_normal(i) for i in range(polytope.nfacets()))
    cones = (vertex.facets() for vertex in polytope.faces(dim=0))
    fan = Fan(cones, rays, lattice=lattice, check=False)
    fan._is_complete = polytope.dim() == polytope.ambient_dim()
    return fan


class Cone_of_fan(ConvexRationalPolyhedralCone):
    r"""
    Construct a cone belonging to a fan.

    .. WARNING::

        This class does not check that the input defines a valid cone of a
        fan. You should not construct objects of this class directly.

    In addition to all of the properties of "regular" :class:`cones
    <sage.geometry.cone.ConvexRationalPolyhedralCone>`, such cones know their
    relation to the fan.

    INPUT:

    - ``ambient`` -- fan whose cone is constructed;

    - ``ambient_ray_indices`` -- increasing list or tuple of integers, indices
      of rays of ``ambient`` generating this cone.

    OUTPUT:

    - cone of ``ambient``.

    EXAMPLES:

    The intended way to get objects of this class is the following::

        sage: P1xP1 = FaceFan(lattice_polytope.octahedron(2))
        sage: cone = P1xP1.generating_cone(0)
        sage: cone
        2-d cone of Rational polyhedral fan in 2-d lattice N
        sage: cone.ambient_ray_indices()
        (0, 3)
        sage: cone.star_generator_indices()
        (0,)
    """

    def __init__(self, ambient, ambient_ray_indices):
        r"""
        See :class:`Cone_of_Fan` for documentation.

        TESTS:

        The following code is likely to construct an invalid object, we just
        test that creation of cones of fans is working::

            sage: P1xP1 = FaceFan(lattice_polytope.octahedron(2))
            sage: cone = sage.geometry.fan.Cone_of_fan(P1xP1, (0,))
            sage: cone
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: TestSuite(cone).run()
        """
        super(Cone_of_fan, self).__init__(
                    ambient=ambient, ambient_ray_indices=ambient_ray_indices)
        self._is_strictly_convex = True
        # Because if not, this cone should not have been constructed

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: P1xP1 = FaceFan(lattice_polytope.octahedron(2))
            sage: cone = P1xP1.generating_cone(0)
            sage: cone._repr_()
            '2-d cone of Rational polyhedral fan in 2-d lattice N'
            sage: cone.facets()[0]._repr_()
            '1-d cone of Rational polyhedral fan in 2-d lattice N'
        """
        # The base class would print "face of" instead of  "cone of"
        return "%d-d cone of %s" % (self.dim(), self.ambient())

    def star_generator_indices(self):
        r"""
        Return indices of generating cones of the "ambient fan" containing
        ``self``.

        OUTPUT:

        - increasing :class:`tuple` of integers.

        EXAMPLES::

            sage: P1xP1 = FaceFan(lattice_polytope.octahedron(2))
            sage: cone = P1xP1.generating_cone(0)
            sage: cone.star_generator_indices()
            (0,)
        """
        if "_star_generator_indices" not in self.__dict__:
            fan = self.ambient()
            sgi = set(range(fan.nrays()))
            for ray in self.ambient_ray_indices():
                sgi.intersection_update(fan._ray_to_cones(ray))
            self._star_generatori_indices = tuple(sorted(sgi))
        return self._star_generator_indices

    def star_generators(self):
        r"""
        Return indices of generating cones of the "ambient fan" containing
        ``self``.

        OUTPUT:

        - increasing :class:`tuple` of integers.

        EXAMPLES::

            sage: P1xP1 = FaceFan(lattice_polytope.octahedron(2))
            sage: cone = P1xP1.generating_cone(0)
            sage: cone.star_generators()
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
        """
        if "_star_generators" not in self.__dict__:
            self._star_generators = tuple(self.ambient().generating_cone(i)
                                    for i in self.star_generator_indices())
        return self._star_generators


class RationalPolyhedralFan(IntegralRayCollection,
                            collections.Callable,
                            collections.Container):
    r"""
    Create a rational polyhedral fan.

    .. WARNING::

        This class does not perform any checks of correctness of input nor
        does it convert input into the standard representation. Use
        :func:`Fan` to construct fans from "raw data" or :func:`FaceFan` and
        :func:`NormalFan` to get fans associated to polytopes.

    Fans are immutable, but they cache most of the returned values.

    INPUT:

    - ``cones`` -- list of generating cones of the fan, each cone given as a
      list of indices of its generating rays in ``rays``;

    - ``rays`` -- list of immutable primitive vectors in ``lattice``
      consisting of exactly the rays of the fan (i.e. no "extra" ones);

    - ``lattice`` -- :class:`ToricLattice
      <sage.geometry.toric_lattice.ToricLatticeFactory>`, `\ZZ^n`, or any
      other object that behaves like these. If ``None``, it will be determined
      as :func:`parent` of the first ray. Of course, this cannot be done if
      there are no rays, so in this case you must give an appropriate
      ``lattice`` directly.

    - ``is_complete`` -- if given, must be ``True`` or ``False`` depending on
      whether this fan is complete or not. By default, it will be determined
      automatically if necessary.

    OUTPUT:

    - rational polyhedral fan.
    """

    def __init__(self, cones, rays, lattice, is_complete=None):
        r"""
        See :class:`RationalPolyhedralFan` for documentation.

        TESTS::

            sage: v = vector([0,1])
            sage: v.set_immutable()
            sage: f = sage.geometry.fan.RationalPolyhedralFan(
            ...                         [(0,)], [v], None)
            sage: f.rays()
            ((0, 1),)
            sage: TestSuite(f).run()
            sage: f = Fan([(0,)], [(0,1)])
            sage: TestSuite(f).run()
        """
        super(RationalPolyhedralFan, self).__init__(rays, lattice)
        self._generating_cones = tuple(Cone_of_fan(self, c) for c in cones)
        for i, cone in enumerate(self._generating_cones):
            cone._star_generator_indices = (i,)
        # Knowing completeness drastically affects the speed of cone lattice
        # computation and containment check, so we have a special way to
        # optimize it.
        if is_complete is not None:
            self._is_complete = is_complete

    def __call__(self, dim=None, codim=None):
        r"""
        Return the specified cones of ``self``.

        .. NOTE::

            "Direct call" syntax is a synonym for :meth:`cones` method.

        INPUT:

        - ``dim`` -- dimension of the requested cones;

        - ``codim`` -- codimension of the requested cones.

        OUTPUT:

        - cones of ``self`` of the specified (co)dimension.

        TESTS::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan(1)
            (1-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N)
            sage: fan(2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(dim=2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(codim=2)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(dim=1, codim=1)
            Traceback (most recent call last):
            ...
            ValueError: dimension and codimension
            cannot be specified together!
        """
        return self.cones(dim, codim)

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is of the same type as ``self``, their rays are
          the same and stored in the same order, and their generating cones
          are the same and stored in the same order. 1 or -1 otherwise.

        TESTS::

            sage: f1 = Fan(cones=[(0,1), (1,2)],
            ...            rays=[(1,0), (0,1), (-1, 0)],
            ...            check=False)
            sage: f2 = Fan(cones=[(1,2), (0,1)],
            ...            rays=[(1,0), (0,1), (-1, 0)],
            ...            check=False)
            sage: f3 = Fan(cones=[(1,2), (0,1)],
            ...            rays=[(1,0), (0,1), (-1, 0)],
            ...            check=False)
            sage: cmp(f1, f2)
            1
            sage: cmp(f2, f1)
            -1
            sage: cmp(f2, f3)
            0
            sage: f2 is f3
            False
            sage: cmp(f1, 1)
            -1
        """
        c = cmp(type(self), type(right))
        if c:
            return c
        return cmp([self.rays(), self.generating_cones()],
                   [right.rays(), right.generating_cones()])

    def __contains__(self, data):
        r"""
        Check if ``data`` is contained in ``self``.

        See :meth:`_contains` (which is called by this function) for
        documentation.

        TESTS::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: cone1 in f
            True
            sage: (1,1) in f
            True
            sage: (-1,-1) in f
            False
        """
        return self._contains(data)

    def __iter__(self):
        r"""
        Return an iterator over generating cones of ``self``.

        OUTPUT:

        -  iterator.

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ...           rays=[(1,0), (0,1), (-1, 0)],
            ...           check=False)
            sage: for cone in f: print cone.rays()
            (N(1, 0), N(0, 1))
            (N(0, 1), N(-1, 0))
         """
        return iter(self.generating_cones())

    def _compute_cone_lattice(self):
        r"""
        Compute the cone lattice of ``self``.

        See :meth:`cone_lattice` for documentation.

        TESTS:

        We use different algorithms depending on available information. One of
        the common cases is a fan which is KNOWN to be complete, i.e. we do
        not even need to check if it is complete.

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.cone_lattice() # indirect doctest
            Finite poset containing 10 elements

        These 10 elements are: 1 origin, 4 rays, 4 generating cones, 1 fan.

        Another common case is the fan of faces of a single cone::

            sage: quadrant = Cone([(1,0), (0,1)])
            sage: fan = Fan([quadrant])
            sage: fan.cone_lattice() # indirect doctest
            Finite poset containing 5 elements

        These 5 elements are: 1 origin, 2 rays, 1 generating cone, 1 fan.

        Finally, we have "intermediate" fans which are incomplete but are
        generated by more than one cone::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.ray_matrix()
            [ 0  1 -1]
            [ 1  0  0]
            sage: for cone in fan: print cone.ambient_ray_indices()
            (0, 1)
            (2,)
            sage: L = fan.cone_lattice() # indirect doctest
            sage: L
            Finite poset containing 6 elements

        Here we got 1 origin, 3 rays (one is a generating cone),
        1 2-dimensional cone (a generating one), and 1 fan.
        """
        # Define a face constructor
        def FanFace(rays, cones):
            if not cones:       # The top face, fan itself
                return self
            if len(cones) == 1: # MAY be a generating cone or NOT!!!
                g_cone = self.generating_cone(cones[0])
                if g_cone.ambient_ray_indices() == rays:
                    return g_cone
            face = Cone_of_fan(ambient=self, ambient_ray_indices=rays)
            face._star_generator_indices=cones
            return face
        # Check directly if we know completeness already, since *determining*
        # completeness relies on this function
        if "_is_complete" in self.__dict__ and self._is_complete:
            # We can use a fast way for complete fans
            self._cone_lattice = hasse_diagram_from_incidences(
                                self._ray_to_cones(),
                                (cone.ambient_ray_indices() for cone in self),
                                FanFace)
        else:
            # For general fans we will "merge" face lattices of generating
            # cones.
            L = DiGraph()
            face_to_rays = dict() # face |---> (indices of fan rays)
            rays_to_index = dict() # (indices of fan rays) |---> face index
            # face index |---> (indices of containing generating cones)
            index_to_cones = []
            # During construction index 0 will correspond to the fan
            # We think of the fan not being in the cone even when there is
            # only one cone
            index_to_cones.append(())
            next_index = 1
            for i, cone in enumerate(self):
                # Set up translation of faces of cone to rays and indices
                # We make a standalone cone to compute its standalone face
                # lattice, since cones of fans get their lattices from fans
                L_cone = Cone(cone.rays(), lattice=self.lattice(),
                              check=False, normalize=False).face_lattice()
                for f in L_cone:
                    f = f.element
                    f_rays = tuple(cone.ambient_ray_indices()[ray]
                                   for ray in f.ambient_ray_indices())
                    face_to_rays[f] = f_rays
                    try:
                        f_index = rays_to_index[f_rays]
                        index_to_cones[f_index].append(i)
                    except KeyError:        # Did not see f before
                        f_index = next_index
                        next_index += 1
                        rays_to_index[f_rays] = f_index
                        index_to_cones.append([i])
                # Add all relations between faces of cone to L
                for f,g in L_cone.cover_relations_iterator():
                    L.add_edge(rays_to_index[face_to_rays[f.element]],
                               rays_to_index[face_to_rays[g.element]])
                # Add the inclusion of cone into the fan itself
                L.add_edge(
                        rays_to_index[face_to_rays[L_cone.top().element]], 0)

            # Enumeration of graph vertices must be a linear extension of the
            # poset
            new_order = L.topological_sort()
            # Make sure that generating cones are in the end in proper order
            tail = [rays_to_index[gc.ambient_ray_indices()] for gc in self]
            tail.append(0) # We know that the fan itself has index 0
            new_order = [n for n in new_order if n not in tail] + tail
            # Make sure that rays are in the beginning in proper order
            head = [rays_to_index[()]] # Empty face
            head.extend(rays_to_index[(n,)] for n in range(self.nrays()))
            new_order = head + [n for n in new_order if n not in head]
            # "Invert" this list to a dictionary
            labels = dict()
            for new, old in enumerate(new_order):
                labels[old] = new
            L.relabel(labels)

            elements = [None] * next_index
            for rays, index in rays_to_index.items():
                elements[labels[index]] = FanFace(
                                           rays, tuple(index_to_cones[index]))
            # We need "special treatment" for the whole fan. If we added its
            # ray incidence information to the total list, it would be
            # confused with the generating cone in the case of a single cone.
            elements[labels[0]] = FanFace(tuple(range(self.nrays())), ())
            self._cone_lattice = FinitePoset(L, elements)

    def _contains(self, data):
        r"""
        Check if ``data`` is contained in ``self``.

        This function is called by :meth:`__contains__` and :meth:`contains`
        to ensure the same call depth for warning messages.

        INPUT:

        - ``data`` -- anything. If it is not a cone, an attempt will be made
          to convert it into a single point of the ambient space of ``self``.
          If it fails, ``False`` will be returned.

        OUTPUT:

        - ``True`` if ``data`` is contained in ``self``, ``False`` otherwise.

        TESTS::

        TESTS::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: f._contains(cone1)
            True
            sage: f._contains((1,1))
            True
        """
        if is_Cone(data):
            if data.lattice() != self.lattice():
                warnings.warn("you have checked if a fan contains a cone "
                              "from another lattice, this is always False!",
                              stacklevel=3)
                # We may need to take into account canonical maps, perhaps...
                return False
            if (not data.is_strictly_convex()
                or not data.ray_set().issubset(self.ray_set())):
                return False
            try:
                return data.is_equivalent(self.cone_containing(data.rays()))
            except ValueError:  # No cone containing all these rays
                return False
        # Now we try to treat data as a point
        try:
            data = self._ambient_space_point(data)
        except TypeError, ex:
            if str(ex).endswith("have incompatible lattices!"):
                warnings.warn("you have checked if a fan contains a point "
                              "from an incompatible lattice, this is False!",
                              stacklevel=3)
            return False
        if self.is_complete():
            return True
        return any(data in cone for cone in self)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ...           rays=[(1,0), (0,1), (-1, 0)],
            ...           check=False)
            sage: f._latex_()
            '\\Sigma^{2}'
        """
        return r"\Sigma^{%s}" % self.lattice_dim()

    def _ray_to_cones(self, i=None):
        r"""
        Return the set of generating cones containing the ``i``-th ray.

        INPUT:

        - ``i`` -- integer, index of a ray of ``self``.

        OUTPUT:

        - :class:`frozenset` of indices of generating cones of ``self``
          containing the ``i``-th ray if ``i`` was given, :class:`tuple` of
          these sets for all rays otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan._ray_to_cones(0)
            frozenset([0, 2])
            sage: fan._ray_to_cones()
            (frozenset([0, 2]), frozenset([2, 3]),
             frozenset([1, 3]), frozenset([0, 1]))
        """
        # This function is close to self(1)[i].star_generator_indices(), but
        # it does not require computation of the cone lattice and is
        # convenient for internal purposes.
        if "_ray_to_cones_tuple" not in self.__dict__:
            ray_to_cones = []
            for ray in self.rays():
                ray_to_cones.append([])
            for k, cone in enumerate(self):
                for j in cone.ambient_ray_indices():
                    ray_to_cones[j].append(k)
            self._ray_to_cones_tuple = tuple(frozenset(rtc)
                                             for rtc in ray_to_cones)
        if i is None:
            return self._ray_to_cones_tuple
        else:
            return self._ray_to_cones_tuple[i]

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: f = Fan(cones=[(0,1), (1,2)],
            ...           rays=[(1,0), (0,1), (-1, 0)],
            ...           check=False)
            sage: f._repr_()
            'Rational polyhedral fan in 2-d lattice N'
            sage: f = Fan(cones=[(0,1), (1,2)],
            ...           rays=[(1,0), (0,1), (-1, 0)],
            ...           lattice=ZZ^2,
            ...           check=False)
            sage: f._repr_()
            'Rational polyhedral fan in 2-d lattice'
        """
        result = "Rational polyhedral fan in"
        if is_ToricLattice(self.lattice()):
            result += " %s" % self.lattice()
        else:
            result += " %d-d lattice" % self.lattice_dim()
        return result

    def _subdivide_palp(self, new_rays, verbose):
        r"""
        Subdivide ``self`` adding ``new_rays`` one by one.

        INPUT:

        - ``new_rays`` -- immutable primitive vectors in the lattice of
          ``self``;

        - ``verbose`` -- if ``True``, some timing information will be printed.

        OUTPUT:

        - rational polyhedral fan.

        .. NOTE::

            All generating cones of ``self`` must be full-dimensional.

        TESTS::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: new_rays = sage.geometry.cone.normalize_rays([(1,1)], None)
            sage: fan = Fan([cone1, cone2])
            sage: fan._subdivide_palp(new_rays, True)
            Traceback (most recent call last):
            ...
            ValueError: PALP-subdividing can be used only for
            fans whose generating cones are full-dimensional!
            sage: fan = Fan([cone1])
            sage: # Timing information will depend on your machine
            sage: new_fan = fan._subdivide_palp(new_rays, True)
            R:1/1  C:2  T:...(ms)  T/new:...(ms)  T/all:...(ms)
            sage: new_fan.ray_matrix()
            [0 1 1]
            [1 0 1]
            sage: for cone in new_fan: print cone.ambient_ray_indices()
            (0, 2)
            (1, 2)
        """
        dim = self.lattice_dim()
        for cone in self:
            if cone.dim() != dim:
                raise ValueError("PALP-subdividing can be used only for fans "
                            "whose generating cones are full-dimensional!")
        # Convert cones to lattice polytopes
        cone_polytopes = [cone.lattice_polytope() for cone in self]
        for cone_polytope in cone_polytopes:
            cone_polytope._dim = dim
        all_faces(cone_polytopes)
        all_facet_equations(cone_polytopes)
        # Iterative subdivision
        for n, ray in enumerate(new_rays):
            start = walltime()
            old_polytopes = []
            new_polytopes = []
            for cone_polytope in cone_polytopes:
                if (cone_polytope.nvertices() == dim + 1 #simplex
                    and ray in cone_polytope.vertices().columns(copy=False)):
                    old_polytopes.append(cone_polytope)
                    continue # Subdivision would not give any new polytopes
                distances = cone_polytope.distances(ray)
                # The origin is the last 0-dimensional face
                cone_facets = cone_polytope.faces(dim=0)[-1].facets()
                if all(distances[fn] >= 0 for fn in cone_facets):
                    # Ray is inside the cone, even if not inside cone_polytope
                    # Do subdivision with cones over each non-containing facet
                    vertices = cone_polytope.vertices().columns(copy=False)
                    for fn in cone_facets:
                        if distances[fn] > 0:
                            new_v = [vertices[v] for v in
                                        cone_polytope.facets()[fn].vertices()]
                            # Add ray keeping the origin last
                            new_v.insert(-1, ray)
                            new_v = matrix(ZZ, new_v).transpose()
                            new_polytope = LatticePolytope(new_v,
                                copy_vertices=False, compute_vertices=False)
                            new_polytope._dim = dim
                            new_polytopes.append(new_polytope)
                else:
                    old_polytopes.append(cone_polytope)
            # Precompute data for new polytopes using single calls to PALP
            all_faces(new_polytopes)
            all_facet_equations(new_polytopes)
            cone_polytopes = old_polytopes + new_polytopes
            if verbose:
                t = walltime(start)
                # Avoid division by zero
                T_new = ("%d" % (t / len(new_polytopes) * 1000)
                         if new_polytopes else "-")
                print("R:%d/%d  C:%d  T:%d(ms)  T/new:%s(ms)  T/all:%d(ms)"
                    % (n + 1, len(new_rays), len(cone_polytopes), t * 1000,
                       T_new, t / len(cone_polytopes) * 1000))
        # Convert lattice polytopes to cones
        new_fan_rays = list(self.rays())
        new_fan_rays.extend(ray for ray in new_rays
                                if ray not in self.ray_set())
        cones = tuple(tuple(new_fan_rays.index(cone_polytope.vertex(v))
                            for v in range(cone_polytope.nvertices() - 1))
                      for cone_polytope in cone_polytopes)
        fan = Fan(cones, new_fan_rays, check=False, normalize=False)
        # Since we already have all lattice polytopes, let's keep them
        for cone, polytope in zip(fan.generating_cones(), cone_polytopes):
            cone._lattice_polytope = polytope
        return fan

    def cone_containing(self, *points):
        r"""
        Return the smallest cone of ``self`` containing all given points.

        INPUT:

        - either one or more indices of rays of ``self``, or one or more
          objects representing points of the ambient space of ``self``, or a
          list of such objects (you CANNOT give a list of indices).

        OUTPUT:

        - :class:`cone of fan <Cone_of_fan>`.

        .. NOTE::

            We think of the origin as of the smallest cone containing no rays
            at all. If there is no ray in ``self`` that contains all ``rays``,
            a ``ValueError`` exception will be raised.

        EXAMPLES::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])
            sage: f.rays()
            (N(0, 1), N(0, -1), N(1, 0))
            sage: f.cone_containing(0)  # ray index
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing(0, 1) # ray indices
            Traceback (most recent call last):
            ...
            ValueError: there is no cone in
            Rational polyhedral fan in 2-d lattice N
            containing all of the given rays! Rays: (N(0, 1), N(0, -1))
            sage: f.cone_containing(0, 2) # ray indices
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((0,1))  # point
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing([(0,1)]) # point
            1-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((1,1))
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((1,1), (1,0))
            2-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing()
            0-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((0,0))
            0-d cone of Rational polyhedral fan in 2-d lattice N
            sage: f.cone_containing((-1,1))
            Traceback (most recent call last):
            ...
            ValueError: there is no cone in
            Rational polyhedral fan in 2-d lattice N
            containing all of the given points! Points: [N(-1, 1)]

        TESTS::

            sage: fan = Fan(cones=[(0,1,2,3), (0,1,4)],
            ...       rays=[(1,1,1), (1,-1,1), (1,-1,-1), (1,1,-1), (0,0,1)])
            sage: fan.cone_containing(0).rays()
            (N(1, 1, 1),)
        """
        if not points:
            return self.cones(dim=0)[0]
        try:
            rays = map(int, points)
            # Got ray indices
            generating_cones = set(range(self.ngenerating_cones()))
            for ray in rays:
                generating_cones.intersection_update(self._ray_to_cones(ray))
            if not generating_cones:
                raise ValueError("there is no cone in %s containing all of "
                        "the given rays! Rays: %s" % (self, self.rays(rays)))
            containing_cone = self.generating_cone(generating_cones.pop())
            for cone in generating_cones:
                containing_cone = containing_cone.intersection(
                                                self.generating_cone(cone))
            if not self.is_complete():
                # This cone may be too big in the case of incomplete fans
                rays = frozenset(rays)
                facets = containing_cone.facets()
                for facet in facets:
                    if rays.issubset(facet._ambient_ray_indices):
                        containing_cone = containing_cone.intersection(facet)
            return containing_cone
        except TypeError:
            # Got points (hopefully)
            try:
                points = map(self._ambient_space_point, points)
            except TypeError:
                if len(points) == 1:
                    points = map(self._ambient_space_point, points[0])
                else:
                    raise
            # If we are still here, points are good
            # First we try to find a generating cone containing all points
            containing_cone = None
            for cone in self:
                contains_all = True
                for point in points:
                    if point not in cone:
                        contains_all = False
                        break
                if contains_all:
                    containing_cone = cone
                    break
            if containing_cone is None:
                raise ValueError("there is no cone in %s containing all of "
                            "the given points! Points: %s" % (self, points))
            # Now we take the intersection of facets that contain all points
            facets = containing_cone.facets()
            for facet in facets:
                contains_all = True
                for point in points:
                    if point not in facet:
                        contains_all = False
                        break
                if contains_all:
                    containing_cone = containing_cone.intersection(facet)
            return containing_cone

    def cone_lattice(self):
        r"""
        Return the cone lattice of ``self``.

        This lattice will have the origin as the bottom (we do not include the
        empty set as a cone) and the fan itself as the top.

        OUTPUT:

        - :class:`finite poset <sage.combinat.posets.posets.FinitePoset` of
          :class:`cones of fan<Cone_of_fan>`, behaving like "regular" cones,
          but also containing the information about their relation to this
          fan, namely, the contained rays and containing generating cones. The
          top of the lattice will be this fan itself (*which is not a*
          :class:`cone of fan<Cone_of_fan>`).

        See also :meth:`cones`.

        EXAMPLES:

        Cone lattices can be computed for arbitrary fans::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.ray_matrix()
            [ 0  1 -1]
            [ 1  0  0]
            sage: for cone in fan: print cone.ambient_ray_indices()
            (0, 1)
            (2,)
            sage: L = fan.cone_lattice()
            sage: L
            Finite poset containing 6 elements

        These 6 elements are the origin, three rays, one two-dimensional
        cone, and the fan itself\ . Since we do add the fan itself as the
        largest face, you should be a little bit careful with this last
        element::

            sage: for face in L: print face.element.ambient_ray_indices()
            Traceback (most recent call last):
            ...
            AttributeError: 'RationalPolyhedralFan'
            object has no attribute 'ambient_ray_indices'
            sage: L.top()
            Rational polyhedral fan in 2-d lattice N

        For example, you can do ::

            sage: for l in L.level_sets()[:-1]:
            ...       print [f.element.ambient_ray_indices() for f in l]
            [()]
            [(0,), (1,), (2,)]
            [(0, 1)]

        If the fan is complete, its cone lattice is atomic and coatomic and
        can (and will!) be computed in a much more efficient way, but the
        interface is exactly the same::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: L = fan.cone_lattice()
            sage: for l in L.level_sets()[:-1]:
            ...       print [f.element.ambient_ray_indices() for f in l]
            [()]
            [(0,), (1,), (2,), (3,)]
            [(0, 1), (1, 2), (0, 3), (2, 3)]

        Let's also consider the cone lattice of a fan generated by a single
        cone::

            sage: fan = Fan([cone1])
            sage: L = fan.cone_lattice()
            sage: L
            Finite poset containing 5 elements

        Here these 5 elements correspond to the origin, two rays, one
        generating cone of dimension two, and the whole fan. While this single
        cone "is" the whole fan, it is consistent and convenient to
        distinguish them in the cone lattice.
        """
        if "_cone_lattice" not in self.__dict__:
            self._compute_cone_lattice()
        return self._cone_lattice

    # Internally we use this name for a uniform behaviour of cones and fans.
    _face_lattice_function = cone_lattice

    def cones(self, dim=None, codim=None):
        r"""
        Return the specified cones of ``self``.

        INPUT:

        - ``dim`` -- dimension of the requested cones;

        - ``codim`` -- codimension of the requested cones.

        .. NOTE::

            You can specify at most one input parameter.

        OUTPUT:

        - :class:`tuple` of cones of ``self`` of the specified (co)dimension,
          if either ``dim`` or ``codim`` is given. Otherwise :class:`tuple` of
          such tuples for all dimensions.

        EXAMPLES::

            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan(dim=0)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(codim=2)
            (0-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: for cone in fan.cones(1): cone.ray(0)
            N(0, 1)
            N(1, 0)
            N(-1, 0)
            sage: fan.cones(2)
            (2-d cone of Rational polyhedral fan in 2-d lattice N,)
            sage: fan(dim=1, codim=1)
            Traceback (most recent call last):
            ...
            ValueError: dimension and codimension
            cannot be specified together!
        """
        if "_cones" not in self.__dict__:
            levels = [(e.element for e in level) # Generators
                      for level in self.cone_lattice().level_sets()]
            levels.pop() # The very last level is this FAN, not cone.
            # It seems that there is no reason to believe that the order of
            # faces in level sets has anything to do with the order of
            # vertices in the Hasse diagram of FinitePoset. So, while
            # hasse_diagram_from_incidences tried to ensure a "good order,"
            # we will sort faces corresponding to rays, as well as faces
            # corresponding to generating cones, if they are all of the same
            # dimension (otherwise it is not very useful).
            if len(levels) >= 3: # There are cones of dimension higher than 1
                top_cones = list(levels[-1])
                if len(top_cones) == self.ngenerating_cones():
                    top_cones.sort(key=lambda cone:
                                            cone.star_generator_indices()[0])
                levels[-1] = top_cones
            if len(levels) >= 2: # We have rays
                rays = list(levels[1])
                rays.sort(key=lambda cone: cone.ambient_ray_indices()[0])
                levels[1] = rays
            self._cones = tuple(tuple(level) for level in levels)
        if dim is None and codim is None:
            return self._cones
        elif dim is None:
            return self._cones[self.dim() - codim]
        elif codim is None:
            return self._cones[dim]
        raise ValueError(
                    "dimension and codimension cannot be specified together!")

    def contains(self, *args):
        r"""
        Check if a given point or cone is contained in ``self``.

        INPUT:

        - anything. If it is not a cone, an attempt will be made to convert
          the input into a point of the ambient space of ``self``. If it
          fails, ``False`` will be returned.

        OUTPUT:

        - ``True`` if the given point or cone is contained in ``self``,
          ``False`` otherwise.

        .. NOTE::

            A point is contained in a fan if it is contained in its support.
            A cone is contained in a fan if it is equivalent to one of the
            cones of the fan, in particular, it is possible that all points of
            the cone are in the fan, but the cone itself is not.

        EXAMPLES:

        We construct a simple fan::

            sage: cone1 = Cone([(0,-1), (1,0)])
            sage: cone2 = Cone([(1,0), (0,1)])
            sage: f = Fan([cone1, cone2])

        We check if some points are in this fan::

            sage: f.contains(f.lattice()(1,0))
            True
            sage: f.contains((1,0))
            True
            sage: f.contains(1,1)
            True
            sage: f.contains((-1,0))
            False
            sage: f.contains(f.lattice().dual()(1,0)) #random output (warning)
            False
            sage: f.contains(f.lattice().dual()(1,0))
            False
            sage: f.contains(1)
            False
            sage: f.contains(1/2, sqrt(3))
            True
            sage: f.contains(-1/2, sqrt(3))
            False

        Now we check if some cones are in this fan. First, we make sure that
        the order of rays of the input cone does not matter (``check=False``
        option ensures that rays of these cones will be listed exactly as they
        are given)::

            sage: f.contains(Cone([(1,0), (0,1)], check=False))
            True
            sage: f.contains(Cone([(0,1), (1,0)], check=False))
            True

        Now we check that a non-generating cone is in our fan::

            sage: f.contains(Cone([(1,0)]))
            True

        Finally, we test some cones which are not in this fan::

            sage: f.contains(Cone([(1,1)]))
            False
            sage: f.contains(Cone([(0,1), (-1,0)]))
            False
        """
        data = flatten(args)
        if len(data) == 1:
           data = data[0]
        return self._contains(data)

    def gale_transform(self):
        r"""
        Return the Gale transform of ``self``.

        OUTPUT:

        - matrix.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.gale_transform()
            [ 1  0  1  0 -2]
            [ 0  1  0  1 -2]
        """
        if "_gale_transform" not in self.__dict__:
            m = self.ray_matrix().augment(matrix(self.lattice_dim(), 1))
            m = m.stack(matrix([1]*m.ncols()))
            self._gale_transform = m.transpose().integer_kernel().matrix()
        return self._gale_transform

    def generating_cone(self, n):
        r"""
        Return the ``n``-th generating cone of ``self``.

        INPUT:

        - ``n`` -- integer, the index of a generating cone.

        OUTPUT:

        - :class:`cone of fan<Cone_of_fan>`.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.generating_cone(0)
            2-d cone of Rational polyhedral fan in 2-d lattice N
        """
        return self._generating_cones[n]

    def generating_cones(self):
        r"""
        Return generating cones of ``self``.

        OUTPUT:

        - :class:`tuple` of :class:`cones of fan<Cone_of_fan>`.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.generating_cones()
            (2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N,
             2-d cone of Rational polyhedral fan in 2-d lattice N)
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.generating_cones()
            (2-d cone of Rational polyhedral fan in 2-d lattice N,
             1-d cone of Rational polyhedral fan in 2-d lattice N)
        """
        return self._generating_cones

    def is_complete(self):
        r"""
        Check if ``self`` is complete.

        A rational polyhedral fan is *complete* if its cones fill the whole
        space.

        OUTPUT:

        - ``True`` if ``self`` is complete and ``False`` otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.is_complete()
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_complete()
            False
        """
        if "_is_complete" in self.__dict__:
            return self._is_complete
        d = self.lattice_dim()
        if self.dim() != d:
            self._is_complete = False
            return False
        for cone in self:
            if cone.dim() != d:
                self._is_complete = False
                return False
        # Now we know that all generating cones are full-dimensional.
        # Then boundary cones are (d-1)-dimensional.
        for cone in self(codim=1):
            if len(cone.star_generator_indices()) != 2:
                self._is_complete = False
                return False
        self._is_complete = True
        return True

    def is_equivalent(self, other):
        r"""
        Check if ``self`` is "mathematically" the same as ``other``.

        INPUT:

        - ``other`` - fan.

        OUTPUT:

        - ``True`` if ``self`` and ``other`` define the same fans as
          collections of equivalent cones in the same lattice, ``False``
          otherwise.

        There are three different equivalences between fans `F_1` and `F_2`
        in the same lattice:

        #. They have the same rays in the same order and the same generating
           cones in the same order.
           This is tested by ``F1 == F2``.
        #. They have the same rays and the same generating cones without
           taking into account any order.
           This is tested by ``F1.is_equivalent(F2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic toric varieties).
           This is tested by ``F1.is_isomorphic(F2)``.

        EXAMPLES::

            sage: fan1 = Fan(cones=[(0,1), (1,2)],
            ...              rays=[(1,0), (0,1), (-1,-1)],
            ...              check=False)
            sage: fan2 = Fan(cones=[(2,1), (0,2)],
            ...              rays=[(1,0), (-1,-1), (0,1)],
            ...              check=False)
            sage: fan3 = Fan(cones=[(0,1), (1,2)],
            ...              rays=[(1,0), (0,1), (-1,1)],
            ...              check=False)
            sage: fan1 == fan2
            False
            sage: fan1.is_equivalent(fan2)
            True
            sage: fan1 == fan3
            False
            sage: fan1.is_equivalent(fan3)
            False
        """
        if (self.lattice() != other.lattice()
            or self.dim() != other.dim()
            or self.ngenerating_cones() != other.ngenerating_cones()
            or self.ray_set() != other.ray_set()):
            return False
        # Now we need to really compare cones, which can take a while
        return sorted(sorted(cone.rays()) for cone in self) \
               == sorted(sorted(cone.rays()) for cone in other)

    def is_isomorphic(self, other):
        r"""
        Check if ``self`` is in the same `GL(n, \ZZ)`-orbit as ``other``.

        INPUT:

        - ``other`` - fan.

        OUTPUT:

        - ``True`` if ``self`` and ``other`` are in the same
          `GL(n, \ZZ)`-orbit, ``False`` otherwise.

        There are three different equivalences between fans `F_1` and `F_2`
        in the same lattice:

        #. They have the same rays in the same order and the same generating
           cones in the same order.
           This is tested by ``F1 == F2``.
        #. They have the same rays and the same generating cones without
           taking into account any order.
           This is tested by ``F1.is_equivalent(F2)``.
        #. They are in the same orbit of `GL(n,\ZZ)` (and, therefore,
           correspond to isomorphic toric varieties).
           This is tested by ``F1.is_isomorphic(F2)``.

        EXAMPLES:

        These fans are "mirrors" of each other::

            sage: fan1 = Fan(cones=[(0,1), (1,2)],
            ...              rays=[(1,0), (0,1), (-1,-1)],
            ...              check=False)
            sage: fan2 = Fan(cones=[(0,1), (1,2)],
            ...              rays=[(1,0), (0,-1), (-1,1)],
            ...              check=False)
            sage: fan1 == fan2
            False
            sage: fan1.is_equivalent(fan2)
            False
            sage: fan1.is_isomorphic(fan2)
            Traceback (most recent call last):
            ...
            NotImplementedError: fan isomorphism is not implemented yet!
        """
        if self.lattice() != other.lattice():
            return False
        raise NotImplementedError("fan isomorphism is not implemented yet!")

    def is_simplicial(self):
        r"""
        Check if ``self`` is simplicial.

        A rational polyhedral fan is **simplicial** if all of its cones are,
        i.e. primitive vectors along generating rays of every cone form a part
        of a *rational* basis of the ambient space.

        OUTPUT:

        - ``True`` if ``self`` is simplicial and ``False`` otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.is_simplicial()
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_simplicial()
            True

        In fact, any fan in a two-dimensional ambient space is simplicial.
        This is no longer the case in dimension three::

            sage: fan = NormalFan(lattice_polytope.octahedron(3))
            sage: fan.is_simplicial()
            False
            sage: fan.generating_cone(0).nrays()
            4
        """
        if "is_simplicial" not in self.__dict__:
            self._is_simplicial = all(cone.is_simplicial() for cone in self)
        return self._is_simplicial

    def is_smooth(self):
        r"""
        Check if ``self`` is smooth.

        A rational polyhedral fan is **smooth** if all of its cones are,
        i.e. primitive vectors along generating rays of every cone form a part
        of an *integral* basis of the ambient space.
        (In this case the corresponding toric variety is smooth.)

        OUTPUT:

        - ``True`` if ``self`` is smooth and ``False`` otherwise.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.is_smooth()
            True
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.is_smooth()
            True
            sage: fan = NormalFan(lattice_polytope.octahedron(2))
            sage: fan.is_smooth()
            False
            sage: fan.generating_cone(0).rays()
            (N(-1, 1), N(-1, -1))
            sage: fan.generating_cone(0).ray_matrix().det()
            2
        """
        if "_is_smooth" not in self.__dict__:
            self._is_smooth = all(cone.is_smooth() for cone in self)
        return self._is_smooth

    def make_simplicial(self, **kwds):
        r"""
        Construct a simplicial fan subdividing ``self``.

        It is a synonym for :meth:`subdivide` with ``make_simplicial=True``
        option.

        INPUT:

        - this functions accepts only keyword arguments. See :meth:`subdivide`
          for documentation.

        OUTPUT:

        - :class:`rational polyhedral fan
          <sage.geometry.fan.RationalPolyhedralFan>`.

        EXAMPLES::

            sage: fan = NormalFan(lattice_polytope.octahedron(3))
            sage: fan.is_simplicial()
            False
            sage: fan.ngenerating_cones()
            6
            sage: new_fan = fan.make_simplicial()
            sage: new_fan.is_simplicial()
            True
            sage: new_fan.ngenerating_cones()
            12
        """
        return self.subdivide(make_simplicial=True, **kwds)

    def ngenerating_cones(self):
        r"""
        Return the number of generating cones of ``self``.

        OUTPUT:

        - integer.

        EXAMPLES::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: fan.ngenerating_cones()
            4
            sage: cone1 = Cone([(1,0), (0,1)])
            sage: cone2 = Cone([(-1,0)])
            sage: fan = Fan([cone1, cone2])
            sage: fan.ngenerating_cones()
            2
        """
        return len(self.generating_cones())

    def subdivide(self, new_rays=None, make_simplicial=False,
                  algorithm="default", verbose=False):
        r"""
        Construct a new fan subdividing ``self``.

        INPUT:

        - ``new_rays`` - list of new rays to be added during subdivision, each
          ray must be a list or a vector. May be empty or ``None`` (default);

        - ``make_simplicial`` - if ``True``, the returned fan is guaranteed to
          be simplicial, default is ``False``;

        - ``algorithm`` - string with the name of the algorithm used for
          subdivision. Currently there is only one available algorithm called
          "default";

        - ``verbose`` - if ``True``, some timing information may be printed
          during the process of subdivision.

        OUTPUT:

        - :class:`rational polyhedral fan
          <sage.geometry.fan.RationalPolyhedralFan>`.

        Currently the "default" algorithm corresponds to iterative stellar
        subdivison for each ray in ``new_rays``.

        EXAMPLES::

            sage: fan = NormalFan(lattice_polytope.octahedron(3))
            sage: fan.is_simplicial()
            False
            sage: fan.ngenerating_cones()
            6
            sage: fan.nrays()
            8
            sage: new_fan = fan.subdivide(new_rays=[(1,0,0)])
            sage: new_fan.is_simplicial()
            False
            sage: new_fan.ngenerating_cones()
            9
            sage: new_fan.nrays()
            9
        """
        # Maybe these decisions should be done inside the algorithms
        # We can figure it out once we have at least two of them.
        if make_simplicial and not self.is_simplicial():
            rays = list(self.rays())
        else:
            rays = []
        rays.extend(ray for ray in normalize_rays(new_rays, self.lattice())
                        if ray not in self.ray_set())
        if not rays:
            return self # Nothing has to be done
        if algorithm == "default":
            algorithm = "palp"
        method_name = "_subdivide_" + algorithm
        if not hasattr(self, method_name):
            raise ValueError('"%s" is an unknown subdivision algorithm!'
                             % algorithm)
        return getattr(self, method_name)(rays, verbose)


class EnhancedCone(Cone_of_fan):
    r"""
    Construct an enhanced cone.

    Enhanced cones are similar to "regular" :class:`convex rational polyhedral
    cones <sage.geometry.cone.ConvexRationalPolyhedralCone>` or, more
    precisely, to :class:`cones of fans <Cone_of_fan>`, but may have some
    extra functionality. :class:`EnhancedCone` class by itself does not
    provide any of such extra functionality, but in conjunction with
    :class:`EnhancedFan` allows convenient derivation of new classes and
    ensures CPU- and memory-efficient interaction between regular cones and
    their enhanced copies.

    INPUT:

    - ``cone`` -- :class:`Cone_of_fan`;

    - ``ambient`` -- ambient fan of this cone. If not given, will be the same
      as for ``cone``;

    OR:

    - ``ambient`` -- ambient fan of this cone;

    - ``ambient_ray_indices`` -- increasing list or tuple of integers, indices
      of rays of ``ambient`` generating this cone.

    OUTPUT:

    - enhanced cone of exactly the same structure as ``cone``.

    EXAMPLES::

        sage: from sage.geometry.fan import EnhancedCone
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: econe = EnhancedCone(fan.generating_cone(0))
        sage: econe == fan.generating_cone(0)
        True
        sage: econe is fan.generating_cone(0)
        False
    """

    def __init__(self, cone=None, ambient=None, ambient_ray_indices=None):
        r"""
        See :class:`EnhancedCone` for documentation.

        TESTS::

            sage: from sage.geometry.fan import EnhancedCone
            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: econe = EnhancedCone(fan.generating_cone(0))
            sage: econe == fan.generating_cone(0)
            True
            sage: econe is fan.generating_cone(0)
            False
            sage: TestSuite(econe).run()
        """
        if cone is None:
            # Then both other attributes must be given and we can make cone
            fan = ambient._base_fan if is_EnhancedFan(ambient) else ambient
            cone = Cone_of_fan(fan, ambient_ray_indices)
        for attribute in [# Cone attributes
                          "_rays",
                          "_lattice",
                          # _ambient will be set later
                          "_ambient_ray_indices",
                          # Cached attributes copied for efficiency
                          "_is_strictly_convex"]:
            setattr(self, attribute, getattr(cone, attribute))
        self._ambient = ambient if ambient is not None else cone.ambient()


def is_EnhancedCone(x):
    r"""
    Check if ``x`` is an enhanced cone.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is an enhanced cone and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.fan import (
        ...     EnhancedCone, is_EnhancedCone)
        sage: is_EnhancedCone(1)
        False
        sage: cone = Cone([(1,0), (0,1)])
        sage: is_EnhancedCone(cone)
        False
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: econe = EnhancedCone(fan.generating_cone(0))
        sage: is_EnhancedCone(econe)
        True
    """
    return isinstance(x, EnhancedCone)


class EnhancedFan(RationalPolyhedralFan):
    r"""
    Construct an enhanced fan.

    Enhanced fans are similar to "regular" :class:`rational polyhedral fans
    <RationalPolyhedralFan>`, but may have some extra functionality on the
    levels of fans themselves as well as their cones. :class:`EnhancedFan`
    class by itself does not provide any extra functionality, but allows
    convenient derivation of new classes and ensures CPU- and memory-efficient
    interaction between regular fans and their enhanced copies.

    INPUT:

    - ``fan`` -- :class:`rational polyhedral fan <RationalPolyhedralFan>`;

    - ``cone_type`` -- class derived from :class:`EnhancedCone`.

    OUTPUT:

    - enhanced fan of exactly the same structure as ``fan``, but all of its
      associated cones are objects of ``cone_type``.

    EXAMPLES::

        sage: from sage.geometry.fan import (
        ...     EnhancedCone, EnhancedFan, is_EnhancedCone)
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: efan = EnhancedFan(fan, EnhancedCone)
        sage: is_EnhancedCone(efan.generating_cone(0))
        True
    """

    def __init__(self, fan, cone_type):
        r"""
        See :class:`EnhancedFan` for documentation.

        TESTS::

            sage: from sage.geometry.fan import (
            ...     EnhancedCone, EnhancedFan, is_EnhancedCone)
            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: efan = EnhancedFan(fan, EnhancedCone)
            sage: is_EnhancedCone(efan.generating_cone(0))
            True
            sage: TestSuite(efan).run()
        """
        for attribute in [# Fan attributes
                          "_rays",
                          "_lattice"]:
            setattr(self, attribute, getattr(fan, attribute))
        self._base_fan = fan
        self._cone_type = cone_type
        self._generating_cones = tuple(cone_type(cone, self) for cone in fan)

    def _compute_cone_lattice(self):
        r"""
        Compute the cone lattice of ``self``.

        See :meth:`~RationalPolyhedralFan.cone_lattice` for documentation.

        TESTS::

            sage: from sage.geometry.fan import (
            ...     EnhancedCone, EnhancedFan, is_EnhancedCone)
            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: efan = EnhancedFan(fan, EnhancedCone)
            sage: L = efan.cone_lattice() # indirect doctest
            sage: L
            Finite poset containing 10 elements
            sage: is_EnhancedCone(L[0].element)
            True
        """
        # Actual computation of the lattice is done in the "usual" fan
        L =  self._base_fan.cone_lattice()
        # We just copy its structure using the appropriate class for cones
        elements = []
        for element in L:
            cone = element.element
            if is_Cone(cone):
                elements.append(self._cone_type(cone))
            else:
                # The last element is the fan itself
                elements.append(self)
        self._cone_lattice = FinitePoset(L, elements)


def is_EnhancedFan(x):
    r"""
    Check if ``x`` is an enhanced fan.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is an enhanced fan and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.fan import (
        ...     EnhancedCone, EnhancedFan, is_EnhancedFan)
        sage: is_EnhancedFan(1)
        False
        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: is_EnhancedFan(fan)
        False
        sage: efan = EnhancedFan(fan, EnhancedCone)
        sage: is_EnhancedFan(efan)
        True
    """
    return isinstance(x, EnhancedFan)
