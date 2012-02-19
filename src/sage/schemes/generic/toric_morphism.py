r"""
Morphisms of Toric Varieties

There are three "obvious" ways to map toric varieties to toric
varieties:

1. Polynomial maps in local coordinates, the usual morphisms in
   algebraic geometry.

2. Polynomial maps in the (global) homogeneous coordinates.

3. Toric morphisms, that is, algebraic morphisms equivariant with
   respect to the torus action on the toric variety.

Both 2 and 3 are special cases of 1, which is just to say that we
always remain within the realm of algebraic geometry. But apart from
that, none is included in one of the other cases. In the examples
below, we will explore some algebraic maps that can or can not be
written as a toric morphism. Often a toric morphism can be written
with polynomial maps in homogeneous coordinates, but sometimes it
cannot.

The toric morphisms are perhaps the most mysterious at the
beginning. Let us quickly review their definition (See Definition
3.3.3 of [CLS]_). Let `\Sigma_1` be a fan in `N_{1,\RR}` and `\Sigma_2` be a
fan in `N_{2,\RR}`. A morphism `\phi: X_{\Sigma_1} \to X_{\Sigma_2}`
of the associated toric varieties is toric if `\phi` maps the maximal
torus `T_{N_1} \subseteq X_{\Sigma_1}` into `T_{N_2} \subseteq
X_{\Sigma_2}` and `\phi|_{T_N}` is a group homomorphism.

The data defining a toric morphism is precisely what defines a fan
morphism (see :mod:`~sage.geometry.fan_morphism`), extending the more
familiar dictionary between toric varieties and fans. Toric geometry
is a functor from the category of fans and fan morphisms to the
category of toric varieties and toric morphisms.

.. note::

    Do not create the toric morphisms (or any morphism of schemes)
    directly from the the ``SchemeMorphism...`` classes. Instead, use the
    :meth:`~sage.schemes.generic.scheme.hom` method common to all
    algebraic schemes to create new homomorphisms.

EXAMPLES:

First, consider the following embedding of `\mathbb{P}^1` into
`\mathbb{P}^2` ::

    sage: P2.<x,y,z> = toric_varieties.P2()
    sage: P1.<u,v> = toric_varieties.P1()
    sage: P1.hom([0,u^2+v^2,u*v], P2)
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to
            [0 : u^2 + v^2 : u*v]

This is a well-defined morphism of algebraic varieties because
homogeneously rescaled coordinates of a point of `\mathbb{P}^1` map to the same
point in `\mathbb{P}^2` up to its homogeneous rescalings. It is not
equivariant with respect to the torus actions

.. math::

    \CC^\times \times \mathbb{P}^1,
    (\mu,[u:v]) \mapsto [u:\mu v]
    \quad\text{and}\quad
    \left(\CC^\times\right)^2 \times \mathbb{P}^2,
    ((\alpha,\beta),[x:y:z]) \mapsto [x:\alpha y:\beta z]
    ,

hence it is not a toric morphism. Clearly, the problem is that
the map in homogeneous coordinates contains summands that transform
differently under the torus action. However, this is not the only
difficulty. For example, consider ::

    sage: phi = P1.hom([0,u,v], P2);  phi
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to
            [0 : u : v]

This map is actually the embedding of the
:meth:`~sage.schemes.generic.toric_variety.ToricVariety_field.orbit_closure`
associated to one of the rays of the fan of `\mathbb{P}^2`. Now the
morphism is equivariant with respect to **some** map `\CC^\times \to
(\CC^\times)^2` of the maximal tori of `\mathbb{P}^1` and
`\mathbb{P}^2`. But this map of the maximal tori cannot be the same as
``phi`` defined above. Indeed, the image of ``phi`` completely misses
the maximal torus `T_{\mathbb{P}^2} = \{ [x:y:z] | x\not=0, y\not=0,
z\not=0 \}` of `\mathbb{P}^2`.

Consider instead the following morphism of fans::

    sage: fm = FanMorphism( matrix(ZZ,[[1,0]]), P1.fan(), P2.fan() );  fm
    Fan morphism defined by the matrix
    [1 0]
    Domain fan: Rational polyhedral fan in 1-d lattice N
    Codomain fan: Rational polyhedral fan in 2-d lattice N

which also defines a morphism of toric varieties::

    sage: P1.hom(fm, P2)
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined by sending the Rational polyhedral fan in 1-d lattice N
            to Rational polyhedral fan in 2-d lattice N.

The fan morphism map is equivalent to the following polynomial map::

    sage: _.as_polynomial_map()
    Scheme morphism:
      From: 1-d CPR-Fano toric variety covered by 2 affine patches
      To:   2-d CPR-Fano toric variety covered by 3 affine patches
      Defn: Defined on coordinates by sending [u : v] to
            [u : v : v]

Finally, here is an example of a fan morphism that cannot be written
using homogeneous polynomials. Consider the blowup `O_{\mathbb{P}^1}(2)
\to \CC^2/\ZZ_2`. In terms of toric data, this blowup is::

    sage: A2_Z2 = toric_varieties.A2_Z2()
    sage: A2_Z2.fan().rays()
    (N(1, 0), N(1, 2))
    sage: O2_P1 = A2_Z2.resolve(new_rays=[(1,1)])
    sage: blowup = O2_P1.hom(identity_matrix(2), A2_Z2)
    sage: blowup.as_polynomial_map()
    Traceback (most recent call last):
    ...
    TypeError: The fan morphism cannot be written in homogeneous polynomials.

If we denote the homogeneous coordinates of `O_{\mathbb{P}^1}(2)` by
`x`, `t`, `y` corresponding to the rays `(1,2)`, `(1,1)`, and `(1,0)`
then the blow-up map is [BB]_:

.. math::

    f: O_{\mathbb{P}^1}(2) \to \CC^2/\ZZ_2, \quad
    (x,t,y) \mapsto \left( x\sqrt{t}, y\sqrt{t} \right)

which requires square roots.

REFERENCES:

.. [BB]
    Gavin Brown, Jaroslaw Buczynski:
    Maps of toric varieties in Cox coordinates,
    http://arxiv.org/abs/1004.4924
"""


#*****************************************************************************
#  Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#  Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sequence  import Sequence
from sage.rings.all import ZZ

from sage.schemes.generic.scheme import is_Scheme
from sage.schemes.generic.morphism import (
    is_SchemeMorphism,
    SchemeMorphism, SchemeMorphism_point, SchemeMorphism_polynomial
)



############################################################################
# A points on a toric variety determined by homogeneous coordinates.
class SchemeMorphism_point_toric_field(SchemeMorphism_point):
    """
    A point of a toric variety determined by homogeneous coordinates
    in a field.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.generic.toric_variety.ToricVariety_field>`
        instead.

    INPUT:

    - ``X`` -- toric variety or subscheme of a toric variety.

    - ``coordinates`` -- list of coordinates in the base field of ``X``.

    - ``check`` -- if ``True`` (default), the input will be checked for
      correctness.

    OUTPUT:

    A :class:`SchemeMorphism_point_toric_field`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1(1,2,3,4)
        [1 : 2 : 3 : 4]
    """
    # Mimicking affine/projective classes
    def __init__(self, X, coordinates, check=True):
        r"""
        See :class:`SchemeMorphism_point_toric_field` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1(1,2,3,4)
            [1 : 2 : 3 : 4]
        """
        # Convert scheme to its set of points over the base ring
        if is_Scheme(X):
            X = X(X.base_ring())
        super(SchemeMorphism_point_toric_field, self).__init__(X)
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



############################################################################
# A morphism of toric varieties determined by homogeneous polynomials.
class SchemeMorphism_polynomial_toric_variety(SchemeMorphism_polynomial):
    """
    A morphism determined by homogeneous polynomials.

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.generic.toric_variety.ToricVariety_field>`
        instead.

    INPUT:

    Same as for
    :class:`~sage.schemes.generic.morphism.SchemeMorphism_polynomial`.

    OUPUT:

    A :class:`~sage.schemes.generic.morphism.SchemeMorphism_polynomial_toric_variety`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1.inject_variables()
        Defining z0, z1, z2, z3
        sage: P1 = P1xP1.subscheme(z0-z2)
        sage: H = P1xP1.Hom(P1)
        sage: import sage.schemes.generic.toric_morphism as MOR
        sage: MOR.SchemeMorphism_polynomial_toric_variety(H, [z0,z1,z0,z3])
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
        See :class:`SchemeMorphism_polynomial_toric_variety` for documentation.

        TESTS::

            sage: fan = FaceFan(lattice_polytope.octahedron(2))
            sage: P1xP1 = ToricVariety(fan)
            sage: P1xP1.inject_variables()
            Defining z0, z1, z2, z3
            sage: P1 = P1xP1.subscheme(z0-z2)
            sage: H = P1xP1.Hom(P1)
            sage: import sage.schemes.generic.toric_morphism as MOR
            sage: MOR.SchemeMorphism_polynomial_toric_variety(H, [z0,z1,z0,z3])
            Scheme morphism:
              From: 2-d toric variety covered by 4 affine patches
              To:   Closed subscheme of 2-d toric variety
                    covered by 4 affine patches defined by:
              z0 - z2
              Defn: Defined on coordinates by sending
                    [z0 : z1 : z2 : z3] to [z0 : z1 : z0 : z3]
        """
        SchemeMorphism_polynomial.__init__(self, parent, polynomials, check)
        if check:
            # Check that defining polynomials are homogeneous (degrees can be
            # different if the target uses weighted coordinates)
            for p in self.defining_polynomials():
                if not self.domain().ambient_space().is_homogeneous(p):
                    raise ValueError("%s is not homogeneous!" % p)

    def as_fan_morphism(self):
        """
        Express the morphism as a map defined by a fan morphism.

        OUTPUT:

        A :class:`SchemeMorphism_polynomial_toric_variety`. Raises a
        ``TypeError`` if the morphism cannot be written in such a way.

        EXAMPLES::

            sage: A1.<z> = toric_varieties.A1()
            sage: P1 = toric_varieties.P1()
            sage: patch = A1.hom([1,z], P1)
            sage: patch.as_fan_morphism()
            Traceback (most recent call last):
            ...
            NotImplementedError: expressing toric morphisms as fan morphisms is
            not implemented yet!
        """
        raise NotImplementedError("expressing toric morphisms as fan "
                                  "morphisms is not implemented yet!")


############################################################################
# A morphism of toric varieties determined by a fan morphism
class SchemeMorphism_fan_toric_variety(SchemeMorphism):
    """
    Construct a morphism determined by a fan morphism

    .. WARNING::

        You should not create objects of this class directly. Use the
        :meth:`~sage.schemes.generic.scheme.hom` method of
        :class:`toric varieties
        <sage.schemes.generic.toric_variety.ToricVariety_field>`
        instead.

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are toric varieties.

    - ``fan_morphism`` -- A morphism of fans whose domain and codomain
      fans equal the fans of the domain and codomain in the ``parent``
      Hom-set.

    - ``check`` -- boolean (optional, default:``True``). Whether to
      check the input for consistency.

    OUPUT:

    A :class:`~sage.schemes.generic.morphism.SchemeMorphism_fan_toric_variety`.

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: dP8 = toric_varieties.dP8()
        sage: f = dP8.hom(identity_matrix(2), P2);  f
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   2-d CPR-Fano toric variety covered by 3 affine patches
          Defn: Defined by sending the Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 2-d lattice N.
        sage: type(f)
        <class 'sage.schemes.generic.toric_morphism.SchemeMorphism_fan_toric_variety'>

    Slightly more explicit construction::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1 = toric_varieties.P1()
        sage: hom_set = P1xP1.Hom(P1)
        sage: fm = FanMorphism( matrix(ZZ,[[1],[0]]), P1xP1.fan(), P1.fan() )
        sage: hom_set(fm)
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending the Rational polyhedral fan in 2-d lattice N
               to Rational polyhedral fan in 1-d lattice N.

        sage: P1xP1.hom(fm, P1)
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending the Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 1-d lattice N.
    """

    def __init__(self, parent, fan_morphism, check=True):
        r"""
        See :class:`SchemeMorphism_polynomial_toric_variety` for documentation.

        TESTS::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: hom_set = P1xP1.Hom(P1)
            sage: fan_morphism = FanMorphism( matrix(ZZ,[[1],[0]]), P1xP1.fan(), P1.fan() )
            sage: from sage.schemes.generic.toric_morphism import SchemeMorphism_fan_toric_variety
            sage: SchemeMorphism_fan_toric_variety(hom_set, fan_morphism)
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined by sending the Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 1-d lattice N.
        """
        SchemeMorphism.__init__(self, parent)
        if check and self.domain().fan()!=fan_morphism.domain_fan():
            raise ValueError('The fan morphism domain must be the fan of the domain.')
        if check and self.codomain().fan()!=fan_morphism.codomain_fan():
            raise ValueError('The fan morphism codomain must be the fan of the codomain.')
        self._fan_morphism = fan_morphism

    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: f = P1xP1.hom(matrix([[1],[0]]), P1)
            sage: f._repr_defn()
            'Defined by sending the Rational polyhedral fan in 2-d lattice N to Rational polyhedral fan in 1-d lattice N.'
        """
        s  = 'Defined by sending the '
        s += str(self.domain().fan())
        s += ' to '
        s += str(self.codomain().fan())
        s += '.'
        return s

    def fan_morphism(self):
        """
        Return the defining fan morphism.

        OUTPUT:

        A :class:`~sage.geometry.fan_morphism.FanMorphism`.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: f = P1xP1.hom(matrix([[1],[0]]), P1)
            sage: f.fan_morphism()
            Fan morphism defined by the matrix
            [1]
            [0]
            Domain fan: Rational polyhedral fan in 2-d lattice N
            Codomain fan: Rational polyhedral fan in 1-d lattice N
        """
        return self._fan_morphism

    def as_polynomial_map(self):
        """
        Express the morphism via homogeneous polynomials.

        OUTPUT:

        A :class:`SchemeMorphism_polynomial_toric_variety`. Raises a
        ``TypeError`` if the morphism cannot be written in terms of
        homogeneous polynomials.

        EXAMPLES::

            sage: A1 = toric_varieties.A1()
            sage: square = A1.hom(matrix([[2]]), A1)
            sage: square.as_polynomial_map()
            Scheme endomorphism of 1-d affine toric variety
              Defn: Defined on coordinates by sending [z] to
                    [z^2]

            sage: P1 = toric_varieties.P1()
            sage: patch = A1.hom(matrix([[1]]), P1)
            sage: patch.as_polynomial_map()
            Scheme morphism:
              From: 1-d affine toric variety
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined on coordinates by sending [z] to
                    [z : 1]
        """
        R = self.domain().coordinate_ring()
        phi = self.fan_morphism()
        polys = [R.one()] * phi.codomain_fan().nrays()
        for rho, x in zip(phi.domain_fan(1), R.gens()):
            ray = rho.ray(0)
            sigma = phi.image_cone(rho)
            degrees = sigma.ray_matrix().solve_right(phi(ray))
            for i, d in zip(sigma.ambient_ray_indices(), degrees):
                try:
                    d = ZZ(d)
                except TypeError:
                    raise TypeError('The fan morphism cannot be written in homogeneous polynomials.')
                polys[i] *= x**d
        return SchemeMorphism_polynomial_toric_variety(self.parent(), polys)

