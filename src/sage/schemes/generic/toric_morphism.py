"""
Morphisms of Toric Varieties


"""


#*****************************************************************************
#  Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#  Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sequence  import Sequence

import scheme
from sage.schemes.generic.morphism import (
    is_SchemeMorphism,
    SchemeMorphism_coordinates, SchemeMorphism_on_points
)



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
        if scheme.is_Scheme(X):
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

    - same as for
      :class:`~sage.schemes.generic.morphism.SchemeMorphism_on_points`.

    OUPUT:

    - :class:`~sage.schemes.generic.morphism.SchemeMorphism_on_points_toric_variety`.

    TESTS::

        sage: fan = FaceFan(lattice_polytope.octahedron(2))
        sage: P1xP1 = ToricVariety(fan)
        sage: P1xP1.inject_variables()
        Defining z0, z1, z2, z3
        sage: P1 = P1xP1.subscheme(z0-z2)
        sage: H = P1xP1.Hom(P1)
        sage: import sage.schemes.generic.toric_morphism as MOR
        sage: MOR.SchemeMorphism_on_points_toric_variety(H, [z0,z1,z0,z3])
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
            sage: import sage.schemes.generic.toric_morphism as MOR
            sage: MOR.SchemeMorphism_on_points_toric_variety(H, [z0,z1,z0,z3])
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



