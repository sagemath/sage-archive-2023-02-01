r"""
Set of homomorphisms between two toric varieties.

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage,
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- Volker Braun (2012-02-18): Initial version

EXAMPLES:

Here is a simple example, the projection of
`\mathbb{P}^1\times\mathbb{P}^1\to \mathbb{P}^1` ::

    sage: P1xP1 = toric_varieties.P1xP1()
    sage: P1 = toric_varieties.P1()
    sage: hom_set = P1xP1.Hom(P1);  hom_set
    Set of morphisms
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches

In terms of the fan, we can define this morphism by the projection
onto the first coordinate. The Hom-set can construct the morphism from
the projection matrix alone::

    sage: hom_set(matrix([[1],[0]]))
    Scheme morphism:
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches
      Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
            to Rational polyhedral fan in 1-d lattice N.
    sage: _.as_polynomial_map()
    Scheme morphism:
      From: 2-d CPR-Fano toric variety covered by 4 affine patches
      To:   1-d CPR-Fano toric variety covered by 2 affine patches
      Defn: Defined on coordinates by sending [s : t : x : y] to
            [s : t]

In the case of toric algebraic schemes (defined by polynomials in
toric varieties), this module defines the underlying morphism of the
ambient toric varieties::

    sage: P1xP1.inject_variables()
    Defining s, t, x, y
    sage: S = P1xP1.subscheme([s*x-t*y])
    sage: type(S.Hom(S))
    <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>
"""



#*****************************************************************************
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import ZZ
from sage.rings.morphism import is_RingHomomorphism

from sage.matrix.matrix import is_Matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.geometry.fan_morphism import FanMorphism

from sage.schemes.generic.homset import (SchemeHomset_generic,
                                         SchemeHomset_points)


class SchemeHomset_toric_variety(SchemeHomset_generic):
    """
    Set of homomorphisms between two toric varieties.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1 = toric_varieties.P1()
        sage: hom_set = P1xP1.Hom(P1);  hom_set
        Set of morphisms
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
        sage: type(hom_set)
        <class 'sage.schemes.toric.homset.SchemeHomset_toric_variety_with_category'>

        sage: hom_set(matrix([[1],[0]]))
        Scheme morphism:
          From: 2-d CPR-Fano toric variety covered by 4 affine patches
          To:   1-d CPR-Fano toric variety covered by 2 affine patches
          Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                to Rational polyhedral fan in 1-d lattice N.
    """

    def __init__(self, X, Y, category=None, check=True, base=ZZ):
        """
        The Python constructor.

        INPUT:

        The same as for any homset, see
        :mod:`~sage.categories.homset`.

        EXAMPLES::

            sage: P1xP1 = toric_varieties.P1xP1()
            sage: P1 = toric_varieties.P1()
            sage: hom_set = P1xP1.Hom(P1);  hom_set
            Set of morphisms
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches

        An integral matrix defines a fan morphism, since we think of
        the matrix as a linear map on the toric lattice. This is why
        we need to ``register_conversion`` in the constructor
        below. The result is::

            sage: hom_set(matrix([[1],[0]]))
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   1-d CPR-Fano toric variety covered by 2 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 1-d lattice N.
        """
        SchemeHomset_generic.__init__(self, X, Y, category=category, check=check, base=base)
        self.register_conversion(MatrixSpace(ZZ, X.fan().dim(), Y.fan().dim()))

    def _element_constructor_(self, x, check=True):
        """
        Construct a scheme morphism.

        INPUT:

        - `x` -- anything that defines a morphism of toric
          varieties. A matrix, fan morphism, or a list or tuple of
          homogeneous polynomials that define a morphism.

        - ``check`` -- boolean (default: ``True``) passed onto
          functions called by this to be more careful about input
          argument type checking

        OUTPUT:

        The morphism of toric varieties determined by ``x``.

        EXAMPLES:

        First, construct from fan morphism::

            sage: dP8.<t,x0,x1,x2> = toric_varieties.dP8()
            sage: P2.<y0,y1,y2> = toric_varieties.P2()
            sage: hom_set = dP8.Hom(P2)

            sage: fm = FanMorphism(identity_matrix(2), dP8.fan(), P2.fan())
            sage: hom_set(fm)     # calls hom_set._element_constructor_()
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 2-d lattice N.

        A matrix will automatically be converted to a fan morphism::

            sage: hom_set(identity_matrix(2))
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N
                    to Rational polyhedral fan in 2-d lattice N.

        Alternatively, one can use homogeneous polynomials to define morphisms::

            sage: P2.inject_variables()
            Defining y0, y1, y2
            sage: dP8.inject_variables()
            Defining t, x0, x1, x2
            sage: hom_set([x0,x1,x2])
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [t : x0 : x1 : x2] to
                    [x0 : x1 : x2]

        A morphism of the coordinate ring will also work::

            sage: ring_hom = P2.coordinate_ring().hom([x0,x1,x2], dP8.coordinate_ring())
            sage: ring_hom
            Ring morphism:
              From: Multivariate Polynomial Ring in y0, y1, y2 over Rational Field
              To:   Multivariate Polynomial Ring in t, x0, x1, x2 over Rational Field
              Defn: y0 |--> x0
                    y1 |--> x1
                    y2 |--> x2
            sage: hom_set(ring_hom)
            Scheme morphism:
              From: 2-d CPR-Fano toric variety covered by 4 affine patches
              To:   2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined on coordinates by sending [t : x0 : x1 : x2] to
                    [x0 : x1 : x2]
        """
        from sage.schemes.toric.morphism import SchemeMorphism_polynomial_toric_variety
        if isinstance(x, (list, tuple)):
            return SchemeMorphism_polynomial_toric_variety(self, x, check=check)

        if is_RingHomomorphism(x):
            assert x.domain() is self.codomain().coordinate_ring()
            assert x.codomain() is self.domain().coordinate_ring()
            return SchemeMorphism_polynomial_toric_variety(self, x.im_gens(), check=check)

        from sage.schemes.toric.morphism import SchemeMorphism_fan_toric_variety
        if isinstance(x, FanMorphism):
            return SchemeMorphism_fan_toric_variety(self, x, check=check)

        if is_Matrix(x):
            fm = FanMorphism(x, self.domain().fan(), self.codomain().fan())
            return SchemeMorphism_fan_toric_variety(self, fm, check=check)

        raise TypeError, "x must be a fan morphism or a list/tuple of polynomials"


    def _an_element_(self):
        """
        Construct a sample morphism.

        OUTPUT:

        An element of the homset.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: homset = P2.Hom(P2)
            sage: homset.an_element()   # indirect doctest
            Scheme endomorphism of 2-d CPR-Fano toric variety covered by 3 affine patches
              Defn: Defined by sending Rational polyhedral fan in 2-d lattice N to
                    Rational polyhedral fan in 2-d lattice N.
        """
        from sage.matrix.constructor import zero_matrix
        zero = zero_matrix(self.domain().dimension_relative(),
                           self.codomain().dimension_relative())
        return self(zero)

class SchemeHomset_points_toric_field(SchemeHomset_points):
    """
    Set of rational points of a toric variety.

    INPUT:

    - same as for :class:`SchemeHomset_points`.

    OUPUT:

    A scheme morphism of type
    :class:`SchemeHomset_points_toric_field`.

    EXAMPLES::

        sage: P1xP1 = toric_varieties.P1xP1()
        sage: P1xP1(QQ)
        Set of rational points of 2-d CPR-Fano toric variety
        covered by 4 affine patches

    TESTS::

        sage: import sage.schemes.toric.homset as HOM
        sage: HOM.SchemeHomset_points_toric_field(Spec(QQ), P1xP1)
        Set of rational points of 2-d CPR-Fano toric variety covered by 4 affine patches
    """
    pass
