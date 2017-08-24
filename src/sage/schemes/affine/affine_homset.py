r"""
Set of homomorphisms between two affine schemes

For schemes `X` and `Y`, this module implements the set of morphisms
`Hom(X,Y)`. This is done by :class:`SchemeHomset_generic`.

As a special case, the Hom-sets can also represent the points of a
scheme. Recall that the `K`-rational points of a scheme `X` over `k`
can be identified with the set of morphisms `Spec(K) \to X`. In Sage
the rational points are implemented by such scheme morphisms. This is
done by :class:`SchemeHomset_points` and its subclasses.

.. note::

    You should not create the Hom-sets manually. Instead, use the
    :meth:`~sage.structure.parent.Hom` method that is inherited by all
    schemes.

AUTHORS:

- William Stein (2006): initial version.
"""


#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.rings.all import ZZ
from sage.rings.rational_field import is_RationalField
from sage.categories.fields import Fields
from sage.categories.number_fields import NumberFields
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
import sage.schemes.generic.homset
from copy import copy

#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeHomset_points_spec(sage.schemes.generic.homset.SchemeHomset_generic):
    """
    Set of rational points of an affine variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
        sage: SchemeHomset_points_spec(Spec(QQ), Spec(QQ))
        Set of rational points of Spectrum of Rational Field
    """

    def _element_constructor_(self, *args, **kwds):
        """
        The element constructor.

        EXAMPLES::

            sage: X = Spec(QQ)
            sage: ring_hom = QQ.hom((1,), QQ);  ring_hom
            Ring endomorphism of Rational Field
              Defn: 1 |--> 1
            sage: H = X.Hom(X)
            sage: H(ring_hom)
            Affine Scheme endomorphism of Spectrum of Rational Field
              Defn: Ring endomorphism of Rational Field
                      Defn: 1 |--> 1

        TESTS::

            sage: H._element_constructor_(ring_hom)
            Affine Scheme endomorphism of Spectrum of Rational Field
              Defn: Ring endomorphism of Rational Field
                      Defn: 1 |--> 1
        """
        return sage.schemes.generic.homset.SchemeHomset_generic._element_constructor_(self, *args, **kwds)

    def _repr_(self):
        """
        Return a string representation of a homset.

        OUTPUT: A string.

        EXAMPLES::

            sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_spec
            sage: S = SchemeHomset_points_spec(Spec(QQ), Spec(QQ))
            sage: S._repr_()
            'Set of rational points of Spectrum of Rational Field'
        """
        return 'Set of rational points of '+str(self.codomain())



#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeHomset_points_affine(sage.schemes.generic.homset.SchemeHomset_points):
    """
    Set of rational points of an affine variety.

    INPUT:

    See :class:`SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.affine.affine_homset import SchemeHomset_points_affine
        sage: SchemeHomset_points_affine(Spec(QQ), AffineSpace(ZZ,2))
        Set of rational points of Affine Space of dimension 2 over Rational Field
    """

    def points(self, B=0):
        r"""
        Return some or all rational points of an affine scheme.

        INPUT:

        - ``B`` -- integer (optional, default: 0). The bound for the
          height of the coordinates.

        OUTPUT:

        - If the base ring is a finite field: all points of the scheme,
          given by coordinate tuples.

        - If the base ring is `\QQ` or `\ZZ`: the subset of points whose
          coordinates have height ``B`` or less.

        EXAMPLES: The bug reported at #11526 is fixed::

            sage: A2 = AffineSpace(ZZ, 2)
            sage: F = GF(3)
            sage: A2(F).points()
            [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)]

            sage: R = ZZ
            sage: A.<x,y> = R[]
            sage: I = A.ideal(x^2-y^2-1)
            sage: V = AffineSpace(R, 2)
            sage: X = V.subscheme(I)
            sage: M = X(R)
            sage: M.points(1)
            [(-1, 0), (1, 0)]

        ::

            sage: u = QQ['u'].0
            sage: K.<v> = NumberField(u^2 + 3)
            sage: A.<x,y> = AffineSpace(K, 2)
            sage: len(A(K).points(9))
            361

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: E = A.subscheme([x^2 + y^2 - 1, y^2 - x^3 + x^2 + x - 1])
            sage: E(A.base_ring()).points()
            [(0, 1), (0, -1), (1, 0), (-1, 0)]

        ::

            sage: A.<x,y> = AffineSpace(CC, 2)
            sage: E = A.subscheme([y^3-x^3-x^2, x*y])
            sage: E(A.base_ring()).points()
            [(-1.00000000000000, 0.000000000000000),
            (0.000000000000000, 0.000000000000000)]
        """
        X = self.codomain()
        from sage.schemes.affine.affine_space import is_AffineSpace
        from sage.rings.all import CC
        if not is_AffineSpace(X) and X.base_ring() in Fields():
            if X.base_ring() == CC:
                complex = True
            else:
                complex = False
            # Then X must be a subscheme
            dim_ideal = X.defining_ideal().dimension()
            if dim_ideal < 0: # no points
                return []
            if dim_ideal == 0: # if X zero-dimensional
                rat_points = []
                PS = X.ambient_space()
                N = PS.dimension_relative()
                BR = X.base_ring()
                #need a lexicographic ordering for elimination
                R = PolynomialRing(BR, N, PS.gens(), order='lex')
                I = R.ideal(X.defining_polynomials())
                I0 = R.ideal(0)
                #Determine the points through elimination
                #This is much faster than using the I.variety() function on each affine chart.
                G = I.groebner_basis()
                if G != [1]:
                    P = {}
                    points = [P]
                    #work backwards from solving each equation for the possible
                    #values of the next coordinate
                    for i in range(len(G) - 1, -1, -1):
                        new_points = []
                        good = 0
                        for P in points:
                            #substitute in our dictionary entry that has the values
                            #of coordinates known so far. This results in a single
                            #variable polynomial (by elimination)
                            L = G[i].substitute(P)
                            if L != 0:
                                if complex:
                                    for pol in L.univariate_polynomial().roots(multiplicities=False):
                                        r = L.variables()[0]
                                        varindex = R.gens().index(r)
                                        P.update({R.gen(varindex):pol})
                                        new_points.append(copy(P))
                                        good = 1
                                else:
                                    L = L.factor()
                                #the linear factors give the possible rational values of
                                #this coordinate
                                    for pol, pow in L:
                                        if pol.degree() == 1 and len(pol.variables()) == 1:
                                            good = 1
                                            r = pol.variables()[0]
                                            varindex = R.gens().index(r)
                                            #add this coordinates information to
                                            #each dictionary entry
                                            P.update({R.gen(varindex):-pol.constant_coefficient() /
                                            pol.monomial_coefficient(r)})
                                            new_points.append(copy(P))
                            else:
                                new_points.append(P)
                                good = 1
                        if good:
                            points = new_points
                    #the dictionary entries now have values for all coordinates
                    #they are the rational solutions to the equations
                    #make them into projective points
                    for i in range(len(points)):
                        if complex:
                            if len(points[i]) == N:
                                S = X.ambient_space()([points[i][R.gen(j)] for j in range(N)])
                                #S.normalize_coordinates()
                                rat_points.append(S)
                        else:
                            if len(points[i]) == N and I.subs(points[i]) == I0:
                                S = X([points[i][R.gen(j)] for j in range(N)])
                                #S.normalize_coordinates()
                                rat_points.append(S)
                #rat_points = sorted(rat_points)
                return rat_points
        R = self.value_ring()
        if is_RationalField(R) or R == ZZ:
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.affine.affine_rational_point import enum_affine_rational_field
            return enum_affine_rational_field(self,B)
        if R in NumberFields():
            if not B > 0:
                raise TypeError("a positive bound B (= %s) must be specified"%B)
            from sage.schemes.affine.affine_rational_point import enum_affine_number_field
            return enum_affine_number_field(self,B)
        elif is_FiniteField(R):
            from sage.schemes.affine.affine_rational_point import enum_affine_finite_field
            return enum_affine_finite_field(self)
        else:
            raise TypeError("unable to enumerate points over %s"%R)
