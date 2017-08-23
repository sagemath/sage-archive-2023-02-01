r"""
Generic dynamical systems on schemes

This is the generic class for dynamical systems and contains the exported
constructor functions. The constructor functions can take either polynomials
(or rational functions in the affine case) or morphisms from which to construct
a dynamical system. If the domain is not specified, it is constructed. However,
if you plan on working with points or subvarieties in the domain, it recommended
to specify the domain. For products of projective spaces the domain must be specified.

The initialization checks are always performed by the constructor functions. It is
possible, but not recommended, to skip these checks by calling the class initialization
directly.

AUTHORS:

- Ben Hutz (July 2017): initial version
"""

#*****************************************************************************
#       Copyright (C) 2017 Ben Hutz <bn4941@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import, print_function
from sage.categories.homset import End
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.rings.finite_rings.finite_field_constructor import is_FiniteField
from sage.rings.fraction_field import is_FractionField
from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
from sage.rings.polynomial.multi_polynomial_ring import is_MPolynomialRing
from sage.rings.quotient_ring import is_QuotientRing
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_affine
from sage.schemes.generic.algebraic_scheme import AlgebraicScheme_subscheme_projective
from sage.schemes.projective.projective_space import is_ProjectiveSpace
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.product_projective.space import is_ProductProjectiveSpaces
from sage.symbolic.ring import SR

from sage.categories.fields import Fields
_Fields = Fields()

def is_DynamicalSystem(f):
    r"""
    Test whether ``f`` is a dynamical system.

    INPUT: ``f`` -- any object.

    OUTPUT: Boolean. Return ``True`` if ``f`` is a dynamical system.

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(QQ,2)
        sage: f = DynamicalSystem_affine([y,x^2+y]); f
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x^2 + y)
        sage: is_DynamicalSystem(f)
        True

    ::

        sage: A.<x,y> = AffineSpace(QQ,2)
        sage: H = End(A)
        sage: f = H([y,x^2+y]); f
        Scheme endomorphism of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x^2 + y)
        sage: is_DynamicalSystem(f)
        False

    ::

        sage: P2.<x,y> = ProjectiveSpace(QQ,1)
        sage: H = End(P2)
        sage: f = H([y^2, x^2 + y^2])
        sage: is_DynamicalSystem(f)
        False
        sage: is_DynamicalSystem(f.as_dynamical_system())
        True

    ::

        sage: PP.<x,y,u,v> = ProductProjectiveSpaces(QQ,[1,1])
        sage: f = DynamicalSystem_projective([u*y^2, v*x^2, u*x, y*v], PP)
        sage: is_DynamicalSystem(f)
        True
        sage: is_DynamicalSystem(f.as_scheme_morphism())
        False
    """
    return isinstance(f, DynamicalSystem_generic)


def DynamicalSystem_affine(morphism_or_polys, domain=None):
    r"""
    Return a dynamical system on an affine scheme.

    INPUT:

    - ``morphism_or_polys`` -- a SchemeMorphism, a polynomial, a
      rational function, or a list or tuple of polynomials or rational
      functions.

    - ``domain`` -- optional affine space or subscheme of such.

      The following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case.

      * ``morphism_or_polys`` is a list of polynomials or rational
        functions that define a rational endomorphism of ``domain``.

      * ``morphism_or_polys`` is a list of polynomials or rational
        functions and ``domain`` is unspecified; ``domain`` is then
        taken to be the affine space of appropriate dimension over the
        base ring of the first element of ``morphism_or_polys``.

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and assumed to be the
        1-dimensional affine space over the base ring of
        ``morphism_or_polys``.

    OUTPUT: :class:`DynamicalSystem_affine`.

    EXAMPLES::

        sage: R.<t> = ZZ[]
        sage: DynamicalSystem_affine(t^2 - 1)
        Dynamical System of Affine Space of dimension 1 over Integer Ring
          Defn: Defined on coordinates by sending (t) to
                (t^2 - 1)

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: X = A.subscheme([x-y^2])
        sage: DynamicalSystem_affine([9/4*x^2, 3/2*y], domain=X)
        Dynamical System of Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              -y^2 + x
              Defn: Defined on coordinates by sending (x, y) to
                    (9/4*x^2, 3/2*y)

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: H = End(A)
        sage: f = H([x^2, y^2])
        sage: DynamicalSystem_affine(f)
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x^2, y^2)

    Notice that `ZZ` becomes `QQ` since the function is rational::

        sage: A.<x,y> = AffineSpace(ZZ, 2)
        sage: DynamicalSystem_affine([3*x^2/(5*y), y^2/(2*x^2)])
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (3*x^2/(5*y), y^2/(2*x^2))

    ::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: DynamicalSystem_affine([3/2*x^2, y^2])
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (3/2*x^2, y^2)

    If you pass in quotient ring elements, they are reduced::

        sage: A.<x,y,z> = AffineSpace(QQ, 3)
        sage: X = A.subscheme([x-y])
        sage: u,v,w = X.coordinate_ring().gens()
        sage: DynamicalSystem_affine([u, v, u+v], domain=X)
        Dynamical System of Closed subscheme of Affine Space of dimension 3
        over Rational Field defined by:
          x - y
          Defn: Defined on coordinates by sending (x, y, z) to
                (y, y, 2*y)

    ::

        sage: R.<t> = PolynomialRing(QQ)
        sage: A.<x,y,z> = AffineSpace(R, 3)
        sage: X = A.subscheme(x^2-y^2)
        sage: H = End(X)
        sage: f = H([x^2/(t*y), t*y^2, x*z])
        sage: DynamicalSystem_affine(f)
        Dynamical System of Closed subscheme of Affine Space of dimension 3
        over Univariate Polynomial Ring in t over Rational Field defined by:
          x^2 - y^2
          Defn: Defined on coordinates by sending (x, y, z) to
                (x^2/(t*y), t*y^2, x*z)

    ::

        sage: x = var('x')
        sage: DynamicalSystem_affine(x^2+1)
        Traceback (most recent call last):
        ...
        TypeError: Symbolic Ring cannot be the base ring

    TESTS::

        sage: A.<x> = AffineSpace(ZZ,1)
        sage: A1.<z> = AffineSpace(CC,1)
        sage: H = End(A1)
        sage: f2 = H([z^2+1])
        sage: f = DynamicalSystem_affine(f2, A)
        sage: f.domain() is A
        False

    ::

        sage: P1.<x,y> = ProjectiveSpace(QQ,1)
        sage: DynamicalSystem_affine([y, 2*x], domain=P1)
        Traceback (most recent call last):
        ...
        ValueError: "domain" must be an affine scheme
        sage: H = End(P1)
        sage: DynamicalSystem_affine(H([y, 2*x]))
        Traceback (most recent call last):
        ...
        ValueError: "domain" must be an affine scheme
    """
    from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_ring
    from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_field
    from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine_finite_field

    if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
        morphism = morphism_or_polys
        R = morphism.base_ring()
        polys = list(morphism)
        domain = morphism.domain()
        if not is_AffineSpace(domain) and not isinstance(domain, AlgebraicScheme_subscheme_affine):
            raise ValueError('"domain" must be an affine scheme')
        if domain != morphism_or_polys.codomain():
            raise ValueError('domain and codomain do not agree')
        if R not in _Fields:
            return DynamicalSystem_affine_ring(polys, domain)
        if is_FiniteField(R):
            return DynamicalSystem_affine_finite_field(polys, domain)
        return DynamicalSystem_affine_field(polys, domain)
    elif isinstance(morphism_or_polys,(list, tuple)):
        polys = list(morphism_or_polys)
    else:
        polys = [morphism_or_polys]

    # We now arrange for all of our list entries to lie in the same ring
    # Fraction field case first
    fraction_field = False
    for poly in polys:
        P = poly.parent()
        if is_FractionField(P):
            fraction_field = True
            break
    if fraction_field:
        K = P.base_ring().fraction_field()
        # Replace base ring with its fraction field
        P = P.ring().change_ring(K).fraction_field()
        polys = [P(poly) for poly in polys]
    else:
        # If any of the list entries lies in a quotient ring, we try
        # to lift all entries to a common polynomial ring.
        quotient_ring = False
        for poly in polys:
            P = poly.parent()
            if is_QuotientRing(P):
                quotient_ring = True
                break
        if quotient_ring:
            polys = [P(poly).lift() for poly in polys]
        else:
            poly_ring = False
            for poly in polys:
                P = poly.parent()
                if is_PolynomialRing(P) or is_MPolynomialRing(P):
                    poly_ring = True
                    break
            if poly_ring:
                polys = [P(poly) for poly in polys]

    if domain is None:
        f = polys[0]
        CR = f.parent()
        if CR is SR:
            raise TypeError("Symbolic Ring cannot be the base ring")
        if fraction_field:
            CR = CR.ring()
        domain = AffineSpace(CR)

    R = domain.base_ring()
    if R is SR:
        raise TypeError("Symbolic Ring cannot be the base ring")
    if not is_AffineSpace(domain) and not isinstance(domain, AlgebraicScheme_subscheme_affine):
        raise ValueError('"domain" must be an affine scheme')
    if R not in _Fields:
        return DynamicalSystem_affine_ring(polys, domain)
    if is_FiniteField(R):
            return DynamicalSystem_affine_finite_field(polys, domain)
    return DynamicalSystem_affine_field(polys, domain)


def DynamicalSystem_projective(morphism_or_polys, domain=None, names=None):
    r"""
    Return a dynamical system on a projective scheme.

    INPUT:

    - ``morphism_or_polys`` -- a SchemeMorphism, a polynomial, a
      rational function, or a list or tuple of homogeneous polynomials.

    - ``domain`` -- optional projective space or projective subscheme.

    - ``names`` -- optional tuple of strings to be used as coordinate
      names for a projective space that is constructed; defaults to ``'X','Y'``.

      The following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case.

      * ``morphism_or_polys`` is a list of homogeneous polynomials
        that define a rational endomorphism of ``domain``.

      * ``morphism_or_polys`` is a list of homogeneous polynomials and
        ``domain`` is unspecified; ``domain`` is then taken to be the
        projective space of appropriate dimension over the base ring of
        the first element of ``morphism_or_polys``.

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and taken to be a
        1-dimensional projective space over the base ring of
        ``morphism_or_polys`` with coordinate names given by ``names``.

    OUTPUT: :class:`DynamicalSystem_projectve`.

    EXAMPLES::

        sage: P1.<x,y> = ProjectiveSpace(QQ,1)
        sage: DynamicalSystem_projective([y, 2*x])
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (y : 2*x)

    We can define dynamical systems on `P^1` by giving a polynomial or
    rational function::

        sage: R.<t> = QQ[]
        sage: DynamicalSystem_projective(t^2 - 3)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 - 3*Y^2 : Y^2)
        sage: DynamicalSystem_projective(1/t^2)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (Y^2 : X^2)

    ::

        sage: R.<x> = PolynomialRing(QQ,1)
        sage: DynamicalSystem_projective(x^2, names=['a','b'])
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (a : b) to
                (a^2 : b^2)

    Symbolic Ring elements are not allows::

        sage: x,y = var('x,y')
        sage: DynamicalSystem_projective([x^2,y^2])
        Traceback (most recent call last):
        ...
        ValueError: [x^2, y^2] must be elements of a polynomial ring

    ::

        sage: R.<x> = PolynomialRing(QQ,1)
        sage: DynamicalSystem_projective(x^2)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 : Y^2)

    ::

        sage: R.<t> = PolynomialRing(QQ)
        sage: P.<x,y,z> = ProjectiveSpace(R, 2)
        sage: X = P.subscheme([x])
        sage: DynamicalSystem_projective([x^2, t*y^2, x*z], domain=X)
        Dynamical System of Closed subscheme of Projective Space of dimension
        2 over Univariate Polynomial Ring in t over Rational Field defined by:
          x
          Defn: Defined on coordinates by sending (x : y : z) to
                (x^2 : t*y^2 : x*z)

    When elements of the quotient ring are used, they are reduced::

        sage: P.<x,y,z> = ProjectiveSpace(CC, 2)
        sage: X = P.subscheme([x-y])
        sage: u,v,w = X.coordinate_ring().gens()
        sage: DynamicalSystem_projective([u^2, v^2, w*u], domain=X)
        Dynamical System of Closed subscheme of Projective Space of dimension
        2 over Complex Field with 53 bits of precision defined by:
          x - y
          Defn: Defined on coordinates by sending (x : y : z) to
                (y^2 : y^2 : y*z)

    We can also compute the forward image of subschemes through
    elimination. In particular, let `X = V(h_1,\ldots, h_t)` and define the ideal
    `I = (h_1,\ldots,h_t,y_0-f_0(\bar{x}), \ldots, y_n-f_n(\bar{x}))`.
    Then the elimination ideal `I_{n+1} = I \cap K[y_0,\ldots,y_n]` is a homogeneous
    ideal and `f(X) = V(I_{n+1})`::

        sage: P.<x,y,z> = ProjectiveSpace(QQ, 2)
        sage: f = DynamicalSystem_projective([(x-2*y)^2, (x-2*z)^2, x^2])
        sage: X = P.subscheme(y-z)
        sage: f(f(f(X)))
        Closed subscheme of Projective Space of dimension 2 over Rational Field
        defined by:
          y - z

    ::

        sage: P.<x,y,z,w> = ProjectiveSpace(QQ, 3)
        sage: f = DynamicalSystem_projective([(x-2*y)^2, (x-2*z)^2, (x-2*w)^2, x^2])
        sage: f(P.subscheme([x,y,z]))
        Closed subscheme of Projective Space of dimension 3 over Rational Field
        defined by:
          w,
          y,
          x

    ::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u^2], domain=T)
        Dynamical System of Product of projective spaces P^2 x P^1 over Rational Field
          Defn: Defined by sending (x : y : z , w : u) to
                (x^2*u : y^2*w : z^2*u , w^2 : u^2).

    TESTS::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: f = DynamicalSystem_projective([x-y, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - y, x*y]) must be of the same degree
        sage: DynamicalSystem_projective([x-1, x*y+x])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

    ::

        sage: DynamicalSystem_projective([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        ValueError: [e^x, e^y] must be elements of a polynomial ring

    ::

        sage: T.<x,y,z,w,u> = ProductProjectiveSpaces([2, 1], QQ)
        sage: DynamicalSystem_projective([x^2*u, y^2*w, z^2*u, w^2, u*z], domain=T)
        Traceback (most recent call last):
        ...
        TypeError: polys (=[x^2*u, y^2*w, z^2*u, w^2, z*u]) must be
        multi-homogeneous of the same degrees (by component)

    ::

        sage: A.<x,y> = AffineSpace(ZZ, 2)
        sage: DynamicalSystem_projective([x^2,y^2], A)
        Traceback (most recent call last):
        ...
        ValueError: "domain" must be a projective scheme
        sage: H = End(A)
        sage: f = H([x,y])
        sage: DynamicalSystem_projective(f)
        Traceback (most recent call last):
        ...
        ValueError: "domain" must be a projective scheme

    ::

        sage: R.<x> = PolynomialRing(QQ,1)
        sage: DynamicalSystem_projective(x^2, names='t')
        Traceback (most recent call last):
        ...
        ValueError: specify 2 variable names

    ::

        sage: P1.<x,y> = ProjectiveSpace(QQ,1)
        sage: DynamicalSystem_projective([y, x, y], domain=P1)
        Traceback (most recent call last):
        ...
        ValueError: polys (=[y, x, y]) do not define a rational endomorphism of the domain

    ::

        sage: A.<x,y> = AffineSpace(QQ,2)
        sage: DynamicalSystem_projective([y,x], domain=A)
        Traceback (most recent call last):
        ...
        ValueError: "domain" must be a projective scheme

    ::

        sage: R.<x> = QQ[]
        sage: DynamicalSystem([x^2])
        Traceback (most recent call last):
        ...
        ValueError: list/tuple must have at least 2 polynomials
    """
    from sage.dynamics.arithmetic_dynamics.product_projective_ds import DynamicalSystem_product_projective_ring
    from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_ring
    from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_field
    from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective_finite_field

    if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
        R = morphism_or_polys.base_ring()
        domain = morphism_or_polys.domain()
        polys = list(morphism_or_polys)
        if domain != morphism_or_polys.codomain():
            raise ValueError('domain and codomain do not agree')
        if not is_ProjectiveSpace(domain) and not isinstance(domain, AlgebraicScheme_subscheme_projective):
            raise ValueError('"domain" must be a projective scheme')
        if R not in _Fields:
            return DynamicalSystem_projective_ring(polys, domain)
        if is_FiniteField(R):
            return DynamicalSystem_projective_finite_field(polys, domain)
        return DynamicalSystem_projective_field(polys, domain)

    if isinstance(morphism_or_polys,(list, tuple)):
        polys = list(morphism_or_polys)
        if len(polys) == 1:
            raise ValueError("list/tuple must have at least 2 polynomials")
        test = lambda x: is_PolynomialRing(x) or is_MPolynomialRing(x)
        if not all(test(poly.parent()) for poly in polys):
            try:
                polys = [poly.lift() for poly in polys]
            except:
                raise ValueError('{} must be elements of a polynomial ring'.format(morphism_or_polys))
    else:
        # homogenize!
        f = morphism_or_polys
        aff_CR = f.parent()
        if not is_PolynomialRing(aff_CR) and not is_FractionField(aff_CR)\
            and not (is_MPolynomialRing(aff_CR) and aff_CR.ngens() == 1):
            msg = '{} is not a single variable polynomial or rational function'
            raise ValueError(msg.format(f))
        if is_FractionField(aff_CR):
            polys = [f.numerator(),f.denominator()]
        else:
            polys = [f, aff_CR(1)]
        d = max(poly.degree() for poly in polys)
        if names is None:
            names = ('X','Y')
        elif len(names) != 2:
            raise ValueError('specify 2 variable names')
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        proj_CR = PolynomialRing(aff_CR.base_ring(), names=names)
        X,Y = proj_CR.gens()
        polys = [proj_CR(Y**d * poly(X/Y)) for poly in polys]
    if domain is None:
        f = polys[0]
        proj_CR = f.parent()
        domain = ProjectiveSpace(proj_CR)
    R = domain.base_ring()
    if R is SR:
        raise TypeError("Symbolic Ring cannot be the base ring")

    if len(polys) != domain.ambient_space().coordinate_ring().ngens():
        msg = 'polys (={}) do not define a rational endomorphism of the domain'
        raise ValueError(msg.format(polys))

    if is_ProductProjectiveSpaces(domain):
        splitpolys = domain._factors(polys)
        for split_poly in splitpolys:
            split_d = domain._degree(split_poly[0])
            if not all(split_d == domain._degree(f) for f in split_poly):
                msg = 'polys (={}) must be multi-homogeneous of the same degrees (by component)'
                raise TypeError(msg.format(polys))
        return DynamicalSystem_product_projective_ring(polys, domain)

    # Now polys define an endomorphism of a scheme in P^n
    if not all(poly.is_homogeneous() for poly in polys):
        msg = 'polys (={}) must be homogeneous'
        raise ValueError(msg.format(polys))
    d = polys[0].degree()
    if not all(poly.degree() == d for poly in polys):
        msg = 'polys (={}) must be of the same degree'
        raise ValueError(msg.format(polys))

    if not is_ProjectiveSpace(domain) and not isinstance(domain, AlgebraicScheme_subscheme_projective):
        raise ValueError('"domain" must be a projective scheme')
    if R not in _Fields:
        return DynamicalSystem_projective_ring(polys, domain)
    if is_FiniteField(R):
            return DynamicalSystem_projective_finite_field(polys, domain)
    return DynamicalSystem_projective_field(polys, domain)


def DynamicalSystem(morphism_or_polys, domain=None, names=None):
    r"""
    Return a dynamical system on a scheme.

    If domain is not specified either through the domain of the morphism
    or the domain parameter, it is assumed to be projective space.

    INPUT:

    - ``morphism_or_polys`` -- a SchemeMorphism, a polynomial, a
      rational function, or a list or tuple of homogeneous polynomials.

    - ``domain`` -- optional scheme space or subscheme.

    - ``names`` -- optional tuple of strings to be used as coordinate
      names for a projective space that is constructed; defaults to ``'X','Y'``.

      The following combinations of ``morphism_or_polys`` and
      ``domain`` are meaningful:

      * ``morphism_or_polys`` is a SchemeMorphism; ``domain`` is
        ignored in this case.

      * ``morphism_or_polys`` is a list of homogeneous polynomials
        that define a rational endomorphism of ``domain``.

      * ``morphism_or_polys`` is a list of homogeneous polynomials and
        ``domain`` is unspecified; ``domain`` is then taken to be the
        projective space of appropriate dimension over the base ring of
        the first element of ``morphism_or_polys``.

      * ``morphism_or_polys`` is a single polynomial or rational
        function; ``domain`` is ignored and taken to be a
        1-dimensional projective space over the base ring of
        ``morphism_or_polys`` with coordinate names given by ``names``.

    OUTPUT: :class:`DynamicalSystem_projective` or :class:`DynamicalSystem_affine`.

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(QQ, 2)
        sage: DynamicalSystem([y, x], domain=A)
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x)
        sage: H = End(A)
        sage: DynamicalSystem(H([y, x]))
        Dynamical System of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x)

        ::

        sage: P1.<x,y> = ProjectiveSpace(CC,1)
        sage: H = End(P1)
        sage: DynamicalSystem(H([y, x]))
        Dynamical System of Projective Space of dimension 1 over Complex Field
        with 53 bits of precision
          Defn: Defined on coordinates by sending (x : y) to
                (y : x)

        DynamicalSystem defaults to projective::

        sage: R.<x,y,z> = QQ[]
        sage: DynamicalSystem([x^2, y^2, z^2])
        Dynamical System of Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x : y : z) to
                (x^2 : y^2 : z^2)

    ::

        sage: R.<t> = QQ[]
        sage: DynamicalSystem(t^2 - 3)
        Dynamical System of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (X : Y) to
                (X^2 - 3*Y^2 : Y^2)
    """
    if isinstance(morphism_or_polys, SchemeMorphism_polynomial):
        domain = morphism_or_polys.domain()
    if not domain is None:
        if is_AffineSpace(domain) or isinstance(domain, AlgebraicScheme_subscheme_affine):
            return DynamicalSystem_affine(morphism_or_polys, domain)
    return DynamicalSystem_projective(morphism_or_polys, domain, names)

class DynamicalSystem_generic(SchemeMorphism_polynomial):
    r"""
    Base class for dynamical systems of schemes.

    INPUT:

    - ``polys_or_rat_fncts`` -- a list of polynomials or rational functions,
      all of which should have the same parent.

    - ``domain`` -- an affine or projective scheme, or product of
      projective schemes, on which ``polys`` defines an endomorphism.
      Subschemes are also ok.

    EXAMPLES::

        sage: A.<x> = AffineSpace(QQ,1)
        sage: f = DynamicalSystem_affine([x^2+1])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.affine_ds.DynamicalSystem_affine_field'>

        ::

        sage: P.<x,y> = ProjectiveSpace(QQ,1)
        sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
        sage: type(f)
        <class 'sage.dynamics.arithmetic_dynamics.projective_ds.DynamicalSystem_projective_field'>
    """

    def __init__(self, polys_or_rat_fncts, domain):
        r"""
        The Python constructor.

        EXAMPLES::

            sage: from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem_generic
            sage: P.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^2+y^2, y^2])
            sage: isinstance(f, DynamicalSystem_generic)
            True
        """
        H = End(domain)
        # All consistency checks are done by the public class constructors,
        # so we can set check=False here.
        SchemeMorphism_polynomial.__init__(self, H, polys_or_rat_fncts, check=False)

    def _repr_type(self):
        r"""
        Return a string representation of the type of a dynamical system.

        OUTPUT: string.

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_type()
            'Dynamical System'
        """
        return "Dynamical System"

    def _repr_(self):
        r"""
        Return a string representation of a dynamical system.

        OUTPUT: string.

        EXAMPLES::

            sage: PS.<x,y> = ProjectiveSpace(QQ,1)
            sage: f = DynamicalSystem_projective([x^3, x*y^2])
            sage: f._repr_()
            'Dynamical System of Projective Space of dimension 1 over Rational Field\n
              Defn: Defined on coordinates by sending (x : y) to\n        (x^3 : x*y^2)'
        """
        s = "%s of %s"%(self._repr_type(), self.domain())
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    def as_scheme_morphism(self):
        """
        Return this dynamical system as :class:`SchemeMorphism_polynomial`.

        OUTPUT: :class:`SchemeMorphism_polynomial`.

        EXAMPLES::

            sage: P.<x,y,z> = ProjectiveSpace(ZZ, 2)
            sage: f = DynamicalSystem_projective([x^2, y^2, z^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space'>

        ::

            sage: P.<x,y> = ProjectiveSpace(QQ, 1)
            sage: f = DynamicalSystem_projective([x^2-y^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_field'>

        ::

            sage: P.<x,y> = ProjectiveSpace(GF(5), 1)
            sage: f = DynamicalSystem_projective([x^2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.projective.projective_morphism.SchemeMorphism_polynomial_projective_space_finite_field'>

        ::

            sage: A.<x,y> = AffineSpace(ZZ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space'>

        ::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_field'>

        ::

            sage: A.<x,y> = AffineSpace(GF(3), 2)
            sage: f = DynamicalSystem_affine([x^2-2, y^2])
            sage: type(f.as_scheme_morphism())
            <class 'sage.schemes.affine.affine_morphism.SchemeMorphism_polynomial_affine_space_finite_field'>
        """
        H = End(self.domain())
        return H(list(self))

    def change_ring(self, R, check=True):
        r"""
        Return a new dynamical system which is this map coerced to ``R``.

        If ``check`` is ``True``, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring or morphism.

        OUTPUT:

        - A new :class:`DynamicalSystem_projective` which is this map coerced to ``R``.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(ZZ, 1)
            sage: f = DynamicalSystem_projective([3*x^2, y^2])
            sage: f.change_ring(GF(5))
            Dynamical System of Projective Space of dimension 1 over Finite Field of size 5
              Defn: Defined on coordinates by sending (x : y) to
                    (-2*x^2 : y^2)
        """
        f = self.as_scheme_morphism()
        F = f.change_ring(R)
        return F.as_dynamical_system()

    def specialization(self, D=None, phi=None, homset=None):
        r"""
        Specialization of this dynamical system.

        Given a family of maps defined over a polynomial ring. A specialization
        is a particular member of that family. The specialization can be specified either
        by a dictionary or a :class:`SpecializationMorphism`.

        INPUT:

        - ``D`` -- dictionary (optional).

        - ``phi`` -- SpecializationMorphism (optional).

        - ``homset`` -- homset of specialized map (optional).

        OUTPUT: :class:`DynamicalSystem`.

        EXAMPLES::

            sage: R.<c> = PolynomialRing(QQ)
            sage: P.<x,y> = ProjectiveSpace(R, 1)
            sage: f = DynamicalSystem_projective([x^2 + c*y^2,y^2], domain=P)
            sage: f.specialization({c:1})
            Dynamical System of Projective Space of dimension 1 over Rational Field
                  Defn: Defined on coordinates by sending (x : y) to
                        (x^2 + y^2 : y^2)
        """
        F = self.as_scheme_morphism().specialization(D, phi, homset)
        return F.as_dynamical_system()
