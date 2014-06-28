r"""
Scheme morphism

.. NOTE::

    You should never create the morphisms directy. Instead, use the
    :meth:`~sage.schemes.generic.scheme.hom` and
    :meth:`~sage.structure.parent.Hom` methods that are inherited by
    all schemes.

If you want to extend the Sage library with some new kind of scheme,
your new class (say, ``myscheme``) should provide a method

* ``myscheme._morphism(*args, **kwds)`` returning a morphism
  between two schemes in your category, usually defined via
  polynomials. Your morphism class should derive from
  :class:`SchemeMorphism_polynomial`. These morphisms will usually be
  elements of the Hom-set
  :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

Optionally, you can also provide a special Hom-set class for your
subcategory of schemes. If you want to do this, you should also
provide a method

* ``myscheme._homset(*args, **kwds)`` returning a
  Hom-set, which must be an element of a derived class of
  `class:`~sage.schemes.generic.homset.SchemeHomset_generic`. If your
  new Hom-set class does not use ``myscheme._morphism`` then you
  do not have to provide it.

Note that points on schemes are morphisms `Spec(K)\to X`, too. But we
typically use a different notation, so they are implemented in a
different derived class. For this, you should implement a method

* ``myscheme._point(*args, **kwds)`` returning a point, that is,
  a morphism `Spec(K)\to X`. Your point class should derive from
  :class:`SchemeMorphism_point`.

Optionally, you can also provide a special Hom-set for the points, for
example the point Hom-set can provide a method to enumerate all
points. If you want to do this, you should also provide a method

* ``myscheme._point_homset(*args, **kwds)`` returning
  the :mod:`~sage.schemes.generic.homset` of points. The Hom-sets of
  points are implemented in classes named ``SchemeHomset_points_...``.
  If your new Hom-set class does not use ``myscheme._point`` then
  you do not have to provide it.

AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.

- Volker Braun (2011-08-08): Renamed classes, more documentation, misc
  cleanups.

- Ben Hutz (June 2012): added support for projective ring

- Simon King (2013-10): copy the changes of :class:`~sage.categories.morphism.Morphism`
  that have been introduced in :trac:`14711`.
"""

# Historical note: in trac #11599, V.B. renamed
# * _point_morphism_class -> _morphism
# * _homset_class -> _point_homset

#*****************************************************************************
#       Copyright (C) 2013 Simon King <simon.king@uni-jena.de>
#       Copyright (C) 2011 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.element   import AdditiveGroupElement, RingElement, Element, generic_power, parent
from sage.structure.sequence  import Sequence
from sage.categories.homset   import Homset, Hom
from sage.rings.all           import Integer
from sage.rings.commutative_ring import is_CommutativeRing
from sage.rings.morphism import is_RingHomomorphism
from point                    import is_SchemeTopologicalPoint
from sage.rings.infinity      import infinity
import scheme

from sage.rings.arith            import gcd, lcm
from sage.categories.gcd_domains import GcdDomains
from sage.rings.quotient_ring    import QuotientRing_generic
from sage.categories.map         import FormalCompositeMap, Map
from sage.misc.constant_function import ConstantFunction
from sage.categories.morphism    import SetMorphism

def is_SchemeMorphism(f):
    """
    Test whether ``f`` is a scheme morphism.

    INPUT:

    - ``f`` -- anything.

    OUTPUT:

    Boolean. Return ``True`` if ``f`` is a scheme morphism or a point
    on an elliptic curve.

    EXAMPLES::

        sage: A.<x,y> = AffineSpace(QQ,2); H = A.Hom(A)
        sage: f = H([y,x^2+y]); f
        Scheme endomorphism of Affine Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (y, x^2 + y)
        sage: from sage.schemes.generic.morphism import is_SchemeMorphism
        sage: is_SchemeMorphism(f)
        True
    """
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
    return isinstance(f, (SchemeMorphism, EllipticCurvePoint_field));


class SchemeMorphism(Element):
    """
    Base class for scheme morphisms

    INPUT:

    - ``parent`` -- the parent of the morphism.

    .. TODO::

        Currently, :class:`SchemeMorphism` copies code from
        :class:`~sage.categories.map.Map` rather than inheriting from it. This
        is to work around a bug in Cython: We want to create a common
        sub-class of :class:`~sage.structure.element.ModuleElement` and
        :class:`SchemeMorphism`, but Cython would currently confuse cpdef
        attributes of the two base classes. Proper inheritance should be used
        as soon as this bug is fixed.

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: Hom = X.Hom(X)
        sage: from sage.schemes.generic.morphism import SchemeMorphism
        sage: f = SchemeMorphism(Hom)
        sage: type(f)
        <class 'sage.schemes.generic.morphism.SchemeMorphism'>

    TESTS::

        sage: A2 = AffineSpace(QQ,2)
        sage: A2.structure_morphism().domain()
        Affine Space of dimension 2 over Rational Field
        sage: A2.structure_morphism().category()
        Category of hom sets in Category of schemes
    """

    def __init__(self, parent, codomain=None):
        """
        The Python constructor.

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: type(f)
            <class 'sage.schemes.generic.morphism.SchemeMorphism'>
        """
        if codomain is not None:
            parent = Hom(parent, codomain)
        if not isinstance(parent, Homset):
            raise TypeError("parent (=%s) must be a Homspace"%parent)
        Element.__init__(self, parent)
        self.domain = ConstantFunction(parent.domain())
        self._codomain = parent.codomain()
        self.codomain = ConstantFunction(self._codomain)

    # We copy methods of sage.categories.map.Map, to make
    # a future transition of SchemeMorphism to a sub-class of Morphism
    # easier.
    def __call__(self, x, *args, **kwds):
        """
        Do not override this method!

        For implementing application of maps, implement a method
        ``_call_(self, x)`` and/or a method ``_call_with_args(x, args, kwds)`.
        In these methods, you can assume that ``x`` belongs to the domain of
        this morphism, ``args`` is a tuple and ``kwds`` is a dict.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y,x^2+y])
            sage: f([2,3])    # indirect doctest
            (3, 7)

        An example with optional arguments::

            sage: PS.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(PS,PS)
            sage: f=H([x^3,x*y^2])
            sage: P=PS(0,1)
            sage: f(P,check=False)     # indirect doctest
            (0 : 0)
        """
        P = parent(x)
        D = self.domain()
        if P is D: # we certainly want to call _call_/with_args
            if not args and not kwds:
                return self._call_(x)
            return self._call_with_args(x, args, kwds)
        # Is there coercion?
        converter = D._internal_coerce_map_from(P)
        if converter is None:
            try:
                return self.pushforward(x,*args,**kwds)
            except (AttributeError, TypeError, NotImplementedError):
                pass # raise TypeError, "%s must be coercible into %s"%(x, self.domain())
            # Here, we would like to do
            ##try:
            ##    x = D(x).
            ##except (TypeError, NotImplementedError):
            ##    raise TypeError, "%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain())
            # However, this would involve a test whether x.codomain() ==
            # self. This would trigger a Groebner basis computation, that
            # (1) could be slow and (2) could involve an even slower toy
            # implementation, resulting in a warning.
            #
            # Contract: If x is a scheme morphism point, then _call_ knows
            # what to do with it (e.g., use the _coords attribute). Otherwise,
            # we can try a conversion into the domain (e.g., if x is a list),
            # WITHOUT to trigger a Groebner basis computation.
            if kwds.get('check', True):
                if not isinstance(x, SchemeMorphism_point):
                    try:
                        x = D(x)
                    except (TypeError, NotImplementedError):
                        raise TypeError, "%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain())
                elif self.domain()!=x.codomain():
                    raise TypeError, "%s fails to convert into the map's domain %s, but a `pushforward` method is not properly implemented"%(x, self.domain())
        else:
            x = converter(x)
        if not args and not kwds:
            return self._call_(x)
        return self._call_with_args(x, args, kwds)

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: f._repr_defn()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def _repr_type(self):
        r"""
        Return a string representation of the type of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A2.structure_morphism()  # indirect doctest
            Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map
        """
        return "Scheme"

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: X = Spec(ZZ)
            sage: Hom = X.Hom(X)
            sage: from sage.schemes.generic.morphism import SchemeMorphism
            sage: f = SchemeMorphism(Hom)
            sage: f._repr_()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if self.is_endomorphism():
            s = "%s endomorphism of %s"%(self._repr_type(), self.domain())
        else:
            s = "%s morphism:"%self._repr_type()
            s += "\n  From: %s"%self.domain()
            s += "\n  To:   %s"%self._codomain
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    def __mul__(self, right):
        """
        We can currently only multiply scheme morphisms.

        If one factor is an identity morphism, the other is returned.
        Otherwise, a formal composition of maps obtained from the scheme
        morphisms is returned.

        EXAMPLES:

        Identity maps do not contribute to the product::

            sage: X = AffineSpace(QQ,2)
            sage: id = X.identity_morphism()
            sage: id^0    # indirect doctest
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
            sage: id^2
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map

        Here, we see a formal composition::

            sage: X = AffineSpace(QQ,2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
             Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring

        Of course, the codomain of the first factor must coincide with the
        domain of the second factor::

            sage: f * g
            Traceback (most recent call last):
            ...
            TypeError: self (=Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map) domain must equal right (=Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map) codomain
        """
        if not isinstance(right, SchemeMorphism):
            raise TypeError, "right (=%s) must be a SchemeMorphism to multiply it by %s"%(right, self)
        if right.codomain() != self.domain():
            raise TypeError, "self (=%s) domain must equal right (=%s) codomain"%(self, right)
        if isinstance(self, SchemeMorphism_id):
            return right
        if isinstance(right, SchemeMorphism_id):
            return self
        return self._composition(right)


    def __pow__(self, n, dummy=None):
        """
        Exponentiate an endomorphism.

        INPUT:

        - ``n`` -- integer. The exponent.

        OUTPUT:

        A composite map that belongs to the same endomorphism set as ``self``.

        EXAMPLES::

            sage: X = AffineSpace(QQ,2)
            sage: id = X.identity_morphism()
            sage: id^0
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
            sage: id^2
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Identity map
        """
        if not self.is_endomorphism():
            raise TypeError, "self must be an endomorphism."
        if n==0:
            return self.domain().identity_morphism()
        return generic_power(self, n)

    def category(self):
        """
        Return the category of the Hom-set.

        OUTPUT:

        A category.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A2.structure_morphism().category()
            Category of hom sets in Category of schemes
        """
        return self.parent().category()

    def category_for(self):
        """
        Return the category which this morphism belongs to.

        EXAMPLES::

            sage: A2 = AffineSpace(QQ,2)
            sage: A2.structure_morphism().category_for()
            Category of schemes
        """
        return self.parent().homset_category()

    def is_endomorphism(self):
        """
        Return wether the morphism is an endomorphism.

        OUTPUT:

        Boolean. Whether the domain and codomain are identical.

        EXAMPLES::

            sage: X = AffineSpace(QQ,2)
            sage: X.structure_morphism().is_endomorphism()
            False
            sage: X.identity_morphism().is_endomorphism()
            True
        """
        return self.parent().is_endomorphism_set()

    def _composition(self, right):
        """
        A helper for multiplying maps by composition.

        .. WARNING::

            Do not override this method! Override :meth:`_composition_`
            instead.

        EXAMPLES::

            sage: X = AffineSpace(QQ,2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
             Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring

            sage: f * g
            Traceback (most recent call last):
            ...
            TypeError: self (=Scheme morphism:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Rational Field
              Defn: Structure map) domain must equal right (=Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map) codomain
        """
        category = self.category_for()._meet_(right.category_for())
        H = Hom(right.domain(), self._codomain, category)
        return self._composition_(right, H)

    def _composition_(self, right, homset):
        """
        Helper to construct the composition of two morphisms.

        Override this if you want to have a different behaviour of composition

        INPUT:

        - ``right`` -- a map or callable
        - ``homset`` -- a homset containing the composed map

        OUTPUT:

        An element of ``homset``. The output is obtained by converting the
        arguments to :class:`~sage.categories.morphism.SetMorphism` if
        necessary, and then forming a :class:`~sage.categories.map.FormalCompositeMap`

        EXAMPLES::

            sage: X = AffineSpace(QQ,2)
            sage: f = X.structure_morphism()
            sage: Y = Spec(QQ)
            sage: g = Y.structure_morphism()
            sage: g * f    # indirect doctest
             Composite map:
              From: Affine Space of dimension 2 over Rational Field
              To:   Spectrum of Integer Ring
              Defn:   Generic morphism:
                      From: Affine Space of dimension 2 over Rational Field
                      To:   Spectrum of Rational Field
                    then
                      Generic morphism:
                      From: Spectrum of Rational Field
                      To:   Spectrum of Integer Ring
        """
        if not isinstance(right, Map):
            right = SetMorphism(right.parent(), right)
        return FormalCompositeMap(homset, right, SetMorphism(self.parent(),self))

    def glue_along_domains(self, other):
        r"""
        Glue two morphism

        INPUT:

        - ``other`` -- a scheme morphism with the same domain.

        OUTPUT:

        Assuming that self and other are open immersions with the same
        domain, return scheme obtained by gluing along the images.

        EXAMPLES:

        We construct a scheme isomorphic to the projective line over
        `\mathrm{Spec}(\QQ)` by gluing two copies of `\mathbb{A}^1`
        minus a point::

            sage: R.<x,y> = PolynomialRing(QQ, 2)
            sage: S.<xbar, ybar> = R.quotient(x*y - 1)
            sage: Rx = PolynomialRing(QQ, 'x')
            sage: i1 = Rx.hom([xbar])
            sage: Ry = PolynomialRing(QQ, 'y')
            sage: i2 = Ry.hom([ybar])
            sage: Sch = Schemes()
            sage: f1 = Sch(i1)
            sage: f2 = Sch(i2)

        Now f1 and f2 have the same domain, which is a
        `\mathbb{A}^1` minus a point. We glue along the domain::

            sage: P1 = f1.glue_along_domains(f2)
            sage: P1
            Scheme obtained by gluing X and Y along U, where
              X: Spectrum of Univariate Polynomial Ring in x over Rational Field
              Y: Spectrum of Univariate Polynomial Ring in y over Rational Field
              U: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
              over Rational Field by the ideal (x*y - 1)

            sage: a, b = P1.gluing_maps()
            sage: a
            Affine Scheme morphism:
             From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
                   over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over
                            Rational Field by the ideal (x*y - 1)
                      Defn: x |--> xbar
            sage: b
            Affine Scheme morphism:
              From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y
                    over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in y over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in y over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over
                            Rational Field by the ideal (x*y - 1)
                      Defn: y |--> ybar
        """
        import glue
        return glue.GluedScheme(self, other)

class SchemeMorphism_id(SchemeMorphism):
    """
    Return the identity morphism from `X` to itself.

    INPUT:

    - ``X`` -- the scheme.

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: X.identity_morphism()  # indirect doctest
        Scheme endomorphism of Spectrum of Integer Ring
          Defn: Identity map
    """
    def __init__(self, X):
        """
        The Python constructor.

        See :class:`SchemeMorphism_id` for details.

        TESTS::

            sage: Spec(ZZ).identity_morphism()
            Scheme endomorphism of Spectrum of Integer Ring
              Defn: Identity map
        """
        SchemeMorphism.__init__(self, X.Hom(X))

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: Spec(ZZ).identity_morphism()._repr_defn()
            'Identity map'
        """
        return 'Identity map'


class SchemeMorphism_structure_map(SchemeMorphism):
    r"""
    The structure morphism

    INPUT:

    - ``parent`` -- Hom-set with codomain equal to the base scheme of
      the domain.

    EXAMPLES::

        sage: Spec(ZZ).structure_morphism()    # indirect doctest
        Scheme morphism:
          From: Spectrum of Integer Ring
          To:   Spectrum of Integer Ring
          Defn: Structure map
    """
    def __init__(self, parent, codomain=None):
        """
        The Python constuctor.

        See :class:`SchemeMorphism_structure_map` for details.

        TESTS::

            sage: from sage.schemes.generic.morphism import SchemeMorphism_structure_map
            sage: SchemeMorphism_structure_map( Spec(QQ).Hom(Spec(ZZ)) )
            Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Structure map
        """
        SchemeMorphism.__init__(self, parent, codomain=None)
        if self.domain().base_scheme() != self._codomain:
            raise ValueError("parent must have codomain equal the base scheme of domain.")

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: Spec(ZZ).structure_morphism()._repr_defn()
            'Structure map'
        """
        return 'Structure map'


class SchemeMorphism_spec(SchemeMorphism):
    """
    Morphism of spectra of rings

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are affine schemes.

    - ``phi`` -- a ring morphism with matching domain and codomain.

    - ``check`` -- boolean (optional, default:``True``). Whether to
      check the input for consistency.

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: phi = R.hom([QQ(7)]); phi
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7

        sage: X = Spec(QQ); Y = Spec(R)
        sage: f = X.hom(phi); f
        Affine Scheme morphism:
          From: Spectrum of Rational Field
          To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
          Defn: Ring morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Rational Field
                  Defn: x |--> 7

        sage: f.ring_homomorphism()
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7
    """
    def __init__(self, parent, phi, check=True):
        """
        The Python constuctor.

        See :class:`SchemeMorphism_structure_map` for details.

        TESTS::

            sage: from sage.schemes.generic.morphism import SchemeMorphism_spec
            sage: SchemeMorphism_spec(Spec(QQ).Hom(Spec(ZZ)), ZZ.hom(QQ))
            Affine Scheme morphism:
              From: Spectrum of Rational Field
              To:   Spectrum of Integer Ring
              Defn: Ring Coercion morphism:
                      From: Integer Ring
                      To:   Rational Field
        """
        SchemeMorphism.__init__(self, parent)
        if check:
            if not is_RingHomomorphism(phi):
                raise TypeError("phi (=%s) must be a ring homomorphism" % phi)
            if phi.domain() != parent.codomain().coordinate_ring():
                raise TypeError("phi (=%s) must have domain %s"
                                % (phi, parent.codomain().coordinate_ring()))
            if phi.codomain() != parent.domain().coordinate_ring():
                raise TypeError("phi (=%s) must have codomain %s"
                                % (phi, parent.domain().coordinate_ring()))
        self.__ring_homomorphism = phi

    def _call_(self, x):
        r"""
        Make morphisms callable.

        INPUT:

        - ``x`` -- a scheme point.

        OUTPUT:

        The image scheme point.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f(X.an_element())    # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: Point on Spectrum of Rational Field defined by the
            Principal ideal (0) of Rational Field fails to convert into the
            map's domain Spectrum of Rational Field, but a `pushforward`
            method is not properly implemented
        """
        # By virtue of argument preprocessing in __call__, we can assume that
        # x is a topological scheme point of self
        S = self.ring_homomorphism().inverse_image(x.prime_ideal())
        return self._codomain(S)

    def _repr_type(self):
        r"""
        Return a string representation of the type of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f._repr_type()
            'Affine Scheme'
        """
        return "Affine Scheme"

    def _repr_defn(self):
        r"""
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: print f._repr_defn()
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
              Defn: x |--> 7
        """
        return repr(self.ring_homomorphism())

    def ring_homomorphism(self):
        """
        Return the underlying ring homomorphism.

        OUTPUT:

        A ring homomorphism.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ)
            sage: phi = R.hom([QQ(7)])
            sage: X = Spec(QQ); Y = Spec(R)
            sage: f = X.hom(phi)
            sage: f.ring_homomorphism()
            Ring morphism:
              From: Univariate Polynomial Ring in x over Rational Field
              To:   Rational Field
              Defn: x |--> 7
        """
        return self.__ring_homomorphism


############################################################################
# Morphisms between schemes given on points
# The _affine and _projective below refer to the CODOMAIN.
# The domain can be either affine or projective regardless
# of the class
############################################################################
class SchemeMorphism_polynomial(SchemeMorphism):
    """
    A morphism of schemes determined by polynomials that define what
    the morphism does on points in the ambient space.

    INPUT:

    - ``parent`` -- Hom-set whose domain and codomain are affine schemes.

    - ``polys`` -- a list/tuple/iterable of polynomials defining the
      scheme morphism.

    - ``check`` -- boolean (optional, default:``True``). Whether to
      check the input for consistency.

    EXAMPLES:

    An example involving the affine plane::

        sage: R.<x,y> = QQ[]
        sage: A2 = AffineSpace(R)
        sage: H = A2.Hom(A2)
        sage: f = H([x-y, x*y])
        sage: f([0,1])
        (-1, 0)

    An example involving the projective line::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x^2+y^2,x*y])
        sage: f([0,1])
        (1 : 0)

    Some checks are performed to make sure the given polynomials
    define a morphism::

        sage: f = H([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of
        Multivariate Polynomial Ring in x, y over Rational Field
    """
    def __init__(self, parent, polys, check=True):
        """
        The Python constructor.

        See :class:`SchemeMorphism_polynomial` for details.

        EXAMPLES::

            sage: A2.<x,y> = AffineSpace(QQ,2)
            sage: H = A2.Hom(A2)
            sage: H([x-y, x*y])
            Scheme endomorphism of Affine Space of dimension 2 over Rational Field
              Defn: Defined on coordinates by sending (x, y) to
                    (x - y, x*y)
        """
        if check:
            if not isinstance(polys, (list, tuple)):
                raise TypeError("polys (=%s) must be a list or tuple"%polys)
            source_ring = parent.domain().coordinate_ring()
            target = parent._codomain.ambient_space()
            if len(polys) != target.ngens():
                raise ValueError("there must be %s polynomials"%target.ngens())
            try:
                polys = [source_ring(poly) for poly in polys]
            except TypeError:
                raise TypeError("polys (=%s) must be elements of %s"%(polys,source_ring))
            if isinstance(source_ring, QuotientRing_generic):
                lift_polys = [f.lift() for f in polys]
            else:
                lift_polys = polys
            polys = Sequence(lift_polys)
        self._polys = polys
        SchemeMorphism.__init__(self, parent)

    def defining_polynomials(self):
        """
        Return the defining polynomials.

        OUTPUT:

        An immutable sequence of polynomials that defines this scheme
        morphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: H([x^3+y, 1-x-y]).defining_polynomials()
            [x^3 + y, -x - y + 1]
        """
        return self._polys

    def _call_(self, x):
        """
        Apply this morphism to a point in the domain.

        INPUT:

        - ``x`` -- a point in the domain or a list or tuple that defines a point in the domain.

        OUTPUT:

        A point in the codomain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y,x^2+y])
            sage: f([2,3])    # indirect doctest
            (3, 7)

        An example with algebraic schemes::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme(x)
            sage: Y = A.subscheme(y)
            sage: Hom_XY = X.Hom(Y)
            sage: f = Hom_XY([y,0])   # (0,y) |-> (y,0)
            sage: f
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y
              Defn: Defined on coordinates by sending (x, y) to
                    (y, 0)
            sage: f([0,3])
            (3, 0)

        The input must be convertible into the map's domain::

            sage: f(0)
            Traceback (most recent call last):
            ...
            TypeError: 0 fails to convert into the map's domain Closed
            subscheme of Affine Space of dimension 2 over Rational Field
            defined by:
              x, but a `pushforward` method is not properly implemented

        ::

        It is possible to avoid the checks on the resulting point which can be
        useful for indeterminacies, but be careful!!

            sage: PS.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(PS,PS)
            sage: f=H([x^3,x*y^2])
            sage: P=PS(0,1)
            sage: f(P,check=False)
            (0 : 0)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2]);
            sage: f([4,4,1])
            (16 : 16 : 1)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2]);
            sage: f(P([4,4,1]))
            Traceback (most recent call last):
            ...
            TypeError: (4 : 4 : 1) fails to convert into the map's domain
            Closed subscheme of Projective Space of dimension 2 over Integer
            Ring defined by:
              x^2 - y^2, but a `pushforward` method is not properly implemented
        """
        # Checks were done in __call__
        P = [f(x._coords) for f in self.defining_polynomials()]
        return self._codomain.point(P,check=True)

    def _call_with_args(self, x, args, kwds):
        """
        Apply this morphism to a point in the domain, with additional arguments

        INPUT:

        - ``x`` -- a point in the domain or a list or tuple that defines a point in the domain.
        - ``check``, a boolean, either provided by position or name.

        OUTPUT:

        A point in the codomain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y,x^2+y])
            sage: f([2,3])
            (3, 7)

        An example with algebraic schemes::

            sage: A.<x,y> = AffineSpace(QQ, 2)
            sage: X = A.subscheme(x)
            sage: Y = A.subscheme(y)
            sage: Hom_XY = X.Hom(Y)
            sage: f = Hom_XY([y,0])   # (0,y) |-> (y,0)
            sage: f
            Scheme morphism:
              From: Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              x
              To:   Closed subscheme of Affine Space of dimension 2 over Rational Field defined by:
              y
              Defn: Defined on coordinates by sending (x, y) to
                    (y, 0)
            sage: f([0,3])
            (3, 0)

        As usual, if the input does not belong to a map's domain, it is first
        attempted to convert it::

            sage: f(0)
            Traceback (most recent call last):
            ...
            TypeError: 0 fails to convert into the map's domain Closed subscheme of
            Affine Space of dimension 2 over Rational Field defined by:
              x, but a `pushforward` method is not properly implemented

        It is possible to avoid the checks on the resulting point which can be
        useful for indeterminacies, but be careful!!

        ::

            sage: PS.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(PS,PS)
            sage: f=H([x^3,x*y^2])
            sage: P=PS(0,1)
            sage: f(P,check=False)     # indirect doctest
            (0 : 0)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2]);
            sage: f([4,4,1])
            (16 : 16 : 1)

        ::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: H=Hom(X,X)
            sage: f=H([x^2,y^2,z^2]);
            sage: f(P([4,4,1]))
            Traceback (most recent call last):
            ...
            TypeError: (4 : 4 : 1) fails to convert into the map's domain
            Closed subscheme of Projective Space of dimension 2 over
            Integer Ring defined by:
              x^2 - y^2, but a `pushforward` method is not properly implemented
        """
        if args:
            check = args[0]
        else:
            check = kwds.get("check", False)
        # containment of x in the domain has already been checked, in __call__
        P = [f(x._coords) for f in self.defining_polynomials()]
        return self._codomain.point(P,check)


    def _repr_defn(self):
        """
        Return a string representation of the definition of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y,x^2+y])
            sage: print f._repr_defn()
            Defined on coordinates by sending (x, y) to
            (y, x^2 + y)
        """
        i = self.domain().ambient_space()._repr_generic_point()
        o = self._codomain.ambient_space()._repr_generic_point(self.defining_polynomials())
        return "Defined on coordinates by sending %s to\n%s"%(i,o)

    def __getitem__(self,i):
        """
        returns the ith poly with self[i]

        INPUT::

        - ``i``- integer

        OTUPUT::

        - element of the coordinate ring of the domain

        Examples::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(P,P)
            sage: f=H([3/5*x^2,6*y^2])
            sage: f[1]
            6*y^2
        """
        return(self._polys[i])

    def __copy__(self):
        r"""
        Returns a copy of ``self``.

        OUTPUT:

        - :class:`SchemeMorphism_polynomial`

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(P,P)
            sage: f=H([3/5*x^2,6*y^2])
            sage: g =copy(f)
            sage: f==g
            True
            sage: f is g
            False

        ::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: X=P.subscheme(x^2-y^2);
            sage: Q=X(23,23,46)
            sage: P=X(1,1,1)
            sage: P!=Q
            True
        """
        return self.parent()(self._polys)

    def base_ring(self):
        r"""
        Return the base ring of ``self``, that is, the ring over which the coefficients
        of ``self`` is given as polynomials.

        OUTPUT:

        - ring

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(P,P)
            sage: f=H([3/5*x^2,6*y^2])
            sage: f.base_ring()
            Rational Field

        ::

            sage: R.<t>=PolynomialRing(ZZ,1)
            sage: P.<x,y>=ProjectiveSpace(R,1)
            sage: H=Hom(P,P)
            sage: f=H([3*x^2,y^2])
            sage: f.base_ring()
            Multivariate Polynomial Ring in t over Integer Ring
        """
        return(self.domain().base_ring())

    def coordinate_ring(self):
        r"""
        Returns the coordinate ring of the ambient projective space
        the multivariable polynomial ring over the base ring

        OUTPUT:

        - ring

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: H=Hom(P,P)
            sage: f=H([3/5*x^2,6*y^2])
            sage: f.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Rational Field

        ::

            sage: R.<t>=PolynomialRing(ZZ,1)
            sage: P.<x,y>=ProjectiveSpace(R,1)
            sage: H=Hom(P,P)
            sage: f=H([3*x^2,y^2])
            sage: f.coordinate_ring()
            Multivariate Polynomial Ring in x, y over Multivariate Polynomial Ring
            in t over Integer Ring
        """
        return(self._polys[0].parent())

    def change_ring(self,R, check=True):
        r"""
        Returns a new :class:`SchemeMorphism_polynomial` which is ``self`` coerced to `R`. If ``check``
        is ``True``, then the initialization checks are performed.

        INPUT:

        - ``R`` -- ring

        - ``check`` -- Boolean

        OUTPUT:

        - element of the coordinate ring of the domain

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(ZZ,1)
            sage: H=Hom(P,P)
            sage: f=H([3*x^2,y^2])
            sage: f.change_ring(GF(3))
            Traceback (most recent call last):
            ...
            ValueError: polys (=[0, y^2]) must be of the same degree

        ::

            sage: P.<x,y,z>=ProjectiveSpace(QQ,2)
            sage: H=Hom(P,P)
            sage: f=H([5/2*x^3 + 3*x*y^2-y^3,3*z^3 + y*x^2, x^3-z^3])
            sage: f.change_ring(GF(3))
            Scheme endomorphism of Projective Space of dimension 2 over Finite Field of size 3
              Defn: Defined on coordinates by sending (x : y : z) to
                    (x^3 - y^3 : x^2*y : x^3 - z^3)

        ::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: X=P.subscheme([5*x^2-y^2])
            sage: H=Hom(X,X)
            sage: f=H([x,y])
            sage: f.change_ring(GF(3))
            Scheme endomorphism of Closed subscheme of Projective Space of dimension
            1 over Finite Field of size 3 defined by:
              -x^2 - y^2
              Defn: Defined on coordinates by sending (x : y) to
                    (x : y)
            """
        F=self._polys
        S=self._codomain.change_ring(R)
        H=Hom(S,S)
        G=[]
        for i in range(self._codomain.ambient_space().dimension_relative()+1):
            G.append(F[i].change_ring(R))
        return(H(G,check))


############################################################################
# Rational points on schemes, which we view as morphisms determined
# by coordinates.
############################################################################

class SchemeMorphism_point(SchemeMorphism):
    """
    Base class for rational points on schemes.

    Recall that the `K`-rational points of a scheme `X` over `k` can
    be identified with the set of morphisms `Spec(K) \to X`. In Sage,
    the rational points are implemented by such scheme morphisms.

    EXAMPLES::

        sage: from sage.schemes.generic.morphism import SchemeMorphism
        sage: f = SchemeMorphism(Spec(ZZ).Hom(Spec(ZZ)))
        sage: type(f)
        <class 'sage.schemes.generic.morphism.SchemeMorphism'>
    """
    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a._repr_()
            '(1, 2)'
        """
        return self._codomain.ambient_space()._repr_generic_point(self._coords)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        OUTPUT:

        String.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: latex(a) == a._latex_()
            True
            sage: a._latex_()
            '\\left(1, 2\\right)'
        """
        return self._codomain.ambient_space()._latex_generic_point(self._coords)

    def __getitem__(self, n):
        """
        Return the ``n``-th coordinate.

        OUTPUT:

        The coordinate values as an element of the base ring.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a[0]
            1
            sage: a[1]
            2
        """
        return self._coords[n]

    def __iter__(self):
        """
        Iterate over the coordinates of the point.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: iter = a.__iter__()
            sage: iter.next()
            1
            sage: iter.next()
            2
            sage: list(a)
            [1, 2]
        """
        return iter(self._coords)

    def __tuple__(self):
        """
        Return the coordinates as a tuple.

        OUTPUT:

        A tuple.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: tuple(a)
            (1, 2)
        """
        return self._coords

    def __len__(self):
        """
        Return the number of coordinates.

        OUTPUT:

        Integer. The number of coordinates used to describe the point.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: len(a)
            2
        """
        return len(self._coords)

    def __cmp__(self, other):
        """
        Compare two scheme morphisms.

        INPUT:

        - ``other`` -- anything. To compare against the scheme
          morphism ``self``.

        OUTPUT:

        ``+1``, ``0``, or ``-1``.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: b = A(3,4)
            sage: a.__cmp__(b)
            -1
            sage: a != b
            True
        """
        if not isinstance(other, SchemeMorphism_point):
            try:
                other = self._codomain.ambient_space()(other)
            except TypeError:
                return -1
        return cmp(self._coords, other._coords)

    def scheme(self):
        """
        Return the scheme whose point is represented.

        OUTPUT:

        A scheme.

        EXAMPLES::

            sage: A = AffineSpace(2, QQ)
            sage: a = A(1,2)
            sage: a.scheme()
            Affine Space of dimension 2 over Rational Field
        """
        return self._codomain

    def change_ring(self,R, check=True):
        r"""
        Returns a new :class:`SchemeMorphism_point` which is self coerced to R. If `check`
        is true, then the initialization checks are performed.

        INPUT:

        - ``R`` -- a ring

        - ``check`` -- Boolean (optional)

        OUTPUT:

        - :class:`SchemeMorphism_point`

        EXAMPLES::

            sage: P.<x,y,z>=ProjectiveSpace(ZZ,2)
            sage: X=P.subscheme(x^2-y^2)
            sage: X(23,23,1).change_ring(GF(13))
            (10 : 10 : 1)

        ::

            sage: P.<x,y>=ProjectiveSpace(QQ,1)
            sage: P(-2/3,1).change_ring(CC)
            (-0.666666666666667 : 1.00000000000000)

        ::

            sage: P.<x,y>=ProjectiveSpace(ZZ,1)
            sage: P(152,113).change_ring(Zp(5))
            (2 + 5^2 + 5^3 + O(5^20) : 3 + 2*5 + 4*5^2 + O(5^20))
        """
        S=self._codomain.change_ring(R)
        Q=[]
        for i in range(len(self._coords)):
            Q.append(R(self._coords[i]))
        return(S.point(Q,check))

    def __copy__(self):
        r"""
        Returns a copy of the :class:`SchemeMorphism_point` self coerced to `R`.

        OUTPUT:

        - :class:`SchemeMorphism_point`

        EXAMPLES::

            sage: P.<x,y>=ProjectiveSpace(ZZ,1)
            sage: Q=P(152,113)
            sage: copy(Q) is Q
            False
            sage: copy(Q) == Q
            True
        """
        return(self._codomain.point(self._coords,False))
