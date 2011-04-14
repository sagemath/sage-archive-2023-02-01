"""
Scheme morphism

AUTHORS:

- David Kohel, William Stein

- William Stein (2006-02-11): fixed bug where P(0,0,0) was allowed as
  a projective point.
"""

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#  Copyright (C) 2006 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import infinity
from sage.structure.element   import AdditiveGroupElement, RingElement, Element, generic_power
from sage.structure.sequence  import Sequence
from sage.categories.morphism import Morphism
from sage.categories.homset   import Homset
from sage.rings.all           import is_RingHomomorphism, is_CommutativeRing, Integer
from point                    import is_SchemeTopologicalPoint
import scheme
import spec

def is_SchemeMorphism(f):
    """
    Return True if f is a scheme morphism or a point on an elliptic curve.

    EXAMPLES::

        sage: R.<x,y> = QQ[]; A.<x,y> = AffineSpace(R); H = A.Hom(A)
        sage: f = H([y,x^2+y])
        sage: sage.schemes.generic.morphism.is_SchemeMorphism(f)
        True
    """
    from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field # TODO: fix circular ref.
    return isinstance(f, (SchemeMorphism, EllipticCurvePoint_field));


class PyMorphism(Element):
    def __init__(self, parent):
        if not isinstance(parent, Homset):
            raise TypeError, "parent (=%s) must be a Homspace"%parent
        Element.__init__(self, parent)
        self._domain = parent.domain()
        self._codomain = parent.codomain()

    def _repr_type(self):
        return "Generic"

    def _repr_defn(self):
        return ""

    def _repr_(self):
        if self.is_endomorphism():
            s = "%s endomorphism of %s"%(self._repr_type(), self.domain())
        else:
            s = "%s morphism:"%self._repr_type()
            s += "\n  From: %s"%self.domain()
            s += "\n  To:   %s"%self.codomain()
        d = self._repr_defn()
        if d != '':
            s += "\n  Defn: %s"%('\n        '.join(self._repr_defn().split('\n')))
        return s

    def domain(self):
        return self._domain

    def codomain(self):
        return self.parent().codomain()

    def category(self):
        return self.parent().category()

    def is_endomorphism(self):
        return self.parent().is_endomorphism_set()

    def _composition_(self, right, homset):
        return FormalCompositeMorphism(homset, right, self)

    def __pow__(self, n, dummy=None):
        if not self.is_endomorphism():
            raise TypeError, "self must be an endomorphism."
        # todo -- what about the case n=0 -- need to specify the identity map somehow.
        return generic_power(self, n)



class SchemeMorphism(PyMorphism):
    """
    A scheme morphism
    """
    def __init__(self, parent):
        PyMorphism.__init__(self, parent)

    def _repr_type(self):
        return "Scheme"

    def glue_along_domains(self, other):
        r"""
        Assuming that self and other are open immersions with the same
        domain, return scheme obtained by gluing along the images.

        EXAMPLES:

        We construct a scheme isomorphic to the projective line over
        `\mathrm{Spec}(\QQ)` by gluing two copies of `\mathbb{A}^1`
        minus a point.

        ::

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
              U: Spectrum of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)

        ::

            sage: a, b = P1.gluing_maps()
            sage: a
            Affine Scheme morphism:
             From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in x over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
                      Defn: x |--> xbar
            sage: b
            Affine Scheme morphism:
              From: Spectrum of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
              To:   Spectrum of Univariate Polynomial Ring in y over Rational Field
              Defn: Ring morphism:
                      From: Univariate Polynomial Ring in y over Rational Field
                      To:   Quotient of Multivariate Polynomial Ring in x, y over Rational Field by the ideal (x*y - 1)
                      Defn: y |--> ybar
        """
        import glue
        return glue.GluedScheme(self, other)

class SchemeMorphism_id(SchemeMorphism):
    """
    Return the identity morphism from X to itself.

    EXAMPLES::

        sage: X = Spec(ZZ)
        sage: X.identity_morphism()
        Scheme endomorphism of Spectrum of Integer Ring
          Defn: Identity map
    """
    def __init__(self, X):
        SchemeMorphism.__init__(self, X.Hom(X))

    def _repr_defn(self):
        return "Identity map"


class SchemeMorphism_structure_map(SchemeMorphism):
    def __init__(self, parent):
        """
        INPUT:

        - ``parent`` - homset with codomain equal to the base
          scheme of the domain.
        """
        SchemeMorphism.__init__(self, parent)
        if self.domain().base_scheme() != self.codomain():
            raise ValueError, "parent must have codomain equal the base scheme of domain."

    def _repr_defn(self):
        return "Structure map"


class SchemeMorphism_spec(SchemeMorphism):
    """
    A morphism of spectrums of rings

    EXAMPLES::

        sage: R.<x> = PolynomialRing(QQ)
        sage: phi = R.hom([QQ(7)]); phi
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7

    ::

        sage: X = Spec(QQ); Y = Spec(R)
        sage: f = X.hom(phi); f
        Affine Scheme morphism:
          From: Spectrum of Rational Field
          To:   Spectrum of Univariate Polynomial Ring in x over Rational Field
          Defn: Ring morphism:
                  From: Univariate Polynomial Ring in x over Rational Field
                  To:   Rational Field
                  Defn: x |--> 7

    ::

        sage: f.ring_homomorphism()
        Ring morphism:
          From: Univariate Polynomial Ring in x over Rational Field
          To:   Rational Field
          Defn: x |--> 7
    """
    def __init__(self, parent, phi, check=True):
        SchemeMorphism.__init__(self, parent)
        if check:
            if not is_RingHomomorphism(phi):
                raise TypeError, "phi (=%s) must be a ring homomorphism"%phi
            if phi.domain() != parent.codomain().coordinate_ring():
                raise TypeError, "phi (=%s) must have domain %s"%(phi,
                                                   parent.codomain().coordinate_ring())
            if phi.codomain() != parent.domain().coordinate_ring():
                raise TypeError, "phi (=%s) must have codomain %s"%(phi,
                                                 parent.domain().coordinate_ring())
        self.__ring_homomorphism = phi

    def __call__(self, P):
        if not is_SchemeTopologicalPoint(P) and P in self.domain():
            raise TypeError, "P (=%s) must be a topological scheme point of %s"%(P, self)
        S = self.ring_homomorphism().inverse_image(P.prime_ideal())
        return self.codomain()(S)

    def _repr_type(self):
        return "Affine Scheme"

    def _repr_defn(self):
        return repr(self.ring_homomorphism())


    def ring_homomorphism(self):
        return self.__ring_homomorphism


############################################################################
# Morphisms between schemes given on points
# The _affine and _projective below refer to the CODOMAIN.
# The domain can be either affine or projective regardless
# of the class
############################################################################

class SchemeMorphism_on_points(SchemeMorphism):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient space.

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

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x^2, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x^2, x*y]) must not have common factors

        sage: f = H([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of Multivariate Polynomial Ring in x, y over Rational Field
    """
    def __init__(self, parent, polys, check=True):
        if check:
            if not isinstance(polys, (list, tuple)):
                raise TypeError, "polys (=%s) must be a list or tuple"%polys
            source_ring = parent.domain().coordinate_ring()
            target = parent.codomain().ambient_space()
            if len(polys) != target.ngens():
                raise ValueError, "there must be %s polynomials"%target.ngens()
            try:
                polys = [source_ring(poly) for poly in polys]
            except TypeError:
                raise TypeError, "polys (=%s) must be elements of %s"%(polys,source_ring)
            from sage.schemes.generic.projective_space import is_ProjectiveSpace
            if is_ProjectiveSpace(target):
                # if the codomain is a subscheme of projective space,
                # then we need to make sure that polys have no common
                # zeros
                from sage.rings.quotient_ring import QuotientRing_generic
                if isinstance(source_ring, QuotientRing_generic):
                    # if the coordinate ring of the domain is a
                    # quotient by an ideal, we need to check that the
                    # gcd of polys and the generators of the ideal is 1
                    lift_polys = [f.lift() for f in polys]
                    lift_polys = lift_polys + list(source_ring.defining_ideal().gens())
                else:
                    # if the domain is affine space, we just need to
                    # check the gcd of polys
                    lift_polys = polys
                from sage.rings.arith import gcd
                if gcd(lift_polys) != 1:
                    raise ValueError, "polys (=%s) must not have common factors"%polys
            polys = Sequence(polys)
            polys.set_immutable()
            # Todo: check that map is well defined (how?)
        self.__polys = polys
        SchemeMorphism.__init__(self, parent)

    def defining_polynomials(self):
        """
        Return the immutable sequence of polynomials that defines this
        scheme morphism.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: H([x^3+y, 1-x-y]).defining_polynomials()
            [x^3 + y, -x - y + 1]
        """
        return self.__polys

    def __call__(self, x):
        """
        Apply this morphism to a point in the domain.

        EXAMPLES::

            sage: R.<x,y> = QQ[]
            sage: A.<x,y> = AffineSpace(R)
            sage: H = A.Hom(A)
            sage: f = H([y,x^2+y])
            sage: f([2,3])
            (3, 7)

        We illustrate type checking of the input::

            sage: f(0)
            Traceback (most recent call last):
            ...
            TypeError: Argument v (=(0,)) must have 2 coordinates.
        """
        dom = self.domain()
        x = dom(x)
        P = [f(x._coords) for f in self.defining_polynomials()]
        return self.codomain()(P)


    def _repr_defn(self):
        """
        This function is used internally for printing.

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
        o = self.codomain().ambient_space()._repr_generic_point(self.defining_polynomials())
        return "Defined on coordinates by sending %s to\n%s"%(i,o)


class SchemeMorphism_on_points_affine_space(SchemeMorphism_on_points):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient affine space.

    EXAMPLES::

        sage: RA.<x,y> = QQ[]
        sage: A2 = AffineSpace(RA)
        sage: RP.<u,v,w> = QQ[]
        sage: P2 = ProjectiveSpace(RP)
        sage: H = A2.Hom(P2)
        sage: f = H([x, y, 1])
        sage: f
        Scheme morphism:
          From: Affine Space of dimension 2 over Rational Field
          To:   Projective Space of dimension 2 over Rational Field
          Defn: Defined on coordinates by sending (x, y) to
                (x : y : 1)
    """

class SchemeMorphism_on_points_projective_space(SchemeMorphism_on_points):
    """
    A morphism of schemes determined by rational functions that define
    what the morphism does on points in the ambient projective space.

    EXAMPLES::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: H([y,2*x])
        Scheme endomorphism of Projective Space of dimension 1 over Rational Field
          Defn: Defined on coordinates by sending (x : y) to
                (y : 2*x)

    An example of a morphism between projective plane curves (see #10297)::

        sage: P2.<x,y,z> = ProjectiveSpace(QQ,2)
        sage: f = x^3+y^3+60*z^3
        sage: g = y^2*z-( x^3 - 6400*z^3/3)
        sage: C = Curve(f)
        sage: E = Curve(g)
        sage: xbar,ybar,zbar = C.coordinate_ring().gens()
        sage: H = C.Hom(E)
        sage: H([zbar,xbar-ybar,-(xbar+ybar)/80])
        Scheme morphism:
          From: Projective Curve over Rational Field defined by x^3 + y^3 + 60*z^3
          To:   Projective Curve over Rational Field defined by -x^3 + y^2*z + 6400/3*z^3
          Defn: Defined on coordinates by sending (x : y : z) to
                (zbar : xbar - ybar : -1/80*xbar - 1/80*ybar)

    We illustrate some error checking::

        sage: R.<x,y> = QQ[]
        sage: P1 = ProjectiveSpace(R)
        sage: H = P1.Hom(P1)
        sage: f = H([x-y, x*y])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - y, x*y]) must be of the same degree

        sage: H([x-1, x*y+x])
        Traceback (most recent call last):
        ...
        ValueError: polys (=[x - 1, x*y + x]) must be homogeneous

        sage: H([exp(x),exp(y)])
        Traceback (most recent call last):
        ...
        TypeError: polys (=[e^x, e^y]) must be elements of Multivariate Polynomial Ring in x, y over Rational Field
    """

    def __init__(self, parent, polys, check=True):
        SchemeMorphism_on_points.__init__(self, parent, polys, check)
        if check:
            # morphisms from projective space are always given by
            # homogeneous polynomials of the same degree
            polys = self.defining_polynomials()
            try:
                d = polys[0].degree()
            except AttributeError:
                polys = [f.lift() for f in polys]
            if not all([f.is_homogeneous() for f in polys]):
                raise  ValueError, "polys (=%s) must be homogeneous"%polys
            degs = [f.degree() for f in polys]
            if not all([d==degs[0] for d in degs[1:]]):
                raise ValueError, "polys (=%s) must be of the same degree"%polys


############################################################################
# Rational points on schemes, which we view as morphisms determined
# by coordinates.
############################################################################

class SchemeMorphism_coordinates(SchemeMorphism):
    def _repr_(self):
        return self.codomain().ambient_space()._repr_generic_point(self._coords)

    def _latex_(self):
        return self.codomain().ambient_space()._latex_generic_point(self._coords)

    def __getitem__(self, n):
        return self._coords[n]

    def __iter__(self):
        return iter(self._coords)

    def __tuple__(self):
        return self._coords

    def __cmp__(self, other):
        if not isinstance(other, SchemeMorphism_coordinates):
            try:
                other = self.codomain().ambient_space()(other)
            except TypeError:
                return -1
        return cmp(self._coords, other._coords)

    def scheme(self):
        return self.codomain()


#*******************************************************************
# Affine varieties
#*******************************************************************
class SchemeMorphism_affine_coordinates(SchemeMorphism_coordinates):
    """
    A morphism determined by giving coordinates in a ring.

    INPUT:


    -  ``X`` - a subscheme of an ambient affine space over
       a ring R.

    -  ``v`` - a list or tuple of coordinates in R


    EXAMPLES::

        sage: A = AffineSpace(2, QQ)
        sage: A(1,2)
        (1, 2)
    """
    def __init__(self, X, v, check=True):
        if scheme.is_Scheme(X):
            X = X(X.base_ring())
        SchemeMorphism.__init__(self, X)
        if check:
            # Verify that there are the right number of coords
            d = X.codomain().ambient_space().ngens()
            if len(v) != d:
                raise TypeError, \
                      "Argument v (=%s) must have %s coordinates."%(v, d)
            if is_SchemeMorphism(v):
                v = list(v)
            if not isinstance(v,(list,tuple)):
                raise TypeError, \
                      "Argument v (= %s) must be a scheme point, list, or tuple."%str(v)
            # Make sure the coordinates all lie in the appropriate ring
            v = Sequence(v, X.value_ring())
            # Verify that the point satisfies the equations of X.
            X.codomain()._check_satisfies_equations(v)
        self._coords = v


#*******************************************************************
# Projective varieties
#*******************************************************************
class SchemeMorphism_projective_coordinates_ring(SchemeMorphism_coordinates):
    """
    A morphism determined by giving coordinates in a ring (how?).
    """
    def __init__(self, X, v, check=True):
        raise NotImplementedError


class SchemeMorphism_projective_coordinates_field(SchemeMorphism_projective_coordinates_ring):
    """
    A morphism determined by giving coordinates in a field.

    INPUT:


    -  ``X`` - a subscheme of an ambient projective space
       over a field K

    -  ``v`` - a list or tuple of coordinates in K


    EXAMPLES::

        sage: P = ProjectiveSpace(3, RR)
        sage: P(2,3,4,5)
        (0.400000000000000 : 0.600000000000000 : 0.800000000000000 : 1.00000000000000)

    ::

        sage: P = ProjectiveSpace(3, QQ)
        sage: P(0,0,0,0)
        Traceback (most recent call last):
        ...
        ValueError: [0, 0, 0, 0] does not define a valid point since all entries are 0
    """
    def __init__(self, X, v, check=True):
        if scheme.is_Scheme(X):
            X = X(X.base_ring())
        SchemeMorphism.__init__(self, X)
        if check:
            from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field # TODO: fix circular ref.
            d = X.codomain().ambient_space().ngens()
            if is_SchemeMorphism(v) or isinstance(v, EllipticCurvePoint_field):
                v = list(v)
            elif v is infinity:
                v = [0] * (d)
                v[1] = 1
            if not isinstance(v,(list,tuple)):
                raise TypeError, \
                      "Argument v (= %s) must be a scheme point, list, or tuple."%str(v)
            if len(v) != d and len(v) != d-1:
                raise TypeError, "v (=%s) must have %s components"%(v, d)
            #v = Sequence(v, X.base_ring())
            v = Sequence(v, X.value_ring())
            if len(v) == d-1:     # very common special case
                v.append(1)

            n = len(v)
            all_zero = True
            for i in range(n):
                if v[n-1-i]:
                    all_zero = False
                    c = v[n-1-i]
                    if c == 1:
                        break
                    for j in range(n-i):
                        v[j] /= c
                    break
            if all_zero:
                raise ValueError, "%s does not define a valid point since all entries are 0"%repr(v)

            X.codomain()._check_satisfies_equations(v)

        self._coords = v



#*******************************************************************
# Abelian varieties
#*******************************************************************
class SchemeMorphism_abelian_variety_coordinates_field(AdditiveGroupElement, SchemeMorphism_projective_coordinates_field):
    def __mul__(self, n):
        if isinstance(n, (RingElement, int, long)):
            # [n]*P - multiplication by n.
            return AdditiveGroupElement._rmul_(self, Integer(n))
        else:
            # Function composition
            return SchemeMorphism_projective_coordinates_field.__mul__(self, n)



#*******************************************************************
# Toric varieties
#*******************************************************************
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
        sage: import sage.schemes.generic.morphism as MOR
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
            sage: import sage.schemes.generic.morphism as MOR
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

