#*****************************************************************************
#      Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.element cimport Element
from sage.categories.morphism cimport Morphism
from sage.categories.map cimport Map

from sage.rings.real_mpfr import RealField_constructor as RealField, mpfr_prec_min
from sage.rings.complex_field import ComplexField
from sage.rings.real_lazy import RLF, CLF


cdef class NumberFieldEmbedding(Morphism):

    cdef _gen_image

    def __init__(self, K, R, gen_embedding):
        """
        If R is a lazy field, the closest root to gen_embedding will be chosen.

        EXAMPLES:
            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^3-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1)
            sage: f(a)^3
            2.00000000000000?
            sage: RealField(200)(f(a)^3)
            2.0000000000000000000000000000000000000000000000000000000000

            sage: sigma_a = K.polynomial().change_ring(CC).roots()[1][0]; sigma_a
            -0.629960524947436 + 1.09112363597172*I
            sage: g = NumberFieldEmbedding(K, CC, sigma_a)
            sage: g(a+1)
            0.370039475052564 + 1.09112363597172*I
        """
        from sage.categories.homset import Hom
        from sage.rings.real_lazy import LazyField, LazyAlgebraic
        Morphism.__init__(self, Hom(K, R))
        if isinstance(R, LazyField) and not isinstance(gen_embedding.parent(), LazyField):
            self._gen_image = LazyAlgebraic(R, K.polynomial(), gen_embedding, prec=0)
        else:
            self._gen_image = R(gen_embedding)

    cpdef Element _call_(self, x):
        """
        EXAMPLES:
            sage: x = polygen(QQ)
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^2-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1.4)
            sage: f(a)
            1.414213562373095?
        """
        return x.polynomial()(self._gen_image)

    def _repr_defn(self):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.number_field_morphisms import NumberFieldEmbedding
            sage: K.<a> = NumberField(x^2-2)
            sage: f = NumberFieldEmbedding(K, RLF, 1.4)
            sage: f
            Generic morphism:
              From: Number Field in a with defining polynomial x^2 - 2
              To:   Real Lazy Field
              Defn: a -> 1.414213562373095?
        """
        return "%s -> %s" % (self._domain.variable_name(), self._gen_image)

    def gen_image(self):
        """
        Returns the image of the generator under this embedding.

        EXAMPLES:
            sage: f = QuadraticField(7, 'a', embedding=2).coerce_embedding()
            sage: f.gen_image()
            2.645751311064591?
        """
        return self._gen_image


cdef class EmbeddedNumberFieldMorphism(NumberFieldEmbedding):
    r"""
    This allows one to go from one number field in another consistantly,
    assuming they both have specified embeddings into an ambient field
    (by default it looks for an embedding into $\C$).
    """
    cdef readonly ambient_field

    def __init__(self, K, L, ambient_field=None):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldMorphism
            sage: K.<a> = NumberField(x^2-17, embedding=4.1)
            sage: L.<b> = NumberField(x^4-17, embedding=2.0)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(a)
            b^2

            sage: K.<zeta12> = CyclotomicField(12)
            sage: L.<zeta36> = CyclotomicField(36)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(zeta12)
            zeta36^3
            sage: f(zeta12^5-zeta12+1)
            zeta36^9 - 2*zeta36^3 + 1
            sage: f
            Generic morphism:
              From: Cyclotomic Field of order 12 and degree 4
              To:   Cyclotomic Field of order 36 and degree 12
              Defn: zeta12 -> zeta36^3

        The embeddings must be compatible:
            sage: F1 = NumberField(x^3 + 2, 'a', embedding=2)
            sage: F2 = NumberField(x^3 + 2, 'a', embedding=CC.0)
            sage: F1.gen() + F2.gen()
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand parent(s) for '+': 'Number Field in a with defining polynomial x^3 + 2' and 'Number Field in a with defining polynomial x^3 + 2'
        """
        if ambient_field is None:
            from sage.rings.complex_double import CDF
            ambient_field = CDF
        gen_image = matching_root(K.polynomial().change_ring(L), K.gen(), ambient_field=ambient_field, margin=2)
        if gen_image is None:
            raise ValueError, "No consistant embedding of all of %s into %s." % (K, L)
        NumberFieldEmbedding.__init__(self, K, L, gen_image)
        self.ambient_field = ambient_field

    def section(self):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldMorphism
            sage: K.<a> = NumberField(x^2-700, embedding=25)
            sage: L.<b> = NumberField(x^6-700, embedding=3)
            sage: f = EmbeddedNumberFieldMorphism(K, L)
            sage: f(2*a-1)
            2*b^3 - 1
            sage: g = f.section()
            sage: g(2*b^3-1)
            2*a - 1
        """
        return EmbeddedNumberFieldConversion(self._codomain, self._domain, self.ambient_field)


cdef class EmbeddedNumberFieldConversion(Map):
    r"""
    This allows one to cast one number field in another consistantly,
    assuming they both have specified embeddings into an ambient field
    (by default it looks for an embedding into $\C$).

    This is done by factoring the minimal polynomial of the input
    in the number field of the codomain. This may fail if the element is
    not actually in the given field.
    """
    cdef _gen_image
    cdef readonly ambient_field

    def __init__(self, K, L, ambient_field=None):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldConversion
            sage: K.<a> = NumberField(x^2-17, embedding=4.1)
            sage: L.<b> = NumberField(x^4-17, embedding=2.0)
            sage: f = EmbeddedNumberFieldConversion(L, K)
            sage: f(b^2)
            a
            sage: f(L(a/2-11))
            1/2*a - 11
        """
        if ambient_field is None:
            from sage.rings.complex_double import CDF
            ambient_field = CDF
        self.ambient_field = ambient_field
        Map.__init__(self, K, L)

    cpdef Element _call_(self, x):
        """
        EXAMPLES:
            sage: from sage.rings.number_field.number_field_morphisms import EmbeddedNumberFieldConversion
            sage: K.<zeta12> = CyclotomicField(12)
            sage: L.<zeta15> = CyclotomicField(15)
            sage: f = EmbeddedNumberFieldConversion(K, L)
            sage: f(zeta12^4)
            zeta15^5
        """
        minpoly = x.minpoly()
        gen_image = matching_root(minpoly.change_ring(self._codomain), x, self.ambient_field, 4)
        if gen_image is None:
            raise ValueError, "No consistant embedding of %s into %s." % (self._domain, self._codomain)
        return gen_image


cpdef matching_root(poly, target, ambient_field=None, margin=1, max_prec=None):
    """
    Given a polynomial and a target, this function chooses the root that
    target best approximates as compared in ambient_field.

    If the parent of target is exact, the equality is required, otherwise
    find closest root (with respect to the \code{abs} function) in the
    ambient field to the target, and return the root of poly (if any) that
    approximates it best.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_morphisms import matching_root
        sage: R.<x> = CC[]
        sage: matching_root(x^2-2, 1.5)
        1.41421356237310
        sage: matching_root(x^2-2, -100.0)
        -1.41421356237310
        sage: matching_root(x^2-2, .00000001)
        1.41421356237310
        sage: matching_root(x^3-1, CDF.0)
        -0.500000000000000 + 0.86602540378443...*I
        sage: matching_root(x^3-x, 2, ambient_field=RR)
        1.00000000000000
    """
    if isinstance(poly, list):
        roots = poly
    else:
        roots = poly.roots()
        if len(roots) == 0:
            return None
        elif isinstance(roots[0], tuple): # as returned from the roots method
            roots = [r for r, e in roots]

    if ambient_field is None:
        ambient_field = target.parent()

    if ambient_field.is_exact():
        target_approx = ambient_field(target)
        for r in roots:
            if ambient_field(r) == target_approx:
                return r

    else:
        # since things are inexact, try and pick the closest one
        if max_prec is None:
            max_prec = ambient_field.prec() * 32
        while ambient_field.prec() < max_prec:
            if isinstance(poly, list):
                ambient_roots = [ambient_field(r) for r in poly]
            else:
                ambient_roots = [r for r, e in poly.change_ring(ambient_field).roots()]
            target_root = closest(ambient_field(target), ambient_roots, margin)
            if target_root is not None:
                for r in roots:
                    if closest(ambient_field(r), ambient_roots, margin) is target_root:
                        return r
            ambient_field = ambient_field.to_prec(ambient_field.prec() * 2)


cpdef closest(target, values, margin=1):
    """
    This is a utility function that returns the item in values closest to
    target (with respect to the \code{abs} function). If margin is greater
    than 1, and x and y are the first and second closest elements to target,
    then only return x if x is margin times closer to target than y, i.e.
    margin * abs(target-x) < abs(target-y).

    TESTS:
        sage: from sage.rings.number_field.number_field_morphisms import closest
        sage: closest(1.2, [0,1,2,3,4])
        1
        sage: closest(1.7, [0,1,2,3,4])
        2
        sage: closest(1.7, [0,1,2,3,4], margin=5)
        sage: closest(1.9, [0,1,2,3,4], margin=5)
        2
        sage: closest(.2, [-1, 1, CDF.0, -CDF.0])
        1
    """
    cdef int i
    if len(values) == 0:
        raise ValueError
    elif len(values) == 1:
        return values[0]
    else:
        dists = [abs(target - r) for r in values]
        sdists = sorted(dists)
        min_dist = sdists[0]
        if margin*min_dist < sdists[1]:
            for i in range(len(values)):
                if dists[i] is min_dist:
                    return values[i]
        else:
            return None

def create_embedding_from_approx(K, gen_image):
    """
    This creates a morphsim into from K into the parent
    of gen_image, choosing as the image of the generator
    the closest root to gen_image in its parent.

    If gen_image is in a real or complex field, then
    it creates an image into a lazy field.

    EXAMPLES:
        sage: from sage.rings.number_field.number_field_morphisms import create_embedding_from_approx
        sage: K.<a> = NumberField(x^3-x+1/10)
        sage: create_embedding_from_approx(K, 1)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> 0.9456492739235915?
        sage: create_embedding_from_approx(K, 0)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> 0.10103125788101081?
        sage: create_embedding_from_approx(K, -1)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Real Lazy Field
          Defn: a -> -1.046680531804603?

    We can define embeddings from one number field to another:
        sage: L.<b> = NumberField(x^6-x^2+1/10)
        sage: create_embedding_from_approx(K, b^2)
        Generic morphism:
          From: Number Field in a with defining polynomial x^3 - x + 1/10
          To:   Number Field in b with defining polynomial x^6 - x^2 + 1/10
          Defn: a -> b^2

    The if the embedding is exact, it must be valid:
        sage: create_embedding_from_approx(K, b)
        Traceback (most recent call last):
        ...
        ValueError: b is not a root of the defining polynomial of Number Field in a with defining polynomial x^3 - x + 1/10
    """
    if gen_image is None:
        return None
    elif isinstance(gen_image, Map):
        return gen_image
    elif isinstance(gen_image, Element):
        f = K.defining_polynomial()
        P = gen_image.parent()
        if not P.is_exact() or f(gen_image) != 0:
            RR = RealField(mpfr_prec_min())
            CC = ComplexField(mpfr_prec_min())
            if RR.has_coerce_map_from(P):
                P = RLF
            elif CC.has_coerce_map_from(P):
                P = CLF
            # padic lazy, when implemented, would go here
            elif f(gen_image) != 0:
                raise ValueError, "%s is not a root of the defining polynomial of %s" % (gen_image, K)
        return NumberFieldEmbedding(K, P, gen_image)
    else:
        raise TypeError, "Embedding (type %s) must be a morphism or element." % type(gen_image)
