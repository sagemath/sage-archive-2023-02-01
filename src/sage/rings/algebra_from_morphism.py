from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.decorators import sage_wraps
from sage.categories.algebras import Algebras
from sage.rings.integer_ring import IntegerRing
from sage.rings.ring import CommutativeRing, CommutativeAlgebra
from sage.rings.morphism import RingHomomorphism

from sage.categories.action import Action
import operator


def coerce_algebras(method):
    @sage_wraps(method)
    def new_method(*algebras, **kwargs):
        other_factors = [ ]
        if base is None:
            base = self.base_ring()
        for i in len(algebras):
            algebra = algebras[i]
            if not isinstance(algebra, CommutativeAlgebra):
                raise TypeError("%s is not a commutative algebra" % algebra)
            if not isinstance(algebra, AlgebraFromMorphism):
                algebra = RingExtension(algebra)
            base_factor = algebra.base_ring()
            if base_factor is base:
                pass
            elif base_factor.has_coerce_map_from(base):
                morphism = base_factor.coerce_map_from(base)
                algebra = self.scalar_restriction(morphism)
            else:
                raise TypeError("No coercion map from %s to %s" % (over, base))
            if i == 0:
                first_factor = algebra
            else:
                other_factors.append(algebra)
        if first_factor is algebras[0]:
            return method(first_factor, *other_factors, **kwargs)
        else:
            return getattr(first_factor, method.__name__)(*other_factors, **kwargs)
    return new_method


# move the next two classes to element.pyx ???
class LeftBaseActionOnRing(Action):
    def __init__(self, algebra):
        if not isinstance(algebra, AlgebraFromMorphism):
            raise TypeError("%s is not an instance of AlgebraFromMorphism" % algebra)
        self._base = algebra.base_ring()
        self._algebra = algebra
        self._defining_morphism = algebra.defining_morphism()
        Action.__init__(self, self._base, algebra, op=operator.mul)

    def _call_(self, scalar, element):
        result = self._defining_morphism(scalar) * element.element_in_ring()
        return self._algebra(result)

class RightBaseActionOnRing(Action):
    def __init__(self, algebra):
        if not isinstance(algebra, AlgebraFromMorphism):
            raise TypeError("%s is not an instance of AlgebraFromMorphism" % algebra)
        self._base = algebra.base_ring()
        self._algebra = algebra
        self._defining_morphism = algebra.defining_morphism()
        Action.__init__(self, self._base, algebra, is_left=0, op=operator.mul)

    def _call_(self, element, scalar):
        result = self._defining_morphism(scalar) * element.element_in_ring()
        return self._algebra(result)


class AlgebraFromMorphism(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, defining_morphism):
        base = defining_morphism.domain()
        self._ring = ring = defining_morphism.codomain()
        self._defining_morphism = defining_morphism

        # To avoid problems...
        if isinstance(base, AlgebraFromMorphism) or isinstance(ring, AlgebraFromMorphism):
            raise NotImplementedError

        # We do not want to have a coercion map base -> self
        # So we first set the base to ring.base() (which indeed coerces to self)
        CommutativeAlgebra.__init__(self, ring.base(), category=Algebras(base)) 
        # and then update the base
        self._base = base

        # Set left and right actions of base on self
        # together with a coercion map self -> ring
        self._unset_coercions_used()
        left = LeftBaseActionOnRing(self)
        right = RightBaseActionOnRing(self)
        from sage.rings.morphism import AlgebraToRing_coercion
        coercion_morphism = AlgebraToRing_coercion(self.Hom(ring), check=False)
        self._populate_coercion_lists_(action_list=[left,right], embedding=coercion_morphism)

        from sage.structure.element import AlgebraFMElement
        self.element_class = AlgebraFMElement

    def _element_constructor_(self, x, *args, **kwargs):
        from sage.structure.element import AlgebraFMElement
        if isinstance(x, AlgebraFMElement):
            x = x.element_in_ring()
        if self._base.has_coerce_map_from(x.parent()):
            x = self._defining_morphism(self._base(x))
        elt = self._ring(x, *args, **kwargs)
        return self.element_class(self, elt)

    def _repr_(self):
        return "%s viewed as an algebra over %s" % (self._ring, self._base)

    def _coerce_map_from_(self, other):
        if isinstance(other, AlgebraFromMorphism):
            return other.base().has_coerce_map_from(self._base) and self._ring.has_coerce_map_from(other.ring())

    def defining_morphism(self):
        return self._defining_morphism

    def ring(self):
        return self._ring

    def _an_element_(self):
        elt = self._ring.an_element()
        return self.element_class(self, elt)

    def random_element(self):
        elt = self._ring.random_element()
        return self.element_class(self, elt)

    def scalar_restriction(self, morphism):
        if isinstance(morphism, AlgebraFromMorphism):
            morphism = morphism.defining_morphism()
        elif isinstance(morphism, CommutativeRing):
            if self._base.has_coerce_map_from(morphism):
                morphism = self._base.coerce_map_from(morphism)
            else:
                raise TypeError("No coercion map from %s to %s" % (morphism, self._base))
        elif morphism.codomain() is not self._base:
            raise TypeError("The codomain of %s is not %s" % (morphism, self._base))
        defining_morphism = self._defining_morphism.pre_compose(morphism)
        return AlgebraFromMorphism(defining_morphism=defining_morphism)

    @coerce_algebras
    def tensor_product(self, other, base=None):
        # Implement this method in subclasses
        raise NotImplementedError("Tensor product of extensions is not implemented")

    def scalar_extension(self, morphism):
        if isinstance(morphism, AlgebraFromMorphism):
            morphism = morphism.defining_morphism()
        elif isinstance(morphism, CommutativeRing):
            if morphism.has_coerce_map_from(self._base):
                morphism = morphism.coerce_map_from(self._base)
            else:
                raise TypeError("No coercion map from %s to %s" % (self._base, morphism))
        elif morphism.domain() is not self._base:
            raise TypeError("The domain of %s is not %s" % (morphism, self._base))
        extension = AlgebraFromMorphism(defining_morphism=morphism)
        algebra, morphisms = self.tensor_product(extension)
        return AlgebraFromMorphism(defining_morphism=morphisms[1].ring())


def RingExtension(arg1, arg2=None):
    if isinstance(arg1, RingHomomorphism):   # fix this
        return AlgebraFromMorphism(arg1)
    elif isinstance(arg1, CommutativeRing):
        ring = arg1
        if arg2 is None:
            base = ring.base_ring()
            if base is None:
                # ZZ is initial in the category of rings
                base = IntegerRing()
        elif isinstance(arg2, CommutativeRing):
            base = arg2
            if not ring.has_coerce_map_from(base):
                raise ValueError("No coercion map from %s to %s" % (base,ring))
        else:
            raise TypeError("The Extension constructor must be called either with a defining morphism or two rings")
        morphism = ring.coerce_map_from(base)
        return AlgebraFromMorphism(morphism)
    else:
        raise TypeError("The Extension constructor must be called either with a defining morphism or two rings")
