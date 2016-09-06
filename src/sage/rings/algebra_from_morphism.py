from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.decorators import sage_wraps
from sage.categories.pushout import pushout
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


# move the next class to element.pyx ???
class BaseActionOnRing(Action):
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


# Constructor for AlgebraFromMorphism
def RingExtension(ring, base=None, defining_morphism=None):
    if defining_morphism is None:
        coerce = True
        if base is None:
            base = ring.base_ring()
            if base is None:
                # ZZ is initial in the category of rings
                base = IntegerRing()
        if ring.has_coerce_map_from(base):
            defining_morphism = ring.coerce_map_from(base)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,ring))
    else:
        domain = defining_morphism.domain()
        codomain = defining_morphism.codomain()
        coerce = (defining_morphism == codomain.coerce_map_from(domain))
        if base is None:
            base = domain
        elif base is domain:
            pass
        elif domain.has_coerce_map_from(base):
            coercion = domain.coerce_map_from(base)
            defining_morphism = defining_morphism.pre_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (base,domain))
        if ring is codomain:
            pass
        elif ring.has_coerce_map_from(codomain):
            coercion = ring.coerce_map_from(codomain)
            defining_morphism = defining_morphism.post_compose(coercion)
        else:
            raise ValueError("No coercion map from %s to %s" % (codomain,ring))
    return AlgebraFromMorphism(defining_morphism, coerce)


class AlgebraFromMorphism(CommutativeAlgebra, UniqueRepresentation):
    def __init__(self, defining_morphism, coerce):
        base = defining_morphism.domain()
        ring = defining_morphism.codomain()

        # To avoid problems...
        if isinstance(base, AlgebraFromMorphism) or isinstance(ring, AlgebraFromMorphism):
            raise NotImplementedError

        # We do not want to have a coercion map base -> self
        # So we first set the base to ring.base() (which indeed coerces to self)
        CommutativeAlgebra.__init__(self, ring.base(), category=Algebras(base)) 
        # and then update the base
        self._base = base
        self._ring = ring
        self._defining_morphism = defining_morphism
        self._coerce = coerce

        # Set left action of base on self
        # and a coercion map self -> ring
        from sage.rings.morphism import AlgebraToRing_coercion
        self._unset_coercions_used()
        self._populate_coercion_lists_(
            action_list = [BaseActionOnRing(self)],
            embedding = AlgebraToRing_coercion(self.Hom(ring), check=False))

        from sage.structure.element import AlgebraFMElement
        self.element_class = AlgebraFMElement

    def _element_constructor_(self, x, *args, **kwargs):
        from sage.structure.element import AlgebraFMElement
        if isinstance(x, AlgebraFMElement):
            x = x.element_in_ring()
        elt = self._ring(x, *args, **kwargs)
        return self.element_class(self, elt)

    def from_base_ring(self, x):
        elt = self._ring(x)
        return self.element_class(self, elt)

    def _pushout_(self, other):
        base = None
        if isinstance(other,AlgebraFromMorphism):
            # TODO: implement pushout when:
            # - defining_morphism are not coercion maps
            #   (this implies to check the equality of two morphisms)
            # - there is no direct coercion map between the bases
            #   but there exists a third ring which coerces to eash base
            #   Question: how can one find the greatest such third ring?
            if self._coerce and other._coerce:
                sbase = self._base
                obase = other._base
                if sbase.has_coerce_map_from(obase):
                    base = obase
                elif obase.has_coerce_map_from(sbase):
                    base = sbase
            ring = pushout(self._ring, other.ring())
        else:
            ring = pushout(self._ring, other)
        if base is None:
            return ring
        else:
            return RingExtension(ring,base)

    def _repr_(self):
        return "%s viewed as an algebra over %s" % (self._ring, self._base)

    def _coerce_map_from_(self, other):
        if isinstance(other, AlgebraFromMorphism):
            if self._coerce and other._coerce:
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
