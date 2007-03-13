import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_generic_element
import sage.rings.padics.padic_extension_element
import sage.rings.polynomial_element
import sage.rings.padics.eisenstein_extension_element
import sage.rings.infinity

infinity = sage.rings.infinity.infinity
Polynomial = sage.rings.polynomial_element.Polynomial
pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pAdicExtensionElement = sage.rings.padics.padic_extension_element.pAdicExtensionElement
EisensteinExtensionElement = sage.rings.padics.eisenstein_extension_element.EisensteinExtensionElement

class EisensteinExtensionCappedRelativeElement(pAdicRingGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        raise NotImplementedError
        if construct:
            (self._ordp, self._unit) = x
            return
        if check:
            from sage.rings.padics.eisenstein_extension_generic import EisensteinExtensionGeneric
            if not isinstance(parent, EisensteinExtensionGeneric):
                raise TypeError, "parent must be an EisensteinExtensionGeneric"
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if absprec is infinity:
            if relprec > parent.precision_cap():
                relprec = parent.precision_cap()

        if isinstance(x, (pAdicGenericElement, pAdicExtensionElement)) and x.parent().is_field() and x.valuation() < 0:
            raise ValueError, "element has negative valuation."
        if isinstance(x, EisensteinRingExtensionCappedRelativeElement) and x.parent() is self.parent():
            (self._ordp, self._unit) = x._val_unit()
            if self._ordp >= absprec:
                self._unit *= 0 #make sure to test for this case
                self._ordp = absprec
            else:
                relprec = min(relprec, absprec - self._ordp)
                self._unit = self._unit.add_bigoh(relprec)
        elif isinstance(x, Polynomial):
            if x.parent() is parent.polynomial_ring():
                x = x.list()
                x = x + [self.base_ring()(0)] * (parent().degree() - len(x))
        if isinstance(x, list):
            if len(x) == parent.degree():
                #should probably be wrapped in a try block, but I'm not sure what types of exceptions I'm worried about.  TypeError?
                absprec = min([v.precision_absolute() for v in x].append(absprec))
                val = min([v.valuation() for v in x])
                absprec = min(absprec, relprec + val)
                x = [parent.base_ring()(v, absprec = absprec).__rshift__(val) for v in x]
                self._unit = EisensteinExtensionElement(parent, x) # the parent here makes my head hurt a bit.
                self._ordp = val
                return
        #We now try casting into the base ring
        x = parent.base_ring()(x, absprec = absprec, relprec = relprec)
        self._ordp = x.valuation()
        self._unit = EisensteinExtensionElement(parent, parent.polynomial_ring()(x.__rshift__(self._ordp)), construct = True)

    def _neg_(self):
        return EisenstenRingExtensionCappedRelativeElement(self.parent(), (self._ordp, -self._unit), construct = True)

    def __pow__(self, right):
        raise NotImplementedError

    def _add_(self, right):
        relprec = min
        #Design flaw, giving up now.
        return EisensteinRingExtensionCappedRelativeElement(self.parent(), (self._ordp * right, -self._unit + right._unit), construct = True)
