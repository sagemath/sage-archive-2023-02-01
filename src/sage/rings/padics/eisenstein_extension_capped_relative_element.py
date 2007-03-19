#import sage.rings.padics.padic_ring_generic_element
import sage.rings.padics.padic_generic_element
import sage.rings.padics.padic_extension_generic_element
import sage.rings.polynomial_element
import sage.rings.padics.eisenstein_extension_generic_element
import sage.rings.infinity
from sage.rings.padics.misc import min

infinity = sage.rings.infinity.infinity
Polynomial = sage.rings.polynomial_element.Polynomial
#pAdicRingGenericElement = sage.rings.padics.padic_ring_generic_element.pAdicRingGenericElement
pAdicGenericElement = sage.rings.padics.padic_generic_element.pAdicGenericElement
pAdicExtensionGenericElement = sage.rings.padics.padic_extension_generic_element.pAdicExtensionGenericElement
EisensteinExtensionGenericElement = sage.rings.padics.eisenstein_extension_generic_element.EisensteinExtensionGenericElement

class EisensteinExtensionCappedRelativeElement(EisensteinExtensionGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        if construct:
            (self._ordp, unit) = x
            EisensteinExtensionGenericElement.__init__(self, parent, unit, construct = True)
            return
        if check:
            from sage.rings.padics.padic_extension_leaves import EisensteinExtensionCappedRelative
            if not isinstance(parent, EisensteinExtensionGeneric):
                raise TypeError, "parent must be an EisensteinExtensionGeneric"
        if not absprec is infinity and not relprec is infinity:
            raise ValueError, "can only specify one of absprec and relprec"
        if absprec is infinity:
            if relprec > parent.precision_cap():
                relprec = parent.precision_cap()

        if isinstance(x, pAdicGenericElement) and x.parent().is_field() and x.valuation() < 0:
            raise ValueError, "element has negative valuation."
        if isinstance(x, EisensteinExtensionCappedRelativeElement):
            if x.parent() is self.parent():
                (self._ordp, unit) = x._val_unit()
                if self._ordp >= absprec:
                    unit *= 0 #make sure to test for this case
                    self._ordp = absprec
                else:
                    relprec = min(relprec, absprec - self._ordp)
                    unit = unit.add_bigoh(relprec)
                EisensteinExtensionGenericElement.__init__(self, parent, unit._value, construct = True)
            elif x.parent().fraction_field() is self.parent():
                (self._ordp, unit) = x._val_unit()
                if self._ordp >= absprec:
                    unit *= 0 #make sure to test for this case
                    self._ordp = absprec
                else:
                    relprec = min(relprec, absprec - self._ordp)
                    unit = unit.add_bigoh(relprec)
                unit = parent._PQR.polynomial_ring()(unit._value)
        EisensteinExtensionGenericElement.__init__(self, parent, x)
        self._ordp = self.valuation()
        self._value = self.__rshift__(self._ordp)._value

    def _normalize(self):
        c = EisensteinExtensionGenericElement.valuation(self)
        self._value = self.__rshift__(c)._value
        self._ordp = self._ordp + c

    def _polynomial(self):
        #Note that self._value._polynomial should have capped-relative coefficients
        return self._value._polynomial * self.parent()._PQR.polynomial_ring().gen()**self._ordp

    def __lshift__(self, shift):
        shift = Integer(shift)
        if shift < 0:
            return self.__rshift__(-shift)
        return EisensteinExtensionCappedRelativeElement(self.parent(), (self.valuation() + shift, self._value._polynomial), construct = True)

    def __rshift__(self, shift):
        shift = Integer(shift)
        if shift < 0:
            return self.__lshift__(-shift)
        if self.parent().is_field() or shift < self._ordp:
            return EisensteinExtensionCappedRelativeElement(self.parent(), (self._ordp - shift, self._value._polynomial), construct = True)
        unit = EisensteinExtensionGenericElement.__rshift__(self, shift - self._ordp)
        curval = unit.valuation()
        if curval > 0:
            return EisensteinExtensionCappedRelativeElement(self.parent(), (curval, unit.__rshift__(curval)._value._polynomial), construct = True)
        else:
            return EisensteinExtensionCappedRelativeElement(self.parent(), (Integer(0), unit._value._polynomial), construct = True)

    def _add_(self, right):
        mordp = min(self._ordp, right._ordp)
        sshift = self.parent()._PQR.polynomial_ring().gen() ** (max(Integer(0), self._ordp - right._ordp))
        rshift = self.parent()._PQR.polynomial_ring().gen() ** (max(Integer(0), right._ordp - self._ordp))
        return EisensteinExtensionCappedRelativeElement(self.parent(), (mordp, self._value._polynomial * sshift + right._value._polynomial * rshift), construct = True)

    def _mul_(self, right):
        return EisensteinExtensionCappedRelativeElement(self.parent(), (self._ordp + right._ordp, self._value._polynomial * right._value._polynomial), construct = True)

    def _neg_(self):
        return EisensteinExtensionCappedRelativeElement(self.parent(), (self._ordp, -self._value._polynomial), construct = True)

    def _sub_(self, right):
        mordp = min(self._ordp, right._ordp)
        sshift = self.parent()._PQR.polynomial_ring().gen() ** (max(Integer(0), self._ordp - right._ordp))
        rshift = self.parent()._PQR.polynomial_ring().gen() ** (max(Integer(0), right._ordp - self._ordp))
        return EisensteinExtensionCappedRelativeElement(self.parent(), (mordp, self._value._polynomial * sshift - right._value._polynomial * rshift), construct = True)

    def copy(self):
        return self.parent()._element_class(self.parent(), (self._ordp, self._value._polynomial), construct = True)

    def list(self):
        if self.parent().is_field():
            self._normalize()
            return EisensteinExtensionGenericElement.list(self)
        else:
            return [self.parent()(0)]*(self._ordp) + EisensteinExtensionGenericElement.list(self)

    def _is_exact_zero(self):
        #This should change once we get new p-adic polynomials
        return self._value._polynomial == 0

    def residue(self, n):
        if self._ordp == 0:
            return EisensteinExtensionGenericElement.residue(self, n)
        else:
            raise NotImplementedError

    def precision_absolute(self):
        return self._ordp + EisensteinExtensionGenericElement.precision_absolute(self)

    def valuation(self):
        self._normalize()
        return self._ordp
