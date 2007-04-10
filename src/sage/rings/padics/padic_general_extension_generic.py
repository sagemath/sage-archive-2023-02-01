import padic_extension_generic

class pAdicGeneralExtensionGeneric(padic_extension_generic.pAdicExtensionGeneric):
    def __init__(self, upoly, epoly, poly, prec, print_mode, names, element_class):
        #This code will change in the future.
        self._PQR = PolynomialQuotientRing(epoly.base_ring(), epoly, names = names)
        pAdicExtensionGeneric.__init__(self, poly, prec, print_mode, names, element_class)
