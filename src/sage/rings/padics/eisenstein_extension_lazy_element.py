import sage.rings.padics.padic_lazy_generic_element

pAdicLazyGenericElement = sage.rings.padics.padic_lazy_generic_element

class EisensteinExtensionLazyElement(pAdicLazyGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        raise NotImplementedError
