import padic_extension_lazy_generic_element
from sage.rings.infinity import infinity

pAdicExtensionLazyGenericElement = padic_extension_lazy_generic_element.pAdicExtensionLazyGenericElement

class EisensteinExtensionLazyElement(pAdicExtensionLazyGenericElement):
    def __init__(self, parent, x, absprec = infinity, relprec = infinity, check = True, construct = False):
        raise NotImplementedError
