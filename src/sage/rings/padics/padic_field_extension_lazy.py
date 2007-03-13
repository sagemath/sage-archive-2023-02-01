import sage.rings.padics.padic_field_extension_generic
import sage.rings.padics.padic_lazy_generic

pAdicFieldExtensionGeneric = sage.rings.padics.padic_field_extension_generic.pAdicFieldExtensionGeneric
pAdicLazyGeneric = sage.rings.padics.padic_lazy_generic.pAdicLazyGeneric

class pAdicFieldExtensionLazy(pAdicFieldExtensionGeneric, pAdicLazyGeneric):
    r"""
    General Extension of a lazy p-adic field.

    You should not create this class directly unless you know what you're doing.  Use ExtensionFactory.
    """

    def __init__(self, base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode):
        raise NotImplementedError
