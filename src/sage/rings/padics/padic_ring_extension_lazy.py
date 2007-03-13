import sage.rings.padics.padic_ring_extension_generic
import sage.rings.padics.padic_lazy_generic

pAdicRingExtensionGeneric = sage.rings.padics.padic_ring_extension_generic.pAdicRingExtensionGeneric
pAdicLazyGeneric = sage.rings.padics.padic_lazy_generic.pAdicLazyGeneric

class pAdicRingExtensionLazy(pAdicRingExtensionGeneric, pAdicLazyGeneric):
    r"""
    General Extension of a lazy p-adic ring.

    You should not create this class directly unless you know what you're doing.  Use ExtensionFactory.
    """

    def __init__(self, base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode):
        raise NotImplementedError
