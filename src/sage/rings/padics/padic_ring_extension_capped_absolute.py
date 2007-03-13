import sage.rings.padics.padic_ring_extension_generic

pAdicRingExtensionGeneric = sage.rings.padics.padic_ring_extension_generic.pAdicRingExtensionGeneric

class pAdicRingExtensionCappedAbsolute(pAdicRingExtensionGeneric):
    r"""
    General Extension of a p-adic ring with capped absolute precision

    You should not create this class directly unless you know what you're doing.  Use ExtensionFactory.
    """

    def __init__(self, base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode):
        raise NotImplementedError
