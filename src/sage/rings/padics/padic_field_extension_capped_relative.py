import sage.rings.padics.padic_field_extension_generic

pAdicFieldExtensionGeneric = sage.rings.padics.padic_field_extension_generic.pAdicFieldExtensionGeneric

class pAdicFieldExtensionCappedRelative(pAdicFieldExtensionGeneric):
    r"""
    General Extension of a p-adic field with capped relative precision

    You should not create this class directly unless you know what you're doing.  Use ExtensionFactory.
    """

    def __init__(self, base, upoly, intermediate_ext, epoly, names, prec, halt, print_mode):
        raise NotImplementedError
