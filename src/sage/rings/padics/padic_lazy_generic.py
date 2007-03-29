from sage.rings.padics.padic_generic import pAdicGeneric

class pAdicLazyGeneric(pAdicGeneric):
    r"""
    Generic Lazy p-adic ring/field
    """
    def __init__(self, p, prec, print_mode, names, halt):
        pAdicGeneric.__init__(self, p, prec, print_mode, names)
        self._halt = halt

    def halting_parameter(self):
        return self._halt
