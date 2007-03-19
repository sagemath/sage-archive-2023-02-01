import lazy_field_generic
import padic_field_generic

class pAdicLazyFieldGeneric(lazy_field_generic.LazyFieldGeneric, padic_field_generic.pAdicFieldGeneric):
    def __init__(self, p, prec, print_mode, names, halt):
        lazy_field_generic.LazyFieldGeneric.__init__(self, prec, names, halt)
        padic_field_generic.pAdicFieldGeneric.__init__(self, p, prec, print_mode, names)
