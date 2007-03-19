import lazy_ring_generic
import padic_ring_generic

class pAdicLazyRingGeneric(lazy_ring_generic.LazyRingGeneric, padic_ring_generic.pAdicRingGeneric):
    def __init__(self, p, prec, print_mode, names, halt):
        lazy_ring_generic.LazyRingGeneric.__init__(self, prec, names, halt)
        padic_ring_generic.pAdicRingGeneric.__init__(self, p, prec, print_mode, names)
