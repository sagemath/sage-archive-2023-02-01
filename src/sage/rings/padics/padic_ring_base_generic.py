import padic_base_generic
import padic_ring_generic

class pAdicRingBaseGeneric(padic_base_generic.pAdicBaseGeneric, padic_ring_generic.pAdicRingGeneric):
    def unit_group(self):
        raise NotImplementedError

    def unit_group_gens(self):
        raise NotImplementedError

    def principal_unit_group(self):
        raise NotImplementedError
