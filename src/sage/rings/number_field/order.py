from sage.rings.ring import DedekindDomain


class Order(DedekindDomain):

    def __init__(self, K):
        self._ambient_global_field = K

    def number_field(self):
        return self._K

    def fraction_field(self):
        return self._K


class AbsoluteOrder(Order):
    def __init__(self, K, basis):
        Order.__init__(K)
        self._basis = basis


class RelativeOrder(Order):

    def __init__(self, K, base, embedding):
        Order.__init__(K)
        self._base = base
        self._embedding = embedding
