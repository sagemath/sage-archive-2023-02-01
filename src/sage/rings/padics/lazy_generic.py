import local_generic

class LazyGeneric(local_generic.LocalGeneric):
    def __init__(self, prec, names, halt):
        self._halt = halt
        local_generic.LocalGeneric.__init__(self, prec, names)

    def halting_parameter(self):
        return self._halt

    def is_lazy(self):
        return True
