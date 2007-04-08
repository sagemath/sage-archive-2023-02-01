import sage.rings.padics.local_generic

class CappedAbsoluteGeneric(sage.rings.padics.local_generic.LocalGeneric):
    def is_capped_absolute(self):
        return True
