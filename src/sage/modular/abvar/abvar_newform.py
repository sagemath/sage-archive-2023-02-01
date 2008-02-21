class ModularAbelianVariety_newform(ModularAbelianVariety):
    """
    A modular abelian variety attached to a specific newform.
    """
    def __init__(self, f):
        """
        Create the modular abelian variety $A_f$ attached to the
        newform $f$.

        INPUT:
            f -- a newform
        """
        self.__f = f

