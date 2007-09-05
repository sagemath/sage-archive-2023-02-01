"""
The class group of a number field.
"""

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class

class ClassGroup(AbelianGroup_class):
    """
    The class group of a number field.

    WARNING/TODO: The map from generators to actual ideals
    of the number field isn't yet implemented.

    """
    def __init__(self, invariants, names, number_field):
        self.__number_field = number_field
        AbelianGroup_class.__init__(self, len(invariants), invariants, names)

    def _repr_(self):
        return '%s as the class group of %s'%(
            AbelianGroup_class._repr_(self),
            self.number_field())

    def number_field(self):
        return self.__number_field

