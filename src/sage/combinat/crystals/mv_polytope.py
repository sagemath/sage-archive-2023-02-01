r"""
A module for implementing ``MVPolytopes``.
"""

from sage.structure.parent import Parent
from sage.structure.element import Element


class MVPolytope(Element):
    def __init__(self, parent, long_word, lusztig_datum):
        Element.__init__(self, parent)
        self._initial_long_word = tuple(long_word)
        self._lusztig_datum = tuple(lusztig_datum)
        self._lusztig_data_dict = {self._initial_long_word: self._lusztig_datum}

class MVPolytopes(Parent):
    def __init__(self, root_system):
        Parent.__init__(self)
        if not root_system.is_finite():
            raise ValueError("{} is not a finite root system".format(root_system))
        self.root_system = root_system

    def _element_constructor_(self, *args):
        return self.element_class(self, *args)

