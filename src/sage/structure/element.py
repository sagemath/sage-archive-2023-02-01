"""
PYREX: sage.ext.element
"""

from sage.ext.element import *
import sage.ext.coerce

class Element_cmp_:
    """
    Class for defining comparisons between elements.
    """
    def __cmp__(self, right):
        if not isinstance(right, Element) or right.parent() != self.parent():
            return sage.ext.coerce.cmp(self, right)
        return self._cmp_(right)

    def _cmp_(self, right):
        raise NotImplementedError
