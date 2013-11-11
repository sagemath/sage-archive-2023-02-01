r"""
Six Vertex Model
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.combinat.combinat import CombinatorialObject

class SixVertexConfiguration(CombinatorialObject, Element):
    """
    A configuration in the six vertex model.
    """
    def __init__(self, parent, lst):
        """
        Initialize ``self``.
        """
        Element.__init__(self, parent)
        CombinatorialObject.__init__(self, lst)

    def to_digraph(self):
        """
        Return a digraph version of ``self``.
        """

    def plot(self):
        """
        Return a plot of ``self``.
        """

class SixVertexModel(Parent, UniqueRepresentation):
    """
    The six vertex model.

    We model a configuration by indicating which configuration by the
    following six configurations which are determined by the two outgoing
    arrows:

    1 - LR
    2 - LU
    3 - LD
    4 - UD
    5 - UR
    6 - RD

    INPUT:

    - ``n`` -- the number of rows
    - ``m`` -- (optional) the number of columns, if not specified, then
      the number of columns is the number of rows
    """
    @staticmethod
    def __classcall_private__(cls, n, m=None, boundary_conditions=None):
        """
        Normalize input to ensure a unique representation.
        """
        if m is None:
            m = n
        if boundary_conditions is not None:
            boundary_conditions = tuple(tuple(x) for x in boundary_conditions)
        else:
            boundary_conditions = (None, None, None, None)
        return super(SixVertexModel, cls).__classcall__(cls, n, m, boundary_conditions)

    def __init__(self, n, m, boundary_conditions):
        """
        Initialize ``self``.
        """
        self._nrows = n
        self._ncols = m
        self._boundary_conditions = boundary_conditions
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def __iter__(self):
        """
        Iterate through ``self``.
        """
        # TODO

    Element = SixVertexConfiguration

