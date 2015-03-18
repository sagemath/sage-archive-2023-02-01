r"""
Fully packed loops
"""
from sage.structure.sage_object import SageObject
from sage.combinat.six_vertex_model import SquareIceModel, SixVertexConfiguration
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix

class FullyPackedLoop(SageObject):
    """
    A class for fully packed loops
    """

    def __init__(self, generator):
        """
        Initialise object: what are we going to use as generators? Perhaps multiple: ASM, and Six Vertex Model

        EXAMPLES:

        We can initiate a fully packed loop using an Alternating Sign Matrix::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.six_vertex_model
            Stuff

        Otherwise we initiate a fully packed loop using a six vertex model::

            sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl.six_vertex_model
                ^    ^    ^
                |    |    |
            --> # -> # <- # <--
                ^    |    ^
                |    V    |
            --> # <- # -> # <--
                |    ^    |
                V    |    V
            --> # -> # <- # <--
                |    |    |
                V    V    V
        """
        if isinstance(generator, AlternatingSignMatrix):
            generator = generator.to_six_vertex_model()
        self.six_vertex_model = generator

    def _repr_():
        """
        Return the ascii key
        """

    def to_signed_matrix():
        """
        Just call from six vertex model
        """

    def plot():
        """
        Tweak the six vertex model stuff
        """
