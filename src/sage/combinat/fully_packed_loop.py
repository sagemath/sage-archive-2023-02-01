r"""
Fully packed loops
"""
from sage.structure.sage_object import SageObject
from sage.combinat.six_vertex_model import SquareIceModel, SixVertexConfiguration
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
from sage.plot.graphics import Graphics
from sage.plot.line import line

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
            sage: fpl.six_vertex_model.to_alternating_sign_matrix()
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
        """
        if isinstance(generator, AlternatingSignMatrix):
            generator = generator.to_six_vertex_model()
        self.six_vertex_model = generator

    def _repr_():
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: M = SixVertexModel(3, boundary_conditions='ice')
            sage: M[0]
                ^    ^    ^
                |    |    |
            --> # <- # <- # <--
                |    ^    ^
                V    |    |  
            --> # -> # <- # <--
                |    |    ^
                V    V    |
            --> # -> # -> # <--
                |    |    |
                V    V    V
        """
        # List are in the order of URDL
        ascii = [[r'  V  ', ' -', r'  ^  ', '- '], # LR
                 [r'  |  ', ' <', r'  ^  ', '- '], # LU
                 [r'  V  ', ' <', r'  |  ', '- '], # LD
                 [r'  |  ', ' <', r'  |  ', '> '], # UD
                 [r'  |  ', ' -', r'  ^  ', '> '], # UR
                 [r'  V  ', ' -', r'  |  ', '> ']] # RD
        ret = '  '
        # Do the top line
        for entry in self[0]:
            if entry == 1 or entry == 3 or entry == 4:
                ret += '  ^  '
            else:
                ret += '  |  '

        # Do the meat of the ascii art
        for row in self:
            ret += '\n  '
            # Do the top row
            for entry in row:
                ret += ascii[entry][0]
            ret += '\n'

            # Do the left-most entry
            if row[0] == 0 or row[0] == 1 or row[0] == 2:
                ret += '<-'
            else:
                ret += '--'

            # Do the middle row
            for entry in row:
                ret += ascii[entry][3] + '#' + ascii[entry][1]

            # Do the right-most entry
            if row[-1] == 0 or row[-1] == 4 or row[-1] == 5:
                ret += '->'
            else:
                ret += '--'

            # Do the bottom row
            ret += '\n  '
            for entry in row:
                ret += ascii[entry][2]

        # Do the bottom line
        ret += '\n  '
        for entry in self[-1]:
            if entry == 2 or entry == 3 or entry == 5:
                ret += '  V  '
            else:
                ret += '  |  '

        return ret

    def to_alternating_sign_matrix(self):
        """

        Returns the alternating sign matrix corresponding to this class.

         EXAMPLES::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl.to_alternating_sign_matrix()
            [ 0  1  0]
            [ 1 -1  1]
            [ 0  1  0]
            sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
            sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl.to_alternating_sign_matrix()
            [ 0  1  0  0]
            [ 0  0  1  0]
            [ 1 -1  0  1]
            [ 0  1  0  0]
        """
        return self.six_vertex_model.to_alternating_sign_matrix()

    def plot(self):
        """
        Return a graphical object of the Fully Packed Loop

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: print fpl.plot().description()
            Line defined by 2 points:       [(-1.0, 1.0), (0.0, 1.0)]
            Line defined by 2 points:       [(0.0, 0.0), (0.0, -1.0)]
            Line defined by 2 points:       [(0.0, 0.0), (1.0, 0.0)]
            Line defined by 2 points:       [(0.0, 2.0), (0.0, 3.0)]
            Line defined by 2 points:       [(0.0, 2.0), (0.0, 3.0)]
            Line defined by 2 points:       [(0.0, 2.0), (1.0, 2.0)]
            Line defined by 2 points:       [(1.0, 1.0), (0.0, 1.0)]
            Line defined by 2 points:       [(1.0, 1.0), (2.0, 1.0)]
            Line defined by 2 points:       [(2.0, 0.0), (1.0, 0.0)]
            Line defined by 2 points:       [(2.0, 0.0), (2.0, -1.0)]
            Line defined by 2 points:       [(2.0, 2.0), (1.0, 2.0)]
            Line defined by 2 points:       [(2.0, 2.0), (2.0, 3.0)]
            Line defined by 2 points:       [(2.0, 2.0), (2.0, 3.0)]
            Line defined by 2 points:       [(3.0, 1.0), (2.0, 1.0)]
            Line defined by 2 points:       [(3.0, 1.0), (2.0, 1.0)]

            sage: A = AlternatingSignMatrix([[0, 1, 0, 0], [1, -1, 0, 1], [0, 1, 0, 0],[0, 0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: print fpl.plot().description()
            Line defined by 2 points:       [(-1.0, 0.0), (0.0, 0.0)]
            Line defined by 2 points:       [(-1.0, 2.0), (0.0, 2.0)]
            Line defined by 2 points:       [(0.0, 1.0), (0.0, 0.0)]
            Line defined by 2 points:       [(0.0, 1.0), (1.0, 1.0)]
            Line defined by 2 points:       [(0.0, 3.0), (0.0, 4.0)]
            Line defined by 2 points:       [(0.0, 3.0), (0.0, 4.0)]
            Line defined by 2 points:       [(0.0, 3.0), (1.0, 3.0)]
            Line defined by 2 points:       [(1.0, 0.0), (1.0, -1.0)]
            Line defined by 2 points:       [(1.0, 0.0), (2.0, 0.0)]
            Line defined by 2 points:       [(1.0, 2.0), (0.0, 2.0)]
            Line defined by 2 points:       [(1.0, 2.0), (2.0, 2.0)]
            Line defined by 2 points:       [(2.0, 1.0), (1.0, 1.0)]
            Line defined by 2 points:       [(2.0, 1.0), (2.0, 2.0)]
            Line defined by 2 points:       [(2.0, 3.0), (1.0, 3.0)]
            Line defined by 2 points:       [(2.0, 3.0), (2.0, 4.0)]
            Line defined by 2 points:       [(2.0, 3.0), (2.0, 4.0)]
            Line defined by 2 points:       [(3.0, 0.0), (2.0, 0.0)]
            Line defined by 2 points:       [(3.0, 0.0), (3.0, -1.0)]
            Line defined by 2 points:       [(3.0, 2.0), (3.0, 1.0)]
            Line defined by 2 points:       [(3.0, 2.0), (3.0, 3.0)]
            Line defined by 2 points:       [(4.0, 1.0), (3.0, 1.0)]
            Line defined by 2 points:       [(4.0, 1.0), (3.0, 1.0)]
            Line defined by 2 points:       [(4.0, 3.0), (3.0, 3.0)]
            Line defined by 2 points:       [(4.0, 3.0), (3.0, 3.0)]

        """
        G = Graphics()
        n=len(self.six_vertex_model)-1
        for j,row in enumerate(reversed(self.six_vertex_model)):
            for i,entry in enumerate(row):
                if i == 0 and (i+j+n+1) % 2 ==0:
                    G+= line([(i-1,j),(i,j)])
                if i == n and (i+j+n+1) % 2 ==0:
                    G+= line([(i+1,j),(i,j)])
                if j == 0 and (i+j+n) % 2 ==0:
                    G+= line([(i,j),(i,j-1)])
                if j == n and (i+j+n) % 2 ==0:
                    G+= line([(i,j),(i,j+1)])
                if entry == 0: # LR
                    if (i+j+n) % 2==0:
                        G += line([(i,j), (i+1,j)])
                    else:
                        G += line([(i,j),(i,j+1)])
                elif entry == 1: # LU
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i,j+1)])
                    else:
                        G += line([(i+1,j), (i,j)])
                elif entry == 2: # LD
                    if (i+j+n) % 2 == 0:
                        pass
                    else:
                        G += line([(i,j+1), (i,j)])
                        G += line([(i+1,j), (i,j)])
                elif entry == 3: # UD
                    if (i+j+n) % 2 == 0:
                        G += line([(i,j), (i,j+1)])
                    else:
                        G += line([(i+1,j), (i,j)])
                elif entry == 4: # UR
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i,j+1)])
                        G += line([(i,j), (i+1,j)])
                    else:
                        pass
                elif entry == 5: # RD
                    if (i+j+n) % 2 ==0:
                        G += line([(i,j), (i+1,j)])
                    else:
                        G += line([(i,j+1), (i,j)])
        G.axes(False)
        return G
