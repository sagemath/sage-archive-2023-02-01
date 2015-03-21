r"""
Fully packed loops

A class for fully packed loops [Propp2001]_.
We can create a fully packed loop using the corresponding alternating sign matrix::

    sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    sage: fpl = FullyPackedLoop(A)
    sage: fpl
        |         |
        |         |
        # -- #    #
             |    |
             |    |
     -- #    #    # --
        |    |
        |    |
        #    # -- #
        |         |
        |         |

The class also has a plot method::

    sage: fpl.plot()
    Graphics object consisting of 15 graphics primitives

which gives:

.. PLOT::
    :width: 200 px

    A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    fpl = FullyPackedLoop(A)
    p = fpl.plot()
    sphinx_plot(p)

Note that we can also create a fully packed loop from a six vertex model configuration::

    sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
    sage: S
        ^    ^    ^
        |    |    |
    --> # -> # -> # <--
        ^    ^    |
        |    |    V
    --> # -> # <- # <--
        ^    |    |
        |    V    V
    --> # <- # <- # <--
        |    |    |
        V    V    V
    sage: fpl = FullyPackedLoop(S)
    sage: fpl
        |         |
        |         |
        # -- #    #
             |    |
             |    |
     -- #    #    # --
        |    |
        |    |
        #    # -- #
        |         |
        |         |

Once we have a fully packed loop we can obtain the corresponding alternating sign matrix::

    sage: fpl.to_alternating_sign_matrix()
    [0 0 1]
    [0 1 0]
    [1 0 0]

REFERENCES:

.. [Propp2001] James Propp.
   *The Many Faces of Alternating Sign Matrices*
   Discrete Mathematics and Theoretical Computer Science 43 (2001): 58
"""
from sage.structure.sage_object import SageObject
from sage.combinat.six_vertex_model import SquareIceModel, SixVertexConfiguration
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
from sage.plot.graphics import Graphics
from sage.plot.line import line
from sage.combinat.perfect_matching import PerfectMatching


class FullyPackedLoop(SageObject):
    """
    A class for fully packed loops
    """

    def __init__(self, generator):
        """
        Initialise object, can take ASM of FPL as generator.

        EXAMPLES:

        We can initiate a fully packed loop using an Alternating Sign Matrix::

            sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl
                |         |
                |         |
                # -- #    #
                     |    |
                     |    |
             -- #    #    # --
                |    |
                |    |
                #    # -- #
                |         |
                |         |

        Otherwise we initiate a fully packed loop using a six vertex model::

            sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl
                |         |
                |         |
                # -- #    #
                     |    |
                     |    |
             -- #    #    # --
                |    |
                |    |
                #    # -- #
                |         |
                |         |
            sage: fpl.six_vertex_model.to_alternating_sign_matrix()
            [0 0 1]
            [0 1 0]
            [1 0 0]

        Note that if anything else is used to generate the fully packed loop an error will occur::

            sage: fpl = FullyPackedLoop(5)
            Traceback (most recent call last):
            ...
            TypeError: The generator for a fully packed loop must either be an AlternatingSignMatrix or a SixVertexConfiguration

            sage: fpl = FullyPackedLoop((1, 2, 3))
            Traceback (most recent call last):
            ...
            TypeError: The generator for a fully packed loop must either be an AlternatingSignMatrix or a SixVertexConfiguration

        """
        if isinstance(generator, AlternatingSignMatrix):
            self.six_vertex_model = generator.to_six_vertex_model()
        elif isinstance(generator, SixVertexConfiguration):
            self.six_vertex_model = generator
        else:
            raise TypeError('The generator for a fully packed loop must either be an AlternatingSignMatrix or a SixVertexConfiguration')

        self.end_points = self._end_point_dictionary()
        self.configuration = list(self.six_vertex_model)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl
                |         |
                |         |
                #    # -- #
                |    |
                |    |
             -- #    #    # --
                     |    |
                     |    |
                # -- #    #
                |         |
                |         |

            sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
            sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl
                |         |
                |         |
                # -- # -- #    # --
                               |
                               |
             -- #    # -- # -- #
                |    |
                |    |
                #    #    # -- # --
                |    |    |
                |    |    |
             -- #    #    # -- #
                     |         |
                     |         |

        """
        # List are in the order of URDL
        # One set of rules for how to draw around even vertex, one set of rules for odd vertex
        n=len(self.six_vertex_model)-1
        ascii1 = [[r'     ', ' -', r'     ', '- '], # LR
                 [r'  |  ', '  ', r'     ', '- '], # LU
                 [r'     ', '  ', r'  |  ', '- '], # LD
                 [r'  |  ', '  ', r'  |  ', '  '], # UD
                 [r'  |  ', ' -', r'     ', '  '], # UR
                 [r'     ', ' -', r'  |  ', '  ']] # RD

        ascii2 = [[r'  |  ', '  ', r'  |  ', '  '], # LR
                 [r'     ', ' -', r'  |  ', '  '], # LU
                 [r'  |  ', ' -', r'     ', '  '], # LD
                 [r'     ', ' -', r'     ', '- '], # UD
                 [r'     ', '  ', r'  |  ', '- '], # UR
                 [r'  |  ', '  ', r'     ', '- ']] # RD
        ret = '  '
        # Do the top line
        for i,entry in enumerate(self.six_vertex_model[0]):
            if i % 2 == 0:
                ret += '  |  '
            else:
                ret += '     '
#            if entry == 1 or entry == 3 or entry == 4:
#                ret += '  ^  '
#            else:
#                ret += '  |  '

        # Do the meat of the ascii art
        for j,row in enumerate(self.six_vertex_model):
            ret += '\n  '
            # Do the top row
            for i,entry in enumerate(row):
                if (i+j) % 2 == 0:
                    ret += ascii1[entry][0]
                else:
                    ret += ascii2[entry][0]
            ret += '\n'

            # Do the left-most entry
            if (j) % 2 == 0:
                ret += '  '
            else:
                ret += ' -'

            # Do the middle row
            for i,entry in enumerate(row):
                if (i+j) % 2 == 0:
                    ret += ascii1[entry][3] + '#' + ascii1[entry][1]
                else:
                    ret += ascii2[entry][3] + '#' + ascii2[entry][1]

            # Do the right-most entry
            if (j+n) % 2 ==0:
                ret += '  '
            else:
                ret += '- '

            # Do the bottom row
            ret += '\n  '
            for i,entry in enumerate(row):
                if (i+j) % 2 ==0:
                    ret += ascii1[entry][2]
                else:
                    ret += ascii2[entry][2]

        # Do the bottom line
        ret += '\n  '
        for i,entry in enumerate(self.six_vertex_model[-1]):
            if (i+n+1) % 2 ==0:
                ret += '     '
            else:
                ret += '  |  '



#            if entry == 2 or entry == 3 or entry == 5:
#                ret += '  V  '
#            else:
#                ret += '  |  '

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

        EXAMPLES:

        Here is the fully packed for :math:`\\begin{pmatrix}0&1&1\\\\1&-1&1\\\\0&1&0\end{pmatrix}`:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        Here is how Sage represents this::

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

        Here are the other 3 by 3 Alternating Sign Matrices and their corresponding fully packed loops:

        .. math::

            A = \\begin{pmatrix}
                1&0&0\\\\
                0&1&0\\\\
                0&0&1\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. math::

            A = \\begin{pmatrix}
                1&0&0\\\\
                0&0&1\\\\
                0&1&0\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[1, 0, 0], [0, 0, 1], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. math::

            A = \\begin{pmatrix}
                0&1&0\\\\
                1&0&0\\\\
                0&0&1\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. math::

            A = \\begin{pmatrix}
                0&1&0\\\\
                0&0&1\\\\
                1&0&0\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. math::

            A = \\begin{pmatrix}
                0&0&1\\\\
                1&0&0\\\\
                0&1&0\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        .. math::

            A = \\begin{pmatrix}
                0&0&1\\\\
                0&1&0\\\\
                1&0&0\\\\
                \end{pmatrix}

        gives:

        .. PLOT::
            :width: 200 px

            A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

        EXAMPLES::

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

        Here is the plot:

        .. PLOT::
            :width: 300 px

            A = AlternatingSignMatrix([[0, 1, 0, 0], [1, -1, 0, 1], [0, 1, 0, 0],[0, 0, 1, 0]])
            fpl = FullyPackedLoop(A)
            p = fpl.plot()
            sphinx_plot(p)

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

    def link_pattern(self):
        """
        Return a :class:`PerfectMatching` class (a non-crossing partition)
        corresponding to a fully packed loop. Note: by convention, we
        choose the top left vertex to be even. See [Propp2001]_.

        EXAMPLES:

        We can extract the underlying link pattern (a non-crossing
        partition) from a fully packed loop::
 
            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.link_pattern()
            [(1, 2), (3, 6), (4, 5)]
        """
        link_pattern = []
        svm = self.six_vertex_model
        n = len(svm)
        boundary_d = self.end_points
        vertices_d = self._vertex_dictionary(n)

        while len(boundary_d) > 2:
            startpoint = boundary_d.keys()[0]
            position = boundary_d[startpoint]

            boundary_d.pop(startpoint)
            vertices_d[position] = 0 # allows us to start

            while not vertices_d[position]:
                vertices_d.pop(position)
                choices = self._get_coordinates(position)
                if choices[0] in vertices_d:
                    position = choices[0]
                else:
                    position = choices[1]

            endpoint = vertices_d[position]
            vertices_d.pop(position)
            link_pattern.append((startpoint, endpoint))
            boundary_d.pop(endpoint)

        link_pattern.append(boundary_d.keys())

        return PerfectMatching(link_pattern)

    def _vertex_dictionary(self, size):
        """
        A function to create a dictionary of all the coordinates.
        """
        n = len(self.six_vertex_model)
        vertices = {}
        for i in range(n):
            for j in range(n):
                vertices[(i, j)] = 0

        for end, vertex in self.end_points:
            vertices[vertex] = end

        return vertices

    def _get_coordinates(self, current):
        # 0 UD, 1 RD, 2 UR, 3 LR, 4 LD, 5 LU
        odd = {0: [[-1, 0], [1, 0]],
               1: [[0, 1], [1, 0]],
               2: [[-1, 0], [0, 1]],
               3: [[0, -1], [0, 1]],
               4: [[0, -1], [1, 0]],
               5: [[0, -1], [-1, 0]]
               }

        # 0 LR, 1 LU, 2 LD, 3 UD, 4 UR, 5 RD
        even = {0: [[0, -1], [0, 1]],
                1: [[0, -1], [-1, 0]],
                2: [[0, -1], [1, 0]],
                3: [[-1, 0], [0, 1]],
                4: [[-1, 0], [1, 0]],
                5: [[0, 1], [1, 0]]
                }

        parity = sum(current)
        c = self.configuration(current)

        if parity % 2 == 0:
            potential_directions = even[parity]
        else:
            potential_directions = odd[parity]

        return [(c[0] + d[0], c[1] + d[1]) for d in potential_directions]

    def _end_point_dictionary(self):
        """
        A function create a dictionary of the endpoints and their
        coordinates.
        """
        n = len(self.six_vertex_model)
        end_points = {}

        for k in range(n):
            if k % 2 == 0:
                # top row
                end_points[1 + k/2] = (0, k)

                # bottom row
                end_points[n + 1 + k/2] = (n-1, n-1-k)

        # sides for even case
        if n % 2 == 0:
            for k in range(n):
                if k % 2 == 0:
                    # left side
                    end_points[((3*n + 2 + k)/2)] = (n-1-k, 0)

                    # right side
                    end_points[(n + 2 + k)/2] = (k, n-1)

        # side for odd case
        if n % 2 == 1:
            for k in range(n):
                if k % 2 == 1:
                    # left side
                    end_points[(3*n + 2 + k)/2] = (n-1-k, 0)

                    # right side
                    end_points[(n + 2 + k)/2] = (k, n-1)

        return end_points
