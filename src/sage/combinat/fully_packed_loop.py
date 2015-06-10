r"""
Fully packed loops

A fully packed loop is a collection of non-intersecting lattice paths on a square grid such that every vertex is part of some path, and the paths are either closed internal loops or have endpoints corresponding to alternate points on the boundary [Propp2001]_. They are known to be in bijection with alternating sign matrices.

To each fully packed loop, we assign a link pattern, which is the non-crossing matching attained by seeing which points on the boundary are connected by open paths in the fully packed loop.

We can create a fully packed loop using the corresponding alternating sign
matrix and also extract the link pattern.::

    sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    sage: fpl = FullyPackedLoop(A)
    sage: fpl.link_pattern()
    [(1, 4), (2, 3), (5, 6)]
    sage: fpl
        |         |
        |         |
        + -- +    +
             |    |
             |    |
     -- +    +    + --
        |    |
        |    |
        +    + -- +
        |         |
        |         |
    sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    sage: fplb = FullyPackedLoop(B)
    sage: fplb.link_pattern()
    [(1, 6), (2, 5), (3, 4)]
    sage: fplb
        |         |
        |         |
        +    + -- +
        |    |
        |    |
     -- +    +    + --
             |    |
             |    |
        + -- +    +
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
        + -- +    +
             |    |
             |    |
     -- +    +    + --
        |    |
        |    |
        +    + -- +
        |         |
        |         |

Once we have a fully packed loop we can obtain the corresponding alternating sign matrix::

    sage: fpl.to_alternating_sign_matrix()
    [0 0 1]
    [0 1 0]
    [1 0 0]

Here are some more examples using bigger ASMs::

    sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
    sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
    sage: fpl = FullyPackedLoop(S)
    sage: fpl.link_pattern()
    [(1, 2), (3, 6), (4, 5), (7, 8)]
    sage: fpl
        |         |
        |         |
        + -- + -- +    + --
                       |
                       |
     -- +    + -- + -- +
        |    |
        |    |
        +    +    + -- + --
        |    |    |
        |    |    |
     -- +    +    + -- +
             |         |
             |         |

    sage: m = AlternatingSignMatrix([[0,0,1,0,0,0],
    ....:                            [1,0,-1,0,1,0],
    ....:                            [0,0,0,1,0,0],
    ....:                            [0,1,0,0,-1,1],
    ....:                            [0,0,0,0,1,0],
    ....:                            [0,0,1,0,0,0]])
    sage: fpl = FullyPackedLoop(m)
    sage: fpl.link_pattern()
    [(1, 12), (2, 7), (3, 4), (5, 6), (8, 9), (10, 11)]
    sage: fpl
        |         |         |
        |         |         |
        + -- +    +    + -- +    + --
             |    |    |         |
             |    |    |         |
     -- + -- +    +    + -- + -- +
                  |
                  |
        + -- +    + -- + -- +    + --
        |    |              |    |
        |    |              |    |
     -- +    +    + -- +    +    +
             |    |    |    |    |
             |    |    |    |    |
        + -- +    + -- +    +    + --
        |                   |
        |                   |
     -- +    + -- + -- +    + -- +
             |         |         |
             |         |         |

    sage: m = AlternatingSignMatrix([[0,1,0,0,0,0,0],
    ....:                            [1,-1,0,0,1,0,0],
    ....:                            [0,0,0,1,0,0,0],
    ....:                            [0,1,0,0,-1,1,0],
    ....:                            [0,0,0,0,1,0,0],
    ....:                            [0,0,1,0,-1,0,1],
    ....:                            [0,0,0,0,1,0,0]])
    sage: fpl = FullyPackedLoop(m)
    sage: fpl.link_pattern()
    [(1, 2), (3, 4), (5, 6), (7, 8), (9, 14), (10, 11), (12, 13)]
    sage: fpl
        |         |         |         |
        |         |         |         |
        + -- + -- +    + -- +    + -- +
                       |         |
                       |         |
     -- + -- + -- +    + -- + -- +    + --
                  |                   |
                  |                   |
        + -- +    + -- + -- +    + -- +
        |    |              |    |
        |    |              |    |
     -- +    +    + -- +    +    +    + --
             |    |    |    |    |    |
             |    |    |    |    |    |
        + -- +    + -- +    +    + -- +
        |                   |
        |                   |
     -- +    + -- + -- +    +    + -- + --
             |         |    |    |
             |         |    |    |
        + -- +    + -- +    +    + -- +
        |         |         |         |
        |         |         |         |


    Gyration on an alternating sign matrix/fully packed loop ``fpl``
    of the link pattern corresponding to ``fpl``::

        sage: ASMs = AlternatingSignMatrices(3).list()
        sage: ncp = FullyPackedLoop(ASMs[1]).link_pattern() # fpl's gyration orbit size is 2
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: fpl = FullyPackedLoop(ASMs[0])
        sage: ncp = fpl.link_pattern() # fpl's gyration size is 3
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,5):
        ....:         a,b=a%6+1,b%6+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],[0,0,1,0,0,0,0],\
        [0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
        sage: fpl = FullyPackedLoop(mat) # n=7
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,13):
        ....:         a,b=a%14+1,b%14+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

        sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0], [0,1,0,0,-1,1], [1,0,-1,1,0,0], \
        [0,0,1,0,0,0], [0,0,0,0,1,0]])
        sage: fpl = FullyPackedLoop(mat) # n =6
        sage: ncp = fpl.link_pattern()
        sage: rotated_ncp=[]
        sage: for (a,b) in ncp:
        ....:     for i in range(0,11):
        ....:         a,b=a%12+1,b%12+1;
        ....:     rotated_ncp.append((a,b))
        sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
        PerfectMatching(rotated_ncp)
        True

REFERENCES:

.. [Propp2001] James Propp.
   *The Many Faces of Alternating Sign Matrices*
   Discrete Mathematics and Theoretical Computer Science 43 (2001): 58

.. [Striker2015] Jessica Striker
   *The toggle group, homomesy, and the Razumov-Stroganov correspondence
   :arxiv:`abs/1503.08898`
"""
from sage.structure.sage_object import SageObject
from sage.combinat.six_vertex_model import SquareIceModel, SixVertexConfiguration
from sage.combinat.alternating_sign_matrix import AlternatingSignMatrix
from sage.plot.graphics import Graphics
from sage.matrix.constructor import matrix
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

        We can initiate a fully packed loop using an alternating sign matrix::

            sage: A = AlternatingSignMatrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl
                |         |
                |         |
                + -- +    +
                     |    |
                     |    |
             -- +    +    + --
                |    |
                |    |
                +    + -- +
                |         |
                |         |

        Otherwise we initiate a fully packed loop using a six vertex model::

            sage: S = SixVertexModel(3, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl
                |         |
                |         |
                + -- +    +
                     |    |
                     |    |
             -- +    +    + --
                |    |
                |    |
                +    + -- +
                |         |
                |         |
            sage: fpl.six_vertex_model().to_alternating_sign_matrix()
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
            self._six_vertex_model = generator.to_six_vertex_model()
        elif isinstance(generator, SixVertexConfiguration):
            self._six_vertex_model = generator
        else:
            raise TypeError('The generator for a fully packed loop must either be an AlternatingSignMatrix or a SixVertexConfiguration')

        self.end_points = self._end_point_dictionary()
        self.configuration = matrix(list(self._six_vertex_model))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: A = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |

            sage: A = AlternatingSignMatrix([[0,1,0,0],[0,0,1,0],[1,-1,0,1],[0,1,0,0]])
            sage: S = SixVertexModel(4, boundary_conditions='ice').from_alternating_sign_matrix(A)
            sage: fpl = FullyPackedLoop(S)
            sage: fpl
                |         |
                |         |
                + -- + -- +    + --
                               |
                               |
             -- +    + -- + -- +
                |    |
                |    |
                +    +    + -- + --
                |    |    |
                |    |    |
             -- +    +    + -- +
                     |         |
                     |         |

        """
        # List are in the order of URDL
        # One set of rules for how to draw around even vertex, one set of rules for odd vertex
        n=len(self._six_vertex_model)-1
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
        for i,entry in enumerate(self._six_vertex_model[0]):
            if i % 2 == 0:
                ret += '  |  '
            else:
                ret += '     '

        plus_sign = '+'

        # Do the meat of the ascii art
        for j,row in enumerate(self._six_vertex_model):
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
                    ret += ascii1[entry][3] + plus_sign + ascii1[entry][1]
                else:
                    ret += ascii2[entry][3] + plus_sign + ascii2[entry][1]

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
        for i,entry in enumerate(self._six_vertex_model[-1]):
            if (i+n+1) % 2 ==0:
                ret += '     '
            else:
                ret += '  |  '

        return ret

    def __eq__(self, other):
        """
        Check equality.

        EXAMPLES::

            sage: A = AlternatingSignMatrices(3)
            sage: M = A.random_element()
            sage: FullyPackedLoop(M) == M.to_fully_packed_loop()
            True

            sage: FullyPackedLoop(A([[1, 0, 0],[0, 1, 0],[0, 0, 1]])) ==\
            FullyPackedLoop(A([[1, 0, 0],[0, 0, 1],[0, 1, 0]]))
            False

            sage: FullyPackedLoop(M) == M
            False
        """
        return repr(self) == repr(other) and self.end_points == self.end_points

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
        return self._six_vertex_model.to_alternating_sign_matrix()

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

        .. MATH::

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

        .. MATH::

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

        .. MATH::

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

        .. MATH::

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

        .. MATH::

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

        .. MATH::

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
        n=len(self._six_vertex_model)-1
        for j,row in enumerate(reversed(self._six_vertex_model)):
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
        Return a list (a non-crossing partition) corresponding to a fully packed loop.
        Note: by convention, we choose the top left vertex to be even. See [Propp2001]_.

        EXAMPLES:

        We can extract the underlying link pattern (a non-crossing
        partition) from a fully packed loop::

            sage: A = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(A)
            sage: fpl.link_pattern()
            [(1, 2), (3, 6), (4, 5)]

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl.link_pattern()
            [(1, 6), (2, 5), (3, 4)]

        Gyration on an alternating sign matrix/fully packed loop ``fpl``
        corresponds to a rotation (i.e. a becomes a-1 mod 2n)
        of the link pattern corresponding to ``fpl``::

            sage: ASMs = AlternatingSignMatrices(3).list()
            sage: ncp = FullyPackedLoop(ASMs[1]).link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[1].gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: fpl = FullyPackedLoop(ASMs[0])
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,5):
            ....:         a,b=a%6+1,b%6+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(ASMs[0].gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,1,0,0,0,0],[1,0,-1,0,1,0,0],[0,0,1,0,0,0,0],\
            [0,1,-1,0,0,1,0],[0,0,1,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,0,1]])
            sage: fpl = FullyPackedLoop(mat) # n=7
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,13):
            ....:         a,b=a%14+1,b%14+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

            sage: mat = AlternatingSignMatrix([[0,0,0,1,0,0], [0,0,1,-1,1,0], [0,1,0,0,-1,1], [1,0,-1,1,0,0], \
            [0,0,1,0,0,0], [0,0,0,0,1,0]])
            sage: fpl = FullyPackedLoop(mat)
            sage: ncp = fpl.link_pattern()
            sage: rotated_ncp=[]
            sage: for (a,b) in ncp:
            ....:     for i in range(0,11):
            ....:         a,b=a%12+1,b%12+1;
            ....:     rotated_ncp.append((a,b))
            sage: PerfectMatching(mat.gyration().to_fully_packed_loop().link_pattern()) ==\
            PerfectMatching(rotated_ncp)
            True

        TESTS:

        We test a previous bug which showed up when this method is called twice::

            sage: A = AlternatingSignMatrices(6)
            sage: B = A.random_element()
            sage: C = FullyPackedLoop(B)
            sage: D = C.link_pattern()
            sage: E = C.link_pattern()
            sage: D == E
            True

        """
        link_pattern = []
        boundary_d = self.end_points.copy()
        vertices_d = self._vertex_dictionary()

        while len(boundary_d) > 2:
            startpoint = boundary_d.keys()[0]
            position = boundary_d[startpoint]

            boundary_d.pop(startpoint)
            vertices_d[position] = False # allows us to start

            while not vertices_d[position]:
                vertices_d.pop(position)
                choices = self._get_coordinates(position)
                if choices[0] in vertices_d:
                    position = choices[0]
                elif choices[1] in vertices_d:
                    position = choices[1]
                else:
                    raise ValueError('No valid choices')

            endpoint = vertices_d[position]
            vertices_d.pop(position)
            link_pattern.append((startpoint, endpoint))
            boundary_d.pop(endpoint)

        link_pattern.append(tuple(boundary_d.keys()))

        return link_pattern

    def six_vertex_model(self):
        """
        Return the underlying six vertex model configuration

        EXAMPLES::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl
                |         |
                |         |
                +    + -- +
                |    |
                |    |
             -- +    +    + --
                     |    |
                     |    |
                + -- +    +
                |         |
                |         |
            sage: fpl.six_vertex_model()
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
        return self._six_vertex_model

    def _vertex_dictionary(self):
        """
        A function to create a dictionary of all the coordinates.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl._vertex_dictionary()
            {(0, 0): 1,
             (0, 1): 0,
             (0, 2): 2,
             (1, 0): 6,
             (1, 1): 0,
             (1, 2): 3,
             (2, 0): 5,
             (2, 1): 0,
             (2, 2): 4}
        """
        n = len(self._six_vertex_model)
        vertices = {}
        for i in range(n):
            for j in range(n):
                vertices[(i, j)] = 0

        for end, vertex in self.end_points.iteritems():
            vertices[vertex] = end

        return vertices

    def _get_coordinates(self, current):
        """
        Returns a list of 2 coordinates that refer to the moves that could
        potentialy be made.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: matrix(list(fpl._six_vertex_model))
            [3 1 1]
            [5 3 1]
            [5 5 3]
            sage: fpl._get_coordinates((0, 1))
            [(1, 1), (0, 2)]
            sage: fpl._get_coordinates((1, 1))
            [(2, 1), (0, 1)]
            sage: fpl._get_coordinates((0, 0))
            [(1, 0), (-1, 0)]
            sage: fpl._get_coordinates((0, 2))
            [(-1, 2), (0, 1)]
            sage: fpl._get_coordinates((2, 1))
            [(1, 1), (2, 0)]


            sage: B = AlternatingSignMatrix([[0, 1, 0], [1, -1, 1], [0, 1, 0]])
            sage: fpl = FullyPackedLoop(B)
            sage: matrix(list(fpl._six_vertex_model))
            [4 3 1]
            [3 0 3]
            [5 3 2]
            sage: fpl._get_coordinates((1, 1))
            [(1, 0), (1, 2)]
            sage: fpl._get_coordinates((0, 0))
            [(-1, 0), (0, 1)]
            sage: fpl._get_coordinates((0, 2))
            [(-1, 2), (0, 1)]
            sage: fpl._get_coordinates((2, 1))
            [(2, 0), (2, 2)]
        """
        # 0 UD, 1 RD, 2 UR, 3 LR, 4 LD, 5 LU
        odd = {0: [(1, 0), (-1, 0)],
                1: [(1, 0), (0, 1)],
                2: [(-1, 0), (0, 1)],
                3: [(0, -1), (0, 1)],
                4: [(0, -1), (1, 0)],
                5: [(-1, 0), (0, -1)]
                }

        # 0 LR, 1 LU, 2 LD, 3 UD, 4 UR, 5 RD
        even = {0: [(0, -1), (0, 1)],
               1: [(-1, 0), (0, -1)],
               2: [(0, -1), (1, 0)],
               3: [(1, 0), (-1, 0)],
               4: [(-1, 0), (0, 1)],
               5: [(1, 0), (0, 1)]
               }

        parity = sum(current)
        conf = self.configuration[current]

        if parity % 2 == 0:
            potential_directions = even[conf]
        else:
            potential_directions = odd[conf]

        return [(current[0] + d[0], current[1] + d[1]) for d in potential_directions]

    def _end_point_dictionary(self):
        r"""
        A function create a dictionary of the endpoints and their
        coordinates.

        TESTS::

            sage: B = AlternatingSignMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: fpl = FullyPackedLoop(B)
            sage: fpl._end_point_dictionary()
            {1: (0, 0), 2: (0, 2), 3: (1, 2), 4: (2, 2), 5: (2, 0), 6: (1, 0)}

        """
        n = len(self._six_vertex_model)
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
