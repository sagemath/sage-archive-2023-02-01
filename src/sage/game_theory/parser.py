class Parser():
    r"""
    A class for parsing the outputs of different algorithms called in other
    software packages.

    At present the only parser included is for the ``'lrs'`` algorithm however
    this is actively being expanded to 'gambit'.
    """

    def __init__(self, raw_string):
        """
        Initialise a Parser instance by storing a raw_string
        (currently only used with H representation of a game).

        TESTS:

        Simply checking that we have the correct string output
        for the H representation::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: raw_string = g._Hrepresentation(A, -A)
            sage: P = Parser(raw_string)
            sage: print P.raw_string[0]
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 1 0 0
            0 0 1 0
            0 1 3  1
            0 2 2  1
            -1 1 1 0
            end
            <BLANKLINE>

            sage: print P.raw_string[1]
            H-representation
            linearity 1 5
            begin
            5 4 rational
            0 -1 -2  1
            0 -3 -2  1
            0 1 0 0
            0 0 1 0
            -1 1 1 0
            end
            <BLANKLINE>

        The specific case of a game with 1 strategy for each player::

            sage: A = matrix([[1]])
            sage: B = matrix([[5]])
            sage: g = NormalFormGame([A,B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print P.raw_string[0]
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 1 0
            0 -5  1
            -1 1 0
            end
            <BLANKLINE>

            sage: print P.raw_string[1]
            H-representation
            linearity 1 3
            begin
            3 3 rational
            0 -1  1
            0 1 0
            -1 1 0
            end
            <BLANKLINE>

        Another test::

            sage: from sage.game_theory.parser import Parser
            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: raw_string = g._Hrepresentation(A, B)
            sage: P = Parser(raw_string)
            sage: print P.raw_string[0]
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            0  9 -6  4  1
            0 -7  2 -6  1
            0 -9  3  10  1
            -1 1 1 1 0
            end
            <BLANKLINE>

            sage: print P.raw_string[1]
            H-representation
            linearity 1 7
            begin
            7 5 rational
            0 7 5 -5  1
            0 -5 -5 -3  1
            0 -1 6 -1  1
            0 1 0 0 0
            0 0 1 0 0
            0 0 0 1 0
            -1 1 1 1 0
            end
            <BLANKLINE>
        """
        self.raw_string = raw_string

    def format_lrs(self):
        """
        Parses the output of lrs so as to return vectors
        corresponding to equilibria.

        TESTS::

            sage: from sage.game_theory.parser import Parser
            sage: from subprocess import Popen, PIPE
            sage: A = matrix([[1, 2], [3, 2]])
            sage: g = NormalFormGame([A])
            sage: game1_str, game2_str = g._Hrepresentation(A, -A)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = file(g1_name, 'w')
            sage: g2_file = file(g2_name, 'w')
            sage: g1_file.write(game1_str)
            sage: g1_file.close()
            sage: g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['nash', g1_name, g2_name], stdout=PIPE)  # optional - lrs
            sage: lrs_output = [row for row in process.stdout]  # optional - lrs

        The above creates a game, writes the H representation to
        temporary files, calls lrs and stores the output in `lrs_output`
        (here slicing to get rid of some system parameters that get returned)::

            sage: lrs_output[5:-4]  # optional - lrs
            ['\n', '***** 4 4 rational\n', '2  0  1  2 \n', '1  1/2  1/2 -2 \n', '\n', '2  0  1  2 \n', '1  0  1 -2 \n', '\n', '*Number of equilibria found: 2\n', '*Player 1: vertices=3 bases=3 pivots=5\n', '*Player 2: vertices=2 bases=1 pivots=6\n']

        The above is pretty messy, here is the output when we put it through
        the parser::

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrs
            sage: nasheq  # optional - lrs
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]

        Another game::

            sage: A = matrix([[-7, -5, 5],
            ....:             [5, 5, 3],
            ....:             [1, -6, 1]])
            sage: B = matrix([[-9, 7, 9],
            ....:             [6, -2, -3],
            ....:             [-4, 6, -10]])
            sage: g = NormalFormGame([A, B])
            sage: game1_str, game2_str = g._Hrepresentation(A, B)
            sage: g1_name = tmp_filename()
            sage: g2_name = tmp_filename()
            sage: g1_file = file(g1_name, 'w')
            sage: g2_file = file(g2_name, 'w')
            sage: g1_file.write(game1_str)
            sage: g1_file.close()
            sage: g2_file.write(game2_str)
            sage: g2_file.close()
            sage: process = Popen(['nash', g1_name, g2_name], stdout=PIPE)  # optional - lrs
            sage: lrs_output = [row for row in process.stdout]  # optional - lrs
            sage: print lrs_output[5:-4]  # optional - lrs
            ['\n', '***** 5 5 rational\n', '2  0  1/6  5/6  10/3 \n', '2  1/7  0  6/7  23/7 \n', '1  1/3  2/3  0  1 \n', '\n', '2  0  0  1  5 \n', '1  1  0  0  9 \n', '\n', '2  1  0  0  5 \n', '1  0  1  0  6 \n', '\n', '*Number of equilibria found: 4\n', '*Player 1: vertices=6 bases=7 pivots=10\n', '*Player 2: vertices=4 bases=2 pivots=14\n']

            sage: nasheq = Parser(lrs_output).format_lrs()  # optional - lrs
            sage: nasheq  # optional - lrs
            [[(1/3, 2/3, 0), (0, 1/6, 5/6)], [(1/3, 2/3, 0), (1/7, 0, 6/7)], [(1, 0, 0), (0, 0, 1)], [(0, 1, 0), (1, 0, 0)]]
        """
        equilibria = []
        from sage.misc.sage_eval import sage_eval
        from itertools import groupby
        for collection in [list(x[1]) for x in groupby(self.raw_string[7:], lambda x: x=='\n')]:
            if collection[0].startswith('2'):
                s1 = tuple([sage_eval(k) for k in collection[-1].split()][1:-1])
                for s2 in collection[:-1]:
                    s2 = tuple([sage_eval(k) for k in s2.split()][1:-1])
                    equilibria.append([s1, s2])

        return equilibria



