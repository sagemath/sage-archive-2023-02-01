class Parser():
    r"""
    A class for parsing the outputs of different algorithms called in other
    software packages.

    At present the only parser included is the for `lrs` algorithm however
    this will be expanded to `gambit` soon.
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
            sage: process = Popen(['nash', g1_name, g2_name], stdout=PIPE)
            sage: lrs_output = [row for row in process.stdout]

        The above creates a game, writes the H representation to
        temporary files, calls lrs and stores the output in `lrs_output`
        (here slicing to get rid of some system parameters that get returned)::

            sage: lrs_output[5:-4]
            ['\n', '***** 4 4 rational\n', '2  0  1  2 \n', '1  1/2  1/2 -2 \n', '\n', '2  0  1  2 \n', '1  0  1 -2 \n', '\n', '*Number of equilibria found: 2\n', '*Player 1: vertices=3 bases=3 pivots=5\n', '*Player 2: vertices=2 bases=1 pivots=6\n']

        The above is pretty messy, here is the output when we put it through
        the parser::

            sage: nasheq = Parser(lrs_output).format_lrs()
            sage: nasheq
            [[(1/2, 1/2), (0, 1)], [(0, 1), (0, 1)]]
        """
        from sage.misc.sage_eval import sage_eval
        p2_strategies = []
        p1_strategies = []
        for i in self.raw_string:
            if i.startswith('2'):
                nums = [sage_eval(k) for k in i.split()]
                p2_strategies.append(tuple(nums[1:-1]))
            elif i.startswith('1'):
                nums = [sage_eval(k) for k in i.split()]
                p1_strategies.append(tuple(nums[1:-1]))

        return [list(a) for a in zip(p1_strategies, p2_strategies)]
