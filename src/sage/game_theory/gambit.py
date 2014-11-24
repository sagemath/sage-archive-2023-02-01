"""
Gambit

This file contains some information and tests for the use of
`gambit<http://www.gambit-project.org/>`_ as a stand alone package.

To install gambit as an optional package run (from root of sage)::

    $ ./sage -i gambit

The `python API documentation for gambit
<http://www.gambit-project.org/gambit14/pyapi.html>_` shows various examples
that can be run easily in an Ipython notebook::

    In [1]: g.players[0].strategies
    Out[1]: [<Strategy [0] '1' for player 'Alphonse' in game 'A
    prisoner's dilemma game'>,
             <Strategy [1] '2' for player 'Alphonse' in game 'A prisoner's dilemma game'>]
    In [2]: len(g.players[0].strategies)
    Out[2]: 2

    In [3]: g.players[0].strategies[0].label = "Cooperate"
    In [4]: g.players[0].strategies[1].label = "Defect"
    In [5]: g.players[0].strategies
    Out[5]: [<Strategy [0] 'Cooperate' for player 'Alphonse' in game 'A
    prisoner's dilemma game'>,
       <Strategy [1] 'Defect' for player 'Alphonse' in game 'A prisoner's dilemma game'>]

    In [6]: g[0,0][0] = 8
    In [7]: g[0,0][1] = 8
    In [8]: g[0,1][0] = 2
    In [9]: g[0,1][1] = 10
    In [10]: g[1,0][0] = 10
    In [11]: g[1,1][1] = 2
    In [12]: g[1,0][1] = 2
    In [13]: g[1,1][0] = 5
    In [14]: g[1,1][1] = 5
    In [15]: solver = gambit.nash.ExternalEnumPureSolver()
    In [16]: solver.solve(g)
    Out[16]: [[1.0, 0.0, 0.0, 1.0, 0.0]]

If one really wants to use gambit directly in Sage then integers must first be
converted to Python integers (due to the preparser)::

Here is an example showing a construction and solution of the Prisoners
Dilemma::

    sage: import gambit  # optional - gambit
    sage: g = gambit.Game.new_table([2,2])  # optional - gambit
    sage: g[int(0), int(0)][int(0)] = int(8)  # optional - gambit
    sage: g[int(0), int(0)][int(1)] = int(8)  # optional - gambit
    sage: g[int(0), int(1)][int(0)] = int(2)  # optional - gambit
    sage: g[int(0), int(1)][int(1)] = int(10)  # optional - gambit
    sage: g[int(1), int(0)][int(0)] = int(10)  # optional - gambit
    sage: g[int(1), int(0)][int(1)] = int(2)  # optional - gambit
    sage: g[int(1), int(1)][int(0)] = int(5)  # optional - gambit
    sage: g[int(1), int(1)][int(1)] = int(5)  # optional - gambit
    sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit
    sage: solver.solve(g)  # optional - gambit
    [<NashProfile for '': [0.0, 1.0, 0.0, 1.0]>]


Here is a list of various other solvers that can be used:

- ExternalEnumPureSolver
- ExternalEnumMixedSolver
- ExternalLPSolver
- ExternalLCPSolver
- ExternalSimpdivSolver
- ExternalGlobalNewtonSolver
- ExternalEnumPolySolver
- ExternalLyapunovSolver
- ExternalIteratedPolymatrixSolver
- ExternalLogitSolver

    sage: solver = gambit.nash.ExternalEnumPureSolver()  # optional - gambit
    sage: solver.solve(g)  # optional - gambit
    [<NashProfile for '': [Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)]>]

Note that the above example used an algorithm that enumerates all pure
strategies. This will fail to find all Nash equilibria in certain games.
For example here is an implementation of Matching Pennies::


    sage: g = gambit.Game.new_table([2,2])  # optional - gambit
    sage: g[int(0), int(0)][int(0)] = int(1)  # optional - gambit
    sage: g[int(0), int(0)][int(1)] = int(-1)  # optional - gambit
    sage: g[int(0), int(1)][int(0)] = int(-1)  # optional - gambit
    sage: g[int(0), int(1)][int(1)] = int(1)  # optional - gambit
    sage: g[int(1), int(0)][int(0)] = int(-1)  # optional - gambit
    sage: g[int(1), int(0)][int(1)] = int(1)  # optional - gambit
    sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
    sage: g[int(1), int(1)][int(1)] = int(-1)  # optional - gambit
    sage: solver = gambit.nash.ExternalEnumPureSolver()  # optional - gambit
    sage: solver.solve(g)  # optional - gambit
    []

If we solve this with the `LCP` solver we get the expected Nash equilibrium::

    sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit
    sage: solver.solve(g)  # optional - gambit
    [<NashProfile for '': [0.5, 0.5, 0.5, 0.5]>]


Here is another example showing the Battle of the Sexes game::

    sage: g = gambit.Game.new_table([2,2])  # optional - gambit
    sage: g[int(0), int(0)][int(0)] = int(2)  # optional - gambit
    sage: g[int(0), int(0)][int(1)] = int(1)  # optional - gambit
    sage: g[int(0), int(1)][int(0)] = int(0)  # optional - gambit
    sage: g[int(0), int(1)][int(1)] = int(0)  # optional - gambit
    sage: g[int(1), int(0)][int(0)] = int(0)  # optional - gambit
    sage: g[int(1), int(0)][int(1)] = int(0)  # optional - gambit
    sage: g[int(1), int(1)][int(0)] = int(1)  # optional - gambit
    sage: g[int(1), int(1)][int(1)] = int(2)  # optional - gambit
    sage: solver = gambit.nash.ExternalLCPSolver()  # optional - gambit
    sage: solver.solve(g)  # optional - gambit
    [<NashProfile for '': [1.0, 0.0, 1.0, 0.0]>,
     <NashProfile for '': [0.6666666667, 0.3333333333, 0.3333333333, 0.6666666667]>,
     <NashProfile for '': [0.0, 1.0, 0.0, 1.0]>]

AUTHOR:

- Vince Knight (11-2014): Original version

"""
