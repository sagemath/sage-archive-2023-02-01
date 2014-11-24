"""
Gambit

This file contains some information and tests for the use of
`gambit<http://www.gambit-project.org/>`_ as a stand alone package.

To install gambit as an optional package run (from root of sage)::

    $ ./sage -i gambit

The `python API documentation for gambit
<http://www.gambit-project.org/gambit14/pyapi.html>_` shows various examples
however due to the preparser in Sage: integers must first be
converted to Python integers.

Here is an example showing a construction and solution of the Prisoners
Dilemma::

    sage: import gambit
    sage: g = gambit.Game.new_table([2,2])
    sage: g[int(0), int(0)][int(0)] = int(8)
    sage: g[int(0), int(0)][int(1)] = int(8)
    sage: g[int(0), int(1)][int(0)] = int(2)
    sage: g[int(0), int(1)][int(1)] = int(10)
    sage: g[int(1), int(0)][int(0)] = int(10)
    sage: g[int(1), int(0)][int(1)] = int(2)
    sage: g[int(1), int(1)][int(0)] = int(5)
    sage: g[int(1), int(1)][int(1)] = int(5)
    sage: solver = gambit.nash.ExternalLCPSolver()
    sage: solver.solve(g)
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

    sage: solver = gambit.nash.ExternalEnumPureSolver()
    sage: solver.solve(g)
    [<NashProfile for '': [Fraction(0, 1), Fraction(1, 1), Fraction(0, 1), Fraction(1, 1)]>]

Note that the above example used an algorithm that enumerates all pure
strategies. This will fail to find all Nash equilibria in certain games.
For example here is an implementation of Matching Pennies::


    sage: g = gambit.Game.new_table([2,2])
    sage: g[int(0), int(0)][int(0)] = int(1)
    sage: g[int(0), int(0)][int(1)] = int(-1)
    sage: g[int(0), int(1)][int(0)] = int(-1)
    sage: g[int(0), int(1)][int(1)] = int(1)
    sage: g[int(1), int(0)][int(0)] = int(-1)
    sage: g[int(1), int(0)][int(1)] = int(1)
    sage: g[int(1), int(1)][int(0)] = int(1)
    sage: g[int(1), int(1)][int(1)] = int(-1)
    sage: solver = gambit.nash.ExternalEnumPureSolver()
    sage: solver.solve(g)
    []

If we solve this with the `LCP` solver we get the expected Nash equilibrium::

    sage: solver = gambit.nash.ExternalLCPSolver()
    sage: solver.solve(g)
    [<NashProfile for '': [0.5, 0.5, 0.5, 0.5]>]


Here is another example showing the Battle of the Sexes game::

    sage: g = gambit.Game.new_table([2,2])
    sage: g[int(0), int(0)][int(0)] = int(2)
    sage: g[int(0), int(0)][int(1)] = int(1)
    sage: g[int(0), int(1)][int(0)] = int(0)
    sage: g[int(0), int(1)][int(1)] = int(0)
    sage: g[int(1), int(0)][int(0)] = int(0)
    sage: g[int(1), int(0)][int(1)] = int(0)
    sage: g[int(1), int(1)][int(0)] = int(1)
    sage: g[int(1), int(1)][int(1)] = int(2)
    sage: solver = gambit.nash.ExternalLCPSolver()
    sage: solver.solve(g)
    [<NashProfile for '': [1.0, 0.0, 1.0, 0.0]>,
     <NashProfile for '': [0.6666666667, 0.3333333333, 0.3333333333, 0.6666666667]>,
     <NashProfile for '': [0.0, 1.0, 0.0, 1.0]>]

Note that all of the above examples can be used in Ipython without the need to
convert integers to instances of `int`.

AUTHOR:

- Vince Knight (11-2014): Original version

"""
