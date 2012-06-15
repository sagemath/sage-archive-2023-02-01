"""
SAT-Solvers via DIMACS Files.

Sage supports calling SAT solvers using the popular DIMACS
format. This module implements infrastructure to make it easy to add
new such interfaces and some example interfaces.

Currently, RSat and Glucose are included.

.. note::

    Our SAT solver interfaces are 1-based, i.e., literals start at
    1. This is consistent with the popular DIMACS format for SAT
    solving but not with Pythion's 0-based convention. However, this
    also allows to construct clauses using simple integers.

AUTHORS:

- Martin Albrecht (2012): first version
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

import os, sys, subprocess, shlex

from sage.sat.solvers.satsolver import SatSolver
from sage.misc.misc import tmp_filename, get_verbose
from time import sleep

class DIMACS(SatSolver):
    """
    Generic DIMACS Solver.

    .. note::

        Usually, users won't have to use this class directly but some
        class which inherits from thsi class.
    """

    command = ""

    def __init__(self, command=None, filename=None, verbosity=0, **kwds):
        """
        Construct a new generic DIMACS solver.

        INPUT:

        - ``command`` - a named format string with the command to
          run. The string must contain {input} and may contain
          {output} if the solvers writes the solution to an output
          file. For example "sat-solver {input}" is a valid
          command. If ``None`` then the class variable ``command`` is
          used. (default: ``None``)

        - ``filename`` - a filename to write clauses to in DIMACS
          format, must be writable. If ``None`` a temporary filename
          is chosen automatically. (default: ``None``)

        - ``verbosity`` - a verbosity level, where zero means silent
          and anything else means verbose output. (default: ``0``)

        - ``**kwds`` - accepted for compatibility with other solves,
          ignored.

        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS()
            DIMACS Solver: ''
        """
        if filename is None:
            filename = tmp_filename()

        self._headname = filename
        self._verbosity = verbosity

        if command is not None:
            self._command = command
        else:
            self._command = self.__class__.command

        self._tail  = open(tmp_filename(),'w')
        self._var = 0
        self._lit = 0

    def __repr__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: DIMACS(command="iliketurtles {input}")
            DIMACS Solver: 'iliketurtles {input}'
        """
        return "DIMACS Solver: '%s'"%(self._command)

    def __del__(self):
        """
        TESTS::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: d = DIMACS(command="iliketurtles {input}")
            sage: del d
        """
        if not self._tail.closed:
            self._tail.close()
        if os.path.exists(self._tail.name):
            os.unlink(self._tail.name)

    def gen(self, decision=None):
        """
        Return a *new* variable.

        INPUT:

        - ``decision`` - accepted for compatibility with other solvers, ignored.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.gen()
            1
        """
        self._var+= 1
        return self._var

    def ngens(self):
        """
        Return the number of variables.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.gen()
            1
            sage: solver.gen(decision=True)
            2
            sage: solver.ngens()
            2
        """
        return self._var

    def add_clause(self, lits):
        """
        Add a new clause to set of clauses.

        INPUT:

        - ``lits`` - a tuple of integers != 0

        .. note::

            If any element ``e`` in ``lits`` has ``abs(e)`` greater
            than the number of variables generated so far, then new
            variables are created automatically.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: solver = DIMACS()
            sage: solver.gen()
            1
            sage: solver.gen(decision=True)
            2
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver
            DIMACS Solver: ''
        """
        l = []
        for lit in lits:
            lit = int(lit)
            while abs(lit) > self.ngens():
                self.gen()
            l.append(str(lit))
        l.append("0\n")
        self._tail.write(" ".join(l) )
        self._lit += 1

    def write(self):
        """
        Write DIMACS file.

        EXAMPLE::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: _ = solver.write()
            sage: for line in open(fn).readlines():
            ...      print line,
            p cnf 3 1
            1 -2 3 0
        """
        head = open(self._headname, "w")
        head.truncate(0)
        head.write("p cnf %d %d\n"%(self._var,self._lit))
        head.close()

        tail = self._tail
        tail.close()

        head = open(self._headname,"a")
        tail = open(self._tail.name,"r")
        head.write(tail.read())
        tail.close()
        head.close()

        self._tail = open(self._tail.name,"a")
        return self._headname

    def __call__(self, assumptions=None):
        """
        Run 'command' and collect output.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        TESTS:

        This class is not meant to be called directly::

            sage: from sage.sat.solvers.dimacs import DIMACS
            sage: fn = tmp_filename()
            sage: solver = DIMACS(filename=fn)
            sage: solver.add_clause( (1, -2 , 3) )
            sage: solver()
            Traceback (most recent call last):
            ...
            ValueError: No SAT solver command selected.
        """
        if assumptions is not None:
            raise NotImplementedError("Assumptions are not supported for DIMACS based solvers.")

        self.write()

        output_filename = None
        self._output = []

        command = self._command.strip()

        if not command:
            raise ValueError("No SAT solver command selected.")


        if "{output}" in command:
            output_filename = tmp_filename()
        command = command.format(input=self._headname, output=output_filename)

        args = shlex.split(command)

        try:
            process = subprocess.Popen(args, stdout=subprocess.PIPE)
        except OSError:
            raise OSError("Could run '%s', perhaps you need to add your SAT solver to $PATH?"%(" ".join(args)))

        try:
            while process.poll() is None:
                for line in iter(process.stdout.readline,''):
                    if get_verbose() or self._verbosity:
                        print line,
                        sys.stdout.flush()
                    self._output.append(line)
                sleep(0.1)
            if output_filename:
                self._output.extend(open(output_filename).readlines())
        except KeyboardInterrupt:
            process.kill()
            raise KeyboardInterrupt

class RSat(DIMACS):
    """
    An instance of the RSat solver.

    For information on RSat see: http://reasoning.cs.ucla.edu/rsat/
    """

    command = "rsat {input} -v -s"

    def __call__(self, assumptions):
        """
        Solve this instance.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``ngens()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLE::

           sage: from sage.sat.boolean_polynomials import solve as solve_sat
           sage: F,s = mq.SR(1,1,1,4,gf2=True,polybori=True).polynomial_system()
           sage: solve_sat(F, solver=sage.sat.solvers.RSat) # not tested, requires RSat in PATH
        """
        DIMACSSolver.__call__(self)

        for line in self._output:
            if line.startswith("c"):
                continue
            if line.startswith("s"):
                if "UNSAT" in line:
                    return False
            if line.startswith("v"):
                return tuple([0] + map(int, line[2:-2].strip().split(" ")))

class Glucose(DIMACS):
    """
    An instance of the Glucose solver.

    For information on Glucose see: http://www.lri.fr/~simon/?page=glucose
    """

    command = "glucose_static -verb=2 {input} {output}"

    def __call__(self, **kwds):
        """
        Solve this instance.

        INPUT:

        - ``assumptions`` - ignored, accepted for compatibility with
          other solvers (default: ``None``)

        OUTPUT:

        - If this instance is SAT: A tuple of length ``ngens()+1``
          where the ``i``-th entry holds an assignment for the
          ``i``-th variables (the ``0``-th entry is always ``None``).

        - If this instance is UNSAT: ``False``

        EXAMPLE::

           sage: from sage.sat.boolean_polynomials import solve as solve_sat
           sage: F,s = mq.SR(1,1,1,4,gf2=True,polybori=True).polynomial_system()
           sage: solve_sat(F, solver=sage.sat.solvers.Glucose) # not tested, requires Glucose in PATH
        """
        DIMACSSolver.__call__(self)

        for line in self._output:

            if line.startswith("c"):
                continue
            if line.startswith("s"):
                if "UNSAT" in line:
                    return False
            try:
                return tuple([0] + map(int, line[2:-2].strip().split(" ")))
            except ValueError:
                pass
        raise ValueError("Could not parse Glucose output.")
