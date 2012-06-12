"""
CryptoMiniSat.

"CryptoMiniSat is an LGPL-licenced SAT solver that aims to become a
premier SAT solver with all the features and speed of successful SAT
solvers, such as MiniSat and PrecoSat. The long-term goals of
CryptoMiniSat are to be an efficient sequential, parallel and
distributed solver. There are solvers that are good at one or the
other, e.g. ManySat (parallel) or PSolver (distributed), but we wish
to excel at all." -- http://www.msoos.org/cryptominisat2/

AUTHORS:

- Martin Albrecht (2012): first version
"""
##############################################################################
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

include "../../../ext/stdsage.pxi"
include "../../../ext/interrupt.pxi"

from libc.stdint cimport uint32_t
from decl cimport lbool, Var, Lit, Clause, l_Undef, l_False
from decl cimport vec, vector
from decl cimport GaussConf
from solverconf cimport SolverConf

from sage.misc.misc import get_verbose

cdef extern from "cryptominisat_helper.h":
     cdef uint32_t** get_sorted_learnts_helper(Solver* solver, uint32_t* num)

cdef class CryptoMiniSat:
    """
    A CryptoMiniSat SAT-solver instance.

    EXAMPLE::

        sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat # optional - cryptominisat
        sage: CMS = CryptoMiniSat()     # optional - cryptominisat
        sage: CMS.add_clause((1,2,-3))  # optional - cryptominisat
        sage: CMS()                     # optional - cryptominisat
        (True, True, False)

    .. note::

        Do not import 'sage.sat.solvers.cryptominisat.cryptominisat'
        directly, but use 'sage.sat.solvers.cryptominisat' which
        throws a friendlier error message if CryptoMiniSat is not
        installed.
    """
    def __cinit__(self, SolverConf sc=None, **kwds):
        cdef SolverConf _sc
        if sc is not None:
             _sc = sc.__copy__()
        else:
             _sc = SolverConf()
        cdef GaussConf gc

        for k,v in kwds.iteritems():
            _sc[k] = v

        self._solver = new Solver(_sc._conf[0], gc)

    def __dealloc__(self):
        del self._solver

    def __repr__(self):
        s = """CryptoMiniSat instance
#vars: %7d, #lits: %7d, #clauses: %7d, #learnt: %7d, #assigns: %7d
"""%(self._solver.nVars(), self._solver.nLiterals(), self._solver.nClauses(), self._solver.nLearnts(), self._solver.nAssigns())
        return s

    def new_gen(self, decision=None):
        cdef Var var
        if decision is None:
            var = self._solver.newVar()
        else:
            var = self._solver.newVar(bool(decision))
        return int(var+1)

    def ngens(self):
        return int(self._solver.nVars())

    def add_clause(self, lits):
        cdef vec[Lit] l
        for lit in lits:
            while abs(lit) > self._solver.nVars():
                self._solver.newVar()
            l.push(Lit(abs(lit)-1,lit<0))
        self._solver.addClause(l)

    def __call__(self, **kwds):
        _sig_on
        cdef lbool r = self._solver.solve()
        _sig_off

        if r == l_False:
            return False
        if r == l_Undef:
            return None

        return tuple([self._solver.model[i].getBool() for i in range(self._solver.model.size())])

    def unitary_learnt_clauses(self):
        cdef vector[Lit] learnt = self._solver.get_unitary_learnts()

        r = []
        for i in range(learnt.size()):
            r.append( (-1)**learnt[i].sign() * (learnt[i].var()+1) )
        return tuple(r)

    def learnt_clauses(self):
        cdef uint32_t num = 0
        cdef uint32_t **learnt = get_sorted_learnts_helper(self._solver,&num)
        cdef uint32_t *clause = NULL

        r = []
        for i in range(num):
            clause = learnt[i]
            C = [(-1)**int(clause[j]&1) * (int(clause[j]>>1)+1) for j in range(1,clause[0]+1)]
            sage_free(clause)
            r.append(tuple(C))
        sage_free(learnt)
        return r



