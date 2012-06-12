"""
SAT Functions for Boolean Polynomial Systems

AUTHOR:

- Martin Albrecht (2012): initial version

"""
from sage.structure.sequence import Sequence

def solve(F, Converter=None, Solver=None, **kwds):
    """
    Solve ``F`` by solving a SAT-problem -- produced by an instance of ``Converter``
    -- using an instance of ``Solver``.

    If ``Converter`` is ``None`` then :cls:`sage.sat.converters.polybori.CNFEncoder`
    is used.

    If ``Solver`` is ``None`` then
    cls:`sage.sat.solvers.cryptominisat.CryptoMiniSat` is used.

    Parameters can be passed to the converter and the solver by prefixing them with
    'c_' and 's_' respectively. For example, to increase CryptoMiniSat's verbosity
    level, pass 's_verbosity=1'.

    INPUT:

    - ``F`` - a system of Boolean polynomials
    - ``Converter`` - an ANF to CNF converter class (default: ``None``)
    - ``Solver`` - a SAT-solver class
    - **kwds - passed to converter and solver
    """
    try:
        m = len(F)
    except AttributeError:
        F = F.gens()
        m = len(F)

    P = iter(F).next().parent()
    K = P.base_ring()

    # instantiate the SAT solver

    if Solver is None:
        from sage.sat.solvers.cryptominisat import CryptoMiniSat as Solver

    solver_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("s_"):
            solver_kwds[k[2:]] = v

    solver = Solver(**solver_kwds)

    # instantiate the ANF to CNF converter

    if Converter is None:
        from sage.sat.converters.polybori import CNFEncoder as Converter

    converter_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("c_"):
            converter_kwds[k[2:]] = v

    converter = Converter(solver, P, **converter_kwds)

    phi = converter(F, **converter_kwds)

    sol = solver(**solver_kwds)

    if sol:
        return dict((phi[i],K(sol[i])) for i in xrange(P.ngens()))
    else:
        try:
            learnt = solver.unitary_learnt_clauses()
            if learnt:
                return dict((phi[abs(i)-1],K(i<0)) for i in learnt)
            else:
                return sol
        except AttributeError:
            # solver does not support recovering learnt clauses
            return sol

def learn(F, Converter=None, Solver=None, max_length=3, **kwds):
    try:
        m = len(F)
    except AttributeError:
        F = F.gens()
        m = len(F)

    P = iter(F).next().parent()
    K = P.base_ring()

    # instantiate the SAT solver

    if Solver is None:
        from sage.sat.solvers.cryptominisat import CryptoMiniSat as Solver

    solver_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("s_"):
            solver_kwds[k[2:]] = v

    solver = Solver(**solver_kwds)

    # instantiate the ANF to CNF converter

    if Converter is None:
        from sage.sat.converters.polybori import CNFEncoder as Converter

    converter_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("c_"):
            converter_kwds[k[2:]] = v

    converter = Converter(solver, P, **converter_kwds)


    phi = converter(F, **converter_kwds)

    sol = solver(**solver_kwds)

    if sol:
        learnt = ((phi[i]+K(sol[i])) for i in xrange(P.ngens()))
    else:

        learnt = []
        for c in solver.unitary_learnt_clauses():
            try:
                f = phi[abs(i)-1]+K(i<0)
            except KeyError:
                # the solver might have learnt clauses that contain CNF
                # variables which have no correspondence to variables in our
                # polynomial ring (XOR chaining variables for example)
                pass
            learnt.append( f )

        for c in solver.learnt_clauses():
            if len(c) <= max_length:
                try:
                    learnt.append( converter.to_polynomial(c) )
                except ValueError:
                    # the solver might have learnt clauses that contain CNF
                    # variables which have no correspondence to variables in our
                    # polynomial ring (XOR chaining variables for example)
                    pass

    return Sequence(learnt)
