from sage.structure.sequence import Sequence

def solve(F, converter_cls=None, solver_cls=None, **kwds):
    try:
        m = len(F)
    except AttributeError:
        F = F.gens()
        m = len(F)

    P = iter(F).next().parent()
    K = P.base_ring()

    # instantiate the SAT solver

    if solver_cls is None:
        from sage.sat.solvers.cryptominisat.cryptominisat import CryptoMiniSat
        solver_cls = CryptoMiniSat

    solver_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("s_"):
            solver_kwds[k[2:]] = v

    solver = solver_cls(**solver_kwds)

    # instantiate the ANF to CNF converter

    if converter_cls is None:
        from sage.sat.converters.polybori import CNFEncoder
        converter_cls = CNFEncoder

    converter_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("c_"):
            converter_kwds[k[2:]] = v

    converter = converter_cls(solver, P, **converter_kwds)

    phi = converter(F, **converter_kwds)

    sol = solver(**solver_kwds)

    if sol:
        sol = dict((phi[i],K(sol[i])) for i in xrange(P.ngens()))
    else:
        learnt = solver.unitary_learnt_clauses()
        sol = dict((phi[abs(i)-1],K(i<0)) for i in learnt)
    return sol

def learn(F, converter_cls=None, solver_cls=None, max_length=3, **kwds):
    try:
        m = len(F)
    except AttributeError:
        F = F.gens()
        m = len(F)

    P = iter(F).next().parent()
    K = P.base_ring()

    # instantiate the SAT solver

    if solver_cls is None:
        from sage.sat.solvers.cryptominisat.cryptominisat import CryptoMiniSat
        solver_cls = CryptoMiniSat

    solver_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("s_"):
            solver_kwds[k[2:]] = v

    solver = solver_cls(**solver_kwds)

    # instantiate the ANF to CNF converter

    if converter_cls is None:
        from sage.sat.converters.polybori import CNFEncoder
        converter_cls = CNFEncoder

    converter_kwds = {}
    for k,v in kwds.iteritems():
        if k.startswith("c_"):
            converter_kwds[k[2:]] = v

    converter = converter_cls(solver, P, **converter_kwds)


    phi = converter(F, **converter_kwds)

    sol = solver(**solver_kwds)

    if sol:
        learnt = Sequence((phi[i]+K(sol[i])) for i in xrange(P.ngens()))
    else:
        learnt1 = solver.unitary_learnt_clauses()
        learnt1 = Sequence((phi[abs(i)-1]+K(i<0)) for i in learnt1)

        learnt2 = []
        for c in solver.learnt_clauses():
            if len(c) <= max_length:
                try:
                    learnt2.append( converter.to_polynomial(c) )
                except ValueError:
                    pass

        learnt = learnt1 + learnt2
    return Sequence(learnt)
