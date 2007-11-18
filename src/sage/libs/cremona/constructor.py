def CremonaModularSymbols(level, sign=0, cuspidal=False, verbose=0):
    """
    Return the space of Cremona modular symbols with given level, sign, etc.

    INPUT:
        level -- a positive integer
        sign -- an integer either 0 (the default) or 1.
        cuspidal -- (default: False); if True, compute only the cuspidal subspace
        verbose -- (default: False): if True, print verbose information while creating space

    EXAMPLES:
    """
    from homspace import ModularSymbols
    return ModularSymbols(level=level, sign=sign, cuspidal=cuspidal, verbose=verbose)
