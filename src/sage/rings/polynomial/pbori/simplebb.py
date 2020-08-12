from .PyPolyBoRi import *
from .interred import interred


def buchberger(l):
    "calculates a (non minimal) Groebner basis"
    l = interred(l)
    #for making sure, that every polynomial has a different leading term
    #needed for add_generator
    if not l:
        return []
    g = GroebnerStrategy(l[0].ring())
    for p in l:
        g.add_generator(p)
    while g.npairs() > 0:
        g.clean_top_by_chain_criterion()
        p = g.next_spoly()
        p = g.nf(p)
        if not p.is_zero():
            g.add_generator(p)
    return list(g)


def less_than_n_solutions(ideal, n):
    l = interred(ideal)
    if not l:
        return False
    g = GroebnerStrategy(l[0].ring())
    all_monomials = Monomial([Variable(i) for i
        in range(number_of_variables())]).divisors()
    monomials_not_in_leading_ideal = all_monomials
    for p in l:
        g.add_generator(p)
    while g.npairs() > 0:
        monomials_not_in_leading_ideal = monomials_not_in_leading_ideal \
            % g.reduction_strategy.minimal_leading_terms
        if len(monomials_not_in_leading_ideal) < n:
            return True
        g.clean_top_by_chain_criterion()
        p = g.next_spoly()
        p = g.nf(p)
        if not p.is_zero():
            g.add_generator(p)
    monomials_not_in_leading_ideal = monomials_not_in_leading_ideal \
        % g.reduction_strategy.minimal_leading_terms
    if len(monomials_not_in_leading_ideal) < n:
        return True
    else:
        return False


def gauss(matrix):
    """Toy Gaussian elimination.
    Example: gauss([[0,1],[1,1]]) """
    from .gbcore import groebner_basis

    def get_num(idx, vars):
        if idx in [var.index() for var in vars.variables()]:
            return 1
        return 0

    nrows = len(matrix)
    ncols = len(matrix[0])
    eqs = [sum([matrix[row][col] * Variable(col) for col in range(ncols)])
            for row in range(nrows)]
    result = groebner_basis(eqs)
    result = result + [BooleConstant(0)] * (nrows - len(result))

    return [[get_num(idx, elt.set().vars()) for idx in range(ncols)]
            for elt in result]

    return result
