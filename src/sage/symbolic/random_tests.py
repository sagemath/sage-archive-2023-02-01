from sage.misc.prandom import randint, random
from sage.ext.fast_callable import ExpressionTreeBuilder
import operator
from sage.rings.all import QQ
import sage.calculus.calculus
import sage.symbolic.pynac
from sage.symbolic.constants import *

fast_binary = [(0.4, operator.add), (0.1, operator.sub), (0.5, operator.mul)]
fast_unary = [(0.8, operator.neg), (0.2, operator.abs)]
fast_nodes = [(0.9, fast_binary, 2), (0.1, fast_unary, 1)]

full_binary = [(0.3, operator.add), (0.1, operator.sub), (0.3, operator.mul), (0.2, operator.div), (0.1, operator.pow)]
full_unary = [(0.8, operator.neg), (0.2, operator.inv)]
full_functions = [(1.0, f, f.number_of_arguments()) for f in sage.symbolic.pynac.symbol_table['functions'].values() if f.number_of_arguments() > 0 and 'elliptic' not in str(f) and 'dickman_rho' not in str(f)]
full_nullary = [(1.0, c) for c in [pi, e]] + [(0.05, c) for c in [golden_ratio, log2, euler_gamma, catalan, khinchin, twinprime, merten, brun]]
full_internal = [(0.6, full_binary, 2), (0.2, full_unary, 1), (0.2, full_functions)]

def normalize_prob_list(pl, extra=()):
    if len(pl) == 0:
        return pl
    result = []
    total = sum([float(p[0]) for p in pl])
    for p in pl:
        prob = p[0]
        val = p[1]
        if len(p) > 2:
            p_extra = p[2:]
        else:
            p_extra = extra
        if isinstance(val, list):
            norm_val = normalize_prob_list(val, extra=p_extra)
            for np in norm_val:
                result.append(((prob/total)*np[0], np[1]) + np[2:])
        else:
            result.append(((prob/total), val) + p_extra)
    return result

def choose_from_prob_list(lst):
    r = random()
    for i in range(len(lst)-1):
        if r < lst[i][0]:
            return lst[i]
        r -= lst[i][0]
    return lst[-1]

def random_integer_vector(n, length):
    r"""
    This is an approximation to IntegerVectors(n, length).random_element().
    That gives values uniformly at random, but might be slow; this
    routine is not uniform, but should always be fast.
    """
    if length == 0:
        return []
    elif length == 1:
        return [n]
    elif length == 2:
        v = randint(0, n)
        return [v, n-v]
    else:
        v = randint(0, 2*n//length)
        return [v] + random_integer_vector(n-v, length-1)

def random_expr_helper(n_nodes, internal, leaves, verbose):
    if n_nodes == 1:
        return choose_from_prob_list(leaves)[1]
    else:
        r = choose_from_prob_list(internal)
        n_nodes -= 1
        n_children = r[2]
        n_spare_nodes = n_nodes - n_children
        if n_spare_nodes <= 0:
            n_spare_nodes = 0
        nodes_per_child = random_integer_vector(n_spare_nodes, n_children)
        children = [random_expr_helper(n+1, internal, leaves, verbose) for n in nodes_per_child]
        if verbose:
            print "About to apply %r to %r" % (r[1], children)
        return r[1](*children)

def random_expr(size, nvars=1, nconsts=None, var_frac=0.5, internal=full_internal, nullary=full_nullary, nullary_frac=0.2, const_generator=QQ.random_element, verbose=False):
    vars = [(1.0, sage.calculus.calculus.var('v%d' % (n+1))) for n in range(nvars)]
    if nconsts is None:
        nconsts = size
    consts = [(1.0, const_generator()) for _ in range(nconsts)]
    leaves = [(var_frac, vars), (1.0 - var_frac - nullary_frac, consts), (nullary_frac, nullary)]
    leaves = normalize_prob_list(leaves)

    internal = normalize_prob_list(internal)

    return random_expr_helper(size, internal, leaves, verbose)
