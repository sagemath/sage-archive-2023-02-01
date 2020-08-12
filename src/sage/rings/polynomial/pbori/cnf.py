from random import Random
from .PyPolyBoRi import (Monomial, BooleSet, Polynomial, if_then_else
    as ite, lp, gauss_on_polys, ll_red_nf_redsb)
from .ll import ll_encode
from .statistics import used_vars_set


class CNFEncoder(object):
    def __init__(self, r, random_seed=16):
        self.random_generator = Random(random_seed)
        self.one_set = r.one().set()
        self.empty_set = r.zero().set()
        self.r = r

    def zero_blocks(self, f):
        """divides the zero set of f into blocks
        >>> from brial import *
        >>> r = declare_ring(["x", "y", "z"], dict())
        >>> e = CNFEncoder(r)
        >>> e.zero_blocks(r.variable(0)*r.variable(1)*r.variable(2))
        [{y: 0}, {z: 0}, {x: 0}]
        """
        f = Polynomial(f)
        variables = f.vars_as_monomial()

        space = variables.divisors()
        variables = list(variables.variables())
        zeros = f.zeros_in(space)
        rest = zeros
        res = list()

        def choose_old(s):
            return next(iter(rest))  # somewhat

            #inefficient compared to polynomials lex_lead
        def choose(s):
            indices = []
            assert not s.empty()
            nav = s.navigation()
            while not nav.constant():
                e = nav.else_branch()
                t = nav.then_branch()
                if e.constant() and not e.terminal_one():
                    indices.append(nav.value())
                    nav = t
                else:
                    if self.random_generator.randint(0, 1):
                        indices.append(nav.value())
                        nav = t

                    else:
                        nav = e
            assert nav.terminal_one()
            res = self.one_set
            for i in reversed(indices):
                res = ite(i, res, self.empty_set)
            return next(iter(res))
        while not rest.empty():
            l = choose(rest)
            l_variables = set(l.variables())

            def get_val(var):
                if var in l_variables:
                    return 1
                return 0
            block_dict = dict([(v, get_val(v)) for v in variables])

            l = l.set()
            self.random_generator.shuffle(variables)
            for v in variables:
                candidate = l.change(v.index())
                if candidate.diff(zeros).empty():
                    l = l.union(candidate)
                    del block_dict[v]
            rest = rest.diff(l)
            res.append(block_dict)
        return res

    def clauses(self, f):
        """
        >>> from brial import *
        >>> r = declare_ring(["x", "y", "z"], dict())
        >>> e = CNFEncoder(r)
        >>> e.clauses(r.variable(0)*r.variable(1)*r.variable(2)) # doctest:+ELLIPSIS
        [{...x: 0...}]
        >>> e.clauses(r.variable(1)+r.variable(0)) # doctest:+ELLIPSIS
        [{...x: 1...}, {...y: 1...}]

        >>> [sorted(c.iteritems()) for c in e.clauses(r.variable(0)*r.variable(1)*r.variable(2))]
        [[(z, 0), (y, 0), (x, 0)]]
        >>> [sorted(c.iteritems()) for c in e.clauses(r.variable(1)+r.variable(0))]
        [[(y, 1), (x, 0)], [(y, 0), (x, 1)]]
        """
        f_plus_one = f + 1
        blocks = self.zero_blocks(f + 1)
        negated_blocks = [dict([(variable, 1 - value) for (variable, value)
            in b.items()]) for b in blocks]
        # we form an expression for a var configuration *not* lying in the
        # block it is evaluated to 0 by f, iff it is not lying in any zero
        # block of f+1
        return negated_blocks

    def polynomial_clauses(self, f):
        """
        >>> from brial import *
        >>> r = declare_ring(["x", "y", "z"], dict())
        >>> e = CNFEncoder(r)
        >>> e.polynomial_clauses(r.variable(0)*r.variable(1)*r.variable(2))
        [x*y*z]
        >>> v = r.variable
        >>> p = v(1)*v(2)+v(2)*v(0)+1
        >>> groebner_basis([p], heuristic = False)==groebner_basis(e.polynomial_clauses(p), heuristic = False)
        True
        """

        def product(l):
            res = l[0]
            for p in l[1:]:
                res = res * p
            # please care about the order of these multiplications for
            # performance
            return res
        return [product([variable + value for (variable, value)
                         in b.items()]) for b in self.clauses(f)]

    def to_dimacs_index(self, v):
        return v.index() + 1

    def dimacs_encode_clause(self, c):
        def get_sign(val):
            if value == 1:
                return 1
            return -1

        items = sorted(c.items(), reverse=True)
        return " ".join(
        [str(v) for v in
            [
            get_sign(value) * self.to_dimacs_index(variable)
            for (variable, value) in items] + [0]])

    def dimacs_encode_polynomial(self, p):
        """
         >>> from brial import *
         >>> d=dict()
         >>> r = declare_ring(["x", "y", "z"], d)
         >>> e = CNFEncoder(r)
         >>> e.dimacs_encode_polynomial(d["x"]+d["y"]+d["z"])
         ['1 2 -3 0', '1 -2 3 0', '-1 -2 -3 0', '-1 2 3 0']
            """
        clauses = self.clauses(p)
        res = []
        for c in clauses:
            res.append(self.dimacs_encode_clause(c))
        return res

    def dimacs_cnf(self, polynomial_system):
        r"""
        >>> from brial import *
        >>> r = declare_ring(["x", "y", "z"], dict())
        >>> e = CNFEncoder(r)
        >>> e.dimacs_cnf([r.variable(0)*r.variable(1)*r.variable(2)])
        'c cnf generated by PolyBoRi\np cnf 3 1\n-1 -2 -3 0'
        >>> e.dimacs_cnf([r.variable(1)+r.variable(0)])
        'c cnf generated by PolyBoRi\np cnf 3 2\n1 -2 0\n-1 2 0'
        >>> e.dimacs_cnf([r.variable(0)*r.variable(1)*r.variable(2), r.variable(1)+r.variable(0)])
        'c cnf generated by PolyBoRi\np cnf 3 3\n-1 -2 -3 0\n-1 2 0\n1 -2 0'
        """
        clauses_list = [c for p in polynomial_system for c in self.
            dimacs_encode_polynomial(p)]
        res = ["c cnf generated by PolyBoRi"]
        r = polynomial_system[0].ring()
        n_variables = r.n_variables()
        res.append("p cnf %s %s" % (n_variables, len(clauses_list)))
        for c in clauses_list:
            res.append(c)
        return "\n".join(res)


class CryptoMiniSatEncoder(CNFEncoder):
    group_counter = 0

    def dimacs_encode_polynomial(self, p):
        r"""
         >>> from brial import *
         >>> d=dict()
         >>> r = declare_ring(["x", "y", "z"], d)
         >>> e = CryptoMiniSatEncoder(r)
         >>> p = d["x"]+d["y"]+d["z"]
         >>> p.deg()
         1
         >>> len(p)
         3
         >>> e.dimacs_encode_polynomial(p)
         ['x1 2 3 0\nc g 1 x + y + z']
         >>> e.dimacs_encode_polynomial(p+1)
         ['x1 2 -3 0\nc g 2 x + y + z + 1']
            """
        if p.deg() != 1 or len(p) <= 1:
            res = super(CryptoMiniSatEncoder, self).dimacs_encode_polynomial(p)
        else:

            if p.has_constant_part():
                invert_last = True
            else:
                invert_last = False
            variables = list(p.vars_as_monomial().variables())
            indices = [self.to_dimacs_index(v) for v in variables]
            if invert_last:
                indices[-1] = -indices[-1]
            indices.append(0)
            res = ["x" + " ".join([str(v) for v in indices])]
        self.group_counter = self.group_counter + 1
        group_comment = "\nc g %s %s" % (self.group_counter, str(p)[:30])
        return [c + group_comment for c in res]

    def dimacs_cnf(self, polynomial_system):
        r"""
            >>> from brial import *
            >>> r = declare_ring(["x", "y", "z"], dict())
            >>> e = CryptoMiniSatEncoder(r)
            >>> e.dimacs_cnf([r.variable(0)*r.variable(1)*r.variable(2)])
            'c cnf generated by PolyBoRi\np cnf 3 1\n-1 -2 -3 0\nc g 1 x*y*z\nc v 1 x\nc v 2 y\nc v 3 z'
            >>> e.dimacs_cnf([r.variable(1)+r.variable(0)])
            'c cnf generated by PolyBoRi\np cnf 3 1\nx1 2 0\nc g 2 x + y\nc v 1 x\nc v 2 y'
            >>> e.dimacs_cnf([r.variable(0)*r.variable(1)*r.variable(2), r.variable(1)+r.variable(0)])
            'c cnf generated by PolyBoRi\np cnf 3 2\n-1 -2 -3 0\nc g 3 x*y*z\nx1 2 0\nc g 4 x + y\nc v 1 x\nc v 2 y\nc v 3 z'
            """
        uv = list(used_vars_set(polynomial_system).variables())
        res = super(CryptoMiniSatEncoder, self).dimacs_cnf(polynomial_system)
        res = res + "\n" + "\n".join(["c v %s %s" % (self.to_dimacs_index(v),
            v) for v in uv])
        return res


def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()
