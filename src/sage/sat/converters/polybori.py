"""
ANF to CNF Converter.

This converter is based on PolyBoRi's internal representation for
polynomials with few variables and XOR chain splitting for polynomials
with many variables.

AUTHORS:

- Martin Albrecht - (2008-09) initial version of 'anf2cnf.py'
- Michael Brickenstein - (2009) 'cnf.py' for PolyBoRi
- Mate Soos - (2010) improved version of 'anf2cnf.py'
- Martin Albrecht - (2012) unified and added to Sage
"""

##############################################################################
#  Copyright (C) 2008-2009 Martin Albrecht <martinralbrecht@googlemail.com>
#  Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#  Copyright (C) 2010 Mate Soos
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from random import Random
from sage.rings.polynomial.pbori import if_then_else as ite
from sage.rings.integer_ring import ZZ
from sage.functions.other import ceil
from sage.misc.cachefunc import cached_method, cached_function
from sage.combinat.permutation import Permutations

class CNFEncoder(object):
    def __init__(self, sat_solver, ring, max_variables=6, use_xor_clauses=None, cutting_number=4, random_seed=16, **kwds):
        """
        INPUT:

        - ``sat_solver`` - a SAT-solver instance
        - ``ring`` - a :cls:`BooleanPolynomialRing`
        - ``max_variables`` - maximum number of variables for direct conversion
        - ``use_xor_clauses`` - use XOR clauses; if ``None`` use if ``sat_solver`` supports it. (default: ``None``)
        - ``cutting_number`` - maximum length of XOR chains
        - ``random_seed`` - (default: 16)
        - ``**kwds`` - ignored
        """
        self.random_generator = Random(random_seed)
        self.one_set = ring.one().set()
        self.empty_set = ring.zero().set()

        self.sat_solver = sat_solver
        self.max_variables = max_variables
        self.cutting_number = cutting_number

        if use_xor_clauses is None:
            use_xor_clauses = hasattr(sat_solver,"add_xor_clause")
        self.use_xor_clauses = use_xor_clauses

        self.ring = ring

        self._phi = [None]
        for x in sorted([x.lm() for x in self.ring.gens()], key=lambda x: x.index()):
            self.gen(x, decision=True)

    def gen(self, m, decision=None):
        self._phi.append(m)
        return self.sat_solver.gen(decision=True)

    def phi(self):
        """
        Map CNF variables to ANF variables.
        """
        return list(self._phi)

    ##################################################
    # Encoding based on roots of polynomial
    ##################################################

    def zero_blocks(self, f):
        """
        Divides the zero set of f into blocks.

        EXAMPLE::

            sage: B.<a,b,c> = BooleanPolynomialRing()
            sage: from sage.sat.converters.polybori import CNFEncoder
            sage: from sage.sat.solvers.cryptominisat import CryptoMiniSat
            sage: cms = CryptoMiniSat()
            sage: e = CNFEncoder(cms, B)
            sage: sorted(e.zero_blocks(a*b*c))
            [{c: 0}, {b: 0}, {a: 0}]
        """
        variables = f.vars_as_monomial()

        space = variables.divisors()
        variables = list(variables.variables())
        zeros = f.zeros_in(space)
        rest = zeros
        res = list()

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
                    if self.random_generator.randint(0,1):
                        indices.append(nav.value())
                        nav = t

                    else:
                        nav = e
            assert nav.terminal_one()
            res = self.one_set
            for i in reversed(indices):
                res = ite(i, res, self.empty_set)
            return iter(res).next()

        while not rest.empty():
            l = choose(rest)
            l_variables = set(l.variables())
            block_dict = dict([(v, 1 if v in l_variables else 0) for v in variables])
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

    def clauses_sparse(self, f):
        """
        INPUT:

        - ``f`` - a :cls:`BooleanPolynomial`
        """
        # we form an expression for a var configuration *not* lying in
        # the block it is evaluated to 0 by f, iff it is not lying in
        # any zero block of f+1

        blocks = self.zero_blocks(f+1)
        C = [dict([(variable, 1-value) for (variable, value) in b.iteritems()]) for b in blocks ]

        def to_dimacs_index(v):
            return v.index()+1

        def clause(c):
            return [to_dimacs_index(variable) if value == 1 else -to_dimacs_index(variable) for (variable, value) in c.iteritems()]

        for c in C:
            self.sat_solver.add_clause(clause(c))

    ###
    # Indirect conversion, may add new variables
    ###

    def clauses_dense(self, f):
        equal_zero = not bool(f.constant_coefficient())

        f = (f - f.constant_coefficient())
        f = [self.monomial(m) for m in f]

        if self.use_xor_clauses:
            self.xor_cnf(f, equal_zero)
        elif f > self.cutting_number:
            for fpart, this_equal_zero in self.split_xor(f, equal_zero):
                self.xor_cnf(fpart, this_equal_zero)
        else:
            self.xor_cnf(f, equal_zero)

    @cached_method
    def monomial(self, m):
        if m.deg() == 1:
            return m.index()+1
        else:
            # we need to encode the relationship between the monomial
            # and its variables
            variables = [self.monomial(v) for v in m.variables()]
            monomial = self.gen(m, decision=False)

            # (a | -w) & (b | -w) & (w | -a | -b) <=> w == a*b
            for v in variables:
                self.sat_solver.add_clause( (v, -monomial) )
            self.sat_solver.add_clause( tuple([monomial] + [-v for v in variables]) )

            return monomial

    @cached_function
    def permutations(length, equal_zero):
        E = []
        for num_negated in range(0, length+1) :
            if (((num_negated % 2) ^ ((length+1) % 2)) == equal_zero) :
                continue
            start = []
            for i in range(num_negated) :
                start.append(1)
            for i in range(length - num_negated) :
                start.append(-1)
            E.extend(Permutations(start))
        return E

    def split_xor(self, monomial_list, equal_zero):
        c = self.cutting_number

        nm = len(monomial_list)
        step = ceil((c-2)/ZZ(nm) * nm)
        M = []

        new_variables = []
        for j in range(0, nm, step):
            m =  new_variables + monomial_list[j:j+step]
            if (j + step) < nm:
                new_variables = [self.gen(None, decision=False)]
                m += new_variables
            M.append([m, equal_zero])
            equal_zero = True
        return M

    def xor_cnf(self, M, equal_zero):
        if self.use_xor_clauses:
            self.sat_solver.add_xor_clause(M, equal_zero)
        else:
            ll = len(M)
            for p in self.permutations(ll, equal_zero):
                self.sat_solver.add_clause([ p[i]*M[i] for i in range(ll) ])

    ####################################################
    # Highlevel Functions
    ###################################################
    def clauses(self, f):
        if f.nvariables() <= self.max_variables:
            self.clauses_sparse(f)
        else:
            self.clauses_dense(f)

    def __call__(self, F, **kwds):
        """
        Encode the polynomial ``F`` .

        INPUT:

        - ``F`` -
        - ``**kwds`` - ignored

        OUTPUT: An inverse map int -> variable
        """
        res = []
        for f in F:
            self.clauses(f)
        return self.phi()

    ##
    # The way back
    ##

    def to_polynomial(self, c):
        def product(l):
            # order of these multiplications for performance
            res = l[0]
            for p in l[1:]:
                res = res*p
            return res

        phi = self.phi()
        product = self.ring(1)
        for v in c:
            if phi[abs(v)] is None:
                raise ValueError("Clause containst an XOR glueing variable.")
            product *= phi[abs(v)] + int(v>0)
        return product
