"""
ANF to CNF Converter based on PolyBoRi's internal representation

AUTHORS:

- Michael Brickenstein (2009) wrote the first version which is
  included in PolyBoRi

- Martin Albrecht (2012) adapted to Sage / updated
"""

##############################################################################
#  Copyright (C) 2009 Michael Brickenstein <brickenstein@mfo.de>
#  Copyright (C) 2012 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from random import Random
from sage.rings.polynomial.pbori import if_then_else as ite

class CNFEncoder(object):
    def __init__(self, R, random_seed = 16, **kwds):
        """
        INPUT:

        - ``R`` - a :cls:`BooleanPolynomialRing`
        - ``random_seed`` - (default: 16)
        - ``**kwds`` - ignored
        """
        self.random_generator = Random(random_seed)
        self.one_set = R.one().set()
        self.empty_set = R.zero().set()
        self.R = R

    def zero_blocks(self, f):
        """
        Divides the zero set of f into blocks.

        EXAMPLE::

            sage: from polybori import *
            sage: r = declare_ring(["x", "y", "z"], dict())
            sage: e = CNFEncoder(r)
            sage: e.zero_blocks(Variable(0)*Variable(1)*Variable(2))
            [{y: 0}, {z: 0}, {x: 0}]
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

    def _clauses(self, f):
        """
        Return clauses representing f as dictionaries of (variable, sign) pairs.

        INPUT:

        - ``f`` - a :cls:BooleanPolynomial

        EXAMPLE::

            sage: B.<x,y,z> = BooleanPolynomialRing()
            sage: from sage.sat.converters.polybori import CNFEncoder
            sage: e = CNFEncoder(R)
            sage: e._clauses(x*y*z)
            [{y: 0, z: 0, x: 0}]
            sage: e._clauses(y + x)
            [{y: 0, x: 1}, {y: 1, x: 0}]
        """
        f_plus_one = f+1
        blocks = self.zero_blocks(f+1)
        negated_blocks = [dict([(variable, 1-value) for (variable, value) in b.iteritems()]) for b in blocks ]

        # we form an expression for a var configuration *not* lying in
        # the block it is evaluated to 0 by f, iff it is not lying in
        # any zero block of f+1

        return negated_blocks

    def clauses(self, f):
        """
        Return clauses representing f as dictionaries of (variable, sign) pairs.

        INPUT:

        - ``f`` - a :cls:BooleanPolynomial

        EXAMPLE::

            sage: B.<x,y,z> = BooleanPolynomialRing()
            sage: from sage.sat.converters.polybori import CNFEncoder
            sage: e = CNFEncoder(R)
            sage: e._clauses(x*y*z)
            [{y: 0, z: 0, x: 0}]
            sage: e._clauses(y + x)
            [{y: 0, x: 1}, {y: 1, x: 0}]
        """
        C = self._clauses(f)

        def to_dimacs_index(v):
            return v.index()+1

        def clause(c):
            return [to_dimacs_index(variable) if value == 1 else -to_dimacs_index(variable) for (variable, value) in c.iteritems()]

        return [clause(c) for c in C]

    def __call__(self, F, satsolver, **kwds):
        """
        Encode the polynomial ``F`` for the SAT-solver ``satsolver''.

        INPUT:

        - ``F`` -
        - ``satsolver`` -
        - ``**kwds`` - ignored

        OUTPUT: An inverse map int -> variable
        """
        res = []
        for f in F:
            clauses = self.clauses(f)
            for c in clauses:
                satsolver.add_clause(c)
        phi = sorted(self.R.gens(), key=lambda x: x.lm().index())
        return phi

    def to_polynomial(self, c):
        def product(l):
            # order of these multiplications for performance
            res = l[0]
            for p in l[1:]:
                res = res*p
            return res
        try:
            return product([variable + value for (variable, value) in c.iteritems()])
        except AttributeError:
            phi = sorted(self.R.gens(), key=lambda x: x.lm().index())

            c = dict((phi[abs(v)-1],int(v>0)) for v in c)
            return product([variable + value for (variable, value) in c.iteritems()])
