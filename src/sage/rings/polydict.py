"""
PolyDict functionality -- an implementation of the underlying
arithmetic for multi-variate polynomial rings using Python dicts.

This class is not meant for end users, but instead for implementing
multivariate polynomial rings over a completely general base.  It does
not do strong type checking or have parents, etc.  For more special
bases, etc., one would implement something similar in Pyrex for speed.

AUTHOR: William Stein and David Joyner
"""

#*****************************************************************************
#       Copyright (C) 2005 William Stein (wstein@ucsd.edu)
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import copy
import sage.rings.arith as arith
from sage.misc.misc import cputime
import sage.misc.latex as latex

"""
SUMMARY:
    The functions in this file use the "dictionary representation"
    of multivariate polynomials

        {(e1,...,er):c1,...} <-> c1*x1^e1*...*xr^er+...,

    which we call a polydict.   This file implements arithmetic
    functionality for polydicts.

AUTHORS:  William Stein
"""


import sage.rings.ring_element as ring_element

class PolyDict:
    def __init__(self, pdict, zero=0, remove_zero=False, force_int_exponents=True):
        """
        INPUT:
            pdict -- list, which represents a multi-variable polynomial
                     with the distribute representation (a copy is not made)

            zero --  (optional) zero in the base ring

            force_int_exponents -- bool (optional) arithmetic with int exponents is much
                      faster than some of the alternatives, so this is True by default.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: PolyDict({(2,3):2, (1,2):3, (2,1):4})
            PolyDict with representation {(1, 2): 3, (2, 3): 2, (2, 1): 4}

            sage: PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1):4}, force_int_exponents=False)
            PolyDict with representation {(2/3, 3, 5): 2, (1, 2, 1): 3, (2, 1): 4}

            sage: PolyDict({(2,3):0, (1,2):3, (2,1):4}, remove_zero=True)
            PolyDict with representation {(1, 2): 3, (2, 1): 4}
        """
        if not isinstance(pdict, dict):
            if isinstance(pdict, list):
                v = {}
                for w in pdict:
                    if w[0] != 0:
                        v[tuple(w[1])] = w[0]
                remove_zero = False
                pdict = v
            else:
                raise TypeError, "pdict must be a list."
        if force_int_exponents:
            new_pdict = {}
            if remove_zero:
                for k, c in pdict.iteritems():
                    if c != zero:
                        new_pdict[tuple([int(x) for x in k])] = c
            else:
                for k, c in pdict.iteritems():
                    new_pdict[tuple([int(x) for x in k])] = c
            pdict = new_pdict
        else:
            if remove_zero:
                for k in pdict.keys():
                    if pdict[k] == zero:
                        del pdict[k]
        self.__repn  = pdict
        self.__zero = zero

    def __cmp__(self, other):
        if not isinstance(other, PolyDict):
            return False
        return cmp(self.__repn, other.__repn)


    def list(self):
        """
        Return a list that defines self. It is safe to change this.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.list()
            [[3, [1, 2]], [2, [2, 3]], [4, [2, 1]]]
        """
        return [[c,list(e)] for e, c in self.__repn.iteritems()]

    def dict(self):
        """
        Return a copy of the dict that defines self.  It is
        safe to change this.  For a reference, use dictref.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.dict()
            {(1, 2): 3, (2, 3): 2, (2, 1): 4}
        """
        return copy.copy(self.__repn)

    def coefficients(self):
        """
        Return the coefficients of self.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.coefficients()
            [3, 2, 4]
        """
        return self.__repn.values()

    def exponents(self):
        """
        Return the exponents of self.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.exponents()
            [(1, 2), (2, 3), (2, 1)]
        """
        return self.__repn.keys()

    def __len__(self):
        """
        Return the number of terms of the polynomial.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: len(f)
            3
        """
        return len(self.exponents())

    def __getitem__(self, e):
        """
        Return a coefficient of the polynomial.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f[1,2]
            3
            sage: f[(2,1)]
            4
        """
        return self.__repn[e]

    def __repr__(self):
        return 'PolyDict with representation %s'%self.__repn

    def degree(self, x=None):
        if x is None:
            return self.total_degree()
        L = x.__repn.keys()
        if len(L) != 1:
            raise TypeError, "x must be one of the generators of the parent."
        L = L[0]
        nonzero_positions = [i for i in range(len(L)) if L[i] != 0]
        if len(nonzero_positions) != 1:
            raise TypeError, "x must be one of the generators of the parent."
        i = nonzero_positions[0]
        if L[i] != 1:
            raise TypeError, "x must be one of the generators of the parent."
        return max([v[i] for v in self.__repn.keys()])

    def total_degree(self):
        return max([-1] + [sum(e) for e in self.__repn.keys()])

    def monomial_coefficient(self, mon):
        K = mon.keys()[0]
        if not self.__repn.has_key(K):
            return 0
        return self.__repn[K]

    def coefficient(self, mon):
        """
        Return a polydict that defines a polynomial in 1 less number
        of variables that gives the coefficient of mon in this
        polynomial.

        The coefficient is defined as follows.  If f is this
        polynomial, then the coefficient is the sum T/mon where the
        sum is over terms T in f that are exactly divisible by mon.
        """
        K = mon.keys()[0]
        nz = set([i for i in range(len(K)) if K[i] != 0])
        ans = {}
        for S in self.__repn.keys():
            exactly_divides = True
            for j in nz:
                if S[j] != K[j]:
                    exactly_divides = False
                    break
            if exactly_divides:
                t = list(S)
                for m in nz:
                    t[m] = 0
                ans[tuple(t)] = self.__repn[S]
        return PolyDict(ans)



    def is_homogeneous(self):
        K = self.__repn.keys()
        if len(K) == 0:
            return True
        # A polynomial is homogeneous if the number of different
        # exponent sums is at most 1.
        return len(set([sum(e) for e in K])) <= 1

    def homogenize(self):
        R = self.__repn
        H = {}
        deg = self.degree()
        for e, val in R.iteritems():
            i = deg - sum(e)
            f = e + (i,)
            H[f] = val
        return PolyDict(H, zero=self.__zero)

    def latex(self, vars, atomic_exponents=True, atomic_coefficients=True):
        """
        Return a nice polynomial latex representation of this PolyDict, where
        the vars are substituted in.

        INPUT:
            vars -- list
            atomic_exponents -- bool (default: True)
            atomic_coefficients -- bool (default: True)

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.latex(['a','WW'])
            '3 aWW^{2} + 4 a^{2}WW + 2 a^{2}WW^{3}'

        When atomic_exponents is False, the exponents are surrounded
        in parenthesis, since ^ has such high precedence.
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: f.latex(['a','b','c'], atomic_exponents=False)
            '2 a^{2/3}b^{3}c^{5} + 3 ab^{2}c + 4 a^{2}bc'
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        E.sort()
        for e in E:
            c = self.__repn[e]
            if c != 0:
                sign_switch = False
                # First determine the multinomial:
                multi = ""
                for j in range(n):
                    if e[j] != 0:
                        multi += vars[j]
                        if e[j] != 1:
                            multi += "^{%s}"%e[j]
                # Next determine coefficient of multinomial
                if len(multi) == 0:
                    multi = latex.latex(c)
                elif c != 1:
                    if not atomic_coefficients:
                        c = latex.latex(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
                    if len(poly) > 0 and c == -1:
                        sign_switch = True
                    else:
                        multi = "%s %s"%(latex.latex(c),multi)

                # Now add on coefficiented multinomials
                if len(poly) > 0:
                    if sign_switch:
                        poly += " - "
                    else:
                        poly += " + "
                poly += multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -","- ")
        if len(poly) == 0:
            return "0"
        return poly


    def poly_repr(self, vars, atomic_exponents=True, atomic_coefficients=True):
        """
        Return a nice polynomial string representation of this PolyDict, where
        the vars are substituted in.

        INPUT:
            vars -- list
            atomic_exponents -- bool (default: True)
            atomic_coefficients -- bool (default: True)

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f.poly_repr(['a','WW'])
            '3*a*WW^2 + 4*a^2*WW + 2*a^2*WW^3'

        When atomic_exponents is False, the exponents are surrounded
        in parenthesis, since ^ has such high precedence.
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: f.poly_repr(['a','b','c'], atomic_exponents=False)
            '2*a^(2/3)*b^(3)*c^(5) + 3*a*b^(2)*c + 4*a^(2)*b*c'
        """
        n = len(vars)
        poly = ""
        E = self.__repn.keys()
        E.sort()
        for e in E:
            c = self.__repn[e]
            if c != 0:
                sign_switch = False
                # First determine the multinomial:
                multi = ""
                for j in range(n):
                    if e[j] != 0:
                        if len(multi) > 0:
                            multi += "*"
                        multi += vars[j]
                        if e[j] != 1:
                            if atomic_exponents:
                                multi += "^%s"%e[j]
                            else:
                                multi += "^(%s)"%e[j]
                # Next determine coefficient of multinomial
                if len(multi) == 0:
                    multi = str(c)
                elif c != 1:
                    if not atomic_coefficients:
                        c = str(c)
                        if c.find("+") != -1 or c.find("-") != -1 or c.find(" ") != -1:
                            c = "(%s)"%c
                    if len(poly) > 0 and c == -1:
                        sign_switch = True
                    else:
                        multi = "%s*%s"%(c,multi)

                # Now add on coefficiented multinomials
                if len(poly) > 0:
                    if sign_switch:
                        poly += " - "
                    else:
                        poly += " + "
                poly += multi
        poly = poly.lstrip().rstrip()
        poly = poly.replace("+ -","- ")
        if len(poly) == 0:
            return "0"
        return poly


    def __add__(self, other):
        """
        Add two PolyDict's in the same number of variables.

        EXAMPLES:
        We add two polynomials in 2 variables:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f+g
            PolyDict with representation {(1, 2): 3, (1, 5): -3, (1, 1): 3, (2, 1): 4}

        Next we add two polynomials with fractional exponents in 3 variables:
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            sage: f+g
            PolyDict with representation {(2/3, 3, 5): 5, (1, 2, 1): 3, (2, 1, 1): 4}
        """
        if not isinstance(other, PolyDict):
            raise TypeError, "other must be a PolyDict."
        zero = self.__zero
        D = copy.copy(self.__repn)
        R = other.__repn
        for e,c in R.iteritems():
            if D.has_key(e):
                D[e] += c
            else:
                D[e] = c
        F = PolyDict(D, zero=zero, remove_zero=True, force_int_exponents=False)
        return F

    def __mul__(self, right):
        """
        Multiply two PolyDict's in the same number of variables.

        EXAMPLES:
        We multiply two polynomials in 2 variables:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(1,5):-3, (2,3):-2, (1,1):3})
            sage: f*g
            PolyDict with representation {(2, 7): -9, (3, 2): 12, (4, 6): -4, (4, 4): -8, (3, 8): -6, (2, 3): 9, (3, 6): -12, (3, 4): 6, (3, 5): -6}

        Next we multiply two polynomials with fractional exponents in 3 variables:
            sage: f = PolyDict({(2/3,3,5):2, (1,2,1):3, (2,1,1):4}, force_int_exponents=False)
            sage: g = PolyDict({(2/3,3,5):3}, force_int_exponents=False)
            sage: f*g
            PolyDict with representation {(8/3, 4, 6): 12, (4/3, 6, 10): 6, (5/3, 5, 6): 9}

        Finally we print the result in a nice format.
            sage: (f*g).poly_repr(['a','b','c'], atomic_exponents = False)
            '6*a^(4/3)*b^(6)*c^(10) + 9*a^(5/3)*b^(5)*c^(6) + 12*a^(8/3)*b^(4)*c^(6)'
        """
        if not isinstance(right, PolyDict):
            raise TypeError, "other must be a PolyDict."
        newpoly = {}
        k = self.__repn.keys()
        if len(k) == 0:   # product is zero anyways
            return self
        n = len(self.__repn.keys()[0])
        r = range(n)
        for e0, c0 in self.__repn.iteritems():
            for e1, c1 in right.__repn.iteritems():
                e = tuple([e0[i] + e1[i] for i in r])
                c = c0*c1
                if e in newpoly:
                    newpoly[e] += c
                else:
                    newpoly[e] = c
        F = PolyDict(newpoly, self.__zero, force_int_exponents=False, remove_zero=True)
        return F

    def scalar_mult(self, s):
        v = {}
        for e, c in self.__repn.iteritems():
            v[e] = c*s
        return PolyDict(v, self.__zero, force_int_exponents=False)

    def __sub__(self, other):
        """
        Subtract two PolyDict's.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: g = PolyDict({(2,3):2, (1,1):-10})
            sage: f - g
            PolyDict with representation {(1, 2): 3, (1, 1): 10, (2, 1): 4}
            sage: g - f
            PolyDict with representation {(1, 2): -3, (1, 1): -10, (2, 1): -4}
        """

        # TOOD: should refactor add, make abstract operator, so can do both +/-; or copy code.
        if not isinstance(other, PolyDict):
            raise TypeError, "other must be a PolyDict."
        return self + other.scalar_mult(-1)

    def __one(self):
        one = self.__zero + 1
        if len(self.__repn.keys()) == 0:
            v = {(0):one}
        else:
            v = {tuple([0]*len(self.__repn.keys()[0])):one}
        return PolyDict(v, self.__zero, force_int_exponents=False)

    def __pow__(self, n):
        """
        Return the n-th nonnegative power of this PolyDict.

        EXAMPLES:
            sage: from sage.rings.polydict import PolyDict
            sage: f = PolyDict({(2,3):2, (1,2):3, (2,1):4})
            sage: f**2
            PolyDict with representation {(4, 6): 4, (4, 4): 16, (4, 2): 16, (3, 5): 12, (2, 4): 9, (3, 3): 24}
        """
        n = int(n)
        if n < 0:
            raise ValueError, "n must be nonnegative."
        if n == 0:
            return self.__one()
        return arith.generic_power(self, n, self.__one())


