"""
Testing modular symbols spaces.
"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@ucsd.edu>
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


import random

import modsym
import sage.modular.dirichlet as dirichlet
import sage.modular.congroup as congroup
from sage.misc.misc import cputime

class Test:
    """
    Modular symbol testing class.
    """
    def __init__(self, levels=20, weights=4, onlyg0=False, onlyg1=False, onlychar=False):
        """
        INPUT:
            levels --  list or int
            weights -- list or int
            onlyg0 -- bool, if True only select Gamma0 spaces for testing
            onlyg1 -- bool, if True only select Gamma1 spaces for testing
            onlychar -- bool, if True only selects spaces with character for testing
        """
        if not isinstance(levels, list):
            levels = range(1,int(levels)+1)
        if not isinstance(weights, list):
            weights = range(2,int(weights)+1)
        self.levels = levels
        self.weights = weights
        if len(levels) < 1:
            raise RuntimeError, "levels must have positive length"
        if len(weights) < 1:
            raise RuntimeError, "weights must have positive length"
        self.current_space = None
        self.onlyg0 = onlyg0
        self.onlyg1 = onlyg1
        self.onlychar = onlychar

    def __repr__(self):
        return "Modular symbols testing class"

    def _modular_symbols_space(self):
        if self.onlyg0:
            which = 0
        elif self.onlyg1:
            which = 1
        elif self.onlychar:
            which = 2
        else:
            which = random.randrange(0,3)
        if which == 0:
            print "gamma0"
            M = self._modular_symbols_space_gamma0()
        elif which == 1:
            print "gamma1"
            M = self._modular_symbols_space_gamma1()
        else:
            print "character"
            M = self._modular_symbols_space_character()
        print "\t",M
        return M

    def _level_weight_sign(self):
        level = random.choice(self.levels)
        weight = random.choice(self.weights)
        sign = random.choice([-1,0,1])
        print "level = %s, weight = %s, sign = %s"%(level,weight,sign)
        return level, weight, sign

    def _modular_symbols_space_gamma0(self):
        level, weight, sign = self._level_weight_sign()
        M = modsym.ModularSymbols(congroup.Gamma0(level), weight, sign)
        self.current_space = M
        return M

    def _modular_symbols_space_gamma1(self):
        level, weight, sign = self._level_weight_sign()
        M = modsym.ModularSymbols(congroup.Gamma1(level), weight, sign)
        self.current_space = M
        return M

    def _modular_symbols_space_character(self):
        level, weight, sign = self._level_weight_sign()
        G = dirichlet.DirichletGroup(level)
        eps = G.random_element()
        M = modsym.ModularSymbols(eps, weight, sign)
        self.current_space = M
        return M

    def _do(self, name):
        print "test_%s"%name
        Test.__dict__["test_%s"%name](self)

    #################################################################
    # The tests
    #################################################################
    def random(self, seconds=0):
        self.test("random", seconds)

    def test(self, name, seconds=0):
        seconds = float(seconds)
        total = cputime()
        n = 1
        while seconds == 0 or cputime(total) < seconds:
            s = "** test_dimension: number %s"%n
            if seconds > 0:
                s += " (will stop after about %s seconds)"%seconds
            t = cputime()
            self._do(name)
            print "\ttime=%s\telapsed=%s"%(cputime(t),cputime(total))
            n += 1

    def test_cs_dimension(self):
        """
        Compute the cuspidal subspace (this implicitly checks that
        the dimension is correct using formulas).
        """
        self._modular_symbols_space().cuspidal_submodule()

    def test_csnew_dimension(self):
        """
        Compute the new cuspidal subspace and verify that the dimension is correct using
        a dimension formula.
        """
        M = self._modular_symbols_space()
        V = M.cuspidal_submodule().new_submodule()
        d = V.dimension()
        d2 = M._cuspidal_new_submodule_dimension_formula()
        if d != d2:
            assert False, "Test failed for M=\"%s\", where computed dimension is %s but formula dimension is %s."%(
                     M, d, d2)

    def test_csns_nscs(self):
        """
        Compute new cuspidal subspace in two ways and verify that the results are the same.
        """
        M = self._modular_symbols_space()
        V1 = M.cuspidal_submodule().new_submodule()
        V2 = M.new_submodule().cuspidal_submodule()
        assert V1 == V2, "Test failed for M=\"%s\", where the new cuspidal and cuspidal new spaces are computed differently."%M
        d = M._cuspidal_new_submodule_dimension_formula()
        if d != V1.dimension():
            assert False, "Test failed for M=\"%s\", where computed dimension is %s but formula dimension is %s."%(
                     M, V1.dimension(), d)


    def test_decomposition(self):
        M = self._modular_symbols_space()
        D = M.decomposition()
        assert M.dimension() == sum([A.dimension() for A in D])

    def test_dimension(self):
        self._modular_symbols_space().dimension()


    def test_random(self):
        """
        Do a random test from all the possible tests.
        """
        tests = [a for a in Test.__dict__.keys() if a[:5] == "test_" and a != "test_random"]
        name = random.choice(tests)
        print "Doing random test %s"%name
        Test.__dict__[name](self)

