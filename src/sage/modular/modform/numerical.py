"""
Numerical computation of newforms
"""

#########################################################################
#       Copyright (C) 2004--2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.structure.sage_object import SageObject
from sage.structure.sequence    import Sequence
from sage.modular.modsym.all    import ModularSymbols
from sage.modular.congroup      import is_CongruenceSubgroup, Gamma0
#from sage.modular.dims          import dimension_new_cusp_forms
from sage.modules.all           import vector
from sage.misc.misc             import verbose
from sage.rings.all             import CDF, Integer, QQ, next_prime, prime_range
from random                     import randint

class NumericalNewforms(SageObject):
    def __init__(self, group, weight=2, eps=1e-20, delta=1e-2, numtp=2):
        if isinstance(group, (int, long, Integer)):
            group = Gamma0(Integer(group))
        self._group  = group
        self._weight = Integer(weight)
        self._numtp = numtp
        if self._weight < 2:
            raise ValueError, "weight must be at least 2"
        self._eps = eps
        self._delta = delta

    def level(self):
        return self._group.level()

    def weight(self):
        return self._weight

    def _repr_(self):
        return "Numerical Hecke eigenvalues for %s of weight %s"%(
            self._group, self._weight)

    def modular_symbols(self):
        try:
            return self.__modular_symbols
        except AttributeError:
            M = ModularSymbols(self._group,
                    self._weight, sign=1)
            if M.base_ring() != QQ:
                raise ValueError, "modular forms space must be defined over QQ"
            self.__modular_symbols = M
            return M

##     def __len__(self):
##         try:
##             return self.__len
##         except AttributeError:
##             pass
##         n = dimension_new_cusp_forms(self._group, self._weight)
##         self.__len = n
##         return n

    def _xxx_eigenvectors(self):
        try:
            return self.__eigenvectors
        except AttributeError:
            pass
        verbose('Finding eigenvector basis')
        M = self.modular_symbols().cuspidal_submodule()
        M_amb = self.modular_symbols()
        N = self.level()
        n = len(self)

        def next_p(p):
            p = next_prime(p)
            while N%p == 0:
                p = next_prime(p)
            return p

        p = next_p(1)

        t = M.T(p).matrix()
        t_amb = M_amb.T(p).matrix()
        eps = self._eps

        w = []
        while len(w) < n:
            tm = verbose('computing complex eigenvectors of %s x %s matrix; found %s of %s so far'%(t.nrows(),t.nrows(), len(w), n))
            evals, B = t.change_ring(CDF).eigen_left()
            verbose('done computing eigenvectors', t=tm)
            # Find the eigenvalues that occur with multiplicity 1 up
            # to the given eps.
            v = list(evals)
            v.sort()
            w = []
            i = 0
            while i < len(v):
                e = v[i]
                if i + 1 == len(v) or abs(e-v[i+1]) >= eps:
                    w.append((e,i))
                i += 1
            if len(w) < n:
                p = next_p(p)
                tt = M.T(p).matrix()
                tt_amb = M_amb.T(p).matrix()
                c = randint(-100,100)
                if c > 1:
                    tt = c*tt
                    tt_amb = c*tt_amb
                t += tt
                t_amb += tt_amb
            #endif
        #end while

        evals, B = t_amb.change_ring(CDF).eigen_left()
        z = []
        # match up the eigenvalues
        for e, _ in w:
            j = -1
            for i in range(len(evals)):
                if abs(e - evals[i]) < eps:
                    print e, evals[i], e-evals[i]
                    if j != -1:
                        raise RuntimeError, "Precision error -- try reducing eps from %s to something smaller"%eps
                    j = i
            if j == -1:
                raise RuntimeError, "Precision error -- try reducing eps from %s to something larger"%eps
            z.append(j)
        self.__eigenvalues = evals
        self.__eigenvectors = B.matrix_from_columns(z)
        return self.__eigenvectors

    def _eigenvectors(self):
        try:
            return self.__eigenvectors
        except AttributeError:
            pass
        verbose('Finding eigenvector basis')
        M = self.modular_symbols()
        N = self.level()

        p = 2
        t = M.T(p).matrix()
        for i in range(self._numtp-1):
            p = next_prime(p)
            t += randint(-50,50)*M.T(p).matrix()

        self._hecke_matrix = t
        evals, B = t.change_ring(CDF).eigen_left()

        # Find the eigenvalues that occur with multiplicity 1 up
        # to the given eps.
        eps = self._eps
        v = list(evals)
        v.sort()
        w = []
        for i in range(len(v)):
            e = v[i]
            uniq = True
            for j in range(len(v)):
                if i != j and abs(e-v[j]) < eps:
                    uniq = False
            if uniq:
                w.append(i)
        self.__eigenvectors = B.matrix_from_columns(w)
        return B

    def _xxx_eigenvectors(self):
        try:
            return self.__eigenvectors
        except AttributeError:
            pass
        verbose('Finding eigenvector basis')
        M = self.modular_symbols().cuspidal_submodule()
        M_amb = self.modular_symbols()
        N = self.level()
        n = len(self)

        def next_p(p):
            p = next_prime(p)
            while N%p == 0:
                p = next_prime(p)
            return p

        p = next_p(1)

        t = M.T(p).matrix()
        t_amb = M_amb.T(p).matrix()
        eps = self._eps

        w = []
        while len(w) < n:
            tm = verbose('computing complex eigenvectors of %s x %s matrix; found %s of %s so far'%(t.nrows(),t.nrows(), len(w), n))
            evals, B = t.change_ring(CDF).eigen_left()
            verbose('done computing eigenvectors', t=tm)
            # Find the eigenvalues that occur with multiplicity 1 up
            # to the given eps.
            v = list(evals)
            v.sort()
            w = []
            i = 0
            while i < len(v):
                e = v[i]
                if i + 1 == len(v) or abs(e-v[i+1]) >= eps:
                    w.append((e,i))
                i += 1
            if len(w) < n:
                p = next_p(p)
                tt = M.T(p).matrix()
                tt_amb = M_amb.T(p).matrix()
                c = randint(-100,100)
                if c > 1:
                    tt = c*tt
                    tt_amb = c*tt_amb
                t += tt
                t_amb += tt_amb
            #endif
        #end while

        evals, B = t_amb.change_ring(CDF).eigen_left()
        z = []
        # match up the eigenvalues
        for e, _ in w:
            j = -1
            for i in range(len(evals)):
                if abs(e - evals[i]) < eps:
                    print e, evals[i], e-evals[i]
                    if j != -1:
                        raise RuntimeError, "Precision error -- try reducing eps from %s to something smaller"%eps
                    j = i
            if j == -1:
                raise RuntimeError, "Precision error -- try reducing eps from %s to something larger"%eps
            z.append(j)
        self.__eigenvalues = evals
        self.__eigenvectors = B.matrix_from_columns(z)
        return self.__eigenvectors

    def _easy_vector(self):
        """
        Return a very sparse vector v such that v times the eigenvector matrix
        has all entries nonzero.

        ALGORITHM:
           1. Choose row with the most nonzero entries.   (put 1 there)
           2. Consider submatrix of columns corresponding
              to zero entries in row chosen in 1.
           3. Find row of submatrix with most nonzero entries,
              and add appropriate multiple.  Repeat.

        EXAMPLES:
            sage: n = numerical_newforms(37)
            sage: n._easy_vector()
            (1.0, 1.0)
            sage: n = numerical_newforms(43)
            sage: n._easy_vector()
            (0, 1.0, 0)
            sage: n = numerical_newforms(125)
            sage: n._easy_vector()
            (1.0, 0, 0, 0, 0, 0, 0, 0)
        """
        try:
            return self.__easy_vector
        except AttributeError:
            pass
        E = self._eigenvectors()
        delta = self._delta
        x = (CDF**E.nrows()).zero_vector()
        if E.nrows() == 0:
            return x



        def best_row(M):
            R = M.rows()
            v = [len(support(r, delta)) for r in R]
            m = max(v)
            i = v.index(m)
            return i, R[i]

        i, e = best_row(E)

        x[i] = 1

        while True:
            s = set(support(e, delta))
            zp = [i for i in range(e.degree()) if not i in s]
            if len(zp) == 0:
                break
            C = E.matrix_from_columns(zp)
            # best row
            i, f = best_row(C)
            x[i] += 1   # simplistic
            e = x*E

        self.__easy_vector = x
        return x

    def _eigendata(self):
        try:
            return self.__eigendata
        except AttributeError:
            pass
        x = self._easy_vector()

        B = self._eigenvectors()
        def phi(y):
            return y.element() * B

        phi_x = phi(x)
        V = phi_x.parent()
        phi_x_inv = V([a**(-1) for a in phi_x])
        eps = self._eps
        nzp = support(x, eps)
        x_nzp = vector(CDF, x.list_from_positions(nzp))
        self.__eigendata = (phi_x, phi_x_inv, nzp, x_nzp)
        return self.__eigendata

    def ap(self, p):
        p = Integer(p)
        try:
            return self._ap[p]
        except AttributeError:
            self._ap = {}
        except KeyError:
            pass
        a = self.eigenvalues([p])[0]
        self._ap[p] = a
        return a

    def eigenvalues(self, primes):
        if isinstance(primes, (int, long, Integer)):
            primes = [Integer(primes)]
        phi_x, phi_x_inv, nzp, x_nzp = self._eigendata()
        B = self._eigenvectors()
        def phi(y):
            return y.element() * B

        ans = []
        m = self.modular_symbols().ambient_module()
        for p in primes:
            t = m._compute_hecke_matrix_prime(p, nzp)
            w = phi(x_nzp*t)
            ans.append([w[i]*phi_x_inv[i] for i in range(w.degree())])
        return ans

    def systems_of_eigenvalues(self, bound):
        P = prime_range(bound)
        e = self.eigenvalues(P)
        v = Sequence([], cr=True)
        if len(e) == 0:
            return v
        for i in range(len(e[0])):
            v.append([e[j][i] for j in range(len(e))])
        v.sort()
        v.set_immutable()
        return v

    def systems_of_abs(self, bound):
        P = prime_range(bound)
        e = self.eigenvalues(P)
        v = Sequence([], cr=True)
        if len(e) == 0:
            return v
        for i in range(len(e[0])):
            v.append([abs(e[j][i]) for j in range(len(e))])
        v.sort()
        v.set_immutable()
        return v

def support(v, eps):
    return [i for i in range(v.degree()) if abs(v[i]) > eps]
