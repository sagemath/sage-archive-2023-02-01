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
from sage.modules.all           import vector
from sage.misc.misc             import verbose
from sage.rings.all             import CDF, Integer, QQ, next_prime, prime_range
from random                     import randint

class NumericalEigenforms(SageObject):
    """
    numerical_eigenforms(group, weight=2, eps=1e-20,
                         delta=1e-2, tp=[2,3,5])

    INPUT:
        group -- a congruence subgroup of a Dirichlet character of
                 order 1 or 2
        weight -- an integer >= 2
        eps -- a small float; abs( ) < eps is what "equal to zero" is
               interpreted as for floating point numbers.
        delta -- a small-ish float; eigenvalues are considered distinct
                 if their difference has absolute value at least delta
        tp -- use the Hecke operators T_p for p in tp when searching
              for a random Hecke operator with distinct Hecke eigenvalues.

    OUTPUT:
        a numerical eigenforms object, with the following useful methods:
            * ap(p) -- return all eigenvalues of $T_p$
            * eigenvalues(primes) -- list of eigenvalues corresponding
                    to the given list of primes, e.g.,:
                        [[eigenvalues of T_2],
                         [eigenvalues of T_3],
                         [eigenvalues of T_5], ...]
            * systems_of_eigenvalues -- a list of the systems of
                 eigenvalues of eigenforms such that the chosen
                 random linear combination of Hecke operators has
                 multiplicity 1 eigenvalues.

    EXAMPLES:
        sage: n = numerical_eigenforms(23)
        sage: n.ap(2)
        [3.0, 0.61803398875, -1.61803398875]
        sage: n.systems_of_eigenvalues(7)
        [
        [-1.61803398875, 2.2360679775, -3.2360679775],
        [0.61803398875, -2.2360679775, 1.2360679775],
        [3.0, 4.0, 6.0]
        ]
        sage: n.systems_of_abs(7)   # random low order bits
        [
        [0.61803398874989446, 2.2360679774997911, 1.2360679774997889],
        [1.6180339887498947, 2.2360679774997894, 3.2360679774997894],
        [3.0, 4.0, 6.0]
        ]
        sage: n.eigenvalues([2,3,5])
        [[3.0, 0.61803398875, -1.61803398875],
         [4.0, -2.2360679775, 2.2360679775],
         [6.0, 1.2360679775, -3.2360679775]]
    """
    def __init__(self, group, weight=2, eps=1e-20,
                 delta=1e-2, tp=[2,3,5]):
        if isinstance(group, (int, long, Integer)):
            group = Gamma0(Integer(group))
        self._group  = group
        self._weight = Integer(weight)
        self._tp = tp
        if self._weight < 2:
            raise ValueError, "weight must be at least 2"
        self._eps = eps
        self._delta = delta

    def level(self):
        """
        Return the level of this set of modular eigenforms.
        """
        return self._group.level()

    def weight(self):
        """
        Return the weight of this set of modular eigenforms.
        """
        return self._weight

    def _repr_(self):
        return "Numerical Hecke eigenvalues for %s of weight %s"%(
            self._group, self._weight)

    def modular_symbols(self):
        """
        Return the space of modular symbols used for computing this
        set of modular eigenforms.
        """
        try:
            return self.__modular_symbols
        except AttributeError:
            M = ModularSymbols(self._group,
                    self._weight, sign=1)
            if M.base_ring() != QQ:
                raise ValueError, "modular forms space must be defined over QQ"
            self.__modular_symbols = M
            return M

    def _eigenvectors(self):
        try:
            return self.__eigenvectors
        except AttributeError:
            pass
        verbose('Finding eigenvector basis')
        M = self.modular_symbols()
        N = self.level()

        tp = self._tp
        p = tp[0]
        t = M.T(p).matrix()
        for p in tp[1:]:
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
            sage: n = numerical_eigenforms(37)
            sage: n._easy_vector()
            (1.0, 1.0, 0)
            sage: n = numerical_eigenforms(43)
            sage: n._easy_vector()
            (1.0, 0, 1.0, 0)
            sage: n = numerical_eigenforms(125)
            sage: n._easy_vector()
            (0, 0, 0, 1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
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
        """
        Return a list of the eigenvalues of the Hecke operator $T_p$
        on all the computed eigenforms.  The eigenvalues match up
        between one prime and the next.

        INPUT:
            p -- integer, a prime number

        OUTPUT:
            list -- a list of double precision complex numbers

        EXAMPLES:
            sage: n = numerical_eigenforms(11,4)
            sage: n.ap(2)  # random low order bits
            [9.0, 9.0, 2.73205080757, -0.732050807569]
            sage: n.ap(3)  # random low order bits
            [28.0, 28.0, -7.92820323028, 5.92820323028]
            sage: m = n.modular_symbols()
            sage: x = polygen(QQ, 'x')
            sage: m.T(2).charpoly(x).factor()
            (x - 9)^2 * (x^2 - 2*x - 2)
            sage: m.T(3).charpoly(x).factor()
            (x - 28)^2 * (x^2 + 2*x - 47)
        """
        p = Integer(p)
        if not p.is_prime():
            raise ValueError, "p must be a prime"
        try:
            return self._ap[p]
        except AttributeError:
            self._ap = {}
        except KeyError:
            pass
        a = Sequence(self.eigenvalues([p])[0], immutable=True)
        self._ap[p] = a
        return a

    def eigenvalues(self, primes):
        """
        Return the eigenvalues of the Hecke operators corresponding
        to the primes in the input list of primes.   The eigenvalues
        match up between one prime and the next.

        INPUT:
            primes -- a list of primes

        OUTPUT:
            list of lists of eigenvalues.

        EXAMPLES:
            sage: n = numerical_eigenforms(1,12)
            sage: n.eigenvalues([3,5,13])            # random low-order bits
            [[177148.0, 252.0], [48828126.0, 4830.0], [1.79216039404e+12, -577737.99999]]
        """
        primes = [Integer(p) for p in primes]
        for p in primes:
            if not p.is_prime():
                raise ValueError, 'each element of primes must be prime.'
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



