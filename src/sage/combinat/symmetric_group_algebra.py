#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

from combinatorial_algebra import CombinatorialAlgebra
import permutation
import partition
from tableau import Tableau, StandardTableaux_n, StandardTableaux_partition, StandardTableaux
from sage.interfaces.all import gap
from sage.rings.all import factorial, QQ

def SymmetricGroupAlgebra(R,n):
    """
    Returns the symmetric group algebra of order n over R.

    EXAMPLES:
        sage: QS3 = SymmetricGroupAlgebra(QQ, 3); QS3
        Symmetric group algebra of order 3 over Rational Field
        sage: QS3(1)
        [1, 2, 3]
        sage: QS3(2)
        2*[1, 2, 3]
        sage: basis = [QS3(p) for p in Permutations(3)]
        sage: a = sum(basis); a
        [3, 1, 2] + [1, 2, 3] + [2, 3, 1] + [2, 1, 3] + [3, 2, 1] + [1, 3, 2]
        sage: a^2
        6*[3, 1, 2] + 6*[1, 2, 3] + 6*[2, 3, 1] + 6*[2, 1, 3] + 6*[3, 2, 1] + 6*[1, 3, 2]
        sage: a^2 == 6*a
        True
        sage: b = QS3([3, 1, 2])
        sage: b
        [3, 1, 2]
        sage: b*a
        [3, 1, 2] + [1, 2, 3] + [2, 3, 1] + [2, 1, 3] + [3, 2, 1] + [1, 3, 2]
        sage: b*a == a
        True
    """
    return SymmetricGroupAlgebra_n(R,n)

class SymmetricGroupAlgebra_n(CombinatorialAlgebra):
    def __init__(self, R, n):
        """
        TESTS:
            #sage: QS3 = SymmetricGroupAlgebra(QQ, 3)
            #sage: QS3 == loads(dumps(QS3))
            #True
        """
        self.n = n
        self._combinatorial_class = permutation.Permutations(n)
        self._name = "Symmetric group algebra of order %s"%self.n
        self._one = permutation.Permutation(range(1,n+1))
        self._prefix = ""
        CombinatorialAlgebra.__init__(self, R)

    def _multiply_basis(self, left, right):
        return left * right

    def _coerce_start(self, x):
        if x == []:
            return self( self._one )
        if len(x) < self.n and x in permutation.Permutations():
            return self( list(x) + range(len(x)+1, self.n+1) )
        raise TypeError

    def cpis(self):
        """
        Returns a list of the centrally primitive idempotents.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: a = QS3.cpis()
            sage: a[0]  # [3]
            1/6*[3, 1, 2] + 1/6*[1, 2, 3] + 1/6*[2, 3, 1] + 1/6*[2, 1, 3] + 1/6*[3, 2, 1] + 1/6*[1, 3, 2]
            sage: a[1]  # [2, 1]
            -1/3*[2, 3, 1] + 2/3*[1, 2, 3] - 1/3*[3, 1, 2]
        """
        return [self.cpi(p) for p in partition.Partitions_n(self.n)]

    def cpi(self, p):
        """
        Returns the centrally primitive idempotent for the symmetric
        group of order n for the irreducible corresponding indexed by
        the partition p.

        EXAMPLES:
            sage: QS3 = SymmetricGroupAlgebra(QQ,3)
            sage: QS3.cpi([2,1])
            -1/3*[2, 3, 1] + 2/3*[1, 2, 3] - 1/3*[3, 1, 2]
            sage: QS3.cpi([3])
            1/6*[3, 1, 2] + 1/6*[1, 2, 3] + 1/6*[2, 3, 1] + 1/6*[2, 1, 3] + 1/6*[3, 2, 1] + 1/6*[1, 3, 2]
            sage: QS3.cpi([1,1,1])
            1/6*[3, 1, 2] + 1/6*[1, 2, 3] + 1/6*[2, 3, 1] - 1/6*[2, 1, 3] - 1/6*[3, 2, 1] - 1/6*[1, 3, 2]
        """
        if p not in partition.Partitions_n(self.n):
            raise TypeError, "p must be a partition of %s"%self.n

        character_table = eval(gap.eval("Display(Irr(SymmetricGroup(%d)));"%self.n))

        cpi = self(0)

        np = partition.Partitions_n(self.n).list()
        np.reverse()
        p_index = np.index(p)

        big_coeff = character_table[p_index][0]/factorial(self.n)

        for g in permutation.StandardPermutations_n(self.n):
            cpi += big_coeff * character_table[p_index][np.index(g.inverse().cycle_type())] * self(g)

        return cpi



    def seminormal_basis(self):
        """
        Returns a list of the seminormal basis elements of self.

        EXAMPLES:

        """
	basis = []
	for part in partition.Partitions_n(self.n):
            stp = StandardTableaux_partition(part)
            for t1 in stp:
                for t2 in stp:
                    basis.append(self.epsilon_ik(t1,t2))
	return basis


    def epsilon_ik(self, itab, ktab, star=0):
        """
        Returns the seminormal basis element of self corresponding
        to the pair of tableaux itab and ktab.

        EXAMPLES:

        """

        it = Tableau(itab)
        kt = Tableau(ktab)

        stn = StandardTableaux_n(self.n)

        if it not in stn:
            raise TypeError, "it must be a standard tableaux of size %s"%self.n

        if kt not in stn:
            raise TypeError, "kt must be a standard tableaux of size %s"%self.n

        return epsilon_ik(it, kt, star=star)





epsilon_ik_cache = {}
def epsilon_ik(itab, ktab, star=0):
    it = Tableau(itab)
    kt = Tableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError, "the two tableaux must be of the same shape"
    elif kt == it:
        return epsilon(itab)
    elif (it, kt) in epsilon_ik_cache:
        return epsilon_ik_cache[(it,kt)]
    else:
        epsilon_ik_cache[(it,kt)] = epsilon(it, star+1)*e_ik(it,kt,star)*epsilon(kt, star+1) * (1/kappa(it.shape()))
        return epsilon_ik_cache[(it,kt)]




epsilon_cache = {}
def epsilon(tab, star=0):
    if star:
        t2 = Tableau(tab)
        t = t2.restrict(t2.size() - star)
    else:
        t = Tableau(tab)

    if t in epsilon_cache:
        return epsilon_cache[t]
    else:
        if t.size() == 2:
            epsilon_cache[t] = e(t)*(1/kappa(t.shape()))
            return epsilon_cache[t]
        #elif t == Tableau([]):
        #	epsilon_cache[t] = SymmetricGroupAlgebraElement(1, [1])
        #	return epsilon_cache[t]
        elif t == Tableau([[1]]):
            epsilon_cache[t] = e(t)
            return epsilon_cache[t]
        else:
            epsilon_cache[t] =  epsilon(t, 1)*e(t)*epsilon(t,1)*( 1 / kappa(t.shape()))
            return epsilon_cache[t]




def pi_ik(itab, ktab, csn=False):
    it = Tableau(itab)
    kt = Tableau(ktab)
    p = [None]*it.size()
    for i in range(len(it)):
        for j in range(len(it[i])):
            p[ kt[i][j] -1 ] = it[i][j]

    QSn = SymmetricGroupAlgebra(QQ, it.size())
    po = permutation.PermutationOptions()
    p = permutation.Permutation(p)
    if po['mult'] == 'l2r':
        return QSn(p.inverse())
    else:
        return QSn(p)


def kappa(alpha):
    try:
        n = alpha.size()
    except:
        n = sum(alpha)
    return factorial(n)/StandardTableaux(alpha).count()


e_cache = {}
def e(tableau, star=0):
    t = Tableau(tableau)
    if star:
        t = t.restrict(t.size()-star)

    if t in e_cache:
        return e_cache[t]
    else:
        rs = t.row_stabilizer()
        cs = t.column_stabilizer()
        n = t.size()

        QSn = SymmetricGroupAlgebra(QQ, n)

        res = 0
        for h in rs:
            for v in cs:
                res += v.sign() * QSn( (h*v).list() )
        e_cache[t] = res
        return res

ehat_cache = {}
def e_hat(tab, star=0):
    t = Tableau(tab)
    if star:
        t = t.restrict(t.size()-star)
    if t in ehat_cache:
        return ehat_cache[t]
    else:
        return (1/kappa(t.shape()))*e(t)

e_ik_cache = {}
def e_ik(itab, ktab, star=0):
    it = Tableau(itab)
    kt = Tableau(ktab)
    if star:
        it = it.restrict(it.size() - star)
        kt = kt.restrict(kt.size() - star)

    if it.shape() != kt.shape():
        raise ValueError, "the two tableaux must be of the same shape"
    elif kt == it:
        return e(itab)
    elif (it, kt) in e_ik_cache:
        return e_ik_cache[(it,kt)]
    else:
        pi = pi_ik(it,kt)
        e_ik_cache[(it,kt)] = e(it)*pi
        return e_ik_cache[(it,kt)]


def seminormal_test(n):
    for part in partition.Partitions_n(n):
        for tab in StandardTableaux(part):
            #3.1.10
            if not e(tab)*(1/kappa(part)) - e_hat(tab) == 0:
                raise ValueError, "3.1.10 - %s"%tab

            #3.2.12.2
            value = e(tab)*epsilon(tab,1)*e(tab) - e(tab)*(kappa(part))
            if value != 0:
                print value
                raise ValueError, "3.2.12.2 - %s"%tab

            for tab2 in StandardTableaux(part):
                #3.2.8 1
                if e_ik(tab, tab2) - e(tab)*pi_ik(tab, tab2)*e(tab2)*(1/kappa(part)) != 0:
                    raise ValueError, "3.2.8.1 - %s, %s"%(tab, tab2)

                #3.2.8.1
                if e(tab)*e_ik(tab, tab2) - e_ik(tab, tab2)*(kappa(part)) != 0:
                    raise ValueError, "3.2.8.2 - %s, %s"%(tab, tab2)

                if tab == tab2:
                    continue

                #3.1.20
                if e_hat(tab)*e_hat(tab2) != 0:
                    raise ValueError, "3.1.20 - %s, %s"%(tab, tab2)
