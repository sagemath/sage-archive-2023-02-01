from sage.combinat.combinat import CombinatorialObject, CombinatorialClass
import sage.combinat.word as word
from sage.combinat.combination import Combinations
from sage.rings.all import QQ, PolynomialRing, prod


class LatticeDiagram(CombinatorialObject):
    def boxes(self):
        """
        EXAMPLES:
            sage: a = LatticeDiagram([3,0,2])
            sage: a.boxes()
            [(1, 1), (1, 2), (1, 3), (3, 1), (3, 2)]
            sage: a = LatticeDiagram([2, 1, 3, 0, 0, 2])
            sage: a.boxes()
            [(1, 1), (1, 2), (2, 1), (3, 1), (3, 2), (3, 3), (6, 1), (6, 2)]

        """
        res = []
        for i in range(1, len(self)+1):
            res += [ (i,j+1) for j in range(self[i]) ]
        return res

    def __getitem__(self, i):
        """
        Returns the $i^{th}$ entry of self. Note that the
        indexing starts for lattice diagrams starts at 1.

        EXAMPLES:
            sage: a = LatticeDiagram([3,0,2])
            sage: a[1]
            3
            sage: a[0]
            Traceback (most recent call last):
            ...
            ValueError: indexing starts at 1
            sage: a[-1]
            2
        """
        if i == 0:
            raise ValueError, "indexing starts at 1"
        elif i < 0:
            i += 1
        return self._list[i-1]

    def leg(self, i, j):
        """
        Returns the leg of the box (i,j) in self.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.leg(5,2)
            [(5, 3)]
        """
        return [(i,j) for j in range(j+1,self[i]+1)]

    def arm_left(self, i, j):
        """
        Returns the left arm of the box (i,j) in self.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm_left(5,2)
            [(1, 2), (3, 2)]
        """
        return [(ip,j) for ip in range(1,i) if self[ip] <= self[i] and j <= self[ip]]

    def arm_right(self, i, j):
        """
        Returns the right arm of the box (i,j) in self.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm_right(5,2)
            [(8, 1)]
        """
        return [(ip,j-1) for ip in range(i+1,len(self)+1) if self[ip] < self[i] and j-1 <= self[ip] ]

    def arm(self, i, j):
        """
        Returns the arm of the box (i,j) in self.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.arm(5,2)
            [(1, 2), (3, 2), (8, 1)]

        """
        return self.arm_left(i,j) + self.arm_right(i,j)

    def l(self, i, j):
        """
        Returns the self[i] - j.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.l(5,2)
            1
        """
        return self[i] - j


    def a(self, i, j):
        """
        Returns len(self.arm(i,j)).

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.a(5,2)
            3
        """
        return len(self.arm(i,j))


    def size(self):
        """
        Returns the number of boxes in self.

        EXAMPLES:
            sage: a = LatticeDiagram([3,1,2,4,3,0,4,2,3])
            sage: a.size()
            22
        """
        return sum(self._list)


    def flip(self):
        """
        Returns the flip of the self where flip is defined as follows.
        Let r = max(self).  Then self.flip()[i] = r - self[i].

        EXAMPLES:
            sage: a = LatticeDiagram([3,0,2])
            sage: a.flip()
            [0, 3, 1]
        """
        r = max(self)
        return LatticeDiagram([r-i for i in self])

class AugmentedLatticeDiagramFilling(CombinatorialObject):
    def __init__(self, l):
        """
        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a == loads(dumps(a))
            True
        """
        self._list = [[i+1]+l[i] for i in range(len(l))]

    def __getitem__(self, i):
        """
        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a[0]
            Traceback (most recent call last):
            ...
            ValueError: indexing starts at 1
            sage: a[1,0]
            1
            sage: a[2,0]
            2
            sage: a[3,2]
            4
            sage: a[3][2]
            4
        """
        if i < 1:
            raise ValueError, "indexing starts at 1"
        if isinstance(i, tuple):
            i,j = i
            return self._list[i-1][j]
        return self._list[i-1]

    def shape(self):
        """
        Returns the shape of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.shape()
            [2, 1, 3, 0, 0, 2]
        """
        return LatticeDiagram([max(0,len(self[i])-1) for i in range(1, len(self)+1)])

    def __contains__(self, ij):
        """
        Returns True if the box (i,j) (= ij) is in self.  Note
        that this does not include the basement row.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: (1,1) in a
            True
            sage: (1,0) in a
            False
        """
        i,j = ij
        if i > 0 and i <= len(self):
            if j > 0 and j <= len(self[i]):
                return True
        return False

    def are_attacking(self, i,j, ii, jj):
        """
        Returns True if the boxes (i,j) and (ii,jj) in self
        are attacking.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: all( a.are_attacking(i,j,ii,jj) for (i,j),(ii,jj) in a.attacking_boxes())
            True
            sage: a.are_attacking(1,1,3,2)
            False
        """
        if j == jj:
            return True

        if jj < j:
            i,j,ii,jj = ii,jj,i,j

        if j == jj - 1 and i > ii:
            return True

        return False

    def boxes(self):
        """
        Returns a list of the coordinates of the boxes of self, including
        the 'basement row'.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.boxes()
            [(1, 1),
             (1, 2),
             (2, 1),
             (3, 1),
             (3, 2),
             (3, 3),
             (6, 1),
             (6, 2),
             (1, 0),
             (2, 0),
             (3, 0),
             (4, 0),
             (5, 0),
             (6, 0)]
        """
        return self.shape().boxes() + [ (i,0) for i in range(1, len(self.shape())+1) ]

    def attacking_boxes(self):
        """
        Returns a list of pairs of boxes in self that are
        attacking.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.attacking_boxes()[:5]
            [((1, 1), (2, 1)),
             ((1, 1), (3, 1)),
             ((1, 1), (6, 1)),
             ((1, 1), (2, 0)),
             ((1, 1), (3, 0))]

        """
        boxes = self.boxes()
        res = []
        for (i,j),(ii,jj) in Combinations(boxes,2):
            if self.are_attacking(i,j,ii,jj):
                res.append( ((i,j),(ii,jj)) )
        return res

    def is_non_attacking(self):
        """
        Returns True if self in non-attacking.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.is_non_attacking()
            True
            sage: a = AugmentedLatticeDiagramFilling([[1, 1, 1], [2, 3], [3]])
            sage: a.is_non_attacking()
            False
        """
        for a,b in self.attacking_boxes():
            if self[a] == self[b]:
                return False
        for i in range(1, len(self)+1):
            if len(self[i]) > 1 and self[i,1] > i:
                return False
        return True

    def weight(self):
        """
        Returns the weight of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.weight()
            [1, 2, 1, 1, 2, 1]
        """
        return word.evaluation(self.reading_word())

    def descents(self):
        """
        Returns a list of the descents of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.descents()
            [(1, 2), (3, 2)]
        """
        res = []
        for i,j in self.shape().boxes():
            if self[i,j] > self[i,j-1]:
                res.append( (i,j) )
        return res

    def maj(self):
        """
        Returns the major index of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.maj()
            3
        """
        res = 0
        shape = self.shape()
        for i,j in self.descents():
            res += shape.l(i,j) + 1
        return res

    def reading_order(self):
        """
        Returns a list of coordinates of the boxes in self, starting
        from the top right, and reading from right to left.  Note that
        this includes the 'basement row' of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.reading_order()
            [(3, 3),
             (6, 2),
             (3, 2),
             (1, 2),
             (6, 1),
             (3, 1),
             (2, 1),
             (1, 1),
             (6, 0),
             (5, 0),
             (4, 0),
             (3, 0),
             (2, 0),
             (1, 0)]

        """
        boxes = self.boxes()
        f = lambda ij: (-ij[1],-ij[0])
        boxes.sort(key=f)
        return boxes

    def reading_word(self):
        """
        Return the reading word of self, obtained by reading the boxes
        entries of self from right to left, starting in the upper right.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.reading_word()
            [2, 5, 4, 6, 5, 3, 2, 1]
        """
        return [self[i,j] for i,j in self.reading_order() if j > 0]


    def inversions(self):
        """
        Returns a list of the inversions of self.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.inversions()[:5]
            [((6, 2), (3, 2)),
             ((1, 2), (6, 1)),
             ((1, 2), (3, 1)),
             ((1, 2), (2, 1)),
             ((6, 1), (3, 1))]
            sage: len(a.inversions())
            25
        """
        atboxes = [set(x) for x in self.attacking_boxes()]
        res = []
        order = self.reading_order()
        for a in range(len(order)):
            i,j = order[a]
            for b in range(a+1,len(order)):
                ii,jj = order[b]
                if self[i,j] > self[ii,jj] and  set( ((i,j),(ii,jj)) ) in atboxes:
                    res.append( ((i,j), (ii,jj)) )
        return res

    def _inv_aux(self):
        """
        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a._inv_aux()
            7
        """
        res = 0
        shape = self.shape()
        for i in range(1, len(self)+1):
            for j in range(i+1, len(self)+1):
                if shape[i] <= shape[j]:
                    res += 1
        return res


    def inv(self):
        """
        Returns self's inversion statistic.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.inv()
            15
        """
        res  = len(self.inversions())
        res -= sum(self.shape().a(i,j) for i,j in self.descents())
        res -= self._inv_aux()
        return res

    def coinv(self):
        """
        Returns self's co-inversion statistic.

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: a.coinv()
            2
        """
        shape = self.shape()
        return sum(shape.a(i,j) for i,j in shape.boxes()) - self.inv()


    def coeff(self, q, t):
        """
        Returns the coefficient in front of self in the HHL formula
        for the expansion of the non-symmetric Macdonald polynomial
        E(self.shape()).

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: q,t = var('q,t')
            sage: a.coeff(q,t)
            (1 - t)^4/((1 - q*t^2)^2*(1 - q^2*t^3)^2)

        """
        res = 1
        shape = self.shape()
        for i,j in shape.boxes():
            if self[i,j] != self[i,j-1]:
                res *= (1-t)/(1-q**(shape.l(i,j)+1)*t**(shape.a(i,j)+1))
        return res

    def coeff_integral(self, q, t):
        """
        Returns the coefficient in front of self in the HHL formula
        for the expansion of the integral non-symmetric Macdonald
        polynomial E(self.shape())

        EXAMPLES:
            sage: a = AugmentedLatticeDiagramFilling([[1,6],[2],[3,4,2],[],[],[5,5]])
            sage: q,t = var('q,t')
            sage: a.coeff_integral(q,t)
            (1 - t)^4*(1 - q*t^2)^2*(1 - q^2*t^3)^2

        """
        res = 1
        shape = self.shape()
        b = []
        for i,j in shape.boxes():
            if self[i,j] != self[i,j-1]:
                res *= (1-q**(shape.l(i,j)+1)*t**(shape.a(i,j)+1))
        for i,j in shape.boxes():
            if self[i,j] == self[i,j-1]:
                res *= (1-t)
        return res

def NonattackingFillings(shape):
    """
    Returning the combinatorial class of nonattacking
    fillings of a given shape.

    EXAMPLES:
        sage: NonattackingFillings([0,1,2])
        Nonattacking fillings of [0, 1, 2]
        sage: NonattackingFillings([0,1,2]).list()
        [[[1], [2, 1], [3, 2, 1]],
         [[1], [2, 1], [3, 2, 2]],
         [[1], [2, 1], [3, 2, 3]],
         [[1], [2, 1], [3, 3, 1]],
         [[1], [2, 1], [3, 3, 2]],
         [[1], [2, 1], [3, 3, 3]],
         [[1], [2, 2], [3, 1, 1]],
         [[1], [2, 2], [3, 1, 2]],
         [[1], [2, 2], [3, 1, 3]],
         [[1], [2, 2], [3, 3, 1]],
         [[1], [2, 2], [3, 3, 2]],
         [[1], [2, 2], [3, 3, 3]]]

    """
    return NonattackingFillings_shape(shape)

class NonattackingFillings_shape(CombinatorialClass):
    def __init__(self, shape):
        """
        EXAMPLES:
            sage: n = NonattackingFillings([0,1,2])
            sage: n == loads(dumps(n))
            True
        """
        self._shape = LatticeDiagram(shape)
        self._name = "Nonattacking fillings of %s"%shape

    def flip(self):
        """
        Returns the nonattacking fillings of the the flipped
        shape.

        EXAMPLES:
            sage: NonattackingFillings([0,1,2]).flip()
            Nonattacking fillings of [2, 1, 0]

        """
        return NonattackingFillings(list(self._shape.flip()))

    def iterator(self):
        """
        EXAMPLES:
            sage: NonattackingFillings([0,1,2]).list() #indirect doctest
            [[[1], [2, 1], [3, 2, 1]],
             [[1], [2, 1], [3, 2, 2]],
             [[1], [2, 1], [3, 2, 3]],
             [[1], [2, 1], [3, 3, 1]],
             [[1], [2, 1], [3, 3, 2]],
             [[1], [2, 1], [3, 3, 3]],
             [[1], [2, 2], [3, 1, 1]],
             [[1], [2, 2], [3, 1, 2]],
             [[1], [2, 2], [3, 1, 3]],
             [[1], [2, 2], [3, 3, 1]],
             [[1], [2, 2], [3, 3, 2]],
             [[1], [2, 2], [3, 3, 3]]]
            sage: len(_)
            12

        """
        k = len(self._shape)
        n = self._shape.size()
        for p in word.Words(k,n):
            c = 0
            res = []
            for i in range(1, len(self._shape)+1):
                res.append(p[c:c+self._shape[i]])
                c += self._shape[i]
            x = AugmentedLatticeDiagramFilling(res)
            if x.is_non_attacking():
                yield x

def _check_muqt(mu, q, t):
    """
    EXAMPLES:
        sage: from sage.combinat.sf.ns_macdonald import _check_muqt
        sage: P, q, t, n, R, x = _check_muqt([0,0,1],None,None)
        sage: P
        Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: q
        q
        sage: t
        t
        sage: n
        Nonattacking fillings of [0, 0, 1]
        sage: R
        Multivariate Polynomial Ring in x0, x1, x2 over Fraction Field of Multivariate Polynomial Ring in q, t over Rational Field
        sage: x
        (x0, x1, x2)

        sage: q,t = var('q,t')
        sage: P, q, t, n, R, x = _check_muqt([0,0,1],q,None)
        Traceback (most recent call last):
        ...
        ValueError: you must specify either both q and t or neither of them

        sage: P, q, t, n, R, x = _check_muqt([0,0,1],q,2)
        Traceback (most recent call last):
        ...
        ValueError: the parents of q and t must be the same

    """
    if q is None and t is None:
        P = PolynomialRing(QQ,'q,t').fraction_field()
        q,t = P.gens()
    elif q is not None and t is not None:
        if q.parent() != t.parent():
            raise ValueError, "the parents of q and t must be the same"
        P = q.parent()
    else:
        raise ValueError, "you must specify either both q and t or neither of them"
    n = NonattackingFillings(mu)
    R = PolynomialRing(P, len(n._shape), 'x')
    x = R.gens()
    return P, q, t, n, R, x

def E(mu, q=None, t=None):
    """
    Returns the non-symmetric Mcadonald polynomial in type A corresponding to
    a shape mu.

    Note that if both q and t are specifed, then they must have the same parent.

    REFERENCE:
        'A combinatorial formula for non-symmetric Macdonald polynomials'. Haiman,
            Haglund, and Loehr. http://arxiv.org/abs/math/0601693

    EXAMPLES:
        sage: from sage.combinat.sf.ns_macdonald import E
        sage: E([0,0,0])
        1
        sage: E([1,0,0])
        x0
        sage: E([0,1,0])
        ((-t + 1)/(-q*t^2 + 1))*x0 + x1
        sage: E([0,0,1])
        ((-t + 1)/(-q*t + 1))*x0 + ((-t + 1)/(-q*t + 1))*x1 + x2
        sage: E([1,1,0])
        x0*x1
        sage: E([1,0,1])
        ((-t + 1)/(-q*t^2 + 1))*x0*x1 + x0*x2
        sage: E([0,1,1])
        ((-t + 1)/(-q*t + 1))*x0*x1 + ((-t + 1)/(-q*t + 1))*x0*x2 + x1*x2
        sage: E([2,0,0])
        x0^2 + ((-q*t + q)/(-q*t + 1))*x0*x1 + ((-q*t + q)/(-q*t + 1))*x0*x2
        sage: E([0,2,0])
        ((-t + 1)/(-q^2*t^2 + 1))*x0^2 + ((-q^2*t^3 + q^2*t^2 - q*t^2 + 2*q*t - q + t - 1)/(-q^3*t^3 + q^2*t^2 + q*t - 1))*x0*x1 + x1^2 + ((q*t^2 - 2*q*t + q)/(q^3*t^3 - q^2*t^2 - q*t + 1))*x0*x2 + ((-q*t + q)/(-q*t + 1))*x1*x2

    """
    P, q, t, n, R, x = _check_muqt(mu, q, t)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj()*t**a.coinv()*a.coeff(q,t)*prod( x[i]**weight[i] for i in range(len(weight)) )
    return res

def E_integral(mu, q=None, t=None):
    """
    Returns the integral form for the non-symmetric Mcadonald polynomial
    in type A corresponding to a shape mu.

    Note that if both q and t are specifed, then they must have
    the same parent.

    REFERENCE:
        'A combinatorial formula for non-symmetric Macdonald polynomials'. Haiman,
            Haglund, and Loehr. http://arxiv.org/abs/math/0601693

    EXAMPLES:
        sage: from sage.combinat.sf.ns_macdonald import E_integral
        sage: E_integral([0,0,0])
        1
        sage: E_integral([1,0,0])
        (-t + 1)*x0
        sage: E_integral([0,1,0])
        (-q*t^2 + 1)*x0 + (-t + 1)*x1
        sage: E_integral([0,0,1])
        (-q*t + 1)*x0 + (-q*t + 1)*x1 + (-t + 1)*x2
        sage: E_integral([1,1,0])
        (t^2 - 2*t + 1)*x0*x1
        sage: E_integral([1,0,1])
        (q*t^3 - q*t^2 - t + 1)*x0*x1 + (t^2 - 2*t + 1)*x0*x2
        sage: E_integral([0,1,1])
        (q^2*t^3 + q*t^4 - q*t^3 - q*t^2 - q*t - t^2 + t + 1)*x0*x1 + (q*t^2 - q*t - t + 1)*x0*x2 + (t^2 - 2*t + 1)*x1*x2
        sage: E_integral([2,0,0])
        (t^2 - 2*t + 1)*x0^2 + (q^2*t^2 - q^2*t - q*t + q)*x0*x1 + (q^2*t^2 - q^2*t - q*t + q)*x0*x2
        sage: E_integral([0,2,0])
        (q^2*t^3 - q^2*t^2 - t + 1)*x0^2 + (q^4*t^3 - q^3*t^2 - q^2*t + q*t^2 - q*t + q - t + 1)*x0*x1 + (t^2 - 2*t + 1)*x1^2 + (q^4*t^3 - q^3*t^2 - q^2*t + q)*x0*x2 + (q^2*t^2 - q^2*t - q*t + q)*x1*x2

    """
    P, q, t, n, R, x = _check_muqt(mu, q, t)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj()*t**a.coinv()*a.coeff_integral(q,t)*prod( x[i]**weight[i] for i in range(len(weight)) )
    return res

def Ht(mu, q=None, t=None):
    """
    Returns the symmetric Macdonald polynomial using the
    Haiman, Haglund, and Loehr formula.

    Note that if both q and t are specifed, then they must have
    the same parent.

    REFERENCE:
        'A combinatorial formula for non-symmetric Macdonald polynomials'. Haiman,
            Haglund, and Loehr. http://arxiv.org/abs/math/0601693

    EXAMPLES:
        sage: from sage.combinat.sf.ns_macdonald import Ht
        sage: HHt = MacdonaldPolynomialsHt(QQ)
        sage: Ht([0,0,1])
        x0 + x1 + x2
        sage: HHt([1]).expand(3)
        x0 + x1 + x2
        sage: Ht([0,0,2])
        x0^2 + (q + 1)*x0*x1 + x1^2 + (q + 1)*x0*x2 + (q + 1)*x1*x2 + x2^2
        sage: HHt([2]).expand(3)
        x0^2 + (q + 1)*x0*x1 + x1^2 + (q + 1)*x0*x2 + (q + 1)*x1*x2 + x2^2

    """
    P, q, t, n, R, x = _check_muqt(mu, q, t)
    res = 0
    for a in n:
        weight = a.weight()
        res += q**a.maj()*t**a.inv()*prod( x[i]**weight[i] for i in range(len(weight)) )
    return res
