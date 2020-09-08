from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.enumerated_sets import EnumeratedSets
from sage.combinat.permutation import Permutations
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat import permutation
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage.monoids.hecke_monoid import HeckeMonoid
from sage.misc.lazy_attribute import lazy_attribute

from sage.rings.integer import Integer

class DecreasingHeckeFactorization:
    """
    Class of decreasing factorizations in the 0-Hecke monoid.
    
    EXAMPLES::

        sage: from sage.combinat.crystals.fully_commutative_stable_grothendieck import DecreasingHeckeFactorization
        sage: t = [[3,2],[],[2,1]]
        sage: h = DecreasingHeckeFactorization(t,3); h
        (3, 2)()(2, 1)
        sage: h.excess
        1
        sage: h.m
        3
        sage: h.k
        3
        sage: h.value
        ((3, 2), (), (2, 1))

        sage: u = [[3,2,1],[3],[2,1]]
        sage: h = HeckeFactorization(u,4); h
        (3, 2, 1)(3)(2, 1)
        sage: h.weight()
        [2, 1, 3]
    """
    def __init__(self, t, k=None):
        """
        Initialize a decreasing factorization for ``self`` given the relevant data.
        """
        if not isinstance(t,(list,tuple)):
            raise ValueError("A list or tuple is expected")
        for factor in t:
            for r, value in enumerate(factor):
                if r<len(factor)-1 and factor[r]<=factor[r+1]:
                    raise ValueError("All the factors need to be decreasing")
        u = [list(factor) for factor in t]
        if k == None:
            k = max([x for factor in t for x in factor])
        self.k = k
        H = HeckeMonoid(SymmetricGroup(k+1))
        w = H.from_reduced_word([x for factor in t for x in factor])
        self.excess = sum(len(l) for l in t)-len(w.reduced_word())
        self.value = tuple([tuple(factors) for factors in t])
    
    @property
    def m(self):
        return len(self.value)

    def __repr__(self):
        return "".join("("+repr(list(factor))[1:-1]+")" for factor in self.value)

    def __key(self):
        return (self.k, self.value)

    def __hash__(self):
        return hash(self.__key())

    def __eq__(self, other):
        return isinstance(self, type(other)) and self.value == other.value

    def __lt__(self,other):
        return (self._weight(),self.value) < (other._weight(),other.value)

    def _latex_(self):
        r"""
        Return LaTeX code for ``self``.

        EXAMPLES::

            sage:

        """
        s = ""
        for factor in self.value:
            if len(factor)>0:
                s += r"\left("+repr(list(factor))[1:-1]+r"\right)"
            else:
                s += r"\left(\;\right)"
        return s

    def weight(self):
        """
        Returns the weight of ``self``

        EXAMPLES:

            sage: t = [[2],[2,1],[],[4,3,1]]
            sage: u = DecreasingHeckeFactorization(t,4)
            sage: t.weight()
            [3, 0, 2, 1]
        """
        return [len(l) for l in self.value][::-1]

    def to_word(self):
        """
        Return the word associated to ``self`` in the 0-Hecke monoid

        EXAMPLES:

            sage: t = [[2],[],[2,1],[4,3,1]]
            sage: u = DecreasingHeckeFactorization(t,4)
            sage: u.to_word()
            [2, 2, 1, 4, 3, 1]
        """
        return [j for factors in self.value for j in factors]

    def to_increasing_hecke_biword(self):
        """
        Return the associated increasing Hecke biword of ``self``

        EXAMPLES:

            sage: t = [[2],[],[2,1],[4,3,1]]
            sage: u = DecreasingHeckeFactorization(t,4)
            sage: u.to_increasing_hecke_biword()
            [[1, 1, 1, 2, 2, 4], [1, 3, 4, 1, 2, 2]]
        """
        L = [[],[]]
        for j in range(len(self.value)):
            L[1] += list(self.value[-j-1][::-1])
            L[0] += [j+1]*len(self.value[-j-1])
        return L


class FullyCommutativeStableGrothendieckCrystal(UniqueRepresentation, Parent):
    """
    The crystal on fully commutative decreasing factorizations in the 0-Hecke 
    monoid, as introduced by [MPPS2020]_.

    INPUT:

    - ``w`` -- an element in a 0-Hecke monoid

    - ``m`` -- the number of factors in the factorization

    - ``excess`` -- the total number of letters in the factorization minus the length of a reduced word for ``w``

    EXAMPLES::

        sage: from sage.monoids.hecke_monoid import HeckeMonoid
        sage: H = HeckeMonoid(SymmetricGroup(3+1))
        sage: w = H.from_reduced_word([1,3,2])
        sage: B = FullyCommutativeStableGrothendieckCrystal(w,3,2); B
        Fully commutative stable Grothendieck crystal of type A_2 associated to [1, 3, 2] with excess 2
    """
    def __init__(self, w, m, excess):
        """
        Initialize a crystal for self given reduced word ``w`` in a 0-Hecke 
        monoid, number of factors ``m`` and``excess`` extra letters.

        EXAMPLES::

            sage: from sage.monoids.hecke_monoid import HeckeMonoid
            sage: H = HeckeMonoid(SymmetricGroup(3+1))
            sage: w = H.from_reduced_word([1,3,2])
            sage: B = FullyCommutativeStableGrothendieckCrystal(w,3,2)
            sage: B.w
            (1, 3, 2)
            sage: B.m
            3
            sage: B.excess
            2
            sage: B.H
            0-Hecke monoid of the Symmetric group of order 4! as a permutation group

            The reduced word ``w`` should be fully commutative, that is, its 
            associated permutation should avoid the pattern 321.

            sage: from sage.monoids.hecke_monoid import HeckeMonoid
            sage: H = HeckeMonoid(SymmetricGroup(3+1))
            sage: w = H.from_reduced_word([1,2,1])
            sage: B = FullyCommutativeStableGrothendieckCrystal(w,4,2)
            Traceback (most recent call last):
            ...
            ValueError: w should be fully commutative

        TESTS::

            sage: from sage.monoids.hecke_monoid import HeckeMonoid
            sage: H = HeckeMonoid(SymmetricGroup(3+1))
            sage: w = H.from_reduced_word([2,1,3,2])
            sage: B = FullyCommutativeStableGrothendieckCrystal(w,4,2)
            sage: TestSuite(B).run()
        """
        # Check if w is fully commutative, ie. if the associated permutation is
        # 321-avoiding 
        word = w.reduced_word()
        p = permutation.from_reduced_word(word)
        if p not in Permutations(avoiding=[3,2,1]):
            raise ValueError("w should be fully commutative")

        Parent.__init__(self, category = ClassicalCrystals())
        self.w = tuple(word)
        self.m = m
        self.H = w.parent()
        self.k = len(self.H.gens())
        self.excess = excess
        self._cartan_type = CartanType(['A', self.m-1])

    @lazy_attribute
    def module_generators(self):
        """
        Return module generators for ``self``.

        EXAMPLES::

            sage: H = HeckeMonoid(SymmetricGroup(3))
            sage: w = H.from_reduced_word([1,3,2])
            sage: B = FullyCommutativeStableGrothendieckCrystal(w,3,2)
            sage: B.module_generators
            ((1)(3, 1)(3, 2), (3)(3, 1)(3, 2))
        """
        return tuple(self(x).to_highest_weight()[0] for x in _lowest_weights(self.h, self.m, self.excess))

    def _repr_(self):
        """
        Return a representation of ``self''.

        EXAMPLES::
        
            sage: H = HeckeMonoid(SymmetricGroup(3+1))
            sage: w = H.from_reduced_word([2,1,3,2])
            sage: FullyCommutativeStableGrothendieckCrystal(w,3,1)
            Fully commutative stable Grothendieck crystal of type A_2 associated to [2, 1, 3, 2] with excess 1
        """
        return "Fully commutative stable Grothendieck crystal of type A_{} associated to {} with excess {}".format(self.m-1, list(self.w), self.excess)

    # temporary workaround while an_element is overriden by Parent
    _an_element_ = EnumeratedSets.ParentMethods._an_element_ 

    class Element(ElementWrapper):

        def e(self, i):
            L = list(self.value[self.m-i-1])
            R = list(self.value[self.m-i])
            b = bracketing_eq(L,R)
            if not b[0]:
                return None
            
            y = b[0][-1]
            if y-1 in L and y-1 in R:
                # special case: (--x+1--)(--x+1,x--) -->> (--x+1,x--)(--x--)
                L.remove(y-1)
            else:
                L.remove(y)
            R.append(y)
            L.sort(reverse=True)
            R.sort(reverse=True)  
            s = DecreasingHeckeFactorization([self.value[j] for j in range(self.m-i-1)] + [L] + [R] + [self.value[j] for j in range(self.m-i+1,self.m)],self.k)
            if not check_local(self,s,i):
                s.bright = True
                print(self,s,i)           
            return self.parent()(s)

        def f(self, i):
            L = list(self.value[self.m-i-1])
            R = list(self.value[self.m-i])
            b = bracketing_eq(L,R)
            if not b[1]:
                return None
            
            x = b[1][0]
            if x+1 in L and x+1 in R:
                # special case: (--x+1--)(--x+1,x--) -->> (--x+1,x--)(--x--)
                R.remove(x+1)
            else:
                R.remove(x)
            L.append(x)
            L.sort(reverse=True)
            R.sort(reverse=True)  
            s = DecreasingHeckeFactorization([self.value[j] for j in range(self.m-i-1)] + [L] + [R] + [self.value[j] for j in range(self.m-i+1,self.m)],self.k)
            if not check_local(self,s,i):
                s.bright = True
                print(self,s,i)  
            return self.parent()(s)

def bracketing_eq(L,R):
    """
        Removes all bracketed letters between `i`-th and `i+1`-th entry.
    """
    right_n = [j for j in R]
    left_n = [j for j in L]
    left_unbracketed = []
    while left_n:
        m = max(left_n)
        left_n.remove(m)
        l = [j for j in right_n if j>=m]
        if l:
            right_n.remove(min(l))
        else:
            left_unbracketed += [m]
    return [[j for j in left_unbracketed],[j for j in right_n]]

def _lowest_weights(w,m,ex):
    """
    Generate all 0-Hecke decreasing factorizations corresponding to 
    some valid semistandard Young tableau.

    The semistandard Young tableaux should have at most m columns 
    and their column reading words should be equivalent to w in some 
    0-Hecke monoid.

    INPUTS:

        w - 321-avoiding reduced word in a 0-Hecke monoid
        
        m - number of factors in decreasing 0-Hecke factorizations
        
        ex - number of extra letters in the decreasing factorizations

    EXAMPLES:

        sage: _lowest_weight_factorizations([1, 2, 1],3,1)
        Traceback (most recent call last):
        ...
        ValueError: The 0-Hecke word w must be 321-avoiding.  

        sage: _lowest_weights([2, 1, 3, 2],4,3)
        [(2, 1)(3, 1)(3, 1)(2), (2, 1)(3, 1)(3, 2)(2)]
        
        sage: _lowest_weights([2, 1, 3, 2],5,3)
        [(2, 1)(3, 1)(3, 1)(2)(),
         (2, 1)(3, 1)(3, 2)(2)(),
         (2, 1)(3, 1)(1)(1)(2),
         (2, 1)(3, 1)(1)(2)(2),
         (2, 1)(3, 1)(2)(2)(2),
         (2, 1)(3, 2)(2)(2)(2)]

        sage: _lowest_weights([1, 3],3,1)
        [(3, 1)(1)(), (3, 1)(3)(), (1)(1)(3), (1)(3)(3)]

        sage: _lowest_weights([3, 2, 1],5,2)
        [(3, 2, 1)(1)(1)()()]
    """
    k = max(w) if w else 1
    P = Permutations(k+1, avoiding=[3,2,1])
    p = permutation.from_reduced_word(w,P)
    if p not in P:
        raise ValueError("The 0-Hecke word w must be 321-avoiding")

    L = list_equivalent_words(canonical_word(w,ex))
    D = {}
    k = max(w)
    for v in L:
        if is_valid_col_word(v, m) == True:
            J = [0] + jumps(v) + [len(v)]
            t = [v[J[i]:J[i+1]] for i in range(len(J)-1)]
            if len(J) < m+1:
                t += [()]*(m+1-len(J))
            h = DecreasingHeckeFactorization(t, k)
            wt = lambda h: [len(l) for l in h.value][::-1]
            weight = tuple(wt(h))
            if weight not in D:
                D[weight] = [h]
            else:
                D[weight] += [h]
    return sorted([h for key in D for h in D[key]], key=lambda h:([i for i in wt(h)],h.value))

########################################################

def canonical_word(w, ex):
    """
    Return a standard Hecke factorization equivalent to w in 
    some 0-Hecke monoid with excess ex

    EXAMPLES:

    sage: w = [1,2,1]
    sage: v = canonical_word(w,ex=2); v
    [1, 1, 1, 2, 1]
    """
    if isinstance(w,list)==True or isinstance(w,tuple)==True:
        L = list(w)
        return [L[0]]*ex+L

def jumps(w):
    """
    Detect positions where indices weakly increases in w
    
    EXAMPLES:

        sage: w = [4, 1, 2, 1, 4, 3, 2, 1, 3, 2, 2]
        sage: jumps(w)
        [2, 4, 8, 10]
    """
    return [i+1 for i in range(len(w)-1) if w[i]<=w[i+1]]

def is_valid_col_word(w, m=None):
    """
    Determine if w is actually a valid column reading word of some 
    semistandard Young tableau with at most m columns

    If m is None, then we determine if w is a valid column reading word 
    of some semistandard Young tableau

    EXAMPLES:

        sage: w = [3, 2, 2, 1, 1]
        sage: is_valid_col_word(w)
        False

        sage: w = [3, 2, 1, 1, 1]
        sage: is_valid_col_word(w,3)
        True

        sage: w = [3, 2, 1, 1, 1]
        sage: is_valid_col_word(w,2)
        False

        sage: w = [3, 2, 1, 3, 1]
        sage: is_valid_col_word(w,2)
        True
    """
    from sage.combinat.tableau import Tableau
    J = [0]+jumps(w)+[len(w)]
    L = [w[J[i]:J[i+1]][::-1] for i in range(len(J)-1)]
    if all(len(L[i])>=len(L[i+1]) for i in range(len(L)-1)):
        if m == None or (len(jumps(w))<=m-1):
            T = Tableau(L)
            return T.conjugate().is_semistandard()
    return False

def list_equivalent_words(w):
    """
    Lists all words v equivalent to w in 0-Hecke monoid

    EXAMPLES:

        sage: w = [1,2,1]
        sage: list_equivalent_words(canonical_word(w,1))
        [(1, 1, 2, 1),
         (1, 2, 1, 1),
         (1, 2, 1, 2),
         (1, 2, 2, 1),
         (2, 1, 1, 2),
         (2, 1, 2, 1),
         (2, 1, 2, 2),
         (2, 2, 1, 2)]
        
        sage: list_equivalent_words([1,1,2,1])
        [(1, 1, 2, 1),
         (1, 2, 1, 1),
         (1, 2, 1, 2),
         (1, 2, 2, 1),
         (2, 1, 1, 2),
         (2, 1, 2, 1),
         (2, 1, 2, 2),
         (2, 2, 1, 2)]
    """
    if all(isinstance(i,(int,Integer)) for i in w):
        u = w
    else:
        raise ValueError("w needs to be a tuple of integers")
    V, queue = set([]), [tuple(u)]
    while len(queue) > 0:
        v = queue.pop(0)
        if tuple(v) not in V:
            V.add(tuple(v))
            L = applicable_relations(v)
            for pair in L:
                position, move = pair
                t = apply_relations(v,position,move)
                queue += [tuple(t)]
    return sorted([v for v in list(V)])

def apply_relations(word,position,move):
    """
    Applies a particular type of move using a relation in 0-Hecke monoid on 
    word at the specified position.

    EXAMPLES:

        sage: w = [2, 1, 3, 4]
        sage: apply_relations(w,position=1,move="pq=qp")
        [2, 3, 1, 4]
        
        sage: w = [1, 3, 2, 1, 2, 4]
        sage: apply_relations(w,position=2,move="pqp=qpq")
        [1, 3, 1, 2, 1, 4]
        
        sage: w = [2, 3, 1, 2, 2, 3]
        sage: apply_relations(w,position=3,move="pp=p")
        [2, 3, 1, 2, 3]

        sage: w = [2, 3, 1, 2, 3]
        sage: apply_relations(w,position=3,move="p=pp")
        [2, 3, 1, 2, 2, 3]
        
        sage: w = [2, 3, 1, 2, 2, 3]
        sage: apply_relations(w,position=2,move="pqq=ppq")
        [2, 3, 1, 1, 2, 3]
        
        sage: w = [2, 3, 1, 1, 2, 3]
        sage: apply_relations(w,position=2,move="ppq=pqq")
        [2, 3, 1, 2, 2, 3]
    """
    w = list(word[:])
    # Type 1
    if move == "pq=qp":
        if position > len(w)-2:
            raise IndexError("position is out of range for relation pq=qp")
        p,q = w[position],w[position+1]
        if abs(p-q)==1:
            raise IndexError("Relation pq=qp does not apply here")
        else:
            w[position] = q
            w[position+1] = p
    # Type 2
    elif move == "pqp=qpq":
        if position > len(w)-3:
            raise IndexError("position is out of range for relation pqp=qpq")
        p,q = w[position],w[position+1]
        if p != w[position+2]:
            raise IndexError("Relation pqp=qpq does not apply here")
        else:
            w[position] = q
            w[position+1] = p
            w[position+2] = q
    # Type 3
    elif move == "pqq=ppq":
        if position > len(w)-3:
            raise IndexError("position is out of range for relation pqq=ppq")
        p,q = w[position],w[position+2]
        if q != w[position+1]:
            raise IndexError("Relation pqq=ppq does not apply here")
        else:
            w[position+1] = p
    # Type 4
    elif move == "ppq=pqq":
        if position > len(w)-3:
            raise IndexError("position is out of range for relation ppq=pqq")
        p,q = w[position],w[position+2]
        if p != w[position+1]:
            raise IndexError("Relation ppq=pqq does not apply here")
        else:
            w[position+1] = q
    # Type 5
    elif move == "pp=p":
        if position > len(w)-2:
            raise IndexError("position is out of range for relation pp=p")
        p = w[position]
        if p != w[position+1]:
            raise IndexError("Relation pp=p does not apply here")
        else:
            w = w[:position+1]+w[position+2:]
    elif move == "p=pp":
        if position > len(w)-1:
            raise IndexError("position is out of range for relation p=pp")
        p = w[position]
        w = w[:position+1]+[p]+w[position+1:]
    return w

def applicable_relations(word):
    L = []
    for i in range(len(word)-1):
        if i < len(word)-2:
            p,q,r = word[i:(i+2)+1]
            if abs(p-q)>1:
                L += [[i,"pq=qp"]]
            if p==r and q!=p and abs(p-q)==1:
                L += [[i,"pqp=qpq"]]
            if q==p and r!=p:
                L += [[i,"ppq=pqq"]]
            if q==r and r!=p:
                L += [[i,"pqq=ppq"]]
        elif i == len(word)-2:
            p,q = word[i:(i+1)+1]
            if abs(p-q)>1:
                L += [[i,"pq=qp"]]
    return L

