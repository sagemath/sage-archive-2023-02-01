# -*- coding: utf-8 -*-

class AssociationScheme:

    def _is_association_scheme(self):

        #check matrix size
        if self._matrix.ncols() != self._nX:
            raise ValueError("matrix has wrong size")
        
        #check R0
        for i in range(self._nX):
            if self._matrix[i][i] != 0:
                raise ValueError("identity is not R_0")
            
        
        #check symmetries
        for i in range(self._nX):
            for j in range(i+1,self._nX):
                if self._matrix[i][j] != self._matrix[j][i]:
                    print("not symmetric")
            

        #check intersection numbers
        self.compute_intersection_numbers()
        r1 = self._r+1
        for i in range(r1):
            Ri = self.R(i)
            for j in range(r1):
                Rj = self.R(j)
                for k in range(r1):
                    Rk = self.R(k)
                    pijk = self.p(i,j,k)
                    for (x,y) in Rk:
                        count = 0
                        for z in self._X:
                            if (x,z) in Ri and (z,y) in Rj:
                                count += 1
                        if pijk != count:
                            raise ValueError("failed p{}{}{} (={}) with pair ({},{})".format(i,j,k,pijk,x,y))
                            return False
        return True
        
        
    def __init__(self, points, matrix, check=True):
        self._X = list(points)
        self._nX = len(self._X)
        self._XtoInt = { x: i for i,x in enumerate(self._X)}
        self._matrix = Matrix(matrix)

        #compute number of classes
        self._r = 0
        for r in self._matrix:
            for c in r:
                if c > self._r: self._r = c

        #empty intersection array
        r1 = self._r+1
        self._P = [ [ [None for k in range(r1) ] for j in range(r1)] for i in range(r1)]

        if check:
            self._is_association_scheme()
            

    def ground_set(self):
        return self._X

    def num_points(self):
        return self._nX

    def matrix(self):
        return self._matrix

    def num_classes(self):
        return self._r

    def has_relation(self,x,y,i):
        return self.which_relation(x,y)  == i

    def which_relation(self,x,y):
        return self._matrix[self._XtoInt[x]][self._XtoInt[y]]
    
    def R(self,k):
        Rk = set()
        for i in range(self._nX):
            for j in range(i,self._nX):
                if self._matrix[i][j] == k:
                    Rk.add((self._X[i],self._X[j]))
                    Rk.add((self._X[j],self._X[i]))

        return Rk

    def p(self,i,j,k):
        if self._P[i][j][k] is not None: return self._P[i][j][k]

        nX = self._nX

        self._P[i][j][k] = 0
        found = False
        for x in range(nX):
            for y in range(nX):
                if self._matrix[x][y] == k:
                    found = True
                    break
            if found:
                break

        if self._matrix[x][y] != k:
            raise RuntimeError("something bad happend in code")

        #now (x,y) in R_k
        for z in range(nX):
            if self._matrix[x][z] == i and self._matrix[z][y] == j:
                self._P[i][j][k] += 1

        return self._P[i][j][k]
        

    def compute_intersection_numbers(self):
        r1 = self._r+1
        for i in range(r1):
            for j in range(r1):
                for k in range(r1):
                    self.p(i,j,k)

    def is_pseudocyclic(self):
        #we need p_{ii}^0 to be a constant f for 1<= i <= r
        #we need sum_{i} p_{ii}^k = f-1 for 1 <= k <= r
        r1 = self._r+1
        f = self.p(1,1,0)
        for i in range(2,r1):
            if self.p(i,i,0) != f:
                return False

        for k in range(1,r1):
            sumP = 0
            for i in range(r1):
                sumP += self.p(i,i,k)

            if sumP != f-1:
                return False
        return True


def cyclotomic_scheme(const int q,const int r,check=True):
    r"""
    Returns an `r`-class association scheme on `q` points.
    The ground set is `GF(q)` and given a primite root `a` `(i,j) \in R_k` iff `j-i \in a^iK`
    where `K` is the subgroup of `r`th powers.

    This is a symmetric association scheme only if `\frac {q-1} r` is even or `q` is a power of 2.

    INPUT:

    - ``q`` -- prime power; number of points of the association scheme

    - ``r`` -- a positive divisor of `q-1`; number of classes of the association scheme

    - ``check`` -- if `True`, then the object constructed is checked to be an association scheme

    EXAMPLES:


    TESTS::


    """
    #for (q-1)/r even or q power of 2, then the association scheme is symmetric
    #and pseudocyclic
    if r <= 0 or (q-1)%r != 0:
        raise ValueError("we need r to be a (positive) divisor of q-1")

    Fq = GF(q)
    X = list(Fq)
    XtoInt = { x: i for i,x in enumerate(X) }
    
    relations = [ [0 for i in range(q)] for j in range(q)] #qxq matrix

    a = Fq.primitive_element()
    ar = a**r
    m = (q-1)//r
    K = [ ar**i for i in range(m)]
    for i in range(1,r+1):
        ai=a**i
        aiK = [ ai*x for x in K]
        for x in X:
            for z in aiK:
                sig_check()
                y = x+z
                relations[XtoInt[x]][XtoInt[y]] = i

    return AssociationScheme(X,relations,check=check)

def distance_regular_association_scheme(const int n, const int r, existence=False, check=True):
    r"""
    Returns an r-class  distance regular association scheme

    We say that an association scheme is distance regular if it is pseudocyclic and
    $p_{i+l mod r, j+l mod r}^{k+l mod r} = p_{i,j}^k for any i,j,k,l \in \{1,...,r\}$
    """
    def result(scheme):
        if check:
            if not scheme.is_pseudocyclic() or scheme.num_classes() != r or scheme.num_points() != n:
                raise RuntimeError("Sage built a wrong association scheme")
        return scheme
    
    if is_prime_power(n) and r > 0 and (n-1)%r == 0 and ( ((n-1)//r)%2 == 0 or n%2 == 0 ):
        #we hve cyclotomic
        if existence: return True
        return result(cyclotomic_scheme(n,r,check))
            
    if existence: return Unknown
    raise RuntimeError("Sage can't build a pseudocyclic association scheme with parameters ({},{})".format(n,r))
