#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#                          Gabriele Nebe <nebe@math.rwth-aachen.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import sage.misc.misc as misc
from sage.rings.arith import LCM
from sage.rings.real_mpfr import RealField
from sage.rings.finite_field import FiniteField
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.integer_ring import IntegerRing
from sage.rings.rational_field import RationalField
from sage.rings.integer import Integer

def Genus(A):
    """
    Given a nonsingular symmetric matrix A, return the genus of A.

    TODO: Implement some datastructure for quadratic forms and lattices.
    """
    return GenusSymbol_global_ring(A)

def LocalGenusSymbol(A,p):
    """
    Given a nonsingular symmetric matrix A, return the local symbol of A.
    """
    val = A.determinant().valuation(p)
    symbol = p_adic_symbol(A, p, val = val)
    return GenusSymbol_p_adic_ring(p, symbol)

def is_GlobalGenus(G):
    """
    Given a genus symbol G (specified by a collection of local symbols), return
    True in G represents the genus of a global quadratic form or lattice.
    """
    D = G.determinant()
    r, s = G.signature()
    oddity = r - s
    for loc in G._local_symbols:
        p = loc._prime
        sym = loc._symbol
        v = sum([ s[0]*s[1] for s in sym ])
        a = D // (p**v)
        b = Integer(misc.prod([ s[2] for s in sym ]))
        if p == 2:
            if not is_2_adic_genus(sym):
   	        # print "False in is_2_adic_genus(sym)"
	        return False
 	    if (a*b).kronecker(p) != 1:
	        # print "False in (%s*%s).kronecker(%s)"%(a,b,p)
                return False
            oddity -= loc.excess()
	else:
 	    if a.kronecker(p) != b:
	        # print "False in %s.kronecker(%s) != *%s"%(a,p,b)
                return False
            oddity += loc.excess()
    if oddity%8 != 0:
  	# print "False in oddity"
        return False
    return True

def is_2_adic_genus(symbol):
    """
    Given a 2-adic local symbols check whether it is the 2-adic symbol
    of a 2-adic form.
    """
    for s in symbol:
        if s[1] == 1:
	    if s[3] == 0 or s[2] != s[4]:
	        return False
        if s[1] == 2 and s[3] == 1:
	    if s[2] in (1,-1):
	       if not s[4] in (0,2,6):
	          return False
	    if s[2] in (3,-3):
	       if not s[4] in (2,4,6):
	          return False
        if (s[1] - s[4])% 2 == 1:
	    return False
	if s[3] == 0 and s[4] != 0:
	    return False
    return True

def canonical_2_adic_compartments(symbol):
    """
    See Conway-Sloane, pp. 382-383.
    """
    compartments = []
    i = 0
    r = len(symbol)
    while i < r:
        s = symbol[i]
        if s[3] == 1:
            v = s[0]
	    c = []
  	    while i < r and symbol[i][3] == 1 and symbol[i][0] == v:
	    	c.append(i)
 	        i += 1
                v += 1
            compartments.append(c)
        else:
	    i += 1
    return compartments

def canonical_2_adic_trains(symbol, compartments):
    """
    See Conway-Sloane, pp. 382-383.
    """
    trains = []
    i = 0
    while i < len(compartments):
        flag = True
        train = [ ]
        while flag:
	    ci = compartments[i]
            j = ci[0]
            if j == 0 or symbol[j-1][0] != symbol[j][0] - 1:
                train += ci
            else:
                train += [j-1] + ci
            act = ci[len(ci)-1]+1
            if i+1 < len(compartments):
                ci_plus = compartments[i+1]
            else:
                if act < len(symbol) and symbol[act][0] == symbol[act-1][0] +1:
                    train += [act]
	        flag = False
            if flag and symbol[ci[len(ci)-1]][0]+2 != symbol[ci_plus[0]][0]:
                if act != ci_plus[0] and symbol[act][0] == symbol[act-1][0] +1:
                    train += [act]
                flag = False
            i += 1
	trains.append(train)
    return trains

def canonical_2_adic_reduction(symbol):
    """
    """
    canonical_symbol = symbol
    # Canonical determinants:
    for i in range(len(symbol)):
        d = symbol[i][2]
	if d in (1,7):
            canonical_symbol[i][2] = 1
	else:
            canonical_symbol[i][2] = -1
    # Oddity fusion:
    compartments = canonical_2_adic_compartments(symbol)
    for compart in compartments:
        oddity = sum([ symbol[i][4] for i in compart ]) % 8
	for i in compart:
	    symbol[i][4] = 0
        symbol[compart[0]][4] = oddity
    #print "End oddity fusion:", canonical_symbol
    # Sign walking:
    trains = canonical_2_adic_trains(symbol, compartments)
    for train in trains:
        t = len(train)
	for i in range(t-1):
	    t1 = train[t-i-1]
	    if canonical_symbol[t1][2] == -1:
	        canonical_symbol[t1][2] = 1
	        canonical_symbol[t1-1][2] *= -1
		for compart in compartments:
		    if t1-1 in compart or t1 in compart:
		        o = canonical_symbol[compart[0]][4]
		        canonical_symbol[compart[0]][4] = (o+4) % 8
    #print "End sign walking:", canonical_symbol
    return canonical_symbol

def basis_complement(B):
    """
    Given an echelonized basis matrix (over a field),
    calculate a matrix whose rows form a basis complement.
    """
    F = B.parent().base_ring()
    m = B.nrows()
    n = B.ncols()
    C = MatrixSpace(F,n-m,n,sparse=True)(0)
    k = 0
    l = 0
    for i in range(m):
        for j in range(k,n):
             if B[i,j] == 0:
                 C[l,j] = 1
	         l += 1
	     else:
                 k = j+1
                 break
    for j in range(k,n):
	C[l+j-k,j] = 1
    return C

def signature(A):
    """
    The signature of a non-degenerate symmetric matrix.

    TODO: Implement a better algorithm.
    """
    r = 0
    s = 0
    e0 = 1
    for i in range(A.nrows()):
        # Argh!...
	e1 = RealField()(A[0:i+1, 0:i+1].determinant()).sign()
	if e0*e1 == 1:
	   r += 1
        else:
	   s += 1
        e0 = e1
    return (r, s)

def p_adic_symbol(A,p,val):
    """
    Given a symmetric matrix A and prime p, return the genus symbol at p.
    val = valuation of the maximal elementary divisor of A
    needed to obtain enough precision
    calculation is modulo p to the val+3
    TODO: Some description of the definition of the genus symbol.
    """
    if p % 2 == 0:
        return two_adic_symbol(A, val)
    m0 = min([ c.valuation(p) for c in A.list() ])
    q = p**m0
    n = A.nrows()
    A = MatrixSpace(IntegerRing(),n,n)([ c // q for c in A.list() ])
    A_p = MatrixSpace(FiniteField(p),n,n)(A)
    B_p = A_p.kernel().echelonized_basis_matrix()
    if B_p.nrows() == 0:
        e0 = Integer(A_p.det()).kronecker(p)
        n0 = A.nrows()
        return [ [m0,n0,e0] ]
    else:
        C_p = basis_complement(B_p)
        e0 = Integer((C_p*A_p*C_p.transpose()).det()).kronecker(p)
        n0 = C_p.nrows()
        sym = [ [0,n0,e0] ]
    r = B_p.nrows()
    B = MatrixSpace(IntegerRing(),r,n)(B_p)
    C = MatrixSpace(IntegerRing(),n-r,n)(C_p)
    # Construct the blocks for the Jordan decomposition [F,X;X,A_new]
    F = MatrixSpace(RationalField(),n-r,n-r)(C*A*C.transpose())
    U = F**-1
    d = LCM([ c.denominator() for c in U.list() ])
    R = IntegerRing().quotient_ring(Integer(p)**(val+3))
    u = R(d)**-1
    MatR = MatrixSpace(R,n-r,n-r)
    MatZ = MatrixSpace(IntegerRing(),n-r,n-r)
    U = MatZ(MatR(MatZ(U*d))*u)
    # X = C*A*B.transpose()
    # A = B*A*B.transpose() - X.transpose()*U*X
    X = C*A
    A = B*(A - X.transpose()*U*X)*B.transpose()
    return [ [s[0]+m0] + s[1:] for s in sym + p_adic_symbol(A, p, val) ]

def is_even(A):
    for i in range(A.nrows()):
        if A[i,i]%2 == 1:
	    return False, i
    return True, -1

def split_odd(A):
    """
    Given a Gram matrix A (mod 8), return a splitting [u] + B
    such that u is odd and B is not even.
    """
    n0 = A.nrows()
    if n0 == 1:
       return A[0,0], MatrixSpace(IntegerRing(),0,A.ncols())([])
    even, i = is_even(A)
    R = A.parent().base_ring()
    C = MatrixSpace(R,n0-1,n0)(0)
    u = A[i,i]
    for j in range(n0-1):
        if j < i:
	    C[j,j] = 1
	    C[j,i] = -A[j,i]*u
        else:
	    C[j,j+1] = 1
	    C[j,i] = -A[j+1,i]*u
        B = C*A*C.transpose()
    even, j = is_even(B)
    if even:
        I = A.parent()(1)
	# TODO: we could manually (re)construct the kernel here...
	if i == 0:
	    I[1,0] = 1 - A[1,0]*u
	    i = 1
	else:
	    I[0,i] = 1 - A[0,i]*u
	    i = 0
	A = I*A*I.transpose()
	u = A[i,i]
        C = MatrixSpace(R,n0-1,n0)(0)
	for j in range(n0-1):
	    if j < i:
	       C[j,j] = 1
	       C[j,i] = -A[j,i]*u
            else:
		C[j,j+1] = 1
		C[j,i] = -A[j+1,i]*u
	    B = C*A*C.transpose()
    even, j = is_even(B)
    if even:
        print "B:"
        print B
	assert False
    return u, B

def trace_diag(A):
    """
    Return the trace of the diagonalised form of A.
    """
    tr = 0
    while A.nrows() > 0:
       u, A = split_odd(A)
       tr += u
    return IntegerRing()(tr)


def two_adic_symbol(A, val):
    """
    Given a symmetric matrix A and prime p, return the genus symbol at p.

    val = valuation of maximal 2-elementary divisor

    The genus symbol of a component 2^m*f is of the form (m,n,s,d[,o]),
    where
        m = valuation of the component
        n = dimension of f
	d = det(f) in {1,3,5,7}
        s = 0 (or 1) if even (or odd)
	o = oddity of f (= 0 if s = 0) in Z/8Z
    """
    m0 = min([ c.valuation(2) for c in A.list() ])
    q = 2**m0
    A = A.parent()([ c // q for c in A.list() ])
    ZZ = IntegerRing()
    n = A.nrows()
    A_2 = MatrixSpace(FiniteField(2),n,n)(A)
    K_2 = A_2.kernel()
    R_8 = ZZ.quotient_ring(Integer(8))
    if K_2.dimension() == 0:
	A_8 = MatrixSpace(R_8,n)(A)
        n0 = A.nrows()
        # d0 = ZZ(A_8.determinant()) # no determinant over Z/8Z
        d0 = ZZ(R_8(MatrixSpace(ZZ,n)(A_8).determinant()))
	if d0 == 0:
	    print "A:"
	    print A
	    assert False
	even, i = is_even(A_2)
	if even:
	    return [ [m0,n0,d0,0,0] ]
	else:
	    tr8 = trace_diag(A_8)
	    return [ [m0,n0,d0,1,tr8] ]
    else:
        B_2 = K_2.echelonized_basis_matrix()
        C_2 = basis_complement(B_2)
	n0 = C_2.nrows()
	C = MatrixSpace(ZZ,n0,n)(C_2)
        A_new = C*A*C.transpose()
	# compute oddity modulo 8:
        A_8 = MatrixSpace(R_8,n0,n0)(A_new)
	# d0 = A_8.det() # no determinant over Z/8Z
        d0 = ZZ(R_8(MatrixSpace(ZZ,n0,n0)(A_8).determinant()))
	if d0 == 0:
	    print "A:"
	    print A_new
	    assert False
	even, i = is_even(A_new)
	if even:
	    sym = [ [0,n0,d0,0,0] ]
	else:
	    tr8 = trace_diag(A_8)
	    sym = [ [0,n0,d0,1,tr8] ]
    r = B_2.nrows()
    B = MatrixSpace(ZZ,r,n)(B_2)
    C = MatrixSpace(IntegerRing(),n-r,n)(C_2)
    F = MatrixSpace(RationalField(),n-r,n-r)(C*A*C.transpose())
    U = F**-1
    d = LCM([ c.denominator() for c in U.list() ])
    R = IntegerRing().quotient_ring(Integer(2)**(val+3))
    u = R(d)**-1
    MatR = MatrixSpace(R,n-r,n-r)
    MatZ = MatrixSpace(IntegerRing(),n-r,n-r)
    U = MatZ(MatR(MatZ(U*d))*u)
    X = C*A
    A = B*(A - X.transpose()*U*X)*B.transpose()
    return [ [s[0]+m0] + s[1:] for s in sym + two_adic_symbol(A, val) ]

def is_trivial_symbol(p, sym):
    if len(sym) != 1:
        return False
    if sym[0] != 0 or sym[2] != 1:
        return False
    if p != 2:
        return True
    return sym[3] == 1 and sym[1] % 8 == sym[4]

class Genus_Symbol_p_adic_ring(object):
    """
    Local genus symbol over a p-adic ring.
    """
    def __init__(self, prime, symbol, check = True):
        """
	Create the local genus symbol of given prime and local invariants.

	INPUT:
	     prime -- the prime
	     symbol -- the list of invariants for Jordan blocks A_t,...,A_t

        The genus symbol of a component p^m*A for odd prime = p is of the
        form (m,n,d), where

            m = valuation of the component
            n = rank of A
            d = det(A) in {1,u} for normalized quadratic non-residue u.

        The genus symbol of a component 2^m*A is of the form (m,n,s,d,o),
        where

            m = valuation of the component
            n = rank of A
            d = det(A) in {1,3,5,7}
            s = 0 (or 1) if even (or odd)
            o = oddity of A (= 0 if s = 0) in Z/8Z
	      = the trace of the diagonalization of A

	The genus symbol is a list of such symbols (ordered by m) for each
        of the Jordan blocks A_1,...,A_t.

	Reference: Conway and Sloane, Chapter 9.
	"""
	if check:
	   pass
	self._prime = prime
	self._symbol = symbol
	self._canonical_symbol = None

    def __repr__(self):
        return "Genus symbol at %s : %s"%(self._prime, self._symbol)

    def __eq__(self,other):
        p = self._prime
	if p != other._prime:
	    return False
        return self.canonical_symbol() == other.canonical_symbol()

    def __ne__(self,other):
        return not self.__eq__(other)

    def canonical_symbol(self):
        symbol = self._symbol
        if self._prime == 2:
	    if self._canonical_symbol is None:
  	        self._canonical_symbol = canonical_2_adic_reduction(symbol)
	    return self._canonical_symbol
	else:
	    return self._symbol

    def symbol(self):
        return self._symbol

    def number_of_blocks(self):
        return len(self._symbol)

    def determinant(self):
        p = self._prime
        return misc.prod([ p**(s[0]*s[1]) for s in self._symbol ])

    def rank(self):
        return sum([ s[1] for s in self._symbol ])

    def dimension(self):
        return self.rank()

    def excess(self):
        """
	The p-excesss in the notation of Conway & Sloane, and the oddity for p = 2.
        """
	p = self._prime
        if self._prime == 2:
	   k = 0
	   for s in self._symbol:
	       if s[0]%2 == 1 and s[2] in (3,5):
	           k += 1
           return Integer(sum([ s[4] for s in self._symbol ]) + 4*k).mod(8)
        else:
	   k = 0
	   for s in self._symbol:
	       if s[0]%2 == 1 and s[2] == -1:
	           k += 1
           return Integer(sum([ s[1]*(p**s[0]-1) for s in self._symbol ]) + 4*k).mod(8)

    def trains(self):
        assert self._prime == 2
	symbol = self._symbol
        compartments = canonical_2_adic_compartments(symbol)
        return canonical_2_adic_trains(symbol, compartments)

    def compartments(self):
        assert self._prime == 2
	symbol = self._symbol
        return canonical_2_adic_compartments(symbol)


class Genus_Symbol_global_ring(object):
    """
    The genus of an integral lattice.
    """
    def __init__(self, A, max_elem_divisors=None):
        """
	Input precision max_elem_divisors for valuation of maximal p-elementary divisor.
        """
        D = A.determinant()
        D = 2*D
        prms = [ p[0] for p in D.factor() ]
	self._representative = A
	self._signature = signature(A)
	self._local_symbols = []
	for p in prms:
	    if max_elem_divisors is None:
	        val = D.valuation(p)
	    symbol = p_adic_symbol(A, p, val = val)
 	    G = Genus_Symbol_p_adic_ring(p, symbol)
      	    self._local_symbols.append(G)

    def __repr__(self):
        return "Genus of %s"%self._representative

    def __eq__(self, other):
        if self is other:
	    return True
	t = len(self._local_symbols)
        if t != len(other._local_symbols):
	    return False
	for i in range(t):
	    if self._local_symbols[i] != other._local_symbols[i]:
	        return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def signature(self):
        return self._signature

    def determinant(self):
        r, s = self.signature()
        return (-1)**s*misc.prod([ G.determinant() for G in self._local_symbols ])
