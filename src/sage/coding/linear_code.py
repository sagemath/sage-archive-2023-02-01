r"""
Linear Codes

VERSION: 0.2

AUTHOR:
    -- David Joyner (2005-11-22, 2006-12-03): written
    -- William Stein (2006-01-23) -- Inclusion in SAGE
    -- David Joyner (2006-01-30, 2006-04): small fixes

This file contains
\begin{enumerate}
\item LinearCode, Codeword class definitions
\item The spectrum (weight distribution) and minimum distance
    programs (calling Steve Linton's C programs)
\item interface with A. Brouwer's online tables
\item wrapped GUAVA's HammingCode, RandomLinearCode,
    Golay codes, binary Reed-Muller code
\item gen_mat, check_mat, decode, dual_code method for LinearCode.
\end{enumerate}


    EXAMPLES:
        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G = MS([[1,1,1,0,0,0,0], [ 1, 0, 0, 1, 1, 0, 0], [ 0, 1, 0, 1, 0, 1, 0], [1, 1, 0, 1, 0, 0, 1]])
        sage: C = LinearCode(G)
        sage: C.basis()
        [(1, 1, 1, 0, 0, 0, 0),
         (1, 0, 0, 1, 1, 0, 0),
         (0, 1, 0, 1, 0, 1, 0),
         (1, 1, 0, 1, 0, 0, 1)]
        sage: c = C.basis()[1]
        sage: c in C
        True
        sage: c.nonzero_positions()
        [0, 3, 4]
        sage: c.support()
        [0, 3, 4]
        sage: c.parent()
        Vector space of dimension 7 over Finite Field of size 2

To be added:
\begin{enumerate}
\item PermutedCode method (with PermutationGroupElement as argument).
\item More wrappers
\item automorphism group.
\item cyclic codes
\item GRS codes and special decoders.
\item $P^1$ Goppa codes and group actions.
\end{enumerate}

"""

#*****************************************************************************
#       Copyright (C) 2005 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.databases.lincodes import *
from sage.interfaces.all import gap
from sage.misc.preparser import *
from sage.matrix.matrix_space import *
from sage.rings.finite_field import *
from sage.misc.sage_eval import sage_eval

VectorSpace = fm.VectorSpace

###################### coding theory functions ##############################

def wtdist(Gmat, F):
    """
    INPUT:
        Gmat -- a string representing a GAP generator matrix G of a linear code.
        F -- a (SAGE) finite field - the base field of the code.

    OUTPUT:
        Returns the spectrum of the associated code.

    EXAMPLES:
        sage: Gstr = 'Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]'
        sage: F = GF(2)
    	sage: sage.coding.linear_code.wtdist(Gstr, F)
    	[1, 0, 0, 7, 7, 0, 0, 1]

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    ALGORITHM:
        Uses C programs written by Steve Linton in the kernel
        of GAP, so is fairly fast.

    AUTHOR: David Joyner (2005-11)
    """
    q = F.order()
    G = gap(Gmat)
    k = gap(F)
    C = G.GeneratorMatCode(k)
    n = int(C.WordLength())
    z = 'Z(%s)*%s'%(q, [0]*n)     # GAP zero vector as a string
    dist = gap.eval("w:=DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    #d = G.DistancesDistributionMatFFEVecFFE(k, z)
    v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,n+2)]
    return v

def min_wt_vec(Gmat,F):
    """
    Uses C programs written by Steve Linton in the kernel of GAP, so is fairly fast.

    INPUT:
        Same as wtdist.

    OUTPUT:
        Returns a minimum weight vector v, the "message" vector m such that m*G = v,
    and the (minimum) weight, as a triple.

    EXAMPLES:
        sage: Gstr = "Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]"
    	sage: sage.coding.linear_code.min_wt_vec(Gstr,GF(2))
    	[[0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 0], 3]

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    AUTHOR: David Joyner (11-2005)
    """
    gap.eval("G:="+Gmat)
    k = int(gap.eval('Length(G)'))
    q = F.order()
    qstr = str(q)
    gap.eval("C:=GeneratorMatCode("+Gmat+",GF("+qstr+"))")
    n = int(gap.eval("WordLength(C)"))
    zerovec = [0 for i in range(n)]
    zerovecstr = "Z("+qstr+")*"+str(zerovec)
    all = []
    for i in range(1,k+1):
        P = gap.eval("P:=AClosestVectorCombinationsMatFFEVecFFECoords("+Gmat+", GF("+qstr+"),"+zerovecstr+","+str(i)+","+str(0)+"); d:=WeightVecFFE(P[1])")
        v = gap.eval("v:=List(P[1], i->i)")
        m = gap.eval("m:=List(P[2], i->i)")
        dist = gap.eval("d")
        #print v,m,dist
        #print [gap.eval("v["+str(i+1)+"]") for i in range(n)]
        all.append([[gap_to_sage(gap.eval("v["+str(i+1)+"]"),F)
                              for i in range(n)],
                    [gap_to_sage(gap.eval("m["+str(i+1)+"]"),F)
                              for i in range(k)],int(dist)])
    ans = all[0]
    for x in all:
        if x[2]<ans[2] and x[2]>0:
            ans = x
    return ans

def minimum_distance_lower_bound(n,k,F):
        """
        Connects to http://www.win.tue.nl/~aeb/voorlincod.html
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven,
        via Steven Sivek's linear_code_bound.

        EXAMPLES:
            sage: sage.coding.linear_code.minimum_distance_upper_bound(7,4,GF(2))     # optional (net connection)
            3

        Obviously requires an internet connection.
        """
        q = F.order()
        bounds = linear_code_bound(q,n,k)
        return bounds[0]

def minimum_distance_upper_bound(n,k,F):
        """
        Connects to http://www.win.tue.nl/~aeb/voorlincod.html
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven
        via Steven Sivek's linear_code_bound.

        EXAMPLES:
            sage: sage.coding.linear_code.minimum_distance_upper_bound(7,4,GF(2))  # optional (net connection)
            3

        Obviously requires an internet connection.
        """
        q = F.order()
        bounds = linear_code_bound(q,n,k)
        return bounds[1]

def minimum_distance_why(n,k,F):
        """
        Connects to http://www.win.tue.nl/~aeb/voorlincod.html
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven
        via Steven Sivek's linear_code_bound.

        EXAMPLES:
            sage: sage.coding.linear_code.minimum_distance_why(7,4,GF(2))  # optional (net connection)
            Lb(7,4) = 3 is found by truncation of:
            Lb(8,4) = 4 is found by the (u|u+v) construction
            applied to [4,3,2] and [4,1,4]-codes
            Ub(7,4) = 3 follows by the Griesmer bound.

        Obviously requires an internet connection.
        """
        q = F.order()
        bounds = linear_code_bound(q,n,k)
        print bounds[2]


########################### linear codes python class #######################

class LinearCode(module.Module):
    """
    A class for linear codes over a finite field or finite ring.

    INPUT:
        G -- A $k\times n$ matrix of (full) rank $k$, $k\leq n$, over
             a finite field $F$.

    OUTPUT:
        The linear code of length $n$ over $F$ having $G$ as a
        generator matrix.

    EXAMPLES:
        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [ 1, 0, 0, 1, 1, 0, 0], [ 0, 1, 0, 1, 0, 1, 0], [1, 1, 0, 1, 0, 0, 1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C.minimum_distance_upper_bound()   # optional (net connection)
        3
        sage: C.base_ring()
        Finite Field of size 2
        sage: C.dimension()
        4
        sage: C.length()
        7
        sage: C.minimum_distance()
        3
        sage: C.spectrum()
        [1, 0, 0, 7, 7, 0, 0, 1]
        sage: C.weight_distribution()
        [1, 0, 0, 7, 7, 0, 0, 1]
        sage: C.minimum_distance_why()     # optional (net connection)
        Ub(7,4) = 3 follows by the Griesmer bound.

    AUTHOR: David Joyner (11-2005)
    """
    def __init__(self, gen_mat):
        self.__gens = gen_mat.rows()
        self.__gen_mat = gen_mat
        self.__base_ring = gen_mat[0][0].parent()
        self.__length = len(gen_mat[1])
        self.__dim = gen_mat.rank()

    def length(self):
        return self.__length

    def dimension(self):
        return self.__dim

    def base_ring(self):
        return self.__base_ring

    def _repr_(self):
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def gen_mat(self):
        return self.__gen_mat

    def gens(self):
        return self.__gens

    def basis(self):
        return self.__gens

    def ambient_space(self):
        return VectorSpace(self.__base_ring,self.__length)

    def __contains__(self,v):
        A = self.ambient_space()
        C = A.subspace(self.gens())
        return C.__contains__(v)

    def characteristic(self):
        return (self.base_ring()).characteristic()

    def minimum_distance_upper_bound(self):
        """
        Connects to http://www.win.tue.nl/~aeb/voorlincod.html
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

        Obviously requires an internet connection
        """
        q = (self.base_ring()).order()
        n = self.length()
        k = self.dimension()
        bounds = linear_code_bound(q,n,k)
        return bounds[1]

    def minimum_distance_why(self):
        """
        Connects to http://www.win.tue.nl/~aeb/voorlincod.html
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

        Obviously requires an internet connection.
        """
        q = (self.base_ring()).order()
        n = self.length()
        k = self.dimension()
        bounds = linear_code_bound(q,n,k)
        lines = bounds[2].split("\n")
        for line in lines:
            if len(line)>0:
                if line[0] == "U":
                    print line

    def minimum_distance(self):
        """
        Uses a GAP kernel function (in C) written by Steve Linton.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3
            sage: C=RandomLinearCode(10,5,GF(4))
            sage: C.gen_mat()                ## random
	    [    1     0     0     0     0 x + 1     1     0     0     0]
	    [x + 1     1     0     1     0 x + 1     1     1     0     0]
	    [    0 x + 1     0 x + 1     0     x x + 1 x + 1 x + 1     0]
	    [    1     0     x     0     1     0     0     0     0     1]
	    [    0     0     1     1     0     0     0     0     x x + 1]
            sage: C.minimum_distance()       ## random
            2
            sage: C.minimum_distance_upper_bound()  # optional (net connection)
            5
            sage: C.minimum_distance_why()          # optional (net connection)
            Ub(10,5) = 5 follows by the Griesmer bound.


        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        Gstr = str(gap(G))+"*Z("+str(q)+")^0"
        return min_wt_vec(Gstr,F)[2]


    def spectrum(self):
        """
        Uses a GAP kernel function (in C) written by Steve Linton.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [ 1, 0, 0, 1, 1, 0, 0], [ 0, 1, 0,1, 0, 1, 0], [1, 1, 0, 1, 0, 0, 1]])
            sage: C = LinearCode(G)
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]

        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        Glist = [list(x) for x in G]
        Gstr = "Z("+str(q)+")*"+str(Glist)
        spec = wtdist(Gstr,F)
        return spec

    def weight_distribution(self):
        #same as spectrum
        return self.spectrum()

    def __cmp__(self, right):
        raise NotImplementedError

    def decode(self, right):
        """
        Wraps GUAVA's Decodeword. Hamming codes have a special
        decoding algorithm. Otherwise, syndrome decoding is used.

        INPUT:
            right must be a vector of length = length(self)

        OUTPUT:
            The codeword c in C closest to r.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: MS = MatrixSpace(GF(2),1,7)
            sage: F=GF(2); a=F.gen()
            sage: v=MS([a,a,F(0),a,a,F(0),a]); v
            [1 1 0 1 1 0 1]
            sage: C.decode(v)
            [1, 1, 0, 1, 0, 0, 1]

        Does not work for very long codes since the syndrome table grows too large.
        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        Gstr = str(gap(G))
        vstr = str(gap(right))
        v = vstr[1:-1]
        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        gap.eval("c:=VectorCodeword(Decodeword( C, Codeword( "+v+" )))")
        ans = [gap_to_sage(gap.eval("c["+str(i)+"]"),F) for i in range(1,n+1)]
        return ans

    def dual_code(self):
        """
        Wraps GUAVA's DualCode.

        OUTPUT:
            The dual code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.dual_code()
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C = HammingCode(3,GF(4))
            sage: C.dual_code()
            Linear code of length 21, dimension 3 over Finite Field in a of size 2^2
        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        Gstr = str(gap(G))
	C = gap.GeneratorMatCode(Gstr, 'GF(%s)'%q)
        H = C.CheckMat()
        A = H._matrix_(GF(q))
        return LinearCode(A)

        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        Hmat = gap.eval("H:=CheckMat( C )")
        H = [[gap_to_sage(gap.eval("H["+str(i)+"]["+str(j)+"]"),F)
              for j in range(1,n+1)] for i in range(1,n-k+1)]
        MS = MatrixSpace(F,n-k,n)
        return LinearCode(MS(H))

    def check_mat(self):
        """
        Returns the check matrix of self.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: Cperp = C.dual_code()
            sage: C; Cperp
            Linear code of length 7, dimension 4 over Finite Field of size 2
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C.gen_mat()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: C.check_mat()
            [0 1 1 1 1 0 0]
            [1 0 1 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: Cperp.check_mat()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
            sage: Cperp.gen_mat()
            [0 1 1 1 1 0 0]
            [1 0 1 1 0 1 0]
            [1 1 0 1 0 0 1]
        """
        Cperp = self.dual_code()
        return Cperp.gen_mat()

######### defining the Codeword class by copying the FreeModuleElement class:
Codeword = fme.FreeModuleElement
Codeword.support = fme.FreeModuleElement.nonzero_positions
is_Codeword = fme.is_FreeModuleElement


##################### wrapped GUAVA functions ############################

def HammingCode(r,F):
    """
    INPUT:
        r, F -- r>1 an integer, F a finite field.

    OUTPUT:
        Returns the $r^{th}$ Hamming code over $F=GF(q)$ of length $n=(q^r-1)/(q-1)$.

    EXAMPLES:
        sage: C = HammingCode(3,GF(3))
        sage: C
        Linear code of length 13, dimension 10 over Finite Field of size 3
        sage: C.minimum_distance()
        3
        sage: C.gen_mat()
        [2 2 1 0 0 0 0 0 0 0 0 0 0]
        [1 2 0 1 0 0 0 0 0 0 0 0 0]
        [2 0 0 0 2 1 0 0 0 0 0 0 0]
        [1 0 0 0 2 0 1 0 0 0 0 0 0]
        [0 2 0 0 2 0 0 1 0 0 0 0 0]
        [2 2 0 0 2 0 0 0 1 0 0 0 0]
        [1 2 0 0 2 0 0 0 0 1 0 0 0]
        [0 1 0 0 2 0 0 0 0 0 1 0 0]
        [2 1 0 0 2 0 0 0 0 0 0 1 0]
        [1 1 0 0 2 0 0 0 0 0 0 0 1]
        sage: C = HammingCode(3,GF(4))
        sage: C
        Linear code of length 21, dimension 18 over Finite Field in a of size 2^2

    AUTHOR: David Joyner (11-2005)
    """
    q = F.order()
    gap.eval("C:=HammingCode("+str(r)+", GF("+str(q)+"))")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def QuadraticResidueCode(n,F):
    """
    INPUT:
        n, F -- n>2 a prime and F a finite prime field F of order q.
                Moreover, q must be a quadratic residue modulo n.

    OUTPUT:
        Returns a quadratic residue code. Its generator polynomial
        is the product of the polynomials $x-\alpha^i$
        ($\alpha$ is a primitive $n^{th}$ root of unity,
        and $i$ is an integer in the set of quadratic residues
        modulo $n$).

    EXAMPLES:
        sage: C = QuadraticResidueCode(7,GF(2))
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: C = QuadraticResidueCode(17,GF(2))
        sage: C
        Linear code of length 17, dimension 9 over Finite Field of size 2

    AUTHOR: David Joyner (11-2005)
    """
    q = F.order()
    gap.eval("C:=QRCode("+str(n)+", GF("+str(q)+"))")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def QuasiQuadraticResidueCode(p):
    """
    INPUT:
        p -- a prime >2.

    OUTPUT:
        Returns a (binary) quasi-quadratic residue code, as defined by
        Proposition 2.2 in Bazzi-Mittel ({\it Some constructions of
        codes from group actions}, (preprint March 2003). Its
        generator matrix has the block form $G=(Q,N)$.
        Here $Q$ is a $p\times p$ circulant matrix whose top row
        is $(0,x_1,...,x_{p-1})$, where $x_i=1$ if and only if $i$
        is a quadratic residue $\mod p$, and $N$ is a $p\times p$
        circulant matrix whose top row is $(0,y_1,...,y_{p-1})$, where
        $x_i+y_i=1$ for all i. (In fact, this matrix can be recovered
        as the component DoublyCirculant of the code.)

    EXAMPLES:
        sage: C = QuasiQuadraticResidueCode(11)
        sage: C
        Linear code of length 22, dimension 11 over Finite Field of size 2

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=QQRCode("+str(p)+")")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def BinaryReedMullerCode(r,k):
    """
    INPUT:
        r, k -- positive integers with $2^k>r$.

    OUTPUT:
        Returns a binary 'Reed-Muller code' with dimension k and
        order r. This is a code with length $2^k$ and minimum
        distance $2^k-r$ (see for example, section 1.10 in
        Huffman-Pless {\it Fundamentals of Coding Theory}).
        By definition, the $r^{th}$ order binary Reed-Muller code of
        length $n=2^m$, for $0 \leq r \leq m$, is the set of
        all vectors $(f(p)\ |\ p in GF(2)^m)$, where $f$ is a
        multivariate polynomial of degree at most $r$ in $m$ variables.

    EXAMPLE:
        sage: C = BinaryReedMullerCode(2,4)
        sage: C
        Linear code of length 16, dimension 11 over Finite Field of size 2
        sage: C.minimum_distance()
        4
        sage: C.gen_mat()
	[1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
	[0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1]
	[0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1]
	[0 0 1 1 0 0 1 1 0 0 1 1 0 0 1 1]
	[0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1]
	[0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1]
	[0 0 0 0 0 0 0 0 0 0 1 1 0 0 1 1]
	[0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1]
	[0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1]
	[0 0 0 0 0 1 0 1 0 0 0 0 0 1 0 1]
	[0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1]

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=ReedMullerCode("+str(r)+", "+str(k)+")")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def BinaryGolayCode():
    """
    OUTPUT:
        BinaryGolayCode returns a binary Golay code. This is a
        perfect [23,12,7] code. It is also cyclic, and has
        generator polynomial $g(x)=1+x^2+x^4+x^5+x^6+x^{10}+x^{11}$.
        Extending it results in an extended Golay code (see
        ExtendedBinaryGolayCode).

    EXAMPLE:
        sage: C = BinaryGolayCode()
        sage: C
        Linear code of length 23, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()               # long
        7

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=BinaryGolayCode()")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def ExtendedBinaryGolayCode():
    """
    OUTPUT:
        BinaryGolayCode returns the extended binary Golay code. This
        is a perfect [24,12,8] code. This code is self-dual.

    EXAMPLES:
        sage: C = ExtendedBinaryGolayCode()
        sage: C
        Linear code of length 24, dimension 12 over Finite Field of size 2
        sage: C.minimum_distance()               # long
        8

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(2)
    gap.eval("C:=ExtendedBinaryGolayCode()")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def TernaryGolayCode():
    """
    OUTPUT:
        TernaryGolayCode returns a ternary Golay code. This is a
        perfect [11,6,5] code. It is also cyclic, and has generator
        polynomial $g(x)=2+x^2+2x^3+x^4+x^5$.

    EXAMPLES:
        sage: C = TernaryGolayCode()
        sage: C
        Linear code of length 11, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        5

    AUTHOR: David Joyner (11-2005)
    """
    F = GF(3)
    gap.eval("C:=TernaryGolayCode()")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def ExtendedTernaryGolayCode():
    """
    OUTPUT:
        ExtendedTernaryGolayCode returns a ternary Golay code.
        This is a self-dual perfect [12,6,6] code.

    EXAMPLES:
        sage: C = ExtendedTernaryGolayCode()
        sage: C
        Linear code of length 12, dimension 6 over Finite Field of size 3
        sage: C.minimum_distance()
        6
        sage: C.gen_mat()
	[1 0 2 1 2 2 0 0 0 0 0 1]
	[0 1 0 2 1 2 2 0 0 0 0 1]
	[0 0 1 0 2 1 2 2 0 0 0 1]
	[0 0 0 1 0 2 1 2 2 0 0 1]
	[0 0 0 0 1 0 2 1 2 2 0 1]
	[0 0 0 0 0 1 0 2 1 2 2 1]

    AUTHOR: David Joyner (11-2005)
    """
    F = FiniteField(3)
    gap.eval("C:=ExtendedTernaryGolayCode()")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def RandomLinearCode(n,k,F):
    """
    INPUT:
        Integers n,k, with n>k>1.

    OUTPUT:
        Returns a "random" linear code with length n,
        dimension k over field F. The method used is to first
        construct a $k\times n$ matrix of the block form $(I,A)$,
        where $I$ is a $k\times k$ identity matrix and $A$ is a
        $k\times (n-k)$ matrix constructed using random elements
        of $F$. Then the columns are permuted
        using a randomly selected element of SymmetricGroup(n).

    EXAMPLES:
        sage: C = RandomLinearCode(30,15,GF(2))
        sage: C                                        # random output
        Linear code of length 30, dimension 15 over Finite Field of size 2
        sage: C = RandomLinearCode(10,5,GF(4))
        sage: C                                       # random output
        Linear code of length 10, dimension 5 over Finite Field in x of size 2^2

    AUTHOR: David Joyner (11-2005)
    """
    q = F.order()
    gap.eval("C:=RandomLinearCode("+str(n)+","+str(k)+", GF("+str(q)+"))")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))





