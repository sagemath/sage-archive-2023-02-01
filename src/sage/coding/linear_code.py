r"""
Linear Codes

VERSION: 0.6

Let $ F$ be a finite field (we denote the finite field with $q$ elements
$GF(q)$ by $\FF_q$). A subspace of $ F^n$ (with the standard basis)
is called a {\it linear code} of length $ n$. If its
dimension is denoted $k$ then we typically store a basis
of $C$ as a $k\\times  n$ matrix (the rows are the basis vectors)
called the {\it generator matrix} of $C$.
The rows of the {\it parity check matrix} of $C$ are a basis
for the code,

\[ C^* = \{ v \in GF(q)^n\ |\ v\cdot c = 0,\ for \ all\ c \in C \},
\]
called the {\it dual space} of  $C$.

If $ F=\FF_2$ then $ C$ is called a {\it binary code}.
If $ F = \FF_q$ then $ C$ is called a {\it $ q$-ary code}.
The elements of a code $ C$ are called {\it codewords}.

Let $ F$ be a finite field with $ q$ elements. Here's a constructive
definition of a cyclic code of length $ n$.

\begin{enumerate}
\item
      Pick a monic polynomial $ g(x)\in F[x]$ dividing $ x^n-1$.
      This is called the {\it generating polynomial} of the code.
\item
      For each polynomial $ p(x)\in F[x]$, compute
      $p(x)g(x)\ ({\rm mod}\ x^n-1). $
      Denote the answer by $ c_0+c_1x+...+c_{n-1}x^{n-1}$.
\item
      $ {\bf c} =(c_0,c_1,...,c_{n-1})$ is a codeword in $ C$. Every
      codeword in $ C$ arises in this way (from some $ p(x)$).
\end{enumerate}
The {\it polynomial notation} for the code is to call
$ c_0+c_1x+...+c_{n-1}x^{n-1}$ the codeword (instead of
$ (c_0,c_1,...,c_{n-1})$). The polynomial $h(x)=(x^n-1)/g(x)$
is called the {\it check polynomial} of $C$.

Let $ n$ be a positive integer relatively prime to $ q$ and
let $ \alpha$ be a primitive $n$-th root of unity. Each generator
polynomial $g$ of a cyclic code $C$ of length $n$ has a factorization
of the form

\[
g(x) = (x - \alpha^{k_1})...(x - \alpha^{k_r}),
\]
where $ \{k_1,...,k_r\} \subset \{0,...,n-1\}$. The numbers
$ \alpha^{k_i}$, $ 1 \leq i \leq r$, are called the {\it zeros}
of the code $ C$. Many families of cyclic codes (such as
the quadratic residue codes) are defined using properties of the
zeros of $C$.

The symmetric group $S_n$ acts on $F^n$ by permuting coordinates.
If an element $p\\in S_n$ sends a code $C$ of length $n$ to itself
(in other words, every codeword of $C$ is sent to some other codeword
of $C$) then $p$ is called a {\it permutation automorphism} of $C$.
The (permutation) automorphism group is denoted $Aut(C)$.

AUTHOR:
    -- David Joyner (2005-11-22, 2006-12-03): written
    -- William Stein (2006-01-23) -- Inclusion in SAGE
    -- David Joyner (2006-01-30, 2006-04): small fixes
    -- David Joyner (2006-07): added documentation, group-theoretical methods,
       ExtendedQuadraticResidueCode, ToricCode
    -- David Joyner (2006-08): hopeful latex fixes to documention, added list
       and __iter__ methods to LinearCode and examples, added hamming_weight
       function, fixed random method to return a vector, TrivialCode,
       fixed subtle bug in dual_code, added galois_closure method,
       fixed mysterious bug in permutation_automorphism_group (GAP
       was over-using "G" somehow?)
    -- David Joyner (2006-08): hopeful latex fixes to documention,
       added CyclicCode, best_known_linear_code, bounds_minimum_distance.

This file contains
\begin{enumerate}
\item
    LinearCode, Codeword class definitions
\item
    The spectrum (weight distribution) and minimum distance
    programs (calling Steve Linton's C programs)
\item
    interface with A. Brouwer's online tables, as well as
    best_known_linear_code, bounds_minimum_distance which call tables
    in GUAVA (updated May 2006) created by Cen Tjhai instead of the online
    internet tables.
\item
    Hamming codes, "random" linear codes, Golay codes, binary Reed-Muller codes,
    binary quadratic and extended quadratic residue codes, ToricCode,
    TrivialCode, cyclic codes.
\item
    gen_mat, list, check_mat, decode, dual_code methods for LinearCode,
\item
    permutation methods:
    is_permutation_automorphism, permutation_automorphism_group,
    permuted_code, standard_form, module_composition_factors.
\end{enumerate}


    EXAMPLES:
        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
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
\item More wrappers
\item GRS codes and special decoders.
\item $P^1$ Goppa codes and group actions.
\end{enumerate}

REFERENCES:
   [HP] W. C. Huffman and V. Pless, Fundamentals of error-correcting codes,
        Cambridge Univ. Press, 2003.
   [Gu] GUAVA manual, http://www.gap-system.org/Packages/guava.html

"""

#*****************************************************************************
#       Copyright (C) 2005 David Joyner <wdj@usna.edu>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.databases.lincodes import linear_code_bound
from sage.interfaces.all import gap
from sage.misc.preparser import *
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_field import *
from sage.groups.perm_gps.permgroup import *
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod
from sage.misc.functional import log

VectorSpace = fm.VectorSpace

###################### coding theory functions ##############################

def hamming_weight(v):
    return len(v.nonzero_positions())

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
        Uses C programs written by Steve Linton in the kernel of GAP, so is fairly fast.

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
    r"""
    Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
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
    r"""
    Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
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
    r"""
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

def best_known_linear_code(n,k,F):
    """
    best_known_linear_code returns the best known (as of 11 May 2006) linear code of
    length n, dimension k over field F. The function uses the tables described in bounds_minimum_distance
    to construct this code.

    This does not require an internet connection.

    EXAMPLES:
        sage: best_known_linear_code(10,5,GF(2))
        'a linear [10,5,4]2..4 shortened code'

    This means that best possible binary linear code of length 10 and dimension 5
    is a code with minimum distance 4 and covering radius somewhere between 2 and 4.
    Use "minimum_distance_why(10,5,GF(2))" or "print bounds_minimum_distance(10,5,GF(2))"
    for further details.
    """
    q = F.order()
    return gap.eval("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))

def bounds_minimum_distance(n,k,F):
    """
    The function bounds_minimum_distance calculates a lower and upper bound for the minimum
    distance of an optimal linear code with word length n, dimension k over field F. The
    function returns a record with the two bounds and an explanation for each bound. The
    function Display can be used to show the explanations.

    The values for the lower and upper bound are obtained from a table constructed by Cen Tjhai for GUAVA,
    derived from the table of Brouwer. (See http://www.win.tue.nl/~aeb/voorlincod.html or use the SAGE function
    minimum_distance_why for the most recent data.) These tables contain lower and upper bounds for
    q=2 (n <= 257), 3 (n <= 243), 4 (n <= 256). (Current as of 11 May 2006.)
    For codes over other fields and for larger word lengths, trivial bounds are used.

    This does not require an internet connection. The format of the output is a
    little non-intuitive. Try print bounds_minimum_distance(10,5,GF(2)) for example.
    """
    q = F.order()
    gap.eval("data := BoundsMinimumDistance(%s,%s,GF(%s))"%(n,k,q))
    Ldata = gap.eval("Display(data)")
    return Ldata

########################### linear codes python class #######################

class LinearCode(module.Module):
    r"""
    A class for linear codes over a finite field or finite ring.
    Each instance is a linear code determined by a generator matrix $G$
    (i.e., a k x n matrix of (full) rank $k$, $k\leq n$ over a finite field $F$.

    INPUT:
        G -- a generator matrix over $F$. (G can be defined over a finite ring but
             the matrices over that ring must have certain attributes, such as "rank".)

    OUTPUT:
        The linear code of length $n$ over $F$ having $G$ as a generator matrix.

    EXAMPLES:
        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
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
        sage: MS = MatrixSpace(IntegerModRing(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Ring of integers modulo 5

    AUTHOR: David Joyner (11-2005)
    """
    def __init__(self, gen_mat):
        self.__gens = gen_mat.rows()
        self.__gen_mat = gen_mat
        self.__base_ring = gen_mat[0][0].parent()
        self.__length = len(gen_mat[0])
        self.__dim = gen_mat.rank()

    def length(self):
        return self.__length

    def dimension(self):
        return self.__dim

    def base_ring(self):
        return self.__base_ring

    def _repr_(self):
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def random(self):
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        Gstr = str(gap(G))
        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        gap.eval("c:=Random( C )")
        ans = [gap_to_sage(gap.eval("c["+str(i)+"]"),F) for i in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(ans)

    def gen_mat(self):
        return self.__gen_mat

    def gens(self):
        return self.__gens

    def basis(self):
        return self.__gens

    def list(self):
        """
        Return list of all elements of this linear code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True

        """
        n = self.length()
        k = self.dimension()
        F = self.base_ring()
        Cs,p = self.standard_form()
        Gs = Cs.gen_mat()
        V = VectorSpace(F,k)
        MS = MatrixSpace(F,n,n)
        ans = []
        perm_mat = MS(p.matrix().rows())**(-1)
        for v in V:
            ans.append((v*Gs)*perm_mat)
        return ans

    def __iter__(self):
        """
        Return an iterator over the elements of this linear code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: [list(c) for c in C if hamming_weight(c) < 4]
            [[0, 0, 0, 0, 0, 0, 0],
             [1, 0, 0, 0, 0, 1, 1],
             [0, 1, 0, 0, 1, 0, 1],
             [0, 0, 1, 0, 1, 1, 0],
             [1, 1, 1, 0, 0, 0, 0],
             [1, 0, 0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0, 1, 0],
             [0, 0, 1, 1, 0, 0, 1]]
        """
        n = self.length()
        k = self.dimension()
        F = self.base_ring()
        Cs,p = self.standard_form()
        Gs = Cs.gen_mat()
        V = VectorSpace(F,k)
        MS = MatrixSpace(F,n,n)
        perm_mat = MS(p.matrix().rows())**(-1)
        for v in V:
            yield (v*Gs)*perm_mat

    def ambient_space(self):
        return VectorSpace(self.__base_ring,self.__length)

    def __contains__(self,v):
        A = self.ambient_space()
        C = A.subspace(self.gens())
        return C.__contains__(v)

    def characteristic(self):
        return (self.base_ring()).characteristic()

    def minimum_distance_lower_bound(self):
        r"""
        Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
        Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

        Obviously requires an internet connection
        """
        q = (self.base_ring()).order()
        n = self.length()
        k = self.dimension()
        bounds = linear_code_bound(q,n,k)
        return bounds[0]

    def minimum_distance_upper_bound(self):
        r"""
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
        r"""
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
        r"""
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
        r"""
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
        if n==k:
            return TrivialCode(F,n)
        Gstr = str(gap(G))
	C = gap.GeneratorMatCode(Gstr, 'GF(%s)'%q)
        #H = C.CheckMat()
        #A = H._matrix_(GF(q))
        #return LinearCode(A)       ## This does not work when k = n-1 for a mysterious reason.
        ##  less pythonic way 1:
        gap.eval("C:=DualCode(GeneratorMatCode("+Gstr+",GF("+str(q)+")))")
        gap.eval("G:=GeneratorMat(C)")
        G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,n-k+1)]
        MS = MatrixSpace(F,n-k,n)
        #print G, MS(G)
        return LinearCode(MS(G))
        ##  less pythonic way 2:
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

    def __eq__(self, right):
        """
        Checks if self == right.

        EXAMPLES:

        """
        slength = self.length()
        rlength = right.length()
        sdim = self.dimension()
        rdim = right.dimension()
        sF = self.base_ring()
        rF = right.base_ring()
        if slength != rlength:
            return 0
        if sdim != rdim:
            return 0
        if sF != rF:
            return 0
        V = VectorSpace(sF,sdim)
        sbasis = self.gens()
        rbasis = right.gens()
        scheck = self.check_mat()
        rcheck = right.check_mat()
        for c in sbasis:
            if rcheck*c != V(0):
                return 0
        for c in rbasis:
            if scheck*c != V(0):
                return 0
        return 1

    def is_permutation_automorphism(self,g):
        """
        Returns 1 if g is an element of S_n (n = length of self) and if
        g is an automorphism of self.

        EXAMPLES:
            sage: C = HammingCode(3,GF(3))
            sage: g = SymmetricGroup(13).random()
            sage: C.is_permutation_automorphism(g)
            0
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: S8 = SymmetricGroup(8)
            sage: g = S8("(2,3)")
            sage: C.is_permutation_automorphism(g)
            1
            sage: g = S8("(1,2,3,4)")
            sage: C.is_permutation_automorphism(g)
            0

        """
        basis = self.gen_mat().rows()
        H = self.check_mat()
        V = H.column_space()
        aut = 1
        for c in basis:
            if H*(g.matrix()*c) != V(0):
                aut = 0
        return aut

    def permutation_automorphism_group(self,mode=None):
        """
        If C is an [n,k,d] code over F this function computes the
        subgroup $Aut(C) \subset S_n$ of all permutation automorphisms of C.
        If mode="verbose" then code-theoretic data is printed out at several
        stages of the computation.

        Combines an idea of mine with an improvement suggested by Cary Huffman.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: C
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            144

        A less easy example involves showing that the permutation automorphism
        group of the extended ternary Golay code is the Mathieu group $M_{11}$.

            sage: C = ExtendedTernaryGolayCode()
            sage: M11 = MathieuGroup(11)
            sage: G = C.permutation_automorphism_group()  ## this should take < 15 seconds
            sage: G.is_isomorphic(M11)
            True

        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())                                 ## G is always full rank
        gap.eval("Gp:=SymmetricGroup(%s)"%n)               ## initializing G in gap
        Sn = SymmetricGroup(n)
        wts = self.spectrum()                                            ## bottleneck 1
        Gstr = str(gap(G))
        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        gap.eval("eltsC:=Elements(C)")
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0]
        if mode=="verbose":
            print "\n Minimum distance: %s \n Weight distribution: \n %s"%(nonzerowts[1],wts)
        stop = 0                                          ## only stop if all gens are autos
        for i in range(1,len(nonzerowts)):
          if stop == 1:
              break
          wt = nonzerowts[i]
          if mode=="verbose":
              size = eval(gap.eval("Size(Gp)"))
              print "\n Using the %s codewords of weight %s \n Supergroup size: \n %s\n "%(wts[wt],wt,size)
          gap.eval("Cwt:=Filtered(eltsC,c->WeightCodeword(c)=%s)"%wt)   ## bottleneck 2 (repeated
          gap.eval("matCwt:=List(Cwt,c->VectorCodeword(c))")            ##        for each i until stop = 1)
          gap.eval("A:=MatrixAutomorphisms(matCwt); GG:=Intersection(Gp,A)")    ## bottleneck 3
          #print i," = i \n",gap.eval("matCwt")," = matCwt\n"
          #print gap.eval("A")," = A \n",gap.eval("GG")," = GG\n\n"
          if eval(gap.eval("Size(GG)"))==0:
              #print gap.eval("GG; Size(GG)")
              #print "GG=0 ", gap.eval("A; Gp")
              return PermutationGroup([()])
          gap.eval("autgp_gens:=GeneratorsOfGroup(GG); Gp:=GG")
          #gap.eval("autgp_gens:=GeneratorsOfGroup(G)")
          N = eval(gap.eval("Length(autgp_gens)"))
          gens = [Sn(gap.eval("autgp_gens[%s]"%i).replace("\n","")) for i in range(1,N+1)]
          stop = 1                                        ## get ready to stop
          for x in gens:                                  ## if one of these gens is not an auto then don't stop
              if not(self.is_permutation_automorphism(x)):
                  stop = 0
                  break
        G = PermutationGroup(gens)
        return G

    def permuted_code(self,p):
        """
        Returns the permuted code - the code C which is equivalenet to self via
        the column permutation p.

        EXAMPLES:
            sage: C = ExtendedQuadraticResidueCode(7,GF(2))
            sage: G = C.permutation_automorphism_group()
            sage: p = G("(1,6,3,5)(2,7,4,8)")
            sage: Cp = C.permuted_code(p)
            sage: C.gen_mat()
            [1 1 0 1 0 0 0 1]
            [0 1 1 0 1 0 0 1]
            [0 0 1 1 0 1 0 1]
            [0 0 0 1 1 0 1 1]
            sage: Cp.gen_mat()
            [0 1 0 0 0 1 1 1]
            [1 1 0 0 1 0 1 0]
            [0 1 1 0 1 0 0 1]
            [1 1 0 1 0 0 0 1]
            sage: Cs1,p1 = C.standard_form(mode="verbose"); p1
            1 . . . 1 1 . 1
            . 1 . . . 1 1 1
             . . 1 . 1 1 1 .
             . . . 1 1 . 1 1
            ()
            sage: Cs2,p2 = Cp.standard_form(mode="verbose"); p2
             1 . . . 1 1 . 1
             . 1 . . . 1 1 1
             . . 1 . 1 1 1 .
             . . . 1 1 . 1 1
            ()

        Therefore you can see that Cs1 and Cs2 are the same, so C and Cp are
        equivalent.

        """
        F = self.base_ring()
        G = self.gen_mat()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        Gp = G*MS(p.matrix().rows())
        return LinearCode(Gp)

    def standard_form(self,mode=None):
        """
        An [n,k]  linear code with generator matrix G is in standard form is the row-reduced
        echelon form of G is (I,A), where I denotes the kxk identity matrix and A is a kx(n-k) block.
        This method returns a pair (C,p) where C is a code permutation equivalent to self
        and p in $S_n$ (n = length of C) is the permutation sending self to C.
        When mode = "verbose" the new generator matrix in the standard form
        (I,A) is "Display"ed.

        EXAMPLES:
            sage: C = ExtendedQuadraticResidueCode(7,GF(2))
            sage: C.gen_mat()
            [1 1 0 1 0 0 0 1]
            [0 1 1 0 1 0 0 1]
            [0 0 1 1 0 1 0 1]
            [0 0 0 1 1 0 1 1]
            sage: Cs,p = C.standard_form()
            sage: Cs.gen_mat()
            [1 0 0 0 1 1 0 1]
            [0 1 0 0 0 1 1 1]
            [0 0 1 0 1 1 1 0]
            [0 0 0 1 1 0 1 1]
            sage: p
            ()
            sage: C.standard_form(mode="verbose")
            <BLANKLINE>
            1 . . . 1 1 . 1
            . 1 . . . 1 1 1
            . . 1 . 1 1 1 .
            . . . 1 1 . 1 1
            <BLANKLINE>
            (Linear code of length 8, dimension 4 over Finite Field of size 2, ())

        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())                                 ## G is always full rank
        Sn = SymmetricGroup(n)
        Gstr = str(gap(G))
        gap.eval( "G:="+Gstr )
        p = Sn(gap.eval("p:=PutStandardForm(G)"))
        if mode=="verbose":
            print "\n",gap.eval("Display(G)"),"\n"
        gap.eval("C:=GeneratorMatCode(G,GF("+str(q)+"))")
        gap.eval("Gp:=GeneratorMat( C )")
        Gp = [[gap_to_sage(gap.eval("Gp["+str(i)+"]["+str(j)+"]"),F)
              for j in range(1,n+1)] for i in range(1,k+1)]
        MS = MatrixSpace(F,k,n)
        return LinearCode(MS(Gp)),p

    def module_composition_factors(self,gp):
        """
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.permutation_automorphism_group()

        Now type "C.module_composition_factors(gp)" to get the record printed.

        """
        F = self.base_ring()
        q = F.order()
        gens = gp.gens()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        MS = MatrixSpace(F,n,n)
        mats = [] # initializing list of mats by which the gens act on self
        Fn = VectorSpace(F, n)
        W = Fn.subspace_with_basis(G.rows()) # this is self
        for g in gens:
            p = MS(g.matrix())
            m = [W.coordinate_vector(r*p) for r in G.rows()]
            mats.append(m)
        mats_str = str(gap([[list(r) for r in m] for m in mats]))
        gap.eval("M:=GModuleByMats("+mats_str+", GF("+str(q)+"))")
        print gap("MTX.CompositionFactors( M )")

    def galois_closure(self, F0):
        """
        If self is a linear code defined over F and F0 is a subfield
        with Galois group G = Gal(F/F0) then this returns the G-module
        $C^-$ containing C.

        EXAMPLES:
            sage: C = HammingCode(3,GF(4))
            sage: C
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            sage: Cc = C.galois_closure(GF(2))
            sage: Cc
            Linear code of length 21, dimension 20 over Finite Field in a of size 2^2
            sage: c = C.random()
            sage: c  ## random output
            (1, 0, 1, 1, 1, a, 0, a, a + 1, a + 1, a + 1, a + 1, a, 1, a + 1, a + 1, 0, a + 1, 0, 1, 1)
            sage: V = VectorSpace(GF(4),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False
            sage: c2 in Cc
            True

        """
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        q0 = F0.order()
        a = log(q,q0)  ## test if F/F0 is a field extension
        if not(a.integer_part() == a):
            raise ValueError,"Base field must be an extension of given field %s"%F0
        n = len(G.columns())
        k = len(G.rows())
        G0 = [[x**q0 for x in g.list()] for g in G.rows()]
        G1 = [[x for x in g.list()] for g in G.rows()]
        G2 = G0+G1
        MS = MatrixSpace(F,2*k,n)
        G3 = MS(G2)
        r = G3.rank()
        MS = MatrixSpace(F,r,n)
        Grref = G3.echelon_form()
        G = MS([Grref[i] for i in range(r)])
        return LinearCode(G)

    def is_galois_closed(self):
        """
        Checks if C is equal to its Galois closure.

        """
        return self == self.galois_closure()

######### defining the Codeword class by copying the FreeModuleElement class:
Codeword = fme.FreeModuleElement
Codeword.support = fme.FreeModuleElement.nonzero_positions
is_Codeword = fme.is_FreeModuleElement


##################### wrapped GUAVA functions ############################

def HammingCode(r,F):
    """
    The $r^{th}$ Hamming code over $F=GF(q)$ is an $[n,k,d]$ code
    with length $n=(q^r-1)/(q-1)$, dimension $k=(q^r-1)/(q-1) - r$ and
    minimum distance $d=3$.
    The parity check matrix of a Hamming code has rows consisting of
    all nonzero vectors of length r in its columns, modulo a scalar factor
    so no parallel columns arise. A Hamming code is a single error-correcting
    code.

    INPUT:
        r -- an integer > 2
        F -- a finite field.

    OUTPUT:
        Returns the r-th q-ary Hamming code.

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
    r"""
    A quadratic residue code (or QR code) is a cyclic code whose
    generator polynomial is the product of the polynomials $x-\alpha^i$
    ($\alpha$ is a primitive $n^{th}$ root of unity; $i$ ranges over
    the set of quadratic residues modulo $n$).

    INPUT:
        n -- an odd prime
        F -- a finite prime field F whose order must be a quadratic
             residue modulo n.

    OUTPUT:
        Returns a quadratic residue code.

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

def ExtendedQuadraticResidueCode(n,F):
    """
    The extended quadratic residue code (or XQR code) is obtained from
    a QR code by adding a check bit to the last coordinate. (These codes
    have very remarkable properties such as large automorphism groups and
    duality properties - see [HP], \S 6.6.3-6.6.4.)

    INPUT:
        n -- an odd prime
        F -- a finite prime field F whose order must be a quadratic
             residue modulo n.

    OUTPUT:
        Returns an extended quadratic residue code.

    EXAMPLES:
        sage: C = ExtendedQuadraticResidueCode(7,GF(2))
        sage: C
        Linear code of length 8, dimension 4 over Finite Field of size 2
        sage: C = ExtendedQuadraticResidueCode(17,GF(2))
        sage: C
        Linear code of length 18, dimension 9 over Finite Field of size 2

    AUTHOR: David Joyner (07-2006)
    """
    q = F.order()
    gap.eval("C:=QRCode("+str(n)+", GF("+str(q)+"))")
    gap.eval("XC:=ExtendedCode(C)")
    gap.eval("G:=GeneratorMat(XC)")
    k = eval(gap.eval("Length(G)"))
    n = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(G))

def QuasiQuadraticResidueCode(p):
    """
    A (binary) quasi-quadratic residue code (or QQR code), as defined by
    Proposition 2.2 in [BM], has a generator matrix in the block form $G=(Q,N)$.
    Here $Q$ is a p x p circulant matrix whose top row
    is $(0,x_1,...,x_{p-1})$, where $x_i=1$ if and only if $i$
    is a quadratic residue $\mod p$, and $N$ is a p x p circulant matrix whose top row
    is $(0,y_1,...,y_{p-1})$, where $x_i+y_i=1$ for all i.

    INPUT:
        p -- a prime >2.

    OUTPUT:
        Returns a QQR code of length 2p.

    EXAMPLES:
        sage: C = QuasiQuadraticResidueCode(11)
        sage: C
        Linear code of length 22, dimension 11 over Finite Field of size 2

    AUTHOR: David Joyner (11-2005)

    REFERENCES:
        [BM] Bazzi-Mittel Some constructions of codes from group actions, (preprint March 2003,
             available on Mittel's MIT website).

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
    The binary 'Reed-Muller code' with dimension k and
    order r is a code with length $2^k$ and minimum distance $2^k-r$
    (see for example, section 1.10 in [HP]). By definition, the
    $r^{th}$ order binary Reed-Muller code of length $n=2^m$, for
    $0 \leq r \leq m$, is the set of all vectors $(f(p)\ |\ p in GF(2)^m)$,
    where $f$ is a multivariate polynomial of degree at most $r$ in $m$ variables.

    INPUT:
        r, k -- positive integers with $2^k>r$.

    OUTPUT:
        Returns the binary 'Reed-Muller code' with dimension k and order r.

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
    BinaryGolayCode() returns a binary Golay code. This is a perfect [23,12,7] code.
    It is also cyclic, and has generator polynomial $g(x)=1+x^2+x^4+x^5+x^6+x^{10}+x^{11}$.
    Extending it yields the extended Golay code (see ExtendedBinaryGolayCode).

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
    ExtendedBinaryGolayCode() returns the extended binary Golay code. This is a
    perfect [24,12,8] code. This code is self-dual.

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
    TernaryGolayCode returns a ternary Golay code. This is a perfect [11,6,5] code.
    It is also cyclic, and has generator polynomial $g(x)=2+x^2+2x^3+x^4+x^5$.

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
    ExtendedTernaryGolayCode returns a ternary Golay code. This is a self-dual perfect [12,6,6] code.

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
    The method used is to first construct a k x n matrix of the block form $(I,A)$,
    where $I$ is a k x k identity matrix and $A$ is a k x (n-k) matrix
    constructed using random elements of $F$. Then the columns are permuted
    using a randomly selected element of the symmetric group $S_n$.

    INPUT:
        Integers n,k, with n>k>1.

    OUTPUT:
        Returns a "random" linear code with length n, dimension k over field F.

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


def ToricCode(P,F):
    """
    Let $P$ denote a list of lattice points in $\Z^d$ and let $T$ denote the
    set of all points in $(F^x )^d$ (ordered in some fixed way). Put $n=|T|$
    and let $k$ denote the dimension of the vector space of functions
    $V = Span \{x^e \ |\ e \\in P\}$. The associated {\it toric code} $C$ is the
    evaluation code which is the image of the evaluation map
    $$
                        eval_T : V \rightarrow F^n,
    $$
    where $x^e$ is the multi-index notation ($x=(x_1,...,x_d)$, $e=(e_1,...,e_d)$, and
    $x^e = x_1^{e_1}...x_d^{e_d}$), where $eval_T (f(x)) = (f(t_1),...,f(t_n))$, and
    where $T=\{t_1,...,t_n\}$. This function returns the toric codes discussed in [J].

    INPUT:
        P -- all the integer lattice points in a polytope defining the toric variety.
        F -- a finite field.

    OUTPUT:
        Returns toric code with length n = , dimension k over field F.

    EXAMPLES:
        sage: C = ToricCode([[0,0],[1,0],[2,0],[0,1],[1,1]],GF(7))
        sage: C
        Linear code of length 36, dimension 5 over Finite Field of size 7
        sage: C.minimum_distance()
        24
        sage: C.minimum_distance_upper_bound()  # optional (requires internet)
        28
        sage: C = ToricCode([[-2,-2],[-1,-2],[-1,-1],[-1,0],[0,-1],[0,0],[0,1],[1,-1],[1,0]],GF(5))
        sage: C
        Linear code of length 16, dimension 9 over Finite Field of size 5
        sage: C.minimum_distance()
        6
        sage: C.minimum_distance_upper_bound()   # optional -- uses internet
        6
        sage: C = ToricCode([ [0,0],[1,1],[1,2],[1,3],[1,4],[2,1],[2,2],[2,3],[3,1],[3,2],[4,1]],GF(8))
        sage: C
        Linear code of length 49, dimension 11 over Finite Field in a of size 2^3
        sage: C.minimum_distance()  ## long time -- very time consuming
        28
        sage: print linear_code_bound(8,49,11)[0]    # optional -- uses internet
        28


    AUTHOR: David Joyner (07-2006)

    REFERENCES:
        [J] D. Joyner, {\it Toric codes over finite fields}, Applicable Algebra in Engineering,
            Communication and Computing, 15, (2004), p. 63--79
    """
    from sage.combinat.combinat import tuples
    mset = [x for x in F if x!=0]
    d = len(P[0])
    pts = tuples(mset,d)
    n = len(pts) ## (q-1)^d
    k = len(P)
    e = P[0]
    B = []
    for e in P:
       tmpvar = [prod([t[i]**e[i] for i in range(d)]) for t in pts]
       B.append(tmpvar)
    ## now B0 *should* be a full rank matrix
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(B))

def TrivialCode(F,n):
    MS = MatrixSpace(F,1,n)
    return LinearCode(MS(0))

def CyclicCode(g,n,F):
    """
    Here $g$ is a polynomial in $F[x]$ which divides $x^n-1$. The
    cyclic code C associated to (g,n,F) is isomorphic as a vector space)
    to the principal ideal $(g)$ in the ring $R = F[x]/(x^n-1)$ generated by g.

    If $g$ does not divide $x^n-1$, instead of returning an error, a
    code generated by $gcd(x^n-1,g)$ is created.

    EXAMPLES:
        sage: x = PolynomialRing(GF(2),"x").gen()
        sage: g = x^3+x+1
        sage: C = CyclicCode(g,7,GF(2))
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
    """
    q = F.order()
    R = gap.PolynomialRing(F)
    I = R.IndeterminatesOfPolynomialRing()
    _ = gap.eval("x := %s[1];"%I.name())
    gap_g = str(gap(str(g)))
    gap.eval("x_1:=X(GF("+str(q)+"))")
    gap.eval("C:=GeneratorPolCode("+gap_g+","+str(n)+", GF("+str(q)+"))")
    gap.eval("G:=GeneratorMat(C)")
    k = eval(gap.eval("Length(G)"))
    n1 = eval(gap.eval("Length(G[1])"))
    G = [[gap_to_sage(gap.eval("G["+str(i)+"]["+str(j)+"]"),F) for j in range(1,n1+1)] for i in range(1,k+1)]
    MS = MatrixSpace(F,k,n1)
    return LinearCode(MS(G))
