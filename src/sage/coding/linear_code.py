r"""
Linear Codes

VERSION: 0.7

Let $ F$ be a finite field (we denote the finite field with $q$ elements
$GF(q)$ by $\FF_q$). A subspace of $ F^n$ (with the standard basis)
is called a {\it linear code} of length $ n$. If its
dimension is denoted $k$ then we typically store a basis
of $C$ as a $k\times  n$ matrix (the rows are the basis vectors)
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
If an element $p\in S_n$ sends a code $C$ of length $n$ to itself
(in other words, every codeword of $C$ is sent to some other codeword
of $C$) then $p$ is called a {\it permutation automorphism} of $C$.
The (permutation) automorphism group is denoted $Aut(C)$.

This file contains
\begin{enumerate}
\item
LinearCode, Codeword class definitions; LinearCode_from_vectorspace
conversion function
\item
The spectrum (weight distribution), minimum distance
programs (calling Steve Linton's C programs), and zeta_function for
the Duursma zeta function.
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
gen_mat, list, check_mat, decode, dual_code, extended_code, genus methods
for LinearCode,
\item
permutation methods:
is_permutation_automorphism, permutation_automorphism_group,
permuted_code, standard_form, module_composition_factors.
\item
design-theoretic methods:
assmus_mattson_designs (implementing Assmus-Mattson Theorem).
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
\item $P^1$ Goppa codes and group actions on $P^1$ RR space codes.
\end{enumerate}

REFERENCES:
   [HP] W. C. Huffman and V. Pless, {\bf Fundamentals of error-correcting codes},
        Cambridge Univ. Press, 2003.

   [Gu] GUAVA manual, http://www.gap-system.org/Packages/guava.html

AUTHOR:
    -- David Joyner (2005-11-22, 2006-12-03): initial version
    -- William Stein (2006-01-23) -- Inclusion in SAGE
    -- David Joyner (2006-01-30, 2006-04): small fixes
    -- DJ (2006-07): added documentation, group-theoretical methods,
ExtendedQuadraticResidueCode, ToricCode
    -- DJ (2006-08): hopeful latex fixes to documention, added list
and __iter__ methods to LinearCode and examples, added hamming_weight
function, fixed random method to return a vector, TrivialCode,
fixed subtle bug in dual_code, added galois_closure method,
fixed mysterious bug in permutation_automorphism_group (GAP
was over-using "G" somehow?)
    -- DJ (2006-08): hopeful latex fixes to documention,
added CyclicCode, best_known_linear_code, bounds_minimum_distance,
assmus_mattson_designs (implementing Assmus-Mattson Theorem).
    -- DJ (2006-09): modified decode syntax, fixed bug in is_galois_closed,
added LinearCode_from_vectorspace, extended_code, zeta_function
    -- Nick Alexander (2006-12-10): factor GUAVA code to guava.py

TESTS:
   sage: MS = MatrixSpace(GF(2),4,7)
   sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
   sage: C  = LinearCode(G)
   sage: C == loads(dumps(C))
   True
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
from sage.databases.lincodes import linear_code_bound
from sage.interfaces.all import gap
from sage.misc.preparser import *
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.finite_field import *
from sage.groups.perm_gps.permgroup import *
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod, add
from sage.misc.functional import log
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens

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

def LinearCode_from_vectorspace(self):
    """
    Converts a VectorSpace over GF(q) into a LinearCode

    EXAMPLES:
        sage: V = VectorSpace(GF(2),7)
        sage: W = V.subspace([[1,0,0,1,1,1,1], [1,1,0,0,0,1,1]])
        sage: C = LinearCode_from_vectorspace(W)
        sage: C
        Linear code of length 7, dimension 2 over Finite Field of size 2
        sage: C.gen_mat()
        [1 0 0 1 1 1 1]
        [0 1 0 1 1 0 0]
    """
    F = self.base_ring()
    B = self.basis()
    n = len(B[0].list())
    k = len(B)
    MS = MatrixSpace(F,k,n)
    G = MS([B[i].list() for i in range(k)])
    return LinearCode(G)


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
        sage: MS = MatrixSpace(GF(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 5

    AUTHOR: David Joyner (11-2005)
    """
    def __init__(self, gen_mat):
        base_ring = gen_mat[0][0].parent()
        ParentWithGens.__init__(self, base_ring)
        self.__gens = gen_mat.rows()
        self.__gen_mat = gen_mat
        self.__length = len(gen_mat[0])
        self.__dim = gen_mat.rank()

    def length(self):
        return self.__length

    def dimension(self):
        return self.__dim

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
        return VectorSpace(self.base_ring(),self.__length)

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
            sage: C=RandomLinearCode(10,5,GF(4,'a'))
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
        gapG = gap(G)
        Gstr = "%s*Z(%s)^0"%(gapG, q)
        return min_wt_vec(Gstr,F)[2]

    def genus(self):
        """
        Returns the "Duursma genus" of the code, gamma_C = n+1-k-d.
        """
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-k-d
        return gammaC

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
        if not isinstance(right, LinearCode):
            return cmp(type(self), type(right))
        return cmp(self.__gen_mat, right.__gen_mat)

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
            sage: F = GF(2); a = F.gen()
            sage: v = [a,a,F(0),a,a,F(0),a]
            sage: C.decode(v)
            (1, 1, 0, 1, 0, 0, 1)

        Does not work for very long codes since the syndrome table grows too large.
        """
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        Gstr = str(gap(G))
        vstr = str(gap(right))
        if vstr[:3] == '[ [':
            vstr = vstr[1:-1]     # added by William Stein so const.tex works 2006-10-01
        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        gap.eval("c:=VectorCodeword(Decodeword( C, Codeword( "+vstr+" )))") # v->vstr, 8-27-2006
        ans = [gap_to_sage(gap.eval("c["+str(i)+"]"),F) for i in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(ans)

    def dual_code(self):
        """
        Wraps GUAVA's DualCode.

        OUTPUT:
            The dual code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.dual_code()
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C = HammingCode(3,GF(4,'a'))
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

    def extended_code(self):
        r"""
        If self is a linear code of length n defined over F then this
        returns the code of length n+1 where the last digit $c_n$
        satisfies the check condition $c_0+...+c_n=0$. If self is an
        $[n,k,d]$ binary code then the extended code $C^{\vee}$ is an
        $[n+1,k,d^{\vee}]$ code, where $d^=d$ (if d is even) and $d^{\vee}=d+1$ (if
        $d$ is odd).

        EXAMPLES:
            sage: C = HammingCode(3,GF(4,'a'))
            sage: C
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            sage: Cx = C.extended_code()
            sage: Cx
            Linear code of length 22, dimension 18 over Finite Field in a of size 2^2

        """
        from sage.rings.finite_field import gap_to_sage
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        n = len(G.columns())
        k = len(G.rows())
        Gstr = str(gap(G))
        gap.eval( "G:="+Gstr )
        gap.eval("C:=GeneratorMatCode(G,GF("+str(q)+"))")
        gap.eval("Gx:=GeneratorMat( ExtendedCode(C) )")
        Gx = [[gap_to_sage(gap.eval("Gx["+str(i)+"]["+str(j)+"]"),F)
              for j in range(1,n+2)] for i in range(1,k+1)]
        MS = MatrixSpace(F,k,n+1)
        return LinearCode(MS(Gx))

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
            return False
        if sdim != rdim:
            return False
        if sF != rF:
            return False
        V = VectorSpace(sF,sdim)
        sbasis = self.gens()
        rbasis = right.gens()
        scheck = self.check_mat()
        rcheck = right.check_mat()
        for c in sbasis:
            if rcheck*c:
                return False
        for c in rbasis:
            if scheck*c:
                return False
        return True

    def is_permutation_automorphism(self,g):
        """
        Returns 1 if g is an element of S_n (n = length of self) and if
        g is an automorphism of self.

        EXAMPLES:
            sage: C = HammingCode(3,GF(3))
            sage: g = SymmetricGroup(13).random_element()
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
        HGm = H*g.matrix()
        # raise TypeError, (type(H), type(V), type(basis[0]), type(Gmc))
        for c in basis:
            if HGm*c != V(0):
                return 0
        return 1

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
            sage: G = C.permutation_automorphism_group()   # long time
            sage: G.order()                                # long time
            144

        A less easy example involves showing that the permutation automorphism
        group of the extended ternary Golay code is the Mathieu group $M_{11}$.

            sage: C = ExtendedTernaryGolayCode()
            sage: M11 = MathieuGroup(11)
            sage: G = C.permutation_automorphism_group()  ## this should take < 15 seconds                                    # long time
            sage: G.is_isomorphic(M11)        # long time
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
            sage: C = HammingCode(3,GF(4,'a'))
            sage: C
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            sage: Cc = C.galois_closure(GF(2))
            sage: Cc
            Linear code of length 21, dimension 20 over Finite Field in a of size 2^2
            sage: c = C.random()
            sage: c  ## random output
            (1, 0, 1, 1, 1, a, 0, a, a + 1, a + 1, a + 1, a + 1, a, 1, a + 1, a + 1, 0, a + 1, 0, 1, 1)
            sage: V = VectorSpace(GF(4,'a'),21)
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
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def assmus_mattson_designs(self,t,mode=None):
        r"""
        Assmus and Mattson Theorem (section 8.4, page 303 of [HP]): Let A0, A1, ..., An
        be the weights of the codewords in a binary linear [n , k, d] code C, and let
        A0*, A1*, ..., An* be the weights of the codewords in its dual [n, n-k, d*] code C*.
        Fix a t, 0<t<d, and let $s = |{ i | Ai* \not= 0, 0<i\leq n-t}|$. Assume $s\leq d-t$.
        (1) If $Ai\not= 0$ and $d\leq i\leq n then Ci = { c in C | wt(c) = i}$
            holds a simple t-design.
        (2) If $Ai*\not= 0$ and $d*\leq i\leq n-t then Ci* = { c in C* | wt(c) = i}$
            holds a simple t-design.

        A {\bf block design} is a pair (X,B), where X is a non-empty finite set of v>0 elements called
        {\bf points}, and B is a non-empty finite multiset of size b whose elements are called {\bf blocks}, such
        that each block is a non-empty finite multiset of k points. A design without repeated blocks is called
        a {\bf simple} block design. If every subset of points of size t is contained in exactly lambda blocks
        the the block design is called a {\bf t-(v,k,lambda) design} (or simply a t-design when the parameters
        are not specfied). When lambda=1 then the block design is called a {\bf S(t,k,v) Steiner system}.

        In the Assmus and Mattson Theorem (1), X is the set {1,2,...,n} of coordinate locations
        and B = {supp(c) | c in Ci} is the set of supports of the codewords of C of weight i. Therefore,
        the parameters of the t-design for Ci are
             t =       given
             v =       n
             k =       i   (k not to be confused with dim(C))
             b =       Ai
             lambda =  b*binomial(k,t)/binomial(v,t)  (by Theorem 8.1.6, p 294, in [HP])

        Setting the mode="verbose" option prints out the values of the parameters.

        The first example below means that the binary [24,12,8]-code C has the property
        that the (support of the) codewords of weight 8 (resp, 12, 16) form a 5-design.
        Similarly for its dual code C* (of course C=C* in this case, so this info is extraneous).
        The test fails to produce 6-designs (ie, the hypotheses of the theorem fail
        to hold, not that the 6-designs definitely don't exist). The command
        assmus_mattson_designs(C,5,mode="verbose") returns the same value but prints out more
        detailed information.
        The second example below illustrates the blocks of the 5-(24, 8, 1) design
        (ie, the S(5,8,24) Steiner system).

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()             #  example 1
            sage: C.assmus_mattson_designs(5)
            ['weights from C: ',
            [8, 12, 16, 24],
            'designs from C: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)], [5, (24, 24, 1)]],
            'weights from C*: ',
            [8, 12, 16],
            'designs from C*: ',
            [[5, (24, 8, 1)], [5, (24, 12, 48)], [5, (24, 16, 78)]]]
            sage: C.assmus_mattson_designs(6)
            0
            sage: X = range(24)                           #  example 2
            sage: blocks = [c.support() for c in C if hamming_weight(c)==8]; len(blocks)  ## long time computation
            759

        REFERENCE:
            [HP] W. C. Huffman and V. Pless, {\bf Fundamentals of ECC}, Cambridge Univ. Press, 2003.
        """
        from sage.rings.arith import binomial
        C = self
        ans = []
        F = C.base_ring()
        q = F.order()
        G = C.gen_mat()
        n = len(G.columns())
        Cp = C.dual_code()
        k = len(G.rows())  ## G is always full rank
        wts = C.spectrum()
        d = min([i for i in range(1,len(wts)) if wts[i]!=0])
        if t>=d:
            return 0
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0 and i<=n and i>=d]
        #print d,t,len(nonzerowts)
        if mode=="verbose":
            #print "\n"
            for w in nonzerowts:
                print "The weight ",w," codewords of C form a t-(v,k,lambda) design, where"
                print "      t = ",t," , v = ",n," , k = ",w," , lambda = ",wts[w]*binomial(w,t)/binomial(n,t)#,"\n"
                print "      There are ",wts[w]," blocks of this design."
        wtsp = Cp.spectrum()
        dp = min([i for i in range(1,len(wtsp)) if wtsp[i]!=0])
        nonzerowtsp = [i for i in range(len(wtsp)) if wtsp[i]!=0 and i<=n-t and i>=dp]
        s = len([i for i in range(1,n) if wtsp[i]!=0 and i<=n-t and i>0])
        if mode=="verbose":
            #print "\n"
            for w in nonzerowtsp:
                print "The weight ",w," codewords of C* form a t-(v,k,lambda) design, where"
                print "      t = ",t," , v = ",n," , k = ",w," , lambda = ",wtsp[w]*binomial(w,t)/binomial(n,t)#,"\n"
                print "      There are ",wts[w]," blocks of this design."
        if s<=d-t:
            des = [[t,(n,w,wts[w]*binomial(w,t)/binomial(n,t))] for w in nonzerowts]
            ans = ans + ["weights from C: ",nonzerowts,"designs from C: ",des]
            desp = [[t,(n,w,wtsp[w]*binomial(w,t)/binomial(n,t))] for w in nonzerowtsp]
            ans = ans + ["weights from C*: ",nonzerowtsp,"designs from C*: ",desp]
            return ans
        return 0

    def zeta_function(self,mode=None):
        """
        Returns the Duursma zeta function of the code.

        mode = "dual" computes both the zeta function of $C$ and that of $C*$
        mode = "normalized" computes the normalized zeta function of $C$, $zeta_C(T)*T(1-genus(C))$
               (NOTE: if xi(T,C) denotes the normalized zeta function then
                $xi(T,C*) = xi(1/(qT),C)$ is equivalent to
                $zeta_{C*}(T) = zeta_C(1/(qT))*q^(gamma-1)T^(gamma+gamma*-2)$, where
                $gamma = gamma(C)$ and $gamma* = gamma(C*)$.)

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.zeta_function()
            (1/5 + 2/5*T + 2/5*T^2)/(1 - 3*T + 2*T^2)
            sage: C = ExtendedTernaryGolayCode()
            sage: C.zeta_function()
            (1/7 + 3/7*T + 3/7*T^2)/(1 - 4*T + 3*T^2)

        Both these examples occur in Duursma's paper below.

        REFERENCES:
             I. Duursma, "Weight distributions of geometric Goppa codes," Trans. A.M.S., 351 (1999)3609-3639.

        NOTE: This is somewhat experimental code. It sometimes returns "fail" for a reason
        I don't fully understand. However, when it does return a polynomial, the answer is
        (as far as I know) correct.
        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing
        from sage.rings.fraction_field import FractionField
        from sage.rings.power_series_ring import PowerSeriesRing
        R = PolynomialRing(QQ,3,"xyT")
        #F = FractionField(R)
        x,y,T = R.gens()
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-d #  C.genus()
        if mode=="dual":
            Cp = self.dual_code()
            #dp = Cp.minimum_distance()
            #kp = Cp.dimension()
            gammaCp = Cp.genus() # = n+1-kp-dp
        q = self.characteristic()
        W = self.weight_distribution()
        A = add([W[i]*y**i*x**(n-i) for i in range(n+1)])
        f = (x*T+y*(1-T))**n*add([T**i for i in range(n+3)])*add([(q*T)**i for i in range(n+3)])
        Mn = [f.coefficient(T**(n-i)) for i in range(d,n)]
        coeffs = [[(Mn[j]+x**(n)).monomial_coefficient(x**i*y**(n-i)) for i in range(n+1)] for j in range(n-d)]
        MS = MatrixSpace(QQ,n-d,n+1)
        M = MS(coeffs)
        F = PowerSeriesRing(PolynomialRing(QQ,2,"xy"),"T")
        for m in range(2,n-d+1):
            #x,y,T = F.gens()
            V = VectorSpace(QQ,m)
            v = V([W[i] for i in range(m)])
            Mmm = M.matrix_from_rows_and_columns(range(m),range(m))
            if Mmm.determinant()!=0:
                Pcoeffs = v*Mmm**(-1)
                P = add([Pcoeffs[i]*T**i for i in range(len(Pcoeffs))])
                #x,y,T = F.gens()
                Z = P*(1-T)**(-1)*(1-q*T)**(-1)
                lhs = P*f
                r = lhs.coefficient(T**(n-d))/((A-x**n)/(q-1))
                if mode=="verbose": print m, P, r
                #if r(x,y,T)==r(1,y,T) and mode=="dual":         ## check if r is a constant in x ...
                #    Z = F(Z*r(1,1,1)**(-1))
                #    x,y,T = F.gens()
                #    Zp = Z(x,y,1/(q*T))*q**(gammaC-1)*T**(gammaC+gammaCp-2)
                #    return Z,Zp
                if r(x,y,T)==r(1,y,T) and mode=="normalized":   ## check if r is a constant in x ...
                    Z = Z*r(1,1,1)**(-1)
                    #x,y,T = F.gens()
                    return Z*T**(1-gammaC)
                if r(x,y,T)==r(1,y,T):                         ## check if r is a constant in x ...
                    Z = Z*r(1,1,1)**(-1)
                    return Z
                if mode=="verbose":
                    print "det = ",Mmm.determinant(),"   m = ",m,"  Z = ",Z
            #if Mmm.determinant()==0:
            #    print "det = ",Mmm.determinant(),"   m = ",m
        return "fails"
