r"""
Linear Codes

VERSION: 0.10

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

If $ F=\FF_2$ then $C$ is called a {\it binary code}.
If $ F = \FF_q$ then $C$ is called a {\it $ q$-ary code}.
The elements of a code $C$ are called {\it codewords}.

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
programs (calling Steve Linton's C programs), characteristic function,
and several implementations of the Duursma zeta function.
\item
interface with best_known_linear_code (interface with A. Brouwer's online tables has
been disabled), bounds_minimum_distance which call tables
in GUAVA (updated May 2006) created by Cen Tjhai instead of the online
internet tables.
\item
ToricCode, TrivialCode (in a separate "guava.py" module, you will find the constructions
Hamming codes, "random" linear codes, Golay codes, binary Reed-Muller codes,
binary quadratic and extended quadratic residue codes, cyclic codes, all
wrapped form the corresponding GUAVA codes),
\item
gen_mat, list, check_mat, decode, dual_code, extended_code, genus, binomial_moment,
and divisor methods for LinearCode,
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
    -- DJ (2007-05): added methods punctured, shortened, divisor,
characteristic_polynomial, binomial_moment, support for LinearCode.
Completely rewritten zeta_function (old version is now zeta_function2)
and a new function, LinearCodeFromVectorSpace.
    -- DJ (2007-11): added zeta_polynomial, weight_enumerator,
                     chinen_polynomial;  improved best_known_code;
                     made some pythonic revisions;
                     added is_equivalent (for binary codes)
    -- dj (2008-01): fixed bug in decode reported by harald schilly,
                     (with M Hansen) added some doctests.
    -- dj (2008-02): translated standard_form, dual_code to Python.

TESTS:
   sage: MS = MatrixSpace(GF(2),4,7)
   sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
   sage: C  = LinearCode(G)
   sage: C == loads(dumps(C))
   True
"""

#*****************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import copy
import sage.modules.free_module as fm
import sage.modules.module as module
import sage.modules.free_module_element as fme
from sage.interfaces.all import gap
from sage.rings.finite_field import FiniteField as GF
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.arith import GCD, rising_factorial, binomial
from sage.groups.all import SymmetricGroup
from sage.misc.sage_eval import sage_eval
from sage.misc.misc import prod, add
from sage.misc.functional import log, is_even, is_odd
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import IntegerRing
from sage.combinat.set_partition import SetPartitions
from sage.rings.real_mpfr import RR      ## RealField

ZZ = IntegerRing()
VectorSpace = fm.VectorSpace

###################### coding theory functions ##############################

def hamming_weight(v):
    return len(v.nonzero_positions())

def wtdist(Gmat, F):
    r"""
    INPUT:
        Gmat -- a string representing a GAP generator matrix G of a
                linear code.
        F -- a (SAGE) finite field - the base field of the code.

    OUTPUT:
        Returns the spectrum of the associated code.

    EXAMPLES:
        sage: Gstr = 'Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]'
        sage: F = GF(2)
    	sage: sage.coding.linear_code.wtdist(Gstr, F)
    	[1, 0, 0, 7, 7, 0, 0, 1]

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    ALGORITHM: Uses C programs written by Steve Linton in the kernel
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
    ## for some reason, this commented code doesn't work:
    #dist0 = gap("DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    #v0 = dist0._matrix_(F)
    #print dist0,v0
    #d = G.DistancesDistributionMatFFEVecFFE(k, z)
    v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,n+2)] ## because GAP returns vectors in compressed form
    return v

def min_wt_vec(Gmat,F):
    r"""
    Uses C programs written by Steve Linton in the kernel of GAP, so
    is fairly fast.

    INPUT:
        Same as wtdist.

    OUTPUT:
        Returns a minimum weight vector v of the code generated by Gmat
         \#\# , the"message" vector m such that m*G = v, and the (minimum) distance, as a triple.

    EXAMPLES:
        sage: Gstr = "Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]"
    	sage: sage.coding.linear_code.min_wt_vec(Gstr,GF(2))
    	(0, 0, 1, 0, 1, 1, 0)

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    AUTHOR: David Joyner (11-2005)
    """
    from sage.interfaces.gap import gfq_gap_to_sage
    gap.eval("G:="+Gmat)
    k = int(gap.eval("Length(G)"))
    q = F.order()
    qstr = str(q)
    C = gap(Gmat).GeneratorMatCode(F)
    n = int(C.WordLength())
    cg = C.MinimumDistanceCodeword()
    c = [gfq_gap_to_sage(cg[j],F) for j in range(1,n+1)]
    V = VectorSpace(F,n)
    return V(c)
    ## this older code returns more info but may be slower:
    #zerovec = [0 for i in range(n)]
    #zerovecstr = "Z("+qstr+")*"+str(zerovec)
    #all = []
    #for i in range(1,k+1):
    #    P = gap.eval("P:=AClosestVectorCombinationsMatFFEVecFFECoords("+Gmat+", GF("+qstr+"),"+zerovecstr+","+str(i)+","+str(0)+"); d:=WeightVecFFE(P[1])")
    #    v = gap("[List(P[1], i->i)]")
    #    m = gap("[List(P[2], i->i)]")
    #    dist = gap.eval("d")
    #    #print v,m,dist
    #    #print [gap.eval("v["+str(i+1)+"]") for i in range(n)]
    #    all.append([v._matrix_(F),m._matrix_(F),int(dist)])
    #ans = all[0]
    #for x in all:
    #    if x[2]<ans[2] and x[2]>0:
    #        ans = x
    #return ans

## def minimum_distance_lower_bound(n,k,F):
##     r"""
##     Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
##     Tables of A. E. Brouwer,   Techn. Univ. Eindhoven,
##     via Steven Sivek's linear_code_bound.

##     EXAMPLES:
##     sage: sage.coding.linear_code.minimum_distance_upper_bound(7,4,GF(2))     # optional (net connection)
##         3

##     Obviously requires an internet connection.
##     """
##     q = F.order()
##     bounds = linear_code_bound(q,n,k)
##     return bounds[0]

## def minimum_distance_upper_bound(n,k,F):
##     r"""
##     Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
##     Tables of A. E. Brouwer,   Techn. Univ. Eindhoven
##     via Steven Sivek's linear_code_bound.

##     EXAMPLES:
##         sage: sage.coding.linear_code.minimum_distance_upper_bound(7,4,GF(2))  # optional (net connection)
##         3

##     Obviously requires an internet connection.
##     """
##     q = F.order()
##     bounds = linear_code_bound(q,n,k)
##     return bounds[1]

## def minimum_distance_why(n,k,F):
##     r"""
##     Connects to http://www.win.tue.nl/~aeb/voorlincod.html
##     Tables of A. E. Brouwer,   Techn. Univ. Eindhoven
##     via Steven Sivek's linear_code_bound.

##     EXAMPLES:
##         sage: sage.coding.linear_code.minimum_distance_why(7,4,GF(2))  # optional (net connection)
##         Lb(7,4) = 3 is found by truncation of:
##         Lb(8,4) = 4 is found by the (u|u+v) construction
##         applied to [4,3,2] and [4,1,4]-codes
##         Ub(7,4) = 3 follows by the Griesmer bound.

##     Obviously requires an internet connection.
##     """
##     q = F.order()
##     bounds = linear_code_bound(q,n,k)
##     print bounds[2]

def best_known_linear_code(n,k,F):
    r"""
    best_known_linear_code returns the best known (as of 11 May 2006)
    linear code of length n, dimension k over field F. The function
    uses the tables described in bounds_minimum_distance to construct
    this code.

    This does not require an internet connection.

    EXAMPLES:
        sage: best_known_linear_code(10,5,GF(2))    # long time
        Linear code of length 10, dimension 5 over Finite Field of size 2
        sage: gap.eval("C:=BestKnownLinearCode(10,5,GF(2))")     # long time
        'a linear [10,5,4]2..4 shortened code'

    This means that best possible binary linear code of length 10 and dimension 5
    is a code with minimum distance 4 and covering radius somewhere between 2 and 4.
    Use "minimum_distance_why(10,5,GF(2))" or "print bounds_minimum_distance(10,5,GF(2))"
    for further details.
    """
    q = F.order()
    C = gap("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))
    G = C.GeneratorMat()
    k = G.Length()
    n = G[1].Length()
    Gs = G._matrix_(F)
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(Gs))
    #return gap.eval("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))

def bounds_minimum_distance(n,k,F):
    r"""
    The function bounds_minimum_distance calculates a lower and upper
    bound for the minimum distance of an optimal linear code with word
    length n, dimension k over field F. The function returns a record
    with the two bounds and an explanation for each bound. The
    function Display can be used to show the explanations.

    The values for the lower and upper bound are obtained from a table
    constructed by Cen Tjhai for GUAVA, derived from the table of
    Brouwer. (See http://www.win.tue.nl/~aeb/voorlincod.html or use
    the SAGE function minimum_distance_why for the most recent data.)
    These tables contain lower and upper bounds for q=2 (n <= 257), 3
    (n <= 243), 4 (n <= 256). (Current as of 11 May 2006.)  For codes
    over other fields and for larger word lengths, trivial bounds are
    used.

    This does not require an internet connection. The format of the
    output is a little non-intuitive. Try print
    bounds_minimum_distance(10,5,GF(2)) for example.
    """
    q = F.order()
    gap.eval("data := BoundsMinimumDistance(%s,%s,GF(%s))"%(n,k,q))
    Ldata = gap.eval("Display(data)")
    return Ldata


########################### linear codes python class #######################

class LinearCode(module.Module):
    r"""
    A class for linear codes over a finite field or finite ring.  Each
    instance is a linear code determined by a generator matrix $G$
    (i.e., a k x n matrix of (full) rank $k$, $k\leq n$ over a finite
    field $F$.

    INPUT:
        G -- a generator matrix over $F$. (G can be defined over a
             finite ring but the matrices over that ring must have
             certain attributes, such as"rank".)

    OUTPUT:
        The linear code of length $n$ over $F$ having $G$ as a
        generator matrix.

    EXAMPLES:
        sage: MS = MatrixSpace(GF(2),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 2
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
        sage: MS = MatrixSpace(GF(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 5

    AUTHOR: David Joyner (11-2005)
    """
    #    sage: C.minimum_distance_upper_bound()   # optional (net connection)
    #    3
    #    sage: C.minimum_distance_why()     # optional (net connection)
    #    Ub(7,4) = 3 follows by the Griesmer bound.
    def __init__(self, gen_mat):
        base_ring = gen_mat[0,0].parent()
        ParentWithGens.__init__(self, base_ring)
        self.__gens = gen_mat.rows()
        self.__gen_mat = gen_mat
        self.__length = len(gen_mat.row(0))
        self.__dim = gen_mat.rank()

    def length(self):
        return self.__length

    def dimension(self):
        return self.__dim

    def _repr_(self):
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def random(self):
        """
        EXAMPLES:
            sage: C = HammingCode(3,GF(4,'a'))
            sage: Cc = C.galois_closure(GF(2))
            sage: c = C.random()
            sage: V = VectorSpace(GF(4,'a'),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False

        """
        from sage.interfaces.gap import gfq_gap_to_sage
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        Cg = gap(G).GeneratorMatCode(F)
        c = Cg.Random()
        ans = [gfq_gap_to_sage(c[i],F) for i in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(ans)

    def gen_mat(self):
        return self.__gen_mat

    def gens(self):
        return self.__gens

    def basis(self):
        return self.__gens

    def list(self):
        r"""
        Return list of all elements of this linear code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 0, 1, 1)
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
            [[0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 1, 0], [1, 1, 0, 0, 0, 0, 1], [0, 0, 1, 1, 0, 0, 1], [0, 1, 1, 0, 0, 1, 0], [0, 0, 0, 0, 1, 1, 1], [0, 1, 0, 1, 1, 0, 0], [1, 0, 1, 0, 1, 0, 0]]

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

##     def minimum_distance_lower_bound(self):
##         r"""
##         Connects to \verb+http://www.win.tue.nl/~aeb/voorlincod.html+
##         Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

##         Obviously requires an internet connection
##         """
##         q = (self.base_ring()).order()
##         n = self.length()
##         k = self.dimension()
##         bounds = linear_code_bound(q,n,k)
##         return bounds[0]

##     def minimum_distance_upper_bound(self):
##         r"""
##         Connects to http://www.win.tue.nl/~aeb/voorlincod.html
##         Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

##         Obviously requires an internet connection
##         """
##         q = (self.base_ring()).order()
##         n = self.length()
##         k = self.dimension()
##         bounds = linear_code_bound(q,n,k)
##         return bounds[1]

##     def minimum_distance_why(self):
##         r"""
##         Connects to http://www.win.tue.nl/~aeb/voorlincod.html
##         Tables of A. E. Brouwer,   Techn. Univ. Eindhoven

##         Obviously requires an internet connection.
##         """
##         q = (self.base_ring()).order()
##         n = self.length()
##         k = self.dimension()
##         bounds = linear_code_bound(q,n,k)
##         lines = bounds[2].split("\n")
##         for line in lines:
##             if len(line)>0:
##                 if line[0] == "U":
##                     print line

    def minimum_distance(self):
        r"""
        Uses a GAP kernel function (in C) written by Steve Linton.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3
        """
        #sage: C.minimum_distance_upper_bound()  # optional (net connection)
        #5
        #    sage: C.minimum_distance_why()          # optional (net connection)
        #    Ub(10,5) = 5 follows by the Griesmer bound.

        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        gapG = gap(G)
        Gstr = "%s*Z(%s)^0"%(gapG, q)
        return hamming_weight(min_wt_vec(Gstr,F))

    def genus(self):
        r"""
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
        r"""
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
            sage: v1 = [a,a,F(0),a,a,F(0),a]
            sage: C.decode(v1)
            (1, 0, 0, 1, 1, 0, 1)
            sage: v2 = matrix([[a,a,F(0),a,a,F(0),a]])
            sage: C.decode(v2)
            (1, 0, 0, 1, 1, 0, 1)
            sage: v3 = vector([a,a,F(0),a,a,F(0),a])
            sage: c = C.decode(v3); c
            (1, 0, 0, 1, 1, 0, 1)
            sage: c in C
            True
            sage: v4 = [[a,a,F(0),a,a,F(0),a]]
            sage: C.decode(v4)
            (1, 0, 0, 1, 1, 0, 1)
            sage: C = HammingCode(2,GF(5))
            sage: v = vector(GF(5),[1,0,0,2,1,0])
            sage: C.decode(v)
 	    (2, 0, 0, 2, 1, 0)

        Does not work for very long codes since the syndrome table grows too large.
        """
        from sage.interfaces.gap import gfq_gap_to_sage
        F = self.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        Gstr = str(gap(G))
        if not(type(right) == list):
            v = right.list()
        else:
            v = right
        vstr = str(gap(v))
        if vstr[:3] == '[ [':
            vstr = vstr[1:-1]     # added by William Stein so const.tex works 2006-10-01
        gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
        gap.eval("c:=VectorCodeword(Decodeword( C, Codeword( "+vstr+" )))") # v->vstr, 8-27-2006
        ans = [gfq_gap_to_sage(gap.eval("c["+str(i)+"]"),F) for i in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(ans)

    def dual_code(self):
        r"""
        This computes the dual code Cd of the code C,
        \[
        Cd = \{ v \in V\ |\ v\cdot c = 0,\ \forall c \in C \}.
        \]
        Does not call GAP.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.dual_code()
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C = HammingCode(3,GF(4,'a'))
            sage: C.dual_code()
            Linear code of length 21, dimension 3 over Finite Field in a of size 2^2
        """
        G = self.gen_mat()
        H = G.transpose().kernel()
        V = H.ambient_vector_space()
        Cd = LinearCodeFromVectorSpace(V.span(H))
        return Cd
        #another way:
        #Gsf, p = standard_form(G)
        #k = len(G.rows())
        #n = len(G.columns())
        #MS = G.parent()
        #sG = G.matrix_from_columns(range(k,n))
        #Inmk = MatrixSpace(F,n-k,n-k).identity_matrix()
        #H = Inmk.augment(sG.transpose())
        #return LinearCode(H)

    def extended_code(self):
        r"""
        If self is a linear code of length n defined over F then this
        returns the code of length n+1 where the last digit $c_n$
        satisfies the check condition $c_0+...+c_n=0$. If self is an
        $[n,k,d]$ binary code then the extended code $C^{\vee}$ is an
        $[n+1,k,d^{\vee}]$ code, where $d^=d$ (if d is even) and
        $d^{\vee}=d+1$ (if $d$ is odd).

        EXAMPLES:
            sage: C = HammingCode(3,GF(4,'a'))
            sage: C
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            sage: Cx = C.extended_code()
            sage: Cx
            Linear code of length 22, dimension 18 over Finite Field in a of size 2^2
        """
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        n = len(G.columns())
        k = len(G.rows())
        g = gap(G)
        gap.eval( "G:="+g.name())
        C = gap("GeneratorMatCode(G,GF("+str(q)+"))")
        Cx = C.ExtendedCode()
        Gx = Cx.GeneratorMat()
        Gxs = Gx._matrix_(F)           # this is the killer
        MS = MatrixSpace(F,k,n+1)
        return LinearCode(MS(Gxs))

    def check_mat(self):
        r"""
        Returns the check matrix of self.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: Cperp = C.dual_code()
            sage: C; Cperp
            Linear code of length 7, dimension 4 over Finite Field of size 2
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C.gen_mat()
            [1 0 0 1 0 1 0]
            [0 1 0 1 0 1 1]
            [0 0 1 1 0 0 1]
            [0 0 0 0 1 1 1]
            sage: C.check_mat()
            [1 0 0 1 1 0 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 1 0]
            sage: Cperp.check_mat()
            [1 0 0 1 0 1 0]
            [0 1 0 1 0 1 1]
            [0 0 1 1 0 0 1]
            [0 0 0 0 1 1 1]
            sage: Cperp.gen_mat()
            [1 0 0 1 1 0 1]
            [0 1 0 1 0 1 1]
            [0 0 1 1 1 1 0]

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
        r"""
        Returns $1$ if $g$ is an element of $S_n$ ($n$ = length of
        self) and if $g$ is an automorphism of self.

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
        r"""
        If $C$ is an $[n,k,d]$ code over $F$, this function computes
        the subgroup $Aut(C) \subset S_n$ of all permutation
        automorphisms of $C$.  If mode="verbose" then
        code-theoretic data is printed out at several stages
        of the computation.

        Combines an idea of mine with an improvement suggested by Cary
        Huffman.

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
            sage: G = C.permutation_automorphism_group()  ## this should take < 15 seconds  # long time
            sage: G.is_isomorphic(M11)        # long time
            True

        WARNING: - *Ugly code, which should be replaced by a call to Robert Miller's nice program.*
                 - Known to mysteriously crash in one example.
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
        r"""
        Returns the permuted code -- the code $C$ which is equivalent
        to self via the column permutation $p$.

        EXAMPLES:
            sage: C = ExtendedQuadraticResidueCode(7,GF(2))
            sage: G = C.permutation_automorphism_group()
            sage: p = G("(1,6,3,5)(2,7,4,8)")
            sage: Cp = C.permuted_code(p)
            sage: C == Cp
            True

        """
        F = self.base_ring()
        G = self.gen_mat()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        Gp = G*MS(p.matrix().rows())
        return LinearCode(Gp)

    def standard_form(self):
        r"""
        An $[n,k]$ linear code with generator matrix $G$ is in
        standard form is the row-reduced echelon form of $G$ is
        $(I,A)$, where $I$ denotes the $k \times k$ identity matrix
        and $A$ is a $k \times (n-k)$ block.  This method returns a
        pair $(C,p)$ where $C$ is a code permutation equivalent to
        self and $p$ in $S_n$ ($n$ = length of $C$) is the permutation
        sending self to $C$. This does not call GAP.

        Thanks to Frank Luebeck for (the GAP version of) this code.

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
            sage: MS = MatrixSpace(GF(3),3,7)
            sage: G = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
            sage: C = LinearCode(G)
            sage: G; C.standard_form()[0].gen_mat()
            [1 0 0 0 1 1 0]
            [0 1 0 1 0 1 0]
            [0 0 0 0 0 0 1]
            [1 0 0 0 1 1 0]
            [0 1 0 1 0 1 0]
            [0 0 1 0 0 0 0]
            sage: C.standard_form()[1]
            (3,7)

        """
        from sage.coding.code_constructions import permutation_action as perm_action
        mat = self.gen_mat()
        MS = mat.parent()
        A = []
        k = len(mat.rows())
        M = mat.echelon_form()
        d = len(mat.columns())
        G = SymmetricGroup(d)
        perm = G([()])
        for i in range(1,k+1):
            r = M.rows()[i-1]
            j = r.nonzero_positions()[0]
            if (j < d and i <> j+1):
                perm = perm *G([(i,j+1)])
        if perm <> G([()]):
            for i in range(k):
                r = M.rows()[i]
                A.append(perm_action(perm,r))
        if perm == G([()]):
            A = M
        return LinearCode(MS(A)), perm

    def module_composition_factors(self,gp):
        r"""
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.permutation_automorphism_group()    # long time

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
        r"""
        If self is a linear code defined over $F$ and $F_0$ is a subfield
        with Galois group $G = Gal(F/F_0)$ then this returns the $G$-module
        $C^-$ containing $C$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(4,'a'))
            sage: Cc = C.galois_closure(GF(2))
            sage: C; Cc
            Linear code of length 21, dimension 18 over Finite Field in a of size 2^2
            Linear code of length 21, dimension 20 over Finite Field in a of size 2^2
            sage: c = C.random()
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
        G = MS([Grref.row(i) for i in range(r)])
        return LinearCode(G)

    def is_galois_closed(self):
        r"""
        Checks if $C$ is equal to its Galois closure.
        """
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def assmus_mattson_designs(self,t,mode=None):
        r"""
        Assmus and Mattson Theorem (section 8.4, page 303 of [HP]):
        Let $A_0, A_1, ..., A_n$ be the weights of the codewords in a
        binary linear $[n , k, d]$ code $C$, and let $A_0^*, A_1^*,
        ..., A_n^*$ be the weights of the codewords in its dual $[n,
        n-k, d^*]$ code $C^*$.  Fix a $t$, $0<t<d$, and let
        $$
             s = |{ i | A_i^* \not= 0, 0<i\leq n-t}|.
        $$
        Assume $s\leq d-t$.

        (1) If $Ai\not= 0$ and $d\leq i\leq n then Ci = { c in C | wt(c) = i}$
            holds a simple t-design.

        (2) If $Ai*\not= 0$ and $d*\leq i\leq n-t then Ci* = { c in C*
            | wt(c) = i}$ holds a simple t-design.

        A {\bf block design} is a pair $(X,B)$, where $X$ is a
        non-empty finite set of $v>0$ elements called {\bf points},
        and $B$ is a non-empty finite multiset of size b whose
        elements are called {\bf blocks}, such that each block is a
        non-empty finite multiset of $k$ points. $A$ design without
        repeated blocks is called a {\bf simple} block design. If
        every subset of points of size $t$ is contained in exactly
        lambda blocks the the block design is called a {\bf
        $t-(v,k,lambda)$ design} (or simply a $t$-design when the
        parameters are not specfied). When $\lambda=1$ then the block
        design is called a {\bf $S(t,k,v)$ Steiner system}.

        In the Assmus and Mattson Theorem (1), $X$ is the set
        $\{1,2,...,n\}$ of coordinate locations and
        $B = \{supp(c) | c in C_i\}$ is the set of
        supports of the codewords of $C$ of weight $i$.
        Therefore, the parameters of the $t$-design for $C_i$ are
        \begin{verbatim}
             t =       given
             v =       n
             k =       i   (k not to be confused with dim(C))
             b =       Ai
             lambda = b*binomial(k,t)/binomial(v,t) (by Theorem 8.1.6,
                      p 294, in [HP])
        \end{verbatim}

        Setting the mode="verbose" option prints out the values of the parameters.

        The first example below means that the binary [24,12,8]-code C
        has the property that the (support of the) codewords of weight
        8 (resp, 12, 16) form a 5-design.  Similarly for its dual code
        C* (of course C=C* in this case, so this info is extraneous).
        The test fails to produce 6-designs (ie, the hypotheses of the
        theorem fail to hold, not that the 6-designs definitely don't
        exist). The command assmus_mattson_designs(C,5,mode="verbose")
        returns the same value but prints out more detailed information.

        The second example below illustrates the blocks of the 5-(24,
        8, 1) design (ie, the S(5,8,24) Steiner system).

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

    def weight_enumerator(self,names="xy"):
        """
        Returns the weight enumerator of the code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.weight_enumerator()
            x^7 + 7*x^4*y^3 + 7*x^3*y^4 + y^7
            sage: C.weight_enumerator(names="st")
            s^7 + 7*s^4*t^3 + 7*s^3*t^4 + t^7
        """
        spec = self.spectrum()
        n = self.length()
        from sage.rings.polynomial.polynomial_ring import PolynomialRing
        R = PolynomialRing(QQ,2,names)
        x,y = R.gens()
        we = sum([spec[i]*x**(n-i)*y**i for i in range(n+1)])
        return we

    def zeta_polynomial(self,name = "T"):
        r"""
        Returns the Duursma zeta polynomial of the code.

        Assumes minimum_distance(C) > 1 and minimum_distance$(C^\perp) > 1$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = best_known_linear_code(6,3,GF(2))  # long time
            sage: C.minimum_distance()                   # long time (because of above)
            3
            sage: C.zeta_polynomial()                    # long time (because of above)
            2/5*T^2 + 2/5*T + 1/5
            sage: C = HammingCode(4,GF(2))
            sage: C.zeta_polynomial()
            16/429*T^6 + 16/143*T^5 + 80/429*T^4 + 32/143*T^3 + 30/143*T^2 + 2/13*T + 1/13

        REFERENCES:

             I. Duursma, "From weight enumerators to zeta functions",
             in {\bf Discrete Applied Mathematics}, vol. 111, no. 1-2,
             pp. 55-73, 2001.
        """
        n = self.length()
        q = (self.base_ring()).characteristic()
        d = self.minimum_distance()
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print "\n WARNING: There is no guarantee this function works when the minimum distance"
            print "            of the code or of the dual code equals 1.\n"
        from sage.rings.polynomial.polynomial_ring import PolynomialRing
        from sage.rings.fraction_field import FractionField
        from sage.rings.power_series_ring import PowerSeriesRing
        RT = PolynomialRing(QQ,"%s"%name)
        R = PolynomialRing(QQ,3,"xy%s"%name)
        x,y,T = R.gens()
        we = self.weight_enumerator()
        A = R(we)
        B = A(x+y,y,T)-(x+y)**n
        Bs = B.coefficients()
        b = [Bs[i]/binomial(n,i+d) for i in range(len(Bs))]
        r = n-d-dperp+2
        #print B,Bs,b,r
        P_coeffs = []
        for i in range(len(b)):
           if i == 0:
              P_coeffs.append(b[0])
           if i == 1:
              P_coeffs.append(b[1] - (q+1)*b[0])
           if i>1:
              P_coeffs.append(b[i] - (q+1)*b[i-1] + q*b[i-2])
        #print P_coeffs
        P = sum([P_coeffs[i]*T**i for i in range(r+1)])
        return RT(P)

    def chinen_polynomial(self):
        """
        Returns the Chinen zeta polynomial of the code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.chinen_polynomial()       # long time
            (2*sqrt(2)*t^3/5 + 2*sqrt(2)*t^2/5 + 2*t^2/5 + sqrt(2)*t/5 + 2*t/5 + 1/5)/(sqrt(2) + 1)
            sage: C = TernaryGolayCode()
            sage: C.chinen_polynomial()       # long time
            (6*sqrt(3)*t^3/7 + 6*sqrt(3)*t^2/7 + 6*t^2/7 + 2*sqrt(3)*t/7 + 6*t/7 + 2/7)/(2*sqrt(3) + 2)

        This last output agrees with the corresponding example given in Chinen's paper below.

        REFERENCES:
            Chinen, K. "An abundance of invariant polynomials
            satisfying the Riemann hypothesis", April 2007 preprint.

        """
        from sage.rings.polynomial.polynomial_ring import PolynomialRing, polygen
        #from sage.calculus.functional import expand
        from sage.calculus.calculus import sqrt, SymbolicExpressionRing, var
        C = self
        n = C.length()
        RT = PolynomialRing(QQ,2,"Ts")
        T,s = FractionField(RT).gens()
        t = PolynomialRing(QQ,"t").gen()
        Cd = C.dual_code()
        k = C.dimension()
        q = (C.base_ring()).characteristic()
        d = C.minimum_distance()
        dperp = Cd.minimum_distance()
        if dperp > d:
            P = RT(C.zeta_polynomial())
            ## SAGE does not find dealing with sqrt(int) *as an algebraic object*
            ## an easy thing to do. Some tricky gymnastics are used to
            ## make SAGE deal with objects over QQ(sqrt(q)) nicely.
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())*T**(dperp - d)
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())*T**(dperp - d)
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp < d:
            P = RT(C.zeta_polynomial())*T**(d - dperp)
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))
        if dperp == d:
            P = RT(C.zeta_polynomial())
            if is_even(n):
                Pd = q**(k-n/2)*RT(Cd.zeta_polynomial())
            if not(is_even(n)):
                Pd = s*q**(k-(n+1)/2)*RT(Cd.zeta_polynomial())
            CP = P+Pd
            f = CP/CP(1,s)
            return f(t,sqrt(q))

    def zeta_function(self,name = "T"):
        r"""
        Returns the Duursma zeta function of the code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.zeta_function()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
        """
        P =  self.zeta_polynomial()
        q = (self.base_ring()).characteristic()
        RT = PolynomialRing(QQ,"%s"%name)
        T = RT.gen()
        return P/((1-T)*(1-q*T))

    def zeta_function2(self,mode=None):
        r"""
        Returns the Duursma zeta function of the code.

        NOTE: This is somewhat experimental code. It sometimes
        returns "fail" for a reason I don't fully understand.
        However, when it does return a polynomial, the answer is
        (as far as I know) correct.  *Experimental code* included to
        study a particular implementation.

        INPUT:

           mode -- string;
                -- "dual" computes both the zeta function
                   of $C$ and that of $C*$,
                -- "normalized" computes the normalized zeta function
                of $C$, $zeta_C(T)*T(1-genus(C))$ (NOTE: if xi(T,C)
                denotes the normalized zeta function then $xi(T,C*) =
                xi(1/(qT),C)$ is equivalent to $zeta_{C*}(T) =
                zeta_C(1/(qT))*q^(gamma-1)T^(gamma+gamma*-2)$, where
                $gamma = gamma(C)$ and $gamma* = gamma(C*)$.)

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.zeta_function2()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
            sage: C = ExtendedTernaryGolayCode()
            sage: C.zeta_function2()
            (3/7*T^2 + 3/7*T + 1/7)/(3*T^2 - 4*T + 1)

        Both these examples occur in Duursma's paper below.

        REFERENCES:

             I. Duursma, "Weight distributions of geometric Goppa codes,"
             Trans. A.M.S., 351 (1999)3609-3639.

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
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print "\n WARNING: There is no guarantee this function works when the minimum distance\n"
            print "            of the code or of the dual code equals 1.\n"
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

    def zeta_function3(self,mode=None):
        r"""
        Returns the Duursma zeta function of the code.

        NOTE: This sometimes returns "fail" for a reason I don't fully
        understand. However, when it does return a polynomial, the answer is
        (as far as I know) correct. *Experimental code* included to
        study a particular implementation.

        INPUT:
            mode -- string
               mode = "dual" computes both the zeta function of C and that of C*
               mode = "normalized" computes the normalized zeta function of C, zeta_C(T)*T(1-genus(C))
            (NOTE: if xi(T,C) denotes the normalized zeta function
            then xi(T,C*) = xi(1/(qT),C) is equivalent to $zeta_{C*}(T)
            = zeta_C(1/(qT))*q^(gamma-1)T^(gamma+gamma*-2)$, where
            gamma = gamma(C) and gamma* = gamma(C*).)

        EXAMPLES:
            sage.: C = HammingCode(3,GF(2))
            sage.: C.zeta_function3()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
            sage.: C = ExtendedTernaryGolayCode()
            sage.: C.zeta_function3()
            (3/7*T^2 + 3/7*T + 1/7)/(3*T^2 - 4*T + 1)
            sage.: C.zeta_function3(mode="dual")
            ((3/7*T^2 + 3/7*T + 1/7)/(3*T^2 - 4*T + 1),
             (729/7*T^2 + 729/7*T + 243/7)/(729*T^2 - 972*T + 243))

        REFERENCES:
            I. Duursma, "Weight distributions of geometric Goppa codes," Trans. A.M.S., 351 (1999)3609-3639.
        """
        R = PolynomialRing(QQ,3,"xyT")
        F = FractionField(R)
        x,y,T = R.gens()
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print "\n WARNING: There is no guarantee this function works when the minimum distance\n"
            print "            of the code or of the dual code equals 1.\n"
        gammaC = n+1-k-d
        if mode=="dual":
            Cp = self.dual_code()
            dp = Cp.minimum_distance()
            kp = Cp.dimension()
            gammaCp = n+1-kp-dp
        q = self.characteristic()
        W = self.weight_distribution()
        A = add([W[i]*y**i*x**(n-i) for i in range(n+1)])
        f = (x*T+y*(1-T))**n*add([T**i for i in range(n+3)])*add([(q*T)**i for i in range(n+3)])
        Mn = [f.coefficient(T**(n-i)) for i in range(d,n)]
        coeffs = [[(Mn[j]+x**(n)).coefficient(x**i*y**(n-i))(1,1,1) for i in range(n+1)] for j in range(n-d)]
        MS = MatrixSpace(QQ,n-d,n+1)
        #print [type(x[0]) for x in coeffs], MS
        M = MS(coeffs)
        for m in range(2,n-d+1):
            V = VectorSpace(QQ,m)
            v = V([W[i] for i in range(m)])
            Mmm = M.matrix_from_rows_and_columns(range(m),range(m))
            if Mmm.determinant()!=0:
                Pcoeffs = v*Mmm**(-1)
                P = add([Pcoeffs[i]*T**i for i in range(len(Pcoeffs))])
                Z = P*(1-T)**(-1)*(1-q*T)**(-1)
                lhs = P*f
                r = lhs.coefficient(T**(n-d))/((A-x**n)/(q-1))
                if mode=="verbose": print m, P, r
                if r(x,y,T)==r(x+T,y,T) and mode=="dual":
                    Z = Z*r(1,1,1)**(-1)
                    x,y,T = F.gens()
                    Zp = Z(x,y,1/(q*T))*q**(gammaC-1)*T**(gammaC+gammaCp-2)
                    return Z,Zp
                if r(x,y,T)==r(x+T,y,T) and mode=="normalized":
                    Z = Z*r(1,1,1)**(-1)
                    x,y,T = F.gens()
                    Zp = Z(x,y,1/(q*T))*q**(gammaC-1)*T**(gammaC+gammaCp-2)
                    return Z*T**(1-gammaC)
                if r(x,y,T)==r(1,y,T): # check if r is a constant in x ...
                    Z = Z*r(1,1,1)**(-1)
                    return Z
                if mode=="verbose":
                    print "det = ",Mmm.determinant(),"   m = ",m,"  Z = ",Z
            if Mmm.determinant()==0:
                pass
                #print "det = ",Mmm.determinant(),"   m = ",m
        return "fails"

    def punctured(self,L):
        r"""
        Returns the code punctured at the positions L, $L \subset \{1,2,...,n\}$.
        If C is a code of length n in GF(q) then the code $C^L$ obtained from C
        by puncturing at the positions in L is the code of length n-|L|
        consisting of codewords of $C$ which have their $i-th$ coordinate
        deleted if $i \in L$ and left alone if $i\notin L$:
        $$
            C^L = \{(c_{i_1},...,c_{i_N})\ |\ (c_1,...,c_n)\in C\},
        $$
        where $\{1,2,...,n\}-T = \{i_1,...,i_N\}$. In particular, if $L=\{j\}$ then
        $C^L$ is simply the code obtained from $C$ by deleting the $j-th$ coordinate
        of each codeword. The code $C^L$ is called the {\it punctured code at $L$}.
        The dimension of $C^L$ can decrease if $|L|>d-1$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.punctured([1,2])
            Linear code of length 5, dimension 4 over Finite Field of size 2

        """
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        n = len(G.columns())
        k = len(G.rows())
        nL = n-len(L)
        Gstr = str(gap(G))
        gap.eval( "G:="+Gstr )
        C = gap("GeneratorMatCode(G,GF("+str(q)+"))")
        CP = C.PuncturedCode(L)
        Gp = CP.GeneratorMat()
        kL = Gp.Length()  ## this is = dim(CL), usually = k
        G2 = Gp._matrix_(F)
        MS = MatrixSpace(F,kL,nL)
        return LinearCode(MS(G2))

    def shortened(self,L):
        r"""
        Returns the code shortened at the positions L, $L \subset \{1,2,...,n\}$.
        Consider the subcode $C(L)$ consisting of all codewords $c\in C$ which
        satisfy $c_i=0$ for all $i\in L$. The punctured code $C(L)^L$ is called
        the {\it shortened code on $L$} and is denoted $C_L$. The code
        constructed is actually only isomorphic to the shortened code defined
        in this way.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.shortened([1,2])
            Linear code of length 5, dimension 2 over Finite Field of size 2

        """
        G = self.gen_mat()
        F = self.base_ring()
        q = F.order()
        n = len(G.columns())
        k = len(G.rows())
        nL = n-len(L)
        Gstr = str(gap(G))
        gap.eval("G:="+Gstr )
        C = gap("GeneratorMatCode(G,GF("+str(q)+"))")
        Cd = C.DualCode()
        Cdp = Cd.PuncturedCode(L)
        Cdpd = Cdp.DualCode()
        Gs = Cdpd.GeneratorMat()
        kL = Gs.Length()  ## this is = dim(CL), usually = k
        Gss = Gs._matrix_(F)
        MS = MatrixSpace(F,kL,nL)
        return LinearCode(MS(Gss))

    def binomial_moment(self,i):
        r"""
        Returns the i-th binomial moment of the $[n,k,d]_q$-code $C$:
        \[
        B_i(C) = \sum_{S, |S|=i} \frac{q^{k_S}-1}{q-1}
        \]
        where $k_S$ is the dimension of the shortened code $C_{J-S}$, $J=[1,2,...,n]$.
        (The normalized binomial moment is $b_i(C) = binomial(n,d+i)^{-1}B_{d+i}(C)$.)
        In other words, $C_{J-S}$ is the isomorphic to the subcode of C of codewords
        supported on S.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.binomial_moment(2)
            0
            sage: C.binomial_moment(3)    # long time
            0
            sage: C.binomial_moment(4)    # long time
            35

        WARNING: This is slow.

        REFERENCE:
            I. Duursma, "Combinatorics of the two-variable zeta function",
            Finite fields and applications, 109--136, Lecture Notes in Comput. Sci.,
            2948, Springer, Berlin, 2004.
        """
        n = self.length()
        #ii = n-i
        k = self.dimension()
        d = self.minimum_distance()
        F = self.base_ring()
        q = F.order()
        J = range(1,n+1)
        Cp = self.dual_code()
        dp = Cp.minimum_distance()
        if i<d:
            return 0
        if i>n-dp and i<=n:
            return binomial(n,i)*(q**(i+k-n) -1)/(q-1)
        P = SetPartitions(J,2).list()
        b = QQ(0)
        for p in P:
            p = list(p)
            S = p[0]
            if len(S)==n-i:
                C_S = self.shortened(S)
                k_S = C_S.dimension()
                b = b + (q**(k_S) -1)/(q-1)
        return b

    def divisor(self):
        r"""
        Returns the divisor of a code (the divisor is the smallest
        integer $d_0>0$ such that each $A_i>0$ iff $i$ is divisible by $d_0$).

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()
            sage: C.divisor()   ## type 2
            4
            sage: C = ExtendedQuadraticResidueCode(17,GF(2))
            sage: C.divisor()   ## type 1
            2

        """
        C = self
        A = C.spectrum()
        n = C.length()
        V = VectorSpace(QQ,n+1)
        S = V(A).nonzero_positions()
        S0 = [S[i] for i in range(1,len(S))]
        #print S0
        if len(S)>1: return GCD(S0)
        return 1

    def support(self):
        r"""
        Returns the set of indices $j$ where $A_j$ is nonzero,
        where spectrum(self) = $[A_0,A_1,...,A_n]$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.support()
            [0, 3, 4, 7]

        """
        n = self.length()
        F = self.base_ring()
        V = VectorSpace(F,n+1)
        return V(self.spectrum()).support()

    def characteristic_polynomial(self):
        r"""
        Returns the characteristic polynomial of a linear code, as
        defined in van Lint's text [vL].

        EXAMPLES:
            sage: C = ExtendedQuadraticResidueCode(7,GF(2))
            sage: C.characteristic_polynomial()
            -2*x + 16

        REFERENCES:
            van Lint, {\it Introduction to coding theory, 3rd ed.}, Springer-Verlag GTM, 86, 1999.

        """
        R = PolynomialRing(QQ,"x")
        x = R.gen()
        C = self
        Cd = C.dual_code()
        Sd = Cd.support()
        k = C.dimension()
        n = C.length()
        q = (C.base_ring()).order()
        return q**(n-k)*prod([1-x/j for j in Sd if j>0])

    def sd_duursma_data(C,i):
        r"""
        INPUT:
           The formally s.d. code C and the type number (1,2,3,4)
           (does not check if C is actually sd)

        RETURN:
           The data v,m as in Duursama [D]

        EXAMPLES:

        REFERENCES:
            [D] I. Duursma, "Extremal weight enumerators and ultraspherical polynomials"

        """
        n = C.length()
        d = C.minimum_distance()
        if i == 1:
            v = (n-4*d)/2 + 4
            m = d-3
        if i == 2:
            v = (n-6*d)/8 + 3
            m = d-5
        if i == 3:
            v = (n-4*d)/4 + 3
            m = d-4
        if i == 4:
            v = (n-3*d)/2 + 3
            m = d-3
        return [v,m]

    def sd_duursma_q(C,i,d0):
        r"""

        INPUT:
           C -- an sd code (does not check if C is actually an sd code),
           i -- the type number, one of 1,2,3,4,
           d0 -- and the divisor d0 (the smallest integer d0>0 such that each
                 A_i>0 iff i is divisible by d0).

        RETURN:
           The coefficients $q_0, q_1, ...,$  of $q(T)$ as in Duursama [D].

        EXAMPLES:

        REFERENCES:
            [D] I. Duursma, "Extremal weight enumerators and ultraspherical polynomials"

        """
        q = (C.base_ring()).order()
        n = C.length()
        d = C.minimum_distance()
        d0 = C.divisor()
        if i==1 or i==2:
            if d>d0:
                c0 = QQ((n-d)*rising_factorial(d-d0,d0+1)*C.spectrum()[d])/rising_factorial(n-d0-1,d0+2)
            else:
                c0 = QQ((n-d)*C.spectrum()[d])/rising_factorial(n-d0-1,d0+2)
        if i==3 or i==4:
	    if d>d0:
                c0 = rising_factorial(d-d0,d0+1)*C.spectrum()[d]/((q-1)*rising_factorial(n-d0,d0+1))
            else:
                c0 = C.spectrum()[d]/((q-1)*rising_factorial(n-d0,d0+1))
        v = ZZ(C.sd_duursma_data(i)[0])
        m = ZZ(C.sd_duursma_data(i)[1])
        #print m,v,d,d0,c0
        if m<0 or v<0:
            raise ValueError("This case not implemented.")
        PR = PolynomialRing(QQ,"T")
        T = PR.gen()
        if i == 1:
            coefs = PR(c0*(1+3*T+2*T**2)**m*(2*T**2+2*T+1)**v).list()
            #print coefs, len(coefs)
            qc = [coefs[j]/binomial(4*m+2*v,m+j) for j in range(2*m+2*v+1)]
            q = PR(qc)
        if i == 2:
            F = ((T+1)**8+14*T**4*(T+1)**4+T**8)**v
            coefs = (c0*(1+T)**m*(1+4*T+6*T**2+4*T**3)**m*F).coeffs()
            qc = [coefs[j]/binomial(6*m+8*v,m+j) for j in range(4*m+8*v+1)]
            q = PR(qc)
        if i == 3:
            F = (3*T^2+4*T+1)**v*(1+3*T^2)**v
            ## Note that: (3*T^2+4*T+1)(1+3*T^2)=(T+1)**4+8*T**3*(T+1)
            coefs = (c0*(1+3*T+3*T**2)**m*F).coeffs()
            qc = [coefs[j]/binomial(4*m+4*v,m+j) for j in range(2*m+4*v+1)]
            q = PR(qc)
        if i == 4:
            coefs = (c0*(1+2*T)**m*(4*T**2+2*T+1)**v).coeffs()
            qc = [coefs[j]/binomial(3*m+2*v,m+j) for j in range(m+2*v+1)]
            q = PR(qc)
        return q/q(1)

    def sd_zeta_polynomial(C,type=1):
        r"""
        Returns the Duursma zeta function of a self-dual code using
        the construction in [D].

        INPUT:
            type -- type of the s.d. code; one of 1,2,3, or 4.

        EXAMPLES:
            sage: C = ExtendedQuadraticResidueCode(7,GF(2))
            sage: C.sd_zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C.zeta_function()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
            sage: C = ExtendedQuadraticResidueCode(17,GF(2))
            sage: P = C.sd_zeta_polynomial(); P
            8/91*T^8 + 16/91*T^7 + 212/1001*T^6 + 28/143*T^5 + 43/286*T^4 + 14/143*T^3 + 53/1001*T^2 + 2/91*T + 1/182
            sage: P(1)
            1
            sage: C.sd_zeta_polynomial()
            8/91*T^8 + 16/91*T^7 + 212/1001*T^6 + 28/143*T^5 + 43/286*T^4 + 14/143*T^3 + 53/1001*T^2 + 2/91*T + 1/182

        It is a general fact about Duursma zeta polynomials that P(1) = 1.

        REFERENCES:
            [D] I. Duursma, "Extremal weight enumerators and ultraspherical polynomials"
        """
        d0 = C.divisor()
        q = (C.base_ring()).order()
        P = C.sd_duursma_q(type,d0)
        PR = P.parent()
        T = FractionField(PR).gen()
        if type == 1: return P
        if type == 2: return P/(1-2*T+2*T**2)
        if type == 3: return P/(1+3*T**2)
        if type == 4: return P/(1+2*T)

def LinearCodeFromVectorSpace(self):
    """
    """
    F = self.base_ring()
    B = self.basis()
    n = len(B[0].list())
    k = len(B)
    MS = MatrixSpace(F,k,n)
    G = MS([B[i].list() for i in range(k)])
    return LinearCode(G)
