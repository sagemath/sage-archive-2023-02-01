r"""
Linear Codes

VERSION: 1.0

Let $ F$ be a finite field (we denote the finite field with $q$ elements
$GF(q)$ by $\FF_q$). A subspace of $ F^n$ (with the standard basis)
is called a {\it linear code} of length $ n$. If its
dimension is denoted $k$ then we typically store a basis
of $C$ as a $k\times  n$ matrix (the rows are the basis vectors)
called the {\it generator matrix} of $C$.
The rows of the {\it parity check matrix} of $C$ are a basis
for the code,

\[
C^* = \{ v \in GF(q)^n\ |\ v\cdot c = 0,\ for \ all\ c \in C \},
\]
called the {\it dual space} of  $C$.

If $ F=\FF_2$ then $C$ is called a {\it binary code}.
If $ F = \FF_q$ then $C$ is called a {\it $ q$-ary code}.
The elements of a code $C$ are called {\it codewords}.

The symmetric group $S_n$ acts on $F^n$ by permuting coordinates.
If an element $p\in S_n$ sends a code $C$ of length $n$ to itself
(in other words, every codeword of $C$ is sent to some other codeword
of $C$) then $p$ is called a {\it permutation automorphism} of $C$.
The (permutation) automorphism group is denoted $Aut(C)$.

This file contains
\begin{enumerate}
\item
LinearCode class definition; LinearCodeFromVectorspace conversion function,
\item
The spectrum (weight distribution), covering_radius, minimum distance
programs (calling Steve Linton's or CJ Tjhal's C programs), characteristic_function,
and several implementations of the Duursma zeta function
(sd_zeta_polynomial, zeta_polynomial, zeta_function, chinen_polynomial,
for example),
\item
interface with best_known_linear_code_www (interface with codetables.de since
A. Brouwer's online tables have been disabled), bounds_minimum_distance
which call tables in GUAVA (updated May 2006) created by Cen Tjhai instead
of the online internet tables,
\item
gen_mat, list, check_mat, decode, dual_code, extended_code,  shortened, punctured,
genus, binomial_moment, and divisor methods for LinearCode,
\item
Boolean-valued functions such as "==", is_self_dual, is_self_orthogonal,
is_subcode, is_permutation_automorphism, is_permutation_equivalent (which
interfaces with Robert Miller's partition refinement code),
\item
permutation methods: automorphism_group_binary_code,
is_permutation_automorphism, (permutation_automorphism_group is deprecated),
permuted_code, standard_form, module_composition_factors,
\item
design-theoretic methods:
assmus_mattson_designs (implementing Assmus-Mattson Theorem),
\item
code constructions, such as HammingCode and ToricCode, are in a separate
\code{code_constructions.py} module; in the separate \code{guava.py} module,
you will find constructions, such as RandomLinearCodeGuava and
BinaryReedMullerCode, wrapped from the corresponding GUAVA codes.
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
    -- DJ (2006-07): added documentation, group-theoretical methods, ToricCode
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
    -- dj (2008-03): translated punctured, shortened, extended_code, random (and renamed
                 random to random_element), deleted zeta_function2,
                 zeta_function3, added wrapper automorphism_group_binary_code to
                 Robert Miller's code), added direct_sum_code, is_subcode,
                 is_self_dual, is_self_orthogonal, redundancy_matrix, did some
                 alphabetical reorganizing to make the file more readable.
                 Fixed a bug in permutation_automorphism_group which caused it to crash.
    -- dj (2008-03): fixed bugs in spectrum and zeta_polynomial (which misbehaved over
                 non-prime base rings.
    -- dj (2008-10): use CJ Tjhal's MinimumWeight if char = 2 or 3 for min_dist;
                     add is_permutation_equivalent and improve permutation_automorphism_group
                     using an interface with Robert Miller's code; added interface with
                     Leon's code for the spectrum method.


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

import urllib
import sage.modules.free_module as fm
import sage.modules.module as module
from sage.interfaces.all import gap
from sage.rings.finite_field import FiniteField as GF
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix
from sage.rings.arith import GCD, rising_factorial, binomial
from sage.groups.all import SymmetricGroup
from sage.misc.misc import prod
from sage.misc.functional import log, is_even
from sage.rings.rational_field import QQ
from sage.structure.parent_gens import ParentWithGens
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer
from sage.combinat.set_partition import SetPartitions

ZZ = IntegerRing()
VectorSpace = fm.VectorSpace

###################### coding theory functions ##############################

def hamming_weight(v):
    return len(v.nonzero_positions())

def code2leon(C,output):
    r"""
    This is the Python/Sage translation of the GuavaToLeon command in Guava's
    codefun.gi file. It takes a code and outputs a file which represents a code
    readable by Leon's wtdist C program.

    INPUT: C -- a linear code (over GF(p), p<11)
           output = a string name
    OUTPUT: A file output.txt in Sage's tmp director.

    EXAMPLES:
        sage: from sage.misc.misc import SAGE_TMP
        sage: tmpdir = SAGE_TMP
        sage: C = HammingCode(3,GF(2)); C
        Linear code of length 7, dimension 4 over Finite Field of size 2
        sage: foo = sage.coding.linear_code.code2leon(C,"output")
        sage: print open(tmpdir+"output").read()
        LIBRARY code;
        code=seq(2,4,7,seq(
        1,0,0,1,0,1,0,
        0,1,0,1,0,1,1,
        0,0,1,1,0,0,1,
        0,0,0,0,1,1,1
        ));
        FINISH;
        sage: open(tmpdir+"output").close()

    This gets added to SAGE_ROOT.
    """
    from sage.misc.misc import SAGE_TMP
    tmpdir = SAGE_TMP
    F = C.base_ring()
    n = C.length()
    k = C.dimension()
    p = F.order()  # must be prine and <11
    s = "LIBRARY code;\n"+"code=seq(%s,%s,%s,seq(\n"%(p,k,n)
    Gr = [str(r) for r in C.gen_mat().rows()]
    for r in Gr:
        r1 = r.replace("(","")
        r2 = r1 # r2 = r1.replace(",","")
        r3 = r2.replace(")",",\n")
        r4 = r3.replace(" ","")
        s = s+r4
    s = s[:-2]+"\n"
    s = s+"));\n"
    s = s+"FINISH;"
    f = open(tmpdir+output,"w")
    f.write(s)
    f.close()
    return f

def wtdist_gap(Gmat, F):
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
    	sage: sage.coding.linear_code.wtdist_gap(Gstr, F)
    	[1, 0, 0, 7, 7, 0, 0, 1]

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    ALGORITHM: Uses C programs written by Steve Linton in the kernel
        of GAP, so is fairly fast.

    AUTHOR: David Joyner (2005-11)
    """
    G = gap(Gmat)
    q = F.order()
    k = gap(F)
    C = G.GeneratorMatCode(k)
    n = int(C.WordLength())
    z = 'Z(%s)*%s'%(q, [0]*n)     # GAP zero vector as a string
    _ = gap.eval("w:=DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    # for some reason, this commented code doesn't work:
    #dist0 = gap("DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    #v0 = dist0._matrix_(F)
    #print dist0,v0
    #d = G.DistancesDistributionMatFFEVecFFE(k, z)
    v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,n+2)] # because GAP returns vectors in compressed form
    return v

def min_wt_vec_gap(Gmat,F):
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
    	sage: sage.coding.linear_code.min_wt_vec_gap(Gstr,GF(2))
    	(0, 0, 1, 0, 1, 1, 0)

    Here Gstr is a generator matrix of the Hamming [7,4,3] binary code.

    AUTHOR: David Joyner (11-2005)
    """
    from sage.interfaces.gap import gfq_gap_to_sage
    gap.eval("G:="+Gmat)
    C = gap(Gmat).GeneratorMatCode(F)
    n = int(C.WordLength())
    cg = C.MinimumDistanceCodeword()
    c = [gfq_gap_to_sage(cg[j],F) for j in range(1,n+1)]
    V = VectorSpace(F,n)
    return V(c)
    # this older code returns more info but may be slower:
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

def best_known_linear_code_www(n, k, F, verbose=False):
    r"""
    Explains the construction of the best known linear code over
    GF(q) with length n and dimension k, courtesy of the www page
         \url{http://www.codetables.de/}.

    INPUT:
        n -- integer, the length of the code
        k -- integer, the dimension of the code
        F -- finite field, whose field order must be in
                      [2, 3, 4, 5, 7, 8, 9]
        verbose -- bool (default=False), print verbose message

    OUTPUT:
        str -- text about why the bounds are as given

    EXAMPLES:
        sage: L = best_known_linear_code_www(72, 36, GF(2)) # requires internet, optional
        sage: print L                                       # requires internet, optional
        Construction of a linear code
        [72,36,15] over GF(2):
        [1]:  [73, 36, 16] Cyclic Linear Code over GF(2)
             CyclicCode of length 73 with generating polynomial x^37 + x^36 + x^34 +
        x^33 + x^32 + x^27 + x^25 + x^24 + x^22 + x^21 + x^19 + x^18 + x^15 + x^11 +
        x^10 + x^8 + x^7 + x^5 + x^3 + 1
        [2]:  [72, 36, 15] Linear Code over GF(2)
             Puncturing of [1] at 1
        last modified: 2002-03-20


    This function raises an IOError if an error occurs downloading
    data or parsing it.  It raises a ValueError if the q input is
    invalid.

    AUTHOR: (2005-11-14) Steven Sivek, (2008-03) modified by David Joyner
    """
    q = F.order()
    if not q in [2, 3, 4, 5, 7, 8, 9]:
        raise ValueError, "q (=%s) must be in [2,3,4,5,7,8,9]"%q
    n = int(n)
    k = int(k)

    param = ("?q=%s&n=%s&k=%s"%(q,n,k)).replace('L','')

    url = "http://iaks-www.ira.uka.de/home/grassl/codetables/BKLC/BKLC.php"+param
    #url = "http://homepages.cwi.nl/htbin/aeb/lincodbd/"+param
    if verbose:
        print "Looking up the bounds at %s"%url
    f = urllib.urlopen(url)
    s = f.read()
    f.close()
    #print s

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError, "Error parsing data (missing pre tags)."
    text = s[i+5:j].strip()
    return text


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

def self_orthogonal_binary_codes(n, k, b=2, parent=None, BC=None, equal=False,
    in_test=None):
    """
    Returns a Python iterator which generates a complete set of representatives
    of all permutation equivalence classes of self-orthogonal binary linear codes
    of length in [1..n] and dimension in [1..k].

    INPUT:
    n -- maximal length
    k -- maximal dimension
    b -- require that the generators all have weight divisible by b (if b=2, all
        self-orthogonal codes are generated, and if b=4, all doubly even codes
        are generated). Must be an even positive integer.
    parent -- dafault None, used in recursion
    BC -- dafault None, used in recursion
    equal -- default False, if True generates only [n, k] codes
    in_test -- default None, used in recursion

    EXAMPLES:
    Generate all self-orthogonal codes of length up to 7 and dimension up to 3:
        sage: for B in self_orthogonal_binary_codes(7,3):
        ...    print B
        ...
        Linear code of length 2, dimension 1 over Finite Field of size 2
        Linear code of length 4, dimension 2 over Finite Field of size 2
        Linear code of length 6, dimension 3 over Finite Field of size 2
        Linear code of length 4, dimension 1 over Finite Field of size 2
        Linear code of length 6, dimension 2 over Finite Field of size 2
        Linear code of length 6, dimension 2 over Finite Field of size 2
        Linear code of length 7, dimension 3 over Finite Field of size 2
        Linear code of length 6, dimension 1 over Finite Field of size 2

    Generate all doubly-even codes of length up to 7 and dimension up to 3:
        sage: for B in self_orthogonal_binary_codes(7,3,4):
        ...    print B; print B.gen_mat()
        ...
        Linear code of length 4, dimension 1 over Finite Field of size 2
        [1 1 1 1]
        Linear code of length 6, dimension 2 over Finite Field of size 2
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]
        Linear code of length 7, dimension 3 over Finite Field of size 2
        [1 0 1 1 0 1 0]
        [0 1 0 1 1 1 0]
        [0 0 1 0 1 1 1]

    Generate all doubly-even codes of length up to 7 and dimension up to 2:
        sage: for B in self_orthogonal_binary_codes(7,2,4):
        ...    print B; print B.gen_mat()
        Linear code of length 4, dimension 1 over Finite Field of size 2
        [1 1 1 1]
        Linear code of length 6, dimension 2 over Finite Field of size 2
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]

    Generate all self-orthogonal codes of length equal to 8 and dimension
    equal to 4:
        sage: for B in self_orthogonal_binary_codes(8, 4, equal=True):
        ...     print B; print B.gen_mat()
        Linear code of length 8, dimension 4 over Finite Field of size 2
        [1 0 0 1 0 0 0 0]
        [0 1 0 0 1 0 0 0]
        [0 0 1 0 0 1 0 0]
        [0 0 0 0 0 0 1 1]
        Linear code of length 8, dimension 4 over Finite Field of size 2
        [1 0 0 1 1 0 1 0]
        [0 1 0 1 1 1 0 0]
        [0 0 1 0 1 1 1 0]
        [0 0 0 1 0 1 1 1]

    Since all the codes will be self-orthogonal, b must be divisible by 2:
        sage: list(self_orthogonal_binary_codes(8, 4, 1, equal=True))
        Traceback (most recent call last):
        ...
        ValueError: b (1) must be a positive even integer.

    """
    d=int(b)
    if d!=b or d%2==1 or d <= 0:
        raise ValueError("b (%s) must be a positive even integer."%b)
    from binary_code import BinaryCode, BinaryCodeClassifier
    if k < 1 or n < 2:
        return
    if equal:
        in_test = lambda M : (M.ncols() - M.nrows()) <= (n-k)
        out_test = lambda C : (C.dimension() == k) and (C.length() == n)
    else:
        in_test = lambda M : True
        out_test = lambda C : True
    if BC is None:
        BC = BinaryCodeClassifier()
    if parent is None:
        for j in xrange(d, n+1, d):
            M = Matrix(GF(2), [[1]*j])
            if in_test(M):
                for N in self_orthogonal_binary_codes(n, k, d, M, BC, in_test=in_test):
                    if out_test(N): yield N
    else:
        C = LinearCode(parent)
        if out_test(C): yield C
        if k == parent.nrows():
            return
        for nn in xrange(parent.ncols()+1, n+1):
            if in_test(parent):
                for child in BC.generate_children(BinaryCode(parent), nn, d):
                    for N in self_orthogonal_binary_codes(n, k, d, child, BC, in_test=in_test):
                        if out_test(N): yield N

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

    def _repr_(self):
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def automorphism_group_binary_code(self):
        r"""
        This only applies to linear binary codes and returns its
        (permutation) automorphism group. In other words, if
        the code C has length $n$ then it returns the subgroup of the
        symmetric group $S_n$:
        \[
        \{ g in S_n\ |\ g(c) \in C, \forall c\in C\},
        \]
        where $S_n$ acts on $GF(2)^n$ by permuting coordinates.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: G = C.automorphism_group_binary_code(); G
            Permutation Group with generators [(3,4)(5,6), (3,5)(4,6), (2,3)(5,7), (1,2)(5,6)]
            sage: G.order()
            168

        """
        C = self
        F = C.base_ring()
        if F!=GF(2):
            raise NotImplementedError, "Only implemented for binary codes."
        from sage.coding.binary_code import BinaryCode, BinaryCodeClassifier
        genmat = C.gen_mat()
        B = BinaryCode(genmat)
        BC = BinaryCodeClassifier()
        autgp = BC._aut_gp_and_can_label(B)
        if autgp[0] == []:
            return PermutationGroup([()])
        Sn = SymmetricGroup(C.length())
        L = [[j+1 for j in autgp[0][i]] for i in range(len(autgp[0]))]
        G = PermutationGroup([Sn(x) for x in L])
        return G

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
        lambda blocks the block design is called a {\bf
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
            sage: blocks = [c.support() for c in C if hamming_weight(c)==8]; len(blocks)  # long time computation
            759

        REFERENCE:
            [HP] W. C. Huffman and V. Pless, {\bf Fundamentals of ECC}, Cambridge Univ. Press, 2003.
        """
        C = self
        ans = []
        G = C.gen_mat()
        n = len(G.columns())
        Cp = C.dual_code()
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

    def basis(self):
        return self.__gens

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

    def __contains__(self,v):
        A = self.ambient_space()
        C = A.subspace(self.gens())
        return C.__contains__(v)

    def characteristic(self):
        return (self.base_ring()).characteristic()

    def characteristic_polynomial(self):
        r"""
        Returns the characteristic polynomial of a linear code, as
        defined in van Lint's text [vL].

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()
            sage: C.characteristic_polynomial()
            -4/3*x^3 + 64*x^2 - 2816/3*x + 4096

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

    def chinen_polynomial(self):
        """
        Returns the Chinen zeta polynomial of the code.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.chinen_polynomial()       # long time
            (2*sqrt(2)*t^3/5 + 2*sqrt(2)*t^2/5 + 2*t^2/5 + sqrt(2)*t/5 + 2*t/5 + 1/5)/(sqrt(2) + 1)
            sage: C = TernaryGolayCode()
            sage: C.chinen_polynomial()       # long time
            (3*sqrt(3)*t^3/7 + 3*sqrt(3)*t^2/7 + 3*t^2/7 + sqrt(3)*t/7 + 3*t/7 + 1/7)/(sqrt(3) + 1)

        This last output agrees with the corresponding example given in Chinen's paper below.

        REFERENCES:
            Chinen, K. "An abundance of invariant polynomials
            satisfying the Riemann hypothesis", April 2007 preprint.

        """
        from sage.calculus.calculus import sqrt
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
            # SAGE does not find dealing with sqrt(int) *as an algebraic object*
            # an easy thing to do. Some tricky gymnastics are used to
            # make SAGE deal with objects over QQ(sqrt(q)) nicely.
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

    def __cmp__(self, right):
        if not isinstance(right, LinearCode):
            return cmp(type(self), type(right))
        return cmp(self.__gen_mat, right.__gen_mat)

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

    def covering_radius(self):
        r"""
        Wraps Guava's CoveringRadius command.

        The {\bf covering radius} of a linear code C is the smallest number r with the
        property that each element v of the ambient vector space of C has at
        most a distance r to the code C. So for each vector v there must be
        an element c of C with $d(v,c) \leq  r$. A binary linear code
        with reasonable small covering radius is often referred to as a {\bf covering code}.

        For example, if C is a perfect code, the covering radius is equal to t, the
        number of errors the code can correct, where d = 2t+1, with d the
        minimum distance of C

        EXAMPLES:
            sage: C = HammingCode(5,GF(2))
            sage: C.covering_radius()
            1

        """
        F = self.base_ring()
        G = self.gen_mat()
        gapG = gap(G)
        C = gapG.GeneratorMatCode(gap(F))
        r = C.CoveringRadius()
        try:
            return ZZ(r)
        except:
            raise ValueError("Sorry, the covering radius of this code cannot be computed by Guava.")

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

    def divisor(self):
        r"""
        Returns the divisor of a code (the divisor is the smallest
        integer $d_0>0$ such that each $A_i>0$ iff $i$ is divisible by $d_0$).

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()
            sage: C.divisor()   # Type II self-dual
            4
            sage: C = QuadraticResidueCodeEvenPair(17,GF(2))[0]
            sage: C.divisor()
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

    def dimension(self):
        return self.__dim

    def direct_sum(self,other):
        """
        C1, C2 must be linear codes defined over the same base ring.
        Returns the (usual vector space) direct sum of the codes.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2))
            sage: C2 = C1.direct_sum(C1); C2
            Linear code of length 14, dimension 8 over Finite Field of size 2
            sage: C3 = C1.direct_sum(C2); C3
            Linear code of length 21, dimension 12 over Finite Field of size 2

        """
        C1 = self; C2 = other
        G1 = C1.gen_mat()
        G2 = C2.gen_mat()
        F = C1.base_ring()
        n1 = len(G1.columns())
        k1 = len(G1.rows())
        n2 = len(G2.columns())
        k2 = len(G2.rows())
        MS1 = MatrixSpace(F,k2,n1)
        MS2 = MatrixSpace(F,k1,n2)
        Z1 = MS1(0)
        Z2 = MS2(0)
        top = G1.augment(Z2)
        bottom = Z1.augment(G2)
        G = top.stack(bottom)
        return LinearCode(G)

    def __eq__(self, right):
        """
        Checks if self == right.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2))
            sage: C2 = HammingCode(3,GF(2))
            sage: C1 == C2
            True
            sage: C2 = C1.extended_code()
            sage: C3 = C2.punctured([7])
            sage: C1 == C3
            True

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
        k = len(G.rows())
        MS1 = MatrixSpace(F,k,1)
        ck_sums = [-sum(G.rows()[i]) for i in range(k)]
        last_col = MS1(ck_sums)
        Gx = G.augment(last_col)
        return LinearCode(Gx)

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
            sage: c = C.basis()[1]
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
        a = log(q,q0)  # test if F/F0 is a field extension
        if not isinstance(a, Integer):
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

    def gen_mat(self):
        r"""
        Returns the generator matrix of the code.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2))
            sage: C1.gen_mat()
            [1 0 0 1 0 1 0]
            [0 1 0 1 0 1 1]
            [0 0 1 1 0 0 1]
            [0 0 0 0 1 1 1]
            sage: C2 = HammingCode(2,GF(4,"a"))
            sage: C2.gen_mat()
            [    1     0     0     1     1]
            [    0     1     0     1 a + 1]
            [    0     0     1     1     a]

        """
        return self.__gen_mat

    def gens(self):
        return self.__gens

    def genus(self):
        r"""
        Returns the "Duursma genus" of the code, gamma_C = n+1-k-d.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2)); C1
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C1.genus()
            1
            sage: C2 = HammingCode(2,GF(4,"a")); C2
            Linear code of length 5, dimension 3 over Finite Field in a of size 2^2
            sage: C2.genus()
            0

        Since all Hamming codes have minimum distance 3, these computations
        agree with the definition, n+1-k-d.

        """
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-k-d
        return gammaC

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
                return False
        return True

    def is_permutation_equivalent(self,other,method=None):
        """
        Returns true if self and other are permutation equivalent codes
        and false otherwise. The method="verbose" option also
        returns a permutation (if true) sending self to other.
        Uses Robert Miller's double coset partition refinement work.

        EXAMPLES:
            sage: P.<x> = PolynomialRing(GF(2),"x")
            sage: g = x^3+x+1
            sage: C1 = CyclicCodeFromGeneratingPolynomial(7,g); C1
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C2 = HammingCode(3,GF(2)); C2
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C1.is_permutation_equivalent(C2)
            True
            sage: C1.is_permutation_equivalent(C2,method="verbose")
            (True, (4,6,5,7))
            sage: C1 = RandomLinearCode(10,5,GF(2))
            sage: C2 = RandomLinearCode(10,5,GF(3))
            sage: C1.is_permutation_equivalent(C2)
            False

        """
        from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct
        F = self.base_ring()
        F_o = other.base_ring()
        q = F.order()
        G = self.gen_mat()
        n = self.length()
        n_o = other.length()
        if F != F_o or n != n_o:
            return False
        k = len(G.rows())
        MS = MatrixSpace(F,q**k,n)
        CW1 = MS(self.list())
        CW2 = MS(other.list())
        B1 = NonlinearBinaryCodeStruct(CW1)
        B2 = NonlinearBinaryCodeStruct(CW2)
        ans = B1.is_isomorphic(B2)
        if ans!=False:
            if method=="verbose":
                Sn = SymmetricGroup(n)
                return True, Sn([i+1 for i in ans])**(-1)
            return True
        return False

    def is_self_dual(self):
        """
        A code C is self-dual if C == C.dual_code() is True.

        Returns True if the code is self-dual (in the
        usual Hamming inner product) and False otherwise.

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()
            sage: C.is_self_dual()
            True
            sage: C = HammingCode(3,GF(2))
            sage: C.is_self_dual()
            False

        """
        C = self
        Cd = C.dual_code()
        return (C == Cd)

    def is_self_orthogonal(C):
        """
        A code C is self-orthogonal if C is a subcode of C.dual_code().

        Returns True if the code is self-dual (in the
        usual Hamming inner product) and False otherwise.

        EXAMPLES:
            sage: C = ExtendedBinaryGolayCode()
            sage: C.is_self_orthogonal()
            True
            sage: C = HammingCode(3,GF(2))
            sage: C.is_self_orthogonal()
            False
            sage: C = QuasiQuadraticResidueCode(11)
            sage: C.is_self_orthogonal()
            True

        """
        Cd = C.dual_code()
        return C.is_subcode(Cd)

    def is_galois_closed(self):
        r"""
        Checks if $C$ is equal to its Galois closure.
        """
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def is_subcode(self,other):
        """
        Returns true if the first is a subcode of the second.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2))
            sage: G1 = C1.gen_mat()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
            sage: C1.is_subcode(C2)
            False
            sage: C3 = C1.extended_code()
            sage: C1.is_subcode(C3)
            False
            sage: C4 = C1.punctured([1])
            sage: C4.is_subcode(C1)
            False
            sage: C5 = C1.shortened([1])
            sage: C5.is_subcode(C1)
            False
            sage: C1 = HammingCode(3,GF(9,"z"))
            sage: G1 = C1.gen_mat()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True

        """
        C1 = self; C2 = other
        G = C1.gen_mat()
        for r in G.rows():
            if not(r in C2):
                return False
        return True

    def length(self):
        return self.__length

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
        return self.gen_mat().row_space().list()

    def minimum_distance(self):
        r"""
        If q is not 2 or 3 then this uses a GAP kernel function (in C) written
        by Steve Linton. If q is 2 or 3 then this uses a very fast program
        written in C written by CJ Tjhal (this is much faster, except in some
        small examples).

        EXAMPLES:
            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3
            sage: C = HammingCode(2,GF(4,"a")); C
            Linear code of length 5, dimension 3 over Finite Field in a of size 2^2
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
        if q == 2 or q == 3:
            C = gapG.GeneratorMatCode(gap(F))
            d = C.MinimumWeight()
            #print "Running Guava's MinimumWeight ...\n"
            return ZZ(d)
        Gstr = "%s*Z(%s)^0"%(gapG, q)
        return hamming_weight(min_wt_vec_gap(Gstr,F))

    def module_composition_factors(self,gp):
        r"""
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.automorphism_group_binary_code()

        Now type "C.module_composition_factors(gp)" to get the record printed.

        """
        F = self.base_ring()
        q = F.order()
        gens = gp.gens()
        G = self.gen_mat()
        n = len(G.columns())
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

    def permutation_automorphism_group(self,method="partition"):
        r"""
        If $C$ is an $[n,k,d]$ code over $F$, this function computes
        the subgroup $Aut(C) \subset S_n$ of all permutation
        automorphisms of $C$. The binary case always uses the (default) partition
        refinement method of Robert Miller.

        Options:
           If method="gap" then GAP's MatrixAutomorphism function (written by
            Thomas Breuer) is used. The implementation combines an idea of mine
            with an improvement suggested by Cary Huffman.
           If method="gap+verbose" then code-theoretic data is printed out at
            several stages of the computation.
           If method="partition" then the (default) partition refinement method
            of Robert Miller is used.

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
            sage: M11.order()
            7920
            sage: G = C.permutation_automorphism_group()  # this should take < 5 seconds
            sage: G.is_isomorphic(M11)                    # this should take < 5 seconds
            True

        In the binary case, uses sage.coding.binary_code:

            sage: C = ExtendedBinaryGolayCode()
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            244823040

        In the non-binary case:

            sage: C = HammingCode(2,GF(3)); C
            Linear code of length 4, dimension 2 over Finite Field of size 3
            sage: C.permutation_automorphism_group(method="partition")
            Permutation Group with generators [(1,2,3)]
            sage: C = HammingCode(2,GF(4,"z")); C
            Linear code of length 5, dimension 3 over Finite Field in z of size 2^2
            sage: C.permutation_automorphism_group(method="partition")
            Permutation Group with generators [(1,2)(3,4), (1,3)(2,4)]
            sage: C.permutation_automorphism_group(method="gap")
            Permutation Group with generators [(1,2)(3,4), (1,3)(2,4)]
            sage: C = TernaryGolayCode()
            sage: C.permutation_automorphism_group(method="gap")
            Permutation Group with generators [(3,4)(5,7)(6,9)(8,11), (3,5,8)(4,11,7)(6,9,10), (2,3)(4,6)(5,8)(7,10), (1,2)(4,11)(5,8)(9,10)]

        However, the option \code{method="gap+verbose"}, will print out

             Minimum distance: 5
             Weight distribution:
             [1, 0, 0, 0, 0, 132, 132, 0, 330, 110, 0, 24]

             Using the 132 codewords of weight 5
             Supergroup size:
             39916800

        in addition to the output of C.permutation_automorphism_group(method="gap").

        """
        F = self.base_ring()
        q = F.order()
        if q == 2 and self.length() <= 64:
            from sage.coding.binary_code import WORD_SIZE
            if self.length() <= WORD_SIZE:
                return self.automorphism_group_binary_code()
        G = self.gen_mat()
        n = len(G.columns())
        k = len(G.rows())
        wts = self.spectrum()                                            # bottleneck 1
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0]
        Sn = SymmetricGroup(n)
        if "gap" in method:
            Gp = gap("SymmetricGroup(%s)"%n)               # initializing G in gap
            Gstr = str(gap(G))
            gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
            gap.eval("eltsC:=Elements(C)")
            if method=="gap+verbose":
                print "\n Minimum distance: %s \n Weight distribution: \n %s"%(nonzerowts[1],wts)
            stop = 0                                          # only stop if all gens are autos
            for i in range(1,len(nonzerowts)):
                if stop == 1:
                    break
                wt = nonzerowts[i]
                if method=="gap+verbose":
                    size = Gp.Size()
                    print "\n Using the %s codewords of weight %s \n Supergroup size: \n %s\n "%(wts[wt],wt,size)
                gap.eval("Cwt:=Filtered(eltsC,c->WeightCodeword(c)=%s)"%wt)   # bottleneck 2 (repeated
                gap.eval("matCwt:=List(Cwt,c->VectorCodeword(c))")            #        for each i until stop = 1)
                A = gap("MatrixAutomorphisms(matCwt)")
                #print "A = ",A, "\n Gp = ", Gp, "\n strGp = ", str(Gp)
                G2 = gap("Intersection2(%s,%s)"%(str(A).replace("\n",""),str(Gp).replace("\n",""))) #  bottleneck 3
                Gp = G2
                if Gp.Size()==1:
                    return PermutationGroup([()])
                autgp_gens = Gp.GeneratorsOfGroup()
                gens = [Sn(str(x).replace("\n","")) for x in autgp_gens]
                stop = 1                         # get ready to stop
                for x in gens:                   # if one of these gens is not an auto then don't stop
                    if not(self.is_permutation_automorphism(x)):
                        stop = 0
                        break
            G = PermutationGroup(gens)
            return G
        if method=="partition":
            from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
            stop = 0                                          # only stop if all gens are autos
            for i in range(1,len(nonzerowts)):
                if stop == 1:
                    break
                wt = nonzerowts[i]
                Cwt = [c for c in self if hamming_weight(c)==wt] # ridiculously slow!!
                MS = MatrixSpace(F,len(Cwt),n)
                Cwords_wt = MS(Cwt)
                M = MatrixStruct(Cwords_wt)
                autgp = M.automorphism_group()
                if autgp[0] == []:
                    return PermutationGroup([()])
                L = [[j+1 for j in autgp[0][i]] for i in range(len(autgp[0]))]
                G = PermutationGroup([Sn(x) for x in L])
                return G
        raise NotImplementedError("The only methods implemented currently are 'gap', 'gap+verbose', and 'partition'.")

    def permuted_code(self,p):
        r"""
        Returns the permuted code -- the code $C$ which is equivalent
        to self via the column permutation $p$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: G = C.automorphism_group_binary_code(); G
            Permutation Group with generators [(3,4)(5,6), (3,5)(4,6), (2,3)(5,7), (1,2)(5,6)]
            sage: g = G("(2,3)(5,7)")
            sage: Cg = C.permuted_code(g)
            sage: Cg
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C == Cg
            True

        """
        F = self.base_ring()
        G = self.gen_mat()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        Gp = G*MS(p.matrix().rows())
        return LinearCode(Gp)

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
        nC = len(G.columns())
        kC = len(G.rows())
        nL = nC-len(L)
        colsGnotL = [G.columns()[i] for i in range(nC) if not(i in L)]
        MS1 = MatrixSpace(F, kC, nL)
        GL1 = MS1(list(colsGnotL))
        VnotL = GL1.row_space()
        B = VnotL.basis()
        k = len(list(B))
        MS2 = MatrixSpace(F, k, nL)
        GL2 = MS2(list(B))
        return LinearCode(GL2)

    def random_element(self):
        """
        Returns a random codeword.

        EXAMPLES:
            sage: C = HammingCode(3,GF(4,'a'))
            sage: Cc = C.galois_closure(GF(2))
            sage: c = C.gen_mat()[1]
            sage: V = VectorSpace(GF(4,'a'),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False
            sage: c2 in Cc
            True

        """
        V = self.ambient_space()
        S = V.subspace(self.basis())
        return S.random_element()

    def redundancy_matrix(C):
        """
        If C is a linear [n,k,d] code then this function returns a
        kx(n-k) matrix A such that G = (I,A) generates a code (in standard
        form) equiv to C. If C is already in standard form and G = (I,A)
        is its gen mat then this function simply returns that A.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.gen_mat()
            [1 0 0 1 0 1 0]
            [0 1 0 1 0 1 1]
            [0 0 1 1 0 0 1]
            [0 0 0 0 1 1 1]
            sage: C.redundancy_matrix()
            [1 1 0]
            [1 1 1]
            [1 0 1]
            [0 1 1]
            sage: C.standard_form()[0].gen_mat()
            [1 0 0 0 1 1 0]
            [0 1 0 0 1 1 1]
            [0 0 1 0 1 0 1]
            [0 0 0 1 0 1 1]
            sage: C = HammingCode(2,GF(3))
            sage: C.gen_mat()
            [1 0 2 2]
            [0 1 2 1]
            sage: C.redundancy_matrix()
            [2 2]
            [2 1]

        """
        n = C.length()
        k = C.dimension()
        C1 = C.standard_form()[0]
        G1 = C1.gen_mat()
        A = G1.matrix_from_columns(range(k,n))
        return A

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
        if m<0 or v<0:
            raise ValueError("This case not implemented.")
        PR = PolynomialRing(QQ,"T")
        T = PR.gen()
        if i == 1:
            coefs = PR(c0*(1+3*T+2*T**2)**m*(2*T**2+2*T+1)**v).list()
            qc = [coefs[j]/binomial(4*m+2*v,m+j) for j in range(2*m+2*v+1)]
            q = PR(qc)
        if i == 2:
            F = ((T+1)**8+14*T**4*(T+1)**4+T**8)**v
            coefs = (c0*(1+T)**m*(1+4*T+6*T**2+4*T**3)**m*F).coeffs()
            qc = [coefs[j]/binomial(6*m+8*v,m+j) for j in range(4*m+8*v+1)]
            q = PR(qc)
        if i == 3:
            F = (3*T^2+4*T+1)**v*(1+3*T^2)**v
            # Note that: (3*T^2+4*T+1)(1+3*T^2)=(T+1)**4+8*T**3*(T+1)
            coefs = (c0*(1+3*T+3*T**2)**m*F).coeffs()
            qc = [coefs[j]/binomial(4*m+4*v,m+j) for j in range(2*m+4*v+1)]
            q = PR(qc)
        if i == 4:
            coefs = (c0*(1+2*T)**m*(4*T**2+2*T+1)**v).coeffs()
            qc = [coefs[j]/binomial(3*m+2*v,m+j) for j in range(m+2*v+1)]
            q = PR(qc)
        return q/q(1)

    def sd_zeta_polynomial(C,typ=1):
        r"""
        Returns the Duursma zeta function of a self-dual code using
        the construction in [D].

        INPUT:
            typ -- type of the s.d. code; one of 1,2,3, or 4.

        EXAMPLES:
            sage: C1 = HammingCode(3,GF(2))
            sage: C2 = C1.extended_code(); C2
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: C2.is_self_dual()
            True
            sage: C2.sd_zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C2.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: P = C2.sd_zeta_polynomial(); P(1)
            1
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.sd_zeta_polynomial(4)
            1

        It is a general fact about Duursma zeta polynomials that P(1) = 1.

        REFERENCES:
            [D] I. Duursma, "Extremal weight enumerators and ultraspherical polynomials"
        """
        d0 = C.divisor()
        P = C.sd_duursma_q(typ,d0)
        PR = P.parent()
        T = FractionField(PR).gen()
        if typ == 1:
            P0 = P
        if typ == 2:
            P0 = P/(1-2*T+2*T**2)
        if typ == 3:
            P0 = P/(1+3*T**2)
        if typ == 4:
            P0 = P/(1+2*T)
        return P0/P0(1)

    def shortened(self,L):
        r"""
        Returns the code shortened at the positions L, $L \subset \{1,2,...,n\}$.
        Consider the subcode $C(L)$ consisting of all codewords $c\in C$ which
        satisfy $c_i=0$ for all $i\in L$. The punctured code $C(L)^L$ is called
        the {\it shortened code on $L$} and is denoted $C_L$. The code
        constructed is actually only isomorphic to the shortened code defined
        in this way.

        By Theorem 1.5.7 in [HP], $C_L$ is $((C^\perp)^L)^\perp$. This is used in the
        construction below.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.shortened([1,2])
            Linear code of length 5, dimension 2 over Finite Field of size 2

        """
        Cd = self.dual_code()
        Cdp = Cd.punctured(L)
        Cdpd = Cdp.dual_code()
        Gs = Cdpd.gen_mat()
        return LinearCode(Gs)

    def spectrum(self, method="gap"):
        r"""
        The default method (gap) uses a GAP kernel function (in C) written by Steve Linton.

        WARNING:
           The optional method (leon) may create a stack smashing error and a
           traceback but should return the correct answer. It appears to run
           much faster than the "gap" method in some ("small") examples
           and much slower than the "gap" method in other ("larger") examples.

        EXAMPLES:
            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: F.<z> = GF(2^2,"z")
            sage: C = HammingCode(2, F); C
            Linear code of length 5, dimension 3 over Finite Field in z of size 2^2
            sage: C.spectrum()
            [1, 0, 0, 30, 15, 18]
            sage: C = HammingCode(3,GF(2)); C
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C.spectrum(method="leon")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C = HammingCode(3,GF(3)); C
            Linear code of length 13, dimension 10 over Finite Field of size 3
            sage: C.spectrum() == C.spectrum(method="leon")
            True
            #[1, 0, 0, 104, 468, 1404, 4056, 8424, 11934, 13442, 11232, 5616, 2080, 288]
            sage: C = HammingCode(2,GF(5)); C
            Linear code of length 6, dimension 4 over Finite Field of size 5
            sage: C.spectrum() == C.spectrum(method="leon")
            True
            #[1, 0, 0, 80, 120, 264, 160]
            sage: C = HammingCode(2,GF(7)); C
            Linear code of length 8, dimension 6 over Finite Field of size 7
            sage: C.spectrum() == C.spectrum(method="leon")
            True
            #[1, 0, 0, 336, 1680, 9072, 26544, 45744, 34272]

        NOTE:
          The GAP command DirectoriesPackageLibrary tells the location of the latest
          version of the Guava libraries, so gives us the location of the Guava binaries too.

        """
        from sage.misc.misc import SAGE_TMP, SAGE_ROOT, tmp_filename
        import commands
        n = self.length()
        F = self.base_ring()
        G = self.gen_mat()
        if method=="gap":
            Gstr = G._gap_init_()
            spec = wtdist_gap(Gstr,F)
            return spec
        if method=="leon":
            if not(F.order() in [2,3,5,7]):
                raise NotImplementedError("The method 'leon' is only implemented for q = 2,3,5,7.")
            wts = []
            code2leon(self,"incode")
            rt = SAGE_ROOT
            #tmp = SAGE_TMP
            tmp_file = tmp_filename()
            #print tmp_file
            #pth = rt+"/local/lib/gap-4.4.10/pkg/guava3.4/bin/"
            spth = gap.eval('DirectoriesPackageLibrary( "guava" )')
            pth = spth[7:][:-8]+"bin/"
            a = commands.getoutput(pth+"wtdist "+rt+"/incode::code > "+tmp_file)
            f = open(tmp_file)
            lines = f.readlines()
            f.close()
            s = 0
            for L in lines:  # find the line where the numbers start
                s = s+1
                if len(L)>2 and L[-2] == "-":
                    break
            for L in lines[s:]:
                N = len(L)
                for m in range(1,N):
                    if L[-m]==" ":
                        break
                if " " in L[-m:]:
                    fs = L[-m:].replace(" ","")
                if " " in L[:N-m]:
                    ws = L[:N-m].replace(" ","")
                wts.append([eval(ws),eval(fs)])
            Wts = [0]*(n+1)
            for x in wts:
                i = int(x[0])
                w = int(x[1])
                if len(x)==2 and (i in range(n+1)):
                    Wts[i] = w
            #print wts, Wts
            return Wts
        raise NotImplementedError("The only methods implemented currently are 'gap' and 'leon'.")

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
            sage: C = HammingCode(3,GF(2))
            sage: C.gen_mat()
            [1 0 0 1 0 1 0]
            [0 1 0 1 0 1 1]
            [0 0 1 1 0 0 1]
            [0 0 0 0 1 1 1]
            sage: Cs,p = C.standard_form()
            sage: p
            (4,5)
            sage: Cs
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: Cs.gen_mat()
            [1 0 0 0 1 1 0]
            [0 1 0 0 1 1 1]
            [0 0 1 0 1 0 1]
            [0 0 0 1 0 1 1]
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
        R = PolynomialRing(QQ,2,names)
        x,y = R.gens()
        we = sum([spec[i]*x**(n-i)*y**i for i in range(n+1)])
        return we

    def zeta_polynomial(self,name = "T"):
        r"""
        Returns the Duursma zeta polynomial of the code C.

        Assumes \code{C.minimum_distance()} > 1 and minimum_distance$(C^\perp) > 1$.

        EXAMPLES:
            sage: C = HammingCode(3,GF(2))
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = best_known_linear_code(6,3,GF(2))
            sage: C.minimum_distance()
            3
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = HammingCode(4,GF(2))
            sage: C.zeta_polynomial()
            16/429*T^6 + 16/143*T^5 + 80/429*T^4 + 32/143*T^3 + 30/143*T^2 + 2/13*T + 1/13
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.zeta_polynomial()
            1

        REFERENCES:

             I. Duursma, "From weight enumerators to zeta functions",
             in {\bf Discrete Applied Mathematics}, vol. 111, no. 1-2,
             pp. 55-73, 2001.
        """
        n = self.length()
        q = (self.base_ring()).order()
        d = self.minimum_distance()
        dperp = (self.dual_code()).minimum_distance()
        if d == 1 or dperp == 1:
            print "\n WARNING: There is no guarantee this function works when the minimum distance"
            print "            of the code or of the dual code equals 1.\n"
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
        return RT(P)/RT(P)(1)

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

    weight_distribution = spectrum


def LinearCodeFromVectorSpace(self):
    """
    Simply converts a vector subspace V of $GF(q)^n$ into a LinearCode.
    """
    F = self.base_ring()
    B = self.basis()
    n = len(B[0].list())
    k = len(B)
    MS = MatrixSpace(F,k,n)
    G = MS([B[i].list() for i in range(k)])
    return LinearCode(G)
