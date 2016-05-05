# -*- coding: utf-8 -*-
r"""
Linear code

VERSION: 1.2

Let `F` be a finite field.  Here, we will denote the finite field with `q`
elements by `\GF{q}`.  A subspace of `F^n` (with the standard basis) is
called a linear code of length `n`.  If its dimension is denoted `k` then we
typically store a basis of `C` as a `k \times n` matrix, with rows the basis
vectors.  It is called the generator matrix of `C`. The rows of the parity
check matrix of `C` are a basis for the code,

.. math::

      C^* = \{ v \in GF(q)^n\ |\ v\cdot c = 0,\ for \ all\ c \in C \},

called the dual space of `C`.

If `F=\GF{2}` then `C` is called a binary code.  If `F = \GF{q}` then `C` is
called a `q`-ary code.  The elements of a code `C` are called codewords.

Let `C`, `D` be linear codes of length `n` and dimension `k`. There are
several notions of equivalence for linear codes:

`C` and `D` are

    - permutational equivalent, if there is some permutation `\pi \in S_n`
      such that `(c_{\pi(0)}, \ldots, c_{\pi(n-1)}) \in D` for all `c \in C`.

    - linear equivalent, if there is some permutation `\pi \in S_n` and a
      vector `\phi` of units of length `n` such that
      `(c_{\pi(0)} \phi_0^{-1}, \ldots, c_{\pi(n-1)} \phi_{n-1}^{-1}) \in D`
      for all `c \in C`.

    - semilinear equivalent, if there is some permutation `\pi \in S_n`, a
      vector `\phi` of units of length `n` and a field automorphism `\alpha`
      such that
      `(\alpha(c_{\pi(0)}) \phi_0^{-1}, \ldots, \alpha( c_{\pi(n-1)}) \phi_{n-1}^{-1} ) \in D`
      for all `c \in C`.

These are group actions. If one of these group elements sends
the linear code `C` to itself, then we will call it an automorphism.
Depending on the group action we will call those groups:

    - permuation automorphism group

    - monomial automorphism group (every linear Hamming isometry is a monomial
      transformation of the ambient space, for `n\geq 3`)

    - automorphism group (every semilinear Hamming isometry is a semimonomial
      transformation of the ambient space, for `n\geq 3`)

This file contains

#. LinearCode class definition; LinearCodeFromVectorspace conversion function,

#. The spectrum (weight distribution), covering_radius, minimum distance
   programs (calling Steve Linton's or CJ Tjhal's C programs),
   characteristic_function, and several implementations of the Duursma zeta
   function (sd_zeta_polynomial, zeta_polynomial, zeta_function,
   chinen_polynomial, for example),

#. interface with best_known_linear_code_www (interface with codetables.de
   since A. Brouwer's online tables have been disabled),
   bounds_minimum_distance which call tables in GUAVA (updated May 2006)
   created by Cen Tjhai instead of the online internet tables,

#. generator_matrix, generator_matrix_systematic, information_set, list, parity_check_matrix,
   decode, dual_code, extended_code, shortened, punctured, genus, binomial_moment,
   and divisor methods for LinearCode,

#. Boolean-valued functions such as "==", is_self_dual, is_self_orthogonal,
   is_subcode, is_permutation_automorphism, is_permutation_equivalent (which
   interfaces with Robert Miller's partition refinement code),

#. permutation methods: is_permutation_automorphism,
   permutation_automorphism_group, permuted_code, standard_form,
   module_composition_factors,

#. design-theoretic methods: assmus_mattson_designs (implementing
   Assmus-Mattson Theorem),

#. code constructions, such as HammingCode and ToricCode, are in a separate
   ``code_constructions.py`` module; in the separate ``guava.py`` module, you
   will find constructions, such as RandomLinearCodeGuava and
   BinaryReedMullerCode, wrapped from the corresponding GUAVA codes.

EXAMPLES::

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

#. More wrappers

#. GRS codes and special decoders.

#. `P^1` Goppa codes and group actions on `P^1` RR space codes.

REFERENCES:

- [HP] W. C. Huffman and V. Pless, Fundamentals of error-correcting codes,
  Cambridge Univ. Press, 2003.

- [Gu] GUAVA manual, http://www.gap-system.org/Packages/guava.html

AUTHORS:

- David Joyner (2005-11-22, 2006-12-03): initial version

- William Stein (2006-01-23): Inclusion in Sage

- David Joyner (2006-01-30, 2006-04): small fixes

- David Joyner (2006-07): added documentation, group-theoretical methods,
  ToricCode

- David Joyner (2006-08): hopeful latex fixes to documentation, added list and
  __iter__ methods to LinearCode and examples, added hamming_weight function,
  fixed random method to return a vector, TrivialCode, fixed subtle bug in
  dual_code, added galois_closure method, fixed mysterious bug in
  permutation_automorphism_group (GAP was over-using "G" somehow?)

- David Joyner (2006-08): hopeful latex fixes to documentation, added
  CyclicCode, best_known_linear_code, bounds_minimum_distance,
  assmus_mattson_designs (implementing Assmus-Mattson Theorem).

- David Joyner (2006-09): modified decode syntax, fixed bug in
  is_galois_closed, added LinearCode_from_vectorspace, extended_code,
  zeta_function

- Nick Alexander (2006-12-10): factor GUAVA code to guava.py

- David Joyner (2007-05): added methods punctured, shortened, divisor,
  characteristic_polynomial, binomial_moment, support for
  LinearCode. Completely rewritten zeta_function (old version is now
  zeta_function2) and a new function, LinearCodeFromVectorSpace.

- David Joyner (2007-11): added zeta_polynomial, weight_enumerator,
  chinen_polynomial; improved best_known_code; made some pythonic revisions;
  added is_equivalent (for binary codes)

- David Joyner (2008-01): fixed bug in decode reported by Harald Schilly,
  (with Mike Hansen) added some doctests.

- David Joyner (2008-02): translated standard_form, dual_code to Python.

- David Joyner (2008-03): translated punctured, shortened, extended_code,
  random (and renamed random to random_element), deleted zeta_function2,
  zeta_function3, added wrapper automorphism_group_binary_code to Robert
  Miller's code), added direct_sum_code, is_subcode, is_self_dual,
  is_self_orthogonal, redundancy_matrix, did some alphabetical reorganizing
  to make the file more readable. Fixed a bug in permutation_automorphism_group
  which caused it to crash.

- David Joyner (2008-03): fixed bugs in spectrum and zeta_polynomial, which
  misbehaved over non-prime base rings.

- David Joyner (2008-10): use CJ Tjhal's MinimumWeight if char = 2 or 3 for
  min_dist; add is_permutation_equivalent and improve
  permutation_automorphism_group using an interface with Robert Miller's code;
  added interface with Leon's code for the spectrum method.

- David Joyner (2009-02): added native decoding methods (see module_decoder.py)

- David Joyner (2009-05): removed dependence on Guava, allowing it to be an
  option. Fixed errors in some docstrings.

- Kwankyu Lee (2010-01): added methods generator_matrix_systematic, information_set, and
  magma interface for linear codes.

- Niles Johnson (2010-08): :trac:`#3893`: ``random_element()`` should pass on ``*args`` and ``**kwds``.

- Thomas Feulner (2012-11): :trac:`13723`: deprecation of ``hamming_weight()``

- Thomas Feulner (2013-10): added methods to compute a canonical representative
  and the automorphism group

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C  = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""

#******************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#******************************************************************************
# python3
from __future__ import division

import sage.modules.module as module
from sage.categories.modules import Modules
from sage.categories.fields import Fields
from copy import copy
from sage.interfaces.all import gap
from sage.rings.finite_rings.finite_field_constructor import FiniteField as GF
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.arith.all import GCD, rising_factorial, binomial
from sage.groups.all import SymmetricGroup
from sage.misc.all import prod
from sage.misc.functional import log, is_even
from sage.rings.rational_field import QQ
from sage.structure.parent import Parent
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField
from sage.rings.integer_ring import IntegerRing
from sage.rings.integer import Integer
from sage.combinat.set_partition import SetPartitions
from sage.modules.free_module import VectorSpace
from sage.misc.randstate import current_randstate
from sage.misc.decorators import rename_keyword
from sage.misc.cachefunc import cached_method
from sage.misc.superseded import deprecated_function_alias
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.combinat.subset import Subsets
from sage.categories.cartesian_product import cartesian_product
# import compatible with py2 and py3
from six.moves.urllib.request import urlopen

ZZ = IntegerRing()

####################### coding theory functions ###############################

def code2leon(C):
    r"""
    Writes a file in Sage's temp directory representing the code C, returning
    the absolute path to the file. This is the Sage translation of the
    GuavaToLeon command in Guava's codefun.gi file.

    INPUT:

    - ``C`` - a linear code (over GF(p), p < 11)

    OUTPUT:

    - Absolute path to the file written

    EXAMPLES::

        sage: C = codes.HammingCode(GF(2), 3); C
        [7, 4] Hamming Code over Finite Field of size 2
        sage: file_loc = sage.coding.linear_code.code2leon(C)
        sage: f = open(file_loc); print f.read()
        LIBRARY code;
        code=seq(2,4,7,seq(
        1,0,0,0,0,1,1,
        0,1,0,0,1,0,1,
        0,0,1,0,1,1,0,
        0,0,0,1,1,1,1
        ));
        FINISH;
        sage: f.close()

    """
    from sage.misc.temporary_file import tmp_filename
    F = C.base_ring()
    p = F.order()  # must be prime and <11
    s = "LIBRARY code;\n"+"code=seq(%s,%s,%s,seq(\n"%(p,C.dimension(),C.length())
    Gr = [str(r)[1:-1].replace(" ","") for r in C.generator_matrix().rows()]
    s += ",\n".join(Gr) + "\n));\nFINISH;"
    file_loc = tmp_filename()
    f = open(file_loc,"w")
    f.write(s)
    f.close()
    return file_loc

def wtdist_gap(Gmat, n, F):
    r"""
    INPUT:

    -  ``Gmat`` - String representing a GAP generator matrix G of a linear code

    -  ``n`` - Integer greater than 1, representing the number of columns of G
       (i.e., the length of the linear code)

    -  ``F`` - Finite field (in Sage), base field the code

    OUTPUT:

    -  Spectrum of the associated code

    EXAMPLES::

        sage: Gstr = 'Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]'
        sage: F = GF(2)
        sage: sage.coding.linear_code.wtdist_gap(Gstr, 7, F)
        [1, 0, 0, 7, 7, 0, 0, 1]

    Here ``Gstr`` is a generator matrix of the Hamming [7,4,3] binary code.

    ALGORITHM:

    Uses C programs written by Steve Linton in the kernel of GAP, so is fairly
    fast.

    AUTHORS:

    - David Joyner (2005-11)
    """
    G = gap(Gmat)
    q = F.order()
    k = gap(F)
    z = 'Z(%s)*%s'%(q, [0]*n)     # GAP zero vector as a string
    _ = gap.eval("w:=DistancesDistributionMatFFEVecFFE("+Gmat+", GF("+str(q)+"),"+z+")")
    v = [eval(gap.eval("w["+str(i)+"]")) for i in range(1,n+2)] # because GAP returns vectors in compressed form
    return v

def min_wt_vec_gap(Gmat, n, k, F, algorithm=None):
    r"""
    Returns a minimum weight vector of the code generated by ``Gmat``.

    Uses C programs written by Steve Linton in the kernel of GAP, so is fairly
    fast. The option ``algorithm="guava"`` requires Guava. The default algorithm
    requires GAP but not Guava.

    INPUT:

    -  ``Gmat`` - String representing a GAP generator matrix G of a linear code
    -  n - Length of the code generated by G
    -  k - Dimension of the code generated by G
    -  F - Base field

    OUTPUT:

    -  Minimum weight vector of the code generated by ``Gmat``

    REMARKS:

    - The code in the default case allows one (for free) to also compute the
      message vector `m` such that `m\*G = v`, and the (minimum) distance, as
      a triple.  however, this output is not implemented.
    - The binary case can presumably be done much faster using Robert Miller's
      code (see the docstring for the spectrum method). This is also not (yet)
      implemented.

    EXAMPLES::

        sage: Gstr = "Z(2)*[[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]]"
        sage: sage.coding.linear_code.min_wt_vec_gap(Gstr,7,4,GF(2))
        (0, 1, 0, 1, 0, 1, 0)

    This output is different but still a minimum weight vector::

        sage: sage.coding.linear_code.min_wt_vec_gap(Gstr,7,4,GF(2),algorithm="guava")    # optional - gap_packages (Guava package)
        (0, 0, 1, 0, 1, 1, 0)

    Here ``Gstr`` is a generator matrix of the Hamming [7,4,3] binary code.

    TESTS:

    We check that :trac:`18480` is fixed::

        sage: codes.HammingCode(GF(2), 2).minimum_distance()
        3

    AUTHORS:

    - David Joyner (11-2005)
    """
    current_randstate().set_seed_gap()

    if algorithm=="guava":
        gap.LoadPackage('"guava"')
        from sage.interfaces.gap import gfq_gap_to_sage
        gap.eval("G:="+Gmat)
        C = gap(Gmat).GeneratorMatCode(F)
        cg = C.MinimumDistanceCodeword()
        c = [gfq_gap_to_sage(cg[j],F) for j in range(1,n+1)]
        V = VectorSpace(F,n)
        return V(c)

    q = F.order()
    ans = None
    dist_min = n + 1
    gap.eval('Gmat:='+Gmat)
    gap.eval('K:=GF({})'.format(q))
    gap.eval('v:=Z({})*{}'.format(q,[0]*n))
    for i in range(1,k+1):
        gap.eval("P:=AClosestVectorCombinationsMatFFEVecFFECoords(Gmat,K,v,{},1)".format(i))
        gap.eval("d:=WeightVecFFE(P[1])")
        v = gap("P[1]")
        dist = gap("d")
        if dist and dist < dist_min:
            dist_min = dist
            ans = list(v)

    if ans is None:
        raise RuntimeError("there is a bug here!")

    # return the result as a vector (and not a 1xn matrix)
    return vector(F, ans)

def best_known_linear_code(n, k, F):
    r"""
    Returns the best known (as of 11 May 2006) linear code of length ``n``,
    dimension ``k`` over field ``F``.  The function uses the tables described
    in ``bounds_minimum_distance`` to construct this code.

    This does not require an internet connection.

    EXAMPLES::

        sage: best_known_linear_code(10,5,GF(2))    # long time; optional - gap_packages (Guava package)
        Linear code of length 10, dimension 5 over Finite Field of size 2
        sage: gap.eval("C:=BestKnownLinearCode(10,5,GF(2))")     # long time; optional - gap_packages (Guava package)
        'a linear [10,5,4]2..4 shortened code'

    This means that best possible binary linear code of length 10 and
    dimension 5 is a code with minimum distance 4 and covering radius
    somewhere between 2 and 4.
    Use ``bounds_minimum_distance(10,5,GF(2))`` for further details.
    """
    q = F.order()
    C = gap("BestKnownLinearCode(%s,%s,GF(%s))"%(n,k,q))
    G = C.GeneratorMat()
    k = G.Length()
    n = G[1].Length()
    Gs = G._matrix_(F)
    MS = MatrixSpace(F,k,n)
    return LinearCode(MS(Gs))

def best_known_linear_code_www(n, k, F, verbose=False):
    r"""
    Explains the construction of the best known linear code over GF(q) with
    length n and dimension k, courtesy of the www page
    http://www.codetables.de/.

    INPUT:

    -  ``n`` - Integer, the length of the code

    -  ``k`` - Integer, the dimension of the code

    -  ``F`` - Finite field, of order 2, 3, 4, 5, 7, 8, or 9

    -  ``verbose`` - Bool (default: ``False``)

    OUTPUT:


    -  Text about why the bounds are as given

    EXAMPLES::

        sage: L = best_known_linear_code_www(72, 36, GF(2)) # optional - internet
        sage: print L                                       # optional - internet
        Construction of a linear code
        [72,36,15] over GF(2):
        [1]:  [73, 36, 16] Cyclic Linear Code over GF(2)
             CyclicCode of length 73 with generating polynomial x^37 + x^36 + x^34 +
        x^33 + x^32 + x^27 + x^25 + x^24 + x^22 + x^21 + x^19 + x^18 + x^15 + x^11 +
        x^10 + x^8 + x^7 + x^5 + x^3 + 1
        [2]:  [72, 36, 15] Linear Code over GF(2)
             Puncturing of [1] at 1
        <BLANKLINE>
        last modified: 2002-03-20

    This function raises an ``IOError`` if an error occurs downloading data or
    parsing it. It raises a ``ValueError`` if the ``q`` input is invalid.

    AUTHORS:

    - Steven Sivek (2005-11-14)
    - David Joyner (2008-03)
    """
    q = F.order()
    if not q in [2, 3, 4, 5, 7, 8, 9]:
        raise ValueError("q (=%s) must be in [2,3,4,5,7,8,9]"%q)
    n = int(n)
    k = int(k)

    param = ("?q=%s&n=%s&k=%s"%(q,n,k)).replace('L','')

    url = "http://iaks-www.ira.uka.de/home/grassl/codetables/BKLC/BKLC.php"+param
    if verbose:
        print "Looking up the bounds at %s"%url
    f = urlopen(url)
    s = f.read()
    f.close()

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError("Error parsing data (missing pre tags).")
    text = s[i+5:j].strip()
    return text

def bounds_minimum_distance(n, k, F):
    r"""
    Calculates a lower and upper bound for the minimum distance of an optimal
    linear code with word length ``n`` and dimension ``k`` over the field
    ``F``.

    The function returns a record with the two bounds and an explanation for
    each bound. The function Display can be used to show the explanations.

    The values for the lower and upper bound are obtained from a table
    constructed by Cen Tjhai for GUAVA, derived from the table of
    Brouwer. See http://www.codetables.de/ for the most recent data.
    These tables contain lower and upper bounds for `q=2` (when ``n <= 257``),
    `q=3` (when ``n <= 243``), `q=4` (``n <= 256``). (Current as of
    11 May 2006.) For codes over other fields and for larger word lengths,
    trivial bounds are used.

    This does not require an internet connection. The format of the output is
    a little non-intuitive. Try ``bounds_minimum_distance(10,5,GF(2))`` for
    an example.

    This function requires optional GAP package (Guava).

    EXAMPLES::

        sage: print bounds_minimum_distance(10,5,GF(2)) # optional - gap_packages (Guava package)
        rec(
          construction :=
           [ <Operation "ShortenedCode">,
              [
                  [ <Operation "UUVCode">,
                      [
                          [ <Operation "DualCode">,
                              [ [ <Operation "RepetitionCode">, [ 8, 2 ] ] ] ],
                          [ <Operation "UUVCode">,
                              [
                                  [ <Operation "DualCode">,
                                      [ [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ]
                                    , [ <Operation "RepetitionCode">, [ 4, 2 ] ] ] ]
                         ] ], [ 1, 2, 3, 4, 5, 6 ] ] ],
          k := 5,
          lowerBound := 4,
          lowerBoundExplanation := ...
          n := 10,
          q := 2,
          references := rec(
               ),
          upperBound := 4,
          upperBoundExplanation := ... )
    """
    q = F.order()
    gap.eval("data := BoundsMinimumDistance(%s,%s,GF(%s))"%(n,k,q))
    Ldata = gap.eval("Display(data)")
    return Ldata

def self_orthogonal_binary_codes(n, k, b=2, parent=None, BC=None, equal=False,
    in_test=None):
    """
    Returns a Python iterator which generates a complete set of
    representatives of all permutation equivalence classes of
    self-orthogonal binary linear codes of length in ``[1..n]`` and
    dimension in ``[1..k]``.

    INPUT:

    -  ``n`` - Integer, maximal length

    -  ``k`` - Integer, maximal dimension

    -  ``b`` - Integer, requires that the generators all have weight divisible
       by ``b`` (if ``b=2``, all self-orthogonal codes are generated, and if
       ``b=4``, all doubly even codes are generated). Must be an even positive
       integer.

    -  ``parent`` - Used in recursion (default: ``None``)

    -  ``BC`` - Used in recursion (default: ``None``)

    -  ``equal`` - If ``True`` generates only [n, k] codes (default: ``False``)

    -  ``in_test`` - Used in recursion (default: ``None``)

    EXAMPLES:

    Generate all self-orthogonal codes of length up to 7 and dimension up
    to 3::

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

    Generate all doubly-even codes of length up to 7 and dimension up
    to 3::

        sage: for B in self_orthogonal_binary_codes(7,3,4):
        ...    print B; print B.generator_matrix()
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

    Generate all doubly-even codes of length up to 7 and dimension up
    to 2::

        sage: for B in self_orthogonal_binary_codes(7,2,4):
        ...    print B; print B.generator_matrix()
        Linear code of length 4, dimension 1 over Finite Field of size 2
        [1 1 1 1]
        Linear code of length 6, dimension 2 over Finite Field of size 2
        [1 1 1 1 0 0]
        [0 1 0 1 1 1]

    Generate all self-orthogonal codes of length equal to 8 and
    dimension equal to 4::

        sage: for B in self_orthogonal_binary_codes(8, 4, equal=True):
        ...     print B; print B.generator_matrix()
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

    Since all the codes will be self-orthogonal, b must be divisible by
    2::

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

class AbstractLinearCode(module.Module):
    """
    Abstract class for linear codes.

    This class contains all methods that can be used on Linear Codes
    and on Linear Codes families.
    So, every Linear Code-related class should inherit from this abstract
    class.

    To implement a linear code, you need to:

    - inherit from AbstractLinearCode

    - call AbstractLinearCode ``__init__`` method in the subclass constructor. Example:
      ``super(SubclassName, self).__init__(base_field, length, "EncoderName", "DecoderName")``.
      By doing that, your subclass will have its ``length`` parameter
      initialized and will be properly set as a member of the category framework.
      You need of course to complete the constructor by adding any additional parameter
      needed to describe properly the code defined in the subclass.

    - fill the dictionary of its encoders in ``sage.coding.__init__.py`` file. Example:
      I want to link the encoder ``MyEncoderClass`` to ``MyNewCodeClass``
      under the name ``MyEncoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_encoders["NameOfMyEncoder"] = MyEncoderClass`` and all instances of
      ``MyNewCodeClass`` will be able to use instances of ``MyEncoderClass``.

    - fill the dictionary of its decoders in ``sage.coding.__init__`` file. Example:
      I want to link the encoder ``MyDecoderClass`` to ``MyNewCodeClass``
      under the name ``MyDecoderName``.
      All I need to do is to write this line in the ``__init__.py`` file:
      ``MyNewCodeClass._registered_decoders["NameOfMyDecoder"] = MyDecoderClass`` and all instances of
      ``MyNewCodeClass`` will be able to use instances of ``MyDecoderClass``.


    As AbstractLinearCode is not designed to be implemented, it does not have any representation
    methods. You should implement ``_repr_`` and ``_latex_`` methods in the subclass.

    .. NOTE::

        :class:`AbstractLinearCode` has generic implementations of the comparison methods ``__cmp__``
        and ``__eq__`` which use the generator matrix and are quite slow. In subclasses you are
        encouraged to override these functions.

    .. WARNING::

        The default encoder should always have `F^{k}` as message space, with `k` the dimension
        of the code and `F` is the base ring of the code.

        A lot of methods of the abstract class rely on the knowledge of a generator matrix.
        It is thus strongly recommended to set an encoder with a generator matrix implemented
        as a default encoder.
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, default_encoder_name, default_decoder_name):
        """
        Initializes mandatory parameters for a Linear Code object.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear code. An abstract
        linear code object should never be created.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``length`` -- the length of ``self``

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        EXAMPLES:

        We first create a new LinearCode subclass::

            sage: class CodeExample(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "GeneratorMatrix", "Syndrome")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "Dummy code of length %d, dimension %d over %s" % (self.length(), self.dimension(), self.base_field())

        We now create a member of our newly made class::

            sage: generator_matrix = matrix(GF(17), 5, 10,
            ....:                           {(i,i):1 for i in range(5)})
            sage: C = CodeExample(GF(17), 10, 5, generator_matrix)

        We can check its existence and parameters::

            sage: C
            Dummy code of length 10, dimension 5 over Finite Field of size 17

        We can check that it is truly a part of the framework category::

            sage: C.parent()
            <class '__main__.CodeExample_with_category'>
            sage: C.category()
            Category of facade finite dimensional vector spaces with basis over Finite Field of size 17

        And any method that works on linear codes works for our new dummy code::

            sage: C.minimum_distance()
            1
            sage: C.is_self_orthogonal()
            False
            sage: print C.divisor() #long time
            1

        TESTS:

        If the length field is neither a Python int nor a Sage Integer, it will
        raise a exception::

            sage: C = CodeExample(GF(17), 10.0, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: length must be a Python int or a Sage Integer

        If the name of the default decoder is not known by the class, it will raise
        a exception::

            sage: class CodeExample(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "GeneratorMatrix", "Fail")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "Dummy code of length %d, dimension %d over %s" % (self.length(), self.dimension(), self.base_field())

            sage: C = CodeExample(GF(17), 10, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders

        If the name of the default encoder is not known by the class, it will raise
        an exception::

            sage: class CodeExample(sage.coding.linear_code.AbstractLinearCode):
            ....:   def __init__(self, field, length, dimension, generator_matrix):
            ....:       sage.coding.linear_code.AbstractLinearCode.__init__(self,field, length, "Fail", "Syndrome")
            ....:       self._dimension = dimension
            ....:       self._generator_matrix = generator_matrix
            ....:   def generator_matrix(self):
            ....:       return self._generator_matrix
            ....:   def _repr_(self):
            ....:       return "Dummy code of length %d, dimension %d over %s" % (self.length(), self.dimension(), self.base_field())

            sage: C = CodeExample(GF(17), 10, 5, generator_matrix)
            Traceback (most recent call last):
            ...
            ValueError: You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders

        A ring instead of a field::

            sage: codes.LinearCode(IntegerModRing(4),matrix.ones(4))
            Traceback (most recent call last):
            ...
            ValueError: 'generator' must be defined on a field (not a ring)
        """
        if not isinstance(length, (int, Integer)):
            raise ValueError("length must be a Python int or a Sage Integer")
        if not base_field.is_field():
            raise ValueError("'base_field' must be a field (and {} is not one)".format(base_field))
        if not default_encoder_name in self._registered_encoders:
            raise ValueError("You must set a valid encoder as default encoder for this code, by filling in the dictionary of registered encoders")
        if not default_decoder_name in self._registered_decoders:
            raise ValueError("You must set a valid decoder as default decoder for this code, by filling in the dictionary of registered decoders")
        self._length = Integer(length)
        self._default_decoder_name = default_decoder_name
        self._default_encoder_name = default_encoder_name
        cat = Modules(base_field).FiniteDimensional().WithBasis().Finite()
        facade_for = VectorSpace(base_field, self._length)
        self.Element = type(facade_for.an_element()) #for when we made this a non-facade parent
        Parent.__init__(self, base=base_field, facade=facade_for, category=cat)

    def _latex_(self):
        """
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: latex(C)
            [7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "[%s, %s]\\textnormal{ Linear code over }%s"\
                % (self.length(), self.dimension(), self.base_ring()._latex_())

    def _an_element_(self):
        r"""
        Return an element of the linear code. Currently, it simply returns
        the first row of the generator matrix.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.an_element()
            (1, 0, 0, 0, 0, 1, 1)
            sage: C2 = C.cartesian_product(C)
            sage: C2.an_element()
            ((1, 0, 0, 0, 0, 1, 1), (1, 0, 0, 0, 0, 1, 1))
        """
        return self.gens()[0]

    def add_decoder(self, name, decoder):
        r"""
        Adds an decoder to the list of registered decoders of ``self``.

        .. NOTE::

            This method only adds ``decoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.decoder.Decoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the decoder

        - ``decoder`` -- the class name of the decoder

        EXAMPLES:

        First of all, we create a (very basic) new decoder::

            sage: class MyDecoder(sage.coding.decoder.Decoder):
            ....:   def __init__(self, code):
            ....:       super(MyDecoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyDecoder decoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new decoder to the list of available decoders of C::

            sage: C.add_decoder("MyDecoder", MyDecoder)
            sage: C.decoders_available()
            ['MyDecoder', 'Syndrome', 'NearestNeighbor']

        We can verify that any new code will not know MyDecoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: C2.decoders_available()
            ['Syndrome', 'NearestNeighbor']

        TESTS:

        It is impossible to use a name which is in the dictionary of available decoders::

            sage: C.add_decoder("Syndrome", MyDecoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered decoder with this name
        """
        if self._registered_decoders == self.__class__._registered_decoders:
            self._registered_decoders = copy(self._registered_decoders)
            reg_dec = self._registered_decoders
            if name in reg_dec:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder
        else:
            if name in self._registered_decoders:
                raise ValueError("There is already a registered decoder with this name")
            reg_dec[name] = decoder

    def add_encoder(self, name, encoder):
        r"""
        Adds an encoder to the list of registered encoders of ``self``.

        .. NOTE::

            This method only adds ``encoder`` to ``self``, and not to any member of the class
            of ``self``. To know how to add an :class:`sage.coding.encoder.Encoder`, please refer
            to the documentation of :class:`AbstractLinearCode`.

        INPUT:

        - ``name`` -- the string name for the encoder

        - ``encoder`` -- the class name of the encoder

        EXAMPLES:

        First of all, we create a (very basic) new encoder::

            sage: class MyEncoder(sage.coding.encoder.Encoder):
            ....:   def __init__(self, code):
            ....:       super(MyEncoder, self).__init__(code)
            ....:   def _repr_(self):
            ....:       return "MyEncoder encoder with associated code %s" % self.code()

        We now create a new code::

            sage: C = codes.HammingCode(GF(2), 3)

        We can add our new encoder to the list of available encoders of C::

            sage: C.add_encoder("MyEncoder", MyEncoder)
            sage: C.encoders_available()
            ['MyEncoder', 'ParityCheck']

        We can verify that any new code will not know MyEncoder::

            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: C2.encoders_available()
            ['ParityCheck']

        TESTS:

        It is impossible to use a name which is in the dictionary of available encoders::

            sage: C.add_encoder("ParityCheck", MyEncoder)
            Traceback (most recent call last):
            ...
            ValueError: There is already a registered encoder with this name
        """
        if self._registered_encoders == self.__class__._registered_encoders:
            self._registered_encoders = copy(self._registered_encoders)
            reg_enc = self._registered_encoders
            if name in reg_enc:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder
        else:
            if name in self._registered_encoders:
                raise ValueError("There is already a registered encoder with this name")
            reg_enc[name] = encoder

    def automorphism_group_gens(self, equivalence="semilinear"):
        r"""
        Return generators of the automorphism group of ``self``.

        INPUT:

        - ``equivalence`` (optional) -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        OUTPUT:

        - generators of the automorphism group of ``self``
        - the order of the automorphism group of ``self``

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4, 'z'), 3)
            sage: C.automorphism_group_gens()
            ([((1, 1, 1, z, z + 1, z + 1, z + 1, z, z, 1, 1, 1, z, z, z + 1, z, z, z + 1, z + 1, z + 1, 1); (1,6,12,17)(2,16,4,5,11,8,14,13)(3,21,19,10,20,18,15,9), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z + 1), ((1, 1, 1, z, z + 1, 1, 1, z, z, z + 1, z, z, z + 1, z + 1, z + 1, 1, z + 1, z, z, 1, 1); (1,6,9,13,15,18)(2,21)(3,16,7)(4,5,11,10,12,14)(17,19), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z); (), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z)], 362880)
            sage: C.automorphism_group_gens(equivalence="linear")
            ([((z, z, 1, 1, z + 1, z, z + 1, z, z, z + 1, 1, 1, 1, z + 1, z, z, z + 1, z + 1, 1, 1, z); (1,5,10,9,4,14,11,16,18,20,6,19,12,15,3,8,2,17,7,13,21), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((z + 1, 1, z, 1, 1, z + 1, z + 1, z, 1, z, z + 1, z, z + 1, z + 1, z, 1, 1, z + 1, z + 1, z + 1, z); (1,17,10)(2,15,13)(4,11,21)(5,18,12)(6,14,19)(7,8,16), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1, z + 1); (), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z)], 181440)
            sage: C.automorphism_group_gens(equivalence="permutational")
            ([((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (1,11)(3,10)(4,9)(5,7)(12,21)(14,20)(15,19)(16,17), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (2,18)(3,19)(4,10)(5,16)(8,13)(9,14)(11,21)(15,20), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (1,19)(3,17)(4,21)(5,20)(7,14)(9,12)(10,16)(11,15), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z), ((1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1); (2,13)(3,14)(4,20)(5,11)(8,18)(9,19)(10,15)(16,21), Ring endomorphism of Finite Field in z of size 2^2
                  Defn: z |--> z)], 64)
        """
        aut_group_can_label = self._canonize(equivalence)
        return aut_group_can_label.get_autom_gens(), \
               aut_group_can_label.get_autom_order()

    def ambient_space(self):
        r"""
        Returns the ambient vector space of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.ambient_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        return VectorSpace(self.base_ring(),self.length())

    def assmus_mattson_designs(self, t, mode=None):
        r"""
        Assmus and Mattson Theorem (section 8.4, page 303 of [HP]): Let
        `A_0, A_1, ..., A_n` be the weights of the codewords in a binary
        linear `[n , k, d]` code `C`, and let `A_0^*, A_1^*, ..., A_n^*` be
        the weights of the codewords in its dual `[n, n-k, d^*]` code `C^*`.
        Fix a `t`, `0<t<d`, and let

        .. math::

           s = |\{ i\ |\ A_i^* \not= 0, 0< i \leq n-t\}|.

        Assume `s\leq d-t`.

        1. If `A_i\not= 0` and `d\leq i\leq n`
           then `C_i = \{ c \in C\ |\ wt(c) = i\}` holds a simple t-design.

        2. If `A_i^*\not= 0` and `d*\leq i\leq n-t` then
           `C_i^* = \{ c \in C^*\ |\ wt(c) = i\}` holds a simple t-design.

        A block design is a pair `(X,B)`, where `X` is a non-empty finite set
        of `v>0` elements called points, and `B` is a non-empty finite
        multiset of size b whose elements are called blocks, such that each
        block is a non-empty finite multiset of `k` points. `A` design without
        repeated blocks is called a simple block design. If every subset of
        points of size `t` is contained in exactly `\lambda` blocks the block
        design is called a `t-(v,k,\lambda)` design (or simply a `t`-design
        when the parameters are not specified). When `\lambda=1` then the
        block design is called a `S(t,k,v)` Steiner system.

        In the Assmus and Mattson Theorem (1), `X` is the set `\{1,2,...,n\}`
        of coordinate locations and `B = \{supp(c)\ |\ c \in C_i\}` is the set
        of supports of the codewords of `C` of weight `i`. Therefore, the
        parameters of the `t`-design for `C_i` are

        ::

            t =       given
            v =       n
            k =       i   (k not to be confused with dim(C))
            b =       Ai
            lambda = b*binomial(k,t)/binomial(v,t) (by Theorem 8.1.6,
                                                       p 294, in [HP])

        Setting the ``mode="verbose"`` option prints out the values of the
        parameters.

        The first example below means that the binary [24,12,8]-code C has
        the property that the (support of the) codewords of weight 8 (resp.,
        12, 16) form a 5-design. Similarly for its dual code `C^*` (of course
        `C=C^*` in this case, so this info is extraneous). The test fails to
        produce 6-designs (ie, the hypotheses of the theorem fail to hold,
        not that the 6-designs definitely don't exist). The command
        assmus_mattson_designs(C,5,mode="verbose") returns the same value
        but prints out more detailed information.

        The second example below illustrates the blocks of the 5-(24, 8, 1)
        design (i.e., the S(5,8,24) Steiner system).

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()             #  example 1
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
            sage: blocks = [c.support() for c in C if c.hamming_weight()==8]; len(blocks)  # long time computation
            759

        REFERENCE:

        - [HP] W. C. Huffman and V. Pless, Fundamentals of ECC,
          Cambridge Univ. Press, 2003.
        """
        C = self
        ans = []
        G = C.generator_matrix()
        n = len(G.columns())
        Cp = C.dual_code()
        wts = C.spectrum()
        d = min([i for i in range(1,len(wts)) if wts[i]!=0])
        if t>=d:
            return 0
        nonzerowts = [i for i in range(len(wts)) if wts[i]!=0 and i<=n and i>=d]
        if mode=="verbose":
            for w in nonzerowts:
                print("The weight w={} codewords of C* form a t-(v,k,lambda) design, where\n \
                        t={}, v={}, k={}, lambda={}. \nThere are {} block of this design.".format(\
                        w,t,n,w,wts[w]*binomial(w,t)//binomial(n,t),wts[w]))
        wtsp = Cp.spectrum()
        dp = min([i for i in range(1,len(wtsp)) if wtsp[i]!=0])
        nonzerowtsp = [i for i in range(len(wtsp)) if wtsp[i]!=0 and i<=n-t and i>=dp]
        s = len([i for i in range(1,n) if wtsp[i]!=0 and i<=n-t and i>0])
        if mode=="verbose":
            for w in nonzerowtsp:
                print("The weight w={} codewords of C* form a t-(v,k,lambda) design, where\n \
                        t={}, v={}, k={}, lambda={}. \nThere are {} block of this design.".format(\
                        w,t,n,w,wts[w]*binomial(w,t)//binomial(n,t),wts[w]))
        if s<=d-t:
            des = [[t,(n,w,wts[w]*binomial(w,t)//binomial(n,t))] for w in nonzerowts]
            ans = ans + ["weights from C: ",nonzerowts,"designs from C: ",des]
            desp = [[t,(n,w,wtsp[w]*binomial(w,t)//binomial(n,t))] for w in nonzerowtsp]
            ans = ans + ["weights from C*: ",nonzerowtsp,"designs from C*: ",desp]
            return ans
        return 0

    def base_field(self):
        r"""
        Return the base field of ``self``.

        EXAMPLES::

            sage: G  = Matrix(GF(2), [[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C.base_field()
            Finite Field of size 2
        """
        return self.base_ring()

    def basis(self):
        r"""
        Returns a basis of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.basis()
            [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.gens()

    # S. Pancratz, 19 Jan 2010:  In the doctests below, I removed the example
    # ``C.binomial_moment(3)``, which was also marked as ``#long``.  This way,
    # we shorten the doctests time while still maintaining a zero and a
    # non-zero example.
    def binomial_moment(self, i):
        r"""
        Returns the i-th binomial moment of the `[n,k,d]_q`-code `C`:

        .. math::

            B_i(C) = \sum_{S, |S|=i} \frac{q^{k_S}-1}{q-1}

        where `k_S` is the dimension of the shortened code `C_{J-S}`,
        `J=[1,2,...,n]`. (The normalized binomial moment is
        `b_i(C) = \binom(n,d+i)^{-1}B_{d+i}(C)`.) In other words, `C_{J-S}`
        is isomorphic to the subcode of C of codewords supported on S.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.binomial_moment(2)
            0
            sage: C.binomial_moment(4)    # long time
            35

        .. warning::

            This is slow.

        REFERENCE:

        - I. Duursma, "Combinatorics of the two-variable zeta function",
          Finite fields and applications, 109-136, Lecture Notes in
          Comput. Sci., 2948, Springer, Berlin, 2004.
        """
        n = self.length()
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
            return binomial(n,i)*(q**(i+k-n) -1)//(q-1)
        P = SetPartitions(J,2).list()
        b = QQ(0)
        for p in P:
            p = list(p)
            S = p[0]
            if len(S)==n-i:
                C_S = self.shortened(S)
                k_S = C_S.dimension()
                b = b + (q**(k_S) -1)//(q-1)
        return b

    @cached_method
    def _canonize(self, equivalence):
        r"""
        Compute a canonical representative and the automorphism group
        under the action of the semimonomial transformation group.

        INPUT:

        - ``equivalence`` -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4, 'z'), 3)
            sage: aut_group_can_label = C._canonize("semilinear")
            sage: C_iso = LinearCode(aut_group_can_label.get_transporter()*C.generator_matrix())
            sage: C_iso == aut_group_can_label.get_canonical_form()
            True
            sage: aut_group_can_label.get_autom_gens()
            [((z, z + 1, 1, z + 1, z, z, z, z + 1, 1, z + 1, z + 1, z, z + 1, 1, z + 1, 1, z, z + 1, 1, z + 1, z + 1); (1,12,21,18,15,20)(2,19,16)(3,4,11,6,13,7)(5,8)(10,14,17), Ring endomorphism of Finite Field in z of size 2^2
            Defn: z |--> z + 1),
            ((z + 1, 1, z, 1, 1, z, 1, z + 1, z, z + 1, z, z + 1, z, 1, z, z + 1, z, z, z + 1, z + 1, 1); (1,20,2,9,13,21,11,17,10,16,3,5,18,8)(4,12,6,15,14,19,7), Ring endomorphism of Finite Field in z of size 2^2
            Defn: z |--> z + 1),
            ((z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z); (), Ring endomorphism of Finite Field in z of size 2^2
            Defn: z |--> z)]
        """
        from sage.coding.codecan.autgroup_can_label import LinearCodeAutGroupCanLabel
        return LinearCodeAutGroupCanLabel(self, algorithm_type=equivalence)

    def canonical_representative(self, equivalence="semilinear"):
        r"""
        Compute a canonical orbit representative under the action of the
        semimonomial transformation group.

        See :mod:`sage.coding.codecan.autgroup_can_label`
        for more details, for example if you would like to compute
        a canonical form under some more restrictive notion of equivalence,
        i.e. if you would like to restrict the permutation group
        to a Young subgroup.

        INPUT:

        - ``equivalence`` (optional) -- which defines the acting group, either

            * ``permutational``

            * ``linear``

            * ``semilinear``

        OUTPUT:

        - a canonical representative of ``self``
        - a semimonomial transformation mapping ``self`` onto its representative

        EXAMPLES::

            sage: F.<z> = GF(4)
            sage: C = codes.HammingCode(F, 3)
            sage: CanRep, transp = C.canonical_representative()

        Check that the transporter element is correct::

            sage: LinearCode(transp*C.generator_matrix()) == CanRep
            True

        Check if an equivalent code has the same canonical representative::

            sage: f = F.hom([z**2])
            sage: C_iso = LinearCode(C.generator_matrix().apply_map(f))
            sage: CanRep_iso, _ = C_iso.canonical_representative()
            sage: CanRep_iso == CanRep
            True

        Since applying the Frobenius automorphism could be extended to an
        automorphism of `C`, the following must also yield ``True``::

            sage: CanRep1, _ = C.canonical_representative("linear")
            sage: CanRep2, _ = C_iso.canonical_representative("linear")
            sage: CanRep2 == CanRep1
            True
        """
        aut_group_can_label = self._canonize(equivalence)
        return aut_group_can_label.get_canonical_form(), \
               aut_group_can_label.get_transporter()

    def __contains__(self, v):
        r"""
        Returns True if `v` can be coerced into `self`. Otherwise, returns False.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: vector((1, 0, 0, 0, 0, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 2, 1, 1)) in C   # indirect doctest
            True
            sage: vector((1, 0, 0, 0, 0, 1/2, 1)) in C # indirect doctest
            False
        """
        if not v in self.ambient_space() or len(v) != self.length():
            return False
        return self.syndrome(v) == 0

    def characteristic(self):
        r"""
        Returns the characteristic of the base ring of `self`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.characteristic()
            2
        """
        return (self.base_ring()).characteristic()

    def characteristic_polynomial(self):
        r"""
        Returns the characteristic polynomial of a linear code, as defined in
        van Lint's text [vL].

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.characteristic_polynomial()
            -4/3*x^3 + 64*x^2 - 2816/3*x + 4096

        REFERENCES:

        - van Lint, Introduction to coding theory, 3rd ed., Springer-Verlag
          GTM, 86, 1999.
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

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.chinen_polynomial()       # long time
            1/5*(2*sqrt(2)*t^3 + 2*sqrt(2)*t^2 + 2*t^2 + sqrt(2)*t + 2*t + 1)/(sqrt(2) + 1)
            sage: C = codes.TernaryGolayCode()
            sage: C.chinen_polynomial()       # long time
            1/7*(3*sqrt(3)*t^3 + 3*sqrt(3)*t^2 + 3*t^2 + sqrt(3)*t + 3*t + 1)/(sqrt(3) + 1)

        This last output agrees with the corresponding example given in
        Chinen's paper below.

        REFERENCES:

        - Chinen, K. "An abundance of invariant polynomials satisfying the
          Riemann hypothesis", April 2007 preprint.
        """
        from sage.functions.all import sqrt
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
            # Sage does not find dealing with sqrt(int) *as an algebraic object*
            # an easy thing to do. Some tricky gymnastics are used to
            # make Sage deal with objects over QQ(sqrt(q)) nicely.
            if is_even(n):
                Pd = q**(k-n//2) * RT(Cd.zeta_polynomial()) * T**(dperp - d)
            else:
                Pd = s * q**(k-(n+1)//2) * RT(Cd.zeta_polynomial()) * T**(dperp - d)
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

    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Cperp = C.dual_code()
            sage: C; Cperp
            [7, 4] Hamming Code over Finite Field of size 2
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C.generator_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.parity_check_matrix()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
            sage: Cperp.parity_check_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: Cperp.generator_matrix()
             [1 0 1 0 1 0 1]
             [0 1 1 0 0 1 1]
             [0 0 0 1 1 1 1]
        """
        G = self.generator_matrix()
        H = G.right_kernel()
        return H.basis_matrix()

    check_mat = deprecated_function_alias(17973, parity_check_matrix)

    @cached_method
    def covering_radius(self):
        r"""
        Wraps Guava's ``CoveringRadius`` command.

        The covering radius of a linear code `C` is the smallest number `r`
        with the property that each element `v` of the ambient vector space
        of `C` has at most a distance `r` to the code `C`. So for each
        vector `v` there must be an element `c` of `C` with `d(v,c) \leq  r`.
        A binary linear code with reasonable small covering radius is often
        referred to as a covering code.

        For example, if `C` is a perfect code, the covering radius is equal
        to `t`, the number of errors the code can correct, where `d = 2t+1`,
        with `d` the minimum distance of `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 5)
            sage: C.covering_radius()  # optional - gap_packages (Guava package)
            1
        """
        F = self.base_ring()
        G = self.generator_matrix()
        gapG = gap(G)
        C = gapG.GeneratorMatCode(gap(F))
        r = C.CoveringRadius()
        try:
            return ZZ(r)
        except TypeError:
            raise RuntimeError("the covering radius of this code cannot be computed by Guava")

    def decode(self, right, algorithm="syndrome"):
        r"""
        Corrects the errors in ``right`` and returns a codeword.

        INPUT:

        - ``right`` -- a vector of the same length as ``self`` over
          the base field of ``self``

        - ``algorithm`` -- (default: ``'syndrome'``) Name of the decoding algorithm which
          will be used to decode ``right``. Can be ``'syndrome'`` or ``'nearest_neighbor'``.

        .. NOTE::

            This is a deprecated method which will soon be removed from Sage. Please use
            :func:`decode_to_code` instead.
        """
        from sage.misc.superseded import deprecation
        deprecation(18813, "decode is deprecated and will be removed from sage soon. Please call decode_to_code instead")
        if algorithm == "syndrome":
            return self.decode_to_code(right, decoder_name="Syndrome")
        elif algorithm == "nearest neighbor":
            return self.decode_to_code(right, decoder_name="NearestNeighbor")
        else:
            return self.decode_to_code(right, decoder_name=algorithm)

    def decode_to_code(self, word, decoder_name=None, **kwargs):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over
          the base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: C.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: C.decoders_available()
            ['Syndrome', 'NearestNeighbor']
            sage: C.decode_to_code(w_err, 'NearestNeighbor')
            (1, 1, 0, 0, 1, 1, 0)
        """
        D = self.decoder(decoder_name, **kwargs)
        return D.decode_to_code(word)

    def decode_to_message(self, word, decoder_name=None, **kwargs):
        r"""
        Correct the errors in word and decodes it to the message space.

        INPUT:

        - ``word`` -- a vector of the same length as ``self`` over the
          base field of ``self``

        - ``decoder_name`` -- (default: ``None``) Name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to :meth:`decoder`

        OUTPUT:

        - A vector of the message space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.decode_to_message(word)
            (0, 1, 1, 0)

        It is possible to manually choose the decoder amongst the list of the available ones::

            sage: C.decoders_available()
            ['Syndrome', 'NearestNeighbor']
            sage: C.decode_to_message(word, 'NearestNeighbor')
            (0, 1, 1, 0)
        """
        return self.unencode(self.decode_to_code(word, decoder_name, **kwargs), **kwargs)

    @cached_method
    def decoder(self, decoder_name=None, **kwargs):
        r"""
        Return a decoder of ``self``.

        INPUT:

        - ``decoder_name`` -- (default: ``None``) name of the decoder which will be
          returned. The default decoder of ``self`` will be used if
          default value is kept.

        - ``kwargs`` -- all additional arguments will be forwarded to the constructor of the decoder
          that will be returned by this method

        OUTPUT:

        - a decoder object

        Besides creating the decoder and returning it, this method also stores
        the decoder in a cache. With this behaviour, each decoder will be created
        at most one time for ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoder()
            Syndrome decoder for Linear code of length 7, dimension 4 over Finite Field of size 2 handling errors of weight up to 1


        If the name of a decoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: C.decoders_available()
            ['Syndrome', 'NearestNeighbor']
            sage: C.decoder('Try')
            Traceback (most recent call last):
            ...
            ValueError: Passed Decoder name not known
        """
        if decoder_name is None:
            decoder_name = self._default_decoder_name
        if decoder_name in self._registered_decoders:
            decClass = self._registered_decoders[decoder_name]
            D = decClass(self, **kwargs)
            return D
        else:
            raise ValueError("Passed Decoder name not known")

    def decoders_available(self, classes=False):
        r"""
        Returns a list of the available decoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``, it also
          returns the decoders' classes associated with the decoders' names.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.decoders_available()
            ['Syndrome', 'NearestNeighbor']

            sage: C.decoders_available(True)
            {'NearestNeighbor': <class 'sage.coding.linear_code.LinearCodeNearestNeighborDecoder'>,
             'Syndrome': <class 'sage.coding.linear_code.LinearCodeSyndromeDecoder'>}
        """
        if classes == True:
            return copy(self._registered_decoders)
        return self._registered_decoders.keys()

    def divisor(self):
        r"""
        Returns the divisor of a code, which is the smallest integer `d_0 > 0`
        such that each `A_i > 0` iff `i` is divisible by `d_0`.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.divisor()   # Type II self-dual
            4
            sage: C = codes.QuadraticResidueCodeEvenPair(17,GF(2))[0]
            sage: C.divisor()
            2
        """
        C = self
        A = C.spectrum()
        n = C.length()
        V = VectorSpace(QQ,n+1)
        S = V(A).nonzero_positions()
        S0 = [S[i] for i in range(1,len(S))]
        if len(S)>1: return GCD(S0)
        return 1

    def is_projective(self):
        r"""
        Test  whether the code is projective.

        A linear code `C` over a field is called *projective* when its dual `Cd`
        has minimum weight `\geq 3`, i.e. when no two coordinate positions of
        `C` are linearly independent (cf. definition 3 from [BS11] or 9.8.1 from
        [BH12]).

        EXAMPLE::

            sage: C = codes.BinaryGolayCode()
            sage: C.is_projective()
            True
            sage: C.dual_code().minimum_distance()
            8

        A non-projective code::

            sage: C = codes.LinearCode(matrix(GF(2),[[1,0,1],[1,1,1]]))
            sage: C.is_projective()
            False

        REFERENCE:

        .. [BS11] \E. Byrne and A. Sneyd,
           On the Parameters of Codes with Two Homogeneous Weights.
           WCC 2011-Workshop on coding and cryptography, pp. 81-90. 2011.
           https://hal.inria.fr/inria-00607341/document
        """
        M = self.generator_matrix().transpose()
        R = self.base_field()

        def projectivize(row):
            if not row.is_zero():
                for i in range(len(row)):
                    if row[i]:
                        break
                row = ~(row[i]) * row
            row.set_immutable()
            return row

        rows = set()
        for row in M.rows():
            row = projectivize(row)
            if row in rows:
                return False
            rows.add(row)

        return True

    def dual_code(self):
        r"""
        Returns the dual code `C^{\perp}` of the code `C`,

        .. math::

            C^{\perp} = \{ v \in V\ |\ v\cdot c = 0,\ \forall c \in C \}.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.dual_code()
            Linear code of length 7, dimension 3 over Finite Field of size 2
            sage: C = codes.HammingCode(GF(4, 'a'), 3)
            sage: C.dual_code()
            Linear code of length 21, dimension 3 over Finite Field in a of size 2^2
        """
        return LinearCode(self.parity_check_matrix())

    def dimension(self):
        r"""
        Returns the dimension of this code.

        EXAMPLES::

            sage: G = matrix(GF(2),[[1,0,0],[1,1,0]])
            sage: C = LinearCode(G)
            sage: C.dimension()
            2
        """
        return self._dimension

    def direct_sum(self, other):
        """
        Returns the code given by the direct sum of the codes ``self`` and
        ``other``, which must be linear codes defined over the same base ring.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = C1.direct_sum(C1); C2
            Linear code of length 14, dimension 8 over Finite Field of size 2
            sage: C3 = C1.direct_sum(C2); C3
            Linear code of length 21, dimension 12 over Finite Field of size 2
        """
        C1 = self; C2 = other
        G1 = C1.generator_matrix()
        G2 = C2.generator_matrix()
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
        Checks if ``self`` is equal to ``right``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = codes.HammingCode(GF(2), 3)
            sage: C1 == C2
            True

        TESTS:

        We check that :trac:`16644` is fixed::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C == ZZ
            False
        """
        if not (isinstance(right, LinearCode)\
                and self.length() == right.length()\
                and self.dimension() == right.dimension()\
                and self.base_ring() == right.base_ring()):
            return False
        Ks = self.parity_check_matrix().right_kernel()
        rbas = right.gens()
        if not all(c in Ks for c in rbas):
            return False
        Kr = right.parity_check_matrix().right_kernel()
        sbas = self.gens()
        if not all(c in Kr for c in sbas):
            return False
        return True

    def __ne__(self, other):
        r"""
        Tests inequality of ``self`` and ``other``.

        This is a generic implementation, which returns the inverse of ``__eq__`` for self.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C1 = LinearCode(G)
            sage: C2 = LinearCode(G)
            sage: C1 != C2
            False
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,1,1]])
            sage: C2 = LinearCode(G)
            sage: C1 != C2
            True
        """
        return not self == other

    def encode(self, word, encoder_name=None, **kwargs):
        r"""
        Transforms an element of a message space into a codeword.

        INPUT:

        - ``word`` -- a vector of a message space of the code.

        - ``encoder_name`` -- (default: ``None``) Name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        OUTPUT:

        - a vector of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word)
            (1, 1, 0, 0, 1, 1, 0)

        It is possible to manually choose the encoder amongst the list of the available ones::

            sage: C.encoders_available()
            ['GeneratorMatrix']
            sage: word = vector((0, 1, 1, 0))
            sage: C.encode(word, 'GeneratorMatrix')
            (1, 1, 0, 0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.encode(word)

    @cached_method
    def encoder(self, encoder_name=None, **kwargs):
        r"""
        Returns an encoder of ``self``.

        The returned encoder provided by this method is cached.

        This methods creates a new instance of the encoder subclass designated by ``encoder_name``.
        While it is also possible to do the same by directly calling the subclass' constructor,
        it is strongly advised to use this method to take advantage of the caching mechanism.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          returned. The default encoder of ``self`` will be used if
          default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the constructor of the encoder
          this method will return.

        OUTPUT:

        - an Encoder object.

        .. NOTE::

            The default encoder always has `F^{k}` as message space, with `k` the dimension
            of ``self`` and `F` the base ring of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoder()
            Generator matrix-based encoder for Linear code of length 7, dimension 4 over Finite Field of size 2

        We check that the returned encoder is cached::

            sage: C.encoder.is_in_cache()
            True

        If the name of an encoder which is not known by ``self`` is passed,
        an exception will be raised::

            sage: C.encoders_available()
            ['GeneratorMatrix']
            sage: C.encoder('NonExistingEncoder')
            Traceback (most recent call last):
            ...
            ValueError: Passed Encoder name not known
        """
        if encoder_name is None:
            encoder_name = self._default_encoder_name
        if encoder_name in self._registered_encoders:
            encClass = self._registered_encoders[encoder_name]
            E = encClass(self, **kwargs)
            return E
        else:
            raise ValueError("Passed Encoder name not known")

    def encoders_available(self, classes=False):
        r"""
        Returns a list of the available encoders' names for ``self``.

        INPUT:

        - ``classes`` -- (default: ``False``) if ``classes`` is set to ``True``, it also
          returns the encoders' classes associated with the encoders' names.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.encoders_available()
            ['GeneratorMatrix']

            sage: C.encoders_available(True)
            {'GeneratorMatrix':
            <class 'sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder'>}
        """
        if classes == True:
            return copy(self._registered_encoders)
        return self._registered_encoders.keys()

    def extended_code(self):
        r"""
        If ``self`` is a linear code of length `n` defined over `F` then this
        returns the code of length `n+1` where the last digit `c_n` satisfies
        the check condition `c_0+...+c_n=0`. If ``self`` is an `[n,k,d]`
        binary code then the extended code `C^{\vee}` is an `[n+1,k,d^{\vee}]`
        code, where `d^=d` (if d is even) and `d^{\vee}=d+1` (if `d` is odd).

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C
            [21, 18] Hamming Code over Finite Field in a of size 2^2
            sage: Cx = C.extended_code()
            sage: Cx
            Linear code of length 22, dimension 18 over Finite Field in a of size 2^2
        """
        G = self.generator_matrix()
        F = self.base_ring()
        k = len(G.rows())
        MS1 = MatrixSpace(F,k,1)
        ck_sums = [-sum(G.rows()[i]) for i in range(k)]
        last_col = MS1(ck_sums)
        Gx = G.augment(last_col)
        return LinearCode(Gx)

    def galois_closure(self, F0):
        r"""
        If ``self`` is a linear code defined over `F` and `F_0` is a subfield
        with Galois group `G = Gal(F/F_0)` then this returns the `G`-module
        `C^-` containing `C`.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: Cc = C.galois_closure(GF(2))
            sage: C; Cc
            [21, 18] Hamming Code over Finite Field in a of size 2^2
            Linear code of length 21, dimension 20 over Finite Field in a of size 2^2
            sage: c = C.basis()[2]
            sage: V = VectorSpace(GF(4,'a'),21)
            sage: c2 = V([x^2 for x in c.list()])
            sage: c2 in C
            False
            sage: c2 in Cc
            True
        """
        G = self.generator_matrix()
        F = self.base_ring()
        q = F.order()
        q0 = F0.order()
        a = log(q,q0)  # test if F/F0 is a field extension
        if not isinstance(a, Integer):
            raise ValueError("Base field must be an extension of given field %s"%F0)
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

    def __getitem__(self, i):
        r"""
        Returns the `i`-th codeword of this code.

        The implementation of this depends on the implementation of the
        :meth:`.__iter__` method.

        The implementation is as follows. Suppose that:

        - the primitive element of the base_ring of this code is `a`,
        - the prime subfield is `p`,
        - the field has order `p^m`,
        - the code has dimension `k`,
        - and the generator matrix is `G`.

        Then the :meth:`.__iter__` method returns the elements in this order:

        1. first, the following ordered list is returned:
           ``[i*a^0 * G[0] for i in range(p)]``
        2. Next, the following ordered list is returned:
           ``[i*a^0 * G[0] + a^1*G[0] for i in range(p)]``
        3. This continues till we get
           ``[(i*a^0 +(p-1)*a^1 +...+ (p-1)*a^(m-1))*G[0] for i in range(p)]``
        4. Then, we move to G[1]:
           ``[i*a^0 * G[0] + a^0*G[1] for i in range(p)]``,
         and so on.
         Hence the `i`-th element can be obtained by the p-adic expansion
         of `i` as ``[i_0, i_1, ...,i_{m-1}, i_m, i_{m+1}, ..., i_{km-1}].``
         The element that is generated is:

        .. math::

             \begin{aligned}
             & (i_0 a^0 + i_1 a^1 + \cdots + i_{m-1} a^{m-1}) G[0] + \\
             & (i_m a^0 + i_{m+1} a^1 + \cdots + i_{2m-1} a^{m-1}) G[1] + \\
             & \vdots\\
             & (i_{(k-1)m} a^0 + \cdots + i_{km-1} a^{m-1}) G[k-1]
             \end{aligned}

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C[24]
            (2, 2, 0, 1, 2, 2, 0)
            sage: C[24] == C.list()[24]
            True

        TESTS::

            sage: C = random_matrix(GF(25,'a'), 2, 7).row_space()
            sage: C = LinearCode(C.basis_matrix())
            sage: Clist = C.list()
            sage: all([C[i]==Clist[i] for i in xrange(len(C))])
            True

        Check that only the indices less than the size of the code are
        allowed::

            sage: C[25**2]
            Traceback (most recent call last):
            ...
            IndexError: The value of the index 'i' (=625) must be between
            0 and 'q^k -1' (=624), inclusive, where 'q' is the size of the
            base field and 'k' is the dimension of the code.

        Check that codewords are immutable. See :trac:`16338`::

            sage: C[0].is_immutable()
            True

        """
        # IMPORTANT: If the __iter__() function implementation is changed
        # then the implementation here must also be changed so that
        # list(self)[i] and self[i] both return the same element.

        from sage.rings.padics.factory import Zp
        F = self.base_ring()
        maxindex = F.order()**self.dimension()-1
        if i < 0 or i > maxindex:
            raise IndexError("The value of the index 'i' (={}) must be between "
                             "0 and 'q^k -1' (={}), inclusive, where 'q' is "
                             "the size of the base field and 'k' is the "
                             "dimension of the code.".format(i, maxindex))

        a = F.primitive_element()
        m = F.degree()
        p = F.prime_subfield().order()
        A = [a**k for k in xrange(m)]
        G = self.generator_matrix()
        N = self.dimension()*F.degree() # the total length of p-adic vector
        Z = Zp(p, N)
        ivec = Z(i).padded_list(N)

        codeword = 0
        row = 0
        for g in G:
            codeword += sum([ivec[j+row*m]*A[j] for j in xrange(m)])*g
            row += 1

        # The codewords for a specific code can not change. So, we set them
        # to be immutable.
        codeword.set_immutable()
        return codeword

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. The default encoder of ``self``
          will be used if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.generator_matrix()

    gen_mat = deprecated_function_alias(17973, generator_matrix)

    def generator_matrix_systematic(self):
        """
        Return a systematic generator matrix of the code.

        A generator matrix of a code is called systematic if it contains
        a set of columns forming an identity matrix.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
            sage: code.generator_matrix_systematic()
            [1 2 0]
            [0 0 1]
        """
        return self.generator_matrix().echelon_form()

    gen_mat_systematic = deprecated_function_alias(17973, generator_matrix_systematic)

    @cached_method
    def gens(self):
        r"""
        Returns the generators of this code as a list of vectors.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.gens()
             [(1, 0, 0, 0, 0, 1, 1), (0, 1, 0, 0, 1, 0, 1), (0, 0, 1, 0, 1, 1, 0), (0, 0, 0, 1, 1, 1, 1)]
        """
        return self.generator_matrix().rows()

    def genus(self):
        r"""
        Returns the "Duursma genus" of the code, `\gamma_C = n+1-k-d`.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3); C1
            [7, 4] Hamming Code over Finite Field of size 2
            sage: C1.genus()
            1
            sage: C2 = codes.HammingCode(GF(4,"a"), 2); C2
            [5, 3] Hamming Code over Finite Field in a of size 2^2
            sage: C2.genus()
            0

        Since all Hamming codes have minimum distance 3, these computations
        agree with the definition, `n+1-k-d`.
        """
        d = self.minimum_distance()
        n = self.length()
        k = self.dimension()
        gammaC = n+1-k-d
        return gammaC

    def __iter__(self):
        """
        Return an iterator over the elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: [list(c) for c in C if c.hamming_weight() < 4]
            [[0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1, 1],
             [0, 1, 0, 0, 1, 0, 1], [0, 0, 1, 0, 1, 1, 0],
             [1, 1, 1, 0, 0, 0, 0], [1, 0, 0, 1, 1, 0, 0],
             [0, 1, 0, 1, 0, 1, 0], [0, 0, 1, 1, 0, 0, 1]]

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: L = list(C)
            sage: L[10].is_immutable()
            True

        """
        from sage.modules.finite_submodule_iter import \
                                                FiniteFieldsubspace_iterator
        return FiniteFieldsubspace_iterator(self.generator_matrix(), immutable=True)

    @cached_method
    def information_set(self):
        """
        Return an information set of the code.

        Return value of this method is cached.

        A set of column positions of a generator matrix of a code
        is called an information set if the corresponding columns
        form a square matrix of full rank.

        OUTPUT:

        - Information set of a systematic generator matrix of the code.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,0,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix_systematic()
            [1 2 0]
            [0 0 1]
            sage: code.information_set()
            (0, 2)
        """
        return self.generator_matrix().transpose().pivot_rows()

    def is_permutation_automorphism(self,g):
        r"""
        Returns `1` if `g` is an element of `S_n` (`n` = length of self) and
        if `g` is an automorphism of self.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(3), 3)
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
        basis = self.generator_matrix().rows()
        H = self.parity_check_matrix()
        V = H.column_space()
        HGm = H*g.matrix()
        for c in basis:
            if HGm*c != V(0):
                return False
        return True

    def is_permutation_equivalent(self,other,algorithm=None):
        """
        Returns ``True`` if ``self`` and ``other`` are permutation equivalent
        codes and ``False`` otherwise.

        The ``algorithm="verbose"`` option also returns a permutation (if
        ``True``) sending ``self`` to ``other``.

        Uses Robert Miller's double coset partition refinement work.

        EXAMPLES::

            sage: P.<x> = PolynomialRing(GF(2),"x")
            sage: g = x^3+x+1
            sage: C1 = codes.CyclicCodeFromGeneratingPolynomial(7,g); C1
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C2 = codes.HammingCode(GF(2), 3); C2
            [7, 4] Hamming Code over Finite Field of size 2
            sage: C1.is_permutation_equivalent(C2)
            True
            sage: C1.is_permutation_equivalent(C2,algorithm="verbose")
            (True, (3,4)(5,7,6))
            sage: C1 = codes.RandomLinearCode(10,5,GF(2))
            sage: C2 = codes.RandomLinearCode(10,5,GF(3))
            sage: C1.is_permutation_equivalent(C2)
            False
        """
        from sage.groups.perm_gps.partn_ref.refinement_binary import NonlinearBinaryCodeStruct
        F = self.base_ring()
        F_o = other.base_ring()
        q = F.order()
        G = self.generator_matrix()
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
            if algorithm=="verbose":
                Sn = SymmetricGroup(n)
                return True, Sn([i+1 for i in ans])**(-1)
            return True
        return False

    def is_self_dual(self):
        """
        Returns ``True`` if the code is self-dual (in the usual Hamming inner
        product) and ``False`` otherwise.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.is_self_dual()
            True
            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.is_self_dual()
            False
        """
        return self == self.dual_code()


    def is_self_orthogonal(self):
        """
        Returns ``True`` if this code is self-orthogonal and ``False``
        otherwise.

        A code is self-orthogonal if it is a subcode of its dual.

        EXAMPLES::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: C.is_self_orthogonal()
            True
            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.is_self_orthogonal()
            False
            sage: C = codes.QuasiQuadraticResidueCode(11)  # optional - gap_packages (Guava package)
            sage: C.is_self_orthogonal()             # optional - gap_packages (Guava package)
            True
        """
        return self.is_subcode(self.dual_code())

    def is_galois_closed(self):
        r"""
        Checks if ``self`` is equal to its Galois closure.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,"a"), 3)
            sage: C.is_galois_closed()
            False
        """
        F = self.base_ring()
        p = F.characteristic()
        return self == self.galois_closure(GF(p))

    def is_subcode(self, other):
        """
        Returns ``True`` if ``self`` is a subcode of ``other``.

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: G1 = C1.generator_matrix()
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
            sage: C1 = codes.HammingCode(GF(9,"z"), 3)
            sage: G1 = C1.generator_matrix()
            sage: G2 = G1.matrix_from_rows([0,1,2])
            sage: C2 = LinearCode(G2)
            sage: C2.is_subcode(C1)
            True
        """
        G = self.generator_matrix()
        for r in G.rows():
            if not(r in other):
                return False
        return True

    def cardinality(self):
        r"""
        Return the size of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.cardinality()
            16
            sage: len(C)
            16
        """
        return self.base_ring().order()**self.dimension()

    __len__ = cardinality

    def length(self):
        r"""
        Returns the length of this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.length()
            7
        """
        return self._length

    def list(self):
        r"""
        Return a list of all elements of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Clist = C.list()
            sage: Clist[5]; Clist[5] in C
            (1, 0, 1, 0, 1, 0, 1)
            True
        """
        return [x for x in self]

    def _magma_init_(self, magma):
        r"""
        Retun a string representation in Magma of this linear code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: Cm = magma(C)                 # optional - magma, indirect doctest
            sage: Cm.MinimumWeight()            # optional - magma
            3

        """
        G = magma(self.generator_matrix())._ref()
        s = 'LinearCode(%s)' % G
        return s

    @cached_method
    def minimum_distance(self, algorithm=None):
        r"""
        Returns the minimum distance of this linear code.

        By default, this uses a GAP kernel function (in C and not part of
        Guava) written by Steve Linton.  If ``algorithm="guava"`` is set  and
        `q` is 2 or 3 then this uses a very fast program written in C written
        by CJ Tjhal. (This is much faster, except in some small examples.)

        Raises a ``ValueError`` in case there is no non-zero vector in this
        linear code.

        The minimum distance of the code is stored once it has been
        computed or provided during the initialization of :class:`LinearCode`.
        If ``algorithm`` is ``None`` and the stored value of minimum
        distance is found, then the stored value will be returned without
        recomputing the minimum distance again.

        INPUT:

        - ``algorithm`` - Method to be used, ``None``, ``"gap"``, or
          ``"guava"`` (default: ``None``).

        OUTPUT:

        - Integer, minimum distance of this code

        EXAMPLES::

            sage: MS = MatrixSpace(GF(3),4,7)
            sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.minimum_distance()
            3

        Once the minimum distance has been computed, it's value is stored.
        Hence the following command will return the value instantly,
        without further computations.::

            sage: C.minimum_distance()
            3

        If ``algorithm`` is provided, then the minimum distance will be
        recomputed even if there is a stored value from a previous run.::

            sage: C.minimum_distance(algorithm="gap")
            3
            sage: C.minimum_distance(algorithm="guava")  # optional - gap_packages (Guava package)
            3

        Another example.::

            sage: C = codes.HammingCode(GF(4,"a"), 2); C
            [5, 3] Hamming Code over Finite Field in a of size 2^2
            sage: C.minimum_distance()
            3

        TESTS::

            sage: C = codes.RandomLinearCode(5, 2, GF(4,"a"))
            sage: C.minimum_distance(algorithm='something')
            Traceback (most recent call last):
            ...
            ValueError: The algorithm argument must be one of None, 'gap' or 'guava'; got 'something'
        """
        # If the minimum distance has already been computed or provided by
        # the user then simply return the stored value.
        # This is done only if algorithm is None.
        if algorithm not in (None, "gap", "guava"):
            raise ValueError("The algorithm argument must be one of None, "
                        "'gap' or 'guava'; got '{0}'".format(algorithm))

        F = self.base_ring()
        q = F.order()
        G = self.generator_matrix()
        n = self.length()
        k = self.dimension()
        gapG = gap(G)
        if (q == 2 or q == 3) and algorithm=="guava":
            C = gapG.GeneratorMatCode(gap(F))
            d = C.MinimumWeight()
            return ZZ(d)
        Gstr = "%s*Z(%s)^0"%(gapG, q)
        return min_wt_vec_gap(Gstr,n,k,F).hamming_weight()

    def module_composition_factors(self, gp):
        r"""
        Prints the GAP record of the Meataxe composition factors module in
        Meataxe notation. This uses GAP but not Guava.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: gp = C.permutation_automorphism_group()

        Now type "C.module_composition_factors(gp)" to get the record printed.
        """
        F = self.base_ring()
        q = F.order()
        gens = gp.gens()
        G = self.generator_matrix()
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

    def permutation_automorphism_group(self, algorithm="partition"):
        r"""
        If `C` is an `[n,k,d]` code over `F`, this function computes the
        subgroup `Aut(C) \subset S_n` of all permutation automorphisms of `C`.
        The binary case always uses the (default) partition refinement
        algorithm of Robert Miller.

        Note that if the base ring of `C` is `GF(2)` then this is the full
        automorphism group. Otherwise, you could use
        :meth:`~sage.coding.linear_code.LinearCode.automorphism_group_gens`
        to compute generators of the full automorphism group.

        INPUT:

        - ``algorithm`` - If ``"gap"`` then GAP's MatrixAutomorphism function
          (written by Thomas Breuer) is used. The implementation combines an
          idea of mine with an improvement suggested by Cary Huffman. If
          ``"gap+verbose"`` then code-theoretic data is printed out at
          several stages of the computation. If ``"partition"`` then the
          (default) partition refinement algorithm of Robert Miller is used.
          Finally, if ``"codecan"`` then the partition refinement algorithm
          of Thomas Feulner is used, which also computes a canonical
          representative of ``self`` (call
          :meth:`~sage.coding.linear_code.LinearCode.canonical_representative`
          to access it).

        OUTPUT:

        - Permutation automorphism group

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,8)
            sage: G  = MS([[1,0,0,0,1,1,1,0],[0,1,1,1,0,0,0,0],[0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0]])
            sage: C  = LinearCode(G)
            sage: C
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            144
            sage: GG = C.permutation_automorphism_group("codecan")
            sage: GG == G
            True

        A less easy example involves showing that the permutation
        automorphism group of the extended ternary Golay code is the
        Mathieu group `M_{11}`.

        ::

            sage: C = codes.ExtendedTernaryGolayCode()
            sage: M11 = MathieuGroup(11)
            sage: M11.order()
            7920
            sage: G = C.permutation_automorphism_group()  # long time (6s on sage.math, 2011)
            sage: G.is_isomorphic(M11)                    # long time
            True
            sage: GG = C.permutation_automorphism_group("codecan") # long time
            sage: GG == G # long time
            True

        Other examples::

            sage: C = codes.ExtendedBinaryGolayCode()
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            244823040
            sage: C = codes.HammingCode(GF(2), 5)
            sage: G = C.permutation_automorphism_group()
            sage: G.order()
            9999360
            sage: C = codes.HammingCode(GF(3), 2); C
            [4, 2] Hamming Code over Finite Field of size 3
            sage: C.permutation_automorphism_group(algorithm="partition")
            Permutation Group with generators [(1,3,4)]
            sage: C = codes.HammingCode(GF(4,"z"), 2); C
            [5, 3] Hamming Code over Finite Field in z of size 2^2
            sage: G = C.permutation_automorphism_group(algorithm="partition"); G
            Permutation Group with generators [(1,3)(4,5), (1,4)(3,5)]
            sage: GG = C.permutation_automorphism_group(algorithm="codecan") # long time
            sage: GG == G # long time
            True
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            Permutation Group with generators [(1,3)(4,5), (1,4)(3,5)]
            sage: C = codes.TernaryGolayCode()
            sage: C.permutation_automorphism_group(algorithm="gap")  # optional - gap_packages (Guava package)
            Permutation Group with generators [(3,4)(5,7)(6,9)(8,11), (3,5,8)(4,11,7)(6,9,10), (2,3)(4,6)(5,8)(7,10), (1,2)(4,11)(5,8)(9,10)]

        However, the option ``algorithm="gap+verbose"``, will print out::

            Minimum distance: 5 Weight distribution: [1, 0, 0, 0, 0, 132, 132,
            0, 330, 110, 0, 24]

            Using the 132 codewords of weight 5 Supergroup size: 39916800

        in addition to the output of
        ``C.permutation_automorphism_group(algorithm="gap")``.
        """
        F = self.base_ring()
        q = F.order()
        G = self.generator_matrix() if 2*self.dimension() <= self.length() else self.dual_code().generator_matrix()
        n = len(G.columns())
        k = len(G.rows())
        if "gap" in algorithm:
            gap.load_package('guava')
            wts = self.spectrum()                                            # bottleneck 1
            nonzerowts = [i for i in range(len(wts)) if wts[i]!=0]
            Sn = SymmetricGroup(n)
            Gp = gap("SymmetricGroup(%s)"%n)               # initializing G in gap
            Gstr = str(gap(G))
            gap.eval("C:=GeneratorMatCode("+Gstr+",GF("+str(q)+"))")
            gap.eval("eltsC:=Elements(C)")
            if algorithm=="gap+verbose":
                print "\n Minimum distance: %s \n Weight distribution: \n %s"%(nonzerowts[1],wts)
            stop = 0                                          # only stop if all gens are autos
            for i in range(1,len(nonzerowts)):
                if stop == 1:
                    break
                wt = nonzerowts[i]
                if algorithm=="gap+verbose":
                    size = Gp.Size()
                    print "\n Using the %s codewords of weight %s \n Supergroup size: \n %s\n "%(wts[wt],wt,size)
                gap.eval("Cwt:=Filtered(eltsC,c->WeightCodeword(c)=%s)"%wt)   # bottleneck 2 (repeated
                gap.eval("matCwt:=List(Cwt,c->VectorCodeword(c))")            # for each i until stop = 1)
                if gap("Length(matCwt)") > 0:
                    A = gap("MatrixAutomorphisms(matCwt)")
                    G2 = gap("Intersection2(%s,%s)"%(str(A).replace("\n",""),str(Gp).replace("\n",""))) #  bottleneck 3
                    Gp = G2
                    if Gp.Size()==1:
                        return PermutationGroup([()])
                    autgp_gens = Gp.GeneratorsOfGroup()
                    gens = [Sn(str(x).replace("\n","")) for x in autgp_gens]
                    stop = 1                    # get ready to stop
                    for x in gens:              # if one of these gens is not an auto then don't stop
                        if not(self.is_permutation_automorphism(x)):
                            stop = 0
                            break
            G = PermutationGroup(gens)
            return G
        if algorithm=="partition":
            if q == 2:
                from sage.groups.perm_gps.partn_ref.refinement_binary import LinearBinaryCodeStruct
                B = LinearBinaryCodeStruct(G)
                autgp = B.automorphism_group()
                L = [[j+1 for j in gen] for gen in autgp[0]]
                AutGp = PermutationGroup(L)
            else:
                from sage.groups.perm_gps.partn_ref.refinement_matrices import MatrixStruct
                from sage.matrix.constructor import matrix
                weights = {}
                for c in self:
                    wt = c.hamming_weight()
                    if wt not in weights:
                        weights[wt] = [c]
                    else:
                        weights[wt].append(c)
                weights.pop(0)
                AutGps = []
                for wt, words in weights.iteritems():
                    M = MatrixStruct(matrix(words))
                    autgp = M.automorphism_group()
                    L = [[j+1 for j in gen] for gen in autgp[0]]
                    G = PermutationGroup(L)
                    AutGps.append(G)
                if len(AutGps) > 0:
                    AutGp = AutGps[0]
                    for G in AutGps[1:]:
                        AutGp = AutGp.intersection(G)
                else:
                    return PermutationGroup([])
            return AutGp
        if algorithm=="codecan":
            gens, _ = self.automorphism_group_gens("permutational")
            return PermutationGroup([x.get_perm() for x in gens])
        raise NotImplementedError("The only algorithms implemented currently are 'gap', 'gap+verbose', and 'partition'.")

    def permuted_code(self, p):
        r"""
        Returns the permuted code, which is equivalent to ``self`` via the
        column permutation ``p``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: G = C.permutation_automorphism_group(); G
            Permutation Group with generators [(4,5)(6,7), (4,6)(5,7), (2,3)(6,7), (2,4)(3,5), (1,2)(5,6)]
            sage: g = G("(2,3)(6,7)")
            sage: Cg = C.permuted_code(g)
            sage: Cg
            Linear code of length 7, dimension 4 over Finite Field of size 2
            sage: C.generator_matrix() == Cg.generator_matrix_systematic()
            True
        """
        F = self.base_ring()
        G = self.generator_matrix()
        n = len(G.columns())
        MS = MatrixSpace(F,n,n)
        Gp = G*MS(p.matrix().rows())
        return LinearCode(Gp)

    def punctured(self, L):
        r"""
        Returns the code punctured at the positions `L`,
        `L \subset \{1,2,...,n\}`. If this code `C` is of length `n` in
        GF(q) then the code `C^L` obtained from `C` by puncturing at the
        positions in `L` is the code of length `n-L` consisting of codewords
        of `C` which have their `i-th` coordinate deleted if `i \in L` and
        left alone if `i\notin L`:

        .. math::

            C^L = \{(c_{i_1},...,c_{i_N})\ |\ (c_1,...,c_n)\in C\},

        where `\{1,2,...,n\}-T = \{i_1,...,i_N\}`. In particular, if `L=\{j\}`
        then `C^L` is simply the code obtainen from `C` by deleting the `j-th`
        coordinate of each codeword. The code `C^L` is called the punctured
        code at `L`. The dimension of `C^L` can decrease if `|L|>d-1`.

        INPUT:

        - ``L`` - Subset of `\{1,...,n\}`, where `n` is the length of ``self``

        OUTPUT:

        - Linear code, the punctured code described above

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.punctured([1,2])
            Linear code of length 5, dimension 4 over Finite Field of size 2
        """
        G = self.generator_matrix()
        GL = G.matrix_from_columns([i for i in range(G.ncols()) if i not in L])
        r = GL.rank()
        if r < GL.nrows():
            GL.echelonize()
            GL = GL[:r]
        return LinearCode(GL)

    def random_element(self, *args, **kwds):
        """
        Returns a random codeword; passes other positional and keyword
        arguments to ``random_element()`` method of vector space.

        OUTPUT:

        - Random element of the vector space of this code

        EXAMPLES::

            sage: C = codes.HammingCode(GF(4,'a'), 3)
            sage: C.random_element() # random test
            (1, 0, 0, a + 1, 1, a, a, a + 1, a + 1, 1, 1, 0, a + 1, a, 0, a, a, 0, a, a, 1)

        Passes extra positional or keyword arguments through::

            sage: C.random_element(prob=.5, distribution='1/n') # random test
            (1, 0, a, 0, 0, 0, 0, a + 1, 0, 0, 0, 0, 0, 0, 0, 0, a + 1, a + 1, 1, 0, 0)

        TESTS:

        Test that the codeword returned is immutable (see :trac:`16469`)::

            sage: c = C.random_element()
            sage: c.is_immutable()
            True

        Test that codeword returned has the same parent as any non-random codeword
        (see :trac:`19653`)::

            sage: C = codes.RandomLinearCode(10, 4, GF(16, 'a'))
            sage: c1 = C.random_element()
            sage: c2 = C[1]
            sage: c1.parent() == c2.parent()
            True
        """
        E = self.encoder()
        M = E.message_space()
        m = M.random_element(*args, **kwds)
        c = E.encode(m)
        c.set_immutable()
        return c

    def redundancy_matrix(C):
        r"""
        If C is a linear [n,k,d] code then this function returns a
        `k \times (n-k)` matrix A such that G = (I,A) generates a code (in
        standard form) equivalent to C. If C is already in standard form and
        G = (I,A) is its generator matrix then this function simply returns
        that A.

        OUTPUT:

        - Matrix, the redundancy matrix

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.generator_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C.redundancy_matrix()
             [0 1 1]
             [1 0 1]
             [1 1 0]
             [1 1 1]
            sage: C.standard_form()[0].generator_matrix()
             [1 0 0 0 0 1 1]
             [0 1 0 0 1 0 1]
             [0 0 1 0 1 1 0]
             [0 0 0 1 1 1 1]
            sage: C = codes.HammingCode(GF(3), 2)
            sage: C.generator_matrix()
            [1 0 1 1]
            [0 1 1 2]
            sage: C.redundancy_matrix()
            [1 1]
            [1 2]
        """
        n = C.length()
        k = C.dimension()
        C1 = C.standard_form()[0]
        G1 = C1.generator_matrix()
        return G1.matrix_from_columns(range(k,n))

    def sd_duursma_data(C, i):
        r"""
        Returns the Duursma data `v` and `m` of this formally s.d. code `C`
        and the type number `i` in (1,2,3,4).  Does *not* check if this code
        is actually sd.

        INPUT:

        - ``i`` - Type number

        OUTPUT:

        - Pair ``(v, m)`` as in Duursma [D]_

        REFERENCES:

        .. [D] \I. Duursma, "Extremal weight enumerators and ultraspherical
           polynomials"

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),2,4)
            sage: G = MS([1,1,0,0,0,0,1,1])
            sage: C = LinearCode(G)
            sage: C == C.dual_code()  # checks that C is self dual
            True
            sage: for i in [1,2,3,4]: print C.sd_duursma_data(i)
            [2, -1]
            [2, -3]
            [2, -2]
            [2, -1]
        """
        n = C.length()
        d = C.minimum_distance()
        if i == 1:
            v = (n-4*d)//2 + 4
            m = d-3
        elif i == 2:
            v = (n-6*d)//8 + 3
            m = d-5
        elif i == 3:
            v = (n-4*d)//4 + 3
            m = d-4
        elif i == 4:
            v = (n-3*d)//2 + 3
            m = d-3
        return [v,m]

    def sd_duursma_q(C,i,d0):
        r"""
        INPUT:

        -  ``C`` - sd code; does *not* check if `C` is actually an sd code
        -  ``i`` - Type number, one of 1,2,3,4
        -  ``d0`` - Divisor, the smallest integer such that each `A_i > 0` iff
           `i` is divisible by `d0`

        OUTPUT:

        - Coefficients `q_0, q_1, ...` of `q(T)` as in Duursma [D]_

        REFERENCES:

        - [D] - I. Duursma, "Extremal weight enumerators and ultraspherical
          polynomials"

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
            sage: C2 = C1.extended_code(); C2
            Linear code of length 8, dimension 4 over Finite Field of size 2
            sage: C2.is_self_dual()
            True
            sage: C2.sd_duursma_q(1,1)
            2/5*T^2 + 2/5*T + 1/5
            sage: C2.sd_duursma_q(3,1)
            3/5*T^4 + 1/5*T^3 + 1/15*T^2 + 1/15*T + 1/15
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
            coefs = (c0*(1+T)**m*(1+4*T+6*T**2+4*T**3)**m*F).coefficients(sparse=False)
            qc = [coefs[j]/binomial(6*m+8*v,m+j) for j in range(4*m+8*v+1)]
            q = PR(qc)
        if i == 3:
            F = (3*T**2+4*T+1)**v*(1+3*T**2)**v
            # Note that: (3*T**2+4*T+1)(1+3*T**2)=(T+1)**4+8*T**3*(T+1)
            coefs = (c0*(1+3*T+3*T**2)**m*F).coefficients(sparse=False)
            qc = [coefs[j]/binomial(4*m+4*v,m+j) for j in range(2*m+4*v+1)]
            q = PR(qc)
        if i == 4:
            coefs = (c0*(1+2*T)**m*(4*T**2+2*T+1)**v).coefficients(sparse=False)
            qc = [coefs[j]/binomial(3*m+2*v,m+j) for j in range(m+2*v+1)]
            q = PR(qc)
        return q/q(1)

    def sd_zeta_polynomial(C, typ=1):
        r"""
        Returns the Duursma zeta function of a self-dual code using the
        construction in [D]_.

        INPUT:

        -  ``typ`` - Integer, type of this s.d. code; one of 1,2,3, or
           4 (default: 1)

        OUTPUT:

        -  Polynomial

        EXAMPLES::

            sage: C1 = codes.HammingCode(GF(2), 3)
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

        It is a general fact about Duursma zeta polynomials that `P(1) = 1`.

        REFERENCES:

        - [D] I. Duursma, "Extremal weight enumerators and ultraspherical
          polynomials"
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

    def shortened(self, L):
        r"""
        Returns the code shortened at the positions ``L``, where
        `L \subset \{1,2,...,n\}`.

        Consider the subcode `C(L)` consisting of all codewords `c\in C` which
        satisfy `c_i=0` for all `i\in L`. The punctured code `C(L)^L` is
        called the shortened code on `L` and is denoted `C_L`. The code
        constructed is actually only isomorphic to the shortened code defined
        in this way.

        By Theorem 1.5.7 in [HP], `C_L` is `((C^\perp)^L)^\perp`. This is used
        in the construction below.

        INPUT:

        - ``L`` - Subset of `\{1,...,n\}`, where `n` is the length of this code

        OUTPUT:

        - Linear code, the shortened code described above

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.shortened([1,2])
            Linear code of length 5, dimension 2 over Finite Field of size 2
        """
        Cd = self.dual_code()
        Cdp = Cd.punctured(L)
        return Cdp.dual_code()

    def spectrum(self, algorithm=None):
        r"""
        Returns the spectrum of ``self`` as a list.

        The default algorithm uses a GAP kernel function (in C) written by
        Steve Linton.

        INPUT:

        - ``algorithm`` - ``None``, ``"gap"``, ``"leon"``, or ``"binary"``;
          defaults to ``"gap"`` except in the binary case.  If ``"gap"`` then
          uses the GAP function, if ``"leon"`` then uses Jeffrey Leon's
          software via Guava, and if ``"binary"`` then uses Sage native Cython
          code

        - List, the spectrum

        The optional algorithm (``"leon"``) may create a stack smashing error
        and a traceback but should return the correct answer. It appears to run
        much faster than the GAP algorithm in some small examples and much
        slower than the GAP algorithm in other larger examples.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G = MS([[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: F.<z> = GF(2^2,"z")
            sage: C = codes.HammingCode(F, 2); C
            [5, 3] Hamming Code over Finite Field in z of size 2^2
            sage: C.spectrum()
            [1, 0, 0, 30, 15, 18]
            sage: C = codes.HammingCode(GF(2), 3); C
            [7, 4] Hamming Code over Finite Field of size 2
            sage: C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.spectrum(algorithm="gap")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.spectrum(algorithm="binary")
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C = codes.HammingCode(GF(3), 3); C
            [13, 10] Hamming Code over Finite Field of size 3
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(GF(5), 2); C
            [6, 4] Hamming Code over Finite Field of size 5
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True
            sage: C = codes.HammingCode(GF(7), 2); C
            [8, 6] Hamming Code over Finite Field of size 7
            sage: C.spectrum() == C.spectrum(algorithm="leon")   # optional - gap_packages (Guava package)
            True

        """
        if algorithm is None:
            if self.base_ring().order() == 2:
                algorithm = "binary"
            else:
                algorithm = "gap"
        n = self.length()
        F = self.base_ring()
        G = self.generator_matrix()
        if algorithm=="gap":
            Gstr = G._gap_init_()
            spec = wtdist_gap(Gstr,n,F)
            return spec
        elif algorithm=="binary":
            from sage.coding.binary_code import weight_dist
            return weight_dist(self.generator_matrix())
        elif algorithm=="leon":
            if not(F.order() in [2,3,5,7]):
                raise NotImplementedError("The algorithm 'leon' is only implemented for q = 2,3,5,7.")
            # The GAP command DirectoriesPackageLibrary tells the location of the latest
            # version of the Guava libraries, so gives us the location of the Guava binaries too.
            guava_bin_dir = gap.eval('DirectoriesPackagePrograms("guava")[1]')
            guava_bin_dir = guava_bin_dir[guava_bin_dir.index('"') + 1:guava_bin_dir.rindex('"')]
            input = code2leon(self) + "::code"
            import os, subprocess
            lines = subprocess.check_output([os.path.join(guava_bin_dir, 'wtdist'), input])
            import StringIO  # to use the already present output parser
            wts = [0]*(n+1)
            s = 0
            for L in StringIO.StringIO(lines).readlines():
                L = L.strip()
                if len(L) > 0:
                    o = ord(L[0])
                    if o >= 48 and o <= 57:
                        wt, num = L.split()
                        wts[eval(wt)] = eval(num)
            return wts
        else:
            raise NotImplementedError("The only algorithms implemented currently are 'gap', 'leon' and 'binary'.")

    def standard_form(self):
        r"""
        Returns the standard form of this linear code.

        An `[n,k]` linear code with generator matrix `G` in standard form is
        the row-reduced echelon form of `G` is `(I,A)`, where `I` denotes the
        `k \times k` identity matrix and `A` is a `k \times (n-k)` block. This
        method returns a pair `(C,p)` where `C` is a code permutation
        equivalent to ``self`` and `p` in `S_n`, with `n` the length of `C`,
        is the permutation sending ``self`` to `C`. This does not call GAP.

        Thanks to Frank Luebeck for (the GAP version of) this code.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
            sage: Cs,p = C.standard_form()
            sage: p
            ()
            sage: MS = MatrixSpace(GF(3),3,7)
            sage: G = MS([[1,0,0,0,1,1,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,1]])
            sage: C = LinearCode(G)
            sage: Cs, p = C.standard_form()
            sage: p
            (3,7)
            sage: Cs.generator_matrix()
             [1 0 0 0 1 1 0]
             [0 1 0 1 0 1 0]
             [0 0 1 0 0 0 0]
        """
        from sage.coding.code_constructions import permutation_action as perm_action
        mat = self.generator_matrix()
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
            if j < d and i != j+1:
                perm = perm *G([(i,j+1)])
        if perm != G([()]):
            for i in range(k):
                r = M.rows()[i]
                A.append(perm_action(perm,r))
        if perm == G([()]):
            A = M
        return LinearCode(MS(A)), perm

    def support(self):
        r"""
        Returns the set of indices `j` where `A_j` is nonzero, where
        spectrum(self) = `[A_0,A_1,...,A_n]`.

        OUTPUT:

        - List of integers

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.spectrum()
            [1, 0, 0, 7, 7, 0, 0, 1]
            sage: C.support()
            [0, 3, 4, 7]
        """
        n = self.length()
        F = self.base_ring()
        V = VectorSpace(F,n+1)
        return V(self.spectrum()).support()

    def syndrome(self, r):
        r"""
        Returns the syndrome of ``r``.

        The syndrome of ``r`` is the result of `H \times r` where `H` is
        the parity check matrix of ``self``. If ``r`` belongs to ``self``,
        its syndrome equals to the zero vector.

        INPUT:

        - ``r`` -- a vector of the same length as ``self``

        OUTPUT:

        - a column vector

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: r = vector(GF(2), (1,0,1,0,1,0,1))
            sage: r in C
            True
            sage: C.syndrome(r)
            (0, 0, 0)

        If ``r`` is not a codeword, its syndrome is not equal to zero::

            sage: r = vector(GF(2), (1,0,1,0,1,1,1))
            sage: r in C
            False
            sage: C.syndrome(r)
            (0, 1, 1)

        Syndrome computation works fine on bigger fields::

            sage: C = codes.RandomLinearCode(12, 4, GF(59))
            sage: r = C.random_element()
            sage: C.syndrome(r)
            (0, 0, 0, 0, 0, 0, 0, 0)
        """
        return self.parity_check_matrix()*r

    def unencode(self, c, encoder_name=None, nocheck=False, **kwargs):
        r"""
        Returns the message corresponding to ``c``.

        This is the inverse of :meth:`encode`.

        INPUT:

        - ``c`` -- a codeword of ``self``.

        - ``encoder_name`` -- (default: ``None``) name of the decoder which will be used
          to decode ``word``. The default decoder of ``self`` will be used if
          default value is kept.

        - ``nocheck`` -- (default: ``False``) checks if ``c`` is in ``self``. You might set
          this to ``True`` to disable the check for saving computation. Note that if ``c`` is
          not in ``self`` and ``nocheck = True``, then the output of :meth:`unencode` is
          not defined (except that it will be in the message space of ``self``).

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        OUTPUT:

        - an element of the message space of ``encoder_name`` of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: c = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: C.unencode(c)
            (0, 1, 1, 0)
        """
        E = self.encoder(encoder_name, **kwargs)
        return E.unencode(c, nocheck)

    def weight_enumerator(self, names="xy", name2=None):
        """
        Returns the weight enumerator of the code.

        INPUT:

        - ``names`` - String of length 2, containing two variable names
          (default: ``"xy"``). Alternatively, it can be a variable name or
          a string, or a tuple of variable names or strings.

        - ``name2`` - string or symbolic variable (default: ``None``).
          If ``name2`` is provided then it is assumed that ``names``
          contains only one variable.

        OUTPUT:

        - Polynomial over `\QQ`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.weight_enumerator()
            x^7 + 7*x^4*y^3 + 7*x^3*y^4 + y^7
            sage: C.weight_enumerator(names="st")
            s^7 + 7*s^4*t^3 + 7*s^3*t^4 + t^7
            sage: (var1, var2) = var('var1, var2')
            sage: C.weight_enumerator((var1, var2))
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7
            sage: C.weight_enumerator(var1, var2)
            var1^7 + 7*var1^4*var2^3 + 7*var1^3*var2^4 + var2^7

        """
        if name2 is not None:
            # We assume that actual variable names or strings are provided
            # for names if names2 is also provided. That is, names is not
            # a tuple or a list. Otherwise, PolynomialRing will return error
            names = (names, name2)
        spec = self.spectrum()
        n = self.length()
        R = PolynomialRing(QQ,2,names)
        x,y = R.gens()
        we = sum([spec[i]*x**(n-i)*y**i for i in range(n+1)])
        return we

    @cached_method
    def zero(self):
        r"""
        Returns the zero vector of ``self``.

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zero()
            (0, 0, 0, 0, 0, 0, 0)
            sage: C.sum(()) # indirect doctest
            (0, 0, 0, 0, 0, 0, 0)
            sage: C.sum((C.gens())) # indirect doctest
            (1, 1, 1, 1, 1, 1, 1)
        """
        return self.ambient_space().zero()

    def zeta_polynomial(self, name="T"):
        r"""
        Returns the Duursma zeta polynomial of this code.

        Assumes that the minimum distances of this code and its dual are
        greater than 1.  Prints a warning to ``stdout`` otherwise.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Polynomial over `\QQ`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zeta_polynomial()
            2/5*T^2 + 2/5*T + 1/5
            sage: C = best_known_linear_code(6,3,GF(2))  # optional - gap_packages (Guava package)
            sage: C.minimum_distance()                   # optional - gap_packages (Guava package)
            3
            sage: C.zeta_polynomial()                    # optional - gap_packages (Guava package)
            2/5*T^2 + 2/5*T + 1/5
            sage: C = codes.HammingCode(GF(2), 4)
            sage: C.zeta_polynomial()
            16/429*T^6 + 16/143*T^5 + 80/429*T^4 + 32/143*T^3 + 30/143*T^2 + 2/13*T + 1/13
            sage: F.<z> = GF(4,"z")
            sage: MS = MatrixSpace(F, 3, 6)
            sage: G = MS([[1,0,0,1,z,z],[0,1,0,z,1,z],[0,0,1,z,z,1]])
            sage: C = LinearCode(G)  # the "hexacode"
            sage: C.zeta_polynomial()
            1

        REFERENCES:

        - I. Duursma, "From weight enumerators to zeta functions", in
          Discrete Applied Mathematics, vol. 111, no. 1-2, pp. 55-73, 2001.
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
        P_coeffs = []
        for i in range(len(b)):
           if i == 0:
              P_coeffs.append(b[0])
           if i == 1:
              P_coeffs.append(b[1] - (q+1)*b[0])
           if i>1:
              P_coeffs.append(b[i] - (q+1)*b[i-1] + q*b[i-2])
        P = sum([P_coeffs[i]*T**i for i in range(r+1)])
        return RT(P)/RT(P)(1)

    def zeta_function(self, name="T"):
        r"""
        Returns the Duursma zeta function of the code.

        INPUT:

        - ``name`` - String, variable name (default: ``"T"``)

        OUTPUT:

        - Element of `\QQ(T)`

        EXAMPLES::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: C.zeta_function()
            (2/5*T^2 + 2/5*T + 1/5)/(2*T^2 - 3*T + 1)
        """
        P =  self.zeta_polynomial()
        q = (self.base_ring()).characteristic()
        RT = PolynomialRing(QQ,"%s"%name)
        T = RT.gen()
        return P/((1-T)*(1-q*T))

    weight_distribution = spectrum

def LinearCodeFromVectorSpace(V, d=None):
    """
    Simply converts a vector subspace `V` of `GF(q)^n` into a `LinearCode`.

    INPUT:

    - ``V`` -- The vector space

    - ``d`` -- (Optional, default: ``None``) the minimum distance of the
      code, if known. This is an optional parameter.

    .. note::
        The veracity of the minimum distance ``d``, if provided, is not
        checked.

    EXAMPLES::

        sage: V = VectorSpace(GF(2), 8)
        sage: L = V.subspace([[1,1,1,1,0,0,0,0],[0,0,0,0,1,1,1,1]])
        sage: C = LinearCodeFromVectorSpace(L)
        sage: C.generator_matrix()
        [1 1 1 1 0 0 0 0]
        [0 0 0 0 1 1 1 1]
        sage: C.minimum_distance()
        4

    Here, we provide the minimum distance of the code.::

        sage: C = LinearCodeFromVectorSpace(L, d=4)
        sage: C.minimum_distance()
        4
    """
    F = V.base_ring()
    B = V.basis()
    n = len(B[0].list())
    k = len(B)
    MS = MatrixSpace(F,k,n)
    G = MS([B[i].list() for i in range(k)])
    return LinearCode(G, d=d)

############################ linear codes python class ########################

class LinearCode(AbstractLinearCode):
    r"""
    Linear codes over a finite field or finite ring, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear codes. This
    means that basic operations on the code, such as the computation of the
    minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating `codes.<tab>`. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.

    See :wikipedia:`Linear_code` for more information on unstructured linear codes.

    INPUT:

    - ``generator`` -- a generator matrix over a finite field (``G`` can be
      defined over a finite ring but the matrices over that ring must have
      certain attributes, such as ``rank``); or a code over a finite field

    - ``d`` -- (optional, default: ``None``) the minimum distance of the code

    .. NOTE::

        The veracity of the minimum distance ``d``, if provided, is not
        checked.

    EXAMPLES::

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

    The minimum distance of the code, if known, can be provided as an
    optional parameter.::

        sage: C  = LinearCode(G, d=3)
        sage: C.minimum_distance()
        3

    Another example.::

        sage: MS = MatrixSpace(GF(5),4,7)
        sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
        sage: C  = LinearCode(G)
        sage: C
        Linear code of length 7, dimension 4 over Finite Field of size 5
        
    Providing a code as the parameter in order to "forget" its structure (see
    :trac:`20198`)::
    
        sage: C = codes.GeneralizedReedSolomonCode(GF(23).list(), 12)
        sage: LinearCode(C)
        Linear code of length 23, dimension 12 over Finite Field of size 23
        
    Another example::
    
        sage: C = codes.HammingCode(GF(7), 3)
        sage: C
        [57, 54] Hamming Code over Finite Field of size 7
        sage: LinearCode(C)
        Linear code of length 57, dimension 54 over Finite Field of size 7

    AUTHORS:

    - David Joyner (11-2005)
    - Charles Prior (03-2016): :trac:`20198`, LinearCode from a code
    """
    def __init__(self, generator, d=None):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)    # indirect doctest
            sage: C
            Linear code of length 7, dimension 4 over Finite Field of size 2

        The minimum distance of the code, if known, can be provided as an
        optional parameter.::

            sage: C  = LinearCode(G, d=3)
            sage: C.minimum_distance()
            3

        TESTS::

            sage: C = codes.HammingCode(GF(2), 3)
            sage: TestSuite(C).run()

        Check that it works even with input matrix with non full rank (see
        :trac:`17452`)::

            sage: K.<a> = GF(4)
            sage: G = matrix([[a, a + 1, 1, a + 1, 1, 0, 0],
            ....:             [0, a, a + 1, 1, a + 1, 1, 0],
            ....:             [0, 0, a, a + 1, 1, a + 1, 1],
            ....:             [a + 1, 0, 1, 0, a + 1, 1, a + 1],
            ....:             [a, a + 1, a + 1, 0, 0, a + 1, 1],
            ....:             [a + 1, a, a, 1, 0, 0, a + 1],
            ....:             [a, a + 1, 1, a + 1, 1, 0, 0]])
            sage: C = LinearCode(G)
            sage: C.basis()
            [(1, 0, 0, a + 1, 0, 1, 0),
             (0, 1, 0, 0, a + 1, 0, 1),
             (0, 0, 1, a, a + 1, a, a + 1)]
            sage: C.minimum_distance()
            3

        Forbid the zero vector space (see :trac:`17452` and :trac:`6486`)::

            sage: G = matrix(GF(2), [[0,0,0]])
            sage: C = LinearCode(G)
            Traceback (most recent call last):
            ...
            ValueError: this linear code contains no non-zero vector
        """
        
        base_ring = generator.base_ring()
        if not base_ring.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")
        
        try:
            basis = generator.row_space().basis() # generator matrix case
            
            # if the matrix does not have full rank we replace it
            if len(basis) != generator.nrows():
                from sage.matrix.constructor import matrix
                generator = matrix(base_ring, basis)
    
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no non-zero vector")
        except AttributeError:
            # Assume input is an AbstractLinearCode, extract its generator matrix
            generator = generator.generator_matrix()

        super(LinearCode, self).__init__(base_ring, generator.ncols(), "GeneratorMatrix", "Syndrome")
        self._generator_matrix = generator
        self._dimension = generator.rank()
        self._minimum_distance = d

    def _repr_(self):
        r"""
        See the docstring for :meth:`LinearCode`.

        EXAMPLES::

            sage: MS = MatrixSpace(GF(2),4,7)
            sage: G  = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
            sage: C  = LinearCode(G)
            sage: C                     # indirect doctest
            Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Linear code of length %s, dimension %s over %s"%(self.length(), self.dimension(), self.base_ring())

    def __hash__(self):
        r"""
        Returns the hash value of ``self``.

        EXAMPLES::

        sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
        sage: C = LinearCode(G)
        sage: hash(C) #random
        9015017528451745710

        If ``C1`` and ``C2`` are two codes which only differ by the coefficients of their
        generator matrices, their hashes are different (we check that the bug found in trac #18813
        is fixed)::

        sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
        sage: C1 = LinearCode(G)
        sage: G = Matrix(GF(2), [[1,0,0,1,0,1,0],[0,1,0,0,1,0,0],[0,0,1,1,0,1,0],[0,0,0,0,0,0,1]])
        sage: C2 = LinearCode(G)
        sage: hash(C1) != hash(C2)
        True
        """
        Str = str(self)
        G = str(self.generator_matrix()) #str because mutable matrices are unhashable
        return hash((Str, G)) ^ hash(Str) ^ hash(G)

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Returns a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = matrix(GF(3),2,[1,-1,1,-1,1,1])
            sage: code = LinearCode(G)
            sage: code.generator_matrix()
            [1 2 1]
            [2 1 1]
        """
        if encoder_name is None or encoder_name is 'GeneratorMatrix':
            return self._generator_matrix
        return super(LinearCode, self).generator_matrix(encoder_name, **kwargs)










####################### encoders ###############################
class LinearCodeGeneratorMatrixEncoder(Encoder):
    r"""
    Encoder based on generator_matrix for Linear codes.

    This is the default encoder of a generic linear code, and should never be used for other codes than
    :class:`LinearCode`.

    INPUT:

    - ``code`` -- The associated :class:`LinearCode` of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        super(LinearCodeGeneratorMatrixEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeGeneratorMatrixEncoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: E1 = LinearCode(G).encoder()
            sage: E2 = LinearCode(G).encoder()
            sage: E1 == E2
            True
        """
        return isinstance(other, LinearCodeGeneratorMatrixEncoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E
            Generator matrix-based encoder for Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Generator matrix-based encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: latex(E)
            \textnormal{Generator matrix-based encoder for }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Generator matrix-based encoder for }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
            sage: E.generator_matrix()
            [1 1 1 0 0 0 0]
            [1 0 0 1 1 0 0]
            [0 1 0 1 0 1 0]
            [1 1 0 1 0 0 1]
        """
        return self.code().generator_matrix()










class LinearCodeParityCheckEncoder(Encoder):
    r"""
    Encoder based on :meth:`parity_check_matrix` for Linear codes.

    It constructs the generator matrix through the parity check matrix.

    INPUT:

    - ``code`` -- The associated code of this encoder.
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E
            Parity check matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        super(LinearCodeParityCheckEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E
            Parity check matrix-based encoder for the Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Parity check matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: latex(E)
            \textnormal{Parity check matrix-based encoder for the }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Parity check matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: E = codes.encoders.LinearCodeParityCheckEncoder(C)
            sage: E.generator_matrix()
            [1 0 0 0 0 1 1]
            [0 1 0 0 1 0 1]
            [0 0 1 0 1 1 0]
            [0 0 0 1 1 1 1]
        """
        return self.code().parity_check_matrix().right_kernel_matrix()

####################### decoders ###############################
class LinearCodeSyndromeDecoder(Decoder):
    r"""
    Constructs a decoder for Linear Codes based on syndrome lookup table.

    The decoding algorithm works as follows:

    - First, a lookup table is built by computing the syndrome of every error
      pattern of weight up to ``maximum_error_weight``.
    - Then, whenever one tries to decode a word ``r``, the syndrome of ``r`` is
      computed. The corresponding error pattern is recovered from the
      pre-computed lookup table.
    - Finally, the recovered error pattern is subtracted from ``r`` to recover
      the original word.

    ``maximum_error_weight`` need never exceed the covering radius of the code,
    since there are then always lower-weight errors with the same syndrome. If
    one sets ``maximum_error_weight`` to a value greater than the covering
    radius, then the covering radius will be determined while building the
    lookup-table. This lower value is then returned if you query
    ``decoding_radius`` after construction.

    If ``maximum_error_weight`` is left unspecified or set to a number at least
    the covering radius of the code, this decoder is complete, i.e. it decodes
    every vector in the ambient space.

    NOTE:

    Constructing the lookup table takes time exponential in the length of the
    code and the size of the code's base field. Afterwards, the individual
    decodings are fast.

    INPUT:

    - ``code`` -- A code associated to this decoder

    - ``maximum_error_weight`` -- (default: ``None``) the maximum number of
      errors to look for when building the table. An error is raised if it is
      set greater than `n-k`, since this is an upper bound on the covering
      radius on any linear code. If ``maximum_error_weight`` is kept
      unspecified, it will be set to `n - k`, where `n` is the length of
      ``code`` and `k` its dimension.

    EXAMPLES::

        sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
        sage: C = LinearCode(G)
        sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
        sage: D
        Syndrome decoder for Linear code of length 9, dimension 3 over Finite Field of size 3 handling errors of weight up to 4

    If one wants to correct up to a lower number of errors, one can do as follows::

        sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight=2)
        sage: D
        Syndrome decoder for Linear code of length 9, dimension 3 over Finite Field of size 3 handling errors of weight up to 2

    If one checks the list of types of this decoder before constructing it,
    one will notice it contains the keyword ``dynamic``.
    Indeed, the behaviour of the syndrome decoder depends on the maximum
    error weight one wants to handle, and how it compares to the minimum
    distance and the covering radius of ``code``.
    In the following examples, we illustrate this property by computing
    different instances of syndrome decoder for the same code.

    We choose the following linear code, whose covering radius equals to 4
    and minimum distance to 5 (half the minimum distance is 2)::

        sage: G = matrix(GF(5), [[1, 0, 0, 0, 0, 4, 3, 0, 3, 1, 0],
        ....:                    [0, 1, 0, 0, 0, 3, 2, 2, 3, 2, 1],
        ....:                    [0, 0, 1, 0, 0, 1, 3, 0, 1, 4, 1],
        ....:                    [0, 0, 0, 1, 0, 3, 4, 2, 2, 3, 3],
        ....:                    [0, 0, 0, 0, 1, 4, 2, 3, 2, 2, 1]])
        sage: C = LinearCode(G)

    In the following examples, we illustrate how the choice of
    ``maximum_error_weight`` influences the types of the instance of
    syndrome decoder, alongside with its decoding radius.

    We build a first syndrome decoder, and pick a ``maximum_error_weight``
    smaller than both the covering radius and half the minimum distance::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 1)
        sage: D.decoder_type()
        {'always-succeed', 'bounded_distance', 'hard-decision', 'unique'}
        sage: D.decoding_radius()
        1

    In that case, we are sure the decoder will always succeed. It is also
    a bounded distance decoder.

    We now build another syndrome decoder, and this time,
    ``maximum_error_weight`` is chosen to be bigger than half the minimum distance,
    but lower than the covering radius::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 3)
        sage: D.decoder_type()
        {'bounded_distance', 'hard-decision', 'might-error', 'unique'}
        sage: D.decoding_radius()
        3

    Here, we still get a bounded distance decoder.
    But because we have a maximum error weight bigger than half the
    minimum distance, we know it might return a codeword which was not
    the original codeword.

    And now, we build a third syndrome decoder, whose ``maximum_error_weight``
    is bigger than both the covering radius and half the minimum distance::

        sage: D = C.decoder("Syndrome", maximum_error_weight = 5)
        sage: D.decoder_type()
        {'complete', 'hard-decision', 'might-error', 'unique'}
        sage: D.decoding_radius()
        4

    In that case, the decoder might still return an unexpected codeword, but
    it is now complete. Note the decoding radius is equal to 4: it was
    determined while building the syndrome lookup table that any error with
    weight more than 4 will be decoded incorrectly. That is because the covering
    radius for the code is 4.

    The minimum distance and the covering radius are both determined while
    computing the syndrome lookup table. They user did not explicitly ask to
    compute these on the code ``C``. The dynamic typing of the syndrome decoder
    might therefore seem slightly surprising, but in the end is quite
    informative.
    """

    def __init__(self, code, maximum_error_weight=None):
        r"""
        TESTS:

        If ``maximum_error_weight`` is greater or equal than `n-k`, where `n`
        is ``code``'s length, and `k` is ``code``'s dimension,
        an error is raised::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, 42)
            Traceback (most recent call last):
            ...
            ValueError: maximum_error_weight has to be less than code's length minus its dimension
        """
        n_minus_k = code.length() - code.dimension()
        if maximum_error_weight == None:
            self._maximum_error_weight = n_minus_k
        elif not isinstance(maximum_error_weight, (Integer, int)):
            raise ValueError("maximum_error_weight has to be a Sage integer or a Python int")
        elif maximum_error_weight > n_minus_k:
            raise ValueError("maximum_error_weight has to be less than code's length minus its dimension")
        else:
            self._maximum_error_weight = maximum_error_weight
        super(LinearCodeSyndromeDecoder, self).__init__(code, code.ambient_space(),\
                code._default_encoder_name)
        self._lookup_table = self._build_lookup_table()

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeSyndromeDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: D1 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: D2 = codes.decoders.LinearCodeSyndromeDecoder(LinearCode(G))
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearCodeSyndromeDecoder)\
                and self.code() == other.code()\
                and self.maximum_error_weight() == other.maximum_error_weight()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D
            Syndrome decoder for Linear code of length 9, dimension 3 over Finite Field of size 3 handling errors of weight up to 4
        """
        return "Syndrome decoder for %s handling errors of weight up to %s" % (self.code(), self.maximum_error_weight())

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: latex(D)
            \textnormal{Syndrome decoder for [9, 3]\textnormal{ Linear code over }\Bold{F}_{3} handling errors of weight up to 4}
        """
        return "\\textnormal{Syndrome decoder for %s handling errors of weight up to %s}" % (self.code()._latex_(), self.maximum_error_weight())

    @cached_method
    def _build_lookup_table(self):
        r"""
        Builds lookup table for all possible error patterns of weight up to :meth:`maximum_error_weight`.

        EXAMPLES::

            sage: G = Matrix(GF(3),[
            ....:   [1, 0, 0, 0, 2, 2, 1, 1],
            ....:   [0, 1, 0, 0, 0, 0, 1, 1],
            ....:   [0, 0, 1, 0, 2, 0, 0, 2],
            ....:   [0, 0, 0, 1, 0, 2, 0, 1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight = 1)
            sage: D._build_lookup_table()
            {(0, 0, 0, 0): (0, 0, 0, 0, 0, 0, 0, 0),
             (0, 0, 0, 1): (0, 0, 0, 0, 1, 0, 0, 0),
             (0, 0, 0, 2): (0, 0, 0, 0, 2, 0, 0, 0),
             (0, 0, 1, 0): (0, 0, 1, 0, 0, 0, 0, 0),
             (0, 0, 1, 2): (0, 0, 0, 0, 0, 0, 0, 1),
             (0, 0, 2, 0): (0, 0, 2, 0, 0, 0, 0, 0),
             (0, 0, 2, 1): (0, 0, 0, 0, 0, 0, 0, 2),
             (0, 1, 0, 0): (0, 1, 0, 0, 0, 0, 0, 0),
             (0, 1, 1, 2): (0, 0, 0, 0, 0, 0, 2, 0),
             (0, 2, 0, 0): (0, 2, 0, 0, 0, 0, 0, 0),
             (0, 2, 2, 1): (0, 0, 0, 0, 0, 0, 1, 0),
             (1, 0, 0, 0): (1, 0, 0, 0, 0, 0, 0, 0),
             (1, 2, 0, 2): (0, 0, 0, 0, 0, 1, 0, 0),
             (1, 2, 2, 0): (0, 0, 0, 1, 0, 0, 0, 0),
             (2, 0, 0, 0): (2, 0, 0, 0, 0, 0, 0, 0),
             (2, 1, 0, 1): (0, 0, 0, 0, 0, 2, 0, 0),
             (2, 1, 1, 0): (0, 0, 0, 2, 0, 0, 0, 0)}
        """
        t = self._maximum_error_weight
        self._code_covering_radius = None
        self._code_minimum_distance = None
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        C = self.code()
        n = C.length()
        k = C.dimension()
        H = C.parity_check_matrix()
        F = C.base_ring()
        l = F.list()
        zero = F.zero()
        #Builds a list of generators of all error positions for all
        #possible error weights
        if zero in l:
            l.remove(zero)
        # Remember to include the no-error-vector to handle codes of minimum
        # distance 1 gracefully
        zero_syndrome = vector(F,[F.zero()]*(n-k))
        zero_syndrome.set_immutable()
        lookup = { zero_syndrome : vector(F,[F.zero()]*n) }
        error_position_tables = [cartesian_product([l]*i) for i in range(1, t+1)]
        first_collision = True
        #Filling the lookup table
        for i in range(1, t+1):
            stop = True
            patterns = Subsets(range(n), i)
            basic = vector(F, n)
            for p in patterns:
                for error in error_position_tables[i-1]:
                    ind = 0
                    e = copy(basic)
                    for pos in p:
                        e[pos] = error[ind]
                        ind += 1
                    s = H * e
                    s.set_immutable()
                    try:
                        e_cur = lookup[s]
                        #if this is the first time we see a collision
                        #we learn the minimum distance of the code
                        if first_collision:
                            self._code_minimum_distance = e.hamming_weight() + e_cur.hamming_weight()
                            first_collision = False
                    except KeyError:
                        stop = False
                        lookup[s] = copy(e)
            #if we reached the early termination condition
            #we learn the covering radius of the code
            if stop:
                self._code_covering_radius = i - 1
                self._maximum_error_weight = self._code_covering_radius
                break
        # Update decoder types depending on whether we are decoding up to covering radius
        if self._code_covering_radius:
            self._decoder_type.add("complete")
        else:
            self._decoder_type.add("bounded_distance")
        # Update decoder types depending on whether we are decoding beyond d/2
        if self._code_minimum_distance:
            if t == (self._code_minimum_distance-1)//2:
                self._decoder_type.add("minimum-distance", "always-succeed")
            else:
                # then t > (d-1)/2
                self._decoder_type.add("might-error")
        else:
            self._decoder_type.add("always-succeed")
        return lookup


    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space

        EXAMPLES::

            sage: G = Matrix(GF(3),[
            ....:   [1, 0, 0, 0, 2, 2, 1, 1],
            ....:   [0, 1, 0, 0, 0, 0, 1, 1],
            ....:   [0, 0, 1, 0, 2, 0, 0, 2],
            ....:   [0, 0, 0, 1, 0, 2, 0, 1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C, maximum_error_weight = 2)
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), 2)
            sage: c = C.random_element()
            sage: r = Chan(c)
            sage: c == D.decode_to_code(r)
            True
        """
        lookup_table = self.syndrome_table()
        s = self.code().parity_check_matrix() * r
        s.set_immutable()
        if s.is_zero():
            return r
        err = lookup_table[s]
        r_corr = copy(r)
        for i in range(self.code().length()):
            r_corr[i] = r[i] - err[i]
        return r_corr

    def maximum_error_weight(self):
        r"""
        Returns the maximal number of errors a received word can have
        and for which ``self`` is guaranteed to return a most likely codeword.

        Same as ``self.decoding_radius``.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.maximum_error_weight()
            4
        """
        return self._maximum_error_weight

    def decoding_radius(self):
        r"""
        Returns the maximal number of errors a received word can have
        and for which ``self`` is guaranteed to return a most likely codeword.

        EXAMPLES::

            sage: G = Matrix(GF(3), [[1,0,0,1,0,1,0,1,2],[0,1,0,2,2,0,1,1,0],[0,0,1,0,2,2,2,1,2]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.decoding_radius()
            4
        """
        return self._maximum_error_weight

    def syndrome_table(self):
        r"""
        Returns the syndrome lookup table of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.syndrome_table()
            {(0, 0, 0): (0, 0, 0, 0, 0, 0, 0),
             (1, 0, 0): (1, 0, 0, 0, 0, 0, 0),
             (0, 1, 0): (0, 1, 0, 0, 0, 0, 0),
             (1, 1, 0): (0, 0, 1, 0, 0, 0, 0),
             (0, 0, 1): (0, 0, 0, 1, 0, 0, 0),
             (1, 0, 1): (0, 0, 0, 0, 1, 0, 0),
             (0, 1, 1): (0, 0, 0, 0, 0, 1, 0),
             (1, 1, 1): (0, 0, 0, 0, 0, 0, 1)}
        """
        return self._lookup_table










class LinearCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Codes. This decoder will decode to the
    nearest codeword found.

    INPUT:

    - ``code`` -- A code associated to this decoder
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        super(LinearCodeNearestNeighborDecoder, self).__init__(code, code.ambient_space(), \
                code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Tests equality between LinearCodeNearestNeighborDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: D1 = codes.decoders.LinearCodeNearestNeighborDecoder(LinearCode(G))
            sage: D2 = codes.decoders.LinearCodeNearestNeighborDecoder(LinearCode(G))
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for Linear code of length 7, dimension 4 over Finite Field of size 2
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: latex(D)
            \textnormal{Nearest neighbor decoder for }[7, 4]\textnormal{ Linear code over }\Bold{F}_{2}
        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: D.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)
        """
        c_min = self.code().zero()
        h_min = r.hamming_weight()
        for c in self.code():
            if (c-r).hamming_weight() < h_min:
                h_min = (c-r).hamming_weight()
                c_min = c
        c_min.set_immutable()
        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeNearestNeighborDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return (self.code().minimum_distance()-1) // 2

####################### registration ###############################

LinearCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder

LinearCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
LinearCodeSyndromeDecoder._decoder_type = {"hard-decision", "unique", "dynamic"}
LinearCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
LinearCodeNearestNeighborDecoder._decoder_type = {"hard-decision", "unique", "always-succeed", "complete"}
