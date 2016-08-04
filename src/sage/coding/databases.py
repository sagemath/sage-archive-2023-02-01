# -*- coding: utf-8 -*-
r"""
Databases and accessors of online databases for coding theory
"""


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
def bounds_minimum_distance(n, k, F):
    r"""
    Explains the construction of the best known linear code over GF(q) with
    length n and dimension k, courtesy of the www page
    http://www.codetables.de/.

    INPUT:

    -  ``n`` - Integer, the length of the code
    Calculates a lower and upper bound for the minimum distance of an optimal
    linear code with word length ``n`` and dimension ``k`` over the field
    ``F``.

    -  ``k`` - Integer, the dimension of the code
    The function returns a record with the two bounds and an explanation for
    each bound. The function Display can be used to show the explanations.

    -  ``F`` - Finite field, of order 2, 3, 4, 5, 7, 8, or 9
    The values for the lower and upper bound are obtained from a table
    constructed by Cen Tjhai for GUAVA, derived from the table of
    Brouwer. See http://www.codetables.de/ for the most recent data.
    These tables contain lower and upper bounds for `q=2` (when ``n <= 257``),
    `q=3` (when ``n <= 243``), `q=4` (``n <= 256``). (Current as of
    11 May 2006.) For codes over other fields and for larger word lengths,
    trivial bounds are used.

    -  ``verbose`` - Bool (default: ``False``)
    This does not require an internet connection. The format of the output is
    a little non-intuitive. Try ``bounds_minimum_distance(10,5,GF(2))`` for
    an example.

    OUTPUT:


    -  Text about why the bounds are as given
    This function requires optional GAP package (Guava).

    EXAMPLES::

        sage: L = best_known_linear_code_www(72, 36, GF(2)) # optional - internet
        sage: print(L)                                      # optional - internet
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
        sage: print(bounds_minimum_distance(10,5,GF(2))) # optional - gap_packages (Guava package)
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
    if not q in [2, 3, 4, 5, 7, 8, 9]:
        raise ValueError("q (=%s) must be in [2,3,4,5,7,8,9]"%q)
    n = int(n)
    k = int(k)

    param = ("?q=%s&n=%s&k=%s"%(q,n,k)).replace('L','')

    url = "http://iaks-www.ira.uka.de/home/grassl/codetables/BKLC/BKLC.php"+param
    if verbose:
        print("Looking up the bounds at %s" % url)
    f = urlopen(url)
    s = f.read()
    f.close()

    i = s.find("<PRE>")
    j = s.find("</PRE>")
    if i == -1 or j == -1:
        raise IOError("Error parsing data (missing pre tags).")
    text = s[i+5:j].strip()
    return text
    gap.eval("data := BoundsMinimumDistance(%s,%s,GF(%s))"%(n,k,q))
    Ldata = gap.eval("Display(data)")
    return Ldata
