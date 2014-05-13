r"""
Database of small combinatorial designs

This module implements known constructions combinatorial designs, and in
particular of :mod:`Orthogonal Arrays
<sage.combinat.designs.orthogonal_arrays>`. Most of them come from the chapter
on :mod:`Mutually Orthogonal Latin Squares
<sage.combinat.designs.latin_squares>` from the Handbook of Combinatorial
Designs.

All the designs returned by these functions can be obtained through the
``designs.<tab>`` functions.

Implemented constructions :

- :func:`OA(6,20) <OA_6_20>`,
  :func:`OA(7,21) <OA_7_21>`,
  :func:`OA(5,22) <OA_5_22>`,
  :func:`OA(9,24) <OA_9_24>`,
  :func:`OA(6,26) <OA_6_26>`,
  :func:`OA(7,28) <OA_7_28>`,
  :func:`OA(6,30) <OA_6_30>`,
  :func:`OA(7,33) <OA_7_33>`,
  :func:`OA(6,34) <OA_6_34>`,
  :func:`OA(7,35) <OA_7_35>`,
  :func:`OA(10,36) <OA_10_36>`,
  :func:`OA(6,38) <OA_6_38>`,
  :func:`OA(7,39) <OA_7_39>`,
  :func:`OA(9,40) <OA_9_40>`,
  :func:`OA(7,42) <OA_7_42>`,
  :func:`OA(7,44) <OA_7_44>`,
  :func:`OA(8,45) <OA_8_45>`,
  :func:`OA(6,46) <OA_6_46>`,
  :func:`OA(10,48) <OA_10_48>`,
  :func:`OA(8,50) <OA_8_50>`,
  :func:`OA(7,51) <OA_7_51>`,
  :func:`OA(7,52) <OA_7_52>`,
  :func:`OA(7,54) <OA_7_54>`,
  :func:`OA(8,55) <OA_8_55>`,
  :func:`OA(9,56) <OA_9_56>`,
  :func:`OA(7,60) <OA_7_60>`,
  :func:`OA(7,62) <OA_7_62>`,
  :func:`OA(9,75) <OA_9_75>`,
  :func:`OA(11,80) <OA_11_80>`,
  :func:`OA(10,82) <OA_10_82>`,
  :func:`OA(10,100) <OA_10_100>`,
  :func:`OA(12,144) <OA_12_144>`,
  :func:`OA(12,210) <OA_12_210>`

- :func:`two MOLS of order 10 <MOLS_10_2>`,
  :func:`five MOLS of order 12 <MOLS_12_5>`,
  :func:`four MOLS of order 14 <MOLS_14_4>`,
  :func:`four MOLS of order 15 <MOLS_15_4>`,
  :func:`three MOLS of order 18 <MOLS_18_3>`


**Dictionaries**

The functions defined here are used by
:func:`~sage.combinat.designs.orthogonal_arrays.orthogonal_array`. Thus, the
functions are indexed by dictionary which associates to every integer ``n`` a
pair ``(k,f)`` where ``f`` is a function such that ``f()`` is a `OA(k,n)`. This
dictionary is defined right after the constructions of OA in the file.

The same goes for the constructions of MOLS, used by
:func:`~sage.combinat.designs.latin_squares.mutually_orthogonal_latin_squares`.

REFERENCES:

.. [DesignHandbook] Handbook of Combinatorial Designs (2ed)
  Charles Colbourn, Jeffrey Dinitz
  Chapman & Hall/CRC
  2012

Functions
---------
"""

from sage.combinat.designs.orthogonal_arrays import (OA_from_quasi_difference_matrix,
                                                     OA_from_Vmt,
                                                     OA_from_wider_OA)

# Cyclic shift of a list
cyclic_shift = lambda l,i : l[-i:]+l[:-i]

def TD_6_12():
    r"""
    Returns a `TD(6,12)` as build in [Hanani75]_.

    This design is Lemma 3.21 from [Hanani75]_.

    EXAMPLE::

        sage: from sage.combinat.designs.database import TD_6_12
        sage: from sage.combinat.designs.orthogonal_arrays import is_transversal_design
        sage: TD = TD_6_12()
        sage: is_transversal_design(TD,6,12)
        True

    The design is available from the general constructor::

        sage: designs.transversal_design(6,12,existence=True)
        True

    REFERENCES:

    .. [Hanani75] Haim Hanani,
      Balanced incomplete block designs and related designs,
      http://dx.doi.org/10.1016/0012-365X(75)90040-0,
      Discrete Mathematics, Volume 11, Issue 3, 1975, Pages 255-369.
    """
    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    G = AdditiveAbelianGroup([2,6])
    d = [[(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)],
         [(0,0),(0,1),(1,0),(0,3),(1,2),(0,4)],
         [(0,0),(0,2),(1,2),(1,0),(0,1),(1,5)],
         [(0,0),(0,3),(0,2),(0,1),(1,5),(1,4)],
         [(0,0),(0,4),(1,1),(1,3),(0,5),(0,2)],
         [(0,0),(0,5),(0,1),(1,5),(1,3),(1,1)],
         [(0,0),(1,0),(1,3),(0,2),(0,3),(1,2)],
         [(0,0),(1,1),(1,5),(1,2),(1,4),(1,0)],
         [(0,0),(1,2),(0,4),(0,5),(0,2),(1,3)],
         [(0,0),(1,3),(1,4),(0,4),(1,1),(0,1)],
         [(0,0),(1,4),(0,5),(1,1),(1,0),(0,3)],
         [(0,0),(1,5),(0,3),(1,4),(0,4),(0,5)]]

    r = lambda x : int(x[0])*6+int(x[1])
    TD = [[i*12+r(G(x)+g) for i,x in enumerate(X)] for X in d for g in G]
    for x in TD: x.sort()

    return TD

def _MOLS_from_string(s,k):
    r"""
    Returns MOLS from a string

    INPUT:

    - ``s`` (string) -- represents the MOLS with entries in a-z. To understand
      how the string should be formatted, read the source code of a constructor
      that uses it.

    - ``k`` (integer) -- the number of MOLS encoded by the string.

    EXAMPLES::

        sage: _ = designs.mutually_orthogonal_latin_squares(10,2) # indirect doctest
    """
    from sage.matrix.constructor import Matrix
    matrices = [[] for _ in range(k)]
    for i,l in enumerate(s.split()):
        l = [ord(x) - 97 for x in l]
        matrices[i%k].append(l)
    return map(Matrix, matrices)

def MOLS_10_2():
    r"""
    Returns a pair of MOLS of order 10

    Data obtained from
    `<http://www.cecm.sfu.ca/organics/papers/lam/paper/html/POLS10/POLS10.html>`_

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: from sage.combinat.designs.database import MOLS_10_2
        sage: MOLS = MOLS_10_2()
        sage: print are_mutually_orthogonal_latin_squares(MOLS)
        True

    The design is available from the general constructor::

        sage: designs.mutually_orthogonal_latin_squares(10,2,existence=True)
        True
    """
    from sage.matrix.constructor import Matrix
    return [Matrix([[1,8,9,0,2,4,6,3,5,7],
                    [7,2,8,9,0,3,5,4,6,1],
                    [6,1,3,8,9,0,4,5,7,2],
                    [5,7,2,4,8,9,0,6,1,3],
                    [0,6,1,3,5,8,9,7,2,4],
                    [9,0,7,2,4,6,8,1,3,5],
                    [8,9,0,1,3,5,7,2,4,6],
                    [2,3,4,5,6,7,1,8,9,0],
                    [3,4,5,6,7,1,2,0,8,9],
                    [4,5,6,7,1,2,3,9,0,8]]),

            Matrix([[1,7,6,5,0,9,8,2,3,4],
                    [8,2,1,7,6,0,9,3,4,5],
                    [9,8,3,2,1,7,0,4,5,6],
                    [0,9,8,4,3,2,1,5,6,7],
                    [2,0,9,8,5,4,3,6,7,1],
                    [4,3,0,9,8,6,5,7,1,2],
                    [6,5,4,0,9,8,7,1,2,3],
                    [3,4,5,6,7,1,2,8,0,9],
                    [5,6,7,1,2,3,4,0,9,8],
                    [7,1,2,3,4,5,6,9,8,0]])]

def MOLS_12_5():
    r"""
    Returns 5 MOLS of order 12

    These MOLS have been found by Brendan McKay.

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: from sage.combinat.designs.database import MOLS_12_5
        sage: MOLS = MOLS_12_5()
        sage: print are_mutually_orthogonal_latin_squares(MOLS)
        True
    """
    M = """
        abcdefghijkl abcdefghijkl abcdefghijkl abcdefghijkl abcdefghijkl
        badcfehgjilk ghefklijcdab dcbahgfelkji jilkbadcfehg klijcdabghef
        cdabghefklij efghijklabcd lkjidcbahgfe ijklabcdefgh fehgjilkbadc
        dcbahgfelkji cdabghefklij ghefklijcdab badcfehgjilk hgfelkjidcba
        ijklabcdefgh klijcdabghef efghijklabcd fehgjilkbadc jilkbadcfehg
        jilkbadcfehg fehgjilkbadc hgfelkjidcba dcbahgfelkji lkjidcbahgfe
        klijcdabghef hgfelkjidcba jilkbadcfehg cdabghefklij dcbahgfelkji
        lkjidcbahgfe ijklabcdefgh badcfehgjilk efghijklabcd ghefklijcdab
        efghijklabcd jilkbadcfehg fehgjilkbadc lkjidcbahgfe cdabghefklij
        fehgjilkbadc dcbahgfelkji cdabghefklij ghefklijcdab badcfehgjilk
        ghefklijcdab badcfehgjilk klijcdabghef hgfelkjidcba ijklabcdefgh
        hgfelkjidcba lkjidcbahgfe ijklabcdefgh klijcdabghef efghijklabcd
        """

    return _MOLS_from_string(M,5)

def MOLS_14_4():
    r"""
    Returns four MOLS of order 14

    These MOLS were shared by Ian Wanless.

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: from sage.combinat.designs.database import MOLS_14_4
        sage: MOLS = MOLS_14_4()
        sage: print are_mutually_orthogonal_latin_squares(MOLS)
        True

    The design is available from the general constructor::

        sage: designs.mutually_orthogonal_latin_squares(14,4,existence=True)
        True
    """
    M = """
        bjihgkecalnfmd  bfmcenidgjhalk  bcdefghijklmna  bcdefghijklmna
        fckjbhledimagn  jcgndfalehkbim  gnkjdmiclbhaef  jflhnkaecmgdib
        mgdlkcbafejnih  ikdhaegnmfblcj  lifhbjemkangcd  emkdjbgfnliahc
        cnhemldbigfkaj  hjlebifkangcmd  dalmgnbjehcfik  anighmflkbdcej
        edabfnmkcjhgli  gbkmfcjeliahdn  njcaeifhbdgkml  kebcajimdgfhln
        nfeicgajldkbhm  khclngdafmjibe  mfbkcdlagnjihe  cgnflembihakjd
        iagfjdhnkmelcb  elbdmahfignkjc  aemnhkjdcifblg  ilabkdnhfcjegm
        dlnkeafimhcjbg  ceabkjnihdmgfl  hdnikbagmcelfj  ljgnihecbamfdk
        gemalfihjnbdkc  adficlkmjbenhg  cgjflhnbiekdam  ndmabcjglfeikh
        jhfnimgdbkacel  liegjdmhnkcfab  fkibmagenldhjc  mbhiefljadkncg
        hkbgajnmeclidf  nmjfhkecbaldgi  imhlneckdfajgb  difjcnkamehgbl
        ablchikgnfdmje  fankgbljdcimeh  klegafdnhjmcbi  ghckmlbdeinjaf
        licmdbjfhagenk  mgialhcbkedjnf  jhadicmlfgbekn  fajlgidkhncbme
        kmjdneclgbihfa  dnhjimbgclfeka  ebgcjlkfamindh  hkemdacngjblfi
        """

    return _MOLS_from_string(M,4)

def MOLS_15_4():
    r"""
    Returns 4 MOLS of order 15.

    These MOLS were shared by Ian Wanless.

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: from sage.combinat.designs.database import MOLS_15_4
        sage: MOLS = MOLS_15_4()
        sage: print are_mutually_orthogonal_latin_squares(MOLS)
        True

    The design is available from the general constructor::

        sage: designs.mutually_orthogonal_latin_squares(15,4,existence=True)
        True
    """
    M = """
        bcdefghijklmnoa  bdgiknfcamehjlo  bhealiofmdjgcnk  blhcmdinejofakg
        abcdefghijklmno  acehjlogdbnfikm  lcifbmjagnekhdo  hcmidnejofkagbl
        oabcdefghijklmn  nbdfikmahecogjl  amdjgcnkbhoflie  midnjeofkaglbhc
        noabcdefghijklm  mocegjlnbifdahk  fbnekhdolciagmj  dnjeokfaglbhmci
        mnoabcdefghijkl  lnadfhkmocjgebi  kgcoflieamdjbhn  jeokfalgbhmcind
        lmnoabcdefghijk  jmobegilnadkhfc  olhdagmjfbnekci  ekfalgbmhcindjo
        klmnoabcdefghij  dknacfhjmobelig  jamiebhnkgcofld  aflgbmhcnidjoek
        jklmnoabcdefghi  helobdgiknacfmj  ekbnjfciolhdagm  lbgmhcnidojekaf
        ijklmnoabcdefgh  kifmacehjlobdgn  nflcokgdjamiebh  gmchnidojeakflb
        hijklmnoabcdefg  oljgnbdfikmaceh  iogmdalhekbnjfc  chndiojeakfblgm
        ghijklmnoabcdef  iamkhocegjlnbdf  djahnebmiflcokg  ndioejakfblgcmh
        fghijklmnoabcde  gjbnliadfhkmoce  hekbiofcnjgmdal  ioejafkblgcmhdn
        efghijklmnoabcd  fhkcomjbegilnad  miflcjagdokhneb  ojafkbglcmhdnie
        defghijklmnoabc  egildankcfhjmob  cnjgmdkbhealiof  fakbglchmdnieoj
        cdefghijklmnoab  cfhjmeboldgikna  gdokhnelcifbmja  kgblchmdineojfa
        """

    return _MOLS_from_string(M,4)

def MOLS_18_3():
    r"""
    Returns 3 MOLS of order 18.

    These MOLS were shared by Ian Wanless.

    EXAMPLES::

        sage: from sage.combinat.designs.latin_squares import are_mutually_orthogonal_latin_squares
        sage: from sage.combinat.designs.database import MOLS_18_3
        sage: MOLS = MOLS_18_3()
        sage: print are_mutually_orthogonal_latin_squares(MOLS)
        True

    The design is available from the general constructor::

        sage: designs.mutually_orthogonal_latin_squares(18,3,existence=True)
        True
    """
    M = """
        bgejhkmodcnarilpfq  beqpodgcflkrjahnim  bcdefghijklmnopqra
        echfbilnprdokajmqg  gcfrqpehdnmlabkioj  rbkamfgdehqjinopcl
        qfdigcjmohaeplkbnr  ehdgarqfibonmkcljp  mlbqhifgajdcrenopk
        prgejhdbnaikfqmlco  jfiehkargqcponldmb  hijbcdefgqraklmnop
        oqahfbiecpkjlgrnmd  hbgjfilkacrdqpomen  gderbkamfpclhqjino
        dprkigcjfeqlbmhaon  kichbgjmlodaerqpnf  fgamlbqhiopkjdcren
        geqaljhdbofrmcnikp  mljdichbngpekfarqo  efghijbcdnopqraklm
        chfrkmbieqpgandojl  onmbejdicphqflgkar  amfgderbkinopclhqj
        fdigalncjmrqhkoepb  dponcfbejaqirgmhlk  qhifgamlbrenopkjdc
        lnbqcogpakmhdrifej  crhapeoqmkbgidnjfl  leqjponacbkidfgmhr
        kmocrdphqblnieajgf  ndaikqfprmlchjeobg  kjmcrponhlbqeafgid
        rlnpdaeqigcmojfkbh  aoekjlrgqhnmdibfpc  dqriklponajbcmhfge
        jamoqekfrihdnpbglc  rkpflbmahdionejcgq  nacleqjpomhrbkidfg
        abknprflgdjieoqchm  ialqgmcnkrejpofbdh  onhkjmcrpgidlbqeaf
        hkcloqagmnebjfprdi  ljkmrhndoiafbqpgce  pondqriklfgeajbcmh
        nildmprkhjofcbgqae  pmblnaioefjkgcrqhd  jponacleqdfgmhrbki
        iojmenqalfbpgdchrk  fqncmokjpegblhdari  crponhkjmeafgidlbq
        mjpbnforklgcqhedia  qgrodnplbjfhcmieka  iklpondqrcmhfgeajb
        """

    return _MOLS_from_string(M,3)

# Index of the MOLS constructions
#
# Associates to n the pair (k,f) where f() is a function that returns k MOLS of order n
#
# This dictionary is used by designs.mutually_orthogonal_latin_squares(n,k).

MOLS_constructions = {
    10 : (2, MOLS_10_2),
    12 : (5, MOLS_12_5),
    14 : (4, MOLS_14_4),
    15 : (4, MOLS_15_4),
    18 : (3, MOLS_18_3)
}

def OA_6_20():
    r"""
    Returns an OA(6,20)

    As explained in the Handbook III.3.49 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_20
        sage: OA = OA_6_20()
        sage: print is_orthogonal_array(OA,6,20,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,20,existence=True)
        True
    """
    M=[[None,   7,  13,   1,  16,   9,   2],
       [   0,   1,  15,   7,  17,   6,  14],
       [   0,  11,  10,  11,   5,   4,   3],
       [   7,None,  13,  16,   1,   2,   9],
       [   1,   0,  15,  17,   7,  14,   6],
       [  11,   0,  10,   5,  11,   3,   4]]

    Mb=[[],[],[],[],[],[]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        Mb[0].extend([a,b,c])
        Mb[1].extend([b,c,a])
        Mb[2].extend([c,a,b])
        Mb[3].extend([d,f,e])
        Mb[4].extend([e,d,f])
        Mb[5].extend([f,e,d])

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    M = OA_from_quasi_difference_matrix(Mb,AdditiveCyclic(19),add_col=False)

    return M

def OA_7_21():
    r"""
    Returns an OA(7,21)

    As explained in the Handbook III.3.50 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_21
        sage: OA = OA_7_21()
        sage: print is_orthogonal_array(OA,7,21,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,21,existence=True)
        True
    """
    M = [[  8,  17,  20,   2],
         [  9,  16,   4,  15],
         [ 11,   5,  10,   6],
         [ 14,   1,   3,  13],
         [ 18,  19,  12,   7]]

    Mb = [[0],[0],[0],[0],[0],[0]]
    for a,b,c,d,e in zip(*M):
        Mb[0].extend([a,b,c,d,e])
        Mb[1].extend([b,c,d,e,a])
        Mb[2].extend([c,d,e,a,b])
        Mb[3].extend([d,e,a,b,c])
        Mb[4].extend([e,a,b,c,d])
        Mb[5].extend([0,0,0,0,0])

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    M = OA_from_quasi_difference_matrix(Mb,AdditiveCyclic(21))
    return M

def OA_5_22():
    r"""
    Returns an OA(5,22)

    As explained in the Handbook III.3.51 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_5_22
        sage: OA = OA_5_22()
        sage: print is_orthogonal_array(OA,5,22,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(5,22,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(21)
    M = [
        [   1,  13,  18,   3,  16,  19,None],
        [  16,  19,   1,  13,  18,   3,   0],
        [  18,   3,  16,  19,   1,  13,   0],
        [   6,  15,   6,  15,   6,  15,   0],
        [  12,   9,  19,  16,   5,   2,   0],
        ]

    Mb=[[],[],[],[],[]]

    for R in zip(*M):
        a,b,c,d,e = [G(x) if x is not None else None for x in R]
        Mb[0].extend([a,16*c,4*b])
        Mb[1].extend([b,None if a is None else 16*a,4*c])
        Mb[2].extend([c,16*b,None if a is None else 4*a])
        Mb[3].extend([d,16*d+7,4*d+14])
        Mb[4].extend([e,16*e+14,4*e+7])

    Mb[0].extend([0,0])
    Mb[1].extend([7,14])
    Mb[2].extend([14,7])
    Mb[3].extend([None,0])
    Mb[4].extend([0,None])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col=False)
    return M

def OA_9_24():
    r"""
    Returns an OA(9,24)

    As explained in the Handbook III.3.52 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_9_24
        sage: OA = OA_9_24()
        sage: print is_orthogonal_array(OA,9,24,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(9,24,existence=True)
        True
    """
    M = ("0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 0000 "+
         "0000 0010 0100 0110 1000 1010 1100 1110 2000 2010 2100 2110 "+
         "0000 0011 1001 2110 0111 2011 2111 1000 0100 1100 1101 2010 "+
         "0000 1010 1011 2000 1101 2110 0001 0101 2100 2001 0111 1100 "+
         "0000 0001 2010 1111 2111 2100 1101 0011 1010 2101 1000 0110 "+
         "0000 1000 2001 1011 0100 1100 0110 2101 2111 0010 1111 2011 "+
         "0000 1001 0111 2100 2000 0010 1110 2011 1100 1011 0101 2111 "+
         "0000 1011 2101 0100 2110 1001 2000 0110 0101 1111 2011 1010 ")

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    G = AdditiveAbelianGroup([3,2,2,2])
    rlabel = {(x%2,x%3):x for x in range(6)}
    M = [G([int(c),int(d),rlabel[int(b),int(a)]]) for a,b,c,d in M.split()]
    M = [M[i*12:(i+1)*12] for i in range(8)]
    Mb = [[] for _ in range(8)]
    for a,b,c,d,e,f,g,h in zip(*M):
        Mb[0].extend([a, a + G([0,0,rlabel[0,0]])])
        Mb[1].extend([b, b + G([0,1,rlabel[0,0]])])
        Mb[2].extend([c, c + G([1,0,rlabel[0,0]])])
        Mb[3].extend([d, d + G([1,1,rlabel[0,0]])])
        Mb[4].extend([e, e + G([0,0,rlabel[1,0]])])
        Mb[5].extend([f, f + G([0,1,rlabel[1,0]])])
        Mb[6].extend([g, g + G([1,0,rlabel[1,0]])])
        Mb[7].extend([h, h + G([1,1,rlabel[1,0]])])

    M = OA_from_quasi_difference_matrix(Mb,G)
    return M

def OA_6_26():
    r"""
    Returns an OA(6,26)

    As explained in the Handbook III.3.53 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_26
        sage: OA = OA_6_26()
        sage: print is_orthogonal_array(OA,6,26,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,26,existence=True)
        True
    """
    M = [
        [None,None,None,None,None],
        [   0,   0,   0,   0,   0],
        [   1,   6,   7,   8,  14],
        [   3,  11,  20,  18,  10],
        [   6,  10,  14,   1,   5],
        [   4,  19,   5,  12,   2],
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(21)
    Mb=[[0],[0],[0],[0],[0],[0]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        Mb[0].extend([a,b,c,d,e,f])
        Mb[1].extend([b,c,d,e,f,a])
        Mb[2].extend([c,d,e,f,a,b])
        Mb[3].extend([d,e,f,a,b,c])
        Mb[4].extend([e,f,a,b,c,d])
        Mb[5].extend([f,a,b,c,d,e])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_7_28():
    r"""
    Returns an OA(7,28)

    As explained in the Handbook III.3.54 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_28
        sage: OA = OA_7_28()
        sage: print is_orthogonal_array(OA,7,28,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,28,existence=True)
        True
    """
    z=2
    M = [
        [(0,0), (z+1,6),(1,1)  ,(1,1)  ,(1,3)  ,(1,4)  ,(0,0)  ,(1,4), (z,5)  ],
        [(z,2), (0,0)  ,(1,5)  ,(z,1)  ,(z,2)  ,(z,6)  ,(z+1,3),(0,0), (z,1)  ],
        [(z,3), (z+1,4),(0,0)  ,(z+1,5),(z+1,2),(z+1,4),(z+1,2),(1,6), (0,0)  ],
        [(0,5), (z,6)  ,(0,5)  ,(0,6)  ,(z,3)  ,(0,0)  ,(0,4)  ,(1,5), (z+1,4)],
        [(0,3), (0,3)  ,(z+1,5),(0,0)  ,(0,5)  ,(z+1,6),(1,1)  ,(0,1), (z,3)  ],
        [(1,3), (0,6)  ,(0,6)  ,(1,5)  ,(0,0)  ,(0,3)  ,(z+1,6),(z,2), (0,2)  ],
        ]

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    from sage.modules.free_module_element import free_module_element as vector
    G = AdditiveAbelianGroup([2,2,7])
    M = [[G(vector([x//2,x%2,y])) for x,y in L] for L in M]

    Mb=[[0],[0],[0],[0],[0],[0]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        Mb[0].extend([a,b,c])
        Mb[1].extend([b,c,a])
        Mb[2].extend([c,a,b])
        Mb[3].extend([d,f,e])
        Mb[4].extend([e,d,f])
        Mb[5].extend([f,e,d])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_6_30():
    r"""
    Returns an OA(6,30)

    As explained in the Handbook III.3.55 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_30
        sage: OA = OA_6_30()
        sage: print is_orthogonal_array(OA,6,30,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,30,existence=True)
        True
    """
    M = [
        [(0,0),None,(0,0),(0,0),(0,0),(0,0),(0,0)],
        [(0,0),(0,0),None,(0,4),(0,2),(0,3),(0,1)],
        [(0,0),(3,1),(3,0),None,(4,0),(1,0),(2,0)],
        [(0,0),(3,0),(0,2),(1,2),None,(0,1),(0,3)],
        [(0,0),(3,3),(1,2),(4,2),(2,0),None,(0,4)],
        [(0,0),(4,2),(2,4),(0,3),(2,3),(3,2),None]
        ]

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    from sage.modules.free_module_element import free_module_element as vector
    G = AdditiveAbelianGroup([5,5])
    M = [[None if x is None else G(vector(x)) for x in L] for L in M]

    Mb=[[],[],[],[],[],[]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        for i in range(5):
            Mb[0].append(None if a is None else a+G(vector((i,i))))
            Mb[1].append(None if b is None else b+G(vector((2*i,i))))
            Mb[2].append(None if c is None else c+G(vector((i,0))))
            Mb[3].append(None if d is None else d+G(vector((4*i,0))))
            Mb[4].append(None if e is None else e+G(vector((3*i,4*i))))
            Mb[5].append(None if f is None else f+G(vector((4*i,4*i))))

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_7_33():
    r"""
    Returns an OA(7,33)

    As explained in the Handbook III.3.56 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_33
        sage: OA = OA_7_33()
        sage: print is_orthogonal_array(OA,7,33,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,33,existence=True)
        True
    """
    M = [
        [   0,   0,   0,   0,   0,   0],
        [  15,  11,  22,   4,  17,   8],
        [  19,   7,  14,  32,  22,  18],
        [  22,  19,   8,  24,  21,   6],
        [   9,  12,  15,   7,  26,  14],
        [  14,  28,  23,   2,  19,   3]
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(33)

    Mb=[[0,1,7],[0,4,28],[0,16,13],[0,31,19],[0,25,10],[0,22,0]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        for i in range(5):
            Mb[0].append(a)
            Mb[1].append(b)
            Mb[2].append(c)
            Mb[3].append(d)
            Mb[4].append(e)
            Mb[5].append(f)
            a,b,c,d,e,f = 4*e,4*a,4*b,4*c,4*d,4*f

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_6_34():
    r"""
    Returns an OA(6,34)

    As explained in the Handbook III.3.57 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_34
        sage: OA = OA_6_34()
        sage: print is_orthogonal_array(OA,6,34,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,34,existence=True)
        True
    """
    M = [
        [None,   0,   0,   0,   0,   0],
        [  30,  17,  10,  25,  23,   8],
        [  22,   4,  32,  29,  28,  22],
        [  25,  10,  20,  15,  21,  16],
        [   0,  12,  15,  16,  32,  23],
        [   6,  11,  18,  14,   9,  20]
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(33)

    Mb=[[0,1,3,10,5],[0,4,12,7,20],[0,16,15,28,14],[0,31,27,13,23],[0,25,9,19,26],[0,11,11,0,None]]

    times4 = lambda x : None if x is None else 4*x
    for R in zip(*M):
        a,b,c,d,e,f = [None if x is None else G(x) for x in R]
        for i in range(5):
            Mb[0].append(a)
            Mb[1].append(b)
            Mb[2].append(c)
            Mb[3].append(d)
            Mb[4].append(e)
            Mb[5].append(f)
            a,b,c,d,e,f = map(times4,[e,a,b,c,d,f])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_7_35():
    r"""
    Returns an OA(7,35)

    As explained in the Handbook III.3.58 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_35
        sage: OA = OA_7_35()
        sage: print is_orthogonal_array(OA,7,35,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,35,existence=True)
        True
    """
    M = [
        [  0, 15, 30, 10, 25,  1, 16, 31, 11, 26,  2, 17, 32, 12,  6,  3, 18, 33, 27, 21,  4, 19, 13,  7, 22,  5, 34, 28,  8, 23, 20, 14, 29,  9, 24],
        [  0, 22, 16,  3,  4,  9, 10, 32, 26, 13, 18,  5, 27, 14, 15, 20,  7,  1, 23, 31, 29,  2, 24, 11, 19, 17, 25, 12,  6, 28, 33, 34, 21,  8, 30],
        [  0, 29,  2, 31, 18, 10, 32, 26, 34, 28, 27, 21, 15,  9, 17, 30,  3,  4,  5, 20, 12,  6, 14, 22, 16,  8, 23, 24, 25, 33, 11, 19, 13,  7,  1],
        [  0,  8,  9, 17, 11, 25, 19, 27, 28,  1, 15, 23, 31,  4, 26, 12,  6, 14, 29, 16,  2,  3, 18, 33, 34, 20,  7, 22, 30, 24, 10, 32,  5, 13, 21],
        [  0,  1, 23, 24, 32, 33,  6,  7, 29, 30, 10, 11, 12, 13, 28,  8,  9, 31,  4,  5, 27, 14, 15, 16,  3, 25, 26, 34, 21, 22,  2, 17, 18, 19, 20],
        [0]*35
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(35)

    M = OA_from_quasi_difference_matrix(M,G,add_col = True)
    return M

def OA_10_36():
    r"""
    Returns an OA(10,36)

    As explained in the Handbook III.3.59 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_10_36
        sage: OA = OA_10_36()
        sage: print is_orthogonal_array(OA,10,36,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(10,36,existence=True)
        True
    """
    M = [
        [(0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0), (0,0,0,0)],
        [(0,0,0,0), (0,1,0,0), (1,0,0,0), (1,1,0,0), (0,0,0,1), (0,1,0,1), (1,0,0,1), (1,1,0,1), (0,0,0,2), (0,1,0,2), (1,0,0,2), (1,1,0,2)],
        [(0,0,0,0), (1,1,1,2), (0,0,2,1), (0,0,1,2), (0,1,2,0), (0,1,0,2), (1,1,1,1), (0,1,1,1), (1,1,1,0), (1,0,2,2), (1,0,0,1), (1,0,1,0)],
        [(0,0,0,0), (0,0,1,0), (1,0,1,0), (0,1,0,0), (1,1,0,0), (1,0,2,0), (1,0,0,0), (0,1,2,0), (1,1,2,0), (0,0,2,0), (1,1,1,0), (0,1,1,0)],
        [(0,0,0,0), (0,1,2,0), (0,0,1,0), (1,1,1,0), (1,0,2,0), (1,0,1,0), (0,1,0,0), (0,0,2,0), (0,1,1,0), (1,1,0,0), (1,1,2,0), (1,0,0,0)],
        [(0,0,0,0), (0,1,1,0), (0,1,2,0), (1,1,2,0), (1,1,0,2), (0,0,1,2), (1,1,2,2), (1,0,0,2), (1,0,0,1), (1,0,1,1), (0,0,2,1), (0,1,1,1)],
        [(0,0,0,0), (1,0,1,0), (1,1,0,1), (1,0,1,2), (1,0,2,2), (0,0,2,1), (0,1,0,1), (0,1,0,0), (1,1,2,2), (0,1,1,0), (0,0,1,2), (1,1,2,1)],
        [(0,0,0,0), (1,1,0,0), (0,1,1,0), (1,0,2,1), (0,1,0,2), (1,0,2,2), (0,0,2,2), (1,1,1,0), (1,0,1,1), (0,1,2,1), (1,1,1,1), (0,0,0,2)],
        [(0,0,0,0), (1,0,0,0), (1,1,1,0), (0,1,1,2), (1,1,2,1), (0,1,1,1), (0,0,1,1), (1,0,2,0), (0,1,2,2), (1,1,0,2), (1,0,2,2), (0,0,0,1)]
        ]

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    from sage.modules.free_module_element import free_module_element as vector
    G = AdditiveAbelianGroup([2,2,3,3])
    M = [[G(vector(x)) for x in L] for L in M]

    Mb=[[],[],[],[],[],[],[],[],[]]

    for R in zip(*M):
        a,b,c,d,e,f,g,h,i = R
        for y in range(3):
            Mb[0].append(a+G(vector([0,0,0,0])))
            Mb[1].append(b+G(vector([0,0,y,0])))
            Mb[2].append(c+G(vector([0,0,2*y,0])))
            Mb[3].append(d+G(vector([0,0,0,y])))
            Mb[4].append(e+G(vector([0,0,0,2*y])))
            Mb[5].append(f+G(vector([0,0,y,y])))
            Mb[6].append(g+G(vector([0,0,2*y,2*y])))
            Mb[7].append(h+G(vector([0,0,y,2*y])))
            Mb[8].append(i+G(vector([0,0,2*y,y])))

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_6_38():
    r"""
    Returns an OA(6,38)

    As explained in the Handbook III.3.60 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_38
        sage: OA = OA_6_38()
        sage: print is_orthogonal_array(OA,6,38,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,38,existence=True)
        True
    """
    M = [
        [None,  10,   1,   2,   6,   3,  22,   5,   7,   9,  14,  18,  28],
        [   0,   1,  10,  20,  23,  30,  35,  13,  33,  16,  29,  32,  21],
        [   0,  26,  26,  15,   8,   4,  17,  19,  34,  12,  31,  24,  25],
        [  10,None,  10,   6,   2,  22,   3,   7,   5,  14,   9,  28,  18],
        [   1,   0,  26,  23,  20,  35,  30,  33,  13,  29,  16,  21,  32],
        [  26,   0,   1,   8,  15,  17,   4,  34,  19,  31,  12,  25,  24]
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(37)

    Mb=[[],[],[],[],[],[]]

    for R in zip(*M):
        a,b,c,d,e,f = R
        Mb[0].extend([a,b,c])
        Mb[1].extend([b,c,a])
        Mb[2].extend([c,a,b])
        Mb[3].extend([d,f,e])
        Mb[4].extend([e,d,f])
        Mb[5].extend([f,e,d])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_7_39():
    r"""
    Returns an OA(7,39)

    As explained in the Handbook III.3.61 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_39
        sage: OA = OA_7_39()
        sage: print is_orthogonal_array(OA,7,39,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,39,existence=True)
        True
    """
    M = [
        [   0,   0,   0,   0,   0,   0],
        [   4,  23,  13,   5,  12,  11],
        [  25,  11,  22,  34,  23,   6],
        [  13,   4,  20,  17,  15,  29],
        [  27,  21,   8,  16,  19,  26],
        [  16,  19,  34,  38,  26,  21]
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(39)

    Mb=[[0,1,-1],[0,16,-16],[0,22,-22],[0,17,-17],[0,38,-38],[0,23,-23]]

    for R in zip(*M):
        a,b,c,d,e,f = [None if x is None else G(x) for x in R]
        for i in range(3):
            Mb[0].extend([a,-a])
            Mb[1].extend([b,-b])
            Mb[2].extend([c,-c])
            Mb[3].extend([d,-d])
            Mb[4].extend([e,-e])
            Mb[5].extend([f,-f])
            a,b,c,d,e,f = [16*x for x in [c,a,b,f,d,e]]

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_9_40():
    r"""
    Returns an OA(9,40)

    As explained in the Handbook III.3.62 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_9_40
        sage: OA = OA_9_40()
        sage: print is_orthogonal_array(OA,9,40,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(9,40,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField

    F8 = FiniteField(8,'x')
    F5 = FiniteField(5)
    F5F8 = F5.cartesian_product(F8)
    w = F8.primitive_element()
    assert w**3 == w+1

    A = [
        [(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None)],
        [(0,None),(1,None),   (2,2),   (3,2),   (4,2),(2,None),(3,None),(4,None),   (0,2),   (1,2)],
        [(0,None),   (2,5),   (4,5),   (1,2),   (3,6),   (3,4),   (0,0),   (2,1),   (4,1),   (1,6)],
        [(0,None),   (3,4),   (1,4),   (4,0),   (2,5),(3,None),   (1,0),   (4,1),   (2,2),   (0,3)],
        [(0,None),   (4,6),(3,None),   (2,3),   (1,4),   (2,1),(1,None),   (0,4),   (4,0),   (3,2)],
        [(0,None),   (1,2),   (4,6),   (4,4),   (1,0),   (0,6),   (2,3),   (3,6),   (3,5),   (2,5)],
        [(1,None),   (0,3),   (1,2),   (4,5),(4,None),   (2,3),   (0,0),   (2,2),   (3,0),(3,None)],
        [(4,None),   (1,3),   (0,0),   (1,1),   (4,0),   (3,1),   (2,5),(0,None),   (2,1),(3,None)]
        ]
    Y = [
        [None, 0, 1, 6, 5, 4, 3, 2],
        [None, 1, 2, 0, 6, 5, 4, 3],
        ]
    r = lambda x : F8(0) if x is None else w**x

    A = [[(F5(a),r(b)) for a,b in L] for L in A]
    Y = [[r(b) for b in L] for L in Y]

    def t(i,(x,y)):
        a,b = A[x][y]
        b = Y[i][x]
        return F5F8((F5(0),b))

    R = {w:(0,1,0),
         w**2:(1,0,0),
         w**3:(0,1,1),
         w**4:(1,1,0),
         w**5:(1,1,1),
         w**6:(1,0,1),
         w**7:(0,0,1),
         F8(0):(0,0,0)
         }
    Mb = [[] for _ in range(8)]
    for y in range(len(A[0])):
        for x in range(len(A)):
            t1,t2 = t(0,(x,y)), t(1,(x,y))
            e = F5F8(A[x][y])
            Mb[x].extend([e,e+t1,e+t2,e+t1+t2])

    M = OA_from_quasi_difference_matrix(Mb,F5F8,add_col = True)
    return M

def OA_7_42():
    r"""
    Returns an OA(7,42)

    As explained in the Handbook III.3.63 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_42
        sage: OA = OA_7_42()
        sage: print is_orthogonal_array(OA,7,42,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,42,existence=True)
        True
    """
    M = [
        [None,None,None,None,None,None,None],
        [   0,   0,   0,   0,   0,   0,   0],
        [  18, -18,  11, -11,   5,  -5,   4],
        [  26, -26,  10, -10,  30, -30,  23],
        [  20, -20,   3,  -3,  33, -33,  23],
        [   5,  -5,  25, -25,  24, -24,   4],
        [  17, -17,   4,  -4,  22, -22,   0]
        ]

    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(35)

    Mb=[[],[],[],[],[],[],[]]

    for R in zip(*M):
        for i in range(7):
            Mb[i].extend(cyclic_shift(R,i))

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_7_44():
    r"""
    Returns an OA(7,44)

    As explained in the Handbook III.3.64 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_44
        sage: OA = OA_7_44()
        sage: print is_orthogonal_array(OA,7,44,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,44,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    from sage.categories.cartesian_product import cartesian_product

    G2 = AdditiveCyclic(2)
    G11 = AdditiveCyclic(11)
    G2211 = cartesian_product((G2,G2,G11))

    M = [
        [(0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0)],
        [(1,1,4), (0,1,4), (1,1,7), (1,0,6), (1,1,9), (0,1,2), (0,1,5), (0,1,1)],
        [(1,0,6), (0,1,3), (1,0,0), (0,1,9), (1,1,1), (0,1,4), (1,1,9), (1,0,9)],
        [(1,1,6), (1,1,9), (0,1,2), (1,1,0), (0,1,0), (1,1,5), (0,0,4), (0,0,9)],
        [(1,0,9), (0,0,2), (0,0,1), (1,0,2), (0,0,7), (1,1,6), (1,1,0), (1,0,7)],
        [(1,0,1), (1,0,6), (1,1,3), (0,1,5), (0,0,5), (0,1,3), (0,1,0), (1,1,0)]
        ]

    M = [[G2211(x) for x in L] for L in M]

    Mb=[[],[],[],[],[],[]]

    for R in zip(*M):
        for c in range(5):
            (x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4),(x5,y5,z5),(x6,y6,z6) = R
            for i,e in enumerate(R):
                Mb[i].append(e)
            R = [(x5,y5,5*z5),
                 (x1,y1,5*z1),
                 (x2,y2,5*z2),
                 (x3,y3,5*z3),
                 (x4,y4,5*z4),
                 (x6,y6,5*z6)]

    for x,y,z in [(0,0,0), (1,0,1),(1,1,2),(0,0,8)]:
        Mb[0].append((x,y,z))
        Mb[1].append((x,y,5*z))
        Mb[2].append((x,y,3*z))
        Mb[3].append((x,y,4*z))
        Mb[4].append((x,y,9*z))
        Mb[5].append((0,0,0))

    M = OA_from_quasi_difference_matrix(Mb,G2211,add_col = True)
    return M

def OA_8_45():
    r"""
    Returns an OA(8,45)

    As explained in the Handbook III.3.65 [DesignHandbook]_.

    ... whose description contained a very deadly typo, kindly fixed by Julian
    R. Abel.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_8_45
        sage: OA = OA_8_45()
        sage: print is_orthogonal_array(OA,8,45,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(8,45,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.categories.cartesian_product import cartesian_product

    G533 = cartesian_product((FiniteField(5),FiniteField(3),FiniteField(3)))

    M = [
        [(0,0,0), (2,2,1), (3,1,1), (4,1,2), (4,0,1), (0,1,1), (0,2,1), (3,2,2)],
        [(0,0,0), (1,2,1), (4,2,2), (1,2,0), (4,1,0), (3,1,1), (3,0,0), (2,1,2)],
        [(0,0,0), (4,1,1), (2,2,1), (3,2,0), (1,2,0), (2,1,0), (1,0,0), (3,2,1)],
        [(0,0,0), (0,1,0), (2,1,1), (4,0,0), (0,0,2), (4,2,2), (3,2,2), (1,2,2)],
        [(0,0,0), (3,1,2), (2,1,0), (0,2,2), (4,2,1), (0,2,1), (2,0,1), (1,1,2)],
        [(0,0,0), (2,1,1), (1,2,2), (3,0,1), (2,0,1), (1,0,0), (4,2,1), (1,1,0)],
        [(0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0), (0,0,0)]
        ]

    for i in range(6):
        M[i].extend(M[5-i][1:8])

    M[6].extend(M[6][1:8])

    Mb=[[],[],[],[],[],[],[]]

    for R in zip(*M):
        (x1,y1,z1),(x2,y2,z2),(x3,y3,z3),(x4,y4,z4),(x5,y5,z5),(x6,y6,z6),(x7,y7,z7) = R
        for i in range(3):
            Mb[0].append((x1, y1    , z1+i  ))
            Mb[1].append((x2, y2+2*i, z2    ))
            Mb[2].append((x3, y3+i  , z3+2*i))
            Mb[3].append((x4, y4+2*i, z4+i  ))
            Mb[4].append((x5, y5+i  , z5    ))
            Mb[5].append((x6, y6    , z6+2*i))
            Mb[6].append((x7, y7    , z7    ))

    M = OA_from_quasi_difference_matrix(Mb,G533,add_col = True)
    return M

def OA_6_46():
    r"""
    Returns an OA(6,46)

    As explained in the Handbook III.3.66 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_6_46
        sage: OA = OA_6_46()
        sage: print is_orthogonal_array(OA,6,46,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(6,46,existence=True)
        True
    """
    M = OA_from_Vmt(4,9,[0, 1, 3, 2, 8])
    return M

def OA_10_48():
    r"""
    Returns an OA(10,48)

    As explained in the Handbook III.3.67 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_10_48
        sage: OA = OA_10_48()
        sage: print is_orthogonal_array(OA,10,48,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(10,48,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    F16 = FiniteField(16,'x')
    F3 = FiniteField(3)
    F3F16 = F3.cartesian_product(F16)
    w = F16.primitive_element()
    assert w**4 == w+1

    A = [
        [ (0, 4), (2, 2),  (2,2), (0,13),  (0,4), (2,13),  (0,1),  (0,7), (1,7) , (2,2) ,  (0,6),  (2,9)],
        [ (2, 7), (0, 9),  (2,7), (2,3) ,  (0,3), (0,9) , (1,12),  (0,6), (0,12), (2,14),  (2,7), (0,11)],
        [ (2,12), (2,12), (0,14), (0,14),  (2,8), (0,8) ,  (0,2),  (1,2), (0,11), (0,1) ,  (2,4), (2,12)],
        [ (1, 3), (0, 2), (0,10), (0,14),  (0,9), (1,3) , (0,12), (2,13), (2,1) , (2,9) ,  (2,0),  (1,7)],
        [ (0, 0), (1, 8),  (0,7), (1,8) ,  (0,4), (0,14),  (2,6),  (0,2), (2,3) , (1,12), (2,14),  (2,5)],
        [ (0,12), (0, 5), (1,13), (0,4) , (1,13), (0,9) ,  (2,8), (2,11), (0,7) , (2,10),  (1,2),  (2,4)],
        [ (1,12), (2, 0), (1,14), (0,6) ,  (1,9), (0,14),  (1,4),  (0,5), (1,8) , (1,3) ,  (2,1),  (1,1)],
        [ (1, 4), (1, 2),  (2,5), (0,4) , (0,11), (1,14), (1,13),  (1,9), (0,10), (1,6) ,  (1,8),  (2,6)],
        [ (2,10), (1, 9),  (1,7), (1,4) ,  (0,9), (0,1) ,  (0,0),  (1,3), (1,14), (2,11), (1,11), (1,13)],
        ]

    A = [[F3F16((F3(a),w**b)) for a,b in L] for L in A]

    Mb = [[] for _ in range(9)]
    for L in zip(*A):
        for i,e in enumerate(L):
            Mb[i].append(e)

        for u in [0,1,4]:
            V = [12,2,7,0,5,10,3,8,13]
            for i,(e,x) in enumerate(zip(L,V)):
                Mb[i].append(e+F3F16((F3(0),w**(x+u))))

    M = OA_from_quasi_difference_matrix(Mb,F3F16,add_col = True)
    return M

def OA_8_50():
    r"""
    Returns an OA(8,50)

    As explained in the Handbook III.3.68 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_8_50
        sage: OA = OA_8_50()
        sage: print is_orthogonal_array(OA,8,50,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(8,50,existence=True)
        True
    """
    M = OA_from_Vmt(6,7,[0, 1, 3, 16, 35, 26, 36])
    return M

def OA_7_51():
    r"""
    Returns an OA(7,51)

    As explained in the Handbook III.3.69 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_51
        sage: OA = OA_7_51()
        sage: print is_orthogonal_array(OA,7,51,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,51,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(51)

    M = [
        [   5,  33,  29,  30,   1],
        [   8,   3,  47,  10,  13],
        [  14,  27,   6,  12,  28],
        [   9,  16,  44,  49,  11],
        [  34,  32,  36,  26,  20]
        ]

    Mb=[[0],[0],[0],[0],[0],[0]*51]

    for R in zip(*M):
        for i in range(5):
            for RR in [R, [-x for x in R]]:
                for i,x in enumerate(RR):
                    Mb[i].append(x)
            R = cyclic_shift(R,1)

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_7_52():
    r"""
    Returns an OA(7,52)

    As explained in the Handbook III.3.70 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_52
        sage: OA = OA_7_52()
        sage: print is_orthogonal_array(OA,7,52,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,52,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    from sage.rings.finite_rings.constructor import FiniteField
    F4  = FiniteField(4,'z')
    G13 = FiniteField(13)
    G = F4.cartesian_product(G13)
    z = F4('z')
    assert z**2 == z+1

    M = [
        [    (0,0),    (0,0),    (0,0),    (0,0),   (0,0)],
        [(z**2,10),    (0,7),   (1,10),   (z,10),(z**2,3)],
        [   (z,10), (z**2,2),   (1,11),    (z,2),(z**2,7)],
        [    (z,8),(z**2,12),   (0,10),(z**2,11),(z**2,6)],
        [    (1,2),    (0,2), (z**2,8),    (z,3),   (z,7)],
        [    (1,6),   (z,12),    (0,7), (z**2,6),   (z,2)]
        ]

    M2 = [
        [    (1,1),(z**2,11)],
        [    (z,3),    (1,7)],
        [ (z**2,9),    (z,8)],
        [    (1,4), (z**2,3)],
        [   (z,12),    (1,9)],
        [(z**2,10),    (z,1)]
        ]

    M = [[G(x) for x in L] for L in M]
    M2= [[G(x) for x in L] for L in M2]

    Mb=[[(0,0)]*6]

    from itertools import product
    p = lambda x,y : G(tuple([x*yy for yy in G(y)]))

    def t1(i,R):
        if i > 1:
            return t1(1,t1(i-1,R))
        ((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6)) = R
        return [(z*x3, 3*y3), (z*x1, 3*y1), (z*x2, 3*y2), (z*x6, 3*y6), (z*x4, 3*y4), (z*x5, 3*y5)]

    def t2(i,R):
        if i > 1:
            return t2(1,t2(i-1,R))
        ((x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6)) = R
        return [(  x3,   y3), (  x1,   y1), (  x2,   y2), (  x5,   y5), (  x6,   y6), (  x4,   y4)]

    for R in zip(*M):
        for c1,c2 in product([1,2,3],repeat=2):
            Mb.append(t2(c2,t1(c1,R)))

    for R in zip(*M2):
        for c2 in [1,2,3]:
            Mb.append(t2(c2,R))

    Mb = zip(*Mb)
    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_7_54():
    r"""
    Returns an OA(7,54)

    As explained in the Handbook III.3.71 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_54
        sage: OA = OA_7_54()
        sage: print is_orthogonal_array(OA,7,54,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,54,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(45)

    M = [
        [None,None,None,None,None,None,None,None,None],
        [   0,   0,   0,   0,   0,   0,   0,   0,   0],
        [   1,  27,  16,   7,  -1, -27, -16,  -7,   3],
        [  24,  40,   1,  35, -24, -40,  -1, -35,   7],
        [  10,  30,  22,  44, -10, -30, -22, -44,   7],
        [   5,  18,  14,  33,  -5, -18, -14, -33,   3],
        [  30,  16,  33,  27, -30, -16, -33, -27,   0],
        ]

    Mb=[[] for _ in range(7)]

    for R in zip(*M):
        for c in range(7):
            for i,x in enumerate(cyclic_shift(R,c)):
                Mb[i].append(x)

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_8_55():
    r"""
    Returns an OA(8,55)

    As explained in the Handbook III.3.72 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_8_55
        sage: OA = OA_8_55()
        sage: print is_orthogonal_array(OA,8,55,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(8,55,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(55)

    M = [
        [ 1  ,  7 , 14 , 19 , 28 , 33 , 40 , 46 , 50],
        [ 2  , 13 , 25 , 38 , 52 , 12 , 20 , 32 , 45],
        [ 39 , 6  ,  8 , 26 , 24 , 51 , 11 , 34 , 37],
        [ 54 , 48 , 41 , 36 , 27 , 22 , 15 , 9  ,  5],
        [ 53 , 42 , 30 , 17 , 3  , 43 , 35 , 23 , 10],
        [ 16 , 49 , 47 , 29 , 31 , 4  , 44 , 21 , 18]
        ]

    Mb=[[0],[0],[0],[0],[0],[0],[0]*55]

    for R in zip(*M):
        for c in range(6):
            for i,x in enumerate(cyclic_shift(R,c)):
                Mb[i].append(x)

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M


def OA_9_56():
    r"""
    Returns an OA(9,56)

    As explained in the Handbook III.3.73 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_9_56
        sage: OA = OA_9_56()
        sage: print is_orthogonal_array(OA,9,56,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(9,56,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    F8  = FiniteField(8,'z')
    F7  = FiniteField(7)
    G   = F8.cartesian_product(F7)

    w = F8.primitive_element()
    assert w**3 == w+1

    M = [
        [(0,0), (w**0,0), (w**1,0), (w**2,0), (w**3,0), (w**4,0), (w**5,0), (w**6,0)],
        [(0,1), (w**1,6), (w**2,1), (w**3,1), (w**4,6), (w**5,1), (w**6,6), (w**0,6)],
        [(0,4), (w**2,3), (w**3,4), (w**4,4), (w**5,3), (w**6,4), (w**0,3), (w**1,3)],
        [(0,2), (w**3,5), (w**4,2), (w**5,2), (w**6,5), (w**0,2), (w**1,5), (w**2,5)],
        [(0,2), (w**4,5), (w**5,2), (w**6,2), (w**0,5), (w**1,2), (w**2,5), (w**3,5)],
        [(0,4), (w**5,3), (w**6,4), (w**0,4), (w**1,3), (w**2,4), (w**3,3), (w**4,3)],
        [(0,1), (w**6,6), (w**0,1), (w**1,1), (w**2,6), (w**3,1), (w**4,6), (w**5,6)],
        [(1,0), (   1,0), (   1,0), (   1,0), (   1,0), (   1,0), (   1,0), (   1,0)]
        ]

    Mb=[[] for _ in range(8)]

    for R in zip(*M):
        for _ in range(7):
            for i,e in enumerate(R):
                Mb[i].append(e)
            (x1,y1),(x2,y2),(x3,y3),(x4,y4),(x5,y5),(x6,y6),(x7,y7),(x8,y8) = R
            R = [(w*x7,y7), (w*x1,y1), (w*x2,y2), (w*x3,y3), (w*x4,y4), (w*x5,y5), (w*x6,y6), (w*x8,y8)]

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_7_60():
    r"""
    Returns an OA(7,60)

    As explained in [JulianAbel13]_.

    REFERENCES:

    .. [JulianAbel13] Existence of Five MOLS of Orders 18 and 60
      R. Julian R. Abel
      Journal of Combinatorial Designs
      2013

    http://onlinelibrary.wiley.com/doi/10.1002/jcd.21384/abstract

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_60
        sage: OA = OA_7_60()
        sage: print is_orthogonal_array(OA,7,60,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,60,existence=True)
        True
    """
    M60 = [[(0,  0), (0, 0), (0,  0), (0,  0), (0,  0), (0,  0), (0,  0), (0,  0), (0,  0), (0,  0)],
           [(1, 10), (1, 6), (0, 17), (0,  7), (1,  5), (0,  9), (0,  3), (1, 13), (1, 17), (0, 13)],
           [(1, 22), (1, 1), (1,  8), (0,  9), (1, 21), (1, 29), (1,  0), (0,  2), (0, 12), (1, 15)],
           [(1, 24), (1, 1), (0, 14), (0,  0), (0, 16), (0, 18), (0,  8), (0, 28), (0, 17), (0,  7)],
           [(0, 17), (0, 7), (0, 20), (0,  1), (1,  4), (0, 26), (0, 19), (0, 28), (1, 21), (0,  6)],
           [(1, 14), (1, 9), (0, 10), (0, 27), (1, 20), (0, 11), (0, 13), (1, 12), (0, 28), (1, 18)]]

    from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
    G = AdditiveAbelianGroup([2,30])
    M60b=[[],[],[],[],[],[]]
    onezero = G((1,0))

    for R in zip(*M60):
        a,b,c,d,e,f = map(G,R)
        M60b[0].extend([a,c,b,-d,-e,-f])
        M60b[1].extend([b,a,c,-e,-f,-d])
        M60b[2].extend([c,b,a,-f,-d,-e])
        M60b[3].extend([d,e,f,-a+onezero,-c+onezero,-b+onezero])
        M60b[4].extend([e,f,d,-b+onezero,-a+onezero,-c+onezero])
        M60b[5].extend([f,d,e,-c+onezero,-b+onezero,-a+onezero])

    M = OA_from_quasi_difference_matrix(M60b,G)
    return M

def OA_7_62():
    r"""
    Returns an OA(7,62)

    As explained in the Handbook III.3.74 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_7_62
        sage: OA = OA_7_62()
        sage: print is_orthogonal_array(OA,7,62,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(7,62,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as AdditiveCyclic
    G = AdditiveCyclic(54)

    M = [
        [ 0 ,None,None,None, 0   ,None ,None ,None,None,None],
        [17 , 0  , 0  , 0  , -17 ,  0  ,  0  ,  0 ,  1 , 11 ],
        [29 , 28 , 35 , 23 , -29 , -28 , -35 , -23,  3 , 19 ],
        [36 , 50 , 5  , 33 , -36 , -50 , -5  , -33,  7 , 33 ],
        [31 ,  2 , 43 , 30 , -31 , - 2 , -43 , -30, 34 , 33 ],
        [16 , 47 , 44 , 51 , -16 , -47 , -44 , -51, 30 , 19 ],
        [41 , 11 ,  1 , 17 , -41 , -11 , - 1 , -17, 28 , 11 ]
        ]

    Mb=[[] for _ in range(7)]

    for R in zip(*M):
        for c in range(7):
            for i,x in enumerate(cyclic_shift(R,c)):
                Mb[i].append(x)

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = False)
    return M

def OA_9_65():
    r"""
    Returns an OA(9,65)

    Construction shared by Julian R. Abel

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_9_65
        sage: OA = OA_9_65()
        sage: print is_orthogonal_array(OA,9,65,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(9,65,existence=True)
        True
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing as G
    from orthogonal_arrays import orthogonal_array

    B = [None,1, 6, 7, 9, 19, 38, 42, 49] # Base block of a (57,8,1)-BIBD
    OA = orthogonal_array(9,9,2)
    M = [R for R in OA if any(R[0] != x for x in R)]

    M = [[B[x] for x in R] for R in M] # replacing [0,..,8] by the elements of B
    M.append([0]*9)

    M = OA_from_quasi_difference_matrix(zip(*M), G(57),add_col=False)
    return M

def OA_9_75():
    r"""
    Returns an OA(9,75)

    As explained in the Handbook III.3.75 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_9_75
        sage: OA = OA_9_75()
        sage: print is_orthogonal_array(OA,9,75,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(9,75,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    from sage.categories.cartesian_product import cartesian_product

    F3 = FiniteField(3)
    F5 = FiniteField(5)
    G  = cartesian_product((F3,F5,F5))

    M = [
        [(2,0,0), (0,0,0), (0,0,0), (1,0,0), (0,0,0), (1,0,0), (1,0,0), (0,0,0)],
        [(0,2,3), (1,4,4), (1,1,3), (1,0,4), (2,4,3), (0,0,3), (1,4,4), (0,0,0)],
        [(1,3,2), (2,1,1), (1,4,0), (0,3,0), (1,0,4), (2,4,1), (0,1,2), (0,0,0)],
        [(0,2,4), (1,3,1), (2,0,2), (0,0,1), (2,4,0), (1,2,2), (0,0,0), (0,0,0)],
        [(1,1,2), (2,2,3), (0,3,1), (1,4,2), (2,1,0), (1,4,3), (2,4,4), (0,0,0)],
        [(0,1,4), (0,4,4), (2,4,1), (1,3,0), (1,3,1), (2,0,0), (2,4,0), (0,0,0)],
        [(0,4,4), (2,0,1), (2,3,3), (2,3,2), (0,0,2), (2,1,2), (1,4,2), (0,0,0)],
        [(2,4,2), (2,4,1), (2,3,1), (1,2,2), (1,3,0), (0,0,2), (2,4,2), (0,0,0)]
        ]

    for i in range(8):
        M[i].extend(M[7-i][:7])

    Mb=[]

    for R in zip(*M):
        for x in range(5):
            V = [(0,0,x), (0,x,0), (0,x,2*x),(0,2*x,2*x), (0,3*x,3*x), (0,4*x,3*x), (0,4*x,0), (0,0,4*x)]
            Mb.append([G(e)+G(ee) for e,ee in zip(R,V)])

    Mb = zip(*Mb)
    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_11_80():
    r"""
    Returns an OA(11,80)

    As explained in the Handbook III.3.76 [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_quasi_difference_matrix`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_11_80
        sage: OA = OA_11_80()
        sage: print is_orthogonal_array(OA,11,80,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(11,80,existence=True)
        True
    """
    from sage.rings.finite_rings.constructor import FiniteField
    F16 = FiniteField(2**4,'w')
    F5  = FiniteField(5)
    G  = F5.cartesian_product(F16)

    w = F16.primitive_element()
    assert w**4 == w+1

    A = [
        [(0,None), (0,None), (0,None), (0,None), (0,None), (0,None), (0,None), (0,None), (0,None), (0,None)],
        [(0,None), (1,None),    (2,3), (3,None),    (4,3), (2,None),    (3,3), (4,None),    (0,3),    (1,3)],
        [(0,None),    (2,8),    (4,6),    (1,3),    (3,3),   (3,13),   (0,13),    (2,6),   (4,14),   (1,12)],
        [(0,None),   (3,11),    (1,0),    (4,9),    (2,0),    (3,7),    (1,8),   (4,10),   (2,10),   (0,11)],
        [(0,None),    (4,8),   (3,14),   (2,14),   (1,12),   (2,10),   (1,10),    (0,3),    (4,5),    (3,8)],
        [(0,None),    (1,8),   (4,14),   (4,12),    (1,1),    (0,1),    (2,8),   (3,12),    (3,6),    (2,1)],
        [(1,None),    (0,6),    (1,1),    (4,4),   (4,13),    (2,6),   (0,14),    (2,9),    (3,0),    (3,3)],
        [(4,None),    (1,9),    (0,7),    (1,1),    (4,8),    (3,5),   (2,14),    (0,0), (2,None),    (3,0)],
        [(4,None),    (4,6),    (1,2), (0,None),   (1,13),    (3,8),    (3,2),    (2,0),   (0,14), (2,None)],
        [(1,None),    (4,9),    (4,1),    (1,0),    (0,4),    (2,5), (3,None),    (3,5), (2,None), (0,None)]
        ]
    Y = [
        [None, 0, 1, 14, 12, 7, 2, 11, 3, 6],
        [None, 1, 2, 0 , 13, 8, 3, 12, 4, 7],
        [None, 2, 3, 1 , 14, 9, 4, 13, 5, 8]
        ]
    r = lambda x : F16(0) if x is None else w**x

    A = [[(F5(a),r(b)) for a,b in L] for L in A]
    Y = [[r(b) for b in L] for L in Y]

    def t(i,(x,y)):
        a,b = A[x][y]
        b = Y[i][x]
        return G((F5(0),b))

    Mb = [[] for _ in range(10)]

    for y in range(len(A[0])):
        for x in range(len(A)):
            t1,t2,t3 = t(0,(x,y)), t(1,(x,y)), t(2,(x,y))
            e = G(A[x][y])
            Mb[x].extend([e,e+t1,e+t2,e+t1+t2,e+t3,e+t3+t1,e+t3+t2,e+t3+t2+t1])

    M = OA_from_quasi_difference_matrix(Mb,G,add_col = True)
    return M

def OA_10_82():
    r"""
    Returns an OA(10,82)

    Given by Julian R. Abel, using a `V(m,t)` from the Handbook
    [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_10_82
        sage: OA = OA_10_82()
        sage: print is_orthogonal_array(OA,10,82,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(10,82,existence=True)
        True
    """
    M = OA_from_Vmt(8,9,[0,1,20,70,23,59,3,8,19])
    return M

def OA_10_100():
    r"""
    Returns an OA(10,100)

    Given by Julian R. Abel, using a `V(m,t)` from the Handbook
    [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_10_100
        sage: OA = OA_10_100()
        sage: print is_orthogonal_array(OA,10,100,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(10,100,existence=True)
        True
    """
    M = OA_from_Vmt(8,11,[0,1,6,56,22,35,47,23,60])
    return M

def OA_12_144():
    r"""
    Returns an OA(12,144)

    Given by Julian R. Abel, using a `V(m,t)` from the Handbook
    [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_12_144
        sage: OA = OA_12_144()
        sage: print is_orthogonal_array(OA,12,144,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(12,144,existence=True)
        True
    """
    M = OA_from_Vmt(10,13,[0, 1, 5, 10, 22, 6, 14, 9, 53, 129, 84])
    return M

def OA_12_210():
    r"""
    Returns an OA(12,210)

    Given by Julian R. Abel, using a `V(m,t)` from the Handbook
    [DesignHandbook]_.

    .. SEEALSO::

        :func:`sage.combinat.designs.orthogonal_arrays.OA_from_Vmt`

    EXAMPLES::

        sage: from sage.combinat.designs.orthogonal_arrays import is_orthogonal_array
        sage: from sage.combinat.designs.database import OA_12_210
        sage: OA = OA_12_210()
        sage: print is_orthogonal_array(OA,12,210,2)
        True

    The design is available from the general constructor::

        sage: designs.orthogonal_array(12,210,existence=True)
        True
    """
    M = OA_from_Vmt(10,19,[0, 1, 3, 96, 143, 156, 182, 142, 4, 189, 25])
    return M

# Index of the OA constructions
#
# Associates to n the pair (k,f) where f() is a function that returns an OA(k,n)
#
# This dictionary is used by designs.orthogonal_array(k,n).

OA_constructions = {
    20  : (6  , OA_6_20),
    21  : (7  , OA_7_21),
    22  : (5  , OA_5_22),
    24  : (9  , OA_9_24),
    26  : (6  , OA_6_26),
    28  : (7  , OA_7_28),
    30  : (6  , OA_6_30),
    33  : (7  , OA_7_33),
    34  : (6  , OA_6_34),
    35  : (7  , OA_7_35),
    36  : (10 , OA_10_36),
    38  : (6  , OA_6_38),
    39  : (7  , OA_7_39),
    40  : (9  , OA_9_40),
    42  : (7  , OA_7_42),
    44  : (7  , OA_7_44),
    45  : (8  , OA_8_45),
    46  : (6  , OA_6_46),
    48  : (10 , OA_10_48),
    50  : (8  , OA_8_50),
    51  : (7  , OA_7_51),
    52  : (7  , OA_7_52),
    54  : (7  , OA_7_54),
    55  : (8  , OA_8_55),
    56  : (9  , OA_9_56),
    60  : (7  , OA_7_60),
    62  : (7  , OA_7_62),
    65  : (9  , OA_9_65),
    75  : (9  , OA_9_75),
    80  : (11 , OA_11_80),
    82  : (10 , OA_10_82),
    100 : (10 , OA_10_100),
    144 : (12 , OA_12_144),
    210 : (12 , OA_12_210)
}

