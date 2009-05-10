r"""
Hexads in S(5,6,12)

  This module completes a 5-element subset of as 12-set X
  into a hexad in a Steiner system S(5,6,12) using Curtis and
  Conway's "kitten method".  The labeling is either the
  "modulo 11" labeling or the "shuffle" labeling.

  The  main functions implemented in this file are
  blackjack_move and find_hexad. Enter "blackjack_move?"
  for help to play blackjack (i.e., the rules of the game), or
  "find_hexad?" for help finding hexads of S(5,6,12) in the
  shuffle labeling.

  This picture is the kitten in the "shuffle" labeling:


                                     6

                                     9

                                10     8

                             7     2      5

                         9     4     11     9

                      10     8     3      10     8

                 1                                      0


The corresponding MINIMOG is

                  _________________
                  |  6  |  3  |  0  |  9  |
                  |---|---|---|---|
                  |  5  |  2  |  7  | 10  |
                  |---|---|---|---|
                  |  4  |  1  |  8  | 11  |
                  |___|____|___|____|

which is specified by the global variable "minimog_shuffle".
See the docstrings for find_hexad and blackjack_move for
further details and examples.

AUTHOR::
  David Joyner (2006-05)

REFERENCES::
    R. Curtis, ``The Steiner system $S(5,6,12)$, the Mathieu group $M_{12}$,
    and the kitten,'' in {\bf Computational group theory}, ed. M. Atkinson,
    Academic Press, 1984.
    J. Conway, ``Hexacode and tetracode - MINIMOG and MOG,'' in {\bf Computational
    group theory}, ed. M. Atkinson, Academic Press, 1984.
    J. Conway and N. Sloane, ``Lexicographic codes: error-correcting codes from
    game theory,'' IEEE Trans. Infor. Theory32(1986)337-348.
    J. Kahane and A. Ryba, ``The hexad game,'' Electronic Journal of Combinatorics, 8 (2001)
    http://www.combinatorics.org/Volume_8/Abstracts/v8i2r11.html

    Some details are also online at:  http://www.permutationpuzzles.org/hexad/

"""
#*****************************************************************************
#       Copyright (C) 2005 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your choice)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.infinity import Infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.rational_field import RationalField
QQ = RationalField()
infinity = Infinity
from sage.rings.finite_field import FiniteField
GF = FiniteField



def setting_up(MINIMOG):
    """
    This initializes the variables in the game.

    EXAMPLES:
        sage: import sage.games.hexad
        sage: XX = sage.games.hexad.setting_up(minimog_shuffle)

    """
    MS34 = MatrixSpace(QQ,3,4)
    A = MS34(MINIMOG)
    MS33 = MatrixSpace(QQ,3,3)
    picture00 = MS33([[A[(1,0)],A[(2,3)],A[(0,1)]],[A[(2,2)],A[(1,1)],A[(2,0)]],[A[(0,3)],A[(1,3)],A[(1,2)]]])
    ####### picture00 is the "picture at 6"
    picture02 = MS33([[A[(1,0)],A[(2,3)],A[(0,1)]],[A[(1,1)],A[(2,0)],A[(2,2)]],[A[(1,2)],A[(0,3)],A[(1,3)]]])
    ####### picture02 is the "picture at 1"
    picture21 = MS33([[A[(2,2)],A[(1,3)],A[(0,1)]],[A[(0,3)],A[(2,3)],A[(2,0)]],[A[(1,0)],A[(1,1)],A[(1,2)]]])
    ####### picture21 is the "picture at 0"

    line = [set([]) for i in range(12)]
    line[0] = set([(0,0),(0,1),(0,2)])
    line[1] = set([(1,0),(1,1),(1,2)])
    line[2] = set([(2,0),(2,1),(2,2)])
    line[3] = set([(0,2),(1,2),(2,2)])
    line[4] = set([(0,1),(1,1),(2,1)])
    line[5] = set([(0,0),(1,0),(2,0)])
    line[6] = set([(0,0),(1,1),(2,2)])
    line[7] = set([(2,0),(0,1),(1,2)])
    line[8] = set([(0,2),(1,0),(2,1)])
    line[9] = set([(2,0),(1,1),(0,2)])
    line[10] = set([(0,0),(1,2),(2,1)])
    line[11] = set([(1,0),(0,1),(2,2)])

    cross = [set([]) for i in range(18)]
    cross[0] = set([(0,0),(0,1),(0,2),(1,0),(2,0)])
    cross[1] = set([(0,0),(0,1),(0,2),(1,2),(2,2)])
    cross[2] = set([(0,0),(1,0),(2,0),(2,1),(2,2)])
    cross[3] = set([(2,0),(2,1),(2,2),(0,2),(1,2)])
    cross[4] = set([(0,0),(0,1),(0,2),(1,1),(2,1)])
    cross[5] = set([(0,0),(1,0),(2,0),(1,1),(1,2)])
    cross[6] = set([(1,0),(1,1),(1,2),(0,2),(2,2)])
    cross[7] = set([(0,1),(1,1),(2,1),(2,0),(2,2)])
    cross[8] = set([(0,0),(0,1),(1,0),(1,1),(2,2)])
    cross[9] = set([(0,0),(1,1),(1,2),(2,1),(2,2)])
    cross[10] = set([(2,0),(2,1),(1,0),(1,1),(0,2)])
    cross[11] = set([(0,1),(0,2),(1,1),(1,2),(2,0)])
    cross[12] = set([(0,0),(1,0),(0,2),(1,2),(2,1)])
    cross[13] = set([(1,0),(0,1),(0,2),(2,1),(2,2)])
    cross[14] = set([(0,1),(1,0),(1,2),(2,0),(2,2)])
    cross[15] = set([(0,0),(0,1),(1,2),(2,0),(2,1)])
    cross[16] = set([(1,0),(1,1),(1,2),(0,1),(2,1)])
    cross[17] = set([(0,0),(0,2),(1,1),(2,0),(2,2)])
    box = set([(i,j) for i in range(3) for j in range(3)])
    square = [set([]) for i in range(18)]
    for i in range(18):
       square[i] = box - cross[i]

    MS34 = MatrixSpace(GF(3),3,4)
    col1 = MS34([[1,0,0,0],[1,0,0,0],[1,0,0,0]])
    col2 = MS34([[0,1,0,0],[0,1,0,0],[0,1,0,0]])
    col3 = MS34([[0,0,1,0],[0,0,1,0],[0,0,1,0]])
    col4 = MS34([[0,0,0,1],[0,0,0,1],[0,0,0,1]])

    tet1 = MS34([[1,1,1,1],[0,0,0,0],[0,0,0,0]])
    tet2 = MS34([[1,0,0,0],[0,1,1,1],[0,0,0,0]])
    tet3 = MS34([[1,0,0,0],[0,0,0,0],[0,1,1,1]])
    tet4 = MS34([[0,1,0,0],[1,0,1,0],[0,0,0,1]])
    tet5 = MS34([[0,0,0,1],[1,1,0,0],[0,0,1,0]])
    tet6 = MS34([[0,0,1,0],[1,0,0,1],[0,1,0,0]])
    tet7 = MS34([[0,1,0,0],[0,0,0,1],[1,0,1,0]])
    tet8 = MS34([[0,0,1,0],[0,1,0,0],[1,0,0,1]])
    tet9 = MS34([[0,0,0,1],[0,0,1,0],[1,1,0,0]])
    col = [col1,col2,col3,col4]
    tet = [tet1,tet2,tet3,tet4,tet5,tet6,tet7,tet8,tet9]
    return picture00,picture02,picture21,line,cross,square,col,tet

############################ initializing "global variables" (for this module) #############

MSq34 = MatrixSpace(QQ,3,4)
minimog_modulo11 = [[0,3,infinity,2],[5,9,8,10],[4,1,6,7]]  ## this is harder to implement,
                                                            ## in SAGE so won't be used
minimog_shuffle = MSq34([[6,3,0,9],[5,2,7,10],[4,1,8,11]])     ## MINIMOG is one of these
XX = setting_up(minimog_shuffle)
picture00 = XX[0]
picture02 = XX[1]
picture21 = XX[2]
line = XX[3]
cross = XX[4]
square = XX[5]
col = XX[6]
tet = XX[7]

#############################################################################

def print_kitten(MINIMOG):
    """
    This simply prints the "kitten" (expressed as a triangular diagram
    of symbols).

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.print_kitten(sage.games.hexad.minimog_shuffle)
                 0
        <BLANKLINE>
                 8
               9  10
              5  11  3
             8  2  4  8
           9  10  7  9  10
        <BLANKLINE>
        6                   1
        sage: sage.games.hexad.print_kitten(sage.games.hexad.minimog_modulo11)
                 +Infinity
        <BLANKLINE>
                 6
               2  10
              5  7  3
             6  9  4  6
           2  10  8  2  10
        <BLANKLINE>
        0                   1

    """
    Kitten1 = [' ',' ',' ',' ',' ',' ',' ',' ',' ',str(MINIMOG[0][2]),' ',' ',' ',' ',' ']
    Kitten2 = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
    Kitten3 = [' ',' ',' ',' ',' ',' ',' ',' ',' ',str(MINIMOG[2][2]),' ',' ',' ',' ',' ']
    Kitten4 = [' ',' ',' ',' ',' ',' ',' ',str(MINIMOG[0][3]),' ',' ',str(MINIMOG[1][3]),' ',' ',' ',' ']
    Kitten5 = [' ',' ',' ',' ',' ',' ',str(MINIMOG[1][0]),' ',' ',str(MINIMOG[2][3]),' ',' ',str(MINIMOG[0][1]),' ',' ',' ']
    Kitten6 = [' ',' ',' ',' ',' ',str(MINIMOG[2][2]),' ',' ',str(MINIMOG[1][1]),' ',' ',str(MINIMOG[2][0]),' ',' ',str(MINIMOG[2][2]),' ',' ']
    Kitten7 = [' ',' ',' ',str(MINIMOG[0][3]),' ',' ',str(MINIMOG[1][3]),' ',' ',str(MINIMOG[1][2]),' ',' ',str(MINIMOG[0][3]),' ',' ',str(MINIMOG[1][3]),' ']
    Kitten8 = [' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
    Kitten9 = [str(MINIMOG[0][0]),' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',str(MINIMOG[2][1])]
    Kitten = [Kitten1, Kitten2, Kitten3, Kitten4, Kitten5, Kitten6, Kitten7, Kitten8, Kitten9]
    kitten = "".join(Kitten1)+"\n"+"".join(Kitten2)+"\n"+"".join(Kitten3)+"\n"+"".join(Kitten4)+"\n"+"".join(Kitten5)+"\n"+"".join(Kitten6)+"\n"+"".join(Kitten7)+"\n"+"".join(Kitten8)+"\n"+"".join(Kitten9)
    print kitten



def view_list(L):
    """
    This provides a printout of the lines, crosses and squares of the
    MINIMOG, as in Curtis' paper.

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.view_list(sage.games.hexad.line[1])
        <BLANKLINE>
        [0 0 0]
        [1 1 1]
        [0 0 0]
        sage: sage.games.hexad.view_list(sage.games.hexad.cross[1])
        <BLANKLINE>
        [1 1 1]
        [0 0 1]
        [0 0 1]
        sage: sage.games.hexad.view_list(sage.games.hexad.square[1])
        <BLANKLINE>
        [0 0 0]
        [1 1 0]
        [1 1 0]

    """
    MS = MatrixSpace(QQ,3,3)
    box = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)]
    A = MS([[0 for i in range(3)] for i in range(3)])
    for x in box:
        if x in L:
            A[x] = 1
    return A

def picture_set(A,L):
    """
    This is needed in the find_hexad function below.

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.picture_set(minimog_shuffle,sage.games.hexad.cross[2])
        set([8, 1, 4, 5, 6])
        sage: sage.games.hexad.picture_set(minimog_shuffle, sage.games.hexad.square[7])
        set([0, 5, 6, 7])

    """
    return set([A[x] for x in L])

def find_hexad0(pts,MINIMOG):
    """
    INPUT::
       pts -- a set of 2 distinct elements of MINIMOG, but not
              including the "points at infinity"
    OUTPUT::
       hexad containing pts and of type 0 (the 3 points "at infinity" union a line)

    NOTE:: 3 points "at infinity" = {MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.find_hexad0(set([2,4]),minimog_shuffle)
        ([0, 1, 2, 4, 6, 8], ['line 1', 'picture 1'])

    """
    L = set(pts)
    H = set([MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]])
    for i in range(12):
        if L <= picture_set(picture02,line[i]):
            WHAT = ["line "+str(i),"picture "+str(1)]
            H = H | picture_set(picture02,line[i])
            return list(H),WHAT
        if L <= picture_set(picture21,line[i]):
            WHAT = ["line "+str(i),"picture "+str(0)]
            H = H | picture_set(picture21,line[i])
            return list(H),WHAT
        if L <= picture_set(picture00,line[i]):
            WHAT = ["line "+str(i),"picture "+str(6)]
            H = H | picture_set(picture00,line[i])
            return list(H),WHAT
    return [],[]


def find_hexad1(pts,MINIMOG):
    """
    INPUT::
       pts -- a set pts of 5 distinct elements of MINIMOG
    OUTPUT::
       hexad containing pts and of type 1 (union of 2 parallel lines -- *no*
       points "at infinity")

    NOTE:: 3 points "at infinity" = {MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.find_hexad1(set([2,3,4,5,8]),minimog_shuffle)
        ([2, 3, 4, 5, 8, 11], ['lines (1, 2)', 'picture 1'])

    """
    H = set(pts)
    L = set(pts)
    linez = [(1,2),(1,3),(2,3),(4,5),(4,6),(5,6),(7,8),(7,9),(8,9),(10,11),(10,12),(11,12)]
    for x in linez:
        x1 = int(x[0]-1)
        x2 = int(x[1]-1)              ## (recall | is union)
        if L <= (picture_set(picture02,line[x1]) | picture_set(picture02,line[x2])):
            WHAT = ["lines "+str(x),"picture "+str(1)]
            H = picture_set(picture02,line[x1]) | picture_set(picture02,line[x2])
            return list(H),WHAT
        if L <= (picture_set(picture21,line[x1]) | picture_set(picture21,line[x2])):
            WHAT = ["lines "+str(x),"picture "+str(0)]
            H = picture_set(picture21,line[x1]) | picture_set(picture21,line[x2])
            return list(H),WHAT
        if L <= (picture_set(picture00,line[x1]) | picture_set(picture00,line[x2])):
            WHAT = ["lines "+str(x),"picture "+str(6)]
            H = picture_set(picture00,line[x1]) | picture_set(picture00,line[x2])
            return list(H),WHAT
    return [],[]


def find_hexad2(pts,x0,MINIMOG):
    """
    INPUT::
       pts -- a list S of 4 elements of MINIMOG, not including
              any "points at infinity"
       x0  -- in {MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}

    OUTPUT::
       hexad containing S \union \{x0\} of type 2

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.find_hexad2([2,3,4,5],1,minimog_shuffle)
        ([], [])

    The above output indicates that there is no hexad of type 2
    containing {2,3,4,5}. However, there is one containing {2,3,4,8}:

        sage: sage.games.hexad.find_hexad2([2,3,4,8],0,minimog_shuffle)
        ([0, 2, 3, 4, 8, 9], ['cross 12', 'picture 0'])

    """
    L = set(pts)
    H = set([x0])
    for i in range(18):
        if (x0 == MINIMOG[2][1] and L <= picture_set(picture02,cross[i])):
            WHAT = ["cross "+str(i),"picture "+str(1)]
            H = H | picture_set(picture02,cross[i])
            return list(H),WHAT
        if (x0 == MINIMOG[0][2] and L <= picture_set(picture21,cross[i])):
            WHAT = ["cross "+str(i),"picture "+str(MINIMOG[0][2])]
            H = H | picture_set(picture21,cross[i])
            return list(H),WHAT
        if (x0 == MINIMOG[0][0] and L <= picture_set(picture00,cross[i])):
            WHAT = ["cross "+str(i),"picture "+str(6)]
            H = H | picture_set(picture00,cross[i])
            return list(H),WHAT
    return [],[]

def find_hexad3(pts,x0,x1,MINIMOG):
    """
    INPUT::
       pts -- a list of 3 elements of MINIMOG, not including
              any "points at infinity"
       x0,x1  -- in {MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}

    OUTPUT::
       hexad containing pts \union \{x0,x1\} of type 3 (square at
       picture of "omitted point at infinity")

    EXAMPLES::
        sage: import sage.games.hexad
        sage: sage.games.hexad.find_hexad3([2,3,4],0,1,minimog_shuffle)
        ([0, 1, 2, 3, 4, 11], ['square 2', 'picture 6'])

    """
    L = set(pts)
    H = set([x0,x1])
    for i in range(18):
        if (not (MINIMOG[0][2] in H) and L <= picture_set(picture21,square[i])):
            WHAT = ["square "+str(i),"picture "+str(MINIMOG[0][2])]
            H = H | picture_set(picture21,square[i])
            return list(H),WHAT
        if (not (MINIMOG[2][1] in H) and L <= picture_set(picture02,square[i])):
            WHAT = ["square "+str(i),"picture "+str(MINIMOG[2][1])]
            H = H | picture_set(picture02,square[i])
            return list(H),WHAT
        if (not (MINIMOG[0][0] in H) and L <= picture_set(picture00,square[i])):
            WHAT = ["square "+str(i),"picture "+str(MINIMOG[0][0])]
            H = H | picture_set(picture00,square[i])
            return list(H),WHAT
    return [],[]

def find_hexad(pts,MINIMOG):
    """
    INPUT::
       pts -- a list S of 5 elements of MINIMOG

    OUTPUT::
       hexad containing S \union \{x0\} of some type

    NOTE:: 3 ``points at infinity" = {MINIMOG[0][2], MINIMOG[2][1], MINIMOG[0][0]}

    Theorem (Curtis, Conway): Each hexads is of exactly one of the following types:
       0. $\{3 ``points at infinity''\}\cup \{{\rm any\ line}\}$,
       1. the union of any two (distinct) parallel lines in the same picture,
       2. one ``point at infinity'' union a cross in the corresponding picture,
       3. two ``points at infinity'' union a square in the picture corresponding
          to the omitted point at infinity.
    More precisely, there are 132 such hexads (12 of type 0, 12 of type 1, 54 of type 2,
    and 54 of type 3). They form a Steiner system of type $(5,6,12)$.


    EXAMPLES::
        sage: find_hexad([0,1,2,3,4],minimog_shuffle)
        ([0, 1, 2, 3, 4, 11], ['square 2', 'picture 6'])
        sage: find_hexad([1,2,3,4,5],minimog_shuffle)
        ([1, 2, 3, 4, 5, 6], ['square 8', 'picture 0'])
        sage: find_hexad([2,3,4,5,8],minimog_shuffle)
        ([2, 3, 4, 5, 8, 11], ['lines (1, 2)', 'picture 1'])
        sage: find_hexad([0,1,2,4,6],minimog_shuffle)
        ([0, 1, 2, 4, 6, 8], ['line 1', 'picture 1'])

    AUTHOR::
       David Joyner (2006-05)

    REFERENCES::
        R. Curtis, ``The Steiner system $S(5,6,12)$, the Mathieu group $M_{12}$,
        and the kitten,'' in {\bf Computational group theory}, ed. M. Atkinson,
        Academic Press, 1984.
        J. Conway, ``Hexacode and tetracode - MINIMOG and MOG,'' in {\bf Computational
        group theory}, ed. M. Atkinson, Academic Press, 1984.

    """
    L = set(pts)
    LL = L.copy()
    ## recall & means intersection
    L2 = LL & set([MINIMOG[0][2],MINIMOG[2][1],MINIMOG[0][0]])
    if len(L2) == 3: ## must be type 0 (line + pts at infty)
        H,WHAT = find_hexad0(LL - set([MINIMOG[0][2],MINIMOG[2][1],MINIMOG[0][0]]),MINIMOG)
        return H,WHAT
    if len(L2) == 2: ## type 0 or 3
        if (MINIMOG[0][2] in LL and MINIMOG[2][1] in LL):
            H,WHAT = find_hexad3(LL - set([MINIMOG[0][2],MINIMOG[2][1]]),MINIMOG[0][2],MINIMOG[2][1],MINIMOG)
            if H != []:   # must be type 3
                return list(H),WHAT
            if H == []:   ## could be type 0
                H,WHAT = find_hexad0(LL - L2,MINIMOG)
                if H != []:   # must be type 0
                    return list(H),WHAT
        if (MINIMOG[2][1] in LL and MINIMOG[0][0] in LL):
            H,WHAT = find_hexad3(LL - set([MINIMOG[2][1],MINIMOG[0][0]]),MINIMOG[2][1],MINIMOG[0][0],MINIMOG)
            if H != []:   # must be type 3
                return list(H),WHAT
            if H == []:  ## could be type 0
                H,WHAT = find_hexad0(LL - L2,MINIMOG)
                if H != []:   # must be type 0
                    return list(H),WHAT
        if (MINIMOG[0][2] in LL and MINIMOG[0][0] in LL):
            H,WHAT = find_hexad3(LL - set([MINIMOG[0][2],MINIMOG[0][0]]),MINIMOG[0][2],MINIMOG[0][0],MINIMOG)
            if H != []:   # must be type 3
                return list(H),WHAT
            if H == []:  ## could be type 0
                H,WHAT = find_hexad0(LL - L2,MINIMOG)
                if H != []:   # must be type 0
                    return list(H),WHAT
    if len(L2) == 1:
        H,WHAT = find_hexad2(LL - L2,list(L2)[0],MINIMOG)
        if H == []:  ## not a cross in picture at infinity
            if list(L2)[0] == MINIMOG[2][1]:
                      L1 = LL - L2
                      H,WHAT = find_hexad3(L1,MINIMOG[0][0],MINIMOG[2][1],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
                      L1 = LL - L2
                      H,WHAT = find_hexad3(L1,MINIMOG[0][2],MINIMOG[2][1],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
            if list(L2)[0] == MINIMOG[0][0]:
                      L1 = (LL - L2)
                      H,WHAT = find_hexad3(L1,MINIMOG[0][0],MINIMOG[2][1],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
                      L1 = (LL - L2)
                      H,WHAT = find_hexad3(L1,MINIMOG[0][0],MINIMOG[0][2],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
            if list(L2)[0] == MINIMOG[0][2]:
                      L1 = (LL - L2)
                      H,WHAT = find_hexad3(L1,MINIMOG[0][0],MINIMOG[0][2],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
                      L1 = (LL - L2)
                      H,WHAT = find_hexad3(L1,MINIMOG[2][1],MINIMOG[0][2],MINIMOG)
                      if H <> []:
                          return list(H),WHAT
        return list(H),WHAT
        ## a cross in a pic at infty
    if len(L2) == 0: ## L is either a union of 2 lines or a cross
        for i in LL:
            for j in [MINIMOG[0][2],MINIMOG[2][1],MINIMOG[0][0]]:
                H,WHAT = find_hexad2(LL - set([i]),j,MINIMOG)
                if (H != [] and i in H):
                    return list(H),WHAT  ## L is in a cross
        H,WHAT = find_hexad1(LL,MINIMOG)  ## L is a union of lines
        return H,WHAT

def blackjack_move(L0,MINIMOG):
    """
    L is a list of cards of length 6, taken from {0,1,...,11}.

         Mathematical blackjack::
    Mathematical blackjack is played with 12 cards, labeled $0,...,11$
    (for example:  king,  ace, $2$, $3$, ..., $10$,  jack, where the  king is
    $0$ and the  jack is $11$). Divide the 12 cards into two piles of $6$ (to be
    fair, this should be done randomly). Each of the $6$ cards of one of these
    piles are to be placed face up on the table. The remaining cards are in a
    stack which is shared and visible to both players. If the sum of the cards
    face up on the table is less than 21 then no legal move is possible so you must
    shuffle the cards and deal a new game. (Conway calls such a game *={0|0},
    where 0={|}; in this game the first player automatically wins.)

    * Players alternate moves.
    * A move consists of exchanging a card on the table with a lower card from the other pile.
    * The player whose move makes the sum of the cards on the table under 21 loses.

    The winning strategy (given below) for this game is due to Conway and Ryba.
    There is a Steiner system $S(5,6,12)$ of hexads in the set $\{0,1,...,11\}$.
    This Steiner system is associated to the MINIMOG of in the "shuffle numbering"
    rather than the "modulo $11$ labeling''.

    Proposition (Ryba, Conway)   For this Steiner system, the winning strategy is to choose a
    move which is a hexad from this system.


    EXAMPLES::
        sage: blackjack_move([0,2,4,6,7,11],minimog_shuffle)
        '4 --> 3. The total went from 30 to 29.'

    Is this really a hexad?

        sage: find_hexad([11,2,3,6,7],minimog_shuffle)
        ([0, 2, 3, 6, 7, 11], ['square 9', 'picture 1'])

    So, yes it is, but here is further confirmation:

        sage: blackjack_move([0,2,3,6,7,11],minimog_shuffle)
        This is a hexad.
        There is no winning move, so make a random legal move.
        [0, 2, 3, 6, 7, 11]

    Now, suppose player 2 replaced the 11 by a 9. Your next move:

        sage: blackjack_move([0,2,3,6,7,9],minimog_shuffle)
        '7 --> 1. The total went from 27 to 21.'

    You have now won. SAGE will even tell you so:

        sage: blackjack_move([0,2,3,6,1,9],minimog_shuffle)
        'No move possible. Shuffle the deck and redeal.'

    AUTHOR::
        David Joyner (2006-05)

    REFERENCES::
        J. Conway and N. Sloane, "Lexicographic codes: error-correcting codes from
        game theory,'' IEEE Trans. Infor. Theory32(1986)337-348.
        J. Kahane and A. Ryba, "The hexad game,'' Electronic Journal of Combinatorics, 8 (2001)
        http://www.combinatorics.org/Volume_8/Abstracts/v8i2r11.html

    """
    total = sum(L0)
    if total <22:
        return "No move possible. Shuffle the deck and redeal."
    L = set(L0)
    for x in L:
        h,WHAT = find_hexad(L - set([x]),MINIMOG)
        if list(L0) == list(h):
            print "      This is a hexad. \n      There is no winning move, so make a random legal move."
            return L0
        y = list(set(h) - (L - set([x])))[0]
        #print x,y,h
        if y < x:
            return str(x) +' --> '+str(y)+". The total went from "+ str(total) +" to "+str(total-x+y)+"."
    print "This is a hexad. \n There is no winning move, so make a random legal move."
    return L0



