r"""
Rubik's cube group functions

NOTE: ``Rubik's cube'' is trademarked.

AUTHOR:
    - David Joyner (2006-10-21): first version


REFERENCES:
    Cameron, P., Permutation Groups. New York: Cambridge University Press, 1999.
    Wielandt, H., Finite Permutation Groups. New York: Academic Press, 1964.
    Dixon, J. and Mortimer, B., Permutation Groups, Springer-Verlag, Berlin/New York, 1996.

TODO:
  Implement 3D plotting (animation?) methods.
  Implement solving methods (see word_problem in permgroup_element.py).


"""

#*****************************************************************************
#       Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.perm_gps.permgroup import PermutationGroup,PermutationGroup_generic, SymmetricGroup, PermutationGroup_subgroup
import random
from types import ListType

import sage.structure.element as element
import sage.groups.group as group

from sage.rings.all      import RationalField, Integer
#from sage.matrix.all     import MatrixSpace
from sage.interfaces.all import gap, is_GapElement, is_ExpectElement
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
import sage.structure.coerce as coerce
from sage.rings.finite_field import GF
from sage.rings.arith import factor
from sage.groups.abelian_gps.abelian_group import AbelianGroup


class CubeGroup(PermutationGroup_generic):
    """
    A python class to help compute Rubik's cube group actions.
    If G denotes the cube group then it may be regarded as a subgroup
    of SymmetricGroup(48), where the 48 facets are labeled as follows.


                             +--------------+
                             |  1    2    3 |
                             |  4  top    5 |
                             |  6    7    8 |
              +--------------+--------------+--------------+--------------+
              |  9   10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
              | 12  left  13 | 20 front  21 | 28 right  29 | 36  rear  37 |
              | 14   15   16 | 22   23   24 | 30   31   32 | 38   39   40 |
              +--------------+--------------+--------------+--------------+
                             | 41   42   43 |
                             | 44 bottom 45 |
                             | 46   47   48 |
                             +--------------+



    EXAMPLES:
        sage: rubik = CubeGroup()
        sage: rubik
        The PermutationGroup of all legal moves of the Rubik's cube.
        sage: print rubik
        The Rubik's cube group with genrators R,L,F,B,U,D in SymmetricGroup(48).

    """
    def __init__(self):
        U = "( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19)" ## U = top
        L = "( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35)" ## L = left
        F = "(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11)" ## F = front
        R = "(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24)" ## R = right
        B = "(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27)" ## B = back or rear
        D = "(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40)" ## D = down or bottom
        self.__gens = [R,L,F,B,U,D]
        self._group = PermutationGroup([R,L,F,B,U,D])
	#H = SymmetricGroup(48)
        #PermutationGroup_subgroup(H,self.__gens)    #### very slow..
	self._group

    def __str__(self):
	return "The Rubik's cube group with genrators R,L,F,B,U,D in SymmetricGroup(48)."

    def __repr__(self):
	return "The PermutationGroup of all legal moves of the Rubik's cube."

    def __call__(self,other):
    	"""
    	EXAMPLES:
    	    sage: rubik = CubeGroup()
            sage: rubik(1)
            ()
    	"""
    	return self._group(other)

    def group(self):
        return self._group

    def gens(self):
        return self.__gens

    def R(self):
        G = self.group()
	g = G(self.gens()[0])
        return g

    def L(self):
        G = self.group()
	g = G(self.gens()[1])
        return g

    def F(self):
        G = self.group()
	g = G(self.gens()[2])
        return g

    def B(self):
        G = self.group()
	g = G(self.gens()[3])
        return g

    def U(self):
        G = self.group()
	g = G(self.gens()[4])
        return g

    def D(self):
	G = self.group()
	g = G(self.gens()[5])
        return g

    def facets(self):
        """
        Returns the set of facets on which the group acts. This function is a "constant".

        EXAMPLES:
            sage: rubik = CubeGroup()
            sage: rubik.facets()==range(1,49)
            True

        """
        fcts = [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 13,14,15,16,17,18,19,20,21,22,23,24,\
          25,26,27,28,29,30,31,32,33,34,35,36, 37,38,39,40,41,42,43,44,45,46,47,48 ]
        return fcts

    def faces(self, mv):
        """
        Returns the set of faces created by the effect of the move mv.

        EXAMPLES:
            sage: rubik = CubeGroup()
            sage: rubik.faces("R*L")
            [[[27, 29, 32], [26, 0, 31], [25, 28, 30]],
             [[11, 13, 16], [10, 0, 15], [9, 12, 14]],
             [[17, 2, 38], [20, 0, 36], [22, 7, 33]],
             [[40, 42, 19], [37, 0, 21], [35, 47, 24]],
             [[41, 18, 3], [44, 0, 5], [46, 23, 8]],
             [[48, 34, 6], [45, 0, 4], [43, 39, 1]]]
            sage: rubik.faces("")
            [[[25, 26, 27], [28, 0, 29], [30, 31, 32]],
             [[9, 10, 11], [12, 0, 13], [14, 15, 16]],
             [[1, 2, 3], [4, 0, 5], [6, 7, 8]],
             [[41, 42, 43], [44, 0, 45], [46, 47, 48]],
             [[17, 18, 19], [20, 0, 21], [22, 23, 24]],
             [[33, 34, 35], [36, 0, 37], [38, 39, 40]]]

        """
        fcts = self.move(mv)[1]
        faceR = [[fcts[24],fcts[25],fcts[26]],[fcts[27],0,fcts[28]],[fcts[29],fcts[30],fcts[31]]]
        faceL = [[fcts[8],fcts[9],fcts[10]],[fcts[11],0,fcts[12]],[fcts[13],fcts[14],fcts[15]]]
        faceU = [[fcts[0],fcts[1],fcts[2]],[fcts[3],0,fcts[4]],[fcts[5],fcts[6],fcts[7]]]
        faceD = [[fcts[40],fcts[41],fcts[42]],[fcts[43],0,fcts[44]],[fcts[45],fcts[46],fcts[47]]]
        faceF = [[fcts[16],fcts[17],fcts[18]],[fcts[19],0,fcts[20]],[fcts[21],fcts[22],fcts[23]]]
        faceB = [[fcts[32],fcts[33],fcts[34]],[fcts[35],0,fcts[36]],[fcts[37],fcts[38],fcts[39]]]
        return [faceR,faceL,faceU,faceD,faceF,faceB]

    def move(self,mv):
        """
       Returns the group element and the reordered list of facets, as moved by
	the list mv (read left-to-right)

	INPUT: mv is a string of the form R^a*Y^b*...",
	       where X, Y, ... are in {R,L,F,B,U,D}
	       and a,b, ... are integers.

	EXAMPLES:
            sage: rubik = CubeGroup()
	    sage: rubik.move("")[0]
	    ()
	    sage: rubik.move("R")[0]
	    (3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)
	    sage: rubik.R()
	    (25,27,32,30)(26,29,31,28)(3,38,43,19)(5,36,45,21)(8,33,48,24)

	"""
	m = mv.split("*")
	M = [x.split("^") for x in m]
	#print M
	n = len(M)
        e = 0
	G = self.group()
	R,L,F,B,U,D = G.gens()
	g = G(1)
	fcts = self.facets()
	for i in range(n):
	    x = M[i][0]
	    if len(M[i])==1:
                if x == "R":   h = self.R()
	        elif x == "L": h = self.L()
	        elif x == "U": h = self.U()
	        elif x == "D": h = self.D()
	        elif x == "F": h = self.F()
	        elif x == "B": h = self.B()
	        else: h = G(1)
	        g = g*h
	    if len(M[i])==2:
                if x == "R":   h = self.R()
	        elif x == "L": h = self.L()
	        elif x == "U": h = self.U()
	        elif x == "D": h = self.D()
	        elif x == "F": h = self.F()
	        elif x == "B": h = self.B()
	        else: h = G(1)
	        e = M[i][1]
		if e==1: g = g*h
		if e==2: g = g*h*h
		if e==3: g = g*h*h*h
	pos = [g(i) for i in fcts]
	return [g,pos]

    def display2d(self,mv):
        """
        Displays a 2d map of the Rubik's cube after the move mv has been made.
        Nothing is returned.

        EXAMPLES:
            sage: rubik = CubeGroup()
            sage: rubik.display2d("")
                                 +--------------+
                                 |  1    2    3 |
                                 |  4  top    5 |
                                 |  6    7    8 |
                  +--------------+--------------+--------------+--------------+
                  |  9   10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
                  | 12  left  13 | 20 front  21 | 28 right  29 | 36  rear  37 |
                  |  14   15   16 | 22   23   24 | 30   31   32 | 38   39   40 |
                  +--------------+--------------+--------------+--------------+
                                 | 41   42   43 |
                                 | 44 bottom 45 |
                                 | 46   47   48 |
                                 +--------------+

            sage: rubik.display2d("R")
                                 +--------------+
                                 |  1    2    38 |
                                 |  4  top    36 |
                                 |  6    7    33 |
                  +--------------+--------------+--------------+--------------+
                  |  9   10   11 | 17   18   3 | 27   29   32 | 48   34   35 |
                  | 12  left  13 | 20 front  5 | 26 right  31 | 45  rear  37 |
                  |  14   15   16 | 22   23   8 | 25   28   30 | 43   39   40 |
                  +--------------+--------------+--------------+--------------+
                                 | 41   42   19 |
                                 | 44 bottom 21 |
                                 | 46   47   24 |
                                 +--------------+

        You can see the right face has been rotated but not the left face.
        """
        lst = self.move(mv)[1]
        line1 = "                     +--------------+\n"
        line2 = "                     |  %s    %s    %s |\n"%(lst[0],lst[1],lst[2])
        line3 = "                     |  %s  top    %s |\n"%(lst[3],lst[4])
        line4 = "                     |  %s    %s    %s |\n"%(lst[5],lst[6],lst[7])
        line5 = "      +--------------+--------------+--------------+--------------+\n"
        line6 = "      |  %s   %s   %s | %s   %s   %s | %s   %s   %s | %s   %s   %s |\n"%(lst[8],lst[9],lst[10],lst[16],lst[17],lst[18],lst[24],lst[25],lst[26],lst[32],lst[33],lst[34])
        line7 = "      | %s  left  %s | %s front  %s | %s right  %s | %s  rear  %s |\n"%(lst[11],lst[12],lst[19],lst[20],lst[27],lst[28],lst[35],lst[36])
        line8 = "      |  %s   %s   %s | %s   %s   %s | %s   %s   %s | %s   %s   %s |\n"%(lst[13],lst[14],lst[15],lst[21],lst[22],lst[23],lst[29],lst[30],lst[31],lst[37],lst[38],lst[39])
        line9 = "      +--------------+--------------+--------------+--------------+\n"
        line10 = "                     | %s   %s   %s |\n"%(lst[40],lst[41],lst[42])
        line11 = "                     | %s bottom %s |\n"%(lst[43],lst[44])
        line12 = "                     | %s   %s   %s |\n"%(lst[45],lst[46],lst[47])
        line13 = "                     +--------------+\n"
        print line1+line2+line3+line4+line5+line6+line7+line8+line9+line10+line11+line12+line13
