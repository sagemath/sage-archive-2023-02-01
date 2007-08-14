r"""
Rubik's cube group functions
\label{sec:rubik}

NOTE: ``Rubik's cube'' is trademarked.

NOTATION: B denotes a clockwise quarter turn of the back face
          D denotes a clockwise quarter turn of the down face
          and similarly for F (front), L (left), R (right), U (up)
          Products of moves are read {\it right to left}, so for example,
          R*U means move U first and then R.

The "Singmaster notation":
  * moves: U, D, R, L, F, B as in the diagram below,
  * corners: xyz means the facet is on face x (in R,F,L,U,D,B)
     and the clockwise rotation of the corner sends x->y->z
  * edges: xy means the facet is on face x and a flip of the edge sends x->y.



                        +-----------------+
                        | 1     2      3 |
                        | 4    up    5 |
                        | 6     7      8 |
     +-----------------+-----------------+--------------------+-------------------+
     | 9    10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
     | 12  left  13 | 20  front 21 | 28  right 29 | 36  back  37 |
     | 14   15   16 | 22    23  24 | 30   31   32 | 38   39   40 |
     +------------------+-----------------+-------------------+-------------------+
                        | 41   42   43 |
                        | 44  down  45 |
                        | 46   47   48 |
                       +------------------+

AUTHOR:
    - David Joyner (2006-10-21): first version
    -      "       (2007-05): changed faces, added legal and solve
    -      "       (2007-06): added plotting functions
    -      "       (2007-08): colors corrected, "solve" rewritten (again),typos fixed.

REFERENCES:
    Cameron, P., Permutation Groups. New York: Cambridge University Press, 1999.
    Wielandt, H., Finite Permutation Groups. New York: Academic Press, 1964.
    Dixon, J. and Mortimer, B., Permutation Groups, Springer-Verlag, Berlin/New York, 1996.
    Joyner, D, Adventures in Group Theory, Johns Hopkins Univ Press, 2002.

"""

#**************************************************************************************
#       Copyright (C) 2006 David Joyner <wdjoyner@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#**************************************************************************************

from sage.groups.perm_gps.permgroup import PermutationGroup,PermutationGroup_generic, PermutationGroup_subgroup
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
import random

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
from sage.plot.plot import PolygonFactory, TextFactory
polygon = PolygonFactory()
text = TextFactory()
from sage.calculus.calculus import Function_sin, Function_cos
sin = Function_sin()
cos = Function_cos()
pi = 3.14159265

####################### predefined colors ##################

red = (1,0,0)                ## F face
green = (0,1,0)            ## R face
blue = (0,0,1)              ## D face
yellow = (1,1,0)           ## L face
white = (1,1,1)             ## none
orange = (1,0.6,0.3)       ## B face
purple = (1,0,1)           ## none
lpurple = (1,0.63,1)       ## U face
lightblue = (0,1,1)        ## none
lgrey = (0.75,0.75,0.75)    ## sagemath.org color

#########################################################
#written by Tom Boothby, placed in the public domain

def xproj(x,y,z,r):
  return (y*r[1] - x*r[3])*r[2]

def yproj(x,y,z,r):
  return z*r[2] - (x*r[1] + y*r[2])*r[0]

def rotation_list(tilt,turn):
    return [ sin(tilt*pi/180.0), sin(turn*pi/180.0), cos(tilt*pi/180.0), cos(turn*pi/180.0) ]

def polygon_plot3d(points, tilt=30, turn=30, **kwargs):
    """
    Plots a polygon viewed from an angle determined by tilt, turn, and vertices
    points.

    WARNING: The ordering of the points is important to get "correct" and
    if you add several of these plots together, the one added first is also
    drawn first (ie, addition of Graphics objects is not commutative).

    The following example produced a green-colored square with vertices
    at the points indicated.

    EXAMPLES:
        sage: from sage.groups.perm_gps.cubegroup import *
        sage: P = polygon_plot3d([[1,3,1],[2,3,1],[2,3,2],[1,3,2],[1,3,1]],rgbcolor=green)
    """
    rot = rotation_list(tilt,turn)
    points2 = [(xproj(x,y,z,rot), yproj(x,y,z,rot)) for (x,y,z) in points ]
    return polygon(points2, **kwargs)

###########################################################

#############  lots of "internal" utility plot functions #########



def inv_list(lst):
    """
    Input a list of ints 1, ..., m (in any order), outputs inverse perm.

    EXAMPLES:
        sage: from sage.groups.perm_gps.cubegroup import *
        sage: L = [2,3,1]
        sage: inv_list(L)
        [3, 1, 2]
    """
    return [lst.index(i)+1 for i in range(1,1+len(lst))]

### bottom layer L, F, R, B

def ldb(rgbclr):
    #square labeled 14
    return polygon([[-3,0],[-2,0], [-2,1], [-3,1]], rgbcolor=rgbclr)

def ld(rgbclr):
    #square labeled 15
    return polygon([[-2,0],[-1,0], [-1,1], [-2,1]], rgbcolor=rgbclr)

def lfd(rgbclr):
    #square labeled 16
    return polygon([[-1,0],[0,0], [0,1], [-1,1]], rgbcolor=rgbclr)

def fdl(rgbclr):
    #square labeled 22
    return polygon([[0,0],[1,0], [1,1], [0,1]], rgbcolor=rgbclr)

def fd(rgbclr):
    #square labeled 23
    return polygon([[1,0],[2,0], [2,1], [1,1]], rgbcolor=rgbclr)

def frd(rgbclr):
    #square labeled 24
    return polygon([[2,0],[3,0], [3,1], [2,1]], rgbcolor=rgbclr)

def rdf(rgbclr):
    #square labeled 30
    return polygon([[3,0],[4,0], [4,1], [3,1]], rgbcolor=rgbclr)

def rd(rgbclr):
    #square labeled 31
    return polygon([[4,0],[5,0], [5,1], [4,1]], rgbcolor=rgbclr)

def rbd(rgbclr):
    #square labeled 32
    return polygon([[5,0],[6,0], [6,1], [5,1]], rgbcolor=rgbclr)

def bdr(rgbclr):
    #square labeled 38
    return polygon([[6,0],[7,0], [7,1], [6,1]], rgbcolor=rgbclr)

def bd(rgbclr):
    #square labeled 39
    return polygon([[7,0],[8,0], [8,1], [7,1]], rgbcolor=rgbclr)

def bld(rgbclr):
    #square labeled 40
    return polygon([[8,0],[9,0], [9,1], [8,1]], rgbcolor=rgbclr)

### middle layer L,F,R, B

def lb(rgbclr):
    #square labeled 12
    return polygon([[-3,1],[-2,1], [-2,2], [-3,2]], rgbcolor=rgbclr)

def l_center(rgbclr):
    return polygon([[-2,1],[-1,1], [-1,2], [-2,2]], rgbcolor=rgbclr)

def lf(rgbclr):
    #square labeled 13
    return polygon([[-1,1],[0,1], [0,2], [-1,2]], rgbcolor=rgbclr)

def fl(rgbclr):
    #square labeled 20
    return polygon([[0,1],[1,1], [1,2], [0,2]], rgbcolor=rgbclr)

def f_center(rgbclr):
    return polygon([[1,1],[2,1], [2,2], [1,2]], rgbcolor=rgbclr)

def fr(rgbclr):
    #square labeled 21
    return polygon([[2,1],[3,1], [3,2], [2,2]], rgbcolor=rgbclr)

def rf(rgbclr):
    #square labeled 28
    return polygon([[3,1],[4,1], [4,2], [3,2]], rgbcolor=rgbclr)

def r_center(rgbclr):
    return polygon([[4,1],[5,1], [5,2], [4,2]], rgbcolor=rgbclr)

def rb(rgbclr):
    #square labeled 29
    return polygon([[5,1],[6,1], [6,2], [5,2]], rgbcolor=rgbclr)

def br(rgbclr):
    #square labeled 36
    return polygon([[6,1],[7,1], [7,2], [6,2]], rgbcolor=rgbclr)

def b_center(rgbclr):
    return polygon([[7,1],[8,1], [8,2], [7,2]], rgbcolor=rgbclr)

def bl(rgbclr):
    #square labeled 37
    return polygon([[8,1],[9,1], [9,2], [8,2]], rgbcolor=rgbclr)

## top layer L, F, R, B

def lbu(rgbclr):
    #square labeled 9
    return polygon([[-3,2],[-2,2], [-2,3], [-3,3]], rgbcolor=rgbclr)

def lu(rgbclr):
    #square labeled 10
    return polygon([[-2,2],[-1,2], [-1,3], [-2,3]], rgbcolor=rgbclr)

def luf(rgbclr):
    #square labeled 11
    return polygon([[-1,2],[0,2], [0,3], [-1,3]], rgbcolor=rgbclr)

def flu(rgbclr):
    #square labeled 17
    return polygon([[0,2],[1,2], [1,3], [0,3]], rgbcolor=rgbclr)

def fu(rgbclr):
    #square labeled 18
    return polygon([[1,2],[2,2], [2,3], [1,3]], rgbcolor=rgbclr)

def fur(rgbclr):
    #square labeled 19
    return polygon([[2,2],[3,2], [3,3], [2,3]], rgbcolor=rgbclr)

def ruf(rgbclr):
    #square labeled 25
    return polygon([[3,2],[4,2], [4,3], [3,3]], rgbcolor=rgbclr)

def ru(rgbclr):
    #square labeled 26
    return polygon([[4,2],[5,2], [5,3], [4,3]], rgbcolor=rgbclr)

def rub(rgbclr):
    #square labeled 27
    return polygon([[5,2],[6,2], [6,3], [5,3]], rgbcolor=rgbclr)

def bur(rgbclr):
    #square labeled 33
    return polygon([[6,2],[7,2], [7,3], [6,3]], rgbcolor=rgbclr)

def bu(rgbclr):
    #square labeled 34
    return polygon([[7,2],[8,2], [8,3], [7,3]], rgbcolor=rgbclr)

def bul(rgbclr):
    #square labeled 35
    return polygon([[8,2],[9,2], [9,3], [8,3]], rgbcolor=rgbclr)

# down face

def dlf(rgbclr):
    #square labeled 41
    return polygon([[0,-1],[1,-1], [1,0], [0,0]], rgbcolor=rgbclr)

def df(rgbclr):
    #square labeled 42
    return polygon([[1,-1],[2,-1], [2,0], [1,0]], rgbcolor=rgbclr)

def dfr(rgbclr):
    #square labeled 43
    return polygon([[2,-1],[3,-1], [3,0], [2,0]], rgbcolor=rgbclr)

def dl(rgbclr):
    #square labeled 44
    return polygon([[0,-2],[1,-2], [1,-1], [0,-1]], rgbcolor=rgbclr)

def d_center(rgbclr):
    return polygon([[1,-2],[2,-2], [2,-1], [1,-1]], rgbcolor=rgbclr)

def dr(rgbclr):
    #square labeled 45
    return polygon([[2,-2],[3,-2], [3,-1], [2,-1]], rgbcolor=rgbclr)

def dlb(rgbclr):
    #square labeled 46
    return polygon([[0,-3],[1,-3], [1,-2], [0,-2]], rgbcolor=rgbclr)

def db(rgbclr):
    #square labeled 47
    return polygon([[1,-3],[2,-3], [2,-2], [1,-2]], rgbcolor=rgbclr)

def drb(rgbclr):
    #square labeled 48
    return polygon([[2,-3],[3,-3], [3,-2], [2,-2]], rgbcolor=rgbclr)

# up face

def ufl(rgbclr):
    #square labeled 6
    return polygon([[0,3],[1,3], [1,4], [0,4]], rgbcolor=rgbclr)

def uf(rgbclr):
    #square labeled 7
    return polygon([[1,3],[2,3], [2,4], [1,4]], rgbcolor=rgbclr)

def urf(rgbclr):
    #square labeled 8
    return polygon([[2,3],[3,3], [3,4], [2,4]], rgbcolor=rgbclr)

def ul(rgbclr):
    #square labeled 4
    return polygon([[0,4],[1,4], [1,5], [0,5]], rgbcolor=rgbclr)

def u_center(rgbclr):
    return polygon([[1,4],[2,4], [2,5], [1,5]], rgbcolor=rgbclr)

def ur(rgbclr):
    #square labeled 5
    return polygon([[2,4],[3,4], [3,5], [2,5]], rgbcolor=rgbclr)

def ulb(rgbclr):
    #square labeled 1
    return polygon([[0,6],[1,6], [1,5], [0,5]], rgbcolor=rgbclr)

def ub(rgbclr):
    #square labeled 2
    return polygon([[1,6],[2,6], [2,5], [1,5]], rgbcolor=rgbclr)

def ubr(rgbclr):
    #square labeled 3
    return polygon([[2,6],[3,6], [3,5], [2,5]], rgbcolor=rgbclr)

####################################################

def index2singmaster(facet):
    """
    Translates index used (eg, 43) to Singmaster facet notation (eg, fdr).

    EXAMPLES:
        sage: from sage.groups.perm_gps.cubegroup import *
        sage: index2singmaster(41)
        'dlf'
    """
    if facet==1: return "ulb"
    if facet==2: return "ub"
    if facet==3: return "ubr"
    if facet==4: return "ul"
    if facet==5: return "ur"
    if facet==6: return "ufl"
    if facet==7: return "uf"
    if facet==8: return "urf"
    if facet==14: return "ldb"
    if facet==15: return "ld"
    if facet==16: return "lfd"
    if facet==12: return "lb"
    if facet==13: return "lf"
    if facet==9: return "lbu"
    if facet==10: return "lu"
    if facet==11: return "luf"
    if facet==17: return "flu"
    if facet==18: return "fu"
    if facet==19: return "fur"
    if facet==20: return "fl"
    if facet==21: return "fr"
    if facet==22: return "fdl"
    if facet==23: return "fd"
    if facet==24: return "frd"
    if facet==41: return "dlf"
    if facet==42: return "df"
    if facet==43: return "dfr"
    if facet==44: return "dl"
    if facet==45: return "dr"
    if facet==46: return "dlb"
    if facet==47: return "db"
    if facet==48: return "drb"
    if facet==33: return "bur"
    if facet==34: return "bu"
    if facet==35: return "bul"
    if facet==36: return "br"
    if facet==37: return "bl"
    if facet==38: return "bdr"
    if facet==39: return "bd"
    if facet==40: return "bld"
    if facet==25: return "ruf"
    if facet==26: return "ru"
    if facet==27: return "rub"
    if facet==28: return "rf"
    if facet==29: return "rb"
    if facet==30: return "rdf"
    if facet==31: return "rd"
    if facet==32: return "rbd"

def color_of_square(facet):
    """
    Returns the color the facet has in the solved state.

    EXAMPLES:
        sage: from sage.groups.perm_gps.cubegroup import *
        sage: color_of_square(41)
        'blue'
    """
    if facet== 1: return "lpurple"
    if facet== 2: return "lpurple"
    if facet== 3: return "lpurple"
    if facet== 4: return "lpurple"
    if facet== 5: return "lpurple"
    if facet== 6: return "lpurple"
    if facet== 7: return "lpurple"
    if facet== 8: return "lpurple"
    if facet== 9: return "yellow"
    if facet==10: return "yellow"
    if facet==11: return "yellow"
    if facet==12: return "yellow"
    if facet==13: return "yellow"
    if facet==14: return "yellow"
    if facet==15: return "yellow"
    if facet==16: return "yellow"
    if facet==41: return "blue"
    if facet==42: return "blue"
    if facet==43: return "blue"
    if facet==44: return "blue"
    if facet==45: return "blue"
    if facet==46: return "blue"
    if facet==47: return "blue"
    if facet==48: return "blue"
    if facet==33: return "orange"
    if facet==34: return "orange"
    if facet==35: return "orange"
    if facet==36: return "orange"
    if facet==37: return "orange"
    if facet==38: return "orange"
    if facet==39: return "orange"
    if facet==40: return "orange"
    if facet==25: return "green"
    if facet==26: return "green"
    if facet==27: return "green"
    if facet==28: return "green"
    if facet==29: return "green"
    if facet==30: return "green"
    if facet==31: return "green"
    if facet==32: return "green"
    if facet==17: return "red"
    if facet==18: return "red"
    if facet==19: return "red"
    if facet==20: return "red"
    if facet==21: return "red"
    if facet==22: return "red"
    if facet==23: return "red"
    if facet==24: return "red"

def cubie_centers(label):
    #  centers of the cubies on the F,U, R faces
    if label == 1: return [1/2,1/2,5/2] #ulb,
    if label == 2: return  [1/2,3/2,5/2] # ub,
    if label == 3: return  [1/2,5/2,5/2] # ubr,
    if label == 4: return  [3/2,1/2,5/2] # ul,
    if label == 5: return [3/2,5/2,5/2] # ur,
    if label == 6: return [5/2,1/2,5/2] # ufl,
    if label == 7: return [5/2,3/2,5/2] # uf,
    if label == 8: return [5/2,5/2,5/2] # urf
    if label == 17: return [5/2,1/2,5/2] # flu
    if label == 18: return [5/2,3/2,5/2] # fu
    if label == 19: return [5/2,5/2,5/2] # fur
    if label == 20: return [5/2,1/2,3/2] # fl
    if label == 21: return [5/2,5/2,3/2] # fr
    if label == 22: return [5/2,1/2,1/2] # fdl
    if label == 23: return [5/2,3/2,1/2] # fd
    if label == 24: return [5/2,5/2,1/2] # frd
    if label == 25: return [5/2,5/2,5/2] #rfu,
    if label == 26: return  [3/2,5/2,5/2] # ru,
    if label == 27: return  [1/2,5/2,5/2] # rub
    if label == 28: return  [5/2,5/2,3/2] # rf,
    if label == 29: return [1/2,5/2,3/2] # rb,
    if label == 30: return [5/2,5/2,1/2] # rdf
    if label == 31: return [3/2,5/2,1/2] # rd,
    if label == 32: return [1/2,5/2,1/2] #rbd,

def cubie_colors(label,state0):
    #  colors of the cubies on the F,U, R faces
    clr_any = white
    state = inv_list(state0)
    if label == 1: return [clr_any, eval(color_of_square(state[1-1])), clr_any] #ulb,
    if label == 2: return  [clr_any,eval(color_of_square(state[2-1])),clr_any] # ub,
    if label == 3: return  [clr_any, eval(color_of_square(state[3-1])), eval(color_of_square(state[27-1]))] # ubr,
    if label == 4: return  [clr_any, eval(color_of_square(state[4-1])), clr_any] # ul,
    if label == 5: return [clr_any, eval(color_of_square(state[5-1])), eval(color_of_square(state[26-1]))] # ur,
    if label == 6: return [eval(color_of_square(state[17-1])), eval(color_of_square(state[6-1])), clr_any] # ufl,
    if label == 7: return [eval(color_of_square(state[18-1])), eval(color_of_square(state[7-1])), clr_any] # uf,
    if label == 8: return [eval(color_of_square(state[19-1])), eval(color_of_square(state[8-1])), eval(color_of_square(state[25-1]))] # urf,
    if label == 17: return [eval(color_of_square(state[17-1])), eval(color_of_square(state[6-1])), clr_any] # flu
    if label == 18: return [eval(color_of_square(state[18-1])), eval(color_of_square(state[7-1])), clr_any] # fu
    if label == 19: return [eval(color_of_square(state[19-1])), eval(color_of_square(state[8-1])), eval(color_of_square(state[25-1]))] # fur
    if label == 20: return [eval(color_of_square(state[20-1])), clr_any, clr_any] # fl
    if label == 21: return [eval(color_of_square(state[21-1])), clr_any, eval(color_of_square(state[28-1]))] # fr
    if label == 22: return [eval(color_of_square(state[22-1])), clr_any, clr_any] # fdl
    if label == 23: return [eval(color_of_square(state[23-1])), clr_any, clr_any] # fd
    if label == 24: return [eval(color_of_square(state[24-1])), clr_any, eval(color_of_square(state[30-1]))] # frd
    if label == 25: return [eval(color_of_square(state[19-1])),eval(color_of_square(state[8-1])),eval(color_of_square(state[25-1]))]  #rfu,
    if label == 26: return [clr_any,eval(color_of_square(state[5-1])),eval(color_of_square(state[26-1]))] # ru,
    if label == 27: return [clr_any,eval(color_of_square(state[3-1])),eval(color_of_square(state[27-1]))] # rub,
    if label == 28: return [eval(color_of_square(state[21-1])),clr_any,eval(color_of_square(state[28-1]))] # rf,
    if label == 29: return [clr_any,clr_any,eval(color_of_square(state[29-1]))] # rb,
    if label == 30: return [eval(color_of_square(state[24-1])),clr_any,eval(color_of_square(state[30-1]))] # rdf,
    if label == 31: return [clr_any,clr_any,eval(color_of_square(state[31-1]))] # rd,
    if label == 32: return [clr_any,clr_any,eval(color_of_square(state[32-1]))] #rbd,

def plot3d_cubie(cnt, clrs):
    """
    Plots the front, up and right face of a cubie centered at cnt and rgbcolors given by clrs (in the
    order FUR).

    Type P.show() to view.

    EXAMPLES:
        sage: from sage.groups.perm_gps.cubegroup import *
        sage: clrF = blue; clrU = red; clrR = green
        sage: P = plot3d_cubie([1/2,1/2,1/2],[clrF,clrU,clrR])

    """
    x = cnt[0]-1/2; y = cnt[1]-1/2; z = cnt[2]-1/2
    #ptsD = [[x+0,y+0,0+z],[x+1,y+0,0+z],[x+1,y+1,0+z],[x+0,y+1,0+z],[x+0,y+0,0+z]]
    ptsF = [[x+1,y+0,0+z],[x+1,y+1,0+z],[x+1,y+1,1+z],[x+1,y+0,1+z],[x+1,y+0,0+z]]
    #ptsB = [[x+0,y+0,0+z],[x+0,y+1,0+z],[x+0,y+1,1+z],[x+0,y+0,1+z],[x+0,y+0,0+z]]
    ptsU = [[x+0,y+0,1+z],[x+1,y+0,1+z],[x+1,y+1,1+z],[x+0,y+1,1+z],[x+0,y+0,1+z]]
    #ptsL = [[x+0,y+0,0+z],[x+1,y+0,0+z],[x+1,y+0,1+z],[x+0,y+0,1+z],[x+0,y+0,0+z]]
    ptsR = [[x+0,y+1,0+z],[x+1,y+1,0+z],[x+1,y+1,1+z],[x+0,y+1,1+z],[x+0,y+1,0+z]]
    PR = polygon_plot3d(ptsR,rgbcolor=clrs[2])
    PU = polygon_plot3d(ptsU,rgbcolor=clrs[1])
    PF = polygon_plot3d(ptsF,rgbcolor=clrs[0])
    P = PR+PF+PU
    P.axes(show=False)
    return P


####################### end of "internal" utility plot functions  #################


class CubeGroup(PermutationGroup_generic):
    """
    A python class to help compute Rubik's cube group actions.
    If G denotes the cube group then it may be regarded as a subgroup
    of SymmetricGroup(48), where the 48 facets are labeled as follows.


                             +--------------+
                             |  1    2    3 |
                             |  4   Up    5 |
                             |  6    7    8 |
              +--------------+--------------+--------------+--------------+
              |  9   10   11 | 17   18   19 | 25   26   27 | 33   34   35 |
              | 12  Left  13 | 20 Front  21 | 28 Right  29 | 36  Back  37 |
              | 14   15   16 | 22   23   24 | 30   31   32 | 38   39   40 |
              +--------------+--------------+--------------+--------------+
                             | 41   42   43 |
                             | 44  Down  45 |
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
        self.__gens = [B,D,F,L,R,U]
        self._group = PermutationGroup([B,D,F,L,R,U])
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

    def B(self):
        G = self.group()
	g = G(self.gens()[0])
        return g

    def D(self):
	G = self.group()
	g = G(self.gens()[1])
        return g

    def F(self):
        G = self.group()
	g = G(self.gens()[2])
        return g

    def L(self):
        G = self.group()
	g = G(self.gens()[3])
        return g

    def R(self):
        G = self.group()
	g = G(self.gens()[4])
        return g

    def U(self):
        G = self.group()
	g = G(self.gens()[5])
        return g


    def facets(self):
        r"""
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
        r"""
        Returns the dictionary of faces created by the effect of the
        move mv, which is a string of the form $X^a*Y^b*...$, where
        $X$, $Y$, ... are in $\{R,L,F,B,U,D\}$ and $a,b, ...$ are
        integers.  We call this ordering of the faces the
        "BDFLRU, L2R, T2B ordering".

        EXAMPLES:
            sage: rubik = CubeGroup()

        Now type \code{rubik.faces("")} for the dictionary of the solved state
        and \code{rubik.faces("R*L")} for the dictionary of the state obtained
        after making the move R followed by L.

        """
        fcts = self.move(mv)[1]
        faceR = [[fcts[24],fcts[25],fcts[26]],[fcts[27],0,fcts[28]],[fcts[29],fcts[30],fcts[31]]]
        faceL = [[fcts[8],fcts[9],fcts[10]],[fcts[11],0,fcts[12]],[fcts[13],fcts[14],fcts[15]]]
        faceU = [[fcts[0],fcts[1],fcts[2]],[fcts[3],0,fcts[4]],[fcts[5],fcts[6],fcts[7]]]
        faceD = [[fcts[40],fcts[41],fcts[42]],[fcts[43],0,fcts[44]],[fcts[45],fcts[46],fcts[47]]]
        faceF = [[fcts[16],fcts[17],fcts[18]],[fcts[19],0,fcts[20]],[fcts[21],fcts[22],fcts[23]]]
        faceB = [[fcts[32],fcts[33],fcts[34]],[fcts[35],0,fcts[36]],[fcts[37],fcts[38],fcts[39]]]
        return {'right':faceR,'left':faceL,'up':faceU,'down':faceD,'front':faceF,'back':faceB}

    def move(self,mv):
        r"""
        Returns the group element and the reordered list of facets, as moved by
	the list mv (read left-to-right)

	INPUT: mv is a string of the form X^a*Y^b*...",
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
	    if len(M[i])==1:
	        M[i] = [M[i][0],"1"]
	#print M
	for i in range(n):
	    x = M[i][0]
            if x == "R":   h = self.R()
	    elif x == "L": h = self.L()
	    elif x == "U": h = self.U()
	    elif x == "D": h = self.D()
	    elif x == "F": h = self.F()
	    elif x == "B": h = self.B()
	    else: h = G(1)
	    e = M[i][1]
	    if e=="1": g = g*h
	    if e=="2": g = g*h*h
	    if e=="3": g = g*h*h*h
            if e=="(-1)": g = g*h*h*h
	pos = [g(i) for i in fcts]
	return [g,pos]

    def display2d(self,mv):
        r"""
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

    def legal(self,state,mode="quiet"):
        r"""
        Returns 1 (true) if the dictionary \code{state} (in the same format as
        returned by the faces method) represents a legal position (or state) of
        the Rubik's cube. Returns 0 (false) otherwise.

        EXAMPLES:
            sage: rubik = CubeGroup()
            sage: G = rubik.group()
            sage: r0 = rubik.faces("")
            sage: r1 = {'back': [[33, 34, 35], [36, 0, 37], [38, 39, 40]], 'down': [[41, 42, 43], [44, 0, 45], [46, 47, 48]],'front': [[17, 18, 19], [20, 0, 21], [22, 23, 24]],'left': [[9, 10, 11], [12, 0, 13], [14, 15, 16]],'right': [[25, 26, 27], [28, 0, 29], [30, 31, 32]],'up': [[1, 2, 3], [4, 0, 5], [6, 8, 7]]}
            sage: rubik.legal(r0)
            1
            sage: rubik.legal(r0,"verbose")
            (1, ())
            sage: rubik.legal(r1)
            0

        """
        state_facets = []
        keyss = state.keys()
        keyss.sort()
        for k in keyss:
            r = state[k][0]+state[k][1]+state[k][2]
            r.remove(0)
            state_facets = state_facets + r
        state0 = self.faces("")
        state0_facets = []
        keyss = state0.keys()
        keyss.sort()
        for k in keyss:
            r = state0[k][0]+state0[k][1]+state0[k][2]
            r.remove(0)
            state0_facets = state0_facets + r
        p1 = [state0_facets.index(x) for x in range(1,49)]
        p2 = [state_facets[j] for j in p1]
        #if mode != "quiet":
        #    print "list (BDFLRU, L2R, T2B ordering):",p2
        G = self.group()
        try:
            g = G(p2)
        except TypeError:
            if mode != "quiet":
                return 0,G([()])
            else:
                return 0
        else:
            if mode != "quiet":
                #print "cube group element:"
                return 1,g
            else:
                return 1

    def solve(self,state):
        r"""
        Solves the cube in the \code{state}, given as a dictionary as
        in \code{legal}.  This uses GAP's \code{EpimorphismFromFreeGroup}
        and \code{PreImagesRepresentative}.

        WARNING: This is currently evidently broken.

        EXAMPLES:
            sage: rubik = CubeGroup()
            sage: state = rubik.faces("R*U")
            sage: rubik.solve(state)  # long time; *omputationally intensive* even for simple moves
            'R*U'

        You can also check this using \code{word_problem} method (eg, G = rubik.group();
        g = G("(3,38,43,19)(5,36,45,21)(8,33,48,24)(25,27,32,30)(26,29,31,28)");
        g.word_problem([b,d,f,l,r,u]), though the output will be less intuitive).

        """
        from sage.groups.perm_gps.permgroup import PermutationGroup
        from sage.interfaces.all import gap
        rubik = self
        leg = rubik.legal(state,"verbose")
        if not(leg[0]):
            return "Illegal or syntactically incorrect state. No solution."
        G = rubik.group()
        g = leg[1]
        hom = G._gap_().EpimorphismFromFreeGroup()
        soln = hom.PreImagesRepresentative(gap(str(g)))
        sol1 = str(soln).replace("x1","B")
        sol2 = sol1.replace("x2","D")
        sol3 = sol2.replace("x3","F")
        sol4 = sol3.replace("x4","L")
        sol5 = sol4.replace("x5","R")
        sol6 = sol5.replace("x6","U")
        return sol6

### default plot method. Goes into global namespace

def plot_cube(mv,title=True):
    """
    Input the move mv, as a string in the Singmaster notation,
    and output the 2-d plot of the cube in that state.

    Type P.show() to display any of the plots below.

    EXAMPLES:
        sage: P = plot_cube("R^2*U^2*R^2*U^2*R^2*U^2", title = False)
        sage: # (R^2U^2)^3  permutes 2 pairs of edges (uf,ub)(fr,br)
        sage: P = plot_cube("R*L*D^2*B^3*L^2*F^2*R^2*U^3*D*R^3*D^2*F^3*B^3*D^3*F^2*D^3*R^2*U^3*F^2*D^3")
        sage: # the superflip (in 20f* moves)
        sage: P = plot_cube("U^2*F*U^2*L*R^(-1)*F^2*U*F^3*B^3*R*L*U^2*R*D^3*U*L^3*R*D*R^3*L^3*D^2")
        sage: # "superflip+4 spot" (in 26q* moves)
    """
    rubik = CubeGroup()
    state = rubik.move(mv)[1]
    #print state
    str_colors = [index2singmaster(state[x])+"("+str(color_of_square(x+1))+")" for x in range(48)]
    colors = [eval(p) for p in str_colors]
    centers = u_center(lpurple)+f_center(red)+l_center(yellow)+r_center(green)+d_center(blue)+b_center(orange)
    clrs = centers+sum(colors)
    clrs.axes(show=False)
    if title == True:
        t = text('sagemath.org', (7.8,-3.5),rgbcolor=lgrey)
        P = clrs+t
        P.axes(show=False)
        return P
    return clrs

####### 3d plot method. Goes into global namespace.

def plot3d_cube(mv,title=True):
    """
    Displays F,U,R faces of the cube after the given move mv, where mv is a string in the Singmaster notation.
    Mostly included for the purpose of drawing pictures and checking moves.

    The first one below is "superflip+4 spot" (in 26q* moves) and the second one is the
    superflip (in 20f* moves). Type show(P) to view them.

    EXAMPLES:
        sage: P = plot3d_cube("U^2*F*U^2*L*R^(-1)*F^2*U*F^3*B^3*R*L*U^2*R*D^3*U*L^3*R*D*R^3*L^3*D^2")
        sage: P = plot3d_cube("R*L*D^2*B^3*L^2*F^2*R^2*U^3*D*R^3*D^2*F^3*B^3*D^3*F^2*D^3*R^2*U^3*F^2*D^3")
    """
    rubik = CubeGroup()
    state = rubik.move(mv)[1]
    clr_any = white
    shown_labels = range(1,9)+range(17,33)
    clr = [color_of_square(state[c-1]) for c in shown_labels]
    cubiesR = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in [32,31,30,29,28,27,26,25]]
    cubeR = sum(cubiesR)
    cubiesU = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in range(1,9)]
    cubeU = sum(cubiesU)
    cubiesF = [plot3d_cubie(cubie_centers(c),cubie_colors(c,state)) for c in [22,23,24,20,21]]
    cubeF = sum(cubiesF)
    centerR =  polygon_plot3d([[1,3,1],[2,3,1],[2,3,2],[1,3,2],[1,3,1]],rgbcolor=green)
    centerF =  polygon_plot3d([[3,1,1],[3,2,1],[3,2,2],[3,1,2],[3,1,1]],rgbcolor=red)
    centerU =  polygon_plot3d([[1,1,3],[1,2,3],[2,2,3],[2,1,3],[1,1,3]],rgbcolor=lpurple)
    centers = centerF+centerR+centerU
    P = cubeR+cubeF+cubeU+centers
    P.axes(show=False)
    if title == True:
        t1 = text('Up, Front, and Right faces. '   , (-0.2,-2.5))
        t2  = text('      sagemath.org', (0.8,-3.1),rgbcolor=lgrey)
        t3 = text("     ",(3.5,0),rgbcolor=white)
        P = P+t1+t2+t3
        P.axes(show=False)
        return P
    return P
