r"""
Plane Partitions

AUTHORS:

- Jang Soo Kim (2016): Initial implementation
- Jessica Striker (2016): Added additional methods
"""
#*****************************************************************************
#       Copyright (C) 2016 Jang Soo Kim <jangsookim@skku.edu>,
#                     2016 Jessica Striker <jessicapalencia@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful, but
#    WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.list_clone import ClonableArray
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.posets.posets import Poset
from sage.rings.integer import Integer
from sage.misc.all import prod
from sage.combinat.tableau import Tableau

class PlanePartition(ClonableArray):
    r"""
    Return a plane partition.

    INPUT:

    - ``A`` -- a list of lists which represents a tableau.

    - ``box_size`` -- a list [A,B,C] of 3 positive integers (default: None)
      If this is not given, it is determined by ``A``.

    OUTPUT:

    The plane partition whose tableau representation is ``A``.

    EXAMPLES:

        sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]]); PP
        The plane partition [[4, 3, 3, 1], [2, 1, 1], [1, 1]].

    """
    
    __metaclass__ = InheritComparisonClasscallMetaclass
    
    @staticmethod
    def __classcall_private__(cls, A, box_size=None):
        if box_size is None:
            box_size = (len(A), len(A[0]), A[0][0])
        return PlanePartitions(box_size)(A)

    def __init__(self, parent, A):
        ClonableArray.__init__(self, parent, A)
        self._max_x = parent._box[0]
        self._max_y = parent._box[1]
        self._max_z = parent._box[2]
        
    def check(self):
        pass

    def _repr_(self):
        return "The plane partition {}.".format(list(self))

    def to_tableau(self):
        r"""
        Return the tableau class of ``self``.
        
        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.to_tableau()
            [[4, 3, 3, 1], [2, 1, 1], [1, 1]]
        
        """
        return Tableau(self)

    def z_tableau(self):
        r"""
        Return the projection of ``self`` in the z direction.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.z_tableau()
            [[4, 3, 3, 1], [2, 1, 1, 0], [1, 1, 0, 0]]

        """
        Z = [[0 for i in range(self._max_y)] for j in range(self._max_x)]
        for C in self.cells():
            Z[C[0]][C[1]] += 1
        return Z

    def y_tableau(self):
        r"""
        Return the projection of ``self`` in the y direction.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.y_tableau()
            [[4, 3, 2], [3, 1, 0], [3, 0, 0], [1, 0, 0]]

        """
        Y = [[0 for i in range(self._max_x)] for j in range(self._max_z)]
        for C in self.cells():
            Y[C[2]][C[0]] += 1
        return Y

    def x_tableau(self):
        r"""
        Return the projection of ``self`` in the x direction.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.x_tableau()
            [[3, 2, 1, 1], [3, 1, 1, 0], [2, 1, 1, 0], [1, 0, 0, 0]]


        """
        X = [[0 for i in range(self._max_z)] for j in range(self._max_y)]
        for C in self.cells():
            X[C[1]][C[2]] += 1
        return X

    def cells(self):
        r"""
        Return the list of cells inside the plane partition.

        EXAMPLES:
            
            sage: PP = PlanePartition([[3,1],[2]])
            sage: PP.cells()
            [[0, 0, 0], [0, 0, 1], [0, 0, 2], [0, 1, 0], [1, 0, 0], [1, 0, 1]]
        """

        L = []
        for r in range(len(self)):
            for c in range(len(self[r])):
                for h in range(self[r][c]):
                    L.append([r,c,h])
        return L

    def _repr_diagram(self, show_box = False):
        x=self._max_x
        y=self._max_y
        z=self._max_z
        drawing = [[" " for i in range(2*x+2*y+1)] for j in range(x+y+z+1)]
        def add_topside(i,j,k):
            drawing[z+1+i+j-k][2*x-1-2*i+2*j] = "/"
            drawing[z+1+i+j-k][2*x+1-2*i+2*j] = "\\"
            drawing[z+2+i+j-k][2*x-1-2*i+2*j] = "\\"
            drawing[z+2+i+j-k][2*x+1-2*i+2*j] = "/"
        def add_rightside(i,j,k):
            drawing[z+1+i+j-k][2*x-2-2*i+2*j] = "|"
            drawing[z+1+i+j-k][2*x-1-2*i+2*j] = "/"
            drawing[z+0+i+j-k][2*x-1-2*i+2*j] = "/"
            drawing[z+0+i+j-k][2*x-0-2*i+2*j] = "|"
        def add_leftside(i,j,k):
            drawing[z+0+i+j-k][2*x-0-2*i+2*j] = "|"
            drawing[z+0+i+j-k][2*x+1-2*i+2*j] = "\\"
            drawing[z+1+i+j-k][2*x+2-2*i+2*j] = "|"
            drawing[z+1+i+j-k][2*x+1-2*i+2*j] = "\\"
        for r in range(len(self.z_tableau())):
            for c in range(len(self.z_tableau()[r])):
                if self.z_tableau()[r][c] > 0 or show_box:
                    add_topside(r,c,self.z_tableau()[r][c])
        for r in range(len(self.y_tableau())):
            for c in range(len(self.y_tableau()[r])):
                if self.y_tableau()[r][c] > 0 or show_box:
                    add_rightside(c,self.y_tableau()[r][c],r)
        for r in range(len(self.x_tableau())):
            for c in range(len(self.x_tableau()[r])):
                if self.x_tableau()[r][c] > 0 or show_box:
                    add_leftside(self.x_tableau()[r][c],r,c)
                    
        ret = ""
        for row in drawing:
            D = ""
            for s in row:
                D = D + s
            ret += D + "\n"
        return ret
        
    def _ascii_art_(self):
        from sage.typeset.ascii_art import AsciiArt
        return AsciiArt(self._repr_diagram().splitlines(), baseline=0)

    def pp(self, show_box=False):
        r"""
        Prints a pretty print of the plane partition.

        INPUT:

        ``show_box`` -- boolean (default:False).
        If True, it also shows the visible tiles on the xy-, yz-, zx-planes.

        OUTPUT:

        A pretty print of the plane partition.

        EXAMPLES:

            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]]).pp()
            <BLANKLINE>              
                 / \       
                |\ /|      
                |\|/ \     
               / \|\ / \   
              |\ /|\|\ /|  
             / \|/ \|\|/|  
            |\ / \ / \|/ \ 
             \|\ /|\ /|\ /|
               \|/ \|/ \|/ 
            <BLANKLINE>  
            <BLANKLINE>                             
            <BLANKLINE>              
            sage: PlanePartition([[4,3,3,1],[2,1,1],[1,1]]).pp(True)
            <BLANKLINE>                             
                 / \       
               /|\ /|\     
             /|/|\|/ \|\   
            |/|/ \|\ / \|\ 
            |/|\ /|\|\ /|\|
            |/ \|/ \|\|/|\|
            |\ / \ / \|/ \|
             \|\ /|\ /|\ /|
               \|/ \|/ \|/ 
                 \ / \ /   
                   \ /     
            <BLANKLINE>              
        """
        print self._repr_diagram(show_box)

    def latex(self, show_box = False, colors = ["white","lightgray","darkgray"]):
        r"""
        Returns a latex code which uses TikZ package to draw the plane partition.

        INPUT:

        - ``show_box`` -- boolean (default:False).
          If True, it also shows the visible tiles on the xy-, yz-, zx-planes.

        - ``colors`` -- a list [A,B,C] of 3 strings representing colors (default:
          ["white","lightgray","darkgray"])

        OUTPUT:
        A latex code for drawing the plane partition.

        EXAMPLES:

            sage: PlanePartition([[1]]).latex()
            <BLANKLINE>  
            \begin{tikzpicture}
            \draw[fill=white,shift={(210:0)},shift={(-30:0)},shift={(90:1)}]
            (0,0)--(-30:1)--(0,-1)--(210:1)--(0,0);
            \draw[fill=darkgray,shift={(210:0)},shift={(-30:1)},shift={(90:0)}]
            (0,0)--(210:1)--(150:1)--(0,1)--(0,0);
            \draw[fill=lightgray,shift={(210:1)},shift={(-30:0)},shift={(90:0)}]
            (0,0)--(0,1)--(30:1)--(-30:1)--(0,0);
            \end{tikzpicture}
            <BLANKLINE>  
        """

        x=self._max_x
        y=self._max_y
        z=self._max_z
        print "\n\\begin{tikzpicture}"
        def add_topside(i,j,k):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(-30:1)--(0,-1)--(210:1)--(0,0);".format(colors[0],i,j,k)
        def add_leftside(i,j,k):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(0,1)--(30:1)--(-30:1)--(0,0);".format(colors[1],i,j,k)
        def add_rightside(i,j,k):
            return "\\draw[fill={},shift={{(210:{})}},shift={{(-30:{})}},shift={{(90:{})}}]\n(0,0)--(210:1)--(150:1)--(0,1)--(0,0);".format(colors[2],i,j,k)
        for r in range(len(self.z_tableau())):
            for c in range(len(self.z_tableau()[r])):
                if self.z_tableau()[r][c] > 0 or show_box:
                    print add_topside(r,c,self.z_tableau()[r][c])
        for r in range(len(self.y_tableau())):
            for c in range(len(self.y_tableau()[r])):
                if self.y_tableau()[r][c] > 0 or show_box:
                    print add_rightside(c,self.y_tableau()[r][c],r)
        for r in range(len(self.x_tableau())):
            for c in range(len(self.x_tableau()[r])):
                if self.x_tableau()[r][c] > 0 or show_box:
                    print add_leftside(self.x_tableau()[r][c],r,c)
        print "\\end{tikzpicture}\n"

    def show(self, show_box = None, colors = ["white","lightgray","darkgray"]):
        r"""
        Plots the plane partition.

        INPUT:

        - ``show_box`` -- boolean (default:False).
          If True, it also shows the visible tiles on the xy-, yz-, zx-planes.

        - ``colors`` -- a list [A,B,C] of 3 strings representing colors 
          (default: ["white","lightgray","darkgray"])       

        OUTPUT:

        A diagram of the plane partition.

        """
        if show_box is None:
            show_box = SHOW_BOX
        x=self._max_x
        y=self._max_y
        z=self._max_z
        Uside = [[0,0], [cos(-pi/6),sin(-pi/6)], [0,-1],[cos(7*pi/6),sin(7*pi/6)]]
        Lside = [[0,0], [cos(-pi/6),sin(-pi/6)], [cos(pi/6),sin(pi/6)],[0,1]]
        Rside = [[0,0], [0,1],[cos(5*pi/6),sin(5*pi/6)], [cos(7*pi/6),sin(7*pi/6)]]
        Xdir = [cos(7*pi/6),sin(7*pi/6)]
        Ydir = [cos(-pi/6),sin(-pi/6)]
        Zdir = [0,1]
        def move(side,i,j,k):
            T = []
            for P in side:
                T.append([P[0]+i*Xdir[0]+j*Ydir[0]+k*Zdir[0],P[1]+i*Xdir[1]+j*Ydir[1]+k*Zdir[1]])
            return T            
        def add_topside(i,j,k):
            return polygon(move(Uside,i,j,k),edgecolor="black",color=colors[0])
        def add_leftside(i,j,k):
            return polygon(move(Lside,i,j,k),edgecolor="black",color=colors[1])
        def add_rightside(i,j,k):
            return polygon(move(Rside,i,j,k),edgecolor="black",color=colors[2])
        TP = plot([])
        for r in range(len(self.z_tableau())):
            for c in range(len(self.z_tableau()[r])):
                if self.z_tableau()[r][c] > 0 or show_box:
                    TP += add_topside(r,c,self.z_tableau()[r][c])
        for r in range(len(self.y_tableau())):
            for c in range(len(self.y_tableau()[r])):
                if self.y_tableau()[r][c] > 0 or show_box:
                    TP += add_rightside(c,self.y_tableau()[r][c],r)
        for r in range(len(self.x_tableau())):
            for c in range(len(self.x_tableau()[r])):
                if self.x_tableau()[r][c] > 0 or show_box:
                    TP += add_leftside(self.x_tableau()[r][c],r,c)
        return TP.show(axes=False)

    def complement(self,tableau_only=False):
        r"""
        Return the complement of ``self``.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.complement()
            The plane partition [[4, 4, 3, 3], [4, 3, 3, 2], [3, 1, 1, 0]].
            sage: PP.complement(True)
            [[4, 4, 3, 3], [4, 3, 3, 2], [3, 1, 1, 0]]

        """
        A = self._max_x
        B = self._max_y
        C = self._max_z
        T = [[C for i in range(B)] for j in range(A)]
        for r in range(A):
            for c in range(B):
                T[A-1-r][B-1-c] = C - self.z_tableau()[r][c]
        if tableau_only:
            return T
        else:
            return type(self)(self.parent(), T)

    def transpose(self,tableau_only=False):
        r"""
        Return the transpose of ``self``.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.transpose()
            The plane partition [[4, 2, 1], [3, 1, 1], [3, 1, 0], [1, 0, 0]].
            sage: PP.transpose(True)
            [[4, 2, 1], [3, 1, 1], [3, 1, 0], [1, 0, 0]]
        """
        T = [[0 for i in range(self._max_x)] for j in range(self._max_y)]
        for r in range(len(self.z_tableau())):
            for c in range(len(self.z_tableau()[r])):
                T[c][r] = self.z_tableau()[r][c]
        if tableau_only:
            return T
        else:
            return type(self)(self.parent(), T)

    def is_SPP(self): 
        r"""
        Check whether ``self`` is a symmetric plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SPP()
            False
            sage: PP = PlanePartition([[3,3,2],[3,3,2],[2,2,2]])
            sage: PP.is_SPP()
            True
        """
        for r in range(len(self.z_tableau())):
            for c in range(r,len(self.z_tableau()[r])):
                if self.z_tableau()[r][c] != self.z_tableau()[c][r]:
                    return False
        return True

    def is_CSPP(self): 
        r"""
        Check whether ``self`` is a cyclically symmetric plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSPP()
            False
            sage: PP = PlanePartition([[3,2,2],[3,1,0],[1,1,0]])
            sage: PP.is_CSPP()
            True
        """
        if self.z_tableau() == self.y_tableau():
            return True
        return False

    def is_TSPP(self):    
        r"""
        Check whether ``self`` is a totally symmetric plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TSPP()
            False
            sage: PP = PlanePartition([[3,3,3],[3,3,2],[3,2,1]])
            sage: PP.is_TSPP()
            True
        """
        if self.is_CSPP() and self.is_SPP():
            return True
        return False

    def is_SCPP(self):  
        r"""
        Check whether ``self`` is a self-complementary plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SCPP()
            False
            sage: PP = PlanePartition([[4,4,4,4],[4,4,2,0],[4,2,0,0],[0,0,0,0]])
            sage: PP.is_SCPP()
            True
        """
        if self.z_tableau() == self.complement(True):
            return True
        return False

    def is_TCPP(self):    
        r"""
        Check whether ``self`` is a transpose-complementary plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,4,2,1],[4,2,0,0],[2,0,0,0]])
            sage: PP.is_TCPP()
            True
        """
        if self.transpose(True) == self.complement(True):
            return True
        return False

    def is_SSCPP(self):  
        r"""
        Check whether ``self`` is a symmetric, self-complementary plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_SSCPP()
            False
            sage: PP = PlanePartition([[4,3,3,2],[3,2,2,1],[3,2,2,1],[2,1,1,0]])
            sage: PP.is_SSCPP()
            True
        """
        if self.is_SPP() and self.is_SCPP():
            return True
        return False

    def is_CSTCPP(self):    
        r"""
        Check whether ``self`` is a cyclically symmetric and transpose-complementary plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSTCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,3,2,1],[3,2,1,0],[2,1,0,0]])
            sage: PP.is_CSTCPP()
            True
        """
        if self.is_CSPP() and self.is_TCPP():
            return True
        return False

    def is_CSSCPP(self):
        r"""
        Check whether ``self`` is a cyclically symmetric and self-complementary 
        plane partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_CSSCPP()
            False
            sage: PP = PlanePartition([[4,4,4,1],[3,3,2,1],[3,2,1,1],[3,0,0,0]])
            sage: PP.is_CSSCPP()
            True
        """
        if self.is_CSPP() and self.is_SCPP():
            return True
        return False

    def is_TSSCPP(self):
        r"""
        Check whether ``self`` is a totally symmetric self-complementary plane 
        partition.

        EXAMPLES:

            sage: PP = PlanePartition([[4,3,3,1],[2,1,1],[1,1]])
            sage: PP.is_TSSCPP()
            False
            sage: PP = PlanePartition([[4,4,3,2],[4,3,2,1],[3,2,1,0],[2,1,0,0]])
            sage: PP.is_TSSCPP()
            True
        """
        if self.is_TSPP() and self.is_SCPP():
            return True
        return False
    
class PlanePartitions(UniqueRepresentation,Parent):
    r"""
    Class of all plane partitions inside a rectangular box of given side lengths.

    A plane partition is a stack of cubes in the positive orthant.

    INPUT:

    - `box_size` -- a triple of positive integers indicating the size of the box 
      containing the plane partition.

    EXAMPLES:

    This will create an instance to manipulate the alternating sign
    matrices of size 3::

        sage: P = PlanePartitions((4,3,2))
        sage: P
        Plane partitions inside a 4 x 3 x 2 box.
        sage: P.cardinality()
        490
    """

    def __init__(self, box_size):
        self._box = box_size
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        return "Plane partitions inside a {} x {} x {} box.".format(self._box[0],self._box[1],self._box[2])

    def __iter__(self):
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]
        for T in SemistandardTableaux([B for i in range(A)],max_entry=C+A):
            PP = [[0 for i in range(B)] for j in range(A)]
            for r in range(A):
                for c in range(B):
                    PP[A-1-r][B-1-c] = T[r][c] - r  - 1          
            yield self.element_class(self, PP)
        
    Element = PlanePartition

    def cardinality(self):
        r"""
        Return the cardinality of ``self``.

        The number of plane partitions inside an a x b x c box is equal to

        .. MATH::

            \prod_{i=1}^{a} \prod_{j=1}^{b} \prod_{k=1}^{c} \frac{i+j+k-1}{i+j+k-2} 

        EXAMPLES::

            sage: PlanePartitions((4,3,5)).cardinality()
            116424
        """
        A = self._box[0]
        B = self._box[1]
        C = self._box[2]        
        return Integer(prod( [ Integer(i+j+k-1)/Integer(i+j+k-2) for i in range(1,A+1) for j in range(1,B+1) for k in range(1,C+1)] ))
    
    def random_element(self):
        r"""
        Return a uniformly random element of ``self``. 
        This uses the random_order_ideal method in posets.

        EXAMPLES::

            sage: PlanePartitions((4,3,5)).random_element()
            The plane partition [[4, 3, 3], [4, 0, 0], [2, 0, 0], [0, 0, 0]].
        """
        def leq(thing1, thing2):
            if all(thing1[i] <= thing2[i] for i in range(len(thing1))):
                return True
            return False
        elem = [(i,j,k) for i in range(self._box[0]) for j in range(self._box[1]) for k in range(self._box[2])]
        myposet = Poset((elem,leq))
        R = myposet.random_order_ideal()
        Z = [[0 for i in range(self._box[1])] for j in range(self._box[0])]
        for C in R:
            Z[C[0]][C[1]] += 1
        return PlanePartition(Z,self._box)    













