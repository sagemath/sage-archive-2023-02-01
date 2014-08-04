r"""
Link class
"""
#*****************************************************************************
#  Copyright (C) 2014
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.groups.free_group import FreeGroupElement
from sage.groups.braid import Braid
from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.groups.braid import BraidGroup
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.integer_mod import Mod
from sage.plot.arrow import arrow2d
from sage.plot.arrow import arrow
from sage.plot.graphics import Graphics
from sage.plot.plot3d.shapes2 import bezier3d
from sage.graphs.digraph import DiGraph
from copy import deepcopy,copy
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import IntegerRing
from sage.combinat.permutation import Permutations
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing

class Link:
    r"""
    The base class for Link, can be initialized by giving in three formats namely
    1. Briadword
    2. Oriented Gauss Code
    3. Planar Diagram Code
    Refer to oriented gauss code, planar diagram for the convention of the code.
    """
    def __init__(self, input = None, oriented_gauss_code = None, PD_code = None):
        if isinstance(input, Braid):
            self._braid = input
            self._PD_code = None
            self._oriented_gauss_code = None

        elif oriented_gauss_code != None:
            self._oriented_gauss_code = oriented_gauss_code
            self._braid = None
            self._PD_code = None

        elif PD_code != None:
            self._PD_code = PD_code
            self._oriented_gauss_code = None
            self._braid = None

        else:
            raise Exception("Invalid input")

    def braidword(self):
        r"""
        Returns the braidword of the link.

        OUTPUT:
            - Braidword representation of the link.

        EXAMPLES::
        """

        if self._braid != None:
            return list(self._braid.Tietze())

        elif self._oriented_gauss_code != None:
            if self.is_knot == True:
                return self.seifert_to_braid()

        elif self._gauss_code != None:
            return "Not implemented Error"

        elif self._dt_code != None:
            return "Not Implemented Error"

    def braid(self):
        r"""
        Returns the braid representation of the link.

        OUTPUT:
            - Braid representation of the link.

        EXAMPLES::
        """
        if self._braid != None:
            return self._braid

        elif self._oriented_gauss_code != None:
            gen = max([abs(i) for i in self.seifert_to_braid()])
            B = BraidGroup(gen + 1)
            return B(self.seifert_to_braid())

        elif self._PD_code != None:
            pd = self._PD_code
            L1 = Link(PD_code = pd)
            ogc = L1.oriented_gauss_code()
            L2 = Link(oriented_gauss_code = ogc)
            return L2.braid()

    def oriented_gauss_code(self):
        r"""
        Returns the oriented gauss code of the input. The oriented gauss
        code has two parts
        a. The gauss code
        b. The orientation of each crossing
        The following orientation was taken into consideration for consturction
        of knots:

        From the outgoing of the overcrossing if we move in the clockwise direction
        to reach the outgoing of the undercrossing then we label that crossing as '-'.

        From the outgoing of the overcrossing if we move in the anticlockwise
        direction to reach the outgoingo the undercrossing then we label that crossing
        as '+'.

        One more consideration we take in while constructing the orientation is:
        The order of the orientation is same as the ordering of the crossings in the
        gauss code.

        OUTPUT:
            - Oriented gauss code of the link.

        EXAMPLES::
        """
        if self._oriented_gauss_code != None:
            return self._oriented_gauss_code

        elif self._PD_code != None:
            pd = self._PD_code
            orient = self.orientation()
            crossing_info = {}
            for i,j in enumerate(pd):
                if orient[i] == '-':
                    crossing_info.update({(j[0],'under',i+1) : j[2]})
                    crossing_info.update({(j[3],'over',i+1) : j[1]})
                elif orient[i] == '+':
                    crossing_info.update({(j[0],'under',i+1) : j[2]})
                    crossing_info.update({(j[1],'over',i+1) : j[3]})
            edges = {}
            cross_number = {}
            for i,j in crossing_info.items():
                edges.update({ i[0] : [j]})
                if i[1] == 'over':
                    cross_number.update({ i[0] : i[2]})
                elif i[1] == 'under':
                    cross_number.update({ i[0] : -i[2]})
            edges_graph = DiGraph(edges)
            d = edges_graph.all_simple_cycles()
            code = []
            for i in d:
                l = []
                for j in i:
                    l.append(cross_number[j])
                del l[-1]
                code.append(l)
            oriented_code = [code, orient]
            return oriented_code

        elif self._braid != None:
            braid = self._braid
            pd = Link(braid).PD_code()
            L = Link(PD_code = pd)
            gc = L.oriented_gauss_code()
            return gc

    def PD_code(self):
        r"""
        Returns the Planar Diagram code of the knot. The Planar Diagram is returned
        in the following format.

        We construct the crossing by starting with the entering component of the
        undercrossing, move in the clockwise direction and then generate the list.
        Suppose if the crossing is given by [a, b, c, d], then we interpret this
        information as
        a is the entering component of the undercrossing
        b, d are the components of the overcrossing
        c is the leaving component of the undercrossing


        OUTPUT:
            - Planar Diagram representation of the link.

        EXAMPLES::
        """
        if self._PD_code != None:
            return self._PD_code

        elif self._oriented_gauss_code != None:
            oriented_gauss_code = self._oriented_gauss_code
            d_dic = {}
            if len(oriented_gauss_code[0]) > 1:
                d = [x for y in oriented_gauss_code[0] for x in y]
                for i,j in enumerate(d):
                    d_dic.update({ j : [i+1, i+2]})
                #here we collect the final component in each gauss code
                last_component = [i[len(i)-1] for i in oriented_gauss_code[0]]
                first_component = [i[0] for i in oriented_gauss_code[0]]
                #here we correct the last_component
                for i,j in zip(last_component,first_component):
                    d_dic[i][1] = d_dic[j][0]
                crossing_dic = {}
                for i in range(len(oriented_gauss_code[1])):
                    if oriented_gauss_code[1][i] == '-':
                        crossing_dic.update({ i + 1 : [d_dic[-(i+1)][0], d_dic[i+1][1], d_dic[-(i+1)][1], d_dic[i+1][0]]})
                    elif oriented_gauss_code[1][i] == '+':
                        crossing_dic.update({ i + 1 : [d_dic[-(i+1)][0], d_dic[i+1][0], d_dic[-(i+1)][1], d_dic[i+1][1]]})
            elif len(oriented_gauss_code[0]) == 1:
                for i,j in enumerate(oriented_gauss_code[0][0]):
                    d_dic.update({ j : [i+1, i+2]})
                d_dic[oriented_gauss_code[0][0][-1]][1] = 1
                crossing_dic = {}
                for i in range(len(oriented_gauss_code[1])):
                    if oriented_gauss_code[1][i] == '-':
                        crossing_dic.update({ i + 1 : [d_dic[-(i+1)][0], d_dic[i+1][1], d_dic[-(i+1)][1], d_dic[i+1][0]]})
                    elif oriented_gauss_code[1][i] == '+':
                        crossing_dic.update({ i + 1 : [d_dic[-(i+1)][0], d_dic[i+1][0], d_dic[-(i+1)][1], d_dic[i+1][1]]})
            pd = [crossing_dic[i] for i in crossing_dic.keys()]
            return pd

        elif self._braid != None:
            b = list(self._braid.Tietze())
            reg = [i for i in range(1, max([abs(i) for i in b]) + 2)]
            regcp = deepcopy(reg)
            pd = [[None for i in range(4)] for i in range(len(b))]
            for i,j in enumerate(b):
                 if cmp(j,0) == -1:
                    pd[i][0] = reg[abs(j)-1]
                    pd[i][1] = max(reg) + 1
                    pd[i][2] = max(reg) + 2
                    pd[i][3] = reg[abs(j)]
                    reg[regcp.index(abs(j))] = max(reg) + 1
                    reg[regcp.index(abs(j)) + 1] = max(reg) + 1
                 elif cmp(j,0) == 1:
                    pd[i][0] = reg[abs(j)]
                    pd[i][1] = reg[abs(j)-1]
                    pd[i][2] = max(reg) + 1
                    pd[i][3] = max(reg) + 2
                    reg[regcp.index(abs(j))] = max(reg) + 1
                    reg[regcp.index(abs(j)) + 1] = max(reg) + 1
            #correcting the last components in the generated pd code
            b_reversed = b[::-1]
            pd_reversed = pd[::-1]
            b_reversed_abs = [abs(i) for i in b_reversed]
            braid_contents = list(set(b_reversed_abs))
            for i in braid_contents:
                for j in b_reversed:
                    if abs(i) <= abs(b_reversed[0]):
                        if j == i:
                            pd_reversed[b_reversed.index(j)][2] = i
                            break
                        elif j == -i:
                            pd_reversed[b_reversed.index(j)][1] = i
                            break
                    elif abs(i) > abs(b_reversed[0]):
                        if j == i:
                            pd_reversed[b_reversed.index(j)][3] = i + 1
                            break
                        elif j == -i:
                            pd_reversed[b_reversed.index(j)][2] = i + 1
                            break
            if b_reversed[0] < 0:
                pd_reversed[0][2] = abs(b_reversed[0]) + 1
            elif b_reversed[0] > 0:
                pd_reversed[0][3] = abs(b_reversed[0]) + 1
            pd_original = pd_reversed[::-1]
            return pd_original

    def gauss_code(self):
        r"""
        Returns the gauss_code of the input. Gauss code is generated by the
        following procedure:

        a. We randomly number the crossings
        b. We select a point on the knot and start moving along the component
        c. At each crossing we take the number of the crossing, along with
           sign, which is '-' if it is a undercrossing and '+' if it is a
           overcrossing.

        OUTPUT:
            - Gauss code representation of the link.

        EXAMPLES::
        """
        if self._braid != None:
            b = self._braid
            L = Link(b).oriented_gauss_code()
            return L[0]

        elif self._PD_code != None:
            pd = self._PD_code
            L = Link(PD_code = pd).oriented_gauss_code()
            return L[0]

        elif self._oriented_gauss_code != None:
            return self._oriented_gauss_code[0]

    def dt_code(self):
        r"""
        Returns the dt_code.DT code is generated by the following way.

        We start moving along the knot, as we encounter the crossings we
        start numbering them, so every crossing has two numbers assigned to
        it once we have traced the entire knot. Now we take the even number
        associated with every knot. The following sign convention is to be
        followed:
        Take the even number with a negative sign if it is an overcrossing
        that we are encountering.

        OUTPUT:
            - DT Code representation of the knot. This is implemented only
              for knots.

        EXAMPLES::
        """
        if self._braid != None:
            b = list(self._braid.Tietze())
            N = len(b)
            label = [0 for i in range(2*N)]
            string = 1
            next_label = 1
            type1 = 0
            crossing = 0
            while(next_label <= 2*N):
                string_found = 0
                for i in range(crossing, N):
                    if(abs(b[i]) == string or abs(b[i]) == string - 1):
                        string_found = 1
                        crossing = i
                        break
                if(string_found == 0):
                    for i in range(0,crossing):
                        if(abs(b[i]) == string or abs(b[i]) == string - 1):
                            string_found = 1
                            crossing = i
                            break
                if(label[2*crossing + next_label%2] == 1):
                    raise Exception("Implemented only for knots")
                else:
                    label[2*crossing + next_label%2] =  next_label
                    next_label = next_label + 1
                if(type1 == 0):
                    if(b[crossing] < 0):
                        type1 = 1
                    else:
                        type1 = -1
                else:
                    type1 = -1 * type1
                    if((abs(b[crossing]) == string and b[crossing] * type1 > 0) or (abs(b[crossing]) != string and b[crossing] * type1 < 0)):
                        if(next_label%2 == 1):
                            label[2*crossing] = label[2*crossing] * -1
                if(abs(b[crossing]) == string):
                    string = string + 1
                else:
                    string = string - 1
                crossing = crossing + 1
            code = [0 for i in range(N)]
            for i in range(N):
                for j in range(N):
                    if label[2*j+1] == 2*i+1:
                        code[i] = label[2*j]
                        break
            return code

        elif self._oriented_gauss_code != None:
            ogc = self._oriented_gauss_code
            L = Link(oriented_gauss_code = ogc)
            b = L.braid()
            L1 = Link(b)
            return L1.dt_code()

        elif self._PD_code != None:
            pdc = self._PD_code
            L = Link(PD_code = pdc)
            b = L.braid()
            L1 = Link(b)
            return L1.dt_code()

    def _braidwordcomponents(self):
        r"""
        Returns the disjoint braid components, if any, else returns the braid itself.
        For example consider the braid [-1, 3, 1, 3] this can be viewed as a braid
        with components as [-1, 1] and [3, 3]. There is no common crossing to these
        two (in sense there is a crossing between strand 1 and 2, crossing between
        3 and 4 but no crossing between strand 2 and 3,so these can be viewed as
        independent components in the braid).

        OUTPUT:
            - List containing the components is returned.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponents()
            [[-1, 1], [3, 3]]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponents()
            [[-1, 1, 1, 1], [3], [5, 7, 6]]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponents()
            [[-2, 1, 1], [4, 4], [6]]
        """
        b = self.braid()
        ml = list(b.Tietze())
        if ml == []:
            raise Exception("The braid remains the same with no components")
        else:
            l = list(set([abs(k) for k in ml]))
            missing1 = list(set(range(min(l),max(l)+1)) - set(l))
            if len(missing1) == 0:
                return [ml]
            else:
                missing = sorted(missing1)
                x = [[] for i in range(len(missing) + 1)]
                for i in range(len(missing)):
                    for j in range(len(ml)):
                        if(ml[j] != 0 and abs(ml[j]) < missing[i]):
                            x[i].append(ml[j])
                            ml[j] = 0
                        elif(ml[j] != 0 and abs(ml[j]) > missing[-1]):
                            x[-1].append(ml[j])
                            ml[j] = 0
                y2 = [x for x in x if x != []]
                return y2

    def _braidwordcomponentsvector(self):
        r"""
        The list from the braidwordcomponents is flattened to give out the vector form.

        OUTPUT:
            - Vector containing braidwordcomponents.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponentsvector()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponentsvector()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponentsvector()
            [-2, 1, 1, 4, 4, 6]
        """
        bc = self._braidwordcomponents()
        return [x for y in bc for x in y]

    def homology_generators(self):
        r"""
        Returns the homology generators of the braidword.
        This method uses the braidwordcomponentsvector to generate the homology generators.
        The position of the repeated element w.r.t the braidwordcomponentvector list is
        compiled into a list.

        OUTPUT:
            - The homology generators relating to the braid word representation

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.homology_generators()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.homology_generators()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.homology_generators()
            [0, 2, 0, 4, 0]
        """
        x4 = self._braidwordcomponentsvector()
        hom_gen = []
        for j in range(len(x4)-1):
            a = abs(x4[j])
            for i in range(j+1, len(x4)):
                if(a == abs(x4[i])):
                    hom_gen.append(i)
                    break
            else:
                hom_gen.append(0)
        return hom_gen

    def Seifert_Matrix(self):
        r"""
        Returns the Seifert Matrix associated with the braidword.
        This is further used to calculate the Alexander Polynomial.

        OUTPUT:
            - Returns the Seifert Matrix of the link.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.Seifert_Matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.Seifert_Matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.Seifert_Matrix()
            [-1  0]
            [ 0 -1]
        """
        x5 = self._braidwordcomponentsvector()
        h = self.homology_generators()
        hl = len(h)
        A = matrix(ZZ, hl, hl)
        for i in range(hl):
            if h[i] != 0:
                for j in range(i,hl):
                        if i == j:
                            A[i,j] = -cmp((x5[i] + x5[h[i]]),0)
                        elif (h[i] > h[j]):
                            A[i,j] = 0
                            A[j,i] = 0
                        elif (h[i] <  j):
                            A[i,j] = 0
                            A[j,i] = 0
                        elif (h[i] == j):
                            if(x5[j] > 0):
                                A[i,j] = 0
                                A[j,i] = 1
                            else:
                                A[i,j] = -1
                                A[j,i] = 0
                        elif abs(abs(x5[i]) - abs(x5[j])) > 1:
                            A[i,j] =  0
                        elif (abs(x5[i]) - abs(x5[j]) == 1):
                            A[i,j] = 0
                            A[j,i] = -1
                        elif (abs(x5[j])- abs(x5[i]) == 1):
                            A[i,j] = 1
                            A[j,i] = 0
                        else: # for debugging
                            A[i,j] = 2
                            A[j,i] = 2
            else:
                for k in range(hl):
                    A[k,i] = 0
                    A[i,k] = 0
        k = []
        for i in range(hl):
                if h[i] == 0:
                    k.append(i)
        for i in reversed(k):
                A = A.delete_rows([i])
                A = A.delete_columns([i])
        return A


    def link_number(self):
        r"""
        Returns the link number

        OUTPUT:
            - Link number of the link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.link_number()
            4
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.link_number()
            5
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.link_number()
            1
        """
        p = self.braid().permutation()
        return len(p.to_cycles())

    def is_knot(self):
        r"""
        Returns true if the link is knot.
        Every knot is a link but the converse is not true.

        OUTPUT:
            - True if knot else False

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([1,3,1,-3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        if self.link_number() == 1:
            return True
        else:
            return False

    def genus(self):
        r"""
        Returns the genus of the link

        OUTPUT:
            - Genus of the Link.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.genus()
            0
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.genus()
            1
        """
        b = self.braidword()
        if b == []:
            return 0
        else:
            B = self.braid().parent()
            x = self._braidwordcomponents()
            q = []
            genus = 0
            s = [Link(B(x[i])).smallest_equivalent() for i in range(len(x))]
            t = [Link(B(s[i])).link_number() for i in range(len(s))]
            for i,j in enumerate(s):
                if j == []:
                    s[i].append(-2)
            for i in s:
                q1 = (abs(k)+1 for k in i)
                q2 = max(q1)
                q.append(q2)
            g = [((2 - t[i]) + len(x[i]) - q[i])/2 for i in range(len(x))]
            for i in range(len(g)):
                genus = genus + g[i]
            return genus

    def smallest_equivalent(self):
        r"""
        Returns the braidword

        OUTPUT:
            - Smallest equivalent of the given braid word representation.

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(5)
            sage: L = link.Link(B([-2, 4, 2, 4]))
            sage: L.smallest_equivalent()
            [-1, 3, 1, 3]
            sage: L = link.Link(B([-1, 1]))
            sage: L.smallest_equivalent()
            []
        """
        b = list(self.braid().Tietze())
        if not b:
            return list(b)
        else:
            b1 = min([abs(k) for k in b])
            for i in range(len(b)):
                if b[i] > 0:
                    b[i] = b[i] - b1 + 1
                else:
                    b[i] = b[i] + b1 - 1
            return b

    def signature(self):
        r"""
        Returns the signature of the link

        OUTPUT:
            - Signature of the Link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.signature()
            -1
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.signature()
            -2
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.signature()
            -2
        """
        m = 2*(self.Seifert_Matrix() + self.Seifert_Matrix().transpose())
        e = m.eigenvalues()
        sum = 0
        s = []
        for i,j in enumerate(e):
            s.append(cmp(j,0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var ='t'):
        r"""
        Returns the alexander polynomial of the link

        OUTPUT:
            - Alexander Polynomial of the Link

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            t^2 - 2*t + 1
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^2 - t + 1
        """
        R = PolynomialRing(ZZ, var)
        t = R.gen()
        m2 = self.Seifert_Matrix() - t* (self.Seifert_Matrix().transpose())
        return m2.determinant()

    def knot_determinant(self):
        r"""
        Returns the determinant of the knot

        OUTPUT:
            - Determinant of the Knot

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 2, 1, 2]))
            sage: L.knot_determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.knot_determinant()
            3
            sage: L = link.Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.knot_determinant()
            65
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            return abs(a(-1))
        else:
            return "is defined for a knot"

    def arf_invariant(self):
        r"""
        Returns the arf invariant only if the link is knot

        OUTPUT:
            - Arf invariant of knot

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, 2, 1, 2]))
            sage: L.arf_invariant()
            0
            sage: B = BraidGroup(8)
            sage: L = link.Link(B([-2, 3, 1, 2, 1, 4]))
            sage: L.arf_invariant()
            0
            sage: L = link.Link(B([1, 2, 1, 2]))
            sage: L.arf_invariant()
            1
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            if ((Mod(a(-1),8) == 1) or (Mod(a(-1),8) == 7)):
                return 0
            else:
                return 1
        else:
            raise Exception("Arf invariant is defined only for knots")

    def is_alternating(self):
        r"""
        Returns True if the given knot diagram is alternating else returns False.
        Alternating diagram implies every over cross is followed by an under cross
        or the vice-versa.

        We look at the gauss code if the sign is alternating, True is returned else
        the knot is not alternating False is returned.

        OUTPUT:
            - True if the knot diagram is alternating else False

        EXAMPLES::

            sage: from sage.knots import link
            sage: B = BraidGroup(4)
            sage: L = link.Link(B([-1, -1, -1, -1]))
            sage: L.is_alternating()
            False
            sage: L = link.Link(B([1, -2, -1, 2]))
            sage: L.is_alternating()
            False
            sage: L = link.Link(B([-1, 3, 1,3, 2]))
            sage: L.is_alternating()
            False
            sage: L = link.Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.is_alternating()
            False
            sage: L = link.Link(B([-1,2,-1,2]))
            sage: L.is_alternating()
            True
        """
        if self.is_knot() == True:
            x = self.gauss_code()
            s = [cmp(i,0) for i in x]
            if s == [(-1)**i+1 for i in range(len(x))] or s == [(-1)**i for i in range(len(x))]:
                return True
            else:
                return False
        else:
            return False

    def knot_diagram(self):
        x = self.PD_code()
        #p = [i for i in range(len(x))]
        p =[[None for i in range(4)] for i in range(len(x))]
        plt = Graphics()
        for i in range(len(x)):
            #print (p[i],p[i] + 1)
            a = x[i][0]
            plt = plt + arrow((i,i,i), (i + 0.4,i,i), legend_color='purple') + arrow((i+0.6,i,i),(i+1,i,i))
            p[i][0] = ((i,i,i)) #((i,i),(i + 0.4, i))
            p[i][2] = ((i+0.6,i,i)) #((i+0.6,i),(i+1,i))
            plt = plt + arrow((i+0.5,i,i-0.5),(i+0.5,i,i-0.1)) + arrow((i+0.5,i,i+0.1),(i+0.5,i,i+0.5))
            p[i][1] = (i+0.5,i,i-0.5) #((i+0.5,i-0.5),(i+0.5,i-0.1))
            p[i][3] = (i+0.5,i,i+0.1) #((i+0.5,i+0.1),(i+0.5,i+0.5))
        #print p
        #plt = plt + arrow2d((0,1),(1,2))
        #plt = plt + arrow((2,1),(3,2))
        q = [x[j][i] for j in range(len(x)) for i in range(4)]
        r = [list(p[j][i]) for j in range(len(p)) for i in range(4)]
        t = []
        print q
        print r
        for i in range(1,len(q)+1):
            for j in range(len(q)):
                if q[j] == i:
                    t.append(j)
                    #plt = plt + bezier_path([[r[j]]])
        print t
        #s = [(-1)*r[t[i]] for i in range(len(t))]
        for i in range(0,len(t),2):
            print r[t[i]], r[t[i+1]]
            path = [[tuple(r[t[i]]),tuple(r[t[i+1]])]]
            b = bezier3d(path, color='green')
            plt = plt + b
            #plt = plt + bezier_path([[(s[i]),(s[i+1])]]).plot3d()
        return plt

    #*****************************************************************************
    # Various methods for the implementation of the Vogel's algorithm
    # start here.The input is the oriented gauss code where in every cross is
    # given a sign for the orientation. The missing information in the last
    # implementation of PD was which part of the over crossing we encounter first.
    # So we took the braid word as the input, now in addition to the gauss code we
    # take the orientation at every crossing and call it the oriented gauss code.
    # Eventually after this algorithmic implementation we would like to work on
    # making the planar diagram code as a standard input.
    #
    #  The crux of the Vogel algorithm is two steps, identify the unoriented
    # Seifert circle and perform a move. The first part deals with identifying
    # the regions and Seifert circles.
    #*****************************************************************************

    #**************************** PART - 1 ***************************************
    def orientation(self):
        r"""
        Returns the orientation of the crossings from the input. We construct the entering, leaving information at
        each crossing to get to the orientation.

        OUTPUT:
            - Orientation  of the crossings.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(PD_code = [[1, 4, 5, 2], [3, 5, 6, 7], [4, 8, 9, 6], [7, 9, 10, 11], [8, 1, 13, 10], [11, 13, 2, 3]])
        sage: L.orientation()
        ['-', '+', '-', '+', '-', '+']
        sage: L = link.Link(PD_code=[[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
        sage: L.orientation()
        ['-', '-', '-', '-', '+', '-', '+']
        """
        y = self.PD_code()
        x = deepcopy(y)
        under = [ [[i[0],"entering"],[i[2],"leaving"]] for i in x]
        under = [a for b in under for a in b]
        over = [ [[i[1],None], [i[3],None]] for i in x]
        over = [a for b in over for a in b]
        for i in over:
            for j in under:
                if i[0] == j[0]:
                    if j[1] == "entering":
                        i[1] = "leaving"
                    elif j[1] == "leaving":
                        i[1] = "entering"
        for i in over:
            if i[1] == None:
                over = rule_1(over)
                over = rule_2(over)
        unfilled = []
        for i in over:
            if i[1] == None:
                unfilled.append(i)
        if len(unfilled) != 0:
            over[over.index(unfilled[0])][1] = "entering"
            for i in over:
                if i[1] == None:
                    over = rule_1(over)
                    over = rule_2(over)
        orientation = []
        for i in range(0, len(over), 2):
            if over[i][1] == "leaving":
                orientation.append('-')
            elif over[i][1] == "entering":
                orientation.append('+')
        return orientation

    def seifert_circles(self):
        r"""
        Returns the seifert circles from the input. Seifert circles are obtained
        by smoothing the crossings in the following way:
        All the crossings are assigned four numbers which form the components of the
        crossings.
        At a crossing the under cross entering component would go to the over cross
        leaving component and the over cross entering component would go to the
        under cross leaving component.
        We start with a component of the crossing and start smoothing as according to
        the above rule until we return to the starting component.

        OUTPUT:
            - Seifert circles of the given link.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[[1, -2, 3, -4, 2, -1, 4, -3]],['+','+','-','-']])
        sage: L.seifert_circles()
        [[6, 2], [8, 4], [7, 5, 3, 1]]
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7]],['-','-','-','-','+','+','-','+']])
        sage: L.seifert_circles()
        [[10, 6, 12, 2], [16, 8, 14, 4], [13, 9, 3, 15, 5, 11, 7, 1]]
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],['-','-','-','-','+','-','+']])
        sage: L.seifert_circles()
        [[13, 9], [12, 10, 4], [8, 14, 6, 2], [7, 3, 11, 5, 1]]
        sage: L = link.Link(PD_code=[[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
        sage: L.seifert_circles()
        [[13, 9], [12, 10, 4], [8, 14, 6, 2], [7, 3, 11, 5, 1]]
        sage: L = link.Link(PD_code = [[1,4,2,3],[4,1,3,2]])
        sage: L.seifert_circles()
        [[4, 1], [3, 2]]
        sage: L = link.Link(PD_code=[[1,11,2,10],[6,2,7,3],[3,12,4,9],[9,5,10,6],[8,1,5,4],[11,8,12,7]])
        sage: L.seifert_circles()
        [[11, 8, 1], [9, 6, 3], [7, 12, 4, 5, 10, 2]]
        """
        pd = self.PD_code()
        orient = self.orientation()
        seifert_pairs = []
        for i,j in enumerate(pd):
            x = [[],[]]
            if orient[i] == '-':
                x[0].append(j[0])
                x[0].append(j[1])
                x[1].append(j[3])
                x[1].append(j[2])
            elif orient[i] == '+':
                x[0].append(j[0])
                x[0].append(j[3])
                x[1].append(j[1])
                x[1].append(j[2])
            seifert_pairs.append(x)
        flatten = [x for y in seifert_pairs for x in y]
        flatten = [x for y in flatten for x in y]
        dic = {}
        dic = {}
        for i in range(0, len(flatten), 2):
            dic.update({ flatten[i] :[flatten[i+1]]})
        D = DiGraph(dic)
        d = D.all_simple_cycles()
        for i in d:
            del i[0]
        return d


    def regions(self):
        r"""
        Returns the regions from the input.Regions are obtained always turning left at the crossing.
        The smoothing used for obtaining the Seifert circles case was specific but here whenever a
        crossing in encountered a left turn is taken until we reach the starting component.
        The following procedure is followed for the development of regions:
        Start at a component of the crossing and keep turning left as the crossings are encountered
        until the starting component is reached.
        Append a negative sign to the component if it were encountered in the opposite direction to
        the direction it is actually in.

        OUTPUT:
            - Regions of the knot.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],['-','-','-','-','+','-','+']])
        sage: L.regions()
        [[4, -11], [2, -7], [6, -1], [13, 9], [-4, -10, -12], [-8, -2, -6, -14], [10, -3, 8, -13], [14, -5, 12, -9], [7, 3, 11, 5, 1]]
        sage: L = link.Link(oriented_gauss_code = [[[1, -2, 3, -4, 2, -1, 4, -3]],['+','+','-','-']])
        sage: L.regions()
        [[-2, -6], [8, 4], [5, 3, -8], [2, -5, -7], [1, 7, -4], [6, -1, -3]]
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7]],['-','-','-','-','+','+','-','+']])
        sage: L.regions()
        [[6, -11], [15, -4], [9, 3, -14], [2, -9, -13], [1, 13, -8], [12, -1, -7], [5, 11, 7, -16], [-3, 10, -5, -15], [-6, -10, -2, -12], [16, 8, 14, 4]]
        """
        pd = self.PD_code()
        orient = self.orientation()
        regions = []
        for i,j in enumerate(pd):
            x = [[],[]]
            if orient[i] == '-':
                x[0].append(j[0])
                x[0].append(j[1])
                x[1].append(j[3])
                x[1].append(-j[0])
            elif orient[i] == '+':
                x[0].append(j[0])
                x[0].append(-j[1])
                x[1].append(j[1])
                x[1].append(j[2])
            regions.append(x)
        pd_opposite = deepcopy(pd)
        for i in range(len(pd_opposite)):
            pd_opposite[i][0] = -pd[i][2]
            pd_opposite[i][1] = -pd[i][3]
            pd_opposite[i][2] = -pd[i][0]
            pd_opposite[i][3] = -pd[i][1]
        regions_op = []
        for i,j in enumerate(pd_opposite):
            x_op = [[],[]]
            if orient[i] == '-':
                x_op[0].append(j[0])
                x_op[0].append(j[1])
                x_op[1].append(j[3])
                x_op[1].append(-j[0])
            elif orient[i] == '+':
                x_op[0].append(j[0])
                x_op[0].append(-j[1])
                x_op[1].append(j[1])
                x_op[1].append(j[2])
            regions_op.append(x_op)
        regions_final = []
        for i in range(len(regions)):
            for j in range(len(regions[i])):
                regions_final.append(regions[i][j])
                regions_final.append(regions_op[i][j])
        dic = {}
        for i in regions_final:
            dic.update({i[0] : [i[1]]})
        D = DiGraph(dic)
        d = D.all_simple_cycles()
        for i in d:
            del i[0]
        return d

    def vogel_move(self):
        r"""
        Returns the Planar Diagram code if there is a vogel's move required, else returns
        "no move required". Whether a move is required or not is decided by the following
        criteria:
        A bad region is one that has two components with the same sign, but that belong
        to different Seifert circles.
        We detect the presence of a bad region and then perform the move. The move is done as follows:
        Perform a Reidemeister move of type 2 by pulling the component with a higher number
        over the component having the lower number. As this is done we have 2 more crossings added to
        the initial system. We renumber everything and return the Planar Diagram Code of the modified
        diagram

        OUTPUT:
            - Planar diagram after the move is performed.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[[1, -2, 3, -4, 2, -1, 4, -3]],['+','+','-','-']])
        sage: L.vogel_move()
        'No move required'
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7]],['-','-','-','-','+','+','-','+']])
        sage: L.vogel_move()
        'No move required'
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],['-','-','-','-','+','-','+']])
        sage: L.vogel_move()
        [[1, 9, 2, 8], [9, 3, 10, 2], [5, 13, 6, 12], [13, 7, 14, 6], [18, 7, 1, 8], [17, 11, 18, 10], [14, 11, 15, 12], [16, 4, 17, 3], [15, 4, 16, 5]]
        sage: L = link.Link(PD_code=[[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
        sage: L.vogel_move()
        [[1, 9, 2, 8], [9, 3, 10, 2], [5, 13, 6, 12], [13, 7, 14, 6], [18, 7, 1, 8], [17, 11, 18, 10], [14, 11, 15, 12], [16, 4, 17, 3], [15, 4, 16, 5]]
        """
        pd = self.PD_code()
        pd_copy = deepcopy(pd)
        sc = self.seifert_circles()
        regions = self.regions()
        regions_copy = deepcopy(regions)
        orient = self.orientation()
        #separating the components into positive and negative
        q = [[[],[]] for i in range(len(regions))]
        for i,j in enumerate(regions):
            for k in j:
                if k < 0:
                    q[i][0].append(k)
                elif k > 0:
                    q[i][1].append(k)
        #making all the components positive
        r = [[] for i in range(len(q))]
        for i in range(len(q)):
            for j in range(len(q[i])):
                r[i].append(None)
                for k in range(len(q[i][j])):
                    if q[i][j][k] < 0:
                        q[i][j][k] = (-1)*q[i][j][k]
        #to find the intersection of regions with seifert circles
        #first clean q that is by removing empty and single length arrays
        clean_q = []
        for i in q:
            for j in i:
                if len(j) >= 2:
                    clean_q.append(j)
        #detecting the bad region
        bad_region = []
        for i in clean_q:
            for j in sc:
                if set(i).intersection(set(j)) == set(i):
                    break
            else:
                bad_region.append(i)
        if bad_region == []:
            return "No Vogel Move"
        else:
            #here there might be many bad regions but we only select one
            #here it is the first one
            bad_region = bad_region[0]
            #finding the max of the pd to develop the new crossings
            pd_max = max([max(i) for i in pd_copy])
            #editing the previous crossings
            #the maximum is corrected
            for i in pd_copy:
                if max(bad_region) == i[0]:
                    i[0] = pd_max + 4
            #editing the contents of the pd code with the minimum
            #the minimum is corrected
            min_cross = [[j,orient[i]] for i,j in enumerate(pd_copy) if min(bad_region) in j]
            y = [[] for i in range(len(min_cross))]
            for i,j in enumerate(min_cross):
                if j[1] == '-':
                    y[i].append((j[0][0],"entering"))
                    y[i].append((j[0][1],"leaving"))
                    y[i].append((j[0][2],"leaving"))
                    y[i].append((j[0][3],"entering"))
                if j[1] == '+':
                    y[i].append((j[0][0],"entering"))
                    y[i].append((j[0][1],"entering"))
                    y[i].append((j[0][2],"leaving"))
                    y[i].append((j[0][3],"leaving"))
            for i in y:
                if (min(bad_region),"entering") in i:
                    pd_copy[pd_copy.index(min_cross[y.index(i)][0])][(pd_copy[pd_copy.index(min_cross[y.index(i)][0])]).index(min(bad_region))] = pd_max + 2
            #sorting the regions in positive and negative ones.
            pos_neg = [[[],[]] for i in range(len(regions_copy))]
            for i,j in enumerate(regions_copy):
                for k in j:
                    if k > 0:
                        pos_neg[i][0].append(k)
                    elif k < 0:
                        pos_neg[i][1].append(k)
            pos = [i[0] for i in pos_neg if i[0] != []]
            neg = [[-j for j in i[1]] for i in pos_neg if i[1] != []]
            #creating the new components, the lesser in the bad region is
            #always over and the greater is always under
            new_component_1 = [min(bad_region), pd_max + 1, pd_max + 2]
            new_component_2 = [max(bad_region), pd_max + 3, pd_max + 4]
            #creating new crossings
            crossing_1 = [None for i in range(4)]
            crossing_2 = [None for i in range(4)]
            if bad_region in pos:
                crossing_1[0] = max(bad_region)
                crossing_1[1] = pd_max + 2
                crossing_1[2] = pd_max + 3
                crossing_1[3] = pd_max + 1
                crossing_2[0] = pd_max + 3
                crossing_2[1] = min(bad_region)
                crossing_2[2] = pd_max + 4
                crossing_2[3] = pd_max + 1
                pd_copy.append(crossing_1)
                pd_copy.append(crossing_2)
                return pd_copy
            elif bad_region in neg:
                crossing_1[0] = max(bad_region)
                crossing_1[1] = pd_max + 1
                crossing_1[2] = pd_max + 3
                crossing_1[3] = pd_max + 2
                crossing_2[0] = pd_max + 3
                crossing_2[1] = pd_max + 1
                crossing_2[2] = pd_max + 4
                crossing_2[3] = min(bad_region)
                pd_copy.append(crossing_1)
                pd_copy.append(crossing_2)
                return pd_copy

    def final(self):
        r"""
        Returns the Planar Diagram code, Seifert circles and regions after all the bad regions
        have been removed by performing the vogel's move.

        OUTPUT:
            - Planar Diagram code, Seifert circle, Regions after the bad regions have been
              removed

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[[1, -2, 3, -4, 2, -1, 4, -3]],['+','+','-','-']])
        sage: L.final()
        [[[6, 2], [8, 4], [7, 5, 3, 1]],
        [[-2, -6], [8, 4], [5, 3, -8], [2, -5, -7], [1, 7, -4], [6, -1, -3]],
        [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]]]
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7]],['-','-','-','-','+','+','-','+']])
        sage: L.final()
        [[[10, 6, 12, 2], [16, 8, 14, 4], [13, 9, 3, 15, 5, 11, 7, 1]],
        [[6, -11], [15, -4], [9, 3, -14], [2, -9, -13], [1, 13, -8], [12, -1, -7], [5, 11, 7, -16], [-3, 10, -5, -15], [-6, -10, -2, -12], [16, 8, 14, 4]],
        [[1, 13, 2, 12], [9, 3, 10, 2], [14, 4, 15, 3], [4, 16, 5, 15], [10, 5, 11, 6], [6, 11, 7, 12], [16, 8, 1, 7], [13, 8, 14, 9]]]
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],['-','-','-','-','+','-','+']])
        sage: L.final()
        [[[18, 4], [21, 15], [9, 3, 19, 11, 17, 5, 13, 7, 1], [10, 20, 16, 12, 6, 14, 22, 8, 2]],
        [[-15, -21], [6, -13], [2, -9], [8, -1], [18, 4], [-3, 10, -19], [12, -5, -17], [20, 16, -11], [14, 22, -7], [19, 11, 17, -4], [21, -14, -6, -12, -16], [15, -20, -10, -2, -8, -22], [5, 13, 7, 1, 9, 3, -18]],
        [[1, 9, 2, 8], [9, 3, 10, 2], [5, 13, 6, 12], [13, 7, 14, 6], [22, 7, 1, 8], [19, 11, 20, 10], [16, 11, 17, 12], [18, 4, 19, 3], [17, 4, 18, 5], [21, 14, 22, 15], [20, 16, 21, 15]]]
        """
        x = self.PD_code()
        while True:
            link = Link(PD_code = x)
            PD_code_old =  x
            x = link.vogel_move()
            if x == "No Vogel Move":
                x = PD_code_old
                break
        L = Link(PD_code = x)
        sc = L.seifert_circles()
        regions = L.regions()
        orientation = L.orientation()
        pd_code = x
        final = [sc, regions, pd_code, orientation]
        return final

    #**************************** PART - 2 ***************************************
    def seifert_to_braid(self):
        r"""
        Returns the braidword of the input. We match the outgoing components to the
        incoming components, in doing so we order the crossings and see to which strand
        they belong, thereby developing the braidword

        OUTPUT:
            - Braidword representation of the link.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[1, -2, 3, -4, 2, -1, 4, -3],['+','+','-','-']])
        sage: L.seifert_to_braid()
        [1, -2, 1, -2]
        sage: L = link.Link(oriented_gauss_code = [[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5],['-','-','-','-','+','-','+']])
        sage: L.seifert_to_braid()
        [-1, -2, -3, 2, 1, -2, -2, 3, 2, -2, -2]
        sage: L = link.Link(oriented_gauss_code = [[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7],['-','-','-','-','+','+','-','+']])
        sage: L.seifert_to_braid()
        [-1, 2, -1, -2, -2, 1, 1, -2]
        """
        #all the data from the previous method
        sc = self.final()[0]
        regions = self.final()[1]
        pd_code = self.final()[2]
        orient = self.final()[3]
        #making the regions positive
        regions_pos = deepcopy(regions)
        for i in range(len(regions_pos)):
            for j in range(len(regions_pos[i])):
                if regions_pos[i][j] < 0:
                    regions_pos[i][j] = (-1)*regions_pos[i][j]
        #finding which sc are same as regions
        #r[0] is the first seifert cirlce and r[1] is the last one
        #which coincides with a region and there are exactly two seifert
        #circles here.
        r = []
        for i in sc:
            for j in regions_pos:
                if set(i) == set(j):
                    r.append(i)
        #here t stores the crossing information required to check which one
        #is the next
        t = [[None for i in range(2*len(sc[j]))] for j in range(len(sc))]
        for i in range(len(sc)):
            for j in range(len(sc[i])):
                t[i][2*j] = sc[i][j] - 1
                t[i][2*j+1] = sc[i][j] + 1
        #here we find the order of the seifert circles
        seifert_order = []
        t1 = deepcopy(t)
        a = r[0]
        b = r[1]
        del t[0]
        i = 0
        while True:
            for i in range(len(t)):
                if len(list(set(a).intersection(set(t[i])))) >= 2:
                    r = t1.index(t[i])
                    seifert_order.append(a)
                    del t[i]
                    a = sc[r]
                    break
            else:
                seifert_order.append(b)
                break
        entering = []
        leaving = []
        for i in range(len(pd_code)):
            t = []
            q = []
            if orient[i] == '-':
                t.append(pd_code[i][0])
                t.append(pd_code[i][3])
                q.append(pd_code[i][1])
                q.append(pd_code[i][2])
            elif orient[i] == '+':
                t.append(pd_code[i][0])
                t.append(pd_code[i][1])
                q.append(pd_code[i][2])
                q.append(pd_code[i][3])
            entering.append(t)
            leaving.append(q)
        #here we correct the leaving and entering components so they belong to the
        #correct seifert circles
        for val in zip(seifert_order, seifert_order[1:]):
            for i in entering:
                if i[1] in val[0] and i[0] in val[1]:
                    a = i[0]
                    i[0] = i[1]
                    i[1] = a
        for val in zip(seifert_order, seifert_order[1:]):
            for i in leaving:
                if i[1] in val[0] and i[0] in val[1]:
                    a = i[0]
                    i[0] = i[1]
                    i[1] = a
        first_seifert = sc[0]
        first_crossing = None
        for i in pd_code:
            if len(list(set(i).intersection(set(first_seifert)))) == 2:
                first_crossing = i
                break
        tmp = []
        if orient[pd_code.index(first_crossing)] == '-':
            tmp.append(first_crossing[1])
            tmp.append(first_crossing[2])
        elif orient[pd_code.index(first_crossing)] == '+':
            tmp.append(first_crossing[2])
            tmp.append(first_crossing[3])
        reg = [0 for i in range(len(sc))]
        if tmp[0] in first_seifert:
            reg[0] = tmp[0]
            reg[1] = tmp[1]
        elif tmp[1] in first_seifert:
            reg[0] = tmp[1]
            reg[1] = tmp[0]
        crossing = []
        crossing.append(first_crossing)
        q = 0
        while q < len(pd_code):
            for val in zip(reg, reg[1:]):
                for i,j in enumerate(entering):
                    if list(val) == j:
                        crossing.append(pd_code[i])
                        reg[reg.index(val[0])] = leaving[i][0]
                        reg[reg.index(val[1])] = leaving[i][1]
                        q = q + 1
                        break
                    if val[0] == j[0] and val[1] == 0:
                        crossing.append(pd_code[i])
                        reg[reg.index(val[0])] = leaving[i][0]
                        reg[reg.index(val[1])] = leaving[i][1]
                        q = q + 1
                        break
        del crossing[-1]
        #before we read the braid words, we make '-' in orient to -1 and similarly
        # '+' to +1
        sign = [None for i in range(len(orient))]
        for i,j in enumerate(orient):
            if j == '-':
                sign[i] = -1
            elif j == '+':
                sign[i] = 1
        #each crossing belongs to two seifert circles, we find the first and break
        braid = []
        for i in crossing:
            for j in seifert_order:
                if len(list(set(i).intersection(set(j)))) == 2:
                    braid.append(sign[pd_code.index(i)]*(seifert_order.index(j) + 1))
                    break
        return braid

    def writhe(self):
        r"""
        Returns the writhe of the knot.

        OUTPUT:
            - Writhe of the knot.

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[1, -2, 3, -4, 2, -1, 4, -3],['+','+','-','-']])
        sage: L.writhe()
        0
        sage: L = link.Link(oriented_gauss_code = [[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5],['-','-','-','-','+','-','+']])
        sage: L.writhe()
        -3
        sage: L = link.Link(oriented_gauss_code = [[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, -4, -7],['-','-','-','-','+','+','-','+']])
        sage: L.writhe()
        -2
        """
        x = self.oriented_gauss_code()
        pos = x[1].count('+')
        neg = (-1)*x[1].count('-')
        return pos + neg

    def jones_polynomial(self, var = 't'):
        r"""
        Returns the jones polynomial of the input. The following procedure is used to
        determine the jones polynomial. We track the trip matrix of the knot and keep
        changing through the various combinations of smoothing the knot.

        OUTPUT:
            - Jones Polynomial of the link.

        Reference:
            - http://katlas.math.toronto.edu/wiki/The_Jones_Polynomial#How_is_the_Jones_polynomial_computed.3F

        EXAMPLES::
        sage: from sage.knots import link
        sage: L = link.Link(oriented_gauss_code = [[[1, -2, 3, -4, 2, -1, 4, -3]],['+','+','-','-']])
        sage: L.jones_polynomial()
        (-t^10 - t^6 - t^2 + t^-6)/(-t^2 - t^-2)
        sage: L = link.Link(oriented_gauss_code = [[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],['-','-','-','-','+','-','+']])
        sage: L.jones_polynomial()
        (-t^6 - 2*t^-2 + t^-18 - t^-22 + t^-26)/(-t^2 - t^-2)
        """
        pd = self.PD_code()
        orient = self.orientation()
        #here we look at the smoothings, either 0 or 1
        label_1 = [0 for i in range(len(pd))]
        label_2 = [1 for i in range(len(pd))]
        label_1.extend(label_2)
        #all the crossings have to be either smoothened by 0 or 1
        P = Permutations(label_1, len(pd))
        R = LaurentPolynomialRing(ZZ, var)
        x = R.gen()
        crossing_to_label = []
        #we record how each crossing is smoothened and also the sign of the crossing
        for i in P.list():
            tmp = {}
            for j,k in enumerate(i):
                tmp.update({tuple(pd[j]) : (k,orient[j])})
            crossing_to_label.append(tmp)
        #we calculate the coefficients
        coeff = []
        for i in crossing_to_label:
            tmp_coeff = []
            for j in pd:
                if i[tuple(j)][0] == 0:
                    tmp_coeff.append(x)
                    #tmp_coeff.append(x**(1/4))
                elif i[tuple(j)][0] == 1:
                    tmp_coeff.append(x**-1)
                    #tmp_coeff.append(x**(-1/4))
            coeff.append(tmp_coeff)
        #the product of the coefficients is calcaluated
        product = []
        for i in coeff:
            p = 1
            for j in i:
                p = p*j
            product.append(p)
        #here we calculate the number of circles after the smoothing has been performed at
        #every crossing.
        pd_edit = []
        for i in crossing_to_label:
            pd_copy = deepcopy(pd)
            tmp = []
            for j in pd_copy:
                if i[tuple(j)][0] == 0:
                    if i[tuple(j)][1] == '-':
                        tmp.append([pd_copy[pd_copy.index(j)][0],pd_copy[pd_copy.index(j)][1]])
                        tmp.append([pd_copy[pd_copy.index(j)][3],pd_copy[pd_copy.index(j)][2]])
                    elif i[tuple(j)][1] == '+':
                        tmp.append([pd_copy[pd_copy.index(j)][1],pd_copy[pd_copy.index(j)][2]])
                        tmp.append([pd_copy[pd_copy.index(j)][0],pd_copy[pd_copy.index(j)][3]])
                elif i[tuple(j)][0] == 1:
                    if i[tuple(j)][1] == '-':
                        tmp.append([pd_copy[pd_copy.index(j)][1],pd_copy[pd_copy.index(j)][2]])
                        tmp.append([pd_copy[pd_copy.index(j)][0],pd_copy[pd_copy.index(j)][3]])
                    elif i[tuple(j)][1] == '+':
                        tmp.append([pd_copy[pd_copy.index(j)][2],pd_copy[pd_copy.index(j)][3]])
                        tmp.append([pd_copy[pd_copy.index(j)][1],pd_copy[pd_copy.index(j)][0]])
            pd_edit.append(tmp)
        #replacing the old edges with the new ones to get the circles and also the number of circles.
        pd_edit = rule_3(pd_edit)
        pd_edit = rule_3(pd_edit)
        for i in pd_edit:
            for j in range(len(i)):
                for k in reversed(range(j+1,len(i))):
                    if i[j] == i[k]:
                        del i[k]
        #the number of circles.
        circle_count = []
        for i in pd_edit:
            circle_count.append(len(i))
        #we calculate the terms of the polynomial
        terms = []
        for i,j in zip(product, circle_count):
            terms.append(i*(((-x**2-(x**(-1))**2))**j))
        #add the terms to generate the polynomial
        poly = 0
        for i in terms:
            poly = i + poly
        wri = self.writhe()
        return ((-x**(3))**wri)*poly / (-x**2-((x**(-1))**2))

def rule_1(over):
    for i in range(0,len(over),2):
        if over[i][1] == None:
            if over[i+1][1] == "entering":
                over[i][1] = "leaving"
            elif over[i+1][1] == "leaving":
                over[i][1] = "entering"
        elif over[i+1][1] == None:
            if over[i][1] == "entering":
                over[i+1][1] = "leaving"
            elif over[i][1] == "leaving":
                over[i+1][1] = "entering"
    return over

def rule_2(over):
    for i in over:
        for j in over:
            if i[0] == j[0] and j[1] == None:
                if i[1] == "entering":
                    j[1] = "leaving"
                    break
                elif i[1] == "leaving":
                    j[1] = "entering"
                    break
    return over

def rule_3(pd):
    for i in pd:
        for j,k in enumerate(i):
            a = min(k)
            b = max(k)
            for l in i[0:]:
                if b in l:
                    i[i.index(l)][l.index(b)] = min(k)
    return pd
