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

class Link:
    r"""
    The base class for Link, taking input in three formats namely Briadword, gauss_code, dt_code
    """
    #here the dt_code and gauss_code should be removed as we do not consider them as standard input.
    def __init__(self, input = None, gauss_code = None, dt_code = None, oriented_gauss_code = None, PD_code = None):
        if type(input) == Braid:
            self._braid = input
            self._gauss_code = None
            self._dt_code = None

        elif gauss_code != None:
            self._braid = None
            self._gauss_code = gauss_code
            self._dt_code = None

        elif dt_code != None:
            self._braid = None
            self._gauss_code = None
            self._dt_code = dt_code

        elif oriented_gauss_code != None:
            self._oriented_gauss_code = oriented_gauss_code
            self._PD_code = None

        elif PD_code != None:
            self._PD_code = PD_code
            self._oriented_gauss_code = None

        else:
            raise Exception("Invalid input")


    def braidword(self):
        r"""
        Returns the braidword

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Braidword representation of the INPUT

        EXAMPLES::
        """

        if self._braid != None:
            return list(self._braid.Tietze())

        elif self._gauss_code != None:
            return "Not implemented Error"

        elif self._dt_code != None:
            return "Not Implemented Error"

    def braid(self):
        r"""
        Returns the braid

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Braid representation of the INPUT

        EXAMPLES::
        """

        if self._braid != None:
            return self._braid

        elif self._gauss_code != None:
            return "Not implemented Error"

        elif self._dt_code != None:
            return "Not Implemented Error"

    def PD_code(self):
        if self._PD_code != None:
            return self._PD_code

        #Whatever is the crossing one should give the sign for that first
        #the order of the sign is as per the ordering of the crossings as in the gauss code
        #so for example if the gauss code is 1 -3 2 -1 3 -2 then the
        #order of the sign of the crossings is sign of crossing 1 then 3 then at 2
        #and so on and so forth.
        elif self._oriented_gauss_code != None:
            gc = self._oriented_gauss_code[0]
            gc_sign = self._oriented_gauss_code[1]
            l = [0 for i in range(len(gc))]
            for i in range(len(gc)):
                k = abs(gc[i])
                if l[2*(k-1)] == 0:
                    l[2*(k-1)] = (i + 1)*(cmp(gc[i],0))
                else:
                    l[2*k-1] = (i + 1)*(cmp(gc[i],0))
            y = [None for i in range(2*len(gc_sign))]
            for i in range(0,2*len(gc_sign),2):
                if l[i] < 0:
                    y[i] = 'under'
                    y[i+1] = 'over'
                elif l[i] > 0:
                    y[i] = 'over'
                    y[i+1] = 'under'
            l1 = [abs(x) for x in l]
            r = []
            p = [[None for i in range(4)] for i in range(len(gc_sign))]
            for i in range(len(gc_sign)):
                if gc_sign[i] == '-':
                    if y[2*i] == 'under':
                        p[i][0] = l1[2*i]
                        p[i][2] = p[i][0] + 1
                        p[i][3] = l1[2*i+1]
                        p[i][1] = p[i][3] + 1
                    elif y[2*i+1] == 'under':
                        p[i][0] = l1[2*i+1]
                        p[i][3] = l1[2*i]
                        p[i][2] = p[i][0] + 1
                        p[i][1] = p[i][3] + 1
                elif gc_sign[i] == '+':
                    if y[2*i] == 'under':
                        p[i][0] = l1[2*i]
                        p[i][2] = p[i][0] + 1
                        p[i][1] = l1[2*i+1]
                        p[i][3] = p[i][1] + 1
                    elif y[2*i+1] == 'under':
                        p[i][0] = l1[2*i+1]
                        p[i][1] = l1[2*i]
                        p[i][2] = p[i][0] + 1
                        p[i][3] = p[i][1] + 1
            for i in range(len(p)):
                for j in range(4):
                    if p[i][j] == max(l1) + 1:
                        p[i][j] = 1
            return p

    def gauss_code(self):
        r"""
        Returns the gauss_code

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Gauss code representation of the INPUT

        EXAMPLES::
        """
        if self._gauss_code != None:
            return self._gauss_code

        elif self._braid != None:
            L = self._braid
            B = L.parent()
            l = self.braidword()
            L = Link(B(l)).dt_code()
            gc = Link(dt_code = L).gauss_code()
            self._gauss_code = gc
            return gc

        elif self._dt_code != None:
            dt = self._dt_code
            gauss = []
            y = [None for i in range(2*len(dt))]
            x = [0 for i in range(2*len(dt))]
            for i in range(len(dt)):
                x[2*i] = 2*i + 1
                x[2*i + 1] = dt[i]
            for i in range(len(dt)):
                if x[2*i+1] > 0:
                    y[2*i+1] = 'under'
                    y[2*i] = 'over'
                elif x[2*i+1] < 0:
                    y[2*i+1] = 'over'
                    y[2*i] = 'under'
            for i in range(1,len(x)+1):
                for j in range(0,len(x)):
                    if abs(x[j]) == i:
                        if y[j] == 'under':
                            gauss.append(-(j//2 + 1))
                        elif y[j] == 'over':
                            gauss.append(j//2 + 1)
            self._gauss_code = gauss
            return gauss

    def dt_code(self):
        r"""
        Returns the dt_code

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - DT Code representation of the INPUT

        EXAMPLES::
        """
        if self._dt_code != None:
            return self._dt_code

        elif self._braid != None:
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

        elif self._gauss_code != None:
            gc = self._gauss_code
            l = [0 for i in range(len(gc))]
            for i in range(len(gc)):
                k = abs(gc[i])
                if l[2*(k-1)] == 0:
                    l[2*(k-1)] = (i + 1)*(cmp(gc[i],0))
                else:
                    l[2*k-1] = (i + 1)*(cmp(gc[i],0))
            y = [l[i] for i in range(len(l)) if abs(l[i])%2 == 0]
            x = [(-1)*y[i] for i in range(len(y))]
            return x

    def _braidwordcomponents(self):
        r"""
        Returns the braid components in an array

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Array containing the components is returned

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
        From braidwordcomponents it is converted to a vector

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Vector containing non-zero values

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
        y2 = self._braidwordcomponents()
        if len(y2) == 1:
            return y2[0]
        else:
            y3 = []
            for i in range(len(y2)):
                y = y2[i]
                for j in range(len(y)):
                    y3.append(y[j])
            return y3

    def homology_generators(self):
        r"""
        Returns the homology generators

        INPUT:
            - Either a braidword, gauss_code, dt_code

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
        Returns the Seifert Matrix

        INPUT:
            - Either a braidword, gauss_code, dt_code

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

        INPUT:
            - Either a braidword, gauss_code, dt_code

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
        Returns true if the link is knot

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - True or False

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

        INPUT:
            - Either a braidword, gauss_code, dt_code

        OUTPUT:
            - Genus of the Link

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
            for i in range(len(s)):
                if s[i] == []:
                    s[i].append(-2)
            for i in range(len(s)):
                q1 = (abs(k)+1 for k in s[i])
                q2 = max(q1)
                q.append(q2)
            g = [((2 - t[i]) + len(x[i]) - q[i])/2 for i in range(len(x))]
            for i in range(len(g)):
                genus = genus + g[i]
            return genus

    def smallest_equivalent(self):
        r"""
        Returns the braidword

        INPUT:
            - Either a braidword, gauss_code, dt_code

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

        INPUT:
            - Either a braidword, gauss_code, dt_code

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
        for i in range(len(e)):
            s.append(cmp(e[i],0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var ='t'):
        r"""
        Returns the alexander polynomial of the link

        INPUT:
            - Either a braidword, gauss_code, dt_code

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

        INPUT:
            - Either a braidword, gauss_code, dt_code

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

        INPUT:
            - Either a braidword, gauss_code, dt_code

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

    #braid to planar
    def pd_code(self):
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
                        type1 = -1
                    else:
                        type1 = 1
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
        s = [None for i in range(len(label))]
        for i in range(N):
            if cmp(label[2*i],0) == -1:
                s[2*i] = 'over'
                s[2*i+1] = 'under'
            if (label[2*i]%2 == 0 and cmp(label[2*i],0) == 1):
                s[2*i] = 'under'
                s[2*i+1] = 'over'
            if (label[2*i+1]%2 == 0 and cmp(label[2*i+1],0) == 1):
                s[2*i+1] = 'under'
                s[2*i] = 'over'
        pd = [[None for i in range(4)] for i in range(N)]
        for j in range(N):
            if s[2*j] == 'under':
                if cmp(b[j],0) == -1:
                    pd[j][0] = abs(label[2*j])
                    pd[j][3] = abs(label[2*j+1])
                    pd[j][2] = pd[j][0] + 1
                    pd[j][1] = pd[j][3] + 1
                elif cmp(b[j],0) == 1:
                    pd[j][0] = abs(label[2*j])
                    pd[j][1] = abs(label[2*j+1])
                    pd[j][2] = pd[j][0] + 1
                    pd[j][3] = pd[j][1] + 1
            elif s[2*j] == 'over':
                if cmp(b[j],0) == -1:
                    pd[j][0] = abs(label[2*j+1])
                    pd[j][2] = pd[j][0] + 1
                    pd[j][3] = abs(label[2*j])
                    pd[j][1] = pd[j][3] + 1
                if cmp(b[j],0) == 1:
                    pd[j][0] = abs(label[2*j+1])
                    pd[j][2] = pd[j][0] + 1
                    pd[j][1] = abs(label[2*j])
                    pd[j][3] = pd[j][1] + 1
        return pd

    def is_alternating(self):
        if self.is_knot() == True:
            x = self.gauss_code()
            s = [cmp(x[i],0) for i in range(len(x))]
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

    def seifert_circles(self):
        y = self.PD_code()
        x = deepcopy(y)
        s = [[None for i in range(4)] for i in range(len(x))]
        for i in range(len(x)):
            s[i][0] = 'entering'
            s[i][2] = 'leaving'
            if x[i][1] < x[i][3]:
                s[i][1] = 'entering'
                s[i][3] = 'leaving'
            if x[i][1] > x[i][3]:
                s[i][1] = 'leaving'
                s[i][3] = 'entering'
        x1 = x #the arrays are not being copied, just a tag
        s1 = s
        r = [[[None,None],[None,None]] for i in range(len(x))]
        for i in range(len(x)):
            for j in [1,3]:
                if s[i][j] == 'leaving':
                    r[i][0][0] = x1[i][0]
                    r[i][0][1] = x1[i][j]
                    del x1[i][j]
                    del x1[i][0]
                    del s1[i][j]
                    del s1[i][0]
                    break
        for i in range(len(x)):
            for j in [0,1]:
                if s1[i][0] == "entering":
                    r[i][1][0] = x1[i][0]
                    r[i][1][1] = x1[i][1]
                elif s1[i][0] == "leaving":
                    r[i][1][0] = x1[i][1]
                    r[i][1][1] = x1[i][0]
        r1 = [a for b in r for a in b]
        dic = dict(sorted(r1))
        for i in dic:
            dic[i] = [dic[i]]
        D = DiGraph(dic)
        d = D.all_simple_cycles()
        for i in d:
            del i[0]
        return d

    def regions(self):
        y = self.PD_code()
        x = deepcopy(y)
        x1 = deepcopy(x)
        s = [[None for i in range(4)] for i in range(len(x))]
        for i in range(len(x)):
            s[i][0] = 'entering'
            s[i][2] = 'leaving'
            if x[i][1] < x[i][3]:
                s[i][1] = 'entering'
                s[i][3] = 'leaving'
            if x[i][1] > x[i][3]:
                s[i][1] = 'leaving'
                s[i][3] = 'entering'
        ou = [["under","over","under","over"] for i in range(len(x))]
        r = [[[None,None],[None,None]] for i in range(len(x))]
        r1 = [[[None,None],[None,None]] for i in range(len(x))]
        for i in range(len(x)):
            r[i][0][0] = x[i][0]
            r1[i][0][0] = "under"
            for j in range(1,4):
                if s[i][j] == "entering":
                    r[i][0][1] = x[i][j]
                    r1[i][0][1] = ou[i][j]
                    del x[i][j]
                    del ou[i][j]
                    s[i][j] = 0
                    del ou[i][0]
                    s[i][0] = 0
                    del x[i][0]
        for i in range(len(x)):
            r[i][1][0] = x[i][0]
            r[i][1][1] = x[i][1]
            r1[i][1][0] = ou[i][0]
            r1[i][1][1] = ou[i][1]
            del x[i][1]
            del x[i][0]
            del ou[i][1]
            del ou[i][0]
        r = [a for b in r for a in b]
        r1 = [a for b in r1 for a in b]
        r3 = deepcopy(r)
        lc = [[None,None] for i in range(len(x))]
        for i in range(len(x)):
            if x1[i][1] > x1[i][3]:
                if r1[2*i+1][0] == 'over':
                    lc[i][0] = r[2*i+1][0]
                lc[i][1] = (-1)*r[2*i][0]
            elif x1[i][1] < x1[i][3]:
                if r1[2*i+1][0] == 'under':
                    lc[i][1] = r[2*i+1][0]
                lc[i][0] = (-1)*r[2*i][1]
        r2 = []
        for i in reversed(range(len(x))):
            r2.append(r[2*i])
            del r[2*i]
        entering = r2[::-1]
        left_component_entering = lc
        leaving = r
        r4 = deepcopy(r)
        for i in r4:
            i[0] = (-1)*i[0]
            i[1] = (-1)*i[1]
        lc1 = [[None,None] for i in range(len(x))]
        for i in range(len(x)):
            if x1[i][1] > x1[i][3]:
                if r1[2*i+1][0] == 'over':
                    lc1[i][0] = r3[2*i+1][1]
                    if r1[2*i][0] == 'over':
                        lc1[i][1] = (-1)*r3[2*i][0]
                    elif r1[2*i][1] == 'over':
                        lc1[i][1] = (-1)*r3[2*i][1]
            elif x1[i][1] < x1[i][3]:
                if r1[2*i+1][1] == 'over':
                    lc1[i][0] = r3[2*i+1][1]
                    if r1[2*i][0] == 'under':
                        lc1[i][1] = (-1)*r3[2*i][0]
                    elif r1[2*i][1] == 'under':
                        lc1[i][1] = (-1)*r3[2*i][1]
        left_component_leaving = lc1
        '''
        sage: L.regions([[-1, +2, -3, 1, -2, +3 ],['-','-','-']])
        [[1, 4], [5, 2], [3, 6]]
        [[5, -1], [3, -5], [1, -3]]
        [[5, 2], [3, 6], [1, 4]]
        read as 1 has left component as 5, 4 has left componet as -1, ...
        read as [1,4], [5,2], [3,6] enter the crossing
        read as [5,2], [3,6], [1,4] leave the crossing
        '''
        w = {}
        w1 = {}
        for i in range(len(x)):
            w.update({entering[i][0] : left_component_entering[i][0]})
            w.update({entering[i][1] : left_component_entering[i][1]})
            w1.update({r4[i][0] : left_component_leaving[i][0]})
            w1.update({r4[i][1] : left_component_leaving[i][1]})
        for i in range(1,len(w)+1):
            w[i] = [w[i]]
        for i in range(-len(w),0):
            w1[i] = [w1[i]]
        w2 = dict(w.items() + w1.items())
        d = DiGraph(w2)
        regions = d.all_simple_cycles()
        for i in regions:
            del i[0]
        check = [i for i in range(-len(w), len(w)+1)]
        del check[len(w)]
        regions1 = [a for b in regions for a in b]
        r5 = sorted(regions1)
        if r5 == check:
            return regions
        else:
            raise Exception("Incorrect Input")

    def remove_regions(self):
        y = self.PD_code()
        x = deepcopy(y)
        sc = self.seifert_circles()
        regions = self.regions()
        q1 = deepcopy(regions)
        q = [[[],[]] for i in range(len(regions))]
        for i in range(len(regions)):
            for j in range(len(regions[i])):
                if regions[i][j] < 0:
                    q[i][0].append(regions[i][j])
                elif regions[i][j] > 0:
                    q[i][1].append(regions[i][j])
        r = [[] for i in range(len(q))]
        for i in range(len(q)):
            for j in range(len(q[i])):
                r[i].append(None)
                for k in range(len(q[i][j])):
                    if q[i][j][k] < 0:
                        q[i][j][k] = (-1)*q[i][j][k]
        for i in range(len(q)):
            for j in range(len(q[i])):
                for k in sc:
                    if set(q[i][j]).intersection(set(k)) == set(q[i][j]):
                        r[i][j] = None
                        break
                else:
                    r[i][j] = q[i][j]
        for i in reversed(range(len(r))):
            for j in reversed(range(len(r[i]))):
                if r[i][j] == None:
                    del r[i][j]
         #we see whether they belong to the same region if it is so we remove the other.
        for i in reversed(range(len(r))):
            if len(r[i]) > 1:
                del r[i][1]
            if len(r[i]) == 0:
                del r[i]
        # r has the components together
        r = [a for b in r for a in b]
        # we flat r to get r1
        r1 = deepcopy(r)
        r1 = sum(r1,[])
        # we sort r1 to get r2 this is in order to computer c
        r2 = sorted(r1)
        # c has the new components after the move has been made
        c = [[None for i in range(3)] for i in range(len(r2))]
        for i in range(len(c)):
            c[i][0] = r2[i] + 2*i
            c[i][1] = c[i][0] + 1
            c[i][2] = c[i][1] + 1
        #here we update c1 after comparing r1 & r2
        c1 = []
        for i in range(len(r1)):
            for j in range(len(r2)):
                if r1[i] == r2[j]:
                    c1.append(c[j])
                    break
        #now using the above data in c1 compute the PD info at the new crossings
        #we always pull the component numbered higher onto the component numbered lower
        npd = [[None for i in range(4)] for i in range(len(r1))]
        for i in range(len(r)):
            if max(r[i]) == r1[2*i]:
                npd[2*i+1][0] = c1[2*i][0]
                npd[2*i+1][2] = npd[2*i][0] + 1
                npd[2*i][0] = npd[2*i][2]
                npd[2*i][2] = npd[2*i][0] + 1
            elif max(r[i]) == r1[2*i+1]:
                npd[2*i+1][0] = c1[2*i+1][0]
                npd[2*i+1][2] = npd[2*i+1][0] + 1
                npd[2*i][0] = npd[2*i+1][2]
                npd[2*i][2] = npd[2*i][0] + 1
        #now we need to decide the over crossings, for this we need to find which crossing has
        #the highest component leaving since that has the new lowest component entering
        for i in range(len(r)):
            if max(r[i]) == r1[2*i+1]:
                if max(c1[2*i+1]) == max(npd[2*i]):
                    npd[2*i][1] = c1[2*i][1]
                    npd[2*i][3] = c1[2*i][0]
                    npd[2*i+1][1] = c1[2*i][1]
                    npd[2*i+1][3] = c1[2*i][2]
                elif max(c1[2*i+1]) == max(npd[2*i+1]):
                    npd[2*i+1][1] = c1[2*i][1]
                    npd[2*i+1][3] = c1[2*i][0]
                    npd[2*i][1] = c1[2*i][1]
                    npd[2*i][3] = c1[2*i][2]
            elif max(r[i]) == r1[2*i]:
                if max(c1[2*i]) == max(npd[2*i]):
                    npd[2*i][1] = c1[2*i][1]
                    npd[2*i][3] = c1[2*i][0]
                    npd[2*i+1][1] = c1[2*i][1]
                    npd[2*i+1][3] = c1[2*i][2]
                elif max(c1[2*i]) == max(npd[2*i+1]):
                    npd[2*i+1][1] = c1[2*i][1]
                    npd[2*i+1][3] = c1[2*i][0]
                    npd[2*i][1] = c1[2*i][1]
                    npd[2*i][3] = c1[2*i][2]
        #here we computer the planar diagram and store it in y with the new components
        #the npd has the PD Code for new crossings, now we need to modify the previous
        #crossings.
        y = [[None for i in range(4)] for i in range(len(x))]
        for i in range(len(r2)):
            if i+1  < len(r2):
                for j in range(len(x)):
                    for k in range(4):
                        if x[j][k] < r2[0]:
                            y[j][k] = x[j][k]
                        elif x[j][k] > r2[i] and x[j][k] < r2[i+1]:
                            y[j][k] = x[j][k] + 2*(i+1)
                        elif x[j][k] > r2[len(r2)-1]:
                            y[j][k] = x[j][k] + 2*(i+2)
        #the above PD code has only few elements and not all, as we are updating the
        #previous crossings we can use the sign information from the input.
        #gc_sign = self.oriented_gauss_code[1]
        for i in range(len(x)):
            if x[i][1] > x[i][3]:
                if y[i][0] == None:
                    y[i][0] = y[i][2] - 1
                if y[i][2] == None:
                    y[i][2] = y[i][0] + 1
                if y[i][3] == None:
                    y[i][3] = y[i][1] - 1
                if y[i][1] == None:
                    y[i][1] = y[i][3] + 1
            elif x[i][1] < x[i][3]:
                if y[i][0] == None:
                    y[i][0] = y[i][2] - 1
                if y[i][2] == None:
                    y[i][2] = y[i][0] + 1
                if y[i][3] == None:
                    y[i][3] = y[i][1] + 1
                if y[i][1] == None:
                    y[i][1] = y[i][3] - 1
        y = y + npd
        return y
