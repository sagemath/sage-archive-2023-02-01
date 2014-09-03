r"""
Link class
"""
#*****************************************************************************
#  Copyright (C) 2014
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.integer_mod import Mod
from sage.graphs.digraph import DiGraph
from copy import deepcopy, copy
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import IntegerRing
from sage.combinat.permutation import Permutations
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.symbolic.ring import SR, var
from sage.rings.integer import Integer


class Link:

    r"""
    The base class for Link.

    INPUT:

    - The different ways in which the input can be provided :

      1. Braid
      2. Oriented Gauss Code
      3. Planar Diagram Code

      EXAMPLES::

          sage: B = BraidGroup(8)
          sage: L = Link(B([1, 2, 1, -2,-1]))
          sage: L
          Link with 2 components represented by 5 crossings
          sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
          sage: L
          Knot represented by 7 crossings
          sage: L = Link([[1,8,2,7],[8,4,9,5],[3,9,4,10],[10,1,7,6],[5,3,6,2]])
          sage: L
          Link with 2 components represented by 5 crossings
    """

    def __init__(self, input_):
        r"""

        The Python constructor.

        A Link can be created by using one of the conventions mentioned below:

        Braid:
        =========

        Generators of the braid group are used to generate the link.

            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, -1, -1, -2,1, -2,3,-2,3]))
            sage: L
            Knot represented by 9 crossings
            sage: L = Link(B([1, 2,1, -2,-1]))
            sage: L
            Link with 2 components represented by 5 crossings

        Oriented Gauss Code:
        ===================

        Randomly number the crossings from 1 to n (where n is the number of
        crossings) and start moving along the link. Trace every component of
        the link, by starting at a particular point on one component of the link and
        taking note of each of the crossings until one returns to the starting
        point. Note each component as a list whose elements are the crossing
        numbers. Compile every component info into a list. We need the orientation
        of every crossing. This is recorded as a list with +1 and -1, +1 is recorded
        if the direction from leaving over-cross to the leaving under-cross is
        anti-clockwise, -1 if the direction from the leaving over-cross to the
        leaving under-cross is clockwise.

            # for knots there is only a single component so the input is as follows
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1,-1,-1,-1,+1,+1,-1,+1]])
            sage: L
            Knot represented by 8 crossings

            # for links there is more than one component and the input is as follows
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L
            Link with 3 components represented by 4 crossings

        Planar Diagram Code:
        ===================

        Select some point on the link. Start numbering the strands in the
        components of the link. For a new component add one to the greatest
        number from the previous component and proceed till all the strands
        are numbered. At every cross contruct the data as follows :
        Start with the strand number of the entering under-cross and move in the
        clockwise direction around the cross and note down the strand numbers.
        Construct this data at every crossing and that would give the PD-Code.

            # there is no particular distinction between knots and links for this input

            # One of the representations of the Trefoil knot
            sage: L = Link([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: L
            Knot represented by 3 crossings

            # One of the representations of the Hopf link
            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L
            Link with 2 components represented by 2 crossings
        """
        if type(input_) == list:
            if len(input_) != 2:
                pd_error = _pd_check_(input_)
                if pd_error == True:
                    raise Exception("Invalid Input")
                elif pd_error == False:
                    self._PD_code = input_
                    self._oriented_gauss_code = None
                    self._braid = None

            elif len(input_) == 2:
                for i in input_[0]:
                    if type(i) == list:
                        ogc = True
                        break
                else:
                    ogc = False
                if ogc == False:
                    pd_error = _pd_check_(input_)
                    if pd_error == True:
                        raise Exception("Invalid Input")
                    elif pd_error == False:
                        self._PD_code = input_
                        self._oriented_gauss_code = None
                        self._braid = None
                elif ogc == True:
                    for i in input_[0]:
                        if type(i) != list:
                            raise Exception("Invalid Input")
                    else:
                        flat = [x for y in input_[0] for x in y]
                        a, b = max(flat), min(flat)
                        if 2 * len(input_[1]) == len(flat) and set(range(b, a + 1)) - set([0]) == set(flat):
                            self._oriented_gauss_code = input_
                            self._PD_code = None
                            self._braid = None
                        else:
                            raise Exception("Invalid Input")
        else:
            from sage.groups.braid import Braid
            if isinstance(input_, Braid):
                self._braid = input_
                self._oriented_gauss_code = None
                self._PD_code = None

            else:
                raise Exception("Invalid Input")

    def __repr__(self):
        r"""
        Return a string representation.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L
            Knot represented by 4 crossings
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L
            Knot represented by 7 crossings
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L
            Link with 3 components represented by 4 crossings
        """
        ncomponents = self.ncomponents()
        pd_len = len(self.PD_code())
        if self.is_knot():
            return 'Knot represented by {} crossings'.format(pd_len)
        else:
            return 'Link with {} components represented by {} crossings'.format(ncomponents, pd_len)

    def braidword(self):
        r"""
        Return the braidword of the link.

        OUTPUT:

            - Braidword representation of the link.

        EXAMPLES::

            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1,-1,-1,-1,+1,+1,-1,+1]])
            sage: L.braidword()
            (-1, 2, -1, -2, -2, 1, 1, -2)
            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L.braidword()
            (-1, -1)
            sage: L = Link([[[1, -2, 3, -4], [-1, 5, -3, 2, -5, 4]], [-1, 1, 1, -1, -1]])
            sage: L.braidword()
            (1, -2, 1, -2, -2)
        """
        return self.braid().Tietze()

    def braid(self):
        r"""
        Return the braid representation of the link.

        OUTPUT:
            - Braid representation of the link.

        EXAMPLES::

            sage: L = Link([[2,3,1,4],[4,1,3,2]])
            sage: L.braid()
            s^2
            sage: L = Link([[[-1, 2, -3, 1, -2, 3]], [-1, -1, -1]])
            sage: L.braid()
            s^-3
            sage: L = Link([[1,8,2,7],[8,4,9,5],[3,9,4,10],[10,1,7,6],[5,3,6,2]])
            sage: L.braid()
            (s0*s1^-1)^2*s1^-1
        """
        if self._braid != None:
            return self._braid

        elif self._oriented_gauss_code != None or self._PD_code != None:
            braid_detection = self._braidword_detection_()
            gen = max([abs(i) for i in braid_detection])
            from sage.groups.braid import BraidGroup
            B = BraidGroup(gen + 1)
            self._braid = B(braid_detection)
            return self._braid

    def oriented_gauss_code(self):
        r"""
        Return the oriented gauss code of the link. The oriented gauss
        code has two parts
        a. The gauss code
        b. The orientation of each crossing
        The following orientation was taken into consideration for consturction
        of knots:

        From the outgoing of the overcrossing if we move in the clockwise direction
        to reach the outgoing of the undercrossing then we label that crossing as '-1'.

        From the outgoing of the overcrossing if we move in the anticlockwise
        direction to reach the outgoingo the undercrossing then we label that crossing
        as '+1'.

        One more consideration we take in while constructing the orientation is:
        The order of the orientation is same as the ordering of the crossings in the
        gauss code.

        Convention : under is denoted by -1, and over by +1 in the crossing info.

        OUTPUT:
            - Oriented gauss code of the link

        EXAMPLES::

            sage: L = Link([[1,11,2,10],[6,2,7,3],[3,12,4,9],[9,5,10,6],[8,1,5,4],[11,8,12,7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]]
            sage: L = Link([[1,4,2,3],[6,1,3,2],[7,4,8,5],[5,8,6,7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]]
            sage: B = BraidGroup(8)
            sage: b=B([1,1,1,1,1])
            sage: L = Link(b)
            sage: L.oriented_gauss_code()
            [[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]]
        """
        if self._oriented_gauss_code != None:
            return self._oriented_gauss_code

        if self._oriented_gauss_code == None:
            self.PD_code()
            pd = self._PD_code
            orient = self.orientation()
            crossing_info = {}
            for i, j in enumerate(pd):
                if orient[i] == -1:
                    crossing_info.update({(j[0], -1, i + 1): j[2]})
                    crossing_info.update({(j[3], 1, i + 1): j[1]})
                elif orient[i] == 1:
                    crossing_info.update({(j[0], -1, i + 1): j[2]})
                    crossing_info.update({(j[1], 1, i + 1): j[3]})
            edges = {}
            cross_number = {}
            for i, j in crossing_info.items():
                edges.update({i[0]: [j]})
                if i[1] == 1:
                    cross_number.update({i[0]: i[2]})
                elif i[1] == -1:
                    cross_number.update({i[0]: -i[2]})
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
            self._oriented_gauss_code = oriented_code
            return self._oriented_gauss_code

    def PD_code(self):
        r"""
        Return the Planar Diagram code of the link. The Planar Diagram is returned
        in the following format.

        We construct the crossing by starting with the entering component of the
        undercrossing, move in the clockwise direction and then generate the list.
        Suppose if the crossing is given by [a, b, c, d], then we interpret this
        information as
        a is the entering component of the undercrossing
        b, d are the components of the overcrossing
        c is the leaving component of the undercrossing

        Convention :
        Orientation of the crossing:
        Leaving over crossing to leaving under crossing in clockwise direction denoted by -1
        Leaving over crossing to leaving under crossing in anticlockwise direction denoted by 1

        OUTPUT:
            - Planar Diagram representation of the link.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1,1,-1,-1]])
            sage: L.PD_code()
            [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]]
            sage: B = BraidGroup(2)
            sage: b=B([1,1,1,1,1])
            sage: L = Link(b)
            sage: L.PD_code()
            [[2, 1, 3, 4], [4, 3, 5, 6], [6, 5, 7, 8], [8, 7, 9, 10], [10, 9, 1, 2]]
            sage: L = Link([[[2, -1], [1, -2]], [1, 1]])
            sage: L.PD_code()
            [[2, 3, 1, 4], [4, 1, 3, 2]]
        """
        if self._PD_code != None:
            return self._PD_code

        elif self._oriented_gauss_code != None:
            oriented_gauss_code = self._oriented_gauss_code
            d_dic = {}
            if len(oriented_gauss_code[0]) > 1:
                d = [x for y in oriented_gauss_code[0] for x in y]
                for i, j in enumerate(d):
                    d_dic.update({j: [i + 1, i + 2]})
                # here we collect the final component in each gauss code
                last_component = [i[len(i) - 1]
                                  for i in oriented_gauss_code[0]]
                first_component = [i[0] for i in oriented_gauss_code[0]]
                # here we correct the last_component
                for i, j in zip(last_component, first_component):
                    d_dic[i][1] = d_dic[j][0]
                crossing_dic = {}
                for i in range(len(oriented_gauss_code[1])):
                    if oriented_gauss_code[1][i] == -1:
                        crossing_dic.update(
                            {i + 1: [d_dic[-(i + 1)][0], d_dic[i + 1][1], d_dic[-(i + 1)][1], d_dic[i + 1][0]]})
                    elif oriented_gauss_code[1][i] == 1:
                        crossing_dic.update(
                            {i + 1: [d_dic[-(i + 1)][0], d_dic[i + 1][0], d_dic[-(i + 1)][1], d_dic[i + 1][1]]})
            elif len(oriented_gauss_code[0]) == 1:
                for i, j in enumerate(oriented_gauss_code[0][0]):
                    d_dic.update({j: [i + 1, i + 2]})
                d_dic[oriented_gauss_code[0][0][-1]][1] = 1
                crossing_dic = {}
                for i in range(len(oriented_gauss_code[1])):
                    if oriented_gauss_code[1][i] == -1:
                        crossing_dic.update(
                            {i + 1: [d_dic[-(i + 1)][0], d_dic[i + 1][1], d_dic[-(i + 1)][1], d_dic[i + 1][0]]})
                    elif oriented_gauss_code[1][i] == 1:
                        crossing_dic.update(
                            {i + 1: [d_dic[-(i + 1)][0], d_dic[i + 1][0], d_dic[-(i + 1)][1], d_dic[i + 1][1]]})
            pd = crossing_dic.values()
            self._PD_code = pd
            return self._PD_code

        elif self._braid != None:
            strings = range(1, self._braid.strands() + 1)
            b = list(self._braid.Tietze())
            pd = []
            strings_max = strings[-1]
            for i in b:
                if i > 0:
                    pd.append([strings[i], strings[i-1], strings_max + 1, strings_max + 2])
                else:
                    pd.append([strings[abs(i)-1], strings_max + 1, strings_max + 2, strings[abs(i)]])
                strings[abs(i)-1] = strings_max + 1
                strings[abs(i)] = strings_max + 2
                strings_max = strings_max + 2
            for i in pd:
                for j in range(4):
                    if i[j] in strings:
                        i[j] = strings.index(i[j]) + 1
            self._PD_code = pd
            return pd

    def gauss_code(self):
        r"""
        Return the gauss_code of the link. Gauss code is generated by the
        following procedure:

        a. We randomly number the crossings
        b. We select a point on the knot and start moving along the component
        c. At each crossing we take the number of the crossing, along with
           sign, which is '-' if it is a undercrossing and '+' if it is a
           overcrossing.

        OUTPUT:
            - Gauss code representation of the link.

        EXAMPLES::

            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L.gauss_code()
            [[-1, 2], [1, -2]]
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, -2, 1, -2, -2]))
            sage: L.gauss_code()
            [[-1, 3, -4, 5], [1, -2, 4, -5, 2, -3]]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L.gauss_code()
            [[-1, 2], [-3, 4], [1, 3, -4, -2]]
        """
        return self.oriented_gauss_code()[0]

    def dt_code(self):
        r"""
        Return the dt_code of the knot.DT code is generated by the following way.

        We start moving along the knot, as we encounter the crossings we
        start numbering them, so every crossing has two numbers assigned to
        it once we have traced the entire knot. Now we take the even number
        associated with every crossing. The following sign convention is to be
        followed:
        Take the even number with a negative sign if it is an overcrossing
        that we are encountering.

        OUTPUT:
            - DT Code representation of the knot. This is implemented only
              for knots.

        EXAMPLES::

            sage: L = Link([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: L.dt_code()
            [4, 6, 2]
            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.dt_code()
            [4, -6, 8, -2]
            sage: L = Link([[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]])
            sage: L.dt_code()
            [6, 8, 10, 2, 4]
        """
        def _dt_internal_(b):
            N = len(b)
            label = [0 for i in range(2 * N)]
            string = 1
            next_label = 1
            type1 = 0
            crossing = 0
            while(next_label <= 2 * N):
                string_found = 0
                for i in range(crossing, N):
                    if(abs(b[i]) == string or abs(b[i]) == string - 1):
                        string_found = 1
                        crossing = i
                        break
                if(string_found == 0):
                    for i in range(0, crossing):
                        if(abs(b[i]) == string or abs(b[i]) == string - 1):
                            string_found = 1
                            crossing = i
                            break
                if(label[2 * crossing + next_label % 2] == 1):
                    raise Exception("Implemented only for knots")
                else:
                    label[2 * crossing + next_label % 2] = next_label
                    next_label = next_label + 1
                if(type1 == 0):
                    if(b[crossing] < 0):
                        type1 = 1
                    else:
                        type1 = -1
                else:
                    type1 = -1 * type1
                    if((abs(b[crossing]) == string and b[crossing] * type1 > 0) or (abs(b[crossing]) != string and b[crossing] * type1 < 0)):
                        if(next_label % 2 == 1):
                            label[2 * crossing] = label[2 * crossing] * -1
                if(abs(b[crossing]) == string):
                    string = string + 1
                else:
                    string = string - 1
                crossing = crossing + 1
            code = [0 for i in range(N)]
            for i in range(N):
                for j in range(N):
                    if label[2 * j + 1] == 2 * i + 1:
                        code[i] = label[2 * j]
                        break
            return code

        if self._braid != None:
            b = self.braidword()
        else:
            b = self._braidword_detection_()
        return _dt_internal_(b)

    def _dowker_notation_(self):
        r"""
        Return the dowker notation of the link. Similar to the PD Code we number the
        components. So every crossing is represented by four numbers. We focus on
        the incoming entites of the under and the over crossing. It is the pair of incoming
        under cross and the incoming over cross. This information at every cross
        gives the dowker notation.

        OUTPUT:
            - List containing the pair of incoming under cross and the incoming
              over cross.

        EXAMPLES::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L._dowker_notation_()
            [(1, 6), (7, 2), (3, 10), (11, 4), (14, 5), (13, 8), (12, 9)]
            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L._dowker_notation_()
            [(2, 1), (3, 5), (6, 4), (7, 9)]
            sage: L = Link([[1,4,2,3],[4,1,3,2]])
            sage: L._dowker_notation_()
            [(1, 3), (4, 2)]
        """
        pd = self.PD_code()
        orient = self.orientation()
        dn = [(i[0], i[3]) if orient[j] == -1 else (i[0], i[1])
              for j, i in enumerate(pd)]
        return dn

    def _braidwordcomponents_(self):
        r"""
        Return the disjoint braid components, if any, else returns the braid itself.
        For example consider the braid [-1, 3, 1, 3] this can be viewed as a braid
        with components as [-1, 1] and [3, 3]. There is no common crossing to these
        two (in sense there is a crossing between strand 1 and 2, crossing between
        3 and 4 but no crossing between strand 2 and 3,so these can be viewed as
        independent components in the braid).

        OUTPUT:
            - List containing the components is returned

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponents_()
            [[-1, 1], [3, 3]]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponents_()
            [[-1, 1, 1, 1], [3], [5, 7, 6]]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponents_()
            [[-2, 1, 1], [4, 4], [6]]
        """
        b = self.braid()
        ml = list(b.Tietze())
        if ml == []:
            raise Exception("The braid remains the same with no components")
        else:
            l = list(set([abs(k) for k in ml]))
            missing1 = list(set(range(min(l), max(l) + 1)) - set(l))
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

    def _braidwordcomponentsvector_(self):
        r"""
        The list from the braidwordcomponents is flattened to give out the vector form.

        OUTPUT:
            - Vector containing braidwordcomponents

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braidwordcomponentsvector_()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braidwordcomponentsvector_()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braidwordcomponentsvector_()
            [-2, 1, 1, 4, 4, 6]
        """
        bc = self._braidwordcomponents_()
        return [x for y in bc for x in y]

    def _homology_generators_(self):
        r"""
        The set of generators for the first homology group of the connected Seifert surface of the given link.
        This method uses the braidwordcomponentsvector to generate the homology generators.
        The position of the repeated element w.r.t the braidwordcomponentvector list is
        compiled into a list.

        OUTPUT:
            - The homology generators relating to the braid word representation

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._homology_generators_()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._homology_generators_()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._homology_generators_()
            [0, 2, 0, 4, 0]
        """
        x4 = self._braidwordcomponentsvector_()
        hom_gen = []
        for j in range(len(x4) - 1):
            a = abs(x4[j])
            for i in range(j + 1, len(x4)):
                if(a == abs(x4[i])):
                    hom_gen.append(i)
                    break
            else:
                hom_gen.append(0)
        return hom_gen

    def Seifert_Matrix(self):
        r"""
        Return the Seifert Matrix associated with the braidword.

        OUTPUT:
            - The intersection matrix of a (not necessarily minimal) Seifert surface of the
              link.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.Seifert_Matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.Seifert_Matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.Seifert_Matrix()
            [-1  0]
            [ 0 -1]
        """
        x5 = self._braidwordcomponentsvector_()
        h = self._homology_generators_()
        hl = len(h)
        A = matrix(ZZ, hl, hl)
        for i in range(hl):
            if h[i] != 0:
                for j in range(i, hl):
                    if i == j:
                        A[i, j] = -cmp((x5[i] + x5[h[i]]), 0)
                    elif (h[i] > h[j]):
                        A[i, j] = 0
                        A[j, i] = 0
                    elif (h[i] < j):
                        A[i, j] = 0
                        A[j, i] = 0
                    elif (h[i] == j):
                        if(x5[j] > 0):
                            A[i, j] = 0
                            A[j, i] = 1
                        else:
                            A[i, j] = -1
                            A[j, i] = 0
                    elif abs(abs(x5[i]) - abs(x5[j])) > 1:
                        A[i, j] = 0
                    elif (abs(x5[i]) - abs(x5[j]) == 1):
                        A[i, j] = 0
                        A[j, i] = -1
                    elif (abs(x5[j]) - abs(x5[i]) == 1):
                        A[i, j] = 1
                        A[j, i] = 0
                    else:  # for debugging
                        A[i, j] = 2
                        A[j, i] = 2
            else:
                for k in range(hl):
                    A[k, i] = 0
                    A[i, k] = 0
        k = []
        for i in range(hl):
            if h[i] == 0:
                k.append(i)
        for i in reversed(k):
            A = A.delete_rows([i])
            A = A.delete_columns([i])
        return A

    def ncomponents(self):
        r"""
        Return the number of connected components of the link.

        OUTPUT:
            - Connected components of the link

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.ncomponents()
            4
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.ncomponents()
            5
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.ncomponents()
            1
        """
        p = self.braid().permutation()
        return len(p.to_cycles())

    def is_knot(self):
        r"""
        Return True if the link is knot.
        Every knot is a link but the converse is not true.

        OUTPUT:
            - True if knot else False

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1,3,1,-3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        if self.ncomponents() == 1:
            return True
        else:
            return False

    def genus(self):
        r"""
        Return the genus of the link

        OUTPUT:
            - Genus of the Link

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.genus()
            0
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.genus()
            1
        """
        b = self.braid().Tietze()
        if b == []:
            return 0
        else:
            B = self.braid().parent()
            x = self._braidwordcomponents_()
            q = []
            genus = 0
            s_tmp = []
            for i in range(len(x)):
                tmp = []
                b1 = min([abs(k) for k in x[i]])
                for j in range(len(x[i])):
                    if x[i][j] > 0:
                        x[i][j] = x[i][j] - b1 + 1
                    else:
                        x[i][j] = x[i][j] + b1 - 1
                    tmp.append(x[i][j])
                s_tmp.append(B(tmp))
            s = []
            for i in s_tmp:
                b = i.Tietze()
                s.append(list(b))
            t = [Link(B(s[i])).ncomponents() for i in range(len(s))]
            for i, j in enumerate(s):
                if j == []:
                    s[i].append(-2)
            for i in s:
                q1 = (abs(k) + 1 for k in i)
                q2 = max(q1)
                q.append(q2)
            g = [((2 - t[i]) + len(x[i]) - q[i]) / 2 for i in range(len(x))]
            for i in range(len(g)):
                genus = genus + g[i]
            return genus

    def signature(self):
        r"""
        Return the signature of the link

        OUTPUT:
            - Signature of the Link

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.signature()
            -1
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.signature()
            -2
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.signature()
            -2
        """
        m = 2 * (self.Seifert_Matrix() + self.Seifert_Matrix().transpose())
        e = m.eigenvalues()
        sum = 0
        s = []
        for i, j in enumerate(e):
            s.append(cmp(j, 0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var='t'):
        r"""
        Return the alexander polynomial of the link

        INPUT:
            - ``var`` -- string (default: ``'t'``); the name of the
                variable in the entries of the matrix

        OUTPUT:
            - Alexander Polynomial of the Link

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            t^-1 - 2 + t
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^-1 - 1 + t
        """
        R = LaurentPolynomialRing(ZZ, var)
        t = R.gen()
        f = (self.Seifert_Matrix() - t *
             (self.Seifert_Matrix().transpose())).determinant()
        return t ** ((-max(f.exponents()) - min(f.exponents())) / 2) * f if f != 0 else f

    def knot_determinant(self):
        r"""
        Return the determinant of the knot

        OUTPUT:
            - Determinant of the Knot

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.knot_determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.knot_determinant()
            3
            sage: L = Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.knot_determinant()
            65
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            return Integer(abs(a(-1)))
        else:
            raise Exception("Determinant implmented only for knots")

    def arf_invariant(self):
        r"""
        Return the arf invariant. Arf invariant is defined only for knots.

        OUTPUT:
            - Arf invariant of knot

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.arf_invariant()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 3, 1, 2, 1, 4]))
            sage: L.arf_invariant()
            0
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.arf_invariant()
            1
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            if ((Mod(a(-1), 8) == 1) or (Mod(a(-1), 8) == 7)):
                return 0
            else:
                return 1
        else:
            raise Exception("Arf invariant is defined only for knots")

    def is_alternating(self):
        r"""
        Return True if the given knot diagram is alternating else returns False.
        Alternating diagram implies every over cross is followed by an under cross
        or the vice-versa.

        We look at the gauss code if the sign is alternating, True is returned else
        the knot is not alternating False is returned.

        OUTPUT:
            - True if the knot diagram is alternating else False

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, -1, -1, -1]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1, -2, -1, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1, 3, 1,3, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1,2,-1,2]))
            sage: L.is_alternating()
            True
        """
        if self.is_knot() == True:
            x = self.gauss_code()
            s = [cmp(i, 0) for i in x[0]]
            if s == [(-1) ** (i + 1) for i in range(len(x[0]))] or s == [(-1) ** i for i in range(len(x[0]))]:
                return True
            else:
                return False
        else:
            return False

    #*************************************************************************
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
    #*************************************************************************

    #**************************** PART - 1 ***********************************
    def orientation(self):
        r"""
        Return the orientation of the crossings from the input. We construct the entering, leaving information at
        each crossing to get to the orientation.

        OUTPUT:
            - Orientation  of the crossings.

        Convention :
        Entering component is denoted by 1
        Leaving component is denoted by -1

        EXAMPLES::

            sage: L = Link([[1, 4, 5, 2], [3, 5, 6, 7], [4, 8, 9, 6], [7, 9, 10, 11], [8, 1, 13, 10], [11, 13, 2, 3]])
            sage: L.orientation()
            [-1, 1, -1, 1, -1, 1]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.orientation()
            [-1, -1, -1, -1, 1, -1, 1]
        """
        y = self.PD_code()
        x = deepcopy(y)
        under = [[[i[0], 1], [i[2], -1]] for i in x]
        under = [a for b in under for a in b]
        over = [[[i[1], None], [i[3], None]] for i in x]
        over = [a for b in over for a in b]
        for i in over:
            for j in under:
                if i[0] == j[0]:
                    if j[1] == 1:
                        i[1] = -1
                    elif j[1] == -1:
                        i[1] = 1
        for i in over:
            if i[1] == None:
                over = _rule_1_(over)
                over = _rule_2_(over)
        unfilled = []
        for i in over:
            if i[1] == None:
                unfilled.append(i)
        if len(unfilled) != 0:
            over[over.index(unfilled[0])][1] = 1
            for i in over:
                if i[1] == None:
                    over = _rule_1_(over)
                    over = _rule_2_(over)
        orientation = []
        for i in range(0, len(over), 2):
            if over[i][1] == -1:
                orientation.append(-1)
            elif over[i][1] == 1:
                orientation.append(1)
        return orientation

    def seifert_circles(self):
        r"""
        Return the seifert circles from the input. Seifert circles are obtained
        by smoothing the crossings in the following way:
        All the crossings are assigned four numbers which form the components of the
        crossings.
        At a crossing the under cross entering component would go to the over cross
        leaving component and the over cross entering component would go to the
        under cross leaving component.
        We start with a component of the crossing and start smoothing as according to
        the above rule until we return to the starting component.

        OUTPUT:
            - Seifert circles of the given knot.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[6, 2], [8, 4], [7, 5, 3, 1]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[10, 6, 12, 2], [16, 8, 14, 4], [13, 9, 3, 15, 5, 11, 7, 1]]
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1,-1,-1,-1,1,-1,1]])
            sage: L.seifert_circles()
            [[13, 9], [12, 10, 4], [8, 14, 6, 2], [7, 3, 11, 5, 1]]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.seifert_circles()
            [[13, 9], [12, 10, 4], [8, 14, 6, 2], [7, 3, 11, 5, 1]]
            sage: L = Link([[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[11, 8, 1], [9, 6, 3], [7, 12, 4, 5, 10, 2]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([1, 1, 1]))
            sage: L.seifert_circles()
            [[3, 5, 1], [4, 6, 2]]
        """
        pd = self.PD_code()
        orient = self.orientation()
        seifert_pairs = []
        for i, j in enumerate(pd):
            x = [[], []]
            if orient[i] == -1:
                x[0].append(j[0])
                x[0].append(j[1])
                x[1].append(j[3])
                x[1].append(j[2])
            elif orient[i] == 1:
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
            dic.update({flatten[i]: [flatten[i + 1]]})
        D = DiGraph(dic)
        d = D.all_simple_cycles()
        for i in d:
            del i[0]
        return d

    def regions(self):
        r"""
        Return the regions from the input.Regions are obtained always turning left at the crossing.
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

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.regions()
            [[4, -11], [2, -7], [6, -1], [13, 9], [-4, -10, -12], [-8, -2, -6, -14], [10, -3, 8, -13], [14, -5, 12, -9], [7, 3, 11, 5, 1]]
            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.regions()
            [[-2, -6], [8, 4], [5, 3, -8], [2, -5, -7], [1, 7, -4], [6, -1, -3]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.regions()
            [[6, -11], [15, -4], [9, 3, -14], [2, -9, -13], [1, 13, -8], [12, -1, -7], [5, 11, 7, -16], [-3, 10, -5, -15], [-6, -10, -2, -12], [16, 8, 14, 4]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([-1, -1, -1]))
            sage: L.regions()
            [[6, -5], [4, -3], [2, -1], [-4, -2, -6], [3, 5, 1]]
            sage: L = Link([[[1, -2, 3, -4], [-1, 5, -3, 2, -5, 4]], [-1, 1, 1, -1, -1]])
            sage: L.regions()
            [[1, -5], [8, 3], [-6, -1, -10], [-2, 6, -9], [10, -4, -7], [9, 7, -3], [4, 5, 2, -8]]
        """
        pd = self.PD_code()
        orient = self.orientation()
        regions = []
        for i, j in enumerate(pd):
            x = [[], []]
            if orient[i] == -1:
                x[0].append(j[0])
                x[0].append(j[1])
                x[1].append(j[3])
                x[1].append(-j[0])
            elif orient[i] == 1:
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
        for i, j in enumerate(pd_opposite):
            x_op = [[], []]
            if orient[i] == -1:
                x_op[0].append(j[0])
                x_op[0].append(j[1])
                x_op[1].append(j[3])
                x_op[1].append(-j[0])
            elif orient[i] == 1:
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
            dic.update({i[0]: [i[1]]})
        D = DiGraph(dic)
        d = D.all_simple_cycles()
        for i in d:
            del i[0]
        return d

    def _vogel_move_(self):
        r"""
        Return the Planar Diagram code if there is a vogel's move required, else returns
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

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L._vogel_move_()
            'No Vogel Move'
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._vogel_move_()
            'No Vogel Move'
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L._vogel_move_()
            [[1, 7, 2, 6], [7, 3, 8, 2], [16, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [18, 9, 14, 8], [12, 9, 13, 10], [13, 15, 17, 16], [17, 15, 18, 3]]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L._vogel_move_()
            [[1, 7, 2, 6], [7, 3, 8, 2], [16, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [18, 9, 14, 8], [12, 9, 13, 10], [13, 15, 17, 16], [17, 15, 18, 3]]
            sage: L = Link([[1,4,2,3],[6,1,3,2],[7,4,8,5],[5,8,6,7]])
            sage: L._vogel_move_()
            [[1, 4, 2, 3], [6, 1, 3, 10], [12, 4, 8, 5], [5, 8, 6, 7], [7, 10, 11, 9], [11, 2, 12, 9]]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L._vogel_move_()
            [[1, 6, 2, 5], [8, 1, 5, 10], [12, 6, 4, 7], [7, 4, 8, 3], [3, 10, 11, 9], [11, 2, 12, 9]]
        """
        pd = self.PD_code()
        pd_copy = deepcopy(pd)
        sc = self.seifert_circles()
        regions = self.regions()
        regions_copy = deepcopy(regions)
        orient = self.orientation()
        # separating the components into positive and negative
        q = [[[], []] for i in range(len(regions))]
        for i, j in enumerate(regions):
            for k in j:
                if k < 0:
                    q[i][0].append(k)
                elif k > 0:
                    q[i][1].append(k)
        # making all the components positive
        r = [[] for i in range(len(q))]
        for i in range(len(q)):
            for j in range(len(q[i])):
                r[i].append(None)
                for k in range(len(q[i][j])):
                    if q[i][j][k] < 0:
                        q[i][j][k] = (-1) * q[i][j][k]
        # to find the intersection of regions with seifert circles
        # first clean q that is by removing empty and single length arrays
        clean_q = []
        for i in q:
            for j in i:
                if len(j) >= 2:
                    clean_q.append(j)
        # detecting the bad region
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
            # here there might be many bad regions but we only select one
            # here it is the first one
            bad_region = bad_region[0]
            # finding the max of the pd to develop the new crossings
            pd_max = max([max(i) for i in pd_copy])
            # editing the previous crossings
            # the maximum is corrected
            for i in pd_copy:
                if max(bad_region) == i[0]:
                    i[0] = pd_max + 4
            # editing the contents of the pd code with the minimum
            # the minimum is corrected
            min_cross = [[j, orient[i]]
                         for i, j in enumerate(pd_copy) if min(bad_region) in j]
            y = [[] for i in range(len(min_cross))]
            for i, j in enumerate(min_cross):
                if j[1] == -1:
                    y[i].append((j[0][0], 1))
                    y[i].append((j[0][1], -1))
                    y[i].append((j[0][2], -1))
                    y[i].append((j[0][3], 1))
                if j[1] == 1:
                    y[i].append((j[0][0], 1))
                    y[i].append((j[0][1], 1))
                    y[i].append((j[0][2], -1))
                    y[i].append((j[0][3], -1))
            for i in y:
                if (min(bad_region), 1) in i:
                    pd_copy[pd_copy.index(min_cross[y.index(i)][0])][
                        (pd_copy[pd_copy.index(min_cross[y.index(i)][0])]).index(min(bad_region))] = pd_max + 2
            # sorting the regions in positive and negative ones.
            pos_neg = [[[], []] for i in range(len(regions_copy))]
            for i, j in enumerate(regions_copy):
                for k in j:
                    if k > 0:
                        pos_neg[i][0].append(k)
                    elif k < 0:
                        pos_neg[i][1].append(k)
            pos = [i[0] for i in pos_neg if i[0] != []]
            neg = [[-j for j in i[1]] for i in pos_neg if i[1] != []]
            # creating the new components, the lesser in the bad region is
            # always over and the greater is always under
            new_component_1 = [min(bad_region), pd_max + 1, pd_max + 2]
            new_component_2 = [max(bad_region), pd_max + 3, pd_max + 4]
            # creating new crossings
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

    def _info_all_moves_(self):
        r"""
        Return the Planar Diagram code, Seifert circles, regions, orientation after all the bad regions
        have been removed by performing the vogel's move.

        OUTPUT:
            - Planar Diagram code, Seifert circle, Regions, orientation after the bad regions have been
              removed

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L._info_all_moves_()
            [[[6, 2], [8, 4], [7, 5, 3, 1]],
            [[-2, -6], [8, 4], [5, 3, -8], [2, -5, -7], [1, 7, -4], [6, -1, -3]],
            [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]],
            [1, 1, -1, -1]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._info_all_moves_()
            [[[10, 6, 12, 2], [16, 8, 14, 4], [13, 9, 3, 15, 5, 11, 7, 1]],
            [[6, -11], [15, -4], [9, 3, -14], [2, -9, -13], [1, 13, -8], [12, -1, -7], [5, 11, 7, -16], [-3, 10, -5, -15], [-6, -10, -2, -12], [16, 8, 14, 4]],
            [[1, 13, 2, 12], [9, 3, 10, 2], [14, 4, 15, 3], [4, 16, 5, 15], [10, 5, 11, 6], [6, 11, 7, 12], [16, 8, 1, 7], [13, 8, 14, 9]],
            [-1, -1, -1, -1, 1, 1, -1, 1]]
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L._info_all_moves_()
            [[[17, 15], [21, 19], [7, 3, 18, 9, 13, 16, 11, 5, 1], [8, 14, 20, 10, 4, 12, 22, 6, 2]],
            [[-19, -21], [4, -11], [2, -7], [6, -1], [17, 15], [-3, 8, -18], [-13, 10, -16], [14, 20, -9], [12, 22, -5], [18, 9, 13, -15], [21, -12, -4, -10, -20], [19, -14, -8, -2, -6, -22], [16, 11, 5, 1, 7, 3, -17]],
            [[1, 7, 2, 6], [7, 3, 8, 2], [16, 11, 4, 10], [11, 5, 12, 4], [22, 5, 1, 6], [18, 9, 14, 8], [20, 9, 13, 10], [13, 15, 17, 16], [17, 15, 18, 3], [14, 20, 21, 19], [21, 12, 22, 19]],
            [-1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1]]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L._info_all_moves_()
            [[[11, 9], [6, 4, 8, 1], [12, 7, 3, 10, 5, 2]],
            [[-9, -11], [7, -4], [5, -1], [3, 10, -8], [2, 12, -6], [9, -3, -7, -12], [11, -2, -5, -10], [6, 4, 8, 1]],
            [[1, 6, 2, 5], [8, 1, 5, 10], [12, 6, 4, 7], [7, 4, 8, 3], [3, 10, 11, 9], [11, 2, 12, 9]],
            [-1, -1, 1, 1, -1, 1]]
        """
        x = self.PD_code()
        while True:
            link = Link(x)
            PD_code_old = x
            x = link._vogel_move_()
            if x == "No Vogel Move":
                x = PD_code_old
                break
        L = Link(x)
        sc = L.seifert_circles()
        regions = L.regions()
        orientation = L.orientation()
        pd_code = x
        final = [sc, regions, pd_code, orientation]
        return final

    #**************************** PART - 2 ***********************************
    def _braidword_detection_(self):
        r"""
        Return the braidword of the input. We match the outgoing components to the
        incoming components, in doing so we order the crossings and see to which strand
        they belong, thereby developing the braidword

        OUTPUT:
            - Braidword representation of the link.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L._braidword_detection_()
            [1, -2, 1, -2]
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L._braidword_detection_()
            [1, -2, -2, 3, 2, -2, -2, -1, -2, -3, 2]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._braidword_detection_()
            [-1, 2, -1, -2, -2, 1, 1, -2]
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L._braidword_detection_()
            [-1, -2, -2, 1, 2, 2]
            sage: B = BraidGroup(8)
            sage: L = Link(B([1,1]))
            sage: L._braidword_detection_()
            [1, 1]
            sage: L = Link(B([1, 2, 1, -2, -1]))
            sage: L._braidword_detection_()
            [1, 2, -1, -2, 2]
        """
        # all the data from the previous method
        sc = self._info_all_moves_()[0]
        regions = self._info_all_moves_()[1]
        pd_code = self._info_all_moves_()[2]
        orient = self._info_all_moves_()[3]
        # making the regions positive
        regions_pos = deepcopy(regions)
        for i in range(len(regions_pos)):
            for j in range(len(regions_pos[i])):
                if regions_pos[i][j] < 0:
                    regions_pos[i][j] = (-1) * regions_pos[i][j]
        # finding which sc are same as regions
        # r[0] is the first seifert cirlce and r[1] is the last one
        # which coincides with a region and there are exactly two seifert
        # circles here.
        r = []
        for i in sc:
            for j in regions_pos:
                if set(i) == set(j):
                    r.append(i)
        # here we find the ordering of the seifert circles
        pd_copy = deepcopy(pd_code)
        seifert_order = []
        seifert_order.append(r[0])
        a = seifert_order[0]
        while True:
            tmp = []
            if a == r[1]:
                break
            for j in pd_copy:
                if len(list(set(a).intersection(set(j)))) != 0:
                    tmp.append(j)
            b = list(set(tmp[0]) - set(a))
            pd_copy = [x for x in pd_copy if x not in tmp]
            for k in sc:
                if len(list(set(b).intersection(set(k)))) != 0:
                    a = sc[sc.index(k)]
                    break
            seifert_order.append(k)
        # here we calculate the entering and leaving of each of the crossing
        entering = []
        leaving = []
        for i in range(len(pd_code)):
            t = []
            q = []
            if orient[i] == -1:
                t.append(pd_code[i][0])
                t.append(pd_code[i][3])
                q.append(pd_code[i][1])
                q.append(pd_code[i][2])
            elif orient[i] == 1:
                t.append(pd_code[i][0])
                t.append(pd_code[i][1])
                q.append(pd_code[i][2])
                q.append(pd_code[i][3])
            entering.append(t)
            leaving.append(q)
        # here we correct the leaving and entering components so they belong to the
        # correct seifert circles
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
        first_seifert = seifert_order[0]
        first_crossing = None
        for i in pd_code:
            if len(list(set(i).intersection(set(first_seifert)))) == 2:
                first_crossing = i
                break
        tmp = []
        if orient[pd_code.index(first_crossing)] == -1:
            tmp.append(first_crossing[1])
            tmp.append(first_crossing[2])
        elif orient[pd_code.index(first_crossing)] == 1:
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
                for i, j in enumerate(entering):
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
        # record the signs
        sign = orient
        # each crossing belongs to two seifert circles, we find the first and
        # break
        braid = []
        for i in crossing:
            for j in seifert_order:
                if len(list(set(i).intersection(set(j)))) == 2:
                    braid.append(
                        sign[pd_code.index(i)] * (seifert_order.index(j) + 1))
                    break
        return braid

    def writhe(self):
        r"""
        Return the writhe of the knot.

        OUTPUT:
            - Writhe of the knot.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.writhe()
            0
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.writhe()
            -3
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.writhe()
            -2
        """
        x = self.oriented_gauss_code()
        pos = x[1].count(1)
        neg = (-1) * x[1].count(-1)
        return pos + neg

    # here we have used the symbolic ring rather than the Laurent Polynomial Ring.
    # The answer has been returned in the symbolic ring. Once the rational powers is
    # functional we can use it and revert back to Laurent Polynomial Ring.
    def jones_polynomial(self):
        r"""
        Return the jones polynomial of the link.
        The following procedure is used to determine the jones polynomial.
        There are two constructions :
        1. Smoothing of crossings in two ways, corresponding coefficients being
           A^-1 and A. The various permutations of the smoothing of the crossings
           is constructed by taking the permutation of the crossings.
        2. After this, the number of circles are calculated after smoothing and the
           polynomial is constructed by summing up in the following way
                polynomial = sigma[(no. of A^-1)*(no. of A)*d**((no. of circles) - 1)]

        OUTPUT:
            - Jones Polynomial of the link.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.jones_polynomial()
            q^2 - q - 1/q + 1/q^2 + 1
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.jones_polynomial()
            1/q + 1/q^3 - 1/q^4
            sage: l1 = [[1,4,2,3],[4,1,3,2]]
            sage: L = Link(l1)
            sage: L.jones_polynomial()
            -1/sqrt(q) - 1/q^(5/2)
            sage: l5 = [[1,8,2,7],[8,4,9,5],[3,9,4,10],[10,1,7,6],[5,3,6,2]]
            sage: L = Link(l5)
            sage: L.jones_polynomial()
            -q^(3/2) + sqrt(q) - 2/sqrt(q) + 1/q^(3/2) - 2/q^(5/2) + 1/q^(7/2)
        """
        pd = self.PD_code()
        x = SR.symbol('q')
        # here we look at the smoothings, either x or x**-1
        label_1 = [x for i in range(len(pd))]
        label_2 = [x ** -1 for i in range(len(pd))]
        label_1.extend(label_2)
        P = Permutations(label_1, len(pd))
        crossing_to_label = []
        # we record how each crossing is smoothened
        for i in P.list():
            tmp = {}
            for j, k in enumerate(i):
                tmp.update({tuple(pd[j]): k})
            crossing_to_label.append(tmp)
        # the product of the coefficients is calcaluated
        product = []
        for i in crossing_to_label:
            p = 1
            for j in i.itervalues():
                p = p * j
            product.append(p)
        # here we calculate the number of circles after the smoothing has been performed at
        # every crossing.
        pd_edit = []
        for i in crossing_to_label:
            pd_copy = deepcopy(pd)
            tmp = []
            for j in pd_copy:
                if i[tuple(j)] == x ** -1:
                    tmp.append(
                        [pd_copy[pd_copy.index(j)][0], pd_copy[pd_copy.index(j)][1]])
                    tmp.append(
                        [pd_copy[pd_copy.index(j)][3], pd_copy[pd_copy.index(j)][2]])
                elif i[tuple(j)] == x:
                    tmp.append(
                        [pd_copy[pd_copy.index(j)][1], pd_copy[pd_copy.index(j)][2]])
                    tmp.append(
                        [pd_copy[pd_copy.index(j)][0], pd_copy[pd_copy.index(j)][3]])
            pd_edit.append(tmp)
        # replacing the old edges with the new ones to get the circles and also
        # the number of circles.
        pd_edit = _rule_3_(pd_edit)
        pd_edit = _rule_3_(pd_edit)
        for i in pd_edit:
            for j in range(len(i)):
                for k in reversed(range(j + 1, len(i))):
                    if i[j] == i[k]:
                        del i[k]
        # the number of circles.
        circle_count = [len(i) for i in pd_edit]
        # we calculate the terms of the polynomial
        terms = [i * ((-x ** 2 - (x ** (-1)) ** 2) ** (j - 1))
                 for i, j in zip(product, circle_count)]
        # add the terms to generate the polynomial
        poly = sum(terms).expand()
        wri = self.writhe()
        f = (((-x ** (3)) ** (-wri)) * poly).expand()
        return f.subs({x: x ** (ZZ(-1) / ZZ(4))})

#********************** Auxillary methods used ********************************
# rule_1 and rule_2 are used in the orientation method looks for entering,
# leaving pairs and fill the gaps where ever necessary.
# rule_3 is used in the jones_polynomial and replace the higer numbered
# components with the lower numbered ones in order to get the number of
# circles.


def _rule_1_(over):
    for i in range(0, len(over), 2):
        if over[i][1] == None:
            if over[i + 1][1] == 1:
                over[i][1] = -1
            elif over[i + 1][1] == -1:
                over[i][1] = 1
        elif over[i + 1][1] == None:
            if over[i][1] == 1:
                over[i + 1][1] = -1
            elif over[i][1] == -1:
                over[i + 1][1] = 1
    return over


def _rule_2_(over):
    for i in over:
        for j in over:
            if i[0] == j[0] and j[1] == None:
                if i[1] == 1:
                    j[1] = -1
                    break
                elif i[1] == -1:
                    j[1] = 1
                    break
    return over


def _rule_3_(pd):
    for i in pd:
        for j, k in enumerate(i):
            a = min(k)
            b = max(k)
            for l in i[0:]:
                if b in l:
                    i[i.index(l)][l.index(b)] = min(k)
    return pd


def _pd_check_(pd):
    pd_error = False
    for i in pd:
        if len(i) != 4:
            pd_error = True
    else:
        flat = [x for y in pd for x in y]
        set_flat = set(flat)
        for i in set_flat:
            if flat.count(i) != 2:
                pd_error = True
    return pd_error
