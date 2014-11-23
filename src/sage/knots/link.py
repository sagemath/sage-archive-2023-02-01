r"""
Link Class

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
"""

##############################################################################
#       Copyright (C) 2014
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.finite_rings.integer_mod import Mod
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from copy import deepcopy, copy
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.integer_ring import IntegerRing
from sage.combinat.permutation import Permutations
from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
from sage.symbolic.ring import SR, var
from sage.rings.integer import Integer
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.functions.generalized import sign
from sage.plot.line import line
from sage.plot.bezier_path import bezier_path
from sage.misc.flatten import flatten


class Link:
    """
    A Link can be created by using one of the conventions mentioned below:

    Braid:

    Generators of the braid group are used to generate the link::

        sage: B = BraidGroup(8)
        sage: L = Link(B([-1, -1, -1, -2,1, -2,3,-2,3]))
        sage: L
        Knot represented by 9 crossings
        sage: L = Link(B([1, 2,1, -2,-1]))
        sage: L
        Link with 2 components represented by 5 crossings

    Oriented Gauss Code:

    Randomly number the crossings from 1 to n (where n is the number of
    crossings) and start moving along the link. Trace every component of
    the link, by starting at a particular point on one component of the link and
    taking note of each of the crossings until one returns to the starting
    point. Note each component as a list whose elements are the crossing
    numbers. Compile every component info into a list. We need the orientation
    of every crossing. This is recorded as a list with +1 and -1, +1 is recorded
    if the direction from leaving over-cross to the leaving under-cross is
    anti-clockwise, -1 if the direction from the leaving over-cross to the
    leaving under-cross is clockwise::

        sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1,-1,-1,-1,+1,+1,-1,+1]])
        sage: L
        Knot represented by 8 crossings

    For links there is more than one component and the input is as follows::

        sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
        sage: L
        Link with 3 components represented by 4 crossings

    Planar Diagram Code:

    Select some point on the link. Start numbering the strands in the
    components of the link. For a new component add one to the greatest
    number from the previous component and proceed till all the strands
    are numbered. At every cross contruct the data as follows :
    Start with the strand number of the entering under-cross and move in the
    clockwise direction around the cross and note down the strand numbers.
    Construct this data at every crossing and that would give the PD-Code.

    There is no particular distinction between knots and links for this input

    One of the representations of the Trefoil knot::

        sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sage: L
        Knot represented by 3 crossings

    One of the representations of the Hopf link::

        sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sage: L
        Link with 2 components represented by 2 crossings
    """
    def __init__(self, data):
        """

        The Python constructor.

        A Link can be created by using one of the conventions mentioned below:
        1. Braid
        2. Oriented Gauss Code
        3. Planar Diagram Code
        """
        if isinstance(data, list):
            if not data:
                raise ValueError("Does not accept empty list as arguement")

            if len(data) != 2 or not all(isinstance(i, list) for i in data[0]):
                for i in data:
                    if len(i) != 4:
                        raise ValueError("Not a valid PD code: crossings must be represented by four segments")
                    else:
                        flat = flatten(data)
                        set_flat = set(flat)
                        for i in set_flat:
                            if flat.count(i) != 2:
                                raise ValueError("Not a valid PD code: each segment must appear twice")
                self._PD_code = data
                self._oriented_gauss_code = None
                self._braid = None

            else:
                flat = flatten(data[0])
                a, b = max(flat), min(flat)
                if 2 * len(data[1]) != len(flat) or set(range(b, a + 1)) - set([0]) != set(flat):
                    raise ValueError("Invalid Input: Data is not a valid Oriented Gauss Code")
                self._oriented_gauss_code = data
                self._PD_code = None
                self._braid = None

        else:
            from sage.groups.braid import Braid
            if isinstance(data, Braid):
                self._braid = data
                self._oriented_gauss_code = None
                self._PD_code = None

            else:
                raise Exception("Invalid Input: Data must be either a list or a Braid")

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT:

        - A string.

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

    def braid(self):
        """
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
        from sage.groups.braid import BraidGroup
        if self._braid != None:
            return self._braid
        # look for possible vogel moves, perform them and call recursively to the modified link
        for region in self.regions():
            n = len(region)
            for i in range(n-1):
                a = region[i]
                seifcirca = filter(lambda x: abs(a) in x, self.seifert_circles())
                for j in range(i+1,n):
                    b = region[j]
                    seifcircb = filter(lambda x: abs(b) in x, self.seifert_circles())
                    if seifcirca != seifcircb and sign(a) == sign(b):
                        tails, heads = self._directions_of_edges_()
                        newedge = max(flatten(self.PD_code())) + 1
                        newPD = deepcopy(self.PD_code())
                        if sign(a) == 1:
                            C1 = newPD[newPD.index(heads[a])]
                            C1[C1.index(a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[b])]
                            C2[C2.index(b)] = newedge + 2
                            newPD.append([newedge + 3, a, b, newedge])
                            newPD.append([newedge + 2, newedge + 1, newedge + 3, newedge])
                            self._braid =  Link(newPD).braid()
                            return self._braid
                        else:
                            C1 = newPD[newPD.index(heads[-a])]
                            C1[C1.index(-a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[-b])]
                            C2[C2.index(-b)] = newedge + 2
                            newPD.append([newedge + 2, newedge, newedge + 3, newedge + 1])
                            newPD.append([newedge + 3, newedge, -b , -a])
                            self._braid = Link(newPD).braid()
                            return self._braid
        # We are in the case where no Vogel moves are necessary.
        G = DiGraph()
        G.add_vertices(map(tuple, self.seifert_circles()))
        for i in range(len(self.PD_code())):
            c = self.PD_code()[i]
            if self.orientation()[i] == 1:
                a  = filter(lambda x: c[1] in x, self.seifert_circles())[0]
                b  = filter(lambda x: c[0] in x, self.seifert_circles())[0]
                G.add_edge(tuple(a), tuple(b))
            else:
                a  = filter(lambda x: c[0] in x, self.seifert_circles())[0]
                b  = filter(lambda x: c[3] in x, self.seifert_circles())[0]
                G.add_edge(tuple(a), tuple(b))
        ordered_cycles = G.all_simple_paths(starting_vertices=G.sources(), ending_vertices=G.sinks())[0]
        B = BraidGroup(len(ordered_cycles))
        available_crossings = copy(self.PD_code())
        crossing = filter(lambda x: set(ordered_cycles[0]).intersection(set(x)), self.PD_code())[0]
        available_crossings.remove(crossing)
        status = [None for i in ordered_cycles]
        if self.orientation()[self.PD_code().index(crossing)] == 1:
            b = B([1])
            status[0] = crossing[2]
            status[1] = crossing[3]
        else:
            b = B([-1])
            status[0] = crossing[1]
            status[1] = crossing[2]
        counter = 0
        while available_crossings:
            possibles = filter(lambda x: status[counter] in x, available_crossings)
            if len(status) < counter + 2 or status[counter + 1] != None:
                possibles = filter(lambda x: status[counter + 1] in x, possibles)
            if possibles:
                added = possibles[0]
                if self.orientation()[self.PD_code().index(added)] == 1:
                    b *= B([counter + 1])
                    status[counter] = added[2]
                    status[counter + 1] = added[3]
                else:
                    b *= B([-counter - 1])
                    status[counter] = added[1]
                    status[counter + 1] = added[2]
                if counter > 0:
                    counter -= 1
                available_crossings.remove(added)
            else:
                counter += 1
        self._braid = b
        return b

    def _directions_of_edges_(self):
        r"""
        Return a tuple of two dictionaries. The first one assigns
        each edge of the PD code to the crossing where it starts.
        The second dictionary assigns it to where it ends.

        EXAMPLES::

            sage: L = Link([[1, 3, 2, 4], [2, 3, 1, 4]])
            sage: L._directions_of_edges_()
            ({1: [2, 3, 1, 4], 2: [1, 3, 2, 4], 3: [1, 3, 2, 4], 4: [2, 3, 1, 4]}, {1: [1, 3, 2, 4], 2: [2, 3, 1, 4], 3: [2, 3, 1, 4], 4: [1, 3, 2, 4]})

        ::

            sage: L = Link([[1,5,2,4],[5,3,6,2],[3,1,4,6]])
            sage: L._directions_of_edges_()
            ({1: [3, 1, 4, 6],
            2: [1, 5, 2, 4],
            3: [5, 3, 6, 2],
            4: [3, 1, 4, 6],
            5: [1, 5, 2, 4],
            6: [5, 3, 6, 2]},
            {1: [1, 5, 2, 4],
            2: [5, 3, 6, 2],
            3: [3, 1, 4, 6],
            4: [1, 5, 2, 4],
            5: [5, 3, 6, 2],
            6: [3, 1, 4, 6]})
        """
        tails = {}
        heads = {}
        for C in self.PD_code():
            tails[C[2]] = C
            a = C[2]
            D = C
            while not a in heads:
                next_crossing = filter(lambda x: a in x and x != D, self.PD_code())
                if len(next_crossing) == 0:
                    heads[a] = D
                    tails[a] = D
                    if D[0] == a:
                        a = D[2]
                    elif D[1] == a:
                        a = D[3]
                    else:
                        a = D[1]
                else:
                    heads[a] = next_crossing[0]
                    tails[a] = D
                    D = next_crossing[0]
                    a = D[(D.index(a)+2) % 4]
        unassigned = set(flatten(self.PD_code())).difference(set(tails.keys()))
        while unassigned:
            a = unassigned.pop()
            D = filter(lambda x: a in x, self.PD_code())[0]
            while not a in heads:
                tails[a] = D
                next_crossing = filter(lambda x: a in x and x != D, self.PD_code())[0]
                heads[a] = next_crossing
                D = next_crossing
                a = D[(D.index(a)+2) % 4]
                if a in unassigned:
                    unassigned.remove(a)
        return tails, heads


    def oriented_gauss_code(self):
        """
        Return the oriented gauss code of the link. The oriented gauss
        code has two parts :

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
        """
        Return the Planar Diagram code of the link. The Planar Diagram is returned
        in the following format.

        We construct the crossing by starting with the entering component of the
        undercrossing, move in the clockwise direction and then generate the list.
        Suppose if the crossing is given by [a, b, c, d], then we interpret this
        information as :

        1. a is the entering component of the undercrossing
        2. b, d are the components of the overcrossing
        3. c is the leaving component of the undercrossing

        Convention :
        Orientation of the crossing :

        1. Leaving over crossing to leaving under crossing in clockwise direction denoted by -1
        2. Leaving over crossing to leaving under crossing in anticlockwise direction denoted by 1

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
                d = flatten(oriented_gauss_code[0])
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
                    pd.append(
                        [strings[i], strings[i - 1], strings_max + 1, strings_max + 2])
                else:
                    pd.append(
                        [strings[abs(i) - 1], strings_max + 1, strings_max + 2, strings[abs(i)]])
                strings[abs(i) - 1] = strings_max + 1
                strings[abs(i)] = strings_max + 2
                strings_max = strings_max + 2
            for i in pd:
                for j in range(4):
                    if i[j] in strings:
                        i[j] = strings.index(i[j]) + 1
            self._PD_code = pd
            return pd

    def gauss_code(self):
        """
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
        """
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
        b = self.braid().Tietze()
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

    def _dowker_notation_(self):
        """
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
            [(2, 1), (3, 6), (7, 5), (8, 10)]
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
        """
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
        """
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
        return flatten(bc)

    def _homology_generators_(self):
        """
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

    def seifert_matrix(self):
        """
        Return the Seifert Matrix associated with the braidword.

        OUTPUT:

        - The intersection matrix of a (not necessarily minimal) Seifert surface of the
          link.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.seifert_matrix()
            [ 0  0]
            [ 0 -1]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L.seifert_matrix()
            [ 0  0  0]
            [ 1 -1  0]
            [ 0  1 -1]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.seifert_matrix()
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
        """
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
        G = Graph()
        G.add_vertices(set(flatten(self.PD_code())))
        for c in self.PD_code():
            G.add_edge(c[0], c[2])
            G.add_edge(c[1], c[3])
        return G.connected_components_number()

    def is_knot(self):
        """
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
        """
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
        """
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
        m = 2 * (self.seifert_matrix() + self.seifert_matrix().transpose())
        e = m.eigenvalues()
        sum = 0
        s = []
        for i, j in enumerate(e):
            s.append(cmp(j, 0))
            sum = sum + s[i]
        return sum

    def alexander_polynomial(self, var='t'):
        """
        Return the alexander polynomial of the link

        INPUT:

        - ``var`` -- (default: ``'t'``); the variable in the polynomial.

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
        f = (self.seifert_matrix() - t *
             (self.seifert_matrix().transpose())).determinant()
        return t ** ((-max(f.exponents()) - min(f.exponents())) / 2) * f if f != 0 else f

    def determinant(self):
        """
        Return the determinant of the knot

        OUTPUT:

        - Determinant of the Knot

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.determinant()
            3
            sage: L = Link(B([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.determinant()
            65
        """
        if self.is_knot() == True:
            a = self.alexander_polynomial()
            return Integer(abs(a(-1)))
        else:
            raise Exception("Determinant implmented only for knots")

    def arf_invariant(self):
        """
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
        """
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


    def orientation(self):
        """
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
        directions = self._directions_of_edges_()[0]
        orientation = []
        for C in self.PD_code():
            if directions[C[1]] == C:
                orientation.append(-1)
            else:
                orientation.append(1)
        return orientation

    def seifert_circles(self):
        """
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
            [[1, 7, 5, 3], [2, 6], [4, 8]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[1, 13, 9, 3, 15, 5, 11, 7], [2, 10, 6, 12], [4, 16, 8, 14]]
            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1,-1,-1,-1,1,-1,1]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.seifert_circles()
            [[1, 7, 3, 11, 5], [2, 8, 14, 6], [4, 12, 10], [9, 13]]
            sage: L = Link([[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 11, 8], [2, 7, 12, 4, 5, 10], [3, 9, 6]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([1, 1, 1]))
            sage: L.seifert_circles()
            [[1, 3, 5], [2, 4, 6]]
        """
        available_segments = set(flatten(self.PD_code()))
        result = []
        tails, heads = self._directions_of_edges_()
        while available_segments:
            a = available_segments.pop()
            if heads[a] == tails[a]:
                result.append([a])
            else:
                C = heads[a]
                par = []
                while not a in par:
                    par.append(a)
                    if tails[C[(C.index(a) + 1) % 4]] == C:
                        a = C[(C.index(a) + 1) % 4]
                    else:
                        a = C[(C.index(a) - 1) % 4]
                    if a in available_segments:
                        available_segments.remove(a)
                    C = heads[a]
                result.append(par)
        return result

    def regions(self):
        """
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

    def writhe(self):
        """
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

    def jones_polynomial(self, var='q'):
        """
        Return the jones polynomial of the link.
        The following procedure is used to determine the jones polynomial.
        There are two constructions :

        1. Smoothing of crossings in two ways, corresponding coefficients being
           A^-1 and A. The various permutations of the smoothing of the crossings
           is constructed by taking the permutation of the crossings.

        2. After this, the number of circles are calculated after smoothing and the
           polynomial is constructed by summing up in the following way:

           polynomial = sigma[(no. of A^-1)*(no. of A)*d**((no. of circles) - 1)]

        Note : here we have used the symbolic ring rather than the Laurent Polynomial Ring.
        The answer has been returned in the symbolic ring. Once the rational powers is
        functional we can use it and revert back to Laurent Polynomial Ring.
        INPUT:

        The name of the variable

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
        t = SR(var)
        poly = self._bracket_(t)
        writhe = self.writhe()
        jones = (poly * (-t) ** (-3 * writhe)).expand()
        return jones.subs({t: t ** (ZZ(-1) / ZZ(4))})

    def _bracket_(self, variable='q'):
        r"""
        Return the Kaufmann bracket polynomial of the diagram.

        INPUT:

        - var (Default = q): the name of the variable of the result.

        Note that this is not an invariant of the link, but of the diagram.
        In particular, it is not invariant under Reidemeister I moves.
        """
        t = SR(variable)
        if len(self.PD_code()) == 1:
            if self.PD_code()[0][0] == self.PD_code()[0][1]:
                return -t**(-3)
            else:
                return -t**3
        cross = self.PD_code()[0]
        rest = deepcopy(self.PD_code()[1:])
        [a, b, c, d] = cross
        if a == b and c == d and len(rest) > 0:
            return ((t ** (-1) + t ** (-5)) * Link(rest)._bracket_(t)).expand()
        elif a == d and c == b and len(rest) > 0:
            return ((t ** (1) + t ** (5)) * Link(rest)._bracket_(t)).expand()
        elif a == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = c
            return ((-t ** (-3)) * Link(rest)._bracket_(t)).expand()
        elif a == d:
            for cross in rest:
                if c in cross:
                    cross[cross.index(c)] = b
            return ((-t ** 3) * Link(rest)._bracket_(t)).expand()
        elif c == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
            return ((-t ** 3) * Link(rest)._bracket_(t)).expand()
        elif c == d:
            for cross in rest:
                if b in cross:
                    cross[cross.index(b)] = a
            return ((-t ** (-3)) * Link(rest)._bracket_(t)).expand()
        else:
            rest_2 = deepcopy(rest)
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
                if c in cross:
                    cross[cross.index(c)] = b
            for cross in rest_2:
                if d in cross:
                    cross[cross.index(d)] = c
                if b in cross:
                    cross[cross.index(b)] = a
            return (t * Link(rest)._bracket_(t) + t ** (-1) * Link(rest_2)._bracket_(t)).expand()

    def _isolated_components_(self):
        r"""
        Return the PD codes of the isolated components of the link.

        Isolated components are links corresponding to subdiagrams that don't
        have any common crossing.

        EXAMPLES::

            sage: L = Link([[1, 1, 2, 2], [3, 3, 4, 4]])
            sage: L._isolated_components_()
            [[[1, 1, 2, 2]], [[3, 3, 4, 4]]]

        """
        G = Graph()
        for c in self.PD_code():
            G.add_vertex(tuple(c))
        for i in range(G.num_verts()-1):
            for j in range(i, G.num_verts()):
                if len(set(G.vertices()[i]).intersection(G.vertices()[j])) > 0:
                    G.add_edge(G.vertices()[i], G.vertices()[j])
        return [[list(i) for i in j] for j in G.connected_components()]

    def plot(self, **kwargs):
        r"""
        Plot the knot or link.
        """
        comp = self._isolated_components_()
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            P1 = L1.plot(**kwargs)
            P2 = L2.plot(**kwargs)
            xtra = P1.get_minmax_data()['xmax'] + P2.get_minmax_data()['xmin'] + 2
            for P in P2:
                if hasattr(P, 'path'):
                    for p in P.path[0]:
                        p[0] += xtra
                    for p in P.vertices:
                        p[0] += xtra
                else:
                    for p in P.xdata:
                        p += xtra
            return P1 + P2
        if not 'color' in kwargs:
            kwargs['color'] = 'blue'
        if not 'axes' in kwargs:
            kwargs['axes'] = False
        if not 'aspect_ratio' in kwargs:
            kwargs['aspect_ratio'] = 1
        # The idea is the same followed in linkplot, but using MLP instead of
        # network flows.
        # We start by computing a way to bend the edges left or right
        # such that the resulting regions are in fact closed regions
        # with straight angles, and using the minimal number of bends.
        regions = sorted(self.regions(), lambda a,b: len(a) < len(b))
        regions = regions[:-1]
        edges = list(set(flatten(self.PD_code())))
        edges.sort()
        MLP = MixedIntegerLinearProgram(maximization = True)
        # v will be the list of variables in the MLP problem. There will be
        # two variables for each edge: number of right bendings and number of
        # left bendings (at the end, since we are minimizing the total, only one
        # of each will be nonzero
        v = MLP.new_variable(nonnegative=True)
        for i in range(2*len(edges)):
            MLP.set_min(v[i], 0)
        # one condition for each region
        for i in range(len(regions)):
            cond = 0
            r = regions[i]
            es = 4 - len(r)
            for e in r:
                if e > 0:
                    cond = cond + v[2*edges.index(e)] - v[2*edges.index(e) + 1]
                else:
                    cond = cond - v[2*edges.index(-e)] + v[2*edges.index(-e) + 1]
            MLP.add_constraint(cond, min=es, max=es)
        MLP.set_objective(-sum(v.values()))
        MLP.solve()
        # we store the result in a vector s packing right bends as negative left ones
        s = range(len(edges))
        values = MLP.get_values(v)
        for i in range(len(edges)):
            s[i] = int(values[2*i] - values[2*i + 1])
        # segments represents the different parts of the previos edges after bending
        segments = {e:[(e,i) for i in range(abs(s[edges.index(e)])+1)] for e in edges}
        pieces = {tuple(i):[i] for j in segments.values() for i in j}
        nregions = []
        for r in regions:
            nregion = []
            for e in r:
                if e>0:
                    rev =  segments[e][:-1]
                    sig = sign(s[edges.index(e)])
                    nregion += [[a, sig] for a in rev]
                    nregion.append([segments[e][-1], 1])
                else:
                    rev = segments[-e][1:]
                    rev.reverse()
                    sig = sign(s[edges.index(-e)])
                    nregion+=[[a, -sig] for a in rev]
                    nregion.append([segments[-e][0], 1])
            nregions.append(nregion)
        N = max(segments.keys()) + 1
        segments = [i for j in segments.values() for i in j]
        badregions = filter(lambda a: -1 in [x[1] for x in a], nregions)
        while len(badregions)>0:
            badregion = badregions[0]
            badturns = []
            a = 0
            while badregion[a][1] != -1:
                a+=1
            c = -1
            b = a
            while c != 2:
                if b == len(badregion)-1:
                    b = 0
                else:
                    b+=1
                c += badregion[b][1]
            otherregion = filter(lambda a: badregion[b][0] in [x[0] for x in a], nregions)
            if len(otherregion) == 1:
                otherregion = None
            elif otherregion[0] == badregion:
                otherregion = otherregion[1]
            else:
                otherregion = otherregion[0]
            N1 = N
            N = N + 2
            N2 = N1 + 1
            segments.append(N1)
            segments.append(N2)
            if type(badregion[b][0]) in (int, Integer):
                segmenttoadd = filter(lambda x:badregion[b][0] in pieces[x], pieces.keys())
                if len(segmenttoadd) > 0:
                    pieces[segmenttoadd[0]].append(N2)
            else:
                pieces[tuple(badregion[b][0])].append(N2)
            if a<b:
                r1 = badregion[:a]+[[badregion[a][0],0], [N1,1]]+badregion[b:]
                r2 = badregion[a+1:b] + [[N2,1],[N1,1]]
            else:
                r1 =  badregion[b:a]+[[badregion[a][0],0], [N1,1]]
                r2 = badregion[:b] + [[N2,1],[N1,1]] + badregion[a+1:]
            if otherregion:
                c = filter(lambda x:badregion[b][0] == x[0], otherregion)
                c = otherregion.index(c[0])
                otherregion.insert(c+1,[N2,otherregion[c][1]])
                otherregion[c][1]=0
            nregions.remove(badregion)
            nregions.append(r1)
            nregions.append(r2)
            badregions = filter(lambda a: -1 in [x[1] for x in a], nregions)
        MLP = MixedIntegerLinearProgram(maximization = True)
        variables = {}
        for e in segments:
            variables[e] = MLP.new_variable(nonnegative=True)
            MLP.set_min(variables[e][0], 1)
        for r in nregions:
            horp = []
            horm = []
            verp = []
            verm = []
            direction = 0
            for se in r:
                if direction % 4 == 0:
                    horp.append(variables[se[0]][0])
                elif direction == 1:
                    verp.append(variables[se[0]][0])
                elif direction == 2:
                    horm.append(variables[se[0]][0])
                elif direction == 3:
                    verm.append(variables[se[0]][0])
                if se[1] == 1:
                    direction += 1
            MLP.add_constraint(sum(horp)-sum(horm), min=0, max=0)
            MLP.add_constraint(sum(verp)-sum(verm), min=0, max=0)
        MLP.set_objective(-sum([x[0] for x in variables.values()]))
        solved = MLP.solve()
        lengths = {piece:sum([MLP.get_values(variables[a])[0] for a in pieces[piece]]) for piece in pieces}
        image = line([], **kwargs)
        crossings = {tuple(self.PD_code()[0]):(0,0,0)}
        availables = self.PD_code()[1:]
        used_edges = []
        horizontal_eq = 0
        vertical_eq = 0
        ims = line([], **kwargs)
        while len(used_edges) < len(edges):
            i = 0
            j = 0
            while crossings.keys()[i][j] in used_edges:
                if j < 3:
                    j += 1
                else:
                    j = 0
                    i+=1
            c = crossings.keys()[i]
            e = c[j]
            used_edges.append(e)
            direction = (crossings[c][2] - c.index(e)) % 4
            orien = self.orientation()[self.PD_code().index(list(c))]
            if s[edges.index(e)] < 0:
                turn = -1
            else:
                turn = 1
            lengthse = [lengths[(e,i)] for i in range(abs(s[edges.index(e)])+1)]
            if c.index(e) == 0 or (c.index(e) == 1 and orien == 1) or (c.index(e) == 3 and orien == -1):
                turn = -turn
                lengthse.reverse()
            if c.index(e) % 2 == 0:
                tailshort=True
            else:
                tailshort = False
            x0 = crossings[c][0]
            y0 = crossings[c][1]
            im = []
            for l1 in range(len(lengthse)):
                l = lengthse[l1]
                if direction == 0:
                    x1 = x0 + l
                    y1 = y0
                elif direction == 1:
                    x1 = x0
                    y1 = y0 + l
                elif direction == 2:
                    x1 = x0 - l
                    y1 = y0
                elif direction == 3:
                    x1 = x0
                    y1 = y0 -l
                im.append(([[x0,y0],[x1,y1]], l, direction))
                direction = (direction + turn) %4
                x0 = x1
                y0 = y1
            direction = (direction - turn) %4
            c2 = [ee for ee in availables if e in ee]
            if len(c2) == 1:
                availables.remove(c2[0])
                crossings[tuple(c2[0])] = (x1,y1,(direction+c2[0].index(e) + 2) % 4)
            c2 = [ee for ee in self.PD_code() if (e in ee and ee != list(c))]
            if c2 == []:
                headshort = not tailshort
            else:
                if c2[0].index(e) % 2 == 0:
                    headshort = True
                else:
                    headshort = False
            a = deepcopy(im[0][0])
            b = deepcopy(im[-1][0])
            if tailshort:
                im[0][0][0][0] += cmp(a[1][0],im[0][0][0][0])*0.1
                im[0][0][0][1] += cmp(a[1][1],im[0][0][0][1])*0.1
            if headshort:
                im[-1][0][1][0] -= cmp(b[1][0],im[-1][0][0][0])*0.1
                im[-1][0][1][1] -= cmp(b[1][1],im[-1][0][0][1])*0.1
            l = line([], **kwargs)
            c = 0
            p = im[0][0][0]
            if len(im) == 4 and max([x[1] for x in im]) == 1:
                l = bezier_path([[im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]], **kwargs)
                p =  im[-1][0][1]
            else:
                while c < len(im)-1:
                    if im[c][1] > 1:
                        (a, b) = im[c][0]
                        if b[0] > a[0]:
                            e = (b[0] - 1, b[1])
                        elif b[0] < a[0]:
                            e = (b[0] + 1, b[1])
                        elif b[1] > a[1]:
                            e = (b[0], b[1] - 1)
                        elif b[1] < a[1]:
                            e = (b[0] , b[1] + 1)
                        l += line((p, e), **kwargs)
                        p = e
                    if im[c+1][1] == 1 and c < len(im) - 2:
                        xr = round(im[c+2][0][1][0])
                        yr = round(im[c+2][0][1][1])
                        xp = xr - im[c+2][0][1][0]
                        yp = yr - im[c+2][0][1][1]
                        q = [p[0] + im[c+1][0][1][0] -   im[c+1][0][0][0] -xp , p[1] + im[c+1][0][1][1] - im[c+1][0][0][1] - yp]
                        l += bezier_path([[p, im[c+1][0][0], im[c+1][0][1], q]], **kwargs)
                        c += 2
                        p = q
                    else:
                        if im[c+1][1] == 1:
                            q = im[c+1][0][1]
                        else:
                            q = [im[c+1][0][0][0] + sign(im[c+1][0][1][0] -  im[c+1][0][0][0]), im[c+1][0][0][1] + sign(im[c+1][0][1][1] -  im[c+1][0][0][1])]
                        l += bezier_path([[p, im[c+1][0][0], q]], **kwargs)
                        p = q
                        c += 1
            l += line([p, im[-1][0][1]], **kwargs)
            image += l
            ims += sum([line(a[0], **kwargs) for a in im])
        return image




