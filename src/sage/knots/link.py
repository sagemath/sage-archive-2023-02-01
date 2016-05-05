r"""
Links

A knot is defined as embedding of the circle `\mathbb{S}^1` in the
3-dimensional sphere `\mathbb{S}^3`, considered up to ambient isotopy.
They represent the physical idea of a knotted rope, but with the
particularity that the rope is closed. That is, the ends of the
rope are joined.

A link is an embedding of one or more copies of `\mathbb{S}^1` in
`\mathbb{S}^3`, considered up to ambient isotopy. That is, a link
represents the idea of one or more tied ropes. Every knot is a link,
but not every link is a knot.

Generically, the projection of a link on `\RR^2` is a curve with
crossings. The crossings are represented to show which strand goes
over the other. This curve is called a planar diagram of the link.
If we remove the crossings, the resulting connected components are
segments. These segments are called the edges of the diagram.

REFERENCES:

- :wikipedia:`Knot_(mathematics)`

.. [Collins13] Julia Collins. *An algorithm for computing the Seifert
   matrix of a link from a braid representation*. (2013).
   http://www.maths.ed.ac.uk/~jcollins/SeifertMatrix/SeifertMatrix.pdf

.. [KnotAtlas] The Knot atlas. http://katlas.org/wiki/Main_Page

.. SEEALSO::

    There are also tables of link and knot invariants at
    http://www.indiana.edu/~knotinfo/
    and http://www.indiana.edu/~linkinfo/.

AUTHORS:

- Miguel Angel Marco Buzunariz
- Amit Jamadagni
"""

#*****************************************************************************
#       Copyright (C) 2014  Miguel Angel Marco Buzunariz
#                           Amit Jamadagni
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ
from sage.graphs.digraph import DiGraph
from sage.graphs.graph import Graph
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.symbolic.ring import SR
from sage.rings.integer import Integer
from sage.numerical.mip import MixedIntegerLinearProgram
from sage.functions.generalized import sign
from sage.misc.flatten import flatten
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method
from copy import deepcopy, copy


class Link(object):
    r"""
    A link.

    A link is an embedding of one or more copies of `\mathbb{S}^1` in
    `\mathbb{S}^3`, considered up to ambient isotopy. That is, a link
    represents the idea of one or more tied ropes. Every knot is a link,
    but not every link is a knot.

    A link can be created by using one of the conventions mentioned below:

    Braid:

    - The closure of a braid is a link::

        sage: B = BraidGroup(8)
        sage: L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 3]))
        sage: L
        Link with 1 component represented by 9 crossings
        sage: L = Link(B([1, 2, 1, -2, -1]))
        sage: L
        Link with 2 components represented by 5 crossings

      .. NOTE::

          The strands of the braid that have no crossings at all
          are removed.

    - Oriented Gauss Code:

      Label the crossings from `1` to `n` (where `n` is the number of
      crossings) and start moving along the link. Trace every component of
      the link, by starting at a particular point on one component of the
      link and writing down each of the crossings that you encounter until
      returning to the starting point. The crossings are written with sign
      depending on whether we cross them as over or undercrossing. Each
      component is then represented as a list whose elements are the
      crossing numbers. A second list of `+1` and `-1`'s keeps track of
      the orientation of each crossing::

        sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
        ....:           [-1, -1, -1, -1, 1, 1, -1, 1]])
        sage: L
        Link with 1 component represented by 8 crossings

      For links there may be more than one component and the input is
      as follows::

        sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
        sage: L
        Link with 3 components represented by 4 crossings

    - Planar Diagram (PD) Code:

      The diagram of the link is formed by segments that are adjacent to
      the crossings. Label each one of this segments with a positive number,
      and for each crossing, write down the four incident segments. The
      order of these segments is clockwise, starting with the incoming
      undercrossing.

      There is no particular distinction between knots and links for
      this input.

    EXAMPLES:

    One of the representations of the trefoil knot::

        sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sage: L
        Link with 1 component represented by 3 crossings

    .. PLOT::
        :width: 300 px

        L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
        sphinx_plot(L.plot())

    One of the representations of the Hopf link::

        sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sage: L
        Link with 2 components represented by 2 crossings

    .. PLOT::
        :width: 300 px

        L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
        sphinx_plot(L.plot())

    We can construct links from from the braid group::

        sage: B = BraidGroup(4)
        sage: L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
        sage: L
        Link with 2 components represented by 8 crossings

    .. PLOT::
        :width: 300 px

        B = BraidGroup(4)
        L = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
        sphinx_plot(L.plot())

    ::

        sage: L = Link(B([1, 2, 1, 3]))
        sage: L
        Link with 2 components represented by 4 crossings

    .. PLOT::
        :width: 300 px

        B = BraidGroup(4)
        L = Link(B([1, 2, 1, 3]))
        sphinx_plot(L.plot())

    We construct the "monster" unknot using a planar code, and
    then construct the oriented Gauss code and braid representation::

        sage: L = Link([[3,1,2,4], [8,9,1,7], [5,6,7,3], [4,18,6,5],
        ....:           [17,19,8,18], [9,10,11,14], [10,12,13,11],
        ....:           [12,19,15,13], [20,16,14,15], [16,20,17,2]])
        sage: L.oriented_gauss_code()
        [[[1, -4, 3, -1, 10, -9, 6, -7, 8, 5, 4, -3, 2, -6, 7, -8, 9, -10, -5, -2]],
         [1, -1, 1, 1, 1, -1, -1, -1, -1, -1]]
        sage: L.braid()
        s0*s1^-1*s2^-1*s3^-1*s2*s1^-1*s0^-1*s1*s2^2*s1^-1*s3*s2*s1^-3

    .. PLOT::
        :width: 300 px

        L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
                  [17,19,8,18],[9,10,11,14],[10,12,13,11],
                  [12,19,15,13],[20,16,14,15],[16,20,17,2]])
        sphinx_plot(L.plot())

    We construct the Ochiai unknot by using an oriented Gauss code::

        sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
        ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
        ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
        sage: L.pd_code()
        [[10, 2, 11, 1], [2, 12, 3, 11], [3, 20, 4, 21], [12, 19, 13, 20],
         [21, 32, 22, 1], [31, 22, 32, 23], [9, 25, 10, 24], [4, 29, 5, 30],
         [23, 30, 24, 31], [28, 14, 29, 13], [17, 14, 18, 15], [5, 17, 6, 16],
         [15, 7, 16, 6], [7, 27, 8, 26], [25, 9, 26, 8], [18, 28, 19, 27]]

    .. PLOT::
        :width: 300 px

        L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
                    -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
                  [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
        sphinx_plot(L.plot())

    We construct the knot `7_1` and compute some invariants::

        sage: B = BraidGroup(2)
        sage: L = Link(B([1]*7))

    .. PLOT::
        :width: 300 px

        B = BraidGroup(2)
        L = Link(B([1]*7))
        sphinx_plot(L.plot())

    ::

        sage: L.alexander_polynomial()
        t^-3 - t^-2 + t^-1 - 1 + t - t^2 + t^3
        sage: L.jones_polynomial()
        -t^10 + t^9 - t^8 + t^7 - t^6 + t^5 + t^3
        sage: L.determinant()
        7
        sage: L.signature()
        -6

    The links here have removed components in which no strand is used::

        sage: B = BraidGroup(8)
        sage: b = B([1])
        sage: L = Link(b)
        sage: b.components_in_closure()
        7
        sage: L.number_of_components()
        1
        sage: L.braid().components_in_closure()
        1
        sage: L.braid().parent()
        Braid group on 2 strands

    .. WARNING::

        Equality of knots is done by comparing the corresponding braids,
        which may give false negatives.

    .. NOTE::

        The behavior of removing unused strands from an element of a
        braid group may change without notice in the future. Do not
        rely on this feature.

    .. TODO::

        Implement methods to creating new links from previously created links.
    """
    def __init__(self, data):
        """
        Initialize ``self``.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, -1, -1, -2,1, -2, 3, -2]))
            sage: TestSuite(L).run()
            sage: L = Link(B([1, 2, 1]))
            sage: TestSuite(L).run()
            sage: L = Link([[1, 1, 2, 2]])
            sage: TestSuite(L).run()

            sage: L = Link(B.one())
            sage: L = Link([])
            sage: L = Link([[], []])

            sage: Link([[[-1, 2, -1, 2]],  [1, 1, 1, 1]])
            Traceback (most recent call last):
            ...
            ValueError: invalid input: data is not a valid oriented Gauss code

            sage: Link([[[-1, 2, 3, 4]]])
            Traceback (most recent call last):
            ...
            ValueError: invalid PD code: crossings must be represented by four segments

            sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 3]])
            Traceback (most recent call last):
            ...
            ValueError: invalid PD code: each segment must appear twice

            sage: L = Link(5)
            Traceback (most recent call last):
            ...
            ValueError: invalid input: data must be either a list or a braid
        """
        if isinstance(data, list):
            if len(data) != 2 or not all(isinstance(i, list) for i in data[0]):
                for i in data:
                    if len(i) != 4:
                        raise ValueError("invalid PD code: crossings must be represented by four segments")
                    else:
                        flat = flatten(data)
                        if any(flat.count(i) != 2 for i in set(flat)):
                            raise ValueError("invalid PD code: each segment must appear twice")
                self._pd_code = data
                self._oriented_gauss_code = None
                self._braid = None

            else:
                flat = flatten(data[0])
                if flat:
                    a, b = max(flat), min(flat)
                    if 2 * len(data[1]) != len(flat) or set(range(b, a + 1)) - set([0]) != set(flat):
                        raise ValueError("invalid input: data is not a valid oriented Gauss code")
                self._oriented_gauss_code = data
                self._pd_code = None
                self._braid = None

        else:
            from sage.groups.braid import Braid, BraidGroup
            if isinstance(data, Braid):
                # Remove all unused strands
                support = sorted(set(abs(x) for x in data.Tietze()))
                i = 0
                cur = 1
                while i < len(support):
                    if support[i] == cur:
                        cur += 1
                        i += 1
                    elif support[i] == cur + 1:
                        support.insert(i, cur+1)
                        cur += 2
                        i += 2
                    else:
                        cur = support[i]
                        i += 1
                d = {}
                for i,s in enumerate(support):
                    d[s] = i+1
                    d[-s] = -i-1
                if not support:
                    B = BraidGroup(2)
                else:
                    B = BraidGroup(len(support)+1)
                self._braid = B([d[x] for x in data.Tietze()])
                self._oriented_gauss_code = None
                self._pd_code = None

            else:
                raise ValueError("invalid input: data must be either a list or a braid")

    def __repr__(self):
        """
        Return a string representation.

        OUTPUT: string representation

        EXAMPLES::

            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L
            Link with 1 component represented by 4 crossings
            sage: L = Link([[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]])
            sage: L
            Link with 3 components represented by 4 crossings
        """
        number_of_components = self.number_of_components()
        if number_of_components > 1:
            plural = 's'
        else:
            plural = ''
        pd_len = len(self.pd_code())
        return 'Link with {} component{} represented by {} crossings'.format(number_of_components, plural, pd_len)

    def __eq__(self, other):
        """
        Check equality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 == L2
            True
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 == L3
            False
        """
        if not isinstance(other, self.__class__):
            return False
        if self._pd_code is not None:
            if self.pd_code() == other.pd_code():
                return True
        if self._oriented_gauss_code is not None:
            if self.oriented_gauss_code() == other.oriented_gauss_code():
                return True
        return self.braid() == other.braid()

    def __ne__(self, other):
        """
        Check inequality.

        TESTS::

            sage: B = BraidGroup(8)
            sage: L1 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L2 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2, 5, 4]))
            sage: L1 != L2
            False
            sage: L3 = Link(B([-1, -1, -1, -2, 1, -2, 3, -2]))
            sage: L1 != L3
            True
        """
        return not self.__eq__(other)

    def braid(self):
        """
        Return a braid representation of ``self``.

        OUTPUT: an element in the braid group

        EXAMPLES::

            sage: L = Link([[2, 3, 1, 4], [4, 1, 3, 2]])
            sage: L.braid()
            s^2
            sage: L = Link([[[-1, 2, -3, 1, -2, 3]], [-1, -1, -1]])
            sage: L.braid()
            s^-3
            sage: L = Link([[1,8,2,7], [8,4,9,5], [3,9,4,10], [10,1,7,6], [5,3,6,2]])
            sage: L.braid()
            (s0*s1^-1)^2*s1^-1

        TESTS::

            sage: L = Link([])
            sage: L.braid()
            1
            sage: L = Link([[], []])
            sage: L.braid()
            1
        """
        if self._braid is not None:
            return self._braid

        from sage.groups.braid import BraidGroup
        comp = self._isolated_components()
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            b1 = L1.braid()
            b2 = L2.braid()
            n1 = b1.parent().strands()
            n2 = b2.parent().strands()
            t1 = list(b1.Tietze())
            t2 = [sign(x)*(abs(x) + n1) for x in b2.Tietze()]
            B = BraidGroup(n1 + n2)
            self._braid = B(t1 + t2)
            return self._braid

        # look for possible Vogel moves, perform them and call recursively to the modified link
        pd_code = self.pd_code()
        if not pd_code:
            B = BraidGroup(2)
            self._braid = B.one()
            return self._braid
        seifert_circles = self.seifert_circles()
        newedge = max(flatten(pd_code)) + 1
        for region in self.regions():
            n = len(region)
            for i in range(n-1):
                a = region[i]
                seifcirca = [x for x in seifert_circles if abs(a) in x]
                for j in range(i+1,n):
                    b = region[j]
                    seifcircb = [x for x in seifert_circles if abs(b) in x]
                    if seifcirca != seifcircb and sign(a) == sign(b):
                        tails, heads = self._directions_of_edges()

                        newPD = deepcopy(pd_code)
                        if sign(a) == 1:
                            C1 = newPD[newPD.index(heads[a])]
                            C1[C1.index(a)] = newedge + 1
                            C2 = newPD[newPD.index(tails[b])]
                            C2[C2.index(b)] = newedge + 2
                            newPD.append([newedge + 3, a, b, newedge])
                            newPD.append([newedge + 2, newedge + 1, newedge + 3, newedge])
                            self._braid = Link(newPD).braid()
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
        G.add_vertices([tuple(c) for c in seifert_circles])
        for i,c in enumerate(pd_code):
            if self.orientation()[i] == 1:
                a  = [x for x in seifert_circles if c[1] in x][0]
                b  = [x for x in seifert_circles if c[0] in x][0]
            else:
                a  = [x for x in seifert_circles if c[0] in x][0]
                b  = [x for x in seifert_circles if c[3] in x][0]
            G.add_edge(tuple(a), tuple(b))

        # Get a simple path from a source to a sink in the digraph
        it = G.all_paths_iterator(starting_vertices=G.sources(), ending_vertices=G.sinks(), simple=True)
        ordered_cycles = next(it)

        B = BraidGroup(len(ordered_cycles))
        available_crossings = copy(pd_code)
        oc_set = set(ordered_cycles[0])
        for i,x in enumerate(pd_code):
            if any(elt in oc_set for elt in x):
                crossing = x
                crossing_index = i
                break
        available_crossings.remove(crossing)
        status = [None for i in ordered_cycles]
        orientation = self.orientation()
        if orientation[crossing_index] == 1:
            b = B([1])
            status[0] = crossing[2]
            status[1] = crossing[3]
        else:
            b = B([-1])
            status[0] = crossing[1]
            status[1] = crossing[2]
        counter = 0
        while available_crossings:
            possibles = [x for x in available_crossings if status[counter] in x]
            if len(status) < counter + 2 or status[counter + 1] is not None:
                possibles = [x for x in possibles if status[counter + 1] in x]
            if possibles:
                added = possibles[0]
                if orientation[pd_code.index(added)] == 1:
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

    def _directions_of_edges(self):
        r"""
        Return the directions of the edges given by the PD code of ``self``.

        OUTPUT:

        A tuple of two dictionaries. The first one assigns
        each edge of the PD code to the crossing where it starts.
        The second dictionary assigns it to where it ends.

        EXAMPLES::

            sage: L = Link([[1, 3, 2, 4], [2, 3, 1, 4]])
            sage: L._directions_of_edges()
            ({1: [2, 3, 1, 4], 2: [1, 3, 2, 4], 3: [1, 3, 2, 4], 4: [2, 3, 1, 4]},
             {1: [1, 3, 2, 4], 2: [2, 3, 1, 4], 3: [2, 3, 1, 4], 4: [1, 3, 2, 4]})

        ::

            sage: L = Link([[1,5,2,4], [5,3,6,2], [3,1,4,6]])
            sage: L._directions_of_edges()
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

        ::

            sage: L = Link([[1,2,3,3], [2,4,5,5], [4,1,7,7]])
            sage: L._directions_of_edges()
            ({1: [4, 1, 7, 7],
              2: [1, 2, 3, 3],
              3: [1, 2, 3, 3],
              4: [2, 4, 5, 5],
              5: [2, 4, 5, 5],
              7: [4, 1, 7, 7]},
             {1: [1, 2, 3, 3],
              2: [2, 4, 5, 5],
              3: [1, 2, 3, 3],
              4: [4, 1, 7, 7],
              5: [2, 4, 5, 5],
              7: [4, 1, 7, 7]})
        """
        tails = {}
        heads = {}
        pd_code = self.pd_code()
        for C in pd_code:
            tails[C[2]] = C
            a = C[2]
            D = C
            while not a in heads:
                next_crossing = [x for x in pd_code if a in x and x != D]
                if not next_crossing:
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

        unassigned = set(flatten(pd_code)).difference(set(tails.keys()))
        while unassigned:
            a = unassigned.pop()
            for x in pd_code:
                if a in x:
                    D = x
                    break
            while not a in heads:
                tails[a] = D
                for x in pd_code:
                    if a in x and x != D:
                        next_crossing = x
                        break
                heads[a] = next_crossing
                D = next_crossing
                a = D[(D.index(a)+2) % 4]
                if a in unassigned:
                    unassigned.remove(a)
        return tails, heads

    def oriented_gauss_code(self):
        """
        Return the oriented Gauss code of ``self``.

        The oriented Gauss code has two parts:

        a. the Gauss code

        b. the orientation of each crossing

        The following orientation was taken into consideration for
        construction of knots:

        From the outgoing of the overcrossing if we move in the clockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `-1`.

        From the outgoing of the overcrossing if we move in the anticlockwise
        direction to reach the outgoing of the undercrossing then we label
        that crossing as `+1`.

        One more consideration we take in while constructing the orientation
        is the order of the orientation is same as the ordering of the
        crossings in the Gauss code.

        .. NOTE::

            Convention: under is denoted by `-1`, and over by `+1` in the
            crossing info.

        EXAMPLES::

            sage: L = Link([[1, 11, 2, 10], [6, 2, 7, 3], [3, 12, 4, 9], [9, 5, 10, 6], [8, 1, 5, 4], [11, 8, 12, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2, -3, 5], [4, -2, 6, -5], [-4, 1, -6, 3]], [-1, 1, 1, 1, -1, -1]]
            sage: L = Link([[1, 4, 2, 3], [6, 1, 3, 2], [7, 4, 8, 5], [5, 8, 6, 7]])
            sage: L.oriented_gauss_code()
            [[[-1, 2], [-3, 4], [1, 3, -4, -2]], [-1, -1, 1, 1]]
            sage: B = BraidGroup(8)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.oriented_gauss_code()
            [[[1, -2, 3, -4, 5, -1, 2, -3, 4, -5]], [1, 1, 1, 1, 1]]

        TESTS::

            sage: L = Link([])
            sage: L.oriented_gauss_code()
            [[], []]
            sage: L = Link(BraidGroup(2).one())
            sage: L.oriented_gauss_code()
            [[], []]
        """
        if self._oriented_gauss_code is not None:
            return self._oriented_gauss_code

        pd = self.pd_code()
        orient = self.orientation()
        crossing_info = {}
        for i, j in enumerate(pd):
            if orient[i] == -1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[3], 1, i + 1)] = j[1]
            elif orient[i] == 1:
                crossing_info[(j[0], -1, i + 1)] = j[2]
                crossing_info[(j[1], 1, i + 1)] = j[3]
        edges = {}
        cross_number = {}
        for i, j in crossing_info.items():
            edges[i[0]] = [j]
            if i[1] == 1:
                cross_number[i[0]] = i[2]
            elif i[1] == -1:
                cross_number[i[0]] = -i[2]
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

    def pd_code(self):
        """
        Return the planar diagram code of ``self``.

        The planar diagram is returned in the following format.

        We construct the crossing by starting with the entering component
        of the undercrossing, move in the clockwise direction and then
        generate the list. If the crossing is given by `[a, b, c, d]`,
        then we interpret this information as:

        1. `a` is the entering component of the undercrossing;
        2. `b, d` are the components of the overcrossing;
        3. `c` is the leaving component of the undercrossing.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]], [1, 1, -1, -1]])
            sage: L.pd_code()
            [[6, 1, 7, 2], [2, 5, 3, 6], [8, 4, 1, 3], [4, 8, 5, 7]]
            sage: B = BraidGroup(2)
            sage: b = B([1, 1, 1, 1, 1])
            sage: L = Link(b)
            sage: L.pd_code()
            [[2, 1, 3, 4], [4, 3, 5, 6], [6, 5, 7, 8], [8, 7, 9, 10], [10, 9, 1, 2]]
            sage: L = Link([[[2, -1], [1, -2]], [1, 1]])
            sage: L.pd_code()
            [[2, 3, 1, 4], [4, 1, 3, 2]]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.pd_code()
            [[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]]

        TESTS::

            sage: L = Link([[], []])
            sage: L.pd_code()
            []
            sage: L = Link(BraidGroup(2).one())
            sage: L.pd_code()
            []
        """
        if self._pd_code is not None:
            return self._pd_code

        if self._oriented_gauss_code is not None:
            oriented_gauss_code = self._oriented_gauss_code
            d_dic = {}
            if len(oriented_gauss_code[0]) > 1:
                d = flatten(oriented_gauss_code[0])
                for i, j in enumerate(d):
                    d_dic[j] = [i + 1, i + 2]
                # here we collect the final component in each Gauss code
                last_component = [i[-1] for i in oriented_gauss_code[0]]
                first_component = [i[0] for i in oriented_gauss_code[0]]
                # here we correct the last_component
                for i, j in zip(last_component, first_component):
                    d_dic[i][1] = d_dic[j][0]
                crossing_dic = {}
                for i,x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]
            elif len(oriented_gauss_code[0]) == 1:
                for i, j in enumerate(oriented_gauss_code[0][0]):
                    d_dic[j] = [i + 1, i + 2]
                d_dic[oriented_gauss_code[0][0][-1]][1] = 1
                crossing_dic = {}
                for i, x in enumerate(oriented_gauss_code[1]):
                    if x == -1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][1],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][0]]
                    elif x == 1:
                        crossing_dic[i + 1] = [d_dic[-(i + 1)][0], d_dic[i + 1][0],
                                               d_dic[-(i + 1)][1], d_dic[i + 1][1]]
            else:
                crossing_dic = {}

            pd = crossing_dic.values()
            self._pd_code = pd
            return self._pd_code

        if self._braid is not None:
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
            self._pd_code = pd
            return pd

        raise AssertionError("invalid state")

    def gauss_code(self):
        """
        Return the Gauss code of ``self``.

        The Gauss code is generated by the following procedure:

        a. Number the crossings from `1` to `n`.
        b. Select a point on the knot and start moving along the component.
        c. At each crossing, take the number of the crossing, along with
           sign, which is `-` if it is an undercrossing and `+` if it is a
           overcrossing.

        EXAMPLES::

            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
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

    def dowker_notation(self):
        """
        Return the Dowker notation of ``self``.

        Similar to the PD code we number the components, so every crossing
        is represented by four numbers. We focus on the incoming entities
        of the under and the overcrossing. It is the pair of incoming
        undercrossing and the incoming overcrossing. This information at
        every crossing gives the Dowker notation.

        OUTPUT:

        A list containing the pair of incoming under cross and the incoming
        over cross.

        EXAMPLES::

            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6,-5]], [-1, -1, -1, -1, 1, -1, 1]])
            sage: L.dowker_notation()
            [(1, 6), (7, 2), (3, 10), (11, 4), (14, 5), (13, 8), (12, 9)]
            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.dowker_notation()
            [(2, 1), (3, 5), (6, 4), (7, 9)]
            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sage: L.dowker_notation()
            [(1, 3), (4, 2)]
        """
        pd = self.pd_code()
        orient = self.orientation()
        dn = [(i[0], i[3]) if orient[j] == -1 else (i[0], i[1])
              for j, i in enumerate(pd)]
        return dn

    def _braid_word_components(self):
        """
        Return the disjoint braid components, if any, else return the braid
        of ``self``.

        For example consider the braid ``[-1, 3, 1, 3]`` this can be viewed
        as a braid with components as ``[-1, 1]`` and ``[3, 3]``. There is no
        common crossing to these two (in sense there is a crossing between
        strand `1` and `2`, crossing between `3` and `4` but no crossing
        between strand `2` and `3`, so these can be viewed as independent
        components in the braid).

        OUTPUT: list containing the components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braid_word_components()
            ([-1, 1], [3, 3])
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braid_word_components()
            ([-1, 1, 1, 1], [3], [5, 7, 6])
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braid_word_components()
            ([-2, 1, 1], [4, 4], [6])
        """
        ml = list(self.braid().Tietze())
        if not ml:
            return tuple()

        l = set(abs(k) for k in ml)
        missing1 = set(range(min(l), max(l) + 1)) - l
        if not missing1:
            return (ml,)

        missing = sorted(missing1)
        x = [[] for i in range(len(missing) + 1)]
        for i,a in enumerate(missing):
            for j, mlj in enumerate(ml):
                if mlj != 0 and abs(mlj) < a:
                    x[i].append(mlj)
                    ml[j] = 0
                elif mlj != 0 and abs(mlj) > missing[-1]:
                    x[-1].append(mlj)
                    ml[j] = 0
        return tuple([a for a in x if a])

    def _braid_word_components_vector(self):
        """
        The list from the :meth:`_braid_word_components` is flattened to
        give out the vector form.

        OUTPUT: list containing braid word components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._braid_word_components_vector()
            [-1, 1, 3, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._braid_word_components_vector()
            [-1, 1, 1, 1, 3, 5, 7, 6]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._braid_word_components_vector()
            [-2, 1, 1, 4, 4, 6]
        """
        return flatten(self._braid_word_components())

    def _homology_generators(self):
        """
        The set of generators for the first homology group of the connected
        Seifert surface of the given link.

        This method uses the :meth:`_braid_word_components_vector` to generate
        the homology generators. The position of the repeated element w.r.t.
        the braid word component vector list is compiled into a list.

        This is based on Lemma 3.1 in [Collins13]_.

        OUTPUT:

        A list of integers `i \in \{1, 2, \ldots, n-1\}` corresponding
        to the simple generators `s_i` that gives a homology generator or
        `0` if the position does not represent a generator.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L._homology_generators()
            [1, 0, 3]
            sage: B = BraidGroup(8)
            sage: L = Link(B([-1, 3, 1, 5, 1, 7, 1, 6]))
            sage: L._homology_generators()
            [1, 2, 3, 0, 0, 0, 0]
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L._homology_generators()
            [0, 2, 0, 4, 0]
        """
        x = self._braid_word_components_vector()
        hom_gen = []
        for j in range(len(x) - 1):
            a = abs(x[j])
            for i in range(j + 1, len(x)):
                if a == abs(x[i]):
                    hom_gen.append(i)
                    break
            else:
                hom_gen.append(0)
        return hom_gen

    @cached_method
    def seifert_matrix(self):
        """
        Return the Seifert matrix associated with ``self``.

        ALGORITHM:

        This is the algorithm presented in Section 3.3 of [Collins13]_.

        OUTPUT:

        The intersection matrix of a (not necessarily minimal) Seifert surface.

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
        x = self._braid_word_components_vector()
        h = self._homology_generators()
        hl = len(h)
        A = matrix(ZZ, hl, hl)
        indices = [i for i,hi in enumerate(h) if hi != 0]
        for i in indices:
            hi = h[i]
            for j in range(i, hl):
                if i == j:
                    A[i, j] = cmp(0, x[i] + x[hi])
                elif hi > h[j]:
                    A[i, j] = 0
                    A[j, i] = 0
                elif hi < j:
                    A[i, j] = 0
                    A[j, i] = 0
                elif hi == j:
                    if x[j] > 0:
                        A[i, j] = 0
                        A[j, i] = 1
                    else:
                        A[i, j] = -1
                        A[j, i] = 0
                elif abs(abs(x[i]) - abs(x[j])) > 1:
                    A[i, j] = 0
                elif abs(x[i]) - abs(x[j]) == 1:
                    A[i, j] = 0
                    A[j, i] = -1
                elif abs(x[j]) - abs(x[i]) == 1:
                    A[i, j] = 1
                    A[j, i] = 0
                else:  # for debugging
                    A[i, j] = 2
                    A[j, i] = 2
        A = A.matrix_from_rows_and_columns(indices, indices)
        A.set_immutable()
        return A

    @cached_method
    def number_of_components(self):
        """
        Return the number of connected components of ``self``.

        OUTPUT: number of connected components

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.number_of_components()
            4
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.number_of_components()
            5
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.number_of_components()
            1
            sage: L = Link(B.one())
            sage: L.number_of_components()
            1
        """
        G = Graph()
        pd = self.pd_code()
        if not pd:
            return ZZ.one()
        G.add_vertices(set(flatten(pd)))
        for c in pd:
            G.add_edge(c[0], c[2])
            G.add_edge(c[1], c[3])
        return G.connected_components_number()

    def is_knot(self):
        """
        Return ``True`` if ``self`` is a knot.

        Every knot is a link but the converse is not true.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1, 3, 1, -3]))
            sage: L.is_knot()
            False
            sage: B = BraidGroup(8)
            sage: L = Link(B([1, 2, 3, 4, 5, 6]))
            sage: L.is_knot()
            True
        """
        return self.number_of_components() == 1

    def genus(self):
        """
        Return the genus of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.genus()
            0
            sage: L = Link(B([1,3]))
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
        if not b:
            return ZZ.zero()

        B = self.braid().parent()
        x = self._braid_word_components()
        q = []
        genus = 0
        s_tmp = []
        for xi in x:
            tmp = []
            b1 = min(abs(k) for k in xi)
            for xij in xi:
                if xij > 0:
                    xij = xij - b1 + 1
                else:
                    xij = xij + b1 - 1
                tmp.append(xij)
            s_tmp.append(B(tmp))
        s = []
        for i in s_tmp:
            b = i.Tietze()
            s.append(list(b))
        t = [Link(B(si)).number_of_components() for si in s]
        for i, j in enumerate(s):
            if not j:
                s[i].append(-2)
        for i in s:
            q2 = max(abs(k) + 1 for k in i)
            q.append(q2)
        g = [((2 - t[i]) + len(x[i]) - q[i]) / 2 for i in range(len(x))]
        for i in range(len(g)):
            genus = genus + g[i]
        return Integer(genus)

    def signature(self):
        """
        Return the signature of ``self``.

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
        tot = ZZ.zero()
        s = []
        for i, j in enumerate(e):
            s.append(cmp(j, 0))
            tot = tot + s[i]
        return tot

    def alexander_polynomial(self, var='t'):
        """
        Return the Alexander polynomial of ``self``.

        INPUT:

        - ``var`` -- (default: ``'t'``) the variable in the polynomial

        EXAMPLES:

        We begin by computing the Alexander polynomial for the
        figure-eight knot::

            sage: B = BraidGroup(3)
            sage: L = Link(B([1, -2, 1, -2]))
            sage: L.alexander_polynomial()
            -t^-1 + 3 - t

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.alexander_polynomial()
            1

        Some additional examples::

            sage: B = BraidGroup(2)
            sage: L = Link(B([1]))
            sage: L.alexander_polynomial()
            1
            sage: L = Link(B.one())
            sage: L.alexander_polynomial()
            1
            sage: B = BraidGroup(3)
            sage: L = Link(B([1, 2, 1, 2]))
            sage: L.alexander_polynomial()
            t^-1 - 1 + t

        When the Seifert surface is disconnected, the Alexander
        polynomial is defined to be `0`::

            sage: B = BraidGroup(4)
            sage: L = Link(B([1,3]))
            sage: L.alexander_polynomial()
            0

        TESTS::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 3, 1, 3]))
            sage: L.alexander_polynomial()
            0
            sage: L = Link(B([1,3,1,1,3,3]))
            sage: L.alexander_polynomial()
            0
            sage: B = BraidGroup(8)
            sage: L = Link(B([-2, 4, 1, 6, 1, 4]))
            sage: L.alexander_polynomial()
            0
        """
        R = LaurentPolynomialRing(ZZ, var)
        # The Alexander polynomial of disjoint links are defined to be 0
        if len(self._braid_word_components()) > 1:
            return R.zero()
        t = R.gen()
        seifert_matrix = self.seifert_matrix()
        f = (seifert_matrix - t * seifert_matrix.transpose()).determinant()
        if f != 0:
            exp = f.exponents()
            return t ** ((-max(exp) - min(exp)) / 2) * f
        return f

    def determinant(self):
        """
        Return the determinant of ``self``.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, 2, 1, 2]))
            sage: L.determinant()
            1
            sage: B = BraidGroup(8)
            sage: L = Link(B([2, 4, 2, 3, 1, 2]))
            sage: L.determinant()
            3
            sage: L = Link(B([1]*16 + [2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.determinant()
            65

        TESTS::

            sage: Link(B([1, 2, 1, -2, -1])).determinant()
            Traceback (most recent call last):
            ...
            NotImplementedError: determinant implemented only for knots
        """
        if self.is_knot():
            a = self.alexander_polynomial()
            return Integer(abs(a(-1)))

        raise NotImplementedError("determinant implemented only for knots")

    def is_alternating(self):
        """
        Return ``True`` if the given knot diagram is alternating else
        returns ``False``.

        Alternating diagram implies every overcross is followed by an
        undercross or the vice-versa.

        We look at the Gauss code if the sign is alternating, ``True``
        is returned else the knot is not alternating ``False`` is returned.

        EXAMPLES::

            sage: B = BraidGroup(4)
            sage: L = Link(B([-1, -1, -1, -1]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1, -2, -1, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1, 3, 1, 3, 2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([1]*16 + [2,1,2,1,2,2,2,2,2,2,2,1,2,1,2,-1,2,-2]))
            sage: L.is_alternating()
            False
            sage: L = Link(B([-1,2,-1,2]))
            sage: L.is_alternating()
            True
        """
        if not self.is_knot():
            return False
        x = self.gauss_code()
        s = [cmp(i, 0) for i in x[0]]
        return (s == [(-1) ** (i + 1) for i in range(len(x[0]))]
                or s == [(-1) ** i for i in range(len(x[0]))])

    def orientation(self):
        r"""
        Return the orientation of the crossings of the link diagram
        of ``self``.

        EXAMPLES::

            sage: L = Link([[1, 4, 5, 2], [3, 5, 6, 7], [4, 8, 9, 6], [7, 9, 10, 11], [8, 1, 13, 10], [11, 13, 2, 3]])
            sage: L.orientation()
            [-1, 1, -1, 1, -1, 1]
            sage: L = Link([[1, 7, 2, 6], [7, 3, 8, 2], [3, 11, 4, 10], [11, 5, 12, 4], [14, 5, 1, 6], [13, 9, 14, 8], [12, 9, 13, 10]])
            sage: L.orientation()
            [-1, -1, -1, -1, 1, -1, 1]
            sage: L = Link([[1, 2, 3, 3], [2, 4, 5, 5], [4, 1, 7, 7]])
            sage: L.orientation()
            [-1, -1, -1]
        """
        directions = self._directions_of_edges()[0]
        orientation = []
        for C in self.pd_code():
            if C[0] == C[1] or C[2] == C[3]:
                orientation.append(-1)
            elif C[1] == C[2] or C[0] == C[3]:
                orientation.append(1)
            elif directions[C[1]] == C:
                orientation.append(-1)
            else:
                orientation.append(1)
        return orientation

    def seifert_circles(self):
        """
        Return the Seifert circles from the link diagram of ``self``.

        Seifert circles are the circles obtained by smoothing all crossings
        respecting the orientation of the segments.

        Each Seifert circle is represented as a list of the segments
        that form it.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]], [1, 1, -1, -1]])
            sage: L.seifert_circles()
            [[1, 7, 5, 3], [2, 6], [4, 8]]
            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]], [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.seifert_circles()
            [[1, 13, 9, 3, 15, 5, 11, 7], [2, 10, 6, 12], [4, 16, 8, 14]]
            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6,-5]], [-1, -1, -1, -1, 1, -1, 1]])
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
        available_segments = set(flatten(self.pd_code()))
        result = []
        tails, heads = self._directions_of_edges()
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
        Return the regions from the link diagram of ``self``.

        Regions are obtained always turning left at each crossing.

        Then the regions are represented as a list with the segments that form
        its boundary, with a sign depending on the orientation of the segment
        as part of the boundary.

        EXAMPLES::

            sage: L = Link([[[-1, +2, -3, 4, +5, +1, -2, +6, +7, 3, -4, -7, -6,-5]],[-1, -1, -1, -1, 1, -1, 1]])
            sage: L.regions()
            [[1, 7, 3, 11, 5], [2, -7], [4, -11], [6, -1], [8, -13, 10, -3], [9, 13], [12, -9, 14, -5], [-14, -8, -2, -6], [-12, -4, -10]]
            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.regions()
            [[1, 7, -4], [2, -5, -7], [3, -8, 5], [4, 8], [6, -1, -3], [-2, -6]]
            sage: L = Link([[[-1, +2, 3, -4, 5, -6, 7, 8, -2, -5, +6, +1, -8, -3, 4, -7]],[-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.regions()
            [[1, 13, -8], [2, -9, -13], [3, -14, 9], [4, 16, 8, 14], [5, 11, 7, -16], [6, -11], [10, -5, -15, -3], [12, -1, -7], [15, -4], [-12, -6, -10, -2]]
            sage: B = BraidGroup(2)
            sage: L = Link(B([-1, -1, -1]))
            sage: L.regions()
            [[1, 3, 5], [2, -1], [4, -3], [6, -5], [-2, -6, -4]]
            sage: L = Link([[[1, -2, 3, -4], [-1, 5, -3, 2, -5, 4]], [-1, 1, 1, -1, -1]])
            sage: L.regions()
            [[1, -5], [2, -8, 4, 5], [3, 8], [6, -9, -2], [7, -3, 9], [10, -4, -7], [-10, -6, -1]]
            sage: L = Link([[1, 2, 3, 3], [2, 5, 4, 4], [5, 7, 6, 6], [7, 1, 8, 8]])
            sage: L.regions()
            [[-3], [-4], [-6], [-8], [1, 2, 5, 7], [-2, 3, -1, 8, -7, 6, -5, 4]]

        .. NOTE::

            The link diagram is assumed to have only one completely isolated
            component. This is because otherwise some regions would have
            disconnected boundary.

        TESTS::

            sage: B = BraidGroup(6)
            sage: L = Link(B([1, 3, 5]))
            sage: L.regions()
            Traceback (most recent call last):
            ...
            NotImplementedError: can only have one isolated component
        """
        if len(self._isolated_components()) != 1:
            raise NotImplementedError("can only have one isolated component")
        pd = self.pd_code()
        tails, heads = self._directions_of_edges()
        available_edges = set(flatten(pd))
        if len(pd) == 1:
            if pd[0][0] == pd[0][1]:
                return [[-pd[0][2]], [pd[0][0]], [pd[0][2], -pd[0][0]]]
            else:
                return [[pd[0][2]], [-pd[0][0]], [-pd[0][2], pd[0][0]]]

        loops = [i for i in available_edges if heads[i] == tails[i]]
        available_edges = available_edges.union({-i for i in available_edges})
        regions = []

        for edge in loops:
            cros = heads[edge]
            if cros[1] == edge:
                regions.append([edge])
            else:
                regions.append([-edge])
            available_edges.remove(edge)
            available_edges.remove(-edge)

        while available_edges:
            edge = available_edges.pop()
            region = []
            while not edge in region:
                region.append(edge)
                if edge > 0 :
                    cros = heads[edge]
                    ind = cros.index(edge)
                else:
                    cros = tails[-edge]
                    ind = cros.index(-edge)
                next_edge = cros[(ind + 1) % 4]
                if [next_edge] in regions:
                    region.append(-next_edge)
                    next_edge = cros[(ind - 1) % 4]
                elif [-next_edge] in regions:
                    region.append(next_edge)
                    next_edge = cros[(ind - 1) % 4]
                if tails[next_edge] == cros:
                    edge = next_edge
                else:
                    edge = -next_edge
                if edge in available_edges:
                    available_edges.remove(edge)
            regions.append(region)
        return regions

    def writhe(self):
        """
        Return the writhe of ``self``.

        EXAMPLES::

            sage: L = Link([[[1, -2, 3, -4, 2, -1, 4, -3]],[1, 1, -1, -1]])
            sage: L.writhe()
            0
            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6,-5]],
            ....:            [-1, -1, -1, -1, 1, -1, 1]])
            sage: L.writhe()
            -3
            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
            ....:            [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L.writhe()
            -2
        """
        x = self.oriented_gauss_code()
        pos = x[1].count(1)
        neg = (-1) * x[1].count(-1)
        return pos + neg

    def jones_polynomial(self, variab=None, skein_normalization=False, algorithm='jonesrep'):
        r"""
        Return the Jones polynomial of ``self``.

        The normalization is so that the unknot has Jones polynomial `1`.
        If ``skein_normalization`` is ``True``, the variable of the result
        is replaced by a itself to the power of `4`, so that the result
        agrees with the conventions of [Lic]_ (which in particular differs
        slightly from the conventions used otherwise in this class), had
        one used the conventional Kauffman bracket variable notation directly.

        If ``variab`` is ``None`` return a polynomial in the variable `A`
        or `t`, depending on the value ``skein_normalization``. In
        particular, if ``skein_normalization`` is ``False``, return the
        result in terms of the variable `t`, also used in [Lic]_.

        ALGORITHM:

        The calculation goes through one of two possible algorithms,
        depending on the value of ``algorithm``. Possible values are
        ``'jonesrep'`` which uses the Jones representation of a braid
        representation of ``self`` to compute the polynomial of the
        trace closure of the braid, and ``statesum`` which recursively
        computes the Kauffman bracket of ``self``. Depending on how the
        link is given, there might be significant time gains in using
        one over the other. When the trace closure of the braid is
        ``self``, the algorithms give the same result.

        INPUT:

        - ``variab`` -- variable (default: ``None``); the variable in the
          resulting polynomial; if unspecified, use either a default variable
          in `\ZZ[A,A^{-1}]` or the variable `t` in the symbolic ring

        - ``skein_normalization`` -- boolean (default: ``False``); determines
          the variable of the resulting polynomial

        - ``algorithm`` -- string (default: ``'jonesrep'``); algorithm to use
          and can be one of the following:

          * ``'jonesrep'`` - use the Jones representation of the braid
            representation
          * ``'statesum'`` - recursively computes the Kauffman bracket

        OUTPUT:

        If ``skein_normalization`` if ``False``, this returns an element
        in the symbolic ring as the Jones polynomial of the link might
        have fractional powers when the link is not a knot. Otherwise the
        result is a Laurant polynomial in ``variab``.

        EXAMPLES:

        The unknot::

            sage: B = BraidGroup(9)
            sage: b = B([1, 2, 3, 4, 5, 6, 7, 8])
            sage: Link(b).jones_polynomial()
            1

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.jones_polynomial()
            1

        The Ochiai unknot::

            sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
            ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
            ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sage: L.jones_polynomial()  # long time
            1

        Two unlinked unknots::

            sage: B = BraidGroup(4)
            sage: b = B([1, 3])
            sage: Link(b).jones_polynomial()
            -sqrt(t) - 1/sqrt(t)

        The Hopf link::

            sage: B = BraidGroup(2)
            sage: b = B([-1,-1])
            sage: Link(b).jones_polynomial()
            -1/sqrt(t) - 1/t^(5/2)

        Different representations of the trefoil and one of its mirror::

            sage: B = BraidGroup(2)
            sage: b = B([-1, -1, -1])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: Link(b).jones_polynomial()
            1/t + 1/t^3 - 1/t^4
            sage: B = BraidGroup(3)
            sage: b = B([-1, -2, -1, -2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            -A^-16 + A^-12 + A^-4
            sage: R.<x> = LaurentPolynomialRing(GF(2))
            sage: Link(b).jones_polynomial(skein_normalization=True, variab=x)
            x^-16 + x^-12 + x^-4
            sage: B = BraidGroup(3)
            sage: b = B([1, 2, 1, 2])
            sage: Link(b).jones_polynomial(skein_normalization=True)
            A^4 + A^12 - A^16

        `K11n42` (the mirror of the "Kinoshita-Terasaka" knot) and `K11n34`
        (the mirror of the "Conway" knot) in [KnotAtlas]_::

            sage: B = BraidGroup(4)
            sage: K11n42 = Link(B([1, -2, 3, -2, 3, -2, -2, -1, 2, -3, -3, 2, 2]))
            sage: K11n34 = Link(B([1, 1, 2, -3, 2, -3, 1, -2, -2, -3, -3]))
            sage: cmp(K11n42.jones_polynomial(), K11n34.jones_polynomial())
            0

        The two algorithms for computation give the same result when the
        trace closure of the braid representation is the link itself::

            sage: L = Link([[[-1, 2, -3, 4, 5, 1, -2, 6, 7, 3, -4, -7, -6, -5]],
            ....:           [-1, -1, -1, -1, 1, -1, 1]])
            sage: jonesrep = L.jones_polynomial(algorithm='jonesrep')
            sage: statesum = L.jones_polynomial(algorithm='statesum')
            sage: cmp(jonesrep, statesum)
            0

        When we have thrown away unknots so that the trace closure of the
        braid is not necessarily the link itself, this is only true up to a
        power of the Jones polynomial of the unknot::

            sage: B = BraidGroup(3)
            sage: b = B([1])
            sage: L = Link(b)
            sage: b.components_in_closure()
            2
            sage: L.number_of_components()
            1
            sage: b.jones_polynomial()
            -sqrt(t) - 1/sqrt(t)
            sage: L.jones_polynomial(algorithm='statesum')
            1

        TESTS::

            sage: L = Link([])
            sage: L.jones_polynomial(algorithm='statesum')
            1

            sage: L.jones_polynomial(algorithm='other')
            Traceback (most recent call last):
            ...
            ValueError: bad value of algorithm
        """
        if algorithm == 'statesum':
            poly = self._bracket()
            t = poly.parent().gens()[0]
            writhe = self.writhe()
            jones = (poly * (-t)**(-3 * writhe))
            # Switch to the variable A to have the result agree with the output
            # of the jonesrep algorithm
            A = LaurentPolynomialRing(ZZ, 'A').gen()
            jones = jones(A**-1)

            if skein_normalization:
                if variab is None:
                    return jones
                else:
                    return jones(variab)
            else:
                if variab is None:
                    variab = 't'
                # We force the result to be in the symbolic ring because of the expand
                return jones(SR(variab)**(ZZ(1)/ZZ(4))).expand()
        elif algorithm == 'jonesrep':
            return self.braid().jones_polynomial(variab, skein_normalization)

        raise ValueError("bad value of algorithm")

    @cached_method
    def _bracket(self):
        r"""
        Return the Kaufmann bracket polynomial of the diagram of ``self``.

        Note that this is not an invariant of the link, but of the diagram.
        In particular, it is not invariant under Reidemeister I moves.

        EXAMPLES::

            sage: L = Link([[[-1, 2, 3, -4, 5, -6, 7, 8, -2, -5, 6, 1, -8, -3, 4, -7]],
            ....:           [-1, -1, -1, -1, 1, 1, -1, 1]])
            sage: L._bracket()
            -t^-10 + 2*t^-6 - t^-2 + 2*t^2 - t^6 + t^10 - t^14
            sage: L = Link([[2, 1, 3, 4], [4, 3, 1, 2]])
            sage: L._bracket()
            -t^-4 - t^4
        """
        t = LaurentPolynomialRing(ZZ, 't').gen()
        pd_code = self.pd_code()
        if not pd_code:
            return t.parent().one()
        if len(pd_code) == 1:
            if pd_code[0][0] == pd_code[0][1]:
                return -t**(-3)
            else:
                return -t**3

        cross = pd_code[0]
        rest = deepcopy(pd_code[1:])
        [a, b, c, d] = cross
        if a == b and c == d and len(rest) > 0:
            return (~t + t**(-5)) * Link(rest)._bracket()
        elif a == d and c == b and len(rest) > 0:
            return (t + t**5) * Link(rest)._bracket()
        elif a == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = c
            return -t**(-3) * Link(rest)._bracket()
        elif a == d:
            for cross in rest:
                if c in cross:
                    cross[cross.index(c)] = b
            return -t**3 * Link(rest)._bracket()
        elif c == b:
            for cross in rest:
                if d in cross:
                    cross[cross.index(d)] = a
            return -t**3 * Link(rest)._bracket()
        elif c == d:
            for cross in rest:
                if b in cross:
                    cross[cross.index(b)] = a
            return -t**(-3) * Link(rest)._bracket()
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
            return t * Link(rest)._bracket() + ~t * Link(rest_2)._bracket()

    def _isolated_components(self):
        r"""
        Return the PD codes of the isolated components of ``self``.

        Isolated components are links corresponding to subdiagrams that
        do not have any common crossing.

        EXAMPLES::

            sage: L = Link([[1, 1, 2, 2], [3, 3, 4, 4]])
            sage: L._isolated_components()
            [[[1, 1, 2, 2]], [[3, 3, 4, 4]]]
        """
        G = Graph()
        for c in self.pd_code():
            G.add_vertex(tuple(c))
        for i in range(G.num_verts()-1):
            for j in range(i, G.num_verts()):
                if len(set(G.vertices()[i]).intersection(G.vertices()[j])) > 0:
                    G.add_edge(G.vertices()[i], G.vertices()[j])
        return [[list(i) for i in j] for j in G.connected_components()]

    def plot(self, gap=0.1, component_gap=0.5, solver=None, **kwargs):
        r"""
        Plot ``self``.

        INPUT:

        - ``gap`` -- (default: 0.1) the size of the blank gap left for
          the crossings

        - ``component_gap`` -- (default: 0.5) the gap between isolated
          components

        - ``solver`` -- the linear solver to use, see
          :class:`~sage.numerical.mip.MixedIntegerLinearProgram`.

        The usual keywords for plots can be used here too.

        EXAMPLES:

        We construct the simplest version of the unknot::

            sage: L = Link([[2, 1, 1, 2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            B = BraidGroup(2)
            L = Link([[2, 1, 1, 2]])
            sphinx_plot(L.plot())

        We construct a more interesting example of the unknot::

            sage: L = Link([[2, 1, 4, 5], [3, 5, 6, 7], [4, 1, 9, 6], [9, 2, 3, 7]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[2,1,4,5], [3,5,6,7], [4,1,9,6], [9,2,3,7]])
            sphinx_plot(L.plot())

        The "monster" unknot::

            sage: L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
            ....:           [17,19,8,18],[9,10,11,14],[10,12,13,11],
            ....:           [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[3,1,2,4],[8,9,1,7],[5,6,7,3],[4,18,6,5],
                      [17,19,8,18],[9,10,11,14],[10,12,13,11],
                      [12,19,15,13],[20,16,14,15],[16,20,17,2]])
            sphinx_plot(L.plot())

        The Ochiai unknot::

            sage: L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
            ....:             -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
            ....:           [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[[1,-2,-3,-8,-12,13,-14,15,-7,-1,2,-4,10,11,-13,12,
                        -11,-16,4,3,-5,6,-9,7,-15,14,16,-10,8,9,-6,5]],
                      [-1,-1,1,1,1,1,-1,1,1,-1,1,-1,-1,-1,-1,-1]])
            sphinx_plot(L.plot())

        One of the representations of the trefoil knot::

            sage: L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
            sage: L.plot()
            Graphics object consisting of 14 graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[1, 5, 2, 4], [5, 3, 6, 2], [3, 1, 4, 6]])
            sphinx_plot(L.plot())

        The figure-eight knot::

            sage: L = Link([[2, 1, 4, 5], [5, 6, 7, 3], [6, 4, 1, 9], [9, 2, 3, 7]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[2,1,4,5], [5,6,7,3], [6,4,1,9], [9,2,3,7]])
            sphinx_plot(L.plot())

        The knot `K11n121` in [KnotAtlas]_::

            sage: L = Link([[4,2,5,1], [10,3,11,4], [5,16,6,17], [7,12,8,13],
            ....:           [18,9,19,10], [2,11,3,12], [13,20,14,21], [15,6,16,7],
            ....:           [22,18,1,17], [8,19,9,20], [21,14,22,15]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[4,2,5,1], [10,3,11,4], [5,16,6,17], [7,12,8,13],
                      [18,9,19,10], [2,11,3,12], [13,20,14,21], [15,6,16,7],
                      [22,18,1,17], [8,19,9,20], [21,14,22,15]])
            sphinx_plot(L.plot())

        One of the representations of the Hopf link::

            sage: L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[1, 4, 2, 3], [4, 1, 3, 2]])
            sphinx_plot(L.plot())

        Plotting links with multiple isolated components::

            sage: L = Link([[[-1, 2, -3, 1, -2, 3], [4, -5, 6, -4, 5, -6]], [1, 1, 1, 1, 1, 1]])
            sage: L.plot()
            Graphics object consisting of ... graphics primitives

        .. PLOT::
            :width: 300 px

            L = Link([[[-1,2,-3,1,-2,3], [4,-5,6,-4,5,-6]], [1,1,1,1,1,1]])
            sphinx_plot(L.plot())

        TESTS:

        Check that :trac:`20315` is fixed::

            sage: L = Link([[2,1,4,5], [5,6,7,3], [6,4,1,9], [9,2,3,7]])
            sage: L.plot(solver='GLPK')
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='Coin')    # optional - cbc
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='CPLEX')   # optional - CPLEX
            Graphics object consisting of ... graphics primitives
            sage: L.plot(solver='Gurobi')  # optional - Gurobi
            Graphics object consisting of ... graphics primitives
        """
        comp = self._isolated_components()
        # Handle isolated components individually
        if len(comp) > 1:
            L1 = Link(comp[0])
            L2 = Link(flatten(comp[1:], max_level=1))
            P1 = L1.plot(gap, **kwargs)
            P2 = L2.plot(gap, **kwargs)
            xtra = P1.get_minmax_data()['xmax'] + component_gap - P2.get_minmax_data()['xmin']
            for P in P2:
                if hasattr(P, 'path'):
                    for p in P.path[0]:
                        p[0] += xtra
                    for p in P.vertices:
                        p[0] += xtra
                else:
                    P.xdata = [p + xtra for p in P.xdata]
            return P1 + P2

        if not 'color' in kwargs:
            kwargs['color'] = 'blue'
        if not 'axes' in kwargs:
            kwargs['axes'] = False
        if not 'aspect_ratio' in kwargs:
            kwargs['aspect_ratio'] = 1

        from sage.plot.line import line
        from sage.plot.bezier_path import bezier_path
        from sage.plot.circle import circle

        # Special case for the unknot
        if not self.pd_code():
            return circle((0,0), ZZ(1)/ZZ(2), **kwargs)

        # The idea is the same followed in spherogram, but using MLP instead of
        # network flows.
        # We start by computing a way to bend the edges left or right
        # such that the resulting regions are in fact closed regions
        # with straight angles, and using the minimal number of bends.
        regions = sorted(self.regions(), key=len)
        regions = regions[:-1]
        edges = list(set(flatten(self.pd_code())))
        edges.sort()
        MLP = MixedIntegerLinearProgram(maximization=False, solver=solver)
        # v will be the list of variables in the MLP problem. There will be
        # two variables for each edge: number of right bendings and number of
        # left bendings (at the end, since we are minimizing the total, only one
        # of each will be nonzero
        v = MLP.new_variable(nonnegative=True, integer=True)
        # one condition for each region
        for i in range(len(regions)):
            cond = 0
            r = regions[i]
            for e in r:
                if e > 0:
                    cond = cond + v[2*edges.index(e)] - v[2*edges.index(e) + 1]
                else:
                    cond = cond - v[2*edges.index(-e)] + v[2*edges.index(-e) + 1]
            MLP.add_constraint(cond == 4 - len(r))
        MLP.set_objective(MLP.sum(v.values()))
        MLP.solve()
        # we store the result in a vector s packing right bends as negative left ones
        s = range(len(edges))
        values = MLP.get_values(v)
        for i in range(len(edges)):
            s[i] = int(values[2*i] - values[2*i + 1])
        # segments represents the different parts of the previos edges after bending
        segments = {e: [(e,i) for i in range(abs(s[edges.index(e)])+1)] for e in edges}
        pieces = {tuple(i): [i] for j in segments.values() for i in j}
        nregions = []
        for r in regions:
            nregion = []
            for e in r:
                if e > 0:
                    rev = segments[e][:-1]
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
        badregions = [nr for nr in nregions if any(-1 == x[1] for x in nr)]
        while len(badregions) > 0:
            badregion = badregions[0]
            badturns = []
            a = 0
            while badregion[a][1] != -1:
                a += 1
            c = -1
            b = a
            while c != 2:
                if b == len(badregion)-1:
                    b = 0
                else:
                    b += 1
                c += badregion[b][1]
            otherregion = [nr for nr in nregions
                           if any(badregion[b][0] == x[0] for x in nr)]
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
                segmenttoadd = [x for x in pieces.keys()
                                if badregion[b][0] in pieces[x]]
                if len(segmenttoadd) > 0:
                    pieces[segmenttoadd[0]].append(N2)
            else:
                pieces[tuple(badregion[b][0])].append(N2)

            if a < b:
                r1 = badregion[:a] + [[badregion[a][0],0], [N1,1]] + badregion[b:]
                r2 = badregion[a+1:b] + [[N2,1],[N1,1]]
            else:
                r1 = badregion[b:a] + [[badregion[a][0],0], [N1,1]]
                r2 = badregion[:b] + [[N2,1],[N1,1]] + badregion[a+1:]

            if otherregion:
                c = [x for x in otherregion if badregion[b][0] == x[0]]
                c = otherregion.index(c[0])
                otherregion.insert(c+1, [N2,otherregion[c][1]])
                otherregion[c][1] = 0
            nregions.remove(badregion)
            nregions.append(r1)
            nregions.append(r2)
            badregions = [nr for nr in nregions if any(x[1] == -1 for x in nr)]
        MLP = MixedIntegerLinearProgram(maximization=False, solver=solver)
        v = MLP.new_variable(nonnegative=True, integer=True)
        for e in segments:
            MLP.set_min(v[e], 1)
        for r in nregions:
            horp = []
            horm = []
            verp = []
            verm = []
            direction = 0
            for se in r:
                if direction % 4 == 0:
                    horp.append(v[se[0]])
                elif direction == 1:
                    verp.append(v[se[0]])
                elif direction == 2:
                    horm.append(v[se[0]])
                elif direction == 3:
                    verm.append(v[se[0]])
                if se[1] == 1:
                    direction += 1
            MLP.add_constraint(MLP.sum(horp) - MLP.sum(horm) == 0)
            MLP.add_constraint(MLP.sum(verp) - MLP.sum(verm) == 0)
        MLP.set_objective(MLP.sum(v.values()))
        solved = MLP.solve()
        v = MLP.get_values(v)
        lengths = {piece: sum(v[a] for a in pieces[piece]) for piece in pieces}
        image = line([], **kwargs)
        crossings = {tuple(self.pd_code()[0]): (0,0,0)}
        availables = self.pd_code()[1:]
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
            orien = self.orientation()[self.pd_code().index(list(c))]
            if s[edges.index(e)] < 0:
                turn = -1
            else:
                turn = 1
            lengthse = [lengths[(e,i)] for i in range(abs(s[edges.index(e)])+1)]
            if c.index(e) == 0 or (c.index(e) == 1 and orien == 1) or (c.index(e) == 3 and orien == -1):
                turn = -turn
                lengthse.reverse()
            tailshort = (c.index(e) % 2 == 0)
            x0 = crossings[c][0]
            y0 = crossings[c][1]
            im = []
            for l in lengthse:
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
                direction = (direction + turn) % 4
                x0 = x1
                y0 = y1
            direction = (direction - turn) % 4
            c2 = [ee for ee in availables if e in ee]
            if len(c2) == 1:
                availables.remove(c2[0])
                crossings[tuple(c2[0])] = (x1, y1, (direction + c2[0].index(e) + 2) % 4)
            c2 = [ee for ee in self.pd_code() if e in ee and ee != list(c)]
            if not c2:
                headshort = not tailshort
            else:
                headshort = (c2[0].index(e) % 2 == 0)
            a = deepcopy(im[0][0])
            b = deepcopy(im[-1][0])
            if tailshort:
                im[0][0][0][0] += cmp(a[1][0], im[0][0][0][0]) * gap
                im[0][0][0][1] += cmp(a[1][1], im[0][0][0][1]) * gap
            if headshort:
                im[-1][0][1][0] -= cmp(b[1][0], im[-1][0][0][0]) * gap
                im[-1][0][1][1] -= cmp(b[1][1], im[-1][0][0][1]) * gap
            l = line([], **kwargs)
            c = 0
            p = im[0][0][0]
            if len(im) == 4 and max([x[1] for x in im]) == 1:
                l = bezier_path([[im[0][0][0], im[0][0][1], im[-1][0][0], im[-1][0][1]]], **kwargs)
                p = im[-1][0][1]
            else:
                while c < len(im)-1:
                    if im[c][1] > 1:
                        (a, b) = im[c][0]
                        if b[0] > a[0]:
                            e = [b[0] - 1, b[1]]
                        elif b[0] < a[0]:
                            e = [b[0] + 1, b[1]]
                        elif b[1] > a[1]:
                            e = [b[0], b[1] - 1]
                        elif b[1] < a[1]:
                            e = [b[0] , b[1] + 1]
                        l += line((p, e), **kwargs)
                        p = e
                    if im[c+1][1] == 1 and c < len(im) - 2:
                        xr = round(im[c+2][0][1][0])
                        yr = round(im[c+2][0][1][1])
                        xp = xr - im[c+2][0][1][0]
                        yp = yr - im[c+2][0][1][1]
                        q = [p[0] + im[c+1][0][1][0] - im[c+1][0][0][0] - xp,
                             p[1] + im[c+1][0][1][1] - im[c+1][0][0][1] - yp]
                        l += bezier_path([[p, im[c+1][0][0], im[c+1][0][1], q]], **kwargs)
                        c += 2
                        p = q
                    else:
                        if im[c+1][1] == 1:
                            q = im[c+1][0][1]
                        else:
                            q = [im[c+1][0][0][0] + sign(im[c+1][0][1][0] - im[c+1][0][0][0]),
                                 im[c+1][0][0][1] + sign(im[c+1][0][1][1] - im[c+1][0][0][1])]
                        l += bezier_path([[p, im[c+1][0][0], q]], **kwargs)
                        p = q
                        c += 1
            l += line([p, im[-1][0][1]], **kwargs)
            image += l
            ims += sum(line(a[0], **kwargs) for a in im)
        return image

