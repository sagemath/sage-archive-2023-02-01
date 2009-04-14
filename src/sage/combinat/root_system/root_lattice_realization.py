"""
Root lattice realization
"""
#*****************************************************************************
#       Copyright (C) 2007 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.combinat.family import Family
from sage.rings.all import ZZ
from sage.misc.all import cached_method

class RootLatticeRealization(object):
    def index_set(self):
        """
        EAMPLES:
           sage: r = RootSystem(['A',4]).root_space()
           sage: r.index_set()
           [1, 2, 3, 4]
        """
        return self.root_system.index_set()

    def dynkin_diagram(self):
        """
        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_space()
            sage: r.dynkin_diagram()
            Dynkin diagram of type ['A', 4]
        """
        return self.root_system.dynkin_diagram()

    def _name_string_helper(self, name, capitalize=True, base_ring=True, type=True):
        """
        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_space()
            sage: r._name_string_helper("root")
            "Root space over the Rational Field of the Root system of type ['A', 4]"
            sage: r._name_string_helper("root", base_ring=False)
            "Root space of the Root system of type ['A', 4]"
            sage: r._name_string_helper("root", base_ring=False, type=False)
            'Root space'
            sage: r._name_string_helper("root", capitalize=False, base_ring=False, type=False)
            'root space'
        """
        s = ""
        if self.root_system.dual_side:
            s += "co"

        s += name + " "

        if self.base_ring() == ZZ:
            s += "lattice "
        else:
            s += "space "
            if base_ring:
                s += "over the %s "%self.base_ring()

        if type:
            s += "of the "
            if self.root_system.dual_side:
                s += repr(self.root_system.dual)
            else:
                s += repr(self.root_system)

        if capitalize:
            s = s[:1].upper() + s[1:]


        return s.strip()

    ##########################################################################
    # checks
    ##########################################################################

    def check(self):
        """
        Checks to make sure that the scalar product between the simple root
        and the simple coroot in this root lattice realization is the
        corresponding entry in the Dynkin diagram.

        EXAMPLES:
            sage: RootSystem(['A',3]).root_lattice().check()
            True
        """
        alpha = self.simple_roots()
        alphacheck = self.simple_coroots()
        dynkin_diagram = self.dynkin_diagram()
        for i in self.index_set():
            for j in self.index_set():
                assert alpha[j].scalar(alphacheck[i]) == dynkin_diagram[i,j]
        return True

    ##########################################################################
    # highest root
    ##########################################################################

    def highest_root(self):
        """
        Returns the highest root (for an irreducible finite root system)

        EXAMPLES:
            sage: RootSystem(['A',4]).ambient_space().highest_root()
            (1, 0, 0, 0, -1)

            sage: RootSystem(['E',6]).weight_space().highest_root()
            Lambda[2]

        """
        assert(self.root_system.is_finite())
        assert(self.root_system.is_irreducible())
        return self.a_long_simple_root().to_positive_chamber()

    def a_long_simple_root(self):
        """
        Returns a long simple root, corresponding to the highest outgoing edge
        in the Dynkin diagram.

        Caveat: This depends on a direct identification between dynkin diagram
        nodes and the simple roots, which is not yet valid in the ambient space

        Caveat: this may be break in affine type A_2n^(2)

        Caveat: meaningful/broken for non irreducible?

        TODO: implement DynkinDiagram.vertices_by_length as in
        MuPAD-Combinat, and use it here

        TESTS:
            sage: X=RootSystem(['A',1]).weight_space()
            sage: X.a_long_simple_root()
            2*Lambda[1]
            sage: X=RootSystem(['A',5]).weight_space()
            sage: X.a_long_simple_root()
            2*Lambda[1] - Lambda[2]
        """
        if self.dynkin_diagram().rank() == 1:
            return self.simple_roots()[self.index_set()[0]]
        longest=self.dynkin_diagram().outgoing_edges()[0]
        for j in self.dynkin_diagram().outgoing_edges():
            if j[2]>longest[2]:
                longest=j
        return self.simple_roots()[longest[0]]


    ##########################################################################
    # simple roots
    ##########################################################################

    def simple_root(self, i):
        """
        Returns the $i$-th simple root.  This should be overridden by any
        subclass.

        EXAMPLES:
            sage: r = RootSystem(["A",3]).root_lattice()
            sage: r.simple_root(1)
            alpha[1]

        TESTS:
            sage: from sage.combinat.root_system.root_lattice_realization import RootLatticeRealization
            sage: RootLatticeRealization.simple_root(r, 1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError

    def simple_roots(self):
        """
        Returns the family $(\alpha_i)_{i\in I}$ of the simple roots.

        EXAMPLES:
            sage: alpha = RootSystem(["A",3]).root_lattice().simple_roots()
            sage: [alpha[i] for i in [1,2,3]]
            [alpha[1], alpha[2], alpha[3]]
        """
        if not hasattr(self,"_simple_roots"):
            self._simple_roots = Family(self.index_set(), self.simple_root)
            # self._simple_roots.rename("alpha")
            # break some doctests
        return self._simple_roots

    def alpha(self):
        """
        Returns the family $(\alpha_i)_{i\in I}$ of the simple roots,
        with the extra feature that, for simple irreducible root
        systems, $\alpha_0$ yields the opposite of the highest root.

        EXAMPLES:
            sage: alpha = RootSystem(["A",2]).root_lattice().alpha()
            sage: alpha[1]
            alpha[1]
            sage: alpha[0]
            alpha[1] + alpha[2]

        """
        if self.root_system.is_finite() and self.root_system.is_irreducible():
            return Family(self.index_set(), self.simple_root, \
                          hidden_keys = [0], hidden_function = lambda i: self.highest_root())
        else:
            return self.simple_roots()

    ##########################################################################
    # roots
    ##########################################################################

    def roots(self):
        """
        Returns the roots of self.

        EXAMPLES:
            sage: RootSystem(['A',2]).ambient_lattice().roots()
            [(1, -1, 0), (1, 0, -1), (0, 1, -1), (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]
        """
        return self.positive_roots() + self.negative_roots()

    ##########################################################################
    # coroots
    ##########################################################################

    def coroot_lattice(self):
        """
        Returns the coroot lattice.

        EXAMPLES:
            sage: RootSystem(['A',2]).root_lattice().coroot_lattice()
            Coroot lattice of the Root system of type ['A', 2]

        """
        return self.root_system.coroot_lattice()

    def simple_coroot(self, i):
        """
        Returns the $i^{th}$ simple coroot.

        EXAMPLES:
            sage: RootSystem(['A',2]).root_lattice().simple_coroot(1)
            alphacheck[1]
        """
        return self.coroot_lattice().simple_root(i)

    def simple_coroots(self):
        """
        Returns the family $(\alpha^\vee_i)_{i\in I}$ of the simple coroots.

        EXAMPLES:
            sage: alphacheck = RootSystem(['A',3]).root_lattice().simple_coroots()
            sage: [alphacheck[i] for i in [1, 2, 3]]
            [alphacheck[1], alphacheck[2], alphacheck[3]]

        TEST:
            sage: alphacheck
            alphacheck
        """
        if not hasattr(self,"cache_simple_coroots"):
            self.cache_simple_coroots = Family(self.index_set(), self.simple_coroot)
            self.cache_simple_coroots.rename("alphacheck")
        return self.cache_simple_coroots

    def alphacheck(self):
        """
        Returns the family $(\alpha^\vee_i)_{i\in I}$ of the simple
        coroots, with the extra feature that for simple irreducible
        root systems, $\alpha^\vee_0$ yields the coroot associated to
        the opposite of the highest root (caveat: this is usually not
        the highest coroot!)

        EXAMPLES:
            sage: alphacheck = RootSystem(["A",3]).ambient_space().alphacheck()
            sage: alphacheck[1]
            (1, -1, 0, 0)
            sage: alphacheck[0]
            (1, 0, 0, -1)
        """
        if self.root_system.is_finite() and self.root_system.is_irreducible():
            return Family(self.index_set(), self.simple_coroot, \
                          hidden_keys = [0], hidden_function = lambda i: self.cohighest_root())
        else:
            return self.simple_roots()

    @cached_method
    def cohighest_root(self):
        """
        Returns the associated coroot of the highest root.  Note that this is
        usually not the highest coroot.

        EXAMPLES:
            sage: RootSystem(['A', 3]).ambient_space().cohighest_root()
            (1, 0, 0, -1)
        """
        return self.highest_root().associated_coroot()


    ##########################################################################
    # reflections
    ##########################################################################

    def reflection(self, root, coroot=None):
        """
        Returns the reflection along the root, and across the
        hyperplane define by coroot, as a function from
        self to self.

        EXAMPLES:
            sage: space = RootSystem(['A',2]).weight_lattice()
            sage: x=space.simple_roots()[1]
            sage: y=space.simple_coroots()[1]
            sage: s = space.reflection(x,y)
            sage: x
            2*Lambda[1] - Lambda[2]
            sage: s(x)
            -2*Lambda[1] + Lambda[2]
            sage: s(-x)
            2*Lambda[1] - Lambda[2]
        """
        if coroot is None:
            coroot = root.associated_coroot()
        return lambda v: v - v.scalar(coroot) * root

    @cached_method
    def simple_reflection(self, i):
        """
        Returns the $i^{th}$ simple reflection, as a function from
        self to self.

        INPUT:
            i -- i is in self's index set

        EXAMPLES:
            sage: space = RootSystem(['A',2]).ambient_lattice()
            sage: s = space.simple_reflection(1)
            sage: x = space.simple_roots()[1]
            sage: x
            (1, -1, 0)
            sage: s(x)
            (-1, 1, 0)
        """
        return self.reflection(self.simple_root(i), self.simple_coroot(i))

    def simple_reflections(self):
        """
        Returns the family $(s_i)_{i\in I}$ of the simple reflections
        of this root system.

        EXAMPLES:
            sage: r = RootSystem(["A", 2]).root_lattice()
            sage: s = r.simple_reflections()
            sage: s[1]( r.simple_root(1) )
            -alpha[1]

        TEST:
            sage: s
            simple reflections
        """
        res =  self.alpha().zip(self.reflection, self.alphacheck())
        res.rename("simple reflections")
        return res

    s = simple_reflections

    ##########################################################################
    # projections
    ##########################################################################

    def projection(self, root, coroot=None, to_negative=True):
        """
        Returns the projection along the root, and across the
        hyperplane define by coroot, as a function $\pi$ from self to
        self. $\pi$ is a half-linear map which stabilizes the negative
        half space, and acts by reflection on the positive half space.

        If to_negative is False, then this project onto the positive
        half space instead.

        EXAMPLES:
            sage: space = RootSystem(['A',2]).weight_lattice()
            sage: x=space.simple_roots()[1]
            sage: y=space.simple_coroots()[1]
            sage: pi = space.projection(x,y)
            sage: x
            2*Lambda[1] - Lambda[2]
            sage: pi(x)
            -2*Lambda[1] + Lambda[2]
            sage: pi(-x)
            -2*Lambda[1] + Lambda[2]
            sage: pi = space.projection(x,y,False)
            sage: pi(-x)
            2*Lambda[1] - Lambda[2]
        """
        if coroot is None:
            coroot = root.associated_coroot()

        return lambda v: v - v.scalar(coroot) * root if ((v.scalar(coroot) > 0) == to_negative) else v

    def simple_projection(self, i, to_negative=True):
        """
        Returns the projection along the $i^{th}$ simple root, and across the
        hyperplane define by the $i^{th}$ simple coroot, as a function from
        self to self.

        INPUT:
            i -- i is in self's index set

        EXAMPLES:
            sage: space = RootSystem(['A',2]).weight_lattice()
            sage: x = space.simple_roots()[1]
            sage: pi = space.simple_projection(1)
            sage: x
            2*Lambda[1] - Lambda[2]
            sage: pi(x)
            -2*Lambda[1] + Lambda[2]
            sage: pi(-x)
            -2*Lambda[1] + Lambda[2]
            sage: pi = space.simple_projection(1,False)
            sage: pi(-x)
            2*Lambda[1] - Lambda[2]
        """
        return self.projection(self.simple_root(i), self.simple_coroot(i), to_negative)

    def simple_projections(self):
        """
        Returns the family $(s_i)_{i\in I}$ of the simple projections
        of this root system

        EXAMPLES:
            sage: space = RootSystem(['A',2]).weight_lattice()
            sage: pi = space.simple_projections()
            sage: x = space.simple_roots()
            sage: pi[1](x[2])
            -Lambda[1] + 2*Lambda[2]

        TESTS:
            sage: pi
            pi
        """
        res = self.alpha().zip(self.projection, self.alphacheck())
        res.rename("pi")
        return res

    pi = simple_projections

    ##########################################################################
    # Weyl group
    ##########################################################################

    def weyl_group(self):
        """
        Returns the Weyl group associated to self.

        EXAMPLES:
            sage: RootSystem(['F',4]).ambient_space().weyl_group()
            Weyl Group of type ['F', 4] (as a matrix group acting on the ambient space)
            sage: RootSystem(['F',4]).root_space().weyl_group()
            Weyl Group of type ['F', 4] (as a matrix group acting on the root space)

        """
        from sage.combinat.root_system.weyl_group import WeylGroup
        return WeylGroup(self)

class RootLatticeRealizationElement(object):
    def scalar(self, lambdacheck):
        """
        The natural pairing between this and the coroot lattice.  This should be overridden
        in subclasses.

        EXAMPLES:
            sage: r = RootSystem(['A',4]).root_lattice()
            sage: cr = RootSystem(['A',4]).coroot_lattice()
            sage: a1 = r.simple_root(1)
            sage: ac1 = cr.simple_root(1)
            sage: a1.scalar(ac1)
            2

        TESTS:
            sage: from sage.combinat.root_system.root_lattice_realization import RootLatticeRealizationElement
            sage: RootLatticeRealizationElement.scalar(a1, ac1)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        raise NotImplementedError


    def simple_reflection(self, i):
        """
        The image of self by the $i$-th simple reflection.

        EXAMPLES:
            sage: alpha = RootSystem(["A", 3]).root_lattice().alpha()
            sage: alpha[1].simple_reflection(2)
            alpha[1] + alpha[2]

        """
        # Subclasses should optimize whenever possible!
        return self.parent().simple_reflections()[i](self)


    def associated_coroot(self):
        """
        Returns the coroot associated to this root.

        EXAMPLES:
            sage: alpha = RootSystem(["A", 3]).root_lattice().alpha()
            sage: alpha[1].associated_coroot()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        #assert(self in self.parent().roots() is not False)
        raise NotImplementedError

    ##########################################################################
    # Descents
    ##########################################################################

    def has_descent(self, i, positive=False):
        """
        Test if self has a descent at position $i$, that is if self is
        on the strict negative side of the $i$-th simple reflection
        hyperplane.

        If positive if True, tests if it is on the strict positive
        side instead.

        EXAMPLES:
            sage: space=RootSystem(['A',5]).weight_space()
            sage: alpha=RootSystem(['A',5]).weight_space().simple_roots()
            sage: [alpha[i].has_descent(1) for i in space.index_set()]
            [False, True, False, False, False]
            sage: [(-alpha[i]).has_descent(1) for i in space.index_set()]
            [True, False, False, False, False]
            sage: [alpha[i].has_descent(1, True) for i in space.index_set()]
            [True, False, False, False, False]
            sage: [(-alpha[i]).has_descent(1, True) for i in space.index_set()]
            [False, True, False, False, False]
            sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(3)
            True
            sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(1)
            False
            sage: (alpha[1]+alpha[2]+alpha[4]).has_descent(1, True)
            True
        """
        s = self.scalar(self.parent().simple_coroots()[i])
        if positive:
            return s > 0
        else:
            return s < 0

    def first_descent(self, index_set=None, positive=False):
        """
        Returns the first descent of pt

        One can use the index_set option to restrict to the parabolic
        subgroup indexed by index_set.

        EXAMPLES:
            sage: space=RootSystem(['A',5]).weight_space()
            sage: alpha=space.simple_roots()
            sage: (alpha[1]+alpha[2]+alpha[4]).first_descent()
            3
            sage: (alpha[1]+alpha[2]+alpha[4]).first_descent([1,2,5])
            5
            sage: (alpha[1]+alpha[2]+alpha[4]).first_descent([1,2,5,3,4])
            5
        """
        if index_set == None:
            index_set = self.parent().index_set()
        for i in index_set:
            if self.has_descent(i, positive):
                return i
        return None

    def descents(self, index_set=None, positive=False):
        """
        Returns the descents of pt

        EXAMPLES:
            sage: space=RootSystem(['A',5]).weight_space()
            sage: alpha=space.simple_roots()
            sage: (alpha[1]+alpha[2]+alpha[4]).descents()
            [3, 5]
        """
        if index_set==None:
            index_set=self.parent().index_set()
        return [ i for i in index_set if self.has_descent(i, positive) ]

    def to_positive_chamber(self, index_set = None, positive = True):
        """
        Returns the unique element of the orbit of pt in the positive
        chamber.

        With the index_set optional parameter, this is done with
        respect to the corresponding parbolic subgroup

        With positive = False, returns the unique element in the
        negative chamber instead


        EXAMPLES:
            sage: space=RootSystem(['A',5]).weight_space()
            sage: alpha=RootSystem(['A',5]).weight_space().simple_roots()
            sage: alpha[1].to_positive_chamber()
            Lambda[1] + Lambda[5]
            sage: alpha[1].to_positive_chamber([1,2])
            Lambda[1] + Lambda[2] - Lambda[3]
        """
        if index_set is None:
            index_set=self.parent().index_set()
        while True:
            # The first index where it is *not* yet on the positive side
            i = self.first_descent(index_set, positive=(not positive))
            if i is None:
                return self
            else:
                self = self.simple_reflection(i)

    def is_dominant(self, index_set = None, positive = True):
        """
        Returns whether self is dominant.

        With positive = False, returns whether self is antidominant

        INPUT:
            v -- an element of the lattice

        EXAMPLES:
            sage: L = RootSystem(['A',2]).ambient_lattice()
            sage: Lambda = L.fundamental_weights()
            sage: [x.is_dominant() for x in Lambda]
            [True, True]
            sage: [x.is_dominant(positive=False) for x in Lambda]
            [False, False]
            sage: (Lambda[1]-Lambda[2]).is_dominant()
            False
            sage: (-Lambda[1]+Lambda[2]).is_dominant()
            False
            sage: (Lambda[1]-Lambda[2]).is_dominant([1])
            True
            sage: (Lambda[1]-Lambda[2]).is_dominant([2])
            False
            sage: [x.is_dominant() for x in L.roots()]
            [False, True, False, False, False, False]
            sage: [x.is_dominant(positive=False) for x in L.roots()]
            [False, False, False, False, True, False]
        """
        return self.first_descent(index_set, not positive) is None
