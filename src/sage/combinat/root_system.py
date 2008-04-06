#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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

import cartan_type
from sage.combinat.combinatorial_algebra import CombinatorialAlgebra
from sage.modules.free_module import FreeModule
from sage.rings.all import ZZ
from sage.misc.misc import prod

def RootSystem(t):
    """
    Returns the root system associated to the Cartan type t.

    EXAMPLES:
        sage: RootSystem(['A',3])
        Root system of type ['A', 3]
    """
    ct = cartan_type.CartanType(t)
    if not ct.affine:
        typ = ct.type()
        if typ == "A":
            return RootSystem_a(ct)
        elif typ == "B":
            return RootSystem_b(ct)
        elif typ == "C":
            return RootSystem_c(ct)
        elif typ == "D":
            return RootSystem_d(ct)
        elif typ == "E":
            return RootSystem_e(ct)
        elif typ == "F":
            return RootSystem_f(ct)
        elif typ == "G":
            return RootSystem_g(ct)
    else:
        return RootSystem_generic(ct)

class RootSystem_generic:
    def __init__(self, ct):
        """
        TESTS:
            sage: R = RootSystem(['A',3])
            sage: loads(dumps(R))
            Root system of type ['A', 3]
        """
        self.ct = ct

    def __repr__(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',3])
            Root system of type ['A', 3]
        """
        return "Root system of type %s"%self.ct

    def cartan_type(self):
        """
        Returns the Cartan type of the root system.

        EXAMPLES:
            sage: R = RootSystem(['A',3])
            sage: R.cartan_type()
            ['A', 3]
        """
        return self.ct

    def cartan_matrix(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',3]).cartan_matrix()
            [ 2 -1  0]
            [-1  2 -1]
            [ 0 -1  2]
        """
        return self.cartan_type().cartan_matrix()

    def index_set(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',3]).index_set()
            [1, 2, 3]
        """
        return self.cartan_type().index_set()

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: r1 = RootSystem(['A',3])
            sage: r2 = RootSystem(['B',3])
            sage: r1 == r1
            True
            sage: r1 == r2
            False
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self.ct != other.ct:
            return cmp(self.ct, other.ct)
        return 0

class RootSystem_a(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',4]).ambient_lattice()
            Ambient lattice of the root system of type ['A', 4]
        """
        return AmbientLattice_a(self.ct)

class RootSystem_b(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',4]).ambient_lattice()
            Ambient lattice of the root system of type ['B', 4]
        """
        return AmbientLattice_b(self.ct)

class RootSystem_c(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',4]).ambient_lattice()
            Ambient lattice of the root system of type ['C', 4]
        """
        return AmbientLattice_c(self.ct)

class RootSystem_d(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_lattice()
            Ambient lattice of the root system of type ['D', 4]
        """
        return AmbientLattice_d(self.ct)

class RootSystem_e(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['E',6]).ambient_lattice()
            Ambient lattice of the root system of type ['E', 6]
        """
        return AmbientLattice_e(self.ct)

class RootSystem_f(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['F',4]).ambient_lattice()
            Ambient lattice of the root system of type ['F', 4]
        """
        return AmbientLattice_f(self.ct)

class RootSystem_g(RootSystem_generic):
    def ambient_lattice(self):
        """
        EXAMPLES:
            sage: RootSystem(['G',2]).ambient_lattice()
            Ambient lattice of the root system of type ['G', 2]
        """
        return AmbientLattice_g(self.ct)



class WeightLatticeRealization_class:
    # Should this be a method or an attribute?
    # same question for the roots, ...
    def rho(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',3]).ambient_lattice().rho()
            (3, 2, 1, 0)
        """
        return sum(self.fundamental_weights())

    # Should it be a method of highest_weight?
    def weyl_dimension(self, highest_weight):
        """
        EXAMPLES:
            sage: hwv = vector([2,1,0,0])
            sage: RootSystem(['A',3]).ambient_lattice().weyl_dimension(hwv)
            20
        """
        # Should assert(highest_weight.is_dominant())
        rho = self.rho()
        n = prod([(rho+highest_weight).dot_product(x) for x in self.positive_roots()])
        d = prod([ rho.dot_product(x) for x in self.positive_roots()])
        return n/d

class AmbientLattice_generic(WeightLatticeRealization_class):
    def __init__(self, ct):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e == loads(dumps(e))
            True
        """
        if not hasattr(self, 'n'):
            self.n  = ct.rank()
        self.ct = ct
        self._free_module = FreeModule(ZZ, self.n)

    def __repr__(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',4]).ambient_lattice()
            Ambient lattice of the root system of type ['B', 4]

        """
        return "Ambient lattice of the root system of type %s"%self.ct

    def __getitem__(self,i):
        """
        Note that indexing starts at 1.

        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_lattice()
            sage: e[1]
            (1, 0, 0)
        """
        return self._term(i-1)

    def roots(self):
        """
        Returns the roots of self.

        EXAMPLES:
            sage: RootSystem(['A',2]).ambient_lattice().roots()
            [(1, -1, 0), (1, 0, -1), (0, 1, -1), (-1, 1, 0), (-1, 0, 1), (0, -1, 1)]
        """
        return self.positive_roots() + self.negative_roots()

    def _term(self, i):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['A',2]).ambient_lattice()
            sage: e._term(0)
            (1, 0, 0)
        """
        return self._free_module.gen(i)

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: e1 = RootSystem(['A',3]).ambient_lattice()
            sage: e2 = RootSystem(['B',3]).ambient_lattice()
            sage: e1 == e1
            True
            sage: e1 == e2
            False
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self.ct != other.ct:
            return cmp(self.ct, other.ct)
        return 0

class AmbientLattice_a(AmbientLattice_generic):
    def __init__(self, ct):
        """
        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.n
            4
        """
        self.n  = ct.rank()+1
        AmbientLattice_generic.__init__(self, ct)

    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['A',3]).ambient_lattice()
            sage: e.root(0,1)
            (1, -1, 0, 0)
        """
        return self._term(i) - self._term(j)

    def simple_roots(self):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.simple_roots()
            [(1, -1, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1)]
        """
        return [ self.root(i, i+1) for i in range(self.n-1) ]

    def negative_roots(self):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 1, 0),
             (0, -1, 0, 1),
             (0, 0, -1, 1)]
        """
        res = []
        for j in range(self.n-1):
            for i in range(j+1,self.n):
                res.append(  self.root(i,j) )
        return res

    def positive_roots(self):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.positive_roots()
            [(1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]

        """
        res = []
        for j in range(self.n):
            for i in range(j):
                res.append(  self.root(i,j) )
        return res

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: e = CartanType(['A',3]).root_system().ambient_lattice()
            sage: e.fundamental_weights()
            [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0)]

        """
        return [ sum([self._term(j) for j in range(i+1)]) for i in range(self.n-1)]

class AmbientLattice_b(AmbientLattice_generic):
    def root(self, i, j):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['B',3]).ambient_lattice()
            sage: e.root(0,1)
            (1, -1, 0)

        """
        return self._term(i) - self._term(j)

    def simple_roots(self):
        """
        EXAMPLES:
            sage: e = RootSystem(['B',4]).ambient_lattice()
            sage: e.simple_roots()
            [(1, -1, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1), (0, 0, 0, 1)]
            sage: e.positive_roots()
            [(1, -1, 0, 0),
            (1, 1, 0, 0),
            (1, 0, -1, 0),
            (1, 0, 1, 0),
            (1, 0, 0, -1),
            (1, 0, 0, 1),
            (0, 1, -1, 0),
            (0, 1, 1, 0),
            (0, 1, 0, -1),
            (0, 1, 0, 1),
            (0, 0, 1, -1),
            (0, 0, 1, 1),
            (1, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 0, 1, 0),
            (0, 0, 0, 1)]
            sage: e.fundamental_weights()
            [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0), (1/2, 1/2, 1/2, 1/2)]

        """
        return [ self.root(i,i+1) for i in range(self.n-1) ] + [ self._term(self.n-1) ]

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_lattice().negative_roots()
            [(-1, 1, 0),
             (-1, -1, 0),
             (-1, 0, 1),
             (-1, 0, -1),
             (0, -1, 1),
             (0, -1, -1),
             (-1, 0, 0),
             (0, -1, 0),
             (0, 0, -1)]

        """
        return [ -a for a in self.positive_roots()]


    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_lattice().positive_roots()
            [(1, -1, 0),
             (1, 1, 0),
             (1, 0, -1),
             (1, 0, 1),
             (0, 1, -1),
             (0, 1, 1),
             (1, 0, 0),
             (0, 1, 0),
             (0, 0, 1)]

        """
        res = []
        for i in range(self.n-1):
            for j in range(i+1,self.n):
                res.append(self._term(i) - self._term(j))
                res.append(self._term(i) + self._term(j))
        for i in range(self.n):
            res.append(self._term(i))
        return res

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: RootSystem(['B',3]).ambient_lattice().fundamental_weights()
            [(1, 0, 0), (1, 1, 0), (1/2, 1/2, 1/2)]

        """
        return [ sum(self._term(j) for j in range(i+1)) for i in range(self.n-1)]\
               + [ sum( self._term(j) for j in range(self.n) ) / 2 ]


class AmbientLattice_c(AmbientLattice_generic):
    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['C',3]).ambient_lattice()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
        """
        return (-1)**p1 * self._term(i) + (-1)**p2 * self._term(j)

    def simple_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_lattice().simple_roots()
            [(1, -1, 0), (0, 1, -1), (0, 0, 2)]

        """
        return [ self.root(i, i+1,0,1) for i in range(self.n-1) ] + [self.root(self.n-1, self.n-1, 0, 0)]

    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_lattice().positive_roots()
            [(1, 1, 0),
             (1, 0, 1),
             (0, 1, 1),
             (1, -1, 0),
             (1, 0, -1),
             (0, 1, -1),
             (2, 0, 0),
             (0, 2, 0),
             (0, 0, 2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,0,p) for i in range(j)])
        res.extend([self.root(i,i,0,0) for i in range(self.n)])
        return res

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_lattice().negative_roots()
            [(-1, 1, 0),
             (-1, 0, 1),
             (0, -1, 1),
             (-1, -1, 0),
             (-1, 0, -1),
             (0, -1, -1),
             (-2, 0, 0),
             (0, -2, 0),
             (0, 0, -2)]
        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend( [self.root(i,j,1,p) for i in range(j) ] )
        res.extend( [ self.root(i,i,1,1) for i in range(self.n) ] )
        return res


    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: RootSystem(['C',3]).ambient_lattice().fundamental_weights()
            [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
        """
        return [ sum(self._term(j) for j in range(i+1)) for i in range(self.n)]


class AmbientLattice_d(AmbientLattice_generic):
    def root(self, i, j, p1, p2):
        """
        Note that indexing starts at 0.

        EXAMPLES:
            sage: e = RootSystem(['D',3]).ambient_lattice()
            sage: e.root(0, 1, 1, 1)
            (-1, -1, 0)
            sage: e.root(0, 0, 1, 1)
            (-1, 0, 0)

        """
        if i != j:
            return (-1)**p1 * self._term(i) + (-1)**p2 * self._term(j)
        else:
            return (-1)**p1 * self._term(i)

    def simple_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_lattice().simple_roots()
            [(1, -1, 0, 0), (0, 1, -1, 0), (0, 0, 1, -1), (0, 0, 1, 1)]

        """
        return [ self.root(i, i+1, 0, 1) for i in range(self.n-1) ] + [self.root(self.n-2, self.n-1, 0, 0)]

    def positive_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_lattice().positive_roots()
            [(1, 1, 0, 0),
             (1, 0, 1, 0),
             (0, 1, 1, 0),
             (1, 0, 0, 1),
             (0, 1, 0, 1),
             (0, 0, 1, 1),
             (1, -1, 0, 0),
             (1, 0, -1, 0),
             (0, 1, -1, 0),
             (1, 0, 0, -1),
             (0, 1, 0, -1),
             (0, 0, 1, -1)]

        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,0,p) for i in range(j)])
        return res

    def negative_roots(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_lattice().negative_roots()
            [(-1, 1, 0, 0),
             (-1, 0, 1, 0),
             (0, -1, 1, 0),
             (-1, 0, 0, 1),
             (0, -1, 0, 1),
             (0, 0, -1, 1),
             (-1, -1, 0, 0),
             (-1, 0, -1, 0),
             (0, -1, -1, 0),
             (-1, 0, 0, -1),
             (0, -1, 0, -1),
             (0, 0, -1, -1)]

        """
        res = []
        for p in [0,1]:
            for j in range(self.n):
                res.extend([self.root(i,j,1,p) for i in range(j)])
        return res


    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: RootSystem(['D',4]).ambient_lattice().fundamental_weights()
            [(1, 0, 0, 0), (1, 1, 0, 0), (1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, 1/2)]

        """
        return [ sum(self._term(j) for j in range(i+1)) for i in range(self.n-2)]+\
               [ sum(self._term(j) for j in range(self.n-1))/2-self._term(self.n-1)/2]+\
               [ sum(self._term(j) for j in range(self.n))/2 ]


class AmbientLattice_e(AmbientLattice_generic):
    """
    The lattice behind E6, E7, or E8.  The computations are based on Bourbaki,
    Groupes et Algebres de Lie, Ch. 4,5,6 (planche V-VII).
    """
    def __init__(self, ct):
        """
        Create the ambient lattice for the root system for E6, E7, E8.
        Specify the Base, i.e., the simple roots w.r. to the canonical
        basis for R^8.

        EXAMPLES:
            sage: e = RootSystem(['E',6]).ambient_lattice()
            sage: e == loads(dumps(e))
            True

        """
        v = ZZ(1)/ZZ(2)
        self.n = 8          # We're always in R^8, but not always the whole space.
        AmbientLattice_generic.__init__(self, ct)

        # Note that the lattices for the root systems E6, E7 have dimensions (ranks) 5, 6;
        #  while that for E8 has dimension 8.
        if ct.n == 6:
            self.dim = 5
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1)]
        elif ct.n == 7:
            self.dim = 6
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1),
                         self.root(4,5,p1=1)]
        elif ct.n == 8:
            self.dim = 8
            self.Base = [v*(self.root(0,7)-self.root(1,2,3,4,5,6)),
                         self.root(0,1),
                         self.root(0,1,p1=1),
                         self.root(1,2,p1=1),
                         self.root(2,3,p1=1),
                         self.root(3,4,p1=1),
                         self.root(4,5,p1=1),
                         self.root(5,6,p1=1)]
        else:
            raise NotImplementedError, "Type \'E\' root systems only come in flavors 6, 7, 8.  Please make another choice"

    def root(self, i1, i2=None, i3=None, i4=None, i5=None, i6=None, i7=None, i8=None, p1=0, p2=0, p3=0, p4=0, p5=0, p6=0, p7=0, p8=0):
        """
        Compute an element of the underlying lattice, using the specified elements of
        the standard basis, with signs dictated by the corresponding 'pi' arguments.
        We rely on the caller to provide the correct arguments.
        This is typically used to generate roots, although the generated elements
        need not be roots themselves.
        We assume that if one of the indices is not given, the rest are not as well.
        This should work for E6, E7, E8.

        EXAMPLES:
            sage: E6 = RootSystem(['E',6])
            sage: LE6 = E6.ambient_lattice()
            sage: [ LE6.root(i,j,p3=1) for i in xrange(LE6.n) for j in xrange(i+1,LE6.n) ]
            [(1, 1, 0, 0, 0, 0, 0, 0), (1, 0, 1, 0, 0, 0, 0, 0), (1, 0, 0, 1, 0, 0, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0), (1, 0, 0, 0, 0, 1, 0, 0), (1, 0, 0, 0, 0, 0, 1, 0), (1, 0, 0, 0, 0, 0, 0, 1), (0, 1, 1, 0, 0, 0, 0, 0), (0, 1, 0, 1, 0, 0, 0, 0), (0, 1, 0, 0, 1, 0, 0, 0), (0, 1, 0, 0, 0, 1, 0, 0), (0, 1, 0, 0, 0, 0, 1, 0), (0, 1, 0, 0, 0, 0, 0, 1), (0, 0, 1, 1, 0, 0, 0, 0), (0, 0, 1, 0, 1, 0, 0, 0), (0, 0, 1, 0, 0, 1, 0, 0), (0, 0, 1, 0, 0, 0, 1, 0), (0, 0, 1, 0, 0, 0, 0, 1), (0, 0, 0, 1, 1, 0, 0, 0), (0, 0, 0, 1, 0, 1, 0, 0), (0, 0, 0, 1, 0, 0, 1, 0), (0, 0, 0, 1, 0, 0, 0, 1), (0, 0, 0, 0, 1, 1, 0, 0), (0, 0, 0, 0, 1, 0, 1, 0), (0, 0, 0, 0, 1, 0, 0, 1), (0, 0, 0, 0, 0, 1, 1, 0), (0, 0, 0, 0, 0, 1, 0, 1), (0, 0, 0, 0, 0, 0, 1, 1)]
        """
        if i1 == i2 or i2 == None:
            return (-1)**p1*self._term(i1)
        if i3 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)
        if i4 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)
        if i5 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)+(-1)**p4*self._term(i4)
        if i6 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)+(-1)**p4*self._term(i4)+(-1)**p5*self._term(i5)
        if i7 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)+(-1)**p4*self._term(i4)+(-1)**p5*self._term(i5)+(-1)**p6*self._term(i6)
        if i8 == None:
            return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)+(-1)**p4*self._term(i4)+(-1)**p5*self._term(i5)+(-1)**p6*self._term(i6)+(-1)**p7*self._term(i7)
        return (-1)**p1*self._term(i1) + (-1)**p2*self._term(i2)+(-1)**p3*self._term(i3)+(-1)**p4*self._term(i4)+(-1)**p5*self._term(i5)+(-1)**p6*self._term(i6)+(-1)**p7*self._term(i7)+(-1)**p8*self._term(i8)

    def simple_roots(self):
        """
        There are computed as what Bourbaki calls the Base:
            a1 = e2-e3, a2 = e3-e4, a3 = e4, a4 = 1/2*(e1-e2-e3-e4)
        EXAMPLES:
            sage: LE6 = RootSystem(['E',6]).ambient_lattice()
            sage: LE6.simple_roots()
            [(1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (1, 1, 0, 0, 0, 0, 0, 0), (-1, 1, 0, 0, 0, 0, 0, 0), (0, -1, 1, 0, 0, 0, 0, 0), (0, 0, -1, 1, 0, 0, 0, 0), (0, 0, 0, -1, 1, 0, 0, 0)]
        """
        return self.Base

    def negative_roots(self):
        """
        The negative postive roots.
        EXAMPLES:
            sage: LE6 =  RootSystem(['E',6]).ambient_lattice()
            sage: LE6.negative_roots()
            [(-1, -1, 0, 0, 0, 0, 0, 0), (-1, 0, -1, 0, 0, 0, 0, 0), (-1, 0, 0, -1, 0, 0, 0, 0), (-1, 0, 0, 0, -1, 0, 0, 0), (0, -1, -1, 0, 0, 0, 0, 0), (0, -1, 0, -1, 0, 0, 0, 0), (0, -1, 0, 0, -1, 0, 0, 0), (0, 0, -1, -1, 0, 0, 0, 0), (0, 0, -1, 0, -1, 0, 0, 0), (0, 0, 0, -1, -1, 0, 0, 0), (1, -1, 0, 0, 0, 0, 0, 0), (1, 0, -1, 0, 0, 0, 0, 0), (1, 0, 0, -1, 0, 0, 0, 0), (1, 0, 0, 0, -1, 0, 0, 0), (0, 1, -1, 0, 0, 0, 0, 0), (0, 1, 0, -1, 0, 0, 0, 0), (0, 1, 0, 0, -1, 0, 0, 0), (0, 0, 1, -1, 0, 0, 0, 0), (0, 0, 1, 0, -1, 0, 0, 0), (0, 0, 0, 1, -1, 0, 0, 0), (-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (1/2, -1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2)]
        """
        return [ -a for a in self.positive_roots()]

    def positive_roots(self):
        """
        These are the roots positive w.r. to lexicographic ordering of the
        basis elements (e1<...<e4).
        EXAMPLES:
            sage: LE6 =  RootSystem(['E',6]).ambient_lattice()
            sage: LE6.positive_roots()
            [(1, 1, 0, 0, 0, 0, 0, 0), (1, 0, 1, 0, 0, 0, 0, 0), (1, 0, 0, 1, 0, 0, 0, 0), (1, 0, 0, 0, 1, 0, 0, 0), (0, 1, 1, 0, 0, 0, 0, 0), (0, 1, 0, 1, 0, 0, 0, 0), (0, 1, 0, 0, 1, 0, 0, 0), (0, 0, 1, 1, 0, 0, 0, 0), (0, 0, 1, 0, 1, 0, 0, 0), (0, 0, 0, 1, 1, 0, 0, 0), (-1, 1, 0, 0, 0, 0, 0, 0), (-1, 0, 1, 0, 0, 0, 0, 0), (-1, 0, 0, 1, 0, 0, 0, 0), (-1, 0, 0, 0, 1, 0, 0, 0), (0, -1, 1, 0, 0, 0, 0, 0), (0, -1, 0, 1, 0, 0, 0, 0), (0, -1, 0, 0, 1, 0, 0, 0), (0, 0, -1, 1, 0, 0, 0, 0), (0, 0, -1, 0, 1, 0, 0, 0), (0, 0, 0, -1, 1, 0, 0, 0), (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2), (1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2), (1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2), (1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2), (1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2), (1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2), (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (-1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2), (-1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2), (-1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2), (-1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2)]
            sage: LE6.rho()
            (0, 1, 2, 3, 4, -4, -4, 4)
            sage: E8=RootSystem(['E',8])
            sage: LE8=E8.ambient_lattice()
            sage: LE8.negative_roots()
            [(-1, -1, 0, 0, 0, 0, 0, 0), (-1, 0, -1, 0, 0, 0, 0, 0), (-1, 0, 0, -1, 0, 0, 0, 0), (-1, 0, 0, 0, -1, 0, 0, 0), (-1, 0, 0, 0, 0, -1, 0, 0), (-1, 0, 0, 0, 0, 0, -1, 0), (-1, 0, 0, 0, 0, 0, 0, -1), (0, -1, -1, 0, 0, 0, 0, 0), (0, -1, 0, -1, 0, 0, 0, 0), (0, -1, 0, 0, -1, 0, 0, 0), (0, -1, 0, 0, 0, -1, 0, 0), (0, -1, 0, 0, 0, 0, -1, 0), (0, -1, 0, 0, 0, 0, 0, -1), (0, 0, -1, -1, 0, 0, 0, 0), (0, 0, -1, 0, -1, 0, 0, 0), (0, 0, -1, 0, 0, -1, 0, 0), (0, 0, -1, 0, 0, 0, -1, 0), (0, 0, -1, 0, 0, 0, 0, -1), (0, 0, 0, -1, -1, 0, 0, 0), (0, 0, 0, -1, 0, -1, 0, 0), (0, 0, 0, -1, 0, 0, -1, 0), (0, 0, 0, -1, 0, 0, 0, -1), (0, 0, 0, 0, -1, -1, 0, 0), (0, 0, 0, 0, -1, 0, -1, 0), (0, 0, 0, 0, -1, 0, 0, -1), (0, 0, 0, 0, 0, -1, -1, 0), (0, 0, 0, 0, 0, -1, 0, -1), (0, 0, 0, 0, 0, 0, -1, -1), (1, -1, 0, 0, 0, 0, 0, 0), (1, 0, -1, 0, 0, 0, 0, 0), (1, 0, 0, -1, 0, 0, 0, 0), (1, 0, 0, 0, -1, 0, 0, 0), (1, 0, 0, 0, 0, -1, 0, 0), (1, 0, 0, 0, 0, 0, -1, 0), (1, 0, 0, 0, 0, 0, 0, -1), (0, 1, -1, 0, 0, 0, 0, 0), (0, 1, 0, -1, 0, 0, 0, 0), (0, 1, 0, 0, -1, 0, 0, 0), (0, 1, 0, 0, 0, -1, 0, 0), (0, 1, 0, 0, 0, 0, -1, 0), (0, 1, 0, 0, 0, 0, 0, -1), (0, 0, 1, -1, 0, 0, 0, 0), (0, 0, 1, 0, -1, 0, 0, 0), (0, 0, 1, 0, 0, -1, 0, 0), (0, 0, 1, 0, 0, 0, -1, 0), (0, 0, 1, 0, 0, 0, 0, -1), (0, 0, 0, 1, -1, 0, 0, 0), (0, 0, 0, 1, 0, -1, 0, 0), (0, 0, 0, 1, 0, 0, -1, 0), (0, 0, 0, 1, 0, 0, 0, -1), (0, 0, 0, 0, 1, -1, 0, 0), (0, 0, 0, 0, 1, 0, -1, 0), (0, 0, 0, 0, 1, 0, 0, -1), (0, 0, 0, 0, 0, 1, -1, 0), (0, 0, 0, 0, 0, 1, 0, -1), (0, 0, 0, 0, 0, 0, 1, -1), (-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2), (-1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2), (-1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2), (-1/2, 1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2), (-1/2, 1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2), (1/2, -1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2), (1/2, -1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2), (1/2, -1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2), (1/2, -1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2), (1/2, -1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2), (1/2, -1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (1/2, -1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2), (1/2, -1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2), (1/2, -1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2), (1/2, -1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2), (1/2, -1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2), (1/2, -1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, -1/2, -1/2, -1/2, -1/2, -1/2), (1/2, 1/2, -1/2, -1/2, -1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, -1/2, 1/2, -1/2, 1/2, -1/2), (1/2, 1/2, -1/2, -1/2, 1/2, 1/2, -1/2, -1/2), (1/2, 1/2, -1/2, 1/2, -1/2, -1/2, 1/2, -1/2), (1/2, 1/2, -1/2, 1/2, -1/2, 1/2, -1/2, -1/2), (1/2, 1/2, -1/2, 1/2, 1/2, -1/2, -1/2, -1/2), (1/2, 1/2, -1/2, 1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, -1/2, -1/2, -1/2, 1/2, -1/2), (1/2, 1/2, 1/2, -1/2, -1/2, 1/2, -1/2, -1/2), (1/2, 1/2, 1/2, -1/2, 1/2, -1/2, -1/2, -1/2), (1/2, 1/2, 1/2, -1/2, 1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, 1/2, -1/2, -1/2, -1/2, -1/2), (1/2, 1/2, 1/2, 1/2, -1/2, 1/2, 1/2, -1/2), (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, 1/2, -1/2), (1/2, 1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2)]
            sage: LE8.rho()
            (0, 1, 2, 3, 4, 5, 6, 23)
        """
        v = ZZ(1)/ZZ(2)
        # Note that
        if not hasattr(self, 'PosRoots'):
            if self.dim == 5:
                self.PosRoots = ( [ self.root(i,j) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ self.root(i,j,p1=1) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ v*(self.root(7)-self.root(6)-self.root(5)+self.root(0,1,2,3,4,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] if (p1+p2+p3+p4+p5)%2 == 0 ])
            elif self.dim == 6:
                self.PosRoots = ( [ self.root(i,j) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ self.root(i,j,p1=1) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ self.root(6,7,p1=1) ] +
                                  [ v*(self.root(7)-self.root(6)+self.root(0,1,2,3,4,5,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] for p6 in [0,1] if (p1+p2+p3+p4+p5+p6)%2 == 1 ])
            elif self.dim == 8:
                self.PosRoots = ( [ self.root(i,j) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ self.root(i,j,p1=1) for i in xrange(self.dim) for j in xrange(i+1,self.dim) ] +
                                  [ v*(self.root(7)+self.root(0,1,2,3,4,5,6,p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7))
                                    for p1 in [0,1] for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] for p5 in [0,1] for p6 in [0,1] for p7 in [0,1] if (p1+p2+p3+p4+p5+p6+p7)%2 == 0 ])

        return self.PosRoots

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: LE6 = RootSystem(['E',6]).ambient_lattice()
            sage: LE6.fundamental_weights()
            [(0, 0, 0, 0, 0, -2/3, -2/3, 2/3), (1/2, 1/2, 1/2, 1/2, 1/2, -1/2, -1/2, 1/2), (-1/2, 1/2, 1/2, 1/2, 1/2, -5/6, -5/6, 5/6), (0, 0, 1, 1, 1, -1, -1, 1), (0, 0, 0, 1, 1, -2/3, -2/3, 2/3), (0, 0, 0, 0, 1, -1/3, -1/3, 1/3)]
        """
        v2 = ZZ(1)/ZZ(2)
        v3 = ZZ(1)/ZZ(3)
        if self.dim == 5:
            return [ 2*v3*self.root(7,6,5,p2=1,p3=1),
                     v2*self.root(0,1,2,3,4,5,6,7,p6=1,p7=1),
                     5*v2*v3*self.root(7,6,5,p2=1,p3=1)+v2*self.root(0,1,2,3,4,p1=1),
                     self.root(2,3,4,5,6,7,p4=1,p5=1),
                     2*v3*self.root(7,6,5,p2=1,p3=1)+self.root(3,4),
                     v3*self.root(7,6,5,p2=1,p3=1)+self.root(4)]
        elif self.dim == 6:
            return [ self.root(7,6,p2=1),
                     v2*self.root(0,1,2,3,4,5)+self.root(6,7,p1=1),
                     v2*(self.root(0,1,2,3,4,5,p1=1)+3*self.root(6,7,p1=1)),
                     self.root(2,3,4,5)+2*self.root(6,7,p1=1),
                     3*v2*self.root(6,7,p1=1)+self.root(3,4,5),
                     self.root(4,5,6,7,p3=1),
                     self.root(5)+v2*self.root(6,7,p1=1)]
        elif self.dim == 8:
            return [ 2*self.root(7),
                     v2*(self.root(0,1,2,3,4,5,6)+5*self.root(7)),
                     v2*(self.root(0,1,2,3,4,5,6,p1=1)+7*self.root(7)),
                     self.root(2,3,4,5,6)+5*self.root(7),
                     self.root(3,4,5,6)+4*self.root(7),
                     self.root(4,5,6)+3*self.root(7),
                     self.root(5,6)+2*self.root(7),
                     self.root(6,7)]


class AmbientLattice_f(AmbientLattice_generic):
    """
    The lattice behind F4.  The computations are based on Bourbaki, Groupes et Algebres de Lie,
    Ch. 4,5,6 (planche VIII).
    """
    def __init__(self, ct):
        """
        Create the ambient lattice for the root system for F4.
        Specify the Base, i.e., the simple roots w.r. to the canonical
        basis for R^4.

        EXAMPLES:
            sage: e = RootSystem(['F',4]).ambient_lattice()
            sage: e == loads(dumps(e))
            True
        """
        v = ZZ(1)/ZZ(2)
        AmbientLattice_generic.__init__(self, ct)
        self.Base = [self.root(1,2,p2=1),
                     self.root(2,3,p2=1),
                     self.root(3),
                     v*(self.root(0)-self.root(1)-self.root(2)-self.root(3))]

    def root(self, i, j=None, k=None, l=None, p1=0, p2=0, p3=0, p4=0):
        """
        Compute a root from base elements of the underlying lattice.
        The arguments specify the basis elements and the signs.
        Sadly, the base elements are indexed zero-based.
        We assume that if one of the indices is not given, the rest are not as well.

        EXAMPLES:
            sage: F4 = RootSystem(['F',4])
            sage: LF4 = F4.ambient_lattice()
            sage: [ LF4.root(i,j,p2=1) for i in xrange(LF4.n) for j in xrange(i+1,LF4.n) ]
            [(1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1), (0, 1, -1, 0), (0, 1, 0, -1), (0, 0, 1, -1)]
        """
        if i == j or j == None:
            return (-1)**p1*self._term(i)
        if k == None:
            return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)
        if l == None:
            return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)+(-1)**p3*self._term(k)
        return (-1)**p1*self._term(i) + (-1)**p2*self._term(j)+(-1)**p3*self._term(k)+(-1)**p4*self._term(l)

    def simple_roots(self):
        """
        There are computed as what Bourbaki calls the Base:
            a1 = e2-e3, a2 = e3-e4, a3 = e4, a4 = 1/2*(e1-e2-e3-e4)

        EXAMPLES:
            sage: LF4 = RootSystem(['F',4]).ambient_lattice()
            sage: LF4.simple_roots()
            [(0, 1, -1, 0), (0, 0, 1, -1), (0, 0, 0, 1), (1/2, -1/2, -1/2, -1/2)]
        """
        return self.Base

    def negative_roots(self):
        """
        Returns the negative roots in self.

        EXAMPLES:
            sage: LF4 =  RootSystem(['F',4]).ambient_lattice()
            sage: LF4.negative_roots()
            [(-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (-1, -1, 0, 0), (-1, 0, -1, 0), (-1, 0, 0, -1), (0, -1, -1, 0), (0, -1, 0, -1), (0, 0, -1, -1), (-1, 1, 0, 0), (-1, 0, 1, 0), (-1, 0, 0, 1), (0, -1, 1, 0), (0, -1, 0, 1), (0, 0, -1, 1), (-1/2, -1/2, -1/2, -1/2), (-1/2, -1/2, -1/2, 1/2), (-1/2, -1/2, 1/2, -1/2), (-1/2, -1/2, 1/2, 1/2), (-1/2, 1/2, -1/2, -1/2), (-1/2, 1/2, -1/2, 1/2), (-1/2, 1/2, 1/2, -1/2), (-1/2, 1/2, 1/2, 1/2)]
        """
        return [ -a for a in self.positive_roots()]

    def positive_roots(self):
        """
        These are the roots positive w.r. to lexicographic ordering of the
        basis elements (e1<...<e4).

        EXAMPLES:
            sage: LF4 = RootSystem(['F',4]).ambient_lattice()
            sage: LF4.positive_roots()
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 1), (0, 1, 1, 0), (0, 1, 0, 1), (0, 0, 1, 1), (1, -1, 0, 0), (1, 0, -1, 0), (1, 0, 0, -1), (0, 1, -1, 0), (0, 1, 0, -1), (0, 0, 1, -1), (1/2, 1/2, 1/2, 1/2), (1/2, 1/2, 1/2, -1/2), (1/2, 1/2, -1/2, 1/2), (1/2, 1/2, -1/2, -1/2), (1/2, -1/2, 1/2, 1/2), (1/2, -1/2, 1/2, -1/2), (1/2, -1/2, -1/2, 1/2), (1/2, -1/2, -1/2, -1/2)]
            sage: LF4.rho()
            (11/2, 5/2, 3/2, 1/2)
        """
        v = ZZ(1)/ZZ(2)
        if not hasattr(self, 'PosRoots'):
            self.PosRoots = ([ self._term(i) for i in xrange(self.n) ] +
                            [ self.root(i,j,p2=0) for i in xrange(self.n) for j in xrange(i+1,self.n) ] +
                            [ self.root(i,j,p2=1) for i in xrange(self.n) for j in xrange(i+1,self.n) ] +
                            [ v*self.root(0,1,2,3,0,p2,p3,p4) for p2 in [0,1] for p3 in [0,1] for p4 in [0,1] ])
        return self.PosRoots

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: LF4 =  RootSystem(['F',4]).ambient_lattice()
            sage: LF4.fundamental_weights()
            [(1, 1, 0, 0), (2, 1, 1, 0), (3/2, 1/2, 1/2, 1/2), (1, 0, 0, 0)]
        """
        v = ZZ(1)/ZZ(2)
        return [ self._term(0)+self._term(1), 2*self._term(0)+self._term(1)+self._term(2), v*(3*self._term(0)+self._term(1)+self._term(2)+self._term(3)), self._term(0)]


class AmbientLattice_g(AmbientLattice_generic):
    """
    TESTS:
        sage: [WeylDim(['G',2],[a,b]) for a,b in [[0,0], [1,0], [0,1], [1,1]]]
        [1, 7, 14, 64]
    """
    def __init__(self, ct):
        """
        EXAMPLES:
            sage: e = RootSystem(['G',2]).ambient_lattice()
            sage: e == loads(dumps(e))
            True
        """
        self.n = 3
        AmbientLattice_generic.__init__(self, ct)

    def simple_roots(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_lattice().simple_roots()
            [(0, 1, -1), (1, -2, 1)]
         """
        return [ self._term(1)-self._term(2),\
                 self._term(0)-2*self._term(1)+self._term(2)]

    def positive_roots(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_lattice().positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return [ c0*self._term(0)+c1*self._term(1)+c2*self._term(2) \
                 for [c0,c1,c2] in
                 [[0,1,-1],[1,-2,1],[1,-1,0],[1,0,-1],[1,1,-2],[2,-1,-1]]]

    def negative_roots(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_lattice().negative_roots()
            [(0, -1, 1), (-1, 2, -1), (-1, 1, 0), (-1, 0, 1), (-1, -1, 2), (-2, 1, 1)]
        """
        return [ c0*self._term(0)+c1*self._term(1)+c2*self._term(2) \
                 for [c0,c1,c2] in
                 [[0,-1,1],[-1,2,-1],[-1,1,0],[-1,0,1],[-1,-1,2],[-2,1,1]]]

    def fundamental_weights(self):
        """
        EXAMPLES:
            sage: CartanType(['G',2]).root_system().ambient_lattice().fundamental_weights()
            [(1, 0, -1), (2, -1, -1)]
        """
        return [ c0*self._term(0)+c1*self._term(1)+c2*self._term(2) \
                 for [c0,c1,c2] in
                 [[1,0,-1],[2,-1,-1]]]


def WeylDim(ct, coeffs):
    """
    The Weyl Dimension Formula. Here type is a Cartan type and coeffs
    are a list of nonnegative integers of length equal to the rank
    type[1]. A dominant weight hwv is constructed by summing the
    fundamental weights with coefficients from this list. The
    dimension of the irreducible representation of the semisimple
    complex Lie algebra with highest weight vector hwv is returned.

    EXAMPLES:
    For SO(7), the Cartan type is B3, so:
        sage: WeylDim(['B',3],[1,0,0]) # standard representation of SO(7)
        7
        sage: WeylDim(['B',3],[0,1,0]) # exterior square
        21
        sage: WeylDim(['B',3],[0,0,1]) # spin representation of spin(7)
        8
        sage: WeylDim(['B',3],[1,0,1]) # sum of the first and third fundamental weights
        48
        sage: [WeylDim(['F',4],x) for x in [1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        [52, 1274, 273, 26]
        sage: [WeylDim(['E', 6], x) for x in [0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 2], [0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0], [1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 1], [2, 0, 0, 0, 0, 0]]
        [1, 78, 27, 351, 351, 351, 27, 650, 351]
    """
    lattice = RootSystem(ct).ambient_lattice()
    rank = ct[1]
    fw = lattice.fundamental_weights()
    hwv = sum(coeffs[i]*fw[i] for i in range(min(rank, len(coeffs))))
    return lattice.weyl_dimension(hwv)
