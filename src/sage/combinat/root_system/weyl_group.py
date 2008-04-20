"""
Weyl Groups
"""
#*****************************************************************************
#       Copyright (C) 2008 Daniel Bump <bump at match.stanford.edu>,
#                          Mike Hansen <mhansen@gmail.com>
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
from sage.misc.functional import transpose
from sage.groups.matrix_gps.matrix_group import MatrixGroup_gens
from sage.groups.matrix_gps.matrix_group_element import MatrixGroupElement
from sage.rings.all import ZZ, QQ
from sage.interfaces.gap import gap
from sage.misc.cache import Cache
from sage.combinat.root_system.cartan_type import CartanType
from sage.matrix.constructor import matrix, diagonal_matrix

def WeylGroup(ct):
    """
    Returns the Weyl group of type type.

    INPUT:
        ct - a Cartan Type.

    EXAMPLES:
      The following constructions yield the same result, namely
      a weight lattice and its corresponding Weyl group:

        sage: G = WeylGroup(['F',4])
        sage: L = G.lattice()

      or alternatively and equivalently:

        sage: L = RootSystem(['F',4]).ambient_lattice()
        sage: G = L.weyl_group()

      Either produces a weight lattice, with access to its
      roots and weights.

        sage: G = WeylGroup(['F',4])
        sage: G.order()
        1152
        sage: [s1,s2,s3,s4] = G.simple_reflections()
        sage: w = s1*s2*s3*s4; w
        [ 1/2  1/2  1/2  1/2]
        [-1/2  1/2  1/2 -1/2]
        [ 1/2  1/2 -1/2 -1/2]
        [ 1/2 -1/2  1/2 -1/2]
        sage: type(w)
        <class 'sage.combinat.root_system.weyl_group.WeylGroupElement'>
        sage: w.order()
        12
        sage: w.length() # length function on Weyl group
        4

        sage: L = G.lattice()
        sage: fw = L.fundamental_weights(); fw
        [(1, 1, 0, 0), (2, 1, 1, 0), (3/2, 1/2, 1/2, 1/2), (1, 0, 0, 0)]
        sage: rho = sum(fw); rho
        (11/2, 5/2, 3/2, 1/2)
        sage: w.action(rho) # action of G on weight lattice
        (5, -1, 3, 2)
    """
    ct = CartanType(ct)
    return weyl_group_cache(ct)

def _WeylGroup(ct):
    """
    Returns the Weyl group of type ct.  This function is wrapped by a cache
    object so that the Weyl group only gets computed once for each type.

    EXAMPLES:
        sage: from sage.combinat.root_system.weyl_group import _WeylGroup
        sage: _WeylGroup(CartanType(['A',2]))
        The Weyl Group of type ['A', 2]
    """
    from sage.combinat.root_system.root_system import RootSystem
    e = RootSystem(ct).ambient_lattice()
    n = e.n
    rank = e.ct[1]
    basis = e._free_module.basis()
    gens = []
    for k in range(rank):
        m = []
        for i in range(n):
            for j in range(n):
                m.append(e.simple_reflection(k+1,basis[j]).inner_product(basis[i]))
        gens.append(matrix(QQ,n,m))
    return WeylGroup_gens(gens,e)

weyl_group_cache = Cache(_WeylGroup)


class WeylGroup_gens(MatrixGroup_gens):
    """Class for a Weyl Group with generators (simple reflections)"""
    def __init__(self, gens, L):
        """
        EXAMPLES:
            sage: G = WeylGroup(['F',4])
            sage: G == loads(dumps(G))
            True
        """
        MatrixGroup_gens.__init__(self, gens)
        self.ambient_lattice = L
        self.reflections = [WeylGroupElement(g, self) for g in gens]
        self.n = L.n

    def __repr__(self):
        """
        EXAMPLES:
            sage: WeylGroup(['A',1]).__repr__()
            "The Weyl Group of type ['A', 1]"
        """
        return "The Weyl Group of type %s"%repr(self.ambient_lattice.ct)

    def __call__(self, x):
        if isinstance(x, WeylGroupElement) and x.parent() is self:
            return x
        M = self.matrix_space()(x)
        g = WeylGroupElement(M, self)
        if not gap(g) in gap(self):
            raise TypeError, "no way to coerce element to self."
        return g

    def list(self):
        """
        Returns a list of the elements of self.

        EXAMPLES:
            sage: G = WeylGroup(['A',1])
            sage: G.list()
            [[0 1]
             [1 0], [1 0]
             [0 1]]

        """
        return [WeylGroupElement(a._matrix_(QQ),self) for a in self._gap_().Elements()]

    def character_table(self):
        """Returns the GAP character table as a string. For larger tables
        you may preface this with a command such as
        # sage: gap.eval("SizeScreen([120,40])")
        in order to widen the screen. Then try
        # print WeylGroup(['F',4]).character_table()
        """
        return gap.eval("Display(CharacterTable(%s))"%gap(self).name())

    def unit(self):
        """Returns the unit element of the Weyl group
        EXAMPLE:
          sage: e = WeylGroup(['A',3]).unit(); e
          [1 0 0 0]
          [0 1 0 0]
          [0 0 1 0]
          [0 0 0 1]
          sage: type(e)
          <class 'sage.combinat.root_system.weyl_group.WeylGroupElement'>
        """
        return WeylGroupElement(matrix(QQ,self.n,self.n,1),self)

    def lattice(self):
        """
        Returns the ambient lattice of self's type.

        EXAMPLES:
            sage: G = WeylGroup(['F',4])
            sage: G.lattice()
            Ambient lattice of the root system of type ['F', 4]
        """
        return self.ambient_lattice

    def simple_reflections(self):
        """
        Returns the list of the simple reflections of self.

        EXAMPLES:
            sage: G = WeylGroup(['F',4])
            sage: [s1,s2,s3,s4] = G.simple_reflections()
            sage: w = s1*s2*s3*s4; w
            [ 1/2  1/2  1/2  1/2]
            [-1/2  1/2  1/2 -1/2]
            [ 1/2  1/2 -1/2 -1/2]
            [ 1/2 -1/2  1/2 -1/2]
            sage: s4^2==G.unit()
            True
            sage: type(w)
            <class 'sage.combinat.root_system.weyl_group.WeylGroupElement'>
        """
        return self.reflections

    def simple_reflection(self, i):
        """
        Returns the $i^th$ simple reflection.

        EXAMPLES:
            sage: G = WeylGroup(['F',4])
            sage: G.simple_reflection(1)
            [1 0 0 0]
            [0 0 1 0]
            [0 1 0 0]
            [0 0 0 1]
        """
        if i not in self.ambient_lattice.ct.index_set():
            raise ValueError, "i must be in the index set"
        return self.reflections[i-1]

    def long_element(self):
        """
        Returns the long Weyl group element.
        EXAMPLES:
           sage: for
        """
        type = self.ambient_lattice.ct
        if type[0] == 'D' and type[1]%2 == 1:
            l = [-1 for i in range(self.n-1)]
            l.append(1)
            m = diagonal_matrix(QQ,l)
        elif type[0] == 'A':
            l = [0 for k in range((self.n)**2)]
            for k in range(self.n-1, (self.n)**2-1, self.n-1):
                l[k] = 1
            m = matrix(QQ, self.n, l)
        elif type[0] == 'E':
            if type[1] == 6:
                half = ZZ(1)/ZZ(2)
                l = [[-half, -half, -half, half, 0, 0, 0, 0], \
                     [-half, -half, half, -half, 0, 0, 0, 0], \
                     [-half, half, -half, -half, 0, 0, 0, 0], \
                     [half, -half, -half, -half, 0, 0, 0, 0], \
                     [0, 0, 0, 0, half, half, half, -half],   \
                     [0, 0, 0, 0, half, half, -half, half],   \
                     [0, 0, 0, 0, half, -half, half, half],   \
                     [0, 0, 0, 0, -half, half, half, half]]
                m = matrix(QQ, 8, l)
            else:
                raise NotImplementedError, "Not implemented yet for this type"
        elif type[0] == 'G':
            third = ZZ(1)/ZZ(3)
            twothirds = ZZ(2)/ZZ(3)
            l = [[-third, twothirds, twothirds], \
                 [twothirds, -third, twothirds], \
                 [twothirds, twothirds, -third]]
            m = matrix(QQ, 3, l)
        else:
            m = diagonal_matrix([-1 for i in range(self.n)])
        print m
        return self.__call__(m)

    def cartan_type(self):
        """
        Returns the CartanType associated to self.

        EXAMPLES:
            sage: G = WeylGroup(['F',4])
            sage: G.cartan_type()
            ['F', 4]

        """
        return self.ambient_lattice.ct

    def __cmp__(self, other):
        """
        TESTS:
            sage: from sage.combinat.root_system.weyl_group import _WeylGroup
            sage: G1 = _WeylGroup(CartanType(['A',2]))
            sage: G2 = _WeylGroup(CartanType(['A',2]))
            sage: G1 == G2
            True
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self.ambient_lattice.ct != other.ambient_lattice.ct:
            return cmp(self.ambient_lattice.ct, other.ambient_lattice.ct)
        return 0

class WeylGroupElement(MatrixGroupElement):
    """Class for a Weyl Group elements"""
    def __init__(self, g, parent):
        """
        EXAMPLES:
            sage: G = WeylGroup(['A',2])
            sage: s1 = G.simple_reflection(1)
            sage: loads(dumps(s1))
            [0 1 0]
            [1 0 0]
            [0 0 1]
        """
        MatrixGroupElement.__init__(self, g, parent)
        self.__matrix = self._MatrixGroupElement__mat
        self._parent = parent

    def lattice(self):
        """
        Returns the ambient lattice associated with self.

        EXAMPLES:
            sage: w = WeylGroup(['A',2])
            sage: s1 = w.simple_reflection(1)
            sage: s1.lattice()
            Ambient lattice of the root system of type ['A', 2]
        """
        return self._parent.lattice()

    def length(self):
        """
        Returns the length of self, which is the smallest number
        of simple reflections into which the element can be
        decomposed.

        EXAMPLES:
            sage: w = WeylGroup(['A',3])
            sage: s1 = w.simple_reflection(1)
            sage: s2 = w.simple_reflection(2)
            sage: s1.length()
            1
            sage: (s1*s2).length()
            2

        """
        ret = 0
        for alph in self.lattice().positive_roots():
            walph = self.__matrix*alph
            for v in self.lattice().fundamental_weights():
                if walph.inner_product(v) < 0:
                    ret += 1
                    break
        return ret

    def matrix(self):
        """
        Returns self as a matrix.

        EXAMPLES:
            sage: w = WeylGroup(['A',2])
            sage: s1 = w.simple_reflection(1)
            sage: m = s1.matrix(); m
            [0 1 0]
            [1 0 0]
            [0 0 1]
            sage: m.parent()
            Full MatrixSpace of 3 by 3 dense matrices over Rational Field
        """
        return self.__matrix

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: w = WeylGroup(['A',3])
            sage: s1 = w.simple_reflection(1)
            sage: s2 = w.simple_reflection(2)
            sage: s1 == s1
            True
            sage: s1 == s2
            False
        """
        if self.__class__ != other.__class__:
            return cmp(self.__class__, other.__class__)
        if self._parent.cartan_type() != other._parent.cartan_type():
            return cmp(self._parent.cartan_type(), other._parent.cartan_type())
        return cmp(self.matrix(), other.matrix())

    def parent(self):
        """
        Returns self's parent.

        EXAMPLES:
            sage: w = WeylGroup(['A',2])
            sage: s1 = w.simple_reflection(1)
            sage: s1.parent()
            The Weyl Group of type ['A', 2]

        """
        return self._parent

    def _mul_(self, other):
        """
        EXAMPLES:
            sage: w = WeylGroup(['A',2])
            sage: s1 = w.simple_reflection(1)
            sage: s1*s1
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return WeylGroupElement(self.__matrix * other.__matrix, self._parent)

    def action(self, v):
        """
        Returns the action of self on the vector v.

        EXAMPLES:
            sage: w = WeylGroup(['A',2])
            sage: s1 = w.simple_reflection(1)
            sage: v = vector(QQ,[1,0,0])
            sage: s1.action(v)
            (0, 1, 0)
        """
        return (self.__matrix*transpose(v)).columns()[0]
