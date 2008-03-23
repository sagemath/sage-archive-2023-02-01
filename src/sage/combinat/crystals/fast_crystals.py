r"""
Fast Rank Two Crystals
"""
#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#                          Ben Brubaker   <brubaker at math.mit.edu>
#                          Daniel Bump    <bump at match.stanford.edu>
#                          Justin Walker  <justin at mac.com>
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
#****************************************************************************

from sage.rings.integer                 import Integer
from sage.structure.element             import Element
from sage.combinat.cartan_type          import CartanType
from sage.combinat.crystals.crystals    import Crystal, ClassicalCrystal, CrystalElement


class FastCrystal(ClassicalCrystal):
    """
    An alternative implementation of rank 2 crystals. The root
    operators are implemented in memory by table lookup.
    This means that in comparison with the Crystals of
    Tableaux, these crystals are slow to instantiate
    but faster for computation. Implemented for types
    A2, B2 and C2.

    Input: CartanType and a shape. The CartanType is
    ['A',2], ['B',2] or ['C',2]. The shape is of
    the form [l1,l2] where l1 and l2 are either integers
    or (in type B) half integers such that l1-l2 is
    integral. It is assumed that l1 >= l2 >= 0.  If l1
    and l2 are integers, this will produce the a crystal
    isomorphic to the one obtained by
    CrystalOfTableaux(type, shape=[l1,l2]).  Furthermore
    FastCrystal(['B', 2], l1+1/2, l2+1/2) produces a
    crystal isomorphic to the following crystal T:

    C = CrystalOfTableaux(['B',2], shape=[l1,l2])
    D = CrystalOfSpins(['B',2])
    T = TensorProductOfCrystals(C,D,C.list()[0],D.list()[0])

    The representation of elements is in term of the
    Berenstein-Zelevinsky-Littelmann strings [a1, a2, ...]
    described under metapost in crystals.py. Alternative
    representations may be obtained by the options format="dual_string"
    or format="simple". In the simple format, the element is
    represented by and integer, and in the dual_string format,
    it is represented by the Berenstein-Zelevinsky-Littelmann
    string, but the underlying decomposition of the long Weyl
    group element into simple reflections is changed.

    TESTS:
        sage: C = FastCrystal(['A',2],shape=[4,1])
        sage: C.count()
        24
        sage: C.check()
        True
        sage: C = FastCrystal(['B',2],shape=[4,1])
        sage: C.count()
        154
        sage: C.check()
        True
        sage: C = FastCrystal(['B',2],shape=[3/2,1/2])
        sage: C.count()
        16
        sage: C.check()
        True
        sage: C = FastCrystal(['C',2],shape=[2,1])
        sage: C.count()
        16
        sage: C = FastCrystal(['C',2],shape=[3,1])
        sage: C.count()
        35
        sage: C.check()
        True
    """
    def __init__(self, type, shape, format="string"):
        if len(shape) == 0:
            l1 = 0
        else:
            l1 = shape[0]
        if len(shape) < 2:
            l2 = 0
        else:
            l2 = shape[1]
        self.delpat = []
        self.gampat = []
        if not type[1] == 2:
            raise NotImplementedError
        if type[0] == 'A':
            for b in range(l2,-1,-1):
                for a in range(l1,l2-1,-1):
                    for c in range(a,b-1,-1):
                        a3 = l1-a
                        a2 = l1+l2-a-b
                        a1 = a-c
                        b1 =  max(a3,a2-a1);
                        b2 = a1+a3;
                        b3 = min(a2-a3,a1);
                        self.delpat.append([a1,a2,a3])
                        self.gampat.append([b1,b2,b3])
        elif type[0] == 'B' or type[0] == 'C':
            if type[0] == 'B':
                [m1, m2] = [l1+l2, l1-l2]
            else:
                [m1, m2] = [l1, l2]
            for b in range(m2,-1,-1):
                for a in range(m1,m2-1,-1):
                    for c in range(b,a+1):
                        for d in range(c,-1,-1):
                            a1 = c-d
                            a2 = m1+m2+c-a-2*b
                            a3 = m1+m2-a-b
                            a4 = m1-a
                            b1 = max(a4,2*a3-a2,a2-2*a1)
                            b2 = max(a3, a1+a4, a1+2*a3-a2)
                            b3 = min(a2, 2*a2-2*a3+a4, 2*a1+a4)
                            b4 = min(a1, a2-a3, a3-a4)
                            if type[0] == 'B':
                                self.delpat.append([a1,a2,a3,a4])
                                self.gampat.append([b1,b2,b3,b4])
                            else:
                                self.gampat.append([a1,a2,a3,a4])
                                self.delpat.append([b1,b2,b3,b4])
        else:
            raise NotImplementedError

        self.format = format
        self.size = len(self.delpat)
        self._rootoperators = []
        for i in range(self.size):
            target = [x for x in self.delpat[i]]
            target[0] = target[0]-1
            if target in self.delpat:
                e1 = self.delpat.index(target)
            else:
                e1 = None
            target = [x for x in self.delpat[i]]
            target[0] = target[0]+1
            if target in self.delpat:
                f1 = self.delpat.index(target)
            else:
                f1 = None

            target = [x for x in self.gampat[i]]
            target[0] = target[0]-1
            if target in self.gampat:
                e2 = self.gampat.index(target)
            else:
                e2 = None
            target = [x for x in self.gampat[i]]
            target[0] = target[0]+1
            if target in self.gampat:
                f2 = self.gampat.index(target)
            else:
                f2 = None
            self._rootoperators.append([e1,f1,e2,f2])

        self.cartanType = CartanType(type)
        if int(2*l1)%2 == 0:
            l1_str = "%d"%l1
            l2_str = "%d"%l2
        else:
            assert self.cartanType[0] == 'B' and int(2*l2)%2 == 1
            l1_str = "%d/2"%int(2*l1)
            l2_str = "%d/2"%int(2*l2)
        self._name = "The fast crystal for %s2 with shape [%s,%s]"%(type[0],l1_str,l2_str)
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(0)]
        self._list = ClassicalCrystal.list(self)
        self._digraph = ClassicalCrystal.digraph(self)
        self._digraph_closure = self.digraph().transitive_closure()

    def __call__(self, value):
        return FastCrystalElement(self, value, self.format)

    def list(self):
        return self._list

    def digraph(self):
        return self._digraph

    def cmp_elements(self, x,y):
        r"""
        Returns whether there is a path from x to y in the crystal graph
        """
        assert x.parent() == self and y.parent() == self;
        if   self._digraph_closure.has_edge(x,y):
            return -1;
        elif self._digraph_closure.has_edge(y,x):
            return 1;
        else:
            return 0;

class FastCrystalElement(CrystalElement):

    def __init__(self, parent, value, format):
        self._parent = parent
        self.value = value
        self.format = format

    def parent(self):
        return self._parent

    def __repr__(self):
        if self.format == "string":
            return str(self._parent.delpat[self.value])
        elif self.format == "dual_string":
            return str(self._parent.gampat[self.value])
        elif self.format == "simple":
            return str(self.value)
        else:
            raise NotImplementedError

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
               self.parent()  == other.parent()   and \
               self.value     == other.value

    def __cmp__(self, other):
        assert self.parent() == other.parent()
        return self.parent().cmp_elements(self, other)

    def e(self, i):
        assert i in self.index_set()
        if i == 1:
            r = self._parent._rootoperators[self.value][0]
        else:
            r = self._parent._rootoperators[self.value][2]
        if r == None:
            return None
        else:
            return self._parent(r)

    def f(self, i):
        assert i in self.index_set()
        if i == 1:
            r = self._parent._rootoperators[self.value][1]
        else:
            r = self._parent._rootoperators[self.value][3]
        if r == None:
            return None
        else:
            return self._parent(r)

