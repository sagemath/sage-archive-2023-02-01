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

    Input: Cartan type, l1, l2. If l1 and l2 are integers,
    this will produce the a crystal isomorphic to the one
    obtained by CrystalOfTableaux(type, shape=[l1,l2]).
    For type ['B',2], l1 and l2 may be half integers but
    l1-l2 is an integer. If l1 and l2 are integers,
    FastCrystal(['B', 2], l1+1/2, l2+1/2) produces a crystal
    isomorphic to the following crystal T:
    C = CrystalOfTableaux(['B',2], shape=[l1,l2])
    D = CrystalOfSpins(['B',2])
    T = TensorProductOfCrystals(C,D,C.list()[0],D.list()[0])

    TESTS:
        sage: C = FastCrystal(['A',2],4,1)
        sage: C.count()
        24
        sage: C.check()
        True
        sage: C = FastCrystal(['B',2],4,1)
        sage: C.count()
        154
        sage: C.check()
        True
        sage: C = FastCrystal(['B',2],3/2,1/2)
        sage: C.count()
        16
        sage: C.check()
        True
        sage: C = FastCrystal(['C',2],2,1)
        sage: C.count()
        16
        sage: C = FastCrystal(['C',2],3,1)
        sage: C.count()
        35
        sage: C.check()
        True
    """
    def __init__(self, type, l1, l2):
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
        self._name = "The fast crystal for %s2 with shape [%d,%d]"%(type[0],l1,l2)
        self.index_set = self.cartanType.index_set()
        self.module_generators = [self(0)]
        self._list = ClassicalCrystal.list(self)
        self._digraph = ClassicalCrystal.digraph(self)
        self._digraph_closure = self.digraph().transitive_closure()

    def __call__(self, value):
        return FastCrystalElement(self, value)

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

    def __init__(self, parent, value):
        self._parent = parent
        self.value = value

    def parent(self):
        return self._parent

    def __repr__(self):
        return str(self._parent.delpat[self.value])

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
               self.parent()  == other.parent()   and \
               self.value     == other.value

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

