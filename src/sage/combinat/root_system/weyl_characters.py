"""
Weyl Characters
"""
#*****************************************************************************
#       Copyright (C) 2008 Daniel Bump <bump at match.stanford.edu>
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
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram
from sage.combinat.root_system.cartan_type import CartanType
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.structure.element import is_Element
from sage.matrix.constructor import matrix
from sage.rings.all import ZZ, QQ
from sage.misc.misc import repr_lincomb
from sage.misc.functional import is_odd, is_even
from sage.misc.flatten import flatten
from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
import sage.structure.parent_base
from sage.misc.cache import Cache


class WeylCharacter(AlgebraElement):
    """
    A class for Weyl Characters. Let K be a compact Lie group, which we
    assume is semisimple and simply-connected. Its complexified Lie
    algebra L is the Lie algebra of a complex analytic Lie group G. The
    following three categories are equivalent: representations of K;
    representations of L; and analytic representations of G. In every
    case, there is a parametrization of the irreducible representations
    by their highest weight vectors. For this theory of Weyl, see (for
    example) J. F. Adams, Lectures on Lie groups; Broecker and Tom
    Dieck, Representations of Compact Lie groups; Bump, Lie Groups,
    Part I; Fulton and Harris, Representation Theory, Part IV; Goodman
    and Wallach, Representations and Invariants of the Classical
    Groups, Chapter 5; Hall, Lie Groups, Lie Algebras and
    Representations; Humphreys, Introduction to Lie Algebras and their
    representations; Procesi, Lie Groups; Samelson, Notes on Lie
    Algebras; Varadarajan, Lie groups, Lie algebras, and their
    representations; or Zhelobenko, Compact Lie Groups and their
    Representations.

    Computations that you can do with these include computing their
    weight multiplicities, products (thus decomposing the tensor
    product of a representation into irreducibles) and branching
    rules (restriction to a smaller group).

    There is associated with K, L or G as above a lattice, the weight
    lattice, whose elements (called weights) are characters of a Cartan
    subgroup or subalgebra. There is an action of the Weyl group W on
    the lattice, and elements of a fixed fundamental domain for W, the
    positive Weyl chamber, are called dominant. There is for each
    representation a unique highest dominant weight that occurs with
    nonzero multiplicity with respect to a certain partial order, and
    it is called the highest weight vector.

    EXAMPLES::

        sage: L = RootSystem("A2").ambient_space()
        sage: [fw1,fw2] = L.fundamental_weights()
        sage: R = WeylCharacterRing(['A',2], prefix="R")
        sage: [R(1),R(fw1),R(fw2)]
        [R(0,0,0), R(1,0,0), R(1,1,0)]

    Here R(1), R(fw1) and R(fw2) are irreducible representations with
    highest weight vectors 0, and the first two fundamental weights.

    A Weyl character is a character (not necessarily irreducible) of a
    reductive Lie group or algebra. It is represented by a pair of
    dictionaries. The mdict is a dictionary of weight multiplicities.
    The hdict is a dictionary of highest weight vectors of the
    irreducible characters that occur in its decomposition into
    irreducibles, with the multiplicities in this decomposition.

    For type A (also `G_2`, `F_4`, `E_6` and `E_7`) we will take as the weight
    lattice not the weight lattice of the semisimple group, but for a
    larger one. For type A, this means we are concerned with the
    representation theory of K=U(n) or G=GL(n,CC) rather than SU(n) or
    SU(n,CC). This is useful since the representation theory of GL(n)
    is ubiquitous, and also since we may then represent the fundamental
    weights (in root_system.py) by vectors with integer entries. If
    you are only interested in SL(3), say, use
    WeylCharacterRing(['A',2]) as above but be aware that R([a,b,c])
    and R([a+1,b+1,c+1]) represent the same character of SL(3) since
    R([1,1,1]) is the determinant.
    """
    def __init__(self, A, hdict, mdict):
        """
        INPUT:


        -  ``A`` - the WeylCharacterRing

        -  ``hdict`` - dictionary of highest weight vector
           coefficients

        -  ``mdict`` - dictionary of weight multiplicities


        EXAMPLES::

            sage: R = WeylCharacterRing("B3", prefix = "R")
            sage: r =  R(1,1,0)
            sage: r == loads(dumps(r))
            True
        """
        AlgebraElement.__init__(self, A)
        self._hdict = hdict
        self._mdict = mdict
        self._parent = A
        self._space = A._space

    def __repr__(self):
        """
        EXAMPLES::

            sage: R = WeylCharacterRing(['B',3], prefix = "R")
            sage: [R(w) for w in R.fundamental_weights()]
            [R(1,0,0), R(1,1,0), R(1/2,1/2,1/2)]
        """
        if self._hdict == {}:
            return "0"
        v = self._hdict.keys()
        v.sort(key = lambda v: tuple(v.to_vector()))
        return repr_lincomb([self._parent.irr_repr(k) for k in v], [self._hdict[k] for k in v])

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: fw = [B3(w) for w in B3.fundamental_weights()]
            sage: sorted(fw)
            [B3(1/2,1/2,1/2), B3(1,0,0), B3(1,1,0)]
            sage: b = B3(1,0,0)
            sage: b == b
            True
        """
        return cmp(self._hdict, right._hdict)

    def _add_(self, y):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: [A2(0,0,0)+A2(2,1,0), A2(2,1,0)+A2(0,0,0), - A2(0,0,0)+2*A2(0,0,0), -2*A2(0,0,0)+A2(0,0,0), -A2(2,1,0)+2*A2(2,1,0)-A2(2,1,0)]
            [A2(0,0,0) + A2(2,1,0), A2(0,0,0) + A2(2,1,0), A2(0,0,0), -A2(0,0,0), 0]
        """
        hdict = {}
        mdict = {}
        for k in self._hdict.keys():
            if k in y._hdict:
                hdict[k] = self._hdict[k] + y._hdict[k]
                if hdict[k] == 0:
                    del hdict[k]
            else:
                hdict[k] = self._hdict[k]
        for k in y._hdict:
            if not k in self._hdict:
                hdict[k] = y._hdict[k]
        for k in self._mdict.keys():
            if k in y._mdict:
                mdict[k] = self._mdict[k] + y._mdict[k]
                if mdict[k] == 0:
                    del mdict[k]
            else:
                mdict[k] = self._mdict[k]
        for k in y._mdict:
            if not k in self._mdict:
                mdict[k] = y._mdict[k]
        return WeylCharacter(self._parent, hdict, mdict)

    def _neg_(self):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: [-x for x in [A2(0,0,0), 2*A2(0,0,0), -A2(0,0,0), -2*A2(0,0,0)]]
            [-A2(0,0,0), -2*A2(0,0,0), A2(0,0,0), 2*A2(0,0,0)]
        """
        hdict = self._hdict.copy()
        mdict = self._mdict.copy()
        for k in self._hdict:
            hdict[k] = - hdict[k]
            mdict[k] = - mdict[k]
        return WeylCharacter(self._parent, hdict, mdict)

    def _sub_(self, y):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: chi = A2(0,0,0)+2*A2(1,0,0)+3*A2(2,0,0)
            sage: mu =  3*A2(0,0,0)+2*A2(1,0,0)+A2(2,0,0)
            sage: chi - mu
            -2*A2(0,0,0) + 2*A2(2,0,0)
        """
        return self._add_(y._neg_())

    def cartan_type(self):
        """
        Returns the Cartan Type.

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: A2([1,0,0]).cartan_type()
            ['A', 2]
        """
        return self._parent._cartan_type

    def degree(self):
        """
        The degree of the character, that is, the dimension of module.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: [B3(x).degree() for x in B3.fundamental_weights()]
            [7, 21, 8]
        """
        return sum(self._mdict[k] for k in self._mdict)

    def check(self, verbose=False):
        """
        To check the correctness of an element, we compare the theoretical
        dimension computed Weyl character formula with the actual one
        obtained by adding up the multiplicities.

        EXAMPLES::

            sage: B4 = WeylCharacterRing("B4")
            sage: [B4(x).check(verbose = true) for x in B4.fundamental_weights()]
            [[9, 9], [36, 36], [84, 84], [16, 16]]
        """
        theoretical = sum(self._hdict[k]*self._space.weyl_dimension(k) for k in self._hdict)
        practical = sum(self._mdict[k] for k in self._mdict)
        if verbose:
            return [theoretical, practical]
        else:
            return theoretical == practical

    def _mul_(self, y):
        """
        EXAMPLES::

            sage: D4 = WeylCharacterRing(['D',4])
            sage: spin_plus = D4(1/2,1/2,1/2,1/2)
            sage: spin_minus = D4(1/2,1/2,1/2,-1/2)
            sage: spin_plus*spin_minus
            D4(1,0,0,0) + D4(1,1,1,0)
            sage: spin_minus*spin_plus
            D4(1,0,0,0) + D4(1,1,1,0)
        """
        mdict = {}
        for k in self._mdict:
            for l in y._mdict:
                m = k+l
                if m in mdict:
                    mdict[m] += self._mdict[k]*y._mdict[l]
                else:
                    mdict[m] = self._mdict[k]*y._mdict[l]

        for k in mdict.keys():
            if mdict[k] == 0:
                del mdict[k]
        hdict = self._parent.char_from_weights(mdict)
        return WeylCharacter(self._parent, hdict, mdict)

    def __pow__(self, n):
        r"""
        Returns `self^n`.

        The coefficients in `chi^k` are the degrees of those
        irreducible representations of the symmetric group `S_k`
        corresponding to partitions of length `\le 3`.

        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: chi = A2(1,0,0)
            sage: [chi^k for k in range(6)]
            [A2(0,0,0),
             A2(1,0,0),
             A2(1,1,0) + A2(2,0,0),
             A2(1,1,1) + 2*A2(2,1,0) + A2(3,0,0),
             3*A2(2,1,1) + 2*A2(2,2,0) + 3*A2(3,1,0) + A2(4,0,0),
             5*A2(2,2,1) + 6*A2(3,1,1) + 5*A2(3,2,0) + 4*A2(4,1,0) + A2(5,0,0)]
        """
        if not n in ZZ:
            raise TypeError, "exponent must be an integer"
        if not n >= 0:
            raise TypeError, "exponent must be nonnegative"
        z = self._parent.__call__(1)
        for i in range(n):
            z *= self
        return z

    def branch(self, S, rule="default"):
        """
        Returns the restriction of the character to the subalgebra. If no
        rule is specified, we will try to specify one.

        INPUT:


        -  ``S`` - a Weyl character ring for a Lie subgroup or
           subalgebra

        -  ``rule`` - a branching rule.


        See branch_weyl_character? for more information about branching
        rules.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: A2 = WeylCharacterRing(['A',2])
            sage: [B3(w).branch(A2,rule="levi") for w in B3.fundamental_weights()]
            [A2(0,0,-1) + A2(0,0,0) + A2(1,0,0),
             A2(0,-1,-1) + A2(0,0,-1) + A2(0,0,0) + A2(1,0,-1) + A2(1,0,0) + A2(1,1,0),
             A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]
        """
        return branch_weyl_character(self, self._parent, S, rule)

    def hlist(self):
        """
        Returns a list of highest weight vectors and multiplicities of the
        irreducible characters in self.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: B3(1/2,1/2,1/2).hlist()
            [[(1/2, 1/2, 1/2), 1]]
        """
        return [[k,m] for k,m in self._hdict.iteritems()]

    def is_irreducible(self):
        """
        Returns True if self is an irreducible character.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: [B3(x).is_irreducible() for x in B3.fundamental_weights()]
             [True, True, True]
            sage: sum(B3(x) for x in B3.fundamental_weights()).is_irreducible()
             False
        """
        h = self.hlist()
        return len(h) is 1 and h[0][1] is 1

    def symmetric_square(self):
        """
        Returns the symmetric square of the character.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2",style="coroots")
            sage: A2(1,0).symmetric_square()
             A2(2,0)
        """

        cmlist = self.mlist()
        mdict = {}
        for j in range(len(cmlist)):
            for i in range(j+1):
                t = cmlist[i][0]+cmlist[j][0]
                if i < j:
                    coef = cmlist[i][1]*cmlist[j][1]
                else:
                    coef = cmlist[i][1]*(cmlist[i][1]+1)/2
                if t in mdict:
                    mdict[t] += coef
                else:
                    mdict[t] = coef
        hdict = self._parent.char_from_weights(mdict)
        return WeylCharacter(self._parent, hdict, mdict)

    def exterior_square(self):
        """
        Returns the exterior square of the character.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2",style="coroots")
            sage: A2(1,0).exterior_square()
             A2(0,1)
        """
        cmlist = self.mlist()
        mdict = {}
        for j in range(len(cmlist)):
            for i in range(j+1):
                t = cmlist[i][0]+cmlist[j][0]
                if i < j:
                    coef = cmlist[i][1]*cmlist[j][1]
                else:
                    coef = cmlist[i][1]*(cmlist[i][1]-1)/2
                if t in mdict:
                    mdict[t] += coef
                else:
                    mdict[t] = coef
        hdict = self._parent.char_from_weights(mdict)
        return WeylCharacter(self._parent, hdict, mdict)

    def frobenius_schur_indicator(self):
        """
        Returns:

        1 if the representation is real (orthogonal)
        -1 if the representation is quaternionic (symplectic)
        0 if the representation is complex (not self dual)

        The Frobenius-Schur indicator of a character 'chi'
        of a compact group G is the Haar integral over the
        group of 'chi(g^2)'. Its value is 1,-1 or 0. This
        method computes it for irreducible characters of
        compact Lie groups by checking whether the symmetric
        and exterior square characters contain the trivial
        character.

        EXAMPLES::

            sage: B2 = WeylCharacterRing("B2",style="coroots")
            sage: B2(1,0).frobenius_schur_indicator()
             1
            sage: B2(0,1).frobenius_schur_indicator()
             -1
        """
        if not self.is_irreducible():
            raise ValueError, "Frobenius-Schur indicator is only valid for irreducible characters"
        z = self._parent.space().zero()
        if z in self.symmetric_square()._hdict:
            return 1
        elif z in self.exterior_square()._hdict:
            return -1
        else:
            return 0

    def mlist(self):
        """
        Returns a list of weights in self with their multiplicities.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: sorted(B3(1/2,1/2,1/2).mlist())
            [[(-1/2, -1/2, -1/2), 1],  [(-1/2, -1/2, 1/2), 1],  [(-1/2, 1/2, -1/2), 1],  [(-1/2, 1/2, 1/2), 1],  [(1/2, -1/2, -1/2), 1],  [(1/2, -1/2, 1/2), 1],  [(1/2, 1/2, -1/2), 1],  [(1/2, 1/2, 1/2), 1]]
        """
        return [[k,m] for k,m in self._mdict.iteritems()]

    def parent(self):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: B3(2).parent()
            The Weyl Character Ring of Type ['B', 3] with Integer Ring coefficients
        """
        return self._parent


def WeylCharacterRing(ct, base_ring=ZZ, prefix=None, cache=False, style="lattice"):
    r"""
    A class for rings of Weyl characters. The Weyl character is a
    character of a semisimple (or reductive) Lie group or algebra. They
    form a ring, in which the addition and multiplication correspond to
    direct sum and tensor product of representations.

    INPUT:

    - ``ct`` - The Cartan Type

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -  (default: `\ZZ`)

    - ``prefix`` (default an automatically generated prefix
      based on Cartan type)

    - ``cache`` -  (default False) setting cache = True is a substantial
      speedup at the expense of some memory.

    - ``style`` - (default "lattice") can be set style = "coroots"
    to obtain an alternative representation of the elements.

    If no prefix specified, one is generated based on the Cartan type.
    It is good to name the ring after the prefix, since then it can
    parse its own output.

    EXAMPLES::

        sage: G2 = WeylCharacterRing(['G',2])
        sage: [fw1,fw2] = G2.fundamental_weights()
        sage: 2*G2(2*fw1+fw2)
        2*G2(4,-1,-3)
        sage: 2*G2(4,-1,-3)
        2*G2(4,-1,-3)
        sage: G2(4,-1,-3).degree()
        189

    Note that since the ring was named `G_2` after its default prefix, it
    was able to parse its own output. You do not have to use the
    default prefix. Thus:

    EXAMPLES::

        sage: R = WeylCharacterRing(['B',3], prefix='R')
        sage: chi = R(R.fundamental_weights()[3]); chi
        R(1/2,1/2,1/2)
        sage: R(1/2,1/2,1/2) == chi
        True

    You may choose an alternative style of labeling the elements.
    If you create the ring with the option style="coroots" then
    the integers in the label are not the components of the
    highest weight vector, but rather the coefficients when
    the highest weight vector is decomposed into a product
    of irreducibles. These coefficients are the values of the
    coroots on the highest weight vector.

    In the coroot style the Lie group or Lie algebra is treated as
    semisimple, so you lose the distinction between GL(n) and
    SL(n). It gives you output that is comparable to that
    in Tables of Dimensions, Indices and Branching Rules for
    Representations of Simple Lie Algebras (Marcel Dekker, 1981).

    EXAMPLES::

        sage: B3 = WeylCharacterRing("B3",style="coroots")
        sage: [fw1,fw2,fw3]=B3.fundamental_weights()
        sage: fw1+fw3
        (3/2, 1/2, 1/2)
        sage: B3(fw1+fw3)
        B3(1,0,1)
        sage: B3(1,0,1)
        B3(1,0,1)

    For type ['A',r], the coroot representation carries
    less information, since elements of the weight lattice
    that are orthogonal to the coroots are represented as
    zero. This means that with the default style you can
    represent the determinant, but not in the coroot style.
    In the coroot style, elements of the Weyl character
    ring represent characters of SL(r+1,CC), while in the default
    style, they represent characters of GL(r+1,CC).

    EXAMPLES:

        sage: A2 = WeylCharacterRing("A2")
        sage: L = A2.space()
        sage: [A2(L.det()), A2(L(0))]
        [A2(1,1,1), A2(0,0,0)]
        sage: A2(L.det()) == A2(L(0))
        False
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: [A2(L.det()), A2(L(0))]
        [A2(0,0), A2(0,0)]
        sage: A2(L.det()) == A2(L(0))
        True

    The multiplication in a Weyl character ring corresponds to the product of
    characters, which you can use to determine the decomposition of tensor
    products into irreducibles. For example, let us compute the tensor product
    of the standard and spin representations of Spin(7).

    EXAMPLES::

        sage: B3 = WeylCharacterRing("B3")
        sage: [fw1,fw2,fw3]=B3.fundamental_weights()
        sage: [B3(fw1).degree(),B3(fw3).degree()]
        [7, 8]
        sage: B3(fw1)*B3(fw3)
        B3(1/2,1/2,1/2) + B3(3/2,1/2,1/2)

    The name of the irreducible representation encodes the highest
    weight vector.

    TESTS::

        sage: F4 = WeylCharacterRing(['F',4], cache = True)
        sage: [fw1,fw2,fw3,fw4] = F4.fundamental_weights()
        sage: chi = F4(fw4); chi, chi.degree()
        (F4(1,0,0,0), 26)
        sage: chi^2
        F4(0,0,0,0) + F4(1,0,0,0) + F4(1,1,0,0) + F4(3/2,1/2,1/2,1/2) + F4(2,0,0,0)
        sage: [x.degree() for x in [F4(0,0,0,0), F4(1,0,0,0), F4(1,1,0,0), F4(3/2,1/2,1/2,1/2), F4(2,0,0,0)]]
        [1, 26, 52, 273, 324]

    You can produce a list of the irreducible elements of an
    irreducible character.

    EXAMPLES::

        sage: R = WeylCharacterRing(['A',2], prefix = R)
        sage: sorted(R([2,1,0]).mlist())
        [[(1, 1, 1), 2],  [(1, 2, 0), 1],  [(1, 0, 2), 1],  [(2, 1, 0), 1],  [(2, 0, 1), 1],  [(0, 1, 2), 1],  [(0, 2, 1), 1]]
    """
    ct = cartan_type.CartanType(ct)
    return cache_wcr(ct, base_ring=base_ring, prefix=prefix, cache=cache, style=style)

class WeylCharacterRing_class(Algebra):
    def __init__(self, ct, base_ring, prefix, cache, style):
        """
        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3])
            sage: R == loads(dumps(R))
            True
        """
        sage.structure.parent_base.ParentWithBase.__init__(self, base_ring)

        self._cartan_type = ct
        self._rank = ct.rank()
        self._base_ring = base_ring
        self._space = RootSystem(self._cartan_type).ambient_space()
        self._origin = self._space.zero()
        if prefix == None:
            if ct.is_atomic():
                prefix = ct[0]+str(ct[1])
            else:
                prefix = ct.__repr__()
        self._prefix = prefix
        self._style = style
        alpha = self._space.simple_roots()
        Lambda = self._space.fundamental_weights()
        self._cache = cache
        if cache:
            self._irreducibles={}

    def __call__(self, *args):
        """
        Coerces the element into the ring. You may pass a vector in the
        ambient space, an element of the base_ring, or an argument list
        of integers (or half-integers for the spin types) which are the
        components of a vector in the ambient space.

        INPUT:


        -  ``x`` - a ring element to be coerced; or

        -  ``*args`` - the components of a vector


        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: [A2(x) for x in [-2,-1,0,1,2]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0)]
            sage: [A2(2,1,0), A2([2,1,0]), A2(2,1,0)== A2([2,1,0])]
            [A2(2,1,0), A2(2,1,0), True]
            sage: A2([2,1,0]) == A2(2,1,0)
            True
            sage: l = -2*A2(0,0,0) - A2(1,0,0) + A2(2,0,0) + 2*A2(3,0,0)
            sage: [l in A2, A2(l) == l]
            [True, True]
            sage: P.<q> = QQ[]
            sage: A2 = WeylCharacterRing(['A',2], base_ring = P)
            sage: [A2(x) for x in [-2,-1,0,1,2,-2*q,-q,q,2*q,(1-q)]]
            [-2*A2(0,0,0), -A2(0,0,0), 0, A2(0,0,0), 2*A2(0,0,0), -2*q*A2(0,0,0), -q*A2(0,0,0),
            q*A2(0,0,0), 2*q*A2(0,0,0), (-q+1)*A2(0,0,0)]
        """
        if len(args) == 1:
            x = args[0]
        else:
            x = args

        if x == 0 and not x in self._space:
            return WeylCharacter(self, {}, {})

        if x in ZZ:
            hdict = {self._origin: x}
            return WeylCharacter(self, hdict, hdict)

        if self._style == "coroots" and all(xv in ZZ for xv in x):
            x = sum(x[i]*list(self.fundamental_weights())[i] for i in range(self._rank))

        if is_Element(x):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return WeylCharacter(self, x._hdict, x._mdict)
            elif x in self.base_ring():
                hdict = {self._origin: x}
                return WeylCharacter(self, hdict, hdict)

        x = self._space(x)

        # if style == "coroots" and type == 'A' subtract a power of det to put self in SL(r+1,CC)
        if self._style == "coroots":
            x = self.coerce_to_sl(x)

        alphacheck = self._space.simple_coroots()
        vp = [x.inner_product(alphacheck[i]) for i in self._cartan_type.index_set()]

        if not all(v in ZZ for v in vp):
            raise ValueError, "not in weight lattice"
        if not all(v >= 0 for v in vp):
            raise ValueError, "the weight%s is not dominant"%x.__repr__()
        if self._cache and x in self._irreducibles:
            return self._irreducibles[x]
        hdict = {x: 1}
        mdict = irreducible_character_freudenthal(x, self._space)
        ret = WeylCharacter(self, hdict, mdict)
        if self._cache:
            self._irreducibles[x] = ret
        return ret

    def coerce_to_sl(self, x):
        """
        For type ['A',r], this coerces an element of the ambient space into SL(r+1,CC)
        by subtracting a (possibly fractional) power of the determinant.
        """
        if self._cartan_type.is_atomic():
            if self._cartan_type[0] == 'A':
                x = x - self._space.det(sum(x.to_vector())/(self._rank+1))
        else:
            xv = x.to_vector()
            shifts = self._cartan_type._shifts
            types = self._cartan_type.component_types()
            for i in range(len(types)):
                if self._cartan_type.component_types()[i][0] == 'A':
                    s = self._space.ambient_spaces()[i].det(sum(xv[shifts[i]:shifts[i+1]])/(types[i][1]+1))
                    x = x - self._space.inject_weights(i, s)
        return x

    def __repr__(self):
        """
        EXAMPLES::

            sage: WeylCharacterRing("A3")
            The Weyl Character Ring of Type ['A', 3] with Integer Ring coefficients
        """
        return "The Weyl Character Ring of Type %s with %s coefficients"%(self._cartan_type.__repr__(), self._base_ring.__repr__())

    def __cmp__(self, x):
        """
        EXAMPLES::

            sage: A3 = WeylCharacterRing(['A',3])
            sage: B3 = WeylCharacterRing(['B',3])
            sage: A3 == A3
            True
            sage: A3 == B3
            False
        """
        return cmp(repr(self), repr(x))

    def base_ring(self):
        """
        Returns the base ring.

        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3], base_ring = CC); R.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def cartan_type(self):
        """
        Returns the Cartan Type.

        EXAMPLES::

            sage: WeylCharacterRing("A2").cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def fundamental_weights(self):
        """
        Returns the fundamental weights.

        EXAMPLES::

            sage: WeylCharacterRing("G2").fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Returns the simple roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def simple_coroots(self):
        """
        Returns the simple coroots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_coroots()

    def positive_roots(self):
        """
        Returns the positive roots.

        EXAMPLES::

            sage: WeylCharacterRing("G2").positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def dynkin_diagram(self):
        """
        Returns the Dynkin diagram.

        EXAMPLES::

            sage: WeylCharacterRing("E7").dynkin_diagram()
                    O 2
                    |
                    |
            O---O---O---O---O---O
            1   3   4   5   6   7
            E7
        """
        return self.space().dynkin_diagram()

    def extended_dynkin_diagram(self):
        """
        Returns the extended Dynkin diagram, which is the Dynkin diagram
        of the corresponding untwisted affine type.

        EXAMPLES::

            sage: WeylCharacterRing("E7").extended_dynkin_diagram()
                        O 2
                        |
                        |
            O---O---O---O---O---O---O
            0   1   3   4   5   6   7
            E7~
        """
        return DynkinDiagram([self.cartan_type()[0],self.cartan_type()[1],1])

    def rank(self):
        """
        Returns the rank.

        EXAMPLES::

            sage: WeylCharacterRing("G2").rank()
            2
        """
        return self._rank

    def _coerce_impl(self, x):
        """
        Coercion from the base ring.

        EXAMPLES::

            sage: R = WeylCharacterRing(['G',2])
            sage: 2 in R
            True
            sage: R._coerce_impl(2)
            2*G2(0,0,0)
        """
        if x in self._base_ring:
            return self.__call__(x)
        raise TypeError, "no canonical coercion of x"

    def space(self):
        """
        Returns the weight space associated to self.

        EXAMPLES::

            sage: WeylCharacterRing(['E',8]).space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def char_from_weights(self, mdict):
        """
        Return the hdict given just a dictionary of weight multiplicities.

        This will not terminate unless the dictionary of weight
        multiplicities is Weyl group invariant.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: m = B3(1,0,0)._mdict
            sage: B3.char_from_weights(m)
            {(1, 0, 0): 1}
        """
        hdict = {}
        ddict = mdict.copy()
        while not ddict == {}:
            highest = max((x.inner_product(self._space.rho()),x) for x in ddict)[1]
            if not highest.is_dominant():
                raise ValueError, "multiplicity dictionary may not be Weyl group invariant"
            if self._cache and highest in self._irreducibles:
                sdict = self._irreducibles[highest]._mdict
            else:
                sdict = irreducible_character_freudenthal(highest, self._space)
            if self._cache and not highest in self._irreducibles:
                self._irreducibles[highest] = WeylCharacter(self, {highest:1}, sdict)
            c = ddict[highest]
            if highest in hdict:
                hdict[highest] += c
            else:
                hdict[highest] = c
            for k in sdict:
                if k in ddict:
                    if ddict[k] == c*sdict[k]:
                        del ddict[k]
                    else:
                        ddict[k] = ddict[k]-c*sdict[k]
                else:
                    ddict[k] = -c*sdict[k]
        return hdict

    def irr_repr(self, hwv):
        """
        Return a string representing the irreducible character with highest
        weight vectr hwv.

        EXAMPLES::

            sage: B3 = WeylCharacterRing("B3")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(1,1,0)', 'B3(1/2,1/2,1/2)']
            sage: B3 = WeylCharacterRing("B3", style="coroots")
            sage: [B3.irr_repr(v) for v in B3.fundamental_weights()]
            ['B3(1,0,0)', 'B3(0,1,0)', 'B3(0,0,1)']
        """
        if self._style == "lattice":
            vec = hwv.to_vector()
        elif self._style == "coroots":
            vec = [hwv.inner_product(x) for x in self.simple_coroots()]
        else:
            raise ValueError, "unknown style"
        hstring = str(vec[0])
        for i in range(1,len(vec)):
            hstring=hstring+","+str(vec[i])
        return self._prefix+"("+hstring+")"

cache_wcr = Cache(WeylCharacterRing_class)


def irreducible_character_freudenthal(hwv, L, debug=False):
    """
    Returns the dictionary of multiplicities for the irreducible
    character with highest weight lamb. The weight multiplicities are
    computed by the Freudenthal multiplicity formula. The algorithm is
    based on recursion relation that is stated, for example, in
    Humphrey's book on Lie Algebras. The multiplicities are invariant
    under the Weyl group, so to compute them it would be sufficient to
    compute them for the weights in the positive Weyl chamber. However
    after some testing it was found to be faster to compute every
    weight using the recursion, since the use of the Weyl group is
    expensive in its current implementation.

    INPUT:

    - ``hwv`` - a dominant weight in a weight lattice.

    - ``L`` - the ambient space
    """

    rho = L.rho()
    mdict = {}
    current_layer = {hwv:1}

    simple_roots = L.simple_roots()
    positive_roots = L.positive_roots()

    while len(current_layer) > 0:
        next_layer = {}
        for mu in current_layer:
            if current_layer[mu] != 0:
                mdict[mu] = current_layer[mu]
                for alpha in simple_roots:
                    next_layer[mu-alpha] = None
        if debug:
            print next_layer

        for mu in next_layer:
            if next_layer[mu] is None:
                accum = 0
                for alpha in positive_roots:
                    mu_plus_i_alpha = mu + alpha
                    while mu_plus_i_alpha in mdict:
                        accum += mdict[mu_plus_i_alpha]*(mu_plus_i_alpha).inner_product(alpha)
                        mu_plus_i_alpha += alpha

                if accum == 0:
                    next_layer[mu] = 0
                else:
                    hwv_plus_rho = hwv + rho
                    mu_plus_rho  = mu  + rho
                    next_layer[mu] = QQ(2*accum)/QQ((hwv_plus_rho).inner_product(hwv_plus_rho)-
                                                    (mu_plus_rho).inner_product(mu_plus_rho))

        current_layer = next_layer
    return mdict

def branch_weyl_character(chi, R, S, rule="default"):
    r"""
    A Branching rule describes the restriction of representations from
    a Lie group or algebra G to a smaller one. See for example, R. C.
    King, Branching rules for classical Lie groups using tensor and
    spinor methods. J. Phys. A 8 (1975), 429-449, Howe, Tan and
    Willenbring, Stable branching rules for classical symmetric pairs,
    Trans. Amer. Math. Soc. 357 (2005), no. 4, 1601-1626, McKay and
    Patera, Tables of Dimensions, Indices and Branching Rules for
    Representations of Simple Lie Algebras (Marcel Dekker, 1981),
    and Fauser, Jarvis, King and Wybourne, New branching rules induced
    by plethysm. J. Phys. A 39 (2006), no. 11, 2611--2655.

    INPUT:

    - ``chi`` - a character of G

    - ``R`` - the Weyl Character Ring of G

    - ``S`` - the Weyl Character Ring of H

    - ``rule`` - a set of r dominant weights in H where r is the rank
      of G.

    You may use a predefined rule by specifying rule = one of"levi",
    "automorphic", "symmetric", "extended", "triality" or
    "miscellaneous". The use of these rules will be explained next.
    After the examples we will explain how to write your own branching
    rules for cases that we have omitted.

    To explain the predefined rules we survey the most important
    branching rules. These may be classified into several cases, and
    once this is understood, the detailed classification can be read
    off from the Dynkin diagrams. Dynkin classified the maximal
    subgroups of Lie groups in Mat. Sbornik N.S. 30(72):349-462
    (1952).

    We will list give predefined rules that cover most cases where the
    branching rule is to a maximal subgroup. For convenience, we
    also give some branching rules to subgroups that are not maximal.
    For example, a Levi subgroup may or may not be maximal.

    LEVI TYPE. These can be read off from the Dynkin diagram. If
    removing a node from the Dynkin diagram produces another Dynkin
    diagram, there is a branching rule. Currently we require that the
    smaller diagram be connected. For these rules use the option
    rule="levi"::

       ['A',r] => ['A',r-1]
       ['B',r] => ['A',r-1]
       ['B',r] => ['B',r-1]
       ['C',r] => ['A',r-1]
       ['C',r] => ['C',r-1]
       ['D',r] => ['A',r-1]
       ['D',r] => ['D',r-1]
       ['E',r] => ['A',r-1] r = 6,7,8 (not implemented yet)
       ['E',r] => ['D',r-1] r = 6,7,8 (not implemented yet)
       ['E',r] => ['E',r-1] r = 6,7 (not implemented yet)
       F4 => B3
       F4 => C3
       G2 => A1 (short root)

    The other Levi branching rule from G2 => A1 corresponding to the
    long root is available by first branching G_2 => A_2 then A2 => A1.

    AUTOMORPHIC TYPE. If the Dynkin diagram has a symmetry, then there
    is an automorphism that is a special case of a branching rule.
    There is also an exotic "triality" automorphism of D4 having order
    3. Use rule="automorphic" or (for D4) rule="triality"

    ['A',r] => ['A',r]
    ['D',r] => ['D',r]
    E6 => E6 (not implemented yet)

    SYMMETRIC TYPE. Related to the automorphic type, when either
    the Dynkin diagram or the extended diagram has a symmetry
    there is a branching rule to the subalgebra (or subgroup) of
    invariants under the automorphism. Use rule="symmetric".
    The last branching rule, D4=>G2 is not to a maximal subgroup
    since D4=>B3=>G2, but it is included for convenience.

    ['A',2r+1] => ['B',r]
    ['A',2r] => ['C',r]
    ['A',2r] => ['D',r]
    ['D',r] => ['B',r-1]
    E6 => F4
    D4 => G2

    EXTENDED TYPE. If removing a node from the extended Dynkin diagram
    results in a Dynkin diagram, then there is a branching rule. Use
    rule="extended" for these. We will also use this classification
    for some rules that are not of this type, mainly involving type B,
    such as D6 => B3xB3.

    G2 => A2
    ['B',r] => ['D',r]
    F4 => B4
    E7 => A7 (not implemented yet)
    E8 => A8 (not implemented yet)

    Here is the extended Dynkin diagram for D6:

        0       6
        O       O
        |       |
        |       |
    O---O---O---O---O
    1   2   3   4   6

    Removing the node 3 results in an embedding D3xD3 -> D6. This
    corresponds to the embedding SO(6)xSO(6) -> SO(12), and is of
    extended type. On the other hand the embedding SO(5)xSO(7)-->SO(12)
    (e.g. B2xB3 -> D6) cannot be explained this way but for
    uniformity is implemented under rule="extended".

    Using rule="extended" you can get any branching rule
    SO(n) => SO(a) x SO(b) x SO(c) x ... where n = a+b+c+ ...
    Sp(2n) => Sp(2a) x Sp(2b) x Sp(2c) x ... where n = a+b+c+ ...
    where O(a) = ['D',r] (a=2r) or ['B',r] (a=2r+1)
    and Sp(2r)=['C',r].

    TENSOR: There are branching rules

    ['A', rs-1] => ['A',r-1] x ['A',s-1]
    ['B',2rs+r+s] => ['B',r] x ['B',s]
    ['D',2rs+s] => ['B',r] x ['D',s]
    ['D',2rs] => ['D',r] x ['D',s]
    ['D',2rs] => ['C',r] x ['C',s]
    ['C',2rs+s] => ['B',r] x ['C',s]
    ['C',2rs] => ['C',r] x ['D',s].

    corresponding to the tensor product homomorphism. For type
    A, the homomorphism is GL(r) x GL(s) -> GL(rs). For the
    classical types, the relevant fact is that if V,W are
    orthogonal or symplectic spaces, that is, spaces endowed
    with symmetric or skew-symmetric bilinear forms, then V
    tensor W is also an orthogonal space (if V and W are both
    orthogonal or both symplectic) or symplectic (if one of
    V and W is orthogonal and the other symplectic).

    The corresponding branching rules are obtained using rule="tensor".


    SYMMETRIC POWER: The k-th symmetric and exterior power homomorphisms
    map GL(n) --> GL(binomial(n+k-1,k)) and GL(binomial(n,k)). The
    corresponding branching rules are not implemented but a special
    case is. The k-th symmetric power homomorphism SL(2) --> GL(k+1)
    has its image inside of SO(2r+1) if k=2r and inside of Sp(2r) if
    k=2r-1. Hence there are branching rules

    ['B',r] => A1
    ['C',r] => A1

    and these may be obtained using the rule "symmetric_power".

    MISCELLANEOUS: Use rule="miscellaneous" for the following rules.

    B3 => G2
    F4 => G2xA1 (not implemented yet)

    BRANCHING RULES FROM PLETHYSMS

    Nearly all branching rules G => H where G is of type A,B,C or D
    are covered by the preceding rules. The function
    branching_rules_from_plethysm covers the remaining cases.

    ISOMORPHIC TYPE: Although not usually referred to as a branching
    rule, the effects of the accidental isomorphisms may be handled
    using rule="isomorphic"

    B2 => C2
    C2 => B2
    A3 => D3
    D3 => A3
    D2 => A1xA1
    B1 => A1
    C1 => A1

    EXAMPLES: (Levi type)

    ::

        sage: A1 = WeylCharacterRing("A1")
        sage: A2 = WeylCharacterRing("A2")
        sage: A3 = WeylCharacterRing("A3")
        sage: A4 = WeylCharacterRing("A4")
        sage: A5 = WeylCharacterRing("A5")
        sage: B2 = WeylCharacterRing("B2")
        sage: B3 = WeylCharacterRing("B3")
        sage: B4 = WeylCharacterRing("B4")
        sage: C2 = WeylCharacterRing("C2")
        sage: C3 = WeylCharacterRing("C3")
        sage: D3 = WeylCharacterRing("D3")
        sage: D4 = WeylCharacterRing("D4")
        sage: D5 = WeylCharacterRing("D5")
        sage: G2 = WeylCharacterRing("G2")
        sage: F4 = WeylCharacterRing("F4") # long time
        sage: [B3(w).branch(A2,rule="levi") for w in B3.fundamental_weights()]
        [A2(0,0,-1) + A2(0,0,0) + A2(1,0,0),
         A2(0,-1,-1) + A2(0,0,-1) + A2(0,0,0) + A2(1,0,-1) + A2(1,0,0) + A2(1,1,0),
         A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]

    The last example must be understood as follows. The representation
    of B3 being branched is spin, which is not a representation of
    SO(7) but of its double cover spin(7). The group A2 is really GL(3)
    and the double cover of SO(7) induces a cover of GL(3) that is
    trivial over SL(3) but not over the center of GL(3). The weight
    lattice for this GL(3) consists of triples (a,b,c) of half integers
    such that a-b and b-c are in `\ZZ`, and this is reflected in the last
    decomposition.

    ::

        sage: [C3(w).branch(A2,rule="levi") for w in C3.fundamental_weights()]
        [A2(0,0,-1) + A2(1,0,0),
         A2(0,-1,-1) + A2(1,0,-1) + A2(1,1,0),
         A2(-1,-1,-1) + A2(1,-1,-1) + A2(1,1,-1) + A2(1,1,1)]
        sage: [D4(w).branch(A3,rule="levi") for w in D4.fundamental_weights()]
        [A3(0,0,0,-1) + A3(1,0,0,0),
         A3(0,0,-1,-1) + A3(0,0,0,0) + A3(1,0,0,-1) + A3(1,1,0,0),
         A3(1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,-1/2),
         A3(-1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,1/2)]
        sage: [B3(w).branch(B2,rule="levi") for w in B3.fundamental_weights()]
        [2*B2(0,0) + B2(1,0), B2(0,0) + 2*B2(1,0) + B2(1,1), 2*B2(1/2,1/2)]
        sage: C3 = WeylCharacterRing(['C',3])
        sage: [C3(w).branch(C2,rule="levi") for w in C3.fundamental_weights()]
        [2*C2(0,0) + C2(1,0),
         C2(0,0) + 2*C2(1,0) + C2(1,1),
         C2(1,0) + 2*C2(1,1)]
        sage: [D5(w).branch(D4,rule="levi") for w in D5.fundamental_weights()]
        [2*D4(0,0,0,0) + D4(1,0,0,0),
         D4(0,0,0,0) + 2*D4(1,0,0,0) + D4(1,1,0,0),
         D4(1,0,0,0) + 2*D4(1,1,0,0) + D4(1,1,1,0),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2)]
        sage: G2(1,0,-1).branch(A1,rule="levi")
         A1(0,-1) + A1(1,-1) + A1(1,0)
        sage: [F4(fw).branch(B3,rule="levi") for fw in F4.fundamental_weights()] # long time
         [B3(0,0,0) + 2*B3(1/2,1/2,1/2) + 2*B3(1,0,0) + B3(1,1,0),
         B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 5*B3(1,0,0) + 7*B3(1,1,0) + 3*B3(1,1,1)
         + 6*B3(3/2,1/2,1/2) + 2*B3(3/2,3/2,1/2) + B3(2,0,0) + 2*B3(2,1,0) + B3(2,1,1),
         3*B3(0,0,0) + 6*B3(1/2,1/2,1/2) + 4*B3(1,0,0) + 3*B3(1,1,0) + B3(1,1,1) + 2*B3(3/2,1/2,1/2),
         3*B3(0,0,0) + 2*B3(1/2,1/2,1/2) + B3(1,0,0)]
        sage: [F4(fw).branch(C3,rule="levi") for fw in F4.fundamental_weights()] # long time
         [3*C3(0,0,0) + 2*C3(1,1,1) + C3(2,0,0),
         3*C3(0,0,0) + 6*C3(1,1,1) + 4*C3(2,0,0) + 2*C3(2,1,0) + 3*C3(2,2,0) + C3(2,2,2) + C3(3,1,0) + 2*C3(3,1,1),
         2*C3(1,0,0) + 3*C3(1,1,0) + C3(2,0,0) + 2*C3(2,1,0) + C3(2,1,1),
         2*C3(1,0,0) + C3(1,1,0)]
        sage: A1xA1 = WeylCharacterRing("A1xA1")
        sage: [A3(hwv).branch(A1xA1,rule="levi") for hwv in A3.fundamental_weights()]
        [A1xA1(0,0,1,0) + A1xA1(1,0,0,0),
         A1xA1(0,0,1,1) + A1xA1(1,0,1,0) + A1xA1(1,1,0,0),
         A1xA1(1,0,1,1) + A1xA1(1,1,1,0)]
        sage: A1xB1=WeylCharacterRing("A1xB1",style="coroots")
        sage: [B3(x).branch(A1xB1,rule="levi") for x in B3.fundamental_weights()]
        [A1xB1(0,2) + 2*A1xB1(1,0),
         3*A1xB1(0,0) + A1xB1(0,2) + 2*A1xB1(1,2) + A1xB1(2,0),
         2*A1xB1(0,1) + A1xB1(1,1)]

    EXAMPLES: (Automorphic type, including D4 triality)

    ::

        sage: [A3(chi).branch(A3,rule="automorphic") for chi in A3.fundamental_weights()]
        [A3(0,0,0,-1), A3(0,0,-1,-1), A3(0,-1,-1,-1)]
        sage: [D4(chi).branch(D4,rule="automorphic") for chi in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1/2,1/2,1/2,-1/2)]
        sage: [D4(chi).branch(D4,rule="triality") for chi in D4.fundamental_weights()]
        [D4(1/2,1/2,1/2,-1/2), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1,0,0,0)]

    EXAMPLES: (Symmetric type)

    ::

        sage: [w.branch(B2,rule="symmetric") for w in [A4(1,0,0,0,0),A4(1,1,0,0,0),A4(1,1,1,0,0),A4(2,0,0,0,0)]]
        [B2(1,0), B2(1,1), B2(1,1), B2(0,0) + B2(2,0)]
        sage: [A5(w).branch(C3,rule="symmetric") for w in A5.fundamental_weights()]
        [C3(1,0,0), C3(0,0,0) + C3(1,1,0), C3(1,0,0) + C3(1,1,1), C3(0,0,0) + C3(1,1,0), C3(1,0,0)]
        sage: [A5(w).branch(D3,rule="symmetric") for w in A5.fundamental_weights()]
        [D3(1,0,0), D3(1,1,0), D3(1,1,-1) + D3(1,1,1), D3(1,1,0), D3(1,0,0)]
        sage: [D4(x).branch(B3,rule="symmetric") for x in D4.fundamental_weights()]
        [B3(0,0,0) + B3(1,0,0), B3(1,0,0) + B3(1,1,0), B3(1/2,1/2,1/2), B3(1/2,1/2,1/2)]
        sage: [D4(x).branch(G2,rule="symmetric") for x in D4.fundamental_weights()]
        [G2(0,0,0) + G2(1,0,-1), 2*G2(1,0,-1) + G2(2,-1,-1), G2(0,0,0) + G2(1,0,-1), G2(0,0,0) + G2(1,0,-1)]
        sage: E6=WeylCharacterRing("E6",style="coroots") # long time
        sage: F4=WeylCharacterRing("F4",style="coroots") # long time
        sage: [E6(fw).branch(F4,rule="symmetric") for fw in E6.fundamental_weights()] # long time
        [F4(0,0,0,0) + F4(0,0,0,1),
         F4(0,0,0,1) + F4(1,0,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(1,0,0,0) + 2*F4(0,0,1,0) + F4(1,0,0,1) + F4(0,1,0,0),
         F4(0,0,0,1) + F4(1,0,0,0) + F4(0,0,1,0),
         F4(0,0,0,0) + F4(0,0,0,1)]

    EXAMPLES: (Extended type)

    ::

        sage: [B3(x).branch(D3,rule="extended") for x in B3.fundamental_weights()]
        [D3(0,0,0) + D3(1,0,0),
         D3(1,0,0) + D3(1,1,0),
         D3(1/2,1/2,-1/2) + D3(1/2,1/2,1/2)]
        sage: [G2(w).branch(A2, rule="extended") for w in G2.fundamental_weights()]
        [A2(0,0,0) + A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3),
         A2(1/3,1/3,-2/3) + A2(2/3,-1/3,-1/3) + A2(1,0,-1)]
        sage: [F4(fw).branch(B4,rule="extended") for fw in F4.fundamental_weights()] # long time
        [B4(1/2,1/2,1/2,1/2) + B4(1,1,0,0),
         B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2) + B4(3/2,3/2,1/2,1/2) + B4(2,1,1,0),
         B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0) + B4(1,1,0,0) + B4(1,1,1,0) + B4(3/2,1/2,1/2,1/2),
         B4(0,0,0,0) + B4(1/2,1/2,1/2,1/2) + B4(1,0,0,0)]

        sage: E6 = WeylCharacterRing("E6", style="coroots")
        sage: A2xA2xA2=WeylCharacterRing("A2xA2xA2",style="coroots")
        sage: A5xA1=WeylCharacterRing("A5xA1",style="coroots")
        sage: G2 = WeylCharacterRing("G2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: F4 = WeylCharacterRing("F4",style="coroots")
        sage: A3xA1 = WeylCharacterRing("A3xA1", style="coroots")
        sage: A2xA2 = WeylCharacterRing("A2xA2", style="coroots")
        sage: A1xC3 = WeylCharacterRing("A1xC3",style="coroots")

        sage: E6(1,0,0,0,0,0).branch(A5xA1,rule="extended") # long time
         A5xA1(0,0,0,1,0,0) + A5xA1(1,0,0,0,0,1)
        sage: E6(1,0,0,0,0,0).branch(A2xA2xA2, rule="extended") # long time
         A2xA2xA2(0,0,0,1,1,0) + A2xA2xA2(0,1,1,0,0,0) + A2xA2xA2(1,0,0,0,0,1)
        sage: F4(1,0,0,0).branch(A1xC3,rule="extended") # long time
         A1xC3(0,2,0,0) + A1xC3(1,0,0,1) + A1xC3(2,0,0,0)
        sage: G2(0,1).branch(A1xA1, rule="extended")
         A1xA1(0,2) + A1xA1(2,0) + A1xA1(3,1)
        sage: F4(0,0,0,1).branch(A2xA2, rule="extended") # long time
         A2xA2(0,0,1,1) + A2xA2(0,1,0,1) + A2xA2(1,0,1,0)
        sage: F4(0,0,0,1).branch(A3xA1,rule="extended") # long time
         A3xA1(0,0,0,0) + A3xA1(0,0,0,2) + A3xA1(0,0,1,1) + A3xA1(0,1,0,0) + A3xA1(1,0,0,1)
        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: D2xD2=WeylCharacterRing("D2xD2",style="coroots") # We get D4 => A1xA1xA1xA1 by remembering that A1xA1 = D2.
        sage: [D4(fw).branch(D2xD2, rule="extended") for fw in D4.fundamental_weights()]
        [D2xD2(0,0,1,1) + D2xD2(1,1,0,0),
         D2xD2(0,0,2,0) + D2xD2(0,0,0,2) + D2xD2(2,0,0,0) + D2xD2(1,1,1,1) + D2xD2(0,2,0,0),
         D2xD2(1,0,0,1) + D2xD2(0,1,1,0),
         D2xD2(1,0,1,0) + D2xD2(0,1,0,1)]


    EXAMPLES: (Tensor type)

    ::
        sage: A5=WeylCharacterRing("A5", style="coroots")
        sage: A2xA1=WeylCharacterRing("A2xA1", style="coroots")
        sage: [A5(hwv).branch(A2xA1, rule="tensor") for hwv in A5.fundamental_weights()]
        [A2xA1(1,0,1),
         A2xA1(0,1,2) + A2xA1(2,0,0),
         A2xA1(0,0,3) + A2xA1(1,1,1),
         A2xA1(1,0,2) + A2xA1(0,2,0),
         A2xA1(0,1,1)]

        sage: B4=WeylCharacterRing("B4",style="coroots")
        sage: B1xB1=WeylCharacterRing("B1xB1",style="coroots")
        sage: [B4(f).branch(B1xB1,rule="tensor") for f in B4.fundamental_weights()]
        [B1xB1(2,2),
         B1xB1(0,2) + B1xB1(2,0) + B1xB1(2,4) + B1xB1(4,2),
         B1xB1(0,2) + B1xB1(0,6) + B1xB1(2,0) + B1xB1(2,2) + B1xB1(2,4) + B1xB1(4,2) + B1xB1(4,4) + B1xB1(6,0),
         B1xB1(1,3) + B1xB1(3,1)]

        sage: D4=WeylCharacterRing("D4",style="coroots")
        sage: C2xC1=WeylCharacterRing("C2xC1",style="coroots")
        sage: [D4(f).branch(C2xC1,rule="tensor") for f in D4.fundamental_weights()]
        [C2xC1(1,0,1),
         C2xC1(0,0,2) + C2xC1(0,1,2) + C2xC1(2,0,0),
         C2xC1(1,0,1),
         C2xC1(0,0,2) + C2xC1(0,1,0)]

        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: B1xC1=WeylCharacterRing("B1xC1",style="coroots")
        sage: [C3(f).branch(B1xC1,rule="tensor") for f in C3.fundamental_weights()]
        [B1xC1(2,1), B1xC1(2,2) + B1xC1(4,0), B1xC1(0,3) + B1xC1(4,1)]

    EXAMPLES: (Symmetric Power)

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: B3=WeylCharacterRing("B3",style="coroots")
        sage: C3=WeylCharacterRing("C3",style="coroots")
        sage: [B3(fw).branch(A1,rule="symmetric_power") for fw in B3.fundamental_weights()]
        [A1(6), A1(2) + A1(6) + A1(10), A1(0) + A1(6)]
        sage: [C3(fw).branch(A1,rule="symmetric_power") for fw in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    EXAMPLES: (Miscellaneous type)

    ::

        sage: G2 = WeylCharacterRing("G2")
        sage: [fw1, fw2, fw3] = B3.fundamental_weights()
        sage: B3(fw1+fw3).branch(G2, rule="miscellaneous")
        G2(1,0,-1) + G2(2,-1,-1) + G2(2,0,-2)

    EXAMPLES: (Isomorphic type)

    ::

        sage: [B2(x).branch(C2, rule="isomorphic") for x in B2.fundamental_weights()]
        [C2(1,1), C2(1,0)]
        sage: [C2(x).branch(B2, rule="isomorphic") for x in C2.fundamental_weights()]
        [B2(1/2,1/2), B2(1,0)]
        sage: [A3(x).branch(D3,rule="isomorphic") for x in A3.fundamental_weights()]
        [D3(1/2,1/2,1/2), D3(1,0,0), D3(1/2,1/2,-1/2)]
        sage: [D3(x).branch(A3,rule="isomorphic") for x in D3.fundamental_weights()]
        [A3(1/2,1/2,-1/2,-1/2), A3(1/4,1/4,1/4,-3/4), A3(3/4,-1/4,-1/4,-1/4)]

    Here A3(x,y,z,w) can be understood as a representation of SL(4).
    The weights x,y,z,w and x+t,y+t,z+t,w+t represent the same
    representation of SL(4) - though not of GL(4) - since
    A3(x+t,y+t,z+t,w+t) is the same as A3(x,y,z,w) tensored with
    `det^t`. So as a representation of SL(4),
    A3(1/4,1/4,1/4,-3/4) is the same as A3(1,1,1,0). The exterior
    square representation SL(4) -> GL(6) admits an invariant symmetric
    bilinear form, so is a representation SL(4) -> SO(6) that lifts to
    an isomorphism SL(4) -> Spin(6). Conversely, there are two
    isomorphisms SO(6) -> SL(4), of which we've selected one.

    In cases like this you might prefer style="coroots".

        sage: A3 = WeylCharacterRing("A3",style="coroots")
        sage: D3 = WeylCharacterRing("D3",style="coroots")
        sage: [D3(fw) for fw in D3.fundamental_weights()]
        [D3(1,0,0), D3(0,1,0), D3(0,0,1)]
        sage: [D3(fw).branch(A3,rule="isomorphic") for fw in D3.fundamental_weights()]
        [A3(0,1,0), A3(0,0,1), A3(1,0,0)]
        sage: D2 = WeylCharacterRing("D2", style="coroots")
        sage: A1xA1 = WeylCharacterRing("A1xA1", style="coroots")
        sage: [D2(fw).branch(A1xA1,rule="isomorphic") for fw in D2.fundamental_weights()]
        [A1xA1(1,0), A1xA1(0,1)]

    EXAMPLES: (Branching rules from plethysms)

    This is a general rule that includes any branching rule
    from types A,B,C or D as a special case. Thus it could be
    used in place of the above rules and would give the same
    results.

    We consider a homomorphism H --> G where G is one of
    SL(r+1), SO(2r+1), Sp(2r) or SO(2r). The function
    branching_rule_from_plethysm produces the corresponding
    branching rule. The main ingredient is the character
    chi of the representation of H that is the homomorphism
    to GL(r+1), GL(2r+1) or GL(2r).

    This rule is so powerful that it contains the other
    rules implemented above as special cases. First let
    us consider the symmetric fifth power representation
    of SL(2).

        sage: A1=WeylCharacterRing("A1",style="coroots")
        sage: chi=A1([5])
        sage: chi.degree()
         6
        sage: chi.frobenius_schur_indicator()
        -1

    This confirms that the character has degree 6 and
    is symplectic, so it corresponds to a homomorphism
    SL(2) --> Sp(6), and there is a corresponding
    branching rule C3 => A1.

        sage: C3 = WeylCharacterRing("C3",style="coroots")
        sage: sym5rule = branching_rule_from_plethysm(chi,"C3")
        sage: [C3(hwv).branch(A1,rule=sym5rule) for hwv in C3.fundamental_weights()]
        [A1(5), A1(4) + A1(8), A1(3) + A1(9)]

    This is identical to the results we would obtain using
    rule="symmetric_power". The next example gives a branching
    not available by other standard rules.

        sage: G2 = WeylCharacterRing("G2",style="coroots")
        sage: D7 = WeylCharacterRing("D7",style="coroots")
        sage: ad=G2(0,1); ad.degree(); ad.frobenius_schur_indicator()
         14
         1
        sage: spin = D7(0,0,0,0,0,1,0); spin.degree()
         64
        sage: spin.branch(G2, rule=branching_rule_from_plethysm(ad, "D7"))
         G2(1,1)

    We have confirmed that the adjoint representation of G2
    gives a homomorphism into SO(14), and that the pullback
    of the one of the two 64 dimensional spin representations
    to SO(14) is an irreducible representation of G2.

    BRANCHING FROM A REDUCIBLE ROOT SYSTEM

    If you are branching from a reducible root system, the rule is
    a list of rules, one for each component type in the root system.
    The rules in the list are given in pairs [type, rule], where
    type is the root system to be branched to, and rule is the
    branching rule.

        sage: D4 = WeylCharacterRing("D4",style="coroots")
        sage: D2xD2 = WeylCharacterRing("D2xD2",style="coroots")
        sage: A1xA1xA1xA1 = WeylCharacterRing("A1xA1xA1xA1",style="coroots")
        sage: rr = [["A1xA1","isomorphic"],["A1xA1","isomorphic"]]
        sage: [D4(fw) for fw in D4.fundamental_weights()]
        [D4(1,0,0,0), D4(0,1,0,0), D4(0,0,1,0), D4(0,0,0,1)]
        sage: [D4(fw).branch(D2xD2,rule="extended").branch(A1xA1xA1xA1,rule=rr) for fw in D4.fundamental_weights()]
        [A1xA1xA1xA1(0,0,1,1) + A1xA1xA1xA1(1,1,0,0),
         A1xA1xA1xA1(0,0,0,2) + A1xA1xA1xA1(0,0,2,0) + A1xA1xA1xA1(0,2,0,0) + A1xA1xA1xA1(1,1,1,1) + A1xA1xA1xA1(2,0,0,0),
         A1xA1xA1xA1(0,1,1,0) + A1xA1xA1xA1(1,0,0,1),
         A1xA1xA1xA1(0,1,0,1) + A1xA1xA1xA1(1,0,1,0)]

    WRITING YOUR OWN RULES

    Suppose you want to branch from a group G to a subgroup H.
    Arrange the embedding so that a Cartan subalgebra U of H is
    contained in a Cartan subalgebra T of G. There is thus
    a mapping from the weight spaces Lie(T)* --> Lie(U)*.
    The embedding must be chosen in such a way that the
    restriction of the image of the positive Weyl chamber
    in Lie(T)* is contained in the positive Weyl chamber
    in Lie(U)*.

    The RULE is this map Lie(T)* = G.space() to Lie(U)* = H.space(),
    which you may implement as a function. As an example, let
    us consider how to implement the branching rule A3 => C2.
    Here H = C2 = Sp(4) embedded as a subgroup in A3 = GL(4). The
    Cartan subalgebra U consists of diagonal matrices with
    eigenvalues u1, u2, -u2, -u1. The C2.space() is the
    two dimensional vector spaces consisting of the linear
    functionals u1 and u2 on U. On the other hand Lie(T) is
    RR^4. A convenient way to see the restriction is to
    think of it as the adjoint of the map [u1,u2] -> [u1,u2,-u2,-u1],
    that is, [x0,x1,x2,x3] -> [x0-x3,x1-x2]. Hence we may
    encode the rule:

    def rule(x):
        return [x[0]-x[3],x[1]-x[2]]

    or simply:

    rule = lambda x : [x[0]-x[3],x[1]-x[2]]

    EXAMPLES::

        sage: A3 = WeylCharacterRing(['A',3])
        sage: C2 = WeylCharacterRing(['C',2])
        sage: rule = lambda x : [x[0]-x[3],x[1]-x[2]]
        sage: branch_weyl_character(A3([1,1,0,0]),A3,C2,rule)
        C2(0,0) + C2(1,1)
        sage: A3(1,1,0,0).branch(C2, rule) == C2(0,0) + C2(1,1)
        True
    """
    if type(rule) == str:
        rule = get_branching_rule(R._cartan_type, S._cartan_type, rule)
    elif R._cartan_type.is_compound():
        Rtypes = R._cartan_type.component_types()
        Stypes = [CartanType(l[0]) for l in rule]
        rules = [l[1] for l in rule]
        ntypes = len(Rtypes)
        rule_list = [get_branching_rule(Rtypes[i], Stypes[i], rules[i]) for i in range(ntypes)]
        shifts = R._cartan_type._shifts
        def rule(x):
            yl = []
            for i in range(ntypes):
                yl.append(rule_list[i](x[shifts[i]:shifts[i+1]]))
            return flatten(yl)
    mdict = {}
    for k in chi._mdict:
        if S._style == "coroots":
            h = S.coerce_to_sl(S._space(rule(list(k.to_vector()))))
        else:
            h = S._space(rule(list(k.to_vector())))
        if h in mdict:
            mdict[h] += chi._mdict[k]
        else:
            mdict[h] = chi._mdict[k]
    hdict = S.char_from_weights(mdict)
    return WeylCharacter(S, hdict, mdict)

def get_branching_rule(Rtype, Stype, rule):
    """
    INPUT:

    - ``R`` - the Weyl Character Ring of G

    - ``S`` - the Weyl Character Ring of H

    - ``rule`` - a string describing the branching rule as a map from
      the weight space of S to the weight space of R.
    """
    r = Rtype.rank()
    s = Stype.rank()
    rdim = Rtype.root_system().ambient_space().dimension()
    sdim = Stype.root_system().ambient_space().dimension()
    if Stype.is_compound():
        stypes = Stype.component_types()
    if rule == "default":
        return lambda x : x
    elif rule == "levi":
        if not s == r-1:
            raise ValueError, "Incompatible ranks"
        if Rtype[0] == 'A':
            if Stype.is_compound():
                if all(ct[0]=='A' for ct in stypes) \
                       and rdim == sdim:
                    return lambda x : x
                else:
                    raise ValueError, "Rule not found"
            elif Stype[0] == 'A':
                return lambda x : list(x)[:r]
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] in ['B', 'C', 'D']:
            if Stype.is_atomic():
                if Stype[0] == 'A':
                    return  lambda x : x
                elif Stype[0] == Rtype[0]:
                    return lambda x : list(x)[1:]
            elif stypes[-1][0] == Rtype[0] and all(t[0] == 'A' for t in stypes[:-1]):
                return lambda x : x
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'E' and Stype[0] in ['A','D','E']:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif Rtype[0] == 'F' and s == 3:
            if Stype[0] == 'B':
                return lambda x : list(x)[1:]
            elif Stype[0] == 'C':
                return lambda x : [x[1]-x[0],x[2]+x[3],x[2]-x[3]]
            else:
                raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif Rtype[0] == 'G' and Stype[0] == 'A':
            return lambda x : list(x)[1:]
        else:
            raise ValueError, "Rule not found"
    elif rule == "automorphic":
        if not Rtype == Stype:
            raise ValueError, "Cartan types must agree for automorphic branching rule"
        elif Rtype[0] == 'E' and r == 6:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif Rtype[0] == 'A':
            def rule(x) : y = [-i for i in x]; y.reverse(); return y
            return rule
        elif Rtype[0] == 'D':
            def rule(x) : x[len(x)-1] = -x[len(x)-1]; return x
            return rule
        elif Rtype[0] == 'E' and r == 6:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "No automorphism found"
    elif rule == "triality":
        if not Rtype == Stype:
            raise ValueError, "Triality is an automorphic type (for D4 only)"
        elif not Rtype[0] == 'D' and r == 4:
            raise ValueError, "Triality is for D4 only"
        else:
            return lambda x : [(x[0]+x[1]+x[2]+x[3])/2,(x[0]+x[1]-x[2]-x[3])/2,(x[0]-x[1]+x[2]-x[3])/2,(-x[0]+x[1]+x[2]-x[3])/2]
    elif rule == "symmetric":
        if Rtype[0] == 'A':
            if (Stype[0] == 'C' or Stype[0] == 'D' and r == 2*s-1) or (Stype[0] == 'B' and r == 2*s):
                return lambda x : [x[i]-x[r-i] for i in range(s)]
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'D' and Stype[0] == 'B' and s == r-1:
            return lambda x : x[:s]
        elif Rtype[0] == 'D' and r == 4 and Stype[0] == 'G':
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        elif Rtype[0] == 'E' and Stype[0] == 'F' and r == 6 and s == 4:
            return lambda x : [(x[4]-3*x[5])/2,(x[0]+x[1]+x[2]+x[3])/2,(-x[0]-x[1]+x[2]+x[3])/2,(-x[0]+x[1]-x[2]+x[3])/2]
        else:
            raise ValueError, "Rule not found"
    elif rule == "extended":
        if Stype.is_compound():
            if Rtype[0] in ['B','D'] and all(t[0] in ['B','D'] for t in stypes):
                if Rtype[0] == 'D':
                    rdeg = 2*r
                else:
                    rdeg = 2*r+1
                sdeg = 0
                for t in stypes:
                    if t[0] == 'D':
                        sdeg += 2*t[1]
                    else:
                        sdeg += 2*t[1]+1
                if rdeg == sdeg:
                    return lambda x : x[:s]
                else:
                    raise ValueError, "Rule not found"
            elif Rtype[0] == 'C'  and s == r:
                if all(t[0] == Rtype[0] for t in stypes):
                    return lambda x : x
            elif Rtype[0] == 'E' and s == r:
                if r == 6:
                    if stypes[0][0] == 'A' and stypes[0][1] == 5:
                        if stypes[1][0] == 'A' and stypes[1][1] == 1:
                            M = matrix(QQ,[(-3, -3, -3, -3, -3, -5, -5, 5), \
                                           (-9, 3, 3, 3, 3, 1, 1, -1), \
                                           (3, -9, 3, 3, 3, 1, 1, -1), \
                                           (3, 3, -9, 3, 3, 1, 1, -1), \
                                           (3, 3, 3, -9, 3, 1, 1, -1), \
                                           (3, 3, 3, 3, -9, 9, -3, 3), \
                                           (-3, -3, -3, -3, -3, -1, 11, 1), \
                                           (3, 3, 3, 3, 3, 1, 1, 11)])/12
                            return lambda x : tuple(M*vector(x))
                    if len(stypes) == 3 and all(x[0] == 'A' and x[1] == 2 for x in stypes):
                        M = matrix(QQ,[(0, 0, -2, -2, -2, -2, -2, 2), \
                                       (-3, 3, 1, 1, 1, 1, 1, -1), \
                                       (3, -3, 1, 1, 1, 1, 1, -1), \
                                       (0, 0, -2, -2, 4, 0, 0, 0), \
                                       (0, 0, -2, 4, -2, 0, 0, 0), \
                                       (0, 0, 4, -2, -2, 0, 0, 0), \
                                       (0, 0, -2, -2, -2, 2, 2, -2), \
                                       (3, 3, 1, 1, 1, -1, -1, 1), \
                                       (-3, -3, 1, 1, 1, -1, -1, 1)])/6
                        return lambda x : tuple(M*vector(x))
            elif Rtype[0] == 'F' and s == r:
                if stypes[0][0] == 'C' and stypes[0][1] == 3:
                    if stypes[1][0] == 'A' and stypes[1][1] == 1:
                        return lambda x : [x[0]-x[1],x[2]+x[3],x[2]-x[3],(-x[0]-x[1])/2,(x[0]+x[1])/2]
                if stypes[0][0] == 'A' and stypes[0][1] == 1:
                    if stypes[1][0] == 'C' and stypes[1][1] == 3:
                        return lambda x : [(-x[0]-x[1])/2,(x[0]+x[1])/2,x[0]-x[1],x[2]+x[3],x[2]-x[3]]
                if stypes[0][0] == 'A' and stypes[1][0] == 'A':
                    if stypes[0][1] == 2 and stypes[1][1] == 2:
                        M = matrix(QQ,[(-2, -1, -1, 0), (1, 2, -1, 0), (1, -1, 2, 0), (1, -1, -1, 3), (1, -1, -1, -3), (-2, 2, 2, 0)])/3
                    elif stypes[0][1] == 3 and stypes[1][1] == 1:
                        M = matrix(QQ,[(-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3), (2, -2, -2, -2), (-2, 2, 2, 2)])/4
                    elif stypes[0][1] == 1 and stypes[1][1] == 3:
                        M = matrix(QQ,[(2, -2, -2, -2), (-2, 2, 2, 2), (-3, -1, -1, -1), (1, 3, -1, -1), (1, -1, 3, -1), (1, -1, -1, 3)])/4
                    return lambda x : tuple(M*vector(x))
                else:
                    raise ValueError, "Rule not found"
            elif Rtype[0] == 'G' and s == r:
                if all(t[0] == 'A' and t[1] == 1 for t in stypes):
                    return lambda x : [(x[1]-x[2])/2,-(x[1]-x[2])/2, x[0]/2, -x[0]/2]
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'B' and Stype[0] == 'D' and s == r:
            return lambda x : x
        elif Rtype[0] == 'G' and Stype[0] == 'A' and s == r:
            return lambda x : [(x[0]-x[2])/3, (-x[1]+x[2])/3, (-x[0]+x[1])/3]
        elif Rtype[0] == 'F' and Stype[0] == 'B' and s == r:
            return lambda x : [-x[0], x[1], x[2], x[3]]
        elif Rtype[0] == 'E' and Stype[0] == 'E' and s == r and r >= 7:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "Rule not found"
    elif rule == "isomorphic":
        if r != s:
            raise ValueError, "Incompatible ranks"
        if Rtype == Stype:
            return lambda x : x
        elif Rtype[0] == 'B' and r == 2 and Stype[0] == 'C':
            def rule(x) : [x1, x2] = x; return [x1+x2, x1-x2]
            return rule
        elif Rtype[0] == 'B' and r == 1 and Stype[0] == 'A':
            return lambda x : [x[0],-x[0]]
        elif Rtype[0] == 'C' and r == 2 and Stype[0] == 'B':
            def rule(x) : [x1, x2] = x; return [(x1+x2)/2, (x1-x2)/2]
            return rule
        elif Rtype[0] == 'C' and r == 1 and Stype[0] == 'A':
            return lambda x : [x[0]/2,-x[0]/2]
        elif Rtype[0] == 'A' and r == 3 and Stype[0] == 'D':
            def rule(x): [x1, x2, x3, x4] = x; return [(x1+x2-x3-x4)/2, (x1-x2+x3-x4)/2, (x1-x2-x3+x4)/2]
            return rule
        elif Rtype[0] == 'D' and r == 2 and Stype.is_compound() and \
                 all(t[0] == 'A' for t in stypes):
            def rule(x): [t1, t2] = x; return [(t1-t2)/2, -(t1-t2)/2, (t1+t2)/2, -(t1+t2)/2]
            return rule
        elif Rtype[0] == 'D' and r == 3 and Stype[0] == 'A':
            def rule(x): [t1, t2, t3] = x; return [(t1+t2+t3)/2, (t1-t2-t3)/2, (-t1+t2-t3)/2, (-t1-t2+t3)/2]
            return rule
        else:
            raise ValueError, "Rule not found"
    elif rule == "tensor" or rule == "tensor-debug":
        if not Stype.is_compound():
            raise ValueError, "Tensor product requires more than one factor"
        if len(stypes) is not 2:
            raise ValueError, "Not implemented"
        if Rtype[0] is 'A':
            nr = Rtype[1]+1
        elif Rtype[0] is 'B':
            nr = 2*Rtype[1]+1
        elif Rtype[0] in ['C', 'D']:
            nr = 2*Rtype[1]
        else:
            raise ValueError, "Rule not found"
        [s1, s2] = [stypes[i][1] for i in range(2)]
        ns = [s1, s2]
        for i in range(2):
            if stypes[i][0] is 'A':
                ns[i] = ns[i]+1
            if stypes[i][0] is 'B':
                ns[i] = 2*ns[i]+1
            if stypes[i][0] in ['C','D']:
                ns[i] = 2*ns[i]
        if nr != ns[0]*ns[1]:
            raise ValueError, "Ranks don't agree with tensor product"
        if Rtype[0] == 'A':
            if all(t[0] == 'A' for t in stypes):
                def rule(x):
                    ret = [sum(x[i*ns[1]:(i+1)*ns[1]]) for i in range(ns[0])]
                    ret.extend([sum(x[ns[1]*j+i] for j in range(ns[0])) for i in range(ns[1])])
                    return ret
                return rule
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'B':
            if not all(t[0] == 'B' for t in stypes):
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'C':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'C':
                pass
            elif stypes[1][0] in ['B','D'] and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError, "Rule not found"
        elif Rtype[0] == 'D':
            if stypes[0][0] in ['B','D'] and stypes[1][0] is 'D':
                pass
            elif stypes[1][0] is 'B' and stypes[0][0] is 'D':
                pass
            elif stypes[1][0] is 'C' and stypes[0][0] is 'C':
                pass
            else:
                raise ValueError, "Rule not found"
        rows = []
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = 1
                rows.append(nextrow)
        if stypes[1][0] == 'B':
            for i in range(s1):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                rows.append(nextrow)
        for i in range(s1):
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[i] = 1
                nextrow[s1+j] = -1
                rows.append(nextrow)
        if stypes[0][0] == 'B':
            for j in range(s2):
                nextrow = (s1+s2)*[0]
                nextrow[s1+j] = 1
                rows.append(nextrow)
        mat = matrix(rows).transpose()
        if rule == "tensor-debug":
            print mat
        return lambda x : tuple(mat*vector(x))
    elif rule == "symmetric_power":
        if Stype[0] == 'A' and s == 1:
            if Rtype[0] == 'B':
                def rule(x):
                    a = sum((r-i)*x[i] for i in range(r))
                    return [a,-a]
                return rule
            elif Rtype[0] == 'C':
                def rule(x):
                    a = sum((2*r-2*i-1)*x[i] for i in range(r))
                    return [a/2,-a/2]
                return rule
            else:
                raise ValueError, "Rule not found"
        else:
            raise ValueError, "Rule not found"
    elif rule == "miscellaneous":
        if Rtype[0] == 'B' and Stype[0] == 'G' and r == 3:
            return lambda x : [x[0]+x[1], -x[1]+x[2], -x[0]-x[2]]
        else:
            raise ValueError, "Rule not found"


def branching_rule_from_plethysm(chi, cartan_type, return_matrix = False):
    """
    INPUT:

    - ``chi`` - the character of an irreducible representation pi of a group G
    - ``cartan_type`` - a classical Cartan type (A,B,C or D).

    It is assumed that the image of the irreducible representation pi
    naturally has its image in the group G.

    Returns a branching rule for this plethysm.

    EXAMPLE:

    The adjoint representation SL(3) --> GL(8) factors
    through SO(8). The branching rule in question will
    describe how representations of SO(8) composed with
    this homomorphism decompose into irreducible characters
    of SL(3).

        sage: A2 = WeylCharacterRing("A2")
        sage: A2 = WeylCharacterRing("A2", style="coroots")
        sage: ad = A2(1,1)
        sage: ad.degree()
         8
        sage: ad.frobenius_schur_indicator()
         1

    This confirms that ad has degree 8 and is orthogonal,
    hence factors through SO(8)=D4.

        sage: br = branching_rule_from_plethysm(ad,"D4")
        sage: D4 = WeylCharacterRing("D4")
        sage: [D4(f).branch(A2,rule = br) for f in D4.fundamental_weights()]
         [A2(1,1), A2(1,1) + A2(0,3) + A2(3,0), A2(1,1), A2(1,1)]
    """
    ct = CartanType(cartan_type)
    if ct[0] not in ["A","B","C","D"]:
        raise ValueError, "not implemented for type %s"%ct[0]
    if ct[0] is "A":
        ret = []
        for [v,n] in chi.mlist():
            ret.extend(n*[v.to_vector()])
        M = matrix(ret).transpose()
        if len(M.columns()) != ct[1] + 1:
            raise ValueError, "representation has wrong degree for type %s"%ct.__repr__()
        return lambda x : tuple(M*vector(x))
    if ct[0] in ["B","D"]:
        if chi.frobenius_schur_indicator() != 1:
            raise ValueError, "character is not orthogonal"
    if ct[0] is "C":
        if chi.frobenius_schur_indicator() != -1:
            raise ValueError, "character is not symplectic"
    if ct[0] is "B":
        if is_even(chi.degree()):
            raise ValueError, "degree is not odd"
    if ct[0] is ["C","D"]:
        if is_odd(chi.degree()):
            raise ValueError, "degree is not even"
    ret = []
    for [v,n] in chi.mlist():
        v = v.to_vector()
        if all(x==0 for x in v):
            if ct[0] is "B":
                n = (n-1)/2
            else:
                n = n/2
        elif [x for x in v if x !=0][0] < 0:
            continue
        ret.extend(n*[v])
    M = matrix(ret).transpose()
    if len(M.columns()) != ct.root_system().ambient_space().dimension():
        raise ValueError, "representation has wrong degree for type %s"%ct.__repr__()
    if return_matrix:
        return M
    else:
        return lambda x : tuple(M*vector(x))


class WeightRingElement(AlgebraElement):
    """
    A class for weights, and linear combinations of weights. See
    WeightRing? for more information.
    """
    def __init__(self, A, mdict):
        """
        INPUT:


        -  ``A`` - the Weight Ring

        -  ``mdict`` - dictionary of weight multiplicities


        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: a2 == loads(dumps(a2))
            True
        """
        AlgebraElement.__init__(self, A)
        self._mdict = mdict
        self._parent = A
        self._prefix = A._prefix
        self._space = A._space

    def _wt_repr(self, k):
        """
        Returns the representation of a single weight.

        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: a2([1,0,0])
            a2(1,0,0)
        """
        hstring = str(k[0])
        for i in range(1,self._space.n):
            hstring = hstring+","+str(k[i])
        return self._prefix+"("+hstring+")"

    def __repr__(self):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: b3 = WeightRing(B3)
            sage: fw = B3.fundamental_weights()
            sage: b3(fw[3])
            b3(1/2,1/2,1/2)
            sage: b3(B3(fw[3]))
            b3(-1/2,-1/2,-1/2) + b3(-1/2,-1/2,1/2) + b3(-1/2,1/2,-1/2) + b3(-1/2,1/2,1/2) + b3(1/2,-1/2,-1/2) + b3(1/2,-1/2,1/2) + b3(1/2,1/2,-1/2) + b3(1/2,1/2,1/2)
        """

        if self._mdict == {}:
            return "0"
        v = self._mdict.keys()
        v.sort(key = lambda v: tuple(v.to_vector()))
        return repr_lincomb([self._wt_repr(k) for k in v], [self._mdict[k] for k in v])

    def __cmp__(left, right):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: b3 = WeightRing(B3)
            sage: fw = [b3(w) for w in B3.fundamental_weights()]
            sage: sorted(fw)
            [b3(1/2,1/2,1/2), b3(1,0,0), b3(1,1,0)]
        """
        return left._mdict.__cmp__(right._mdict)

    def _add_(self, y):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: [a2(0,0,0)+a2(2,1,0), a2(2,1,0)+a2(0,0,0), - a2(0,0,0)+2*a2(0,0,0), -2*a2(0,0,0)+a2(0,0,0), -a2(2,1,0)+2*a2(2,1,0)-a2(2,1,0)]
            [a2(0,0,0) + a2(2,1,0), a2(0,0,0) + a2(2,1,0), a2(0,0,0), -a2(0,0,0), 0]
        """
        mdict = {}
        for k in self._mdict.keys():
            if k in y._mdict:
                mdict[k] = self._mdict[k] + y._mdict[k]
                if mdict[k] == 0:
                    del mdict[k]
            else:
                mdict[k] = self._mdict[k]
        for k in y._mdict:
            if not k in self._mdict:
                mdict[k] = y._mdict[k]
        return WeightRingElement(self._parent, mdict)

    def _neg_(self):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: [-x for x in [a2(0,0,0), 2*a2(0,0,0), -a2(0,0,0), -2*a2(0,0,0)]]
            [-a2(0,0,0), -2*a2(0,0,0), a2(0,0,0), 2*a2(0,0,0)]
        """
        mdict = self._mdict.copy()
        for k in self._mdict:
            mdict[k] = - mdict[k]
        return WeightRingElement(self._parent, mdict)

    def _sub_(self, y):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: chi = a2(0,0,0)+2*a2(1,0,0)+3*a2(2,0,0)
            sage: mu =  3*a2(0,0,0)+2*a2(1,0,0)+a2(2,0,0)
            sage: chi - mu
            -2*a2(0,0,0) + 2*a2(2,0,0)
        """
        return self._add_(y._neg_())

    def _mul_(self, y):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: [chi, mu] = [A2(1,0,0), A2(1,1,0)]
            sage: chi*mu
            A2(1,1,1) + A2(2,1,0)
            sage: a2(chi)*a2(mu)
            a2(0,1,2) + a2(0,2,1) + a2(1,0,2) + 3*a2(1,1,1) + a2(1,2,0) + a2(2,0,1) + a2(2,1,0)
            sage: (a2(chi)*a2(mu)).character() == chi*mu
            True
        """
        mdict = {}
        for k in self._mdict:
            for l in y._mdict:
                m = k+l
                if m in mdict:
                    mdict[m] += self._mdict[k]*y._mdict[l]
                else:
                    mdict[m] = self._mdict[k]*y._mdict[l]

        for k in mdict.keys():
            if mdict[k] == 0:
                del mdict[k]
        return WeightRingElement(self._parent, mdict)

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: chi = A2([2,1,0])
            sage: chi^4
            8*A2(4,4,4) + 32*A2(5,4,3) + 20*A2(5,5,2) + 20*A2(6,3,3) + 33*A2(6,4,2) + 15*A2(6,5,1) + 2*A2(6,6,0) + 15*A2(7,3,2) + 12*A2(7,4,1) + 3*A2(7,5,0) + 2*A2(8,2,2) + 3*A2(8,3,1) + A2(8,4,0)
            sage: a2(chi)^4
            a2(0,4,8) + 4*a2(0,5,7) + 6*a2(0,6,6) + 4*a2(0,7,5) + a2(0,8,4) + 4*a2(1,3,8) + 20*a2(1,4,7) + 40*a2(1,5,6) + 40*a2(1,6,5) + 20*a2(1,7,4) + 4*a2(1,8,3) + 6*a2(2,2,8) + 40*a2(2,3,7) + 106*a2(2,4,6) + 144*a2(2,5,5) + 106*a2(2,6,4) + 40*a2(2,7,3) + 6*a2(2,8,2) + 4*a2(3,1,8) + 40*a2(3,2,7) + 144*a2(3,3,6) + 260*a2(3,4,5) + 260*a2(3,5,4) + 144*a2(3,6,3) + 40*a2(3,7,2) + 4*a2(3,8,1) + a2(4,0,8) + 20*a2(4,1,7) + 106*a2(4,2,6) + 260*a2(4,3,5) + 346*a2(4,4,4) + 260*a2(4,5,3) + 106*a2(4,6,2) + 20*a2(4,7,1) + a2(4,8,0) + 4*a2(5,0,7) + 40*a2(5,1,6) + 144*a2(5,2,5) + 260*a2(5,3,4) + 260*a2(5,4,3) + 144*a2(5,5,2) + 40*a2(5,6,1) + 4*a2(5,7,0) + 6*a2(6,0,6) + 40*a2(6,1,5) + 106*a2(6,2,4) + 144*a2(6,3,3) + 106*a2(6,4,2) + 40*a2(6,5,1) + 6*a2(6,6,0) + 4*a2(7,0,5) + 20*a2(7,1,4) + 40*a2(7,2,3) + 40*a2(7,3,2) + 20*a2(7,4,1) + 4*a2(7,5,0) + a2(8,0,4) + 4*a2(8,1,3) + 6*a2(8,2,2) + 4*a2(8,3,1) + a2(8,4,0)
            sage: (a2(chi)^4).character() == chi^4
            True
        """
        if not n in ZZ:
            raise TypeError, "exponent must be an integer"
        if not n >= 0:
            raise TypeError, "exponent must be nonnegative"
        z = self._parent.__call__(1)
        for _ in range(n):
            z *= self
        return z

    def cartan_type(self):
        """
        Returns the Cartan Type.

        EXAMPLES::

            sage: A2=WeylCharacterRing("A2")
            sage: a2 = WeightRing(A2)
            sage: a2([0,1,0]).cartan_type()
            ['A', 2]
        """
        return self._parent._cartan_type

    def mlist(self):
        """
        Returns a list of weights in self with their multiplicities.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: pr = sum(g2(a) for a in G2.positive_roots())
            sage: sorted(pr.mlist())
            [[(1, -2, 1), 1],  [(1, -1, 0), 1],  [(1, 1, -2), 1],  [(1, 0, -1), 1],  [(2, -1, -1), 1],  [(0, 1, -1), 1]]
        """
        return [[k,m] for k,m in self._mdict.iteritems()]

    def weyl_group_action(self, w):
        """
        Returns the actionof the Weyl group element w on self.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: L = g2.space()
            sage: [fw1, fw2] = L.fundamental_weights()
            sage: sum(g2(fw2).weyl_group_action(w) for w in L.weyl_group())
            2*g2(-2,1,1) + 2*g2(-1,-1,2) + 2*g2(-1,2,-1) + 2*g2(1,-2,1) + 2*g2(1,1,-2) + 2*g2(2,-1,-1)
        """
        return WeightRingElement(self._parent, dict([[w.action(x),m] for x,m in self._mdict.iteritems()]))

    def character(self):
        """
        Assuming that self is invariant under the Weyl group, this will
        express it as a linear combination of characters. If self is not
        Weyl group invariant, this method will not terminate.

        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: W = a2.space().weyl_group()
            sage: mu = a2(2,1,0)
            sage: nu = sum(mu.weyl_group_action(w) for w in W)
            sage: nu
            a2(0,1,2) + a2(0,2,1) + a2(1,0,2) + a2(1,2,0) + a2(2,0,1) + a2(2,1,0)
            sage: nu.character()
            -2*A2(1,1,1) + A2(2,1,0)
        """
        return WeylCharacter(self._parent._parent, self._parent._parent.char_from_weights(self._mdict), self._mdict)


class WeightRing(Algebra):
    def __init__(self, A, prefix=None):
        """
        Creates an auxiliary ring for the weights of a representation. This
        ring is associated with a WeylCharacterRing. The weights of a
        representation of a Lie group G or algebra are the characters of
        the characters of a Cartan subgroup T or subalgebra that occur when
        the representation is restricted to that Cartan. If a linear
        combination of weights is invariant under the Weyl group, then it
        is a linear combination of characters of G.

        As with the WeylCharacterRing, you may want to make sure that the
        prefix matches the name that you assign the ring, so that the
        __call__ method can parse the ring's own output. If you do not
        assign a prefix, one is automatically generated by changing the
        case of the prefix to the associated WeylCharacterRing, from upper
        case to lower case. However you may assign your own prefixes to
        both rings.

        INPUT:


        -  ``A`` - a WeylCharacterRing.


        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: wd = prod(a2(x/2)-a2(-x/2) for x in a2.space().positive_roots()); wd
            -a2(-1,0,1) + a2(-1,1,0) + a2(0,-1,1) - a2(0,1,-1) - a2(1,-1,0) + a2(1,0,-1)
            sage: chi = A2([5,3,0]); chi
            A2(5,3,0)
            sage: a2(chi)
            a2(0,3,5) + a2(0,4,4) + a2(0,5,3) + a2(1,2,5) + 2*a2(1,3,4) + 2*a2(1,4,3) + a2(1,5,2) +
            a2(2,1,5) + 2*a2(2,2,4) + 3*a2(2,3,3) + 2*a2(2,4,2) + a2(2,5,1) + a2(3,0,5) + 2*a2(3,1,4) +
            3*a2(3,2,3) + 3*a2(3,3,2) + 2*a2(3,4,1) + a2(3,5,0) + a2(4,0,4) + 2*a2(4,1,3) + 2*a2(4,2,2) +
            2*a2(4,3,1) + a2(4,4,0) + a2(5,0,3) + a2(5,1,2) + a2(5,2,1) + a2(5,3,0)
            sage: a2(chi)*wd
            -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
            sage: sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
            -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
            sage: a2(chi)*wd == sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.space().weyl_group())
            True

        In the above example, we create a WeylCharacterRing A2 whose
        objects are the characters. Attached to this is a second WeightRing
        called a2. The Weyl denominator is created and labeled wd. This is
        the product of factors a2(alpha/2)-a2(-alpha/2) where alpha runs
        through the positive roots. Then character chi with highest weight
        [5,3,0] is created and coerced into a2. It has many terms. It is
        multiplied by wd and compared with the alternating sum of its
        images under Weyl group elements. These are equal, illustrating the
        Weyl character formula.

        EXAMPLES::

            sage: R = WeylCharacterRing(['G',2], prefix = "R", base_ring = QQ)
            sage: S = WeightRing(R, prefix = "S")
            sage: L = S.space()
            sage: vc = [sum(S(1/2)*S(f).weyl_group_action(w) for w in L.weyl_group()) for f in L.fundamental_weights()]; vc
            [S(-1,0,1) + S(-1,1,0) + S(0,-1,1) + S(0,1,-1) + S(1,-1,0) + S(1,0,-1),
            S(-2,1,1) + S(-1,-1,2) + S(-1,2,-1) + S(1,-2,1) + S(1,1,-2) + S(2,-1,-1)]
            sage: [v.character() for v in vc]
            [-R(0,0,0) + R(1,0,-1), -R(0,0,0) - R(1,0,-1) + R(2,-1,-1)]

        TESTS::

            sage: R = WeylCharacterRing(['G',2], prefix = "R", base_ring = QQ)
            sage: S = WeightRing(R, prefix = "S")
            sage: S == loads(dumps(S))
            True
        """
        self._parent = A
        self._cartan_type = self._parent._cartan_type
        if prefix is None:
            if self._parent._prefix.isupper():
                prefix = self._parent._prefix.lower()
            elif self._parent._prefix.islower():
                prefix = self._parent._prefix.upper()
            else:
                prefix = (self._cartan_type[0].lower()+str(self._rank))
        self._base_ring = self._parent._base_ring
        self._space = self._parent._space
        self._origin = self._parent._origin
        self._prefix = prefix

    def __call__(self, *args):
        """
        Coerces the element into the ring.

        INPUT:


        -  ``x`` - a ring element


        EXAMPLES::

            sage: a2 = WeightRing(WeylCharacterRing(['A',2]))
            sage: a2(-1)
            -a2(0,0,0)
        """
        if len(args) == 1:
            x = args[0]
        else:
            x = args
        if x == 0:
            return WeightRingElement(self, {})
        if x in ZZ:
            mdict = {self._origin: x}
            return WeightRingElement(self, mdict)
        if is_Element(x):
            P = x.parent()
            if P is self:
                return x
            elif x in self.base_ring():
                mdict = {self._origin: x}
                return WeightRingElement(self, mdict)
        try:
            if x.parent() == self._parent:
                return WeightRingElement(self, x._mdict)
        except AttributeError:
            pass
        x = self._space(x)
        mdict = {x: 1}
        return WeightRingElement(self, mdict)

    def __repr__(self):
        """
        EXAMPLES::

            sage: P.<q>=QQ[]
            sage: G2 = WeylCharacterRing(['G',2], base_ring = P)
            sage: WeightRing(G2)
            The Weight ring attached to The Weyl Character Ring of Type ['G', 2] with Univariate Polynomial Ring in q over Rational Field coefficients
        """
        return "The Weight ring attached to %s"%self._parent.__repr__()

    def base_ring(self):
        """
        Returns the base ring.

        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3], base_ring = CC); R.base_ring()
            Complex Field with 53 bits of precision
        """
        return self._base_ring

    def _coerce_impl(self, x):
        """
        Coercion from the base ring.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: 2 in g2
            True
        """
        if x in self._base_ring:
            return self.__call__(x)
        raise TypeError, "no canonical coercion of x"

    def __cmp__(self, x):
        """
        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: A3 = WeylCharacterRing(['A', 3])
            sage: a3 = WeightRing(A3)
            sage: e8 == e8
            True
            sage: e8 == a3
            False
        """
        return cmp(repr(self), repr(x))

    def cartan_type(self):
        """
        Returns the Cartan Type.

        EXAMPLES::

            sage: A2 = WeylCharacterRing("A2")
            sage: WeightRing(A2).cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def space(self):
        """
        Returns the weight space realization associated to self.

        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: e8.space()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._space

    def fundamental_weights(self):
        """
        Returns the fundamental weights.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).fundamental_weights()
            Finite family {1: (1, 0, -1), 2: (2, -1, -1)}
        """
        return self._space.fundamental_weights()

    def simple_roots(self):
        """
        Returns the simple roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).simple_roots()
            Finite family {1: (0, 1, -1), 2: (1, -2, 1)}
        """
        return self._space.simple_roots()

    def positive_roots(self):
        """
        Returns the positive roots.

        EXAMPLES::

            sage: WeightRing(WeylCharacterRing("G2")).positive_roots()
            [(0, 1, -1), (1, -2, 1), (1, -1, 0), (1, 0, -1), (1, 1, -2), (2, -1, -1)]
        """
        return self._space.positive_roots()

    def wt_repr(self, wt):
        """
        Returns a string representing the irreducible character with
        highest weight vectr wt.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: g2.wt_repr([1,0,0])
            'g2(1,0,0)'
        """
        hstring = str(wt[0])
        for i in range(1,self._space.n):
            hstring=hstring+","+str(wt[i])
        return self._prefix+"("+hstring+")"


