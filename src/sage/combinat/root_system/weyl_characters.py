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
from sage.modules.free_module import VectorSpace
from sage.structure.element import is_Element
from sage.rings.all import ZZ, QQ
from sage.misc.misc import repr_lincomb
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

    There is associated with K, L or G as above a lattice, the weight
    lattice, whose elements (called weights) are characters of a Cartan
    subgroup or subalgebra. There is an action of the Weyl group W on
    the lattice, and elements of a fixed fundamental domain for W, the
    positive Weyl chamber, are called dominant. There is for each
    representation a unique highest dominant weight that occurs with
    nonzero multiplicity with respect to a certain partial order, and
    it is called the highest weight vector.

    EXAMPLES::

        sage: L = RootSystem(['A',2]).ambient_space()
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

            sage: R = WeylCharacterRing(['B',3], prefix = "R")
            sage: r =  R(1,1,0)
            sage: r == loads(dumps(r))
            True
        """
        AlgebraElement.__init__(self, A)
        self._hdict = hdict
        self._mdict = mdict
        self._parent = A
        self._lattice = A._lattice

    def __repr__(self):
        """
        EXAMPLES::

            sage: R = WeylCharacterRing(['B',3], prefix = "R")
            sage: [R(w) for w in R.lattice().fundamental_weights()]
            [R(1,0,0), R(1,1,0), R(1/2,1/2,1/2)]
        """
        if self._hdict == {}:
            return "0"
        v = self._hdict.keys()
        # Just a workaround to keep the same sorting as before when
        # the dictionary was indexed by tuples
        v.sort(key = lambda v: tuple(v.to_vector()))
        return repr_lincomb([self._parent.irr_repr(k) for k in v], [self._hdict[k] for k in v])

    def __cmp__(self, right):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: fw = [B3(w) for w in B3.lattice().fundamental_weights()]
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

    def degree(self):
        """
        The degree of the character, that is, the dimension of module.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: [B3(x).degree() for x in B3.lattice().fundamental_weights()]
            [7, 21, 8]
        """
        return sum(self._mdict[k] for k in self._mdict)

    def check(self, verbose=False):
        """
        To check the correctness of an element, we compare the theoretical
        dimension computed Weyl character formula with the actual one
        obtained by adding up the multiplicities.

        EXAMPLES::

            sage: F4 = WeylCharacterRing(['F',4])
            sage: [F4(x).check(verbose = true) for x in F4.lattice().fundamental_weights()]
            [[52, 52], [1274, 1274], [273, 273], [26, 26]]
        """
        theoretical = sum(self._hdict[k]*self._lattice.weyl_dimension(k) for k in self._hdict)
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
            sage: [B3(w).branch(A2,rule="levi") for w in B3.lattice().fundamental_weights()]
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

    def mlist(self):
        """
        Returns a list of weights in self with their multiplicities.

        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: B3(1/2,1/2,1/2).mlist()
            [[(1/2, -1/2, -1/2), 1], [(-1/2, 1/2, -1/2), 1], [(1/2, 1/2, 1/2), 1], [(1/2, 1/2, -1/2), 1], [(-1/2, -1/2, 1/2), 1], [(-1/2, -1/2, -1/2), 1], [(1/2, -1/2, 1/2), 1], [(-1/2, 1/2, 1/2), 1]]
        """
# Why did the test not pass with the following indentation?
#             [[( 1/2, -1/2, -1/2), 1],
#              [(-1/2,  1/2, -1/2), 1],
#              [( 1/2,  1/2,  1/2), 1],
#              [( 1/2,  1/2, -1/2), 1],
#              [(-1/2, -1/2,  1/2), 1],
#              [(-1/2, -1/2, -1/2), 1],
#              [( 1/2, -1/2,  1/2), 1],
#              [(-1/2,  1/2,  1/2), 1]]
        return [[k,m] for k,m in self._mdict.iteritems()]

    def parent(self):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: B3(2).parent()
            The Weyl Character Ring of Type [B,3] with Integer Ring coefficients
        """
        return self._parent


def WeylCharacterRing(ct, base_ring=ZZ, prefix=None, cache=False):
    r"""
    A class for rings of Weyl characters. The Weyl character is a
    character of a semisimple (or reductive) Lie group or algebra. They
    form a ring, in which the addition and multiplication correspond to
    direct sum and tensor product of representations.

    INPUT:

    - ``ct`` - The Cartan Type

    OPTIONAL ARGUMENTS:

    - ``base_ring`` -  (default: `\mathbb{Z}`)

    - ``prefix`` (default an automatically generated prefix
      based on Cartan type)

    - ``cache`` -  (default False) setting cache = True is a substantial
      speedup at the expense of some memory.

    If no prefix specified, one is generated based on the Cartan type.
    It is good to name the ring after the prefix, since then it can
    parse its own output.

    EXAMPLES::

        sage: G2 = WeylCharacterRing(['G',2])
        sage: [fw1,fw2] = G2.lattice().fundamental_weights()
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
        sage: chi = R(R.lattice().fundamental_weights()[3]); chi
        R(1/2,1/2,1/2)
        sage: R(1/2,1/2,1/2) == chi
        True

    The multiplication in R corresponds to the product of characters,
    which you can use to determine the decomposition of tensor products
    into irreducibles. For example, let us compute the tensor product
    of the standard and spin representations of Spin(7).

    EXAMPLES::

        sage: B3 = WeylCharacterRing(['B',3])
        sage: [fw1,fw2,fw3]=B3.lattice().fundamental_weights()
        sage: [B3(fw1).degree(),B3(fw3).degree()]
        [7, 8]
        sage: B3(fw1)*B3(fw3)
        B3(1/2,1/2,1/2) + B3(3/2,1/2,1/2)

    The name of the irreducible representation encodes the highest
    weight vector.

    TESTS::

        sage: F4 = WeylCharacterRing(['F',4], cache = True)
        sage: [fw1,fw2,fw3,fw4] = F4.lattice().fundamental_weights()
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
        sage: R([2,1,0]).mlist()
        [[(1, 0, 2), 1],
         [(1, 1, 1), 2],
         [(2, 1, 0), 1],
         [(1, 2, 0), 1],
         [(2, 0, 1), 1],
         [(0, 1, 2), 1],
         [(0, 2, 1), 1]]
    """
    ct = cartan_type.CartanType(ct)
    return cache_wcr(ct, base_ring=base_ring, prefix=prefix, cache=cache)

# TODO: inherit all the data structure from CombinatorialFreeModule(base_ring, self._lattice)

class WeylCharacterRing_class(Algebra):
    def __init__(self, ct, base_ring, prefix, cache):
        """
        EXAMPLES::

            sage: R = WeylCharacterRing(['A',3])
            sage: R == loads(dumps(R))
            True
        """
        sage.structure.parent_base.ParentWithBase.__init__(self, base_ring)

        self.cartan_type = ct
        self._base_ring = base_ring
        self._lattice = RootSystem(self.cartan_type).ambient_space()
        self._origin = self._lattice.zero()
        if prefix == None:
            prefix = ct[0]+str(ct[1])
        self._prefix = prefix
        alpha = self._lattice.simple_roots()
        Lambda = self._lattice.fundamental_weights()
        # FIXME: indexing of fundamental weights
        self._ip = [Lambda[i].inner_product(alpha[i])
                    for i in ct.index_set()]
        self._cache = cache
        if cache:
            self._irreducibles={}


    def __call__(self, *args):
        """
        Coerces the element into the ring. You may pass a vector in the
        ambient lattice, an element of the base_ring, or an argument list
        of integers (or half-integers for the spin types) which are the
        components of a vector in the ambient lattice.

        INPUT:


        -  ``x`` - a ring element to be coerced; or

        -  ``*args`` - the components of a vector


        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
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

        if x == 0:
            return WeylCharacter(self, {}, {})

        if x in ZZ:
            hdict = {self._origin: x}
            return WeylCharacter(self, hdict, hdict)

        if is_Element(x):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return WeylCharacter(self, x._hdict, x._mdict)
            elif x in self.base_ring():
                hdict = {self._origin: x}
                return WeylCharacter(self, hdict, hdict)

        x = self._lattice(x)

        alphacheck = self._lattice.simple_coroots()
        vp = [x.inner_product(alphacheck[i])
              for i in self.cartan_type.index_set()]
        if not all(v in ZZ for v in vp):
            raise ValueError, "not in weight lattice"
        if not all(v >= 0 for v in vp):
            raise ValueError, "the weight%s is not dominant"%x.__repr__()
        if self._cache and x in self._irreducibles:
            return self._irreducibles[x]
        hdict = {x: 1}
        mdict = irreducible_character_freudenthal(x, self._lattice)
        ret = WeylCharacter(self, hdict, mdict)
        if self._cache:
            self._irreducibles[x] = ret
        return ret

    def __repr__(self):
        """
        EXAMPLES::

            sage: WeylCharacterRing(['A',3])
            The Weyl Character Ring of Type [A,3] with Integer Ring coefficients
        """
        return "The Weyl Character Ring of Type [%s,%d] with %s coefficients"%(self.cartan_type[0], self.cartan_type[1], self._base_ring.__repr__())

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

    # FIXME: should be something like weight_lattice_realization?
    def lattice(self):
        """
        Returns the weight lattice associated to self.

        EXAMPLES::

            sage: WeylCharacterRing(['E',8]).lattice()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._lattice

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
            highest = max((x.inner_product(self._lattice.rho()),x) for x in ddict)[1]
            if not highest.is_dominant():
                raise ValueError, "multiplicity dictionary may not be Weyl group invariant"
            if self._cache and highest in self._irreducibles:
                sdict = self._irreducibles[highest]._mdict
            else:
                sdict = irreducible_character_freudenthal(highest, self._lattice)
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

            sage: A2 = WeylCharacterRing(['A',2])
            sage: A2.irr_repr([2,1,0])
            'A2(2,1,0)'
        """
        hstring = str(hwv[0])
        for i in range(1,self._lattice.n):
            hstring=hstring+","+str(hwv[i])
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

    - ``L`` - the ambient lattice
    """

    rho = L.rho()
    mdict = {}
    current_layer = {hwv:1}
    while len(current_layer) > 0:
        next_layer = {}
        for mu in current_layer:
            if not current_layer[mu] == 0:
                mdict[mu] = current_layer[mu]
                for alpha in L.simple_roots():
                    next_layer[mu-alpha] = None
        if debug:
            print next_layer
        for mu in next_layer:
            if next_layer[mu] == None:
                if debug:
                    print "  mu:", mu
                accum = 0
                for alpha in L.positive_roots():
                    i = 1
                    while mu+i*alpha in mdict:
                        if debug:
                            print "    ", mu+i*alpha,
                            print mdict[mu + i*alpha]*(mu + i*alpha).inner_product(alpha)
                        accum += mdict[mu + i*alpha]*(mu + i*alpha).inner_product(alpha)
                        i += 1
                if accum == 0:
                    next_layer[mu] = 0
                else:
                    next_layer[mu] = QQ(2*accum)/QQ((hwv+rho).inner_product(hwv+rho)-(mu+rho).inner_product(mu+rho))
        current_layer = next_layer
    return mdict

def branch_weyl_character(chi, R, S, rule="default"):
    r"""
    A Branching rule describes the restriction of representations from
    a Lie group or algebra G to a smaller one. See for example, R. C.
    King, Branching rules for classical Lie groups using tensor and
    spinor methods. J. Phys. A 8 (1975), 429-449 or Howe, Tan and
    Willenbring, Stable branching rules for classical symmetric pairs,
    Trans. Amer. Math. Soc. 357 (2005), no. 4, 1601-1626.

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
    rules for cases that we have omitted. (Rules for the exceptional
    groups have not yet been coded.)

    To explain the predefined rules we survey the most important
    branching rules. These may be classified into several cases, and
    once this is understood, the detailed classification can be read
    off from the Dynkin diagrams. Dynkin classified the maximal
    subgroups of Lie groups in Mat. Sbornik N.S. 30(72):349-462
    (1952).

    We will list these for the cases where the Dynkin diagram of S is
    connected. This excludes branching rules such as A3 -> A1 x A1,
    which are not yet implemented.

    LEVI TYPE. These can be read off from the Dynkin diagram. If
    removing a node from the Dynkin diagram produces another Dynkin
    diagram, there is a branching rule. Currently we require that the
    smaller diagram be connected. For these rules use the option
    rule="levi"::

       ['A',r] -> ['A',r-1]
       ['B',r] -> ['A',r-1]
       ['B',r] -> ['B',r-1]
       ['C',r] -> ['A',r-1]
       ['C',r] -> ['C',r-1]
       ['D',r] -> ['A',r-1]
       ['D',r] -> ['D',r-1]
       ['E',r] -> ['A',r-1] r = 6,7,8 (not implemented yet)
       ['E',r] -> ['D',r-1] r = 6,7,8 (not implemented yet)
       ['E',r] -> ['E',r-1] r = 6,7 (not implemented yet)
       ['F',4] -> ['B',3] (not implemented yet)
       ['F',4] -> ['C',3] (not implemented yet)
       ['G',2] -> ['A',1] (short root) (not implemented yet)

    The other Levi branching rule from `G_2` -> `A_1` corresponding to the
    long root is available by first branching `G_2` -> `A_2` then branching to
    A1.

    AUTOMORPHIC TYPE. If the Dynkin diagram has a symmetry, then there
    is an automorphism that is a special case of a branching rule.
    There is also an exotic"triality" automorphism of D4 having order
    3. Use rule="automorphic" or rule="triality"

    ['A',r] - ['A',r] ['D',r] - ['D',r] ['E',6] - ['E',6] (not
    implemented yet)

    SYMMETRIC TYPE. Related to the automorphic type, when the Dynkin
    diagram has a symmetry there is a branching rule to the subalgebra
    (or subgroup) of invariants under the automorphism. Use
    rule="symmetric".

    ['A',2r+1] - ['B',r] ['A',2r] - ['C',r] ['D',r] - ['B',r-1] ['E',6]
    - ['F',4] (not implemented yet) ['D',4] - ['G',2] (not implemented
    yet)

    EXTENDED TYPE. If removing a node from the extended Dynkin diagram
    results in a Dynkin diagram, then there is a branching rule. Use
    rule="extended" for these.

    ['G',2] - ['A',2] (not implemented yet) ['B',r] - ['D',r] ['F',4] -
    ['B',4] (not implemented yet) ['E',7] - ['A',7] (not implemented
    yet) ['E',8] - ['A',8] (not implemented yet)

    MISCELLANEOUS: Use rule="miscellaneous" for the following rule,
    which does not fit into the above framework.

    ['B',3] - ['G',2] (Not implemented yet.)

    ISOMORPHIC TYPE: Although not usually referred to as a branching
    rule, the effects of the accidental isomorphisms may be handled
    using rule="isomorphic"

    ['B',2] - ['C',2] ['C',2] - ['B',2] ['A',3] - ['D',3] ['D',3] -
    ['A',3]

    EXAMPLES: (Levi type)

    ::

        sage: A2 = WeylCharacterRing(['A',2])
        sage: B2 = WeylCharacterRing(['B',2])
        sage: C2 = WeylCharacterRing(['C',2])
        sage: A3 = WeylCharacterRing(['A',3])
        sage: B3 = WeylCharacterRing(['B',3])
        sage: C3 = WeylCharacterRing(['C',3])
        sage: D3 = WeylCharacterRing(['D',3])
        sage: A4 = WeylCharacterRing(['A',4])
        sage: D4 = WeylCharacterRing(['D',4])
        sage: A5 = WeylCharacterRing(['A',5])
        sage: D5 = WeylCharacterRing(['D',5])
        sage: [B3(w).branch(A2,rule="levi") for w in B3.lattice().fundamental_weights()]
        [A2(0,0,-1) + A2(0,0,0) + A2(1,0,0),
         A2(0,-1,-1) + A2(0,0,-1) + A2(0,0,0) + A2(1,0,-1) + A2(1,0,0) + A2(1,1,0),
         A2(-1/2,-1/2,-1/2) + A2(1/2,-1/2,-1/2) + A2(1/2,1/2,-1/2) + A2(1/2,1/2,1/2)]

    The last example must be understood as follows. The representation
    of B3 being branched is spin, which is not a representation of
    SO(7) but of its double cover spin(7). The group A2 is really GL(3)
    and the double cover of SO(7) induces a cover of GL(3) that is
    trivial over SL(3) but not over the center of GL(3). The weight
    lattice for this GL(3) consists of triples (a,b,c) of half integers
    such that a-b and b-c are in `\mathbb{Z}`, and this is reflected in the last
    decomposition.

    ::

        sage: [C3(w).branch(A2,rule="levi") for w in C3.lattice().fundamental_weights()]
        [A2(0,0,-1) + A2(1,0,0),
         A2(0,-1,-1) + A2(1,0,-1) + A2(1,1,0),
         A2(-1,-1,-1) + A2(1,-1,-1) + A2(1,1,-1) + A2(1,1,1)]
        sage: [D4(w).branch(A3,rule="levi") for w in D4.lattice().fundamental_weights()]
        [A3(0,0,0,-1) + A3(1,0,0,0),
         A3(0,0,-1,-1) + A3(0,0,0,0) + A3(1,0,0,-1) + A3(1,1,0,0),
         A3(1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,-1/2),
         A3(-1/2,-1/2,-1/2,-1/2) + A3(1/2,1/2,-1/2,-1/2) + A3(1/2,1/2,1/2,1/2)]
        sage: [B3(w).branch(B2,rule="levi") for w in B3.lattice().fundamental_weights()]
        [2*B2(0,0) + B2(1,0), B2(0,0) + 2*B2(1,0) + B2(1,1), 2*B2(1/2,1/2)]
        sage: C3 = WeylCharacterRing(['C',3])
        sage: [C3(w).branch(C2,rule="levi") for w in C3.lattice().fundamental_weights()]
        [2*C2(0,0) + C2(1,0),
         C2(0,0) + 2*C2(1,0) + C2(1,1),
         C2(1,0) + 2*C2(1,1)]
        sage: [D5(w).branch(D4,rule="levi") for w in D5.lattice().fundamental_weights()]
        [2*D4(0,0,0,0) + D4(1,0,0,0),
         D4(0,0,0,0) + 2*D4(1,0,0,0) + D4(1,1,0,0),
         D4(1,0,0,0) + 2*D4(1,1,0,0) + D4(1,1,1,0),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2),
         D4(1/2,1/2,1/2,-1/2) + D4(1/2,1/2,1/2,1/2)]

    EXAMPLES: (Automorphic type, including D4 triality)

    ::

        sage: [A3(chi).branch(A3,rule="automorphic") for chi in A3.lattice().fundamental_weights()]
        [A3(0,0,0,-1), A3(0,0,-1,-1), A3(0,-1,-1,-1)]
        sage: [D4(chi).branch(D4,rule="automorphic") for chi in D4.lattice().fundamental_weights()]
        [D4(1,0,0,0), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1/2,1/2,1/2,-1/2)]
        sage: [D4(chi).branch(D4,rule="triality") for chi in D4.lattice().fundamental_weights()]
        [D4(1/2,1/2,1/2,-1/2), D4(1,1,0,0), D4(1/2,1/2,1/2,1/2), D4(1,0,0,0)]

    EXAMPLES: (Symmetric type)

    ::

        sage: [w.branch(B2,rule="symmetric") for w in [A4(1,0,0,0,0),A4(1,1,0,0,0),A4(1,1,1,0,0),A4(2,0,0,0,0)]]
        [B2(1,0), B2(1,1), B2(1,1), B2(0,0) + B2(2,0)]
        sage: [A5(w).branch(C3,rule="symmetric") for w in A5.lattice().fundamental_weights()]
        [C3(1,0,0),
         C3(0,0,0) + C3(1,1,0),
         C3(1,0,0) + C3(1,1,1),
         C3(0,0,0) + C3(1,1,0),
         C3(1,0,0)]
        sage: [A5(w).branch(D3,rule="symmetric") for w in A5.lattice().fundamental_weights()]
        [D3(1,0,0), D3(1,1,0), D3(1,1,-1) + D3(1,1,1), D3(1,1,0), D3(1,0,0)]
        sage: [D4(x).branch(B3,rule="symmetric") for x in D4.lattice().fundamental_weights()]
        [B3(0,0,0) + B3(1,0,0),
         B3(1,0,0) + B3(1,1,0),
         B3(1/2,1/2,1/2),
         B3(1/2,1/2,1/2)]

    EXAMPLES: (Extended type)

    ::

        sage: [B3(x).branch(D3,rule="extended") for x in B3.lattice().fundamental_weights()]
        [D3(0,0,0) + D3(1,0,0),
         D3(1,0,0) + D3(1,1,0),
         D3(1/2,1/2,-1/2) + D3(1/2,1/2,1/2)]

    EXAMPLES: (Isomorphic type)

    ::

        sage: [B2(x).branch(C2, rule="isomorphic") for x in B2.lattice().fundamental_weights()]
        [C2(1,1), C2(1,0)]
        sage: [C2(x).branch(B2, rule="isomorphic") for x in C2.lattice().fundamental_weights()]
        [B2(1/2,1/2), B2(1,0)]
        sage: [A3(x).branch(D3,rule="isomorphic") for x in A3.lattice().fundamental_weights()]
        [D3(1/2,1/2,1/2), D3(1,0,0), D3(1/2,1/2,-1/2)]
        sage: [D3(x).branch(A3,rule="isomorphic") for x in D3.lattice().fundamental_weights()]
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

    You may also write your own rules. We may arrange a Cartan
    subalgebra U of H to be contained in a Cartan subalgebra T of G.
    This embedding must be chosen in a particular way - the restriction
    of the positive Weyl chamber in G.lattice() must be contained in
    the positive Weyl chamber in H.lattice(). The RULE is the adjoint
    of this embedding U - T, which is a mapping from the dual space of
    T, which is the weight space of G, to the weight space of H.

    EXAMPLES::

        sage: A3 = WeylCharacterRing(['A',3])
        sage: C2 = WeylCharacterRing(['C',2])
        sage: rule = lambda x : [x[0]-x[3],x[1]-x[2]]
        sage: branch_weyl_character(A3([1,1,0,0]),A3,C2,rule)
        C2(0,0) + C2(1,1)
        sage: A3(1,1,0,0).branch(C2, rule) == C2(0,0) + C2(1,1)
        True
    """
    r = R.cartan_type[1]
    s = S.cartan_type[1]
    # Each rule takes a tuple or list and returns a list
    if rule == "levi":
        if not s == r-1:
            raise ValueError, "Rank is wrong"
        if R.cartan_type[0] == 'A':
            if S.cartan_type[0] == 'A':
                rule = lambda x : list(x)[:r]
            else:
                raise ValueError, "Rule not found"
        elif R.cartan_type[0] in ['B', 'C', 'D']:
            if S.cartan_type[0] == 'A':
                rule = lambda x : x
            elif S.cartan_type[0] == R.cartan_type[0]:
                rule = lambda x : list(x)[1:]
            else:
                raise ValueError, "Rule not found"
        elif S.cartan_type[0] == 'E' and R.cartan_type[0] in ['A','D','E']:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif S.cartan_type[0] == 'F' and R.cartan_type[0] in ['B','C']:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif S.cartan_type[0] == 'G' and R.cartan_type[0] == 'A':
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "Rule not found"
    elif rule == "automorphic":
        if not R.cartan_type == S.cartan_type:
            raise ValueError, "Cartan types must agree for automorphic branching rule"
        elif R.cartan_type[0] == 'E' and R.cartan_type[1] == 6:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif R.cartan_type[0] == 'A':
            def rule(x) : y = [-i for i in x]; y.reverse(); return y
        elif R.cartan_type[0] == 'D':
            def rule(x) : x[len(x)-1] = -x[len(x)-1]; return x
        elif R.cartan_type[0] == 'E' and R.cartan_type[1] == 6:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "No automorphism found"
    elif rule == "triality":
        if not R.cartan_type == S.cartan_type:
            raise ValueError, "Triality is an automorphic type (for D4 only)"
        elif not R.cartan_type[0] == 'D' and r == 4:
            raise ValueError, "Triality is for D4 only"
        else:
            def rule(x):
                [x1,x2,x3,x4] = x
                return [(x1+x2+x3+x4)/2,(x1+x2-x3-x4)/2,(x1-x2+x3-x4)/2,(-x1+x2+x3-x4)/2]
    elif rule == "symmetric":
        if R.cartan_type[0] == 'A':
            if (S.cartan_type[0] == 'C' or S.cartan_type[0] == 'D' and r == 2*s-1) or (S.cartan_type[0] == 'B' and r == 2*s):
                def rule(x):
                    return [x[i]-x[r-i] for i in range(s)]
            else:
                print S.cartan_type, r, s
                raise ValueError, "Rule not found"

        elif R.cartan_type[0] == 'D' and S.cartan_type[0] == 'B' and s == r-1:
            rule = lambda x : x[:s]

        elif R.cartan_type[0] == 'E' and S.cartan_type[0] == '6':
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "Rule not found"
    elif rule == "extended":
        if R.cartan_type[0] == 'B' and S.cartan_type[0] == 'D' and s == r:
            rule = lambda x : x
        elif R.cartan_type[0] == 'G' and S.cartan_type[0] == 'A' and s == r:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif R.cartan_type[0] == 'F' and S.cartan_type[0] == 'B' and s == r:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        elif R.cartan_type[0] == 'E' and S.cartan_type[0] == 'E' and s == r and r >= 7:
            raise NotImplementedError, "Exceptional branching rules are yet to be implemented"
        else:
            raise ValueError, "Rule not found"
    elif rule == "isomorphic":
        if R.cartan_type[0] == 'B' and S.cartan_type[0] == 'C' and s == 2 and r == 2:
            def rule(x):
                [x1, x2] = x
                return [x1+x2, x1-x2]
        elif R.cartan_type[0] == 'C' and S.cartan_type[0] == 'B' and s == 2 and r == 2:
            def rule(x):
                [x1, x2] = x
                return [(x1+x2)/2, (x1-x2)/2]
        elif R.cartan_type[0] == 'A' and S.cartan_type[0] == 'D' and s == 3 and r == 3:
            def rule(x):
                [x1, x2, x3, x4] = x
                return [(x1+x2-x3-x4)/2, (x1-x2+x3-x4)/2, (x1-x2-x3+x4)/2]
        elif R.cartan_type[0] == 'D' and S.cartan_type[0] == 'A' and s == 3 and r == 3:
            def rule(x):
                [t1, t2, t3] = x
                return [(t1+t2+t3)/2, (t1-t2-t3)/2, (-t1+t2-t3)/2, (-t1-t2+t3)/2]
        else:
            raise ValueError, "Rule not found"

    mdict = {}
    for k in chi._mdict:
        h = S._lattice(rule(list(k.to_vector())))
        if h in mdict:
            mdict[h] += chi._mdict[k]
        else:
            mdict[h] = chi._mdict[k]
    hdict = S.char_from_weights(mdict)
    return WeylCharacter(S, hdict, mdict)



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
        self._lattice = A._lattice

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
        for i in range(1,self._lattice.n):
            hstring = hstring+","+str(k[i])
        return self._prefix+"("+hstring+")"

    def __repr__(self):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: b3 = WeightRing(B3)
            sage: fw = b3.lattice().fundamental_weights()
            sage: b3(fw[3])
            b3(1/2,1/2,1/2)
            sage: b3(B3(fw[3]))
            b3(-1/2,-1/2,-1/2) + b3(-1/2,-1/2,1/2) + b3(-1/2,1/2,-1/2) + b3(-1/2,1/2,1/2) + b3(1/2,-1/2,-1/2) + b3(1/2,-1/2,1/2) + b3(1/2,1/2,-1/2) + b3(1/2,1/2,1/2)
        """

        if self._mdict == {}:
            return "0"
        v = self._mdict.keys()
        # Just a workaround to keep the same sorting as before when
        # the dictionary was indexed by tuples
        v.sort(key = lambda v: tuple(v.to_vector()))
        return repr_lincomb([self._wt_repr(k) for k in v], [self._mdict[k] for k in v])

    def __cmp__(left, right):
        """
        EXAMPLES::

            sage: B3 = WeylCharacterRing(['B',3])
            sage: b3 = WeightRing(B3)
            sage: fw = [b3(w) for w in b3.lattice().fundamental_weights()]
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

    def mlist(self):
        """
        Returns a list of weights in self with their multiplicities.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: pr = sum(g2(a) for a in g2.lattice().positive_roots())
            sage: pr.mlist()
            [[(1, -2,  1), 1],
             [(1, -1,  0), 1],
             [(1,  0, -1), 1],
             [(2, -1, -1), 1],
             [(0,  1, -1), 1],
             [(1,  1, -2), 1]]
        """
        return [[k,m] for k,m in self._mdict.iteritems()]

    def weyl_group_action(self, w):
        """
        Returns the actionof the Weyl group element w on self.

        EXAMPLES::

            sage: G2 = WeylCharacterRing(['G',2])
            sage: g2 = WeightRing(G2)
            sage: L = g2.lattice()
            sage: [fw1, fw2] = L.fundamental_weights()
            sage: sum(g2(fw2).weyl_group_action(w) for w in L.weyl_group())
            2*g2(-2,1,1) + 2*g2(-1,-1,2) + 2*g2(-1,2,-1) + 2*g2(1,-2,1) + 2*g2(1,1,-2) + 2*g2(2,-1,-1)
        """
        return WeightRingElement(self._parent, dict([[w.action(x),self._mdict[x]] for x in self._mdict.keys()]))

    def character(self):
        """
        Assuming that self is invariant under the Weyl group, this will
        express it as a linear combination of characters. If self is not
        Weyl group invariant, this method will not terminate.

        EXAMPLES::

            sage: A2 = WeylCharacterRing(['A',2])
            sage: a2 = WeightRing(A2)
            sage: W = a2.lattice().weyl_group()
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
            sage: wd = prod(a2(x/2)-a2(-x/2) for x in a2.lattice().positive_roots()); wd
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
            sage: sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.lattice().weyl_group())
            -a2(-1,3,6) + a2(-1,6,3) + a2(3,-1,6) - a2(3,6,-1) - a2(6,-1,3) + a2(6,3,-1)
            sage: a2(chi)*wd == sum((-1)^w.length()*a2([6,3,-1]).weyl_group_action(w) for w in a2.lattice().weyl_group())
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
            sage: L = S.lattice()
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
        self.cartan_type = self._parent.cartan_type
        if prefix == None:
            if self._parent._prefix.isupper():
                prefix = self._parent._prefix.lower()
            elif self._parent._prefix.islower():
                prefix = self._parent._prefix.upper()
            else:
                prefix = (self.cartan_type[0].lower()+str(self.cartan_type[1]))
        self._base_ring = self._parent._base_ring
        self._lattice = self._parent._lattice
        self._origin = self._parent._origin
        self._prefix = prefix
        self._ip = self._parent._ip

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
        x = self._lattice(x)
        mdict = {x: 1}
        return WeightRingElement(self, mdict)

    def __repr__(self):
        """
        EXAMPLES::

            sage: P.<q>=QQ[]
            sage: G2 = WeylCharacterRing(['G',2], base_ring = P)
            sage: WeightRing(G2)
            The Weight ring attached to The Weyl Character Ring of Type [G,2] with Univariate Polynomial Ring in q over Rational Field coefficients
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

    def lattice(self):
        """
        Returns the weight lattice realization associated to self.

        EXAMPLES::

            sage: E8 = WeylCharacterRing(['E',8])
            sage: e8 = WeightRing(E8)
            sage: e8.lattice()
            Ambient space of the Root system of type ['E', 8]
        """
        return self._lattice

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
        for i in range(1,self._lattice.n):
            hstring=hstring+","+str(wt[i])
        return self._prefix+"("+hstring+")"


