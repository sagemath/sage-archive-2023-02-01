r"""
Discrete dynamical systems
==========================

A *discrete dynamical system* (short: *DDS*) is a pair `(S, \phi)` of
a set `S` and a map `\phi : S \to S`.
(This is one of several things known as a "discrete dynamical system"
in mathematics.)
Thus, a DDS is the same as an endomorphism of a set.
The DDS is said to be *finite* if `S` is finite;
this code is mostly concerned with finite DDSs.

Given a DDS `(S, \phi)`, we can study its orbits, its invariants
(i.e., maps `f : S \to X` satisfying `f \circ \phi = f`),
its homomesies (i.e., maps `h : S \to A` to a `\QQ`-vector space `A`
such that the average of the values of `h` on each `\phi`-orbit is
the same), and various other phenomena.
(Some of these require `S` to be finite or at least to have finite
orbits.)

- :class:`DiscreteDynamicalSystem`: general discrete dynamical
  system, as above.

- :class:`FiniteDynamicalSystem`: finite discrete dynamical
  system.

.. TODO ::

- DiscreteDynamicalSystem class, and orbit method.

- Implement basic functionality for homomesy and invariance testing
  when the endo is an auto:
  is_homomesic, is_invariant, orbits, orbit, orbit_lengths,
  orbit_lengths_lcm.

- Examples:
  - finite set with permutation;
  - rowmotion on order ideals;
  - promotion on std YTs.
  - rotation of bitstrings.

- Implement caching for orbits.

- non-auto functionality: is_recurrent, recurrent_entries,
  first_recurrent_image, lollipop.

- Examples for non-auto functionality:
  - finite set with non-permutation;
  - infection on a chessboard;
  - Conway's GoL.

- Wrap (some of) the cyclic_sieving_phenomenon.py methods.

"""
#*****************************************************************************
#       Copyright (C) 2018 Darij Grinberg <darijgrinberg@gmail.com>,
#                     2018 Tom Roby <tomrobyuconn@gmail.com>
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
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.abstract_method import abstract_method
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_generic

class DiscreteDynamicalSystem(object):
    r"""
    A discrete dynamical system.

    A *discrete dynamical system* (henceforth *DDS*) is a
    pair `(S, \phi)` of a set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.

    At the moment, support exists only for the case when `\phi`
    is invertible.
    """
    def __init__(self, X, phi, invertible=True, cache_orbits=False):
        r"""
        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D.ground_set()
            Non negative integer semiring
            sage: D.evolution()(5)
            7
        """
        # what if X is a list? tuple it?
        self._X = X
        self._phi = phi
        self._invertible = invertible
        self._cache_orbits = cache_orbits

    def ground_set(self):
        # wrap it to avoid clobbering?
        return self._X

    def evolution(self):
        return self._phi

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D[:3]
            [0, 1, 2]
        """
        return iter(self._X)

    def __getitem__(self, i):
        return self._X[i]

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D # indirect doctest
            Discrete dynamical system with ground set Non negative integer semiring
            and evolution <function <lambda> at ...
        """
        return "Discrete dynamical system with ground set " \
               + repr(self._X) + " and evolution " + \
               repr(self._phi)

    def orbit(self, x):
        r"""
        Return the orbit of the element ``x`` of the ground set
        of ``self``.

        This orbit is a list beginning with ``x`` and ending
        with the last element until ``x`` reappears.
        If ``x`` never reappears, then this will not terminate!

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(range(8), lambda x : (x + 2) % 8)
            sage: D.ground_set()
            [0, 1, 2, 3, 4, 5, 6, 7]
            sage: D.orbit(2)
            [2, 4, 6, 0]
            sage: D.orbit(5)
            [5, 7, 1, 3]
        """
        orb = [x]
        phi = self._phi
        curr = phi(x)
        while curr != x:
            orb.append(curr)
            curr = phi(curr)
        return orb

class FiniteDynamicalSystem(DiscreteDynamicalSystem):
    r"""
    A finite discrete dynamical system.

    A *finite discrete dynamical system* (henceforth *FDDS*) is a
    pair `(S, \phi)` of a finite set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.

    At the moment, support exists only for the case when `\phi`
    is invertible.

    EXAMPLES::

        sage: D = FiniteDynamicalSystem(range(5), lambda x : (x + 2) % 5)
        sage: D.ground_set()
        [0, 1, 2, 3, 4]
        sage: D.evolution()(4)
        1
        sage: D.orbits()
        [[4, 1, 3, 0, 2]]

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (s[0],)
        sage: D = FiniteDynamicalSystem(Y, rot)
        sage: D.evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
        sage: sorted(D.orbit_lengths())
        [2, 4, 8, 8, 8, 8, 8, 8, 8, 8]
    """
    def __repr__(self):
        r"""
        """
        return "Finite discrete dynamical system with ground set " \
               + repr(self._X) + " and evolution " + \
               repr(self._phi)

    def orbits(self):
        r"""
        Return a list of all orbits of ``self``.

        If the evolution of ``self`` is not invertible, then
        this will not terminate!
        """
        phi = self._phi
        l = list(self)
        orbs = []
        while l:
            start = l.pop()
            orb = [start]
            curr = phi(start)
            while curr != start:
                l.remove(curr)
                orb.append(curr)
                curr = phi(curr)
            orbs.append(orb)
        return orbs

    def orbit_lengths(self):
        return [len(orb) for orb in self.orbits()]

    def is_homomesic(self, h, average=None):
        r"""
        Check if ``h`` is homomesic with respect to ``self``,
        with average ``average`` along each orbit (if not provided,
        then just checks that the averages are equal).
        """
        orbavgs = []
        for orb in self.orbits():
            l = len(orb)
            avg = QQ(1) / QQ(l) * sum(h(i) for i in orb)
            if avg not in orbavgs:
                orbavgs.append(avg)
            if len(orbavgs) > 1:
                return False
        if not orbavgs:
            return True
        if average is None:
            return True
        return orbavgs[0] == average

    def is_invariant(self, f):
        r"""
        Check if ``f`` is an invariant of ``self``.
        """
        phi = self._phi
        for i in self._X:
            if f(phi(i)) != f(i):
                return False
        return True

def FDDSFromPermutation(pi):
    r"""
    Return the finite discrete dynamical system induced by
    the permutation ``pi``.
    """
    from sage.combinat.permutation import Permutation
    pi = Permutation(pi)
    n = len(pi)
    X = range(1, n+1)
    return FiniteDynamicalSystem(X, pi)

def FDDSFromSYTPromotion(lam):
    r"""
    Return the finite discrete dynamical system consisting
    of all standard tableaux of shape ``lam`` and evolving
    according to promotion.

    EXAMPLES::

        sage: from sage.combinat.finite_dynamical_system import FDDSFromSYTPromotion
        sage: F = FDDSFromSYTPromotion([4, 4, 4])
        sage: sorted(F.orbit_lengths())
        [3, 3, 4, 4, 4, 6, 6, 6, 6, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
         12, 12, 12, 12, 12]
        sage: G = FDDSFromSYTPromotion([4, 3, 1])
        sage: sorted(G.orbit_lengths())
        [16, 22, 32]
    """
    from sage.combinat.partition import Partition
    lam = Partition(lam)
    from sage.combinat.tableau import StandardTableaux
    X = StandardTableaux(lam)
    return FiniteDynamicalSystem(X, lambda T : T.promotion())

def FDDSFromRowmotion(P):
    r"""
    Return the finite discrete dynamical system consisting
    of all order ideals of the poset ``P``.

    EXAMPLES::

        sage: P = RootSystem(["A", 6]).root_poset()
        sage: from sage.combinat.finite_dynamical_system import FDDSFromRowmotion
        sage: F = FDDSFromRowmotion(P)
        sage: sorted(F.orbit_lengths())
        [2, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
         14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
        sage: F.is_homomesic(lambda I: len(I))
        False
        sage: F.is_homomesic(lambda I: sum((-1)**(P.rank(i)) for i in I))
        True

    (These are perhaps long-time.)
    """
    from sage.sets.set import Set
    X = [Set(P.order_ideal(A)) for A in P.antichains()]
    # Using P.order_ideals_lattice() instead causes intransparency issues:
    # sage can't always do P.rowmotion(I) when I is in P.order_ideals_lattice().
    # Bug in P.order_ideals_lattice() when P is facade?
    phi = lambda I : P.rowmotion(I)
    return FiniteDynamicalSystem(X, phi)


