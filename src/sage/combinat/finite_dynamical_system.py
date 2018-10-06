r"""
Discrete dynamical systems
==========================

A *discrete dynamical system* (short: *DDS*) is a pair `(S, \phi)`
of a set `S` and a map `\phi : S \to S`.
(This is one of several things known as a "discrete dynamical
system" in mathematics.)
Thus, a DDS is the same as an endomorphism of a set.
The DDS is said to be *finite* if `S` is finite.
The DDS is said to be *invertible* if the map `\phi` is
invertible.
The set `S` is called the *ground set* of the DDS;
the map `\phi` is called the *evolution* of the DDS;
the inverse map `\phi^{-1}` is called the *inverse evolution*
of the DDS.

Given a DDS `(S, \phi)`, we can study its orbits, its invariants
(i.e., maps `f : S \to X` satisfying `f \circ \phi = f`),
its homomesies (i.e., maps `h : S \to A` to a `\QQ`-vector space `A`
such that the average of the values of `h` on each `\phi`-orbit is
the same), and various other phenomena.
(Some of these require `S` to be finite or at least to have finite
orbits.)

This file implements the following classes for discrete
dynamical systems:

- :class:`DiscreteDynamicalSystem`: general discrete dynamical
  system, as above.
  Inherit from this class if the ground set of your DDS is
  infinite or large enough that you want to avoid it getting
  stored as a list.

- :class:`FiniteDynamicalSystem`: finite discrete dynamical
  system.

- :class:`InvertibleDiscreteDynamicalSystem`: invertible
  discrete dynamical system.
  This implements an ``inverse_evolution`` method for `\phi^{-1}`
  (the default implementation simply applies `\phi` over and
  over until the original value is revisited; the last value
  before that is then taken to be the result).

- :class:`InvertibleFiniteDynamicalSystem`: invertible
  finite discrete dynamical system.

.. TODO::

    - Implement basic functionality for homomesy and invariance testing
    when the endo is an auto:
    is_homomesic, is_invariant, orbits, orbit, orbit_lengths,
    orbit_lengths_lcm.

    - General classcall that delegates to subclasses?

    - Subclasses for DDSes whose ground set is an enumerated set?

    - Implement caching for orbits?

    - Possibly implement lazy ground set -- or should this just be
    what you get when you use DiscreteDynamicalSystem rather than
    FiniteDynamicalSystem?

    - non-auto functionality: is_recurrent, recurrent_entries,
    first_recurrent_image, rho, cycles, homomesy?

    - Examples for non-auto functionality:
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
#from sage.categories.enumerated_sets import EnumeratedSets
#from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
#from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
#from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
from sage.misc.abstract_method import abstract_method
from sage.rings.rational_field import QQ
from sage.sets.recursively_enumerated_set import RecursivelyEnumeratedSet_generic
from sage.structure.sage_object import SageObject

class DiscreteDynamicalSystem(SageObject):
    r"""
    A discrete dynamical system.

    A *discrete dynamical system* (henceforth *DDS*) is a
    pair `(S, \phi)` of a set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the DDS, while
    the map `\phi` is called the *evolution* of the DDS.

    See the module-level doc for details.
    """
    def __init__(self, X, phi, cache_orbits=False, create_tuple=False):
        r"""
        INPUT:

        - ``X`` (set, list, tuple, or another iterable, or
          ``None``) -- the ground set for the DDS; this can be
          ``None`` in case of a
          :class:`DiscreteDynamicalSystem` or a
          :class:`InvertibleDiscreteDynamicalSystem`.
          Make sure to set the ``create_tuple`` argument to
          ``True`` if you provide an iterator or a list for
          ``X``, as otherwise the input would be exposed.

        - ``phi`` (function, or callable that acts like a
          function) -- the evolution of the DDS.

        - ``cache_orbits`` (boolean) -- (default: ``False``)
          whether or not the orbits should be cached once they
          are computed.

        - ``create_tuple`` (boolean) -- (default: ``False``)
          whether or not the input ``X`` should be translated
          into a tuple (set this to ``True`` to prevent
          mutation if ``X`` is a list, and to prevent
          exhaustion if ``X`` is an iterator).

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D.ground_set()
            Non negative integer semiring
            sage: D.evolution()(5)
            7

        The necessity of ``create_tuple=True``::

            sage: X = [0, 1, 2, 3, 4]
            sage: D_wrong = DiscreteDynamicalSystem(X, lambda x : (x**3) % 5)
            sage: D_right = DiscreteDynamicalSystem(X, lambda x : (x**3) % 5, create_tuple=True)
            sage: X[4] = 666 # evil
            sage: D_wrong.ground_set()
            [0, 1, 2, 3, 666]
            sage: D_right.ground_set()
            (0, 1, 2, 3, 4)
        """
        if create_tuple:
            X = tuple(X)
        self._X = X
        self._phi = phi
        self._cache_orbits = cache_orbits

    def ground_set(self):
        r"""
        Return the ground set of ``self``.

        This will return ``None`` if no ground set was
        provided in the construction of ``self``.

        .. WARNING::

            Unless ``self`` has been constructed with the
            ``create_tuple`` parameter set to ``True``,
            this method will return whatever ground set was
            provided to the constructor.
            In particular, if a list was provided, then this
            precise list will be returned; mutating this list
            will then corrupt ``self``.

        EXAMPLES::
        
            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x : (3 if x == 4 else 1), create_tuple=True)
            sage: D.ground_set()
            (1, 3, 4)
        """
        return self._X

    def evolution(self):
        r"""
        Return the evolution of ``self``.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem([1, 3, 4], lambda x : (3 if x == 4 else 1), create_tuple=True)
            sage: ev = D.evolution()
            sage: ev(1)
            1
            sage: ev(4)
            3
        """
        return self._phi

    def __iter__(self):
        r"""
        Iterate over the ground set of ``self``.

        This assumes that an iterable ground set of ``self``
        has been provided.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D[:3]
            [0, 1, 2]
        """
        return iter(self._X)

    def __getitem__(self, i):
        r"""
        Accessor for the ``i``-th entry of the ground set of
        ``self``, assuming that the latter ground set was
        provided as a list or tuple when ``self`` was constructed.

        EXAMPLES::

            sage: X = (5, 6, 7, 8)
            sage: D = DiscreteDynamicalSystem(X, lambda x : (x**3) % 4 + 5)
            sage: X[3]
            8
            sage: X[2]
            7
        """
        return self._X[i]

    def _repr_(self):
        r"""
        String representation of ``self``.
        
        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = DiscreteDynamicalSystem(NN, lambda x : x + 2)
            sage: D # indirect doctest
            A discrete dynamical system with ground set
             Non negative integer semiring
            sage: D = DiscreteDynamicalSystem(None, lambda x : x + 2)
            sage: D # indirect doctest
            A discrete dynamical system with unspecified ground set
        """
        if self._X is None:
            return "A discrete dynamical system with unspecified ground set"
        return "A discrete dynamical system with ground set " \
               + repr(self._X)

    def orbit(self, x, preperiod=False):
        r"""
        Return the orbit of the element ``x`` of the ground set
        of ``self`` under the evolution `\phi` of ``self``.

        This orbit is a list beginning with ``x`` and ending
        with the last element that is not a repetition of a
        previous element.
        If the orbit is infinite, then this method does not
        terminate!

        If the optional argument ``preperiod`` is set to
        ``True``, then this method returns a pair ``(o, k)``,
        where ``o`` is the orbit of ``self``, while `k` is the
        smallest nonnegative integer such that
        `\phi^k(x) \in \left\{ \phi^i(x) \mid i > k \right\}`.

        EXAMPLES::

            sage: D = DiscreteDynamicalSystem(tuple(range(11)), lambda x : (x ** 2) % 11)
            sage: D.orbit(6)
            [6, 3, 9, 4, 5]
            sage: D.orbit(6, preperiod=True)
            ([6, 3, 9, 4, 5], 1)
            sage: D.orbit(3)
            [3, 9, 4, 5]
            sage: D.orbit(3, preperiod=True)
            ([3, 9, 4, 5], 0)
            sage: D.orbit(9)
            [9, 4, 5, 3]
            sage: D.orbit(0)
            [0]
        """
        orb = [x]
        phi = self._phi
        curr = phi(x)
        while curr not in orb:
            orb.append(curr)
            curr = phi(curr)
        if not preperiod:
            return orb
        return (orb, orb.index(curr))

class InvertibleDiscreteDynamicalSystem(DiscreteDynamicalSystem):
    r"""
    An invertible discrete dynamical system.

    A *discrete dynamical system* (henceforth *DDS*) is a
    pair `(S, \phi)` of a set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the DDS, while
    the map `\phi` is called the *evolution* of the DDS.
    An *invertible DDS* is a DDS `(S, \phi)` whose evolution
    `\phi` is invertible.
    In that case, `\phi^{-1}` is called the *inverse evolution*
    of the DDS.

    See the module-level doc for details.
    """
    def __init__(self, X, phi, inverse=None, cache_orbits=False, create_tuple=False):
        r"""
        INPUT:

        - ``X`` (set, list, tuple, or another iterable, or
          ``None``) -- the ground set for the DDS; this can be
          ``None`` in case of a
          :class:`DiscreteDynamicalSystem` or a
          :class:`InvertibleDiscreteDynamicalSystem`.
          Make sure to set the ``create_tuple`` argument to
          ``True`` if you provide an iterator or a list for
          ``X``, as otherwise the input would be exposed.

        - ``phi`` (function, or callable that acts like a
          function) -- the evolution of the DDS.

        - ``inverse`` (function, or callable that acts like a
          function) -- the inverse evolution
          of the DDS. (A default implementation is implemented
          when this argument is not provided; but it assumes
          the orbits to be finite.)

        - ``cache_orbits`` (boolean) -- (default: ``False``)
          whether or not the orbits should be cached once they
          are computed.

        - ``create_tuple`` (boolean) -- (default: ``False``)
          whether or not the input ``X`` should be translated
          into a tuple (set this to ``True`` to prevent
          mutation if ``X`` is a list, and to prevent
          exhaustion if ``X`` is an iterator).

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = InvertibleDiscreteDynamicalSystem(NN, lambda x : (x + 2 if x % 4 < 2 else x - 2))
            sage: D.ground_set()
            Non negative integer semiring
            sage: D.evolution()(5)
            7
            sage: D.evolution()(6)
            4
            sage: D.evolution()(4)
            6
            sage: D.inverse_evolution()(4)
            6

        The necessity of ``create_tuple=True``::

            sage: X = [0, 1, 2, 3, 4]
            sage: D_wrong = InvertibleDiscreteDynamicalSystem(X, lambda x : (x**3) % 5)
            sage: D_right = InvertibleDiscreteDynamicalSystem(X, lambda x : (x**3) % 5, create_tuple=True)
            sage: X[4] = 666 # evil
            sage: D_wrong.ground_set()
            [0, 1, 2, 3, 666]
            sage: D_right.ground_set()
            (0, 1, 2, 3, 4)
        """
        if create_tuple:
            X = tuple(X)
        self._X = X
        self._phi = phi
        if inverse is None:
            self._inverse = self.inverse_evolution_default
        else:
            self._inverse = inverse
        self._cache_orbits = cache_orbits

    def _repr_(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            sage: D = InvertibleDiscreteDynamicalSystem(NN, lambda x : (x + 2 if x % 4 < 2 else x - 2))
            sage: D # indirect doctest
            An invertible discrete dynamical system with ground set
             Non negative integer semiring
            sage: D = InvertibleDiscreteDynamicalSystem(None, lambda x : x + 2)
            sage: D # indirect doctest
            An invertible discrete dynamical system with unspecified ground set
        """
        if self._X is None:
            return "An invertible discrete dynamical system with unspecified ground set"
        return "An invertible discrete dynamical system with ground set " \
               + repr(self._X)

    def inverse_evolution(self):
        r"""
        Return the inverse evolution of ``self`` (as a map
        from the ground set of ``self`` to itself).

        EXAMPLES::

            sage: D = InvertibleDiscreteDynamicalSystem(tuple(range(8)), lambda x : (x + 2) % 8)
            sage: D.inverse_evolution()(1)
            7
            sage: D.inverse_evolution()(3)
            1
        """
        return self._inverse

    def verify_inverse_evolution(self, x=None):
        r"""
        Verify that the composition of evolution and
        inverse evolution on ``self`` is the identity
        (both ways).

        The optional argument ``x``, if provided, restricts
        the testing to the element ``x`` only.
        Otherwise, all elements of the ground set are
        tested (if they can be enumerated).

        This is mostly used to check the correctness of
        self-implemented inverse evolution methods.

        EXAMPLES::

            sage: D = InvertibleDiscreteDynamicalSystem(tuple(range(8)), lambda x : (x + 2) % 8)
            sage: D.verify_inverse_evolution()
            True
            sage: D.verify_inverse_evolution(3)
            True
            sage: fake_inverse = lambda x : x
            sage: D = InvertibleDiscreteDynamicalSystem(tuple(range(8)), lambda x : (x + 2) % 8, inverse=fake_inverse)
            sage: D.verify_inverse_evolution()
            False
            sage: D.verify_inverse_evolution(3)
            False
        """
        if x is None:
            els = self.ground_set()
        else:
            els = [x]
        for y in els:
            if self.evolution()(self.inverse_evolution()(y)) != y:
                return False
            if self.inverse_evolution()(self.evolution()(y)) != y:
                return False
        return True

    def orbit(self, x, preperiod=False):
        r"""
        Return the orbit of the element ``x`` of the ground set
        of ``self``.

        This orbit is a list beginning with ``x`` and ending
        with the last element until ``x`` reappears.
        If ``x`` never reappears, then this will not terminate!

        If the optional argument ``preperiod`` is set to
        ``True``, then this method returns a pair ``(o, k)``,
        where ``o`` is the orbit of ``self``, while `k` is the
        smallest nonnegative integer such that
        `\phi^k(x) \in \left\{ \phi^i(x) \mid i > k \right\}`.
        Note that `k` is necessarily `0`, since the DDS
        ``self`` is invertible!

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import DiscreteDynamicalSystem
            sage: D = InvertibleDiscreteDynamicalSystem(tuple(range(8)), lambda x : (x + 2) % 8)
            sage: D.ground_set()
            (0, 1, 2, 3, 4, 5, 6, 7)
            sage: D.orbit(2)
            [2, 4, 6, 0]
            sage: D.orbit(5)
            [5, 7, 1, 3]
            sage: D.orbit(5, preperiod=True)
            ([5, 7, 1, 3], 0)
        """
        orb = [x]
        phi = self._phi
        curr = phi(x)
        while curr != x:
            orb.append(curr)
            curr = phi(curr)
        if not preperiod:
            return orb
        return (orb, 0)

    def inverse_evolution_default(self, x):
        r"""
        Return the inverse evolution of ``self``, applied
        to the element ``x`` of the ground set of ``self``.

        This is the default implementation, assuming that
        the orbit of ``x`` is finite.

        EXAMPLES::

            sage: D = InvertibleDiscreteDynamicalSystem(tuple(range(8)), lambda x : (x + 2) % 8)
            sage: D.inverse_evolution_default(1)
            7
            sage: D.inverse_evolution_default(3)
            1
        """
        return self.orbit(x)[-1]

class FiniteDynamicalSystem(DiscreteDynamicalSystem):
    r"""
    A finite discrete dynamical system.

    A *finite discrete dynamical system* (henceforth *FDDS*) is a
    pair `(S, \phi)` of a finite set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.

    The ground set `S` should always be provided as an
    iterable when defining a :class:`FiniteDynamicalSystem`.

    EXAMPLES::

        sage: D = FiniteDynamicalSystem(tuple(range(11)), lambda x : (x**2) % 11)
        sage: D.ground_set()
        (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        sage: D.evolution()(4)
        5
        sage: D.orbit(4)
        [4, 5, 3, 9]
        sage: D.orbit(1)
        [1]
        sage: D.orbit(2)
        [2, 4, 5, 3, 9]

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (0,)
        sage: D = FiniteDynamicalSystem(Y, rot)
        sage: D.evolution()((1, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
    """
    def _repr_(self):
        r"""
        String representation of ``self``.
        
        EXAMPLES::
        
            sage: D = FiniteDynamicalSystem(tuple(range(11)), lambda x : (x**2) % 11)
            sage: D
            A finite discrete dynamical system with ground set
            (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
        """
        return "A finite discrete dynamical system with ground set " \
               + repr(self._X)

    def is_invariant(self, f):
        r"""
        Check if ``f`` is an invariant of ``self``.

        EXAMPLES::

            sage: W = Words(2, 5)
            sage: F = FiniteDynamicalSystem(W, lambda x : x[1:] + Word([x[0]]))
            sage: F.is_invariant(lambda w: sum(w))
            True
            sage: F.is_invariant(lambda w: 1)
            True
            sage: F.is_invariant(lambda w: w[0] - w[1])
            False
            sage: F.is_invariant(lambda w: sum(i**2 for i in w))
            True

        An invariant of a permutation::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.permutation([3, 4, 5, 6, 1, 2])
            sage: F.is_invariant(lambda i: i % 2)
            True
            sage: F.is_invariant(lambda i: i % 3)
            False
        """
        phi = self._phi
        for i in self._X:
            if f(phi(i)) != f(i):
                return False
        return True

class InvertibleFiniteDynamicalSystem(InvertibleDiscreteDynamicalSystem, FiniteDynamicalSystem):
    r"""
    An invertible finite discrete dynamical system.

    A *finite discrete dynamical system* (henceforth *FDDS*) is a
    pair `(S, \phi)` of a finite set `S` and a map `\phi : S \to S`.
    This set `S` is called the *ground set* of the FDDS, while
    the map `\phi` is called the *evolution* of the FDDS.
    An FDDS `(S, \phi)` is called *invertible* if the map `\phi`
    is invertible; in this case, `\phi^{-1}` is called the
    *inverse evolution* of the FDDS.

    The ground set `S` should always be provided as an
    iterable when defining a :class:`FiniteDynamicalSystem`.

    EXAMPLES::

        sage: D = InvertibleFiniteDynamicalSystem(tuple(range(5)), lambda x : (x + 2) % 5)
        sage: D.ground_set()
        (0, 1, 2, 3, 4)
        sage: D.evolution()(4)
        1
        sage: D.orbits()
        [[4, 1, 3, 0, 2]]
        sage: D.inverse_evolution()(2)
        0
        sage: D.inverse_evolution()(1)
        4

        sage: X = cartesian_product([[0, 1]]*8)
        sage: Y = [s for s in X if sum(s) == 4]
        sage: rot = lambda s : s[1:] + (s[0],)
        sage: D = InvertibleFiniteDynamicalSystem(Y, rot)
        sage: D.evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 1, 0, 1, 0, 0, 1, 0)
        sage: D.inverse_evolution()((0, 1, 1, 0, 1, 0, 0, 1))
        (1, 0, 1, 1, 0, 1, 0, 0)
        sage: sorted(D.orbit_lengths())
        [2, 4, 8, 8, 8, 8, 8, 8, 8, 8]
    """
    def _repr_(self):
        r"""
        String representation of ``self``.
        
        EXAMPLES::
        
            sage: D = InvertibleFiniteDynamicalSystem(tuple(range(5)), lambda x : (x + 2) % 5)
            sage: D
            An invertible finite discrete dynamical system with ground set
            (0, 1, 2, 3, 4)
        """
        return "An invertible finite discrete dynamical system with ground set " \
               + repr(self._X)

    def orbits(self):
        r"""
        Return a list of all orbits of ``self``, up to
        cyclic rotation.
        
        EXAMPLES::
        
            sage: D = InvertibleFiniteDynamicalSystem(tuple(range(6)), lambda x : (x + 2) % 6)
            sage: D.orbits()
            [[5, 1, 3], [4, 0, 2]]
            sage: D = InvertibleFiniteDynamicalSystem(tuple(range(6)), lambda x : (x + 3) % 6)
            sage: D.orbits()
            [[5, 2], [4, 1], [3, 0]]
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
        r"""
        Return a list of the lengths of all orbits of
        ``self``.
        
        EXAMPLES::
        
            sage: D = InvertibleFiniteDynamicalSystem(tuple(range(6)), lambda x : (x + 2) % 6)
            sage: D.orbit_lengths()
            [3, 3]
            sage: D = InvertibleFiniteDynamicalSystem(tuple(range(6)), lambda x : (x + 3) % 6)
            sage: D.orbit_lengths()
            [2, 2, 2]
        """
        return [len(orb) for orb in self.orbits()]

    def is_homomesic(self, h, average=None, find_average=False):
        r"""
        Check if ``h`` (a map from the ground set of ``self`` to
        a `\QQ`-vector space) is homomesic with respect to ``self``.

        If the optional argument ``average`` is provided, then
        this also checks that the averages are equal to ``average``.

        If the optional argument ``find_average`` is set to
        ``True``, then this method returns the average of ``h``
        in case ``h`` is homomesic (instead of returning ``True``).

        EXAMPLES::

            sage: W = Words(2, 5)
            sage: F = InvertibleFiniteDynamicalSystem(W, lambda x : x[1:] + Word([x[0]]))
            sage: F.is_homomesic(lambda w: sum(w))
            False
            sage: F.is_homomesic(lambda w: 1, average=1)
            True
            sage: F.is_homomesic(lambda w: 1, average=0)
            False
            sage: F.is_homomesic(lambda w: 1)
            True
            sage: F.is_homomesic(lambda w: 1, find_average=True)
            1
            sage: F.is_homomesic(lambda w: w[0] - w[1], average=0)
            True
            sage: F.is_homomesic(lambda w: w[0] - w[1], find_average=True)
            0
        """
        orbavgs = []
        for orb in self.orbits():
            l = len(orb)
            avg = ~(QQ(l)) * sum(h(i) for i in orb)
            if avg not in orbavgs:
                if orbavgs:
                    return False
                orbavgs.append(avg)
        if not orbavgs:
            return True
        if average is None:
            if find_average:
                return orbavgs[0]
            return True
        return orbavgs[0] == average

class discrete_dynamical_systems():
    r"""
    A class consisting of constructors for several specific
    discrete dynamical systems.
    """

    @staticmethod
    def permutation(pi, invertible=True):
        r"""
        Return the invertible finite discrete dynamical system
        induced by the permutation ``pi`` of the set
        `\{1, 2, \ldots, n\}`.

        The permutation ``pi`` can be provided as an instance
        of :class:`Permutation`, but also as anything that can
        be cast into that class.

        See also :meth:`one_line` for a more general
        construction, which takes any map from
        `\{1, 2, \ldots, n\}` to `\{1, 2, \ldots, n\}` given in
        one-line notation, and builds a (not necessarily
        invertible) discrete dynamical system out of it.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.permutation([3, 5, 4, 1, 2])
            sage: F.verify_inverse_evolution()
            True
            sage: sorted(F.orbit_lengths())
            [2, 3]
            sage: F.orbit(3)
            [3, 4, 1]
            sage: F.is_homomesic(lambda x: 1)
            True
            sage: F.is_homomesic(lambda x: x)
            False
        """
        from sage.combinat.permutation import Permutation
        pi = Permutation(pi)
        n = len(pi)
        X = range(1, n+1)
        return InvertibleFiniteDynamicalSystem(X, pi, inverse=pi.inverse(), create_tuple=True)

    @staticmethod
    def one_line(xs):
        r"""
        Return the finite discrete dynamical system
        with ground set `\{1, 2, \ldots, n\}` and evolution
        sending each `i` to `x_i`,
        where `(x_1, x_2, \ldots, x_n)` is the argument ``xs``
        provided.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.one_line([2, 2, 1, 2, 3])
            sage: F.orbit(3)
            [3, 1, 2]
            sage: F.orbit(5)
            [5, 3, 1, 2]
            sage: F.orbit(2)
            [2]
        """
        n = len(xs)
        X = range(1, n+1)
        xs2 = tuple(xs)
        def pi(i):
            return xs2[i - 1]
        return FiniteDynamicalSystem(X, pi, create_tuple=True)

    @staticmethod
    def bitstring_rotation(n, ones=None):
        r"""
        Return the invertible finite discrete dynamical system
        consisting of all bitstrings of size `n` (that is, of
        all `n`-tuples `(i_1, \ldots, i_n) \in \{0, 1\}^n`),
        evolving by cyclic rotation.

        If the optional parameter ``ones`` is provided, the
        system is restricted only to those bitstrings whose
        number of ones is the value of this parameter.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.bitstring_rotation(5)
            sage: sorted(F.orbit_lengths())
            [1, 1, 5, 5, 5, 5, 5, 5]
            sage: F.orbit((0, 1, 1, 0, 1))
            [(0, 1, 1, 0, 1),
             (1, 1, 0, 1, 0),
             (1, 0, 1, 0, 1),
             (0, 1, 0, 1, 1),
             (1, 0, 1, 1, 0)]
            sage: F.is_homomesic(lambda x: sum(1 for i in range(5) for j in range(i) if x[j] > x[i]))
            False
            sage: F.is_homomesic(lambda x: x[0])
            False
            sage: F = discrete_dynamical_systems.bitstring_rotation(5, ones=3)
            sage: F.verify_inverse_evolution()
            True
            sage: sorted(F.orbit_lengths())
            [5, 5]
            sage: F.orbit((0, 1, 1, 0, 1))
            [(0, 1, 1, 0, 1),
             (1, 1, 0, 1, 0),
             (1, 0, 1, 0, 1),
             (0, 1, 0, 1, 1),
             (1, 0, 1, 1, 0)]
            sage: F.is_homomesic(lambda x: sum(1 for i in range(5) for j in range(i) if x[j] > x[i]))
            True
            sage: F.is_homomesic(lambda x: x[0])
            True
        """
        from sage.categories.cartesian_product import cartesian_product
        X = cartesian_product([[0,1]] * n)
        if ones is not None:
            # Not the best method...
            X = [x for x in X if sum(x) == ones]
        if n == 0:
            phi = lambda x : x
            psi = phi
        else:
            phi = lambda x : x[1:] + (x[0],)
            psi = lambda x : (x[-1],) + x[:-1]
        return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)

    @staticmethod
    def striker_sweep(E, predicate, elements, lazy=False):
        r"""
        Return the invertible finite discrete dynamical system
        on all subsets of a finite set ``E`` satisfying a
        boolean predicate ``predicate``, where evolution is
        the "Striker sweep" -- i.e., the composition of
        Striker toggles corresponding to the elements in the
        iterable ``elements`` (from first to last).

        Let `E` be a finite set.
        Let `\mathcal{L}` be a subset of the powerset of `E`.
        (In this implementation, `\mathcal{L}` should be
        specified via the boolean predicate ``predicate``,
        which takes a subset `F` of `E` and returns the truth
        value of `F \in \mathcal{L}`.)
        For any `e \in E`, the *Striker toggle* `t_e` is
        the involution of the set `\mathcal{L}` that sends
        each `F \in \mathcal{L}` to the symmetric difference
        `F \triangle \{ e \}` if this symmetric difference
        is in `\mathcal{L}`, and otherwise to `F` itself.
        If `(e_1, e_2, \ldots, e_k)` is a finite sequence of
        elements of `E` (to be provided as a tuple, via the
        argument ``elements``), then the *Striker sweep*
        corresponding to this sequence is the composition of
        maps
        `t_{e_k} \circ t_{e_{k-1}} \circ \cdots \circ t_{e_1}`.

        This generalizes classical constructions such as
        rowmotion on order ideals.

        The optional argument ``lazy`` can be set to
        ``True``; in that case, the ground set of the
        dynamical system will not be explicitly computed.

        .. WARNING::

            The sequence ``elements`` should be provided as
            a tuple, or as a list that is guaranteed not to
            be mutated (as otherwise, mutation will corrupt
            this DDS).

        .. TODO::

            Implement the ``lazy=True`` case.
            This should differ in that ``X`` is no longer a
            list, but an enumerated set.
            But how to build an enumerated set filtering a
            given enumerated set according to a predicate?

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: StS = discrete_dynamical_systems.striker_sweep
            sage: E = range(1, 5)
            sage: lac = lambda S: all(s + 1 not in S for s in S) # lacunarity predicate
            sage: F = StS(E, lac, [1, 2, 3, 4])
            sage: F.ground_set()
            [{}, {1}, {2}, {3}, {4}, {1, 3}, {1, 4}, {2, 4}]
            sage: F.evolution()(Set([2, 4]))
            {}
            sage: F.evolution()(Set([]))
            {1, 3}
            sage: F.evolution()(Set([1, 3]))
            {4}
            sage: F.evolution()(Set([4]))
            {1}
            sage: F.inverse_evolution()(Set([1]))
            {4}
            sage: F.verify_inverse_evolution()
            True
            sage: sorted(F.orbit_lengths())
            [3, 5]
            sage: F.orbit(Set([2, 4]))
            [{2, 4}, {}, {1, 3}, {4}, {1}]
            sage: F.is_homomesic(lambda S: S.cardinality())
            False
            sage: F.is_homomesic(lambda S: bool(1 in S) - bool(4 in S), find_average=True)
            0
            sage: F.is_homomesic(lambda S: bool(2 in S) - bool(3 in S), find_average=True)
            0
            sage: F.is_homomesic(lambda S: bool(1 in S), find_average=True)
            False
        """
        from sage.combinat.subset import Subsets
        from sage.sets.set import Set
        X = [F for F in Subsets(E) if predicate(F)]
        def phi(F):
            for e in elements:
                G = F.symmetric_difference(Set([e]))
                if predicate(G):
                    F = G
            return F
        def psi(F):
            for e in reversed(elements):
                G = F.symmetric_difference(Set([e]))
                if predicate(G):
                    F = G
            return F
        return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)

    @staticmethod
    def syt_promotion(lam):
        r"""
        Return the invertible finite discrete dynamical system
        consisting of all standard tableaux of shape ``lam`` (a
        given partition) and evolving according to promotion.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.syt_promotion([4, 4, 4])
            sage: sorted(F.orbit_lengths())
            [3, 3, 4, 4, 4, 6, 6, 6, 6, 12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
            12, 12, 12, 12, 12]
            sage: G = discrete_dynamical_systems.syt_promotion([4, 3, 1])
            sage: sorted(G.orbit_lengths())
            [16, 22, 32]
            sage: G.verify_inverse_evolution()
            True
        """
        from sage.combinat.partition import Partition
        lam = Partition(lam)
        from sage.combinat.tableau import StandardTableaux
        X = StandardTableaux(lam)
        return InvertibleFiniteDynamicalSystem(X, lambda T : T.promotion(), inverse=lambda T : T.promotion_inverse())

    @staticmethod
    def order_ideal_rowmotion(P):
        r"""
        Return the invertible finite discrete dynamical system
        consisting of all order ideals of the poset ``P``,
        evolving according to rowmotion.

        EXAMPLES::

            sage: P = RootSystem(["A", 6]).root_poset()
            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.order_ideal_rowmotion(P)
            sage: sorted(F.orbit_lengths())
            [2, 7, 7, 7, 7, 7, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
            14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
            sage: F.is_homomesic(lambda I: len(I))
            False
            sage: F.is_homomesic(lambda I: sum((-1)**(P.rank(i)) for i in I))
            True

            sage: P = RootSystem(["A", 3]).root_poset()
            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: F = discrete_dynamical_systems.order_ideal_rowmotion(P)
            sage: F.verify_inverse_evolution()
            True
        """
        from sage.sets.set import Set
        X = [Set(P.order_ideal(A)) for A in P.antichains()]
        # Using P.order_ideals_lattice() instead causes intransparency issues:
        # sage can't always do P.rowmotion(I) when I is in P.order_ideals_lattice().
        # Bug in P.order_ideals_lattice() when P is facade?
        phi = lambda I : P.rowmotion(I)
        def psi(I): # inverse of rowmotion
            result = I
            for i in P.linear_extension():
                result = P.order_ideal_toggle(result, i)
            return result
        return InvertibleFiniteDynamicalSystem(X, phi, inverse=psi)

    @staticmethod
    def bulgarian_solitaire(n):
        r"""
        Return the finite discrete dynamical system defined
        by Bulgarian solitaire on partitions of size `n`.

        Let `n` be a nonnegative integer.
        Let `P` be the set of all integer partitions of
        size `n`.
        Let `B : P \to P` be the map that sends each
        partition
        `\lambda = (\lambda_1, \lambda_2, \ldots, \lambda_k)`
        of `n` (with all the `\lambda_i` positive) to
        `(\lambda_1 - 1, \lambda_2 - 1, \ldots, \lambda_k - 1, k)`,
        where zero entries have been removed and the remaining
        entries sorted into decreasing order.
        (For example,
        `B(5, 4, 2, 2, 1, 1) = (6, 4, 3, 1, 1)`.)
        This method yields the finite DDS whose ground set
        is `P` and whose evolution is `B`.

        EXAMPLES::

            sage: from sage.combinat.finite_dynamical_system import discrete_dynamical_systems
            sage: BS = discrete_dynamical_systems.bulgarian_solitaire
            sage: BS(3).evolution()(Partition([3]))
            [2, 1]
            sage: BS(3).evolution()(Partition([2, 1]))
            [2, 1]
            sage: BS(3).evolution()(Partition([1, 1, 1]))
            [3]
            sage: BS(4).evolution()(Partition([4]))
            [3, 1]
            sage: BS(4).orbit(Partition([4]))
            [[4], [3, 1], [2, 2], [2, 1, 1]]
            sage: BS(4).orbit(Partition([3, 1]))
            [[3, 1], [2, 2], [2, 1, 1]]
            sage: BS(7).orbit(Partition([6, 1]), preperiod=True)
            ([[6, 1], [5, 2], [4, 2, 1], [3, 3, 1], [3, 2, 2], [3, 2, 1, 1]], 2)
        """
        from sage.combinat.partition import Partition, Partitions
        X = Partitions(n)
        def phi(lam):
            mu = [p - 1 for p in lam if p > 0]
            nu = list(reversed(sorted(mu + [len(lam)])))
            return Partition(nu)
        return FiniteDynamicalSystem(X, phi)

