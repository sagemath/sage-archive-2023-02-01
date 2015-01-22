r"""
Evenly distributed sets in finite fields

This module consists of a simple class :class:`EvenlyDistributedSetsBacktracker`. Its main
purpose is to iterate through the evenly distributed sets of a given finite
field.

The naive backtracker implemented here is not directly used to generate
difference family as even for small parameters it already takes time to run.
Instead, its output has been stored into a database
:mod:`sage.combinat.designs.database`. If the backtracker is improved, then one
might want to update this database with more values.
"""

include "sage/ext/stdsage.pxi"

cimport cython

from libc.limits cimport UINT_MAX
from libc.string cimport memset
from libc.stdlib cimport qsort

from sage.ext.memory cimport check_malloc
from sage.libs.gmp.types cimport mpz_t

from sage.rings.integer cimport Integer,smallInteger

cdef class EvenlyDistributedSetsBacktracker:
    r"""
    Set of evenly distributed subsets in finite fields.

    Let `K` be a finite field of cardinality `q` and `k` an integer so that
    `k(k-1)` divides `q-1`. Let `H = K^*` be the multiplicative group of
    invertible elements in `K`. A `k`-*evenly distributed set* in `K` is a set
    `B = \{b_1, b_2, \ldots, b_k\}` of `k` elements of `K` so that the `k(k-1)`
    differences `\Delta B = \{b_i - b_j; i \not= j\}` fill each coset modulo
    `H^{2(q-1)/(k(k-1))}` exactly twice.

    Evenly distributed sets were introduced by Wilson [Wi72]_. He proved that
    for any `k`, and for any prime power `q` large enough such that `k(k-1)`
    divides `k` there exists a `k`-evenly distributed set in the field of
    cardinality `q`. This existence result based on a counting argument does
    not provide a simple method to generate them.

    From a `k`-evenly distributed set, it is straightforward to build a
    `(q,k,1)`-difference family. Another approach to generate difference family,
    somehow dual to this one, is via radical difference family (see in
    particular
    :func:`~sage.combinat.designs.difference_family.radical_difference_family`
    from the module :mod:`~sage.combinat.designs.difference_family`).

    INPUT:

    - ``K`` -- a finite field of cardinality `q`

    - ``k`` -- a positive integer such that `k(k-1)` divides `q-1`

    - ``up_to_isomorphism`` - (boolean, default ``True``) whether only consider
      evenly distributed set up to automorphisms of the field of the form
      `x \mapsto ax + b`.

    - ``check`` -- boolean (default is ``False``). Whether you want to check
      intermediate steps of the iterator. This is mainly intended for debugging
      purpose. Setting it to ``True`` will considerably slow the iteration.

    EXAMPLES::

        sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

    The main part of the code is contained in the iterator. To get one evenly
    distributed set just do::

        sage: E = EvenlyDistributedSetsBacktracker(Zmod(151),6)
        sage: B = E.an_element()
        sage: B
        [0, 1, 69, 36, 57, 89]

    The class has a method to convert it to a difference family::

        sage: E.to_difference_family(B)
        [[0, 1, 69, 36, 57, 89],
         [0, 132, 48, 71, 125, 121],
         [0, 59, 145, 10, 41, 117],
         [0, 87, 114, 112, 127, 42],
         [0, 8, 99, 137, 3, 108]]

    It is also possible to run over all evenly distributed sets::

        sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=False)
        sage: for B in E: print B
        [0, 1, 11, 5]
        [0, 1, 4, 6]
        [0, 1, 9, 3]
        [0, 1, 8, 10]


        sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=True)
        sage: for B in E: print B
        [0, 1, 11, 5]
    """
    # PYTHON DATA
    cdef K                       # the underlying field
    cdef list list_K             # the elements of K  (i -> x)
    cdef dict K_to_int           # inverse of list_K  (x -> i)

    # FLAGS
    cdef int count               # do we count or do we iterate
    cdef int check               # do we need to check (debug)
    cdef int up_to_isom          # do we care only about isomorphisms

    # STATIC DATA
    cdef unsigned int q          # cardinality of the field
    cdef unsigned int k          # size of the subsets
    cdef unsigned int e          # k(k-1)/2
    cdef unsigned int ** diff    # qxq array: diff[x][y]  = x - y
    cdef unsigned int ** ratio   # qxq array: ratio[x][y] = x / y

    # DYNAMIC DATA
    cdef unsigned int * B        # current stack of elements of {0,...,q-1}
    cdef unsigned int * cosets   # cosets of differences of elts in B
    cdef unsigned int * t        # temporary variable for updates

    def __dealloc__(self):
        if self.diff != NULL:
            sage_free(self.diff[0])
            sage_free(self.diff)
        if self.ratio != NULL:
            sage_free(self.ratio[0])
            sage_free(self.ratio)
        sage_free(self.B)
        sage_free(self.cosets)
        sage_free(self.t)

    def __init__(self, K, k, up_to_isomorphism=True, check=False):
        r"""
        TESTS::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: EvenlyDistributedSetsBacktracker(Zmod(4),2)
            Traceback (most recent call last):
            ...
            ValueError: Ring of integers modulo 4 is not a field

            sage: EvenlyDistributedSetsBacktracker(Zmod(71),7)
            Traceback (most recent call last):
            ...
            ValueError: k(k-1)=42 does not divide q-1=71
        """
        self.check      = 1 if check else 0
        self.up_to_isom = 1 if up_to_isomorphism else 0
        self.count      = 0

        cdef unsigned int i,j,ell

        if not K.is_field():
            raise ValueError("{} is not a field".format(K))
        cdef unsigned int q = K.cardinality()
        cdef unsigned int e = k*(k-1)/2
        if (q-1) % (2*e) != 0:
            raise ValueError("k(k-1)={} does not divide q-1={}".format(k*(k-1),q))
        cdef unsigned int m = (q-1)/e

        self.q = q
        self.e = e
        self.k = k
        self.K = K

        self.diff = <unsigned int **> check_malloc(q*sizeof(unsigned int *))
        self.diff[0] = <unsigned int *> check_malloc(q*q*sizeof(unsigned int))
        for i in range(1,self.q):
            self.diff[i] = self.diff[i-1] + q
        self.ratio = <unsigned int **> check_malloc(q*sizeof(unsigned int *))
        self.ratio[0] = <unsigned int *> check_malloc(q*q*sizeof(unsigned int))
        for i in range(1,self.q):
            self.ratio[i] = self.ratio[i-1] + q

        self.B  = <unsigned int *> check_malloc(k*sizeof(unsigned int))
        self.cosets = <unsigned int *> check_malloc(e*sizeof(unsigned int))
        self.t = <unsigned int *> check_malloc(e*sizeof(unsigned int))

        x = K.multiplicative_generator()
        list_K = []
        for i in range(e):
            list_K.extend(sorted(x**(j*e+i) for j in range(m)))
        list_K.append(K.zero())
        self.list_K = list_K
        K_to_int = self.K_to_int = {y:i for i,y in enumerate(list_K)}
        assert self.K_to_int[K.zero()] == q-1
        assert self.K_to_int[K.one()] == 0
        assert set(K) == set(list_K)

        for i,x in enumerate(self.list_K):
            for j,y in enumerate(self.list_K):
                self.diff[i][j]  = K_to_int[x-y]
                if y:
                    self.ratio[i][j] = K_to_int[x/y]
                else:
                    self.ratio[i][j] = UINT_MAX

    def to_difference_family(self, B):
        r"""
        Given an evenly distributed set ``B`` convert it to a difference family.

        This is useful if you want to obtain the difference family from the
        output of the iterator.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(41),5)
            sage: B = E.an_element(); B
            [0, 1, 13, 38, 31]
            sage: D = E.to_difference_family(B); D
            [[0, 1, 13, 38, 31], [0, 32, 6, 27, 8]]

            sage: from sage.combinat.designs.difference_family import is_difference_family
            sage: is_difference_family(Zmod(41),D,41,5,1)
            True
        """
        xe = self.K.multiplicative_generator() ** (self.e)
        return [[xe**j*b for b in B] for j in range((self.q-1)/(2*self.e))]

    def an_element(self):
        r"""
        Return an evenly distributed set.

        If there is no such subset raise an ``EmptySetError``.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(41),5)
            sage: E.an_element()
            [0, 1, 13, 38, 31]

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(61),6)
            sage: E.an_element()
            Traceback (most recent call last):
            ...
            EmptySetError: no 6-evenly distributed set in Ring of integers modulo 61
        """
        from sage.categories.sets_cat import EmptySetError
        it = iter(self)
        try:
            B = it.next()
        except StopIteration:
            raise EmptySetError("no {}-evenly distributed set in {}".format(self.k,self.K))
        return B

    def __repr__(self):
        r"""
        A string representative.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: EvenlyDistributedSetsBacktracker(GF(25,'a'), 4)
            4-evenly distributed sets (up to isomorphism) in Finite Field in a of size 5^2
            sage: EvenlyDistributedSetsBacktracker(GF(25,'a'), 4, up_to_isomorphism=False)
            4-evenly distributed sets in Finite Field in a of size 5^2
        """
        return "{}-evenly distributed sets {} in {}".format(
                self.k,
                '(up to isomorphism)' if self.up_to_isom else '',
                self.K)

    def cardinality(self):
        r"""
        Return the number of evenly distributed sets.

        Use with precaution as there can be a lot of such sets and this method
        might be very long to answer!

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(GF(25,'a'),4)
            sage: E
            4-evenly distributed sets (up to isomorphism) in Finite Field in a of size 5^2
            sage: E.cardinality()
            4

            sage: E = EvenlyDistributedSetsBacktracker(GF(25,'a'), 4, up_to_isomorphism=False)
            sage: E.cardinality()
            40
        """
        cdef n = 0
        self.count = 1
        for a in self:
            n += a
        self.count = 0
        return smallInteger(n)

    def _B_automorphisms(self):
        r"""
        Check whether B is minimal among its relabelization.

        This is an internal function and should only be call by the backtracker
        implemented in the method `__iter__`.

        OUTPUT:

        - ``False`` if ``self.B`` is not minimal

        - ``True`` if ``self.B`` is minimal and ``self.up_to_isom = 1``

        - the list of relabeled versions if ``self.B`` is minimal and ``self.up_to_isom = 0``

        TESTS::

            sage: from sage.combinat.designs.evenly_distributed_sets import \
            ....:     EvenlyDistributedSetsBacktracker
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=True)
            sage: E.cardinality()   # indirect doctest
            1
            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13), 4, up_to_isomorphism=False)
            sage: E.cardinality()   # indirect doctest
            4
        """
        cdef unsigned int i,j,k,tmp1,tmp2,verify
        cdef list B = [self.B[i] for i in range(1,self.k)]
        B.append(self.q-1)
        cdef list BB = [None]*self.k
        cdef set relabs = set([tuple(B)])

        #  z -> (z - B[i]) / (B[j] - B[i])
        for i in range(self.k):
            for j in range(self.k):
                if i == j:
                    continue
                tmp1 = self.diff[self.B[j]][self.B[i]]

                verify = 0
                for k in range(self.k):
                    if k == i:
                        BB[k] = self.q-1
                    elif k == j:
                        BB[k] = 0
                    else:
                        tmp2 = self.ratio[self.diff[self.B[k]][self.B[i]]][tmp1]
                        if tmp2 == 0 or tmp2 == self.q-1 or tmp2 < self.B[2]:
                            raise RuntimeError("there is a problem got tmp2={}".format(tmp2,self.B[2]))
                        elif tmp2 == self.B[2]:
                            verify = 1
                        BB[k] = tmp2

                if verify:
                    BB.sort()
                    if BB < B:
                        return False

                if not self.up_to_isom:
                    if not verify:
                        BB.sort()
                    if BB > B:
                        relabs.add(tuple(BB))

        return sorted(relabs)

    @cython.cdivision(True)
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __iter__(self):
        r"""
        Iterator through all evenly distributed sets that start with `[0,1]`.

        EXAMPLES::

            sage: from sage.combinat.designs.evenly_distributed_sets import EvenlyDistributedSetsBacktracker

            sage: E = EvenlyDistributedSetsBacktracker(Zmod(13),4)
            sage: for B in E:
            ....:     print B
            [0, 1, 11, 5]
        """
        cdef unsigned int k = self.k
        cdef unsigned int q = self.q
        cdef unsigned int e = self.e
        cdef unsigned int m = (q-1) / e

        # in the list B we store the candidate for being an e.d.s.
        # we always have B[0] = 0 and B[1] = 1
        # because 0 is in B, the cosets of the elements of B must be
        # disjoint.
        cdef unsigned int kk = 2
        cdef unsigned int * B  = self.B
        B[0] = q-1  # the element 0 in K
        B[1] = 0    # the element 1 in K

        memset(self.cosets, 0, e * sizeof(unsigned int))
        memset(self.t, 0, e * sizeof(unsigned int))

        self.cosets[0] = 1  # the coset of 1

        cdef unsigned int x = m
        while True:
            if self.check:
                self._check(kk)
                assert m < x < q-1, "got x < m or x > q where x={}".format(x)

            # try to append x
            if self._add_element(x,kk):
                # note: the element x is already added to B in ._add_element()
                if kk == k-1:
                    ans = self._B_automorphisms()
                    if ans is False:
                        continue
                    elif self.count:
                        yield len(ans)
                    else:
                        for a in ans:
                            yield [self.list_K[q-1]] + [self.list_K[a[r]] for r in range(k-1)]
                    self._pop(kk)
                else:
                    # we need to choose the next element from a new coset. In
                    # order to jump we artificially set x to the last element of
                    # the coset.
                    x += m - x%m - 1
                    kk += 1

            # now we determine the next element x to be tested
            if x == q-2:
                kk -= 1
                x = B[kk]
                self._pop(kk)
                if x == q-2:
                    kk -= 1
                    x = B[kk]
                    self._pop(kk)
                    if x == q-2:
                        raise RuntimeError("this is impossible!")

            if kk == 1:
                return

            x += 1

    @cython.cdivision(True)
    cdef inline int _add_element(self, unsigned int x, unsigned int kk) except -1:
        r"""
        Add the element ``x`` to ``B`` in position kk if the resulting set is
        still evenly distributed.
        """
        cdef unsigned int i,j
        cdef unsigned int tmp
        cdef unsigned int m = (self.q-1)/self.e

        # we check that we do not get a smaller evenly distributed set by a
        # relabeling of the form
        #   1. z -> (z - x) / (B[i] - x)    that maps x -> 0 and B[i] -> 1
        #   2. z -> (z - B[i]) / (x - B[i]) that maps B[i] -> 0 and x -> 1
        if kk > 2:
            for i in range(kk):
                # isom of form 1
                tmp = self.diff[self.B[i]][x]
                for j in range(kk):
                    if j == i:
                        continue
                    if self.ratio[tmp][self.diff[self.B[j]][x]] < self.B[2]:
                        return 0

                # isom of form 2
                tmp = self.diff[x][self.B[i]]
                for j in range(kk):
                    if j == i:
                        continue
                    if self.ratio[self.diff[self.B[j]][self.B[i]]][tmp] < self.B[2]:
                        return 0

        for j in range(kk):
            i = self.diff[x][self.B[j]] / m
            if self.cosets[i] or self.t[i]:
                memset(self.t, 0, self.e*sizeof(unsigned int))
                return 0
            self.t[i] = 1

        self.B[kk] = x
        for i in range(self.e):
            if self.t[i]:
                self.cosets[i] = 1
                self.t[i] = 0
        return 1

    cdef inline void _pop(self, unsigned int kk):
        r"""
        Pop the element of ``self.B`` at position ``kk`` and updates
        ``self.cosets``
        """
        cdef unsigned int i,j,x
        cdef unsigned int m = (self.q-1)/self.e
        x = self.B[kk]
        for j in range(kk):
            i = self.diff[x][self.B[j]] / m
            self.cosets[i] = 0

    cdef int _check(self, unsigned int kk) except -1:
        r"""
        Sanity check (only for debug purposes).
        """
        cdef unsigned int i,j,x,y
        cdef set s = set()
        for i in range(kk):
            x = self.B[i]
            for j in range(i):
                y = self.B[j]
                s.add(self.diff_to_coset[x][y])

        if set(j for j in range(self.e) if self.cosets[j]) != s:
            raise RuntimeError("self.cosets is not synchronized with self.B!")
