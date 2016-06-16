# -*- coding: utf-8 -*-
r"""
Finite simplicial sets

AUTHORS:

- John H. Palmieri (2016-05)

This module implements finite simplicial sets: simplicial sets with
finitely many non-degenerate simplices.

A *simplicial set* `X` is a collection of sets `X_n` indexed by the
non-negative integers; the set `X_n` is called the set of
`n`-simplices. These sets are connected by maps

.. MATH::

    d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
    s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

satisfying the *simplicial identities*:

.. MATH::

    d_i d_j &= d_{j-1} d_i \ \  \text{if } i<j \\
    d_i s_j &= s_{j-1} d_i \ \  \text{if } i<j \\
    d_j s_j &= 1 = d_{j+1} s_j \\
    d_i s_j &= s_{j} d_{i-1} \ \  \text{if } i>j+1 \\
    s_i s_j &= s_{j+1} s_{i} \ \  \text{if } i<j+1

See :wikipedia:`Simplicial_set`, Peter May's seminal book [May]_, or
Greg Friedman's "Illustrated introduction" :arxiv:`0809.4221` for more
information.

Several simplicial sets are predefined, and users can construct others
either by hand (using :class:`SimplicialSet`) or from existing ones
using pushouts, pullbacks, etc.

EXAMPLES:

Some of the predefined simplicial sets::

    sage: simplicial_sets.Torus()
    Torus
    sage: simplicial_sets.RealProjectiveSpace(7)
    RP^7
    sage: S5 = simplicial_sets.Sphere(5)
    sage: S5
    S^5
    sage: S5.nondegenerate_simplices()
    [v_0, sigma_5]

Type ``simplicial_sets.`` and hit the ``TAB`` key to get a full
list. To define a simplicial set by hand, first define some simplices,
then use them to construct a simplicial set by specifying their
faces::

    sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
    sage: v = AbstractSimplex(0, name='v')
    sage: w = AbstractSimplex(0, name='w')
    sage: e = AbstractSimplex(1, name='e')
    sage: f = AbstractSimplex(1, name='f')
    sage: X = SimplicialSet({e: (v,w), f: (w,w)})

Now `e` is an edge from `v` to `w` and `f` is an edge starting and
ending at `w`. Therefore the first homology group of `X` should be a
copy of the integers::

    sage: X.homology(1)
    Z

For most simplicial sets (the ``Point`` is the main exception), each
time it is constructed, it gives a distinct copy, and two distinct
simplicial sets are never equal::

    sage: T = simplicial_sets.Torus()
    sage: T == simplicial_sets.Torus()
    False
    sage: T == T
    True
    sage: simplicial_sets.Torus() == simplicial_sets.Torus()
    False
    sage: simplicial_sets.Point() == simplicial_sets.Point()
    True

You can construct subsimplical sets by specifying a list of simplices,
and then you can define the quotient simplicial set::

    sage: v = AbstractSimplex(0, name='v')
    sage: w = AbstractSimplex(0, name='w')
    sage: e = AbstractSimplex(1, name='e')
    sage: f = AbstractSimplex(1, name='f')
    sage: X = SimplicialSet({e: (v,w), f: (w,w)})
    sage: Y = X.subsimplicial_set([v,w])
    sage: Z = X.quotient(Y)  # identify v and w
    sage: Z.homology(1)
    Z x Z

    sage: K = simplicial_sets.Simplex(3)
    sage: L = K.n_skeleton(2) # a subsimplicial set of K
    sage: M = K.quotient(L)
    sage: M
    Simplicial set with 2 non-degenerate simplices

Note that subsimplicial sets and quotients come equipped with
inclusion and quotient morphisms::

    sage: inc = L.inclusion_map()
    sage: inc.domain() == L and inc.codomain() == K
    True
    sage: q = M.quotient_map()
    sage: q.domain()
    3-simplex
    sage: q.codomain() == M
    True

You can construct new simplicial sets from old by taking disjoint
unions, wedges (if they are pointed), products, pushouts, pullbacks,
cones, and suspensions, most of which also have maps associated with
them. Wedges, for example::

    sage: T = simplicial_sets.Torus()
    sage: S3 = simplicial_sets.Sphere(3)
    sage: T.is_pointed() and S3.is_pointed()
    True
    sage: T.wedge(S3)
    Simplicial set with 7 non-degenerate simplices
    sage: T.disjoint_union(S3) == T.coproduct(S3)
    False

    sage: W = T.wedge(S3)
    sage: W.inclusion(0).domain()
    Torus
    sage: W.projection(1).codomain()
    Simplicial set with 2 non-degenerate simplices

If the `5`-sphere were not already available via
``simplicial_sets.Sphere(5)``, you could construct it as follows::

    sage: pt = simplicial_sets.Simplex(0)
    sage: edge = pt.cone()
    sage: S1 = edge.quotient(edge.n_skeleton(0))

At this point, ``S1`` is pointed: every quotient is automatically
given a base point, namely the image of the subcomplex. So its
suspension is the reduced suspension, and therefore is small::

    sage: S5 = S1.suspension(4)
    sage: S5
    Simplicial set with 2 non-degenerate simplices

If we forget about the base point in ``S1``, we would get the
unreduced suspension instead. The unreduced suspension is constructed
as the quotient `CX/X` and hence is again pointed, so if we want to
keep getting unreduced suspensions, we need to forget the base point
each time::

    sage: Z1 = S1.unset_base_point()
    sage: Z2 = Z1.suspension().unset_base_point()
    sage: Z3 = Z2.suspension().unset_base_point()
    sage: Z4 = Z3.suspension().unset_base_point()
    sage: Z5 = Z4.suspension().unset_base_point()
    sage: Z5
    Simplicial set with 10 non-degenerate simplices

The cone on a pointed simplicial set is the reduced cone. The
`n`-simplex in Sage is not pointed, but the simplicial set ``Point``
is. ::

    sage: simplicial_sets.Simplex(0).cone()
    Simplicial set with 3 non-degenerate simplices
    sage: simplicial_sets.Point().cone()
    Simplicial set with 1 non-degenerate simplex

Finally, you can compute homology groups and the fundamental group of
any simplicial set::

    sage: S1 = simplicial_sets.Sphere(1)
    sage: eight = S1.wedge(S1)
    sage: eight.fundamental_group()
    Finitely presented group < e0, e1 | >

    sage: Sigma3 = groups.permutation.Symmetric(3)
    sage: BSigma3 = Sigma3.nerve_n_skeleton(3)
    sage: pi = BSigma3.fundamental_group(); pi
    Finitely presented group < e0, e1 | e0^2, e1^3, (e0*e1^-1)^2 >
    sage: pi.order()
    6
    sage: pi.is_abelian()
    False

    sage: RP6 = simplicial_sets.RealProjectiveSpace(6)
    sage: RP6.homology(reduced=False, base_ring=GF(2))
    {0: Vector space of dimension 1 over Finite Field of size 2,
     1: Vector space of dimension 1 over Finite Field of size 2,
     2: Vector space of dimension 1 over Finite Field of size 2,
     3: Vector space of dimension 1 over Finite Field of size 2,
     4: Vector space of dimension 1 over Finite Field of size 2,
     5: Vector space of dimension 1 over Finite Field of size 2,
     6: Vector space of dimension 1 over Finite Field of size 2}
    sage: RP6.homology(reduced=False, base_ring=QQ)
    {0: Vector space of dimension 1 over Rational Field,
     1: Vector space of dimension 0 over Rational Field,
     2: Vector space of dimension 0 over Rational Field,
     3: Vector space of dimension 0 over Rational Field,
     4: Vector space of dimension 0 over Rational Field,
     5: Vector space of dimension 0 over Rational Field,
     6: Vector space of dimension 0 over Rational Field}

REFERENCES:

.. [May] \J. P. May, Simplicial Objects in Algebraic Topology,
   University of Chicago Press (1967)
"""
#*****************************************************************************
#  Copyright (C) 2016 John H. Palmieri <palmieri at math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#
#*****************************************************************************

from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from sage.homology.cell_complex import GenericCellComplex
from sage.homology.simplicial_complex import SimplicialComplex
import sage.homology.simplicial_complexes_catalog as simplicial_complexes
from sage.homology.delta_complex import DeltaComplex, delta_complexes
from sage.rings.integer import Integer
from sage.rings.infinity import Infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.matrix.constructor import matrix
from sage.homology.chain_complex import ChainComplex
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.sets.set import Set
from sage.misc.cachefunc import cached_method, cached_function
from pyparsing import OneOrMore, nestedExpr
from sage.graphs.graph import Graph
from sage.env import SAGE_ENV
import copy
import itertools
import re
import os
from functools import total_ordering

from sage.misc.lazy_import import lazy_import
lazy_import('sage.categories.simplicial_sets', 'SimplicialSets')

@total_ordering
class AbstractSimplex(SageObject):
    """
    An abstract simplex, a building block of a simplicial set.

    In a simplicial set, a simplex either is non-degenerate or is
    obtained by applying degeneracy maps to a non-degenerate simplex.

    INPUT:

    - ``dim`` -- a non-negative integer, the dimension of the
      underlying non-degenerate simplex.

    - ``degeneracies`` (optional, default `()`) -- a list or tuple of
      non-negative integers, the degeneracies to be applied.

    - ``underlying`` (optional) -- a non-degenerate simplex to which
      the degeneracies are being applied.

    - ``name`` (optional) -- string, a name for this simplex.

    So to define a simplex formed by applying the degeneracy maps `s_2
    s_1` to a 1-simplex, call ``AbstractSimplex(1, (2, 1))``.

    Specify ``underlying`` if you need to keep explicit track of the
    underlying non-degenerate simplex, for example when computing
    faces of another simplex. This is mainly for use by the method
    :meth:`apply_degeneracies`.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex
        sage: AbstractSimplex(3, (3, 1))
        Simplex obtained by applying degeneracies s_3 s_1 to Non-degenerate simplex of dimension 3
        sage: AbstractSimplex(3, None)
        Non-degenerate simplex of dimension 3
        sage: AbstractSimplex(3)
        Non-degenerate simplex of dimension 3

    Simplices may be named (or renamed), affecting how they are printed::

        sage: AbstractSimplex(0)
        Non-degenerate simplex of dimension 0
        sage: v = AbstractSimplex(0, name='v')
        sage: v
        v
        sage: v.rename('w_0')
        sage: v
        w_0

    The simplicial identities are used to put the degeneracies in
    standard decreasing form::

        sage: x = AbstractSimplex(0, (0, 0, 0))
        sage: x
        Simplex obtained by applying degeneracies s_2 s_1 s_0 to Non-degenerate simplex of dimension 0
        sage: x.degeneracies()
        [2, 1, 0]

    Use of the ``underlying`` argument::

        sage: v = AbstractSimplex(0, name='v')
        sage: e = AbstractSimplex(0, (0,), underlying=v)
        sage: e
        Simplex obtained by applying degeneracy s_0 to v
        sage: e.nondegenerate() is v
        True

        sage: e.dimension()
        1
        sage: e.is_degenerate()
        True

    Distinct simplices are never equal::

        sage: AbstractSimplex(0, None) == AbstractSimplex(0, None)
        False
        sage: AbstractSimplex(0, (2,1,0)) == AbstractSimplex(0, (2,1,0))
        False

        sage: e = AbstractSimplex(0, ((0,)))
        sage: f = AbstractSimplex(0, ((0,)))
        sage: e == f
        False
        sage: e.nondegenerate() == f.nondegenerate()
        False

    This means that if, when defining a simplicial set, you specify
    the faces of a 2-simplex as::

        (e, e, e)

    then the faces are the same degenerate vertex, but if you specify
    the faces as::

        (AbstractSimplex(0, ((0,))), AbstractSimplex(0, ((0,))), AbstractSimplex(0, ((0,))))

    then the faces are three different degenerate vertices.
    """
    def __init__(self, dim, degeneracies=(), underlying=None, name=None):
        """
        A simplex of dimension ``dim``.

        INPUT:

        - ``dim`` -- integer, the dimension
        - ``degeneracies`` (optional) -- iterable, the indices of the degeneracy maps
        - ``underlying`` (optional) -- a non-degenerate simplex
        - ``name`` (optional) -- string

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (1,2))
            Simplex obtained by applying degeneracies s_3 s_1 to Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, None)
            Non-degenerate simplex of dimension 3

            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(0, (0,), underlying=v)
            sage: e
            Simplex obtained by applying degeneracy s_0 to v
            sage: e.nondegenerate() is v
            True

            sage: AbstractSimplex(3.2, None)
            Traceback (most recent call last):
            ...
            ValueError: the dimension must be an integer
            sage: AbstractSimplex(-3, None)
            Traceback (most recent call last):
            ...
            ValueError: the dimension must be non-negative

            sage: AbstractSimplex(0, (1,))
            Traceback (most recent call last):
            ...
            ValueError: invalid list of degeneracy maps on 0-simplex
        """
        try:
            Integer(dim)
        except TypeError:
            raise ValueError('the dimension must be an integer')
        if dim < 0:
            raise ValueError('the dimension must be non-negative')
        self._dim = dim
        if degeneracies:
            self._degens = standardize_degeneracies(*degeneracies)
            for (d, s) in enumerate(reversed(self._degens)):
                if d + dim < s:
                    raise ValueError('invalid list of degeneracy maps '
                                     'on {}-simplex'.format(dim))
            if underlying is None:
                self._underlying = AbstractSimplex(dim, None)
            else:
                self._underlying = underlying
        else:
            self._degens = ()
            if underlying is None:
                self._underlying = self
            else:
                self._underlying = underlying
        if name is not None:
            self.rename(name)
        # _faces: a dictionary storing the faces of this simplex in any
        # given simplicial set. Each key: a simplicial set. The
        # corresponding value: the tuple of faces of this simplex, or
        # None if it is a vertex. Thus you can tell whether a simplex
        # sigma is in a simplicial set X by testing
        # 'X in sigma._faces'. This dictionary is not copied when
        # copying a simplex.
        self._faces = {}

    def __hash__(self):
        """
        If nondegenerate: return the id of this simplex.

        Otherwise, combine the id of its underlying nondegenerate
        simplex with the tuple of indeterminacies.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: hash(v) == hash(w)
            False
            sage: x = v.apply_degeneracies(2,1,0)
            sage: id(x) == id(v.apply_degeneracies(2,1,0))
            False
            sage: hash(x) == hash(v.apply_degeneracies(2,1,0))
            True
        """
        if self.is_nondegenerate():
            return id(self)
        return hash(self.nondegenerate()) ^ hash(self._degens)

    def __eq__(self, other):
        """
        Two nondegenerate simplices are equal if they are identical.
        Two degenerate simplices are equal if their underlying
        nondegenerate simplices are identical and their tuples of
        degeneracies are equal.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: v == w
            False
            sage: v.apply_degeneracies(2,1,0) is v.apply_degeneracies(2,1,0)
            False
            sage: v.apply_degeneracies(2,1,0) == v.apply_degeneracies(2,1,0)
            True
        """
        if self.is_nondegenerate() and other.is_nondegenerate():
            return self is other
        return (self._degens == other._degens
                and self.nondegenerate() is other.nondegenerate())

    def __ne__(self, other):
        """
        This returns the negation of `__eq__`.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: v != w
            True
            sage: x = v.apply_degeneracies(1, 0)
            sage: y = v.apply_degeneracies(1, 0)
            sage: x != y
            False
        """
        return not self == other

    def __lt__(self, other):
        """
        We implement sorting in the hopes that sorted lists of simplices,
        for example as defining data for a simplicial set, will be
        well-defined invariants.

        Sort by dimension first. If dimensions are equal, if only one
        has a custom name (as set by specifying ``name=NAME`` upon
        creation or by calling ``object.rename('NAME')``, put it
        first. If both have custom names, sort by the names. As a last
        resort, sort by their id.

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)

        At this point, comparision between v and w is random, based on
        their location in memory. ::

            sage: v < w and w < v
            False
            sage: v < w or w < v
            True
            sage: (v < w and w > v) or (w < v and v > w)
            True

        Now we add names to force an ordering::

            sage: w.rename('w')
            sage: v < w
            False
            sage: v > w
            True
            sage: v >= w
            True
            sage: v.rename('v')
            sage: v < w
            True
            sage: v <= w
            True
            sage: v > w
            False

        Test other sorting. Dimensions::

            sage: AbstractSimplex(0) < AbstractSimplex(3)
            True
            sage: AbstractSimplex(0, ((0,0,0))) <= AbstractSimplex(2)
            False

        Degenerate comes after non-degenerate, and if both are
        degenerate, sort on the degeneracies::

            sage: AbstractSimplex(0, ((0,0))) <= AbstractSimplex(2)
            False
            sage: AbstractSimplex(1, ((0,))) < AbstractSimplex(1, ((1,)))
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False
            sage: w.rename('a')
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            True

        Testing `<=`, `>`, `>=`::

            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: v <= v
            True
            sage: w <= v
            False
            sage: v.apply_degeneracies(1,0) <= w.apply_degeneracies(1,0)
            True

            sage: v > v
            False
            sage: w > v
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False

            sage: v >= v
            True
            sage: w >= v
            True
            sage: v.apply_degeneracies(1,0) >= w.apply_degeneracies(1,0)
            False
        """
        if self.dimension() < other.dimension():
            return True
        if self.dimension() > other.dimension():
            return False
        if self.degeneracies() and not other.degeneracies():
            return False
        if other.degeneracies() and not self.degeneracies():
            return True
        if self.degeneracies() and other.degeneracies() and self.degeneracies() != other.degeneracies():
            return self.degeneracies() < other.degeneracies()
        if hasattr(self.nondegenerate(), '__custom_name'):
            if hasattr(other.nondegenerate(), '__custom_name'):
                return str(self) < str(other)
            return True
        else:
            if hasattr(other, '__custom_name'):
                return False
        return id(self) < id(other)

    def nondegenerate(self):
        """
        The non-degenerate simplex underlying this one.

        Therefore return itself if this simplex is non-degenerate.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0, name='v')
            sage: sigma = v.apply_degeneracies(1, 0)
            sage: sigma.nondegenerate()
            v
            sage: tau = AbstractSimplex(1, (3,2,1))
            sage: x = tau.nondegenerate(); x
            Non-degenerate simplex of dimension 1
            sage: x == tau.nondegenerate()
            True

            sage: AbstractSimplex(1, None)
            Non-degenerate simplex of dimension 1
            sage: AbstractSimplex(1, None) == x
            False
            sage: AbstractSimplex(1, None) == tau.nondegenerate()
            False
        """
        return self._underlying

    def degeneracies(self):
        """
        Return the list of indices for the degeneracy maps for this
        simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(4, (0,0,0)).degeneracies()
            [2, 1, 0]
            sage: AbstractSimplex(4, None).degeneracies()
            []
        """
        return list(self._degens)

    def is_degenerate(self):
        """
        True if this simplex is degenerate.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).is_degenerate()
            True
            sage: AbstractSimplex(3, None).is_degenerate()
            False
        """
        return bool(self.degeneracies())

    def is_nondegenerate(self):
        """
        True if this simplex is non-degenerate.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).is_nondegenerate()
            False
            sage: AbstractSimplex(3, None).is_nondegenerate()
            True
            sage: AbstractSimplex(5).is_nondegenerate()
            True
        """
        return not self.is_degenerate()

    def dimension(self):
        """
        The dimension of this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, (2,1)).dimension()
            5
            sage: AbstractSimplex(3, None).dimension()
            3
            sage: AbstractSimplex(7).dimension()
            7
        """
        return self._dim + len(self.degeneracies())

    def apply_degeneracies(self, *args):
        """
        Apply the degeneracies given by the arguments ``args`` to this simplex.

        INPUT:

        - ``args`` -- integers

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: e = v.apply_degeneracies(0)
            sage: e.nondegenerate() == v
            True
            sage: f = e.apply_degeneracies(0)
            sage: f
            Simplex obtained by applying degeneracies s_1 s_0 to Non-degenerate simplex of dimension 0
            sage: f.degeneracies()
            [1, 0]
            sage: f.nondegenerate() == v
            True
            sage: v.apply_degeneracies(1, 0)
            Simplex obtained by applying degeneracies s_1 s_0 to Non-degenerate simplex of dimension 0

        TESTS::

            sage: e.apply_degeneracies() == e
            True

        Do not pass an explicit list or tuple as the argument: call
        this with the syntax ``x.apply_degeneracies(1,0)``, not
        ``x.apply_degeneracies([1,0])``::

            sage: e.apply_degeneracies([1,0])
            Traceback (most recent call last):
            ...
            TypeError: degeneracies are indexed by non-negative integers; do not use an explicit list or tuple
        """
        if not args:
            return self
        underlying = self.nondegenerate()
        return AbstractSimplex(underlying.dimension(),
                               degeneracies= list(args) + self.degeneracies(),
                               underlying=underlying)

    def __copy__(self):
        """
        Return a copy of this simplex.

        Forget the "underlying" non-degenerate simplex. If this
        simplex has a name, then its copy's name is obtained by adding
        a prime ``'`` at the end.

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: v = AbstractSimplex(0)
            sage: copy(v) == v
            False
            sage: copy(v).nondegenerate() == v
            False
            sage: x = v.apply_degeneracies(1, 0)
            sage: y = copy(v).apply_degeneracies(1, 0)
            sage: z = copy(x)
            sage: x == y or x == z or y == z
            False
            sage: x.nondegenerate() == copy(v)
            False
            sage: y.nondegenerate() == v
            False

            sage: v.rename('v')
            sage: copy(v)
            v'
            sage: copy(copy(v))
            v''
        """
        # Don't preserve the underlying simplex when copying, just the
        # dimension and the degeneracies.
        sigma = AbstractSimplex(self._dim, degeneracies=self.degeneracies())
        if hasattr(self, '__custom_name'):
            sigma.rename(str(self) + "'")
        return sigma

    def _repr_(self):
        """
        Print representation.

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: AbstractSimplex(3, None)
            Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, (0,))
            Simplex obtained by applying degeneracy s_0 to Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, (0, 0))
            Simplex obtained by applying degeneracies s_1 s_0 to Non-degenerate simplex of dimension 3

        Test renaming::

            sage: v = AbstractSimplex(0)
            sage: v
            Non-degenerate simplex of dimension 0
            sage: v.rename('v')
            sage: v
            v
            sage: v.apply_degeneracies(1, 0)
            Simplex obtained by applying degeneracies s_1 s_0 to v
        """
        if self.degeneracies():
            degens = ' '.join(['s_{}'.format(i) for i in self.degeneracies()])
            if len(self.degeneracies()) == 1:
                return 'Simplex obtained by applying degeneracy {} to {!r}'.format(degens,
                                        self.nondegenerate())
            return 'Simplex obtained by applying degeneracies {} to {!r}'.format(degens,
                                        self.nondegenerate())
        return 'Non-degenerate simplex of dimension {}'.format(self._dim)

    def _latex_(self):
        r"""
        LaTeX representation.

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: latex(AbstractSimplex(18, None))
            \Delta^{18}
            sage: latex(AbstractSimplex(3, (0, 0,)))
            s_{1} s_{0} \Delta^{3}
        """
        simplex = "\\Delta^{{{}}}".format(self._dim)
        if self.degeneracies():
            degens = ' '.join(['s_{{{}}}'.format(i) for i in self.degeneracies()])
            return degens + ' ' + simplex
        return simplex


########################################################################
# The main class

class SimplicialSet(GenericCellComplex, Parent):
    r"""
    A simplicial set.

    A simplicial set `X` is a collection of sets `X_n`, the
    *n-simplices*, indexed by the non-negative integers, together with
    maps

    .. MATH::

        d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. MATH::

        d_i d_j &= d_{j-1} d_i \ \  \text{if } i<j \\
        d_i s_j &= s_{j-1} d_i \ \  \text{if } i<j \\
        d_j s_j &= 1 = d_{j+1} s_j \\
        d_i s_j &= s_{j} d_{i-1} \ \  \text{if } i>j+1 \\
        s_i s_j &= s_{j+1} s_{i} \ \  \text{if } i<j+1

    INPUT:

    - ``data`` -- a dictionary with entries of the form ``simplex:
      tuple``, where ``simplex`` is a non-degenerate simplex of
      dimension `d` and ``tuple`` has length `d+1`. The `i`-th entry of
      ``tuple`` is the `i`-th face of the simplex. This face is
      recorded as a simplex (degenerate or not) of dimension `d-1`.

    - Alternatively, ``data`` may be a simplicial complex, a
      `\Delta`-complex, or a simplicial set. In the first two cases, it
      is converted to a simplicial set: the simplices in the original
      complex become the non-degenerate simplices in the simplicial
      set. If ``data`` is a simplicial set, return a copy.

    - ``base_point`` (optional) -- one may specify a base point (a
      0-simplex). This is used for some constructions, e.g., the
      wedge. If the simplicial set is reduced (has only one
      0-simplex), then that point is used as the base point.

    - ``name`` (optional) -- string, a name for this simplicial set.

    - ``check`` (optional, default ``True``) -- boolean. If ``True``,
      check the simplicial identity `d_i d_j = d_{j-1} d_i`.

    A note on degeneracies: simplicial sets in Sage are defined by
    specifying the non-degenerate simplices and allowing the
    degeneracies to act freely, subject only to the simplicial
    identities. This is valid for any simplicial set: see (8.3) in
    [EZ50]_, for example.

    REFERENCES:

    .. [EZ50] \S. Eilenberg and J. A. Zilber, *Semi-simplicial
       complexes and singular homology*, Ann. Math. 51 (1950),
       499-513.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: v = AbstractSimplex(0)
        sage: e = AbstractSimplex(1)
        sage: S1 = SimplicialSet({e: (v, v)})  # the circle
        sage: S1
        Simplicial set with 2 non-degenerate simplices

    Simplicial sets may be given names::

        sage: SimplicialSet({e: (v, v)}, name='Simplicial circle')
        Simplicial circle

        sage: s_0_v = v.apply_degeneracies(0) # s_0 applied to v
        sage: sigma = AbstractSimplex(2)
        sage: S2 = SimplicialSet({sigma: (s_0_v, s_0_v, s_0_v)})  # the 2-sphere
        sage: S2
        Simplicial set with 2 non-degenerate simplices
        sage: S2.rename('S^2')
        sage: S2
        S^2

    You don't have to specify the faces of a vertex, but you may,
    using ``None``. This is necessary if you wish to specify a vertex
    which is not contained in the face of any non-degenerate simplex.
    For example, to form a 0-sphere::

        sage: v = AbstractSimplex(0)
        sage: w = AbstractSimplex(0)
        sage: S0 = SimplicialSet({v: None, w: None})
        sage: S0.nondegenerate_simplices()
        [Non-degenerate simplex of dimension 0, Non-degenerate simplex of dimension 0]

    If you want more precise information about the simplices, you can
    name them::

        sage: v.rename('v')
        sage: w.rename('w')
        sage: sorted(S0.nondegenerate_simplices())
        [v, w]

    You can specify a base point::

        sage: S0.is_pointed()
        False
        sage: S0 = S0.set_base_point(v)
        sage: S0.is_pointed()
        True
        sage: S0.base_point()
        v

    You can also specify the base point when you construct the simplicial set::

        sage: S0 = SimplicialSet({v: None, w: None}, base_point=w)
        sage: S0.base_point()
        w

    Converting simplicial complexes and `\Delta`-complexes to
    simplicial sets::

        sage: RP3 = simplicial_complexes.RealProjectiveSpace(3)
        sage: X = SimplicialSet(RP3)
        sage: X
        Simplicial set with 182 non-degenerate simplices

    Simplicial complexes include the "empty simplex" in dimension
    `-1`, but simplicial sets do not::

        sage: RP3.f_vector()
        [1, 11, 51, 80, 40]
        sage: X.f_vector()
        [11, 51, 80, 40]
        sage: RP3.homology(base_ring=GF(2)) == X.homology(base_ring=GF(2))
        True

    The non-degenerate simplices in a converted simplicial complex
    automatically have appropriate names::

        sage: Y = simplicial_complexes.Sphere(1)
        sage: Y.n_cells(1)
        [(0, 1), (0, 2), (1, 2)]
        sage: SimplicialSet(Y).n_cells(1)
        [(0, 1), (0, 2), (1, 2)]

    `\Delta`-complexes: as with simplicial complexes,
    `\Delta`-complexes include the empty simplex in dimension `-1`,
    while simplicial sets do not::

        sage: K = delta_complexes.KleinBottle()
        sage: K
        Delta complex with 1 vertex and 7 simplices
        sage: SimplicialSet(K)
        Simplicial set with 6 non-degenerate simplices

    In a converted `\Delta`-complex, the non-degenerate `d`-simplices
    have names ``Delta_{d,0}``, ``Delta_{d,1}``, ``Delta_{d,2}``, ...::

        sage: SimplicialSet(K).n_cells(0)
        [Delta_{0,0}]
        sage: SimplicialSet(K).n_cells(1)
        [Delta_{1,0}, Delta_{1,1}, Delta_{1,2}]
        sage: SimplicialSet(K).n_cells(2)
        [Delta_{2,0}, Delta_{2,1}]
    """
    def __init__(self, data, base_point=None, name=None, check=True, category=None):
        r"""
        See :class:`SimplicialSet` for documentation.

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: SimplicialSet({e: (v, v, v)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1
            sage: SimplicialSet({e: (v,)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1

        Base points::

            sage: SimplicialSet({e: (v,v)}, base_point=AbstractSimplex(0))
            Traceback (most recent call last):
            ...
            ValueError: the base point is not a simplex in this simplicial set
            sage: SimplicialSet({e: (v,v)}, base_point=e)
            Traceback (most recent call last):
            ...
            ValueError: the base "point" is not a zero-simplex

        Simplicial identity::

            sage: sigma = AbstractSimplex(2)
            sage: w = AbstractSimplex(0)
            sage: K = SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0))})
            sage: SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                        v.apply_degeneracies(0),
            ....:                        w.apply_degeneracies(0))})
            Traceback (most recent call last):
            ...
            ValueError: simplicial identity d_i d_j = d_{j-1} d_i fails in dimension 2

        Returning a copy of the orignal::

            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: SimplicialSet(S1) == S1
            False
        """
        def face(sigma, i):
            """
            i-th face of sigma, a simplex in this simplicial set.

            Once the simplicial set has been fully initialized, you
            can use the :meth:`face` method instead.
            """
            if sigma.is_nondegenerate():
                return data[sigma][i]
            else:
                underlying = sigma.nondegenerate()
                J, t = face_degeneracies(i, sigma.degeneracies())
                if t is None:
                    return underlying.apply_degeneracies(*J)
                else:
                    return data[underlying][t].apply_degeneracies(*J)

        if isinstance(data, GenericCellComplex):
            # Construct new data appropriately.
            if isinstance(data, SimplicialComplex):
                simplices = {}
                faces = {}
                for d in range(data.dimension()+1):
                    old_faces = faces
                    faces = {}
                    for idx, sigma in enumerate(data.n_faces(d)):
                        new_sigma = AbstractSimplex(d)
                        new_sigma.rename(str(sigma))
                        if d > 0:
                            simplices[new_sigma] = [old_faces[_] for _ in sigma.faces()]
                        else:
                            simplices[new_sigma] = None
                        faces[sigma] = new_sigma
                data = simplices

            elif isinstance(data, DeltaComplex):
                simplices = {}
                current = []
                for d in range(data.dimension()+1):
                    faces = tuple(current)
                    current = []
                    for idx, sigma in enumerate(data.n_cells(d)):
                        new_sigma = AbstractSimplex(d)
                        # Name: Delta_{d,idx} where d is dimension,
                        # idx is its index in the list of d-simplices.
                        new_sigma.rename('Delta_{{{},{}}}'.format(d, idx))
                        if d > 0:
                            simplices[new_sigma] = [faces[_] for _ in sigma]
                        else:
                            simplices[new_sigma] = None
                        current.append(new_sigma)
                data = simplices
            elif isinstance(data, SimplicialSet):
                data = dict(copy.deepcopy(data._data))
            else:
                raise NotImplementedError('I do not know how to convert this '
                                          'to a simplicial set')
        # Convert each value in data to a tuple, and then convert all
        # of data to a tuple, so that it is hashable.
        for x in data:
            if data[x]:
                if x.dimension() != len(data[x]) - 1:
                    raise ValueError('wrong number of faces for simplex '
                                     'in dimension {}'.format(x.dimension()))
                if not all(y.dimension() == x.dimension() - 1 for y in data[x]):
                    raise ValueError('faces of a {}-simplex have the wrong '
                                     'dimension'.format(x.dimension()))
                data[x] = tuple(data[x])

        # To obtain the non-degenerate simplices, look at both the
        # keys for data and also the underlying non-degenerate
        # simplices in its values.
        simplices = set(data.keys())
        for t in data.values():
            if t:
                simplices.update([_.nondegenerate() for _ in t])

        for x in simplices:
            if x not in data:
                # x had better be a vertex.
                assert(x.dimension() == 0)
                data[x] = None

        # Check the simplicial identity d_i d_j = d_{j-1} d_i.
        if check:
            for sigma in simplices:
                d = sigma.dimension()
                if d >= 2:
                    for j in range(d+1):
                        for i in range(j):
                            if face(face(sigma, j), i) != face(face(sigma, i), j-1):
                                raise ValueError('simplicial identity d_i d_j '
                                                 '= d_{{j-1}} d_i fails '
                                                 'in dimension {}'.format(d))

        # Now define the attributes for an instance of this class.
        # self._data: a tuple representing the defining data of the
        # simplicial set.
        self._data = tuple(data.items())
        # self._simplices: a sorted tuple of non-degenerate simplices.
        self._simplices = sorted(tuple(simplices))
        # self._basepoint: the base point, or None.
        if base_point is not None:
            if base_point not in simplices:
                raise ValueError('the base point is not a simplex in '
                                 'this simplicial set')
            if base_point.dimension() != 0:
                raise ValueError('the base "point" is not a zero-simplex')
            self._basepoint = base_point
        if category is None:
            if base_point is None:
                category = category=SimplicialSets().Finite()
            else:
                category = category=SimplicialSets().Finite().Pointed()
        Parent.__init__(self, category=category)
        if name:
            self.rename(name)
        # Finally, store the faces of each nondegenerate simplex as an
        # attribute of the simplex itself.
        for x in simplices:
            x._faces[self] = data[x]

    def __eq__(self, other):
        """
        True if ``self`` and ``other`` are equal as simplicial sets.

        Two simplicial sets are equal if they have the same defining
        data.  This means that they have *the same* simplices in each
        dimension, not just that they have the same numbers of
        `n`-simplices for each `n` with the face maps.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({e: (w, w)})
            sage: X == X
            True
            sage: X == SimplicialSet({e: (v, v)})
            True
            sage: X == Y
            False
        """
        if self.is_pointed():
            return (isinstance(other, SimplicialSet)
                    and other.is_pointed()
                    and sorted(self._data) == sorted(other._data)
                    and self.base_point() == other.base_point())
        else:
            return (isinstance(other, SimplicialSet)
                    and not other.is_pointed()
                    and sorted(self._data) == sorted(other._data))

    def __ne__(self, other):
        """
        True if ``self`` and ``other`` are not equal as simplicial sets.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None})
            sage: Y = SimplicialSet({w: None})
            sage: X != X
            False
            sage: X != SimplicialSet({v: None})
            False
            sage: X != Y
            True
        """
        return not (self == other)

    # This is cached because it is used frequently in simplicial set
    # construction: the last two lines in the __init__ method access
    # dictionaries which use instances of SimplicialSet as keys and so
    # computes their hash. If the tuple self._data is long, this can
    # take a long time.
    @cached_method
    def __hash__(self):
        """
        The hash is formed from that of the tuple ``self._data``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: X = SimplicialSet({v: None})
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2)
            sage: Y = SimplicialSet({tau: (degen, degen, degen)})

            sage: hash(X) # random
            17
            sage: hash(X) != hash(Y)
            True
        """
        if self.is_pointed():
            return hash(self._data) ^ hash(self.base_point())
        else:
            return hash(self._data)

    def face_data(self):
        """
        The face-map data defining this simplicial set.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, w)})
            sage: X.face_data()[e]
            (v, w)

            sage: Y = SimplicialSet({v: None, w: None})
            sage: v in Y.face_data()
            True
            sage: Y.face_data()[v] is None
            True
        """
        return dict(self._data)

    # This is cached because it is used frequently in morphism
    # construction when verifying that the morphism commutes with the
    # face maps.
    @cached_method
    def faces(self, simplex):
        """
        The list of faces of ``simplex`` in this simplicial set.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial set, either
          degenerate or not

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: S2.faces(sigma)
            (Simplex obtained by applying degeneracy s_0 to v_0,
             Simplex obtained by applying degeneracy s_0 to v_0,
             Simplex obtained by applying degeneracy s_0 to v_0)
            sage: S2.faces(sigma.apply_degeneracies(0))
            [sigma_2,
             sigma_2,
             Simplex obtained by applying degeneracies s_1 s_0 to v_0,
             Simplex obtained by applying degeneracies s_1 s_0 to v_0]

        TESTS::

            sage: v_0 = S2.n_cells(0)[0]
            sage: S2.faces(v_0) is None
            True

            sage: from sage.homology.simplicial_set import AbstractSimplex
            sage: w = AbstractSimplex(0)
            sage: S2.faces(w)
            Traceback (most recent call last):
            ...
            ValueError: this simplex is not in this simplicial set
        """
        if simplex.nondegenerate() not in self._simplices:
            raise ValueError('this simplex is not in this simplicial set')
        if simplex in self._simplices:
            return simplex._faces[self]
        underlying = simplex.nondegenerate()
        faces = []
        for J, t in [face_degeneracies(m, simplex.degeneracies())
                     for m in range(simplex.dimension()+1)]:
            if t is None:
                faces.append(underlying.apply_degeneracies(*J))
            else:
                faces.append(self.face(underlying, t).apply_degeneracies(*J))
        return faces

    def face(self, simplex, i):
        """
        The `i`-th face of ``simplex`` in this simplicial set.

        INPUT:

        - ``simplex`` -- a simplex in this simplicial set
        - ``i`` -- integer

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: v_0 = S2.n_cells(0)[0]
            sage: S2.face(sigma, 0)
            Simplex obtained by applying degeneracy s_0 to v_0
            sage: S2.face(sigma, 0) == v_0.apply_degeneracies(0)
            True
            sage: S2.face(S2.face(sigma, 0), 0) == v_0
            True
        """
        if i < 0 or i > simplex.dimension():
            raise ValueError('cannot compute face {} of {}-dimensional '
                             'simplex'.format(i, simplex.dimension()))
        faces = self.faces(simplex)
        if faces is not None:
            return self.faces(simplex)[i]
        return None

    def __contains__(self, x):
        """
        True if ``x`` is a simplex which is contained in this complex.

        EXAMPLES::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: v0 = S0.n_cells(0)[0]
            sage: v0 in S0
            True
            sage: v0 in S1
            False

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: e = AbstractSimplex(1)
            sage: K = SimplicialSet({e: (v, v)})  # the circle
            sage: v in K
            True
            sage: v0 in K
            False
            sage: S1.n_cells(1)[0] in K
            False
        """
        return x.nondegenerate() in self._simplices

    def alexander_whitney(self, simplex, dim_left):
        r"""
        'Subdivide' ``simplex`` in this simplicial set into a pair of
        simplices.

        The left factor should have dimension ``dim_left``, so the
        right factor should have dimension ``dim - dim_left``, if
        ``dim`` is the dimension of the starting simplex. The results
        are obtained by applying iterated face maps to
        ``simplex``. Writing `d` for ``dim`` and `j` for ``dim_left``:
        apply `d_{j+1} d_{j+2} ... d_{d}` to get the left factor,
        `d_0 ... d_0` to get the right factor.

        INPUT:

        - ``dim_left`` -- integer, the dimension of the left-hand factor

        OUTPUT: a list containing the triple ``(c, left, right)``,
        where ``left`` and ``right`` are the two simplices described
        above. If either ``left`` or ``right`` is degenerate, ``c`` is
        0; otherwise, ``c`` is 1. This is so that, when used to
        compute cup products, it is easy to ignore terms which have
        degenerate factors.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: sigma = S2.n_cells(2)[0]
            sage: S2.alexander_whitney(sigma, 0)
            [(1, v_0, sigma_2)]
            sage: S2.alexander_whitney(sigma, 1)
            [(0,
             Simplex obtained by applying degeneracy s_0 to v_0,
             Simplex obtained by applying degeneracy s_0 to v_0)]
        """
        dim = simplex.dimension()
        if dim_left < 0 or dim_left > dim:
            raise ValueError('alexander_whitney is only valid if dim_left '
                             'is between 0 and the dimension of the simplex')
        left = simplex
        for i in range(dim, dim_left, -1):
            left = self.face(left, i)
        right = simplex
        for i in range(dim_left):
            right = self.face(right, 0)
        if left.is_degenerate() or right.is_degenerate():
            c = ZZ.zero()
        else:
            c = ZZ.one()
        return [(c, left, right)]

    def nondegenerate_simplices(self):
        """
        The sorted list of non-degenerate simplices in this simplicial set.

        The sorting is in increasing order of dimension, and within
        each dimension, by the name (if present) of each simplex.

        .. NOTE::

            The sorting is done when the simplicial set is
            constructed, so changing the name of a simplex after
            construction will not affect the ordering.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.nondegenerate_simplices()
            [Non-degenerate simplex of dimension 0, Non-degenerate simplex of dimension 0]

        Name the vertices and reconstruct the simplicial set: they
        should be ordered alphabetically::

            sage: v.rename('v')
            sage: w.rename('w')
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.nondegenerate_simplices()
            [v, w]

        Rename but do not reconstruct the set; the ordering does not
        take the new names into account::

            sage: v.rename('z')
            sage: S0.nondegenerate_simplices() # old ordering is used
            [z, w]

            sage: X0 = SimplicialSet({v: None, w: None})
            sage: X0.nondegenerate_simplices() # new ordering is used
            [w, z]
        """
        return list(self._simplices)

    def cells(self, subcomplex=None):
        """
        Dictionary of all non-degenerate simplices.

        INPUT:

        - ``subcomplex`` (optional) -- a subsimplicial set of this
          simplicial set

        Each key is a dimension, and the corresponding value is the
        list of simplices in that dimension. If ``subcomplex`` is
        specified, then the corresponding value is the list of
        simplices in the quotient by the subcomplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.cells()
            {0: [Non-degenerate simplex of dimension 0, Non-degenerate simplex of dimension 0]}

            sage: v.rename('v')
            sage: w.rename('w')
            sage: S0.cells()
            {0: [v, w]}

            sage: e = AbstractSimplex(1, name='e')
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: S1.cells()
            {0: [v], 1: [e]}

            sage: S0.cells(S0.subsimplicial_set([v, w]))
            {0: [*]}

            sage: X = SimplicialSet({e: (v,w)})
            sage: X.cells(X.subsimplicial_set([v, w]))
            {0: [*], 1: [e]}
        """
        if subcomplex is None:
            simplices = {}
            for sigma in self.nondegenerate_simplices():
                if sigma.dimension() in simplices:
                    simplices[sigma.dimension()].append(sigma)
                else:
                    simplices[sigma.dimension()] = [sigma]
            return {d: sorted(simplices[d]) for d in simplices}
        return self.quotient(subcomplex).cells()

    def f_vector(self):
        """
        Return the list of the number of non-degenerate simplices in each
        dimension.

        Unlike for some other cell complexes in Sage, this does not
        include the empty simplex in dimension `-1`; thus its `i`-th
        entry is the number of `i`-dimensional simplices.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.f_vector()
            [2]

            sage: e = AbstractSimplex(1)
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: S1.f_vector()
            [1, 1]
            sage: simplicial_sets.Sphere(3).f_vector()
            [1, 0, 0, 1]
        """
        return [len(self.n_cells(_)) for _ in range(self.dimension()+1)]

    def euler_characteristic(self):
        r"""
        The Euler characteristic of this simplicial set: the
        alternating sum over `n \geq 0` of the number of
        nondegenerate `n`-simplices.

        EXAMPLES::

            sage: simplicial_sets.RealProjectiveSpace(4).euler_characteristic()
            1
            sage: simplicial_sets.Sphere(6).euler_characteristic()
            2
            sage: simplicial_sets.KleinBottle().euler_characteristic()
            0
        """
        return sum([(-1)**n * num for (n, num) in enumerate(self.f_vector())])

    def all_simplices(self, dim):
        """
        All simplices, non-degenerate and degenerate, in dimension ``dim``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2, name='tau')
            sage: Y = SimplicialSet({tau: (degen, degen, degen), w: None})

        ``Y`` is the disjoint union of a 2-sphere, with vertex ``v``
        and non-degenerate 2-simplex ``tau``, and a point ``w``. ::

            sage: Y.all_simplices(0)
            [v, w]
            sage: Y.all_simplices(1)
            [Simplex obtained by applying degeneracy s_0 to v,
             Simplex obtained by applying degeneracy s_0 to w]
            sage: Y.all_simplices(2)
            [tau,
             Simplex obtained by applying degeneracies s_1 s_0 to v,
             Simplex obtained by applying degeneracies s_1 s_0 to w]
        """
        non_degen = [_ for _ in self.nondegenerate_simplices() if _.dimension() <= dim]
        ans = set([_ for _ in non_degen if _.dimension() == dim])
        for sigma in non_degen:
            d = sigma.dimension()
            ans.update([sigma.apply_degeneracies(*_)
                        for _ in all_degeneracies(d, dim - d)])
        return sorted(list(ans))

    def _map_from_empty_set(self):
        """
        The unique map from the empty set to this simplicial set.

        This is used to in the method :meth:`disjoint_union` to
        construct disjoint unions as pushouts.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: T._map_from_empty_set()
            Simplicial set morphism:
              From: Empty simplicial set
              To:   Torus
              Defn: [] --> []
        """
        return Empty().Hom(self)({})

    def constant_map(self, codomain=None):
        """
        Return a constant map with this space as its domain.

        INPUT:

        - ``codomain`` -- optional, default ``None``. If ``None``, the
          codomain is the standard one-point space constructed by
          :func:`Point`. Otherwise, the codomain must be a pointed
          simplicial set, in which case the map is constant at the
          base point.

        EXAMPLES::

            sage: S4 = simplicial_sets.Sphere(4)
            sage: S4.constant_map()
            Simplicial set morphism:
              From: S^4
              To:   Point
              Defn: [v_0, sigma_4] --> [*, Simplex obtained by applying degeneracies s_3 s_2 s_1 s_0 to *]
            sage: S0 = simplicial_sets.Sphere(0)
            sage: S4.constant_map(codomain=S0)
            Simplicial set morphism:
              From: S^4
              To:   S^0
              Defn: [v_0, sigma_4] --> [v_0, Simplex obtained by applying degeneracies s_3 s_2 s_1 s_0 to v_0]

        TESTS::

            sage: S0 = S0.unset_base_point()
            sage: S4.constant_map(codomain=S0)
            Traceback (most recent call last):
            ...
            ValueError: the codomain must have a base point
        """
        if codomain is not None:
            if not codomain.is_pointed():
                raise ValueError('the codomain must have a base point')
        else:
            codomain = Point()
        base_pt = codomain.base_point()
        d = {}
        cells = self.cells()
        for dim in cells:
            d.update({sigma:base_pt.apply_degeneracies(*range(dim-1, -1, -1))
                      for sigma in cells[dim]})
        return self.Hom(codomain)(d)

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this simplicial set.

        That is, the subsimplicial set generated by all nondegenerate
        simplices of dimension at most `n`.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2, name='tau')
            sage: Y = SimplicialSet({tau: (degen, degen, degen), w: None})

        ``Y`` is the disjoint union of a 2-sphere, with vertex ``v``
        and non-degenerate 2-simplex ``tau``, and a point ``w``. ::

            sage: Y.nondegenerate_simplices()
            [v, w, tau]
            sage: Y.n_skeleton(1).nondegenerate_simplices()
            [v, w]
            sage: Y.n_skeleton(2).nondegenerate_simplices()
            [v, w, tau]
        """
        data = [x for x in self.nondegenerate_simplices()
                if x.dimension() <= n]
        return self.subsimplicial_set(data)

    def chain_complex(self, dimensions=None, base_ring=ZZ, augmented=False,
                      cochain=False, verbose=False, subcomplex=None,
                      check=False):
        r"""
        The normalized chain complex.

        INPUT:

        - ``dimensions`` -- if ``None``, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.

        - ``base_ring`` (optional, default ``ZZ``) -- commutative ring

        - ``augmented`` (optional, default ``False``) -- if ``True``,
          return the augmented chain complex (that is, include a class
          in dimension `-1` corresponding to the empty cell).

        - ``cochain`` (optional, default ``False``) -- if ``True``,
          return the cochain complex (that is, the dual of the chain
          complex).

        - ``verbose`` (optional, default ``False``) -- ignored.

        - ``subcomplex`` (optional, default ``None``) -- if present,
          compute the chain complex relative to this subcomplex.

        - ``check`` (optional, default ``False``) -- If ``True``, make
           sure that the chain complex is actually a chain complex:
           the differentials are composable and their product is zero.

        The normalized chain complex of a simplicial set is isomorphic
        to the chain complex obtained by modding out by degenerate
        simplices, and the latter is what is actually constructed
        here.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: degen = v.apply_degeneracies(1, 0) # s_1 s_0 applied to v
            sage: sigma = AbstractSimplex(3)
            sage: S3 = SimplicialSet({sigma: (degen, degen, degen, degen)}) # the 3-sphere
            sage: S3.chain_complex().homology()
            {0: Z, 3: Z}
            sage: S3.chain_complex(augmented=True).homology()
            {-1: 0, 0: 0, 3: Z}
            sage: S3.chain_complex(dimensions=range(3), base_ring=QQ).homology()
            {0: Vector space of dimension 1 over Rational Field}

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5.chain_complex(subcomplex=RP2).homology()
            {0: Z, 3: C2, 4: 0, 5: Z}

        TESTS:

        Convert some simplicial complexes and `\Delta`-complexes to
        simplicial sets, and compare homology calculations::

            sage: T = simplicial_complexes.Torus()
            sage: T.homology() == SimplicialSet(T).homology()
            True
            sage: RP2 = delta_complexes.RealProjectivePlane()
            sage: RP2.homology() == SimplicialSet(RP2).homology()
            True
            sage: RP2.cohomology(base_ring=GF(2)) == SimplicialSet(RP2).cohomology(base_ring=GF(2))
            True
        """
        if dimensions is None:
            if not self.cells(): # Empty
                if cochain:
                    return ChainComplex({-1: matrix(base_ring, 0, 0)},
                                        degree_of_differential=1)
                return ChainComplex({0: matrix(base_ring, 0, 0)},
                                    degree_of_differential=-1)
            dimensions = range(self.dimension()+1)
        else:
            if not isinstance(dimensions, (list, tuple)):
                dimensions = range(dimensions-1, dimensions+2)
            augmented = False

        differentials = {}
        # Convert the tuple self._data to a dictionary indexed by the
        # non-degenerate simplices.
        if subcomplex:
            X = self.quotient(subcomplex)
            face_data = X.face_data()
            nondegens = X.nondegenerate_simplices()
        else:
            face_data = self.face_data()
            nondegens = self.nondegenerate_simplices()
        # simplices: dictionary indexed by dimension, values the list
        # of non-degenerate simplices in that dimension.
        simplices = {}
        for sigma in nondegens:
            if sigma.dimension() in simplices:
                simplices[sigma.dimension()].append(sigma)
            else:
                simplices[sigma.dimension()] = [sigma]
        first = dimensions.pop(0)
        if first in simplices:
            rank = len(simplices[first])
            current = sorted(simplices[first])
        else:
            rank = 0
            current = []
        if augmented:
            differentials[first-1] = matrix(base_ring, 0, 1)
            differentials[first] = matrix(base_ring, 1, rank,
                                      [1] * rank)
        else:
            differentials[first] = matrix(base_ring, 0, rank)

        for d in dimensions:
            old_rank = rank
            faces = {_[1]:_[0] for _ in enumerate(current)}
            if d in simplices:
                current = sorted(simplices[d])
                rank = len(current)
                # old_rank: number of simplices in dimension d-1.
                # faces: list of simplices in dimension d-1.
                # rank: number of simplices in dimension d.
                # current: list of simplices in dimension d.
                if not faces:
                    differentials[d] = matrix(base_ring, old_rank, rank)
                else:
                    matrix_data = {}
                    for col, sigma in enumerate(current):
                        sign = 1
                        for tau in face_data[sigma]:
                            if tau.is_nondegenerate():
                                row = faces[tau]
                                if (row, col) in matrix_data:
                                    matrix_data[(row, col)] += sign
                                else:
                                    matrix_data[(row, col)] = sign
                            sign *= -1

                    differentials[d] = matrix(base_ring, old_rank,
                                              rank, matrix_data)

            else:
                rank = 0
                current = []
                differentials[d] = matrix(base_ring, old_rank, rank)

        if cochain:
            new_diffs = {}
            for d in differentials:
                new_diffs[d-1] = differentials[d].transpose()
            return ChainComplex(new_diffs, degree_of_differential=1,
                                check=check)
        return ChainComplex(differentials, degree_of_differential=-1,
                            check=check)

    @cached_method
    def algebraic_topological_model(self, base_ring=None):
        r"""
        Algebraic topological model for this simplicial set with
        coefficients in ``base_ring``.

        The term "algebraic topological model" is defined by Pilarczyk
        and Réal [PR]_.

        INPUT:

        - ``base_ring`` - coefficient ring (optional, default
          ``QQ``). Must be a field.

        Denote by `C` the chain complex associated to this simplicial
        set. The algebraic topological model is a chain complex `M`
        with zero differential, with the same homology as `C`, along
        with chain maps `\pi: C \to M` and `\iota: M \to C` satisfying
        `\iota \pi = 1_M` and `\pi \iota` chain homotopic to
        `1_C`. The chain homotopy `\phi` must satisfy

        - `\phi \phi = 0`,
        - `\pi \phi = 0`,
        - `\phi \iota = 0`.

        Such a chain homotopy is called a *chain contraction*.

        OUTPUT: a pair consisting of

        - chain contraction ``phi`` associated to `C`, `M`, `\pi`, and
          `\iota`
        - the chain complex `M`

        Note that from the chain contraction ``phi``, one can recover the
        chain maps `\pi` and `\iota` via ``phi.pi()`` and
        ``phi.iota()``. Then one can recover `C` and `M` from, for
        example, ``phi.pi().domain()`` and ``phi.pi().codomain()``,
        respectively.

        EXAMPLES::

            sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
            sage: phi, M = RP2.algebraic_topological_model(GF(2))
            sage: M.homology()
            {0: Vector space of dimension 1 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}
            sage: T = simplicial_sets.Torus()
            sage: phi, M = T.algebraic_topological_model(QQ)
            sage: M.homology()
            {0: Vector space of dimension 1 over Rational Field,
             1: Vector space of dimension 2 over Rational Field,
             2: Vector space of dimension 1 over Rational Field}
        """
        from algebraic_topological_model import algebraic_topological_model_delta_complex
        if base_ring is None:
            base_ring = QQ
        return algebraic_topological_model_delta_complex(self, base_ring)

    def __copy__(self):
        """
        Return a copy of this simplicial set.

        All simplices in the new version are copies of, not identical
        to, the old simplices.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v,v)})
            sage: Y = copy(X)
            sage: Y
            Simplicial set with 2 non-degenerate simplices

        Because the simplices in ``Y`` are copies of those in ``X``,
        they have been renamed::

            sage: Y.nondegenerate_simplices()
            [v', e']

            sage: Y.is_pointed()
            False
            sage: Y = X.set_base_point(v)
            sage: Y.is_pointed()
            True
        """
        translate = {x: copy.copy(x) for x in self.nondegenerate_simplices()}
        old_data = self.face_data()
        new_data = {}
        for x in old_data:
            faces = old_data[x]
            if faces is None: # isolated vertex
                new_data[translate[x]] = None
            else:
                new_faces = []
                for y in faces:
                    if y.is_nondegenerate():
                        new_faces.append(translate[y])
                    else:
                        z = y.nondegenerate()
                        new_faces.append(translate[z].apply_degeneracies(*y.degeneracies()))
                new_data[translate[x]] = tuple(new_faces)

        if self.is_pointed():
            base_point=translate[self.base_point()]
            C = SimplicialSet(new_data, base_point=base_point)
        else:
            C = SimplicialSet(new_data)
        return C

    def subsimplicial_set(self, simplices):
        """
        The sub-simplicial set of this simplicial set determined by
        ``simplices``, a set of nondegenerate simplices.

        INPUT:

        - ``simplices`` -- set, list, or tuple of nondegenerate
          simplices in this simplicial set, or a simplicial
          complex -- see below.

        If ``simplices`` is a simplicial complex, then the original
        simplicial set should itself have been converted from a
        simplicial complex, and ``simplices`` should be a subcomplex
        of that.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')

            sage: X = SimplicialSet({e: (v, w), f: (w, v)})
            sage: Y = X.subsimplicial_set([e])
            sage: Y
            Simplicial set with 3 non-degenerate simplices
            sage: Y.nondegenerate_simplices()
            [v, w, e]

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: K = SimplicialSet(S3)
            sage: tau = K.n_cells(3)[0]
            sage: tau.dimension()
            3
            sage: K.subsimplicial_set([tau])
            Simplicial set with 15 non-degenerate simplices

        TESTS:

        Make sure vertices are treated properly::

            sage: X.subsimplicial_set([v]).nondegenerate_simplices()
            [v]
            sage: X.subsimplicial_set([v, w]).nondegenerate_simplices()
            [v, w]
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.subsimplicial_set([w]).nondegenerate_simplices()
            [w]

        Raise an error if an element of ``simplices`` is not actually
        in the original simplicial set::

            sage: sigma = AbstractSimplex(2, name='sigma_2')
            sage: Z = X.subsimplicial_set([e, sigma])
            Traceback (most recent call last):
            ...
            ValueError: not all simplices are in the original simplicial set

        Simplicial complexes::

            sage: X = simplicial_complexes.ComplexProjectivePlane()
            sage: Y = X._contractible_subcomplex()
            sage: CP2 = SimplicialSet(X)
            sage: sub = CP2.subsimplicial_set(Y)
            sage: CP2.f_vector()
            [9, 36, 84, 90, 36]
            sage: K = CP2.quotient(sub)
            sage: K.f_vector()
            [1, 0, 23, 45, 24]
            sage: K.homology()
            {0: 0, 1: 0, 2: Z, 3: 0, 4: Z}

        Try to construct a subcomplex from a simplicial complex which
        is not actually contained in ``self``::

            sage: Z = SimplicialComplex([[0,1,2,3,4]])
            sage: CP2.subsimplicial_set(Z)
            Traceback (most recent call last):
            ...
            ValueError: not all simplices are in the original simplicial set
        """
        # If simplices is a simplicial complex, turn it into a list of
        # nondegenerate simplices.
        if isinstance(simplices, SimplicialComplex):
            new = []
            for f in simplices.facets():
                d = f.dimension()
                found = False
                for x in self.n_cells(d):
                    if str(x) == str(f):
                        new.append(x)
                        found = True
                        break
                if not found:
                    raise ValueError('not all simplices are in the original simplicial set')
            simplices = new

        data = self.face_data()
        vertices = set([])
        keep = set(simplices)
        old_keep = set([])
        while keep != old_keep:
            old_keep = copy.copy(keep)
            for x in old_keep:
                underlying = x.nondegenerate()
                if underlying not in self._simplices:
                    raise ValueError('not all simplices are in the original simplicial set')
                keep.add(underlying)
                if underlying in data and data[underlying]:
                    keep.update([f.nondegenerate() for f in data[underlying]])
                else:
                    # x is a vertex
                    assert(underlying.dimension() == 0)
                    vertices.add(underlying)
        missing = set(self.nondegenerate_simplices()).difference(keep)
        for x in missing:
            if x in data:
                del data[x]
        for x in vertices:
            data[x] = None
        return SubSimplicialSet(data, self)

    def is_subset(self, other):
        """
        Return ``True`` iff ``self`` is a subsimplicial set of ``other``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: S1 = simplicial_sets.Sphere(1)
            sage: pt = S1.n_cells(0)[0]
            sage: other_pt = AbstractSimplex(0)

            sage: SimplicialSet({pt: None}).is_subset(S1)
            True
            sage: SimplicialSet({other_pt: None}).is_subset(S1)
            False
        """
        return all(x in other._simplices for x in self.nondegenerate_simplices())

    def quotient(self, subcomplex, vertex_name='*'):
        """
        The quotient of this simplicial set by ``subcomplex``.

        That is, ``subcomplex`` is replaced by a vertex.

        INPUT:

        - ``subcomplex`` -- subsimplicial set of this simplicial set,
          or a list, tuple, or set of simplices defining a
          subsimplicial set.

        - ``vertex_name`` (optional) -- string, name to be given to the new
          vertex. By default, use ``'*'``.

        Base points: if the original simplicial set has a base point
        contained in ``subcomplex``, then ``*`` is the base point in
        the quotient. If the original simplicial set has a base point
        not contained in ``subcomplex``, then use that base point. If
        the original simplicial set had no base point, then use
        ``*``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, w), f: (v, w)})
            sage: Y = X.quotient([f])
            sage: Y.nondegenerate_simplices()
            [*, e]
            sage: Y.homology(1)
            Z

            sage: E = SimplicialSet({e: (v, w)})
            sage: Z = E.quotient([v, w])
            sage: Z.nondegenerate_simplices()
            [*, e]
            sage: Z.homology(1)
            Z

            sage: F = E.quotient([v])
            sage: F.nondegenerate_simplices()
            [*, w, e]
            sage: F.base_point()
            *

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 0 over Finite Field of size 2,
             3: Vector space of dimension 1 over Finite Field of size 2,
             4: Vector space of dimension 1 over Finite Field of size 2,
             5: Vector space of dimension 1 over Finite Field of size 2}

        TESTS::

            sage: pt = RP5.quotient(RP5.n_skeleton(5))
            sage: pt
            Simplicial set with 1 non-degenerate simplex
        """
        if not isinstance(subcomplex, SimplicialSet):
            # If it's not a simplicial set, subcomplex should be a
            # list, tuple, or set of simplices, so form the actual
            # subcomplex:
            subcomplex = self.subsimplicial_set(subcomplex)
        else:
            # Test whether subcomplex is actually a subcomplex of
            # self.
            if not subcomplex.is_subset(self):
                raise ValueError('the "subcomplex" is not actually a subcomplex')
            pass
        return QuotientOfSimplicialSet(subcomplex, vertex_name=vertex_name)

    def disjoint_union(self, *others):
        """
        The disjoint union of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        The simplicial sets must have no simplices in common.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({f: (v, w)})
            sage: Z = X.disjoint_union(Y)

        Since ``X`` and ``Y`` have simplices in common, Sage uses a
        copy of ``Y`` when constructing the disjoint union. Note the
        name conflict in the list of simplices: ``v`` appears twice::

            sage: Z = X.disjoint_union(Y)
            sage: Z.nondegenerate_simplices()
            [v, v, w, e, f]
        """
        maps = [self._map_from_empty_set()] + [X._map_from_empty_set() for X in others]
        return PushoutOfSimplicialSets(maps)

    def coproduct(self, *others):
        """
        The coproduct of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        If these simplicial sets are pointed, return their wedge sum;
        if they are not, return their disjoint union. If some are
        pointed and some are not, raise an error: it is not clear in
        which category to work.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: D3 = simplicial_sets.Simplex(3)
            sage: Y = S2.unset_base_point()
            sage: Z = K.unset_base_point()

            sage: S2.coproduct(K).is_pointed()
            True
            sage: S2.coproduct(K)
            Simplicial set with 7 non-degenerate simplices
            sage: D3.coproduct(Y, Z).is_pointed()
            False
            sage: D3.coproduct(Y, Z)
            Simplicial set with 23 non-degenerate simplices

            sage: D3.coproduct(S2, Z)
            Traceback (most recent call last):
            ...
            ValueError: some, but not all, of the simplicial sets are pointed, so the categorical coproduct is not defined: the category is ambiguous
        """
        if self.is_pointed() and all(X.is_pointed() for X in others):
            return self.wedge(*others)
        if self.is_pointed() or any(X.is_pointed() for X in others):
            raise ValueError('some, but not all, of the simplicial sets are pointed, '
                             'so the categorical coproduct is not defined: the '
                             'category is ambiguous')
        return self.disjoint_union(*others)

    def product(self, *others):
        r"""
        The product of this simplicial set with ``others``.

        INPUT:

        - ``others`` -- one or several simplicial sets

        If `X` and `Y` are simplicial sets, then their product `X
        \times Y` is defined to be the simplicial set with
        `n`-simplices `X_n \times Y_n`. See
        :class:`ProductOfSimplicialSets` for more information.

        If a simplicial set is constructed as a product, the factors
        are recorded and are accessible via the method
        :meth:`~ProductOfSimplicialSets.factors`. If it is constructed
        as a product and then copied, this information is lost.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, w)})
            sage: square = X.product(X)

        ``square`` is now the standard triangulation of the square: 4
        vertices, 5 edges (the four on the border and the diagonal), 2
        triangles::

            sage: square.f_vector()
            [4, 5, 2]

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: T.homology(reduced=False)
            {0: Z, 1: Z x Z, 2: Z}

        Since ``S1`` is pointed, so is ``T``::

            sage: S1.is_pointed()
            True
            sage: S1.base_point()
            v_0
            sage: T.is_pointed()
            True
            sage: T.base_point()
            (v_0, v_0)

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: Z = S2.product(S3)
            sage: Z.homology(reduced=False)
            {0: Z, 1: 0, 2: Z, 3: Z, 4: 0, 5: Z}

            sage: Z.factors() == (S2, S3)
            True
            sage: Z.factors() == (S3, S2)
            False

        More than two factors::

            sage: W = S1.product(S1, S1)
            sage: W.homology(reduced=False)
            {0: Z, 1: Z x Z x Z, 2: Z x Z x Z, 3: Z}
            sage: W.factors()[1] == S1
            True
            sage: W.factors()
            (S^1, S^1, S^1)
        """
        return ProductOfSimplicialSets((self,) + others)

    def pushout(self, *maps):
        r"""
        The pushout obtained from given ``maps``.

        INPUT:

        - ``maps`` -- several maps of simplicial sets, each of which
          has this simplicial set as its domain

        If only a single map `f: X \to Y` is given, then return
        `Y`. If more than one map is given, say `f_i: X \to Y_i` for
        `0 \leq i \leq m`, then return the pushout defined by those
        maps. If no maps are given, return the empty simplicial set.

        EXAMPLES:

        Construct the 4-sphere as a quotient of a 4-simplex::

            sage: K = simplicial_sets.Simplex(4)
            sage: L = K.n_skeleton(3)
            sage: S4 = L.pushout(L.constant_map(), L.inclusion_map())
            sage: S4
            Simplicial set with 2 non-degenerate simplices
            sage: S4.homology(4)
            Z

        TESTS::

            sage: K = simplicial_sets.Simplex(5)
            sage: K.pushout()
            Simplicial set with 0 non-degenerate simplices

            sage: S0 = simplicial_sets.Sphere(0)
            sage: pt_map = S0.base_point_map()
            sage: pt_map.domain().pushout(pt_map) == S0
            True

            sage: K.pushout(K.constant_map(), pt_map)
            Traceback (most recent call last):
            ...
            ValueError: the domains of the maps must be equal
        """
        if any(self != f.domain() for f in maps):
            raise ValueError('the domains of the maps must be equal')
        return PushoutOfSimplicialSets(maps)

    def pullback(self, *maps):
        r"""
        The pullback obtained from given ``maps``.

        INPUT:

        - ``maps`` -- several maps of simplicial sets, each of which
          has this simplicial set as its codomain

        If only a single map `f: X \to Y` is given, then return
        `X`. If more than one map is given, say `f_i: X_i \to Y` for
        `0 \leq i \leq m`, then return the pullback defined by those
        maps. If no maps are given, return the one-point simplicial
        set.

        EXAMPLES:

        Construct a product as a pullback::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: pt = simplicial_sets.Point()
            sage: P = pt.pullback(S2.constant_map(), S2.constant_map())
            sage: P.homology(2)
            Z x Z

        TESTS::

            sage: pt.pullback(S2.constant_map(), S2.base_point_map())
            Traceback (most recent call last):
            ...
            ValueError: the codomains of the maps must be equal
        """
        if any(self != f.codomain() for f in maps):
            raise ValueError('the codomains of the maps must be equal')
        return PullbackOfSimplicialSets(maps)

    # Ideally, this would be defined at the category level and only
    # for pointed simplicial sets, but the abstract_method "wedge" in
    # cell_complex.py shadows that.
    def wedge(self, *others):
        r"""
        The wedge sum of this pointed simplicial set with ``others``.

        - ``others`` -- one or several simplicial sets

        This constructs the quotient of the disjoint union in which
        the base points of all of the simplicial sets have been
        identified. This is the coproduct in the category of pointed
        simplicial sets.

        This raises an error if any of the factors is not pointed.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: W = X.wedge(Y)
            sage: W.nondegenerate_simplices()
            [*, e, f]
            sage: W.homology()
            {0: 0, 1: Z x Z}
            sage: S2 = simplicial_sets.Sphere(2)
            sage: X.wedge(S2).homology(reduced=False)
            {0: Z, 1: Z, 2: Z}
            sage: X.wedge(X).nondegenerate_simplices()
            [*, e, e]

        TESTS::

            sage: Z = SimplicialSet({e: (v,w)})
            sage: X.wedge(Z)
            Traceback (most recent call last):
            ...
            ValueError: the simplicial sets must be pointed
        """
        from sage.homology.simplicial_set import WedgeOfSimplicialSets
        return WedgeOfSimplicialSets((self,) + others)

    def cone(self):
        r"""
        The (reduced) cone on this simplicial set.

        If this simplicial set `X` is not pointed, construct the
        ordinary cone: add a point `v` (which will become the base
        point) and for each simplex `\sigma` in `X`, add both `\sigma`
        and a simplex made up of `v` and `\sigma` (topologically, form
        the join of `v` and `\sigma`).

        If this simplicial set is pointed, then construct the reduced
        cone: take the quotient of the unreduced cone by the 1-simplex
        connecting the old base point to the new one.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone()  # unreduced cone, since X not pointed
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *
            sage: X = X.set_base_point(v)
            sage: CX = X.cone()  # reduced cone
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]
        """
        if self.is_pointed():
            return ReducedConeOfSimplicialSet(self)
        return ConeOfSimplicialSet(self)

    def suspension(self, n=1):
        """
        The (reduced) `n`-th suspension of this simplicial set.

        INPUT:

        - ``n`` (optional, default 1) -- integer, suspend this many
          times.

        If this simplicial set `X` is not pointed, return the
        suspension: the quotient `CX/X`, where `CX` is the (ordinary,
        unreduced) cone on `X`. If `X` is pointed, then use the
        reduced cone instead, and so return the reduced suspension.

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: SigmaRP4 = RP4.suspension()
            sage: S1_smash_RP4 = S1.smash_product(RP4)
            sage: SigmaRP4.homology() == S1_smash_RP4.homology()
            True

        The version of the suspension obtained by the smash product is
        typically less efficient than the reduced suspension produced
        here::

            sage: SigmaRP4
            Simplicial set with 5 non-degenerate simplices
            sage: S1_smash_RP4
            Simplicial set with 25 non-degenerate simplices

        TESTS::

            sage: RP4.suspension(-3)
            Traceback (most recent call last):
            ...
            ValueError: n must be non-negative
        """
        if n == 0:
            return self
        if n < 0:
            raise ValueError('n must be non-negative')
        CX = self.cone()
        if self.is_pointed():
            X = CX.map_from_X().image()
        else:
            X = CX.X_as_subset()
        Sigma = CX.quotient(X)
        return Sigma.suspension(n-1)

    def is_reduced(self):
        """
        Return ``True`` if this simplicial set has only one vertex.

        EXAMPLES::

            sage: simplicial_sets.Sphere(0).is_reduced()
            False
            sage: simplicial_sets.Sphere(3).is_reduced()
            True
        """
        return len(self.n_cells(0)) == 1

    def graph(self):
        """
        The 1-skeleton of this simplicial set, as a graph.

        EXAMPLES::

            sage: Delta3 = simplicial_sets.Simplex(3)
            sage: G = Delta3.graph()
            sage: G.edges()
            [((0,), (1,), (0, 1)),
             ((0,), (2,), (0, 2)),
             ((0,), (3,), (0, 3)),
             ((1,), (2,), (1, 2)),
             ((1,), (3,), (1, 3)),
             ((2,), (3,), (2, 3))]

            sage: T = simplicial_sets.Torus()
            sage: T.graph()
            Looped multi-graph on 1 vertex
            sage: len(T.graph().edges())
            3

            sage: CP3 = simplicial_sets.ComplexProjectiveSpace(3)
            sage: G = CP3.graph()
            sage: len(G.vertices())
            1
            sage: len(G.edges())
            0
        """
        edges = self.n_cells(1)
        vertices = self.n_cells(0)
        used_vertices = set([])  # vertices which are in an edge
        d = {}
        for e in edges:
            v = self.face(e, 0)
            w = self.face(e, 1)
            if v in d:
                if w in d[v]:
                    d[v][w] = d[v][w] + [e]
                else:
                    d[v][w] = [e]
            else:
                d[v] = {w: [e]}
            used_vertices.update([v, w])
        for v in vertices:
            if v not in used_vertices:
                d[v] = {}
        return Graph(d, format='dict_of_dicts')

    def reduce(self):
        """
        Reduce this simplicial set.

        That is, take the quotient by a spanning tree of the
        1-skeleton, so that the resulting simplicial set has only one
        vertex. This only makes sense if the simplicial set is
        connected, so raise an error if not. If already reduced,
        return itself.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(2)
            sage: K.is_reduced()
            False
            sage: X = K.reduce()
            sage: X.is_reduced()
            True

        ``X`` is reduced, so calling ``reduce`` on it again
        returns ``X`` itself::

            sage: X is X.reduce()
            True
            sage: K is K.reduce()
            False

        Raise an error for disconnected simplicial sets::

            sage: S0 = simplicial_sets.Sphere(0)
            sage: S0.reduce()
            Traceback (most recent call last):
            ...
            ValueError: this simplicial set is not connected
        """
        if self.is_reduced():
            return self
        if not self.is_connected():
            raise ValueError("this simplicial set is not connected")
        graph = self.graph()
        spanning_tree = [e[2] for e in graph.min_spanning_tree()]
        return self.quotient(spanning_tree)

    def is_connected(self):
        """
        True if this simplicial set is connected

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: X = T.disjoint_union(K)
            sage: T.is_connected()
            True
            sage: K.is_connected()
            True
            sage: X.is_connected()
            False
            sage: simplicial_sets.Sphere(0).is_connected()
            False
        """
        return self.graph().is_connected()

    def fundamental_group(self, simplify=True):
        r"""
        Return the fundamental group of this simplicial set.

        INPUT:

        - ``simplify`` (bool, optional True) -- if False, then return a
          presentation of the group in terms of generators and
          relations. If True, the default, simplify as much as GAP is
          able to.

        Algorithm: we compute the edge-path group -- see Section 19 of
        [Kan]_ and :wikipedia:`Fundamental_group`. Choose a spanning
        tree for the connected component of the 1-skeleton containing
        the base point, and then the group's generators are given by
        the non-degenerate edges. There are two types of relations:
        `e=1` if `e` is in the spanning tree, and for every 2-simplex,
        if its edges are `e_0`, `e_1`, and `e_2`, then we impose the
        relation `e_0 e_1^{-1} e_2 = 1`, where we first set `e_i=1` if
        `e_i` is degenerate.

        REFERENCES:

        .. [Kan] \D. M. Kan, *A combinatorial definition of homotopy
           groups*, Ann. Math. (2) 67 (1958), 282-312.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: eight = S1.wedge(copy(S1))
            sage: eight.fundamental_group() # free group on 2 generators
            Finitely presented group < e0, e1 |  >

        The fundamental group of a disjoint union of course depends on
        the choice of base point::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: X = T.disjoint_union(K)

            sage: X_0 = X.set_base_point(X.n_cells(0)[0])
            sage: X_0.fundamental_group().is_abelian()
            True
            sage: X_1 = X.set_base_point(X.n_cells(0)[1])
            sage: X_1.fundamental_group().is_abelian()
            False

        Compute the fundamental group of some classifying spaces::

            sage: RP3 = simplicial_sets.RealProjectiveSpace(3)
            sage: RP3.fundamental_group()
            Finitely presented group < e | e^2 >

            sage: C5 = groups.misc.MultiplicativeAbelian([5])
            sage: BC5 = C5.nerve_n_skeleton(3)
            sage: BC5.fundamental_group()
            Finitely presented group < e0 | e0^5 >

            sage: Sigma3 = groups.permutation.Symmetric(3)
            sage: BSigma3 = Sigma3.nerve_n_skeleton(3)
            sage: pi = BSigma3.fundamental_group(); pi
            Finitely presented group < e0, e1 | e0^2, e1^3, (e0*e1^-1)^2 >
            sage: pi.order()
            6
            sage: pi.is_abelian()
            False
        """
        from sage.groups.free_group import FreeGroup
        from sage.interfaces.gap import gap

        if not self.is_connected() and not self.is_pointed():
            raise ValueError("this simplicial set is not connected, so you must set a base point")
        graph = self.graph()
        if not self.is_connected():
            graph = graph.subgraph(self.base_point())

        edges = [e[2] for e in graph.edges()]
        spanning_tree = [e[2] for e in graph.min_spanning_tree()]
        gens = [e for e in edges if e not in spanning_tree]

        if not gens:
            return gap.TrivialGroup()

        gens_dict = dict(zip(gens, range(len(gens))))
        FG = FreeGroup(len(gens), 'e')
        rels = []

        for f in self.n_cells(2):
            z = dict()
            for i, sigma in enumerate(self.faces(f)):
                if sigma in spanning_tree:
                    z[i] = FG.one()
                elif sigma.is_degenerate():
                    z[i] = FG.one()
                elif sigma in edges:
                    z[i] = FG.gen(gens_dict[sigma])
                else:
                    # sigma is not in the correct connected component.
                    z[i] = FG.one()
            rels.append(z[0]*z[1].inverse()*z[2])
        if simplify:
            return FG.quotient(rels).simplified()
        else:
            return FG.quotient(rels)

    def is_simply_connected(self):
        """
        True if this simplicial set is simply connected.

        .. WARNING::

           Determining simple connectivity is not always possible,
           because it requires determining when a group, as given by
           generators and relations, is trivial. So this conceivably
           may give a false negative in some cases.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: T.is_simply_connected()
            False
            sage: T.suspension().is_simply_connected()
            True
            sage: simplicial_sets.KleinBottle().is_simply_connected()
            False

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: (S2.wedge(S3)).is_simply_connected()
            True
            sage: (S2.disjoint_union(S3)).is_simply_connected()
            False

            sage: C3 = groups.misc.MultiplicativeAbelian([3])
            sage: BC3 = simplicial_sets.ClassifyingSpace(C3, 4)
            sage: BC3.is_simply_connected()
            False
        """
        if not self.is_connected():
            return False
        try:
            return bool(self.fundamental_group().IsTrivial())
        except AttributeError:
            try:
                return self.fundamental_group().order() == 1
            except (NotImplementedError, RuntimeError):
                # I don't know of any simplicial sets for which the
                # code reaches this point, but there are certainly
                # groups for which these errors are raised. 'IsTrivial'
                # works for all of the examples I've seen, though.
                raise ValueError('unable to determine if the fundamental group is trival')

    def connectivity(self):
        """
        The connectivity of this simplicial set.

        The dimension of the first nonzero homotopy group. If simply
        connected, this is the same as the dimension of the first
        nonzero homology group.

        .. WARNING::

           See the warning for the :meth:`is_simply_connected` method.

        The connectivity of a contractible space is ``+Infinity``.

        EXAMPLES::

            sage: simplicial_sets.Sphere(3).connectivity()
            2
            sage: simplicial_sets.Sphere(0).connectivity()
            -1
            sage: simplicial_sets.Simplex(4).connectivity()
            +Infinity
            sage: X = simplicial_sets.Torus().suspension(2)
            sage: X.connectivity()
            2
        """
        if not self.is_connected():
            return Integer(-1)
        if not self.is_simply_connected():
            return Integer(0)
        H = self.homology()
        for i in range(2, self.dimension()+1):
            if i in H and H[i].order() != 1:
                return i-1
        return Infinity

    def _Hom_(self, other, category=None):
        """
        Return the set of simplicial maps between simplicial sets
        ``self`` and ``other``.

        INPUT:

        - ``other`` -- another simplicial set
        - ``category`` -- optional, the category in which to compute
          the maps. By default this is ``SimplicialSets``, and it must
          be a subcategory of this or else an error is raised.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3._Hom_(S2)
            Set of Morphisms from S^3 to S^2 in Category of finite pointed simplicial sets
            sage: Hom(S3, S2)
            Set of Morphisms from S^3 to S^2 in Category of finite pointed simplicial sets
            sage: K4 = simplicial_sets.Simplex(4)
            sage: S3._Hom_(K4)
            Set of Morphisms from S^3 to 4-simplex in Category of finite simplicial sets
        """
        # Error-checking on the ``category`` argument is done when
        # calling Hom(X,Y), so no need to do it again here.
        if category is None:
            if self.is_pointed() and other.is_pointed():
                category = SimplicialSets().Finite().Pointed()
            else:
                category = SimplicialSets().Finite()
        from sage.homology.simplicial_set_morphism import SimplicialSetHomset
        return SimplicialSetHomset(self, other, category=category)

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0)
            sage: w = AbstractSimplex(0)
            sage: degen = v.apply_degeneracies(0)
            sage: tau = AbstractSimplex(2)
            sage: SimplicialSet({tau: (degen, degen, degen), w: None})
            Simplicial set with 3 non-degenerate simplices
            sage: SimplicialSet({w: None})
            Simplicial set with 1 non-degenerate simplex

        Test names and renaming::

            sage: SimplicialSet({w: None}, name='pt')
            pt
            sage: K = SimplicialSet({w: None}, name='pt')
            sage: K.rename('point')
            sage: K
            point
        """
        num = len(self.nondegenerate_simplices())
        if num == 1:
            return "Simplicial set with 1 non-degenerate simplex"
        return "Simplicial set with {} non-degenerate simplices".format(num)


########################################################################
# classes which inherit from SimplicialSet

class SubSimplicialSet(SimplicialSet):
    def __init__(self, data, ambient=None):
        r"""
        A simplicial set as a subsimplicial set of another simplicial set.

        This keeps track of the ambient simplicial set and the
        inclusion map from the subcomplex into it.

        INPUT:

        - ``data`` -- the data defining the subset: a list of
          nondegenerate simplices in the ambient simplicial set

        - ``ambient`` -- the ambient simplicial set. If omitted, use
          the same simplicial set as the subset and the ambient
          complex.

        This is intended to be accessed via the
        :meth:`SimplicialSet.subsimplicial_set` method for the class
        :class:`SimplicialSet`, rather than be called directly.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: K = simplicial_sets.KleinBottle()
            sage: X = S3.disjoint_union(K)
            sage: Y = X.induced_map(0).image() # the S3 summand
            sage: Y.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   Simplicial set with 8 non-degenerate simplices
              Defn: [v_0, sigma_3] --> [v_0, sigma_3]
            sage: Y.ambient_space()
            Simplicial set with 8 non-degenerate simplices
        """
        if ambient is None:
            ambient = self
        if ambient.is_pointed() and ambient.base_point() in data:
            SimplicialSet.__init__(self, data, base_point=ambient.base_point())
        else:
            SimplicialSet.__init__(self, data)
        self._inclusion = self.Hom(ambient)({x:x for x in data})

    def inclusion_map(self):
        """
        The inclusion map from this subsimplicial set into its ambient space.

        EXAMPLES::

            sage: RP6 = simplicial_sets.RealProjectiveSpace(6)
            sage: K = RP6.n_skeleton(2)
            sage: K.inclusion_map()
            Simplicial set morphism:
              From: Simplicial set with 3 non-degenerate simplices
              To:   RP^6
              Defn: [1, f, f * f] --> [1, f, f * f]
        """
        return self._inclusion

    def ambient_space(self):
        """
        The simplicial set of which this is a subsimplicial set.

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: eight = T.wedge_as_subset()
            sage: eight
            Simplicial set with 3 non-degenerate simplices
            sage: eight.fundamental_group()
            Finitely presented group < e0, e1 |  >
            sage: eight.ambient_space()
            Torus
        """
        return self._inclusion.codomain()


class PullbackOfSimplicialSets(SimplicialSet):
    def __init__(self, maps=None):
        r"""
        The pullback obtained from the morphisms ``maps``.

        INPUTS:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets

        If only a single map `f: X \to Y` is given, then return
        `X`. If no maps are given, return the one-point simplicial
        set. Otherwise, given a simplicial set `Y` and maps `f_i: X_i
        \to Y` for `0 \leq i \leq m`, construct the pullback `P`: see
        :wikipedia:`Pullback_(category_theory)`. This is constructed
        as pullbacks of sets for each set of `n`-simplices, so `P_n`
        is the subset of the product `\prod (X_i)_n` consisting of
        those elements `(x_i)` for which `f_i(x_i) = f_j(x_j)` for all
        `i`, `j`.

        This is pointed if the maps `f_i` are.

        EXAMPLES:

        The pullback of a quotient map by a subsimplicial set and the
        base point map gives a simplicial set isomorphic to the
        original subcomplex::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: X = K.pullback(K.quotient_map(), K.base_point_map())
            sage: X.homology() == RP5.n_skeleton(2).homology()
            True

        Pullbacks of identity maps::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: one = S2.Hom(S2).identity()
            sage: P = S2.pullback(one, one)
            sage: P.homology()
            {0: 0, 1: 0, 2: Z}

        The pullback is constructed in terms of the product -- of
        course, the product is a special case of the pullback -- and
        the simplices are named appropriately::

            sage: P.nondegenerate_simplices()
            [(v_0, v_0), (sigma_2, sigma_2)]

        There are defining maps `f_i: X_i \to Y` and induced maps
        `\bar{f}_i: Y_i \to P`. ::

            sage: P.defining_map(0) == one
            True
            sage: P.induced_map(1)
            Simplicial set morphism:
              From: Simplicial set with 2 non-degenerate simplices
              To:   S^2
              Defn: [(v_0, v_0), (sigma_2, sigma_2)] --> [v_0, sigma_2]
            sage: P.induced_map(0).domain() == P
            True
            sage: P.induced_map(0).codomain() == S2
            True

        The pullback `P` also has a universal property: given a
        simplicial set `Z` and maps `g_i: Z \to X_i` such that `f_i
        g_i = f_j g_j` for all `i`, `j`, there is an induced map `g: Z
        \to P` making the appropriate diagram commute: that is,
        `\bar{f}_i g = g_i` for all `i`. For example, given maps `f: X
        \to Y` and `g: X \to Z`, there is an induced map `g: X \to Y
        \times Z`::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: D = S1.cone()      # the cone C(S^1)
            sage: g = D.map_from_X() # map from S^1 to C(S^1)
            sage: P = T.product(D)
            sage: h = P.universal_property(f, g)
            sage: h.domain() == S1
            True
            sage: h.codomain() == P
            True
        """
        from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
        if any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')
        if not maps:
            star = AbstractSimplex(0, name='*')
            SimplicialSet.__init__(self, {star: None}, base_point=star, name='Point')
            self._maps = ()
            self._induced = ()
            self._translation = {}
            return
        if len(maps) == 1:
            f = maps[0]
            if f.is_pointed():
                SimplicialSet.__init__(self, f.domain().face_data(),
                                       base_point=f.domain().base_point())
            else:
                SimplicialSet.__init__(self, f.domain().face_data())
            self._maps = (f,)
            self._induced = (f,)
            self._translation = {(((sigma,()),), sigma)
                                 for sigma in self.nondegenerate_simplices()}
            return
        codomain = maps[0].codomain()
        if any(codomain != f.codomain() for f in maps[1:]):
            raise ValueError('the codomains of the maps must be equal')
        # Now construct the pullback by constructing the product and only
        # keeping the appropriate simplices.
        domains = [f.domain() for f in maps]
        nondegen = [X.nondegenerate_simplices() for X in domains]
        data_factors = [X.face_data() for X in domains]
        # data: dictionary to construct the new simplicial set.
        data = {}
        # translate: keep track of the nondegenerate simplices in the
        # new simplicial set for computing faces: keys are tuples of
        # pairs (sigma, degens), with sigma a nondegenerate simplex in
        # one of the factors, degens the tuple of applied
        # degeneracies. The associated value is the actual simplex in
        # the product.
        translate = {}
        for simplices in itertools.product(*nondegen):
            dims =  [_.dimension() for _ in simplices]
            dim_min = min(dims)
            dim_max = max(dims)
            sum_dims = sum(dims)
            for d in range(dim_max, sum_dims + 1):
                for I in itertools.product(*[Set(range(d)).subsets(d - _.dimension()) for _ in simplices]):
                    if len(set.intersection(*[set(_) for _ in I])) > 0:
                        # To get a nondegenerate face, can't have a
                        # degeneracy in common for all the factors.
                        continue
                    degens = [tuple(sorted(_, reverse=True)) for _ in I]


                    sigma = simplices[0].apply_degeneracies(*degens[0])
                    target = maps[0](sigma)
                    if any(target != f(tau.apply_degeneracies(*degen))
                               for (f, tau, degen) in zip(maps[1:], simplices[1:], degens[1:])):
                        continue

                    simplex_factors = tuple(zip(simplices, tuple(degens)))
                    s = '(' + ', '.join(['{}'.format(_[0].apply_degeneracies(*_[1]))
                                         for _ in simplex_factors]) + ')'
                    simplex = AbstractSimplex(d, name=s)
                    translate[simplex_factors] = simplex
                    # Now compute the faces of simplex.
                    if d == 0:
                        # It's a vertex, so it has no faces.
                        faces = None
                    else:
                        faces = []
                        for i in range(d+1):
                            # Compute d_i on simplex.
                            #
                            # face_degens: tuple of pairs (J, t): J is the
                            # list of degeneracies to apply to the
                            # corresponding entry in simplex_factors, t is
                            # the face map to apply.
                            face_degens = [face_degeneracies(i, _) for _ in degens]
                            face_factors = []
                            new_degens = []
                            for x, Face, face_dict in zip(simplices, face_degens, data_factors):
                                J = Face[0]
                                t = Face[1]
                                if t is None:
                                    face_factors.append(x.nondegenerate())
                                else:
                                    underlying = face_dict[x][t]
                                    temp_degens = underlying.degeneracies()
                                    underlying = underlying.nondegenerate()
                                    J = standardize_degeneracies(*(J + list(temp_degens)))
                                    face_factors.append(underlying)
                                new_degens.append(J)

                            # By the simplicial identities, s_{i_1}
                            # s_{i_2} ... s_{i_n} z (if decreasing) is in
                            # the image of s_{i_k} for each k.
                            #
                            # So find the intersection K of each J, the
                            # degeneracies applied to left_face and
                            # right_face. Then the face will be s_{K}
                            # (s_{J'_L} left_face, s_{J'_R} right_face),
                            # where you get J'_L from J_L by pulling out K
                            # from J_L.
                            #
                            # J'_L is obtained as follows: for each j in
                            # J_L, decrease j by q if q = #{x in K: x < j}
                            K = set.intersection(*[set(J) for J in new_degens])

                            face_degens = []
                            for J in new_degens:
                                new_J = []
                                for j in J:
                                    if j not in K:
                                        q = len([x for x in K if x < j])
                                        new_J.append(j - q)
                                face_degens.append(tuple(new_J))
                            K = sorted(K, reverse=True)
                            underlying_face = translate[tuple(zip(tuple(face_factors), tuple(face_degens)))]
                            faces.append(underlying_face.apply_degeneracies(*K))
                        data[simplex] = faces

        if all(f.is_pointed() for f in maps):
            basept = translate[tuple([(sset.base_point(), ()) for sset in domains])]
            SimplicialSet.__init__(self, data, base_point=basept)
        else:
            SimplicialSet.__init__(self, data)
        # self._factors = tuple(domains)
        self._maps = maps
        # self._translation: tuple converted from dict. keys: tuples
        # of pairs (sigma, degens), with sigma a nondegenerate simplex
        # in one of the factors, degens the tuple of applied
        # degeneracies. The associated value is the actual simplex in
        # the product.
        self._translation = tuple(translate.items())

    def defining_map(self, i):
        r"""
        The `i`-th map defining the pullback.

        INPUT:

        - ``i`` -- integer

        If this pullback was constructed as ``Y.pullback(f_0, f_1, ...)``,
        this returns `f_i`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: Y = K.pullback(K.quotient_map(), K.base_point_map())
            sage: Y.defining_map(1)
            Simplicial set morphism:
              From: Point
              To:   Simplicial set with 4 non-degenerate simplices
              Defn: [*] --> [*]
            sage: Y.defining_map(0).domain()
            RP^5
        """
        return self._maps[i]

    def induced_map(self, i):
        r"""
        The `i`-th induced map from the pullback.

        INPUT:

        - ``i`` -- integer

        If this pullback `P` was constructed as ``Y.pullback(f_0, f_1,
        ...)``, where `f_i: X_i \to Y`, then there are induced maps
        `\bar{f}_i: P \to X_i`. This method constructs `\bar{f}_i`.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: K = RP5.quotient(RP5.n_skeleton(2))
            sage: Y = K.pullback(K.quotient_map(), K.base_point_map())
            sage: Y.induced_map(0)
            Simplicial set morphism:
              From: Simplicial set with 3 non-degenerate simplices
              To:   RP^5
              Defn: [(1, *), (f, Simplex obtained by applying degeneracy s_0 to *), (f * f, Simplex obtained by applying degeneracies s_1 s_0 to *)] --> [1, f, f * f]
            sage: Y.induced_map(1).codomain()
            Point
        """
        f = {}
        for x in self._translation:
            f[x[1]] = x[0][i][0].apply_degeneracies(*x[0][i][1])
        codomain = self.defining_map(i).domain()
        return self.Hom(codomain)(f)

    def universal_property(self, *maps):
        r"""
        Construct the map induced by ``maps``

        INPUT:

        - ``maps`` -- maps from a simplicial set `Z` to the "factors"
          `X_i` forming the pullback.

        If the pullback `P` is formed by maps `f_i: X_i \to Y`, then
        given maps `g_i: Z \to X_i` such that `f_i g_i = f_j g_j` for
        all `i`, `j`, then there is a unique map `g: Z \to P` making
        the appropriate diagram commute. This constructs that map.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: P = S1.product(T)
            sage: P.universal_property(S1.Hom(S1).identity(), f)
            Simplicial set morphism:
              From: S^1
              To:   Simplicial set with 26 non-degenerate simplices
              Defn: [v_0, sigma_1] --> [(v_0, (v_0, v_0)), (sigma_1, (sigma_1, Simplex obtained by applying degeneracy s_0 to v_0))]
        """
        domain = maps[0].domain()
        if any(g.domain() != domain for g in maps[1:]):
            raise ValueError('the maps do not all have the same codomain')
        composite = self._maps[0] * maps[0]
        if any(f*g != composite for f,g in zip(self._maps[1:], maps[1:])):
            raise ValueError('the maps are not compatible')
        data = {}
        translate = dict(self._translation)
        for sigma in domain.nondegenerate_simplices():
            target = tuple([(f(sigma).nondegenerate(), tuple(f(sigma).degeneracies()))
                               for f in maps])
            data[sigma] = translate[target]
        return domain.Hom(self)(data)


class ProductOfSimplicialSets(PullbackOfSimplicialSets):
    r"""
    The product of simplicial sets.

    INPUT:

    - ``factors`` -- a list or tuple of simplicial sets

    Return the product of the simplicial sets in ``factors``.

    If `X` and `Y` are simplicial sets, then their product `X \times
    Y` is defined to be the simplicial set with `n`-simplices `X_n
    \times Y_n`.  Therefore the simplices in the product have the form
    `(s_I \sigma, s_J \tau)`, where `s_I = s_{i_1} ... s_{i_p}` and
    `s_J = s_{j_1} ... s_{j_q}` are composites of degeneracy maps,
    written in decreasing order.  Such a simplex is nondegenerate if
    the indices `I` and `J` are disjoint. Therefore if `\sigma` and
    `\tau` are nondegenerate simplices of dimensions `m` and `n`, in
    the product they will lead to nondegenerate simplices up to
    dimension `m+n`, and no further.

    This extends in the more or less obvious way to products with more
    than two factors: with three factors, a simplex `(s_I \sigma, s_J
    \tau, s_K \rho)` is nondegenerate if `I \cap J \cap K` is empty,
    etc.

    If a simplicial set is constructed as a product, the factors are
    recorded and are accessible via the method :meth:`factors`. If it
    is constructed as a product and then copied, this information is
    lost.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: v = AbstractSimplex(0, name='v')
        sage: w = AbstractSimplex(0, name='w')
        sage: e = AbstractSimplex(1, name='e')
        sage: X = SimplicialSet({e: (v, w)})
        sage: square = X.product(X)

    ``square`` is now the standard triangulation of the square: 4
    vertices, 5 edges (the four on the border plus the diagonal), 2
    triangles::

        sage: square.f_vector()
        [4, 5, 2]

        sage: S1 = simplicial_sets.Sphere(1)
        sage: T = S1.product(S1)
        sage: T.homology(reduced=False)
        {0: Z, 1: Z x Z, 2: Z}

    Since ``S1`` is pointed, so is ``T``::

        sage: S1.is_pointed()
        True
        sage: S1.base_point()
        v_0
        sage: T.is_pointed()
        True
        sage: T.base_point()
        (v_0, v_0)

        sage: S2 = simplicial_sets.Sphere(2)
        sage: S3 = simplicial_sets.Sphere(3)
        sage: Z = S2.product(S3)
        sage: Z.homology()
        {0: 0, 1: 0, 2: Z, 3: Z, 4: 0, 5: Z}

    Factors and copying::

        sage: Z.factors() == (S2, S3)
        True
        sage: copy(Z).factors() # copying forgets the product structure
        Traceback (most recent call last):
        ...
        AttributeError: 'SimplicialSet_with_category' object has no attribute 'factors'
    """
    def __init__(self, factors=None):
        r"""
        The product of simplicial sets.

        See :class:`ProductOfSimplicialComplexes` for more
        information.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet, ProductOfSimplicialSets
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: W = ProductOfSimplicialSets([X, X, X])
            sage: W.homology()
            {0: 0, 1: Z x Z x Z, 2: Z x Z x Z, 3: Z}
            sage: W.is_pointed()
            False

            sage: X = X.set_base_point(v)
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1)
            sage: Y = SimplicialSet({f: (v,w)}, base_point=w)
            sage: Z = ProductOfSimplicialSets([Y, X])
            sage: Z.is_pointed()
            True
            sage: Z.base_point()
            (w, v)
        """
        PullbackOfSimplicialSets.__init__(self, [space.constant_map() for space in factors])
        self._factors = tuple([f.domain() for f in self._maps])

    def factors(self):
        """
        The factors of this product of simplicial sets.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2.product(S3).factors() == (S2, S3)
            True
            sage: S2.product(S3).factors()[0]
            S^2
        """
        return self._factors

    def factor(self, i, as_subset=False):
        r"""
        The $i$-th factor of the product.

        INPUT:

        - ``i`` -- integer, the index of the factor

        - ``as_subset`` -- boolean, optional (default ``False``)

        If ``as_subset`` is ``True``, return the $i$-th factor as a
        subsimplicial set of the product, identifying it with its
        product with the base point in each other factor. As a
        subsimplicial set, it comes equipped with an inclusion
        map. This option will raise an error if any factor does not
        have a base point.

        If ``as_subset`` is ``False``, return the $i$-th factor in
        its original form as a simplicial set.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: K = S2.product(S3)
            sage: K.factor(0)
            S^2

            sage: K.factor(0, as_subset=True)
            Simplicial set with 2 non-degenerate simplices
            sage: K.factor(0, as_subset=True).homology()
            {0: 0, 1: 0, 2: Z}

            sage: K.factor(0) is S2
            True
            sage: K.factor(0, as_subset=True) is S2
            False
        """
        if as_subset:
            if any(not _.is_pointed() for _ in self.factors()):
                raise ValueError('"as_subset=True" is only valid '
                                 'if each factor is pointed')

            basept_factors = [sset.base_point() for sset in self.factors()]
            basept_factors = basept_factors[:i] + basept_factors[i+1:]
            to_factors = dict((v,k) for k,v in self._translation)
            simps = []
            for x in self.nondegenerate_simplices():
                simplices = [sigma[0] for sigma in to_factors[x]]
                if simplices[:i] + simplices[i+1:] == basept_factors:
                    simps.append(x)
            return self.subsimplicial_set(simps)
        return self.factors()[i]

    def projection_map(self, i):
        """
        The map projecting onto the $i$-th factor.

        INPUT:

        - ``i`` -- integer, the index of the projection map

        EXAMPLES::

            sage: T = simplicial_sets.Torus()
            sage: f_0 = T.projection_map(0)
            sage: f_1 = T.projection_map(1)
            sage: m_0 = f_0.induced_homology_morphism().to_matrix(1) # matrix in dim 1
            sage: m_1 = f_1.induced_homology_morphism().to_matrix(1)
            sage: m_0.rank()
            1
            sage: m_0 == m_1
            False
        """
        return self.induced_map(i)

    def wedge_as_subset(self):
        """
        The wedge as a subsimplicial set of this product of pointed
        simplicial sets.

        This will raise an error if any factor is not pointed.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: w = AbstractSimplex(0, name='w')
            sage: f = AbstractSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: P = X.product(Y)
            sage: W = P.wedge_as_subset()
            sage: W.nondegenerate_simplices()
            [(v, w),
             (Simplex obtained by applying degeneracy s_0 to v, f),
             (e, Simplex obtained by applying degeneracy s_0 to w)]
            sage: W.homology()
            {0: 0, 1: Z x Z}
        """
        basept_factors = [sset.base_point() for sset in self.factors()]
        to_factors = dict((v,k) for k,v in self._translation)
        simps = []
        for x in self.nondegenerate_simplices():
            simplices = to_factors[x]
            not_base_pt = 0
            for sigma, star in zip(simplices, basept_factors):
                if not_base_pt > 1:
                    continue
                if sigma[0].nondegenerate() != star:
                    not_base_pt += 1
            if not_base_pt <= 1:
                simps.append(x)
        return self.subsimplicial_set(simps)

    def fat_wedge(self):
        """
        The fat wedge as a subsimplicial set of this product of pointed
        simplicial sets.

        The fat wedge consists of those terms where at least one
        factor is the base point. Thus with two factors this is the
        ordinary wedge, but with more factors, it is larger.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: X = S1.product(S1, S1)
            sage: W = X.fat_wedge()
            sage: W.homology()
            {0: 0, 1: Z x Z x Z, 2: Z x Z x Z}
        """
        basept_factors = [sset.base_point() for sset in self.factors()]
        to_factors = dict((v,k) for k,v in self._translation)
        N = len(basept_factors)
        simps = []
        for x in self.nondegenerate_simplices():
            simplices = to_factors[x]
            combined = zip(simplices, basept_factors)
            if any(sigma[0] == pt for (sigma,pt) in combined):
                simps.append(x)
        return self.subsimplicial_set(simps)


class PushoutOfSimplicialSets(SimplicialSet):
    def __init__(self, maps=None, vertex_name=None):
        r"""
        The pushout obtained from the morphisms ``maps``.

        INPUTS:

        - ``maps`` -- a list or tuple of morphisms of simplicial sets
        - ``vertex_name`` -- optional, default ``None``

        If only a single map `f: X \to Y` is given, then return
        `Y`. If no maps are given, return the empty simplicial
        set. Otherwise, given a simplicial set `X` and maps `f_i: X
        \to Y_i` for `0 \leq i \leq m`, construct the pushout `P`: see
        :wikipedia:`Pushout_(category_theory)`. This is constructed as
        pushouts of sets for each set of `n`-simplices, so `P_n` is
        the disjoint union of the sets `(Y_i)_n`, with elements
        `f_i(x)` identified for `n`-simplex `x` in `X`.

        Simplices in the pushout are given names as follows: if a
        simplex comes from a single `Y_i`, it inherits its
        name. Otherwise it must come from a simplex (or several) in
        `X`, and then it inherits one of those names, and it should be
        the first alphabetically. For examnple, if vertices `v`, `w`,
        and `z` in `X` are glued together, then the resulting vertex
        in the pushout will be called `v`.

        Base points are taken care of automatically: if each of the
        maps `f_i` is pointed, so is the pushout. If `X` is a point or
        if `X` is nonempty and any of the spaces `Y_i` is a point, use
        those for the base point. In all of these cases, if
        ``vertex_name`` is ``None``, generate the name of the base
        point automatically; otherwise, use ``vertex_name`` for its
        name.

        In all other cases, the pushout is not pointed.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: a = AbstractSimplex(0, name='a')
            sage: b = AbstractSimplex(0, name='b')
            sage: c = AbstractSimplex(0, name='c')
            sage: e0 = AbstractSimplex(1, name='e_0')
            sage: e1 = AbstractSimplex(1, name='e_1')
            sage: e2 = AbstractSimplex(1, name='e_2')
            sage: X = SimplicialSet({e2: (b, a)})
            sage: Y0 = SimplicialSet({e2: (b,a), e0: (c,b), e1: (c,a)})
            sage: Y1 = simplicial_sets.Simplex(0)
            sage: f0_data = {a:a, b:b, e2: e2}
            sage: v = Y1.n_cells(0)[0]
            sage: f1_data = {a:v, b:v, e2:v.apply_degeneracies(0)}
            sage: f0 = X.Hom(Y0)(f0_data)
            sage: f1 = X.Hom(Y1)(f1_data)
            sage: P = X.pushout(f0, f1)
            sage: P.nondegenerate_simplices()
            [a, c, e_0, e_1]

        There are defining maps `f_i: X \to Y_i` and induced maps
        `\bar{f}_i: Y_i \to P`. ::

            sage: P.defining_map(0) == f0
            True
            sage: P.induced_map(1)
            Simplicial set morphism:
              From: 0-simplex
              To:   Simplicial set with 4 non-degenerate simplices
              Defn: [(0,)] --> [a]
            sage: P.induced_map(0).domain() == Y0
            True
            sage: P.induced_map(0).codomain() == P
            True

        An inefficient way of constructing a suspension for an
        unpointed set: take the pushout of two copies of the inclusion
        map `X \to CX`::

            sage: T = simplicial_sets.Torus()
            sage: T = T.unset_base_point()
            sage: CT = T.cone()
            sage: inc = CT.X_as_subset().inclusion_map()
            sage: P = T.pushout(inc, inc)
            sage: P
            Simplicial set with 20 non-degenerate simplices
            sage: P.homology()
            {0: 0, 1: 0, 2: Z x Z, 3: Z}

        It is more efficient to construct the suspension as the
        quotient `CX/X`::

            sage: CT.quotient(CT.X_as_subset())
            Simplicial set with 8 non-degenerate simplices

        It is more efficient still if the original simplicial set has
        a base point::

            sage: T = simplicial_sets.Torus()
            sage: T.suspension()
            Simplicial set with 6 non-degenerate simplices

        The pushout `P` also has a universal property: given a
        simplicial set `Z` and maps `g_i: Y_i \to Z` such that `g_i
        f_i = g_j f_j` for all `i`, `j`, there is an induced map `g: P
        \to Z` making the appropriate diagram commute: that is, `g
        \bar{f}_i = g_i` for all `i`. ::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = S1.product(S1)
            sage: K = T.factor(0, as_subset=True)
            sage: f = S1.Hom(T)({S1.n_cells(0)[0]:K.n_cells(0)[0], S1.n_cells(1)[0]:K.n_cells(1)[0]})
            sage: W = S1.wedge(T) # wedge, constructed as a pushout

        The maps `f: S^1 \to T` and `1: T \to T` should induce a map `S^1 \vee T \to T`::

            sage: g = W.universal_property(f, Hom(T,T).identity())
            sage: g.domain() == W
            True
            sage: g.codomain() == T
            True

            sage: S1 = simplicial_sets.Sphere(1)
            sage: pt = simplicial_sets.Point()
            sage: bouquet = pt.pushout(S1.base_point_map(), S1.base_point_map(), S1.base_point_map())
            sage: bouquet.homology(1)
            Z x Z x Z
        """
        from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
        if any(not isinstance(f, SimplicialSetMorphism) for f in maps):
            raise ValueError('the maps must be morphisms of simplicial sets')
        if not maps:
            SimplicialSet.__init__(self, {})
            self._maps = ()
            self._induced = ()
            return
        domain = maps[0].domain()
        if len(maps) == 1:
            # f: X --> Y
            f = maps[0]
            codomain = f.codomain()
            if f.is_pointed():
                base_point=codomain.base_point()
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            elif len(domain.nondegenerate_simplices()) == 1:
                # X is a point.
                base_point = f(domain().n_cells(0)[0])
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            elif len(codomain.nondegenerate_simplices()) == 1:
                # Y is a point.
                base_point = codomain.n_cells(0)[0]
                if vertex_name is not None:
                    base_point.rename(vertex_name)
                SimplicialSet.__init__(self, codomain.face_data(),
                                       base_point=base_point)
            else:
                SimplicialSet.__init__(self, codomain.face_data())
            self._maps = (f,)
            self._induced = (f,)
            return
        if any(domain != f.domain() for f in maps[1:]):
            raise ValueError('the domains of the maps must be equal')
        # Data to define the pushout:
        data = {}
        codomains = [f.codomain() for f in maps]
        # spaces: indexed list of spaces. Entries are of the form
        # (space, int) where int=-1 for the domain, and for the
        # codomains, int is the corresponding index.
        spaces = [(Y,i-1) for (i,Y) in enumerate([domain] + codomains)]
        # Dictionaries to translate from simplices in domain,
        # codomains to simplices in the pushout. The keys are of the
        # form (space, int). int=-1 for the domain, and for the
        # codomains, int is the corresponding index.
        _to_P = {Y:{} for Y in spaces}
        max_dim = max(Y.dimension() for Y in codomains)
        for n in range(1 + max_dim):
            # Now we impose an equivalence relation on the simplices,
            # setting x equivalent to f_i(x) for each simplex x in X
            # and each defining map f_i. We do this by constructing a
            # graph and finding its connected components: the vertices
            # of the graph are be the n-cells of X and the Y_i, and
            # there are edges from x to f_i(x).
            vertices = []
            for (Y,i) in spaces:
                vertices.extend([(cell,i) for cell in Y.n_cells(n)])
            edges = []
            for x in domain.n_cells(n):
                edges.extend([[(x,-1), (f(x),i)] for (i,f) in enumerate(maps)])
            G = Graph([vertices, edges], format='vertices_and_edges')
            data[n] = [set(_) for _ in G.connected_components()]
        # data is now a dictionary indexed by dimension, and data[n]
        # consists of sets of n-simplices of the domain and the
        # codomains, each set an equivalence class of n-simplices
        # under the gluing. So if any element of one of those sets is
        # degenerate, we can throw the whole thing away. Otherwise, we
        # can choose a representative to compute the faces.
        simplices = {}
        for dim in sorted(data.keys()):
            for s in data[dim]:
                degenerate = any(sigma[0].is_degenerate() for sigma in s)
                if degenerate:
                    # Identify the degeneracies involved.
                    degens = []
                    for (sigma,j) in s:
                        if len(sigma.degeneracies()) > len(degens):
                            degens = sigma.degeneracies()
                            underlying = sigma
                            space = spaces[j+1]
                            old = _to_P[space][sigma.nondegenerate()]
                    for (sigma,j) in s:
                        # Now update the _to_P[space] dictionaries.
                        space = spaces[j+1]
                        _to_P[space][sigma] = old.apply_degeneracies(*degens)
                else: # nondegenerate
                    if len(s) == 1:
                        name = str(list(s)[0][0])
                    else:
                        # Choose a name from a simplex in domain.
                        for (sigma,j) in s:
                            if j == -1:
                                name = str(sigma)
                                break
                    new = AbstractSimplex(dim, name=name)
                    if dim == 0:
                        faces = None
                    for (sigma,j) in s:
                        space = spaces[j+1]
                        _to_P[space][sigma] = new
                        if dim > 0:
                            faces = [_to_P[space][tau.nondegenerate()].apply_degeneracies(*tau.degeneracies())
                                     for tau in space[0].faces(sigma)]
                    simplices[new] = faces

        some_Y_is_pt = False
        if len(domain.nondegenerate_simplices()) > 1:
            # Only investigate this if X is not empty and not a point.
            for (Y,i) in spaces:
                if len(Y.nondegenerate_simplices()) == 1:
                    some_Y_is_pt = True
                    break
        if len(domain.nondegenerate_simplices()) == 1:
            # X is a point.
            base_point = _to_P[(domain,-1)][domain.n_cells(0)[0]]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet.__init__(self, simplices, base_point=base_point)
        elif some_Y_is_pt:
            # We found (Y,i) above.
            base_point = _to_P[(Y,i)][Y.n_cells(0)[0]]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet.__init__(self, simplices, base_point=base_point)
        elif all(f.is_pointed() for f in maps):
            pt = _to_P[(codomains[0],0)][codomains[0].base_point()]
            if any(_to_P[(Y,i)][Y.base_point()] != pt for (Y,i) in spaces[2:]):
                raise ValueError('something unexpected went wrong with base points')
            base_point = _to_P[(domain,-1)][domain.base_point()]
            if vertex_name is not None:
                base_point.rename(vertex_name)
            SimplicialSet.__init__(self, simplices, base_point=base_point)
        else:
            SimplicialSet.__init__(self, simplices)
        # The relevant maps:
        self._maps = maps
        self._induced = tuple([Y.Hom(self)(_to_P[(Y,i)]) for (Y,i) in spaces[1:]])

    def defining_map(self, i):
        r"""
        The `i`-th map defining the pushout.

        INPUT:

        - ``i`` -- integer

        If this pushout was constructed as ``X.pushout(f_0, f_1, ...)``,
        this returns `f_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: X = S1.wedge(T) # a pushout
            sage: X.defining_map(0)
            Simplicial set morphism:
              From: Point
              To:   S^1
              Defn: [*] --> [v_0]
            sage: X.defining_map(1).domain()
            Point
            sage: X.defining_map(1).codomain()
            Torus
        """
        return self._maps[i]

    def induced_map(self, i):
        r"""
        The $i$-th induced map to the pushout.

        INPUT:

        - ``i`` -- integer

        If this pushout `Z` was constructed as ``X.pushout(f_0, f_1, ...)``,
        where `f_i: X \to Y_i`, then there are induced maps
        `\bar{f}_i: Y_i \to Z`. This method constructs `\bar{f}_i`.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: T = simplicial_sets.Torus()
            sage: X = S1.disjoint_union(T) # a pushout
            sage: X.induced_map(0)
            Simplicial set morphism:
              From: S^1
              To:   Simplicial set with 8 non-degenerate simplices
              Defn: [v_0, sigma_1] --> [v_0, sigma_1]
            sage: X.induced_map(1).domain()
            Torus
            sage: X.induced_map(1).codomain()
            Simplicial set with 8 non-degenerate simplices
        """
        return self._induced[i]

    def universal_property(self, *maps):
        r"""
        Construct the map induced by ``maps``

        INPUT:

        - ``maps`` -- maps "factors" `Y_i` forming the pushout to a
          fixed simplicial set `Z`.

        If the pushout `P` is formed by maps `f_i: X \to Y_i`, then
        given maps `g_i: Y_i \to Z` such that `g_i f_i = g_j f_j` for
        all `i`, `j`, then there is a unique map `g: P \to Z` making
        the appropriate diagram commute. This constructs that map.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: w = AbstractSimplex(0, name='w')
            sage: x = AbstractSimplex(0, name='x')
            sage: evw = AbstractSimplex(1, name='vw')
            sage: evx = AbstractSimplex(1, name='vx')
            sage: ewx = AbstractSimplex(1, name='wx')
            sage: X = SimplicialSet({evw: (w, v), evx: (x, v)})
            sage: Y_0 = SimplicialSet({evw: (w, v), evx: (x, v), ewx: (x, w)})
            sage: Y_1 = SimplicialSet({evx: (x, v)})

            sage: f_0 = Hom(X, Y_0)({v:v, w:w, x:x, evw:evw, evx:evx})
            sage: f_1 = Hom(X, Y_1)({v:v, w:v, x:x, evw:v.apply_degeneracies(0), evx:evx})
            sage: P = X.pushout(f_0, f_1)

            sage: one = Hom(Y_1, Y_1).identity()
            sage: g = Hom(Y_0, Y_1)({v:v, w:v, x:x, evw:v.apply_degeneracies(0), evx:evx, ewx:evx})
            sage: P.universal_property(g, one)
            Simplicial set morphism:
              From: Simplicial set with 4 non-degenerate simplices
              To:   Simplicial set with 3 non-degenerate simplices
              Defn: [v, x, vx, wx] --> [v, x, vx, vx]
        """
        codomain = maps[0].codomain()
        if any(g.codomain() != codomain for g in maps[1:]):
            raise ValueError('the maps do not all have the same codomain')
        composite = maps[0] * self._maps[0]
        if any(g*f != composite for g,f in zip(maps[1:], self._maps[1:])):
            raise ValueError('the maps are not compatible')
        data = {}
        for i,g in enumerate(maps):
            f_i_dict = self.induced_map(i)._dictionary
            for sigma in f_i_dict:
                tau = f_i_dict[sigma]
                # For sigma_i in Y_i, define the map G by
                # G(\bar{f}_i)(sigma_i) = g_i(sigma_i).
                if tau not in data:
                    data[tau] = g(sigma)
        return self.Hom(codomain)(data)


class QuotientOfSimplicialSet(PushoutOfSimplicialSets):
    def __init__(self, subcomplex, vertex_name='*'):
        r"""
        The quotient of a simplicial set by a subsimplicial set

        INPUT:

        - ``subcomplex`` -- a subcomplex (= subsimplicial set) of a
          simplicial set
        - ``vertex_name`` -- optional, default ``'*'``

        A subcomplex `A` comes equipped with the inclusion map `A \to
        X` to its ambient complex `X`, so the ambient complex is not
        specified explicitly. This constructs the quotient `X/A`,
        collapsing `A` to a point.

        The resulting point is called ``vertex_name``, which is
        ``'*'`` by default.

        Quotients come equipped with a quotient map, available using
        the :meth:`quotient_map` method.

        EXAMPLES::

            sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
            sage: RP2 = RP5.n_skeleton(2)
            sage: RP5_2 = RP5.quotient(RP2)
            sage: RP5_2
            Simplicial set with 4 non-degenerate simplices
            sage: RP5_2.quotient_map()
            Simplicial set morphism:
              From: RP^5
              To:   Simplicial set with 4 non-degenerate simplices
              Defn: [1, f, f * f, f * f * f, f * f * f * f, f * f * f * f * f] --> [*, Simplex obtained by applying degeneracy s_0 to *, Simplex obtained by applying degeneracies s_1 s_0 to *, f * f * f, f * f * f * f, f * f * f * f * f]
        """
        PushoutOfSimplicialSets.__init__(self, [subcomplex.inclusion_map(),
                                                subcomplex.constant_map()],
                                         vertex_name=vertex_name)

    def quotient_map(self):
        """
        The quotient map from the original simplicial set to the quotient.

        EXAMPLES::

            sage: K = simplicial_sets.Simplex(1)
            sage: S1 = K.quotient(K.n_skeleton(0))
            sage: q = S1.quotient_map()
            sage: q
            Simplicial set morphism:
              From: 1-simplex
              To:   Simplicial set with 2 non-degenerate simplices
              Defn: [(0,), (1,), (0, 1)] --> [*, *, (0, 1)]
            sage: q.domain() == K
            True
            sage: q.codomain() == S1
            True
        """
        return self.induced_map(0)


class WedgeOfSimplicialSets(PushoutOfSimplicialSets):
    def __init__(self, factors=None):
        r"""
        The wedge sum of pointed simplicial sets.

        INPUT:

        - ``factors`` -- a list or tuple of simplicial sets

        Return the wedge of the simplicial sets in ``factors``.

        All of the simplicial sets must be pointed. Their base points
        need not be the same, but during the construction, the base
        points will become identified and renamed '*'.

        The wedge comes equipped with maps to and from each factor, or
        actually, maps from each factor, and maps to simplicial sets
        isomorphic to each factor. The codomains of the latter maps
        are quotients of the wedge, not identical to the original
        factors.

        EXAMPLES::

            sage: CP2 = simplicial_sets.ComplexProjectiveSpace(2)
            sage: K = simplicial_sets.KleinBottle()
            sage: W = CP2.wedge(K)
            sage: W.homology()
            {0: 0, 1: Z x C2, 2: Z, 3: 0, 4: Z}

            sage: W.inclusion(1)
            Simplicial set morphism:
              From: Simplicial set with 6 non-degenerate simplices
              To:   Simplicial set with 14 non-degenerate simplices
              Defn: [Delta_{0,0}, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}] --> [*, Delta_{1,0}, Delta_{1,1}, Delta_{1,2}, Delta_{2,0}, Delta_{2,1}]

            sage: W.projection(0).domain()
            Simplicial set with 14 non-degenerate simplices
            sage: W.projection(0).codomain() # copy of CP^2
            Simplicial set with 9 non-degenerate simplices
            sage: W.projection(0).codomain().homology()
            {0: 0, 1: 0, 2: Z, 3: 0, 4: Z}

        An error occurs if any of the factors is not pointed::

            sage: CP2.wedge(simplicial_sets.Simplex(1))
            Traceback (most recent call last):
            ...
            ValueError: the simplicial sets must be pointed
        """
        if any(not space.is_pointed() for space in factors):
            raise ValueError('the simplicial sets must be pointed')
        PushoutOfSimplicialSets.__init__(self, [space.base_point_map()
                                                for space in factors])
        self.base_point().rename('*')
        self._factors = factors

    def inclusion(self, i):
        """
        Inclusion map of the $i$-th factor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: K1 = copy(S1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: W = S1.wedge(S2, K1)
            sage: W.inclusion(1)
            Simplicial set morphism:
              From: S^2
              To:   Simplicial set with 4 non-degenerate simplices
              Defn: [v_0, sigma_2] --> [*, sigma_2]
            sage: W.inclusion(0).domain()
            S^1
            sage: W.inclusion(2).domain()
            Simplicial set with 2 non-degenerate simplices
        """
        return self.induced_map(i)

    def projection(self, i):
        """
        Projection map onto the $i$-th factor.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: K1 = copy(S1)
            sage: S2 = simplicial_sets.Sphere(2)
            sage: W = S1.wedge(S2, K1)
            sage: W.projection(1)
            Simplicial set morphism:
              From: Simplicial set with 4 non-degenerate simplices
              To:   Simplicial set with 2 non-degenerate simplices
              Defn: [*, sigma_1, sigma_1', sigma_2] --> [*, Simplex obtained by applying degeneracy s_0 to *, Simplex obtained by applying degeneracy s_0 to *, sigma_2]
            sage: W.projection(1).image().homology(1)
            0
            sage: W.projection(1).image().homology(2)
            Z
        """
        m = len(self._factors)
        simplices = ([self.inclusion(j).image().nondegenerate_simplices() for j in range(i)]
                     + [self.inclusion(j).image().nondegenerate_simplices() for j in range(i+1,m)])
        return self.quotient(list(itertools.chain(*simplices))).quotient_map()


class ConeOfSimplicialSet(SimplicialSet):
    def __init__(self, X):
        r"""
        The unreduced cone on a simplicial set.

        INPUT:

        - ``X`` -- return the cone on this simplicial set.

        Add a point `*` (which will become the base point) and for
        each simplex `\sigma` in `X`, add both `\sigma` and a simplex
        made up of `*` and `\sigma` (topologically, form the join of
        `*` and `\sigma`).

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone() # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *
        """
        star = AbstractSimplex(0, name='*')
        data = {}
        data[star] = None
        # Dictionary for translating old simplices to new: keys are
        # old simplices, corresponding value is the new simplex
        # (sigma, *).
        new_simplices = {None: star}
        for sigma in X.nondegenerate_simplices():
            new = AbstractSimplex(sigma.dimension()+1,
                                       name='({},*)'.format(sigma))
            if sigma.dimension() == 0:
                data[sigma] = None
                data[new] = (star, sigma)
            else:
                sigma_faces = X.face_data()[sigma]
                data[sigma] = sigma_faces
                new_faces = [new_simplices[face.nondegenerate()].apply_degeneracies(*face.degeneracies())
                             for face in sigma_faces]
                data[new] = (new_faces + [sigma])
            new_simplices[sigma] = new
        SimplicialSet.__init__(self, data, base_point=star)
        self._X = X

    def X_as_subset(self):
        """
        If this is the cone `CX` on `X`, return `X` as a subsimplicial set.

        EXAMPLES::

            sage: X = simplicial_sets.RealProjectiveSpace(4).unset_base_point()
            sage: Y = X.cone()
            sage: Y.X_as_subset()
            Simplicial set with 5 non-degenerate simplices
            sage: Y.X_as_subset() == X
            True
        """
        X = self._X
        return self.subsimplicial_set(X.nondegenerate_simplices())


class ReducedConeOfSimplicialSet(QuotientOfSimplicialSet):
    def __init__(self, X):
        r"""
        The reduced cone on a simplicial set.

        INPUT:

        - ``X`` -- return the cone on this simplicial set.

        Start with the unreduced cone: take `X` and add a point `*`
        (which will become the base point) and for each simplex
        `\sigma` in `X`, add both `\sigma` and a simplex made up of ``
        and `\sigma` (topologically, form the join of `*` and
        `\sigma`).

        Now reduce: take the quotient by the 1-simplex connecting the
        old base point to the new one.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
            sage: v = AbstractSimplex(0, name='v')
            sage: e = AbstractSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: X = X.set_base_point(v)
            sage: CX = X.cone()  # indirect doctest
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]
        """
        CX = ConeOfSimplicialSet(X)
        edge_faces = sorted([CX.base_point(), X.base_point()])
        for t in CX.n_cells(1):
            if sorted(CX.faces(t)) == edge_faces:
                edge = t
                break
        QuotientOfSimplicialSet.__init__(self, CX.subsimplicial_set([edge]))
        self._X = X

    def map_from_X(self):
        r"""
        If this is the cone `\tilde{C}X` on `X`, return the map from `X`.

        The map is defined to be the composite `X \to CX \to
        \tilde{C}X`.  This is used by the
        :meth:`SimplicialSet.suspension` method to construct the
        reduced suspension: take the quotient of the reduced cone by
        the image of `X` therein.

        EXAMPLES::

            sage: S3 = simplicial_sets.Sphere(3)
            sage: CS3 = S3.cone()
            sage: CS3.map_from_X()
            Simplicial set morphism:
              From: S^3
              To:   Simplicial set with 3 non-degenerate simplices
              Defn: [v_0, sigma_3] --> [*, sigma_3]
        """
        quotient_map = self.quotient_map()
        unreduced = quotient_map.domain()
        temp_map = unreduced.X_as_subset().inclusion_map()
        X = self._X
        incl = X.Hom(unreduced)(temp_map._dictionary)
        return quotient_map * incl


########################################################################
# Functions for manipulating face and degeneracy maps.

def standardize_degeneracies(*L):
    r"""
    INPUT:

    - ``L`` -- list of integers, representing a composition of
      degeneracies in a simplicial set.

    OUTPUT: an equivalent list of degeneracies, standardized to be
    written in decreasing order, using the simplicial identity

    .. MATH::

        s_i s_j = s_{j+1} s_i \ \   \text{if }   i \leq j.

    For example, `s_0 s_2 = s_3 s_0` and `s_0 s_0 = s_1 s_0`.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import standardize_degeneracies
        sage: standardize_degeneracies(0, 0)
        (1, 0)
        sage: standardize_degeneracies(0, 0, 0, 0)
        (3, 2, 1, 0)
        sage: standardize_degeneracies(1, 2)
        (3, 1)

    TESTS::

        sage: standardize_degeneracies()
        ()
        sage: standardize_degeneracies(2, -1)
        Traceback (most recent call last):
        ...
        ValueError: degeneracies are indexed by non-negative integers
        sage: standardize_degeneracies([2, 1])
        Traceback (most recent call last):
        ...
        TypeError: degeneracies are indexed by non-negative integers; do not use an explicit list or tuple
    """
    J = list(L)
    for m in J:
        try:
            if Integer(m) < 0:
                raise ValueError('degeneracies are indexed by non-negative integers')
        except TypeError:
            # Likely if called via standard_degeneracies([1,2,3])
            # rather than          standard_degeneracies(1,2,3).
            raise TypeError('degeneracies are indexed by non-negative integers; do not use an explicit list or tuple')
    inadmissible = True
    while inadmissible:
        inadmissible = False
        for idx in range(len(J)-1):
            if J[idx] <= J[idx + 1]:
                inadmissible = True
                tmp = J[idx]
                J[idx] = J[idx + 1] + 1
                J[idx + 1] = tmp
    return tuple(J)

def all_degeneracies(n, l=1):
    r"""
    List of all composites of degeneracies (written in "admissible"
    form, i.e., as a strictly decreasing sequence) of length `l` on an
    `n`-simplex.

    INPUT:

    - ``n``, ``l`` -- integers

    On an `n`-simplex, one may apply the degeneracies `s_i` for `0
    \leq i \leq n`. Then on the resulting `n+1`-simplex, one may apply
    `s_i` for `0 \leq i \leq n+1`, and so on. But one also has to take
    into account the simplicial identity

    .. MATH::

        s_i s_j = s_{j+1} s_i \ \   \text{if }   i \leq j.

    There are `\binom{l+n}{n}` such composites: each non-degenerate
    `n`-simplex leads to `\binom{l+n}{n}` degenerate `l+n` simplices.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import all_degeneracies
        sage: all_degeneracies(0, 3)
        {(2, 1, 0)}
        sage: all_degeneracies(1, 1)
        {(0,), (1,)}
        sage: all_degeneracies(1, 3)
        {(2, 1, 0), (3, 1, 0), (3, 2, 0), (3, 2, 1)}
    """
    if l == 0:
        return set(())
    if l == 1:
        return set([tuple([_]) for _ in range(n+1)])
    ans = set()
    for i in range(n+l):
        ans.update(set([tuple(standardize_degeneracies(*([i] + list(_)))) for _ in all_degeneracies(n, l-1)]))
    return ans

def standardize_face_maps(*L):
    r"""
    INPUT:

    - ``L`` -- list of integers, representing a composition of
      face maps in a simplicial set.

    OUTPUT: an equivalent list of face maps, standardized to be
    written in non-increasing order, using the simplicial identity

    .. MATH::

        d_i d_j = d_{j-1} d_i \ \   \text{if }   i<j.

    For example, `d_0 d_2 = d_1 d_0` and `d_0 d_1 = d_0 d_0`.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import standardize_face_maps
        sage: standardize_face_maps(0, 1)
        (0, 0)
        sage: standardize_face_maps(0, 2)
        (1, 0)
        sage: standardize_face_maps(1, 3, 5)
        (3, 2, 1)
    """
    J = list(L)
    for m in J:
        if Integer(m) < 0:
            raise ValueError('faces are indexed by non-negative integers')
    inadmissible = True
    while inadmissible:
        inadmissible = False
        for idx in range(len(J)-1):
            if J[idx] < J[idx + 1]:
                inadmissible = True
                tmp = J[idx]
                J[idx] = J[idx + 1] - 1
                J[idx + 1] = tmp
    return tuple(J)

def face_degeneracies(m, I):
    r"""
    Result of applying the face map `d_m` to the iterated degeneracy
    `s_I = s_{i_1} s_{i_2} ... s_{i_n}`.

    INPUT:

    - ``m`` -- integer
    - ``I`` -- tuple ``(i_1, i_2, ..., i_n)`` of integers. We assume
      that this sequence is strictly decreasing.

    Using the simplicial identities (see :class:`SimplicialSet`), we
    can rewrite

    .. MATH::

        d_m s_{i_1} s_{i_2} ... s_{i_n}

    in one of the forms

    .. MATH::

        s_{j_1} s_{j_2} ... s_{j_n} d_t, \quad
        s_{j_1} s_{j_2} ... s_{j_{n-1}}.


    OUTPUT: the pair ``(J, t)`` or ``(J, None)``. ``J`` is returned as
    a list.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import face_degeneracies
        sage: face_degeneracies(0, (1, 0))
        ([0], None)
        sage: face_degeneracies(1, (1, 0))
        ([0], None)
        sage: face_degeneracies(2, (1, 0))
        ([0], None)
        sage: face_degeneracies(3, (1, 0))
        ([1, 0], 1)
        sage: face_degeneracies(3, ())
        ([], 3)
    """
    if not I:
        return ([], m)
    J = []
    t = m
    for i in I:
        if t is None:
            J.append(i)
        elif t<i:
            J.append(i-1)
        elif t == i or t == i+1:
            t = None
        else:
            J.append(i)
            t -= 1
    return (J, t)


########################################################################

def shrink_simplicial_complex(K):
    """
    Convert the simplicial complex ``K`` to a "small" simplicial set.

    First convert ``K`` naively, then mod out by a large contractible
    subcomplex, as found by
    :meth:`~sage.homology.simplicial_complex.SimplicialComplex._contractible_subcomplex`.
    This will produce a simplicial set no larger than, and sometimes
    much smaller than, the initial simplicial complex.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import shrink_simplicial_complex
        sage: K = simplicial_complexes.Simplex(3)
        sage: X = shrink_simplicial_complex(K)
        sage: X.f_vector()
        [1]

        sage: Y = simplicial_complexes.Sphere(2)
        sage: S2 = shrink_simplicial_complex(Y)
        sage: S2
        Simplicial set with 2 non-degenerate simplices
        sage: S2.f_vector()
        [1, 0, 1]
        sage: S2.homology()
        {0: 0, 1: 0, 2: Z}

        sage: Z = simplicial_complexes.SurfaceOfGenus(3)
        sage: Z.f_vector()
        [1, 15, 57, 38]
        sage: Z.homology()
        {0: 0, 1: Z^6, 2: Z}
        sage: M = shrink_simplicial_complex(Z)
        sage: M.f_vector()
        [1, 32, 27]
        sage: M.homology()
        {0: 0, 1: Z^6, 2: Z}
    """
    L = K._contractible_subcomplex()
    return SimplicialSet(K).quotient(L)


########################################################################
# Catalog of examples. These are accessed via simplicial_set_catalog.py.

def Sphere(n):
    r"""
    The `n`-sphere as a simplicial set.

    It is constructed with two non-degenerate simplices: a vertex
    `v_0` (which is the base point) and an `n`-simplex `\sigma_n`.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: S0 = simplicial_sets.Sphere(0)
        sage: S0
        S^0
        sage: S0.nondegenerate_simplices()
        [v_0, w_0]
        sage: S0.is_pointed()
        True
        sage: simplicial_sets.Sphere(4)
        S^4
        sage: simplicial_sets.Sphere(4).nondegenerate_simplices()
        [v_0, sigma_4]
    """
    v_0 = AbstractSimplex(0, name='v_0')
    if n == 0:
        w_0 = AbstractSimplex(0, name='w_0')
        return SimplicialSet({v_0: None, w_0: None}, base_point=v_0, name='S^0')
    degens = range(n-2, -1, -1)
    degen_v = v_0.apply_degeneracies(*degens)
    sigma = AbstractSimplex(n, name='sigma_{}'.format(n))
    return SimplicialSet({sigma: [degen_v] * (n+1)}, base_point=v_0, name='S^{}'.format(n))


def ClassifyingSpace(group, n):
    r"""
    The `n`-skeleton of the classifying space of ``group``,
    as a simplicial set.

    INPUT:

    - ``group`` -- a finite group or finite monoid
    - ``n`` -- a non-negative integer

    See
    :meth:`~sage.categories.finite_monoids.FiniteMonoids.ParentMethods.nerve_n_skeleton`
    for more details and more examples.

    EXAMPLES::

        sage: C2 = groups.misc.MultiplicativeAbelian([2])
        sage: BC2 = simplicial_sets.ClassifyingSpace(C2, 8)
        sage: H = BC2.homology(base_ring=GF(2))
        sage: [H[i].dimension() for i in range(9)]
        [0, 1, 1, 1, 1, 1, 1, 1, 1]

        sage: Klein4 = groups.misc.MultiplicativeAbelian([2, 2])
        sage: BK = simplicial_sets.ClassifyingSpace(Klein4, 5)
        sage: BK
        5-skeleton of classifying space of Multiplicative Abelian group isomorphic to C2 x C2
        sage: BK.homology(base_ring=GF(2))
        {0: Vector space of dimension 0 over Finite Field of size 2,
         1: Vector space of dimension 2 over Finite Field of size 2,
         2: Vector space of dimension 3 over Finite Field of size 2,
         3: Vector space of dimension 4 over Finite Field of size 2,
         4: Vector space of dimension 5 over Finite Field of size 2,
         5: Vector space of dimension 185 over Finite Field of size 2}

    The 5-dimensional homology differs from that for `C_2 \times C_2`
    because this is the 5-skeleton, not the whole classifying space.
    """
    X = group.nerve_n_skeleton(n)
    X.rename('{}-skeleton of classifying space of {}'.format(n, group))
    return X


def RealProjectiveSpace(n):
    """
    Real `n`-dimensional projective space, as a simplicial set.

    This is constructed as the `n`-skeleton of the nerve of the group
    of order 2, and therefore has a single non-degenerate simplex in
    each dimension up to `n`.

    EXAMPLES::

        sage: simplicial_sets.RealProjectiveSpace(7)
        RP^7
        sage: RP5 = simplicial_sets.RealProjectiveSpace(5)
        sage: RP5.homology()
        {0: 0, 1: C2, 2: 0, 3: C2, 4: 0, 5: Z}
    """
    X = AbelianGroup([2]).nerve_n_skeleton(n)
    X.rename('RP^{}'.format(n))
    return X


def KleinBottle():
    r"""
    The Klein bottle as a simplicial set.

    This converts the `\Delta`-complex version to a simplicial set. It
    has one 0-simplex, three 1-simplices, and two 2-simplices.

    EXAMPLES::

        sage: K = simplicial_sets.KleinBottle()
        sage: K.f_vector()
        [1, 3, 2]
        sage: K.homology(reduced=False)
        {0: Z, 1: Z x C2, 2: 0}
    """
    temp = SimplicialSet(delta_complexes.KleinBottle())
    pt = temp.n_cells(0)[0]
    return SimplicialSet(temp.face_data(), base_point=pt)


def Torus():
    r"""
    The torus as a simplicial set.

    This computes the product of the circle with itself, where the
    circle is represented using a single 0-simplex and a single
    1-simplex. Thus it has one 0-simplex, three 1-simplices, and two
    2-simplices.

    EXAMPLES::

        sage: T = simplicial_sets.Torus()
        sage: T.f_vector()
        [1, 3, 2]
        sage: T.homology(reduced=False)
        {0: Z, 1: Z x Z, 2: Z}
    """
    S1 = Sphere(1)
    T = S1.product(S1)
    T.rename('Torus')
    return T


def Simplex(n):
    r"""
    The `n`-simplex as a simplicial set.

    EXAMPLES::

        sage: K = simplicial_sets.Simplex(2)
        sage: K
        2-simplex
        sage: K.n_cells(0)
        [(0,), (1,), (2,)]
        sage: K.n_cells(1)
        [(0, 1), (0, 2), (1, 2)]
        sage: K.n_cells(2)
        [(0, 1, 2)]
    """
    return SimplicialSet(simplicial_complexes.Simplex(n), name='{}-simplex'.format(n))


@cached_function
def Empty():
    """
    The empty simplicial set.

    This should return the same simplicial set each time it is called.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import Empty
        sage: E = Empty()
        sage: E
        Empty simplicial set
        sage: E.nondegenerate_simplices()
        []
        sage: E is Empty()
        True
    """
    return SimplicialSet({}, name='Empty simplicial set')


@cached_function
def Point():
    """
    A single point called "*" as a simplicial set.

    This should return the same simplicial set each time it is called.

    EXAMPLES::

        sage: P = simplicial_sets.Point()
        sage: P.is_pointed()
        True
        sage: P.nondegenerate_simplices()
        [*]

        sage: Q = simplicial_sets.Point()
        sage: P is Q
        True
        sage: P == Q
        True
    """
    star = AbstractSimplex(0, name='*')
    return SimplicialSet({star: None}, base_point=star, name='Point')


def Horn(n, k):
    r"""
    The horn $\Lambda^n_k$.

    This is the subsimplicial set of the $n$-simplex obtained by
    removing its $k$-th face.

    EXAMPLES::

        sage: L = simplicial_sets.Horn(3, 0)
        sage: L.n_cells(3)
        []
        sage: L.n_cells(2)
        [(0, 1, 2), (0, 1, 3), (0, 2, 3)]
    """
    K = Simplex(n)
    sigma = K.n_cells(n)[0]
    return K.subsimplicial_set(K.faces(sigma)[:k] + K.faces(sigma)[k+1:])


def ComplexProjectiveSpace(n):
    r"""
    Complex `n`-dimensional projective space, as a simplicial set.

    This is only defined when `n` is at most 4. It is constructed
    using the simplicial set decomposition provided by Kenzo, as
    described by Sergeraert [Ser]_

    REFERENCES:

    .. [Ser] \F. Sergeraert, *Triangulations of complex projective
       spaces* in Scientific contributions in honor of Mirian
       Andrés Gómez, pp 507-519, Univ. La Rioja Serv. Publ., Logroño (2010).

    EXAMPLES::

        sage: simplicial_sets.ComplexProjectiveSpace(2).homology(reduced=False)
        {0: Z, 1: 0, 2: Z, 3: 0, 4: Z}
        sage: CP3 = simplicial_sets.ComplexProjectiveSpace(3)
        sage: CP3
        CP^3
        sage: CP3.f_vector()
        [1, 0, 3, 10, 25, 30, 15]

        sage: K = CP3.suspension()
        sage: R = K.cohomology_ring(GF(2))
        sage: R.gens()
        (h^{0,0}, h^{3,0}, h^{5,0}, h^{7,0})
        sage: x = R.gens()[1]
        sage: x.Sq(2)
        h^{5,0}

        sage: simplicial_sets.ComplexProjectiveSpace(4).f_vector()
        [1, 0, 4, 22, 97, 255, 390, 315, 105]

        sage: simplicial_sets.ComplexProjectiveSpace(5)
        Traceback (most recent call last):
        ...
        ValueError: complex projective spaces are only available in dimensions between 0 and 4
    """
    if n < 0 or n > 4:
        raise ValueError('complex projective spaces are only available in dimensions between 0 and 4')
    if n == 0:
        return Point()
    if n == 1:
        return Sphere(2)
    if n == 2:
        v = AbstractSimplex(0, name='<<GBar>>')
        f2_1 = AbstractSimplex(2, name='<<GBar<- (1)><- NIL>>>')
        f2_2 = AbstractSimplex(2, name='<<GBar<- (2)><- NIL>>>')
        f3_110 = AbstractSimplex(3, name='<<GBar<- (1 1)><0 NIL><- NIL>>>')
        f3_011 = AbstractSimplex(3, name='<<GBar<0 (1)><- (1)><- NIL>>>')
        f3_111 = AbstractSimplex(3, name='<<GBar<1 (1)><- (1)><- NIL>>>')
        f4_101101 = AbstractSimplex(4, name='<<GBar<1-0 (1)><1-0 NIL><- (1)><- NIL>>>')
        f4_201110 = AbstractSimplex(4, name='<<GBar<2-0 (1)><1 (1)><0 NIL><- NIL>>>')
        f4_211010 = AbstractSimplex(4, name='<<GBar<2-1 (1)><0 (1)><0 NIL><- NIL>>>')
        return SimplicialSet({f2_1: (v.apply_degeneracies(0),
                                     v.apply_degeneracies(0),
                                     v.apply_degeneracies(0)),
                              f2_2: (v.apply_degeneracies(0),
                                     v.apply_degeneracies(0),
                                     v.apply_degeneracies(0)),
                              f3_110: (f2_1, f2_2, f2_1, v.apply_degeneracies(1, 0)),
                              f3_011: (f2_1, f2_1, f2_1, f2_1),
                              f3_111: (v.apply_degeneracies(1, 0), f2_1, f2_2, f2_1),
                              f4_101101: (f2_1.apply_degeneracies(0),
                                          f2_1.apply_degeneracies(0),
                                          f3_011,
                                          f2_1.apply_degeneracies(2),
                                          f2_1.apply_degeneracies(2)),
                              f4_201110: (f2_1.apply_degeneracies(1),
                                          f3_111,
                                          f3_011,
                                          f3_110,
                                          f2_1.apply_degeneracies(1)),
                              f4_211010: (f2_1.apply_degeneracies(2),
                                          f3_111,
                                          f2_1.apply_degeneracies(1),
                                          f3_110,
                                          f2_1.apply_degeneracies(0))},
                             name='CP^2', base_point=v)
    if n == 3:
        file = os.path.join(SAGE_ENV['SAGE_SRC'], 'ext', 'kenzo', 'CP3.txt')
        data = simplicial_data_from_kenzo_output(file)
        v = [_ for _ in data.keys() if _.dimension() == 0][0]
        return SimplicialSet(data, name='CP^3', base_point=v)
    if n == 4:
        file = os.path.join(SAGE_ENV['SAGE_SRC'], 'ext', 'kenzo', 'CP4.txt')
        data = simplicial_data_from_kenzo_output(file)
        v = [_ for _ in data.keys() if _.dimension() == 0][0]
        return SimplicialSet(data, name='CP^4', base_point=v)


def simplicial_data_from_kenzo_output(filename):
    """
    INPUT:

    - ``filename`` -- name of file containing the output from Kenzo's
      ``show-structure`` function

    OUTPUT: data to construct a simplicial set from the Kenzo output

    Several files with Kenzo output are in the directory
    ``SAGE_ROOT/src/ext/kenzo/``.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import simplicial_data_from_kenzo_output, SimplicialSet
        sage: sphere = os.path.join(SAGE_ENV['SAGE_SRC'], 'ext', 'kenzo', 'S4.txt')
        sage: S4 = SimplicialSet(simplicial_data_from_kenzo_output(sphere))
        sage: S4.homology(reduced=False)
        {0: Z, 1: 0, 2: 0, 3: 0, 4: Z}
    """
    with open(filename, 'r') as f:
        data = f.read()
    dim = 0
    start = 0
    # simplex_data: data for constructing the simplicial set.
    simplex_data = {}
    # simplex_names: simplices indexed by their names
    simplex_names = {}
    dim_idx = data.find('Dimension = {}:'.format(dim), start)
    while dim_idx != -1:
        start = dim_idx + len('Dimension = {}:'.format(dim))
        new_dim_idx = data.find('Dimension = {}:'.format(dim+1), start)
        if new_dim_idx == -1:
            end = len(data)
        else:
            end = new_dim_idx
        if dim == 0:
            simplex_string = data[data.find('Vertices :') + len('Vertices :'):end]
            vertices = OneOrMore(nestedExpr()).parseString(simplex_string).asList()[0]
            for v in vertices:
                vertex = AbstractSimplex(0, name=v)
                simplex_data[vertex] = None
                simplex_names[v] = vertex
        else:
            simplex_string = data[start:end].strip()

            for s in [_.strip() for _ in simplex_string.split('Simplex : ')]:
                if s:
                    name, face_str = [_.strip() for _ in s.split('Faces : ')]
                    face_str = face_str.strip('()')
                    face_str = face_str.split('<AbSm ')
                    faces = []
                    for f in face_str[1:]:
                        # f has the form 'DEGENS NAME>', possibly with a trailing space.

                        # DEGENS is a hyphen-separated list, like
                        # '3-2-1-0' or '0' or '-'.
                        m = re.match('[-[0-9]+', f)
                        degen_str = m.group(0)
                        if degen_str.find('-') != -1:
                            if degen_str == '-':
                                degens = []
                            else:
                                degens = [Integer(_) for _ in degen_str.split('-')]
                        else:
                            degens = [Integer(degen_str)]

                        face_name = f[m.end(0):].strip()[:-1]
                        nondegen = simplex_names[face_name]
                        faces.append(nondegen.apply_degeneracies(*degens))

                    simplex = AbstractSimplex(dim, name=name)
                    simplex_data[simplex] = faces
                    simplex_names[name] = simplex
        dim += 1
        dim_idx = new_dim_idx
    return simplex_data
