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

.. math::

    d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
    s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

satisfying the *simplicial identities*:

.. math::

    d_i d_j &= d_{j-1} d_i \ \  \text{if } i<j \\
    d_i s_j &= s_{j-1} d_i \ \  \text{if } i<j \\
    d_j s_j &= 1 = d_{j+1} s_j \\
    d_i s_j &= s_{j} d_{i-1} \ \  \text{if } i>j+1 \\
    s_i s_j &= s_{j+1} s_{i} \ \  \text{if } i<j+1

See :wikipedia:`Simplicial_set`, Peter May's seminal book [May]_, or
Greg Friedman's "Illustrated introduction" :arxiv:`0809.4221` for more
information.

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
from sage.misc.cachefunc import cached_method
from pyparsing import OneOrMore, nestedExpr
from sage.graphs.graph import Graph
from sage.env import SAGE_ENV
import copy
import itertools
import re
import os

from sage.misc.lazy_import import lazy_import
lazy_import('sage.categories.simplicial_sets', 'SimplicialSets')

class AbstractSimplex(SageObject):
    """
    An abstract simplex, a building block of a simplicial set.

    In a simplicial set, a simplex either is non-degenerate or is
    obtained by applying degeneracy maps to a non-degenerate simplex.

    INPUT:

    - ``dim`` -- a non-negative integer, the dimension of the
      underlying non-degenerate simplex.

    - ``degeneracies`` -- a list or tuple of non-negative integers,
      the degeneracies to be applied.

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

        sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
        sage: AbstractSimplex(3, (3, 1))
        Simplex obtained by applying degeneracies s_3 s_1 to Non-degenerate simplex of dimension 3
        sage: AbstractSimplex(3, None)
        Non-degenerate simplex of dimension 3

    Simplices may be named (or renamed), affecting how they are printed::

        sage: NonDegenerateSimplex(0)
        Non-degenerate simplex of dimension 0
        sage: v = NonDegenerateSimplex(0, name='v')
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

        sage: v = NonDegenerateSimplex(0, name='v')
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
    def __init__(self, dim, degeneracies, **kwds):
        """
        A simplex of dimension ``dim``.

        INPUT:

        - ``dim`` -- integer, the dimension
        - ``degeneracies`` -- iterable, the indices of the degeneracy maps
        - ``underlying`` (optional) -- a non-degenerate simplex
        - ``name`` (optional) -- string

        TESTS::

            sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
            sage: AbstractSimplex(3, (1,2))
            Simplex obtained by applying degeneracies s_3 s_1 to Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, None)
            Non-degenerate simplex of dimension 3

            sage: v = NonDegenerateSimplex(0, name='v')
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
            self._underlying = kwds.get('underlying', AbstractSimplex(dim, None))
        else:
            self._degens = ()
            self._underlying = kwds.get('underlying', self)
        if 'name' in kwds:
            self.rename(kwds['name'])
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, AbstractSimplex
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)

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

            sage: NonDegenerateSimplex(0) < NonDegenerateSimplex(3)
            True
            sage: AbstractSimplex(0, ((0,0,0))) <= NonDegenerateSimplex(2)
            False

        Degenerate comes after non-degenerate, and if both are
        degenerate, sort on the degeneracies::

            sage: AbstractSimplex(0, ((0,0))) <= NonDegenerateSimplex(2)
            False
            sage: AbstractSimplex(1, ((0,))) < AbstractSimplex(1, ((1,)))
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False
            sage: w.rename('a')
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            True
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

    def __le__(self, other):
        """
        True if ``self`` is less than or equal to ``other``.

        See :meth:`__lt__` for details on the ordering.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, AbstractSimplex
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: v <= v
            True
            sage: w <= v
            False
            sage: v.apply_degeneracies(1,0) <= w.apply_degeneracies(1,0)
            True
        """
        return self == other or self < other

    def __gt__(self, other):
        """
        True if ``self`` is greater than ``other``.

        See :meth:`__lt__` for details on the ordering.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, AbstractSimplex
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: v > v
            False
            sage: w > v
            True
            sage: v.apply_degeneracies(1,0) > w.apply_degeneracies(1,0)
            False
        """
        return not self <= other

    def __ge__(self, other):
        """
        True if ``self`` is greater than or equal to ``other``.

        See :meth:`__lt__` for details on the ordering.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, AbstractSimplex
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: v >= v
            True
            sage: w >= v
            True
            sage: v.apply_degeneracies(1,0) >= w.apply_degeneracies(1,0)
            False
        """
        return not self < other

    def nondegenerate(self):
        """
        The non-degenerate simplex underlying this one.

        Therefore return itself if this simplex is non-degenerate.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, AbstractSimplex
            sage: v = NonDegenerateSimplex(0, name='v')
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

            sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
            sage: AbstractSimplex(3, (2,1)).is_nondegenerate()
            False
            sage: AbstractSimplex(3, None).is_nondegenerate()
            True
            sage: NonDegenerateSimplex(5).is_nondegenerate()
            True
        """
        return not self.is_degenerate()

    def dimension(self):
        """
        The dimension of this simplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
            sage: AbstractSimplex(3, (2,1)).dimension()
            5
            sage: AbstractSimplex(3, None).dimension()
            3
            sage: NonDegenerateSimplex(7).dimension()
            7
        """
        return self._dim + len(self.degeneracies())


    def apply_degeneracies(self, *args):
        """
        Apply the degeneracies given by the arguments ``args`` to this simplex.

        INPUT:

        - ``args`` -- integers

        EXAMPLES::

            sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
            sage: v = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: v = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import AbstractSimplex, NonDegenerateSimplex
            sage: AbstractSimplex(3, None)
            Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, (0,))
            Simplex obtained by applying degeneracy s_0 to Non-degenerate simplex of dimension 3
            sage: AbstractSimplex(3, (0, 0))
            Simplex obtained by applying degeneracies s_1 s_0 to Non-degenerate simplex of dimension 3

        Test renaming::

            sage: v = NonDegenerateSimplex(0)
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


class NonDegenerateSimplex(AbstractSimplex):
    """
    A non-degenerate simplex of dimension ``dim``.

    INPUT:

    - ``dim`` -- a non-negative integer, the dimension
    - ``name`` (optional) -- string, a name for the simplex

    EXAMPLES::

        sage: from sage.homology.simplicial_set import NonDegenerateSimplex
        sage: e = NonDegenerateSimplex(1, name='e')
        sage: e
        e
        sage: e.rename('e_1')
        sage: e
        e_1

    You can apply degeneracy maps to a non-degenerate simplex to get a
    degenerate simplex::

        sage: sigma = e.apply_degeneracies(2,0); sigma
        Simplex obtained by applying degeneracies s_2 s_0 to e_1
        sage: sigma.dimension()
        3
    """
    def __init__(self, dim, name=None):
        """
        A non-degenerate simplex of dimension ``dim``.

        This is just a shortcut to a special case of :class:`AbstractSimplex`.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: NonDegenerateSimplex(0)
            Non-degenerate simplex of dimension 0

            sage: NonDegenerateSimplex(pi)
            Traceback (most recent call last):
            ...
            ValueError: the dimension must be an integer
        """
        AbstractSimplex.__init__(self, dim, degeneracies=None, name=name)


class SimplicialSet(GenericCellComplex, Parent):
    r"""
    A simplicial set.

    A simplicial set `X` is a collection of sets `X_n`, the
    *n-simplices*, indexed by the non-negative integers, together with
    maps

    .. math::

        d_i: X_n \to X_{n-1}, \ \ 0 \leq i \leq n \ \  \text{(face maps)} \\
        s_j: X_n \to X_{n+1}, \ \ 0 \leq j \leq n \ \  \text{(degeneracy maps)}

    satisfying the *simplicial identities*:

    .. math::

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

        sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
        sage: v = NonDegenerateSimplex(0)
        sage: e = NonDegenerateSimplex(1)
        sage: S1 = SimplicialSet({e: (v, v)})  # the circle
        sage: S1
        Simplicial set with 2 non-degenerate simplices

    Simplicial sets may be given names::

        sage: SimplicialSet({e: (v, v)}, name='Simplicial circle')
        Simplicial circle

        sage: s_0_v = v.apply_degeneracies(0) # s_0 applied to v
        sage: sigma = NonDegenerateSimplex(2)
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

        sage: v = NonDegenerateSimplex(0)
        sage: w = NonDegenerateSimplex(0)
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
        sage: S0.base_point()   # No base point, so no output
        sage: S0.set_base_point(v)
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
        [(1, 2), (0, 2), (0, 1)]
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
    def __init__(self, data, base_point=None, name=None, check=True):
        r"""
        See :class:`SimplicialSet` for documentation.

        TESTS::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: e = NonDegenerateSimplex(1)
            sage: SimplicialSet({e: (v, v, v)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1
            sage: SimplicialSet({e: (v,)})
            Traceback (most recent call last):
            ...
            ValueError: wrong number of faces for simplex in dimension 1

        Base points::

            sage: SimplicialSet({e: (v,v)}, base_point=NonDegenerateSimplex(0))
            Traceback (most recent call last):
            ...
            ValueError: the base point is not a simplex in this simplicial set
            sage: SimplicialSet({e: (v,v)}, base_point=e)
            Traceback (most recent call last):
            ...
            ValueError: the base "point" is not a zero-simplex

        Simplicial identity::

            sage: sigma = NonDegenerateSimplex(2)
            sage: w = NonDegenerateSimplex(0)
            sage: K = SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0),
            ....:                            v.apply_degeneracies(0))})
            sage: SimplicialSet({sigma: (v.apply_degeneracies(0),
            ....:                        v.apply_degeneracies(0),
            ....:                        w.apply_degeneracies(0))})
            Traceback (most recent call last):
            ...
            ValueError: simplicial identity d_i d_j = d_{j-1} d_i fails in dimension 2
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
                        new_sigma = NonDegenerateSimplex(d)
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
                        new_sigma = NonDegenerateSimplex(d)
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
                return copy.copy(self)
            else:
                raise NotImplementedError('I do not know how to convert this to a simplicial set')
        # Convert each value in data to a tuple, and then convert all
        # of data to a tuple, so that it is hashable.
        for x in data:
            if data[x]:
                if x.dimension() != len(data[x]) - 1:
                    raise ValueError('wrong number of faces for simplex in dimension {}'.format(x.dimension()))
                if not all(y.dimension() == x.dimension() - 1 for y in data[x]):
                    raise ValueError('faces of a {}-simplex have the wrong dimension'.format(x.dimension()))
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
                                raise ValueError('simplicial identity d_i d_j = d_{{j-1}} d_i'
                                                 ' fails in dimension {}'.format(d))

        # Now define the attributes for an instance of this class.
        # self._data: a tuple representing the defining data of the
        # simplicial set.
        self._data = tuple(data.items())
        # self._simplices: a sorted tuple of non-degenerate simplices.
        self._simplices = sorted(tuple(simplices))
        # self._basepoint: the base point, or None.
        if base_point is not None:
            if base_point not in simplices:
                raise ValueError('the base point is not a simplex in this simplicial set')
            if base_point.dimension() != 0:
                raise ValueError('the base "point" is not a zero-simplex')
        else:
            # base_point is None. If there is only one vertex, make
            # that the base point.
            if (len(self._simplices) == 1
                or (len(self._simplices) > 1 and self._simplices[1].dimension() > 1)):
                base_point = self._simplices[0]
        self._basepoint = base_point
        Parent.__init__(self, category=SimplicialSets().Finite())
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
            sage: e = NonDegenerateSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({e: (w, w)})
            sage: X == X
            True
            sage: X == SimplicialSet({e: (v, v)})
            True
            sage: X == Y
            False
        """
        return isinstance(other, SimplicialSet) and sorted(self._data) == sorted(other._data)

    def __ne__(self, other):
        """
        True if ``self`` and ``other`` are not equal as simplicial sets.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: X = SimplicialSet({v: None})
            sage: degen = v.apply_degeneracies(0)
            sage: tau = NonDegenerateSimplex(2)
            sage: Y = SimplicialSet({tau: (degen, degen, degen)})

            sage: hash(X) # random
            17
            sage: hash(X) != hash(Y)
            True
        """
        return hash(self._data)

    def face_data(self):
        """
        The face-map data defining this simplicial set.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex
            sage: w = NonDegenerateSimplex(0)
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

        .. note::

            The sorting is done when the simplicial set is
            constructed, so changing the name of a simplex after
            construction will not affect the ordering.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
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

        - ``subcomplex`` (optional) -- a subcomplex of this simplicial
          set

        Each key is a dimension, and the corresponding value is the
        list of simplices in that dimension. If ``subcomplex`` is
        specified, then the corresponding value is the list of
        simplices in the quotient by the subcomplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.cells()
            {0: [Non-degenerate simplex of dimension 0, Non-degenerate simplex of dimension 0]}

            sage: v.rename('v')
            sage: w.rename('w')
            sage: S0.cells()
            {0: [v, w]}

            sage: e = NonDegenerateSimplex(1, name='e')
            sage: S1 = SimplicialSet({e: (v, v)})
            sage: S1.cells()
            {0: [v], 1: [e]}

            sage: S0.cells(S0.subcomplex([v, w]))
            {0: [*]}

            sage: X = SimplicialSet({e: (v,w)})
            sage: X.cells(X.subcomplex([v, w]))
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.f_vector()
            [2]

            sage: e = NonDegenerateSimplex(1)
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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = NonDegenerateSimplex(2, name='tau')
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

    def base_point(self):
        """
        Return this simplicial set's base point, or return ``None`` if no
        base point has been defined.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='*')
            sage: e = NonDegenerateSimplex(1)
            sage: S1 = SimplicialSet({e: (v, v)})  # the circle
            sage: S1.base_point()  # no base point, so no output
            sage: S1 = SimplicialSet({e: (v, v)}, base_point=v)
            sage: S1.base_point()
            *
        """
        return self._basepoint

    def is_pointed(self):
        """
        True if this simplicial set is pointed, i.e., has a base point.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
            sage: e = NonDegenerateSimplex(1)
            sage: X = SimplicialSet({e: (v, w)})
            sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
            sage: X.is_pointed()
            False
            sage: Y.is_pointed()
            True
        """
        return self._basepoint is not None

    def set_base_point(self, point):
        """
        Define ``point`` to be the base point of this simplicial set.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v_0')
            sage: w = NonDegenerateSimplex(0, name='w_0')
            sage: e = NonDegenerateSimplex(1)
            sage: X = SimplicialSet({e: (v, w)})
            sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
            sage: Y.base_point()
            w_0
            sage: X.set_base_point(w)
            sage: X.base_point()
            w_0
            sage: Y.set_base_point(v)
            sage: Y.base_point()
            v_0

        TESTS::

            sage: X.set_base_point(e)
            Traceback (most recent call last):
            ...
            ValueError: the "point" is not a zero-simplex
            sage: pt = NonDegenerateSimplex(0)
            sage: X.set_base_point(pt)
            Traceback (most recent call last):
            ...
            ValueError: the point is not a simplex in this simplicial set
        """
        if point is not None:
            if point.dimension() != 0:
                raise ValueError('the "point" is not a zero-simplex')
            if point not in self._simplices:
                raise ValueError('the point is not a simplex in this simplicial set')
        self._basepoint = point

    def unset_base_point(self):
        """
        Unset the base point for this simplicial set.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v_0')
            sage: w = NonDegenerateSimplex(0, name='w_0')
            sage: e = NonDegenerateSimplex(1)
            sage: Y = SimplicialSet({e: (v, w)}, base_point=w)
            sage: Y.is_pointed()
            True
            sage: Y.base_point()
            w_0
            sage: Y.unset_base_point()
            sage: Y.is_pointed()
            False
            sage: Y.unset_base_point()
            sage: Y.is_pointed()
            False
        """
        self._basepoint = None

    def n_skeleton(self, n):
        """
        The `n`-skeleton of this simplicial set, obtained by removing all
        non-degenerate simplices of dimension larger than `n`.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: degen = v.apply_degeneracies(0)
            sage: tau = NonDegenerateSimplex(2, name='tau')
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
        vertices = self.n_cells(0)
        data = [_ for _ in self._data if _[0].dimension() <= n]
        for v in vertices:
            if all(_[0] != v for _ in data):
                data.append((v, None))
        return SimplicialSet(dict(data))

    def chain_complex(self, **kwds):
        r"""
        The normalized chain complex.

        INPUT:

        - ``max_dimension`` -- if ``None`` (the default), compute the
          chain complex up to the maximum dimension of the
          non-degenerate simplices. Otherwise, compute the chain
          complex through this dimension.

        - ``base_ring`` (optional, default ``ZZ``) -- commutative ring

        - ``augmented`` (optional, default ``False``) -- if ``True``,
          return the augmented chain complex (that is, include a class
          in dimension `-1` corresponding to the empty cell).

        - ``cochain`` (optional, default ``False``) -- if ``True``,
          return the cochain complex (that is, the dual of the chain
          complex).

        The normalized chain complex of a simplicial set is isomorphic
        to the chain complex obtained by modding out by degenerate
        simplices, and the latter is what is actually constructed
        here.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: degen = v.apply_degeneracies(1, 0) # s_1 s_0 applied to v
            sage: sigma = NonDegenerateSimplex(3)
            sage: S3 = SimplicialSet({sigma: (degen, degen, degen, degen)}) # the 3-sphere
            sage: S3.chain_complex().homology()
            {0: Z, 3: Z}
            sage: S3.chain_complex(augmented=True).homology()
            {-1: 0, 0: 0, 3: Z}
            sage: S3.chain_complex(max_dimension=2, base_ring=QQ).homology()
            {0: Vector space of dimension 1 over Rational Field}

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
        augmented = kwds.get('augmented', False)
        cochain = kwds.get('cochain', False)
        base_ring = kwds.get('base_ring', ZZ)
        max_dim = kwds.get('max_dimension', None)

        if not max_dim:
            max_dim = self.dimension()

        differentials = {}
        # Convert the tuple self._data to a dictionary indexed by the
        # non-degenerate simplices.
        face_data = self.face_data()
        # simplices: dictionary indexed by dimension, values the list
        # of non-degenerate simplices in that dimension.
        simplices = {}
        for sigma in self.nondegenerate_simplices():
            if sigma.dimension() in simplices:
                simplices[sigma.dimension()].append(sigma)
            else:
                simplices[sigma.dimension()] = [sigma]
        rank = len(simplices[0])
        if augmented:
            differentials[-1] = matrix(base_ring, 0, 1)
            differentials[0] = matrix(base_ring, 1, rank,
                                      [1] * rank)
        else:
            differentials[0] = matrix(base_ring, 0, rank)
        current = sorted(simplices[0])

        for d in range(1, max_dim + 1):
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
            return ChainComplex(new_diffs, degree_of_differential=1)

        try:
            return ChainComplex(differentials, degree_of_differential=-1)
        except ValueError:
            return differentials

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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: e = NonDegenerateSimplex(1, name='e')
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
            sage: X.set_base_point(v)
            sage: Y = copy(X)
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

        base_point=translate.get(self.base_point(), None)
        C = SimplicialSet(new_data, base_point=base_point)
        return C

    def subcomplex(self, simplices):
        """
        The subcomplex of this simplicial set determined by ``simplices``,
        a set of nondegenerate simplices.

        INPUT:

        - ``simplices`` -- set, list, or tuple of nondegenerate
          simplices in this simplicial set, or a simplicial 
          complex -- see below.

        If ``simplices`` is a simplicial complex, then the original
        simplicial set should itself have been converted from a
        simplicial complex, and ``simplices`` should be a subcomplex
        of that.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: f = NonDegenerateSimplex(1, name='f')

            sage: X = SimplicialSet({e: (v, w), f: (w, v)})
            sage: Y = X.subcomplex([e])
            sage: Y
            Simplicial set with 3 non-degenerate simplices
            sage: Y.nondegenerate_simplices()
            [v, w, e]

            sage: S3 = simplicial_complexes.Sphere(3)
            sage: K = SimplicialSet(S3)
            sage: tau = K.n_cells(3)[0]
            sage: tau.dimension()
            3
            sage: K.subcomplex([tau])
            Simplicial set with 15 non-degenerate simplices

        TESTS:

        Make sure vertices are treated properly::

            sage: X.subcomplex([v]).nondegenerate_simplices()
            [v]
            sage: X.subcomplex([v, w]).nondegenerate_simplices()
            [v, w]
            sage: S0 = SimplicialSet({v: None, w: None})
            sage: S0.subcomplex([w]).nondegenerate_simplices()
            [w]

        Raise an error if an element of ``simplices`` is not actually
        in the original simplicial set::

            sage: sigma = NonDegenerateSimplex(2, name='sigma_2')
            sage: Z = X.subcomplex([e, sigma])
            Traceback (most recent call last):
            ...
            ValueError: not all simplices are in the original simplicial set

        Simplicial complexes::

            sage: X = simplicial_complexes.ComplexProjectivePlane()
            sage: Y = X._contractible_subcomplex()
            sage: CP2 = SimplicialSet(X)
            sage: sub = CP2.subcomplex(Y)
            sage: CP2.f_vector()
            [9, 36, 84, 90, 36]
            sage: CP2.quotient(sub).f_vector()
            [1, 0, 23, 45, 24]
            sage: CP2.quotient(sub).homology()
            {0: 0, 1: 0, 2: Z, 3: 0, 4: Z}

        Try to construct a subcomplex from a simplicial complex which
        is not actually contained in ``self``::

            sage: Z = SimplicialComplex([[0,1,2,3,4]])
            sage: CP2.subcomplex(Z)
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
                if x not in data and x.dimension() > 0:
                    raise ValueError('not all simplices are in the original simplicial set')
                keep.add(x)
                if x in data and data[x]:
                    keep.update([f.nondegenerate() for f in data[x]])
                else:
                    # x is a vertex
                    assert(x.dimension() == 0)
                    vertices.add(x)
        missing = set(self.nondegenerate_simplices()).difference(keep)
        for x in missing:
            if x in data:
                del data[x]
        for x in vertices:
            data[x] = None
        return SimplicialSet(data)

    def is_subcomplex(self, other):
        """
        Return ``True`` iff ``self`` is a subcomplex of ``other``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: S1 = simplicial_sets.Sphere(1)
            sage: pt = S1.n_cells(0)[0]
            sage: other_pt = NonDegenerateSimplex(0)

            sage: SimplicialSet({pt: None}).is_subcomplex(S1)
            True
            sage: SimplicialSet({other_pt: None}).is_subcomplex(S1)
            False
            sage: S1.is_subcomplex(S1.wedge(S1))
            True
        """
        return all(x in other._simplices for x in self.nondegenerate_simplices())

    def quotient(self, subcomplex, vertex_name=None):
        """
        The quotient of this simplicial set by ``subcomplex``.

        That is, ``subcomplex`` is replaced by a vertex.

        INPUT:

        - ``subcomplex`` -- subcomplex of this simplicial set, or a
          list, tuple, or set of simplices defining a subcomplex.

        - ``vertex_name`` (optional) -- string, name to be given to the new
          vertex. If omitted, use '*'.

        Base points: if the original simplicial set has a base point
        contained in ``subcomplex``, then ``*`` is the base point in
        the quotient. If the original simplicial set has a base point
        not contained in ``subcomplex``, then use that base point. If
        the original simplicial set had no base point, then use
        ``*``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: f = NonDegenerateSimplex(1, name='f')
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
        # v: new vertex, to replace ``subcomplex``.
        if vertex_name:
            v = NonDegenerateSimplex(0, name=vertex_name)
        else:
            v = NonDegenerateSimplex(0, name='*')
        if not isinstance(subcomplex, SimplicialSet):
            # If it's not a simplicial set, subcomplex should be a
            # list, tuple, or set of simplices, so form the actual
            # subcomplex:
            subcomplex = self.subcomplex(subcomplex)
        else:
            # Test whether subcomplex is actually a subcomplex of
            # self.
            if not subcomplex.is_subcomplex(self):
                raise ValueError('the "subcomplex" is not actually a subcomplex')
            pass
        removed = subcomplex.nondegenerate_simplices()
        data = self.face_data()
        # Remove the simplices in the subcomplex as keys:
        for x in removed:
            if x in data:
                del data[x]
        # Remove the simplices in the subcomplex as faces, replacing
        # them with v:
        for y in data:
            if data[y]:
                d = y.dimension()
                faces = []
                for x in data[y]:
                    if x.nondegenerate() in removed:
                        faces.append(v.apply_degeneracies(*range(d-2, -1, -1)))
                    else:
                        faces.append(x)
                data[y] = faces

        if self.is_pointed():
            if self.base_point() in removed:
                base_pt = v
            else:
                base_pt = self.base_point()
        else:
            base_pt = v
        data[v] = None
        return SimplicialSet(data, base_point = base_pt)

    def disjoint_union(self, other):
        """
        The disjoint union of this simplicial set with ``other``.

        If there are any simplices in common between the two
        simplicial sets, use a copy of ``other``.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: f = NonDegenerateSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)})
            sage: Y = SimplicialSet({f: (v, w)})
            sage: Z = X.disjoint_union(Y)

        Since ``X`` and ``Y`` have simplices in common, ``Y`` has been
        copied and so its simplices have been renamed::

            sage: Z.nondegenerate_simplices()
            [v, v', w', e, f']

        Test that copying works properly::

            sage: Y = copy(X)
            sage: Y.nondegenerate_simplices()
            [v', e']
            sage: Z = Y.disjoint_union(Y)
            sage: Z.nondegenerate_simplices()
            [v', v'', e', e'']
        """
        X = other
        simps = set(self.nondegenerate_simplices())
        while not simps.isdisjoint(set(X.nondegenerate_simplices())):
            X = copy.copy(X)
        data = self.face_data()
        data.update(X._data)
        return SimplicialSet(data)

    def wedge(self, other):
        r"""
        The wedge product of this simplicial set with ``other``.

        If there are any simplices in common between the two
        simplicial sets, use a copy of ``other``, unless of course the
        base point is the only simplex in common.

        Both this simplicial set and ``other`` must be pointed. Their
        base points need not be the same, but during the construction,
        the base points will become identified.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: f = NonDegenerateSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: W = X.wedge(Y)
            sage: W.nondegenerate_simplices()
            [v, e, f]
            sage: W.homology()
            {0: 0, 1: Z x Z}
            sage: S2 = simplicial_sets.Sphere(2)
            sage: X.wedge(S2).homology(reduced=False)
            {0: Z, 1: Z, 2: Z}
            sage: X.wedge(X).nondegenerate_simplices()
            [v, e, e']

        TESTS::

            sage: Z = SimplicialSet({e: (v,w)})
            sage: X.wedge(Z)
            Traceback (most recent call last):
            ...
            ValueError: both simplicial sets must be pointed
        """
        if not self.is_pointed() or not other.is_pointed():
            raise ValueError('both simplicial sets must be pointed')
        basept = self.base_point()
        X = other
        simps = set(self.nondegenerate_simplices())
        points = set([basept]).intersection(set([X.base_point()]))
        while not simps.intersection(set(X.nondegenerate_simplices())).issubset(points):
            X = copy.copy(X)
            points = set([basept]).intersection(set([X.base_point()]))
        # Now X and self have no simplices in common, except perhaps
        # for the base point. If the base point is not in common,
        # replace the base point in X with basept.
        data = X.face_data()
        y = X.base_point()
        if basept != y:
            if basept in data:
                data[basept] = None
            try:
                del data[y]
            except KeyError:
                pass
            for x in data:
                faces = data[x]
                new_faces = []
                if faces:
                    for f in faces:
                        if f.nondegenerate() == y:
                            new_faces.append(basept.apply_degeneracies(*f.degeneracies()))
                        else:
                            new_faces.append(f)
                    data[x] = tuple(new_faces)
        # Now that the base points agree, we just have to merge the
        # data for the two simplicial sets.
        data.update(self._data)
        return SimplicialSet(data, base_point=basept)

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

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: e = NonDegenerateSimplex(1, name='e')
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
        return ProductOfSimplicialSets([self] + list(others))

    def pushout(self, f, g):
        return PushoutOfSimplicialSets(self, f, g)

    def smash_product(self, other):
        """
        The smash product of this simplicial set with ``other``.

        EXAMPLES::

            sage: S1 = simplicial_sets.Sphere(1)
            sage: RP2 = simplicial_sets.RealProjectiveSpace(2)
            sage: X = S1.smash_product(RP2)
            sage: X.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 0 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2,
             3: Vector space of dimension 1 over Finite Field of size 2}

            sage: T = S1.product(S1)
            sage: X = T.smash_product(S1)
            sage: X.homology(reduced=False)
            {0: Z, 1: 0, 2: Z x Z, 3: Z}
        """
        prod = self.product(other)
        wedge = prod._wedge_as_subcomplex()
        return prod.quotient(wedge)

    def cone(self, reduced=None):
        r"""
        The (reduced) cone on this simplicial set.

        INPUT:

        - ``reduced`` (optional) -- boolean, default ``None``.

        If this simplicial set `X` is not pointed or ``reduced`` is
        ``False``, construct the ordinary cone: add a point `v` (which
        will become the base point) and for each simplex `\sigma` in
        `X`, add both `\sigma` and a simplex made up of `v` and
        `\sigma` (topologically, form the join of `v` and `\sigma`).

        If this simplicial set is pointed and ``reduced`` is ``True``
        or ``None``, then construct the reduced cone: take the
        quotient by the 1-simplex connecting the old base point to the
        new one.

        (If the simplicial set is not pointed and ``reduced`` is
        ``True``, raise an error.)

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: X = SimplicialSet({e: (v, v)})
            sage: CX = X.cone()  # unreduced cone, since X not pointed
            sage: CX.nondegenerate_simplices()
            [*, v, (v,*), e, (e,*)]
            sage: CX.base_point()
            *
            sage: X.set_base_point(v)
            sage: CX = X.cone()  # reduced cone
            sage: CX.nondegenerate_simplices()
            [*, e, (e,*)]

        TESTS::

            sage: X = SimplicialSet({e: (v,v)})
            sage: X.cone(reduced=True)
            Traceback (most recent call last):
            ...
            ValueError: this simplicial set is not pointed, so cannot construct the reduced cone
        """
        if self.is_pointed():
            if reduced is None:
                reduced = True
        else:
            if reduced:
                raise ValueError('this simplicial set is not pointed, so '
                                 'cannot construct the reduced cone')
        star = NonDegenerateSimplex(0, name='*')
        data = {}
        data[star] = None
        # Dictionary for translating old simplices to new: keys are
        # old simplices, corresponding value is the new simplex
        # (sigma, *).
        new_simplices = {None: star}
        for sigma in self.nondegenerate_simplices():
            new = NonDegenerateSimplex(sigma.dimension()+1,
                                       name='({},*)'.format(sigma))
            if sigma.dimension() == 0:
                data[sigma] = None
                if reduced:
                    if sigma == self.base_point():
                        edge = new
                data[new] = (star, sigma)
            else:
                sigma_faces = self.face_data()[sigma]
                data[sigma] = sigma_faces
                new_faces = [new_simplices[face.nondegenerate()].apply_degeneracies(*face.degeneracies())
                             for face in sigma_faces]
                data[new] = (new_faces + [sigma])
            new_simplices[sigma] = new
        Cone = SimplicialSet(data, base_point=star)
        if reduced:
            return Cone.quotient([edge])
        return Cone

    def suspension(self, n=1, reduced=None):
        """
        The (reduced) `n`-th suspension of this simplicial set.

        INPUT:

        - ``n`` (optional, default 1) -- integer, suspend this many
          times.

        - ``reduced`` (optional, default None) -- boolean

        If this simplicial set `X` is not pointed, return the
        suspension: the quotient `CX/X`, where `CX` is the (ordinary,
        unreduced) cone on `X`. If `X` is pointed, then by default,
        use the reduced cone instead, and so return the reduced
        suspension. If ``reduced`` is ``False``, then return the
        unreduced suspension.

        (If the simplicial set is not pointed and ``reduced`` is
        ``True``, raise an error.)

        EXAMPLES::

            sage: RP4 = simplicial_sets.RealProjectiveSpace(4)
            sage: S1 = simplicial_sets.Sphere(1)
            sage: SigmaRP4 = RP4.suspension()
            sage: S1_smash_RP4 = S1.smash_product(RP4)
            sage: SigmaRP4.homology() == S1_smash_RP4.homology()
            True

        The version of the suspension obtained by the smash product is
        typically less efficient than the model produced here::

            sage: SigmaRP4
            Simplicial set with 5 non-degenerate simplices
            sage: S1_smash_RP4
            Simplicial set with 25 non-degenerate simplices

        TESTS::

            sage: RP4.suspension(-3)
            Traceback (most recent call last):
            ...
            ValueError: n must be positive
        """
        if n == 0:
            return self
        if self.is_pointed():
            if reduced is None:
                reduced = True
        else:
            if reduced:
                raise ValueError('this simplicial set is not pointed, so '
                                 'cannot construct the reduced suspension')
        if n < 0:
            raise ValueError('n must be positive')
        CX = self.cone(reduced=reduced)
        X = set(self.nondegenerate_simplices()).intersection(CX.nondegenerate_simplices())
        Sigma = CX.quotient(CX.subcomplex(X))
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
            sage: eight = S1.wedge(S1)
            sage: eight.fundamental_group() # free group on 2 generators
            Finitely presented group < e0, e1 |  >

        The fundamental group of a disjoint union of course depends on
        the choice of base point::

            sage: T = simplicial_sets.Torus()
            sage: K = simplicial_sets.KleinBottle()
            sage: X = T.disjoint_union(K)

            sage: X.set_base_point(X.n_cells(0)[0])
            sage: X.fundamental_group().is_abelian()
            True
            sage: X.set_base_point(X.n_cells(0)[1])
            sage: X.fundamental_group().is_abelian()
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

        if len(gens) == 0:
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

        .. warning::

           Determining simple connectivity is not always possible,
           because it requires determining when a group, as given by
           generators and relations, is trivial. So this may give a
           false negative in some cases.

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

        .. warning::

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
            Set of Morphisms from S^3 to S^2 in Category of finite simplicial sets
            sage: Hom(S3, S2)
            Set of Morphisms from S^3 to S^2 in Category of finite simplicial sets
        """
        # Error-checking on the ``category`` argument is done when
        # calling Hom(X,Y), so no need to do it again here.
        if category is None:
            category = SimplicialSets().Finite()
        from sage.homology.simplicial_set_morphism import SimplicialSetHomset
        return SimplicialSetHomset(self, other, category=category)

    def _repr_(self):
        """
        Print representation.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0)
            sage: w = NonDegenerateSimplex(0)
            sage: degen = v.apply_degeneracies(0)
            sage: tau = NonDegenerateSimplex(2)
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


class ProductOfSimplicialSets(SimplicialSet):
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

        sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
        sage: v = NonDegenerateSimplex(0, name='v')
        sage: w = NonDegenerateSimplex(0, name='w')
        sage: e = NonDegenerateSimplex(1, name='e')
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
        The product of this simplicial set with ``other``.

        See :class:`ProductOfSimplicialComplexes` for more
        information.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet, ProductOfSimplicialSets
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: e = NonDegenerateSimplex(1)
            sage: X = SimplicialSet({e: (v, v)})
            sage: W = ProductOfSimplicialSets([X, X, X])
            sage: W.homology()
            {0: 0, 1: Z x Z x Z, 2: Z x Z x Z, 3: Z}
            sage: W.is_pointed()
            False

            sage: X.set_base_point(v)
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: f = NonDegenerateSimplex(1)
            sage: Y = SimplicialSet({f: (v,w)}, base_point=w)
            sage: Z = ProductOfSimplicialSets([Y, X])
            sage: Z.is_pointed()
            True
            sage: Z.base_point()
            (w, v)
        """
        nondegen = [factor.nondegenerate_simplices() for factor in factors]
        data_factors = [factor.face_data() for factor in factors]
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
                    simplex_factors = tuple(zip(simplices, tuple(degens)))
                    s = '(' + ', '.join(['{}'.format(_[0].apply_degeneracies(*_[1]))
                                         for _ in simplex_factors]) + ')'
                    simplex = NonDegenerateSimplex(d, name=s)
                    translate[simplex_factors] = simplex
                    # Now compute the faces of simplex.
                    if d == 0:
                        # It's a vertex, so it has no faces.
                        faces = None
                        break
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

        if all(sset.is_pointed() for sset in factors):
            basept = translate[tuple([(sset.base_point(), ()) for sset in factors])]
            SimplicialSet.__init__(self, data, base_point=basept)
        else:
             SimplicialSet.__init__(self, data)
        self._factors = tuple(factors)
        self._translation = tuple(translate.items())

    def factors(self):
        """
        The factors of this product of simplicial sets.

        EXAMPLES::

            sage: S2 = simplicial_sets.Sphere(2)
            sage: S3 = simplicial_sets.Sphere(3)
            sage: S2.product(S3).factors() == (S2, S3)
            True
            sage: S2.product(S3).factors() == (S3, S2)
            False
        """
        return self._factors

    def _wedge_as_subcomplex(self):
        """
        The wedge of this product of simplicial sets as a subcomplex.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, name='v')
            sage: e = NonDegenerateSimplex(1, name='e')
            sage: w = NonDegenerateSimplex(0, name='w')
            sage: f = NonDegenerateSimplex(1, name='f')
            sage: X = SimplicialSet({e: (v, v)}, base_point=v)
            sage: Y = SimplicialSet({f: (w, w)}, base_point=w)
            sage: P = X.product(Y)
            sage: W = P._wedge_as_subcomplex()
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
        return self.subcomplex(simps)


class PushoutOfSimplicialSets(SimplicialSet):
    def __init__(self, domain, f, g):
        r"""
        The pushout obtained from the morphisms `f` and `g`.

        Given a simplicial set `X` (= ``domain``) and maps `f: X \to Y`,
        `g: X \to Z`, construct the pushout `P`: see
        :wikipedia:`Pushout_(category_theory)`. This is constructed as
        pushouts of sets for each set of `n`-simplices, so `P_n` is
        the disjoint union of `Y_n` and `Z_n`, with elements `f(x)`
        and `g(x)` identified for each `x \in X_n`.

        Names: if simplices have names and are identified, choose
        their names from ``Y``.  Decorate somehow?

        Pointed if `f` and `g` are.

        EXAMPLES::

            sage: from sage.homology.simplicial_set import NonDegenerateSimplex, SimplicialSet
            sage: v = NonDegenerateSimplex(0, 'v')
            sage: a = NonDegenerateSimplex(0, 'a')
            sage: b = NonDegenerateSimplex(0, 'b')
            sage: c = NonDegenerateSimplex(0, 'c')
            sage: e0 = NonDegenerateSimplex(1, 'e_0')
            sage: e1 = NonDegenerateSimplex(1, 'e_1')
            sage: e2 = NonDegenerateSimplex(1, 'e_2')
            sage: X = SimplicialSet({e2: (b, a)})
            sage: Y = SimplicialSet({e2: (b,a), e0: (c,b), e1: (c,a)})
            sage: Z = simplicial_sets.Simplex(0)
            sage: f_data = {a:a, b:b, e2: e2}
            sage: v = Z.n_cells(0)[0]
            sage: g_data = {a:v, b:v, e2:v.apply_degeneracies(0)}
            sage: f = X.Hom(Y)(f_data)
            sage: g = X.Hom(Z)(g_data)
            sage: P = X.pushout(f,g)
        """
        from sage.homology.simplicial_set_morphism import SimplicialSetMorphism
        if (not isinstance(f, SimplicialSetMorphism) or
            not isinstance(g, SimplicialSetMorphism)):
            raise ValueError('f and g must be morphisms of simplicial sets')
        if not f.domain() == g.domain() == domain:
            raise ValueError('the domain of f and g must both equal "self"')
        # Data to define the pushout:
        data = {}
        X = f.domain()
        Y = f.codomain()
        Z = g.codomain()
        # Dictionaries to translate from simplices in Y, Z to simplices in P.
        Y_to_P = {}
        Z_to_P = {}
        for n in range(1 + max(Y.dimension(), Z.dimension())):
            # This isn't doing the right thing.



            Y_cells = {y:y for y in Y.n_cells(n)}
            Z_cells = {z:z for z in Z.n_cells(n)}
            for x in X.n_cells(n):
                # Identify f(x) with g(x)
                Y_cells[f(x)] = x
                Z_cells[g(x)] = x


            # Don't take a union here, in case Y and Z have cells in
            # common. Maybe make copies of one of them.
            cells = set(Y_cells.values()).union(Z_cells.values())
            data[n] = {'Y': Y_cells, 'Z': Z_cells}
            # data[n] = cells

        self._data = data


        if f.is_pointed() and g.is_pointed():
            if Y_to_P[Y.base_point()] != Z_to_P[Z.base_point()]:
                raise ValueError('something unexpected went wrong with base points')
            pass
            # SimplicialSet.__init__(self, data, base_point=Y_to_P[Y.base_point()])
        else:
            pass
            # SimplicialSet.__init__(self, data)
        # Would be good to keep track of the maps Y -> P, Z -> P. 
        # Also the original maps X -> Y, X -> Z?
        
        # self._f = f
        # self._g = g
        # self._fbar = fbar
        # self._gbar = gbar


########################################################################
# Functions for manipulating face and degeneracy maps.

def standardize_degeneracies(*L):
    r"""
    INPUT:

    - ``L`` -- list of integers, representing a composition of
      degeneracies in a simplicial set.

    OUTPUT: an equivalent list of degeneracies, standardized to be
    written in decreasing order, using the simplicial identity

    .. math::

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

    .. math::

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

    .. math::

        d_i d_j &= d_{j-1} d_i \ \   \text{if }   i<j.

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

    .. math::

        d_m s_{i_1} s_{i_2} ... s_{i_n}

    in one of the forms

    .. math::

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
    The `n`-sphere as a simplicial set,
    constructed with two non-degenerate simplices:
    a vertex `v_0` and an `n`-simplex `\sigma_n`.

    INPUT:

    - ``n`` -- integer

    EXAMPLES::

        sage: S0 = simplicial_sets.Sphere(0)
        sage: S0
        Simplicial set with 2 non-degenerate simplices
        sage: S0.nondegenerate_simplices()
        [v_0, w_0]
        sage: simplicial_sets.Sphere(4)
        S^4
        sage: simplicial_sets.Sphere(4).nondegenerate_simplices()
        [v_0, sigma_4]
    """
    v_0 = NonDegenerateSimplex(0, name='v_0')
    if n == 0:
        w_0 = NonDegenerateSimplex(0, name='w_0')
        return SimplicialSet({v_0: None, w_0: None})
    degens = range(n-2, -1, -1)
    degen_v = v_0.apply_degeneracies(*degens)
    sigma = NonDegenerateSimplex(n, name='sigma_{}'.format(n))
    return SimplicialSet({sigma: [degen_v] * (n+1)}, base_point=v_0, name='S^{}'.format(n))


def ClassifyingSpace(group, n):
    r"""
    The `n` skeleton of the classifying space of ``group``,
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
    return SimplicialSet(delta_complexes.KleinBottle())


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


def Horn(n, k):
    r"""
    The horn $\Lambda^n_k$.

    This is the subcomplex of the $n$-simplex obtained by removing its
    $k$th face.

    EXAMPLES::

        sage: L = simplicial_sets.Horn(3, 0)
        sage: L.n_cells(3)
        []
        sage: L.n_cells(2)
        [(0, 1, 2), (0, 1, 3), (0, 2, 3)]
    """
    K = Simplex(n)
    sigma = K.n_cells(n)[0]
    return K.subcomplex(K.faces(sigma)[:k] + K.faces(sigma)[k+1:])


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

        sage: simplicial_sets.ComplexProjectiveSpace(5)
        Traceback (most recent call last):
        ...
        ValueError: complex projective spaces are only available in dimensions between 0 and 2
    """
    if n < 0 or n > 2:
        raise ValueError('complex projective spaces are only available in dimensions between 0 and 2')
    if n == 0:
        K = Simplex(0)
        v = K.n_cells(0)[0]
        K.set_base_point(v)
        return K
    if n == 1:
        return Sphere(2)
    if n == 2:
        v = NonDegenerateSimplex(0, name='<<GBar>>')
        f2_1 = NonDegenerateSimplex(2, name='<<GBar<- (1)><- NIL>>>')
        f2_2 = NonDegenerateSimplex(2, name='<<GBar<- (2)><- NIL>>>')
        f3_110 = NonDegenerateSimplex(3, name='<<GBar<- (1 1)><0 NIL><- NIL>>>')
        f3_011 = NonDegenerateSimplex(3, name='<<GBar<0 (1)><- (1)><- NIL>>>')
        f3_111 = NonDegenerateSimplex(3, name='<<GBar<1 (1)><- (1)><- NIL>>>')
        f4_101101 = NonDegenerateSimplex(4, name='<<GBar<1-0 (1)><1-0 NIL><- (1)><- NIL>>>')
        f4_201110 = NonDegenerateSimplex(4, name='<<GBar<2-0 (1)><1 (1)><0 NIL><- NIL>>>')
        f4_211010 = NonDegenerateSimplex(4, name='<<GBar<2-1 (1)><0 (1)><0 NIL><- NIL>>>')
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
