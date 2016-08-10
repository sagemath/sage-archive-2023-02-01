r"""
Cluster algebras

Implementation of cluster algebras as an algebra using mainly structural theorems from CA IV

TODO: We should write a nice paragraph here.

AUTHORS:

- Dylan Rupel (2015-06-15): initial version

- Salvatore Stella (2015-06-15): initial version

EXAMPLES::

    sage: A = ClusterAlgebra(['A',2])
    sage: A.b_matrix()
    [ 0  1]
    [-1  0]
    sage: A.explore_to_depth(2)
    sage: A.g_vectors_so_far()
    [(0, 1), (0, -1), (1, 0), (-1, 1), (-1, 0)]
    sage: A.F_polynomials_so_far()
    [1, u1 + 1, 1, u0 + 1, u0*u1 + u0 + 1]
    sage: A.cluster_variables_so_far()
    [x1, (x0 + 1)/x1, x0, (x1 + 1)/x0, (x0 + x1 + 1)/(x0*x1)]
    sage: B = A.mutate_initial(0)
    sage: B.b_matrix()
    [ 0 -1]
    [ 1  0]
    sage: B.g_vectors_so_far()
    [(0, 1), (0, -1), (1, 0), (1, -1), (-1, 0)]
    sage: B.F_polynomials_so_far()
    [1, u0*u1 + u1 + 1, 1, u1 + 1, u0 + 1]
    sage: B.cluster_variables_so_far()
    [x1, (x0 + x1 + 1)/(x0*x1), x0, (x0 + 1)/x1, (x1 + 1)/x0]
"""

#*****************************************************************************
#       Copyright (C) 2015 Dylan Rupel and Salvatore Stella
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from functools import wraps
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.rings import Rings
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.combinat.permutation import Permutation
from sage.functions.generalized import sign
from sage.functions.other import binomial
from sage.geometry.cone import Cone
from sage.geometry.fan import Fan
from sage.matrix.constructor import identity_matrix, matrix
from sage.matrix.special import block_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing, LaurentPolynomialRing_generic
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from types import MethodType

##############################################################################
# Helper functions
##############################################################################
def _mutation_parse(mutate):
    r"""
    Preparse input for mutation functions.

    This wrapper provides:
        - inplace (only for seeds)
        - mutate along sequence
        - mutate at all sinks/sources

    Possible things to implement later include:
        - mutate at a cluster variariable
        - mutate at a g-vector (it is hard to distinguish this case from a generic sequence)
        - urban renewals
        - other?
    """
    doc = mutate.__doc__.split("INPUT:")
    doc[0] += "INPUT:"
    if mutate.__name__ == "mutate":
        doc[0] += r"""

        - ``inplace`` -- bool (default True) whether to mutate in place or to return a new object
    """
    doc[0] += r"""

        - ``direction`` -- in which direction(s) to mutate. It can be
            - an integer in ``range(self.rk())`` to mutate in one direction only;
            - an iterable of such integers to mutate along a sequence;
            - a string "sinks" or "sources" to mutate at all sinks or sources simultaneously.
    """
    mutate.__doc__ = doc[0] + doc[1]

    @wraps(mutate)
    def mutate_wrapper(self, direction, **kwargs):
        inplace = kwargs.pop('inplace', True) and mutate.__name__ != "mutate_initial"
        if inplace:
            to_mutate = self
        else:
            to_mutate = copy(self)

        if direction == "sinks":
            B = self.b_matrix()
            seq = [ i for i in range(B.ncols()) if all( x<=0 for x in B.column(i) ) ]
        elif direction == "sources":
            B = self.b_matrix()
            seq = [ i for i in range(B.ncols()) if all( x>=0 for x in B.column(i) ) ]
        else:
            try:
                seq = iter(direction)
            except TypeError:
                seq = iter((direction,))

        for k in seq:
            mutate(to_mutate, k, **kwargs)

        if not inplace:
            return to_mutate

    return mutate_wrapper

##############################################################################
# Elements of a cluster algebra
##############################################################################

class ClusterAlgebraElement(ElementWrapper):

    def __init__(self, parent, value):
        r"""
        An element of a Cluster Algebra.

        INPUT:

        - ``parent`` -- a :class:`ClusterAlgebra`: the algebra to which the element belongs;

        - ``value`` -- the value of the element.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: from sage.algebras.cluster_algebra import ClusterAlgebraElement
            sage: ClusterAlgebraElement(A,1)
            1
        """
        ElementWrapper.__init__(self, parent=parent, value=value)

        # setup methods defined only in special cases
        if parent._is_principal:
            self.g_vector = MethodType(g_vector, self, self.__class__)
            self.is_homogeneous = MethodType(is_homogeneous, self, self.__class__)
            self.homogeneous_components = MethodType(homogeneous_components, self, self.__class__)

    # AdditiveMagmas.Subobjects currently does not implements _add_
    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        INPUT:

        - ``other`` - an element of ``self.parent()``

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: A.an_element() + A.an_element()
            2*x0
        """
        return self.parent().retract(self.lift() + other.lift())

    def _div_(self, other):
        r"""
        Return the quotient of ``self`` and ``other``.

        WARNING::

            This method returns an element of ``self.parent().ambient()``
            rather than an element of ``self.parent()`` because, a priori,
            we cannot guarantee membership.

            You can force the result to be an element of ``self.parent()``
            by feeding it into ``self.parent().retract``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: A.an_element() / A.an_element()
            1
        """
        return self.lift()/other.lift()

    def d_vector(self):
        r"""
        Return the denominator vector of ``self`` as a tuple of integers.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4], principal_coefficients=True)
            sage: A.current_seed().mutate([0, 2, 1])
            sage: x = A.cluster_variable((-1, 2, -2, 2)) * A.cluster_variable((0,0,0,1))**2
            sage: x.d_vector()
            (1, 1, 2, -2)
        """
        monomials = self.lift()._dict().keys()
        minimal = map(min, zip(*monomials))
        return tuple(-vector(minimal))[:self.parent().rk()]

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4], principal_coefficients=True)
            sage: A.current_seed().mutate([0, 2, 1])
            sage: A.cluster_variable((-1, 2, -2, 2))
            (x0*x2^2*y0*y1*y2^2 + x1^3*x3^2 + x1^2*x3^2*y0 + 2*x1^2*x3*y2 + 2*x1*x3*y0*y2 + x1*y2^2 + y0*y2^2)/(x0*x1*x2^2)
        """
        numer, denom = self.lift()._fraction_pair()
        return repr(numer/denom)

####
# Methods not always defined
####

def g_vector(self):   #READY
    r"""
    Return the g-vector of ``self``.

    EXAMPLES::

        sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
        sage: A.cluster_variable((1,0)).g_vector() == (1,0)
        True
    """
    components = self.homogeneous_components()
    if len(components) == 1:
        return components.keys()[0]
    else:
        raise ValueError("This element is not homogeneous.")

def is_homogeneous(self):
    r"""
    Return ``True`` if ``self`` is an homogeneous element of ``self.parent()``.

    EXAMPLES::

        sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
        sage: A.cluster_variable((1,0)).is_homogeneous()
        True
        sage: x = A.cluster_variable((1,0)) + A.cluster_variable((0,1))
        sage: x.is_homogeneous()
        False
    """
    return len(self.homogeneous_components()) == 1

def homogeneous_components(self):
    r"""
    Return a dictionary of the homogeneous components of ``self``.

    OUTPUT:

    A dictionary whose keys are homogeneous degrees and whose values are the
    summands of ``self`` of the given degree.

    EXAMPLES::

        sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
        sage: x = A.cluster_variable((1,0)) + A.cluster_variable((0,1))
        sage: x.homogeneous_components()
        {(0, 1): x1, (1, 0): x0}
    """
    deg_matrix = block_matrix([[identity_matrix(self.parent().rk()),-self.parent().b_matrix()]])
    components = dict()
    x = self.lift()
    monomials = x.monomials()
    for m in monomials:
        g_vect = tuple(deg_matrix*vector(m.exponents()[0]))
        if g_vect in components:
            components[g_vect] += self.parent().retract(x.monomial_coefficient(m)*m)
        else:
            components[g_vect] = self.parent().retract(x.monomial_coefficient(m)*m)
    return components


##############################################################################
# Seeds
##############################################################################

class ClusterAlgebraSeed(SageObject):

    def __init__(self, B, C, G, parent, **kwargs):
        r"""
        A seed in a Cluster Algebra.

        INPUT:

        - ``B`` -- a skew-symmetrizable integer matrix;

        - ``C`` -- the matrix of c-vectors of ``self``;

        - ``G`` -- the matrix of g-vectors of ``self``;

        - ``parent`` -- a :class:`ClusterAlgebra`: the algebra to which the
          seed belongs;

        - ``path`` -- list (default []) the mutation sequence from the initial
          seed of ``parent`` to `self``

        WARNING:

            Seeds should **not** be created manually: no test is performed to
            assert that they are built from consistent data nor that they
            really are seeds of ``parent``. If you create seeds with
            inconsistent data all sort of things can go wrong, even
            :meth:`__eq__` is no longer guaranteed to give correct answers.
            Use at your own risk.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: from sage.algebras.cluster_algebra import ClusterAlgebraSeed
            sage: ClusterAlgebraSeed(A.b_matrix(),identity_matrix(4),identity_matrix(4),A,path=[1,2,3])
            The seed of a Cluster Algebra with cluster variables x0, x1, x2, x3 and no coefficients over Integer Ring obtained from the initial by mutating along the sequence [1, 2, 3]
        """
        self._B = copy(B)
        self._C = copy(C)
        self._G = copy(G)
        self._parent = parent
        self._path = kwargs.get('path', [])

    def __copy__(self):
        r"""
        Return a copy of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = copy(A.current_seed())
            sage: S == A.current_seed()
            True
            sage: S is not A.current_seed()
            True
        """
        other = type(self).__new__(type(self))
        other._B = copy(self._B)
        other._C = copy(self._C)
        other._G = copy(self._G)
        other._parent = self._parent
        other._path = copy(self._path)
        return other

    def __eq__(self, other):
        r"""
        Test equality of two seeds.

        INPUT:

        - ``other`` -- a :class:`ClusterAlgebraSeed`

        ALGORITHM:

            ``self`` and ``other`` are deemed to be equal if they have the same
            parent and their set of g-vectors coincide, i.e. this tests
            equality of unlabelled seeds.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = copy(A.current_seed())
            sage: S.mutate([0,2,0])
            sage: S == A.current_seed()
            False
            sage: S.mutate(2)
            sage: S == A.current_seed()
            True

            sage: B = ClusterAlgebra(['B',2],principal_coefficients=True)
            sage: S = B.current_seed()
            sage: S.mutate(0)
            sage: S == B.current_seed()
            True
        """
        return type(self) == type(other) and self.parent() == other.parent() and frozenset(self.g_vectors()) == frozenset(other.g_vectors())

    def __contains__(self, element):
        r"""
        Test whether ``element`` belong to ``self``

        INPUT:

        - ``element`` -- either a g-vector or an element of :meth:`parent`

        EXAMPLES::
            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: (1,0,0) in S
            True
            sage: (1,1,0) in S
            False
            sage: A.cluster_variable((1,0,0)) in S
            True
        """
        if isinstance(element, ClusterAlgebraElement):
            cluster = self.cluster_variables()
        else:
            element = tuple(element)
            cluster = self.g_vectors()
        return element in cluster

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.current_seed(); S
            The initial seed of a Cluster Algebra with cluster variables x0, x1, x2 and no coefficients over Integer Ring
            sage: S.mutate(0); S
            The seed of a Cluster Algebra with cluster variables x0, x1, x2 and no coefficients over Integer Ring obtained from the initial by mutating in direction 0
            sage: S.mutate(1); S
            The seed of a Cluster Algebra with cluster variables x0, x1, x2 and no coefficients over Integer Ring obtained from the initial by mutating along the sequence [0, 1]
        """
        if self._path == []:
            return "The initial seed of a %s"%str(self.parent())[2:]
        elif len(self._path) == 1:
            return "The seed of a %s obtained from the initial by mutating in direction %s"%(str(self.parent())[2:],str(self._path[0]))
        else:
            return "The seed of a %s obtained from the initial by mutating along the sequence %s"%(str(self.parent())[2:],str(self._path))

    def parent(self):
        r"""
        Return the parent of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B',3])
            sage: A.current_seed().parent() == A
            True
        """
        return self._parent

    def depth(self):
        r"""
        Return the length of a mutation sequence from the initial seed of :meth:`parent` to ``self``.

        WARNING:
            This is the length of the mutation sequence returned by
            :meth:`path_from_initial_seed` which need not be the shortest
            possible.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: S1 = A.initial_seed()
            sage: S1.mutate([0,1,0,1])
            sage: S1.depth()
            4
            sage: S2 = A.initial_seed()
            sage: S2.mutate(1)
            sage: S2.depth()
            1
            sage: S1 == S2
            True
        """
        return len(self._path)

    def path_from_initial_seed(self):
        r"""
        Return a mutation sequence from the initial seed of :meth:`parent` to ``self``.

        WARNING:

            This is the path used to compute ``self`` and it does not have to
            be the shortest possible.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: S1 = A.initial_seed()
            sage: S1.mutate([0,1,0,1])
            sage: S1.path_from_initial_seed()
            [0, 1, 0, 1]
            sage: S2 = A.initial_seed()
            sage: S2.mutate(1)
            sage: S2.path_from_initial_seed()
            [1]
            sage: S1 == S2
            True
        """
        return copy(self._path)

    def b_matrix(self):
        r"""
        Return the exchange matrix of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.b_matrix()
            [ 0  1  0]
            [-1  0 -1]
            [ 0  1  0]
        """
        return copy(self._B)

    def c_matrix(self):
        r"""
        Return the matrix whose columns are the c-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.c_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return copy(self._C)

    def c_vector(self, j):
        r"""
        Return the j-th c-vector of ``self``.

        INPUT:

        - ``j`` -- an integer in ``range(self.parent().rk())``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.c_vector(0)
            (1, 0, 0)
            sage: S.mutate(0)
            sage: S.c_vector(0)
            (-1, 0, 0)
            sage: S.c_vector(1)
            (1, 1, 0)
        """
        return tuple(self._C.column(j))

    def c_vectors(self):
        r"""
        Return all the c-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.c_vectors()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return map(tuple, self._C.columns())

    def g_matrix(self):
        r"""
        Return the matrix whose columns are the g-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.g_matrix()
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        return copy(self._G)

    def g_vector(self, j):
        r"""
        Return the j-th g-vector of ``self``.

        INPUT:

        - ``j`` -- an integer in ``range(self.parent().rk())``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.g_vector(0)
            (1, 0, 0)
        """
        return tuple(self._G.column(j))

    def g_vectors(self):
        r"""
        Return all the g-vectors of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.g_vectors()
            [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        """
        return map(tuple, self._G.columns())

    def F_polynomial(self, j):
        r"""
        Return the j-th F-polynomial of ``self``.

        INPUT:

        - ``j`` -- an integer in ``range(self.parent().rk())``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.F_polynomial(0)
            1
        """
        return self.parent().F_polynomial(self.g_vector(j))

    def F_polynomials(self):
        r"""
        Return all the F-polynomials of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.F_polynomials()
            [1, 1, 1]
        """
        return [self.parent().F_polynomial(g) for g in self.g_vectors()]

    def cluster_variable(self, j):
        r"""
        Return the j-th cluster variable of ``self``.

        INPUT:

        - ``j`` -- an integer in ``range(self.parent().rk())``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.cluster_variable(0)
            x0
            sage: S.mutate(0)
            sage: S.cluster_variable(0)
            (x1 + 1)/x0
        """
        return self.parent().cluster_variable(self.g_vector(j))

    def cluster_variables(self):
        r"""
        Return all the cluster variables of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.initial_seed()
            sage: S.cluster_variables()
            [x0, x1, x2]
        """
        return [self.parent().cluster_variable(g) for g in self.g_vectors()]

    @_mutation_parse
    def mutate(self, k, mutating_F=True):
        r"""
        Mutate ``self``.

        INPUT:

        - ``k`` --  an integer in ``range(self.parent().rk())``: the direction
          in which we are mutating

        - ``mutating_F`` -- bool (default True) whether to compute F-polynomials
          also. While knowing F-polynomials is essential to computing
          cluster variables, the process of mutating them is quite slow. If you
          care only about combinatorial data like g-vectors and c-vectors,
          setting ``mutating_F=False`` yields significant benefits in terms of
          speed.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: S = A.initial_seed()
            sage: S.mutate(0); S
            The seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring obtained from the initial by mutating in direction 0
            sage: S.mutate(5)
            Traceback (most recent call last):
            ...
            ValueError: Cannot mutate in direction 5.
        """
        n = self.parent().rk()

        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction ' + str(k) + '.')

        # store new mutation path
        if self._path != [] and self._path[-1] == k:
            self._path.pop()
        else:
            self._path.append(k)

        # find sign of k-th c-vector
        if any(x > 0 for x in self._C.column(k)):
            eps = +1
        else:
            eps = -1

        # store the g-vector to be mutated in case we are mutating F-polynomials also
        old_g_vector = self.g_vector(k)

        # compute new G-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[j,k] += max(0, -eps*self._B[j,k])
        J[k,k] = -1
        self._G = self._G*J

        # path to new g-vector (we store the shortest encountered so far)
        g_vector = self.g_vector(k)
        if not self.parent()._path_dict.has_key(g_vector) or len(self.parent()._path_dict[g_vector]) > len(self._path):
            self.parent()._path_dict[g_vector] = copy(self._path)

        # compute F-polynomials
        if mutating_F and not self.parent()._F_poly_dict.has_key(g_vector):
            self.parent()._F_poly_dict[g_vector] = self._mutated_F(k, old_g_vector)

        # compute new C-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[k,j] += max(0, eps*self._B[k,j])
        J[k,k] = -1
        self._C = self._C*J

        # compute new B-matrix
        self._B.mutate(k)

    def _mutated_F(self, k, old_g_vector):
        r"""
        Compute new F-polynomial obtained by mutating in direction ``k``.

        INPUT:

        - ``k`` --  an integer in ``range(self.parent().rk())``: the direction
          in which we are mutating

        - ``old_g_vector`` -- tuple: the k-th g-vector of ``self`` before
          mutating

        NOTE:

            This function is the bottleneck of :meth:`mutate`. The problem is
            that operations on polynomials are slow. One can get a significant
            speed boost by disabling this method calling :meth:`mutate` with
            ``mutating_F=False``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: S = A.initial_seed()
            sage: S.mutate(0)
            sage: S._mutated_F(0,(1,0))
            u0 + 1
        """
        alg = self.parent()
        pos = alg._U(1)
        neg = alg._U(1)
        for j in xrange(alg.rk()):
            if self._C[j,k] > 0:
                pos *= alg._U.gen(j)**self._C[j,k]
            else:
                neg *= alg._U.gen(j)**(-self._C[j,k])
            if self._B[j,k] > 0:
                pos *= self.F_polynomial(j)**self._B[j,k]
            elif self._B[j,k] <0:
                neg *= self.F_polynomial(j)**(-self._B[j,k])
        return (pos+neg)/alg.F_polynomial(old_g_vector)

##############################################################################
# Cluster algebras
##############################################################################

class ClusterAlgebra(Parent):

    Element = ClusterAlgebraElement

    def __init__(self, data, **kwargs): # READY
        r"""
        A Cluster Algebra.

        INPUT:

        - ``data`` -- some data defining a cluster algebra. It can be anything
          that can be parsed by :class:`ClusterQuiver`.

        - ``scalars`` -- (default ZZ) the scalars on which the cluster algebra
          is defined.

        - ``cluster_variable_prefix`` -- string (default 'x'); it needs to be
          a valid variable name.

        - ``cluster_variable_names`` -- a list of strings.  Supersedes
          ``cluster_variable_prefix``. Each element needs to be a valid
          variable name.

        - ``coefficient_prefix`` -- string (default 'y'); it needs to be
          a valid variable name.

        - ``coefficient_names`` -- a list of strings. Supersedes
          ``cluster_variable_prefix``. Each element needs to be a valid
          variable name.

        - ``principal_coefficients`` -- bool (default False). Supersedes any
          coefficient defined by ``data``.

        EXAMPLES::

            sage: B = matrix([(0, 1, 0, 0),(-1, 0, -1, 0),(0, 1, 0, 1),(0, 0, -2, 0),(-1, 0, 0, 0),(0, -1, 0, 0)])
            sage: A = ClusterAlgebra(B); A
            A Cluster Algebra with cluster variables x0, x1, x2, x3 and coefficients y0, y1 over Integer Ring
            sage: A.gens()
            [x0, x1, x2, x3, y0, y1]
            sage: A = ClusterAlgebra(['A',2]); A
            A Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
            sage: A = ClusterAlgebra(['A',2], principal_coefficients=True); A.gens()
            [x0, x1, y0, y1]
            sage: A = ClusterAlgebra(['A',2], principal_coefficients=True, coefficient_prefix='x'); A.gens()
            [x0, x1, x2, x3]
            sage: A = ClusterAlgebra(['A',3], principal_coefficients=True, cluster_variable_names=['a','b','c']); A.gens()
            [a, b, c, y0, y1, y2]
            sage: A = ClusterAlgebra(['A',3], principal_coefficients=True, cluster_variable_names=['a','b'])
            Traceback (most recent call last):
            ...
            ValueError: cluster_variable_names should be a list of 3 valid variable names

        ALGORITHM:

            The implementation is mainly based on [FZ07]_ and [NZ12]_.

        REFERENCES:

        .. [FZ07] \S. Fomin and \A. Zelevinsky, "Cluster algebras IV.
           Coefficients", Compos. Math. 143 (2007), no. 1, 112-164.

        .. [NZ12] \T. Nakanishi and \A. Zelevinsky, "On tropical dualities in
           cluster algebras', Algebraic groups and quantum groups, Contemp.
           Math., vol. 565, Amer. Math. Soc., Providence, RI, 2012, pp.
           217-226.
        """
        # Temporary variables
        Q = ClusterQuiver(data)
        n = Q.n()
        I = identity_matrix(n)
        if kwargs.get('principal_coefficients', False):
            M0 = I
        else:
            M0 = Q.b_matrix()[n:,:]
        B0 = block_matrix([[Q.b_matrix()[:n,:]],[M0]])
        m = M0.nrows()

        # Ambient space for F-polynomials
        # NOTE: for speed purposes we need to have QQ here instead of the more
        # natural ZZ. The reason is that _mutated_F is faster if we do not cast
        # the result to polynomials but then we get "rational" coefficients
        self._U = PolynomialRing(QQ, ['u%s'%i for i in xrange(n)])

        # Storage for computed data
        self._path_dict = dict([ (v, []) for v in map(tuple,I.columns()) ])
        self._F_poly_dict = dict([ (v, self._U(1)) for v in self._path_dict ])

        # Determine the names of the initial cluster variables
        variables_prefix = kwargs.get('cluster_variable_prefix','x')
        variables = list(kwargs.get('cluster_variable_names', [variables_prefix+str(i) for i in xrange(n)]))
        if len(variables) != n:
             raise ValueError("cluster_variable_names should be a list of %d valid variable names"%n)

        # Determine scalars
        scalars = kwargs.get('scalars', ZZ)

        # Determine coefficients and setup self._base
        if m>0:
            coefficient_prefix = kwargs.get('coefficient_prefix', 'y')
            if coefficient_prefix == variables_prefix:
                offset = n
            else:
                offset = 0
            coefficients = list(kwargs.get('coefficient_names', [coefficient_prefix+str(i) for i in xrange(offset,m+offset)]))
            if len(coefficients) != m:
                raise ValueError("coefficient_names should be a list of %d valid variable names"%m)
            base = LaurentPolynomialRing(scalars, coefficients)
        else:
            base = scalars
            coefficients = []

        # setup Parent and ambient
        self._ambient = LaurentPolynomialRing(scalars, variables+coefficients)
        Parent.__init__(self, base=base, category=Rings(scalars).Commutative().Subobjects(), names=variables+coefficients)

        # Data to compute cluster variables using separation of additions
        self._y = dict([ (self._U.gen(j), prod([self._base.gen(i)**M0[i,j] for i in xrange(m)])) for j in xrange(n)])
        self._yhat = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(n+m)])) for j in xrange(n)])

        # Have we got principal coefficients?
        self._is_principal = (M0 == I)

        # Store initial data
        # NOTE: storing both _B0 as rectangular matrix and _yhat is redundant.
        # We keep both around for speed purposes.
        self._B0 = copy(B0)
        self._n = n
        self.reset_current_seed()

        # Internal data for exploring the exchange graph
        self.reset_exploring_iterator()

        # Add methods that are defined only for special cases
        if n == 2 and m == 0:
            self.greedy_element = MethodType(greedy_element, self, self.__class__)
            self._greedy_coefficient = MethodType(_greedy_coefficient, self, self.__class__)

        # Register embedding into self.ambient()
        embedding = SetMorphism(Hom(self,self.ambient()), lambda x: x.lift())
        self._populate_coercion_lists_(embedding=embedding)

    def __copy__(self): # READY
        r"""
        Return a copy of ``self``.

        EXAMPLES::

            sage: A1 = ClusterAlgebra(['A',3])
            sage: A2 = copy(A1)
            sage: A2 == A1
            True
            sage: A2 is not A1
            True
        """
        n = self.rk()
        cv_names = self.initial_cluster_variable_names()
        coeff_names = self.coefficient_names()
        other = ClusterAlgebra(self._B0, cluster_variable_names=cv_names,
                coefficient_names=coeff_names, scalars=self.scalars())
        other._F_poly_dict = copy(self._F_poly_dict)
        other._path_dict = copy(self._path_dict)
        S = self.current_seed()
        S._parent = other
        other.set_current_seed(S)
        return other

    def __eq__(self, other):    # READY
        r"""
        Test equality of two cluster algebras.

        INPUT:

        - ``other`` -- a :class:`ClusterAlgebra`

        ALGORITHM:

            ``self`` and ``other`` are deemed to be equal if they have the same
            initial exchange matrix and their ambients coincide. In
            particular we do not keep track of how much each algebra has been
            explored.

        EXAMPLES::

            sage: A1 = ClusterAlgebra(['A',3])
            sage: A2 = copy(A1)
            sage: A1 is not A2
            True
            sage: A1 == A2
            True
            sage: A1.current_seed().mutate([0,1,2])
            sage: A1 == A2
            True
        """
        return type(self) == type(other) and self._B0 == other._B0 and self.ambient() == other.ambient()

    def _repr_(self):   # READY
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(matrix(1),principal_coefficients=True); A
            A Cluster Algebra with cluster variable x0 and coefficient y0 over Integer Ring
            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
        """
        var_names = self.initial_cluster_variable_names()
        var_names = (" " if len(var_names)==1 else "s ") + ", ".join(var_names)
        coeff_names = self.coefficient_names()
        coeff_prefix = " and" +(" " if len(coeff_names) >0 else " no ") + "coefficient"
        coeff = coeff_prefix + (" " if len(coeff_names)==1 else "s ") + ", ".join(coeff_names) + (" " if len(coeff_names)>0 else "")
        return "A Cluster Algebra with cluster variable" + var_names + coeff + "over " + repr(self.scalars())

    def _an_element_(self): # READY
        r"""
        Return an element of ``self``.

        EXAMPLES::
            sage: A = ClusterAlgebra(['A',2])
            sage: A.an_element()
            x0
        """
        return self.current_seed().cluster_variable(0)

    def _coerce_map_from_(self, other): # READY
        r"""
        Test whether there is a coercion from ``other`` to ``self``.

        ALGORITHM:

            If ``other`` is an instance of :class:`ClusterAlgebra` then allow
            coercion if ``other.ambient()`` can be coerced into
            ``self.ambient()`` and other can be obtained from ``self`` by
            permuting variables and coefficients and/or freezing some initial
            cluster variables.

            Otherwise allow anything that coerces into ``self.base()`` to coerce
            into ``self``.

        EXAMPLES::

            sage: B1 = matrix([(0, 1, 0, 0),(-1, 0, -1, 0),(0, 1, 0, 1),(0, 0, -2, 0),(-1, 0, 0, 0),(0, -1, 0, 0)])
            sage: B2 = B1.matrix_from_columns([0,1,2])
            sage: A1 = ClusterAlgebra(B1, coefficient_prefix='x')
            sage: A2 = ClusterAlgebra(B2, coefficient_prefix='x')
            sage: A1.has_coerce_map_from(A2)
            True
            sage: A2.has_coerce_map_from(A1)
            False
            sage: f = A1.coerce_map_from(A2)
            sage: A2.find_g_vector((-1, 1, -1))
            [0, 2, 1]
            sage: S = A1.initial_seed(); S.mutate([0, 2, 1])
            sage: S.cluster_variable(1) == f(A2.cluster_variable((-1, 1, -1)))
            True
        """
        if isinstance(other, ClusterAlgebra):
            gen_s = self.gens()
            gen_o = other.gens()
            if len(gen_s) == len(gen_o):
                f = self.ambient().coerce_map_from(other.ambient())
            if f is not None:
                perm = Permutation([ gen_s.index(self(f(v)))+1 for v in gen_o ]).inverse()
                n = self.rk()
                m = len(perm) - n
                M = self._B0[n:,:]
                B = block_matrix([[self.b_matrix(),-M.transpose()],[M,matrix(m)]])
                B.permute_rows_and_columns(perm,perm)
                return B.matrix_from_columns(range(other.rk())) == other._B0

        # everything that is in the base can be coerced to self
        return self.base().has_coerce_map_from(other)

    def rk(self):   # READY
        r"""
        Return the rank of ``self`` i.e. the number of cluster variables in any seed.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
            sage: A.rk()
            2
        """
        return self._n

    def current_seed(self): # READY
        r"""
        Return the current seed of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.current_seed()
            The initial seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
        """
        return self._seed

    def set_current_seed(self, seed):   # READY
        r"""
        Set the value reported by :meth:`current_seed`  to ``seed`` if it makes sense.

        INPUT:

        - ``seed`` -- an instance of :class:`ClusterAlgebraSeed`

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: S = copy(A.current_seed())
            sage: S.mutate([0,1,0])
            sage: A.current_seed() == S
            False
            sage: A.set_current_seed(S)
            sage: A.current_seed() == S
            True
        """
        if self.contains_seed(seed):
            self._seed = seed
        else:
            raise ValueError("This is not a seed in this cluster algebra.")

    def contains_seed(self, seed):  # READY
        r"""
        Test if ``seed`` is a seed in ``self``.

        INPUT:

        - ``seed`` -- an instance of :class:`ClusterAlgebraSeed`

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True); A
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
            sage: S = copy(A.current_seed())
            sage: A.contains_seed(S)
            True
        """
        computed_sd = self.initial_seed()
        computed_sd.mutate(seed._path, mutating_F=False)
        return computed_sd == seed

    def reset_current_seed(self):   # READY
        r"""
        Reset the value reported by :meth:`current_seed` to :meth:`initial_seed`.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.current_seed().mutate([1,0])
            sage: A.current_seed() == A.initial_seed()
            False
            sage: A.reset_current_seed()
            sage: A.current_seed() == A.initial_seed()
            True
        """
        self._seed = self.initial_seed()

    def initial_seed(self): # READY
        r"""
        Return the initial seed of ``self``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.initial_seed()
            The initial seed of a Cluster Algebra with cluster variables x0, x1 and no coefficients over Integer Ring
        """
        n = self.rk()
        I = identity_matrix(n)
        return ClusterAlgebraSeed(self.b_matrix(), I, I, self)

    def b_matrix(self): # READY
        r"""
        Return the initial exchange matrix of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.b_matrix()
            [ 0  1]
            [-1  0]
        """
        n = self.rk()
        return copy(self._B0[:n,:])

    def g_vectors_so_far(self): # READY
        r"""
        Return a list of the g-vectors of cluster variables encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.current_seed().mutate(0)
            sage: A.g_vectors_so_far()
            [(0, 1), (1, 0), (-1, 1)]
        """
        return self._path_dict.keys()

    def cluster_variables_so_far(self): # READY
        r"""
        Return a list of the cluster variables encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.current_seed().mutate(0)
            sage: A.cluster_variables_so_far()
            [x1, x0, (x1 + 1)/x0]
        """
        return map(self.cluster_variable, self.g_vectors_so_far())

    def F_polynomials_so_far(self): # READY
        r"""
        Return a list of the cluster variables encountered so far.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.current_seed().mutate(0)
            sage: A.F_polynomials_so_far()
            [1, 1, u0 + 1]
        """
        return self._F_poly_dict.values()

    def F_polynomial(self, g_vector):   # READY
        r"""
        Return the F-polynomial with g-vector ``g_vector`` if it has been found.

        INPUT:

        - ``g_vector`` -- a tuple: the g-vector of the F-polynomial to return.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.F_polynomial((-1, 1))
            Traceback (most recent call last):
            ...
            KeyError: 'The g-vector (-1, 1) has not been found yet.'
            sage: A.initial_seed().mutate(0,mutating_F=False)
            sage: A.F_polynomial((-1, 1))
            Traceback (most recent call last):
            ...
            KeyError: 'The F-polynomial with g-vector (-1, 1) has not been computed yet.
            You can compute it by mutating from the initial seed along the sequence [0].'
            sage: A.initial_seed().mutate(0)
            sage: A.F_polynomial((-1, 1))
            u0 + 1
        """
        g_vector = tuple(g_vector)
        try:
            return self._F_poly_dict[g_vector]
        except KeyError:
            if g_vector in self._path_dict:
                msg = "The F-polynomial with g-vector %s has not been computed yet. "%str(g_vector)
                msg += "You can compute it by mutating from the initial seed along the sequence "
                msg += str(self._path_dict[g_vector]) + "."
                raise KeyError(msg)
            else:
                raise KeyError("The g-vector %s has not been found yet."%str(g_vector))

    @cached_method(key=lambda a,b: tuple(b) )
    def cluster_variable(self, g_vector):   # READY
        r"""
        Return the cluster variable with g-vector ``g_vector`` if it has been found.

        INPUT:

        - ``g_vector`` -- a tuple: the g-vector of the cluster variable to return.

        ALGORITHM:

            This function computes cluster variables from their g-vectors and
            and F-polynomials using the "separation of additions" formula of
            Theorem 3.7 in [FZ07]_.

        EXAMPLE::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.initial_seed().mutate(0)
            sage: A.cluster_variable((-1,1))
            (x1 + 1)/x0
        """
        g_vector = tuple(g_vector)
        F = self.F_polynomial(g_vector)
        F_std = F.subs(self._yhat)
        g_mon = prod([self.ambient().gen(i)**g_vector[i] for i in xrange(self.rk())])
        F_trop = self.ambient()(F.subs(self._y))._fraction_pair()[1]
        return self.retract(g_mon*F_std*F_trop)

    def find_g_vector(self, g_vector, depth=infinity):  # READY
        r"""
        Return a mutation sequence to obtain a seed containing the g-vector ``g_vector`` from the initial seed.

        INPUT:

        - ``g_vector`` -- a tuple: the g-vector to find.

        - ``depth`` -- a positive integer: the maximum distance from ``self.current_seed`` to reach.

        OUTPUT:

            This function returns a list of integers if it can find
            ``g_vector`` otherwise it returns ``None``.  If the exploring
            iterator stops it means that the algebra is of finite type and
            ``g_vector`` is not the g-vector of any cluster variable. In this
            case the fuction resets the iterator and raises an error.

        EXAMPLES::

            sage: A = ClusterAlgebra(['G',2],principal_coefficients=True)
            sage: A.find_g_vector((-2, 3), depth=2)
            sage: A.find_g_vector((-2, 3), depth=3)
            [0, 1, 0]
        """
        g_vector = tuple(g_vector)
        while g_vector not in self.g_vectors_so_far() and self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth()
            except StopIteration:
                # Unless self._sd_iter has been manually altered, we checked
                # all the seeds of self and did not find g_vector.
                # Do some house cleaning before failing
                self.reset_exploring_iterator()
                raise ValueError("%s is not the g-vector of any cluster variable of %s."%(str(g_vector),str(self)))
        return copy(self._path_dict.get(g_vector,None))

    def ambient(self):  # READY
        r"""
        Return the Laurent polynomial ring containing ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.ambient()
            Multivariate Laurent Polynomial Ring in x0, x1, y0, y1 over Integer Ring
        """
        return self._ambient

    def scalars(self):  # READY
        r"""
        Return the scalars on which ``self`` is defined.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.scalars()
            Integer Ring
        """
        return self.base().base()

    def lift(self, x):  # READY
        r"""
        Return ``x`` as an element of :meth:`ambient`.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: x = A.cluster_variable((1,0))
            sage: A.lift(x).parent()
            Multivariate Laurent Polynomial Ring in x0, x1, y0, y1 over Integer Ring
        """
        return self.ambient()(x.value)

    def retract(self, x):   # READY
        r"""
        Return ``x`` as an element of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: L = A.ambient()
            sage: x = L.gen(0)
            sage: A.retract(x).parent()
            A Cluster Algebra with cluster variables x0, x1 and coefficients y0, y1 over Integer Ring
        """
        return self(x)

    def gens(self): # READY
        r"""
        Return the list of initial cluster variables and coefficients of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.gens()
            [x0, x1, y0, y1]
        """
        return map(self.retract, self.ambient().gens())

    def coefficient(self, j):   # READY
        r"""
        Return the j-th coefficient of ``self``.

        INPUT:

        - ``j`` -- an integer: the index of the coefficient to return.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.coefficient(0)
            y0
        """
        if isinstance(self.base(), LaurentPolynomialRing_generic):
            return self.retract(self.base().gen(j))
        else:
            raise ValueError("generator not defined")

    def coefficients(self): # READY
        r"""
        Return the list of coefficients of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.coefficients()
            [y0, y1]
        """
        if isinstance(self.base(), LaurentPolynomialRing_generic):
            return map(self.retract, self.base().gens())
        else:
            return []

    def coefficient_names(self):    # READY
        r"""
        Return the list of coefficient names.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: A.coefficient_names()
            ()
            sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
            sage: A.coefficient_names()
            ('y0', 'y1')
        """
        return self.variable_names()[self.rk():]

    def initial_cluster_variable(self, j):
        r"""
        Return the j-th initial cluster variable of ``self``.

        INPUT:

        - ``j`` -- an integer in ``range(self.parent().rk())``

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.initial_cluster_variable(0)
            x0
        """
        return self.retract(self.ambient().gen(j))

    def initial_cluster_variables(self):    # READY
        r"""
        Return the list of initial cluster variables of ``self``.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2],principal_coefficients=True)
            sage: A.initial_cluster_variables()
            [x0, x1]
        """
        return map(self.retract, self.ambient().gens()[:self.rk()])

    def initial_cluster_variable_names(self):   # READY
        r"""
        Return the list of initial variable names.

        EXAMPLES::

            sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
            sage: A.initial_cluster_variable_names()
            ('x0', 'x1')
        """
        return self.variable_names()[:self.rk()]

    def seeds(self, **kwargs):  # READY
        r"""
        Return an iterator running over seeds of ``self``.

        INPUT:

        - ``from_current_seed`` -- bool (default False): whether to start the
          iterator from :meth:`current_seed` or :meth:`initial_seed`.

        - ``mutating_F`` -- bool (default True): wheter to compute also
          F-polynomials; for speed considerations you may want to disable this.

        - ``allowed_directions`` -- a tuple of integers (default
          ``range(self.rk())``): the directions in which to mutate.

        - ``depth`` -- (defaulf ``infinity``):  the maximum depth at which to stop searching.

        ALGORITHM:

            This function traverses the exchange graph in a breadth-first search.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',4])
            sage: seeds = A.seeds(allowed_directions=[3,0,1])
            sage: _ = list(seeds)
            sage: A.g_vectors_so_far()
            [(-1, 0, 0, 0),
             (1, 0, 0, 0),
             (0, 0, 0, 1),
             (0, -1, 0, 0),
             (0, 0, 1, 0),
             (0, 1, 0, 0),
             (-1, 1, 0, 0),
             (0, 0, 0, -1)]
        """
        # should we begin from the current seed?
        if kwargs.get('from_current_seed', False):
            seed = copy(self.current_seed())
        else:
            seed = self.initial_seed()

        # yield first seed
        yield seed

        # some initialization
        depth_counter = 0
        n = self.rk()

        # do we mutate F-polynomials?
        mutating_F = kwargs.get('mutating_F', True)

        # which directions are we allowed to mutate into
        allowed_dirs = list(sorted(kwargs.get('allowed_directions', range(n))))

        # setup seeds storage
        cl = frozenset(seed.g_vectors())
        clusters = {}
        clusters[cl] = [ seed, copy(allowed_dirs) ]

        # ready, set, go!
        gets_bigger = True
        while gets_bigger and depth_counter < kwargs.get('depth', infinity):
            # remember if we got a new seed
            gets_bigger = False

            for key in clusters.keys():
                sd, directions = clusters[key]
                while directions:
                    # we can mutate in some direction
                    i = directions.pop()
                    new_sd  = sd.mutate(i, inplace=False, mutating_F=mutating_F)
                    new_cl = frozenset(new_sd.g_vectors())
                    if new_cl in clusters:
                        # we already had new_sd, make sure it does not mutate to sd during next round
                        j = map(tuple,clusters[new_cl][0].g_vectors()).index(new_sd.g_vector(i))
                        try:
                            clusters[new_cl][1].remove(j)
                        except ValueError:
                            pass
                    else:
                        # we got a new seed
                        gets_bigger = True
                        # next round do not mutate back to sd and do commuting mutations only in directions j > i
                        new_directions = [ j for j in allowed_dirs if j > i or new_sd.b_matrix()[j,i] != 0 ]
                        clusters[new_cl] = [ new_sd, new_directions ]
                        yield new_sd
            # we went one step deeper
            depth_counter += 1

    def reset_exploring_iterator(self, mutating_F=True):    # READY
        r"""
        Reset the iterator used to explore ``self``.

        INPUT:

        - ``mutating_F`` -- bool (default True): whether to also compute
          F-polynomials; for speed considerations you may want to disable this.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',4])
            sage: A.reset_exploring_iterator(mutating_F=False)
            sage: A.explore_to_depth(infinity)
            sage: len(A.g_vectors_so_far())
            14
            sage: len(A.F_polynomials_so_far())
            4
        """
        self._sd_iter = self.seeds(mutating_F=mutating_F)
        self._explored_depth = 0

    def explore_to_depth(self, depth):  # READY
        r"""
        Explore the exchange graph of ``self`` up to distance ``depth`` from the initial seed.

        INPUT:

        - ``depth`` -- the maximum depth at which to stop searching.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',4])
            sage: A.explore_to_depth(infinity)
            sage: len(A.g_vectors_so_far())
            14
        """
        while self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth()
            except:
                break

    def cluster_fan(self, depth=infinity):
        r"""
        Return the cluster fan (the fan of g-vectors) of ``self``.

        INPUT:

        - ``depth`` -- (default ``infinity``): the maximum depth at which to compute.

        EXAMPLES::

            sage: A = ClusterAlgebra(['A',2])
            sage: A.cluster_fan()
            Rational polyhedral fan in 2-d lattice N
        """
        seeds = self.seeds(depth=depth, mutating_F=False)
        cones = map(lambda s: Cone(s.g_vectors()), seeds)
        return Fan(cones)

    @_mutation_parse
    def mutate_initial(self, k):
        r"""
        Return the cluster algebra obtained by mutating ``self`` at the initial seed.

        INPUT:

        - ``k`` --  an integer in ``range(self.parent().rk())``: the direction
          in which we are mutating

        ALGORITHM:

            This function computes data for the new algebra from known data for
            the old algebra using [NZ12]_ equation (4.2) for g-vectors, and
            [FZ07]_ equation (6.21) for F-polynomials. The exponent ``h`` in the
            formula for F-polynomials is ``-min(0,old_g_vect[k])`` due to [NZ12]_
            Proposition 4.2.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: A.explore_to_depth(infinity)
            sage: B = A.b_matrix()
            sage: B.mutate(0)
            sage: A1 = ClusterAlgebra(B)
            sage: A1.explore_to_depth(infinity)
            sage: A2 = A1.mutate_initial(0)
            sage: A2._F_poly_dict == A._F_poly_dict
            True
        """
        n = self.rk()
        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction ' + str(k) + '.')

        # save computed data
        old_F_poly_dict = copy(self._F_poly_dict)
        old_path_dict = copy(self._path_dict)
        old_path_to_current =  copy(self.current_seed().path_from_initial_seed())

        # mutate initial exchange matrix
        B0 = copy(self._B0)
        B0.mutate(k)

        # HACK: pretend there is no embedding of self into self.parent(): we
        # will recreate this in a moment.  The problem is that there can be
        # only one embedding per object and it is set up in __init__.
        self._unset_embedding()

        # update algebra
        cv_names = self.initial_cluster_variable_names()
        coeff_names = self.coefficient_names()
        scalars = self.scalars()
        self.__init__(B0, cluster_variable_names=cv_names,
                coefficient_names=coeff_names, scalars=scalars)

        # substitution data to compute new F-polynomials
        Ugen = self._U.gens()
        # here we have \mp B0 rather then \pm B0 because we want the k-th row of the old B0
        F_subs = [Ugen[k]**(-1) if j==k else Ugen[j]*Ugen[k]**max(B0[k,j],0)*(1+Ugen[k])**(-B0[k,j]) for j in xrange(n)]

        # restore computed data
        for old_g_vect in old_path_dict:
            # compute new g-vector
            J = identity_matrix(n)
            eps = sign(old_g_vect[k])
            for j in xrange(n):
                # here we have -eps*B0 rather than eps*B0 because we want the k-th column of the old B0
                J[j,k] += max(0, -eps*B0[j,k])
            J[k,k] = -1
            new_g_vect = tuple(J*vector(old_g_vect))

            #compute new path
            new_path = old_path_dict[old_g_vect]
            new_path = ([k]+new_path[:1] if new_path[:1] != [k] else []) + new_path[1:]
            self._path_dict[new_g_vect] = new_path

            #compute new F-polynomial
            if old_F_poly_dict.has_key(old_g_vect):
                h =  -min(0,old_g_vect[k])
                new_F_poly = old_F_poly_dict[old_g_vect](F_subs)*Ugen[k]**h*(Ugen[k]+1)**old_g_vect[k]
                self._F_poly_dict[new_g_vect] = new_F_poly

        # reset self.current_seed() to the previous location
        S = self.initial_seed()
        S.mutate([k]+old_path_to_current, mutating_F=False)
        self.set_current_seed(S)

    # DESIDERATA
    # Some of these are probably unrealistic
    def upper_cluster_algebra(self):
        r"""
        Return the upper cluster algebra associated to ``self``.
        """
        raise NotImplementedError("Not implemented yet.")

    def upper_bound(self):
        r"""
        Return the upper bound associated to ``self``.
        """
        raise NotImplementedError("Not implemented yet.")

    def lower_bound(self):
        r"""
        Return the lower bound associated to ``self``.
        """
        raise NotImplementedError("Not implemented yet.")

    def theta_basis_element(self, g_vector):
        r"""
        Return the element of the theta basis with g-vector ``g_vetor``.
        """
        raise NotImplementedError("Not implemented yet.")

####
# Methods only defined for special cases
####

def greedy_element(self, d_vector):
    r"""
    Return the greedy element with d-vector ``d_vector``.

    INPUT

    - ``d_vector`` -- tuple of 2 integers: the d-vector of the element to compute.

    ALGORITHM:

        This implements greedy elements of a rank 2 cluster algebra from [LLZ14]_ equation (1.5).

    REFERENCES:

    .. [LLZ14] \K. Lee, \L. Li, and \A. Zelevinsky, "Greedy elements in rank 2
       cluster algebras", Selecta Math. 20 (2014), 57-82.

    EXAMPLES::

        sage: A = ClusterAlgebra(['A',[1,1],1])
        sage: A.greedy_element((1,1))
        (x0^2 + x1^2 + 1)/(x0*x1)
    """
    b = abs(self.b_matrix()[0,1])
    c = abs(self.b_matrix()[1,0])
    (a1, a2) = d_vector
    # here we use the generators of self.ambient() because cluster variables do not have an inverse.
    (x1, x2) = self.ambient().gens()
    if a1 < 0:
        if a2 < 0:
            return self.retract(x1**(-a1) * x2**(-a2))
        else:
            return self.retract(x1**(-a1) * ((1+x2**c)/x1)**a2)
    elif a2 < 0:
        return self.retract(((1+x1**b)/x2)**a1 * x2**(-a2))
    output = 0
    for p in xrange(0, a2+1):
        for q in xrange(0, a1+1):
            output += self._greedy_coefficient(d_vector, p, q) * x1**(b*p) * x2**(c*q)
    return self.retract(x1**(-a1) * x2**(-a2) * output)

def _greedy_coefficient(self,d_vector,p,q): # READY
    r"""
    Return the coefficient of the monomial ``x1**(b*p) * x2**(c*q)`` in the numerator of the greedy element with denominator vector ``d_vector``.

    EXAMPLES::

        sage: A = ClusterAlgebra(['A',[1,1],1])
        sage: A.greedy_element((1,1))
        (x0^2 + x1^2 + 1)/(x0*x1)
        sage: A._greedy_coefficient((1,1),0,0)
        1
        sage: A._greedy_coefficient((1,1),1,0)
        1
    """
    b = abs(self.b_matrix()[0,1])
    c = abs(self.b_matrix()[1,0])
    (a1, a2) = d_vector
    p = Integer(p)
    q = Integer(q)
    if p == 0 and q == 0:
        return Integer(1)
    sum1 = 0
    for k in range(1,p+1):
        bino = 0
        if a2 - c*q + k - 1 >= k:
            bino = binomial(a2 - c*q + k - 1, k)
        sum1 += (-1)**(k-1) * self._greedy_coefficient(d_vector, p-k, q) * bino
    sum2 = 0
    for l in range(1, q+1):
        bino = 0
        if a1 - b*p + l - 1 >= l:
            bino = binomial(a1 - b*p + l - 1, l)
        sum2 += (-1)**(l-1) * self._greedy_coefficient(d_vector, p, q-l) * bino
    return Integer(max(sum1, sum2))
