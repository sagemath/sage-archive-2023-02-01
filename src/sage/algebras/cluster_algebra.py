r"""
Cluster algebras

Implementation of cluster algebras as an algebra using mainly structural theorems from CA IV

TODO: We should write a nice paragraph here.

AUTHORS:

- Dylan Rupel (2015-06-15): initial version

- Salvatore Stella (2015-06-15): initial version

EXAMPLES::

    TODO: we need to write a complete example of usage here
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

# TODO: check that we import all we need and possibly move some import used
# rarely close to where needed
from copy import copy
from functools import wraps
from sage.categories.homset import Hom
from sage.categories.morphism import SetMorphism
from sage.categories.quotient_fields import QuotientFields
from sage.categories.rings import Rings
from sage.combinat.cluster_algebra_quiver.quiver import ClusterQuiver
from sage.functions.generalized import sign
from sage.functions.other import binomial
from sage.matrix.constructor import identity_matrix
from sage.matrix.special import block_matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.modules.free_module_element import vector
from sage.rings.infinity import infinity
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject
from types import MethodType

##############################################################################
# Helper functions
##############################################################################
def mutation_parse(mutate): # READY
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
    def mutate_wrapper(self, direction, *args, **kwargs):
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
            mutate(to_mutate, k, *args, **kwargs)

        if not inplace:
            return to_mutate

    return mutate_wrapper

##############################################################################
# Elements of a cluster algebra
##############################################################################

class ClusterAlgebraElement(ElementWrapper):    # READY

    def __init__(self, parent, value):  # READY
        r"""
        An element of a cluster algebra.

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
    def _add_(self, other): # READY
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

    def _div_(self, other): # READY
        r"""
        Return the quotient of ``self`` and ``other``.

        WARNING::
            
            This method returns an element of ``self.parent().ambient()``
            rather than an element of ``self.parent()`` because, a priori,
            we cannot guarantee membership.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: A.an_element() / A.an_element()
            1
        """
        return self.lift()/other.lift()

    def d_vector(self): # READY
        r"""
        Return the d-vector of ``self`` as a tuple of integers.

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

    def _repr_(self):   # READY
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

def g_vector(self): #READY
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

def is_homogeneous(self):   #READY
    r"""
    Return ``True`` if ``self`` is an homogeneous element of ``self.parent()``.
    
    EXAMPLES::
    
        sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
        sage: x = A.cluster_variable((1,0)) + A.cluster_variable((0,1))
        sage: x.is_homogeneous()
        False
    """
    return len(self.homogeneous_components()) == 1

def homogeneous_components(self):   #READY
    r"""
    Return a dictionary of the homogeneous components of ``self''.

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

    def __init__(self, B, C, G, parent, **kwargs):  # READY
        r"""
        A seed in a cluster algebra.

        INPUT:
        
        - ``B`` -- a skew-symmetrizable integer matrix;

        - ``C`` -- the matrix of c-vectors of ``self``;

        - ``G`` -- the matrix of g-vectors of ``self``;

        - ``parent`` -- a :class:`ClusterAlgebra`: the algebra to which the
          seed belongs;

        - ``path`` -- list (default: []) the mutation sequence from the initial
          seed of ``parent`` to `self``

        WARNING:
            
            Seeds should **not** be created manually: no test is performed to
            assert that they are built from consistent data nor that they
            really are seeds of ``parent``. If you create seeds with
            inconsistent data all sort of things can go wrong, even
            :meth:`__eq__` is no longer guaranteed to give correct answers. 
            Use at your ouwn risk.

        EXAMPLES::

            sage: A = ClusterAlgebra(['F',4])
            sage: from sage.algebras.cluster_algebra import ClusterAlgebraSeed
            sage: ClusterAlgebraSeed(A.b_matrix(),identity_matrix(4),identity_matrix(4),A,path=[1,2,3])
            The seed of Cluster Algebra of rank 4 obtained from the initial by mutating along the sequence [1, 2, 3]

        """
        self._B = copy(B)
        self._C = copy(C)
        self._G = copy(G)
        self._parent = parent
        self._path = kwargs.get('path', [])

    def __copy__(self): # READY
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

    def __eq__(self, other):    # READY
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
        """
        return self.parent() == other.parent() and frozenset(self.g_vectors()) == frozenset(other.g_vectors())

    def _repr_(self):   # READY
        r"""
        Return the string representation of ``self``.
 
        EXAMPLES::

            sage: A = ClusterAlgebra(['A',3])
            sage: S = A.current_seed(); S
            The initial seed of Cluster Algebra of rank 3
            sage: S.mutate(0); S
            The seed of Cluster Algebra of rank 3 obtained from the initial by mutating in direction 0
            sage: S.mutate(1); S
            The seed of Cluster Algebra of rank 3 obtained from the initial by mutating along the sequence [0, 1]
        """
        if self._path == []:
            return "The initial seed of %s"%str(self.parent())
        elif len(self._path) == 1:
            return "The seed of %s obtained from the initial by mutating in direction %s"%(str(self.parent()),str(self._path[0]))
        else:
            return "The seed of %s obtained from the initial by mutating along the sequence %s"%(str(self.parent()),str(self._path))

    def depth(self):    # READY
        r"""
        Retun the length of a path from the initial seed of :meth:`parent` to ``self``.

        WARNING::
            This is the length of the path returned by
            :meth:`path_from_initial_seed` which needs not be the shortest
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

    def parent(self):
        r"""
        Return the parent of ``self``.

        EXAMPLES::
            sage: A = ClusterAlgebra(['B',3])
            sage: A.current_seed().parent() == A
            True
        """
        return self._parent

    def b_matrix(self):
        return copy(self._B)

    def c_matrix(self):
        return copy(self._C)

    def c_vector(self, j):
        return tuple(self._C.column(j))

    def c_vectors(self):
        return map(tuple, self._C.columns())

    def g_matrix(self):
        return copy(self._G)

    def g_vector(self, j):
        return tuple(self._G.column(j))

    def g_vectors(self):
        return map(tuple, self._G.columns())

    def F_polynomial(self, j):
        return self.parent().F_polynomial(self.g_vector(j))

    def F_polynomials(self):
        return (self.parent().F_polynomial(g) for g in self.g_vectors())

    def cluster_variable(self, j):
        return self.parent().cluster_variable(self.g_vector(j))

    def cluster_variables(self):
        return (self.parent().cluster_variable(g) for g in self.g_vectors())

    @mutation_parse
    def mutate(self, k, mutating_F=True):
        r"""
        mutate seed
        
        INPUT:

        bla bla ba
        """
        n = self.parent().rk()

        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction ' + str(k) + '.')

        # store mutation path
        if self._path != [] and self._path[-1] == k:
            self._path.pop()
        else:
            self._path.append(k)

        # find sign of k-th c-vector
        # Will this be used enough that it should be a built-in function?
        if any(x > 0 for x in self._C.column(k)):
            eps = +1
        else:
            eps = -1

        # store the g-vector to be mutated in case we are mutating F-polynomials also
        old_g_vector = self.g_vector(k)

        # G-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[j,k] += max(0, -eps*self._B[j,k])
        J[k,k] = -1
        self._G = self._G*J

        # g-vector path list
        g_vector = self.g_vector(k)
        if g_vector not in self.parent().g_vectors_so_far():
            self.parent()._path_dict[g_vector] = copy(self._path)
            # F-polynomials
            if mutating_F:
                self.parent()._F_poly_dict[g_vector] = self._mutated_F(k, old_g_vector)

        # C-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[k,j] += max(0, eps*self._B[k,j])
        J[k,k] = -1
        self._C = self._C*J

        # B-matrix
        self._B.mutate(k)

    def _mutated_F(self, k, old_g_vector):
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
        # TODO: understand why using // instead of / here slows the code down by
        # a factor of 3 but in the original experiments we made at sage days it
        # was much faster with // (we were working with cluter variables at the
        # time).
        # By the way, as of August 28th 2015 we split in half the running time
        # compared to sage days. With my laptop plugged in I get
        # sage: A = ClusterAlgebra(['E',8])
        # sage: seeds = A.seeds()
        # sage: %time void = list(seeds)
        # CPU times: user 26.8 s, sys: 21 ms, total: 26.9 s
        # Wall time: 26.8 s
        #####
        # Bad news: as of 19/10/2015 we got a huge slowdown:
        # right now it takes 150s with / and 100s with //
        # what did we do wrong?
        ##
        # I figured it out: the problem is that casting the result to alg._U is
        # quite slow: it amounts to run // instead of / :(
        # do we need to perform this or can we be less precise here and allow
        # F-polynomials to be rational funtions?
        # I am partucularly unhappy about this, for the moment the correct and
        # slow code is commented
        #return alg._U((pos+neg)/alg.F_polynomial(old_g_vector))
        ##
        # One more comment: apparently even without casting the result is a
        # polynomial! This is really weird but I am not going to complain. I
        # suppose we should not do the casting then

        # DR: Now I get the same computation time for / and //, 49.7s while simultaneously rebuiling sage
        return (pos+neg)/alg.F_polynomial(old_g_vector)

    def path_from_initial_seed(self):
        return copy(self._path)

    # TODO: ideally we should allow to mutate in direction "this g-vector" or
    # "this cluster variable" or "sink", "urban renewal" and all the other
    # options provided by Gregg et al. To do so I guess the best option is to
    # have a generic function transforming all these into an index and use it as
    # a decorator. In this way we can also use it in this __contains__ even
    # though one may write weird things like "sink" in A.current_seed and get
    # True as an answer.
    def __contains__(self, element):
        if isinstance(element, ClusterAlgebraElement):
            cluster = [ self.cluster_variable(i) for i in xrange(self.parent().rk()) ]
        else:
            element = tuple(element)
            cluster = self.g_vectors()
        return element in cluster

##############################################################################
# Cluster algebras
##############################################################################

class ClusterAlgebra(Parent):
    r"""
    INPUT:

    - ``data`` -- some data defining a cluster algebra.

    - ``scalars`` -- (default ZZ) the scalars on which the cluster algebra
      is defined.

    - ``cluster_variables_prefix`` -- string (default 'x').

    - ``cluster_variables_names`` -- a list of strings.  Superseedes
      ``cluster_variables_prefix``.

    - ``coefficients_prefix`` -- string (default 'y').

    - ``coefficients_names`` -- a list of strings. Superseedes
      ``cluster_variables_prefix``.

    - ``principal_coefficients`` -- bool (default: False). Superseedes any
      coefficient defined by ``data``.
    """

    Element = ClusterAlgebraElement

    def __init__(self, data, **kwargs):
        r"""
        See :class:`ClusterAlgebra` for full documentation.
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
        variables_prefix = kwargs.get('cluster_variables_prefix','x')
        variables = kwargs.get('cluster_variables_names', [variables_prefix+str(i) for i in xrange(n)])
        if len(variables) != n:
             raise ValueError("cluster_variables_names should be a list of %d valid variable names"%n)

        # Determine scalars
        scalars = kwargs.get('scalars', ZZ)

        # Determine coefficients and setup self._base
        if m>0:
            coefficients_prefix = kwargs.get('coefficients_prefix', 'y')
            if coefficients_prefix == variables_prefix:
                offset = n
            else:
                offset = 0
            coefficients = kwargs.get('coefficients_names', [coefficients_prefix+str(i) for i in xrange(offset,m+offset)])
            if len(coefficients) != m:
                raise ValueError("coefficients_names should be a list of %d valid variable names"%m)
            base = LaurentPolynomialRing(scalars, coefficients)
        else:
            base = scalars
            coefficients = []

        # setup Parent and ambient
        self._ambient = LaurentPolynomialRing(scalars, variables+coefficients)
        Parent.__init__(self, base=base, category=Rings(scalars).Commutative().Subobjects(), names=variables+coefficients)

        # Data to compute cluster variables using separation of additions
        # BUG WORKAROUND: if your sage installation does not have trac:`19538`
        # merged uncomment the following line and comment the next
        #self._y = dict([ (self._U.gen(j), prod([self._ambient.gen(n+i)**M0[i,j] for i in xrange(m)])) for j in xrange(n)])
        self._y = dict([ (self._U.gen(j), prod([self._base.gen(i)**M0[i,j] for i in xrange(m)])) for j in xrange(n)])
        self._yhat = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(n+m)])) for j in xrange(n)])

        # Have we got principal coefficients?
        self._is_principal = (M0 == I)

        # Store initial data
        self._B0 = copy(B0)
        self._n = n
        self.reset_current_seed()

        # Internal data for exploring the exchange graph
        self.reset_exploring_iterator()

        # Add methods that are defined only for special cases
        if n == 2 and m == 0:
            self.greedy_element = MethodType(greedy_element, self, self.__class__)
            self.greedy_coefficient = MethodType(greedy_coefficient, self, self.__class__)

        # Register embedding into self.ambient()
        embedding = SetMorphism(Hom(self,self.ambient()), lambda x: x.lift())
        self._populate_coercion_lists_(embedding=embedding)

    def __copy__(self):
        other = type(self).__new__(type(self))
        other._U = self._U
        other._path_dict = copy(self._path_dict)
        other._F_poly_dict = copy(self._F_poly_dict)
        other._ambient = self._ambient
        other._y = copy(self._y)
        other._yhat = copy(self._yhat)
        other._is_principal = self._is_principal
        other._B0 = copy(self._B0)
        other._n = self._n
        other._seed = copy(self._seed)
        # We probably need to put n=2 initializations here also
        # TODO: we may want to use __init__ to make the initialization somewhat easier (say to enable special cases) This might require a better written __init__
        return other

    def __eq__(self, other):
        return type(self) == type(other) and self._B0 == other._B0 and  self._yhat == other._yhat and self.base() == other.base()

    def _coerce_map_from_(self, other):
        # if other is a cluster algebra allow inherit coercions from ambients
        if isinstance(other, ClusterAlgebra):
            return self.ambient().has_coerce_map_from(other.ambient())
        # everything that is in the base can be coerced to self
        return self.base().has_coerce_map_from(other)

    def _repr_(self):
        return "Cluster Algebra of rank %s"%self.rk()

    def _an_element_(self):
        return self.current_seed().cluster_variable(0)

    def rk(self):
        r"""
        The rank of ``self`` i.e. the number of cluster variables in any seed of
        ``self``.
        """
        return self._n

    def current_seed(self):
        r"""
        The current seed of ``self``.
        """
        return self._seed

    def set_current_seed(self, seed):
        r"""
        Set ``self._seed`` to ``seed`` if it makes sense.
        """
        if self.contains_seed(seed):
            self._seed = seed
        else:
            raise ValueError("This is not a seed in this cluster algebra.")

    def contains_seed(self, seed):
        computed_sd = self.initial_seed()
        computed_sd.mutate(seed._path, mutating_F=False)
        return computed_sd == seed

    def reset_current_seed(self):
        r"""
        Reset the current seed to the initial one
        """
        self._seed = self.initial_seed()

    def initial_seed(self):
        r"""
        Return the initial seed
        """
        n = self.rk()
        I = identity_matrix(n)
        return ClusterAlgebraSeed(self.b_matrix(), I, I, self)

    def b_matrix(self): 
        r"""
        Return the initial exchange matrix of ``self``.
        """
        n = self.rk()
        return copy(self._B0[:n,:])

    def g_vectors_so_far(self):
        r"""
        Return the g-vectors of cluster variables encountered so far.
        """
        return self._path_dict.keys()

    def F_polynomial(self, g_vector):
        g_vector= tuple(g_vector)
        try:
            return self._F_poly_dict[g_vector]
        except:
            # If the path is known, should this method perform that sequence of mutations to compute the desired F-polynomial?
            # Yes, perhaps with the a prompt first, something like:
            #comp = raw_input("This F-polynomial has not been computed yet.  It can be found using %s mutations.  Continue? (y or n):"%str(directions.__len__()))
            #if comp == 'y':
            #    ...compute the F-polynomial...
            if g_vector in self._path_dict:
                raise ValueError("The F-polynomial with g-vector %s has not been computed yet.  You probably explored the exchange tree with compute_F=False.  You can compute this F-polynomial by mutating from the initial seed along the sequence %s."%(str(g_vector),str(self._path_dict[g_vector])))
            else:
                raise ValueError("The F-polynomial with g-vector %s has not been computed yet."%str(g_vector))

    @cached_method(key=lambda a,b: tuple(b) )
    def cluster_variable(self, g_vector):
        g_vector = tuple(g_vector)
        if not g_vector in self.g_vectors_so_far():
            # Should we let the self.F_polynomial below handle raising the exception?
            raise ValueError("This Cluster Variable has not been computed yet.")
        F_std = self.F_polynomial(g_vector).subs(self._yhat)
        g_mon = prod([self.ambient().gen(i)**g_vector[i] for i in xrange(self.rk())])
        # LaurentPolynomial_mpair does not know how to compute denominators, we need to lift to its fraction field
        F_trop = self.ambient()(self.F_polynomial(g_vector).subs(self._y))._fraction_pair()[1]
        return self.retract(g_mon*F_std*F_trop)

    def find_cluster_variable(self, g_vector, depth=infinity):
        r"""
        Returns the shortest mutation path to obtain the cluster variable with
        g-vector ``g_vector`` from the initial seed.

        ``depth``: maximum distance from ``self.current_seed`` to reach.

        WARNING: if this method is interrupted then ``self._sd_iter`` is left in
        an unusable state. To use again this method it is then necessary to
        reset ``self._sd_iter`` via self.reset_exploring_iterator()
        """
        g_vector = tuple(g_vector)
        mutation_counter = 0
        while g_vector not in self.g_vectors_so_far() and self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth()
            except:
                raise ValueError("Could not find a cluster variable with g-vector %s up to mutation depth %s after performing %s mutations."%(str(g_vector),str(depth),str(mutation_counter)))

            # If there was a way to have the seeds iterator continue after the depth_counter reaches depth,
            # the following code would allow the user to continue searching the exchange graph
            #cont = raw_input("Could not find a cluster variable with g-vector %s up to mutation depth %s."%(str(g_vector),str(depth))+"  Continue searching? (y or n):")
            #if cont == 'y':
            #    new_depth = 0
            #    while new_depth <= depth:
            #        new_depth = raw_input("Please enter a new mutation search depth greater than %s:"%str(depth))
            #    seeds.send(new_depth)
            #else:
            #    raise ValueError("Could not find a cluster variable with g-vector %s after %s mutations."%(str(g_vector),str(mutation_counter)))

            mutation_counter += 1
        return copy(self._path_dict[g_vector])

    def ambient(self):
        return self._ambient

    def lift(self, x):
        r"""
        Return x as an element of self._ambient
        """
        return self.ambient()(x.value)

    def retract(self, x):
        return self(x)

    def gens(self):
        r"""
        Return the generators of :meth:`self.ambient`
        """
        return map(self.retract, self.ambient().gens())

    def seeds(self, depth=infinity, mutating_F=True, from_current_seed=False):
        r"""
        Return an iterator producing all seeds of ``self`` up to distance
        ``depth`` from ``self.initial_seed`` or ``self.current_seed``.

        If ``mutating_F`` is set to false it does not compute F_polynomials
        """
        if from_current_seed:
            seed = self.current_seed()
        else:
            seed = self.initial_seed()
        # add allowed_directions

        yield seed
        depth_counter = 0
        n = self.rk()
        cl = frozenset(seed.g_vectors())
        clusters = {}
        clusters[cl] = [ seed, range(n) ]
        gets_bigger = True
        while gets_bigger and depth_counter < depth:
            gets_bigger = False
            keys = clusters.keys()
            for key in keys:
                sd, directions = clusters[key]
                while directions:
                    i = directions.pop()
                    new_sd  = sd.mutate(i, inplace=False, mutating_F=mutating_F)
                    new_cl = frozenset(new_sd.g_vectors())
                    if new_cl in clusters:
                        j = map(tuple,clusters[new_cl][0].g_vectors()).index(new_sd.g_vector(i))
                        try:
                            clusters[new_cl][1].remove(j)
                        except:
                            pass
                    else:
                        gets_bigger = True
                        # doublecheck this way of producing directions for the new seed: it is taken almost verbatim fom ClusterSeed
                        new_directions = [ j for j in xrange(n) if j > i or new_sd.b_matrix()[j,i] != 0 ]
                        clusters[new_cl] = [ new_sd, new_directions ]
                        # Use this if we want to have the user pass info to the
                        # iterator
                        #new_depth = yield new_sd
                        #if new_depth > depth:
                        #    depth = new_depth
                        yield new_sd
            depth_counter += 1

    def reset_exploring_iterator(self, mutating_F=True):
        self._sd_iter = self.seeds(mutating_F=mutating_F)
        self._explored_depth = 0

    @mutation_parse
    def mutate_initial(self, k):
        # WARNING: at the moment this function does not behave well with respect to coefficients:
        # sage: A = ClusterAlgebra(['B',2],principal_coefficients=True)
        # sage: A.initial_seed().b_matrix()
        # [ 0  1]
        # [-2  0]
        # sage: A.initial_seed().c_matrix()
        # [1 0]
        # [0 1]
        # sage: A.mutate_initial(0)
        # sage: A.initial_seed().b_matrix()
        # [ 0 -1]
        # [ 2  0]
        # sage: A.initial_seed().c_matrix()
        # [1 0]
        # [0 1]
        # this creates several issues 
        # BTW: mutating the initial seed in place creates several issues with elements of the algebra: for example:
        # sage: A = ClusterAlgebra(['B',4],principal_coefficients=True)
        # sage: A.explore_to_depth(infinity)
        # sage: x = A.cluster_variable((-1, 1, -2, 2))
        # sage: x.is_homogeneous()
        # True
        # sage: A.mutate_initial([0,3])
        # sage: x in A
        # True
        # sage: (-1, 1, -2, 2) in A.g_vectors_so_far()
        # False
        # sage: x.is_homogeneous()
        # False

        r"""
        Mutate ``self`` in direction `k` at the initial cluster.

        INPUT:
        - ``k`` -- integer in between 0 and ``self.rk``
        """
        n = self.rk()

        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction %s, please try a value between 0 and %s.'%(str(k),str(n-1)))
        
        #store current seed location
        path_to_current = self.current_seed().path_from_initial_seed()
        #modify self._path_dict using Nakanishi-Zelevinsky (4.1) and self._F_poly_dict using CA-IV (6.21)
        new_path_dict = dict()
        new_F_dict = dict()
        new_path_dict[tuple(identity_matrix(n).column(k))] = []
        new_F_dict[tuple(identity_matrix(n).column(k))] = self._U(1)

        poly_ring = PolynomialRing(ZZ,'u')
        F_subs_tuple = tuple([self._U.gen(k)**(-1) if j==k else self._U.gen(j)*self._U.gen(k)**max(-self._B0[k][j],0)*(1+self._U.gen(k))**(self._B0[k][j]) for j in xrange(n)])

        for g_vect in self._path_dict:
            #compute new path
            path = self._path_dict[g_vect]
            if g_vect == tuple(identity_matrix(n).column(k)):
                new_path = [k]
            elif path != []:
                if path[0] != k:
                    new_path = [k] + path
                else:
                    new_path = path[1:]
            else:
                new_path = []

            #compute new g-vector
            new_g_vect = vector(g_vect) - 2*g_vect[k]*identity_matrix(n).column(k)
            for i in xrange(n):
                new_g_vect += max(sign(g_vect[k])*self._B0[i,k],0)*g_vect[k]*identity_matrix(n).column(i)
            new_path_dict[tuple(new_g_vect)] = new_path

            #compute new F-polynomial
            h =  -min(0,g_vect[k])
            new_F_dict[tuple(new_g_vect)] = self._F_poly_dict[g_vect](F_subs_tuple)*self._U.gen(k)**h*(self._U.gen(k)+1)**g_vect[k]

        self._path_dict = new_path_dict
        self._F_poly_dict = new_F_dict

        self._B0.mutate(k)
       
        # keep the current seed were it was on the exchange graph
        self._seed = self.initial_seed().mutate([k]+self.current_seed().path_from_initial_seed(), mutating_F=False, inplace=False)

    def explore_to_depth(self, depth):
        while self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth()
            except:
                break

    def cluster_fan(self, depth=infinity):
        from sage.geometry.cone import Cone
        from sage.geometry.fan import Fan
        seeds = self.seeds(depth=depth, mutating_F=False)
        cones = map(lambda s: Cone(s.g_vectors()), seeds)
        return Fan(cones)

    # DESIDERATA. Some of these are probably unrealistic
    def upper_cluster_algebra(self):
        pass

    def upper_bound(self):
        pass

    def lower_bound(self):
        pass

    def theta_basis_element(self, g_vector):
        pass

####
# Methods only defined for special cases
####

# Greedy elements exist only in rank 2
# Does not yet take into account coefficients, this can probably be done by
# using the greedy coefficients to write down the F-polynomials
def greedy_element(self, d_vector):
    b = abs(self._B0[0,1])
    c = abs(self._B0[1,0])
    a1 = d_vector[0]
    a2 = d_vector[1]
    # TODO: we need to have something like initial_cluster_variables so that we
    # do not have to use the generators of the ambient field. (this would also
    # make it better behaved when allowing different names)
    # Warning: there might be issues with coercions, make sure there are not
    x1 = self._ambient.gens()[0]
    x2 = self._ambient.gens()[1]
    if a1 < 0:
        if a2 < 0:
            return self.retract(x1**(-a1)*x2**(-a2))
        else:
            return self.retract(x1**(-a1)*((1+x2**c)/x1)**a2)
    elif a2 < 0:
        return self.retract(((1+x1**b)/x2)**a1*x2**(-a2))
    output = 0
    for p in xrange(0,a2+1):
        for q in xrange(0,a1+1):
            output += self.greedy_coefficient(d_vector,p,q)*x1**(b*p)*x2**(c*q)
    return self.retract(x1**(-a1)*x2**(-a2)*output)

# Is this function something we want to make public or do we want to make this a
# private method changing it to _greedy_coefficient ?
def greedy_coefficient(self,d_vector,p,q):
    b = abs(self._B0[0,1])
    c = abs(self._B0[1,0])
    a1 = d_vector[0]
    a2 = d_vector[1]
    p = Integer(p)
    q = Integer(q)
    if p == 0 and q == 0:
        return 1
    sum1 = 0
    for k in range(1,p+1):
        bin = 0
        if a2-c*q+k-1 >= k:
            bin = binomial(a2-c*q+k-1,k)
        sum1 += (-1)**(k-1)*self.greedy_coefficient(d_vector,p-k,q)*bin
    sum2 = 0
    for l in range(1,q+1):
        bin = 0
        if a1-b*p+l-1 >= l:
            bin = binomial(a1-b*p+l-1,l)
        sum2 += (-1)**(l-1)*self.greedy_coefficient(d_vector,p,q-l)*bin
    #print "sum1=",sum1,"sum2=",sum2
    return max(sum1,sum2)



