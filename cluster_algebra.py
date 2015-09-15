import time
from types import MethodType
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ

####
# Helper functions.
####
# These my have to go in a more standard place or just outside this file

####
# Elements of a cluster algebra
####

class ClusterAlgebraElement(ElementWrapper):

    # this is to store extra information like g-vector: what I am thinking is to
    # store the g_vector whenever possible and pass it along when doing sums of
    # elements with the same degree or multiplications of any two elements.
    # This can potentially slow things down and make life harder. We need to
    # redefine _add_ _mul_ and _lmul_ _rmul_ accordingly if we decide that there
    # is no other way to compute g-vectors
    def __init__(self, parent, value, g_vector=None):
        ElementWrapper.__init__(self, parent=parent, value=value)
        self._g_vector = g_vector

    # This function needs to be removed once AdditiveMagmas.Subobjects
    # implements _add_
    def _add_(self, other):
        return self.parent().retract(self.lift() + other.lift())

    # HACK: LaurentPolynomial_mpair does not know how to compute denominators, we need to lift to its fraction field
    def lift_to_field(self):
        return self.parent().lift_to_field(self)

    # this function is quite disgusting but at least it works for any element of
    # the algebra, can we do better?
    def d_vector(self):
        n = self.parent().rk
        one = self.parent().ambient_field()(1)
        factors = self.lift_to_field().factor()
        initial = []
        non_initial = []
        [(initial if x[1] > 0 and len(x[0].monomials()) == 1 else non_initial).append(x[0]**x[1]) for x in factors]
        initial = prod(initial+[one]).numerator()
        non_initial = prod(non_initial+[one]).denominator()
        v1 = vector(non_initial.exponents()[0][:n])
        v2 = vector(initial.exponents()[0][:n])
        return tuple(v1-v2)

    def g_vector(self):
        # at the moment it is not immediately clear to me how to compute this
        # assuming this is a generic element of the cluster algebra it is not
        # going to be homogeneous, expecially if we are not using principal
        # coefficients. I am not sure it can be done if the bottom part of the
        # exchange matrix is not invertible.
        raise NotImplementederror("This should be computed by substitution")

####
# Seeds
####

class ClusterAlgebraSeed(SageObject):

    def __init__(self, B, C, G, parent):
        self._B = copy(B)
        self._C = copy(C)
        self._G = copy(G)
        self._parent = parent
        self._path = []

    def __copy__(self):
        other = type(self).__new__(type(self))
        other._B = copy(self._B)
        other._C = copy(self._C)
        other._G = copy(self._G)
        other._parent = self._parent
        other._path = copy(self._path)
        return other

    def _repr_(self):
        if self._path == []:
            return "The initial seed of %s"%str(self.parent())
        elif self._path.__len__() == 1:
            return "The seed of %s obtained from the initial by mutating in direction %s"%(str(self.parent()),str(self._path[0]))
        else:
            return "The seed of %s obtained from the initial by mutating along the sequence %s"%(str(self.parent()),str(self._path))
    
    @property
    def depth(self):
        return len(self._path)

    def parent(self):
        return self._parent

    def b_matrix(self):
        return copy(self._B)

    def c_matrix(self):
        return copy(self._C)

    def c_vector(self, j):
        return tuple(self._C.column(j))

    def g_matrix(self):
        return copy(self._G)

    def g_vector(self, j):
        return tuple(self._G.column(j))

    def F_polynomial(self, j):
        return self.parent().F_polynomial(self.g_vector(j))

    def cluster_variable(self, j):
        return self.parent().cluster_variable(self.g_vector(j))
    
    def mutate(self, k, inplace=True, mutating_F=True):
        if inplace:
            seed = self
        else:
            seed = copy(self)

        n = seed.parent().rk

        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction ' + str(k) + '.')

        # store mutation path
        if seed._path != [] and seed._path[-1] == k:
            seed._path.pop()
        else:
            seed._path.append(k)

        # find sign of k-th c-vector
        # Will this be used enough that it should be a built-in function?
        if any(x > 0 for x in seed._C.column(k)):
            eps = +1
        else:
            eps = -1

        # store the g-vector to be mutated in case we are mutating F-polynomials also
        old_g_vector = seed.g_vector(k)

        # G-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[j,k] += max(0, -eps*seed._B[j,k])
        J[k,k] = -1
        seed._G = seed._G*J

        # g-vector path list
        g_vector = seed.g_vector(k)
        if g_vector not in seed.parent().g_vectors_so_far():
            seed.parent()._g_vect_set.add(g_vector)
            seed.parent()._path_dict[g_vector] = copy(seed._path)
            # F-polynomials
            if mutating_F:
                seed.parent()._F_poly_dict[g_vector] = seed._mutated_F(k, old_g_vector)

        # C-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[k,j] += max(0, eps*seed._B[k,j])
        J[k,k] = -1
        seed._C = seed._C*J

        # B-matrix
        seed._B.mutate(k)

        # wrap up
        if not inplace:
            return seed

    def _mutated_F(self, k, old_g_vector):
        alg = self.parent()
        pos = alg._U(1)
        neg = alg._U(1)
        for j in xrange(alg.rk):
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
        return alg._U((pos+neg)/alg.F_polynomial(old_g_vector))

    def mutation_sequence(self, sequence, inplace=True, mutating_F=True):
        seq = iter(sequence)

        if inplace:
            seed = self
        else:
            seed = self.mutate(next(seq), inplace=False, mutating_F=mutating_F)

        for k in seq:
            seed.mutate(k, inplace=True, mutating_F=mutating_F)

        if not inplace:
            return seed

    def path_form_initial_seed(self):
        return copy(self._path)

    # TODO: ideally we should allow to mutate in direction "this g-vector" or
    # "this cluster variable" or "sink", "urban renewal" and all the other
    # options provided by Gregg et al. To do so I guess the best option is to
    # have a generic function transforming all these into an index and use it as
    # a decorator. In this way we can also use it in this __contains__ even
    # though one may write weird things like "sink" in A.current_seed and get
    # True as an answer.
    def __contains__(self, element):
        if isinstance(element, ClusterAlgebraElement ):
            cluster = [ self.cluster_variable(i) for i in xrange(self.parent().rk) ]
        else:
            element = tuple(element)
            cluster = map( tuple, self.g_matrix().columns() )
        return element in cluster


####
# Cluster algebras
####

class ClusterAlgebra(Parent):
    # it would be nice to have inject_variables() to allow the user to
    # automagically export the initial cluster variables into the sage shell.
    # Unfortunately to do this we need to inherit form ParentWithGens and, as a
    # drawback, we get several functions that are meaningless/misleading in a
    # cluster algebra like ngens() or gens()
    # If we only care about initial cluster variables we can always do the following:
    # def inject_variables(self, scope=None, verbose=True):
    #     self._ambient.inject_variables(scope=scope,verbose=verbose)
    # for labeled non-initial cluster variables we can also manually add them to the scope.

    Element = ClusterAlgebraElement

    def __init__(self, data, scalars=ZZ):
        # Temporary variables
        # TODO: right now  we use ClusterQuiver to parse input data. It looks
        # like a good idea but we should make sure it is.
        Q = ClusterQuiver(data)
        B0 = Q.b_matrix()
        n = B0.ncols()
        # We use a different m than ClusterSeed and ClusterQuiver: their m is our m-n.
        # Should we merge this behaviour? what is the notation in CA I-IV? 
        # They use the m-n convention, the m convention comes out of Fock-Goncharov or Gekhtman-Shapiro-Vainshtein
        m = B0.nrows()
        I = identity_matrix(n)

        # add methods that are defined only for special cases
        if n == 2:
            self.greedy_element = MethodType(greedy_element, self, self.__class__)
            self.greedy_coeff = MethodType(greedy_coeff, self, self.__class__)
            self.theta_basis_element = MethodType(theta_basis_element, self, self.__class__)

        # TODO: is ZZ the correct ambient here?
        # ambient space for F-polynomials
        self._U = PolynomialRing(ZZ, ['u%s'%i for i in xrange(n)])

        # already computed data
        # TODO: I am unhappy because _g_vect_set is slightly redundant (we could
        # use _path_dict.keys() instead) but it is faster to check membership in
        # sets than in lists and _path_dict.keys() returns a list. Depending on
        # the number of cluster variables this may be relevant or not.
        self._g_vect_set = set([ tuple(v) for v in I.columns() ])
        self._F_poly_dict = dict([ (v, self._U(1)) for v in self._g_vect_set ])
        self._path_dict = dict([ (v, []) for v in self._g_vect_set ])
        
        # setup Parent and ambient
        # TODO: at the moment self.base() is not a subobject of self.ambient()
        # unfortunately if we change `scalars` to `base` in *** it becomes
        # harder to list generators of ambient 
        if m>n: 
            # do we want change this to its fraction field (so that we can
            # invert coefficients)? 
            base = LaurentPolynomialRing(scalars, ['x%s'%i for i in xrange(n,m)])
        else:
            base = scalars
        # TODO: understand why using CommutativeAlgebras() instead of Rings() makes A(1) complain of missing _lmul_
        Parent.__init__(self, base=base, category=Rings(scalars).Subobjects())
        # *** here we might want to replace scalars with base and m with n
        self._ambient = LaurentPolynomialRing(scalars, 'x', m)
        self._ambient_field = self._ambient.fraction_field()
        # TODO: understand if we need this
        #self._populate_coercion_lists_()

        # these are used for computing cluster variables using separation of
        # additions
        self._y = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(n,m)])) for j in xrange(n)])
        self._yhat = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(m)])) for j in xrange(n)])

        self._B0 = copy(B0)
        self._n = n
        self._m = m
        self.reset_current_seed()

        # internal data for exploring the exchange graph
        self.reset_exploring_iterator()

    # enable standard cohercions: everything that is in the base can be coherced
    def _coerce_map_from_(self, other):
        return self.base().has_coerce_map_from(other)

    def _repr_(self):
        return "Cluster Algebra of rank %s"%self.rk

    def _an_element_(self):
        return self.current_seed.cluster_variable(0)

    @property
    def rk(self):
        r"""
        The rank of ``self`` i.e. the number of cluster variables in any seed of
        ``self``.
        """
        return self._n

    @property
    def current_seed(self):
        r"""
        The current seed of ``self``.
        """
        return self._seed

    def reset_current_seed(self):
        r"""
        Reset the current seed to the initial one
        """
        self._seed = self.initial_seed

    @property
    def initial_seed(self):
        r"""
        Return the initial seed
        """
        n = self.rk
        I = identity_matrix(n)
        return ClusterAlgebraSeed(self._B0[:n,:n], I, I, self)


    def g_vectors_so_far(self):
        r"""
        Return the g-vectors of cluster variables encountered so far.
        """
        return self._g_vect_set

    def F_polynomial(self, g_vector):
        g_vector= tuple(g_vector)
        try:
            return self._F_poly_dict[g_vector]
        except:
            # TODO: improve this error message to include the case in which we
            # already know the path
            # If the path is known, should this method perform that sequence of mutations to compute the desired F-polynomial?
            # Yes, perhaps with the a prompt first, something like:
            #comp = raw_input("This F-polynomial has not been computed yet.  It can be found using %s mutations.  Continue? (y or n):"%str(directions.__len__()))
            #if comp == 'y':
            #    ...compute the F-polynomial...
            raise ValueError("This F-polynomial has not been computed yet. Did you explore the tree with compute_F=False ?")

    @cached_method(key=lambda a,b: tuple(b) )
    def cluster_variable(self, g_vector):
        g_vector = tuple(g_vector)
        if not g_vector in self.g_vectors_so_far():
            raise ValueError("This Cluster Variable has not been computed yet.")
        g_mon = prod([self.ambient().gen(i)**g_vector[i] for i in xrange(self.rk)])
        F_std = self.F_polynomial(g_vector).subs(self._yhat)
        # LaurentPolynomial_mpair does not know how to compute denominators, we need to lift to its fraction field
        F_trop = self.ambient_field()(self.F_polynomial(g_vector).subs(self._y)).denominator()
        return self.retract(g_mon*F_std*F_trop)

    def find_cluster_variable(self, g_vector, depth=infinity):
        r"""
        Returns the shortest mutation path to obtain the cluster variable with
        g-vector ``g_vector`` from the initial seed.

        ``depth``: maximum distance from ``self.current_seed`` to reach.
        """
        g_vector = tuple(g_vector)
        mutation_counter = 0
        while g_vector not in self.g_vectors_so_far() and self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth
            except:
                raise ValueError("Could not find a cluster variable with g-vector %s after %s mutations."%(str(g_vector),str(mutation_counter)))

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

    def ambient_field(self):
        return self._ambient_field

    def lift_to_field(self, x):
        return self.ambient_field()(1)*x.value

    def lift(self, x):
        return x.value

    def retract(self, x):
        return self(x)
    
    def seeds(self, depth=infinity, mutating_F=True, from_current_seed=False):
        r"""
        Return an iterator producing all seeds of ``self`` up to distance
        ``depth`` from ``self.initial_seed`` or ``self.current_seed``.

        If ``mutating_F`` is set to false it does not compute F_polynomials
        """
        if from_current_seed:
            seed = self.current_seed
        else:
            seed = self.initial_seed

        yield seed
        depth_counter = 0
        n = self.rk
        cl = frozenset(seed.g_matrix().columns())
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
                    new_cl = frozenset(new_sd.g_matrix().columns())
                    if new_cl in clusters:
                        if i in clusters[new_cl][1]:
                            clusters[new_cl][1].remove(i)
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

    def reset_exploring_iterator(self):
        self._sd_iter = self.seeds()
        self._explored_depth = 0 

    def explore_to_depth(self, depth):
        while self._explored_depth <= depth:
            try:
                seed = next(self._sd_iter)
                self._explored_depth = seed.depth

    # DESIDERATA. Some of these are probably unrealistic
    def upper_cluster_algebra(self):
        pass

    def upper_bound(self):
        pass

    def lower_bound(self):
        pass




####
# Methods that are only defined for special cases
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
            output += self.greedy_coeff(d_vector,p,q)*x1**(b*p)*x2**(c*q)
    return self.retract(x1**(-a1)*x2**(-a2)*output)

# Is this function something we want to make public or do we want to make this a
# private method changing it to _greedy_coeff ?
# Since we are giving long names to things we might want to change this into
# greedy_coefficient
def greedy_coeff(self,d_vector,p,q):
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
        sum1 += (-1)**(k-1)*self.greedy_coeff(d_vector,p-k,q)*bin
    sum2 = 0
    for l in range(1,q+1):
        bin = 0
        if a1-b*p+l-1 >= l:
            bin = binomial(a1-b*p+l-1,l)
        sum2 += (-1)**(l-1)*self.greedy_coeff(d_vector,p,q-l)*bin
    #print "sum1=",sum1,"sum2=",sum2
    return max(sum1,sum2)


# At the moment I know only how to compute theta basis in rank 2
# maybe we should let ClusterAlgebra have this method for any rank and have a
# NotImplementedError to encourage someone (read: Greg) to code this
#I think Greg already has some code to do this
def theta_basis_element(self, g_vector):
    pass

####
# Scratchpad
####

# Shall we use properties with setters and getters? This is the example
# maybe it is not a great idea but it saves on parenthesis and makes quantities immutable at the same time
# I am not sure I can give an opinion yet, maybe after I see it in action
#class C(object):
#def __init__(self):
#    self._x = None

#@property
#def x(self):
#    """I'm the 'x' property."""
#    return self._x

#@x.setter
#def x(self, value):
#    self._x = value

#@x.deleter
#def x(self):
#    del self._x


