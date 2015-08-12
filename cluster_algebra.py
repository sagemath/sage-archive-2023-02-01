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

# Define decorator to automatically transform vectors and lists to tuples
# is there not a standard function wrapper we could use instead of this?
# TODO: this messes up the function signature. In particular ? and ?? do not
# give the desired results. It is a known issue that is fized in python 3.4!!!!
# do we ever hope to have sage running on python 3?
from functools import wraps
def make_hashable(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        hashable_args = []
        for x in args:
            try:
                hashable_args.append(tuple(x))
            except:
                hashable_args.append(x)

        hashable_kwargs = {}
        for x in kwargs:
            try:
                hashable_kwargs[x] = tuple(kwargs[x])
            except:
                hashable_kwargs[x] = kwargs[x]

        return func(*hashable_args, **hashable_kwargs)
    return wrapper

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

    def parent(self):
        return self._parent

    def b_matrix(self):
        return copy(self._B)

    def F_polynomial(self, j):
        g_vector = tuple(self._G.column(j))
        return self.parent().F_polynomial(g_vector)

    def cluster_variable(self, j):
        g_vector = tuple(self._G.column(j))
        return self.parent().cluster_variable(g_vector)

    def g_vector(self, j):
        return tuple(self._G.column(j))

    def g_matrix(self):
        return copy(self._G)

    def c_vector(self, j):
        return tuple(self._C.column(j))

    def c_matrix(self):
        return copy(self._C)

    def mutate(self, k, inplace=True, mutating_F=True):
        if inplace:
            seed = self
        else:
            seed = copy(self)

        n = self.parent().rk

        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction ' + str(k) + '.')

        # find sign of k-th c-vector
        if any(x > 0 for x in seed._C.column(k)):
            eps = +1
        else:
            eps = -1

        # store the g-vector to be mutated in case we are mutating also F-polynomials
        old_g_vector = tuple(seed._G.column(k))

        # mutate G-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[j,k] += max(0, -eps*seed._B[j,k])
        J[k,k] = -1
        seed._G = seed._G*J

        # F-polynomials
        if mutating_F:
            g_vector = tuple(seed._G.column(k))
            if g_vector not in self.parent().g_vectors_so_far():
                self.parent()._data_dict.setdefault(g_vector, (self._mutated_F(k, old_g_vector),self._path+[k]))

        # C-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[k,j] += max(0, eps*seed._B[k,j])
        J[k,k] = -1
        seed._C = seed._C*J

        # B-matrix
        seed._B.mutate(k)

        # store mutation path
        if seed._path != [] and seed._path[-1] == k:
            seed._path.pop()
        else:
            seed._path.append(k)

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
        return (pos+neg)//alg.F_polynomial(old_g_vector)

    def mutation_sequence(self, sequence, inplace=True, mutating_F=True):
        seq = iter(sequence)

        if inplace:
            seed = self
        else:
            seed = self.mutate(seq.next(), inplace=False, mutating_F=mutating_F)

        for k in seq:
            seed.mutate(k, inplace=True, mutating_F=mutating_F)

        if not inplace:
            return seed

    def path_form_initial_seed(self):
        return copy(self._path)


####
# Cluster algebras
####

class ClusterAlgebra(Parent):
    # it would be nice to have inject_variables() to allow the user to
    # automagically export the initial cluster variables into the sage shell.
    # Unfortunately to do this we need to inherit form ParentWithGens and, as a
    # drawback, we get several functions that are meaningless/misleading in a
    # cluster algebra like ngens() or gens()

    Element = ClusterAlgebraElement

    def __init__(self, B0, scalars=ZZ):
        # Temporary variables
        n = B0.ncols()
        m = B0.nrows()
        I = identity_matrix(n)

        # add methods that are defined only for special cases
        if n == 2:
            self.greedy_element = MethodType(greedy_element, self, self.__class__)
            self.theta_basis_element = MethodType(theta_basis_element, self, self.__class__)

        # TODO: maybe here scalars can be replaced with just ZZ
        # ambient space for F-polynomials
        self._U = PolynomialRing(scalars,['u%s'%i for i in xrange(n)])

        # dictionary of already computed data:
        # index is the g-vector, first entry is the F-polynomial, second entry is path from initial seed
        self._data_dict = dict([ (tuple(v), (self._U(1),[])) for v in I.columns() ])

        self._seed = ClusterAlgebraSeed(B0[:n,:n], I, I, self)

        base = LaurentPolynomialRing(scalars, 'x', n)
        # TODO: understand why using CommutativeAlgebras() instead of Rings() makes A(1) complain of missing _lmul_
        Parent.__init__(self, base=base, category=Rings(scalars).Subobjects())

        self._ambient = LaurentPolynomialRing(scalars, 'x', m)
        self._ambient_field = self._ambient.fraction_field()

        self._y = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(n,m)])) for j in xrange(n)])
        self._yhat = dict([ (self._U.gen(j), prod([self._ambient.gen(i)**B0[i,j] for i in xrange(m)])) for j in xrange(n)])

        self._B0 = copy(B0)
        self._n = n
        self._m = m

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

    def g_vectors_so_far(self):
        return self._data_dict.keys()

    @make_hashable
    def F_polynomial(self, g_vector):
        try:
            return self._data_dict[g_vector][0]
        except:
            # TODO: improve this error message
            raise ValueError("This F-polynomial has not been computed yet. Did you explore the tree with compute_F=False ?")

    @make_hashable
    @cached_method
    def cluster_variable(self, g_vector):
        if not g_vector in self.g_vectors_so_far():
            raise ValueError("This Cluster Variable has not been computed yet.")
        g_mon = prod([self.ambient().gen(i)**g_vector[i] for i in xrange(self.rk)])
        F_std = self.F_polynomial(g_vector).subs(self._yhat)
        # LaurentPolynomial_mpair does not know how to compute denominators, we need to lift to its fraction field
        F_trop = self.ambient_field()(self.F_polynomial(g_vector).subs(self._y)).denominator()
        return self.retract(g_mon*F_std*F_trop)

    @make_hashable
    def find_cluster_variable(self, g_vector, depth=infinity):
        r"""
        Returns the shortest mutation path to obtain the cluster variable with
        given g-vector from the initial seed
        """
        seeds = self.seeds(depth=depth)
        mutation_counter = 0
        while g_vector not in self.g_vectors_so_far():
            try:
                seeds.next()
            except:
                raise ValueError("Could not find a cluster variable with g-vector %s after %s mutations."%(str(g_vector),str(mutation_counter)))
            mutation_counter += 1
        return copy(self._data_dict[g_vector][1])

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

    def seeds(self, depth=infinity):
        yield self.current_seed
        depth_counter = 0
        n = self.rk
        cl = Set(self.current_seed.g_matrix().columns())
        clusters = {}
        clusters[ cl ] = [ self.current_seed, range(n) ]
        gets_bigger = True
        while gets_bigger and depth_counter < depth:
            gets_bigger = False
            keys = clusters.keys()
            for key in keys:
                sd, directions = clusters[key]
                while directions:
                    i = directions.pop()
                    new_sd  = sd.mutate( i, inplace=False )
                    new_cl = Set(new_sd.g_matrix().columns())
                    if new_cl in clusters:
                        if i in clusters[new_cl][1]:
                            clusters[new_cl][1].remove(i)
                    else:
                        gets_bigger = True
                        # doublecheck this way of producing directions for the new seed: it is taken almost verbatim fom ClusterSeed
                        new_directions = [ j for j in xrange(n) if j > i or new_sd.b_matrix()[j,i] != 0 ]
                        clusters[ new_cl ] = [ new_sd, new_directions ]
                        yield new_sd
            depth_counter += 1

    # DESIDERATA. Some of these are probably irrealistic
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
def greedy_element(self, d_vector):
    pass

# At the moment I know only how to compute theta basis in rank 2
def theta_basis_element(self, g_vector):
    pass

####
# Scratchpad
####

# Shall we use properties with setters and getters? This is the example
# maybe it is not a great idea but it saves on parenthesis and makes quantities immutable at the same time
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

