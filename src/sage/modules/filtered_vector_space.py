"""
`\ZZ`-filtered vector spaces

This module implements filtered vector spaces, that is, a descending
sequence of vector spaces

.. math::

    \cdots \supset F_d \supset F_{d+1} \supset F_{d+2} \supset \cdots 

To construct a filtered vector space, use the
:func:`FilteredVectorSpace` command. It supports easy creation of
simple filtrations, for example the trivial one::

    sage: FilteredVectorSpace(2, base_ring=RDF)
    RDF^2

The next-simplest filtration has a single non-trivial inclusion
between `V_d` and `V_{d+1}`::
  
    sage: d = 1
    sage: V = FilteredVectorSpace(2, d);  V
    QQ^2 >= 0
    sage: [V.get_degree(i).dimension() for i in range(0,4)]
    [2, 2, 0, 0]

To construct general filtrations, you need to pass generating vectors
for the nested subspaces. For example::

    sage: r1 = (1, 0, 0)
    sage: r2 = (0, 1, 0)
    sage: r3 = (1, 2, 0)
    sage: V = FilteredVectorSpace({0:[r1, r2], 1:[r3], 3:[]});  V
    QQ^3 >= QQ^2 >= QQ^1 >= QQ^1 >= 0

For degrees `d` that are not specified, the associated vector subspace
is the same as the next-lower degree, that is, `V_d \simeq
V_{d-1}`. In the above example, this means that

* `V_d \simeq \QQ^3` for `d<0`
* `V_0 = \mathop{span}(r_1, r_2) \simeq \QQ^2` 
* `V_1 = V_2 = \mathop{span}(r_3) \simeq \QQ` 
* `V_d = 0` for `d \geq 3` 

That is::

    sage: QQ3 = V.ambient_vector_space();  QQ3
    Vector space of dimension 3 over Rational Field
    sage: V.get_degree(-1) == QQ3
    True
    sage: V.get_degree(0) == QQ3.span([r1, r2])
    True
    sage: V.get_degree(1) == V.get_degree(2) == QQ3.span([r3])
    True
    sage: V.get_degree(3) == V.get_degree(4) == QQ3.span([])
    True

Alternatively, you can just pass the rays once and then refer to them
by index::

    sage: FilteredVectorSpace([r1, r2, r3], {0:[0,1], 1:[2], 3:[]})
    QQ^3 >= QQ^2 >= QQ^1 >= QQ^1 >= 0


"""

#*****************************************************************************
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of 
#  the License, or (at your option) any later version.  
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.all import QQ, ZZ, RDF, RR, Integer
from sage.categories.fields import Fields
from sage.modules.free_module import FreeModule_ambient_field, VectorSpace
from sage.matrix.constructor import vector, matrix, block_matrix, zero_matrix, identity_matrix
from sage.misc.all import uniq, cached_method, prod


def is_FilteredVectorSpace(X):
    """
    Test whether ``X`` is a filtered vector space.
    
    This function is for library use only.

    INPUT:
    
    - ``X`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::
    
        sage: from sage.modules.filtered_vector_space import is_FilteredVectorSpace
        sage: V = FilteredVectorSpace(2,1)
        sage: is_FilteredVectorSpace(V)
        True
        sage: is_FilteredVectorSpace('ceci n\'est pas une pipe')
        False
    """
    return isinstance(X,FilteredVectorSpace_class)


def FilteredVectorSpace(arg1, arg2=None, base_ring=QQ, check=True):
    """
    Construct a filtered vector space.

    INPUT:

    This function accepcts various input that determines the vector space and filtration.

    - Just the dimensionFilteredVectorSpace(dimension): Return the trivial filtration
      (where all vector spaces are isomorphic).

    - FilteredVectorSpace(rays, filtration, base_ring=QQ)`

    - Dimension and cutoff, see
      :func:`FilteredVectorSpace_from_dim_cutoff` for
      arguments. Construct a filtration with only one non-trivial step
      `V\supset 0` at the given cutoff degree.

    In addition, the following keyword arguments are supported:

    - ``base_ring`` -- a field (optional, default `\QQ`). The base
      field of the vector space. Must be a field.

    EXAMPLES::

        sage: FilteredVectorSpace(2,1)
        QQ^2 >= 0
        sage: FilteredVectorSpace(2)
        QQ^2
        sage: FilteredVectorSpace([(1,0),(0,1)], [(1,[1]), (3,)])
        QQ^2 >= QQ^1 >= QQ^1 >= 0
    """
    if base_ring not in Fields():
        raise ValueError('the base_ring argument must be a field')
    try: 
        return FilteredVectorSpace_from_dim_cutoff(arg1, arg2, base_ring, check)
    except (ValueError, TypeError):
        pass
    return FilteredVectorSpace_from_rays_filtration(arg1, arg2, base_ring, check)
    

def FilteredVectorSpace_from_dim_cutoff(dim, cutoff=None, base_ring=QQ, check=True):
    """
    Construct a filtered vector space.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import FilteredVectorSpace_from_dim_cutoff
        sage: V = FilteredVectorSpace_from_dim_cutoff(2,5);  V
        QQ^2 >= 0
        sage: V.get_degree_indices(5)
        (0, 1)
        sage: V.get_degree_indices(6)
        ()
    """
    if dim not in ZZ:
        raise ValueError('dimension must be an integer')
    dim = ZZ(dim)
    rays = identity_matrix(base_ring, dim).columns()
    filtration = dict()
    if cutoff is not None:
        filtration[ZZ(cutoff+1)] = tuple()
    return FilteredVectorSpace_class(rays, filtration, base_ring, dim, check=check)


def FilteredVectorSpace_from_rays_filtration(rays, filtration, base_ring=QQ, check=True):
    """
    Construct a filtered vector space.
    
    INPUT:

    - ``rays`` -- a list of vectors, or something convertible to
      them. The rays spanning various subspaces. The set of all rays
      must be of full rank.

    - ``filtration`` -- a list or iterable of filtration steps. Each
      filtration step is a pair ``(degree, ray_indices)``. The
      ``ray_indices`` are a list or iterable of ray indices, which
      span a subspace of the vector space. The integer ``degree``
      stipulates that all filtration steps of degree higher or equal
      than ``degree`` (up to the next filtration step) are said
      subspace.

    Instead of the ``rays``, you can just specify the filtration using
    rays instead of ray indices.

    EXAMPLES::
    
        sage: from sage.modules.filtered_vector_space import FilteredVectorSpace_from_rays_filtration
        sage: V = FilteredVectorSpace_from_rays_filtration(
        ....:      [(1,0),(0,1),(-1,-1)], [(1,[1]), (3,)]);  V
        QQ^2 >= QQ^1 >= QQ^1 >= 0
        sage: V.get_degree_indices(0)
        (0, 1, 2)

    Alternative notation: Use rays instead of ray indices::
    
        sage: r = (1,0)
        sage: FilteredVectorSpace({1:[r], 3:[]}).is_isomorphic(V)
        True
    """
    if filtration is None:
        filtration = rays
        rays = None

    # convert filtration list/tuple to dict
    if not isinstance(filtration, dict):
        filtration_dict = dict()
        for filt in filtration:
            d = ZZ(filt[0])
            filtration_dict[d] = [] if len(filt) == 1 else list(filt[1])
        filtration = filtration_dict

    # convert ray notation to ray indices
    if rays is None:
        if len(filtration) == 0:
            raise ValueError('you need to specify at least one ray to deduce the dimension')
        rays = []
        for r in filtration.values():
            rays += list(r)
        rays = tuple(uniq(rays))
    else:
        rays = tuple(rays)
    for d in filtration.keys():
        rays_d = filtration[d]
        try:
            rays_d = [ZZ(rays.index(r)) for r in rays_d]
        except ValueError:
            rays_d = map(ZZ, list(rays_d))
        filtration[d] = tuple(rays_d)

    # verify filtration data
    for d, rays_d in filtration.iteritems():
        if any(r not in ZZ for r in rays_d):
            raise ValueError('ray index must be integer')
        if any(r < 0 or r >= len(rays) for r in rays_d):
            raise ValueError('ray index out of bounds')

    # deduce dimension
    if len(rays)==0:
        dim = ZZ(0)
    else:
        dim = ZZ(len(rays[0]))

    # complete rays to a generating set
    ray_matrix = matrix(base_ring, rays)
    if ray_matrix.rank() < dim:
        complement = VectorSpace(base_ring, dim).span(rays).complement()
        rays = rays + complement.gens()
        
    return FilteredVectorSpace_class(rays, filtration, base_ring, dim, check=check)


##############################################################################
class RayCollection(FreeModule_ambient_field):
    """
    A collection of rays in an ambient vector space.

    .. warning::
    
        This class is only used as a base class for filtered vector
        spaces. You should not use it yourself.

    INPUT:

    - ``dim`` -- integer. The dimension of the ambient vector space.

    - ``base_ring`` -- a field. The base field of the ambient vector space.

    - ``rays`` -- any list/iterable of things than can be converted
      into vectors of the ambient vector space. These will be used to
      span the subspaces of the filtration. Must span the ambient
      vector space.

    EXAMPLES::
    
        sage: from sage.modules.filtered_vector_space import RayCollection
        sage: R = RayCollection([(1,0), (0,1), (1,2)], QQ, 2);  R
        Vector space of dimension 2 over Rational Field
    
    TESTS::

        sage: R._rays
        ((1, 0), (0, 1), (1, 2))
        sage: r = R._rays[0]
        sage: type(r)
        <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        sage: r.parent() is R
        True
        sage: r.is_immutable()
        True
    """

    def __init__(self, rays, base_ring, dim):
        super(RayCollection, self).__init__(base_ring, dim)
        self._n_rays = len(rays)
        self._rays = tuple(self(r) for r in rays)
        for r in self._rays:
            r.set_immutable()
        if matrix(base_ring, self._rays).rank() != self.degree():
            raise ValueError('the rays must span the ambient vector space')
        self._all_indices = tuple(map(ZZ, range(0, self._n_rays)))


##############################################################################
class RayCollection_tensor_operation(SageObject):
    """
    Auxiliary class to compute the tensor product of two
    :class:`RayCollection` objects.

    .. warning::
    
        This class is only used as a base class for filtered vector
        spaces. You should not use it yourself.

    INPUT:

    - ``P``, ``Q`` -- two :class:`RayCollection` objects.

    EXAMPLES::

        sage: from sage.modules.filtered_vector_space import \
        ...      RayCollection, RayCollection_tensor_operation
        sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
        sage: S = RayCollection([(1,), (-1,)], QQ, 1)
        sage: R_tensor_S = RayCollection_tensor_operation([R,S])
        sage: R_tensor_S(0,0)
        0
        sage: matrix(ZZ, 3,2, lambda i,j:R_tensor_S(i,j))
        [0 1]
        [2 3]
        [3 2]
        sage: R_tensor_S._rays
        [(1, 0), (-1, 0), (1, 2), (-1, -2)]
    """
    def __init__(self, V_list, operation='product'):
        assert all(isinstance(V, RayCollection) for V in V_list)
        self._V_list = V_list
        self._base_ring = self._V_list[0].base_ring()
        self._rays = []
        self._result_dict = dict()
        if operation=='product':
            self._product()
            self._symmetrize_result = False
        elif operation=='symmetric':
            assert all(V is V_list[0] for V in V_list)
            self._symmetric()
            self._symmetrize_result = True
        elif operation=='antisymmetric':
            assert all(V is V_list[0] for V in V_list)
            self._antisymmetric()
            self._symmetrize_result = True
        else:
            raise ValueError('I don\'t understand operation='+str(operation)+'.')

    @cached_method
    def _symmetrized_coordinate_sums(self):
        """
        EXAMPLES::
        
            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: R_tensor_S = RayCollection_tensor_operation([R]*2, operation='symmetric')
            sage: R_tensor_S._symmetrized_coordinate_sums()
            ((0, 1) + (1, 0), (0, 0), (1, 1))
        """
        from sage.structure.formal_sum import FormalSum
        d = self._V_list[0].dimension()
        coordinates = [ range(0,d) for V in self._V_list ]
        table = dict()
        from sage.combinat.cartesian_product import CartesianProduct
        for i in CartesianProduct(*coordinates):
            sort_i = tuple(sorted(i))
            x = table.get(sort_i,[])
            x.append([+1,tuple(i)])
            table[sort_i] = x
        table = tuple( FormalSum(x) for x in table.values() )
        return table

    @cached_method
    def _antisymmetrized_coordinate_sums(self):
        """
        EXAMPLES::
        
            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0,0), (1,2,0), (-1,-2,1)], QQ, 3)
            sage: R_tensor_S = RayCollection_tensor_operation([R]*2, operation='symmetric')
            sage: R_tensor_S._antisymmetrized_coordinate_sums()
            ((0, 1) - (1, 0), (0, 2) - (2, 0), (1, 2) - (2, 1))
        """
        from sage.structure.formal_sum import FormalSum
        dim = self._V_list[0].dimension()
        n = len(self._V_list)
        table = []
        from sage.groups.perm_gps.permgroup_named import SymmetricGroup
        S_d = SymmetricGroup(n)
        from sage.combinat.combination import Combinations
        for i in Combinations(range(0,dim), n):
            i = tuple(i)
            x = []
            for g in S_d:
                x.append([g.sign(), g(i)])
            x = FormalSum(x)
            table.append(x)
        return tuple(table)

    def _product_ray(self, i):
        """
        EXAMPLES::

            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: S = RayCollection([(1,), (-1,)], QQ, 1)
            sage: R_tensor_S = RayCollection_tensor_operation([R,S])
            sage: R_tensor_S(1,1)
            3
            sage: R_tensor_S(2,0)
            3
            sage: R_tensor_S._rays[3]
            (-1, -2)
        """
        rays = [ self._V_list[j]._rays[k] for j,k in enumerate(i) ]
        v = []
        from sage.combinat.cartesian_product import CartesianProduct
        for r in CartesianProduct(*rays):
            v.append(prod(r))
        v = vector(self._base_ring, v)
        try:
            result = self._rays.index(v)
        except ValueError:
            self._rays.append(v)
            result = len(self._rays)-1
        return result

    def _power_operation_ray(self, i, linear_combinations):
        """
        EXAMPLES::

            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym2_R = RayCollection_tensor_operation([R,R], operation='symmetric')
            sage: Sym2_R._rays    # indirect doctest
            [(0, 1, 0), (2, 1, 0), (-2, -1, 0), (4, 1, 4), (-4, -1, -4)]
            sage: Alt2_R = RayCollection_tensor_operation([R,R], operation='antisymmetric')
            sage: Alt2_R._rays    # indirect doctest
            [(2), (-2)]
        """
        rays = [ self._V_list[j]._rays[k] for j,k in enumerate(i) ]
        v = []
        for coordinate_linear_combination in linear_combinations:
            v_entry = self._base_ring(0)
            for coeff, index in coordinate_linear_combination:
                v_entry += coeff * prod(rays[j][k] for  j,k in enumerate(index))
            v.append(v_entry)
        v = vector(self._base_ring, v)
        if v.is_zero():
            return None
        try:
            result = self._rays.index(v)
        except ValueError:
            self._rays.append(v)
            result = len(self._rays)-1
        return result

    def _product(self):
        """
        EXAMPLES::
        
            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: S = RayCollection([(1,), (-1,)], QQ, 1)
            sage: R_tensor_S = RayCollection_tensor_operation([R,S], operation='product')
            sage: R_tensor_S._result_dict
            {(0, 1): 1, (0, 0): 0, (2, 1): 2, (2, 0): 3, (1, 0): 2, (1, 1): 3}
        """
        V_list_indices = [ range(0,V._n_rays) for V in self._V_list ]
        from sage.combinat.cartesian_product import CartesianProduct
        for i in CartesianProduct(*V_list_indices):
            self._result_dict[tuple(i)] = self._product_ray(i)

    def _symmetric(self):
        """
        EXAMPLES::

            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym2_R = RayCollection_tensor_operation([R,R], operation='symmetric')
            sage: Sym2_R._result_dict
            {(0, 1): 1, (1, 2): 4, (0, 0): 0, (1, 1): 3, (2, 2): 3, (0, 2): 2}
        """
        V_list_indices = [ range(0,V._n_rays) for V in self._V_list ]
        Sym = self._symmetrized_coordinate_sums()
        from sage.combinat.cartesian_product import CartesianProduct
        for i in CartesianProduct(*V_list_indices):
            if any(i[j-1]>i[j] for j in range(1,len(i))):
                continue
            self._result_dict[tuple(i)] = self._power_operation_ray(i, Sym)

    def _antisymmetric(self):
        """
        EXAMPLES::

            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Alt2_R = RayCollection_tensor_operation([R,R], operation='antisymmetric')
            sage: Alt2_R._result_dict
            {(0, 1): 0, (0, 2): 1}
        """
        k = len(self._V_list)
        n = self._V_list[0]._n_rays
        Alt = self._antisymmetrized_coordinate_sums()
        from sage.combinat.combination import Combinations
        for i in Combinations(range(0, n), k):
            ray = self._power_operation_ray(i, Alt)
            if ray is not None:
                self._result_dict[tuple(i)] = ray

    def __call__(self, *i):
        """
        Return the result of the tensor operation on the rays given by
        the multiindex ``i``.

        EXAMPLES::

            sage: from sage.modules.filtered_vector_space import \
            ...      RayCollection, RayCollection_tensor_operation
            sage: R = RayCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym3_R = RayCollection_tensor_operation([R]*3, 'symmetric')
            sage: Sym3_R.__call__(0,1,2)
            4
            sage: Sym3_R.__call__(2,0,1)
            4
            sage: Sym3_R.__call__(2,1,0)
            4
        """
        if len(i)==1 and isinstance(i[0],(list,tuple)):
            i = tuple(i[0])
        if self._symmetrize_result:
            i = tuple(sorted(i))
        return self._result_dict[i]


##############################################################################
class FilteredVectorSpace_class(RayCollection):

    def __init__(self, rays, filtration, base_ring, dim, check=True):
        r"""
        A descending filtration of a vector space
    
        INPUT:

        - ``dim`` -- integer. The dimension of the ambient vector space.

        - ``base_ring`` -- a field. The base field of the ambient vector space.

        - ``rays`` -- any list of things than can be converted into
          vectors of the ambient vector space. These will be used to
          span the subspaces of the filtration.

        - ``filtration`` -- a dictionary of filtration steps in ray
          index notation. See
          :func:`FilteredVectorSpace_from_rays_filtration` for
          details.

        TESTS::

            sage: from sage.modules.filtered_vector_space import FilteredVectorSpace_class
            sage: rays = [(1,0,0), (1,1,0), (1,2,0), (-1,-1, 0), (0,0,1)]
            sage: FilteredVectorSpace_class(rays, {2:(4,)}, QQ, 3)
            QQ^3 >= QQ^1
            sage: FilteredVectorSpace_class(rays, {2:(4,), 3:()}, QQ, 3)
            QQ^3 >= QQ^1 >= 0

        The trivial filtration::
            
            sage: FilteredVectorSpace_class(rays, {}, QQ, 3)
            QQ^3

        The empty vector space::

            sage: FilteredVectorSpace_class([], {}, QQ, 0)
            0
        """
        assert isinstance(dim, Integer)
        assert base_ring in Fields()
        super(FilteredVectorSpace_class, self).__init__(rays, base_ring, dim)

        assert isinstance(filtration, dict)
        for degree, generators in filtration.iteritems():
            assert isinstance(degree, Integer)
            assert isinstance(generators, tuple)
            assert all(isinstance(r, Integer) for r in generators)
        filtration = sorted(list(filtration.iteritems()))  # sort by degree
        self._filt = tuple(filtration)
        self._trivial = (len(self._filt) == 0)

        if not check: return
        assert matrix(self._rays).rank() == self.dimension()
        # check that it is an increasing filtration
        #prev_gens_set = set(self._all_indices)
        #for deg, gens in self._filt:
        #    gens_set = set(gens)
        #    assert gens_set.issubset(prev_gens_set)
        #    prev_gens_set = gens_set

    def is_trivial(self):
        """
        Return whether the filtration is trivial.

        OUTPUT:
        
        Boolean.
        
        EXAMPLES::

            sage: V = FilteredVectorSpace(1,0);  V
            QQ^1 >= 0
            sage: V.is_trivial()
            False

            sage: V = FilteredVectorSpace([[1]], []);  V
            QQ^1
            sage: V.is_trivial()
            True
        """
        return self._trivial

    @cached_method
    def min_degree(self):
        r"""
        Return the lowest degree of the filtration.

        OUTPUT:

        Integer. The largest degree of the (descending) filtration
        such that the filtered vector space is still the whole ambient
        space.
        
        0 is returned for the trivial filtration.
        
        EXAMPLES::

            sage: FilteredVectorSpace(1,3).min_degree()
            3
            sage: FilteredVectorSpace([[1]], []).min_degree()
            0
        """
        if self._trivial: 
            return 0
        return self._filt[0][0]-1

    @cached_method
    def max_degree(self):
        r"""
        Return the highest degree of the filtration.

        OUTPUT:

        Integer. The smallest degree of the filtration such that the
        filtration is constant to the right.

        0 is returned for the trivial filtration.
        
        EXAMPLES::

            sage: FilteredVectorSpace(1,3).max_degree()
            4
            sage: FilteredVectorSpace([[1]], []).max_degree()
            0
        """
        if self._trivial: 
            return 0
        return self._filt[-1][0]

    @cached_method
    def is_finite(self):
        """
        Return whether the filtration eventually becomes the trivial vector space.

        EXAMPLES::

            sage: FilteredVectorSpace(1,4).is_finite()
            True
            sage: FilteredVectorSpace(1).is_finite()
            False
        """
        gens_oo = self.get_degree_indices(self.max_degree())
        return len(gens_oo) == 0

    def changes_in_degree(self, d):
        r"""
        Whether the filtration step at degree `d` is non-trivial.
        
        That is, the quotient `F_{d-1}/F_d \not= 0`. In other words,
        the when reading the inclusion sequence from the left the
        filtration just switched to a smaller subspace at position
        `d`.

        EXAMPLES:

        The following filtration is `F_d=\QQ` for `d\leq 4` and
        `F_d=0` for d\geq 5`::

            sage: F = FilteredVectorSpace(1,4);  F
            QQ^1 >= 0
            sage: F.changes_in_degree(4)
            False
            sage: F.changes_in_degree(5)
            True
            sage: F.changes_in_degree(6)
            False
        """
        for F_i in self._filt:
            if F_i[0]==d:
                return True
        return False

    def get_degree_indices(self, d):
        r"""
        Return the indices of rays generating the degree-``d`` entry
        of the filtration.

        INPUT:

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        A tuple of ray indices spanning the degree-``d`` vector space
        in the filtration.

        EXAMPLES::
        
            sage: rays = [(1,0),(1,1),(1,2),(-1,-1)]
            sage: F = FilteredVectorSpace(rays, [(1,[1]), (3,)])
            sage: for d in range(F.min_degree(), F.max_degree()+1):
            ...       print d, F.get_degree_indices(d)
             0 (0, 1, 2, 3)
             1 (1,)
             2 (1,)
             3 ()
        """
        if self._trivial:
            return self._all_indices
        generators = self._all_indices
        for deg, gens in self._filt:
            if deg>d:
                return generators
            generators = gens
        return generators

    def get_degree_gens(self, d):
        r""" 
        Return rays generating the degree-``d`` entry of the
        filtration.

        INPUT:

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        A tuple of rays spanning the degree-``d`` vector space in the
        filtration.

        EXAMPLES::
        
            sage: rays = [(1,0,0),(1,1,0),(1,2,0),(-1,-1,-1)]
            sage: F = FilteredVectorSpace(rays, [[3,[1]], [1,[1,2]]])
            sage: for d in range(F.min_degree(), F.max_degree()+1):
            ...       print 'd=', d, '  indices=', F.get_degree_indices(d)
            d= 0   indices= (0, 1, 2, 3)
            d= 1   indices= (1, 2)
            d= 2   indices= (1, 2)
            d= 3   indices= (1,)
        """
        return tuple( self._rays[i] for i in self.get_degree_indices(d) )

    def get_degree(self, d):
        r"""
        Return the degree-``d`` entry of the filtration.

        INPUT:

        - ``d`` -- Integer. The desired degree of the filtration.

        OUTPUT:

        The degree-``d`` vector space in the filtration as subspace of
        the ambient space.

        EXAMPLES::
        
            sage: rays = [(1,0),(1,1),(1,2),(-1,-1)]
            sage: F = FilteredVectorSpace(rays, [[3,[1]], [1,[1,2]]])
            sage: F.get_degree(2)
            Vector space of degree 2 and dimension 2 over Rational Field
            Basis matrix:
            [1 0]
            [0 1]
        """
        try:
            cache = self._get_degree
        except AttributeError:
            cache = dict()
            self._get_degree = cache
        d_equiv = self.equivalent_degree(d)
        try:
            return cache[d_equiv]
        except KeyError:
            result = self.span( self.get_degree_gens(d) )
            cache[d_equiv] = result
            return result

    def _repr_(self):
        r"""
        Return as string representation of ``self``.

        OUTPUT: 
        
        A string.

        EXAMPLES::
        
            sage: rays = [(1,0),(1,1),(1,2),(-1,-1)]
            sage: FilteredVectorSpace(rays, [[1,[1]], [3,[]]])._repr_()
            'QQ^2 >= QQ^1 >= QQ^1 >= 0'
            sage: FilteredVectorSpace(rays, [[1,[1]]])._repr_()
            'QQ^2 >= QQ^1'
        """
        if self.dimension()==0:
            return '0'
        field_name = str(self.base_ring())
        if self.base_ring() == QQ:
            field_name = 'QQ'
        elif self.base_ring() == RDF:
            field_name = 'RDF'
        elif self.base_ring() == RR:
            field_name = 'RR'
        dims = [ field_name + '^' + str(self.get_degree(i).dimension())
                 for i in range(self.min_degree(), self.max_degree()) ] 
        final = self.get_degree(self.max_degree()).dimension()
        if final == 0:
            dims.append('0')
        else:
            dims.append(field_name + '^' + str(final))
        return (' >= ').join(dims)
        
    def __cmp__(self, other):
        """
        Compare two filtered vector spaces.

        .. warning::
        
            This method tests whether the spanning rays are the
            same. Use :meth:`is_isomorphic` to test for mathematical
            equivalence.

        EXAMPLES::

            sage: V = FilteredVectorSpace(2,0)
            sage: W = FilteredVectorSpace([(1,0),(0,1)], [[1]])
            sage: V == W
            True
            sage: V is W
            False

            sage: W = FilteredVectorSpace([(1,0),(1,1)], [[1]])
            sage: V == W
            False

        TESTS::

            sage: P = toric_varieties.P2()
            sage: T_P = P.sheaves.tangent_bundle()
            sage: O_P = P.sheaves.trivial_bundle(1)
            sage: S1 = T_P + O_P
            sage: S2 = O_P + T_P
            sage: S1._filt[0].is_isomorphic(S2._filt[0])  # known bug
            True
            sage: FilteredVectorSpace(2, base_ring=QQ) == FilteredVectorSpace(2, base_ring=GF(5))
            False
        """
        c = super(FilteredVectorSpace_class, self).__cmp__(other)
        if c!=0: return c
        c = cmp(type(self), type(other))
        if c!=0: return c
        for r_self, r_other in zip(self._rays, other._rays):
            # avoid implicit comparison of parents
            c = cmp(tuple(r_self), tuple(r_other))
            if c!=0: return c
        return cmp(self._filt, other._filt)

    def chomp(self):
        """
        Return the same filtration with trivial filtration steps
        removed.

        EXAMPLES::
          
            sage: F = FilteredVectorSpace([(1,0),(0,1),(-1,-1)], [(1,[1,2]), (3,[1]), (4,[]), (6,)]);  F
            QQ^2 >= QQ^2 >= QQ^2 >= QQ^1 >= QQ^0 >= QQ^0 >= 0
            sage: F.chomp()
            QQ^2 >= QQ^1 >= 0
        """
        filt = []
        dim = self.dimension()
        for deg, gens in self._filt:
            new_dim = self.get_degree(deg).dimension()
            if new_dim == dim:
                continue
            dim = new_dim
            filt.append([deg, gens])
        return FilteredVectorSpace(self._rays, filt, self.base_ring(), self.dimension())

    def is_isomorphic(self, other):
        """
        Test whether ``self`` and ``other`` are mathematically equivalent.

        INPUT:
        
        - ``other`` -- anything.

        OUTPUT:
        
        Boolean.

        EXAMPLES::
        
            sage: V = FilteredVectorSpace(2,0)
            sage: W = FilteredVectorSpace([(1,-1),(2,1)], [[1]])
            sage: V == W
            False
            sage: V.is_isomorphic(W)
            True
        """
        if not is_FilteredVectorSpace(other):
            return False
        if self.ambient_vector_space() != other.ambient_vector_space():
            return False
        min_deg = self.min_degree()
        if not min_deg==other.min_degree():
            return False
        max_deg = self.max_degree()
        if not max_deg==other.max_degree():
            return False
        # two vector space are isomorphic iff their dimension is equal
        return all(self.get_degree(i).dimension() == 
                   other.get_degree(i).dimension()
                   for i in range(min_deg, max_deg+1))
        
    def ambient_vector_space(self):
        """
        Return the ambient vector space.
        
        EXAMPLES::
        
            sage: M = QQ^3
            sage: M.ambient_vector_space()
            Vector space of dimension 3 over Rational Field
        """
        return VectorSpace(self.base_ring(), self.dimension())
    
    ambient_module = ambient_vector_space

    def direct_sum(self, other):
        """
        Return the direct sum of ``self`` and ``other``.
        
        EXAMPLES::
        
            sage: V = FilteredVectorSpace(2, 0)
            sage: W = FilteredVectorSpace([(1,-1),(2,1)], [[1]])
            sage: V.direct_sum(W)
            QQ^4 >= 0
            sage: V+W
            QQ^4 >= 0
            sage: _.is_isomorphic(FilteredVectorSpace(4,0))
            True

            sage: W = FilteredVectorSpace([(1,-1),(2,1)], [[2,[1]],[3]])
            sage: V+W
            QQ^4 >= QQ^2 >= QQ^1 >= 0
        """
        def shift(lst):
            return tuple(x+self._n_rays for x in lst)
        rays = \
            [ list(r) + [0]*other.dimension() for r in self._rays  ] + \
            [ [0]*self.dimension() + list(r)  for r in other._rays ]
        min_deg = min(self.min_degree(), other.min_degree())
        max_deg = max(self.max_degree(), other.max_degree())
        filtration = dict()
        for deg in range(min_deg+1, max_deg+1):
            if self.changes_in_degree(deg) or other.changes_in_degree(deg):
                gens = self.get_degree_indices(deg) + shift(other.get_degree_indices(deg))
                filtration[ZZ(deg)] = gens
        dim = self.dimension() + other.dimension()
        return self.__class__(rays, filtration, self.base_ring(), dim)

    __add__ = direct_sum
    
    def steps_right_iter(self):
        """
        Iterate over "right edges" of the filtration steps.

        The ray indices returned at degree `d` generate `F_d` and all
        `F_k` with `k>=d` up to the next higher filtration step.

        EXAMPLES::

            sage: W = FilteredVectorSpace([(1,-1),(2,1)], [[0,[1]],[2]])
            sage: for deg in range(W.min_degree(), W.max_degree()+1):
            ...       print 'deg =', deg, '  gens =', W.get_degree_indices(deg)
            deg = -1   gens = (0, 1)
            deg = 0   gens = (1,)
            deg = 1   gens = (1,)
            deg = 2   gens = ()

            sage: for deg, gens in W.steps_right_iter():
            ...       print 'deg =', deg, '  gens =', gens
            deg = -1   gens = (0, 1)
            deg = 0   gens = (1,)
            deg = 2   gens = ()

            sage: for deg, gens in W.steps_left_iter():
            ...       print 'deg =', deg, '  gens =', gens
            deg = -1   gens = (0, 1)
            deg = 1   gens = (1,)
            deg = 2   gens = ()
        """
        try:
            last_full_degree = self._filt[0][0] - 1
        except IndexError:
            last_full_degree = 0
        yield (last_full_degree, self._all_indices)
        for step in self._filt:
            yield step

    def steps_left_iter(self):
        """
        Iterate over "left edges" of the filtration steps.

        The ray indices returned at degree `d` generate `F_d` and all
        `F_k` with `k<=d` up to the next lower filtration step.

        EXAMPLES::

            sage: W = FilteredVectorSpace([(1,-1),(2,1)], [[0,[1]],[2]])
            sage: for deg in range(W.min_degree(), W.max_degree()+1):
            ...       print 'deg =', deg, '  gens =', W.get_degree_indices(deg)
            deg = -1   gens = (0, 1)
            deg = 0   gens = (1,)
            deg = 1   gens = (1,)
            deg = 2   gens = ()

            sage: for deg, gens in W.steps_right_iter():
            ...       print 'deg =', deg, '  gens =', gens
            deg = -1   gens = (0, 1)
            deg = 0   gens = (1,)
            deg = 2   gens = ()

            sage: for deg, gens in W.steps_left_iter():
            ...       print 'deg =', deg, '  gens =', gens
            deg = -1   gens = (0, 1)
            deg = 1   gens = (1,)
            deg = 2   gens = ()
        """
        deg = 0
        gens = self._all_indices
        for step in self._filt:
            deg = step[0]-1
            yield deg, gens
            gens = step[1]
        yield deg+1, gens
        
    @cached_method
    def equivalent_degree(self, d):
        r"""
        Return a choice of degree `d'` with the same filtered
        subspace.

        This method is intended to be used if you are caching
        quantities derived from the filtered vector space `F_d`. As
        long as the quantity is only dependent on the vector subspace,
        there is no need to cache the result multiple times.

        OUTPUT:

        An integer `d'` such that `F_d = F_{d'}`.

        EXAMPLES::

            sage: rays = [(1,0),(0,1)]
            sage: F = FilteredVectorSpace(rays, [(1,[0]), (3,[])])
            sage: for i in range(-1,5):
            ...       i_equiv = F.equivalent_degree(i) 
            ...       print i, i_equiv, F.get_degree(i)==F.get_degree(i_equiv), F.get_degree_gens(i)
            -1 0 True ((1, 0), (0, 1))
            0 0 True ((1, 0), (0, 1))
            1 1 True ((1, 0),)
            2 1 True ((1, 0),)
            3 3 True ()
            4 3 True ()
        """
        deg = self.min_degree()
        for F in self._filt:
            if F[0]>d:
                return deg
            deg = F[0]
        return deg

    def tensor_product(self, other):
        """
        Return the tensor product of ``self`` and ``other``

        EXAMPLES::

            sage: F1 = FilteredVectorSpace(1,1)
            sage: F2 = FilteredVectorSpace(1,2)
            sage: F1.tensor_product(F2)
            QQ^1 >= 0
            sage: F1 * F2
            QQ^1 >= 0
            
            sage: F1.min_degree()
            1
            sage: F2.min_degree()
            2
            sage: (F1*F2).min_degree()
            3
        """
        ray_tensor_index = RayCollection_tensor_operation([self, other], 'product')

        filtration = dict()
        for self_step in self.steps_right_iter():
            for other_step in other.steps_right_iter():
                deg = self_step[0]+other_step[0]
                filt_step = filtration.get(deg,[])
                for i in self_step[1]:
                    for j in other_step[1]:
                        i_tensor_j = ray_tensor_index(i,j)
                        filt_step.append(ZZ(i_tensor_j))
                filtration[ZZ(deg)] = filt_step
        return FilteredVectorSpace(ray_tensor_index._rays, filtration, base_ring=self.base_ring()).chomp()
        
    __mul__ = tensor_product

    def _power_operation(self, n, operation):
        """
        Return tensor power operation.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1,1) + FilteredVectorSpace(1,2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F._power_operation(2, 'symmetric')
            QQ^3 >= QQ^3 >= QQ^2 >= QQ^1 >= QQ^0 >= 0
            sage: F._power_operation(2, 'antisymmetric')
            QQ^1 >= QQ^1 >= QQ^1 >= QQ^0 >= QQ^0 >= 0
        """
        ray_power_index = RayCollection_tensor_operation([self]*n, operation)
        iters = [ list(self.steps_right_iter()) ] * n
        filt_dict = dict()
        from sage.combinat.cartesian_product import CartesianProduct
        for steps in CartesianProduct(*iters):
            deg = sum(s[0] for s in steps)
            filt_step = filt_dict.get(deg, set())
            for i in CartesianProduct(*[s[1] for s in steps]):
                try:
                    pow_i = ray_power_index(i)
                    filt_step.update([pow_i])
                except KeyError:
                    pass
            filt_dict[deg] = filt_step
        filt = [ (deg, list(gens.difference([None]))) 
                 for deg, gens in filt_dict.iteritems() ]
        return FilteredVectorSpace(ray_power_index._rays, filt, base_ring=self.base_ring())

    def exterior_power(self, n):
        """
        Return the `n`-th graded exterior power.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1,1) + FilteredVectorSpace(1,2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F.exterior_power(1)
            QQ^2 >= QQ^1 >= 0
            sage: F.exterior_power(2)
            QQ^1 >= 0
            sage: F.exterior_power(3)
            0
            sage: F.wedge(2)
            QQ^1 >= 0
        """
        return self._power_operation(n, 'antisymmetric').chomp()
    
    wedge = exterior_power

    def symmetric_power(self, n):
        """
        Return the `n`-th graded symmetric power.

        EXAMPLES::

            sage: F = FilteredVectorSpace(1,1) + FilteredVectorSpace(1,2);  F
            QQ^2 >= QQ^1 >= 0
            sage: F.symmetric_power(2)
            QQ^3 >= QQ^2 >= QQ^1 >= 0
        """
        return self._power_operation(n, 'symmetric').chomp()

