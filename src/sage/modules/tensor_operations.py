"""
Helper Classes to implement Tensor Operations

.. warning::

   This module is not meant to be used directly. It just provides
   functionality for other classes to implement tensor operations.

The :class:`VectorCollection` constructs the basis of tensor products
(and symmetric/exterior powers) in terms of a chosen collection of
vectors that generate the vector space(s).

EXAMPLES::

    sage: from sage.modules.tensor_operations import VectorCollection, TensorOperation
    sage: V = VectorCollection([(1,0), (-1, 0), (1,2)], QQ, 2)
    sage: W = VectorCollection([(1,1), (1,-1), (-1, 1)], QQ, 2)
    sage: VW = TensorOperation([V, W], operation='product')

Here is the tensor product of two vectors:

    sage: V.vectors()[0]    
    (1, 0)
    sage: W.vectors()[1]
    (1, -1)
    sage: VW.vectors()[VW.index_map(0, 1)]
    (1, -1, 0, 0)

    sage: rows = []
    sage: for i, j in CartesianProduct(range(3), range(3)): 
    ....:     v = V.vectors()[i]
    ....:     w = W.vectors()[j]
    ....:     i_tensor_j = VW.index_map(i, j)
    ....:     vw = VW.vectors()[i_tensor_j]
    ....:     rows.append([i, v, j, w, i_tensor_j, vw]) 
    sage: table(rows)
          0   (1, 0)    0   (1, 1)    0   (1, 1, 0, 0)
          0   (1, 0)    1   (1, -1)   1   (1, -1, 0, 0)
          0   (1, 0)    2   (-1, 1)   2   (-1, 1, 0, 0)
          1   (-1, 0)   0   (1, 1)    3   (-1, -1, 0, 0)
          1   (-1, 0)   1   (1, -1)   2   (-1, 1, 0, 0)
          1   (-1, 0)   2   (-1, 1)   1   (1, -1, 0, 0)
          2   (1, 2)    0   (1, 1)    4   (1, 1, 2, 2)
          2   (1, 2)    1   (1, -1)   5   (1, -1, 2, -2)
          2   (1, 2)    2   (-1, 1)   6   (-1, 1, -2, 2)
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
from sage.modules.free_module import FreeModule_ambient_field, VectorSpace
from sage.misc.all import cached_method, prod
from sage.matrix.constructor import vector, matrix
from sage.rings.all import ZZ



def symmetrized_coordinate_sums(dim, n):
    """
    Return formal symmetrized sum of multi-indices

    INPUT:

    - ``dim`` -- integer. The dimension (range of each index).

    - ``n`` -- integer. The total number of indices.

    OUTPUT:

    A symmetrized formal sum of multi-indices (tuples of integers)

    EXAMPLES::
    
        sage: from sage.modules.tensor_operations import symmetrized_coordinate_sums
        sage: symmetrized_coordinate_sums(2, 2)
        ((0, 1) + (1, 0), (0, 0), (1, 1))
    """
    from sage.structure.formal_sum import FormalSum
    coordinates = [range(dim) for i in range(n)]
    table = dict()
    from sage.combinat.cartesian_product import CartesianProduct
    for i in CartesianProduct(*coordinates):
        sort_i = tuple(sorted(i))
        x = table.get(sort_i, [])
        x.append([+1, tuple(i)])
        table[sort_i] = x
    return tuple(FormalSum(x) for x in table.values())


def antisymmetrized_coordinate_sums(dim, n):
    """
    Return formal anti-symmetrized sum of multi-indices

    INPUT:

    - ``dim`` -- integer. The dimension (range of each index).

    - ``n`` -- integer. The total number of indices.

    OUTPUT:

    An anti-symmetrized formal sum of multi-indices (tuples of integers)

    EXAMPLES::
    
        sage: from sage.modules.tensor_operations import antisymmetrized_coordinate_sums
        sage: antisymmetrized_coordinate_sums(3, 2)
        ((0, 1) - (1, 0), (0, 2) - (2, 0), (1, 2) - (2, 1))
    """
    from sage.structure.formal_sum import FormalSum
    table = []
    from sage.groups.perm_gps.permgroup_named import SymmetricGroup
    S_d = SymmetricGroup(n)
    from sage.combinat.combination import Combinations
    for i in Combinations(range(dim), n):
        i = tuple(i)
        x = []
        for g in S_d:
            x.append([g.sign(), g(i)])
        x = FormalSum(x)
        table.append(x)
    return tuple(table)


class VectorCollection(FreeModule_ambient_field):
    """
    A collection of rays in a vector space.

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
    
        sage: from sage.modules.tensor_operations import VectorCollection
        sage: R = VectorCollection([(1,0), (0,1), (1,2)], QQ, 2);  R
        Vector space of dimension 2 over Rational Field
    
    TESTS::

        sage: R.vectors()
        ((1, 0), (0, 1), (1, 2))
        sage: r = R._vectors[0]
        sage: type(r)
        <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>
        sage: r.parent() is R
        True
        sage: r.is_immutable()
        True
    """
    def __init__(self, vector_collection, base_ring, dim):
        super(VectorCollection, self).__init__(base_ring, dim)
        self._n_vectors = len(vector_collection)
        self._vectors = tuple(self(r) for r in vector_collection)
        for r in self._vectors:
            r.set_immutable()
        if matrix(base_ring, self._vectors).rank() != self.degree():
            raise ValueError('the vectors must span the ambient vector space')
        self._all_indices = tuple(map(ZZ, range(0, self._n_vectors)))

    def vectors(self):
        """
        Return the collection of vectors

        OUTPUT:

        A tuple of vectors. The vectors that were specified in the
        constructor, in the same order.
        
        EXAMPLES::

            sage: from sage.modules.tensor_operations import VectorCollection
            sage: V = VectorCollection([(1,0), (0,1), (1,2)], QQ, 2)
            sage: V.vectors()
            ((1, 0), (0, 1), (1, 2))
        """
        return self._vectors
    
    def n_vectors(self):
        """
        Return the number of vectors

        OUTPUT:

        Integer.


            sage: from sage.modules.tensor_operations import VectorCollection
            sage: V = VectorCollection([(1,0), (0,1), (1,2)], QQ, 2)
            sage: V.n_vectors()
            3
        """
        return len(self._vectors)


class TensorOperation(VectorCollection):
    """
    Auxiliary class to compute the tensor product of two
    :class:`VectorCollection` objects.

    .. warning::
    
        This class is only used as a base class for filtered vector
        spaces. You should not use it yourself.

    INPUT:

    - ``vector_collections`` -- a nonempty list/tuple/iterable of
      :class:`VectorCollection` objects.

    - ``operation`` -- string. The tensor operation. Currently allowed
      values are ``product``, ``symmetric``, and ``antisymmetric``.
    
    .. todo::

        More general tensor operations (specified by Young tableaux)
        should be implemented.

    EXAMPLES::

        sage: from sage.modules.tensor_operations import \
        ....:      VectorCollection, TensorOperation
        sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
        sage: S = VectorCollection([(1,), (-1,)], QQ, 1)
        sage: R_tensor_S = TensorOperation([R, S])
        sage: R_tensor_S.index_map(0, 0)
        0
        sage: matrix(ZZ, 3, 2, lambda i,j: R_tensor_S.index_map(i, j))
        [0 1]
        [2 3]
        [3 2]
        sage: R_tensor_S.vectors()
        ((1, 0), (-1, 0), (1, 2), (-1, -2))
    """
    def __init__(self, vector_collections, operation='product'):
        assert all(isinstance(V, VectorCollection) for V in vector_collections)
        self._base_ring = base_ring = vector_collections[0].base_ring()
        assert all(V.base_ring() is base_ring for V in vector_collections)
        self._V = tuple(vector_collections)
        self._vectors = []
        self._index_map = dict()
        if operation=='product':
            self._init_product()
        elif operation=='symmetric':
            assert all(V is self._V[0] for V in self._V)
            self._init_symmetric()
        elif operation=='antisymmetric':
            assert all(V is self._V[0] for V in self._V)
            self._init_antisymmetric()
        else:
            raise ValueError('invalid operation')
        vectors = self._vectors
        dim = 0 if len(vectors) == 0 else len(vectors[0])
        del self._vectors
        del self._base_ring
        super(TensorOperation, self).__init__(vectors, base_ring, dim)

    def _init_product_vectors(self, i):
        """
        Helper to build up ``self._vectors`` incrementally during the constructor.

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: S = VectorCollection([(1,), (-1,)], QQ, 1)
            sage: R_tensor_S = TensorOperation([R,S])
            sage: R_tensor_S.index_map(1, 1)
            3
            sage: R_tensor_S.index_map(2, 0)
            3
            sage: R_tensor_S.vectors()   # indirect doctest
            ((1, 0), (-1, 0), (1, 2), (-1, -2))
        """
        rays = [ self._V[j].vectors()[k] for j, k in enumerate(i) ]
        v = []
        from sage.combinat.cartesian_product import CartesianProduct
        for r in CartesianProduct(*rays):
            v.append(prod(r))
        v = tuple(v)
        try:
            result = self._vectors.index(v)
        except ValueError:
            self._vectors.append(v)
            result = len(self._vectors)-1
        return result

    def _init_power_operation_vectors(self, i, linear_combinations):
        """
        Helper to build up ``self._vectors`` incrementally during the constructor.

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym2_R = TensorOperation([R,R], operation='symmetric')
            sage: Sym2_R.vectors()    # indirect doctest
            ((0, 1, 0), (2, 1, 0), (-2, -1, 0), (4, 1, 4), (-4, -1, -4))
            sage: Alt2_R = TensorOperation([R, R], operation='antisymmetric')
            sage: Alt2_R.vectors()    # indirect doctest
            ((2), (-2))
        """
        rays = [ self._V[j].vectors()[k] for j, k in enumerate(i) ]
        v = []
        for coordinate_linear_combination in linear_combinations:
            v_entry = self._base_ring.zero()
            for coeff, index in coordinate_linear_combination:
                v_entry += coeff * prod(rays[j][k] for  j,k in enumerate(index))
            v.append(v_entry)
        v = tuple(v)
        if all(vi == 0 for vi in v):
            return None
        try:
            result = self._vectors.index(v)
        except ValueError:
            self._vectors.append(v)
            result = len(self._vectors)-1
        return result

    def _init_product(self):
        """
        Initialization for the tensor product
        
        EXAMPLES::
        
            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: S = VectorCollection([(1,), (-1,)], QQ, 1)
            sage: R_tensor_S = TensorOperation([R,S], operation='product')
            sage: sorted(R_tensor_S._index_map.iteritems())
            [((0, 0), 0), ((0, 1), 1), ((1, 0), 2), ((1, 1), 3), ((2, 0), 3), ((2, 1), 2)]
        """
        V_list_indices = [range(V.n_vectors()) for V in self._V]
        from sage.combinat.cartesian_product import CartesianProduct
        for i in CartesianProduct(*V_list_indices):
            self._index_map[tuple(i)] = self._init_product_vectors(i)
        self._symmetrize_indices = False

    def _init_symmetric(self):
        """
        Initialization for the symmetric product.

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym2_R = TensorOperation([R,R], operation='symmetric')
            sage: sorted(Sym2_R._index_map.iteritems())
            [((0, 0), 0), ((0, 1), 1), ((0, 2), 2), ((1, 1), 3), ((1, 2), 4), ((2, 2), 3)]
        """
        V_list_indices = [range(V.n_vectors()) for V in self._V]
        Sym = symmetrized_coordinate_sums(self._V[0].dimension(), len(self._V))
        from sage.combinat.cartesian_product import CartesianProduct
        for i in CartesianProduct(*V_list_indices):
            if any(i[j-1]>i[j] for j in range(1,len(i))):
                continue
            self._index_map[tuple(i)] = self._init_power_operation_vectors(i, Sym)
        self._symmetrize_indices = True

    def _init_antisymmetric(self):
        """
        Initialization for the antisymmetric product

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Alt2_R = TensorOperation([R, R], operation='antisymmetric')
            sage: sorted(Alt2_R._index_map.iteritems())
            [((0, 1), 0), ((0, 2), 1)]
        """
        n = len(self._V)
        dim = self._V[0].degree()
        Alt = antisymmetrized_coordinate_sums(dim, n)
        from sage.combinat.combination import Combinations
        for i in Combinations(range(self._V[0].n_vectors()), n):
            ray = self._init_power_operation_vectors(i, Alt)
            if ray is not None:
                self._index_map[tuple(i)] = ray
        self._symmetrize_indices = True

    def index_map(self, *i):
        """
        Return the result of the tensor operation.

        INPUT: 
        
        - ``*i`` -- list of integers. The indices (in the
        corresponding factor of the tensor operation) of the domain
        vector.

        OUTPUT::

        The index (in :meth:`vectors`) of the image of the tensor
        product/operation acting on the domain vectors indexed by `i`.

        ``None`` is returned if the tensor operation maps the
        generators to zero (usually because of antisymmetry).

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (1,2), (-1,-2)], QQ, 2)
            sage: Sym3_R = TensorOperation([R]*3, 'symmetric')
            sage: Sym3_R.index_map(0, 1, 2)
            4
            sage: Sym3_R.index_map(2, 0, 1)
            4
            sage: Sym3_R.index_map(2, 1, 0)
            4

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (0,1), (-2,-3)], QQ, 2)
            sage: detR = TensorOperation([R]*2, 'antisymmetric') 
            sage: detR.index_map(1, 0)
            0
            sage: detR.index_map(0, 1)
            0

        TESTS::

            sage: sorted(detR._index_map.iteritems())
            [((0, 1), 0), ((0, 2), 1), ((1, 2), 2)]
            sage: detR.vectors()
            ((1), (-3), (2))
        """
        if len(i) == 1 and isinstance(i[0], (list, tuple)):
            i = tuple(i[0])
        if self._symmetrize_indices:
            i = tuple(sorted(i))
        try:
            return self._index_map[i]
        except KeyError:
            return None

    def preimage(self):
        """
        A choice of pre-image multi-indices.

        OUTPUT:
        
        A list of multi-indices (tuples of integers) whose image is
        the entire image under the :meth:`index_map`.

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (0,1), (-2,-3)], QQ, 2)
            sage: detR = TensorOperation([R]*2, 'antisymmetric') 
            sage: sorted(detR.preimage())
            [(0, 1), (0, 2), (1, 2)]
            sage: sorted(detR.codomain())
            [0, 1, 2]
        """
        return self._index_map.keys()

    def codomain(self):
        """
        The codomain of the index map.

        OUTPUT:

        A list of integers. The image of :meth:`index_map`.

        EXAMPLES::

            sage: from sage.modules.tensor_operations import \
            ....:      VectorCollection, TensorOperation
            sage: R = VectorCollection([(1,0), (0,1), (-2,-3)], QQ, 2)
            sage: detR = TensorOperation([R]*2, 'antisymmetric') 
            sage: sorted(detR.preimage())
            [(0, 1), (0, 2), (1, 2)]
            sage: sorted(detR.codomain())
            [0, 1, 2]
        """
        return self._index_map.values()






