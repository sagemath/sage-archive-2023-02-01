"""
Tensor Products of R-Linear Functions and Free R-Modules

In Sage, matrices assume that the base is a ring. Hence, we cannot
construct matrices whose entries are linear functions in Sage. Really,
they should be thought of as the tensor product of the R-module of
linear functions and the R-module of vector/matrix spaces (`R` is
``QQ`` or ``RDF`` for our purposes).

You should not construct any tensor products by calling the parent
directly. This is also why none of the classes are imported in the
global namespace. The come into play whenever you have vector or
matrix MIP linear expressions/constraints. The intented way to
construct them is implinicly by acting with vectors or matrices on
linear functions. For example::

    sage: mip.<x> = MixedIntegerLinearProgram('ppl')   # base ring is QQ
    sage: 3 + x[0] + 2*x[1]            # a linear function
    3 + x_0 + 2*x_1
    sage: x[0] * vector([3,4]) + 1     # vector linear function
    (1, 1) + (3, 4)*x_0
    sage: x[0] * matrix([[3,1],[4,0]]) + 1   # matrix linear function
    [1 + 3*x_0 x_0]
    [4*x_0     1  ]

Just like :mod:`sage.numerical.linear_functions`, (in)equalities
become symbolic inequalities. See
:mod:`~sage.numerical.linear_tensor_constraints` for detais.

.. NOTE::

    For breverity, we just use ``LinearTensor`` in class names. It is
    understood that this refers to the above tensor product
    construction.
"""
from copy import copy

from sage.structure.parent import Parent
from sage.misc.cachefunc import cached_function
from sage.matrix.matrix_space import is_MatrixSpace
from sage.modules.free_module import is_FreeModule
from sage.numerical.linear_functions import is_LinearFunction, LinearFunctionsParent_class
from sage.numerical.linear_tensor_element import LinearTensor

#*****************************************************************************
#
# Utility functions to test that something is a linear function / constraint
#
#*****************************************************************************

def is_LinearTensor(x):
    """
    Test whether ``x`` is a tensor product of linear functions with a
    free module.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: x = p.new_variable(nonnegative=False)
        sage: from sage.numerical.linear_tensor import is_LinearTensor
        sage: is_LinearTensor(x[0] - 2*x[2])
        False
        sage: is_LinearTensor('a string')
        False
    """
    return isinstance(x, LinearTensor)


#*****************************************************************************
#
# Factory functions for the parents to ensure uniqueness
#
#*****************************************************************************

@cached_function
def LinearTensorParent(free_module_parent, linear_functions_parent):
    """
    Return the parent for the tensor product over the common ``base_ring``.

    The output is cached, so only a single parent is ever constructed
    for a given base ring.

    INPUT:
    
    - ``free_module_parent`` -- module. A free module, like vector or
      matrix space.

    - ``linear_functions_parent`` -- linear functions. The linear
      functions parent.

    OUTPUT:

    The parent of the tensor product of a free module and linear
    functions over a common base ring.

    EXAMPLES::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent
        sage: from sage.numerical.linear_tensor import LinearTensorParent
        sage: LinearTensorParent(QQ^3, LinearFunctionsParent(QQ))
        Tensor product of Vector space of dimension 3 over Rational Field and Linear functions over Rational Field

        sage: LinearTensorParent(ZZ^3, LinearFunctionsParent(QQ))
        Traceback (most recent call last):
        ...
        ValueError: base rings must match
    """
    if free_module_parent.base_ring() != linear_functions_parent.base_ring():
        raise ValueError('base rings must match')
    if not isinstance(linear_functions_parent, LinearFunctionsParent_class):
        raise TypeError('linear_functions_parent must be a parent of linear functions')
    return LinearTensorParent_class(free_module_parent, linear_functions_parent)



#*****************************************************************************
#
# Parent of linear functions tensored with a free module
#
#*****************************************************************************

class LinearTensorParent_class(Parent):
    r"""
    The parent for all linear functions over a fixed base ring.

    .. warning::

        You should use :func:`LinearTensorParent` to construct
        instances of this class.

    INPUT/OUTPUT:

    See :func:`LinearTensorParent`

    EXAMPLES::

        sage: from sage.numerical.linear_tensor import LinearTensorParent_class
        sage: LinearTensorParent_class
        <class 'sage.numerical.linear_tensor.LinearTensorParent_class'>
    """
    Element = LinearTensor

    def __init__(self, free_module, linear_functions):
        """
        The Python constructor

        INPUT/OUTPUT:

        See :func:`LinearTensorParent`

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LinearFunctionsParent(RDF).tensor(RDF^2)
            Tensor product of Vector space of dimension 2 over Real Double
            Field and Linear functions over Real Double Field
        """
        self._free_module = free_module
        self._linear_functions = linear_functions
        base_ring = linear_functions.base_ring()
        from sage.categories.modules_with_basis import ModulesWithBasis
        Parent.__init__(self, base=base_ring, category=ModulesWithBasis(base_ring))

    def free_module(self):
        return self._free_module
        
    def linear_functions(self):
        return self._linear_functions

    def _repr_(self):
        """
        Return a string representation
        
        OUTPUT:

        String.

        EXAMPLES::

            sage: MixedIntegerLinearProgram().linear_functions_parent()
            Linear functions over Real Double Field
        """
        return 'Tensor product of {0} and {1}'.format(self.free_module(), self.linear_functions())

    def _convert_constant(self, m):
        """
        Convert ``m`` to a constant free module element.

        OUTPUT:

        A :meth:`free_module`` element.

        EXAMPLES::

            sage: mip = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LF.tensor(RDF^2)._convert_constant(42)
            (42.0, 42.0)
            sage: LF.tensor(RDF^(2,2))._convert_constant(42)
            [42.0  0.0]
            [ 0.0 42.0]
        """
        M = self.free_module()
        m = M.base_ring()(m)
        if is_MatrixSpace(M):
            # Turn constants into diagonal matrices
            m_matrix = copy(M.zero_matrix())
            for i in range(min(M.ncols(), M.nrows())):
                m_matrix[i, i] = m
            m_matrix.set_immutable()
            return m_matrix
        elif is_FreeModule(M):
            # Turn constants into vectors with all entries equal
            m_vector = M([m] * M.degree())
            return m_vector
        else:
            return M(m)
        
    def _element_constructor_(self, x):
        """
        Construt a :class:`LinearTensor` from ``x``.

        INPUT:

        - ``x`` -- anything that defines a
          :class:`~sage.numerical.linear_tensor_element.LinearTensor`. See
          examples.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: LT = p.linear_functions_parent().tensor(RDF^2)
            sage: LT._element_constructor_(123)
            (123.0, 123.0)

        Construct from dict with type conversion to RDF vector::

            sage: LT({1:[1, 2]})    # indirect doctest
            (1.0, 2.0)*x_1
            sage: type(_)
            <type 'sage.numerical.linear_tensor_element.LinearTensor'>

        Construct from scalar:

            sage: LT(123)    # indirect doctest
            (123.0, 123.0)
        
        Similar, over ``QQ`` and with matrices instead of vectors::

            sage: p_QQ = MixedIntegerLinearProgram(solver='ppl')
            sage: LT_QQ = p_QQ.linear_functions_parent().tensor(QQ^(2, 2))
            sage: LT_QQ({-1:[[1/2, 1/3], [2, 3]], 2:[[3/4, 1/4], [0, 0]]})
            [1/2 + 3/4*x_2 1/3 + 1/4*x_2]
            [2             3            ]
            sage: LT_QQ(42.1)
            [421/10 0     ]
            [0      421/10]

        Construct from a linear function::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LF_ZZ = LinearFunctionsParent(ZZ)
            sage: lf = LF_ZZ({-1:3, 1:2, 3:1})
            sage: LT(lf)
            (3.0, 3.0) + (2.0, 2.0)*x_1 + (1.0, 1.0)*x_3
        """
        M = self.free_module()
        if is_LinearTensor(x):
            if x.parent() is self:
                return x
            else:
                x = x.dict()
        elif is_LinearFunction(x):
            x = dict([key, self._convert_constant(value)] for key, value in x.dict().items())
        elif isinstance(x, dict):
            x = dict([int(key), M(value)] for key, value in x.items())
        else:
            try:
                x = {-1: M(x)}
            except (TypeError, ValueError):
                x_R = M.base_ring()(x)
                x = {-1: self._convert_constant(x_R)}
        return self.element_class(self, x)

    def _coerce_map_from_(self, R):
        """
        Allow coercion of scalars into tensors.

        INPUT:

        - ``R`` -- a ring.

        OUTPUT:

        Boolean. Whether there is a coercion map.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: parent = p.linear_functions_parent()
            sage: parent.coerce(int(2))
            2
            sage: parent._coerce_map_from_(int)
            True

            sage: mip.<x> = MixedIntegerLinearProgram()
            sage: LF = mip.linear_functions_parent()
            sage: LT = LF.tensor(RDF^3)
            sage: LT.has_coerce_map_from(LF)
            True
        """
        if self.free_module().has_coerce_map_from(R):
            return True
        if self.linear_functions().has_coerce_map_from(R):
            return True
        return False

    def _an_element_(self):
        """
        Returns an element

        OUTPUT:

        A linear function tensored with a free module.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: p._an_element_()
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
            sage: p.an_element()   # indirect doctest
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
        """
        m = self.free_module().an_element()
        return self._element_constructor_({-1:m, 2:5*m, 5:7*m})



