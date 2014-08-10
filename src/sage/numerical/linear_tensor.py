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

    sage: p = MixedIntegerLinearProgram('ppl')   # base ring is QQ
    sage: x = p.new_variable(nonnegative=False)
    sage: 3 + x[0] + 2*x[1]            # a linear function
    3 + x_0 + 2*x_1
    sage: x[0] * vector([3,4]) + 1     # vector linear function
    (1, 1) + (3, 4)*x_0
    sage: x[0] * matrix([[3,1],[4,0]]) + 1   # matrix linear function
    [1 + 3*x_0 x_0]
    [4*x_0     1  ]

.. NOTE::

    For breverity, we just use ``LinearTensor`` in class names. It is
    understood that this refers to the above tensor product
    construction.
"""

from sage.structure.parent import Parent
from sage.structure.element import ModuleElement, Element
from sage.misc.cachefunc import cached_function
from sage.matrix.matrix_space import is_MatrixSpace
from sage.modules.free_module import is_FreeModule
from sage.numerical.linear_functions import is_LinearFunction, LinearFunctionsParent_class

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
# Elements of linear functions tensored with a free module
#
#*****************************************************************************

class LinearTensor(ModuleElement):
    r"""
    A linear function tensored with a free module

    .. warning::

        You should never instantiate :class:`LinearTensor`
        manually. Use the element constructor in the parent
        instead.

    EXAMPLES::

        sage: parent = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
        sage: parent({0: [1,2], 3: [-7,-8]})
        (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
    """

    def __init__(self, parent, f):
        r"""
        Constructor taking a dictionary as its argument.

        INPUT:
        
        - ``parent`` -- the parent :class:`LinearTensorParent_class`.
        
        - ``f`` -- A linear function tensored by a free module is
          represented as a dictionary. The values are the coefficient
          (free module elements) of the variable represented by the
          keys. The key ``-1`` corresponds to the constant term.

        EXAMPLES:

        With a dictionary::

            sage: LT = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]})
            (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
        
            sage: TestSuite(LT).run(skip=['_test_additive_associativity', '_test_elements', '_test_pickling', '_test_zero'])
        """
        ModuleElement.__init__(self, parent)
        assert isinstance(f, dict)
        self._f = f

    def __getitem__(self, *indices):
        """
        Return the linear function component with given tensor indices.

        INPUT:

        - ``*indices`` -- one or more integers. The basis indices of
          the free module. E.g. a single integer for vectors, two for
          matrices.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: lt = p({0:[1,2], 3:[4,5]});  lt
            (1.0, 2.0)*x_0 + (4.0, 5.0)*x_3
            sage: lt[0]
            x_0 + 4*x_3
            sage: lt[1]
            2*x_0 + 5*x_3
        """
        f = dict([key, value.__getitem__(*indices)] for key, value in self._f.items())
        LF = self.parent().linear_functions()
        return LF(f)

    def dict(self):
        r"""
        Return the dictionary corresponding to the tensor product.

        OUTPUT:

        The linear function tensor product is represented as a
        dictionary. The value are the coefficient (free module
        elements) of the variable represented by the keys (which are
        integers). The key ``-1`` corresponds to the constant term.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent().tensor(RDF^2)
            sage: lt = p({0:[1,2], 3:[4,5]})
            sage: lt.dict()
            {0: (1.0, 2.0), 3: (4.0, 5.0)}
        """
        return dict(self._f)

    def _repr_(self):
        """
        Return a string representation.

        OUTPUT:

        String.
        
        EXAMPLES::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: R.<s,t> = RDF[]
            sage: LT = LinearFunctionsParent(RDF).tensor(R)
            sage: LT.an_element()  # indirect doctest
            (s) + (5.0*s)*x_2 + (7.0*s)*x_5

            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT.an_element()  # indirect doctest
            (1.0, 0.0) + (5.0, 0.0)*x_2 + (7.0, 0.0)*x_5
        """
        if is_MatrixSpace(self.parent().free_module()):
            return self._repr_matrix()
        terms = []
        for key in sorted(self._f.keys()):
            coeff = self._f[key]
            if coeff._is_atomic():
                if key == -1:
                    term = '({1})'.format(key, coeff)
                else:
                    term = '({1})*x_{0}'.format(key, coeff)
            else:
                if key == -1:
                    term = '{1}'.format(key, coeff)
                else:
                    term = '{1}*x_{0}'.format(key, coeff)
            terms.append(term)
        return ' + '.join(terms)
            
    def _repr_matrix(self):
        """
        Return a matrix-like string representation.

        OUTPUT:

        String.

        EXAMPLES::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^(2,2))
            sage: LT.an_element()  # indirect doctest
            [1 + 5*x_2 + 7*x_5 0]
            [0                 0]
        """
        MS = self.parent().free_module()
        assert is_MatrixSpace(MS)
        col_lengths = []
        columns = []
        for c in range(MS.ncols()):
            column = []
            for r in range(MS.nrows()):
                cell = repr(self[r, c])
                column.append(cell)
            columns.append(column)
            col_lengths.append(max(map(len, column)))
        s = ''
        for r in range(MS.nrows()):
            if r > 0:
                s += '\n'
            s += '['
            for c in range(MS.ncols()):
                if c > 0:
                    s += ' '
                s += columns[c][r].ljust(col_lengths[c])
            s += ']'
        return s

    def _add_(self, b):
        r"""
        Return sum.

        INPUT:
        
        - ``b`` -- a :class:`LinearTensor`.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]}) + LT({2: [5,6], 3: [2,-2]}) + 16
            (16.0, 16.0) + (1.0, 2.0)*x_0 + (5.0, 6.0)*x_2 + (-5.0, -10.0)*x_3
        """
        result = dict(self._f)
        for key, coeff in b.dict().iteritems():
            result[key] = self._f.get(key, 0) + coeff
        return self.parent()(result)

    def _neg_(self):
        r"""
        Return the negative.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: -LT({0: [1,2], 3: [-7,-8]})
            (-1.0, -2.0)*x_0 + (7.0, 8.0)*x_3
        """
        result = dict([key, -coeff] for key, coeff in self._f.items())
        return self.parent()(result)

    def _sub_(self, b):
        r"""
        Return difference.

        INPUT:
        
        - ``b`` -- a :class:`LinearTensor`.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: LT({0: [1,2], 3: [-7,-8]}) - LT({1: [1,2]})
            (1.0, 2.0)*x_0 + (-1.0, -2.0)*x_1 + (-7.0, -8.0)*x_3
            sage: LT({0: [1,2], 3: [-7,-8]}) - 16
            (-16.0, -16.0) + (1.0, 2.0)*x_0 + (-7.0, -8.0)*x_3
        """
        result = dict(self._f)
        for key, coeff in b.dict().iteritems():
            result[key] = self._f.get(key, 0) - coeff
        return self.parent()(result)

    def _rmul_(self, b):
        r"""
        Return right multiplication by scalar.

        INPUT:

        - ``b`` -- base ring element. The scalar to multiply by.

        OUTPUT:

        A :class:`LinearTensor`.

        EXAMPLE::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LT = LinearFunctionsParent(RDF).tensor(RDF^2)
            sage: 10 * LT({0: [1,2], 3: [-7,-8]})
            (10.0, 20.0)*x_0 + (-70.0, -80.0)*x_3
        """
        result = dict([key, b*coeff] for key, coeff in self._f.items())
        return self.parent()(result)


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
        
    def _element_constructor_(self, x):
        """
        Construt a :class:`LinearTensor` from ``x``.

        INPUT:

        - ``x`` -- anything that defines a :class:`LinearTensor`. See
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
            <class 'sage.numerical.linear_tensor.LinearTensorParent_class_with_category.element_class'>

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
        """
        if is_LinearTensor(x):
            if x.parent() is self:
                return x
            else:
                x = x.dict()
        M = self.free_module()
        if isinstance(x, dict):
            x = dict( (int(key), M(value)) for key, value in x.items() )
        else:
            try:
                x = {-1: M(x)}
            except (TypeError, ValueError):
                x_R = M.base_ring()(x)
                if is_MatrixSpace(M):
                    # Turn constants into diagonal matrices
                    x_matrix = M.zero_matrix()
                    for i in range(min(M.ncols(), M.nrows())):
                        x_matrix[i, i] = x_R
                    x = {-1: x_matrix}
                elif is_FreeModule(M):
                    # Turn constants into vectors with all entries equal
                    x = {-1: M([x_R] * M.degree())}
                else:
                    raise
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
        """
        if R in [int, float, long, complex]:
            return True
        from sage.rings.real_mpfr import mpfr_prec_min
        from sage.rings.all import RealField
        if RealField(mpfr_prec_min()).has_coerce_map_from(R):
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

