"""
Linear Functions and Constraints

This module implements linear functions (see :class:`LinearFunction`)
in formal variables and chained (in)equalities between them (see
:class:`LinearConstraint`). By convention, these are always written as
either equalities or less-or-equal. For example::

    sage: p = MixedIntegerLinearProgram()
    sage: x = p.new_variable()
    sage: f = 1 + x[1] + 2*x[2];  f     #  a linear function
    1 + x_0 + 2*x_1
    sage: type(f)
    <type 'sage.numerical.linear_functions.LinearFunction'>

    sage: c = (0 <= f);  c    # a constraint
    0 <= 1 + x_0 + 2*x_1
    sage: type(c)
    <type 'sage.numerical.linear_functions.LinearConstraint'>

Note that you can use this module without any reference to linear
programming, it only implements linear functions over a base ring and
constraints. However, for ease of demonstration we will always
construct them out of linear programs (see
:mod:`~sage.numerical.mip`).

Constraints can be equations or (non-strict) inequalities. They can be
chained::

    sage: p = MixedIntegerLinearProgram()
    sage: x = p.new_variable()
    sage: x[0] == x[1] == x[2] == x[3]
    x_0 == x_1 == x_2 == x_3

    sage: ieq_01234 = x[0] <= x[1] <= x[2] <= x[3] <= x[4]
    sage: ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4

If necessary, the direction of inequality is flipped to always write
inequalities as less or equal::

    sage: x[5] >= ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4 <= x_5

    sage: (x[5] <= x[6]) >= ieq_01234
    x_0 <= x_1 <= x_2 <= x_3 <= x_4 <= x_5 <= x_6
    sage: (x[5] <= x[6]) <= ieq_01234
    x_5 <= x_6 <= x_0 <= x_1 <= x_2 <= x_3 <= x_4

.. WARNING::

    The implementation of chained inequalities uses a Python hack to
    make it work, so it is not completely robust. In particular, while
    constants are allowed, no two constants can appear next to
    eachother. The following does not work for example::

        sage: x[0] <= 3 <= 4
        True

    If you really need this for some reason, you can explicitly convert
    the constants to a :class:`LinearFunction`::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent
        sage: LF = LinearFunctionsParent(QQ)
        sage: x[0] <= LF(3) <= LF(4)
        x_0 <= 3 <= 4

TESTS:

See :trac:`12091`::

    sage: p = MixedIntegerLinearProgram()
    sage: b = p.new_variable()
    sage: b[0] <= b[1] <= 2
    x_0 <= x_1 <= 2
    sage: list(b[0] <= b[1] <= 2)
    [x_0, x_1, 2]
    sage: 1 >= b[1] >= 2*b[0]
    2*x_0 <= x_1 <= 1
    sage: b[2] >= b[1] >= 2*b[0]
    2*x_0 <= x_1 <= x_2
"""

#*****************************************************************************
#       Copyright (C) 2012 Nathann Cohen <nathann.cohen@gmail.com>
#       Copyright (C) 2012 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2016 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from cpython.object cimport Py_EQ, Py_GE, Py_LE, Py_GT, Py_LT

from sage.misc.fast_methods cimport hash_by_id
from sage.structure.parent cimport Parent
from sage.structure.element cimport ModuleElement, Element
from sage.misc.cachefunc import cached_function

#*****************************************************************************
#
# Utility functions to test that something is a linear function / constraint
#
#*****************************************************************************

cpdef is_LinearFunction(x):
    """
    Test whether ``x`` is a linear function

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: x = p.new_variable()
        sage: from sage.numerical.linear_functions import is_LinearFunction
        sage: is_LinearFunction(x[0] - 2*x[2])
        True
        sage: is_LinearFunction('a string')
        False
    """
    return isinstance(x, LinearFunction)

def is_LinearConstraint(x):
    """
    Test whether ``x`` is a linear constraint

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: x = p.new_variable()
        sage: ieq = (x[0] <= x[1])
        sage: from sage.numerical.linear_functions import is_LinearConstraint
        sage: is_LinearConstraint(ieq)
        True
        sage: is_LinearConstraint('a string')
        False
    """
    return isinstance(x, LinearConstraint)

#*****************************************************************************
#
# Factory functions for the parents to ensure uniqueness
#
#*****************************************************************************

@cached_function
def LinearFunctionsParent(base_ring):
    """
    Return the parent for linear functions over ``base_ring``.

    The output is cached, so only a single parent is ever constructed
    for a given base ring.

    INPUT:

    - ``base_ring`` -- a ring. The coefficient ring for the linear
      functions.

    OUTPUT:

    The parent of the linear functions over ``base_ring``.

    EXAMPLES::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent
        sage: LinearFunctionsParent(QQ)
        Linear functions over Rational Field
    """
    return LinearFunctionsParent_class(base_ring)

@cached_function
def LinearConstraintsParent(linear_functions_parent):
   """
   Return the parent for linear functions over ``base_ring``.

   The output is cached, so only a single parent is ever constructed
   for a given base ring.

    INPUT:

    - ``linear_functions_parent`` -- a
      :class:`LinearFunctionsParent_class`. The type of linear
      functions that the constraints are made out of.

    OUTPUT:

    The parent of the linear constraints with the given linear functions.

    EXAMPLES::

        sage: from sage.numerical.linear_functions import (
        ....:     LinearFunctionsParent, LinearConstraintsParent)
        sage: LF = LinearFunctionsParent(QQ)
        sage: LinearConstraintsParent(LF)
        Linear constraints over Rational Field
   """
   return LinearConstraintsParent_class(linear_functions_parent)


#*****************************************************************************
#
# Elements of linear functions or constraints
#
#*****************************************************************************

cdef chained_comparator_left = None
cdef chained_comparator_right = None
cdef LinearConstraint chained_comparator_replace = None

cdef class LinearFunctionOrConstraint(ModuleElement):
    """
    Base class for :class:`LinearFunction` and :class:`LinearConstraint`.

    This class exists solely to implement chaining of inequalities
    in constraints.
    """
    def __richcmp__(py_left, py_right, int op):
        """
        Create an inequality or equality object, possibly chaining
        several together.

        EXAMPLES::

            sage: p.<x> = MixedIntegerLinearProgram()
            sage: x[0].__le__(x[1])
            x_0 <= x_1

        ::

            sage: p = MixedIntegerLinearProgram()
            sage: from sage.numerical.linear_functions import LinearFunction
            sage: p({2 : 5, 3 : 2}) <= p({2 : 3, 9 : 2})
            5*x_2 + 2*x_3 <= 3*x_2 + 2*x_9

            sage: p({2 : 5, 3 : 2}) >= p({2 : 3, 9 : 2})
            3*x_2 + 2*x_9 <= 5*x_2 + 2*x_3

            sage: p({2 : 5, 3 : 2}) == p({2 : 3, 9 : 2})
            5*x_2 + 2*x_3 == 3*x_2 + 2*x_9

        We can chain multiple (in)equalities::

            sage: p.<b> = MixedIntegerLinearProgram()
            sage: b[0] == 1 == b[1] == 2 == b[2] == 3
            x_0 == 1 == x_1 == 2 == x_2 == 3
            sage: b[0] <= 1 <= b[1] <= 2 <= b[2] <= 3
            x_0 <= 1 <= x_1 <= 2 <= x_2 <= 3
            sage: b[0] <= b[1] <= b[2] <= b[3]
            x_0 <= x_1 <= x_2 <= x_3

        Other comparison operators are not allowed::

            sage: b[0] < b[1]
            Traceback (most recent call last):
            ...
            ValueError: strict < is not allowed, use <= instead
            sage: b[0] > b[1]
            Traceback (most recent call last):
            ...
            ValueError: strict > is not allowed, use >= instead
            sage: b[0] != b[1]
            Traceback (most recent call last):
            ...
            ValueError: inequality != is not allowed, use one of <=, ==, >=

        Mixing operators is also not allowed::

            sage: 1 <= b[1] >= 2
            Traceback (most recent call last):
            ...
            ValueError: incorrectly chained inequality
            sage: 1 >= b[1] <= 2
            Traceback (most recent call last):
            ...
            ValueError: incorrectly chained inequality
            sage: 1 == b[1] <= 2
            Traceback (most recent call last):
            ...
            ValueError: cannot mix equations and inequalities
            sage: 1 >= b[1] == 2
            Traceback (most recent call last):
            ...
            ValueError: cannot mix equations and inequalities

        TESTS::

            sage: p.<x> = MixedIntegerLinearProgram()
            sage: cm = sage.structure.element.get_coercion_model()
            sage: cm.explain(10, p(1), operator.le)
            Coercion on left operand via
                Conversion map:
                  From: Integer Ring
                  To:   Linear functions over Real Double Field
            Arithmetic performed after coercions.
            Result lives in Linear functions over Real Double Field
            Linear functions over Real Double Field

            sage: operator.le(10, x[0])
            10 <= x_0
            sage: x[0] <= 1
            x_0 <= 1
            sage: x[0] >= 1
            1 <= x_0
            sage: 1 <= x[0]
            1 <= x_0
            sage: 1 >= x[0]
            x_0 <= 1

        This works with non-Sage types too, see :trac:`14540`::

            sage: p.<b> = MixedIntegerLinearProgram()
            sage: int(1) <= b[0] <= int(2)
            1 <= x_0 <= 2
            sage: int(1) >= b[0] >= int(2)
            2 <= x_0 <= 1
            sage: int(1) == b[0] == int(2)
            x_0 == 1 == 2
            sage: float(0) <= b[0] <= int(1) <= b[1] <= float(2)
            0 <= x_0 <= 1 <= x_1 <= 2
        """
        global chained_comparator_left
        global chained_comparator_right
        global chained_comparator_replace

        # Ensure that we use chained_comparator_replace only for this
        # call to __richcmp__.
        cdef LinearConstraint replace = chained_comparator_replace
        chained_comparator_replace = None

        # At least one of py_left or py_right must be of type
        # LinearFunctionOrConstraint: figure out its parent
        try:
            LC = (<LinearFunctionOrConstraint?>py_left)._parent
        except TypeError:
            LC = (<LinearFunctionOrConstraint>py_right)._parent

        # We want the parent to be a LinearConstraintsParent
        if not isinstance(LC, LinearConstraintsParent_class):
            LC = LinearConstraintsParent(LC)

        # Check op and change >= to <=
        cdef bint equality = False
        if op == Py_LE:
            pass
        elif op == Py_GE:
            py_left, py_right = py_right, py_left
        elif op == Py_EQ:
            equality = True
        elif op == Py_LT:
            raise ValueError("strict < is not allowed, use <= instead")
        elif op == Py_GT:
            raise ValueError("strict > is not allowed, use >= instead")
        else:
            raise ValueError("inequality != is not allowed, use one of <=, ==, >=")

        # Convert py_left and py_right to constraints left and right
        cdef LinearConstraint left
        cdef LinearConstraint right

        try:
            left = <LinearConstraint?>py_left
        except TypeError:
            left = <LinearConstraint>LC(py_left, equality=equality)
        else:
            if left._parent is not LC:
                left = <LinearConstraint>LC(left.constraints, equality=left.equality)

        try:
            right = <LinearConstraint?>py_right
        except TypeError:
            right = <LinearConstraint>LC(py_right, equality=equality)
        else:
            if right._parent is not LC:
                right = <LinearConstraint>LC(right.constraints, equality=right.equality)

        # HACK to allow chained constraints: Python translates
        # x <= y <= z into:
        #
        #     temp = x <= y      # calls x.__richcmp__(y)
        #     if temp:           # calls temp.__nonzero__()
        #         return y <= z  # calls y.__richcmp__(z)
        #     else:
        #         return temp
        #
        # or, if x <= y is not implemented (for example, if x is a
        # non-Sage type):
        #
        #     temp = y >= x      # calls y.__richcmp__(x)
        #     if temp:           # calls temp.__nonzero__()
        #         return y <= z  # calls y.__richcmp__(z)
        #     else:
        #         return temp
        #
        # but we would like x <= y <= z as output. The trick to make it
        # work is to store x and y in the first call to __richcmp__
        # and temp in the call to __nonzero__. Then we can replace x
        # or y by x <= y in the second call to __richcmp__.
        if replace is not None:
            # First, check for correctly-chained inequalities
            # x <= y <= z or z <= y <= x.
            if py_left is chained_comparator_right:
                left = replace
            elif py_right is chained_comparator_left:
                right = replace
            # Next, check for wrongly-chained inequalities like
            # x <= y >= z. In case of equality, these are accepted
            # anyway.
            elif py_left is chained_comparator_left:
                if not equality:
                    raise ValueError("incorrectly chained inequality")
                left = replace
            elif py_right is chained_comparator_right:
                if not equality:
                    raise ValueError("incorrectly chained inequality")
                right = replace

        chained_comparator_left = py_left
        chained_comparator_right = py_right

        # Finally, chain the (in)equalities
        if left.equality != equality or right.equality != equality:
            raise ValueError("cannot mix equations and inequalities")
        return LC(left.constraints + right.constraints, equality=equality)

    cdef void chained_comparator_hack_nonzero(self):
        """
        Hack to allow chained (in)equalities, see ``__richcmp__``.
        """
        global chained_comparator_replace
        chained_comparator_replace = self

    def __hash__(self):
        r"""
        Return a hash from the ``id()``.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: f = p({2 : 5, 3 : 2})
            sage: f.__hash__()   # random output
            103987752
            sage: d = {}
            sage: d[f] = 3
        """
        # see __cmp__() if you want to change the hash function
        return hash_by_id(<void*>self)

    def __cmp__(left, right):
        """
        Implement comparison of two linear functions or constraints.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: f = p({2 : 5, 3 : 2})
            sage: cmp(f, f)
            0
            sage: abs(cmp(f, f+0))     # since we are comparing by id()
            1
            sage: abs(cmp(f, f+1))
            1
            sage: len(set([f, f]))
            1
            sage: len(set([f, f+0]))
            2
            sage: len(set([f, f+1]))
            2
            sage: abs(cmp(f <= 0, f <= 0))
            1
        """
        # Note: if you want to implement smarter comparison, you also
        # need to change __hash__(). The comparison function must
        # satisfy cmp(x,y)==0 => hash(x)==hash(y)
        if left is right:
            return 0
        if <size_t><void*>left < <size_t><void*>right:
            return -1
        else:
            return 1


#*****************************************************************************
#
# Parent of linear functions
#
#*****************************************************************************

cdef class LinearFunctionsParent_class(Parent):
    r"""
    The parent for all linear functions over a fixed base ring.

    .. warning::

        You should use :func:`LinearFunctionsParent` to construct
        instances of this class.

    INPUT/OUTPUT:

    See :func:`LinearFunctionsParent`

    EXAMPLES::

        sage: from sage.numerical.linear_functions import LinearFunctionsParent_class
        sage: LinearFunctionsParent_class
        <type 'sage.numerical.linear_functions.LinearFunctionsParent_class'>
    """
    def __cinit__(self):
        """
        Cython initializer

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent_class
            sage: LF = LinearFunctionsParent_class.__new__(LinearFunctionsParent_class)
            sage: LF._multiplication_symbol
            '*'
        """
        # Do not use coercion framework for __richcmp__
        self.flags |= Parent_richcmp_element_without_coercion
        self._multiplication_symbol = '*'

    def __init__(self, base_ring):
        """
        The Python constructor

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LinearFunctionsParent(RDF)
            Linear functions over Real Double Field
        """
        from sage.categories.modules_with_basis import ModulesWithBasis
        Parent.__init__(self, base=base_ring, category=ModulesWithBasis(base_ring))

    def set_multiplication_symbol(self, symbol='*'):
        """
        Set the multiplication symbol when pretty-printing linear functions.

        INPUT:

        - ``symbol`` -- string, default: ``'*'``. The multiplication
          symbol to be used.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: f = -1-2*x[0]-3*x[1]
            sage: LF = f.parent()
            sage: LF._get_multiplication_symbol()
            '*'
            sage: f
            -1 - 2*x_0 - 3*x_1
            sage: LF.set_multiplication_symbol(' ')
            sage: f
            -1 - 2 x_0 - 3 x_1
            sage: LF.set_multiplication_symbol()
            sage: f
            -1 - 2*x_0 - 3*x_1
        """
        self._multiplication_symbol = symbol

    def _get_multiplication_symbol(self):
        """
        Return the multiplication symbol.

        OUTPUT:

        String. By default, this is ``'*'``.

        EXAMPLES::

            sage: LF = MixedIntegerLinearProgram().linear_functions_parent()
            sage: LF._get_multiplication_symbol()
            '*'
        """
        return self._multiplication_symbol

    def tensor(self, free_module):
        """
        Return the tensor product with ``free_module``.

        INPUT:

        - ``free_module`` -- vector space or matrix space over the
          same base ring.

        OUTPUT:

        Instance of
        :class:`sage.numerical.linear_tensor.LinearTensorParent_class`.

        EXAMPLES::

            sage: LF = MixedIntegerLinearProgram().linear_functions_parent()
            sage: LF.tensor(RDF^3)
            Tensor product of Vector space of dimension 3 over Real Double Field
            and Linear functions over Real Double Field
            sage: LF.tensor(QQ^2)
            Traceback (most recent call last):
            ...
            ValueError: base rings must match
        """
        from sage.numerical.linear_tensor import LinearTensorParent
        return LinearTensorParent(free_module, self)

    def gen(self, i):
        """
        Return the linear variable `x_i`.

        INPUT:

        - ``i`` -- non-negative integer.

        OUTPUT:

        The linear function `x_i`.

        EXAMPLES::

            sage: LF = MixedIntegerLinearProgram().linear_functions_parent()
            sage: LF.gen(23)
            x_23
        """
        return LinearFunction(self, {i:1})

    def _repr_(self):
        """
        Return as string representation

        EXAMPLES::

            sage: MixedIntegerLinearProgram().linear_functions_parent()
            Linear functions over Real Double Field
        """
        return 'Linear functions over '+str(self.base_ring())

    cpdef _element_constructor_(self, x):
        """
        Construt a :class:`LinearFunction` from ``x``.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: LF = p.linear_functions_parent()
            sage: LF._element_constructor_(123)
            123
            sage: p(123)    # indirect doctest
            123
            sage: type(_)
            <type 'sage.numerical.linear_functions.LinearFunction'>

            sage: p_QQ = MixedIntegerLinearProgram(solver='ppl')
            sage: LF_QQ = p_QQ.linear_functions_parent()
            sage: f = LF({-1:1/2, 2:3/4});  f
            0.5 + 0.75*x_2
            sage: LF(f) is f
            True
            sage: LF_QQ(f)
            1/2 + 3/4*x_2
            sage: LF_QQ(f) is f
            False
        """
        if is_LinearFunction(x):
            if x.parent() is self:
                return x
            else:
                return LinearFunction(self, (<LinearFunction>x)._f)
        return LinearFunction(self, x)

    cpdef _coerce_map_from_(self, R):
        """
        Allow coercion of scalars into linear functions.

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

        A linear function.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent()
            sage: p._an_element_()
            5*x_2 + 7*x_5
            sage: p.an_element()   # indirect doctest
            5*x_2 + 7*x_5
        """
        return self._element_constructor_({2:5, 5:7})


#*****************************************************************************
#
# Elements of linear functions
#
#*****************************************************************************

cdef class LinearFunction(LinearFunctionOrConstraint):
    r"""
    An elementary algebra to represent symbolic linear functions.

    .. warning::

        You should never instantiate :class:`LinearFunction`
        manually. Use the element constructor in the parent
        instead. For convenience, you can also call the
        :class:`MixedIntegerLinearProgram` instance directly.

    EXAMPLES:

    For example, do this::

        sage: p = MixedIntegerLinearProgram()
        sage: p({0 : 1, 3 : -8})
        x_0 - 8*x_3

    or this::

        sage: parent = p.linear_functions_parent()
        sage: parent({0 : 1, 3 : -8})
        x_0 - 8*x_3

    instead of this::

        sage: from sage.numerical.linear_functions import LinearFunction
        sage: LinearFunction(p.linear_functions_parent(), {0 : 1, 3 : -8})
        x_0 - 8*x_3
    """

    def __init__(self, parent, f):
        r"""
        Constructor taking a dictionary or a numerical value as its argument.

        A linear function is represented as a dictionary. The
        values are the coefficient of the variable represented
        by the keys ( which are integers ). The key ``-1``
        corresponds to the constant term.

        EXAMPLES:

        With a dictionary::

            sage: p = MixedIntegerLinearProgram()
            sage: p({0 : 1, 3 : -8})
            x_0 - 8*x_3

        Using the constructor with a numerical value::

            sage: p = MixedIntegerLinearProgram()
            sage: p(25)
            25
        """
        ModuleElement.__init__(self, parent)
        R = self.base_ring()
        if isinstance(f, dict):
            self._f = dict( (int(key), R(value)) for key, value in f.iteritems() )
        else:
            self._f = {-1: R(f)}

    cpdef iteritems(self):
        """
        Iterate over the index, coefficient pairs

        OUTPUT:

        An iterator over the ``(key, coefficient)`` pairs. The keys
        are integers indexing the variables. The key ``-1``
        corresponds to the constant term.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram(solver = 'ppl')
            sage: x = p.new_variable()
            sage: f = 0.5 + 3/2*x[1] + 0.6*x[3]
            sage: for id, coeff in f.iteritems():
            ....:    print 'id =', id, '  coeff =', coeff
            id = 0   coeff = 3/2
            id = 1   coeff = 3/5
            id = -1   coeff = 1/2
        """
        return self._f.iteritems()

    def dict(self):
        r"""
        Return the dictionary corresponding to the Linear Function.

        OUTPUT:

        The linear function is represented as a dictionary. The value
        are the coefficient of the variable represented by the keys (
        which are integers ). The key ``-1`` corresponds to the
        constant term.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: lf = p({0 : 1, 3 : -8})
            sage: lf.dict()
            {0: 1.0, 3: -8.0}
        """
        return dict(self._f)

    def coefficient(self, x):
        r"""
        Return one of the the coefficients.

        INPUT:

        - ``x`` -- a linear variable or an integer. If an integer `i`
          is passed, then `x_i` is used as linear variable.

        OUTPUT:

        A base ring element. The coefficient of ``x`` in the linear
        function. Pass ``-1`` for the constant term.

        EXAMPLE::

            sage: mip.<b> = MixedIntegerLinearProgram()
            sage: lf = -8 * b[3] + b[0] - 5;  lf
            -5 - 8*x_0 + x_1
            sage: lf.coefficient(b[3])
            -8.0
            sage: lf.coefficient(0)      # x_0 is b[3]
            -8.0
            sage: lf.coefficient(4)
            0.0
            sage: lf.coefficient(-1)
            -5.0

        TESTS::

            sage: lf.coefficient(b[3] + b[4])
            Traceback (most recent call last):
            ...
            ValueError: x is a sum, must be a single variable
            sage: lf.coefficient(2*b[3])
            Traceback (most recent call last):
            ...
            ValueError: x must have a unit coefficient
            sage: mip.<q> = MixedIntegerLinearProgram(solver='ppl')
            sage: lf.coefficient(q[0])
            Traceback (most recent call last):
            ...
            ValueError: x is from a different linear functions module
        """
        if is_LinearFunction(x):
            if self.parent() != x.parent():
                raise ValueError('x is from a different linear functions module')
            if len((<LinearFunction>x)._f) != 1:
                raise ValueError('x is a sum, must be a single variable')
            i = (<LinearFunction>x)._f.keys()[0]
            if (<LinearFunction>x)._f[i] != 1:
                raise ValueError('x must have a unit coefficient')
        else:
            i = int(x)
        try:
            return self._f[i]
        except KeyError:
            return self.parent().base_ring().zero()

    cpdef ModuleElement _add_(self, ModuleElement b):
        r"""
        Defining the + operator

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p({0 : 1, 3 : -8}) + p({2 : 5, 3 : 2}) - 16
            -16 + x_0 + 5*x_2 - 6*x_3
        """
        e = dict(self._f)
        for (id,coeff) in b.dict().iteritems():
            e[id] = self._f.get(id,0) + coeff
        P = self.parent()
        return P(e)

    cpdef ModuleElement _neg_(self):
        r"""
        Defining the - operator (opposite).

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: - p({0 : 1, 3 : -8})
            -1*x_0 + 8*x_3
        """
        P = self.parent()
        return P(dict([(id,-coeff) for (id, coeff) in self._f.iteritems()]))

    cpdef ModuleElement _sub_(self, ModuleElement b):
        r"""
        Defining the - operator (substraction).

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p({2 : 5, 3 : 2}) - 3
            -3 + 5*x_2 + 2*x_3
            sage: p({0 : 1, 3 : -8}) - p({2 : 5, 3 : 2}) - 16
            -16 + x_0 - 5*x_2 - 10*x_3
        """
        e = dict(self._f)
        for (id,coeff) in b.dict().iteritems():
            e[id] = self._f.get(id,0) - coeff
        P = self.parent()
        return P(e)

    cpdef ModuleElement _rmul_(self, RingElement b):
        r"""
        Left multiplication by scalars

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: p({2 : 5, 3 : 2}) * 3
            15*x_2 + 6*x_3
        """
        P = self.parent()
        return P(dict([(id,b*coeff) for (id, coeff) in self._f.iteritems()]))

    cpdef ModuleElement _lmul_(self, RingElement b):
        r"""
        Right multiplication by scalars

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: 3 * p({2 : 5, 3 : 2})
            15*x_2 + 6*x_3
        """
        return self._rmul_(b)

    cpdef _acted_upon_(self, x, bint self_on_left):
       """
       Act with scalars that do not have a natural coercion into
       ``self.base_ring()``

       EXAMPLES:

       Note that there is no coercion from ``RR`` to ``QQ``, but you
       can explicitly convert them::

           sage: 1/2 * 0.6
           0.300000000000000

       If there were a coercion, then 0.6 would have been coerced into
       6/10 and the result would have been rational. This is
       undesirable, which is why there cannot be a coercion on the
       level of the coefficient rings.

       By declaring an action of ``RR`` on the linear functions over
       ``QQ`` we make the following possible::

           sage: p = MixedIntegerLinearProgram(solver='ppl')
           sage: p.base_ring()
           Rational Field
           sage: x = p.new_variable()
           sage: x[0] * 0.6
           3/5*x_0

           sage: vf = (2 + x[0]) * vector(ZZ, [3,4]);  vf
           (6, 8) + (3, 4)*x_0
           sage: vf.parent()
           Tensor product of Vector space of dimension 2 over Rational Field
           and Linear functions over Rational Field

           sage: tf = x[0] * identity_matrix(2);  tf
           [x_0 0  ]
           [0   x_0]
           sage: tf.parent()
           Tensor product of Full MatrixSpace of 2 by 2 dense matrices over
           Rational Field and Linear functions over Rational Field
       """
       R = self.base_ring()
       try:
           x_R = R(x)
       except TypeError:
           M = x.parent().base_extend(R)
           x_M = M(x)
           from sage.numerical.linear_tensor import LinearTensorParent
           P = LinearTensorParent(M, self.parent())
           tensor = dict()
           for k, v in self._f.items():
               tensor[k] = x_M * v
           return P(tensor)
       return self._rmul_(x_R)

    def _coeff_formatter(self, coeff, constant_term=False):
        """
        Pretty-print multiplicative coefficient ``x``

        OUTPUT:

        String, including a trailing space if the coefficient is not
        one. Empty string otherwise.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: f = p(1);  type(f)
            <type 'sage.numerical.linear_functions.LinearFunction'>
            sage: f._coeff_formatter(1)
            ''
            sage: f._coeff_formatter(1, constant_term=True)
            '1'
            sage: f._coeff_formatter(RDF(12.0))
            '12*'
            sage: f._coeff_formatter(RDF(12.3))
            '12.3*'

            sage: p = MixedIntegerLinearProgram(solver='ppl')
            sage: f = p(1)
            sage: f._coeff_formatter(13/45)
            '13/45*'

            sage: from sage.rings.number_field.number_field import QuadraticField
            sage: K.<sqrt5> = QuadraticField(5, 'sqrt5')
            sage: p = MixedIntegerLinearProgram(solver='interactivelp', base_ring=K)
            sage: f = p(1)
            sage: f._coeff_formatter(sqrt5)
            'sqrt5*'

            sage: from sage.rings.all import AA
            sage: sqrt5 = AA(5).sqrt()
            sage: p = MixedIntegerLinearProgram(solver='interactivelp', base_ring=AA)
            sage: f = p(1)
            sage: f._coeff_formatter(sqrt5)
            '2.236067977499790?*'
        """
        R = self.base_ring()
        if coeff == R.one() and not constant_term:
            return ''
        try:
            from sage.rings.all import ZZ
            coeff = ZZ(coeff)    # print as integer if possible
        except (TypeError, ValueError):
            pass
        if constant_term:
            return str(coeff)
        else:
            return str(coeff) + self.parent()._multiplication_symbol

    def _repr_(self):
        r"""
        Returns a string version of the linear function.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram(solver='GLPK')
            sage: p({-1: -15, 2 : -5.1, 3 : 2/3})
            -15 - 5.1*x_2 + 0.666666666667*x_3
            sage: p = MixedIntegerLinearProgram(solver='ppl')
            sage: p({-1: -15, 2 : -5.1, 3 : 2/3})
            -15 - 51/10*x_2 + 2/3*x_3
        """
        cdef dict d = dict(self._f)
        cdef bint first = True
        t = ""

        if -1 in d:
            coeff = d.pop(-1)
            if coeff!=0:
                t = self._coeff_formatter(coeff, constant_term=True)
                first = False

        cdef list l = sorted(d.items())
        for id,coeff in l:
            sign = cmp(coeff,0)
            if sign == 0:
                continue
            if not first:
                if sign == -1:
                    t += ' - '
                if sign == +1:
                    t += ' + '
                t += self._coeff_formatter(abs(coeff)) + 'x_' + str(id)
            else:
                t += self._coeff_formatter(coeff) + 'x_' + str(id)
            first = False

        if first:
            return '0'
        else:
            return t

    cpdef is_zero(self):
        """
        Test whether ``self`` is zero.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: (x[1] - x[1] + 0*x[2]).is_zero()
            True
        """
        for coeff in self._f.values():
            if not coeff.is_zero():
                return False
        return True

    cpdef equals(LinearFunction left, LinearFunction right):
        """
        Logically compare ``left`` and ``right``.

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: (x[1] + 1).equals(3/3 + 1*x[1] + 0*x[2])
            True
        """
        return (left-right).is_zero()


#*****************************************************************************
#
# Parent of linear constraints
#
#*****************************************************************************

cdef class LinearConstraintsParent_class(Parent):
    """
    Parent for :class:`LinearConstraint`

    .. warning::

        This class has no reason to be instantiated by the user, and
        is meant to be used by instances of
        :class:`MixedIntegerLinearProgram`. Also, use the
        :func:`LinearConstraintsParent` factory function.

    INPUT/OUTPUT:

        See :func:`LinearFunctionsParent`

    EXAMPLES::

        sage: p = MixedIntegerLinearProgram()
        sage: LC = p.linear_constraints_parent();  LC
        Linear constraints over Real Double Field
        sage: from sage.numerical.linear_functions import LinearConstraintsParent
        sage: LinearConstraintsParent(p.linear_functions_parent()) is LC
        True
    """
    def __cinit__(self, linear_functions_parent):
        """
        Cython initializer

        TESTS::

            sage: from sage.numerical.linear_functions import LinearConstraintsParent_class
            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LF = LinearFunctionsParent(RDF)
            sage: LinearConstraintsParent_class.__new__(LinearConstraintsParent_class, LF)
            Linear constraints over Real Double Field
            sage: LinearConstraintsParent_class.__new__(LinearConstraintsParent_class, None)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert NoneType to sage.numerical.linear_functions.LinearFunctionsParent_class
        """
        self._LF = <LinearFunctionsParent_class?>linear_functions_parent
        # Do not use coercion framework for __richcmp__
        self.flags |= Parent_richcmp_element_without_coercion

    def __init__(self, linear_functions_parent):
        """
        The Python constructor

        INPUT/OUTPUT:

        See :func:`LinearFunctionsParent`

        TESTS::

            sage: from sage.numerical.linear_functions import LinearFunctionsParent
            sage: LF = LinearFunctionsParent(RDF)
            sage: from sage.numerical.linear_functions import LinearConstraintsParent
            sage: LinearConstraintsParent(LF)
            Linear constraints over Real Double Field
            sage: LinearConstraintsParent(None)
            Traceback (most recent call last):
            ...
            TypeError: Cannot convert NoneType to sage.numerical.linear_functions.LinearFunctionsParent_class
        """
        Parent.__init__(self)

    def linear_functions_parent(self):
        """
        Return the parent for the linear functions

        EXAMPLES::

            sage: LC = MixedIntegerLinearProgram().linear_constraints_parent()
            sage: LC.linear_functions_parent()
            Linear functions over Real Double Field
        """
        return self._LF

    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String.

        EXAMPLES::

            sage: MixedIntegerLinearProgram().linear_constraints_parent()
            Linear constraints over Real Double Field
        """
        return 'Linear constraints over '+str(self.linear_functions_parent().base_ring())

    cpdef _element_constructor_(self, left, right=None, equality=False):
        """
        Construt a :class:`LinearConstraint`.

        INPUT:

        - ``left`` -- a :class:`LinearFunction`, or something that can
          be converted into one, a list/tuple of
          :class:`LinearFunction`, or an existing
          :class:`LinearConstraint`.

        - ``right`` -- a :class:`LinearFunction` or ``None``
          (default).

        - ``equality`` -- boolean (default: ``True``). Whether to
          construct an equation or an inequality.

        OUTPUT:

        The :class:`LinearConstraint` constructed from the input data.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: LC = p.linear_constraints_parent()
            sage: LC._element_constructor_(1, 2)
            1 <= 2

            sage: x = p.new_variable()
            sage: LC([x[0], x[1], x[2]])
            x_0 <= x_1 <= x_2

            sage: LC([x[0], x[1], x[2]], equality=True)
            x_0 == x_1 == x_2

            sage: type(_)
            <type 'sage.numerical.linear_functions.LinearConstraint'>

        TESTS::

            sage: inequality = LC([x[0], 1/2*x[1], 3/4*x[2]]); inequality
            x_0 <= 0.5*x_1 <= 0.75*x_2
            sage: LC(inequality) is inequality
            True
            sage: p_QQ = MixedIntegerLinearProgram(solver='ppl')
            sage: LC_QQ = p_QQ.linear_constraints_parent()
            sage: LC_QQ(inequality)
            x_0 <= 1/2*x_1 <= 3/4*x_2
            sage: LC_QQ(inequality) is inequality
            False
        """
        if right is None and is_LinearConstraint(left):
            if (left.parent() is self) and (left.is_equation() == equality):
                return left
            else:
                return LinearConstraint(self, (<LinearConstraint>left).constraints,
                                        equality=equality)
        if right is None:
            if isinstance(left, (list,tuple)):
                return LinearConstraint(self, left, equality=equality)
            else:
                return LinearConstraint(self, [left], equality=equality)
        else:
            return LinearConstraint(self, [left, right], equality=equality)

    cpdef _coerce_map_from_(self, R):
        """
        Allow coercion of scalars into linear functions.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: parent = p.linear_constraints_parent()
            sage: parent.coerce(int(2))
            trivial constraint starting with 2
            sage: parent._coerce_map_from_(int)
            True
        """
        return self.linear_functions_parent().has_coerce_map_from(R)

    def _an_element_(self):
        """
        Returns an element

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram().linear_functions_parent()
            sage: p._an_element_()
            5*x_2 + 7*x_5
            sage: p.an_element()   # indirect doctest
            5*x_2 + 7*x_5
        """
        LF = self.linear_functions_parent()
        return self(0) <= LF.an_element()


#*****************************************************************************
#
# Elements of linear constraints
#
#*****************************************************************************

cdef class LinearConstraint(LinearFunctionOrConstraint):
    """
    A class to represent formal Linear Constraints.

    A Linear Constraint being an inequality between
    two linear functions, this class lets the user
    write ``LinearFunction1 <= LinearFunction2``
    to define the corresponding constraint, which
    can potentially involve several layers of such
    inequalities (``A <= B <= C``), or even equalities
    like ``A == B == C``.

    Trivial constraints (meaning that they have only one term and no
    relation) are also allowed. They are required for the coercion
    system to work.

    .. warning::

        This class has no reason to be instantiated by the user, and
        is meant to be used by instances of
        :class:`MixedIntegerLinearProgram`.

    INPUT:

    - ``parent`` -- the parent, a :class:`LinearConstraintsParent_class`

    - ``terms`` -- a list/tuple/iterable of two or more linear
      functions (or things that can be converted into linear
      functions).

    - ``equality`` -- boolean (default: ``False``). Whether the terms
      are the entries of a chained less-or-equal (``<=``) inequality
      or a chained equality.

    EXAMPLE::

        sage: p = MixedIntegerLinearProgram()
        sage: b = p.new_variable()
        sage: b[2]+2*b[3] <= b[8]-5
        x_0 + 2*x_1 <= -5 + x_2
    """

    def __init__(self, parent, terms, equality=False):
        r"""
        Constructor for ``LinearConstraint``

        INPUT:

        See :class:`LinearConstraint`.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: b[2]+2*b[3] <= b[8]-5
            x_0 + 2*x_1 <= -5 + x_2
        """
        assert len(terms) > 0
        super(LinearConstraint, self).__init__(parent)
        self.equality = equality
        LF = parent.linear_functions_parent()
        self.constraints = [ LF(term) for term in terms ]

    cpdef equals(LinearConstraint left, LinearConstraint right):
        """
        Compare ``left`` and ``right``.

        OUTPUT:

        Boolean. Whether all terms of ``left`` and ``right`` are
        equal. Note that this is stronger than mathematical
        equivalence of the relations.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: x = p.new_variable()
            sage: (x[1] + 1 >= 2).equals(3/3 + 1*x[1] + 0*x[2] >= 8/4)
            True
            sage: (x[1] + 1 >= 2).equals(x[1] + 1-1 >= 1-1)
            False
        """
        if len(left.constraints) != len(right.constraints):
            return False
        if left.equality != right.equality:
            return False
        cdef LinearFunction l, r
        for i in range(len(left.constraints)):
            l = <LinearFunction>(left.constraints[i])
            r = <LinearFunction>(right.constraints[i])
            if not l.equals(r):
                return False
        return True

    def is_equation(self):
        """
        Whether the constraint is a chained equation

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: (b[0] == b[1]).is_equation()
            True
            sage: (b[0] <= b[1]).is_equation()
            False
        """
        return self.equality

    def is_less_or_equal(self):
        """
        Whether the constraint is a chained less-or_equal inequality

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: (b[0] == b[1]).is_less_or_equal()
            False
            sage: (b[0] <= b[1]).is_less_or_equal()
            True
        """
        return not self.equality

    def is_trivial(self):
        """
        Test whether the constraint is trivial.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: LC = p.linear_constraints_parent()
            sage: ieq = LC(1,2);  ieq
            1 <= 2
            sage: ieq.is_trivial()
            False

            sage: ieq = LC(1);  ieq
            trivial constraint starting with 1
            sage: ieq.is_trivial()
            True
        """
        return len(self.constraints) < 2

    def __iter__(self):
        """
        Iterate over the terms of the chained (in)-equality

        OUTPUT:

        A generator yielding the individual terms of the constraint in
        left-to-right order.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: ieq = 1 <= b[0] <= b[2] <= 3 <= b[3];  ieq
            1 <= x_0 <= x_1 <= 3 <= x_2
            sage: list(ieq)
            [1, x_0, x_1, 3, x_2]
            sage: for term in ieq:
            ....:     print term
            1
            x_0
            x_1
            3
            x_2
        """
        for term in self.constraints:
            yield term

    def equations(self):
        """
        Iterate over the unchained(!) equations

        OUTPUT:

        An iterator over pairs ``(lhs, rhs)`` such that the individual
        equations are ``lhs == rhs``.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: eqns = 1 == b[0] == b[2] == 3 == b[3];  eqns
            1 == x_0 == x_1 == 3 == x_2
            sage: for lhs, rhs in eqns.equations():
            ....:     print str(lhs) + ' == ' + str(rhs)
            1 == x_0
            x_0 == x_1
            x_1 == 3
            3 == x_2
        """
        if not self.is_equation() or self.is_trivial():
            raise StopIteration
        term_iter = iter(self)
        lhs = next(term_iter)
        rhs = next(term_iter)
        while True:
            yield (lhs, rhs)
            lhs = rhs
            rhs = next(term_iter)

    def inequalities(self):
        """
        Iterate over the unchained(!) inequalities

        OUTPUT:

        An iterator over pairs ``(lhs, rhs)`` such that the individual
        equations are ``lhs <= rhs``.

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: ieq = 1 <= b[0] <= b[2] <= 3 <= b[3]; ieq
            1 <= x_0 <= x_1 <= 3 <= x_2

            sage: for lhs, rhs in ieq.inequalities():
            ....:     print str(lhs) + ' <= ' + str(rhs)
            1 <= x_0
            x_0 <= x_1
            x_1 <= 3
            3 <= x_2
        """
        if not self.is_less_or_equal() or self.is_trivial():
            raise StopIteration
        term_iter = iter(self)
        lhs = next(term_iter)
        rhs = next(term_iter)
        while True:
            yield (lhs, rhs)
            lhs = rhs
            rhs = next(term_iter)

    def _repr_(self):
        r"""
        Returns a string representation of the constraint.

        OUTPUT:

        String.

        EXAMPLE::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: b[3] <= b[8] + 9
            x_0 <= 9 + x_1

            sage: LC = p.linear_constraints_parent()
            sage: LC(b[3], b[8] + 9)
            x_0 <= 9 + x_1
            sage: LC(b[3])
            trivial constraint starting with x_0
        """
        comparator = ( ' == ' if self.equality else ' <= ' )
        result = comparator.join(map(str, self))
        if self.is_trivial():
            return 'trivial constraint starting with '+result
        return result

    def __nonzero__(self):
        """
        Part of the hack to allow chained (in)equalities

        EXAMPLES::

            sage: p = MixedIntegerLinearProgram()
            sage: b = p.new_variable()
            sage: ieq = (b[3] <= b[8] + 9)
            sage: ieq <= ieq <= ieq
            x_0 <= 9 + x_1 <= x_0 <= 9 + x_1 <= x_0 <= 9 + x_1
        """
        self.chained_comparator_hack_nonzero()
        return True
