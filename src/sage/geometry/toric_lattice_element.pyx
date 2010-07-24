r"""
Toric lattice elements

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.generic.toric_variety`,
:mod:`~sage.schemes.generic.fano_toric_variety`).

AUTHORS:

- Andrey Novoseltsev (2010-05-27): initial version.

TESTS:

Let's create some lattices first::

    sage: N = ToricLattice(3)
    sage: M = N.dual()

Now we are ready to create elements of toric lattices::

    sage: n = N([1,2,3])
    sage: n
    N(1, 2, 3)
    sage: m = M(1,2,3)
    sage: m
    M(1, 2, 3)

Dual lattices can act on each other::

    sage: n * m
    14
    sage: m * n
    14

You can also add elements of the same lattice or scale them::

    sage: 2 * n
    N(2, 4, 6)
    sage: n * 2
    N(2, 4, 6)
    sage: n + n
    N(2, 4, 6)

However, you cannot "mix wrong lattices" in your expressions::

    sage: n + m
    Traceback (most recent call last):
    ...
    TypeError: unsupported operand parent(s) for '+':
    '3-d lattice N' and '3-d lattice M'
    sage: n * n
    Traceback (most recent call last):
    ...
    TypeError: elements of the same toric lattice cannot be multiplied!
    sage: n == m
    False

Note that ``n`` and ``m`` are not equal to each other even though they are
both "just (1,2,3)." Moreover, you cannot easily convert elements between
toric lattices::

    sage: M(n)
    Traceback (most recent call last):
    ...
    TypeError: N(1, 2, 3) cannot be converted to 3-d lattice M!

If you really need to consider elements of one lattice as elements of another,
you can either use intermediate conversion to "just a vector"::

    sage: ZZ3 = ZZ^3
    sage: n_in_M = M(ZZ3(n))
    sage: n_in_M
    M(1, 2, 3)
    sage: n == n_in_M
    False
    sage: n_in_M == m
    True

Or you can create a homomorphism from one lattice to any other::

    sage: h = N.hom(identity_matrix(3), M)
    sage: h(n)
    M(1, 2, 3)
"""
# The "tutorial" above is a truncated version of one in toric_lattice.py.


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


include '../ext/stdsage.pxi' # Needed for PY_NEW


from sage.modules.vector_integer_dense cimport Vector_integer_dense
from sage.structure.coerce_exceptions import CoercionException
from sage.structure.element cimport Element, Vector


def is_ToricLatticeElement(x):
    r"""
    Check if ``x`` is an element of a toric lattice.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is an element of a toric lattice, ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.toric_lattice_element import (
        ...     is_ToricLatticeElement)
        sage: is_ToricLatticeElement(1)
        False
        sage: e = ToricLattice(3).an_element()
        sage: e
        N(1, 0, 0)
        sage: is_ToricLatticeElement(e)
        True
    """
    return isinstance(x, ToricLatticeElement)


# Why do we need a special class:
# - customize output to include lattice name
# - prohibit operations mixing "wrong" lattices
cdef class ToricLatticeElement(Vector_integer_dense):
    r"""
    Create an element of a toric lattice.

    .. WARNING::

        You probably should not construct such elements explicitly.

    INPUT:

    - same as for
      :class:`~sage.modules.vector_integer_dense.Vector_integer_dense`.

    OUTPUT:

    - element of a toric lattice.

    TESTS::

        sage: N = ToricLattice(3)
        sage: from sage.geometry.toric_lattice_element import (
        ...             ToricLatticeElement)
        sage: e = ToricLatticeElement(N, [1,2,3])
        sage: e
        N(1, 2, 3)
        sage: TestSuite(e).run()
    """

    # We do not add any new functionality, we actually limit the existing one
    # instead. In particular, there is no need in __init__, but the following
    # function ensures that _add_ etc. return ToricLatticeElement, rather
    # than Vector_integer_dense. Its code is copied from the base class with
    # the type name replacements.
    # It is not detected by doctest coverage, the original has no
    # documentation, and I don't feel I can clearly define the exact part of
    # the initialization process which is done by it. So I will leave it
    # without any further documentation as well...
    cdef _new_c(self):
        cdef ToricLatticeElement y
        y = PY_NEW(ToricLatticeElement)
        y._init(self._degree, self._parent)
        return y

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is an equal element of the same toric lattice as
          ``self``, 1 or -1 otherwise.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: cmp(n, m)
            1
            sage: n2 = N(1,2,3)
            sage: cmp(n, n2)
            0
            sage: n is n2
            False
            sage: n == 1
            False
        """
        c = cmp(type(self), type(right))
        if c:
            return c
        c = cmp(self.parent(), right.parent())
        if c:
            return c
        # Now use the real comparison of vectors
        return self._cmp_c_impl(right)

    # For some reason, vectors work just fine without redefining this function
    # from the base class, but if it is not here, we get "unhashable type"...
    def __hash__(self):
        r"""
        Return the hash of ``self``.

        OUTPUT:

        - integer.

        TESTS::

            sage: N = ToricLattice(3)
            sage: n = N(1,2,3)
            sage: hash(n)
            Traceback (most recent call last):
            ...
            TypeError: mutable vectors are unhashable
            sage: n.set_immutable()
            sage: hash(n) == hash(n)
            True
        """
        return Vector_integer_dense.__hash__(self)

    cpdef _act_on_(self, other, bint self_on_left):
        """
        Act on ``other``.

        INPUT:

        - ``other`` - :class:`ToricLatticeElement`.

        OUTPUT:

        - integer, if ``other`` is an element of the dual lattice of ``self``;

        - ``CoercionException`` is raised if ``other`` is an element of
          an incompatible toric lattice;

        - standard output for ``self`` acting as an integral vector on
          ``other`` if the latter one is not an element of a toric lattice.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: n * m # indirect doctest
            14

        Now we test behaviour with other types::

            sage: v = vector([1, 2, 3])
            sage: v * n == n * v
            True
            sage: v = vector([1, 1/2, 3/4])
            sage: v * n == n * v
            True
            sage: A = matrix(3, range(9))
            sage: A * n
            (8, 26, 44)
            sage: n * A
            (24, 30, 36)
            sage: B = A / 3
            sage: B * n
            (8/3, 26/3, 44/3)
            sage: n * B
            (8, 10, 12)
        """
        # We try to deal only with the case of two lattice elements...
        if is_ToricLatticeElement(other):
            if other.parent() is self.parent().dual():
                # Our own _dot_product_ is disabled
                return Vector_integer_dense._dot_product_(self, other)
            raise CoercionException("only elements of dual toric lattices "
                                    "can act on each other!")
        # ... however we also need to treat the case when other is an integral
        # vector, since otherwise it will be coerced to the parent of self and
        # then the dot product will be called for elements of the same lattice
        if isinstance(other, Vector_integer_dense):
            return Vector_integer_dense._dot_product_(self, other)
        # Now let the standard framework work...
        return Vector_integer_dense._act_on_(self, other, self_on_left)

    # We need to override this function to prohibit default behaviour.
    # It seems to be called when right is in the same lattice as self, which
    # is wrong from our point of view.
    cpdef Element _dot_product_(self, Vector right):
        """
        Raise a ``TypeError`` exception.

        Dot product is not defined on toric lattices (there are actions of
        dual lattices on each other instead).

        INPUT:

        - ``right`` - vector.

        OUTPUT:

        - ``TypeError`` exception is raised.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: n = N(1,2,3)
            sage: m = M(1,2,3)
            sage: n * n # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: elements of the same
            toric lattice cannot be multiplied!
        """
        raise TypeError("elements of the same toric lattice cannot be "
                        "multiplied!")

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: Ld = ToricLattice(3, "L").dual()
            sage: e = Ld(1,2,3)
            sage: e._latex_()
            '\\left(1,2,3\\right)_{L^*}'
        """
        return "%s_{%s}" % (super(ToricLatticeElement, self)._latex_(),
                            self.parent()._latex_name)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: Ld = ToricLattice(3, "L").dual()
            sage: e = Ld(1,2,3)
            sage: e._repr_()
            'L*(1, 2, 3)'
        """
        return self.parent()._name + super(ToricLatticeElement, self)._repr_()
