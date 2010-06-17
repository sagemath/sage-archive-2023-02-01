r"""
Toric lattices

This module was designed as a part of the framework for toric varieties
(:mod:`~sage.schemes.generic.toric_variety`,
:mod:`~sage.schemes.generic.fano_toric_variety`).

All toric lattices are isomorphic to `\ZZ^n` for some `n`, but will prevent
you from doing "wrong" operations with objects from different lattices.

AUTHORS:

- Andrey Novoseltsev (2010-05-27): initial version.

EXAMPLES:

The simplest way to create a toric lattice is to specify its dimension only::

    sage: N = ToricLattice(3)
    sage: N
    3-d lattice N

While our lattice ``N`` is called exactly "N" it is a coincidence: all
lattices are called "N" by default::

    sage: another_name = ToricLattice(3)
    sage: another_name
    3-d lattice N

If fact, the above lattice is exactly the same as before as an object in
memory::

    sage: N is another_name
    True

There are actually four names associated to a toric lattice and they all must
be the same for two lattices to coincide::

    sage: N, N.dual(), latex(N), latex(N.dual())
    (3-d lattice N, 3-d lattice M, N, M)

Notice that the lattice dual to ``N`` is called "M" which is standard in toric
geometry. This happens only if you allow completely automatic handling of
names::

    sage: another_N = ToricLattice(3, "N")
    sage: another_N.dual()
    3-d lattice N*
    sage: N is another_N
    False

What can you do with toric lattices? Well, their main purpose is to allow
creation of elements of toric lattices::

    sage: n = N([1,2,3])
    sage: n
    N(1, 2, 3)
    sage: M = N.dual()
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

.. WARNING::

    While integer vectors (elements of `\ZZ^n`) are printed as ``(1,2,3)``,
    in the code ``(1,2,3)`` is a :class:`tuple`, which has nothing to do
    neither with vectors, nor with toric lattices, so the following is
    probably not what you want while working with toric geometry objects::

        sage: (1,2,3) + (1,2,3)
        (1, 2, 3, 1, 2, 3)

    Instead, use syntax like ::

        sage: N(1,2,3) + N(1,2,3)
        N(2, 4, 6)
"""
# Parts of the "tutorial" above are also in toric_lattice_element.pyx.


#*****************************************************************************
#       Copyright (C) 2010 Andrey Novoseltsev <novoselt@gmail.com>
#       Copyright (C) 2010 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.geometry.toric_lattice_element import (ToricLatticeElement,
                                                 is_ToricLatticeElement)
from sage.modules.free_module import FreeModule_ambient_pid
from sage.rings.all import ZZ
from sage.structure.factory import UniqueFactory


def is_ToricLattice(x):
    r"""
    Check if ``x`` is a toric lattice.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is a toric lattice and ``False`` otherwise.

    EXAMPLES::

        sage: from sage.geometry.toric_lattice import (
        ...     is_ToricLattice)
        sage: is_ToricLattice(1)
        False
        sage: N = ToricLattice(3)
        sage: N
        3-d lattice N
        sage: is_ToricLattice(N)
        True
    """
    return isinstance(x, ToricLatticeClass)


class ToricLatticeFactory(UniqueFactory):
    r"""
    Create a lattice for toric geometry objects.

    INPUT:

    - ``rank`` -- nonnegative integer, the only mandatory parameter;

    - ``name`` -- string;

    - ``dual_name`` -- string;

    - ``latex_name`` -- string;

    - ``latex_dual_name`` -- string.

    OUTPUT:

    - lattice.

    A toric lattice is uniquely determined by its rank and associated names.
    There are four such "associated names" whose meaning should be clear from
    the names of the corresponding parameters, but the choice of default
    values is a little bit involved. So here is the full description of the
    "naming algorithm":

    #. If no names were given at all, then this lattice will be called "N" and
       the dual one "M". These are the standard choices in toric geometry.

    #. If ``name`` was given and ``dual_name`` was not, then ``dual_name``
       will be ``name`` followed by "*".

    #. If LaTeX names were not given, they will coincide with the "usual"
       names, but if ``dual_name`` was constructed automatically, the trailing
       star will be typeset as a superscript.

    EXAMPLES:

    Let's start with no names at all and see how automatic names are given::

        sage: L1 = ToricLattice(3)
        sage: L1
        3-d lattice N
        sage: L1.dual()
        3-d lattice M

    If we give the name "N" explicitly, the dual lattice will be called "N*"::

        sage: L2 = ToricLattice(3, "N")
        sage: L2
        3-d lattice N
        sage: L2.dual()
        3-d lattice N*

    However, we can give an explicit name for it too::

        sage: L3 = ToricLattice(3, "N", "M")
        sage: L3
        3-d lattice N
        sage: L3.dual()
        3-d lattice M

    If you want, you may also give explicit LaTeX names::

        sage: L4 = ToricLattice(3, "N", "M", r"\mathbb{N}", r"\mathbb{M}")
        sage: latex(L4)
        \mathbb{N}
        sage: latex(L4.dual())
        \mathbb{M}

    While all four lattices above are called "N", only two of them are equal
    (and are actually the same)::

        sage: L1 == L2
        False
        sage: L1 == L3
        True
        sage: L1 is L3
        True
        sage: L1 == L4
        False

    The reason for this is that ``L2`` and ``L4`` have different names either
    for dual lattices or for LaTeX typesetting.
    """

    def create_key(self, rank, name=None, dual_name=None,
                   latex_name=None, latex_dual_name=None):
        """
        Create a key that uniquely identifies this toric lattice.

        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        .. WARNING::

            You probably should not use this function directly.

        TESTS::

            sage: ToricLattice.create_key(3)
            (3, 'N', 'M', 'N', 'M')
            sage: N = ToricLattice(3)
            sage: loads(dumps(N)) is N
            True
            sage: TestSuite(N).run()
        """
        rank = int(rank)
        # Should we use standard M and N lattices?
        if name is None:
            if dual_name is not None:
                raise ValueError("you can name the dual lattice only if you "
                                 "also name the original one!")
            name = "N"
            dual_name = "M"
        if latex_name is None:
            latex_name = name
        # Now name and latex_name are set
        # The default for latex_dual_name depends on whether dual_name was
        # given or constructed, so we determine it before dual_name
        if latex_dual_name is None:
            latex_dual_name = (dual_name if dual_name is not None
                                         else latex_name + "^*")
        if dual_name is None:
            dual_name = name + "*"
        return (rank, name, dual_name, latex_name, latex_dual_name)

    def create_object(self, version, key):
        r"""
        Create the toric lattice described by ``key``.

        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        .. WARNING::

            You probably should not use this function directly.

        TESTS::

            sage: key = ToricLattice.create_key(3)
            sage: ToricLattice.create_object(1, key)
            3-d lattice N
        """
        return ToricLatticeClass(*key)


ToricLattice = ToricLatticeFactory("ToricLattice")


# Possible TODO's:
# - implement a better construction() method, which still will prohibit
#   operations mixing lattices by conversion to ZZ^n
# - coersion between lattices and sublattices
# - maybe __call__ is not the right place to prohibit conversion between
#   lattices (we need it now so that morphisms behave nicely)
class ToricLatticeClass(FreeModule_ambient_pid):
    r"""
    Create a toric lattice.

    See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

    .. WARNING::

        There should be only one toric lattice with the given rank and
        associated names. Using this class directly to create toric lattices
        may lead to unexpected results. Please, use :class:`ToricLattice
        <ToricLatticeFactory>` to create toric lattices.

    TESTS::

        sage: from sage.geometry.toric_lattice import (
        ...         ToricLatticeClass)
        sage: N = ToricLatticeClass(3, "N", "M", "N", "M")
        sage: N
        3-d lattice N
        sage: TestSuite(N).run()
    """

    def __init__(self, rank, name, dual_name, latex_name, latex_dual_name):
        r"""
        See :class:`ToricLattice <ToricLatticeFactory>` for documentation.

        TESTS::

            sage: from sage.geometry.toric_lattice import (
            ...         ToricLatticeClass)
            sage: ToricLatticeClass(3, "N", "M", "N", "M")
            3-d lattice N
        """
        super(ToricLatticeClass, self).__init__(ZZ, rank)
        self._name = name
        self._dual_name = dual_name
        self._latex_name = latex_name
        self._latex_dual_name = latex_dual_name
        # This is how other free modules work now, but it seems that things
        # should be a bit different in the new coersion model
        self._element_class = ToricLatticeElement

    # It is not recommended to override __call__ in Parent-derived objects
    # since it may interfere with the coersion model. We do it here to allow
    # N(1,2,3) to be interpreted as N([1,2,3]). We also prohibit N(m) where
    # m is an element of another lattice. Otherwise morphisms will care only
    # about dimension of lattices.
    def __call__(self, *args, **kwds):
        r"""
        Construct a new element of ``self``.

        INPUT:

        - anything that can be interpreted as coordinates, except for elements
          of other lattices.

        OUTPUT:

        - :class:`ToricLatticeElement`.

        TESTS::

            sage: N = ToricLattice(3)
            sage: N([1,2,3])
            N(1, 2, 3)

        The point of overriding this function was to allow writing the above
        command as::

            sage: N(1,2,3)
            N(1, 2, 3)

        And to prohibit conversion between different lattices::

            sage: M = N.dual()
            sage: M(N(1,2,3))
            Traceback (most recent call last):
            ...
            TypeError: N(1, 2, 3) cannot be converted to 3-d lattice M!

        We also test that the special treatment of zero still works::

            sage: N(0)
            N(0, 0, 0)
        """
        supercall = super(ToricLatticeClass, self).__call__
        if args == (0, ):
            # Special treatment for N(0) to return (0,...,0)
            return supercall(*args, **kwds)
        try:
            coordinates = map(ZZ, args)
        except TypeError:
            # Prohibit conversion of elements of other lattices
            if (is_ToricLatticeElement(args[0])
                and args[0].parent() is not self):
                raise TypeError("%s cannot be converted to %s!"
                                % (args[0], self))
            # "Standard call"
            return supercall(*args, **kwds)
        # Coordinates were given without packing them into a list or a tuple
        return supercall(coordinates, **kwds)

    def __cmp__(self, right):
        r"""
        Compare ``self`` and ``right``.

        INPUT:

        - ``right`` -- anything.

        OUTPUT:

        - 0 if ``right`` is a toric lattice of the same dimension as ``self``
          and their associated names are the same, 1 or -1 otherwise.

        TESTS::

            sage: N3 = ToricLattice(3)
            sage: N4 = ToricLattice(4)
            sage: M3 = N3.dual()
            sage: cmp(N3, N4)
            -1
            sage: cmp(N3, M3)
            1
            sage: abs( cmp(N3, 3) )
            1
            sage: cmp(N3, ToricLattice(3))
            0
        """
        if self is right:
            return 0
        c = cmp(type(self), type(right))
        if c:
            return c
        c = cmp(self.rank(), right.rank())
        if c:
            return c
        # If lattices are the same as ZZ-modules, compare associated names
        return cmp([self._name, self._dual_name,
                    self._latex_name, self._latex_dual_name],
                   [right._name, right._dual_name,
                    right._latex_name, right._latex_dual_name])

    def __contains__(self, point):
        r"""
        Check if ``point`` is an element of ``self``.

        INPUT:

        - ``point`` -- anything.

        OUTPUT:

        - ``True`` if ``point`` is an element of ``self``, ``False``
          otherwise.

        TESTS::

            sage: N = ToricLattice(3)
            sage: M = N.dual()
            sage: L = ToricLattice(3, "L")
            sage: 1 in N
            False
            sage: (1,0) in N
            False
            sage: (1,0,0) in N
            True
            sage: N(1,0,0) in N
            True
            sage: M(1,0,0) in N
            False
            sage: L(1,0,0) in N
            False
            sage: (1/2,0,0) in N
            False
            sage: (2/2,0,0) in N
            True
        """
        try:
            self(point)
        except TypeError:
            return False
        return True

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.dual()._latex_()
            'L^*'
        """
        return self._latex_name

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: L = ToricLattice(3, "L")
            sage: L.dual()._repr_()
            '3-d lattice L*'
        """
        return "%d-d lattice %s" % (self.dimension(), self._name)

    # We need to override this function, otherwise e.g. the sum of elements of
    # different lattices of the same dimension will live in ZZ^n.
    def construction(self):
        r"""
        Return the functorial construction of ``self``.

        OUTPUT:

        - ``None``, we do not think of toric lattices as constructed from
          simpler objects since we do not want to perform arithmetic involving
          different lattices.

        TESTS::

            sage: print ToricLattice(3).construction()
            None
        """
        return None

    def dual(self):
        r"""
        Return the lattice dual to ``self``.

        OUTPUT:

        - :class:`toric lattice <ToricLatticeFactory>`.

        EXAMPLES::

            sage: N = ToricLattice(3)
            sage: N
            3-d lattice N
            sage: M = N.dual()
            sage: M
            3-d lattice M
            sage: M.dual() is N
            True

        Elements of dual lattices can act on each other::

            sage: n = N(1,2,3)
            sage: m = M(4,5,6)
            sage: n * m
            32
            sage: m * n
            32
        """
        if "_dual" not in self.__dict__:
            self._dual = ToricLattice(self.rank(), self._dual_name,
                          self._name, self._latex_dual_name, self._latex_name)
        return self._dual
