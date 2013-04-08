"""
Free Groups

Free groups and finitely presented groups are implemented as a wrapper
over the corresponing GAP objects.

A free group can be created by giving the number of generators, or their names.
It is also possible to create indexed generators::

    sage: G.<x,y,z> = FreeGroup();  G
    Free Group on generators {x, y, z}
    sage: FreeGroup(3)
    Free Group on generators {x0, x1, x2}
    sage: FreeGroup('a,b,c')
    Free Group on generators {a, b, c}
    sage: FreeGroup(3,'t')
    Free Group on generators {t0, t1, t2}

The elements can be created by operating with the generators, or by passing a list
with the indices of the letters to the group:

EXAMPLES::

    sage: G.<a,b,c> = FreeGroup()
    sage: a*b*c*a
    a*b*c*a
    sage: G([1,2,3,1])
    a*b*c*a
    sage: a * b / c * b^2
    a*b*c^-1*b^2
    sage: G([1,1,2,-1,-3,2])
    a^2*b*a^-1*c^-1*b

You can use call syntax to replace the generators with a set of
arbitrary ring elements::

    sage: g =  a * b / c * b^2
    sage: g(1,2,3)
    8/3
    sage: M1 = identity_matrix(2)
    sage: M2 = matrix([[1,1],[0,1]])
    sage: M3 = matrix([[0,1],[1,0]])
    sage: g([M1, M2, M3])
    [1 3]
    [1 2]

AUTHORS:

- Miguel Angel Marco Buzunariz
- Volker Braun
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.groups.group import Group
from sage.groups.libgap_wrapper import ParentLibGAP, ElementLibGAP
from sage.structure.unique_representation import UniqueRepresentation
from sage.libs.gap.libgap import libgap
from sage.libs.gap.element import GapElement
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.misc.cachefunc import cached_method


def is_FreeGroup(x):
    """
    Test whether ``x`` is a :class:`FreeGroup_class`.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    Boolean.

    EXAMPLES::

        sage: from sage.groups.free_group import is_FreeGroup
        sage: is_FreeGroup('a string')
        False
        sage: is_FreeGroup(FreeGroup(0))
        True
    """
    return isinstance(x, FreeGroup_class)


class FreeGroupElement(ElementLibGAP):
    """
    A wrapper of GAP's Free Group elements.

    INPUT:

    - ``x`` -- something that determines the group element. Either a
      :class:`~sage.libs.gap.element.GapElement` or the Tietze list
      (see :meth:`Tietze`) of the group element.

    - ``parent`` -- the parent :class:`FreeGroup`.

    EXAMPLES::

        sage: G = FreeGroup('a, b')
        sage: x = G([1, 2, -1, -2])
        sage: x
        a*b*a^-1*b^-1
        sage: y = G([2, 2, 2, 1, -2, -2, -2])
        sage: y
        b^3*a*b^-3
        sage: x*y
        a*b*a^-1*b^2*a*b^-3
        sage: y*x
        b^3*a*b^-3*a*b*a^-1*b^-1
        sage: x^(-1)
        b*a*b^-1*a^-1
        sage: x == x*y*y^(-1)
        True
    """

    def __init__(self, x, parent):
        """
        The Python constructor.

        See :class:`FreeGroupElement` for details.

        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: x = G([1, 2, -1, -1])
            sage: x # indirect doctest
            a*b*a^-2
            sage: y = G([2, 2, 2, 1, -2, -2, -1])
            sage: y # indirect doctest
            b^3*a*b^-2*a^-1

            sage: TestSuite(G).run()
            sage: TestSuite(x).run()
        """
        if not isinstance(x, GapElement):
            try:
                l = x.Tietze()
            except AttributeError:
                l = list(x)
            if len(l)>0:
                if min(l) < -parent.ngens() or parent.ngens() < max(l):
                    raise ValueError('generators not in the group')
            if 0 in l:
                raise ValueError('zero does not denote a generator')
            i=0
            while i<len(l)-1:
                if l[i]==-l[i+1]:
                    l.pop(i)
                    l.pop(i)
                    if i>0:
                        i=i-1
                else:
                    i=i+1
            AbstractWordTietzeWord = libgap.eval('AbstractWordTietzeWord')
            x = AbstractWordTietzeWord(l, parent._gap_gens())
        ElementLibGAP.__init__(self, x, parent)

    def __reduce__(self):
        """
        Implement pickling.

        TESTS::

            sage: F.<a,b> = FreeGroup()
            sage: a.__reduce__()
            (Free Group on generators {a, b}, ((1,),))
            sage: (a*b*a^-1).__reduce__()
            (Free Group on generators {a, b}, ((1, 2, -1),))
        """
        return (self.parent(), (self.Tietze(),))

    @cached_method
    def Tietze(self):
        """
        Return the Tietze list of the element.

        The Tietze list of a word is a list of integers that represent
        the letters in the word.  A positive integer `i` represents
        the letter corresponding to the `i`-th generator of the group.
        Negative integers represent the inverses of generators.

        OUTPUT:

        A tuple of integers.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: a.Tietze()
            (1,)
            sage: x = a^2 * b^(-3) * a^(-2)
            sage: x.Tietze()
            (1, 1, -2, -2, -2, -1, -1)

        TESTS::

            sage: type(a.Tietze())
            <type 'tuple'>
            sage: type(a.Tietze()[0])
            <type 'sage.rings.integer.Integer'>
        """
        tl = self.gap().TietzeWordAbstractWord()
        return tuple(tl.sage())

    def fox_derivative(self, gen):
        """
        Return the Fox derivative of self with respect to a given generator.

        INPUT:

        - ``gen`` : the generator with respect to which the derivative will be computed.

        OUTPUT:

        An element of the group algebra with integer coefficients.

        EXAMPLES::

            sage: G = FreeGroup(5)
            sage: G.inject_variables()
            Defining x0, x1, x2, x3, x4
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x0)
            -B[x0^-1] + B[x0^-1*x1] - B[x0^-1*x1*x0*x2*x0^-1]
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x1)
            B[x0^-1]
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x2)
            B[x0^-1*x1*x0]
            sage: (~x0*x1*x0*x2*~x0).fox_derivative(x3)
            0
        """
        if not gen in self.parent().generators():
            raise ValueError("Fox derivative can only be computed with respect to generators of the group")
        l = list(self.Tietze())
        R = self.parent().algebra(IntegerRing())
        i = gen.Tietze()[0]
        a = R.zero()
        while len(l)>0:
            b = l.pop(-1)
            if b == i:
                p = R(self.parent()(l))
                a += p
            elif b == -i:
                p = R(self.parent()(l+[b]))
                a -= p
        return a

    @cached_method
    def syllables(self):
        r"""
        Return the syllables of the word.

        Consider a free group element `g = x_1^{n_1} x_2^{n_2} \cdots
        x_k^{n_k}`. The uniquely-determined subwords `x_i^{e_i}`
        consisting only of powers of a single generator are called the
        syllables of `g`.

        OUTPUT:

        The tuple of syllables. Each syllable is given as a pair
        `(x_i, e_i)` consisting of a generator and a non-zero integer.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: w = a^2 * b^-1 * a^3
            sage: w.syllables()
            ((a, 2), (b, -1), (a, 3))
        """
        g = self.gap().UnderlyingElement()
        k = g.NumberSyllables().sage()
        gen = self.parent().gen
        exponent_syllable  = libgap.eval('ExponentSyllable')
        generator_syllable = libgap.eval('GeneratorSyllable')
        result = []
        gen = self.parent().gen
        for i in range(k):
            exponent  = exponent_syllable(g, i+1).sage()
            generator = gen(generator_syllable(g, i+1).sage() - 1)
            result.append( (generator, exponent) )
        return tuple(result)

    def __call__(self, *values):
        """
        Replace the generators of the free group.

        INPUT:

        - ``*values`` -- a sequence of values, or a
          list/tuple/iterable of the same length as the number of
          generators.

        OUTPUT:

        The product of ``values`` in the order and with exponents
        specified by ``self``.

        EXAMPLES::

            sage: G.<a,b> = FreeGroup()
            sage: w = a^2 * b^-1 * a^3
            sage: w(1, 2)
            1/2
            sage: w(2, 1)
            32
            sage: w.subs(b=1, a=2)   # indirect doctest
            32

        TESTS::

            sage: w([1, 2])
            1/2
            sage: w((1, 2))
            1/2
            sage: w(i+1 for i in range(2))
            1/2
        """
        if len(values) == 1:
            try:
                values = list(values[0])
            except TypeError:
                pass
        G = self.parent()
        if len(values) != G.ngens():
            raise ValueError('number of values has to match the number of generators')
        replace = dict(zip(G.gens(), values))
        from sage.misc.all import prod
        return prod( replace[gen] ** power for gen, power in self.syllables() )


def FreeGroup(n=None, names='x'):
    """
    Construct a Free Group

    INPUT:

    - ``n`` -- integer or ``None`` (default). The nnumber of
      generators. If not specified the ``names`` are counted.

    - ``names`` -- string or list/tuple/iterable of strings (default:
      ``'x'``). The generator names or name prefix.

    EXAMPLES::

        sage: G.<a,b> = FreeGroup();  G
        Free Group on generators {a, b}
        sage: H = FreeGroup('a, b')
        sage: G is H
        True
        sage: FreeGroup(0)
        Free Group on generators {}

    The entry can be either a string with the names of the generators,
    or the number of generators and the prefix of the names to be
    given. The default prefix is ``'x'`` ::

        sage: FreeGroup(3)
        Free Group on generators {x0, x1, x2}
        sage: FreeGroup(3, 'g')
        Free Group on generators {g0, g1, g2}
        sage: FreeGroup()
        Free Group on generators {x}

    TESTS::

        sage: G1 = FreeGroup(2, 'a,b')
        sage: G2 = FreeGroup('a,b')
        sage: G3.<a,b> = FreeGroup()
        sage: G1 is G2, G2 is G3
        (True, True)
    """
    # Support Freegroup('a,b') syntax
    if n is not None:
        try:
            n = Integer(n)
        except TypeError:
            names = n
            n = None
    # derive n from counting names
    if n is None:
        if isinstance(names, basestring):
            n = len(names.split(','))
        else:
            names = list(names)
            n = len(names)
    from sage.structure.parent import normalize_names
    names = tuple(normalize_names(n, names))
    return FreeGroup_class(names)


def wrap_FreeGroup(libgap_free_group):
    """
    Wrap a LibGAP free group.

    This function changes the comparison method of
    ``libgap_free_group`` to comparison by Python ``id``. If you want
    to put the LibGAP free group into a container (set, dict) then you
    should understand the implications of
    :meth:`~sage.libs.gap.element.GapElement._set_compare_by_id`. To
    be safe, it is recommended that you just work with the resulting
    Sage :class:`FreeGroup_class`.

    INPUT:

    - ``libgap_free_group`` -- a LibGAP free group.

    OUTPUT:

    A Sage :class:`FreeGroup_class`.

    EXAMPLES:

    First construct a LibGAP free group::

        sage: F = libgap.FreeGroup(['a', 'b'])
        sage: type(F)
        <type 'sage.libs.gap.element.GapElement'>

    Now wrap it::

        sage: from sage.groups.free_group import wrap_FreeGroup
        sage: wrap_FreeGroup(F)
        Free Group on generators {a, b}

    TESTS:

    Check that we can do it twice (see :trac:`12339`) ::

        sage: G = libgap.FreeGroup(['a', 'b'])
        sage: wrap_FreeGroup(G)
        Free Group on generators {a, b}
    """
    assert libgap_free_group.IsFreeGroup()
    libgap_free_group._set_compare_by_id()
    names = tuple( str(g) for g in libgap_free_group.GeneratorsOfGroup() )
    return FreeGroup_class(names, libgap_free_group)


class FreeGroup_class(UniqueRepresentation, Group, ParentLibGAP):
    """
    A class that wraps GAP's FreeGroup

    See :func:`FreeGroup` for details.

    TESTS::

        sage: G = FreeGroup('a, b')
        sage: TestSuite(G).run()
    """
    Element = FreeGroupElement

    def __init__(self, generator_names, libgap_free_group=None):
        """
        Python constructor.

        INPUT:

        - ``generator_names`` -- a tuple of strings. The names of the
          generators.

        - ``libgap_free_group`` -- a LibGAP free group or ``None``
          (default). The LibGAP free group to wrap. If ``None``, a
          suitable group will be constructed.

        TESTS::

            sage: G.<a,b> = FreeGroup() # indirect doctest
            sage: G
            Free Group on generators {a, b}
            sage: G.variable_names()
            ('a', 'b')
        """
        n = len(generator_names)
        self._assign_names(generator_names)
        if libgap_free_group is None:
            libgap_free_group = libgap.FreeGroup(generator_names)
        ParentLibGAP.__init__(self, libgap_free_group)
        Group.__init__(self)

    def _repr_(self):
        """
        TESTS::

            sage: G = FreeGroup('a, b')
            sage: G # indirect doctest
            Free Group on generators {a, b}
            sage: G._repr_()
            'Free Group on generators {a, b}'
        """
        return 'Free Group on generators {'+ ', '.join(self.variable_names()) + '}'

    def rank(self):
        """
        Return the number of generators of self.

        Alias for :meth:`ngens`.

        OUTPUT:

        Integer.

        EXAMPLES::

            sage: G = FreeGroup('a, b');  G
            Free Group on generators {a, b}
            sage: G.rank()
            2
            sage: H = FreeGroup(3, 'x')
            sage: H
            Free Group on generators {x0, x1, x2}
            sage: H.rank()
            3
        """
        return self.ngens()

    def _gap_init_(self):
        """
        Return the string used to construct the object in gap.

        EXAMPLES::

            sage: G = FreeGroup(3)
            sage: G._gap_init_()
            'FreeGroup(["x0", "x1", "x2"])'
        """
        gap_names = [ '"' + s + '"' for s in self.variable_names() ]
        gen_str = ', '.join(gap_names)
        return 'FreeGroup(['+gen_str+'])'

    def _element_constructor_(self, *args, **kwds):
        """
        TESTS::

            sage: G.<a,b> = FreeGroup()
            sage: G([1, 2, 1]) # indirect doctest
            a*b*a
            sage: G([1, 2, -2, 1, 1, -2]) # indirect doctest
            a^3*b^-1

            sage: G( G._gap_gens()[0] )
            a
            sage: type(_)
            <class 'sage.groups.free_group.FreeGroup_class_with_category.element_class'>
        """
        if len(args)!=1:
            return self.element_class(*args, parent=self, **kwds)
        x = args[0]
        if x==1:
            return self.one()
        try:
            P = x.parent()
        except AttributeError:
            return self.element_class(x, parent=self, **kwds)
        if hasattr(P, '_freegroup_'):
            if P.FreeGroup() is self:
                return self.element_class(x.Tietze(), parent=self, **kwds)
        return self.element_class(x, parent=self, **kwds)

    def abelian_invariants(self):
        r"""
        Return the Abelian invariants of self.

        The Abelian invariants are given by a list of integers `i_1 \dots i_j`, such that the
        abelianization of the group is isomorphic to

        .. math::

            \ZZ / (i_1) \times \dots \times \ZZ / (i_j)

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: F.abelian_invariants()
            (0, 0)
        """
        inv = self.gap().AbelianInvariants()
        return tuple(inv.sage())

    def quotient(self, relations):
        """
        Return the quotient of self by the normal subgroup generated
        by the given elements.

        This quotient is a finitely presented groups with the same
        generators as ``self``, and relations given by the elements of
        ``relations``.

        INPUT:

        - ``relations`` -- A list/tuple/iterable with the elements of
          the free group.

        OUTPUT:

        A finitely presented group, with generators corresponding to
        the generators of the free group, and relations corresponding
        to the elements in ``relations``.

        EXAMPLES::

            sage: F.<a,b> = FreeGroup()
            sage: F.quotient([a*b^2*a, b^3])
            Finitely presented group < a, b | a*b^2*a, b^3 >

        Division is shorthand for :meth:`quotient` ::

            sage: F /  [a*b^2*a, b^3]
            Finitely presented group < a, b | a*b^2*a, b^3 >
        """
        from sage.groups.finitely_presented import FinitelyPresentedGroup
        return FinitelyPresentedGroup(self, tuple(relations))

    __div__ = quotient
