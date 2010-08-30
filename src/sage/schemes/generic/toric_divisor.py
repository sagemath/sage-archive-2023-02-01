r"""
Toric divisors and divisor classes

Let `X` be a :class:`toric variety
<sage.schemes.generic.toric_variety.ToricVariety_field>` corresponding to a
:class:`rational polyhedral fan <sage.geometry.fan.RationalPolyhedralFan>`
`\Sigma`. A :class:`toric divisor <ToricDivisor_generic>` `D` is a T-Weil
divisor over a given coefficient ring (usually `\ZZ` or `\QQ`), i.e. a formal
linear combination of torus-invariant subvarieties of `X` of codimension one.
In homogeneous coordinates `[z_0:\cdots:z_k]`, these are the subvarieties
`\{z_i=0\}`. Note that there is a finite number of such subvarieties, one for
each ray of `\Sigma`. We generally identify

  * Toric divisor `D`,

  * Sheaf `\mathcal{O}(D)` (if `D` is Cartier, it is a line bundle),

  * Support function `\phi_D` (if `D` is `\QQ`-Cartier, it is a function
    linear on each cone of `\Sigma`).

EXAMPLES:

We start with an illustration of basic divisor arithmetic::

    sage: dP6 = toric_varieties.dP6()
    sage: Dx,Du,Dy,Dv,Dz,Dw = dP6.toric_divisor_group().gens()
    sage: Dx
    V(x)
    sage: -Dx
    -V(x)
    sage: 2*Dx
    2*V(x)
    sage: Dx*2
    2*V(x)
    sage: (1/2)*Dx + Dy/3 - Dz
    1/2*V(x) + 1/3*V(y) - V(z)
    sage: Dx.parent()
    Group of toric ZZ-Weil divisors
    on 2-d CPR-Fano toric variety covered by 6 affine patches
    sage: (Dx/2).parent()
    Group of toric QQ-Weil divisors
    on 2-d CPR-Fano toric variety covered by 6 affine patches

Now we create a more complicated variety to demonstrate divisors of different
types::

    sage: F = Fan(cones=[(0,1,2,3), (0,1,4)],
    ...       rays=[(1,1,1), (1,-1,1), (1,-1,-1), (1,1,-1), (0,0,1)])
    sage: X = ToricVariety(F)
    sage: QQ_Cartier = X.divisor([2,2,1,1,1])
    sage: Cartier = 2 * QQ_Cartier
    sage: Weil = X.divisor([1,1,1,0,0])
    sage: QQ_Weil = 1/2 * Weil
    sage: [QQ_Weil.is_QQ_Weil(),
    ...    QQ_Weil.is_Weil(),
    ...    QQ_Weil.is_QQ_Cartier(),
    ...    QQ_Weil.is_Cartier()]
    [True, False, False, False]
    sage: [Weil.is_QQ_Weil(),
    ...    Weil.is_Weil(),
    ...    Weil.is_QQ_Cartier(),
    ...    Weil.is_Cartier()]
    [True, True, False, False]
    sage: [QQ_Cartier.is_QQ_Weil(),
    ...    QQ_Cartier.is_Weil(),
    ...    QQ_Cartier.is_QQ_Cartier(),
    ...    QQ_Cartier.is_Cartier()]
    [True, True, True, False]
    sage: [Cartier.is_QQ_Weil(),
    ...    Cartier.is_Weil(),
    ...    Cartier.is_QQ_Cartier(),
    ...    Cartier.is_Cartier()]
    [True, True, True, True]

The toric (`\QQ`-Weil) divisors modulo linear equivalence form the
toric divisor class group :class:`ToricRationalDivisorClassGroup`.  If
the toric variety `X` is smooth, this equals the **Picard group**
`\mathop{\mathrm{Pic}}(X)`::

    sage: Cl = dP6.divisor_class_group(); Cl
    The toric QQ-divisor class group of a 2-d CPR-Fano toric variety covered by 6 affine patches
    sage: Dx.divisor_class()
    Divisor class [1,0,0,0]
    sage: Dx.divisor_class() in Cl
    True

The (rational) divisor class group is where the Kahler cone lives::

    sage: Kc = dP6.Kaehler_cone(); Kc
    4-d cone in 4-d lattice
    sage: Kc.rays()
    (Divisor class [1,1,0,0], Divisor class [1,1,1,0], Divisor class [0,1,1,1],
    Divisor class [0,0,1,1], Divisor class [0,1,1,0])
    sage: Kc.ray(1).lift()
    V(x) + V(u) + V(y)
"""


########################################################################
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################


from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import is_Vector
from sage.geometry.cone import is_Cone
from sage.geometry.toric_lattice_element import is_ToricLatticeElement
from sage.misc.all import latex
from sage.modules.all import vector
from sage.modules.free_module import FreeModule_ambient_field
from sage.modules.vector_rational_dense import Vector_rational_dense
from sage.rings.all import QQ, ZZ
from sage.matrix.constructor import identity_matrix
from sage.schemes.generic.divisor import Divisor_generic
from sage.schemes.generic.divisor_group import DivisorGroup_generic
from sage.schemes.generic.toric_variety import CohomologyRing, is_ToricVariety


#********************************************************
class ToricDivisorGroup(DivisorGroup_generic):
    r"""
    The group of (`\QQ`-T-Weil) divisors on a toric variety.

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: P2.toric_divisor_group()
        Group of toric ZZ-Weil divisors
        on 2-d CPR-Fano toric variety covered by 3 affine patches
    """

    def __init__(self, toric_variety, base_ring):
        r"""
        Construct an instance of :class:`ToricDivisorGroup`.

        INPUT:

        - ``toric_variety`` -- a
          :class:`toric variety
          <sage.schemes.generic.toric_variety.ToricVariety_field>``;

        - ``base_ring`` -- the coefficient ring of this divisor group,
          usually `\ZZ` (default) or `\QQ`.

        Implementation note: :meth:`__classcall__` sets the default
        value for ``base_ring``.

        OUTPUT:

        Divisor group of the toric variety.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: from sage.schemes.generic.toric_divisor import ToricDivisorGroup
            sage: ToricDivisorGroup(P2, base_ring=ZZ)
            Group of toric ZZ-Weil divisors
            on 2-d CPR-Fano toric variety covered by 3 affine patches

        Note that :class:`UniqueRepresentation` correctly distinguishes the
        parent classes even if the schemes are the same::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivisorGroup(P2,ZZ) is ToricDivisorGroup(P2,ZZ)
            False
        """
        assert is_ToricVariety(toric_variety), str(toric_variety)+' is not a toric variety!'
        super(ToricDivisorGroup, self).__init__(toric_variety, base_ring)

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        OUTPUT:

        - string.

        TESTS::

            sage: toric_varieties.P2().toric_divisor_group()._latex_()
            '\\mathrm{Div_T}\\left(\\mathbb{P}_{\\Delta^{2}}, \\Bold{Z}\\right)'
        """
        return (r"\mathrm{Div_T}\left(%s, %s\right)"
                % (latex(self.scheme()), latex(self.base_ring())))

    def _repr_(self):
        """
        Return a string representation of the toric divisor group.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: toric_varieties.P2().toric_divisor_group()._repr_()
            'Group of toric ZZ-Weil divisors
            on 2-d CPR-Fano toric variety covered by 3 affine patches'
        """
        ring = self.base_ring()
        if ring == ZZ:
            base_ring_str = 'ZZ'
        elif ring == QQ:
            base_ring_str = 'QQ'
        else:
            base_ring_str = '('+str(ring)+')'
        return 'Group of toric '+base_ring_str+'-Weil divisors on '+str(self.scheme())

    def ngens(self):
        r"""
        Return the number of generators.

        OUTPUT:

        The number of generators of ``self``, which equals the number of
        rays in the fan of the toric variety.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: TDiv = P2.toric_divisor_group()
            sage: TDiv.ngens()
            3
        """
        return self.scheme().fan().nrays()

    def gens(self):
        r"""
        Return the generators of the divisor group.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: TDiv = P2.toric_divisor_group()
            sage: TDiv.gens()
            (V(x), V(y), V(z))
        """
        # Note: self._gens is originally incorrectly set by the parent class
        if self._gens is None:
            one = self.base_ring().one()
            self._gens = tuple(ToricDivisor_generic([(one, c)], self)
                               for c in self.scheme().gens())
        return self._gens

    def gen(self,i):
        r"""
        Return the ``i``-th generator of the divisor group.

        INPUT:

        - ``i`` -- integer.

        OUTPUT:

        The divisor `z_i=0`, where `z_i` is the `i`-th homogeneous
        coordinate.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: TDiv = P2.toric_divisor_group()
            sage: TDiv.gen(2)
            V(z)
        """
        return self.gens()[i]

    def _element_constructor_(self, x, check=True, reduce=True):
        r"""
        Construct a :class:`ToricDivisor_generic`

        INPUT:

        - ``x`` -- something defining a toric divisor, see
          :func:`ToricDivisor`.

        - ``check``, ``reduce`` -- boolean. See
          :meth:`ToricDivisor_generic.__init__`.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: TDiv = P2.toric_divisor_group()
            sage: TDiv._element_constructor_([ (1,P2.gen(2)) ])
            V(z)
        """
        if is_ToricDivisor(x):
            if x.parent() is self:
                return x
            else:
                x = x._data
        return ToricDivisor(self.scheme(), x, self.base_ring(), check, reduce)

    def base_extend(self, R):
        """
        Extend the scalars of ``self`` to ``R``.

        INPUT:

        - ``R`` -- ring.

        OUTPUT:

        - toric divisor group.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: DivZZ = P2.toric_divisor_group()
            sage: DivQQ = P2.toric_divisor_group(base_ring=QQ)
            sage: DivZZ.base_extend(QQ) is DivQQ
            True
        """
        # This check prevents extension to cohomology rings via coercion
        if isinstance(R,CohomologyRing):
            raise TypeError, 'Coefficient ring cannot be a cohomology ring.'
        if self.base_ring().has_coerce_map_from(R):
            return self
        elif R.has_coerce_map_from(self.base_ring()):
            return ToricDivisorGroup(self.scheme(), base_ring=R)
        else:
            raise ValueError("the base of %s cannot be extended to %s!"
                             % ( self, R))

#********************************************************
def is_ToricDivisor(x):
    r"""
    Test whether ``x`` is a toric divisor.

    INPUT:

    - ``x`` -- anything.

    OUTPUT:

    - ``True`` if ``x`` is an instance of :class:`ToricDivisor_generic` and
      ``False`` otherwise.

    EXAMPLES::

        sage: from sage.schemes.generic.toric_divisor import is_ToricDivisor
        sage: is_ToricDivisor(1)
        False
        sage: P2 = toric_varieties.P2()
        sage: D = P2.divisor(0); D
        V(x)
        sage: is_ToricDivisor(D)
        True
    """
    return isinstance(x, ToricDivisor_generic)


#********************************************************
def ToricDivisor(toric_variety, arg=None, ring=None, check=False, reduce=False):
    r"""
    Construct a divisor of ``toric_variety``.

    INPUT:

    - ``toric_variety`` -- a :class:`toric variety
      <sage.schemes.generic.toric_variety.ToricVariety_field>`;

    - ``arg`` -- one of the following description of the toric divisor to be
      constructed:

      * ``None`` (the trivial divisor);

      * integer (divisor corresponding to the ``arg``-th ray of the fan, order
        is the same as in ``toric_variety.fan().rays()``, which is also the
        same as for the one-cones ``toric_variety.fan(dim=1)``);

      * monomial in the homogeneous coordinates;

      * one-dimensional cone of the fan of ``toric_variety`` or a lattice
        point generating such a cone;

      * sequence of rational numbers, specifying multiplicities for each of
        the toric divisors.

    - ``ring`` -- usually either `\ZZ` or `\QQ`. The base ring of the
      divisor group. If ``ring`` is not specified, a coefficient ring
      suitable for ``arg`` is derived.

    - for ``check`` and ``reduce`` see :meth:`ToricDivisor_generic`.

    OUTPUT:

    - A :class:`sage.schemes.generic.toric_divisor.ToricDivisor_generic`

    EXAMPLES::

        sage: from sage.schemes.generic.toric_divisor import ToricDivisor
        sage: dP6 = toric_varieties.dP6()
        sage: ToricDivisor(dP6, [(1,dP6.gen(2)), (1,dP6.gen(1))])
        V(u) + V(y)
        sage: ToricDivisor(dP6, 1) + ToricDivisor(dP6, 2)
        V(u) + V(y)
        sage: ToricDivisor(dP6, (0,1,1,0,0,0), ring=QQ)
        V(u) + V(y)
        sage: dP6.inject_variables()
        Defining x, u, y, v, z, w
        sage: ToricDivisor(dP6, u+y)
        Traceback (most recent call last):
        ...
        ValueError: u + y is not a monomial!
        sage: ToricDivisor(dP6, u*y)
        V(u) + V(y)
        sage: ToricDivisor(dP6, dP6.fan(dim=1)[0] )
        V(x)
        sage: N = dP6.fan().lattice()
        sage: ToricDivisor(dP6, N(1,1) )
        V(w)
    """
    assert is_ToricVariety(toric_variety)
    # Zero divisor
    if arg is None:
        arg = []
    # Divisor by lattice point (corresponding to a ray)
    if is_ToricLatticeElement(arg):
        if arg not in toric_variety.fan().lattice():
            raise ValueError("%s is not in the ambient lattice of %s!"
                             % (arg, toric_variety.fan()))
        arg = toric_variety.fan().cone_containing(arg)
    # Divisor by a one-cone
    if is_Cone(arg):
        if arg.dim() != 1:
            raise ValueError("Only 1-dimensional cones of the toric variety "
                             "define divisors.")
        arg = arg.ambient_ray_indices()[0]
    # Divisor by a ray index
    if arg in ZZ:
        arg = [(1, toric_variety.gen(arg))]
    # Divisor by monomial
    if arg in toric_variety.coordinate_ring():
        if len(list(arg)) != 1:
            raise ValueError("%s is not a monomial!" % arg)
        arg = arg.exponents()[0]
    # By now either we have converted arg to a list, or it is something else
    # which should be convertible to a list
    if not isinstance(arg, list):
        try:
            arg = list(arg)
        except TypeError:
            raise TypeError("%s does not define a divisor!" % arg)
    # Now we have a list of multiplicities or pairs multiplicity-coordinate,
    # if the coefficient ring was not given, try to use the most common ones.
    if ring is None:
        try:
            return ToricDivisor(toric_variety, arg, ring=ZZ,
                                check=check, reduce=reduce)
        except TypeError:
            pass
        try:
            return ToricDivisor(toric_variety, arg, ring=QQ,
                                check=check, reduce=reduce)
        except TypeError:
            raise TypeError("cannot deduce coefficient ring for %s!" % arg)
    # The ring was given, if we are still here
    n_rays = toric_variety.fan().nrays()
    if len(arg) == n_rays and arg[0] in ring:
        arg = zip(arg, toric_variety.gens())
    # Now we MUST have a LIST of multiplicity-coordinate pairs
    try:
        for i, pair in enumerate(arg):
            pair = tuple(pair)
            if len(pair) != 2 or pair[1] not in toric_variety.gens():
                raise TypeError
            arg[i] = (ring(pair[0]), pair[1])
    except TypeError:
        raise TypeError("%s does not define a divisor!" % arg)
    TDiv = ToricDivisorGroup(toric_variety, ring)
    return ToricDivisor_generic(arg, TDiv, check, reduce)


#********************************************************
class ToricDivisor_generic(Divisor_generic):
    """
    Construct a :class:`(toric Weil) divisor <ToricDivisor_generic>` on the
    given toric variety.

    INPUT:

    - ``v`` -- a list of tuples (multiplicity, coordinate).

    - ``parent`` -- :class:`ToricDivisorGroup`. The parent divisor group.

    - ``check`` -- boolean. Type-check the entries of ``v``, see
      :meth:`sage.schemes.generic.divisor_group.DivisorGroup_generic.__init__`.

    - ``reduce`` -- boolean. Combine coefficients in ``v``, see
      :meth:`sage.schemes.generic.divisor_group.DivisorGroup_generic.__init__`.

    .. WARNING::

        Do not construct :class:`ToricDivisor_generic` objects manually.
        Instead, use either the function :func:`ToricDivisor` or the method
        :meth:`~sage.schemes.generic.toric_variety.ToricVariety_field.divisor`
        of toric varieties.

    EXAMPLES::

        sage: dP6 = toric_varieties.dP6()
        sage: ray = dP6.fan().ray(0)
        sage: ray
        N(0, 1)
        sage: D = dP6.divisor(ray); D
        V(x)
        sage: D.parent()
        Group of toric ZZ-Weil divisors
        on 2-d CPR-Fano toric variety covered by 6 affine patches
    """

    def __init__(self, v, parent, check=True, reduce=True):
        """
        See :class:`ToricDivisor_generic` for documentation.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: from sage.schemes.generic.toric_divisor import ToricDivisor_generic
            sage: TDiv = dP6.toric_divisor_group()
            sage: ToricDivisor_generic([], TDiv)
            0
            sage: ToricDivisor_generic([(2,dP6.gen(1))], TDiv)
            2*V(u)
        """
        self._variety = parent.scheme()
        self._fan = self._variety.fan()
        super(ToricDivisor_generic,self).__init__(v, parent, check, reduce)

    def _vector_(self, ring=None):
        r"""
        Return a vector representation.

        INPUT:

        - ``ring`` -- a ring (usually `\ZZ` or `\QQ`) for the
          coefficients to live in). This is an optional argument, by
          default a suitable ring is chosen automatically.

        OUTPUT:

        A vector whose ``self.scheme().fan().nrays()`` components are
        the coefficients of the divisor.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: D = dP6.divisor((0,1,1,0,0,0)); D
            V(u) + V(y)
            sage: D._vector_()
            (0, 1, 1, 0, 0, 0)
            sage: vector(D)        # syntactic sugar
            (0, 1, 1, 0, 0, 0)
            sage: type( vector(D) )
            <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>
            sage: D_QQ = dP6.divisor((0,1,1,0,0,0), base_ring=QQ);
            sage: vector(D_QQ)
            (0, 1, 1, 0, 0, 0)
            sage: type( vector(D_QQ) )
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>

        The vector representation is a suitable input for :func:`ToricDivisor` ::

            sage: dP6.divisor(vector(D)) == D
            True
        """
        if ring is None:
            ring = self.base_ring()
        v = vector(ring, [0]*self._fan.nrays() )
        for coeff, variable in self:
            v[ self._variety.gens().index(variable) ] += coeff
        return v

    def coefficient(self, x):
        r"""
        Return the coefficient of ``x``.

        INPUT:

        - ``x`` -- one of the homogeneous coordinates, either given by
          the variable or its index.

        OUTPUT:

        The coefficient of ``x``.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: D = P2.divisor((11,12,13)); D
            11*V(x) + 12*V(y) + 13*V(z)
            sage: D.coefficient(1)
            12
            sage: P2.inject_variables()
            Defining x, y, z
            sage: D.coefficient(y)
            12
        """
        try:
            index = ZZ(x)
            variable = self._variety.gen(index)
        except TypeError:
            variable = x
        total = self.base_ring().zero()
        for coeff, var in self:
            if var == variable:
                total += coeff
        return total

    def function_value(self, point):
        r"""
        Return the value of the support function at ``point``.

        Let `X` be the ambient toric variety of ``self``, `\Sigma` the fan
        associated to `X`, and `N` the ambient lattice of `\Sigma`.

        INPUT:

        - ``point`` -- either an integer, interpreted as the index of a ray of
          `\Sigma`, or a point of the lattice `N`.

        OUTPUT:

        - an interger or a rational number.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: D = P2.divisor([11,22,44])   # total degree 77
            sage: D.function_value(0)
            11
            sage: N = P2.fan().lattice()
            sage: D.function_value( N(1,1) )
            33
            sage: D.function_value( P2.fan().ray(0) )
            11
        """
        if not self.is_QQ_Cartier():
            raise ValueError("support functions are associated to QQ-Cartier "
                             "divisors only, %s is not QQ-Cartier!" % self)
        try:
            index = ZZ(point)
            return self.coefficient(index)
        except TypeError:
            pass

        assert point in self._fan.lattice(), 'The point '+str(point)+' is not in the N-lattice.'
        cone = self._fan.cone_containing(point)
        return point * self.m(cone)

    def m(self, cone):
        r"""
        Return `m_\sigma` representing `\phi_D` on ``cone``.

        Let `X` be the ambient toric variety of this divisor `D` associated to
        the fan `\Sigma` in lattice `N`. Let `M` be the lattice dual to `N`.
        Given the cone `\sigma =\langle v_1, \dots, v_k \rangle` in `\Sigma`,
        this method searches for a vector `m_\sigma \in M_\QQ` such that
        `\phi_D(v_i) = \langle m_\sigma, v_i \rangle` for all `i=1, \dots, k`,
        where `\phi_D` is the support function of `D`.

        INPUT:

        - ``cone`` -- A cone in the fan of the toric variety.

        OUTPUT:

        - If possible, a point of lattice `M`.

        - If the dual vector cannot be chosen integral, a rational vector is
          returned.

        - If there is no such vector (i.e. ``self`` is not even a
          `\QQ`-Cartier divisor), a ``ValueError`` is raised.

        EXAMPLES::

            sage: F = Fan(cones=[(0,1,2,3), (0,1,4)],
            ...       rays=[(1,1,1), (1,-1,1), (1,-1,-1), (1,1,-1), (0,0,1)])
            sage: X = ToricVariety(F)
            sage: square_cone = X.fan().cone_containing(0,1,2,3)
            sage: triangle_cone = X.fan().cone_containing(0,1,4)
            sage: ray = X.fan().cone_containing(0)
            sage: QQ_Cartier = X.divisor([2,2,1,1,1])
            sage: QQ_Cartier.m(ray)
            M(0, 2, 0)
            sage: QQ_Cartier.m(square_cone)
            (3/2, 0, 1/2)
            sage: QQ_Cartier.m(triangle_cone)
            M(1, 0, 1)
            sage: Weil = X.divisor([1,1,1,0,0])
            sage: Weil.m(square_cone)
            Traceback (most recent call last):
            ...
            ValueError: V(z0) + V(z1) + V(z2) is not QQ-Cartier,
            cannot choose a dual vector on 3-d cone
            of Rational polyhedral fan in 3-d lattice N!
            sage: Weil.m(triangle_cone)
            M(1, 0, 0)
        """
        try:
            return self._m[cone]
        except AttributeError:
            self._m = {}
        except KeyError:
            pass

        M = self._fan.lattice().dual()
        if cone.is_trivial():
            m = M(0)
            self._m[cone] = m
            return m

        b = vector(self.coefficient(i) for i in cone.ambient_ray_indices())
        A = cone.ray_matrix()
        try:
            if cone.dim() == self._variety.dimension():
                # either unique solution or ValueError (if not QQ-Cartier)
                m = A.solve_left(b)  # A m = b
            else:
                # under-determined system; try to find integral solution
                D,U,V = A.smith_form()   # D = U*A*V
                bV = b*V
                m = D.solve_left(bV) * U
        except ValueError:
            raise ValueError("%s is not QQ-Cartier, cannot choose a dual "
                             "vector on %s!" % (self, cone))

        try:
            m = M(m)
        except TypeError:  # not integral
            pass
        self._m[cone] = m
        return m

    def is_Weil(self):
        """
        Return whether the divisor is a Weil-divisor.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: D = P2.divisor([1,2,3])
            sage: D.is_Weil()
            True
            sage: (D/2).is_Weil()
            False
        """
        if self.base_ring() == ZZ:
            return True
        try:
            vector(ZZ, vector(self))
            return True
        except TypeError:
            return False

    def is_QQ_Weil(self):
        r"""
        Return whether the divisor is a `\QQ`-Weil-divisor.

        .. NOTE::

            This function returns always ``True`` since
            :class:`ToricDivisor <ToricDivisor_generic>` can only
            describe `\QQ`-Weil divisors.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: D = P2.divisor([1,2,3])
            sage: D.is_QQ_Weil()
            True
            sage: (D/2).is_QQ_Weil()
            True
        """
        return True

    def is_Cartier(self):
        r"""
        Return whether the divisor is a Cartier-divisor.

        .. NOTE::

            The sheaf `\mathcal{O}(D)` associated to the given divisor
            `D` is a line bundle if and only if the divisor is
            Cartier.

        EXAMPLES::

            sage: X = toric_varieties.P4_11169()
            sage: D = X.divisor(3)
            sage: D.is_Cartier()
            False
            sage: D.is_QQ_Cartier()
            True
        """
        try:
            return self._is_Cartier
        except AttributeError:
            pass

        self._is_Cartier = self.is_QQ_Cartier()
        if self._is_Cartier:
            M = self._fan.lattice().dual()
            self._is_Cartier = all(self.m(c) in M for c in self._fan)
        return self._is_Cartier

    def is_QQ_Cartier(self):
        """
        Return whether the divisor is a `\QQ`-Cartier divisor.

        A `\QQ`-Cartier divisor is a divisor such that some multiple
        of it is Cartier.

        EXAMPLES::

            sage: X = toric_varieties.P4_11169()
            sage: D = X.divisor(3)
            sage: D.is_QQ_Cartier()
            True

            sage: X = toric_varieties.Cube_face_fan()
            sage: D = X.divisor(3)
            sage: D.is_QQ_Cartier()
            False
        """
        try:
            return self._is_QQ_Cartier
        except AttributeError:
            pass

        try:
            [self.m(c) for c in self._fan]
            self._is_QQ_Cartier = True
        except ValueError:
            self._is_QQ_Cartier = False
        return self._is_QQ_Cartier

    def is_integral(self):
        r"""
        Return whether the support function is integral on the rays.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: DZZ = P2.toric_divisor_group(base_ring=ZZ).gen(0); DZZ
            V(x)
            sage: DQQ = P2.toric_divisor_group(base_ring=QQ).gen(0); DQQ
            V(x)
            sage: DZZ.is_integral()
            True
            sage: DQQ.is_integral()
            True
        """
        return all( self.function_value(r) in ZZ
                    for r in self._fan.rays() )

    def move_away_from(self, cone):
        """
        Move the divisor away from the orbit closure of ``cone``.

        INPUT:

        - A ``cone`` of the fan of the toric variety.

        OUTPUT:

        A (rationally equivalent) divisor that is moved off the
        orbit closure of the given cone.

        .. NOTE::

            A divisor that is Weil but not Cartier might be impossible
            to move away. In this case, a ``ValueError`` is raised.

        EXAMPLES::

            sage: F = Fan(cones=[(0,1,2,3), (0,1,4)],
            ...       rays=[(1,1,1), (1,-1,1), (1,-1,-1), (1,1,-1), (0,0,1)])
            sage: X = ToricVariety(F)
            sage: square_cone = X.fan().cone_containing(0,1,2,3)
            sage: triangle_cone = X.fan().cone_containing(0,1,4)
            sage: line_cone = square_cone.intersection(triangle_cone)
            sage: Cartier = X.divisor([2,2,1,1,1])
            sage: Cartier
            2*V(z0) + 2*V(z1) + V(z2) + V(z3) + V(z4)
            sage: Cartier.move_away_from(line_cone)
            -V(z2) - V(z3) + V(z4)
            sage: QQ_Weil = X.divisor([1,0,1,1,0])
            sage: QQ_Weil.move_away_from(line_cone)
            V(z2)
        """
        m = self.m(cone)
        if m in self._fan.lattice():
            ring = self._ring
        else:
            ring = m.base_ring()
        divisor = vector(self)
        values = [ divisor[i] - m*self._fan.ray(i)
                   for i in range(len(divisor)) ]
        return ToricDivisor(self._variety, values, ring=ring)

    def cohomology_class(self):
        r"""
        Return the degree-2 cohomology class associated to the divisor.

        OUTPUT:

        Returns the corresponding cohomology class as an instance of
        :class:`~sage.schemes.generic.toric_variety.CohomologyClass`.
        The cohomology class is the first Chern class of the
        associated line bundle `\mathcal{O}(D)`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: D = dP6.divisor(dP6.fan().ray(0) )
            sage: D.cohomology_class()
            [y + v - w]
        """
        divisor = vector(self)
        return sum([ divisor[i] * cone.cohomology_class()
                     for [i,cone] in enumerate(self._fan(dim=1)) ])

    def Chern_character(self):
        r"""
        Return the Chern character of the sheaf `\mathcal{O}(D)`
        defined by the divisor `D`.

        You can also use a shortcut :meth:`ch`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: N = dP6.fan().lattice()
            sage: D3 = dP6.divisor(dP6.fan().cone_containing( N(0,1)   ))
            sage: D5 = dP6.divisor(dP6.fan().cone_containing( N(-1,-1) ))
            sage: D6 = dP6.divisor(dP6.fan().cone_containing( N(0,-1)  ))
            sage: D = -D3 + 2*D5 - D6
            sage: D.Chern_character()
            [5*w^2 + y - 2*v + w + 1]
            sage: dP6.integrate( D.ch() * dP6.Td() )
            -4
        """
        return self.cohomology_class().exp()

    ch = Chern_character

    def divisor_class(self):
        r"""
        Return the linear equivalence class of the divisor.

        OUTPUT:

        Returns the class of the divisor in `\mathop{Cl}(X)
        \otimes_\ZZ \QQ` as an instance of
        :class:`ToricRationalDivisorClass`.

        EXAMPLE:

            sage: dP6 = toric_varieties.dP6()
            sage: D = dP6.divisor(0)
            sage: D.divisor_class()
            Divisor class [1,0,0,0]
        """
        if '_divisor_class' not in self.__dict__:
            Cl = self._variety.divisor_class_group()(self)
            self._divisor_class = Cl
        return self._divisor_class


#********************************************************
class ToricRationalDivisorClassGroup(FreeModule_ambient_field, UniqueRepresentation):
    r"""
    The rational divisor class group of a toric variety.

    The **T-Weil divisor class group** `\mathop{Cl}(X)` of a toric
    variety `X` is a finitely generated abelian group and can contain
    torsion. Its rank equals the number of rays in the fan of `X`
    minus the dimension of `X`.

    The **rational divisor class group** is `\mathop{Cl}(X)
    \otimes_\ZZ \QQ` and never includes torsion. If `X` is *smooth*,
    this equals the **Picard group** `\mathop{\mathrm{Pic}}(X)`, whose
    elements are the isomorphism classes of line bundles on `X`. The
    group law (which we write as addition) is the tensor product of
    the line bundles. The Picard group of a toric variety is always
    torsion-free.

    .. WARNING::

        You don't need to instantiate this class yourself. Use
        :meth:`sage.schemes.generic.ToricVariety_field.divisor_class_group`
        if you need the divisor class group. Or you can obtain it as
        the parent of any divisor class constructed, for example, with
        :meth:`ToricDivisor_generic.divisor_class()`

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: P2.divisor_class_group()
        The toric QQ-divisor class group of a 2-d CPR-Fano
        toric variety covered by 3 affine patches
        sage: D = P2.divisor(1); D
        V(y)
        sage: Dclass = D.divisor_class(); Dclass
        Divisor class [1]
        sage: Dclass.lift()
        V(x)
        sage: Dclass.parent()
        The toric QQ-divisor class group of a 2-d CPR-Fano
        toric variety covered by 3 affine patches
    """

    def __init__(self, toric_variety):
        r"""
        Construct the toric rational divisor class group.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: from sage.schemes.generic.toric_divisor import ToricRationalDivisorClassGroup
            sage: ToricRationalDivisorClassGroup(P2)
            The toric QQ-divisor class group of a 2-d CPR-Fano
            toric variety covered by 3 affine patches
        """
        self._variety = toric_variety
        nrays = toric_variety.fan().nrays()
        rk = nrays - toric_variety.fan().lattice_dim()
        super(ToricRationalDivisorClassGroup,self).__init__(base_field=QQ, dimension=rk, sparse=False)
        gale = toric_variety.fan().gale_transform()
        self._projection_matrix = gale.matrix_from_columns(range(0,nrays))
        self._lift_matrix = self._projection_matrix.solve_right( identity_matrix(rk) )

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: from sage.schemes.generic.toric_divisor import ToricRationalDivisorClassGroup
            sage: ToricRationalDivisorClassGroup(P2)._repr_()
            'The toric QQ-divisor class group of a 2-d CPR-Fano toric variety covered by 3 affine patches'
        """
        return 'The toric QQ-divisor class group of a '+self._variety._repr_()

    def _latex_(self):
        r"""
        Return a LaTeX representation of ``self``.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: from sage.schemes.generic.toric_divisor import ToricRationalDivisorClassGroup
            sage: ToricRationalDivisorClassGroup(P2)._latex_()
            '\\mathop{Cl}_{\\QQ}\\left(\\mathbb{P}_{\\Delta^{2}}\\right)'
        """
        return '\\mathop{Cl}_{\\QQ}\\left('+self._variety._latex_()+'\\right)'

    def gen(self, i):
        r"""
        EXAMPLES::

            sage: toric_varieties.dP8().divisor_class_group().gen(0)
            Divisor class [1,0]
        """
        return self._element_constructor_( super(ToricRationalDivisorClassGroup,self).gen(i) )

    def gens(self):
        r"""
        EXAMPLES::

            sage: toric_varieties.dP8().divisor_class_group().gens()
            (Divisor class [1,0], Divisor class [0,1])
        """
        return tuple( self._element_constructor_(x)
                      for x in  super(ToricRationalDivisorClassGroup,self).gens() )

    def _element_constructor_(self, x):
        r"""
        Construct a :class:`ToricRationalDivisorClass`

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: Cl = dP6.divisor_class_group()
            sage: D = dP6.divisor(2)
            sage: Cl._element_constructor_(D)
            Divisor class [0,0,1,0]
            sage: Cl(D)
            Divisor class [0,0,1,0]
        """
        if is_ToricDivisor(x):
            coefficients = list(self._projection_matrix * vector(x))
            return ToricRationalDivisorClass(self, coefficients)
        if is_Vector(x):
            return ToricRationalDivisorClass(self, list(x))
        return ToricRationalDivisorClass(self, x)

    # parent does not conform to the new-style coercion model
    __call__ = _element_constructor_


#********************************************************
class ToricRationalDivisorClass(Vector_rational_dense):
    r"""
    An element of :class:`ToricRationalDivisorClassGroup`

    EXAMPLES::

        sage: dP6 = toric_varieties.dP6()
        sage: Cl = dP6.divisor_class_group()
        sage: D = dP6.divisor(2)
        sage: Cl(D)
        Divisor class [0,0,1,0]
    """
    def __init__(self, parent, x, coerce=True, copy=True):
        r"""
        Construct a :class:`ToricRationalDivisorClass`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: Cl = dP6.divisor_class_group()
            sage: from sage.schemes.generic.toric_divisor import ToricRationalDivisorClass
            sage: ToricRationalDivisorClass(Cl, [1,2,3,4])
            Divisor class [1,2,3,4]
        """
        assert isinstance(parent, ToricRationalDivisorClassGroup)
        super(ToricRationalDivisorClass,self).__init__(parent, x, coerce, copy)
        self.set_immutable()

    def _vector_(self, ring=QQ):
        r"""
        Return a vector representation.

        INPUT:

        - ``ring`` -- the coefficient ring.

        OUTPUT:

        A vector with coefficients in ``ring``.

        EXAMPLES::

            sage: dP8 = toric_varieties.dP8()
            sage: c = sum( dP8.divisor_class_group().gens() )
            sage: c._vector_()
            (1, 1)
            sage: vector(c)    # indirect test
            (1, 1)
        """
        if ring==None:
            ring = QQ
        return vector(ring, list(self))

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: toric_varieties.P2().divisor(0).divisor_class()._repr_()
            'Divisor class [1]'
        """
        return 'Divisor class [' + ','.join(map(str,list(self))) + ']'

    def _add_(self, right):
        r"""
        Return the sum.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._add_(c[1])
            Divisor class [1,1]
            sage: c[0] + c[1]      # indirect test
            Divisor class [1,1]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._add_(right) )

    def _sub_(self, right):
        r"""
        Return the difference.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._sub_(c[1])
            Divisor class [1,-1]
            sage: c[0] - c[1]      # indirect test
            Divisor class [1,-1]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._sub_(right) )

    def _dot_product_(self, right):
        r"""
        Return ``ValueError``.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._dot_product_(c[1])
            Traceback (most recent call last):
            ...
            ValueError: Cannot multiply two divisor classes.
            sage: c[0] * c[1]      # indirect test
            Traceback (most recent call last):
            ...
            ValueError: Cannot multiply two divisor classes.
        """
        raise ValueError, 'Cannot multiply two divisor classes.'

    def _pairwise_product_(self, right):
        r"""
        Return the pairwise product.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._pairwise_product_(c[1])
            Divisor class [0,0]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._pairwise_product_(right) )

    def _rmul_(self, left):
        r"""
        Return left scalar multiple.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._rmul_(2)
            Divisor class [2,0]
            sage: 2 * c[0]
            Divisor class [2,0]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._rmul_(left) )

    def _lmul_(self, right):
        r"""
        Return right scalar multiple.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._lmul_(2)
            Divisor class [2,0]
            sage: c[0] * 2
            Divisor class [2,0]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._lmul_(right) )

    def _neg_(self):
        r"""
        Return the negative.

        EXAMPLES::

            sage: c = toric_varieties.dP8().divisor_class_group().gens()
            sage: c[0]._neg_()
            Divisor class [-1,0]
            sage: - c[0]      # indirect test
            Divisor class [-1,0]
        """
        return self.parent()( super(ToricRationalDivisorClass,self)._neg_() )

    def lift(self):
        r"""
        Return a divisor in the given divisor class.

        OUTPUT:

        An instance of :class:`ToricDivisor` in the same divisor
        class.

        EXAMPLES::

            sage: X = toric_varieties.Cube_nonpolyhedral()
            sage: D = X.divisor([0,1,2,3,4,5,6,7]); D
            V(z1) + 2*V(z2) + 3*V(z3) + 4*V(z4) + 5*V(z5) + 6*V(z6) + 7*V(z7)
            sage: D.divisor_class()
            Divisor class [29,6,8,10,0]
            sage: Dequiv = D.divisor_class().lift(); Dequiv
            29/2*V(z0) + 6*V(z1) + 8*V(z2) + 10*V(z3)
            sage: Dequiv == D
            False
            sage: Dequiv.divisor_class() == D.divisor_class()
            True
        """
        if '_lift' not in self.__dict__:
            lift = self.parent()._lift_matrix * vector(self)
            self._lift = self.parent()._variety.divisor(lift)
        return self._lift


