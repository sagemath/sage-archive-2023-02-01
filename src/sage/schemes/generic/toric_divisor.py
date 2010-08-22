r"""
Toric Divisors

A ``ToricDivisor`` is, in general, a `\QQ` T-Weil divisor, that is, a
formal `\QQ`-linear combination of torus-invariant subvarieties. In
homogeneous coordinates `[z_0:\cdots:z_k]`, these are the subvarieties
`\{z_i=0\}`. Note that there is a finite number of such subvarieties,
one for each ray of the fan. We generally identify

  * Toric divisors `D`

  * Piecewise linear support functions on the fan

  * The sheaf `\mathcal{O}(D)`

If the divisor is actually a Cartier divisor (see
:meth:`ToricDivisor_generic.is_Cartier`), then the support function is linear on
each cone and `\mathcal{O}(D)` is a line bundle.

EXAMPLES::

    sage: dP6 = toric_varieties.dP6()
    sage: Dx,Du,Dy,Dv,Dz,Dw = dP6.toric_divisor_group().gens()
    sage: Dx
    Divisor x
    sage: +Dx
    Divisor x
    sage: -Dx
    Divisor -x
    sage: 2*Dx
    Divisor 2*x
    sage: Dx*2
    Divisor 2*x
    sage: Dx + Dy
    Divisor y + x
    sage: Dx - Dy
    Divisor -y + x
    sage: (1/2)*Dx
    Divisor 1/2*x
    sage: Dx / 2
    Divisor 1/2*x
    sage: Dx.parent()
    Group of toric ZZ-Weil divisors on 2-d CPR-Fano toric variety covered by 6 affine patches
    sage: (Dx/2).parent()
    Group of toric QQ-Weil divisors on 2-d CPR-Fano toric variety covered by 6 affine patches
"""


########################################################################
#       Copyright (C) 2010 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.modules.all import vector
from sage.rings.all import ZZ, QQ, Integer, is_Ideal
from sage.matrix.constructor import matrix
from sage.misc.all import flatten, prod, subsets

from sage.geometry.toric_lattice import ToricLattice, is_ToricLatticeElement
from sage.geometry.cone import Cone, is_Cone
from sage.geometry.fan import Fan
from sage.schemes.generic.divisor_group import DivisorGroup_generic
from sage.schemes.generic.divisor import Divisor_generic
from sage.schemes.generic.toric_variety import ToricVariety, is_ToricVariety, is_CohomologyClass, CohomologyRing
from sage.schemes.generic.toric_variety_library import toric_varieties



#********************************************************
class ToricDivisorGroup(DivisorGroup_generic):
    r"""
    The group of (`\QQ` T-Weil) divisors on a toric variety.

    .. WARNING::

        You should not create instances by hand, but either use
        :func:`ToricDivisorGroup` or
        :meth:`sage.schemes.generic.toric_variety.ToricVariety_field.divisor_group`.

    EXAMPLES::

        sage: P2 = toric_varieties.P2()
        sage: P2.toric_divisor_group()
        Group of toric ZZ-Weil divisors on 2-d CPR-Fano toric variety covered by 3 affine patches
    """

    def __init__(self, toric_variety, base_ring):
        r"""
        Construct an instance of :class:`ToricDivisorGroup`

        INPUT:

        - ``toric_variety`` -- a
          :class:`sage.schemes.generic.toric_variety.ToricVariety_generic``.

        - ``base_ring`` -- the coefficient ring of the toric variety,
          usually `\ZZ` or `\QQ`.

        Implementation note: :meth:`__classcall__` sets the default
        value for ``base_ring``.

        OUTPUT:

        Divisor group of the toric variety.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: from sage.schemes.generic.toric_divisor import ToricDivisorGroup
            sage: ToricDivisorGroup(P2, base_ring=ZZ)
            Group of toric ZZ-Weil divisors on 2-d CPR-Fano toric variety covered by 3 affine patches

        Note that :class:`UniqueRepresentation` correctly
        distinguishes the parent classes even if the schemes are the
        same::

            sage: from sage.schemes.generic.divisor_group import DivisorGroup
            sage: DivisorGroup(P2,ZZ) is ToricDivisorGroup(P2,ZZ)
            False
        """
        assert is_ToricVariety(toric_variety), str(toric_variety)+' is not a toric variety!'
        super(ToricDivisorGroup,self).__init__(toric_variety, base_ring)

    def _repr_(self):
        """
        Return a string representation of the toric divisor group.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: toric_varieties.P2().toric_divisor_group()._repr_()
            'Group of toric ZZ-Weil divisors on 2-d CPR-Fano toric variety covered by 3 affine patches'
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

        The number of generators equals, which equals the number of
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
            [Divisor x, Divisor y, Divisor z]
        """
        # Note: self._gens is already defined by parent class
        try:
            return self._toric_gens
        except AttributeError:
            pass
        self._toric_gens = [ ToricDivisor_generic([(self.base_ring().one(),c)],self)
                             for c in self.scheme().gens() ]
        return self._toric_gens

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
            Divisor z
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
            Divisor z
        """
        if isinstance(x, ToricDivisor_generic):
            P = x.parent()
            if P is self:
                return x
            elif P == self:
                return ToricDivisor_generic(x._data, check=False, reduce=False, parent=self)
            else:
                x = x._data

        if isinstance(x, list):
            all( len(xi)==2 for xi in x )
            all( xi[0] in self.base_ring() for xi in x )
            all( xi[1] in self.scheme().gens() for xi in x )
            return ToricDivisor_generic(x, check=check, reduce=reduce, parent=self)
        elif x in self.scheme().gens():
            return ToricDivisor_generic([(self.base_ring()(1), x)],
                                        check=False, reduce=False, parent=self)
        elif isinstance(x, ZZ) and x == 0:
            return ToricDivisor_generic([], check=False, reduce=False, parent=self)

        return ToricDivisor(self.scheme(), x)


    def base_extend(self, R):
        """
        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: DivZZ = P2.toric_divisor_group()
            sage: DivQQ = P2.toric_divisor_group(base_ring=QQ)
            sage: DivZZ.base_extend(QQ) is DivQQ
            True
        """
        if isinstance(R,CohomologyRing):
            raise TypeError, 'Coefficient ring cannot be a cohomology ring.'
        if self.base_ring().has_coerce_map_from(R):
            return self
        elif R.has_coerce_map_from(self.base_ring()):
            return ToricDivisorGroup(self.scheme(), base_ring=R)


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
        Divisor x
        sage: is_ToricDivisor(D)
        True
    """
    return isinstance(x, ToricDivisor_generic)


#********************************************************
def ToricDivisor(toric_variety, arg=None, ring=None, check=False, reduce=False):
    r"""
    Construct a divisor of a toric variety.

    INPUT:

    - ``toric_variety`` must be a :class:`ToricVariety <ToricVariety_field>`

    - if ``arg`` equals ``None``, the trivial divisor is constructed.

    - if ``arg`` is an integer, the divisor corresponding to the
      ``arg``-th ray of the fan is constructed. The divisors are in
      the same order as ``toric_variety.fan().rays()``, which is also
      the same order as the one-cones ``toric_variety.fan(dim=1)``

    - if ``arg`` is a monomial in the homogeneous coordinates, the
      corresponding toric divisor is constructed.

    - if ``arg`` is a one-dimensional cone of the fan of the toric
      variety or a N-lattice point generating such a ray, the
      corresponding divisor is constructed.

    - if ``arg`` is a sequence of rational numbers, one for each ray
      of the fan, the corresponding linear combination of toric
      divisors is returned.

    - if ``arg`` is a valid argument for the constructor of
      :class:`ToricDivisor_generic`, then it is simply passed on.

    - ``ring`` -- usually either `\ZZ` or `\QQ`. The base ring of the
      divisor group. If ``ring`` is not specified, a coefficient ring
      suitable for ``arg`` is derived.

    - for ``check`` and ``reduce`` see
      :meth:`ToricDivisor_generic.__init__`.

    OUTPUT:

    - A :class:`sage.schemes.generic.toric_divisor.ToricDivisor_generic`

    EXAMPLES::

        sage: dP6 = toric_varieties.dP6()
        sage: ToricDivisor(dP6, [(1,dP6.gen(2)), (1,dP6.gen(1))])
        Divisor y + u
        sage: ToricDivisor(dP6, 1) + ToricDivisor(dP6, 2)
        Divisor y + u
        sage: ToricDivisor(dP6, (0,1,1,0,0,0), ring=QQ)
        Divisor u + y
        sage: dP6.inject_variables()
        Defining x, u, y, v, z, w
        sage: ToricDivisor(dP6, u+y)
        Traceback (most recent call last):
        ...
        ValueError: The polynomial u + y must consist of a single monomial.
        sage: ToricDivisor(dP6, u*y)
        Divisor u + y
        sage: ToricDivisor(dP6, dP6.fan(dim=1)[0] )
        Divisor x
        sage: N = dP6.fan().lattice()
        sage: ToricDivisor(dP6, N(1,1) )
        Divisor w
    """
    assert is_ToricVariety(toric_variety);
    n_rays = toric_variety.fan().nrays()

    if ring==None:
        R = ZZ
    else:
        R = ring

    if arg==None:
        divisor = []
    elif arg in ZZ:
        cone_index = arg
        divisor = [ (R(1),toric_variety.gen(cone_index)) ]
    elif is_Cone(arg):
        if arg.dim() != 1:
            raise ValueError, 'Only 1-dimensional cones of the toric variety define divisors.'
        cone_index = arg.ambient_ray_indices()[0]
        divisor = [ (R(1),toric_variety.gen(cone_index)) ]
    elif is_ToricLatticeElement(arg):
        if arg not in toric_variety.fan().lattice():
            raise ValueError, 'The argument must be in the N-lattice of the toric variety.'
        cone = toric_variety.fan().cone_containing(arg)
        if cone.dim() != 1:
            raise ValueError, 'The point '+str(arg)+' does not lie on a ray of the fan.'
        cone_index = cone.ambient_ray_indices()[0]
        divisor = [ (R(1),toric_variety.gen(cone_index)) ]
    elif arg in toric_variety.coordinate_ring():
        if len(list(arg)) != 1:
            raise ValueError, \
                'The polynomial '+str(arg)+' must consist of a single monomial.'
        exponent = arg.exponents()[0]
        divisor = [ (R(exponent[i]),toric_variety.gen(i))
                    for i in range(0,n_rays) ]
    else:
        # we tried everything else, so assume arg is iterable
        divisor=list(arg)
        try:
            # is arg a valid argument for ToricDivisor_generic.__init__?
            all( len(x)==2 for x in divisor )
            all( x[1] in toric_variety.gens() for x in divisor )
        except TypeError:
            # if all fails, try to convert arg to a vector over the ring
            if ring!=None:  # the
                divisor = vector(ring, arg)
            else:
                # Try to deduce a suitable ring
                try:
                    divisor = vector(ZZ, arg)
                    R = ZZ
                except TypeError:
                    divisor = vector(QQ, arg)
                    R = QQ
            assert len(divisor) == n_rays
            divisor = [ (R(divisor[i]),toric_variety.gen(i))
                        for i in range(0,n_rays) ]

    TDiv = ToricDivisorGroup(toric_variety, R)
    return ToricDivisor_generic(divisor, TDiv, check, reduce)


#********************************************************
class ToricDivisor_generic(Divisor_generic):
    """
    A (`\QQ`-Weil) divisor of a toric variety.

    EXAMPLES::

        sage: dP6 = toric_varieties.dP6()
        sage: ray = dP6.fan().ray(0)
        sage: ray
        N(0, 1)
        sage: D = ToricDivisor(dP6, ray); D
        Divisor x
        sage: D.parent()
        Group of toric ZZ-Weil divisors on 2-d CPR-Fano toric variety covered by 6 affine patches
    """

    _desc = "A divisor"

    def __init__(self, v, parent, check=True, reduce=True):
        """
        Construct a :class:`ToricDivisor_generic` object on the given
        toric variety.

        INPUT:

        - ``v`` -- a list of tuples (multiplicity, equation).

        - ``parent`` -- :class:`ToricDivisorGroup`. The parent divisor
          group.

        - ``check`` -- boolean. Type check the entries of ``v``, see
          :meth:'sage.schemes.generic.divisor_group.DivisorGroup_generic.__init__'.

        - ``reduce`` -- boolean. Combine coefficients in ``v``, see
          :meth:'sage.schemes.generic.divisor_group.DivisorGroup_generic.__init__'.

        .. WARNING::

            Do not construct :class:`ToricDivisor_generic` objects
            manually. Instead, use either the function
            :func:`ToricDivisor` or the method
            :meth:`sage.schemes.generic.toric_variety.divisor` to
            construct toric divisor.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: from sage.schemes.generic.toric_divisor import ToricDivisorGroup, ToricDivisor_generic
            sage: TDiv = ToricDivisorGroup(dP6,ZZ)
            sage: ToricDivisor_generic([], TDiv)
            Divisor 0
            sage: ToricDivisor_generic([(2,dP6.gen(1))], TDiv)
            Divisor 2*u
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
            sage: D = ToricDivisor(dP6, (0,1,1,0,0,0)); D
            Divisor u + y
            sage: D._vector_()
            (0, 1, 1, 0, 0, 0)
            sage: vector(D)        # syntactic sugar
            (0, 1, 1, 0, 0, 0)
            sage: type( vector(D) )
            <type 'sage.modules.vector_integer_dense.Vector_integer_dense'>
            sage: D_QQ = ToricDivisor(dP6, (0,1,1,0,0,0), ring=QQ);
            sage: vector(D_QQ)
            (0, 1, 1, 0, 0, 0)
            sage: type( vector(D_QQ) )
            <type 'sage.modules.vector_rational_dense.Vector_rational_dense'>

        The vector representation is a suitable input for :func:`ToricDivisor` ::

            sage: ToricDivisor(dP6, vector(D)) == D
            True
        """
        if ring==None:
            ring = self.base_ring()
        v = vector(ring, [0]*self._fan.nrays() )
        for coeff, variable in self:
            v[ self._variety.gens().index(variable) ] = coeff
        return v;

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
            sage: D = ToricDivisor(P2, (11,12,13)); D
            Divisor 11*x + 12*y + 13*z
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
        for coeff, var in self:
            if var==variable:
                return coeff
        return self.base_ring().zero()

    def function_value(self, point):
        r"""
        Return the value of the piecewise-linear function at the point.

        INPUT:

        - Either the index of a ray; Note that indexing starts at 0 and
          ends at ``ToricVariety.fan().n_rays()-1``,

        - Or a N-lattice point.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: D = ToricDivisor(P2, [11,22,44] )   # total degree 77
            sage: D.function_value(0)
            11
            sage: N = P2.fan().lattice()
            sage: D.function_value( N(1,1) )
            33
            sage: D.function_value( P2.fan().ray(0) )
            11
        """
        try:
            index = ZZ(point)
            return self.coefficient(index)
        except TypeError:
            pass

        assert point in self._fan.lattice(), 'The point '+str(point)+' is not in the N-lattice.'
        cone = self._fan.cone_containing(point)
        try:
           m = self.m(cone)
        except ValueError: # Not Cartier divisor
           raise NotImplementedError
        return point * m

    def _repr_(self):
        r"""
        Return the string representation of a divisor.
        EXAMPLES::
           sage: dP6 = toric_varieties.dP6()
           sage: ToricDivisor(dP6, dP6.fan().ray(0) )._repr_()
           'Divisor x'
        """
        return 'Divisor '+super(ToricDivisor_generic,self)._repr_()

    def m(self, cone):
        r"""
        Return a choice of dual vector that represents the linear
        function on the cone.

        Given the cone `\sigma =\langle v_1, \dots \rangle` in a
        lattice `N` with dual lattice `M`, this method searches for a
        vector `m_\sigma \in M_\QQ` such that `f(v_i) = m_\sigma \cdot
        v_i`.

        INPUT:

        - ``cone`` -- A cone in the fan of the toric variety.

        OUTPUT:

        - If possible, a M-lattice point.

        - If the dual vector cannot be chosen integral, a QQ-vector is
          returned.

        - If there is no such vector, a ``ValueError`` is raised.

        EXAMPLES::

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: square_cone = X.fan().cone_containing(0,1,2,3)
            sage: triangle_cone = X.fan().cone_containing(0,1,4)
            sage: ray = X.fan().cone_containing(0)
            sage: QQ_Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1])
            sage: QQ_Cartier.m(ray)
            M(0, 2, 0)
            sage: QQ_Cartier.m(square_cone)
            (3/2, 0, 1/2)
            sage: QQ_Cartier.m(triangle_cone)
            M(1, 0, 1)
            sage: Weil = ToricDivisor(X,[1,1,1,0,0])
            sage: try:
            ...       Weil.m(square_cone)
            ... except ValueError:
            ...       print "ValueError raised"
            ...
            ValueError raised
            sage: Weil.m(triangle_cone)
            M(1, 0, 0)
        """
        try:
            cached_values = self._m
        except AttributeError:
            self._m = {}
            cached_values = self._m

        try:
            return cached_values[cone]
        except KeyError:
            pass

        M = self._fan.lattice().dual()

        if cone.is_trivial():
            m = M( [0]*self._variety.dimension() )
            cached_values[cone]=m
            return m

        b = vector([ self.coefficient(i) for i in cone.ambient_ray_indices() ])
        A = cone.ray_matrix()

        if cone.dim()==self._variety.dimension():
            # either unique solution or ValueError (if not QQ-Cartier)
            m=A.solve_left(b)  # A m = b
        else:
            # under-determined system; try to find integral solution
            D,U,V = A.smith_form()   # D = U*A*V
            bV = b*V
            m = D.solve_left(bV) * U

        try:
            m = M(m)
        except TypeError:  # not integral
            pass

        cached_values[cone]=m
        return m

    def is_Weil(self):
        """
        Return whether the divisor is a Weil-divisor.

        EXAMPLES::

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: QQ_Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1])
            sage: [ QQ_Cartier.is_Cartier(), QQ_Cartier.is_QQ_Cartier(), QQ_Cartier.is_Weil(), QQ_Cartier.is_QQ_Weil() ]
            [False, True, True, True]
            sage: Cartier = 2 * QQ_Cartier
            sage: [ Cartier.is_Cartier(), Cartier.is_QQ_Cartier(), Cartier.is_Weil(), Cartier.is_QQ_Weil() ]
            [True, True, True, True]
            sage: Weil = ToricDivisor(X,[1,1,1,0,0])
            sage: [ Weil.is_Cartier(), Weil.is_QQ_Cartier(), Weil.is_Weil(), Weil.is_QQ_Weil() ]
            [False, False, True, True]
            sage: QQ_Weil = 1/2 * Weil
            sage: [ QQ_Weil.is_Cartier(), QQ_Weil.is_QQ_Cartier(), QQ_Weil.is_Weil(), QQ_Weil.is_QQ_Weil() ]
            [False, False, False, True]
        """
        if self.base_ring()==ZZ:
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

            This function returns always true since ``ToricDivisor`` can only
            describe `\QQ`-Weil divisors.

        EXAMPLES::

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: QQ_Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1])
            sage: [ QQ_Cartier.is_Cartier(), QQ_Cartier.is_QQ_Cartier(), QQ_Cartier.is_Weil(), QQ_Cartier.is_QQ_Weil() ]
            [False, True, True, True]
            sage: Cartier = 2 * QQ_Cartier
            sage: [ Cartier.is_Cartier(), Cartier.is_QQ_Cartier(), Cartier.is_Weil(), Cartier.is_QQ_Weil() ]
            [True, True, True, True]
            sage: Weil = ToricDivisor(X,[1,1,1,0,0])
            sage: [ Weil.is_Cartier(), Weil.is_QQ_Cartier(), Weil.is_Weil(), Weil.is_QQ_Weil() ]
            [False, False, True, True]
            sage: QQ_Weil = 1/2 * Weil
            sage: [ QQ_Weil.is_Cartier(), QQ_Weil.is_QQ_Cartier(), QQ_Weil.is_Weil(), QQ_Weil.is_QQ_Weil() ]
            [False, False, False, True]
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

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: QQ_Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1])
            sage: [ QQ_Cartier.is_Cartier(), QQ_Cartier.is_QQ_Cartier(), QQ_Cartier.is_Weil(), QQ_Cartier.is_QQ_Weil() ]
            [False, True, True, True]
            sage: Cartier = 2 * QQ_Cartier
            sage: [ Cartier.is_Cartier(), Cartier.is_QQ_Cartier(), Cartier.is_Weil(), Cartier.is_QQ_Weil() ]
            [True, True, True, True]
            sage: Weil = ToricDivisor(X,[1,1,1,0,0])
            sage: [ Weil.is_Cartier(), Weil.is_QQ_Cartier(), Weil.is_Weil(), Weil.is_QQ_Weil() ]
            [False, False, True, True]
            sage: QQ_Weil = 1/2 * Weil
            sage: [ QQ_Weil.is_Cartier(), QQ_Weil.is_QQ_Cartier(), QQ_Weil.is_Weil(), QQ_Weil.is_QQ_Weil() ]
            [False, False, False, True]
        """
        try:
            return self._is_Cartier
        except AttributeError:
            pass

        self._is_Cartier = True
        for c in self._fan.generating_cones():
            try:
                m = self.m(c)
            except ValueError:
                self._is_Cartier = False
                break

            self._is_Cartier = ( m in self._fan.lattice().dual() )

            if not self._is_Cartier:
                break

        return self._is_Cartier

    def is_QQ_Cartier(self):
        """
        Return whether the divisor is a `\QQ`-Cartier divisor.

        A `\QQ`-Cartier divisor is a divisor such that some multiple
        of it is Cartier.

        EXAMPLES::

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: QQ_Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1])
            sage: [ QQ_Cartier.is_Cartier(), QQ_Cartier.is_QQ_Cartier(), QQ_Cartier.is_Weil(), QQ_Cartier.is_QQ_Weil() ]
            [False, True, True, True]
            sage: Cartier = 2 * QQ_Cartier
            sage: [ Cartier.is_Cartier(), Cartier.is_QQ_Cartier(), Cartier.is_Weil(), Cartier.is_QQ_Weil() ]
            [True, True, True, True]
            sage: Weil = ToricDivisor(X,[1,1,1,0,0])
            sage: [ Weil.is_Cartier(), Weil.is_QQ_Cartier(), Weil.is_Weil(), Weil.is_QQ_Weil() ]
            [False, False, True, True]
            sage: QQ_Weil = 1/2 * Weil
            sage: [ QQ_Weil.is_Cartier(), QQ_Weil.is_QQ_Cartier(), QQ_Weil.is_Weil(), QQ_Weil.is_QQ_Weil() ]
            [False, False, False, True]
        """
        try:
            return self._is_QQ_Cartier
        except AttributeError:
            pass

        try:
            if self._is_Cartier:
                return True
        except AttributeError:
            pass

        self._is_QQ_Cartier = True
        for c in self._fan.generating_cones():
            try:
                self.m(c)
            except ValueError:
                self._is_QQ_Cartier = False
                break
        return self._is_QQ_Cartier

    def is_integral(self):
        r"""
        Return whether the support function is integral on the rays.

        EXAMPLES::

            sage: P2 = toric_varieties.P2()
            sage: DZZ = P2.toric_divisor_group(base_ring=ZZ).gen(0); DZZ
            Divisor x
            sage: DQQ = P2.toric_divisor_group(base_ring=QQ).gen(0); DQQ
            Divisor x
            sage: DZZ.is_integral()
            True
            sage: DQQ.is_integral()
            True
        """
        return all( self.function_value(r) in ZZ
                    for r in self._fan.rays() )

    def move_away_from(self, cone):
        """
        Move the divisor away from the orbit closure of the cone.

        INPUT:

        - A ``cone`` of the fan of the toric variety.

        OUTPUT:

        A (rationally equivalent) divisor that is moved off the
        orbit closure of the given cone.

        .. NOTE::

            A divisor that is Weil but not Cartier might be impossible
            to move away. In this case, a ``ValueError`` is raised.

        EXAMPLES::

            sage: X = ToricVariety(Fan([[0,1,2,3],[0,1,4]], [(1, 1, 1), (1, -1, 1), (1, -1, -1), (1, 1, -1), (0, 0, 1)]))
            sage: square_cone = X.fan().cone_containing(0,1,2,3)
            sage: triangle_cone = X.fan().cone_containing(0,1,4)
            sage: line_cone = square_cone.intersection(triangle_cone)
            sage: Cartier = ToricDivisor(X, [1,1,1,1,0]) + ToricDivisor(X, [1,1,0,0,1]); Cartier
            Divisor z4 + z3 + z2 + 2*z1 + 2*z0
            sage: Cartier.move_away_from(line_cone)
            Divisor -z2 - z3 + z4
            sage: QQ_Weil = ToricDivisor(X, [1,0,1,1,0])
            sage: QQ_Weil.move_away_from(line_cone)
            Divisor z2
        """
        m = self.m(cone)
        if m in self._fan.lattice():
            ring = self._ring
        else:
            ring = m.base_ring()

        divisor = vector(self)
        values = [ divisor[i] - m*self._fan.ray(i)
                   for i in range(0,len(divisor)) ]

        return ToricDivisor(self._variety, values, ring=ring)

    def c1(self):
        r"""
        Return the degree-2 cohomology class associated to the divisor.

        OUTPUT:

        Returns the corresponding cohomology class as an instance of
        :class:`sage.schemes.generic.toric_variety.CohomologyClass`.
        The cohomology class is the first Chern class of the
        associated line bundle `\mathcal{O}(D)`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: D = ToricDivisor(dP6, dP6.fan().ray(0) )
            sage: D.c1()
            [y + v - w]
        """
        divisor = vector(self)
        return sum([ divisor[i] * cone.cohomology_class()
                     for [i,cone] in enumerate(self._fan(dim=1)) ])

    def ch(self):
        r"""
        Return the Chern character of the sheaf `\mathcal{O}(D)`
        defined by the divisor `D`.

        EXAMPLES::

            sage: dP6 = toric_varieties.dP6()
            sage: N = dP6.fan().lattice()
            sage: D3 = ToricDivisor(dP6, dP6.fan().cone_containing( N(0,1)   ))
            sage: D5 = ToricDivisor(dP6, dP6.fan().cone_containing( N(-1,-1) ))
            sage: D6 = ToricDivisor(dP6, dP6.fan().cone_containing( N(0,-1)  ))
            sage: D = -D3 + 2*D5 - D6
            sage: D.ch()
            [5*w^2 + y - 2*v + w + 1]
            sage: dP6.integrate( D.ch() * dP6.Td() )
            -4
        """
        return self.c1().exp()

     # def A(self, ring=ZZ):
     #    r"""
     #    Returns the Chow homology class of the divisor.
     #    """
     #    return sum([ self._div[i] * f.A(ring)
     #                 for [i,f] in enumerate(self._variety.fan(1)) ])

