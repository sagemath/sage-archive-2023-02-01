r"""
Dynamical systmes on Berkovich space over \CC_p.

A dynamical system on Berkovich space over \CC_p is 
determined by a dynamical system on A^1(\CC_p) or P^1(\CC_p),
which naturally induces a dynamical system on affine or
projective Berkovich space.
"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element
from sage.dynamics.arithmetic_dynamics.generic_ds import DynamicalSystem
from six import add_metaclass
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.classcall_metaclass import typecall
from sage.schemes.generic.morphism import SchemeMorphism_polynomial
from sage.rings.padics.generic_nodes import is_pAdicField
from sage.schemes.berkovich.berkovich_space import (Berkovich_Cp_Affine,
                                Berkovich_Cp_Projective, is_Berkovich_Cp,
                                Berkovich_Element_Cp_Affine, Berkovich_Element_Cp_Projective)
from sage.rings.padics.factory import Qp
from sage.structure.element import get_coercion_model
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine
from sage.rings.rational_field import QQ

@add_metaclass(InheritComparisonClasscallMetaclass)
class DynamicalSystem_Berkovich(Element):
    r"""
    A dynamical system on Berkovich space over `\CC_p`.

    A dynamical system on Berkovich space over `\CC_p` is 
    determined by a dynamical system on `A^1(\CC_p)` or `P^1(\CC_p)`,
    which naturally induces a dynamical system on affine or
    projective Berkovich space.

    INPUT::

    - ``dynamical_system`` -- any input which defines a dynamical
        system over affine or projective space, or a dynamical system
        over affine or projective space. If the input is a list
        of homogenous polynomials, then ``domain`` is taken to
        be projective Berkovich space, unless specified.
        If this input is not defined over a padic field,
        then ``domain`` MUST be specified.

    - ``domain`` -- (optional) affine or projective Berkovich space 
        over `\CC_p`. If the input to ``dynamical_system`` is 
        not defined over `\QQ_p` or a finite extension, ``domain``
        must be specified.

    EXAMPLES:

    We can easily create a dynamical system on Berkovich space
    using a dynamical system on projective space over `\QQ_p`::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
        sage: DynamicalSystem_Berkovich(f)
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map
            Defn: Defined on coordinates by sending (x : y) to
                ((2 + O(3^20))*x^2 + (1 + 3 + O(3^20))*y^2 : (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)

    Or from a morphism::

        sage: P1.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: H = End(P1)
        sage: DynamicalSystem_Berkovich(H([y, x]))
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map 
            Defn: Defined on coordinates by sending (x : y) to
                (y : x)

    Or from polynomials::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: DynamicalSystem_Berkovich([x^2+y^2, y^2])
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map
            Defn: Defined on coordinates by sending (x : y) to
                (x^2 + y^2 : y^2)

    Note that the default behavior on polynomial input is to construct
    a dynamical system on the projective line. To construct a dynamical
    system on on the affine line, specify ``domain``::

        sage: A.<x> = AffineSpace(Qp(3),1)
        sage: B = Berkovich_Cp_Affine(3)
        sage: DynamicalSystem_Berkovich([x^2+1], B)
        Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map 
            Defn: Defined on coordinates by sending (x) to
                (x^2 + 1 + O(3^20))

    ``domain`` is ignored if a dynamical system or an endomorphism is
    passsed in, unless that morphism is not defined over a padic ring/field::

        sage: f = DynamicalSystem_affine(x^2+1)
        sage: C = Berkovich_Cp_Projective(3)
        sage: DynamicalSystem_Berkovich(f, C)
        Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (x) to
                (x^2 + 1 + O(3^20))

    Creating a map on Berkovich space creates the Berkovich space
    it acts on::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: f = DynamicalSystem_projective([x^2, y^2])
        sage: g = DynamicalSystem_Berkovich(f)
        sage: B = g.domain(); B
        Projective Berkovich line over Cp(3) of precision 20

    We can take the image of points of the domain::

        sage: Q1 = B(2)
        sage: g(Q1)
        Type I point centered at (1 + 3 + O(3^20) : 1 + O(3^20))
    """

    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None):
        """
        Returns the appropriate dynamical system on Berkovich space.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: B = Berkovich_Cp_Affine(3)
            sage: DynamicalSystem_Berkovich(t^2 - 3,B)
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
              Defn: Defined on coordinates by sending (t) to
                    (t^2 + 2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^7 + 2*3^8 + 2*3^9 + 2*3^10 + 2*3^11 + 2*3^12 + 2*3^13 + 2*3^14 + 2*3^15 + 2*3^16 + 2*3^17 + 2*3^18 + 2*3^19 + 2*3^20 + O(3^21))
        """
        if not(domain is None or is_Berkovich_Cp(domain)):
            raise ValueError('domain must be a Berkovich space over Cp')
        morphism_domain = None
        if isinstance(dynamical_system, SchemeMorphism_polynomial):
            morphism_domain = dynamical_system.domain()

        if not domain is None:
            if isinstance(domain, Berkovich_Cp_Affine):
                return DynamicalSystem_Berkovich_affine(dynamical_system,domain)
        if not morphism_domain is None:
            from sage.schemes.affine.affine_subscheme import AlgebraicScheme_subscheme_affine
            if is_AffineSpace(morphism_domain) or isinstance(domain, AlgebraicScheme_subscheme_affine):
                return DynamicalSystem_Berkovich_affine(dynamical_system,domain)

        return DynamicalSystem_Berkovich_projective(dynamical_system,domain)

    def __init__(self, dynamical_system, domain):
        r"""
        The Python constructor

        TESTS::

            sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: isinstance(g, DynamicalSystem_Berkovich)
            True
        """
        self._system = dynamical_system
        self._polys = dynamical_system._polys
        self._domain = domain

    def domain(self):
        """
        Returns the domain of this dynamical system.

        OUTPUT: A Berkovich space over ``Cp``

        EXAMPLES::

            sage: Q.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([3*x^2,2*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.domain()
            Projective Berkovich line over Cp(3) of precision 20
        """
        return self._domain

    def defining_polynomials(self):
        """
        Return the defining polynomials.

        OUTPUT:

        A tuple of polynomials that defines the
        dynamical system.

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([2*x^2 + 4*y^2, 3*x^2 + 9*y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.defining_polynomials()
            ((2 + O(3^20))*x^2 + (1 + 3 + O(3^20))*y^2,
            (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)
        """
        return self._polys

    def _repr_(self):
        r"""
        Return a string representation of this dynamical system.

        OUTPUT: a string

        EXAMPLES::

            sage: Q.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([3*x^2,2*y^2])
            sage: f = DynamicalSystem_Berkovich(f)
            sage: f._repr_()
            'Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map\n
              Defn: Defined on coordinates by sending (x : y) to\n        ((3 + O(3^21))*x^2 : (2 + O(3^20))*y^2)'
        """
        domain_str = self._domain._repr_()
        return "Dynamical system of " + domain_str + " induced by the map" + \
            "\n  Defn: %s"%('\n        '.join(self._system._repr_defn().split('\n')))

class DynamicalSystem_Berkovich_projective(DynamicalSystem_Berkovich):
    r"""
    A dynamical system on projective Berkovich space over `\CC_p`.

    A dynamical system on projective Berkovich space over `\CC_p` is 
    determined by a dynamical system on `A^1(\CC_p)` or `P^1(\CC_p)`,
    which naturally induces a dynamical system on affine or
    projective Berkovich space.

    INPUT::

    - ``dynamical_system`` -- any input which defines a dynamical
        system over affine or projective space, or a dynamical system
        over affine or projective space. If the input is a list
        of homogenous polynomials, then ``domain`` is taken to
        be projective Berkovich space, unless specified.
        If this input is not defined over a padic field,
        then ``domain`` MUST be specified.

    - ``domain`` -- (optional) affine or projective Berkovich space 
        over `\CC_p`. If the input to ``dynamical_system`` is 
        not defined over `\QQ_p` or a finite extension, ``domain``
        must be specified.

    EXAMPLES:

    We can easily create a dynamical system on Berkovich space
    using a dynamical system on projective space over `\QQ_p`::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: f = DynamicalSystem_projective([1/2*x^2 + x*y + 3*y^2, 3*x^2 + 9*y^2])
        sage: DynamicalSystem_Berkovich_projective(f)
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map
          Defn: Defined on coordinates by sending (x : y) to
                ((2 + 3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 3^7 + 3^8 + 3^9 + 3^10 + 3^11 + 3^12 + 3^13 + 3^14 + 3^15 + 3^16 + 3^17 + 3^18 + 3^19 + O(3^20))*x^2 + x*y + (3 + O(3^21))*y^2 : (3 + O(3^21))*x^2 + (3^2 + O(3^22))*y^2)

    Or from a morphism::

        sage: P1.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: H = End(P1)
        sage: DynamicalSystem_Berkovich_projective(H([y, x]))
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map 
            Defn: Defined on coordinates by sending (x : y) to
                (y : x)

    Or from polynomials::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: DynamicalSystem_Berkovich_projective([x^2+y^2, y^2])
        Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map
            Defn: Defined on coordinates by sending (x : y) to
                (x^2 + y^2 : y^2)

    Creating a map on Berkovich space creates the Berkovich space
    it acts on::

        sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
        sage: f = DynamicalSystem_projective([x^2, y^2])
        sage: g = DynamicalSystem_Berkovich(f)
        sage: B = g.domain(); B
        Projective Berkovich line over Cp(3) of precision 20

    We can take the image of points of the domain::

        sage: Q1 = B(2)
        sage: g(Q1)
        Type I point centered at (1 + 3 + O(3^20) : 1 + O(3^20))
    """
    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None):
        """
        Returns the approapriate dynamical system on projective Berkovich space over ``Cp``.
        """
        if not isinstance(dynamical_system, DynamicalSystem):
            if not isinstance(dynamical_system, DynamicalSystem_projective):
                dynamical_system = DynamicalSystem_projective(dynamical_system)
            else:
                raise ValueError('affine dynamical system passed to projective constructor')
        R = dynamical_system.base_ring()
        morphism_domain = dynamical_system.domain()
        if not isinstance(R, pAdicBaseGeneric):
            if domain is None:
                raise ValueError('dynamical system not defined over padic field and domain is None')
            try:
                #TODO change to Qpbar
                dynamical_system = dynamical_system.change_ring(Qp(domain.prime()))
                morphism_domain = dynamical_system.domain()
                R = dynamical_system.base_ring()
                flag = False
            except:
                flag = True
            if flag:
                raise ValueError('dynamical_system could not be converted to Qp(%s)' %domain.prime())
        if morphism_domain != morphism_domain.ambient_space():
            raise ValueError('morphism must be defined on the ambient space')
        if morphism_domain.dimension_absolute() != 1:
            raise ValueError('domain not dimension 1')
        domain = Berkovich_Cp_Projective(ProjectiveSpace(R, 1))
        return typecall(cls,dynamical_system,domain)

    def __init__(self, dynamical_system, domain=None):
        """
        Python constructor.
        """
        DynamicalSystem_Berkovich.__init__(self, dynamical_system, domain)

    def dehomogenize(self, n):
        """
        Returns the map induced by the standard dehomogenization.

        The dehomogenization is done at the ``n[0]`` coordinate
        of the domain and the ``n[1]`` coordinate of the codomain.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers; if ``n`` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: A dynamical system on affine Berkovich space

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([x^2 + y^2, x*y + y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: g.dehomogenize(1)
            Dynamical system of Affine Berkovich line over Cp(3) of precision 20 induced by the map
              Defn: Defined on coordinates by sending (x) to
                    ((x^2 + 1 + O(3^20))/(x + 1 + O(3^20)))
        """
        new_system = self._system.dehomogenize(n)
        base_ring = self.domain().base_ring().base_ring() #2 base rings since this is projective berkovich space
        new_domain = Berkovich_Cp_Affine(base_ring)
        return DynamicalSystem_Berkovich_affine(new_system,new_domain)

    def __call__(self, x):
        """
        Makes dynamical systems on Berkovich space over ``Cp`` callable.

        INPUT:

        - ``x`` -- a point of Berkovich space over ``Cp``

        EXAMPLES::

            sage: P.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([x^2, y^2])
            sage: g = DynamicalSystem_Berkovich(f)
            sage: B = g.domain()
            sage: Q1 = B(2)
            sage: g(Q1)
            Type I point centered at (1 + 3 + O(3^20) : 1 + O(3^20))
        """
        if not isinstance(x.parent(), type(self._domain)):
            try:
                x = self.domain()(x)
            except:
                raise ValueError('action of dynamical system not defined on %s' %x.parent())
        if x.type_of_point() == 1:
            return self.domain()(self._system(x.center()))
        if x.type_of_point() == 4:
            raise NotImplementedError('action on Type IV points not implemented')
        f = self._system
        if x.type_of_point() == 2:
            from sage.matrix.constructor import Matrix
            from sage.modules.free_module_element import vector
            y,w = f.domain().gens()[0],f.domain().gens()[1]
            field = f.domain().base_ring()
            M = Matrix([[field(x.prime()**(-1*x.power())),x.center()[0]],[field(0),field(1)]])
            X = M * vector(f[0].parent().gens())
            F = vector(f._polys)
            F = list(F(list(X)))
            R = field['z']
            z = R.gens()[0]
            for i in range(len(F)):
                F[i] = F[i].subs({y:z,w:1})
            print(F)
            lcm = field(1)
            for poly in F:
                for i in poly:
                    if i != 0:
                        lcm = i.denominator().lcm(lcm)
            for i in range(len(F)):
                F[i] *= lcm
            gcd = [i for i in F[0] if i != 0][0]
            for poly in F:
                for i in poly:
                    if i != 0:
                        gcd = gcd*i*gcd.lcm(i).inverse_of_unit()
            for i in range(len(F)):
                F[i] *= gcd.inverse_of_unit()
            gcd = F[0].gcd(F[1])
            F[0] = F[0].quo_rem(gcd)[0]
            F[1] = F[1].quo_rem(gcd)[0]
            fraction = []
            for poly in F:
                new_poly = []
                for i in poly:
                    new_poly.append((i).residue())
                new_poly = R(new_poly)
                fraction.append((new_poly))
            gcd = fraction[0].gcd(fraction[1])
            num = fraction[0].quo_rem(gcd)[0]
            dem = fraction[1].quo_rem(gcd)[0]
            if dem.is_zero():
                f = DynamicalSystem_affine(F[0]/F[1]).homogenize(1)
                f = f.conjugate(Matrix([[0, 1], [1 , 0]]))
                g = DynamicalSystem_Berkovich(f)
                return g(self.domain()(QQ(0),QQ(1))).involution_map()
            #if the reduction is not constant, the image
            #is the Gauss point
            if not(num.is_constant() and dem.is_constant()):
                return self.domain()(QQ(0),QQ(1))
            reduced_value = field(num*dem.inverse_of_unit()).lift_to_precision(field.precision_cap())
            new_num = F[0]-reduced_value*F[1]
            power_of_p = min([i.valuation() for i in new_num])
            inverse_map = (field(x.prime()**power_of_p)*z + reduced_value)
            return self.domain()(inverse_map(0),(inverse_map(1)-inverse_map(0)).abs())
        #TODO write a better check for zeros in disk
        P = f.domain()
        if P.gens()[0] in f.defining_polynomials()[1].variables():
            raise ValueError('action on Type II/III points only implemented for polynomials')
        nth_derivative = f.dehomogenize(1).defining_polynomials()[0]
        variable = nth_derivative.parent().gens()[0]
        a = x.center()[0]
        Taylor_expansion = []
        from sage.functions.other import factorial
        for i in range(f.degree()+1):
            Taylor_expansion.append(nth_derivative(a)*1/factorial(i))
            nth_derivative = nth_derivative.derivative(variable)
        r = x.radius()
        new_center = f(a)
        new_radius = max([Taylor_expansion[i].abs()*r**i for i in range(1,len(Taylor_expansion))])
        return self.domain()(new_center, new_radius)

class DynamicalSystem_Berkovich_affine(DynamicalSystem_Berkovich):
    @staticmethod
    def __classcall_private__(cls, dynamical_system, domain=None):
        if not isinstance(dynamical_system, DynamicalSystem_affine):
            dynamical_system = DynamicalSystem_affine(dynamical_system)
        R = dynamical_system.base_ring()
        morphism_domain = dynamical_system.domain()
        if not isinstance(R, pAdicBaseGeneric):
            if domain is None:
                raise ValueError('dynamical system not defined over padic field and domain is None')
            try:
                #TODO change to Qpbar
                dynamical_system = dynamical_system.change_ring(Qp(domain.prime()))
                morphism_domain = dynamical_system.domain()
                R = dynamical_system.base_ring()
                flag = False
            except:
                flag = True
            if flag:
                raise ValueError('dynamical_system could not be converted to Qp(%s)' %domain.prime())
        if morphism_domain != morphism_domain.ambient_space():
            raise ValueError('morphism must be defined on the ambient space')
        if morphism_domain.dimension_absolute() != 1:
            raise ValueError('domain not dimension 1')
        domain = Berkovich_Cp_Affine(R)
        return typecall(cls,dynamical_system,domain)

    def __init__(self, dynamical_system, domain):
        """
        Python constructor.
        """
        DynamicalSystem_Berkovich.__init__(self, dynamical_system, domain)

    def homogenize(self, n):
        """
        Returns the homogenization of this dynamical system.

        For dynamical systems on Berkovich space, this is the dynamical
        system on projective space induced by the homogenization of
        the dynamical system.

        INPUT:

        - ``n`` -- a tuple of nonnegative integers. If ``n`` is an integer,
          then the two values of the tuple are assumed to be the same

        OUTPUT: a dynamical system on projective Berkovich space

        EXAMPLES::

            sage: A.<x> = AffineSpace(Qp(3),1)
            sage: f = DynamicalSystem_affine(1/x)
            sage: f = DynamicalSystem_Berkovich(f)
            sage: f.homogenize(1)
            Dynamical system of Projective Berkovich line over Cp(3) of precision 20 induced by the map
                  Defn: Defined on coordinates by sending (x0 : x1) to
                        (x1 : x0)

        """
        new_system = self._system.homogenize(n)
        base_ring = self.domain().base_ring()
        new_domain = Berkovich_Cp_Affine(base_ring)
        return DynamicalSystem_Berkovich_projective(new_system,new_domain)

    def __call__(self, x):
        """
        Makes this dynamical system callable.

        EXAMPLES::

            sage: P.<x> = AffineSpace(Qp(3),1)
            sage: f = DynamicalSystem_affine(x^2)
            sage: g = DynamicalSystem_Berkovich(f)
            sage: B = g.domain()
            sage: Q1 = B(2)
            sage: g(Q1)
            Type I point centered at 1 + 3 + O(3^20)
        """
        if not isinstance(x, Berkovich_Element_Cp_Affine):
            try:
                x = self.domain()(x)
            except:
                raise ValueError('action of dynamical system not defined on %s' %x)
        proj_system = self.homogenize(1)
        return self.domain()(proj_system(x))