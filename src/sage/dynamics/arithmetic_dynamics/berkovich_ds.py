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
                                Berkovich_Cp_Projective, is_Berkovich_Cp, Berkovich_Element_Cp)
from sage.rings.padics.factory import Qp
from sage.structure.element import get_coercion_model
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.affine.affine_space import is_AffineSpace
from sage.rings.padics.padic_base_generic import pAdicBaseGeneric
from sage.dynamics.arithmetic_dynamics.projective_ds import DynamicalSystem_projective
from sage.dynamics.arithmetic_dynamics.affine_ds import DynamicalSystem_affine

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
        if isinstance(dynamical_system, SchemeMorphism_polynomial):
            morphism_domain = dynamical_system.domain()

        if not domain is None:
            if isinstance(domain, Berkovich_Cp_Affine):
                return DynamicalSystem_Berkovich_affine(dynamical_system,domain)
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
            raise ValueError('action of dynamical system not defined on %s' %x.parent())
        if x.type_of_point() == 1:
            return self.domain()(self._system(x.center()))
        if x.type_of_point() == 4:
            raise NotImplementedError('action on Type IV points not implemented')
        #TODO write a better check for zeros in disk

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
        if not isinstance(dynamical_system, DynamicalSystem_projective):
            dynamical_system = DynamicalSystem_projective(dynamical_system)
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