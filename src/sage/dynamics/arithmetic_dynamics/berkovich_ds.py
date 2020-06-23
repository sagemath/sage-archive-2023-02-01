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

    @staticmethod
    def __classcall_private__(cls, system_morphism_polys, domain=None):
        r"""
        Base class for dynamical systems on Berkovich space over \CC_p.

        INPUT::

        - ``system_morphism_polys`` -- a list of polynomials or rational functions,
          or a dynamical system. In any case, this input should be
          defined over affine or projective space of a finite extension of \QQ_p.

        - ``domain`` -- (optional) affine or projective Berkovich space over \CC_p
        """
        if not(domain is None or is_Berkovich_Cp(domain)):
            raise ValueError('domain must be a Berkovich space over Cp')
        if isinstance(system_morphism_polys,(list,tuple)):
            if len(system_morphism_polys) not in [1,2]:
                raise ValueError('list of polynomials too long, must be length 1 or 2')
            from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
            from sage.rings.polynomial.multi_polynomial_ring_base import is_MPolynomialRing
            test = lambda x: is_PolynomialRing(x) or is_MPolynomialRing(x)
            polys = list(system_morphism_polys)
            if not all(test(poly.parent()) for poly in polys):
                try:
                    polys = [poly.lift() for poly in polys]
                except AttributeError:
                    raise ValueError('{} must be elements of a polynomial ring'\
                        .format(system_morphism_polys))
            for poly in polys:
                if not isinstance(poly.base_ring(), pAdicBaseGeneric):
                    if domain is None:
                        raise ValueError('polynomials not defined over a padic with no specified domain')
                    try:
                        poly = poly.change_ring(Qp(domain.prime())) #TODO change to Qpbar
                        flag = False
                    except:
                        flag = True
                    if flag:
                        raise ValueError('{} does not convert to Qp'.format(poly))
            PR = get_coercion_model().common_parent(*polys)
            if isinstance(domain, Berkovich_Cp_Affine):
                if domain.prime() != PR.prime():
                    raise ValueError('specified domain has an incorrect residue characteristic')
                if len(polys) != 1:
                    raise ValueError('list of polynomials too long for affine Berkovich space')
                if len(polys[0].gens()) != 1:
                    raise ValueError('too many variables for dynamical system' + \
                        'on affine Berkovich space')
                system = DynamicalSystem_affine(polys)
                return typecall(cls,system,domain)
            else:
                if domain is None:
                    P = ProjectiveSpace(PR, 1)
                    domain = Berkovich_Cp_Projective(P)
                if domain.prime() != PR.prime():
                    raise ValueError('specified domain has an incorrect residue characteristic')
                if not (len(polys) in [1,2]):
                    raise ValueError('list of polynomials too long for affine Berkovich space')
                if not(len(polys[0].gens()) in [1,2]):
                    raise ValueError('too many variables for dynamical system' + \
                        'on projective Berkovich space')
                system = DynamicalSystem_projective(polys)
                return typecall(cls,system,domain)

        if isinstance(system_morphism_polys, SchemeMorphism_polynomial):
            R = system_morphism_polys.base_ring()
            morphism_domain = system_morphism_polys.domain()
            if morphism_domain != system_morphism_polys.codomain():
                raise ValueError('domain and codomain do not agree')
            if not isinstance(R, pAdicBaseGeneric):
                if domain is None:
                    raise ValueError('system_morphism_polys not defined over padic field and domain is None')
                try:
                    #TODO change to Qpbar
                    system_morphism_polys = system_morphism_polys.change_ring(Qp(domain.prime()))
                    morphism_domain = system_morphism_polys.domain()
                    R = system_morphism_polys.base_ring()
                    flag = False
                except:
                    flag = True
                if flag:
                    raise ValueError('system_morphism_polys could not be converted to Qp(%s)' %domain.prime())
            if morphism_domain != morphism_domain.ambient_space():
                raise ValueError('morphism must be defined on the ambient space')
            if morphism_domain.dimension_absolute() != 1:
                raise ValueError('domain not dimension 1')
            if is_AffineSpace(morphism_domain):
                domain = Berkovich_Cp_Affine(R)
            else:
                domain = Berkovich_Cp_Projective(ProjectiveSpace(R, 1))
            if not isinstance(system_morphism_polys, DynamicalSystem):
                system_morphism_polys = DynamicalSystem(system_morphism_polys)
            return typecall(cls,system_morphism_polys,domain)
        
        raise ValueError('system_morphism_polys was not a dynamical system, a morphism, or a polynomial')

    def __init__(self, dynamical_system, domain):
        r"""
        The Python constructor
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
        """
        if not isinstance(x.parent(), type(self._domain)):
            raise ValueError('action of dynamical system not defined on %s' %x)
        if x.type_of_point() != 1:
            raise NotImplementedError('action on Type II, III, and IV points not implemented')
        return self._system(x.center())

    def _repr_(self):
        """
        Return a string representation of this dynamical system.

        OUTPUT: a string

        EXAMPLES::

            sage: Q.<x,y> = ProjectiveSpace(Qp(3),1)
            sage: f = DynamicalSystem_projective([3*x^2,2*y^2])
            sage: f = DynamicalSystem_Berkovich(f)
            sage: f._repr_()
            'Dynamical system of Projective Space of dimension 1 over 
            3-adic Field with capped relative precision 20 induced by the map
             \n  Defn: Defined on coordinates by sending (x : y) to\n        
             ((3 + O(3^21))*x^2 : (2 + O(3^20))*y^2)'
        """
        domain_str = self._domain._repr_()
        return "Dynamical system of " + domain_str + " induced by the map " + \
            "\n  Defn: %s"%('\n        '.join(self._system._repr_defn().split('\n')))