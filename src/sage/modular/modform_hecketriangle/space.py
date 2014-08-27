r"""
Modular forms for Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ, infinity

from sage.modules.module import Module
from sage.categories.all import Modules
from sage.modules.free_module import FreeModule
from sage.modules.free_module_element import vector
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from hecke_triangle_groups import HeckeTriangleGroup
from abstract_space import FormsSpace_abstract

def canonical_parameters(group, base_ring, k, ep, n=None):
    r"""
    Return a canonical version of the parameters.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.space import canonical_parameters
        sage: canonical_parameters(5, ZZ, 20/3, int(1))
        (Hecke triangle group for n = 5, Integer Ring, 20/3, 1, 5)
    """

    if not (n is None):
        group = n

    if (group == infinity):
        group = HeckeTriangleGroup(infinity)
    else:
        try: 
            group = HeckeTriangleGroup(ZZ(group))
        except TypeError:
            group = HeckeTriangleGroup(group.n())

    n = group.n()
    k = QQ(k)
    if (ep == None):   
        ep = (-1)**(k*ZZ(n-2)/ZZ(2))
    ep = ZZ(ep)
    num = (k-(1-ep)*n/(n-2))*(n-2)/4
    try:
        num = ZZ(num)
    except TypeError:
        raise ValueError("Invalid or non-occuring weight!")

    return (group, base_ring, k, ep, n)


class QuasiMeromorphicModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiMeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 20/3, int(1))
            sage: QuasiMeromorphicModularForms(5, ZZ, 20/3, int(1)) == QuasiMeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiMeromorphicModularForms
            sage: MF = QuasiMeromorphicModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiMeromorphicModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi meromorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_space() == MF
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["quasi", "mero"])

class QuasiWeakModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiWeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(4, ZZ, 8, -1)
            sage: QuasiWeakModularForms(4, ZZ, 8, -1) == QuasiWeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiWeakModularForms
            sage: MF = QuasiWeakModularForms(4, ZZ, 8, 1)
            sage: MF
            QuasiWeakModularForms(n=4, k=8, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi weakly holomorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["quasi", "weak"])

class QuasiModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, ZZ, 10/3, -1)
            sage: QuasiModularForms(5, ZZ, 10/3) == QuasiModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(5, ZZ, 20/3, 1)
            sage: MF
            QuasiModularForms(n=5, k=20/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["quasi", "holo"])

    def quasi_part_gens(self, r=0):
        r"""
        Return a basis of ``self`` for the submodule
        of quasi modular forms of the form ``E2^r*f``,
        where ``f`` is a modular form.

        INPUT:

        - `r` -- An integer.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=5, k=6, ep=-1)
            sage: MF.default_prec(2)
            sage: MF.dimension()
            3
            sage: MF.quasi_part_gens(r=0)
            [1 - 37/(200*d)*q + O(q^2)]
            sage: MF.quasi_part_gens(r=0)[0] == MF.E6()
            True
            sage: MF.quasi_part_gens(r=1)
            [1 + 33/(200*d)*q + O(q^2)]
            sage: MF.quasi_part_gens(r=1)[0] == MF.E2()*MF.E4()
            True
            sage: MF.quasi_part_gens(r=2)
            []
            sage: MF.quasi_part_gens(r=3)
            [1 - 27/(200*d)*q + O(q^2)]
            sage: MF.quasi_part_gens(r=3)[0] == MF.E2()^3
            True
        """

        r = ZZ(r)
        if (r < 0 or 2*r > self._weight):
            return []

        gens = self.graded_ring().reduce_type("holo", degree=(self._weight-QQ(2*r), self._ep*(-1)**r)).gens()
        if (len(gens)>0):
            (x,y,z,d) = self.rat_field().gens()
            #gens = [ self.graded_ring().E2()**r*gen for gen in gens ]
            return [ self(z**r*gen._rat) for gen in gens ]
        else:
            return []

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=5, k=6, ep=-1)
            sage: MF.default_prec(2)
            sage: MF.gens()
            [1 - 37/(200*d)*q + O(q^2),
             1 + 33/(200*d)*q + O(q^2),
             1 - 27/(200*d)*q + O(q^2)]
        """

        gens = []
        for r in range(ZZ(0), QQ(self._weight/ZZ(2)).floor()+1):
            gens.extend(self.quasi_part_gens(r))

        return gens

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiModularForms
            sage: MF = QuasiModularForms(n=5, k=6, ep=-1)
            sage: MF.dimension()
            3
            sage: len(MF.gens()) == MF.dimension()
            True
        """
        n  = self.hecke_n()
        k  = self.weight()
        ep = self.ep()
        return sum([ 
            max(QQ((k-2*r)*(n-2)/(4*n) - (1-ep*(-1)**r)/4).floor() + 1, 0)\
            for r in range(ZZ(0), QQ(k / ZZ(2)).floor() + 1)])

    # TODO: it is possible to define coordinate_vector!
    # for this a routine needs to be written to additively decompose
    # the polynomial according to the degree of z
    # the basis vector with respect to the above basis is then given
    # by concatenating the coordinate vectors for each (decomposed) part
    # However the above gens() has no nice relation to the Fourier coefficients
    # and it is expected to be hard(er) to write an algorithm to determine
    # the form by its fourier coefficients

class QuasiCuspForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) quasi cusp forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, QuasiCuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(8, ZZ, 16/3, None)
            sage: QuasiCuspForms(8, ZZ, 16/3) == QuasiCuspForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) quasi cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(8, ZZ, 16/3)
            sage: MF
            QuasiCuspForms(n=8, k=16/3, ep=1) over Integer Ring
            sage: MF.analytic_type()
            quasi cuspidal
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["quasi", "cusp"])

    def quasi_part_gens(self, r=0):
        r"""
        Return a basis of ``self`` for the submodule
        of quasi cusp forms of the form ``E2^r*f``,
        where ``f`` is a cusp form.

        INPUT:

        - `r` -- An integer.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms, CuspForms
            sage: MF = QuasiCuspForms(n=5, k=18, ep=-1)
            sage: MF.default_prec(4)
            sage: MF.dimension()
            8
            sage: MF.quasi_part_gens(r=0)
            [q - 34743/(640000*d^2)*q^3 + O(q^4), q^2 - 69/(200*d)*q^3 + O(q^4)]
            sage: MF.quasi_part_gens(r=1)
            [q - 9/(200*d)*q^2 + 37633/(640000*d^2)*q^3 + O(q^4),
             q^2 + 1/(200*d)*q^3 + O(q^4)]
            sage: MF.quasi_part_gens(r=2)
            [q - 1/(4*d)*q^2 - 24903/(640000*d^2)*q^3 + O(q^4)]
            sage: MF.quasi_part_gens(r=3)
            [q + 1/(10*d)*q^2 - 7263/(640000*d^2)*q^3 + O(q^4)]
            sage: MF.quasi_part_gens(r=4)
            [q - 11/(20*d)*q^2 + 53577/(640000*d^2)*q^3 + O(q^4)]
            sage: MF.quasi_part_gens(r=5)
            [q - 1/(5*d)*q^2 + 4017/(640000*d^2)*q^3 + O(q^4)]

            sage: MF.quasi_part_gens(r=1)[0] == MF.E2() * CuspForms(n=5, k=16, ep=1).gen(0)
            True
            sage: MF.quasi_part_gens(r=1)[1] == MF.E2() * CuspForms(n=5, k=16, ep=1).gen(1)
            True
            sage: MF.quasi_part_gens(r=3)[0] == MF.E2()^3 * MF.Delta()
            True
        """

        r = ZZ(r)
        if (r < 0 or 2*r > self._weight):
            return []

        gens = self.graded_ring().reduce_type("cusp", degree=(self._weight-QQ(2*r), self._ep*(-1)**r)).gens()
        if (len(gens)>0):
            (x,y,z,d) = self.rat_field().gens()
            #gens = [ self.graded_ring().E2()**r*gen for gen in gens ]
            return [ self(z**r*gen._rat) for gen in gens ]
        else:
            return []

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(n=8, k=46/3, ep=-1)
            sage: MF.default_prec(4)
            sage: MF.dimension()
            7
            sage: MF.gens()
            [q - 17535/(262144*d^2)*q^3 + O(q^4),
             q^2 - 47/(128*d)*q^3 + O(q^4),
             q - 9/(128*d)*q^2 + 15633/(262144*d^2)*q^3 + O(q^4),
             q^2 - 7/(128*d)*q^3 + O(q^4),
             q - 23/(64*d)*q^2 - 3103/(262144*d^2)*q^3 + O(q^4),
             q - 3/(64*d)*q^2 - 4863/(262144*d^2)*q^3 + O(q^4),
             q - 27/(64*d)*q^2 + 17217/(262144*d^2)*q^3 + O(q^4)]
        """

        gens = []
        for r in range(ZZ(0), QQ(self._weight/ZZ(2)).floor()+1):
            gens.extend(self.quasi_part_gens(r))

        return gens

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import QuasiCuspForms
            sage: MF = QuasiCuspForms(n=8, k=46/3, ep=-1)
            sage: MF.default_prec(3)
            sage: MF.dimension()
            7
            sage: len(MF.gens()) == MF.dimension()
            True
        """

        n  = self.hecke_n()
        k  = self.weight()
        ep = self.ep()
        return sum([ 
            max( QQ( (k - 2*r)*(n - 2)/(4*n) - (1 - ep*(-1)**r)/4 ).floor() + 0, 0)\
            for r in range(ZZ(0), QQ(k/ZZ(2)).floor() + 1) ])

    # TODO: it is possible to define coordinate_vector! (see above)

class MeromorphicModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) meromorphic modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, MeromorphicModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, 1)
            sage: MeromorphicModularForms() == MeromorphicModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) meromorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import MeromorphicModularForms
            sage: MF = MeromorphicModularForms()
            sage: MF
            MeromorphicModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            meromorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["mero"])

class WeakModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) weakly holomorphic modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, WeakModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(5, CC, 20/3, None)
            sage: WeakModularForms(5, CC, 20/3) == WeakModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) weakly holomorphic modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import WeakModularForms
            sage: MF = WeakModularForms(5, CC, 20/3)
            sage: MF
            WeakModularForms(n=5, k=20/3, ep=1) over Complex Field with 53 bits of precision
            sage: MF.analytic_type()
            weakly holomorphic modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["weak"])

class ModularForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) modular forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, ModularForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(3, ZZ, 0, None)
            sage: ModularForms() == ModularForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) modular forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms()
            sage: MF
            ModularForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF.analytic_type()
            modular
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type = self.AT(["holo"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: MF.gens()
            [1 + 360360*q^4 + O(q^5),
             q + 21742*q^4 + O(q^5),
             q^2 + 702*q^4 + O(q^5),
             q^3 - 6*q^4 + O(q^5)]
        """

        return [ self.F_basis(m) for m in range(ZZ(0), -(self._l1 + 1), -1)]

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: len(MF.gens()) == MF.dimension()
            True
        """

        return max(self._l1+1, ZZ(0))

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.
        
        INPUT:

        - ``v``    - An element of ``self``.

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ModularForms
            sage: MF = ModularForms(n=6, k=20, ep=1)
            sage: MF.dimension()
            4
            sage: el = MF.E4()^2*MF.Delta()
            sage: el
            q + 78*q^2 + 2781*q^3 + 59812*q^4 + O(q^5)
            sage: vec = el.coordinate_vector()
            sage: vec
            (0, 1, 13/(18*d), 103/(432*d^2))
            sage: vec.parent()
            Vector space of dimension 4 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()
            True
            sage: el == vec[0]*MF.gen(0) + vec[1]*MF.gen(1) + vec[2]*MF.gen(2) + vec[3]*MF.gen(3)
            True
            sage: el == MF.element_from_coordinates(vec)
            True
        """

        vec = v.q_expansion(prec=self.degree()).add_bigoh(self.degree()).padded_list()
        return self._module(vector(self.coeff_ring(), vec))

class CuspForms(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Module of (Hecke) cusp forms
    for the given group, base ring, weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, CuspForms)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(6, ZZ, 6, 1)
            sage: CuspForms(6, ZZ, 6, 1) == CuspForms(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the Module of (Hecke) cusp forms
        of weight ``k`` with multiplier ``ep`` for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(6, ZZ, 6, 1)
            sage: MF
            CuspForms(n=6, k=6, ep=1) over Integer Ring
            sage: MF.analytic_type()
            cuspidal
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.module()
            Vector space of dimension 1 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT(["cusp"])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF=CuspForms(n=12, k=72/5, ep=1)
            sage: MF
            CuspForms(n=12, k=72/5, ep=1) over Integer Ring
            sage: MF.dimension()
            3
            sage: MF.gens()
            [q + 296888795/(10319560704*d^3)*q^4 + O(q^5),
             q^2 + 6629/(221184*d^2)*q^4 + O(q^5),
             q^3 - 25/(96*d)*q^4 + O(q^5)]
          """

        return [ self.F_basis(m) for m in range(ZZ(-1), -(self._l1 + 1), -1)]

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=12, k=72/5, ep=1)
            sage: MF.dimension()
            3
            sage: len(MF.gens()) == MF.dimension()
            True
        """

        return max(self._l1, ZZ(0))

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.
        
        INPUT:

        - ``v``    - An element of ``self``.

        OUTPUT:

        An element of ``self.module()``, namely the
        corresponding coordinate vector of ``v`` with respect
        to the basis ``self.gens()``.

        The module is the free module over the coefficient
        ring of ``self`` with the dimension of ``self``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import CuspForms
            sage: MF = CuspForms(n=12, k=72/5, ep=-1)
            sage: MF.default_prec(4)
            sage: MF.dimension()
            2
            sage: el = MF(MF.f_i()*MF.Delta())
            sage: el
            q - 1/(288*d)*q^2 - 96605/(1327104*d^2)*q^3 + O(q^4)
            sage: vec = el.coordinate_vector()
            sage: vec
            (1, -1/(288*d))
            sage: vec.parent()
            Vector space of dimension 2 over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
            sage: vec.parent() == MF.module()
            True
            sage: el == vec[0]*MF.gen(0) + vec[1]*MF.gen(1)
            True
            sage: el == MF.element_from_coordinates(vec)
            True
        """

        vec = v.q_expansion(prec=self.degree()+1).add_bigoh(self.degree()+1).padded_list()
        vec.pop(0)
        return self._module(vector(self.coeff_ring(), vec))

class ZeroForm(FormsSpace_abstract, Module, UniqueRepresentation):
    r"""
    Zero Module for the zero form for the given group, base ring
    weight and multiplier
    """
            
    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, k=QQ(0), ep=None, n=None):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import (canonical_parameters, ZeroForm)
            sage: (group, base_ring, k, ep, n) = canonical_parameters(6, CC, 3, -1)
            sage: ZeroForm(6, CC, 3, -1) == ZeroForm(group, base_ring, k, ep, n)
            True
        """

        (group, base_ring, k, ep, n) = canonical_parameters(group, base_ring, k, ep, n)
        return super(FormsSpace_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, k=k, ep=ep, n=n)

    def __init__(self, group, base_ring, k, ep, n):
        r"""
        Return the zero Module for the zero form of weight ``k`` with multiplier ``ep``
        for the given ``group`` and ``base_ring``.

        The ZeroForm space coerces into any forms space or ring with a compatible group.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm(6, CC, 3, -1)
            sage: MF
            ZeroForms(n=6, k=3, ep=-1) over Complex Field with 53 bits of precision
            sage: MF.analytic_type()
            zero
            sage: MF.category()
            Category of vector spaces over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
            sage: MF.module()
            Vector space of dimension 0 over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
            sage: MF.ambient_module() == MF.module()
            True
            sage: MF.is_ambient()
            True
        """

        FormsSpace_abstract.__init__(self, group=group, base_ring=base_ring, k=k, ep=ep, n=n)
        Module.__init__(self, base=self.coeff_ring())
        self._analytic_type=self.AT([])
        self._module = FreeModule(self.coeff_ring(), self.dimension())

    def _change_degree(self, k, ep):
        r"""
        Return the Zeroform space with a different weight ``k`` and multiplier ``ep``
        for the same group and base_ring.

        All those spaces coerce into each other.

        INPUT:

        - ``k``   -- A rational number, the weight.

        - ``ep``  -- ``1`` or ``-1``, the multiplier.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm()
            sage: MF
            ZeroForms(n=3, k=0, ep=1) over Integer Ring
            sage: MF._change_degree(14, -1)
            ZeroForms(n=3, k=14, ep=-1) over Integer Ring
        """
        return ZeroForm(group=self.group(), base_ring=self.base_ring(), k=k, ep=ep)

    @cached_method
    def gens(self):
        r"""
        Return a basis of ``self`` as a list of basis elements.
        Since this is the zero module an empty list is returned.

        EXAMPLES::
 
            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: ZeroForm(6, CC, 3, -1).gens()
            []
        """

        return []

    @cached_method
    def dimension(self):
        r"""
        Return the dimension of ``self``.
        Since this is the zero module ``0`` is returned.

        EXAMPLES::
 
            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: ZeroForm(6, CC, 3, -1).dimension()
            0
        """

        return 0

    @cached_method
    def coordinate_vector(self, v):
        r"""
        Return the coordinate vector of ``v`` with respect to
        the basis ``self.gens()``.

        Since this is the zero module which only contains
        the zero form the trivial vector in the trivial
        module of dimension ``0`` is returned.

        INPUT:

        - ``v``    - An element of ``self``, i.e. in this case the zero vector.

        EXAMPLES::
 
            sage: from sage.modular.modform_hecketriangle.space import ZeroForm
            sage: MF = ZeroForm(6, CC, 3, -1)
            sage: el = MF(0)
            sage: el
            O(q^5)
            sage: vec = el.coordinate_vector()
            sage: vec
            ()
            sage: vec.parent()
            Vector space of dimension 0 over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
            sage: vec.parent() == MF.module()
            True
        """

        vec = []
        return self._module(vector(self.coeff_ring(), vec))
