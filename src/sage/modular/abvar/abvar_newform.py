"""
TODO
"""


###########################################################################
#       Copyright (C) 2008 William Stein <wstein@gmail.com>               #
#  Distributed under the terms of the GNU General Public License (GPL)    #
#                  http://www.gnu.org/licenses/                           #
###########################################################################

from sage.databases.cremona import cremona_letter_code

from sage.rings.all  import QQ, ZZ

from sage.modular.modform.element import Newform
from sage.modular.congroup import is_Gamma0, is_Gamma1, is_GammaH

from abvar import ModularAbelianVariety_modsym_abstract
import homspace

class ModularAbelianVariety_newform(ModularAbelianVariety_modsym_abstract):
    """
    A modular abelian variety attached to a specific newform.
    """
    def __init__(self, f):
        """
        Create the modular abelian variety $A_f$ attached to the
        newform $f$.

        INPUT:
            f -- a newform

        EXAMPLES:
            sage: f = CuspForms(37).newforms('a')[0]
            sage: f.abelian_variety()
            Modular abelian variety attached to the newform q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
        """
        if not isinstance(f, Newform):
            raise TypeError, "f must be a newform"
        self.__f = f
        self._is_hecke_stable = True
        ModularAbelianVariety_modsym_abstract.__init__(self, QQ)

    def _modular_symbols(self,sign=0):
        """
        EXAMPLES:
            sage: f = CuspForms(52).newforms('a')[0]
            sage: A = f.abelian_variety()
            sage: A._modular_symbols()
            Modular Symbols subspace of dimension 2 of Modular Symbols space of dimension 15 for Gamma_0(52) of weight 2 with sign 0 over Rational Field
            sage: A._modular_symbols(1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 10 for Gamma_0(52) of weight 2 with sign 1 over Rational Field
            sage: A._modular_symbols(-1)
            Modular Symbols subspace of dimension 1 of Modular Symbols space of dimension 5 for Gamma_0(52) of weight 2 with sign -1 over Rational Field
        """
        return self.__f.modular_symbols(sign=sign)

    def newform(self):
        """
        Return the newform that this modular abelian variety is attached to.

        EXAMPLES:
            sage: f = CuspForms(37).newforms('a')[0]
            sage: A = f.abelian_variety()
            sage: A.newform()
            q - 2*q^2 - 3*q^3 + 2*q^4 - 2*q^5 + O(q^6)
            sage: A.newform() is f
            True
        """
        return self.__f

    def label(self):
        """
        Return canonical label that defines this newform modular
        abelian variety.

        OUTPUT:
            string

        EXAMPLES:
            sage: f = CuspForms(43).newforms('a')[1]
            sage: A = f.abelian_variety()
            sage: A.label()
            '43b'
        """
        G = self.__f.group()
        if is_Gamma0(G):
            group = ''
        elif is_Gamma1(G):
            group = 'G1'
        elif is_GammaH(G):
            group = 'GH[' + ','.join([str(z) for z in G._generators_for_H()]) + ']'
        return '%s%s%s'%(self.level(), cremona_letter_code(self.factor_number()), group)

    def factor_number(self):
        """
        Return factor number.

        OUTPUT:
            int

        EXAMPLES:
            sage: f = CuspForms(43).newforms('a')[1]
            sage: A = f.abelian_variety()
            sage: A.factor_number()
            1
        """
        try:
            return self.__factor_number
        except AttributeError:
            self.__factor_number = self.__f.number()
            return self.__factor_number

    def _repr_(self):
        """
        String representation of this modular abelian variety.

        EXAMPLES:
            sage: f = CuspForms(43).newforms('a')[1]
            sage: A = f.abelian_variety()
            sage: A._repr_()
            'Modular abelian variety attached to the newform q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)'
        """
        return "Modular abelian variety attached to the newform %s"%self.newform()

    def endomorphism_ring(self):
        """
        EXAMPLES:
            sage: f = CuspForms(43).newforms('a')[1]
            sage: f
            q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)
            sage: A = f.abelian_variety()
            sage: B = A.endomorphism_ring(); B
            Subring of endomorphism ring of Modular abelian variety attached to the newform q + a1*q^2 - a1*q^3 + (-a1 + 2)*q^5 + O(q^6)
            sage: B.gens()
            ([1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1],
             [ 0  1  0  0]
            [ 1 -2  0  0]
            [ 1  0 -2 -1]
            [ 0  1 -1  0])
        """
        try:
            return self.__endomorphism_ring
        except AttributeError:
            pass

        M = self.modular_symbols()
        bound = M.sturm_bound()

        d = self.dimension()
        T1list = self.hecke_operator(1).matrix().list()
        EndVecZ = ZZ**(len(T1list))
        V = EndVecZ.submodule([T1list])
        n = 2

        while V.dimension() < d:
            W = EndVecZ.submodule([((self.hecke_operator(n).matrix())**i).list()
                                   for i in range(1,d+1)])
            V = V+W
            n += 1

        E = homspace.EndomorphismSubring(self)
        E._set_generators(V.saturation().basis())
        self.__endomorphism_ring = E
        return self.__endomorphism_ring



