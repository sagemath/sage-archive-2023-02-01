"""
The class group of a number field.

EXAMPLES:
    sage: K.<a> = NumberField(x^2 + 23)
    sage: I = K.class_group().gen(); I
    Fractional ideal class (2, 1/2*a - 1/2) of Number Field in a with defining polynomial x^2 + 23
    sage: J = I * I; J
    Fractional ideal class (2, 1/2*a + 1/2) of Number Field in a with defining polynomial x^2 + 23
    sage: O = K.OK(); O
    Order with module basis 1/2*a + 1/2, a in Number Field in a with defining polynomial x^2 + 23
    sage: O*(2, 1/2*a + 1/2)
    Fractional ideal (2, 1/2*a + 1/2) of Number Field in a with defining polynomial x^2 + 23
    sage: (O*(2, 1/2*a + 1/2)).is_principal()
    False
    sage: (O*(2, 1/2*a + 1/2))^3
    Fractional ideal (1/2*a - 3/2) of Number Field in a with defining polynomial x^2 + 23
"""

from sage.groups.abelian_gps.abelian_group import AbelianGroup_class

from sage.structure.sequence import Sequence

from sage.structure.element import MultiplicativeGroupElement

class ClassGroup(AbelianGroup_class):
    """
    The class group of a number field.
    """
    def __init__(self, invariants, names, number_field, gens):
        self.__number_field = number_field
        self.__gens = Sequence([FractionalIdealClass(x, self) for x in gens], immutable=True,
                               universe=self, check=False)
        AbelianGroup_class.__init__(self, len(invariants), invariants, names)

    def __call__(self, x):
        return FractionalIdealClass(self.__number_field.ideal(x), self)

    def _coerce_impl(self, x):
        return self(x)

    def gens(self):
        return self.__gens

    def ngens(self):
        return len(self.__gens)

    def gen(self, i=0):
        if i < 0 or i >= len(self.__gens):
            raise IndexError
        return self.__gens[i]

    def _repr_(self):
        return 'Class group of order %s with structure %s of %s'%(
            self.order(),
            self._group_notation(self.invariants()),
            self.number_field())

    def number_field(self):
        return self.__number_field


class FractionalIdealClass(MultiplicativeGroupElement):
    """
    A fractional ideal in a number field.

    EXAMPLES:
        sage: G = NumberField(x^2 + 23,'a').class_group(); G
        Class group of order 3 with structure C3 of Number Field in a with defining polynomial x^2 + 23
        sage: I = G.0; I
        Fractional ideal class (2, 1/2*a - 1/2) of Number Field in a with defining polynomial x^2 + 23
        sage: I*I
        Fractional ideal class (2, 1/2*a + 1/2) of Number Field in a with defining polynomial x^2 + 23
        sage: I*I*I
        Trivial principal fractional ideal class of Number Field in a with defining polynomial x^2 + 23
    """
    def __init__(self, ideal, class_group):
        self.__ideal = ideal
        MultiplicativeGroupElement.__init__(self, class_group)

    def _repr_(self):
        if self.is_principal():
            return 'Trivial principal fractional ideal class of %s'%self.__ideal.number_field()
        return 'Fractional ideal class %s of %s'%(self.__ideal._repr_short(), self.__ideal.number_field())

    def __cmp__(self, other):
        q = self.__ideal / other.__ideal
        if q.is_principal():
            return 0
        return cmp(self.__ideal, other.__ideal)

    def _mul_(self, other):
        return self.parent()((self.__ideal * other.__ideal).reduce_equiv())

    def is_principal(self):
        return self.__ideal.is_principal()

    def reduce(self):
        """
        Return representative for this ideal class that has been
        reduced using PARI's idealred.

        EXAMPLES:
            sage: k.<a> = NumberField(x^2 + 20072); G = k.class_group(); G
            Class group of order 76 with structure C38 x C2 of Number Field in a with defining polynomial x^2 + 20072
            sage: I = G.0; I
            Fractional ideal class (41, a + 10) of Number Field in a with defining polynomial x^2 + 20072
            sage: J = G(I.ideal()^5); J
            Fractional ideal class (115856201, 1/2*a + 40407883) of Number Field in a with defining polynomial x^2 + 20072
            sage: J.reduce()
            Fractional ideal class (57, 1/2*a + 44) of Number Field in a with defining polynomial x^2 + 20072
        """
        return self.parent()(self.__ideal.reduce_equiv())

    #def multiplicative_order(self):
    #    try:
    #        return self.__multiplicative_order
    #    except AttributeError:

    def ideal(self):
        """
        Return a representative ideal in this ideal class.
        """
        return self.__ideal

    def gens(self):
        """
        Return generators for a representative ideal in this
        ideal class.
        """
        return self.ideal().gens()
