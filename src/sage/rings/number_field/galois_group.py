"""
Galois Groups of Number Fields

TESTS:

Standard test of pickleability:
    sage: G = NumberField(x^3 + 2, 'alpha').galois_group(); G
    Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in alpha with defining polynomial x^3 + 2
    sage: G == loads(dumps(G))
    True
"""

from sage.structure.sage_object import SageObject

class GaloisGroup(SageObject):
    r"""
    The Galois group of a number field.

    This is just a fairly minimal object at present.  To get the
    underlying group, do \code{G.group()}, and to get the
    corresponding number field do \code{G.number_field()}.  Galois
    groups are mainly useful in Sage right now for getting their
    structure and order, but not much more.  Of course much more
    general functionality is planned.

    EXAMPLES:
        sage: K = QQ[2^(1/3)]
        sage: G = K.galois_group(); G
        Galois group PARI group [6, -1, 2, "S3"] of degree 3 of the Number Field in a with defining polynomial x^3 - 2
        sage: G.order()
        6
        sage: G.group()
        PARI group [6, -1, 2, "S3"] of degree 3
        sage: G.number_field()
        Number Field in a with defining polynomial x^3 - 2
    """
    def __init__(self, group, number_field):
        """
        Create a Galois group.

        EXAMPLES:
            sage: NumberField([x^2 + 1, x^2 + 2],'a').galois_group()
            Galois group PARI group [4, 1, 2, "E(4) = 2[x]2"] of degree 4 of the Number Field in a0 with defining polynomial x^2 + 1 over its base field
        """
        self.__group = group
        self.__number_field = number_field

    def __cmp__(self, other):
        """
        Compare two number field Galois groups.  First the number
        fields are compared, then the Galois groups if the number
        fields are equal.  (Of course, if the number fields are the
        same, the Galois groups are automatically equal.)

        EXAMPLES:
            sage: G = NumberField(x^3 + 2, 'alpha').galois_group()
            sage: H = QQ[sqrt(2)].galois_group()
            sage: cmp(G,H)
            -1
            sage: H == H
            True
            sage: G == G
            True
        """
        if not isinstance(other, GaloisGroup):
            return cmp(type(self), type(other))
        return cmp( (self.__number_field, self.__group),
                    (other.__number_field, other.__group) )

    def __repr__(self):
        """
        Display print representation of a Galois group.

        EXAMPLES:
            sage: G = NumberField(x^4 + 2*x + 2, 'a').galois_group()
            sage: G.__repr__()
            'Galois group PARI group [24, -1, 5, "S4"] of degree 4 of the Number Field in a with defining polynomial x^4 + 2*x + 2'
        """
        return "Galois group %s of the %s"%(
            self.__group, self.__number_field)

    def group(self):
        """
        Return the underlying abstract group.

        EXAMPLES:
            sage: G = NumberField(x^3 + 2*x + 2, 'theta').galois_group()
            sage: H = G.group(); H
            PARI group [6, -1, 2, "S3"] of degree 3
            sage: P = H.permutation_group(); P  # optional -- requires Gap optional databases
            Transitive group number 2 of degree 3
            sage: list(P)                       # optional
            [(), (2,3), (1,2), (1,2,3), (1,3,2), (1,3)]
        """
        return self.__group

    def order(self):
        """
        Return the order of this Galois group.

        EXAMPLES:
            sage: G = NumberField(x^5 + 2, 'theta_1').galois_group(); G
            Galois group PARI group [20, -1, 3, "F(5) = 5:4"] of degree 5 of the Number Field in theta_1 with defining polynomial x^5 + 2
            sage: G.order()
            20
        """
        return self.__group.order()

    def number_field(self):
        """
        Return the number field of which this is the Galois group.

        EXAMPLES:
            sage: G = NumberField(x^6 + 2, 't').galois_group(); G
            Galois group PARI group [12, -1, 3, "D(6) = S(3)[x]2"] of degree 6 of the Number Field in t with defining polynomial x^6 + 2
            sage: G.number_field()
            Number Field in t with defining polynomial x^6 + 2
        """
        return self.__number_field
