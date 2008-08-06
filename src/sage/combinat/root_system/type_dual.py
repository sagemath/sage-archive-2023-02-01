from sage.combinat.root_system.cartan_type import CartanType_abstract

class CartanType(CartanType_abstract):
    r"""
    A class for dual Cartan types

    The dual of a cartan type is a cartan type with the same index
    set, but all arrows reversed in the Dynkin diagram (otherwise
    said, the Cartan matrix is transposed). It shares a lot of
    properties in common with its dual. In particular, the Weyl group
    is isomorphic to that of the dual as a Coxeter group.
    """

    def __init__(self, type):
        """
        INPUT:
           type -- a Cartan type

        EXAMPLES:
           sage: ct = CartanType(['F',4]).dual()
           sage: ct == loads(dumps(ct))
           True
        """
        self._dual = type

    def __repr__(self):
        """
        EXAMPLES:
           sage: ct = CartanType(['F', 4]).dual()
           sage: ct
           ['F', 4]^*
        """
        return self.dual().__repr__()+"^*"

    def __cmp__(self, other):
        """
        EXAMPLES:
            sage: ct1 = CartanType(['A',1],['B',2])
            sage: ct2 = CartanType(['B',2],['A',1])
            sage: ct3 = CartanType(['A',4])
            sage: ct = [ct1, ct2, ct3]
            sage: [[x.__cmp__(y) for x in ct] for y in ct]
            [[0, 1, -1], [-1, 0, -1], [1, 1, 0]]
            sage: sorted(ct)
            [['A', 4], A1xB2, B2xA1]
        """
        if other.__class__ != self.__class__:
            return cmp(self.__class__, other.__class__)
        return cmp(self._dual, other._dual)

    def dual(self):
        return self._dual

    def is_irreducible(self):
        return self._dual.is_irreducible()

    def is_finite(self):
        return self._dual.is_finite()

    def is_crystalographic(self):
        return self._dual.is_crystalographic

    def is_affine(self):
        return self._dual.is_affine()

    def rank(self):
        return self._dual.rank()

    def index_set(self):
        return self._dual.index_set()

    def dynkin_diagram(self):
        return self._dual.dynkin_diagram().dual()

#class ambient_space(AmbientSpace):
# todo?
