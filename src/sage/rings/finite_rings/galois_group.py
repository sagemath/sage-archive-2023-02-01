r"""
Galois groups of Finite Fields
"""

from sage.groups.galois_group import GaloisGroup as GaloisGroup_base

class GaloisGroup_GF(GaloisGroup_base):
    r"""
    The Galois group of a finite field.
    """
    def __init__(self, finite_field):
        r"""
        Create a Galois group.

        TESTS::

            sage: TestSuite(GF(9).galois_group()).run()
        """
        super().__init__(finite_field, algorithm=None, names=None, gc_numbering=False)

    def order(self, algorithm=None, recompute=False):
        r"""
        Return the order of this Galois group, which is just the degree of the extension since finite fields are Galois.

        EXAMPLES::

            sage: GF(9).galois_group().order()
            2
        """
        return self._field.degree()
