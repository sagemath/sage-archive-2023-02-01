"""
This file contains test functions that can be used to search
bugs by testing random objects.

As an examples: if a lattice is distributive, then it must be also
modular, and if a lattice is upper semimodular, then the dual
of the lattice must be lower semimodular.
"""
def test_finite_lattice(L):
    """
    Test several functions on a given finite lattice.

    EXAMPLES::

        sage: from sage.tests.test import test_finite_lattice
        sage: L = Posets.RandomLattice(10, 0.98)
        sage: test_finite_lattice(L) is None
        True
    """
    if L.is_distributive() and not L.is_modular():
        raise ValueError("bug: distributive lattice should be modular")
    if L.is_upper_semimodular() and not L.dual().is_lower_semimodular():
        raise ValueError("bug: dual of upper semimodular lattice should be lower semimodular")
