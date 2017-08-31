"""
This file contains test functions that can be used to search
bugs by testing random finite posets and lattices.

As an examples: if a lattice is distributive, then it must be also
modular, and if a poset is ranked, then the dual poset must also
be ranked.

TODO: Currently only finite lattices have test function. Add some to
general posets too.
"""
def test_finite_lattice(L):
    """
    Test several functions on a given finite lattice.

    EXAMPLES::

        sage: from sage.tests.finite_poset import test_finite_lattice
        sage: L = Posets.RandomLattice(10, 0.98)
        sage: test_finite_lattice(L) is None
        True
    """
    if L.cardinality() < 4:
        # Special cases should be tested in specific TESTS-sections.
        return None
    if L.is_distributive() and not L.is_modular():
        raise ValueError("bug: distributive lattice should be modular")
    if L.is_upper_semimodular() and not L.dual().is_lower_semimodular():
        raise ValueError("bug: dual of upper semimodular lattice should be lower semimodular")
