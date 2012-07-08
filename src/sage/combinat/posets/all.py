from posets import Poset

from lattices import LatticePoset
from lattices import MeetSemilattice
from lattices import JoinSemilattice

from poset_examples import posets, Posets

###########################################################################
##### DEPRECATION WARNINGS ################################################
##### Added 28 April 2009 #################################################
###########################################################################

from sage.misc.superseded import deprecation

def BooleanLattice(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.BooleanLattice`` instead.

    TESTS::

        sage: BooleanLattice(3)
        doctest:...: DeprecationWarning: BooleanLattice is deprecated, use Posets.BooleanLattice instead!
        See http://trac.sagemath.org/10998 for details.
        Finite lattice containing 8 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("BooleanLattice", "BooleanLattice"))
    return Posets.BooleanLattice(*args, **kwds)

def ChainPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.ChainPoset`` instead.

    TESTS::

        sage: ChainPoset(3)
        doctest:...: DeprecationWarning: ChainPoset is deprecated, use Posets.ChainPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite lattice containing 3 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("ChainPoset","ChainPoset"))
    return Posets.ChainPoset(*args, **kwds)

def AntichainPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.AntichainPoset`` instead.

    TESTS::

        sage: AntichainPoset(3)
        doctest:...: DeprecationWarning: AntichainPoset is deprecated, use Posets.AntichainPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 3 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("AntichainPoset","AntichainPoset"))
    return Posets.AntichainPoset(*args, **kwds)

def PentagonPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.PentagonPoset`` instead.

    TESTS::

        sage: PentagonPoset()
        doctest:...: DeprecationWarning: PentagonPoset is deprecated, use Posets.PentagonPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite lattice containing 5 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("PentagonPoset","PentagonPoset"))
    return Posets.PentagonPoset(*args, **kwds)

def DiamondPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.DiamondPoset`` instead.

    TESTS::

        sage: DiamondPoset(3)
        doctest:...: DeprecationWarning: DiamondPoset is deprecated, use Posets.DiamondPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite lattice containing 3 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("DiamondPoset","DiamondPoset"))
    return Posets.DiamondPoset(*args, **kwds)

def PosetOfIntegerCompositions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.IntegerCompositions`` instead.

    TESTS::

        sage: PosetOfIntegerCompositions(3)
        doctest:...: DeprecationWarning: PosetOfIntegerCompositions is deprecated, use Posets.IntegerCompositions instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 4 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfIntegerCompositions","IntegerCompositions"))
    return Posets.IntegerCompositions(*args, **kwds)

def PosetOfIntegerPartitions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.IntegerPartitions`` instead.

    TESTS::

        sage: PosetOfIntegerPartitions(3)
        doctest:...: DeprecationWarning: PosetOfIntegerPartitions is deprecated, use Posets.IntegerPartitions instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 3 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfIntegerPartitions","IntegerPartitions"))
    return Posets.IntegerPartitions(*args, **kwds)

def PosetOfRestrictedIntegerPartitions(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.RestrictedIntegerPartitions`` instead.

    TESTS::

        sage: PosetOfRestrictedIntegerPartitions(3)
        doctest:...: DeprecationWarning: PosetOfRestrictedIntegerPartitions is deprecated, use Posets.RestrictedIntegerPartitions instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 3 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("PosetOfRestrictedIntegerPartitions","RestrictedIntegerPartitions"))
    return Posets.RestrictedIntegerPartitions(*args, **kwds)

def RandomPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.RandomPoset`` instead.

    TESTS::

        sage: RandomPoset(17,.15)
        doctest:...: DeprecationWarning: RandomPoset is deprecated, use Posets.RandomPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 17 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("RandomPoset","RandomPoset"))
    return Posets.RandomPoset(*args, **kwds)


def SymmetricGroupBruhatOrderPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.SymmetricGroupBruhatOrderPoset`` instead.

    TESTS::

        sage: SymmetricGroupBruhatOrderPoset(3)
        doctest:...: DeprecationWarning: SymmetricGroupBruhatOrderPoset is deprecated, use Posets.SymmetricGroupBruhatOrderPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 6 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("SymmetricGroupBruhatOrderPoset","SymmetricGroupBruhatOrderPoset"))
    return Posets.SymmetricGroupBruhatOrderPoset(*args, **kwds)

def SymmetricGroupWeakOrderPoset(*args, **kwds):
    r"""
    This function is deprecated and will be removed in a future
    version of Sage. Please use ``Posets.SymmetricGroupWeakOrderPoset`` instead.

    TESTS::

        sage: SymmetricGroupWeakOrderPoset(3)
        doctest:...: DeprecationWarning: SymmetricGroupWeakOrderPoset is deprecated, use Posets.SymmetricGroupWeakOrderPoset instead!
        See http://trac.sagemath.org/10998 for details.
        Finite poset containing 6 elements
    """
    deprecation(10998, "%s is deprecated, use Posets.%s instead!" % \
           ("SymmetricGroupWeakOrderPoset","SymmetricGroupWeakOrderPoset"))
    return Posets.SymmetricGroupWeakOrderPoset(*args, **kwds)
