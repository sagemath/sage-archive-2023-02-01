############################################################
# Backwards compatible unpickling
############################################################

def unpickle_QuaternionAlgebra_v0(*key):
    """
    The 0th version of pickling for quaternion algebras.

    EXAMPLES::

        sage: Q = QuaternionAlgebra(-5,-19)
        sage: f, t = Q.__reduce__()
        sage: import sage.algebras.quaternion_algebra
        sage: sage.algebras.quaternion_algebra.unpickle_QuaternionAlgebra_v0(*t)
        Quaternion Algebra (-5, -19) with base ring Rational Field
    """
    from quatalg.quaternion_algebra import QuaternionAlgebra
    return QuaternionAlgebra(*key)

