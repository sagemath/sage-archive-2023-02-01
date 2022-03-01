r"""
Features for testing the presence of Singular
"""
from . import Executable


class Singular(Executable):
    r"""
    A :class:`~sage.features.Feature` describing the presence of the Singular executable.

    EXAMPLES::

        sage: from sage.features.singular import Singular
        sage: Singular().is_present()
        True
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.singular import Singular
            sage: isinstance(Singular(), Singular)
            True
        """
        Executable.__init__(self, "singular", "Singular",
                            spkg='singular')
