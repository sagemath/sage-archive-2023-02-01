r"""
Check for pynormaliz
"""
from . import PythonModule
from .join_feature import JoinFeature


class pynormaliz(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of the
    Python package ``PyNormaliz``.

    EXAMPLES::

        sage: from sage.features.normaliz import pynormaliz
        sage: pynormaliz().is_present()                    # optional - pynormaliz
        FeatureTestResult('pynormaliz', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.normaliz import pynormaliz
            sage: isinstance(pynormaliz(), pynormaliz)
            True
        """
        JoinFeature.__init__(self, 'pynormaliz',
                             [PythonModule('PyNormaliz', spkg="pynormaliz")])
