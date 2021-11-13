r"""
Check for working interpreter interfaces
"""
import importlib

from . import Feature, FeatureTestResult, PythonModule


class InterfaceFeature(Feature):

    @staticmethod
    def __classcall__(cls, name, module, description=None):
        if isinstance(module, str):
            module = PythonModule(module)
        return Feature.__classcall__(cls, name, module, description)

    def __init__(self, name, module, description):
        super().__init__(name, description=description)
        self.module = module

    def _is_present(self):
        result = self.module.is_present()
        if not result:
            return result
        m = importlib.import_module(self.module.name)
        try:
            interface = getattr(m, self.name)
            interface('2+3')
            return FeatureTestResult(self, True)
        except Exception as exception:
            return FeatureTestResult(self, False,
                                     reason=f"Interface {interface} is not functional: {exception}")

# The following are provided by external software only (no SPKG)

class Magma(InterfaceFeature):
    r"""
    A :class:`sage.features.Feature` describing whether :class:`sage.interfaces.magma.Magma`
    is present and functional.

    EXAMPLES::

        sage: from sage.features.interfaces import Magma
        sage: Magma().is_present()  # random
        FeatureTestResult('jupymake', False)
    """

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'magma', 'sage.interfaces.magma')


class Matlab(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'matlab', 'sage.interfaces.matlab')


class Mathematica(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'mathematica', 'sage.interfaces.mathematica')


class Maple(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'maple', 'sage.interfaces.maple')


class Macaulay2(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'macaulay2', 'sage.interfaces.macaulay2')


class Octave(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'octave', 'sage.interfaces.octave')


class Scilab(InterfaceFeature):

    @staticmethod
    def __classcall__(cls):
        return InterfaceFeature(cls, 'scilab', 'sage.interfaces.scilab')
