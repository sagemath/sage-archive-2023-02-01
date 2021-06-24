from sage.env import SAGE_ENV

from sage.features import Executable

class FourTi2Executable(Executable):
    r"""
    Feature for the 4ti2 executables.
    """
    def __init__(self, name):
        Executable.__init__(self,
                            name="4ti2-" + name,
                            executable=SAGE_ENV.get("FOURTITWO_" + name.upper(), None) or name,
                            spkg="4ti2")
