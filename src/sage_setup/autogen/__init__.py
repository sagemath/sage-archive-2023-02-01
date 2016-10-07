import os
from sage.env import SAGE_SRC

def autogen_all():
    """
    Regenerate the automatically generated files of the Sage library.
    """
    from sage_setup.autogen import pari
    pari.rebuild()

    from sage_setup.autogen import interpreters
    interpreters.rebuild(os.path.join(SAGE_SRC, "sage", "ext", "interpreters"))
