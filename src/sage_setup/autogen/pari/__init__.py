from sage_setup.autogen.pari.generator import PariFunctionGenerator
from sage_setup.autogen.pari.parser import pari_share
import os

def rebuild(force=False):
    stamp = os.path.join(os.path.dirname(__file__), 'timestamp')
    desc = os.path.join(pari_share(), 'pari.desc')
    try:
        if not force and os.stat(stamp).st_mtime >= os.stat(desc).st_mtime:
            # No need to rebuild
            return
    except OSError:
        pass

    print("Generating PARI functions.")
    G = PariFunctionGenerator()
    G()

    open(stamp, 'w').close()
