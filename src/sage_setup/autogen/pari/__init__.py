from sage_setup.autogen.pari.generator import PariFunctionGenerator
import os

def rebuild(force=False):
    stamp = os.path.join(os.path.dirname(__file__), 'timestamp')
    desc = os.path.join(os.environ['SAGE_LOCAL'], 'share', 'pari', 'pari.desc')
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
