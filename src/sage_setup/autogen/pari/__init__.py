import glob
from os.path import join, getmtime, exists

from sage_setup.autogen.pari.generator import PariFunctionGenerator
from sage_setup.autogen.pari.parser import pari_share

def rebuild(force=False):
    src_files = [join(pari_share(), 'pari.desc'),
                 join('sage', 'libs', 'pari', 'decl.pxi')] + \
                glob.glob(join('sage_setup', 'autogen', 'pari', '*.py'))
    gen_files = [join('sage', 'libs', 'pari', 'auto_gen.pxi')]

    if all(exists(f) for f in gen_files):
        src_mtime = max(getmtime(f) for f in src_files)
        gen_mtime = min(getmtime(f) for f in gen_files)

        if gen_mtime > src_mtime and not force:
            return

    G = PariFunctionGenerator()
    G()
