from . import main
main()

# Remove old documentation
import os
from sage.env import SAGE_DOC_SRC, SAGE_DOC
old_doc = os.path.join(SAGE_DOC_SRC, "output")

if os.path.exists(old_doc):
    try:
        os.unlink(old_doc)
    except OSError:
        from shutil import rmtree
        rmtree(old_doc)

# Add symlink ouput -> . inside SAGE_DOC (for backwards compatibility)
try:
    os.symlink(".", os.path.join(SAGE_DOC, "output"))
except OSError:
    pass
