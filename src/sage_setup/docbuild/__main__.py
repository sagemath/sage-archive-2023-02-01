from . import main
main()

# Remove old documentation output
import os
from sage.env import SAGE_DOC
old_doc = os.path.join(SAGE_DOC, "output")
if os.path.exists(old_doc):
    from shutil import rmtree
    rmtree(old_doc)
