from . import main
main()

# Replace old documentation output by symbolic link to SAGE_DOC_OUTPUT
import os
from sage.env import SAGE_DOC, SAGE_DOC_OUTPUT
old_doc = os.path.join(SAGE_DOC, "output")

if not os.path.islink(old_doc):
    if os.path.exists(old_doc):
        from shutil import rmtree
        rmtree(old_doc)
    os.symlink(SAGE_DOC_OUTPUT, old_doc)
