###############################################
# Options for building the Sage documentation #
###############################################

import os, re
SAGE_DOC = os.environ['SAGE_DOC']
LANGUAGES = [d for d in os.listdir(SAGE_DOC) if re.match('^[a-z][a-z]$', d)]
SPHINXOPTS = ""
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = "-D latex_paper_size=" + PAPER
else:
    PAPEROPTS = ""

#Note that this needs to have the doctrees dir
ALLSPHINXOPTS   = SPHINXOPTS + " " + PAPEROPTS + " "
