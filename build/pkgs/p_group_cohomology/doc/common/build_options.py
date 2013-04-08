###############################################
# Options for building the Sage documentation #
###############################################

import os
SAGE_DOC = os.path.realpath('.')
LANGUAGES = ['en']
SPHINXOPTS = ""
PAPER = ""
OMIT = ["introspect"]  # docs/dirs to omit when listing and building 'all'

if PAPER:
    PAPEROPTS = "-D latex_paper_size=" + PAPER
else:
    PAPEROPTS = ""

#Note that this needs to have the doctrees dir
ALLSPHINXOPTS   = SPHINXOPTS + " " + PAPEROPTS + " "
