# -*- coding: utf-8 -*-
#
# Sage documentation build configuration file, created by
# sphinx-quickstart on Thu Aug 21 20:15:55 2008.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os
from sage.env import SAGE_DOC
sys.path.append(SAGE_DOC)
from common.conf import *

# settings for the intersphinx extension:

ref_src = os.path.join(SAGE_DOC, 'en', 'reference')
ref_out = os.path.join(SAGE_DOC, 'output', 'html', 'en', 'reference')
intersphinx_mapping[ref_out] = None

for doc in os.listdir(ref_src):
    if os.path.exists(os.path.join(ref_src, doc, 'index.rst')):
        intersphinx_mapping[os.path.join(ref_out, doc)] = None

# General information about the project.
project = u"Sage Reference Manual"
name = "reference"

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = project + " v"+release

html_short_title = u'Sage Reference v' + release

# Output file base name for HTML help builder.
htmlhelp_basename = name

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
  ('index', name + '.tex', u'Sage Reference Manual',
   u'The Sage Development Team', 'manual'),
]

latex_elements['preamble'] += r'''
% One-column index
\makeatletter
\renewenvironment{theindex}{
  \chapter*{\indexname}
  \markboth{\MakeUppercase\indexname}{\MakeUppercase\indexname}
  \setlength{\parskip}{0.1em}
  \relax
  \let\item\@idxitem
}{}
\makeatother
'''

#Ignore all .rst in the _sage subdirectory
exclude_trees = exclude_trees + ['_sage']

multidocs_is_master = True

# List of subdocs
multidocs_subdoc_list = [
    'algebras',
    'arithgroup',
    'calculus',
    'categories',
    'cmd',
    'coding',
    'coercion',
    'combinat',
    'constants',
    'cryptography',
    'databases',
    'finance',
    'finite_rings',
    'function_fields',
    'functions',
    'games',
    'geometry',
    'graphs',
    'groups',
    'hecke',
    'history_and_license',
    'homology',
    'interfaces',
    'lfunctions',
    'libs',
    'logic',
    'matrices',
    'misc',
    'modabvar',
    'modfrm',
    'modmisc',
    'modsym',
    'modules',
    'monoids',
    'notebook',
    'number_fields',
    'numerical',
    'padics',
    'parallel',
    'plane_curves',
    'plot3d',
    'plotting',
    'polynomial_rings',
    'power_series',
    'probability',
    'quadratic_forms',
    'quat_algebras',
    'rings',
    'rings_numerical',
    'rings_standard',
    'sat',
    'schemes',
    'semirings',
    'stats',
    'structure',
    'tensor'
    ]

# List of directories, relative to source directory, that shouldn't be
# searched for source files.
exclude_trees += multidocs_subdoc_list + [
    'sage', 'sagenb', 'options'
    ]
