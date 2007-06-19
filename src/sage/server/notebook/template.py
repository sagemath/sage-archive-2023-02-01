#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

"""
HTML templating for the notebook

AUTHOR:
    -- Bobby Moretti
"""

from string import Template
from sage.misc.misc import SAGE_EXTCODE

import os

pjoin = os.path.join
path = pjoin(SAGE_EXTCODE, "notebook/templates")


class PageTemplate:
    def __init__(self, filename):
        file = open(filename, 'r')
        self.__template = Template(file.read())
        file.close()

    def __call__(self, **kwds):
        return self.__template.substitute(kwds)

# Define variables for each template
G = globals()
templates = ['login', 'yes_no', 'failed_login']
for name in templates:
    G[name + '_template'] =  PageTemplate(pjoin(path, '%s.template'%name))
