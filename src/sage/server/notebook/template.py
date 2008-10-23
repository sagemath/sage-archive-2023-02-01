"""nodoctest
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

"""
HTML templating for the notebook

AUTHORS:
    -- Bobby Moretti
    -- Timothy Clemans
"""

from jinja import Environment, FileSystemLoader
from jinja.exceptions import TemplateNotFound
from os.path import join, exists, getmtime
from sage.misc.misc import SAGE_ROOT

env = Environment(loader=FileSystemLoader(SAGE_ROOT+"/devel/sage/sage/server/notebook/templates"))

class PageTemplate:
    def __init__(self, filename):
        self.__template = env.get_template(filename)

    def __call__(self, **kwds):
        return str(self.__template.render(kwds))

# Define variables for each template
G = globals()
templates = ['login', 'yes_no', 'registration', 'account_settings']
for name in templates:
    G[name] =  PageTemplate('%s.html'%name)
