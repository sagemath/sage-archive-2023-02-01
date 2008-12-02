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
    -- Bobby Moretti (2007-07-18): initial version
    -- Timothy Clemans and Mike Hansen (2008-10-27): major update
"""

import jinja
import sage.misc.misc as misc

env = jinja.Environment(loader=jinja.FileSystemLoader(misc.SAGE_ROOT + '/devel/sage/sage/server/notebook/templates'))

def contained_in(container):
    def wrapped(env, context, value):
        return value in container
    return wrapped

env.tests['contained_in'] = contained_in

default_context = {'sitename': 'Sage Notebook'}

def template(filename, **context):
    """
    Returns a compiled template.
    """
    try:
        tmpl = env.get_template(filename)
    except jinja.exceptions.TemplateNotFound:
        return template('template_error.html')
    context.update(default_context)
    return str(tmpl.render(**context))