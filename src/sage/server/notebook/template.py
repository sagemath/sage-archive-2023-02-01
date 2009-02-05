"""
HTML templating for the notebook

AUTHORS:
    -- Bobby Moretti (2007-07-18): initial version
    -- Timothy Clemans and Mike Hansen (2008-10-27): major update

"""
#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################
import jinja
import sage.misc.misc
from sage.version import version

TEMPLATE_PATH = sage.misc.misc.SAGE_ROOT + '/devel/sage/sage/server/notebook/templates'
env = jinja.Environment(loader=jinja.FileSystemLoader(TEMPLATE_PATH))

def contained_in(container):
    """
    Returns a function which takes in an environment, context, and value
    and returns True if that value is in the container and False
    otherwise.  This is registered and used as a test in the templates.

    EXAMPLES:
        sage: from sage.server.notebook.template import contained_in
        sage: f = contained_in([1,2,3])
        sage: f(None, None, 2)
        True
        sage: f(None, None, 4)
        False
    """
    def wrapped(env, context, value):
        return value in container
    return wrapped

env.tests['contained_in'] = contained_in

#A dictionary containing the default context
#The values in this dictionary will be updated
#by the
default_context = {'sitename': 'Sage Notebook',
                   'sage_version': version}

def template(filename, **user_context):
    """
    Returns a rendered template as a string.

    INPUT:
        filename -- the filename of the template relative to
                    $SAGE_ROOT/devel/sage/sage/server/notebook/templates

    EXAMPLES:
        sage: from sage.server.notebook.template import template
        sage: s = template('yes_no.html'); type(s)
        <type 'str'>
        sage: 'Yes' in s
        True

        sage: from sage.server.notebook.template import template
        sage: u = unicode('Are Gröbner bases awesome?','utf-8')
        sage: s = template('yes_no.html',message=u)
        sage: 'Gr\xc3\xb6bner' in s
        True
    """
    try:
        tmpl = env.get_template(filename)
    except jinja.exceptions.TemplateNotFound:
        return template('template_error.html', template=filename)
    context = dict(default_context)
    context.update(user_context)
    r = tmpl.render(**context)
    return r.encode('utf-8')
