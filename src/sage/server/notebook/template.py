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
templates = ['login', 'yes_no', 'failed_login', 'register']
for name in templates:
    G[name + '_template'] =  PageTemplate(pjoin(path, '%s.template'%name))

def login_page_template(accounts, default_user, is_username_error=False, is_password_error=False, welcome=None):
    if accounts:
        reg = "<a href='/register'><b>Sign up for a new SAGE Notebook account</b></a>"
    else:
        reg = ""
    if is_username_error:
        u_e = '<tr><td align=right><span style="color:red">Error:</span></td><td>Username is not in the system</td></tr>'
    else:
        u_e = ''
    if is_password_error:
        p_e = '<tr><td align=right><span style="color:red">Error:</span></td><td>Wrong password</td></tr>'
    else:
        p_e = ''
    if welcome:
        welcome = '<h2>Congratulations %s! You can now sign into the Sage Notebook.</h2>' % welcome
    else:
        welcome = ''
    return login_template(register = reg, default=default_user, username_error=u_e, password_error=p_e, welcome=welcome)

def registration_page_template(error=None):
    if error:
        def error_html(msg):
            return '<p><span class="error">Error:</span> ' + msg + '</p>'

        error_msg = '<h2 class="error_found">Error found</h2>'
        username_error = ''
        password_error = ''
        confirm_pass_error = ''
        email_error = ''
        if error == 'username_taken':
            username_error = error_html("The username given is not available.")
        if error == 'username_invalid':
            username_error = error_html("The username given contains characters that are not allowed.")
        if error == 'username_missing':
            username_error = error_html("There was no username given.")
        if error == 'password_too_short':
            password_error = error_html("The password given was too short.")
        if error == 'password_missing':
            password_error = error_html("There was no password given.")
        if error == 'email_invalid':
            email_error = error_html('The e-mail address given is invalid.')
        return register_template(error=error_msg, username_error=username_error, password_error=password_error, confirm_pass_error=confirm_pass_error, email_error=email_error)
    else:
        return register_template(error='', username_error='', password_error='', confirm_pass_error='', email_error='')
