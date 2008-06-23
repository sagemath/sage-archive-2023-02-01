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

def registration_page_template(error=None, input=None):
    keywords = dict([(i, '') for i in ['error', 'username', 'username_error', 'password_error', 'confirm_pass_error', 'email', 'email_error']])
    def error_html(msg):
        return '<p><span class="error">Error:</span> ' + msg + '</p>'

    if error:
        error_msgs = {'error_msg': '<h2 class="error_found">Error%s found</h2>' % 's' if len(error) > 1 else '',
                      'username_taken': error_html("Username taken"),
                      'username_invalid': error_html("Invalid username"),
                      'username_missing': error_html("No username given"),
                      'password_invalid': error_html("Bad password"),
                      'password_missing': error_html("No password given"),
                      'passwords_dont_match': error_html("Passwords didn't match"),
                      'retype_password_missing': error_html("Passwords didn't match"),
                      'email_missing': error_html('No email address given'),
                      'email_invalid': error_html('Invalid email address')}

        keywords['error'] = error_msgs['error_msg']

        for i in error:
            if i != 'passwords_dont_match' and i != 'retype_password_missing':
                keywords[i.split('_')[0] + '_error'] = error_msgs[i]
            else:
                keywords['confirm_pass_error'] = error_msgs[i]

    file = open(pjoin(path, 'register.template'), 'r')
    if input:
        if 'password' in input:
            del input['password']
        if 'retype_password' in input:
            del input['retype_password']
        keywords.update(input)
    template = Template(file.read()).substitute(keywords)
    file.close()
    return template
