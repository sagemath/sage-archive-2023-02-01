"""nodoctest
"""

#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

"""
Helper functions dealing with the verification of user
"""

def build_msg(key, username, addr, port, secure):
    url_prefix = "https" if secure else "http"
    s  = "Hi %s!\n\n" % username
    s += """\
Thank you for registering for the SAGE notebook. To complete your registration,
copy and paste the following link into your browser:

%s://%s:%s/confirm?key=%s

You will be taken to a page which will confirm that you have indeed
registered.""" % (url_prefix, addr, port, key)
    return s

def build_password_msg(key, username, addr, port, secure):
    url_prefix = "https" if secure else "http"
    s  = "Hi %s!\n\n" % username
    s += """\
Your new password is %s

Sign in at %s://%s:%s/

Make sure to reset your password by going to Settings in the upper right bar.""" % (key, url_prefix, addr, port)
    return s

def make_key():
    from random import randint
    key = randint(0,2**128-1)
    return key
