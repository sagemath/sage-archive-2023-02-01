#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


"""
Server the SAGE Notebook.
"""

##########################################################
# This actually serves up the notebook.
##########################################################

from sage.misc.misc import DOT_SAGE
from   sage.server.misc import print_open_msg
import os, shutil, socket

conf_path       = os.path.join(DOT_SAGE, 'notebook')

private_pem = conf_path + '/private.pem'
public_pem  = conf_path + '/public.pem'

def notebook_setup(self=None):
    if not os.path.exists(conf_path):
        os.makedirs(conf_path)
    print "Using dsage certificates."
    dsage = os.path.join(DOT_SAGE, 'dsage')
    import sage.dsage.all
    sage.dsage.all.dsage.setup()
    shutil.copyfile(dsage + '/cacert.pem', private_pem)
    shutil.copyfile(dsage + '/pubcert.pem', public_pem)
    print "Successfully configured notebook."

def notebook_twisted(self,
             directory   = 'sage_notebook',
             port        = 8000,
             address     = 'localhost',
             port_tries  = 0,
             secure      = True,
             server_pool = None,
             ulimit      = None):
    r"""
    Experimental twisted version of the SAGE Notebook.

    INPUT:
        directory  -- (default: 'sage_notebook') directory that contains
                      the SAGE notebook files
        port       -- (default: 8000), port to serve the notebook on
        address    -- (default: 'localhost'), address to listen on
        port_tries -- (default: 0), number of additional ports to try if the
                      first one doesn't work (*not* implemented)
        secure     -- (default: True) if True use https so all
                      communication, e.g., logins and passwords,
                      between web browsers and the SAGE notebook is
                      encrypted (via GNU TLS).
    ADVANCED OPTIONS:
        server_pool -- (default: None), if given, should be a list like
                      ['sage1@localhost', 'sage2@localhost'], where
                      you have setup ssh keys so that typing
                         ssh sage1@localhost
                      logs in without requiring a password, e.g., by typing
                      as the notebook server user
                          cd; ssh-keygen -t rsa
                      then putting ~/.ssh/id_rsa.pub as the file .ssh/authorized_keys2.
        ulimit      -- (default: None -- leave as is), if given and server_pool is also given,
                      the worksheet processes are run with these constraints.
                      See the ulimit documentation. Common options include:
                           -f   The maximum size of files created by the shell
                           -t   The maximum amount of cpu time in seconds.
                           -u   The maximum number of processes available to a single user.
                           -v   The maximum amount of virtual memory available to the process.
                      Values are in 1024-byte increments, except for `-t', which is in seconds.
                      Example:  ulimit="-v 400000 -t 30"
    """
    if not os.path.exists(directory):
        os.makedirs(directory)
    port = int(port)
    conf = '%s/twistedconf.py'%directory

    # We load the notebook to make sure it is created with the
    # given options, then delete it.  The notebook is later
    # loaded by the *other* Twisted process below.
    if not server_pool is None or not ulimit is None:
        from sage.server.notebook.notebook import load_notebook
        nb = load_notebook(directory)
        nb.set_server_pool(server_pool)
        nb.set_ulimit(ulimit)
        nb.save()
        del nb

    def run(port):
        ## Create the config file
        if secure:
            if not os.path.exists(private_pem) or not os.path.exists(public_pem):
                print "In order to use an SECURE encrypted notebook, you must first run notebook.setup()."
                print "Now running notebook.setup()"
                notebook_setup()
            if not os.path.exists(private_pem) or not os.path.exists(public_pem):
                print "Failed to setup notebook.  Please try notebook.setup() again manually."
            strport = 'tls:%s:privateKey=%s:certKey=%s'%(port, private_pem, public_pem)
        else:
            strport = 'tcp:%s'%port

        notebook_opts = '"%s",address="%s",port=%s,secure=%s' % (os.path.abspath(directory),
                address, port, secure)
        config = open(conf, 'w')
        config.write("""
import sage.server.notebook.notebook
sage.server.notebook.notebook.JSMATH=True
import sage.server.notebook.notebook as notebook
import sage.server.notebook.twist as twist
twist.notebook = notebook.load_notebook(%s)
import sage.server.notebook.worksheet as worksheet
worksheet.init_sage_prestart(twist.notebook.get_server(), twist.notebook.get_ulimit())

import signal, sys
def my_sigint(x, n):
    twist.notebook.save()
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    print "(Notebook cleanly saved. Press control-C again to exit.)"

signal.signal(signal.SIGINT, my_sigint)

## Use Knoboo's authentication framework
from twisted.web2 import log, server, channel
from twisted.cred import portal, checkers, credentials
import sage.server.notebook.guard as guard
import sage.server.notebook.avatars as avatars

from twisted.cred import portal

password_file = 'passwords.txt'
realm = avatars.LoginSystem(password_file)
p = portal.Portal(realm)
# p.registerChecker(avatars.PasswordDataBaseChecker(DBCONNECTION))
p.registerChecker(avatars.PasswordFileChecker(password_file))
# p.registerChecker(checkers.AllowAnonymousAccess(), credentials.IAnonymous)
p.registerChecker(checkers.AllowAnonymousAccess())
rsrc = guard.MySessionWrapper(p)
log.DefaultCommonAccessLoggingObserver().start()
site = server.Site(rsrc)
factory = channel.HTTPFactory(site)

from twisted.web2 import channel
from twisted.application import service, strports
application = service.Application("SAGE Notebook")
s = strports.service('%s', factory)
s.setServiceParent(application)
"""%(notebook_opts, strport))


        config.close()

        ## Start up twisted
        print_open_msg(address, port, secure=secure)
        e = os.system('cd "%s" && sage -twistd -ny twistedconf.py'%directory)
        if e == 256:
            raise socket.error


    for i in range(int(port_tries)+1):
        try:
            run(port + i)
        except socket.error:
            print "Port %s is already in use.  Trying next port..."%port
        else:
            break

    return True
