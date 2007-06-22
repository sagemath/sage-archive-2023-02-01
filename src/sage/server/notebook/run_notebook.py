#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


"""
Server the SAGE Notebook.
"""

import getpass

##########################################################
# This actually serves up the notebook.
##########################################################

from sage.misc.misc import DOT_SAGE
from   sage.server.misc import print_open_msg
import os, shutil, socket

import notebook

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
             reset       = False,

             server_pool = None,
             ulimit      = None):
    r"""
    The SAGE Notebook.

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
        reset      -- (default: False) if True allows you to set the
                      admin password.  Use this if you forget your
                      admin password.
    \begin{verbatim}
    EXAMPLES:

    1. I want to run the SAGE notebook server on a remote machine
       and be the only person allowed to log in.  Type

                   notebook(address="")

       the first time you do this you'll be prompted to set
       an administrator password.  Use this to login.

    2. I just want to run the server locally and do not want to be
       bothered with SSL, accounts, etc., and I am the only
       user so I do not have to worry about somebody else using
       the notebook on localhost and deleting my files.  Use

             notebook(secure = False)

    3. I want to create a SAGE notebook server that is open to anybody
       in the world to create new accounts, etc.  To run the SAGE
       notebook publically (1) at a minimu run it from a chroot jail
       (see the SAGE install guide), and (2) use a command like

    notebook(secure=True, server_pool=['sage1@localhost'], ulimit='-v 500000')

       The secure option enables enccryption between all users and the
       notebook server.  The server_pool option specifies that
       worksheet processes run as a separate user.  The ulimit option
       restricts the memory available to each worksheet processes to
       500MB.  You will have to explicitly login as the admin user
       and click "Server" in the upper right home page to configure the
       server to allow users to create new accounts.



    INPUT:  (more advanced)

      NOTE: The values of these two properties default to what they were
     last time the notebook command was called.

        server_pool -- (default: None), if given, should be a list like
                      ['sage1@localhost', 'sage2@localhost'], where
                      you have setup ssh keys so that typing
                         ssh sage1@localhost
                      logs in without requiring a password, e.g., by typing
                      as the notebook server user
                          cd; ssh-keygen -t rsa
                      then put ~/.ssh/id_rsa.pub as the file .ssh/authorized_keys2.
                      Note -- you have to get the permissions of files
                      and directories just right -- do a web search
                      for more details.

        ulimit      -- (default: None -- leave as is), if given and server_pool is also given,
                      the worksheet processes are run with these constraints.
                      See the ulimit documentation. Common options include:
                           -f   The maximum size of files created by the shell
                           -t   The maximum amount of cpu time in seconds.
                           -u   The maximum number of processes available to a single user.
                           -v   The maximum amount of virtual memory available to the process.
                      Values are in 1024-byte increments, except for `-t', which is in seconds.
                      Example:  ulimit="-v 400000 -t 30"

    \end{verbatim}
    """

    if not os.path.exists(directory):
        os.makedirs(directory)
    port = int(port)
    conf = '%s/twistedconf.tac'%directory

    nb = notebook.load_notebook(directory)
    if reset or not nb.user_exists('admin'):
        if not reset:
            print '\n' + '-'*70
            print "Creating a new SAGE notebook, which will be stored in the directory"
            print os.path.abspath(directory)
        while True:
            print "Setting password for admin user."
            passwd = getpass.getpass("Enter new password: ")
            passwd2 = getpass.getpass("Retype new password: ")
            if passwd != passwd2:
                print "Sorry, passwords do not match."
            else:
                break
        if reset:
            nb.user('admin').set_password(passwd)
            print "Password changed for user 'admin'."
        else:
            nb.create_default_users(passwd)
            print "User admin created with the password you specified."

    if not server_pool is None:
        nb.set_server_pool(server_pool)

    if not ulimit is None:
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
####################################################################
# WARNING -- Do not edit this file!   It is autogenerated each time
# the notebook(...) command is executed.
####################################################################
import sage.server.notebook.notebook
sage.server.notebook.notebook.JSMATH=True
import sage.server.notebook.notebook as notebook
import sage.server.notebook.twist as twist
twist.notebook = notebook.load_notebook(%s)
twist.OPEN_MODE = %s
twist.init_updates()
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

realm = avatars.LoginSystem()
p = portal.Portal(realm)
password_checker = avatars.PasswordChecker()
p.registerChecker(password_checker)
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
"""%(notebook_opts, not secure, strport))


        config.close()

        ## Start up twisted
        print_open_msg(address, port, secure=secure)
        e = os.system('cd "%s" && sage -twistd -ny twistedconf.tac'%directory)
        if e == 256:
            raise socket.error


    for i in range(int(port_tries)+1):
        try:
            run(port + i)
        except socket.error:
            print "Port %s is already in use." #  Trying next port..."%port
        else:
            break

    return True










#######

##     if address != 'localhost':
##         if not secure:
##             print "*"*70
##             print "WARNING: Insecure notebook server listening on external address."
##             print "Unless you are running this via ssh port forwarding, you are"
##             print "**crazy**!  You should run the notebook with the option secure=True."
##             print "*"*70

##         if secure and not server_pool:
##             print "*"*70
##             print "You are running an ecrypted secure server, but without an external"
##             print "server pool of worksheet users.  This is a bad idea, since anybody"
##             print "can trivially ** kill ** your server at any time."
##             print "*"*70

