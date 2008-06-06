"""nodoctest
"""

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
from   sage.server.misc import print_open_msg, find_next_available_port
import os, shutil, socket, pexpect

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
             directory   = None,
             port        = 8000,
             address     = 'localhost',
             port_tries  = 50,
             secure      = False,
             reset       = False,
             accounts    = False,
             require_login = True,

             server_pool = None,
             ulimit      = None,

             timeout     = 0,

             open_viewer = True,

             sagetex_path = "",
             start_path = "",
             fork = False,
             quiet = False):

    if directory is None:
        directory = '%s/sage_notebook'%DOT_SAGE
    else:
        if isinstance(directory, basestring) and len(directory) > 0 and directory[-1] == "/":
            directory = directory[:-1]

    if not os.path.exists(directory):
        os.makedirs(directory)

    if not quiet:
        print "The notebook files are stored in:", directory
    # First change to the directory that contains the notebook directory
    wd = os.path.split(directory)
    if wd[0]: os.chdir(wd[0])
    directory = wd[1]


    port = int(port)
    conf = '%s/twistedconf.tac'%directory

    if not secure and address != 'localhost':
        print '*'*70
        print "WARNING: Running the notebook insecurely not on localhost is dangerous"
        print "because its possible for people to sniff passwords and gain access to"
        print "your account. Make sure you know what you are doing."
        print '*'*70

    nb = notebook.load_notebook(directory)

    nb.conf()['idle_timeout'] = int(timeout)

    if nb.user_exists('root') and not nb.user_exists('admin'):
        # This is here only for backward compatibility with one
        # version of the notebook.
        s = nb.create_user_with_same_password('admin', 'root')
        # It would be a security risk to leave an escalated account around.

    if not nb.user_exists('admin'):
        reset = True

    if reset:
        passwd = get_admin_passwd()
        if reset:
            nb.user('admin').set_password(passwd)
            print "Password changed for user 'admin'."
        else:
            nb.create_default_users(passwd)
            print "User admin created with the password you specified."
            print "\n\n"
            print "*"*70
            print "\n"
            if secure:
                print "Login to the SAGE notebook as admin with the password you specified above."
        #nb.del_user('root')

    if not server_pool is None:
        nb.set_server_pool(server_pool)

    if not ulimit is None:
        nb.set_ulimit(ulimit)

    nb.set_accounts(accounts)

    if os.path.exists('%s/nb-older-backup.sobj'%directory):
        nb._migrate_worksheets()
        os.unlink('%s/nb-older-backup.sobj'%directory)
        print "Updating to new format complete."

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
            strport = 'tls:%s:interface=%s:privateKey=%s:certKey=%s'%(port, address, private_pem, public_pem)
        else:
            strport = 'tcp:%s:interface=%s'%(port, address)

        notebook_opts = '"%s",address="%s",port=%s,secure=%s' % (os.path.abspath(directory),
                address, port, secure)

        if open_viewer:
            if require_login:
                start_path = "'/?startup_token=%s' % startup_token"
            else:
                start_path = "'/'"
            open_page = "from sage.server.misc import open_page; open_page('%s', %s, %s, %s)"%(address, port, secure, start_path)
        else:
            open_page = ''

        config = open(conf, 'w')

        config.write("""
####################################################################
# WARNING -- Do not edit this file!   It is autogenerated each time
# the notebook(...) command is executed.
####################################################################
from twisted.internet import reactor

# Now set things up and start the notebook
import sage.server.notebook.notebook
sage.server.notebook.notebook.JSMATH=True
import sage.server.notebook.notebook as notebook
import sage.server.notebook.twist as twist
twist.notebook = notebook.load_notebook(%s)
twist.SAGETEX_PATH = "%s"
twist.OPEN_MODE = %s
twist.SID_COOKIE = str(hash("%s"))
twist.init_updates()
import sage.server.notebook.worksheet as worksheet
worksheet.init_sage_prestart(twist.notebook.get_server(), twist.notebook.get_ulimit())

import signal, sys, random
def save_notebook():
    print "Saving notebook..."
    twist.notebook.save()
    reactor.stop()
    print "Notebook cleanly saved."

def my_sigint(x, n):
    save_notebook()
    signal.signal(signal.SIGINT, signal.SIG_DFL)


signal.signal(signal.SIGINT, my_sigint)

## Disable client-side certificate request for gnutls
try:
    import gnutls.connection
    gnutls.connection.CERT_REQUEST = 0
except OSError:
    print "Note: GNUTLS not available."


## Authentication framework (ported from Knooboo)
from twisted.web2 import log, server, channel
from twisted.cred import portal, checkers, credentials
import sage.server.notebook.guard as guard
import sage.server.notebook.avatars as avatars

from twisted.cred import portal

realm = avatars.LoginSystem()
p = portal.Portal(realm)
startup_token = '%%x' %% random.randint(0, 2**128)
startup_checker = avatars.OneTimeTokenChecker()
startup_checker.token = startup_token
p.registerChecker(startup_checker)
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
%s
s.setServiceParent(application)

reactor.addSystemEventTrigger('before', 'shutdown', save_notebook)

"""%(notebook_opts, sagetex_path, not require_login,
     os.path.abspath(directory), strport, open_page))


        config.close()

        ## Start up twisted
        if not quiet:
            print_open_msg('localhost' if not address else address, port, secure=secure)
        if secure and not quiet:
            print "There is an admin account.  If you do not remember the password,"
            print "quit the notebook and type notebook(reset=True)."
        cmd = 'sage -twistd --pidfile="%s"/twistd.pid -ny "%s"/twistedconf.tac'%(directory, directory)
        if fork:
            return pexpect.spawn(cmd)
        else:
            e = os.system(cmd)
        if e == 256:
            raise socket.error
        return True
        # end of inner function run

    if address != 'localhost' and not secure:
            print "*"*70
            print "WARNING: Insecure notebook server listening on external address."
            print "Unless you are running this via ssh port forwarding, you are"
            print "**crazy**!  You should run the notebook with the option secure=True."
            print "*"*70

    port = find_next_available_port(port, port_tries)
    #if open_viewer:
    #    open_page(address, port, secure, pause=PAUSE)
    if open_viewer:
        "Open viewer automatically isn't fully implemented.  You have to manually open your web browser to the above URL."
    return run(port)






#######


def get_admin_passwd():
    print "\n"*2
    print "Please choose a new password for the SAGE Notebook 'admin' user."
    print "Do _not_ choose a stupid password, since anybody who could guess your password"
    print "and connect to your machine could access or delete your files."
    print "NOTE: Only the md5 hash of the password you type is stored by SAGE."
    print "You can change your password by typing notebook(reset=True)."
    print "\n"*2
    while True:
        passwd = getpass.getpass("Enter new password: ")
        if len(passwd) < 6:
            print "That password is way too short. Enter a password with at least 6 characters."
            continue
        passwd2 = getpass.getpass("Retype new password: ")
        if passwd != passwd2:
            print "Sorry, passwords do not match."
        else:
            break

    print "Please login to the notebook with the username 'admin' and the above password."
    return passwd
