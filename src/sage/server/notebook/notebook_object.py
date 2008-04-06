"""nodoctest
"""
#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################


"""
Notebook control object

This is used for configuring and starting the SAGE notebook server.
"""

import time, os, shutil, signal, tempfile

import notebook as _notebook

import run_notebook

class NotebookObject:
    r"""
    Start the SAGE Notebook server.

    INPUT:
        directory  -- directory that contains the SAGE notebook files;
                      The default is .sage/sage_notebook, in your home directory.
        port       -- (default: 8000), port to serve the notebook on
        address    -- (default: 'localhost'), address to listen on
        port_tries -- (default: 0), number of additional ports to try if the
                      first one doesn't work (*not* implemented)
        secure     -- (default: True) if True use https so all
                      communication, e.g., logins and passwords,
                      between web browsers and the SAGE notebook is
                      encrypted (via GNU TLS).  *Highly recommended!*
        reset      -- (default: False) if True allows you to set the
                      admin password.  Use this if you forget your
                      admin password.
        accounts   -- (default: False) if True, any visitor to the website
                      will be able to create a new account.  If False,
                      only the admin can create accounts (currently, this
                      can only be done by running with accounts=True for
                      a few minutes, or on the command line with, e.g.,
                          nb = load('sage_notebook/nb.sobj')
                          nb.set_accounts(True)
                          nb.add_user("username", "password", "email@place", "user")
                          nb.save()
        open_viewer -- (default: True) whether to pop up a web browser.
                      You can override the default browser by setting
                      the SAGE_BROWSER environment variable, e.g., by putting
                         export SAGE_BROWSER="firefox"
                      in the file .bashrc in your home directory.
        timeout    -- (default: 0) seconds until idle worksheet sessions
                      automatically timeout, i.e., the corresponding
                      Sage session terminates.  0 means 'never timeout'.
        server_pool -- list;   The server_pool option specifies that worksheet processes run
                      as a separate user (chosen from the list in the server_pool -- see below).

    \begin{verbatim}

    NOTE: If you have problems with the server certificate hostname not
    matching, do \code{notebook.setup()}.

    EXAMPLES:

    1. I just want to run the SAGE notebook.  Type

             notebook()

    2. I want to run the SAGE notebook server on a remote machine
       and be the only person allowed to log in.  Type

         notebook(address="address.of.remote.machine")

       the first time you do this you'll be prompted to set
       an administrator password.  Use this to login.
       NOTE: You may have to run notebook.setup() again and change
       the hostname.

    3. I just want to run the server locally on my laptop at a coffee
       shop with no wifi and do not want to be bothered with SSL,
       accounts, etc., and I am *absolutely certain* I am the only
       user logged into my laptop so I do not have to worry about
       somebody else using the notebook on localhost and deleting my
       files.  Use

                  notebook(secure=False)

    4. I want to create a SAGE notebook server that is open to anybody
       in the world to create new accounts, etc.  To run the SAGE
       notebook publically (1) at a minimu run it from a chroot jail
       (see the SAGE install guide), and (2) use a command like

         notebook(address="address.of.remote.machine",
            server_pool=['sage1@localhost'], ulimit='-v 500000', accounts=True)

       The server_pool option specifies that worksheet processes run
       as a separate user.  The ulimit option restricts the memory
       available to each worksheet processes to 500MB.  See help on
       the accounts option above.

       Be sure to make that the sage_notebook/nb.sobj and contents
       of sage_notebook/backups is chmod og-rwx, i.e., only readable
       by the notebook process, since otherwise any user can read
       nb.sobj, which contains user email addresses and account
       information (password are stored hashed, so less worries there).


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
    def __call__(self, *args, **kwds):
        return self.notebook(*args, **kwds)

    notebook = run_notebook.notebook_twisted
    setup    = run_notebook.notebook_setup

notebook = NotebookObject()


def inotebook(*args, **kwds):
    """
    Exactly the same as notebook(...) but with secure=False.
    """
    kwds['secure'] = False
    notebook(*args, **kwds)


def test_notebook(admin_passwd, directory=None, port=8050, address='localhost'):
    """
    This function is used to test notebook server functions.

    EXAMPLE:
        sage: from sage.server.notebook.notebook_object import test_notebook
        sage: passwd = str(randint(1,1<<128))
        sage: nb = test_notebook(passwd, address='localhost', port=8060)
        sage: import urllib
        sage: h = urllib.urlopen('https://localhost:8060')
        sage: homepage = h.read()
        sage: h.close()
        sage: 'html' in homepage
        True
        sage: nb.dispose()
        """
    import socket

    if directory is None:
        directory = tmp_dir = tempfile.mkdtemp()
    else:
        tmp_dir = None

    if not os.path.exists(directory):
        os.makedirs(directory)

    nb = _notebook.load_notebook(directory)
    nb.set_accounts(True)
    nb.add_user('admin', admin_passwd, '')
    nb.set_accounts(False)
    nb.save()

    p = notebook(directory=directory, port=port, address=address, open_viewer=False, fork=True, quiet=True)
    p.expect("Starting factory")
    def dispose():
        p.send('\x03') # control-C
        p.expect("Press control-C again to exit")
        p.send('\x03')
        p.close(force=True)
        shutil.rmtree(nb.directory())
    p.dispose = dispose
    return p
