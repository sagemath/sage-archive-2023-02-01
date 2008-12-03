"""
The Sage Notebook object
"""

#############################################################################
#       Copyright (C) 2006, 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import os
import random
import re
import shutil
import socket
import time
import bz2
import cPickle

# Sage libraries
from   sage.structure.sage_object import SageObject, load
from   sage.misc.misc       import (alarm, cancel_alarm,
                                    tmp_dir, pad_zeros, cputime)
from   sage.misc.package   import is_package_installed
from   sage.version        import version
# Sage Notebook
import css          # style
import js           # javascript
import worksheet    # individual worksheets (which make up a notebook)
import config       # internal configuration stuff (currently, just keycodes)
import keyboards    # keyboard layouts
import server_conf  # server configuration
import user_conf    # user configuration
import user         # users


SYSTEMS = ['sage', 'gap', 'gp', 'jsmath', 'html', 'latex', 'maxima', 'python', 'r', 'sage', 'sh', 'singular', 'axiom (optional)', 'kash (optional)', 'macaulay2 (optional)', 'magma (optional)', 'maple (optional)', 'mathematica (optional)', 'matlab (optional)', 'mupad (optional)', 'octave (optional)']

# We also record the system names without (optional) since they are
# used in some of the html menus, etc.
SYSTEM_NAMES = [v.split()[0] for v in SYSTEMS]

JSMATH = True

JQUERY = True

if is_package_installed("jsmath-image-fonts"):
    JSMATH_IMAGE_FONTS = True
else:
    JSMATH_IMAGE_FONTS = False

vbar = '<span class="vbar"></span>'

DOC_TIMEOUT = 120

class Notebook(SageObject):
    def __init__(self,
                 dir,
                 system=None,
                 pretty_print=False,
                 show_debug = False,
                 address='localhost',
                 port=8000,
                 secure=True,
                 server_pool = []):
        if isinstance(dir, basestring) and len(dir) > 0 and dir[-1] == "/":
            dir = dir[:-1]
        self.__dir = dir
        self.__absdir = os.path.abspath(dir)

        self.__server_pool = server_pool
        self.set_system(system)
        self.set_pretty_print(pretty_print)
        self.__worksheets = {}
        self.__filename      = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir    = '%s/objects'%dir
        self.__makedirs()
        self.__history = []
        self.__history_count = 0
        self.__server_log = [] #server log list
        self.__show_debug = show_debug
        self.__admins = []
        self.__conf = server_conf.ServerConfiguration()

        # Install this copy of the notebook in twist.py as *the*
        # global notebook object used for computations.  This is
        # mainly to avoid circular references, etc.  This also means
        # only one notebook can actually be used at any point.
        import sage.server.notebook.twist
        sage.server.notebook.twist.notebook = self

        # This must happen after twist.notebook is set.
        self.save()

    def _migrate_worksheets(self):
        v = []
        for key, W in self.__worksheets.iteritems():
            if not '/' in W.filename():
                v.append((key, W))
        if len(v) > 0:
            print "Migrating from old to new worksheet format"
            D = self.directory()
##             if os.path.exists('%s/worksheets'%D):
##                 import shutil
##                 target = '%s/../old_worksheets.tar.bz2'%D
##                 print "First archiving old worksheets and objects directory to '%s'"%target
##                 os.system('tar jcf "%s" "%s/worksheets" "%s/objects"'%(target, D, D))
##                 ws_tree = "%s/worksheets"%D
##                 print "Now removing ", ws_tree
##                 shutil.rmtree(ws_tree)
##                 obj_tree = "%s/objects"%D
##                 if os.path.exists(obj_tree):
##                     shutil.rmtree(obj_tree)
            for key, W in v:
                print W.name()
                txt = W.edit_text()
                N = self.create_new_worksheet(W.name(), 'pub')
                N.edit_save(txt, ignore_ids=True)
                del self.__worksheets[key]
            print "Your old worksheets are all available by clicking the published link"
            print "in the upper right corner."
            print "If you want to save disk space, you could immediately remove"
            print "the objects and worksheets directories in your Sage notebook, as"
            print "they are no longer used.  Do this now or never."

    def delete(self):
        """
        Delete all files related to this notebook.

        This is used for doctesting mainly.  This command
        is obviously *VERY* dangerous to use on a notebook
        you actually care about.  You could easily lose
        all data.

        EXAMPLES:
            sage: tmp = tmp_dir()
            sage: nb = sage.server.notebook.notebook.Notebook(tmp)
            sage: sorted(os.listdir(tmp))
            ['backups', 'nb.sobj', 'objects', 'worksheets']
            sage: nb.delete()

        Now the directory is gone.
            sage: os.listdir(tmp)
            Traceback (most recent call last):
            ...
            OSError: [Errno 2] No such file or directory: '...
        """
        try:
            dir = self.__absdir
        except AttributeErrro:
            dir = self.__dir
        import shutil
        # We ignore_errors because in rare parallel doctesting
        # situations sometimes the directory gets cleaned up too
        # quickly, etc.
        shutil.rmtree(dir, ignore_errors=True)

    ##########################################################
    # Users
    ##########################################################
    def create_default_users(self, passwd):
        """
        Create the default users for a notebook.

        INPUT:
            passwd -- a string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: list(sorted(nb.users().iteritems()))
            [('_sage_', _sage_), ('admin', admin), ('guest', guest), ('pub', pub)]
            sage: list(sorted(nb.passwords().iteritems()))
            [('_sage_', 'aaQSqAReePlq6'), ('admin', 'aajfMKNH1hTm2'), ('guest', 'aaQSqAReePlq6'), ('pub', 'aaQSqAReePlq6')]
            sage: nb.create_default_users('newpassword')
            Creating default users.
            WARNING: User 'pub' already exists -- and is now being replaced.
            WARNING: User '_sage_' already exists -- and is now being replaced.
            WARNING: User 'guest' already exists -- and is now being replaced.
            WARNING: User 'admin' already exists -- and is now being replaced.
            sage: list(sorted(nb.passwords().iteritems()))
            [('_sage_', 'aaQSqAReePlq6'), ('admin', 'aajH86zjeUSDY'), ('guest', 'aaQSqAReePlq6'), ('pub', 'aaQSqAReePlq6')]
        """
        print "Creating default users."
        self.add_user('pub', '', '', account_type='user', force=True)
        self.add_user('_sage_', '', '', account_type='user', force=True)
        self.add_user('guest', '', '', account_type='guest', force=True)
        self.add_user('admin', passwd, '', account_type='admin', force=True)

    def user_exists(self, username):
        """
        Return whether or not a user exists given a username.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.user_exists('admin')
            True
            sage: nb.user_exists('pub')
            True
            sage: nb.user_exists('mark')
            False
            sage: nb.user_exists('guest')
            True
        """
        return username in self.users()

    def users(self):
        """
        Return dictionary of users in a notebook.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: list(sorted(nb.users().iteritems()))
            [('_sage_', _sage_), ('admin', admin), ('guest', guest), ('pub', pub)]
        """
        try:
            return self.__users
        except AttributeError:
            self.__users = {}
            return self.__users

    def user(self, username):
        """
        Return an instance of the User class given the username of a user in
        a notebook.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.user('admin')
            admin
            sage: nb.user('admin')._User__email
            ''
            sage: nb.user('admin')._User__password
            'aajfMKNH1hTm2'
        """
        if not isinstance(username, str) or '/' in username:
            raise KeyError, "no user '%s'"%username
        try:
            return self.users()[username]
        except KeyError:
            if username in ['pub', '_sage_']:
                self.add_user(username, '', '', account_type='user', force=True)
                return self.users()[username]
            elif username == 'admin':
                self.add_user(username, '', '', account_type='admin', force=True)
                return self.users()[username]
            elif username == 'guest':
                self.add_user('guest', '', '', account_type='guest', force=True)
                return self.users()[username]
            raise KeyError, "no user '%s'"%username

    def create_user_with_same_password(self, user, other_user):
        r"""
        INPUT:
            user -- a string
            other_user -- a string
        OUTPUT:
            Changes password of \var{user} to that of \var{other_user}.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('bob', 'an**d', 'bob@gmail.com', force=True)
            sage: nb.user('bob').password()
            'aa4Q6Jbx/MiUs'
            sage: nb.add_user('mary', 'ccd', 'mary@gmail.com', force=True)
            sage: nb.user('mary').password()
            'aaxr0gcWJMXKU'
            sage: nb.create_user_with_same_password('bob', 'mary')
            sage: nb.user('bob').password() == nb.user('mary').password()
            True
        """
        U = self.user(user)
        O = self.user(other_user)
        passwd = O.password()
        U.set_hashed_password(passwd)

    def user_is_admin(self, user):
        """
        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Administrator', 'password', '', 'admin', True)
            sage: nb.add_user('RegularUser', 'password', '', 'user', True)
            sage: nb.user_is_admin('Administrator')
            True
            sage: nb.user_is_admin('RegularUser')
            False
        """
        return self.user(user).is_admin()

    def user_is_guest(self, username):
        """
        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.user_is_guest('guest')
            True
            sage: nb.user_is_guest('admin')
            False
        """
        try:
            return self.user(username).is_guest()
        except KeyError:
            return False

    def user_list(self):
        """
        Return list of user objects.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: sorted(nb.user_list(), key=lambda k: k.username())
            [_sage_, admin, guest, pub]
        """
        return list(self.users().itervalues())

    def usernames(self):
        """
        Return list of usernames.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: sorted(nb.usernames())
            ['_sage_', 'admin', 'guest', 'pub']
        """
        U = self.users()
        return U.keys()

    def valid_login_names(self):
        """
        Return list of users that can be signed in.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.valid_login_names()
            ['admin']
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: nb.add_user('Sarah', 'password', '', force=True)
            sage: nb.add_user('David', 'password', '', force=True)
            sage: sorted(nb.valid_login_names())
            ['David', 'Mark', 'Sarah', 'admin']
        """
        return [x for x in self.usernames() if not x in ['guest', '_sage_', 'pub']]

    def default_user(self):
        """
        Return a default login name that the user will see when confronted with the
        Sage notebook login page.

        OUTPUT:
            string

        Currently this returns 'admin' if that is the *only* user.  Otherwise it
        returns the string ''.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.default_user()
            'admin'
            sage: nb.add_user('AnotherUser', 'password', '', force=True)
            sage: nb.default_user()
            ''
        """
        if self.valid_login_names() == ['admin']:
            return 'admin'
        else:
            return ''

    def set_accounts(self, value):
        r"""
        Changes \var{__accounts} to \var{value}

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.get_accounts()
            False
            sage: nb.set_accounts(True)
            sage: nb.get_accounts()
            True
            sage: nb.set_accounts(False)
            sage: nb.get_accounts()
            False
        """
        self.__accounts = value

    def get_accounts(self):
        r"""
        Return \var{__accounts}

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.get_accounts()
            False
            sage: nb.set_accounts(True)
            sage: nb.get_accounts()
            True
        """
        try:
            return self.__accounts
        except AttributeError:
            self.__accounts = False
            return False

    def add_user(self, username, password, email, account_type="user", force=False):
        """
        INPUT:
            username -- the username
            password -- the password
            email -- the email address
            account_type -- one of 'user', 'admin', or 'guest'
            force -- bool

        If the method get_accounts return False then user can only be added if force=True

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: nb.user('Mark')
            Mark
            sage: nb.add_user('Sarah', 'password', '')
            Traceback (most recent call last):
            ValueError: creating new accounts disabled.
            sage: nb.set_accounts(True)
            sage: nb.add_user('Sarah', 'password', '')
            sage: nb.user('Sarah')
            Sarah
        """
        if not self.get_accounts() and not force:
            raise ValueError, "creating new accounts disabled."

        us = self.users()
        if us.has_key(username):
            print "WARNING: User '%s' already exists -- and is now being replaced."%username
        U = user.User(username, password, email, account_type)
        us[username] = U

    def change_password(self, username, password):
        """
        INPUT:
            username -- the username
            password -- the password to change the user's password to

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: nb.user('Mark').password()
            'aajfMKNH1hTm2'
            sage: nb.change_password('Mark', 'different_password')
            sage: nb.user('Mark').password()
            'aaTlXok5npQME'
        """
        self.user(username).set_password(password)

    def del_user(self, username):
        """
        Deletes the given user

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: nb.user('Mark')
            Mark
            sage: nb.del_user('Mark')
            sage: nb.user('Mark')
            Traceback (most recent call last):
            KeyError: "no user 'Mark'"
        """
        us = self.users()
        if us.has_key(username):
            del us[username]

    def passwords(self):
        """
        Return the username:password dictionary.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: list(sorted(nb.passwords().iteritems()))
            [('Mark', 'aajfMKNH1hTm2'), ('_sage_', 'aaQSqAReePlq6'), ('admin', 'aajfMKNH1hTm2'), ('guest', 'aaQSqAReePlq6'), ('pub', 'aaQSqAReePlq6')]
        """
        return dict([(user.username(), user.password()) for user in self.user_list()])

    def user_conf(self, username):
        """
        Return a user's configuration.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: config = nb.user_conf('admin')
            sage: config['max_history_length']
            500
            sage: config['default_system']
            'sage'
            sage: config['autosave_interval']
            180
            sage: config['default_pretty_print']
            False
        """
        return self.users()[username].conf()

    ##########################################################
    # Publishing worksheets
    ##########################################################
    def _initialize_worksheet(self, src, W):
        r"""
        \var{src} and \var{W} are worksheets and \var{W} is brand new.
        """
        # Copy over images and other files
        data = src.data_directory()
        if os.path.exists(data):
            shutil.copytree(data, W.directory() + '/data')
        cells = src.cells_directory()
        if os.path.exists(cells):
            shutil.copytree(cells, W.directory() + '/cells')
        W.edit_save(src.edit_text())

    def publish_worksheet(self, worksheet, username):
        r"""
        Publish the given worksheet.

        This creates a new worksheet in the \file{pub} directory with the
        same contents as \var{worksheet}.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Mark','password','',force=True)
            sage: W = nb.new_worksheet_with_title_from_text('First steps', owner='Mark')
            sage: nb.worksheet_names()
            ['Mark/0']
            sage: nb.publish_worksheet(nb.get_worksheet_with_filename('Mark/0'), 'Mark')
            <BLANKLINE>
            [Cell 0; in=, out=]
            sage: sorted(nb.worksheet_names())
            ['Mark/0', 'pub/0']
        """
        for X in self.__worksheets.itervalues():
            if X.is_published() and X.worksheet_that_was_published() == worksheet:
                # Update X based on worksheet instead of creating something new
                # 1. delete cells and data directories
                # 2. copy them over
                # 3. update worksheet text
                if os.path.exists(X.data_directory()):
                    shutil.rmtree(X.data_directory(), ignore_errors=True)
                if os.path.exists(X.cells_directory()):
                    shutil.rmtree(X.cells_directory(), ignore_errors=True)
                if os.path.exists(X.snapshot_directory()):
                    shutil.rmtree(X.snapshot_directory(), ignore_errors=True)
                self._initialize_worksheet(worksheet, X)
                X.set_worksheet_that_was_published(worksheet)
                X.move_to_archive(username)
                worksheet.set_published_version(X.filename())
                X.record_edit(username)
                return X

        # Have to create a new worksheet
        W = self.create_new_worksheet(worksheet.name(), 'pub')
        self._initialize_worksheet(worksheet, W)
        W.set_worksheet_that_was_published(worksheet)
        W.move_to_archive(username)
        worksheet.set_published_version(W.filename())
        return W

    ##########################################################
    # Moving, copying, creating, renaming, and listing worksheets
    ##########################################################

    def scratch_worksheet(self):
        try:
            return self.__scratch_worksheet
        except AttributeError:
            W = self.create_new_worksheet('scratch', '_sage_', add_to_list=False)
            self.__scratch_worksheet = W
            return W

    def create_new_worksheet(self, worksheet_name, username, docbrowser=False, add_to_list=True):
        if username!='pub' and self.user_is_guest(username):
            raise ValueError, "guests cannot create new worksheets"

        filename = worksheet.worksheet_filename(worksheet_name, username)
        if self.__worksheets.has_key(filename):
            return self.__worksheets[filename]
        i = 0
        dir = self.worksheet_directory() + '/' + username
        if os.path.exists(dir):
            D = os.listdir(dir)
            D.sort()
            dirname = str(i)
            while dirname in D:
                i += 1
                dirname = str(i)
        else:
            dirname = '0'

        W = worksheet.Worksheet(worksheet_name, dirname, self.worksheet_directory(),
                                system = self.system(username),
                                owner=username,
                                docbrowser = docbrowser,
                                auto_publish = False)

        if add_to_list:
            self.__worksheets[W.filename()] = W
        return W

    def copy_worksheet(self, ws, owner):
        W = self.create_new_worksheet('default', owner)
        self._initialize_worksheet(ws, W)
        name = "Copy of %s"%ws.name()
        W.set_name(name)
        return W

    def delete_worksheet(self, filename):
        """
        Delete the given worksheet and remove its name from the
        worksheet list.
        """
        if not (filename in self.__worksheets.keys()):
            print self.__worksheets.keys()
            raise KeyError, "Attempt to delete missing worksheet '%s'"%filename
        W = self.__worksheets[filename]
        W.quit()
        shutil.rmtree(W.directory(), ignore_errors=True)
        self.deleted_worksheets()[filename] = W
        del self.__worksheets[filename]

    def deleted_worksheets(self):
        try:
            return self.__deleted_worksheets
        except AttributeError:
            self.__deleted_worksheets = {}
            return self.__deleted_worksheets

    def empty_trash(self, username):
        """
        Empty the trash for the given user.

        INPUT:
            username -- a string

        This empties the trash for the given user and cleans up all
        files associated with the worksheets that are in the trash.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.new_worksheet_with_title_from_text('Sage', owner='sage')
            sage: W.move_to_trash('sage')
            sage: nb.worksheet_names()
            ['sage/0']
            sage: nb.empty_trash('sage')
            sage: nb.worksheet_names()
            []
        """
        X = self.get_worksheets_with_viewer(username)
        X = [W for W in X if W.is_trashed(username)]
        for W in X:
            W.delete_user(username)
            if W.owner() is None:
                self.delete_worksheet(W.filename())

    def worksheet_names(self):
        """
        Return a list of all the names of worksheets in this notebook.

        OUTPUT:
            list of strings.

        EXAMPLES:
        We make a new notebook with two users and two worksheets, then list their names:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.new_worksheet_with_title_from_text('Sage', owner='sage')
            sage: nb.add_user('wstein','sage','wstein@sagemath.org',force=True)
            sage: W2 = nb.new_worksheet_with_title_from_text('Elliptic Curves', owner='wstein')
            sage: nb.worksheet_names()
            ['sage/0', 'wstein/0']
        """
        W = self.__worksheets.keys()
        W.sort()
        return W

    def migrate_old(self):
        """
        Migrate all old worksheets, i.e., ones with no owner,
        to \file{/pub}.
        """
        raise NotImplementedError
        for w in self.__worksheets.itervalues():
            if not '/' in w.filename():
                print "Moving worksheet ", w.name()
                w.set_owner('old')
                self.rename_worksheet_filename(w, w.filename())

    ##########################################################
    # Information about the pool of worksheet compute servers
    ##########################################################

    def server_pool(self):
        try:
            return self.__server_pool
        except AttributeError:
            self.__server_pool = []
            return []

    def set_server_pool(self, servers):
        self.__server_pool = servers

    def get_ulimit(self):
        try:
            return self.__ulimit
        except AttributeError:
            self.__ulimit = ''
            return ''

    def set_ulimit(self, ulimit):
        self.__ulimit = ulimit

    def get_server(self):
        P = self.server_pool()
        if P is None or len(P) == 0:
            return None
        try:
            self.__server_number = (self.__server_number + 1)%len(P)
            i = self.__server_number
        except AttributeError:
            self.__server_number = 0
            i = 0
        return P[i]

    ##########################################################
    # The default math software system for new worksheets for
    # a given user or the whole notebook (if username is None).
    ##########################################################

    # TODO -- only implemented for the notebook right now
    def system(self, username=None):
        return self.user(username).conf()['default_system']

    def set_system(self, system):
        self.__system = system

    ##########################################################
    # The default typeset setting for new worksheets for
    # a given user or the whole notebook (if username is None).
    ##########################################################

    # TODO -- only implemented for the notebook right now
    def pretty_print(self, username=None):
        return self.user(username).conf()['default_pretty_print']

    def set_pretty_print(self, pretty_print):
        self.__pretty_print = pretty_print

    ##########################################################
    # The default color scheme for the notebook.
    ##########################################################
    def color(self):
        try:
            return self.__color
        except AttributeError:
            self.__color = 'default'
            return self.__color

    def set_color(self,color):
        self.__color = color

    ##########################################################
    # The directory the notebook runs in.
    ##########################################################
    def set_directory(self, dir):
        if dir == self.__dir:
            return
        if isinstance(dir, basestring) and len(dir) > 0 and dir[-1] == "/":
            dir = dir[:-1]
        self.__dir = dir
        self.__filename = '%s/nb.sobj'%dir
        self.__worksheet_dir = '%s/worksheets'%dir
        self.__object_dir = '%s/objects'%dir

    ##########################################################
    # The notebook history.
    ##########################################################
    def user_history(self, username):
        U = self.user(username)
        return U.history_list()

    def create_new_worksheet_from_history(self, name, username, maxlen=None):
        W = self.create_new_worksheet(name, username)
        W.edit_save('Log Worksheet\n' + self.user_history_text(username, maxlen=None))
        return W

    def user_history_text(self, username, maxlen=None):
        H = self.user_history(username)
        if maxlen:
            H = H[-maxlen:]
        return '\n\n'.join([L.strip() for L in H])

    def user_history_html(self, username):
        t = self.user_history_text(username)
        t = t.replace('<','&lt;')
        s = """
        <html>
        <head>
           <link rel=stylesheet href="/css/main.css">
           <title>Sage: History for %s</title>
        </head>
        <body>
        %s
        <pre>
        %s
        </pre>
        <hr class="usercontrol">
        <a title="Click here to turn the above into a Sage worksheet" href="/live_history">Create a new Sage worksheet version of the last 100 commands in the above log.</a>
        <a name="bottom"></a>
        <script type="text/javascript"> window.location="#bottom"</script>
        </body>
        </html>
        """%(username, self.html_worksheet_list_top(username, actions=False), t)
        return s


    def add_to_user_history(self, entry, username):
        H = self.user_history(username)
        H.append(entry)
        maxlen = self.user_conf(username)['max_history_length']
        while len(H) > maxlen:
            del H[0]

    def add_to_history(self, input_text):
        H = self.history()
        H.append(input_text)
        while len(H) > self.max_history_length():
            del H[0]

    def history_count_inc(self):
        self.__history_count += 1

    def history_count(self):
        return self.__history_count

    def history(self):
        try:
            s = self.__history
        except AttributeError:
            self.__history = []
            s = self.__history
        return s

    def history_text(self):
        return '\n\n'.join([H.strip() for H in self.history()])

    def max_history_length(self):
        try:
            return self.conf()['max_history_length']
        except KeyError:
            return MAX_HISTORY_LENGTH

    def history_html(self):
        t = self.history_text()
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>Command History</title>\n'
        s += '</head>\n'
        s += '<body>\n'
        s += '<pre>' + t + '</pre>\n'
        s += '<a name="bottom"></a>\n'
        s += '<script type="text/javascript"> window.location="#bottom"</script>\n'
        s += '</body>\n'
        return s


    def history_with_start(self, start):
        n = len(start)
        return [x for x in self.history() if x[:n] == start]

    ##########################################################
    # Importing and exporting worksheets to files
    ##########################################################
    def export_worksheet(self, worksheet_filename, output_filename, verbose=True):
        """
        Export a worksheet with given directory filenmae to output_filename.

        INPUT:
            worksheet_filename -- string
            output_filename -- string
            verbose -- bool (default: True) if True print some the tar
                       command used to extract the sws file.

        OUTPUT:
            creates a file on the filesystem
        """
        W = self.get_worksheet_with_filename(worksheet_filename)
        W.save()
        path = W.filename_without_owner()
        cmd = 'cd "%s/%s/" && tar -jcf "%s" "%s"'%(
            self.__worksheet_dir, W.owner(),
            os.path.abspath(output_filename), path)
        if verbose:
            print cmd
        e = os.system(cmd)
        if e:
            print "Failed to execute command to export worksheet:\n'%s'"%cmd

    def new_worksheet_with_title_from_text(self, text, owner):
        name, _ = worksheet.extract_name(text)
        W = self.create_new_worksheet(name, owner)
        return W

    def change_worksheet_key(self, old_key, new_key):
        ws = self.__worksheets
        W = ws[old_key]
        ws[new_key] = W
        del ws[old_key]

    def import_worksheet(self, filename, owner):
        r"""
        Upload the worksheet with name \var{filename} and make it have the
        given owner.

        INPUT:
            filename -- a string
            owner -- a string

        OUTPUT:
            worksheet -- a newly created worksheet

        EXAMPLES:
        We create a notebook and import a plain text worksheet into it.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: name = tmp_filename() + '.txt'
            sage: open(name,'w').write('foo\n{{{\n2+3\n}}}')
            sage: W = nb.import_worksheet(name, 'admin')

        W is our newly-created worksheet, with the 2+3 cell in it:
            sage: W.name()
            'foo'
            sage: W.cell_list()
            [Cell 0; in=2+3, out=]
        """
        if not os.path.exists(filename):
            raise ValueError, "no file %s"%filename

        # Figure out the file extension
        ext = os.path.splitext(filename)[1]
        if ext.lower() == '.txt':
            # A plain text file with {{{'s that defines a worksheet (not graphics).
            return self._import_worksheet_txt(filename, owner)
        elif ext.lower() == '.sws':
            # An sws file (really a tar.bz2) which defines a worksheet with graphics,
            # revisions, etc.
            return self._import_worksheet_sws(filename, owner)
        else:
            # We only support txt or sws files.
            raise ValueError, "unknown extension '%s'"%ext

    def _import_worksheet_txt(self, filename, owner):
        """
        Import a plain text file as a new worksheet.

        INPUT:
            filename -- string; a filename that ends in .txt
            owner -- string; who will own this worksheet when imported

        OUTPUT:
            a new worksheet

        EXAMPLES:
            We write a plain text worksheet to a file and import it using this function.

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: name = tmp_filename() + '.txt'
            sage: open(name,'w').write('foo\n{{{\na = 10\n}}}')
            sage: W = nb._import_worksheet_txt(name, 'admin'); W
            [Cell 0; in=a = 10, out=]
        """
        # Open the worksheet txt file and load it in.
        worksheet_txt = open(filename).read()
        # Create a new worksheet with the write title and owner.
        worksheet = self.new_worksheet_with_title_from_text(worksheet_txt, owner)
        # Set the new worksheet to have the contents specified by that file.
        worksheet.edit_save(worksheet_txt)
        return worksheet

    def _import_worksheet_sws(self, filename, owner, verbose=True):
        """
        Import an sws format worksheet into this notebook as a new worksheet.

        INPUT:
            filename -- string; a filename that ends in .sws; internally
                        it must be a tar'd bz2'd file.
            owner -- string
            verbose -- bool (default: True) if True print some the tar
                       command used to extract the sws file.

        OUTPUT:
            a new worksheet

        EXAMPLES:
        We create a notebook, then make a worksheet from a plain text file first.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: name = tmp_filename() + '.txt'
            sage: open(name,'w').write('foo\n{{{\n2+3\n}}}')
            sage: W = nb.import_worksheet(name, 'admin')
            sage: W.filename()
            'admin/0'


        We then export the worksheet to an sws file.
            sage: nb.export_worksheet(W.filename(),  'tmp.sws', verbose=False)

        Now we import the sws.
            sage: nb._import_worksheet_sws('tmp.sws', 'admin', verbose=False)
            [Cell 0; in=2+3, out=]

        Yep, it's there now (as admin/2):
            sage: nb.worksheet_names()
            ['admin/0', 'admin/2']
        """
        # Decompress the worksheet to a temporary directory.
        tmp = tmp_dir()
        cmd = 'cd "%s"; tar -jxf "%s"'%(tmp, os.path.abspath(filename))
        if verbose:
            print cmd
        e = os.system(cmd)
        if e:
            raise ValueError, "Error decompressing saved worksheet."

        # Find the worksheet text representation and load it into memory.
        try:
            D = os.listdir(tmp)[0]
        except IndexError:
            raise ValueError, "invalid worksheet"
        text_filename = '%s/%s/worksheet.txt'%(tmp,D)
        worksheet_txt = open(text_filename).read()
        worksheet = self.new_worksheet_with_title_from_text(worksheet_txt, owner)
        worksheet.set_owner(owner)
        name = worksheet.filename_without_owner()

        # Change the filename of the worksheet, if necessary
        names = [w.filename_without_owner() for w in self.get_worksheets_with_owner(owner)]
        if name in names:
            name = 0
            while str(name) in names:
                name += 1
            name = str(name)
            worksheet.set_filename_without_owner(name)

        # Change the display name of the worksheet if necessary
        name = worksheet.name()
        display_names = [w.name() for w in self.get_worksheets_with_owner(owner)]
        if name in display_names:
            j = name.rfind('(')
            if j != -1:
                name = name[:j].rstrip()
            i = 2
            while name + " (%s)"%i in display_names:
                i += 1
            name = name + " (%s)"%i
            worksheet.set_name(name)


        # Put the worksheet files in the target directory.
        S = self.__worksheet_dir
        target = '%s/%s'%(os.path.abspath(S), worksheet.filename())
        if not os.path.exists(target):
            os.makedirs(target)
        cmd = 'rm -rf "%s"/*; mv "%s/%s/"* "%s/"'%(target, tmp, D, target)
        #print cmd
        if os.system(cmd):
            raise ValueError, "Error moving over files when loading worksheet."

        worksheet.edit_save(worksheet_txt)

        shutil.rmtree(tmp, ignore_errors=True)

        return worksheet


    ##########################################################
    # Importing and exporting worksheets to a plain text format
    ##########################################################

    def plain_text_worksheet_html(self, name, prompts=True):
        W = self.get_worksheet_with_filename(name)
        t = W.plain_text(prompts = prompts)
        t = t.replace('<','&lt;')
        s = '<head>\n'
        s += '<title>Sage Worksheet: %s</title>\n'%W.name()
        s += '</head>\n'
        s += '<body>\n'
        s += '<h1><a href=".">Sage Worksheet: %s</a></h1>\n'%W.name()
        s += '<pre>' + t + '</pre>'
        s += '</body>\n'
        return s

    ##########################################################
    # Directories for worksheets, etc.
    ##########################################################
    def directory(self):
        if not os.path.exists(self.__dir):
            # prevent "rm -rf" accidents.
            os.makedirs(self.__dir)
        return self.__dir

    def DIR(self):
        """
        Return the absolute path to the directory that contains
        the Sage Notebook directory.
        """
        P = os.path.abspath('%s/..'%self.__dir)
        if not os.path.exists(P):
            # prevent "rm -rf" accidents.
            os.makedirs(P)
        return P

    def worksheet_directory(self):
        return self.__worksheet_dir

    def __makedirs(self):
        if not os.path.exists(self.__dir):
            os.makedirs(self.__dir)
        if not os.path.exists(self.__worksheet_dir):
            os.makedirs(self.__worksheet_dir)
        if not os.path.exists(self.__object_dir):
            os.makedirs(self.__object_dir)

    ##########################################################
    # Server configuration
    ##########################################################
    def conf(self):
        try:
            return self.__conf
        except AttributeError:
            C = server_conf.ServerConfiguration()
            self.__conf = C
            return C


    def set_debug(self,show_debug):
        self.__show_debug = show_debug

    def number_of_backups(self):
        return self.conf()['number_of_backups']

    def backup_directory(self):
        try:
            D = self.__backup_dir
        except AttributeError:
            D = self.__dir + "/backups/"
            self.__backup_dir = D
        if not os.path.exists(D):
            os.makedirs(D)
        return D


    ##########################################################
    # The object store for the notebook.
    ##########################################################
    # Todo: like with worksheets, objects should belong to
    # users, some should be published, rateable, etc.
    def object_directory(self):
        O = self.__object_dir
        if not os.path.exists(O):
            os.makedirs(O)
        return O

    def objects(self):
        L = [x[:-5] for x in os.listdir(self.object_directory())]
        L.sort()
        return L

    def object_list_html(self):
        m = max([len(x) for x in self.objects()] + [30])
        s = []
        a = '<a href="/%s.sobj" class="object_name">\n'
        for name in self.objects():
            s.append(a%name + name + '</a>\n')  # '&nbsp;'*(m-len(name)) +
        return '<br>\n'.join(s)

    ##########################################################
    # Computing control
    ##########################################################
    def set_not_computing(self):
        # unpickled, no worksheets will think they are
        # being computed, since they clearly aren't (since
        # the server just started).
        for W in self.__worksheets.values():
            W.set_not_computing()

    def quit(self):
        for W in self.__worksheets.itervalues():
            W.quit()

    def quit_idle_worksheet_processes(self):
        timeout = self.conf()['idle_timeout']
        if timeout == 0:
            # Quit only the doc browser worksheets
            for W in self.__worksheets.itervalues():
                if W.docbrowser() and W.compute_process_has_been_started():
                    W.quit_if_idle(DOC_TIMEOUT)
            return

        for W in self.__worksheets.itervalues():
            if W.compute_process_has_been_started():
                W.quit_if_idle(timeout)


    ##########################################################
    # Worksheet HTML generation
    ##########################################################
    def list_window_javascript(self, worksheet_filenames):
        s = """
           <script type="text/javascript" src="/javascript/main.js"></script>
           <script type="text/javascript">
           var worksheet_filenames = %s;
           </script>
        """%(worksheet_filenames)

        return s

    def worksheet_html(self, filename, do_print=False):
        W = self.get_worksheet_with_filename(filename)
        s = '<head>\n'
        s += '<title>Sage Worksheet: %s</title>\n'%W.name()
        s += '<script type="text/javascript" src="/javascript/main.js"></script>\n'
        if do_print:
            s += '<script type="text/javascript" src="/javascript/jsmath/jsMath.js"></script>\n'
        s += '<link rel=stylesheet href="/css/main.css">\n'
        s += '</head>\n'
        if do_print:
            s += '<body>\n'
            s += '<div class="worksheet_print_title">%s</div>'%W.name()
        else:
            s += '<body onLoad="initialize_the_notebook();">\n'
        s += W.html(include_title=False, do_print=do_print)
        if do_print:
            s += '<script type="text/javascript">jsMath.Process();</script>\n'
        s += '\n</body>\n'
        return s



    def worksheet_list_for_public(self, username, sort='last_edited', reverse=False, search=None):
        W = [x for x in self.__worksheets.itervalues() if x.is_published() and not x.is_trashed(user)]

        if search:
            W = [x for x in W if x.satisfies_search(search)]

        sort_worksheet_list(W, sort, reverse)  # changed W in place
        return W

    def worksheet_list_for_user(self, user, typ="active", sort='last_edited', reverse=False, search=None):
        X = self.get_worksheets_with_viewer(user)
        if typ == "trash":
            W = [x for x in X if x.is_trashed(user)]
        elif typ == "active":
            W = [x for x in X if x.is_active(user)]
        else: # typ must be archived or "all"
            W = [x for x in X if not x.is_trashed(user)]
        if search:
            W = [x for x in W if x.satisfies_search(search)]
        sort_worksheet_list(W, sort, reverse)  # changed W in place
        return W

    def html_topbar(self, user, pub=False):
        s = ''
        entries = []

        if self.user_is_guest(user):
            entries.append(('/', 'Log in', 'Please log in to the Sage notebook'))
        else:
            entries.append(('/home/%s'%user, 'Home', 'Back to your personal worksheet list'))
            entries.append(('/pub', 'Published', 'Browse the published worksheets'))
            entries.append(('help()', 'Help', 'Documentation'))

        ## TODO -- settings
        #if self.user(user).is_admin():
        #    entries.insert(1, ('/notebook_settings', 'Server', 'Change general Sage notebook server configuration'))
        if not pub:
            entries.insert(2, ('history_window()', 'Log', 'View a log of recent computations'))
        if not self.user_is_guest(user):
            entries.append(('/settings', 'Settings', 'Change account settings including password'))
            entries.append(('/logout', 'Sign out', 'Log out of the Sage notebook'))

        s += self.html_banner_and_control(user, entries)
        s += '<hr class="usercontrol">'
        return s

    def html_worksheet_list_top(self, user, actions=True, typ='active', pub=False, search=None):
        s = self.html_topbar(user, pub)
        if not pub:
            s += self.html_new_or_upload()
        s += self.html_search(search, typ)
        s += '<br>'
        s += '<hr class="usercontrol">'
        if actions:
            s += self.html_worksheet_actions(user, typ=typ)
        return s

    def html_banner_and_control(self, user, entries):
        return """
        <table width="100%%"><tr><td>
        %s
        </td><td align=right>
        %s
        </td></tr>
        </table>
        """%(self.html_banner(),
             self.html_user_control(user, entries))


    def html_user_control(self, user, entries):
        s = ''
        s += '<span class="username">%s</span>'%user
        for href, name, title in entries:
            if '(' in href:
                action = 'onClick="%s"'%href
            else:
                action = 'href="%s"'%href
            x = '<a title="%s" class="usercontrol" %s>%s</a>\n'%(title, action, name)
            s += vbar + x
        return s

    def html_banner(self):
        ver=version
        s = """
        <div class="banner">
        <table width="100%%"><tr><td>
        <a class="banner" href="http://www.sagemath.org"><img align="top" src="/images/sagelogo.png" alt="Sage"> Notebook</a></td><td><span class="ping" id="ping">Searching for Sage server...</span></td>
        </tr><tr><td style="font-size:xx-small; text-indent:13px; color:black">Version %s</td><td></td></tr></table>
        </div>
        """%ver
        return s


    def html_worksheet_actions(self, user, typ):
        s = ''

        if not self.user_is_guest(user):
            if typ == 'archive':
                s += '<button onClick="make_active_button();" title="Unarchive selected worksheets so it appears in the default worksheet list">Unarchive</button>'
            else:
                s += '<button onClick="archive_button();" title="Archive selected worksheets so they do not appear in the default worksheet list">Archive</button>'

            if typ != 'trash':
                s += '&nbsp;&nbsp;<button onClick="delete_button();" title="Move the selected worksheets to the trash">Delete</button>'
            else:
                s += '&nbsp;&nbsp;<button onClick="make_active_button();" title="Move the selected worksheets out of the trash">Undelete</button>'

            s += '&nbsp;&nbsp;<button onClick="stop_worksheets_button();" title="Stop selected worksheets">Stop</button>'

            s += '<span>'
            s += '&nbsp;'*10
            #s += '<a class="control" href="/pub" title="Browse everyone\'s published worksheets">Published Worksheets</a>'
            s += '&nbsp;'*10
            s += "Current Folder: "
            s += '&nbsp;<a class="%susercontrol" href=".">Active</a>'%('bold' if typ=='active' else '')
            s += '&nbsp;<a class="%susercontrol" href=".?typ=archive">Archived</a>'%('bold' if typ=='archive' else '')
            s += '&nbsp;<a class="%susercontrol" href=".?typ=trash">Trash</a>&nbsp;&nbsp;'%('bold' if typ=='trash' else '')
            if typ == 'trash':
                s += '&nbsp;<a class="boldusercontrol" onClick="empty_trash();return false;" href="">(Empty Trash)</a>'
            s += '</span>'
        return s


    ##########################################################
    # Revision history for a worksheet
    ##########################################################
    def html_worksheet_revision_list(self, username, worksheet):
        head, body = self.html_worksheet_page_template(worksheet, username, "Revision history", select="revisions")
        data = worksheet.snapshot_data()  # pairs ('how long ago', key)
        rows = []
        i = 0
        for i in range(len(data)):
            desc, key = data[i]
            rows.append('<tr><td></td><td><a href="revisions?rev=%s">Revision %s</a></td><td><span class="revs">%s</span></td></tr>'%
                        (key, i, desc))

        rows = list(reversed(rows))
        rows = '\n'.join(rows)
        body += """
        <hr class="usercontrol">
<table width="100%%">
<tr><td width="1%%"></td><td width="20%%"><b>Revision</b></td> <td width="20%%"><b>Last Edited</b></td><td width="30%%"></td>
%s
</table>
"""%rows

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)


    def html_specific_revision(self, username, ws, rev):
        t = time.time() - float(rev[:-4])
        when = worksheet.convert_seconds_to_meaningful_time_span(t)
        head, body = self.html_worksheet_page_template(ws, username,
                                       "Revision from %s ago&nbsp;&nbsp;&nbsp;&nbsp;<a href='revisions'>Revision List</a>"%when, select="revisions")

        filename = ws.get_snapshot_text_filename(rev)
        txt = bz2.decompress(open(filename).read())
        W = self.scratch_worksheet()
        W.delete_cells_directory()
        W.edit_save(txt)
        html = W.html_worksheet_body(do_print=True, publish=True)

        data = ws.snapshot_data()  # pairs ('how long ago', key)
        prev_rev = None
        next_rev = None
        for i in range(len(data)):
            if data[i][1] == rev:
                if i > 0:
                    prev_rev = data[i-1][1]
                if i < len(data)-1:
                    next_rev = data[i+1][1]
                break

        if prev_rev:
            prev = '<a class="listcontrol" href="revisions?rev=%s">Older</a>&nbsp;&nbsp;'%prev_rev
        else:
            prev = 'Oldest'

        if next_rev:
            next = '<a class="listcontrol" href="revisions?rev=%s">Newer</a>&nbsp;&nbsp;'%next_rev
        else:
            next = 'Newest'

        actions = """
        %s
        %s
        <a class="listcontrol" href="revisions?rev=%s&action=revert">Revert to this one</a>  <span class="lastedit">(note that images are note recorded)</span>&nbsp;&nbsp;
        <a class="listcontrol" href="revisions?rev=%s&action=publish">Publish this one</a>&nbsp;&nbsp;
        """%(prev, next, rev, rev)

        s = """
        %s
        <hr class="usercontrol">
<table width="100%%">
%s
        <hr class="usercontrol">
%s
</table>
"""%(actions, html, actions)
        body += s

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_worksheet_page_template(self, worksheet, username, title, select=None, backwards=False):
        head = self._html_head(worksheet_filename=worksheet.filename(), username=username)
        head += '<script  type="text/javascript">worksheet_filename="%s"; worksheet_name="%s"; server_ping_while_alive(); </script>'%(worksheet.filename(), worksheet.name())
        body = self._html_body(worksheet.filename(), top_only=True, username=username)
        body += self.html_worksheet_topbar(worksheet, select=select, username=username, backwards=backwards)
        body += '<hr class="usercontrol">'
        body += '<span class="sharebar">%s</span>'%title
        body += '<br>'*3
        return head, body


    def html_share(self, worksheet, username):
        head, body = self.html_worksheet_page_template(worksheet, username, "Share this document", select="share")

        if not (self.user(username).is_admin() or username == worksheet.owner()):
            body += "Only the owner of a worksheet is allowed to share it."
            body += 'You can do whatever you want if you <a href="copy">make your own copy</a>.'
        else:
            body += 'This Sage Worksheet is currently shared with the people listed in the box below.<br>'
            body += 'You may add or remove collaborators (separate user names by commas).<br><br>'

            collabs = ', '.join(worksheet.collaborators())
            body += '<form width=70% method="post" action="invite_collab">\n'
            body += '<textarea name="collaborators" rows=5 cols=70 class="edit" id="collaborators">%s</textarea><br><br>'%collabs
            body += '<input type="submit" title="Give access to your worksheet to the above collaborators" value="Invite Collaborators">'
            body += '</form>'

            body += '<br>'*2
            body += '<hr class="usercontrol">'
            body += '<span class="username">Sage Users:</span>'
            U = self.users()
            K = [x for x, u in U.iteritems() if not u.is_guest() and not u.username() in [username, 'pub', '_sage_']]
            def mycmp(x,y):
                return cmp(x.lower(), y.lower())
            K.sort(mycmp)
            body += '<span class="users">%s</span>'%(', '.join(K))


        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)



    def html_download_or_delete_datafile(self, ws, username, filename):
        head, body = self.html_worksheet_page_template(ws, username, "Data file: %s"%filename)
        path = "/home/%s/data/%s"%(ws.filename(), filename)
        body += 'You may download <a href="%s">%s</a>'%(path, filename)

        X = self.get_worksheets_with_viewer(username)
        v = [x for x in X if x.is_active(username)]
        sort_worksheet_list(v, 'name', False)
        ws_form = ['<option selected>select worksheet</option>'] + \
                  ["""<option value='link_datafile("%s","%s")'>%s</option>"""%(
                           x.filename(), filename, x.name()) for x in v]
        ws_form = '\n'.join(ws_form)
        ws_form = "<select onchange='go_option(this);' class='worksheet'>%s</select>"%ws_form
        body += ' or create a linked copy to the worksheet %s,'%ws_form
        body += ' or <a href="/home/%s/datafile?name=%s&action=delete">delete %s.</a>'%(ws.filename(),filename, filename)

        body += "<br><br>Access %s in this worksheet by typing <tt>DATA+'%s'</tt>.  Here DATA is a special variable that gives the exact path to all data files uploaded to this worksheet.<br><br>"%(filename, filename)

        body += '<hr class="usercontrol">'
        ext = os.path.splitext(filename)[1].lower()
        if ext in ['.png', '.jpg', '.gif']:
            body += '<div align=center><img src="%s"></div>'%path
        elif ext in ['.txt', '.tex', '.sage', '.spyx', '.py', '.f', '.f90', '.c']:
            body += '<form method="post" action="savedatafile" enctype="multipart/form-data">'
            body += '<input type="submit" value="Save Changes" name="button_save"> <input type="submit" value="Cancel" name="button_cancel"><br>'
            body += '<textarea class="edit" name="textfield" rows=17 cols=70 id="textfield">%s</textarea>'%open('%s/%s'%(ws.data_directory(), filename)).read()
            body += '<input type="hidden" name="filename" value="%s" id="filename">'%filename
            body += '</form>'

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)




    ##########################################################
    # Accessing all worksheets with certain properties.
    ##########################################################
    def get_all_worksheets(self):
        return [x for x in self.__worksheets.itervalues() if not x.owner() in ['_sage_', 'pub']]

    def get_worksheets_with_collaborator(self, user):
        if self.user_is_admin(user): return self.get_all_worksheets()
        return [w for w in self.__worksheets.itervalues() if w.user_is_collaborator(user)]

    def get_worksheet_names_with_collaborator(self, user):
        if self.user_is_admin(user): return [W.name() for W in self.get_all_worksheets()]
        return [W.name() for W in self.get_worksheets_with_collaborator(user)]

    def get_worksheets_with_viewer(self, user):
        if self.user_is_admin(user): return self.get_all_worksheets()
        return [w for w in self.__worksheets.itervalues() if w.user_is_viewer(user)]

    def get_worksheets_with_owner(self, owner):
        return [w for w in self.__worksheets.itervalues() if w.owner() == owner]

    def get_worksheets_with_owner_that_are_viewable_by_user(self, owner, user):
        return [w for w in self.get_worksheets_with_owner(owner) if w.user_is_viewer(user)]

    def get_worksheet_names_with_viewer(self, user):
        if self.user_is_admin(user): return [W.name() for W in self.get_all_worksheets()]
        return [W.name() for W in self.get_worksheets_with_viewer(user) if not W.docbrowser()]

    def get_worksheet_with_name(self, name):
        for W in self.__worksheets.itervalues():
            if W.name() == name:
                return W
        raise KeyError, "No worksheet with name '%s'"%name

    def get_worksheet_with_filename(self, filename):
        """
        Get the worksheet with given filename.  If there is no such
        worksheet, raise a \exception{KeyError}.

        INPUT:
            string
        OUTPUT:
            a worksheet or KeyError
        """
        if self.__worksheets.has_key(filename):
            return self.__worksheets[filename]
        raise KeyError, "No worksheet with filename '%s'"%filename

    ###########################################################
    # Saving the whole notebook
    ###########################################################

    def save(self, filename=None, verbose=False):

        if filename is None:
            F = os.path.abspath(self.__filename)
            backup_dir = self.backup_directory()
            backup = backup_dir + '/nb-backup-'
            for i in range(self.number_of_backups()-1,0,-1):
                a = pad_zeros(i-1); b = pad_zeros(i)
                try:
                    shutil.move(backup + '%s.sobj'%a, backup + '%s.sobj'%b)
                except IOError, msg:
                    pass
            a = '%s.sobj'%pad_zeros(0)
            try:
                shutil.copy(F, backup + a)
            except Exception, msg:
                pass
            F = os.path.abspath(self.__filename)
        else:
            F = os.path.abspath(filename)

        D, _ = os.path.split(F)
        if not os.path.exists(D):
            os.makedirs(D)

        t = cputime()
        out = cPickle.dumps(self, 2)
        if verbose: print "Dumped notebook to pickle (%s seconds)"%cputime(t)

        t = cputime()
        # Assuming an exception wasn't raised during pickling we write to the file.
        # This is vastly superior to writing to a file immediately, which can easily
        # result in a poor empty file.
        open(F,'w').write(out)
        if verbose: print "Wrote notebook pickle to file '%s' (%s seconds)"%(F,cputime(t))

    def delete_doc_browser_worksheets(self):
        names = self.worksheet_names()
        for n in self.__worksheets.keys():
            if n.startswith('doc_browser'):
                self.delete_worksheet(n)

    ###########################################################
    # HTML -- generate most html related to the whole notebook page
    ###########################################################
    def html_slide_controls(self):
        return """
          <div class="hidden" id="slide_controls">
            <div class="slideshow_control">
             <a class="slide_arrow" onClick="slide_next()">&gt;</a>
              <a class="slide_arrow" onClick="slide_last()">&gt;&gt;</a> %s
              <a class="cell_mode" onClick="cell_mode()">Exit</a>
            </div>
            <div class="slideshow_progress" id="slideshow_progress" onClick="slide_next()">
              <div class="slideshow_progress_bar" id="slideshow_progress_bar">&nbsp;</div>
              <div class="slideshow_progress_text" id="slideshow_progress_text">&nbsp;</div>
            </div>
            <div class="slideshow_control">
              <a class="slide_arrow" onClick="slide_first()">&lt;&lt;</a>
              <a class="slide_arrow" onClick="slide_prev()">&lt;</a>
            </div>
          </div>
          """%vbar

    def html_debug_window(self):
        return """
        <div class='debug_window'>
        <div class='debug_output'><pre id='debug_output'></pre></div>
        <textarea rows=5 id='debug_input' class='debug_input'
         onKeyPress='return debug_keypress(event);'
         onFocus='debug_focus();' onBlur='debug_blur();'></textarea>
        </div>"""


    def _html_head(self, worksheet_filename, username):
        if worksheet_filename is not None:
            worksheet = self.get_worksheet_with_filename(worksheet_filename)
            head = '\n<title>%s (Sage)</title>'%(worksheet.name())
        else:
            head = '\n<title>Sage Notebook | Welcome</title>'

        # Load the Sage javascript libray.
        head += '\n<script type="text/javascript" src="/javascript/main.js"></script>\n'
        head += '\n<link rel=stylesheet href="/css/main.css" type="text/css">\n'

        if JSMATH:
            # turn off the ugly scary font warning.
            head += '\n <STYLE> #jsMath_Warning {display: none} </STYLE>\n'
            head += '<script type="text/javascript">jsMath = {Controls: {cookie: {scale: 115}}}</script>\n'
            if not JSMATH_IMAGE_FONTS:
                head +=' <script type="text/javascript" src="/javascript/jsmath/plugins/noImageFonts.js"></script>\n'

            # Move the jsMath button 20 pixels from the right edge
            # (apparently in some browsers, it covers up the scroll
            # bar)
            head += """<script type="text/javascript">
jsMath = {styles: {
        '#jsMath_button':   'position:fixed; bottom:1px; right:20px; background-color:white; '
                                + 'border: solid 1px #959595; margin:0px; padding: 0px 3px 1px 3px; '
                                + 'z-index:102; color:black; text-decoration:none; font-size:x-small; '
                                + 'width:auto; cursor:hand;',
      }};
    </script>
"""
            head += '<script type="text/javascript" src="/javascript/jsmath/jsMath.js"></script>\n'

        if JQUERY:
            # Load the jquery and ui-jquery javascript library.
            # This is used for interact functionality in the notebook, and will be used
            # to enable drag and drop, image zoom, etc.
            head += '''
<script type="text/javascript" src="/javascript/jquery/jquery.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/jquery.dimensions.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.mouse.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.slider.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.draggable.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.draggable.ext.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.resizable.js"></script>
<script type="text/javascript" src="/javascript/jqueryui/ui.dialog.js"></script>
<link rel="stylesheet" href="/javascript/jqueryui/themes/flora/flora.all.css">

<script type="text/javascript" src="/javascript/farbtastic/farbtastic.js"></script>
<link rel="stylesheet" href="/javascript/farbtastic/farbtastic.css" type="text/css" />
         '''

        # This was for syntax hilighting
#        head +=' <script type="text/javascript" src="/javascript/highlight/prettify.js"></script>\n'
#        head += '<link rel=stylesheet href="/css/highlight/prettify.css" type="text/css">\n'

        head +=' <script type="text/javascript" src="/javascript/sage3d.js"></script>\n'

        # Jmol -- embedded 3d graphics.
        head +=' <script type="text/javascript" src="/java/jmol/appletweb/Jmol.js"></script>\n'

        head +=' <script>jmolInitialize("/java/jmol");</script>\n' # this must stay in the <head>
        return head

    def html_worksheet_topbar(self, worksheet, select=None, username='guest', backwards=False):
        body = ''
        body += """
<table width="100%%" id="topbar">
<tr>
  <td align=left> %s </td>   <td align=right> %s </td>
</tr>
<tr>
  <td align=left> %s </td>   <td align=right> %s </td>
</tr>
</table>
"""%(worksheet.html_title(username), worksheet.html_save_discard_buttons(),
     worksheet.html_menu(), worksheet.html_share_publish_buttons(select=select, backwards=backwards))

        body += self.html_slide_controls()
        return body


    def _html_body(self, worksheet_filename, show_debug=False, username='', top_only=False):
        worksheet = self.get_worksheet_with_filename(worksheet_filename)
        worksheet_html = worksheet.html()

        body = ''

        if worksheet.is_published() or self.user_is_guest(username):
            original_worksheet = worksheet.worksheet_that_was_published()
            if original_worksheet.user_is_collaborator(username) or original_worksheet.is_owner(username):
                s = "Edit this."
                url = 'edit_published_page'
            elif self.user_is_guest(username):
                s = 'Log in to edit a copy.'
                url = '/'
            else:
                s = 'Edit a copy.'
                url = 'edit_published_page'
            r = worksheet.rating()
            if r == -1:
                rating = ''
            else:
                rating = '<a class="usercontrol" href="rating_info">This page is rated %.1f.</a>'%r
            if not self.user_is_guest(username) \
                   and not worksheet.is_publisher(username):
                if worksheet.is_rater(username):
                    action = "Rerate"
                else:
                    action = "Rate"
                rating += '&nbsp;&nbsp; <span class="usercontrol">%s it: </span>'%action
                rating += '  '.join(['<a class="usercontrol" onClick="rate_worksheet(%s)">&nbsp;%s&nbsp;</a>'%(i,i) for
                                   i in range(5)])
                rating += '&nbsp;&nbsp; <input name="rating_comment" id="rating_comment"></input>'

            download_name = os.path.split(worksheet.name())[-1]
            edit_line = '<a class="usercontrol" href="%s">%s</a>'%(url, s) + \
                        '  <a class="usercontrol" href="download/%s.sws">Download.</a>'%download_name + \
                        '  <span class="ratingmsg">%s</span>'%rating

            body += edit_line
            #This document was published using <a href="/">Sage</a>.'
            body += '<span class="pubmsg">'
            body += '<a href="/pub/">Other published documents...</a></span>'
            body += '<hr class="usercontrol">'
            body += '<h1 align=center>%s</h1>'%original_worksheet.name()
            body += '<h2 align=center>%s</h2>'%worksheet.html_time_since_last_edited()
            body += worksheet_html
            body += '<hr class="usercontrol">'
            body += '&nbsp;'*10



        else:

            entries = [('toggle_top()', 'Toggle', 'Toggle the top bar'),
                       ('/', 'Home', 'Back to your personal worksheet list'),
                       ('/pub', 'Published', 'Browse the published worksheets'),
                       ('history_window()', 'Log', 'View a log of recent computations'),
                       ('/settings', 'Settings', 'Account Settings'),
                       ('bugreport()', 'Report a Problem', 'Report a problem or submit a bug to improve Sage'),
                       ('help()', 'Help', 'Documentation')]

            if not self.user_is_guest(username):
                entries.append(('/logout', 'Sign out', 'Log out of the Sage notebook'))

            body += self.html_banner_and_control(username, entries)
            if top_only:
                return body

            if worksheet_filename:
                body += self.html_worksheet_topbar(worksheet, select="use", username=username)

            if self.__show_debug or show_debug:
                body += self.html_debug_window()


            body += '<div class="worksheet" id="worksheet">%s</div>'%worksheet_html

        endpanespan = '</td></tr></table></span>\n'


        if worksheet is None:
             return body + endpanespan

        if worksheet.user_is_only_viewer(username):
            body += '<script type="text/javascript">worksheet_locked=true;</script>'
        else:
            body += '<script type="text/javascript">worksheet_locked=false;</script>'

        if worksheet.computing():
            # Set the update checking back in motion.
            body += '<script type="text/javascript"> active_cell_list = %r; \n'%worksheet.queue_id_list()
            body += 'for(var i = 0; i < active_cell_list.length; i++)'
            body += '    cell_set_running(active_cell_list[i]); \n'
            body += 'start_update_check();\n'
            body +=' </script>\n'

        return body

    def html_plain_text_window(self, worksheet, username):
        """
        Return a window that displays a plain text version of the worksheet

        INPUT:
            worksheet -- a worksheet
            username -- name of the user
        """
        head, body = self.html_worksheet_page_template(worksheet, username, 'View plain text', select="text")

        t = worksheet.plain_text(prompts=True, banner=False)
        t = t.replace('<','&lt;')
        body += """
        <pre class="plaintext" id="cell_intext" name="textfield">%s
        </pre>
        """%t.strip()

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_edit_window(self, worksheet, username):
        r"""
        Return a window for editing \var{worksheet}.

        INPUT:
            worksheet -- a worksheet
        """
        head, body = self.html_worksheet_page_template(worksheet, username, 'Edit plain text &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type="submit" value="Save Changes" name="button_save" id="button_save"> <input type="submit" value="Cancel" name="button_cancel">', select="edit")


        body += """<script type="text/javascript">
function save_worksheet() {
}
function save_worksheet_and_close() {
}
        </script>
        """
        t = worksheet.edit_text()
        t = t.replace('<','&lt;')
        body = '<form method="post" action="save" enctype="multipart/form-data">' + body
        body += """
        <textarea class="plaintextedit" id="cell_intext" name="textfield" rows="%s">%s</textarea>
        </form>
        """%(t.count("\n")+1,t)

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_beforepublish_window(self, worksheet, username):
        """
        Return the html code for a page dedicated to worksheet publishing prior to the
        publication of the given worksheet.

        INPUT:
            worksheet - instance of Worksheet
            username - string
        """
        msg = """You can publish your worksheet to the Internet, where anyone will be able to access and view it online.
        Your worksheet will be assigned a unique address (URL) that you can send to your friends and colleagues.<br/><br/>
        Do you want to publish this worksheet?<br/><br/>
        <form method="get" action=".">
        <input type="hidden" name="yes" value="" />
        <input type="submit" value="Yes" style="margin-left:10px" />
        <input type="button" value="No" style="margin-left:5px" onClick="parent.location=\'../'"><br/><br/>
        <input type="checkbox" name="auto" style="margin-left:13px" /> Automatically re-publish when changes are made
        </form>
        """
        head, body = self.html_worksheet_page_template(worksheet, username, msg, select="publish", backwards=True)

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_afterpublish_window(self, worksheet, username, addr, dtime):
        """
        Return the html code for a page dedicated to worksheet publishing after the
        publication of the given worksheet.

        INPUT:
            worksheet - instance of Worksheet
            username - string
            addr - string
            dtime - instance of time.struct_time
        """
        from time import strftime
        time = strftime("%B %d, %Y %I:%M %p", dtime)
        msg = """Worksheet is publicly viewable at <a href="%s" style="color:#FFF" target="_blank">%s</a><br />
        Published on %s<br/><br />
        <input type="button" value="Re-publish worksheet" onClick="parent.location=\'?re'"><input type="button" value="Stop publishing" style="margin-left:5px" onClick="parent.location=\'?stop'"><br /><br />
<input type="checkbox" name="auto"%s onchange="parent.location=\'?auto'"/> Automatically re-publish when changes are made
        """ % (addr, addr, time, ' checked="true" ' if worksheet.is_auto_publish() else '')
        head, body = self.html_worksheet_page_template(worksheet, username, msg, select="publish", backwards=True)

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_notebook_help_window(self, username):
        top = self._html_head(None, username) + self.html_topbar(username)

        from tutorial import notebook_help
        s = """
        <html>
        <title>Sage Documentation</title>

        <body>
        """ + top + \
        """
        <style>

        div.help_window {
            font-family: sans-serif;
            background-color:white;
            border: 3px solid #3d86d0;
            top: 10ex;
            bottom:10%;
            left:25%;
            right:15%;
            padding:2ex;
            width:80%;
        }


        table.help_window {
            background-color:white;
            width:95%;
        }

        td.help_window_cmd {
            background-color: #f5e0aa;
            width:30%;
            padding:1ex;
            font-weight:bold;
        }

        td.help_window_how {
            padding:1ex;
            width:70%;
        }
        </style>

        <center>
        <br>
        <a class="control" title="To quickly try out Sage start here" href="/doc/live/tut/index.html">Tutorial</a>
        &nbsp;&nbsp;
        <a class="control" title="View a 2000 page reference manual about Sage" href="/doc/live/ref/index.html">Reference Manual</a>
        &nbsp;&nbsp;
        <a class="control" title="Learn to write Sage programs" href="/doc/live/prog/index.html">Programming Guide</a>
        &nbsp;&nbsp;
        <a class="control" title="How do I construct ... in Sage?" href="/doc/live/const/const.html">Constructions</a>
        <br><br>
        <center>
        <a class="control" title="Static version..." href="/doc/static/">Fast Static Versions of the Documentation</a>
        </center>
        <hr class="usercontrol">
        <br>

        <div class="help_window">
        <h2>How to use the Sage Notebook</h2>

        <br><br>

        A <i>worksheet</i> is an ordered list of Sage calculations with output. <br>
        A <i>session</i> is a worksheet and a set of variables in some state.<br>
        The <i>Sage notebook</i> is a collection of worksheets, saved objects, and user information.<br>
        <br>
        <br>
        To get started with Sage, <a href="/doc/live/tut/tut.html">work through the tutorial</a> (if
        you have trouble with it, view the <a href="/doc/static/tut/tut.html">static version</a>).
        <br><br>

        <table class="help_window">
        """

        for x, y in notebook_help:
            s += '<tr><td class="help_window_cmd">%s</td><td class="help_window_how">%s</td></tr>\n'%(x,y)
        s += '</table></div>'

        s +="""
        <br>        <br>
        The Sage Notebook was primarily written by William Stein with substantial contributions from Tom Boothby, Timothy Clemans, Alex Clemesha, Bobby Moretti, Yi Qiang, and Dorian Raymer.
        </center>
        </body>
        </html>
        """
        return s

    def upload_window(self):
        return """
          <html>
            <head>
              <title>Upload File</title>
              <style>%s</style>
            </head>
            <body>
              <div class="upload_worksheet_menu" id="upload_worksheet_menu">
              %s
              <h1><font size=+1>Upload worksheet (an sws or txt file) to the Sage Notebook</font></h1>
              <hr>
              <form method="POST" action="upload_worksheet"
                    name="upload" enctype="multipart/form-data">
              <table><tr>
              <td>
              Browse your computer to select a worksheet file to upload:<br>
              <input class="upload_worksheet_menu" size="50" type="file" name="fileField" id="upload_worksheet_filename"></input><br><br>
              Or enter the url of a worksheet file on the web:<br>

              <input class="upload_worksheet_menu" size="50" type="text" name="urlField" id="upload_worksheet_url"></input>
              <br><br>
              What do you want to call it? (if different than the original name)<br>
              <input class="upload_worksheet_menu" size="50" type="text" name="nameField" id="upload_worksheet_name"></input></br>
              </td>
              </tr>
              <tr>
              <td><br><input type="button" class="upload_worksheet_menu" value="Upload Worksheet" onClick="form.submit();"></td>
              </tr>
              </form><br>
              </div>
            </body>
          </html>
         """%(css.css(self.color()),self.html_banner())

    def html_upload_data_window(self, ws, username):
        head, body = self.html_worksheet_page_template(ws, username, "Upload or Create Data File")

        body += """
              <div class="upload_worksheet_menu" id="upload_worksheet_menu">
              <h1><font size=+1>Upload or create data file attached to the worksheet '%s'</font></h1>
              <hr>
              <form method="POST" action="do_upload_data"
                    name="upload" enctype="multipart/form-data">
              <table><tr>
              <td>
              Browse your computer to select a file to upload:<br>
              <input class="upload_worksheet_menu" size="50" type="file" name="fileField" value="" id="upload_filename"></input><br><br>
              Or enter the url of a file on the web:<br>

              <input class="upload_worksheet_menu" size="50" type="text" name="urlField" value="" id="upload_url"></input></br>
              <br><br>
              Or enter the name of a new file, which will be created:<br>
              <input class="upload_worksheet_menu" size="50" type="text" name="newField" value="" id="upload_filename"></input><br><br>

              What do you want to call it? (if different than the original name)<br>
              <input class="upload_worksheet_menu" size="50" type="text" name="nameField" value="" id="upload_name"></input></br>
              </td>
              </tr>
              <tr>
              <td><br><input type="button" class="upload_worksheet_menu" value="Upload File" onClick="form.submit();"></td>
              </tr>
              </form><br>
              </div>
            </body>
          </html>
         """%(ws.name())

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html(self, worksheet_filename=None, username='guest', show_debug=False, admin=False):
        if worksheet_filename is None or worksheet_filename == '':
            worksheet_filename = None
            W = None
        else:
            try:
                W = self.get_worksheet_with_filename(worksheet_filename)
            except KeyError:
                W = None

        head = self._html_head(worksheet_filename=worksheet_filename, username=username)
        body = self._html_body(worksheet_filename=worksheet_filename, username=username, show_debug=show_debug)

        head += '<script type="text/javascript">user_name="%s"; </script>'%username

        if worksheet_filename is not None:
            head += '<script  type="text/javascript">worksheet_filename="%s"; worksheet_name="%s"; server_ping_while_alive(); </script>'%(worksheet_filename, W.name())

            # Uncomment this to force rename when the worksheet is opened (annoying!)
            #if W and W.name() == "Untitled":
            #    head += '<script  type="text/javascript">setTimeout("rename_worksheet()",1)</script>'

        return """
        <!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
        <html>
        <head>%s</head>
        <body onLoad="initialize_the_notebook();">%s</body>
        </html>
        """%(head, body)

    def _html_authorize(self):
        return """
        <h1>Sage Notebook Server</h1>
        <div id="mainbody" class="login">Sign in to the Sage Notebook<br>
        <form>
        <table>
        <tr><td>
          <span class="username">Username:</span></td>
          <td><input name="username" class="username"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <tr><td>
           <span class="password">Password:</span></td>
           <td><input name="password" class="username" type="password"
                      onKeyPress="if(is_submit(event)) login(username.value, password.value)"></td>
        </tr>
        <td>&nbsp</td>
        <td>
           <input type='button' onClick="login(username.value,password.value);" value="Sign in">
           </td></table>
                   </form></div>

        """


    ####################################################################
    # Configuration html.
    # In each case the settings html is a form that when submitted
    # pulls up another web page and sets the corresponding options.
    ####################################################################
    def html_system_select_form_element(self, ws):
        system = ws.system()
        options = ''
        i = SYSTEM_NAMES.index(system)
        for j, S in enumerate(SYSTEMS):
            if i == j:
                selected = "selected"
            else:
                selected = ''
            T = SYSTEM_NAMES[j]
            options += '<option title="Evaluate all input cells using %s" %s value="%s">%s</option>\n'%(T, selected, T,S)
        s = """<select  onchange="go_system_select(this, %s);" class="worksheet">
            %s
            </select>"""%(i, options)
        return s

    def html_pretty_print_check_form_element(self, ws):
        pretty_print = ws.pretty_print()
        if pretty_print:
            check='checked="checked"'
        else:
            check=''
        s = """<input type="checkbox" title="Enable/disable pretty_printing"
        onchange="pretty_print_check(this.checked);"
        class="worksheet" value="pretty_print" %s>&nbsp;Typeset"""%(check)
        return s


    def html_worksheet_settings(self, ws, username):
        head, body = self.html_worksheet_page_template(ws, username, 'Worksheet Settings &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<button name="button_save">Save Settings</button>  <input type="submit" value="Cancel" name="button_cancel"/>')

        body = '<form width=70%% method="post" action="input_settings"  enctype="multipart/form-data">' + body
        body += '</form>'

        return """
        <html>
        <head>%s</head>
        <body>%s</body>
        </html>
        """%(head, body)

    def html_settings(self):
        s = """
        <h1>Settings</h1>
        """
        return s

    def html_user_settings(self, username):
        s = self.html_settings()
        return s

    def html_notebook_settings(self):
        s = self.html_settings()
        return s

    def html_doc(self, username):
        top = self._html_head(None, username) + self.html_topbar(username)
        body = """
        <br>
        <div class="docidx">
        <h1>Sage Documentation</h1>
        <br>
        <hr class="usercontrol">
        <br><br>
        <font size=+2>
        <a href="/doc/live/">Live Documentation</a><br><br>
        <a href="/doc/static/">Static Documentation</a><br><br>
        <a href="/help/">Sage Notebook Howto</a><br><br>
        <br><br>
        <br>
        <hr class="usercontrol">
        </font>
        </div>
        """
        #(<a href="/doc/static/">static</a>)

        s = """
        <html>
        %s
        <body>
        %s
        </body>
        </html>
        """%(top, body)

        return s

    def html_src(self, filename, username):
        top = self._html_head(None, username) + self.html_topbar(username)

        if not os.path.exists(filename):
            file = "No such file '%s'"%filename
        else:
            file = open(filename).read()
        file = file.replace('<','&lt;')
        s = """
<html>
<head>
"""
        s += '<title>%s | Sage Source Code</title>' % filename
#        s += '<link rel=stylesheet href="/highlight/prettify.css" type="text/css" />\n'
        s += """
</head>
<body>
"""
        s += '<h1 align=center>Sage Source Browser</h1>\n'
        s += '<h2 align=center><tt>%s  <a href="..">(browse directory)</a></tt></h2>\n'%filename
        s += '<br><hr><br>\n'
        s += '<font size=+1><pre id="code">%s</pre></font>\n'%file
        s += '<br><hr><br>\n'
#        s += '<script src="/highlight/prettify.js" type="text/javascript"></script>\n'
        s += """<script type="text/javascript">
function get_element(id) {
  if(document.getElementById)
    return document.getElementById(id);
  if(document.all)
    return document.all[id];
  if(document.layers)
    return document.layers[id];
}

var x = get_element("code");
x.innerHTML = prettyPrintOne(x.innerHTML);
</script>
"""

        s = """
        <html>
        %s
        <body>
        %s
        </body>
        </html>
        """%(top, s)

        return s


####################################################################

def load_notebook(dir, address=None, port=None, secure=None):
    """
    Load the notebook from the given directory, or create one in that directory
    if one isn't already there.

    INPUT:
        dir -- a string that defines a directory name
        address -- the address that the notebook server will listen on
        port -- the port the server listens on
        secure -- whether or not the notebook is secure
    """
    sobj = '%s/nb.sobj'%dir
    nb = None
    if os.path.exists(dir):
        try:
            nb = load(sobj, compress=False)
        except:
            backup = '%s/backups/'%dir
            if os.path.exists(backup):
                print "****************************************************************"
                print "  * * * WARNING   * * * WARNING   * * * WARNING   * * * "
                print "WARNING -- failed to load notebook object. Trying backup files."
                print "****************************************************************"
                for F in os.listdir(backup):
                    file = backup + '/' + F
                    try:
                        nb = load(file, compress=False)
                    except Exception, msg:
                        print "Failed to load backup '%s'"%file
                    else:
                        print "Successfully loaded backup '%s'"%file
                        break
                if nb is None:
                    print "Unable to restore notebook from *any* auto-saved backups."
                    print "This is a serious problem."
    if nb is None:
        nb = Notebook(dir)
    dir = make_path_relative(dir)
    nb.delete_doc_browser_worksheets()
    nb.set_directory(dir)
    nb.set_not_computing()
    nb.address = address
    nb.port = port
    nb.secure = secure

    return nb


def make_path_relative(dir):
    r"""
    If easy, replace the absolute path \var{dir} by a relative one.
    """
    base, file = os.path.split(dir)
    if os.path.exists(file):
        return file
    return dir

##########################################################
# Misc
##########################################################

def clean_name(name):
    return ''.join([x if (x.isalnum() or x == '_') else '_' for x in name])

def sort_worksheet_list(v, sort, reverse):
    """
    INPUT:
        sort -- 'last_edited', 'owner', 'rating', or 'name'
        reverse -- if True, reverse the order of the sort.
    """
    f = None
    if sort == 'last_edited':
        def c(a, b):
            return -cmp(a.last_edited(), b.last_edited())
        f = c
    elif sort == 'name':
        def c(a,b):
            return cmp((a.name().lower(), -a.last_edited()), (b.name().lower(), -b.last_edited()))
        f = c
    elif sort == 'owner':
        def c(a,b):
            return cmp((a.owner().lower(), -a.last_edited()), (b.owner().lower(), -b.last_edited()))
        f = c
    elif sort == "rating":
        def c(a,b):
            return -cmp((a.rating(), -a.last_edited()), (b.rating(), -b.last_edited()))
        f = c
    else:
        raise ValueError, "invalid sort key '%s'"%sort
    v.sort(cmp = f, reverse=reverse)
