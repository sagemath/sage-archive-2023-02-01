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
from template import template

from cgi import escape

# latex macros
from sage.misc.latex_macros import sage_jsmath_macros

SYSTEMS = ['sage', 'gap', 'gp', 'jsmath', 'html', 'latex', 'maxima', 'python', 'r', 'sage', 'sh', 'singular', 'axiom (optional)', 'kash (optional)', 'macaulay2 (optional)', 'magma (optional)', 'maple (optional)', 'mathematica (optional)', 'matlab (optional)', 'mupad (optional)', 'octave (optional)']

# We also record the system names without (optional) since they are
# used in some of the html menus, etc.
SYSTEM_NAMES = [v.split()[0] for v in SYSTEMS]

JSMATH = True

if is_package_installed("jsmath-image-fonts"):
    JSMATH_IMAGE_FONTS = True
else:
    JSMATH_IMAGE_FONTS = False

if is_package_installed("tinyMCE"):
    JEDITABLE_TINYMCE = True
else:
    JEDITABLE_TINYMCE = False

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

        This is used for doctesting mainly. This command is obviously
        *VERY* dangerous to use on a notebook you actually care about.
        You could easily lose all data.

        EXAMPLES::

            sage: tmp = tmp_dir()
            sage: nb = sage.server.notebook.notebook.Notebook(tmp)
            sage: sorted(os.listdir(tmp))
            ['backups', 'nb.sobj', 'objects', 'worksheets']
            sage: nb.delete()

        Now the directory is gone.

        ::

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

    def systems(self):
        return SYSTEMS

    def system_names(self):
        return SYSTEM_NAMES
    ##########################################################
    # Users
    ##########################################################
    def create_default_users(self, passwd):
        """
        Create the default users for a notebook.

        INPUT:

        -  ``passwd`` - a string

        EXAMPLES::

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
        Return whether a user with the given ``username`` exists.

        INPUT:

        - ``username`` - a string

        OUTPUT:

        - a bool

        EXAMPLES::

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
        Return a dictionary of users in a notebook.

        OUTPUT:

        - a string:User instance dictionary

        EXAMPLES::

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
        Return an instance of the User class given the ``username`` of a user
        in a notebook.

        INPUT:

        - ``username`` - a string

        OUTPUT:

        - an instance of User

        EXAMPLES::

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
        Change the password of ``user`` to that of ``other_user``.

        INPUT:

        -  ``user`` - a string

        -  ``other_user`` - a string

        EXAMPLES::

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
        Return true if ``user`` is an admin.

        INPUT:

        - ``user`` - an instance of User

        OUTPUT:

        - a bool

        EXAMPLES::

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
        Return true if ``username`` is a guest.

        INPUT:

        - ``username`` - a string

        OUTPUT:

        - a bool

        EXAMPLES::

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
        Return a list of user objects.

        OUTPUT:

        - a list of User instances

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: sorted(nb.user_list(), key=lambda k: k.username())
            [_sage_, admin, guest, pub]
        """
        return list(self.users().itervalues())

    def usernames(self):
        """
        Return a list of usernames.

        OUTPUT:

        - a list of strings

        EXAMPLES::

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
        Return a list of users that can log in.

        OUTPUT:

        - a list of strings

        EXAMPLES::

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
        r"""
        Return a default login name that the user will see when
        confronted with the Sage notebook login page.  Currently, this
        returns 'admin' if that is the *only* user.  Otherwise it
        returns an empty string ('').

        OUTPUT:

        - a string - the default username.

        EXAMPLES::

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
        Set the attribute ``__accounts`` to ``value``.

        INPUT:

        - ``value`` - a bool

        EXAMPLES::

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
        Return ``__accounts``.

        OUTPUT:

        - a bool

        EXAMPLES::

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
        Add a user with the given credentials.

        INPUT:

        -  ``username`` - the username

        -  ``password`` - the password

        -  ``email`` - the email address

        -  ``account_type`` - one of 'user', 'admin', or 'guest'

        -  ``force`` - a bool (default: False)

        If the method :meth:`get_accounts` returns False then user can
        only be added if ``force`` is True.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('Mark', 'password', '', force=True)
            sage: nb.user('Mark')
            Mark
            sage: nb.add_user('Sarah', 'password', ")
            Traceback (most recent call last):
            ValueError: creating new accounts disabled.
            sage: nb.set_accounts(True)
            sage: nb.add_user('Sarah', 'password', ")
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
        Change a user's password.

        INPUT:

        - ``username`` - a string, the username

        - ``password`` - a string, the user's new password

        EXAMPLES::

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
        Delete the given user.

        INPUT:

        - ``username`` - a string

        EXAMPLES::

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
        Return a username:password dictionary.

        OUTPUT:

        - a string:string dictionary

        EXAMPLES::

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
        Return a user's configuration object.

        OUTPUT:

        - an instance of Configuration.

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.create_default_users('password')
            Creating default users.
            sage: config = nb.user_conf('admin')
            sage: config['max_history_length']
            1000
            sage: config['default_system']
            'sage'
            sage: config['autosave_interval']
            3600
            sage: config['default_pretty_print']
            False
        """
        return self.users()[username].conf()

    ##########################################################
    # Publishing worksheets
    ##########################################################
    def _initialize_worksheet(self, src, W):
        r"""
        Initialize a new worksheet from a source worksheet.

        INPUT:

        - ``src`` - a Worksheet instance; the source

        - ``W`` - a new Worksheet instance; the target
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
        Publish a user's worksheet.  This creates a new worksheet in
        the 'pub' directory with the same contents as ``worksheet``.

        INPUT:

        - ``worksheet`` - an instance of Worksheet

        - ``username`` - a string

        OUTPUT:

        - a new or existing published instance of Worksheet

        EXAMPLES::

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
        Delete the given worksheet and remove its name from the worksheet
        list.  Raise a KeyError, if it is missing.

        INPUT:

        - ``filename`` - a string
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

        -  ``username`` - a string

        This empties the trash for the given user and cleans up all files
        associated with the worksheets that are in the trash.

        EXAMPLES::

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

        - a list of strings.

        EXAMPLES: We make a new notebook with two users and two worksheets,
        then list their names::

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
        Migrate all old worksheets, i.e., those with no owner, to
        ``/pub``.  Currently, this always raises a
        NotImplementedError.
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
        return template("notebook/command_history.html", history_text = escape(self.history_text()))


    def history_with_start(self, start):
        n = len(start)
        return [x for x in self.history() if x[:n] == start]

    ##########################################################
    # Importing and exporting worksheets to files
    ##########################################################
    def export_worksheet(self, worksheet_filename, output_filename, verbose=True):
        """
        Export a worksheet, creating a sws file on the file system.

        INPUT:

        -  ``worksheet_filename`` - a string

        -  ``output_filename`` - a string

        - ``verbose`` - a bool (default: True); if True, print the tar
           command used to create the sws file.
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
        Import a worksheet with the given ``filename`` and set its
        ``owner``.  If the file extension is not txt or sws, raise a
        ValueError.

        INPUT:

        -  ``filename`` - a string

        -  ``owner`` - a string

        OUTPUT:

        -  ``worksheet`` - a newly created Worksheet instance

        EXAMPLES: We create a notebook and import a plain text worksheet
        into it.

        ::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: name = tmp_filename() + '.txt'
            sage: open(name,'w').write('foo\n{{{\n2+3\n}}}')
            sage: W = nb.import_worksheet(name, 'admin')

        W is our newly-created worksheet, with the 2+3 cell in it::

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
        r"""
        Import a plain text file as a new worksheet.

        INPUT:

        -  ``filename`` - a string; a filename that ends in .txt

        -  ``owner`` - a string; the imported worksheet's owner

        OUTPUT:

        -  a new instance of Worksheet

        EXAMPLES: We write a plain text worksheet to a file and import it
        using this function.

        ::

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
        r"""
        Import an sws format worksheet into this notebook as a new
        worksheet.  If the worksheet cannot be read, raise a
        ValueError.

        INPUT:

        - ``filename`` - a string; a filename that ends in .sws;
           internally it must be a tar'd bz2'd file.

        - ``owner`` - a string

        - ``verbose`` - a bool (default: True) if True print some the
           tar command used to extract the sws file.

        OUTPUT:

        - a new Worksheet instance

        EXAMPLES: We create a notebook, then make a worksheet from a plain
        text file first.

        ::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: name = tmp_filename() + '.txt'
            sage: open(name,'w').write('foo\n{{{\n2+3\n}}}')
            sage: W = nb.import_worksheet(name, 'admin')
            sage: W.filename()
            'admin/0'

        We then export the worksheet to an sws file.

        ::

            sage: nb.export_worksheet(W.filename(),  'tmp.sws', verbose=False)

        Now we import the sws.

        ::

            sage: nb._import_worksheet_sws('tmp.sws', 'admin', verbose=False)
            [Cell 0; in=2+3, out=]

        Yes, it's there now (as admin/2)::

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

    def plain_text_worksheet_html(self, filename, prompts=True):
        """
        Return HTML containing the plain text version of a worksheet.

        INPUT:

        - ``filename`` - a string; filename of a worksheet

        - ``prompts`` - a bool (default: True); whether to format the
          text for inclusion in docstrings

        OUTPUT:

        - a string - the worksheet's HTML representation
        """
        worksheet = self.get_worksheet_with_filename(filename)
        text = escape(worksheet.plain_text(prompts = prompts))
        return template("notebook/plain_text_worksheet.html",
                        worksheet_name = worksheet.name(),
                        worksheet_plain_text = text)

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
        Return the absolute path to the parent of this Notebook
        instance's home directory.

        OUTPUT:

        - a string
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
           <script type="text/javascript" src="/javascript_local/jquery/jquery.js"></script>
           <script type="text/javascript" src="/javascript/main.js"></script>
           <script type="text/javascript">
           var worksheet_filenames = %s;
           </script>
        """%(worksheet_filenames)

        return s

    def worksheet_html(self, filename, do_print=False):
        r"""
        Return the HTML for a given worksheet.

        INPUT:

        - ``filename`` - a string; the worksheet's filename

        - ``do_print`` - a bool (default: False); whether this is a
          printed worksheet

        OUTPUT:

        - a string - the worksheet rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.worksheet_html(W.filename())
            '\n<!D...ript type="text/javascript">cell_id_list=[0];</script>\n\n\n\n\n\n    </body>\n</html>'
        """
        worksheet = self.get_worksheet_with_filename(filename)
        return template("notebook/worksheet.html", worksheet_name = worksheet.name(),
                 worksheet_html = worksheet.html(include_title=False, do_print=do_print),
                        do_print = do_print)



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
        else: # typ must be archived
            W = [x for x in X if not (x.is_trashed(user) or x.is_active(user))]
        if search:
            W = [x for x in W if x.satisfies_search(search)]
        sort_worksheet_list(W, sort, reverse)  # changed W in place
        return W

    ##########################################################
    # Revision history for a worksheet
    ##########################################################
    def html_worksheet_revision_list(self, username, worksheet):
        r"""
        Return HTML for the revision list of a worksheet.

        INPUT:

        - ``username`` - a string

        - ``worksheet`` - an instance of Worksheet

        OUTPUT:

        - a string - the HTML for the revision list

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_worksheet_revision_list('admin', W)
            '\n<!D...seconds ago</span></td>\n    </tr>\n\n</table>\n\n\n    </body>\n</html>'
        """
        data = worksheet.snapshot_data()  # pairs ('how long ago', key)

        return template("notebook/worksheet_revision_list.html", data = data,
                        worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username,
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)


    def html_specific_revision(self, username, ws, rev):
        r"""
        Return the HTML for a specific revision of a worksheet.

        INPUT:

        - ``username`` - a string

        - ``ws`` - an instance of Worksheet

        - ``rev`` - a string containing the key of the revision

        OUTPUT:

        - a string - the revision rendered as HTML
        """
        t = time.time() - float(rev[:-4])
        time_ago = worksheet.convert_seconds_to_meaningful_time_span(t)

        filename = ws.get_snapshot_text_filename(rev)
        txt = bz2.decompress(open(filename).read())
        W = self.scratch_worksheet()
        W.delete_cells_directory()
        W.edit_save(txt)
        body_worksheet_html = W.html_worksheet_body(do_print=True, publish=True)

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

        return template("notebook/specific_revision.html", worksheet = ws,
                        worksheet_filename = ws.filename(),
                        username = username, rev = rev,
                        prev_rev = prev_rev, next_rev = next_rev,
                        time_ago = time_ago,
                        body_worksheet_html = body_worksheet_html,
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)


    def html_share(self, worksheet, username):
        r"""
        Return the HTML for the "share" page of a worksheet.

        INPUT:

        - ``username`` - a string

        - ``worksheet`` - an instance of Worksheet

        OUTPUT:

        - string - the share page's HTML representation

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_share(W, 'admin')
            '\n<!D...span class="username">Sage Users:</span>\n<span class="users">\n    \n</span>\n\n\n\n    </body>\n</html>'
        """
        U = self.users()
        other_users = [x for x, u in U.iteritems() if not u.is_guest() and not u.username() in [username, 'pub', '_sage_']]
        other_users.sort(lambda x,y: cmp(x.lower(), y.lower()))

        return template("notebook/worksheet_share.html", worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username, other_users = other_users,
                        user_is_admin = self.user(username).is_admin(),
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)


    def html_download_or_delete_datafile(self, ws, username, filename):
        r"""
        Return the HTML for the download or delete datafile page.

        INPUT:

        - ``username`` - a string

        - ``ws`` - an instance of Worksheet

        - ``filename`` - a string; the name of the file

        OUTPUT:

        - a string - the page rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_download_or_delete_datafile(W, 'admin', 'bar')
            '\n<!D...ploaded to this worksheet.</p>\n\n<hr class="usercontrol" />\n\n\n\n\n    </body>\n</html>'
        """
        path = "/home/%s/data/%s"%(ws.filename(), filename)

        worksheets = self.get_worksheets_with_viewer(username)
        active_worksheets = [worksheet for worksheet in worksheets if worksheet.is_active(username)]
        sort_worksheet_list(active_worksheets, 'name', False)

        ext = os.path.splitext(filename)[1].lower()
        file_is_image, file_is_text = False, False
        text_file_content = ""

        if ext in ['.png', '.jpg', '.gif']:
            file_is_image = True
        if ext in ['.txt', '.tex', '.sage', '.spyx', '.py', '.f', '.f90', '.c']:
            file_is_text = True
            text_file_content = open('%s/%s'%(ws.data_directory(), filename)).read()

        return template("notebook/download_or_delete_datafile.html",
                        worksheet = ws,
                        worksheet_filename = ws.filename(),
                        username = username,
                        active_worksheets = active_worksheets,
                        filename_ = filename,
                        path = path,
                        file_is_image = file_is_image,
                        file_is_text = file_is_text,
                        text_file_content = text_file_content,
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)




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
        Get the worksheet with the given filename.  If there is no
        such worksheet, raise a ``KeyError``.

        INPUT:

        - ``filename`` - a string

        OUTPUT:

        - a Worksheet instance
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
    def html_debug_window(self):
        r"""
        Return the HTML for the debug window.

        OUTPUT:

        - a string - the debug window rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.html_debug_window()
            "\n<div class='debug_window'>...</div>"
        """
        return template("notebook/debug_window.html")


    def html_plain_text_window(self, worksheet, username):
        r"""
        Return HTML for the window that displays a plain text version
        of the worksheet.

        INPUT:

        -  ``worksheet`` - a Worksheet instance

        -  ``username`` - a string

        OUTPUT:

        - a string - the plain text window rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_plain_text_window(W, 'admin')
            '\n<!D...>\n\n<pre class="plaintext" id="cell_intext" name="textfield"></pre>\n\n\n    </body>\n</html>'
        """
        plain_text = worksheet.plain_text(prompts=True, banner=False)
        plain_text = escape(plain_text).strip()

        return template("notebook/plain_text_window.html", worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username,
                        plain_text = plain_text, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

    def html_edit_window(self, worksheet, username):
        r"""
        Return HTML for a window for editing ``worksheet``.

        INPUT:

        - ``username`` - a string containing the username

        - ``worksheet`` - a Worksheet instance

        OUTPUT:

        - a string - the editing window's HTML representation

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_edit_window(W, 'admin')
            '\n<!D...Test\nsystem:sage\n\n{{{id=0|\n\n///\n}}}</textarea>\n</form>\n\n\n    </body>\n</html>'
        """
        text = worksheet.edit_text()
        text = escape(text)
        n_lines = text.count("\n")+1

        return template("notebook/edit_window.html", worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username, text = text,
                        n_lines = n_lines, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

    def html_beforepublish_window(self, worksheet, username):
        r"""
        Return HTML for the warning and decision page displayed prior
        to publishing the given worksheet.

        INPUT:

        - ``worksheet`` - an instance of Worksheet

        - ``username`` - a string

        OUTPUT:

        - a string - the pre-publication page rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_beforepublish_window(W, 'admin')
            '\n<!D...publish when changes are made</form></span>\n<br /><br /><br />\n\n\n    </body>\n</html>'
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
        return template("notebook/beforepublish_window.html", worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

    def html_afterpublish_window(self, worksheet, username, url, dtime):
        r"""
        Return HTML for a given worksheet's post-publication page.

        INPUT:

        - ``worksheet`` - an instance of Worksheet

        - ``username`` - a string

        - ``url`` - a string representing the URL of the published
          worksheet

        - ``dtime`` - an instance of time.struct_time representing the
          publishing time

        OUTPUT:

        - a string - the post-publication page rendered as HTML

        """
        from time import strftime
        time = strftime("%B %d, %Y %I:%M %p", dtime)

        return template("notebook/afterpublish_window.html", worksheet = worksheet,
                        worksheet_filename = worksheet.filename(),
                        username = username, url = url,
                        time = time, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

    def html_upload_data_window(self, ws, username):
        r"""
        Return HTML for the "Upload Data" window.

        INPUT:

        - ``worksheet`` - an instance of Worksheet

        - ``username`` - a string

        OUTPUT:

        - a string - the HTML representation of the data upload window

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_upload_data_window(W, 'admin')
            '\n<!D...orksheet_menu" value="Upload File" onClick="form.submit()...r />\n</div>\n\n\n    </body>\n</html>'
        """
        return template("notebook/upload_data_window.html", worksheet = worksheet,
                        worksheet_filename = ws.filename(),
                        username = username, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)


    def html(self, worksheet_filename=None, username='guest', show_debug=False, admin=False):
        r"""
        Return the HTML for a worksheet's index page.

        INPUT:

        - ``worksheet_filename`` - a string (default: None)

        - ``username`` - a string (default: 'guest')

        - ``show_debug`` - a bool (default: False)

        - ``admin`` - a bool (default: False)

        OUTPUT:

        - a string - the worksheet rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html(W.filename(), 'admin')
            '\n<!D...ipt type="text/javascript">worksheet_locked=false;</script>\n\n\n\n    </body>\n</html>'
        """
        if worksheet_filename is None or worksheet_filename == '':
            worksheet_filename = None
            W = None
        else:
            try:
                W = self.get_worksheet_with_filename(worksheet_filename)
            except KeyError:
                W = None

        template_page = "notebook/index.html"
        if W.docbrowser():
            template_page = "notebook/doc_page.html"

        return template(template_page, worksheet = W,
                        worksheet_filename = W.filename(),
                        worksheet_html = W.html(),
                        notebook = self, username = username,
                        show_debug = show_debug,
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

    ####################################################################
    # Configuration html.
    # In each case the settings html is a form that when submitted
    # pulls up another web page and sets the corresponding options.
    ####################################################################


    def html_worksheet_settings(self, ws, username):
        r"""
        Return the HTML for a worksheet's settings page.

        INPUT:

        - ``ws`` - an instance of Worksheet

        - ``username`` - a string

        OUTPUT:

        - a string - HTML representation of the settings page

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_worksheet_settings(W, 'admin')
            '\n<!D...lue="Cancel" name="button_cancel"/></span>\n<br /><br /><br />\n\n</form>\n\n\n    </body>\n</html>'
        """
        return template("notebook/worksheet_settings.html", worksheet = ws,
                        worksheet_filename = ws.filename(),
                        username = username, JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)

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
        r"""
        Return the HTML for the a documentation page.

        INPUT:

        - ``username`` - a string

        OUTPUT:

        - a string - the doc page rendered as HTML

        EXAMPLES::

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: nb.html_doc('admin')
            '\n<!D...c Documentation</a><br /><br />\n        <a href="/help/">Sage Notebook Howto...   </body>\n</html>'
        """
        return template("notebook/doc.html", username = username,
                        JSMATH = JSMATH,
                        JSMATH_IMAGE_FONTS = JSMATH_IMAGE_FONTS,
                        JEDITABLE_TINYMCE = JEDITABLE_TINYMCE,
                        sage_jsmath_macros = sage_jsmath_macros)


####################################################################

def load_notebook(dir, address=None, port=None, secure=None):
    """
    Load and return a notebook from a given directory.  Create a new
    one in that directory, if one isn't already there.

    INPUT:

    -  ``dir`` - a string that defines a directory name

    -  ``address`` - the address the server listens at

    -  ``port`` - the port the server listens on

    -  ``secure`` - whether the notebook is secure

    OUTPUT:

    - a Notebook instance
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
    Replace an absolute path with a relative path, if possible.
    Otherwise, return the given path.

    INPUT:

    - ``dir`` - a string containing, e.g., a directory name

    OUTPUT:

    - a string
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
    Sort a given list on a given key, in a given order.

    INPUT:

    - ``sort`` - a string; 'last_edited', 'owner', 'rating', or 'name'

    - ``reverse`` - a bool; if True, reverse the order of the sort.

    OUTPUT:

    - the sorted list
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
