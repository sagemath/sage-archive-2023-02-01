r"""
A Worksheet.

A worksheet is embedded in a webpage that is served by the Sage
server.  It is a linearly-ordered collections of numbered cells, where
a cell is a single input/output block.

The worksheet module is responsible for running calculations in a
worksheet, spawning Sage processes that do all of the actual work and
are controlled via pexpect, and reporting on results of calculations.
The state of the cells in a worksheet is stored on the filesystem (not
in the notebook pickle sobj).

AUTHOR:
    -- William Stein
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

# Import standard Python libraries that we will use below
import os
import copy
import shutil
import re
import string
import traceback
import time
import crypt
import bz2
import re

# A library that we ship with sage
import pexpect

# General sage library code
import sage.misc.remote_file as remote_file
import sage.misc.cython as cython
from   sage.structure.sage_object  import load, save
from   sage.interfaces.sage0 import Sage
from   sage.misc.preparser   import preparse_file
import sage.misc.interpreter
from   sage.misc.misc        import alarm, cancel_alarm, verbose, DOT_SAGE, walltime
import sage.server.support   as support

# Imports specifically relevant to the sage notebook
import worksheet_conf
import twist
from   cell import Cell, TextCell

# Set some constants that will be used for regular expressions below.
whitespace = re.compile('\s')  # Match any whitespace character
non_whitespace = re.compile('\S')

# Constants that control the behavior of the worksheet.
INTERRUPT_TRIES = 3    # number of times to send control-c to subprocess before giving up
INITIAL_NUM_CELLS = 1  # number of empty cells in new worksheets

WARN_THRESHOLD = 100   # The number of seconds, so if there was no activity on
                       # this worksheet for this many seconds, then editing
                       # is considered safe.  Used when multiple people are editing
                       # the same worksheet.

# The strings used to synchronized the compute subprocesses.
# WARNING:  If you make any changes to this, be sure to change the
# error line below that looks like this:
#         cmd += 'print "\\x01r\\x01e%s"'%self.synchro()
SC         = '\x01'
SAGE_BEGIN = SC + 'b'
SAGE_END   = SC + 'e'
SAGE_ERROR = SC + 'r'

# Integers that define which folder this worksheet is in
# relative to a given user.
ARCHIVED = 0
ACTIVE   = 1
TRASH    = 2

# The default is for there to be one sage session for each worksheet.
# If this is False, then there is just one global Sage session, like
# with Mathematica. The multisessin variable gets possibly changed
# when the notebook function in notebook.py is called.

multisession = True

def initialized_sage(server, ulimit):
    """
    Return one copy of a Sage compute process that has initialization
    code run.

    INPUT:
       server -- if sessions will be run via ssh on a remote account then
                 this string specifies that account (passed on to the Sage
                 pexpect interface).
       ulimit -- string; passed to the ulimit command before running
                 the subprocess

    OUTPUT:
        a pexpect interface to a local or remote copy of Sage

    EXAMPLES:
        sage: S = sage.server.notebook.worksheet.initialized_sage(None,None)
        sage: S
        Sage
    """
    # Create new pexpect interface to a Python instance
    S = Sage(server=server, ulimit=ulimit, maxread = 1, python=True, verbose_start=False)

    # Start up the subprocess but do not block when starting
    S._start(block_during_init=False)

    # Get at the underlying pexpect interface
    E = S.expect()

    # Send some code to it (nonblocking)
    E.sendline('\n')
    cmd = 'from sage.all_notebook import *;'
    cmd += 'import sage.server.support as _support_; '
    E.sendline(cmd)

    # Return our new Sage instance.
    return S

_a_sage = None
def init_sage_prestart(server, ulimit):
    """
    Set the module-scope variable _a_sage to
    an initialized sage server.

    INPUT:
        server, ulimit -- strings that are passed
        to the Sage pexpect interface constructor

    EXAMPLES:
    The _a_sage variable is initially set to None:
        sage: sage.server.notebook.worksheet._a_sage

    We call init_sage_prestart and now _a_sage is a Sage instance:
        sage: sage.server.notebook.worksheet.init_sage_prestart(None,None)
        sage: sage.server.notebook.worksheet._a_sage
        Sage
    """
    global _a_sage
    _a_sage = initialized_sage(server, ulimit)

def one_prestarted_sage(server, ulimit):
    """
    Return a Sage interface that has been initialized.

    INPUT:
        server, ulimit -- strings that are passed
        to the Sage pexpect interface constructor
    OUTPUT:
        -- an interface to a running copy of Sage

    If the global variable multisession is true, each call to
    one_prestarted_sage returns a new Sage compute instance.
    Otherwise it always returns the same instance.

    EXAMPLES:
        sage: sage.server.notebook.worksheet.one_prestarted_sage(None,None)
        Sage
        sage: sage.server.notebook.worksheet.multisession=False
        sage: sage.server.notebook.worksheet.one_prestarted_sage(None,None) is sage.server.notebook.worksheet._a_sage
        True
        sage: sage.server.notebook.worksheet.multisession=True
    """
    global _a_sage
    X = _a_sage
    if multisession:
        init_sage_prestart(server, ulimit)
    return X

import notebook as _notebook
def worksheet_filename(name, owner):
    """
    Return the relative directory name of this worksheet
    with given name and owner.

    INPUT:
        name -- string, which may have spaces and funny characters, which
                are replaced by underscores.
        owner -- string, with no spaces or funny characters

    OUTPUT:
        string

    EXAMPLES:
        sage: sage.server.notebook.worksheet.worksheet_filename('Example worksheet 3', 'sage10')
        'sage10/Example_worksheet_3'
        sage: sage.server.notebook.worksheet.worksheet_filename('Example#%&! work\\sheet 3', 'sage10')
        'sage10/Example_____work_sheet_3'
    """
    return owner + '/' + _notebook.clean_name(name)

class Worksheet:
    def __init__(self, name, dirname, notebook_worksheet_directory, system, owner,
                 docbrowser=False, pretty_print=False, auto_publish=False):
        """
        Create and initialize a new worksheet.

        INPUT:
            name    -- string; the name of this worksheet
            dirname -- string; name of the directory in which the worksheet's
                       data is stored
            notebook_worksheet_directory -- string; the directory in which the
                       notebook object that contains this worksheet
                       stores worksheets, i.e., nb.worksheet_directory().
            system -- string; 'sage', 'gp', 'singular', etc. -- the math software
                      system in which all code is evaluated by default
            owner  -- string; username of the owner of this worksheet
            docbrowser -- bool (default: False); whether this is a docbrowser worksheet
            pretty_print -- bool (default: False); whether all output is pretty printed
                            by default.

        EXAMPLES:
        We test the constructor via an indirect doctest:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test', 'admin')
            sage: W
            [Cell 0; in=, out=]
        """
        # Record the basic properties of the worksheet
        self.__system   = system
        self.__pretty_print = pretty_print
        self.__owner         = owner
        self.__viewers       = []
        self.__collaborators = []
        self.__docbrowser = docbrowser
        self.__autopublish = auto_publish

        # Initialize the cell id counter.
        self.__next_id = 0

        self.set_name(name)

        # set the directory in which the worksheet files will be stored.
        # We also add the hash of the name, since the cleaned name loses info, e.g.,
        # it could be all _'s if all characters are funny.
        filename = '%s/%s'%(owner, dirname)
        self.__filename = filename
        self.__dir = '%s/%s'%(notebook_worksheet_directory, filename)

        self.clear()
        self.save_snapshot(owner)

    def __cmp__(self, other):
        """
        We compare two worksheets.

        INPUT:
            self, other -- worksheets
        OUTPUT:
            -1,0,1 -- comparison is on the underlying filenames.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W2 = nb.create_new_worksheet('test2', 'admin')
            sage: W1 = nb.create_new_worksheet('test1', 'admin')
            sage: cmp(W1, W2)
            1
            sage: cmp(W2, W1)
            -1
        """
        try:
            return cmp(self.filename(), other.filename())
        except AttributeError:
            return cmp(type(self), type(other))

    def __repr__(self):
        """
        Return string representation of this worksheet, which is
        simply the string representation of the underlying list of
        cells.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.__repr__()
            '[Cell 0; in=, out=]'
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n5\n}}}\n{{{id=10|\n2+8\n///\n10\n}}}')
            sage: W.__repr__()
            '[Cell 0; in=2+3, out=5, Cell 10; in=2+8, out=10]'
        """
        return str(self.cell_list())

    def __len__(self):
        """
        Return the number of cells in this worksheet.

        OUTPUT:
            int

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: len(W)
            1
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n5\n}}}\n{{{id=10|\n2+8\n///\n10\n}}}')
            sage: len(W)
            2
        """
        return len(self.cell_list())

    def docbrowser(self):
        """
        Return True if this is a docbrowser worksheet.

        OUTPUT:
            bool

        EXAMPLES:
        We first create a standard worksheet for which docbrowser is of course False:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.docbrowser()
            False

        We create a worksheet for which docbrowser is True:
            sage: W = nb.create_new_worksheet('docs', 'admin', docbrowser=True)
            sage: W.docbrowser()
            True
        """
        try:
            return self.__docbrowser
        except AttributeError:
            return False

    ##########################################################
    # Configuration
    ##########################################################
    def conf(self):
        """
        Return the configuration object for this worksheet, which is
        stored in an sobj in the worksheet directory.

        OUTPUT:
            worksheet configuration object.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.conf()
            Configuration: {}
        """
        try:
            # after unpickling the worksheet this __conf attribute
            # will not be set because we explicitly delete it in
            # the __getstate__ method.
            return self.__conf
        except AttributeError:
            file = '%s/conf.sobj'%self.directory()
            if os.path.exists(file):
                C = load(file)
            else:
                C = worksheet_conf.WorksheetConfiguration()
            self.__conf = C
            return C

    ##########################################################
    # Basic properties
    ##########################################################
    def collaborators(self):
        """
        Return a (reference to the) list of the collaborators who can
        also view and modify this worksheet.

        OUTPUT:
            list

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: C = W.collaborators(); C
            []
            sage: C.append('sage')
            sage: W.collaborators()
            ['sage']
        """
        try:
            return self.__collaborators
        except AttributeError:
            self.__collaborators = []
            return self.__collaborators

    def collab(self):
        collab = [x for x in self.collaborators() if x != self.owner()]
        collaborators = ', '.join([x for x in collab])
        if len(collaborators) > 21:
            collaborators = collaborators[:21] + '...'
        return collaborators

    def set_collaborators(self, v):
        """
        Set the list of collaborators to those listed in the
        list v of strings.

        INPUT:
            v -- a list of strings

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: nb.add_user('hilbert','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.set_collaborators(['sage', 'admin', 'hilbert', 'sage'])

        Note that repeats are not added multiple times and admin --
        the owner -- isn't added:
            sage: W.collaborators()
            ['hilbert', 'sage']
        """
        n = self.notebook()
        U = n.users().keys()
        L = [x.lower() for x in U]
        owner = self.owner()
        self.__collaborators = []
        for x in v:
            y = x.lower()
            try:
                i = L.index(y)
                z = U[i]
                if z != owner and z not in self.__collaborators:
                    self.__collaborators.append(z)
            except ValueError:
                pass
        self.__collaborators.sort()

    def viewers(self):
        """
        Return list of viewers of this worksheet.

        OUTPUT:
           list -- of string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: nb.add_user('hilbert','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.add_viewer('hilbert')
            sage: W.viewers()
            ['hilbert']
            sage: W.add_viewer('sage')
            sage: W.viewers()
            ['hilbert', 'sage']
        """
        try:
            return self.__viewers
        except AttributeError:
            self.__viewers = []
            return self.__viewers

    def delete_notebook_specific_data(self):
        """
        Delete data from this worksheet this is specific to a certain
        notebook.  This means deleting the attached files,
        collaborators, and viewers.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('hilbert','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('test1', 'admin')
            sage: W.add_viewer('hilbert')
            sage: W.delete_notebook_specific_data()
            sage: W.viewers()
            []
            sage: W.add_collaborator('hilbert')
            sage: W.collaborators()
            ['admin', 'hilbert']
            sage: W.delete_notebook_specific_data()
            sage: W.collaborators()
            ['admin']
        """
        self.__attached = {}
        self.__collaborators = [self.owner()]
        self.__viewers = []

    def name(self):
        """
        Return the name of this worksheet.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.name()
            'A Test Worksheet'
            sage: W.set_name('Title<script>alert("Uh oh");</script>')
            sage: W.name()
            'Title&amp;lt;script&amp;gt;alert&#40;&amp;quot;Uh oh&amp;quot;&#41;;&amp;lt;/script&amp;gt;'

        """
        name = self.__name
        for chara, replacement in [('<', '&lt;'), ('>', '&gt;'), ('"', '&quot;'),
                                   ('&', '&amp;'), ('\'', '&apos;'), ('(', '&#40;'),
                                   (')', '&#41;'), ('{', '&#123;'), ('\'', '&#124;')]:
            name = name.replace(chara, replacement)
        return name

    def set_name(self, name):
        """
        Set the name of this worksheet.

        INPUT:
            name -- string

        EXAMPLES:
        We create a worksheet and change the name:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.set_name('A renamed worksheet')
            sage: W.name()
            'A renamed worksheet'
        """
        if len(name.strip()) == 0:
            name = 'Untitled'
        self.__name = name

    def set_filename_without_owner(self, nm):
        r"""
        Set this worksheet filename (actually directory) by getting
        the owner from the pre-stored owner via \code{self.owner()}.

        INPUT:
            nm -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.filename()
            'admin/0'
            sage: W.set_filename_without_owner('5')
            sage: W.filename()
            'admin/5'
        """
        filename = '%s/%s'%(self.owner(), nm)
        self.set_filename(filename)

    def set_filename(self, filename):
        """
        Set the worksheet filename  (actually directory).

        INPUT:
           filename -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.filename()
            'admin/0'
            sage: W.set_filename('admin/10')
            sage: W.filename()
            'admin/10'
        """
        old_filename = self.__filename
        self.__filename = filename
        self.__dir = '%s/%s'%(self.notebook().worksheet_directory(), filename)
        self.notebook().change_worksheet_key(old_filename, filename)

    def filename(self):
        """
        Return the filename (really directory) where the files associated
        to this worksheet are stored.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.filename()
            'admin/0'
            sage: sorted(os.listdir(nb.directory() + '/worksheets/' + W.filename()))
            ['snapshots', 'worksheet.txt']
        """
        return self.__filename

    def filename_without_owner(self):
        """
        Return the part of the worksheet filename after the last /,
        i.e., without any information about the owner of this
        worksheet.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.filename_without_owner()
            '0'
            sage: W.filename()
            'admin/0'
        """
        return os.path.split(self.__filename)[-1]

    def directory(self):
        """
        Return the full path to the directory where this
        worksheet is stored.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.directory()
            '.../worksheets/admin/0'
        """
        return self.__dir

    def data_directory(self):
        """
        Return path to directory where worksheet data is stored.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.data_directory()
            '.../worksheets/admin/0/data/'
        """
        d = self.directory() + '/data/'
        if not os.path.exists(d):
            os.makedirs(d)
        return d

    def attached_data_files(self):
        """
        Return a list of the filenames of files in the worksheet data directory.

        OUTPUT:
            list of strings

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.attached_data_files()
            []
            sage: open('%s/foo.data'%W.data_directory(),'w').close()
            sage: W.attached_data_files()
            ['foo.data']
        """
        D = self.data_directory()
        if not os.path.exists(D):
            return []
        return os.listdir(D)

    def cells_directory(self):
        """
        Return the directory in which the cells of this worksheet
        are evaluated.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.cells_directory()
            '.../worksheets/admin/0/cells/'
        """
        return self.directory() + '/cells/'

    def notebook(self):
        """
        Return the notebook that contains this worksheet.

        OUTPUT:
            a Notebook object.

        EXAMPLES:
        This really returns the Notebook object that is set as a global
        variable of the twist module.

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.notebook()
            <class 'sage.server.notebook.notebook.Notebook'>
            sage: W.notebook() is sage.server.notebook.twist.notebook
            True
        """
        return twist.notebook

    def DIR(self):
        """
        Return the absolute path to the directory that contains
        the Sage Notebook directory for the notebook that contains
        this worksheet.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.DIR()   # random output
            '/Users/was/.sage/temp/teragon_2.local/19129'
        """
        return self.notebook().DIR()

    def system(self):
        """
        Return the math software system in which by default all input
        to the notebook is evaluated.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.system()
            'sage'
            sage: W.set_system('mathematica')
            sage: W.system()
            'mathematica'
        """
        try:
            return self.__system
        except AttributeError:
            self.__system = 'sage'
            return 'sage'

    def set_system(self, system='sage'):
        """
        Set the math software system in which input is evaluated by
        default.

        INPUT:
            sysem -- string (default: 'sage')

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.set_system('magma')
            sage: W.system()
            'magma'
        """
        self.__system = system.strip()

    def pretty_print(self):
        """
        Return True if output shold be pretty printed by default.

        OUTPUT:
            bool -- True of False

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.pretty_print()
            False
            sage: W.set_pretty_print('true')
            sage: W.pretty_print()
            True
        """
        try:
            return self.__pretty_print
        except AttributeError:
            self.__pretty_print = False
            return self.__pretty_print

    def set_pretty_print(self, check='false'):
        """
        Set whether or not output should be pretty printed by default.

        INPUT:
            check -- string (default: 'false'); either 'true' or 'false'.

        NOTE: The reason the input is a string and lower case instead
        of a Python bool is because this gets called indirectory from
        javascript.  (And, Jason Grout wrote this and didn't realize
        how unpythonic this design is -- it should be redone to use
        True/False.)

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('A Test Worksheet', 'admin')
            sage: W.set_pretty_print('false')
            sage: W.pretty_print()
            False
            sage: W.set_pretty_print('true')
            sage: W.pretty_print()
            True
        """
        if check == 'false':
            check=False
        else:
            check=True
        self.__pretty_print = check
        self.eval_asap_no_output("pretty_print_default(%r)"%(check))

    ##########################################################
    # Publication
    ##########################################################
    def is_auto_publish(self):
        """
        Returns boolean of "Is this worksheet set to be published automatically when saved?"
        if private variable "autopublish" is set otherwise False is returned and the variable
        is set to False.
        """
        try:
            return self.__autopublish
        except AttributeError:
            self.__autopublish = False
            return False

    def set_auto_publish(self):
        """
        Sets the worksheet to be published automatically when the worksheet is saved if the worksheet
        isn't already set to this otherwise it is set not to.
        """
        try:
            self.__autopublish = False if self.__autopublish else True
        except AttributeError:
            self.__autopublish = True

    def is_published(self):
        """
        Return True if this worksheet is a published worksheet.

        OUTPUT:
            bool -- whether or not owner is 'pub'

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.is_published()
            False
            sage: W.set_owner('pub')
            sage: W.is_published()
            True
        """
        return self.owner() == 'pub'

    def worksheet_that_was_published(self):
        """
        Return the worksheet that was published to get this worksheet,
        if this worksheet was published.  Otherwise just return this
        worksheet.

        OUTPUT:
            Worksheet

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.worksheet_that_was_published() is W
            True

            sage: S = nb.publish_worksheet(W, 'admin')
            sage: S.worksheet_that_was_published() is S
            False
            sage: S.worksheet_that_was_published() is W
            True
        """
        try:
            return self.__worksheet_came_from
        except AttributeError:
            return self

    def publisher(self):
        """
        Return username of user that published this worksheet.

        OUTPUT:
            string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: S = nb.publish_worksheet(W, 'admin')
            sage: S.publisher()
            'admin'
        """
        return self.worksheet_that_was_published().owner()

    def is_publisher(self, username):
        """
        Return True if username is the username of the publisher
        of this worksheet, assuming this worksheet was published.

        INPUT:
            username -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: P = nb.publish_worksheet(W, 'admin')
            sage: P.is_publisher('hearst')
            False
            sage: P.is_publisher('admin')
            True
        """
        return self.publisher() == username

    def has_published_version(self):
        """
        Return True if there is a published version of this worksheet.

        OUTPUT:
            bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: P = nb.publish_worksheet(W, 'admin')
            sage: P.has_published_version()
            False
            sage: W.has_published_version()
            True
        """
        try:
            self.published_version()
            return True
        except ValueError:
            return False

    def set_published_version(self, filename):
        """
        Set the published version of this worksheet to be the
        worksheet with given filename.

        INPUT:
            filename -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: P = nb.publish_worksheet(W, 'admin')  # indirect test
            sage: W._Worksheet__published_version
            'pub/0'
            sage: W.set_published_version('pub/0')
        """
        self.__published_version = filename

    def published_version(self):
        """
        If this worksheet was published, return the published version
        of this worksheet.  Otherwise, raise a ValueError.

        OUTPUT:
            a worksheet (or raise a ValueError)

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: P = nb.publish_worksheet(W, 'admin')
            sage: W.published_version() is P
            True
        """
        try:
            filename =self.__published_version
            try:
                W = self.notebook().get_worksheet_with_filename(filename)
                return W
            except KeyError:
                del self.__published_version
                raise ValueError
        except AttributeError:
            raise ValueError, "no published version"

    def set_worksheet_that_was_published(self, W):
        """
        Set the worksheet that was published to get self to W.

        INPUT:
            W -- a Worksheet

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: P = nb.publish_worksheet(W, 'admin')
            sage: P.worksheet_that_was_published() is W
            True

        We fake things and make it look like P published itself:
            sage: P.set_worksheet_that_was_published(P)
            sage: P.worksheet_that_was_published() is P
            True
        """
        if not isinstance(W, Worksheet):
            raise TypeError, "W must be a worksheet"
        self.__worksheet_came_from = W

    def rate(self, x, comment, username):
        """
        Set the rating on this worksheet by the given user to x and
        also set the given comment.

        INPUT:
            x -- integer
            comment -- string
            usename -- string

        EXAMPLES:
        We create a worksheet and rate it, then look at the ratings.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.rate(3, 'this is great', 'hilbert')
            sage: W.ratings()
            [('hilbert', 3, 'this is great')]

        Note that only the last rating by a user counts:
            sage: W.rate(1, 'this lacks content', 'riemann')
            sage: W.rate(0, 'this lacks content', 'riemann')
            sage: W.ratings()
            [('hilbert', 3, 'this is great'), ('riemann', 0, 'this lacks content')]
        """
        r = self.ratings()
        for i in range(len(r)):
            if r[i][0] == username:
                r[i] = (username, x, comment)
                return
        else:
            r.append((username, x, comment))

    def is_rater(self, username):
        """
        Return True is the user with given username has rated this worksheet.

        INPUT:
            username -- string

        OUTPUT:
            bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.rate(0, 'this lacks content', 'riemann')
            sage: W.is_rater('admin')
            False
            sage: W.is_rater('riemann')
            True
        """
        try:
            return username in [x[0] for x in self.ratings()]
        except TypeError:
            return False

    def ratings(self):
        """
        Return all the ratings of this worksheet.

        OUTPUT:
            list -- a reference to the list of ratings.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.ratings()
            []
            sage: W.rate(0, 'this lacks content', 'riemann')
            sage: W.rate(3, 'this is great', 'hilbert')
            sage: W.ratings()
            [('riemann', 0, 'this lacks content'), ('hilbert', 3, 'this is great')]
        """
        try:
            return self.__ratings
        except AttributeError:
            v = []
            self.__ratings = v
            return v

    def html_ratings_info(self):
        """
        Return html that renders to give a summary of how this
        worksheet has been rated.

        OUTPUT:
           string -- a string of HTML as a bunch of table rows.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.rate(0, 'this lacks content', 'riemann')
            sage: W.rate(3, 'this is great', 'hilbert')
            sage: W.html_ratings_info()
            '<tr><td>hilbert</td><td align=center>3</td><td>this is great</td></tr>\n<tr><td>riemann</td><td align=center>0</td><td>this lacks content</td></tr>'
        """
        ratings = self.ratings()
        lines = []
        for z in sorted(ratings):
            if len(z) == 2:
                person, rating = z
                comment = ''
            else:
                person, rating, comment = z
            lines.append('<tr><td>%s</td><td align=center>%s</td><td>%s</td></tr>'%(
                person, rating, '&nbsp;' if not comment else comment))
        return '\n'.join(lines)

    def rating(self):
        """
        Return overall aerage rating of self.

        OUTPUT:
            float or the int -1 to mean "not rated"

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.rating()
            -1
            sage: W.rate(0, 'this lacks content', 'riemann')
            sage: W.rate(3, 'this is great', 'hilbert')
            sage: W.rating()
            1.5
        """
        r = [x[1] for x in self.ratings()]
        if len(r) == 0:
            rating = -1    # means "not rated"
        else:
            rating = float(sum(r))/float(len(r))
        return rating

    ##########################################################
    # Active, trash can and archive
    ##########################################################
    def everyone_has_deleted_this_worksheet(self):
        """
        Return True if all users have deleted this worksheet,
        so we know we can safely purge it from disk.

        OUTPUT:
           bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.everyone_has_deleted_this_worksheet()
            False
            sage: W.move_to_trash('admin')
            sage: W.everyone_has_deleted_this_worksheet()
            True
        """
        for user in self.__collaborators + [self.owner()]:
            # When the worksheet has been deleted by the owner,
            # self.owner() returns None, so we have to be careful
            # about that case.
            if user is not None and not self.is_trashed(user):
                return False
        return True

    def user_view(self, user):
        """
        Return the view that the given user has of this worksheet.  If
        the user currently doesn't have a view set it to ACTIVE and
        return ACTIVE.

        INPUT:
            user -- a string

        OUTPUT:
            Python int -- one of ACTIVE, ARCHIVED, TRASH, which are
            defined in worksheet.py

        EXAMPLES:
        We create a new worksheet and get the view, which is ACTIVE:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.user_view('admin')
            1
            sage: sage.server.notebook.worksheet.ACTIVE
            1

        Now for the admin user we move W to the archive:
            sage: W.move_to_archive('admin')

        The view is now archive.
            sage: W.user_view('admin')
            0
            sage: sage.server.notebook.worksheet.ARCHIVED
            0

        For any other random viewer the view is set by default
        to ACTIVE.
            sage: W.user_view('foo')
            1
        """
        try:
            return self.__user_view[user]
        except AttributeError:
            self.__user_view = {}
        except KeyError:
            pass
        self.__user_view[user] = ACTIVE
        return ACTIVE

    def set_user_view(self, user, x):
        """
        Set the view on this worksheet for the given user.

        INPUT:
            user -- a string
            x -- int, one of the variables ACTIVE, ARCHIVED, TRASH in worksheet.py

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.set_user_view('admin', sage.server.notebook.worksheet.ARCHIVED)
            sage: W.user_view('admin') == sage.server.notebook.worksheet.ARCHIVED
            True
        """
        try:
            self.__user_view[user] = x
        except (KeyError, AttributeError):
            self.user_view(user)
            self.__user_view[user] = x

    def user_view_is(self, user, x):
        """
        Return True if the user view of user is x.

        INPUT:
            user -- a string
            x -- int, one of the variables ACTIVE, ARCHIVED, TRASH in worksheet.py

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Publish Test', 'admin')
            sage: W.user_view_is('admin', sage.server.notebook.worksheet.ARCHIVED)
            False
            sage: W.user_view_is('admin', sage.server.notebook.worksheet.ACTIVE)
            True
        """
        return self.user_view(user) == x

    def is_archived(self, user):
        """
        Return True if this worksheet is archived for the given user.

        INPUT:
            user -- string
        OUTPUT:
            bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Archived Test', 'admin')
            sage: W.is_archived('admin')
            False
            sage: W.move_to_archive('admin')
            sage: W.is_archived('admin')
            True
        """
        return self.user_view_is(user, ARCHIVED)

    def is_active(self, user):
        """
        Return True if this worksheet is active for the given user.

        INPUT:
            user -- string
        OUTPUT:
            bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Active Test', 'admin')
            sage: W.is_active('admin')
            True
            sage: W.move_to_archive('admin')
            sage: W.is_active('admin')
            False
        """
        return self.user_view_is(user, ACTIVE)

    def is_trashed(self, user):
        """
        Return True if this worksheet is in the trash for the given user.

        INPUT:
            user -- string
        OUTPUT:
            bool

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Trash Test', 'admin')
            sage: W.is_trashed('admin')
            False
            sage: W.move_to_trash('admin')
            sage: W.is_trashed('admin')
            True
        """
        return self.user_view_is(user, TRASH)

    def move_to_archive(self, user):
        """
        Move this worksheet to be archived for the given user.

        INPUT:
            user -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Archive Test', 'admin')
            sage: W.move_to_archive('admin')
            sage: W.is_archived('admin')
            True
        """
        self.set_user_view(user, ARCHIVED)

    def set_active(self, user):
        """
        Set his worksheet to be active for the given user.

        INPUT:
            user -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Active Test', 'admin')
            sage: W.move_to_archive('admin')
            sage: W.is_active('admin')
            False
            sage: W.set_active('admin')
            sage: W.is_active('admin')
            True
        """
        self.set_user_view(user, ACTIVE)

    def move_to_trash(self, user):
        """
        Move this worksheet to the trash for the given user.

        INPUT:
            user -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Trash Test', 'admin')
            sage: W.move_to_trash('admin')
            sage: W.is_trashed('admin')
            True
        """
        self.set_user_view(user, TRASH)

    def move_out_of_trash(self, user):
        """
        Exactly the same as set_active(user).

        INPUT:
           user -- string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Active Test', 'admin')
            sage: W.move_to_trash('admin')
            sage: W.is_active('admin')
            False
            sage: W.move_out_of_trash('admin')
            sage: W.is_active('admin')
            True
        """
        self.set_active(user)

    #############

    def delete_cells_directory(self):
        """
        Delete the directory in which all the cell computations occur.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\n{{{\n3^20\n}}}')
            sage: sorted(os.listdir(W.directory()))
            ['snapshots', 'worksheet.txt']
            sage: W.cell_list()[0].evaluate()
            sage: sorted(os.listdir(W.directory()))
            ['cells', 'code', 'data', 'snapshots', 'worksheet.txt']
            sage: W.delete_cells_directory()
            sage: sorted(os.listdir(W.directory()))
            ['code', 'data', 'snapshots', 'worksheet.txt']
        """
        dir = self.cells_directory()
        if os.path.exists(dir):
            shutil.rmtree(dir)


    ##########################################################
    # Owner/viewer/user management
    ##########################################################

    def owner(self):
        try:
            return self.__owner
        except AttributeError:
            self.__owner = 'pub'
            return 'pub'

    def is_owner(self, username):
        return self.owner() == username

    def set_owner(self, owner):
        self.__owner = owner
        if not owner in self.__collaborators:
            self.__collaborators.append(owner)

    def user_is_only_viewer(self, user):
        try:
            return user in self.__viewers
        except AttributeError:
            return False

    def user_is_viewer(self, user):
        try:
            return user in self.__viewers or user in self.__collaborators or user == self.publisher()
        except AttributeError:
            return True

    def user_is_collaborator(self, user):
        try:
            return user in self.__collaborators
        except AttributeError:
            return True

    def user_can_edit(self, user):
        """
        Return True if the user with given name is allowed to edit this worksheet.

        INPUT:
            user -- string
        OUTPUT:
            bool

        EXAMPLES:
        We create a notebook with one worksheet and two users.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: nb.add_user('william', 'william', 'wstein@sagemath.org', force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.user_can_edit('sage')
            True

        At first the user 'william' can't edit this worksheet:
            sage: W.user_can_edit('william')
            False

        After adding 'william' as a collaborator he can edit the worksheet.
            sage: W.add_collaborator('william')
            sage: W.user_can_edit('william')
            True

        Clean up:
            sage: nb.delete()
        """
        return self.user_is_collaborator(user) or self.is_owner(user)

    def delete_user(self, user):
        """
        Delete a user from having any view or ownership of this worksheet.

        INPUT:
            user -- string; the name of a user

        EXAMPLES:
        We create a notebook with 2 users and 1 worksheet that both view.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('wstein','sage','wstein@sagemath.org',force=True)
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.new_worksheet_with_title_from_text('Sage', owner='sage')
            sage: W.add_viewer('wstein')
            sage: W.owner()
            'sage'
            sage: W.viewers()
            ['wstein']

        We delete the sage user from the worksheet W.   This makes wstein the
        new owner.
            sage: W.delete_user('sage')
            sage: W.viewers()
            ['wstein']
            sage: W.owner()
            'wstein'

        Then we delete wstein from W, which makes the owner None:
            sage: W.delete_user('wstein')
            sage: W.owner() is None
            True
            sage: W.viewers()
            []

        Finally, we clean up.
            sage: nb.delete()
        """
        if user in self.__collaborators:
            self.__collaborators.remove(user)
        if user in self.__viewers:
            self.__viewers.remove(user)
        if self.__owner == user:
            if len(self.__collaborators) > 0:
                self.__owner = self.__collaborators[0]
            elif len(self.__viewers) > 0:
                self.__owner = self.__viewers[0]
            else:
                # Now there is nobody to take over ownership.  We
                # assign the owner None, which means nobody owns it.
                # It will get purged elsewhere.
                self.__owner = None


    def add_viewer(self, user):
        """
        Add the given user as an allowed viewer of this worksheet.

        INPUT:
            user -- string (username)

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('diophantus','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Viewer test', 'admin')
            sage: W.add_viewer('diophantus')
            sage: W.viewers()
            ['diophantus']
        """
        try:
            if not user in self.__viewers:
                self.__viewers.append(user)
        except AttributeError:
            self.__viewers = [user]

    def add_collaborator(self, user):
        """
        Add the given user as a collaborator on this worksheet.

        INPUT:
            user -- a string

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('diophantus','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Collaborator test', 'admin')
            sage: W.collaborators()
            []
            sage: W.add_collaborator('diophantus')
            sage: W.collaborators()
            ['diophantus']
        """
        try:
            if not user in self.__collaborators:
                self.__collaborators.append(user)
        except AttributeError:
            self.__collaborators = [user]


    ##########################################################
    # Searching
    ##########################################################
    def satisfies_search(self, search):
        """
        Return True if all words in search are in the saved text of
        the worksheet.

        INPUT:
            search is a string that describes a search query, i.e.,
            a space-separated collections of words.

        OUTPUT:
            True if the search is satisfied by self, i.e., all
            the words appear in the text version of self.
        """
        # Load the worksheet data file from disk.
        r = open('%s/worksheet.txt'%self.__dir).read().lower()
        # Check that every single word is in the file from disk.
        for W in split_search_string_into_keywords(search):
            if W.lower() not in r:
                # Some word from the text is not in the search list, so
                # we return False.
                return False
        # Every single word is there.
        return True


    ##########################################################
    # Saving
    ##########################################################

    def save(self):
        path = self.__dir
        E = self.edit_text()
        self.save_snapshot(self.owner(), E)
        save(self.conf(), path + '/conf.sobj')

    def save_snapshot(self, user, E=None):
        self.uncache_snapshot_data()
        path = self.snapshot_directory()
        basename = str(int(time.time()))
        filename = '%s/%s.bz2'%(path, basename)
        if E is None:
            E = self.edit_text()
        open(filename, 'w').write(bz2.compress(E))
        open('%s/worksheet.txt'%self.__dir, 'w').write(E)
        try:
            X = self.__saved_by_info
        except AttributeError:
            X = {}
            self.__saved_by_info = X
        X[basename] = user
        if self.is_auto_publish():
            self.notebook().publish_worksheet(self, user)

    def get_snapshot_text_filename(self, name):
        path = self.snapshot_directory()
        return '%s/%s'%(path, name)

    def user_autosave_interval(self, username):
        return self.notebook().user(username)['autosave_interval']

    def autosave(self, username):
        try:
            last = self.__last_autosave
        except AttributeError:
            self.__last_autosave = time.time()
            return
        t = time.time()
        if t - last >= self.user_autosave_interval(username):
            self.__last_autosave = t
            self.save_snapshot(username)

    def revert_to_snapshot(self, name):
        path = self.snapshot_directory()
        filename = '%s/%s.txt'%(path, name)
        E = bz2.decompress(open(filename).read())
        self.edit_save(E)

    def _saved_by_info(self, x):
        try:
            u = self.__saved_by_info[x]
            return ' ago by %s'%u
        except (KeyError,AttributeError):
            return ' ago'

    def snapshot_data(self):
        try:
            return self.__snapshot_data
        except AttributeError:
            pass
        filenames = os.listdir(self.snapshot_directory())
        filenames.sort()
        t = time.time()
        v = [(convert_seconds_to_meaningful_time_span(t - float(os.path.splitext(x)[0]))+ self._saved_by_info(x), x)  \
             for x in filenames]
        self.__snapshot_data = v
        return v

    def uncache_snapshot_data(self):
        try:
            del self.__snapshot_data
        except AttributeError:
            pass

    def revert_to_last_saved_state(self):
        filename = '%s/worksheet.txt'%(self.__dir)
        E = open(filename).read()
        self.edit_save(E)

    def snapshot_directory(self):
        path = os.path.abspath(self.__dir) + '/snapshots/'
        if not os.path.exists(path):
            os.makedirs(path)
        return path


    ##########################################################
    # Stuff to customize how pickling works slightly.
    ##########################################################

    def __getstate__(self):
        """
        The getstate method makes sure that the self.__cells
        dictionary is not saved in the pickle since it could be huge.

        OUTPUT:
            a dictionary; same as self.__dict__ but with some fields deleted.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test Edit Save', 'admin')
            sage: v = W.__getstate__().keys(); v.sort(); v
            ['_Worksheet__autopublish', '_Worksheet__collaborators', '_Worksheet__comp_is_running', '_Worksheet__dir', '_Worksheet__docbrowser', '_Worksheet__filename', '_Worksheet__name', '_Worksheet__next_id', '_Worksheet__owner', '_Worksheet__pretty_print', '_Worksheet__queue', '_Worksheet__saved_by_info', '_Worksheet__system', '_Worksheet__viewers']
        """
        d = copy.copy(self.__dict__)

        # These attributes can take a while too and there is no need to cache them
        for attr in ['html', 'notebook', 'conf']:
            mangled = '_Worksheet__%s'%attr
            if d.has_key(mangled):
                del d[mangled]

        if d.has_key('_Worksheet__cells'):
            try:
                #print "Saving ", self.directory()
                self.save()  # make sure the worksheet.txt file is up to date.
                del d['_Worksheet__cells']
            except:
                # It is important to catch all exceptions.  If
                # *anything* goes wrong here we must catch it or the
                # whole notebook sobj could get messed up,
                # potentially, in theory, maybe.

                # There is one possible easy fix related to absolute paths in old old versions
                # of the notebook.  Check for this.
                try:
                    dir = self.directory()
                    i = dir.find('/worksheets/')
                    if i != -1:
                        i = dir[:i].rfind('/')
                        if i != -1:
                            self.__dir = dir[i+1:]
                            d['_Worksheet__dir'] = dir[i+1:]
                            self.save()
                            del d['_Worksheet__cells']
                            #print "Saving worksheet %s -- getting rid of absolute path."%self.directory()
                            return d
                except Exception, msg:
                    print msg

                    print "Unable to save worksheet %s; you may have a permissions or other problem that could result in data loss."%self.directory()
        return d

    # The following setstate method is here
    # so that when this object is pickled and
    # unpickled, the self.__sage attribute
    # will not be set, so it will properly initialized.
    def __setstate__(self, state):
        self.__dict__ = state
        try:
            del self.__sage
            self.__queue = []
        except AttributeError:
            pass

    ##########################################################
    # Exporting the worksheet in plain text command-line format
    ##########################################################
    def plain_text(self, prompts=False, banner=True):
        """
        Return a plain-text version of the worksheet.

        INPUT:
            prompts -- if True format for inclusion in docstrings.
        """
        s = ''
        if banner:
            s += "#"*80 + '\n'
            s += "# Worksheet: %s"%self.name() + '\n'
            s += "#"*80+ '\n\n'

        for C in self.cell_list():
            t = C.plain_text(prompts=prompts).strip('\n')
            if t != '':
                s += '\n' + t
        return s

    def input_text(self):
        """
        Return text version of the input to the worksheet.
        """
        return '\n\n---\n\n'.join([C.input_text() for C in self.cell_list()])

    ##########################################################
    # Editing the worksheet in plain text format (export and import)
    ##########################################################
    def edit_text(self):
        """
        Returns a plain-text version of the worksheet with \{\{\{\}\}\} wiki-formatting,
        suitable for hand editing.
        """
        s = self.name() + '\n'
        s += 'system:%s'%self.system()

        for C in self.cell_list():
            t = C.edit_text().strip()
            if t != '':
                s += '\n\n' + t
        return s

    def reset_interact_state(self):
        """
        Reset the interact state of this worksheet.
        """
        try:
            S = self.__sage
        except AttributeError:
            return
        try:
            S._send('sage.server.notebook.interact.reset_state()')
        except OSError:
            # Dosn't matter, since if S is not running, no need
            # to zero out the state dictionary.
            return

    def edit_save(self, text, ignore_ids=False):
        """
        Set the contents of this worksheet to the worksheet defined by
        the plain text string text, which should be a sequence of html
        and {{{}}}'s code blocks.

        INPUT:
            text -- a string
            ignore_ids -- bool (default: False); if True ignore all the
                          id's in the {{{}}} code block.

        EXAMPLES:
        We create a new test notebook and a worksheet.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test Edit Save', 'sage')

        We set the contents of the worksheet using the edit_save command.
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n5\n}}}\n{{{\n2+8\n///\n10\n}}}')
            sage: W
            [Cell 0; in=2+3, out=5, Cell 1; in=2+8, out=10]
            sage: W.name()
            'Sage'
        """
        # Clear any caching.
        try:
            del self.__html
        except AttributeError:
            pass

        self.reset_interact_state()

        text.replace('\r\n','\n')
        name, i = extract_name(text)
        self.set_name(name)
        text = text[i:]

        system, i = extract_system(text)
        if system == "None":
            system = "sage"
        self.set_system(system)
        text = text[i:]

        data = []
        while True:
            plain_text = extract_text_before_first_compute_cell(text).strip()
            if len(plain_text) > 0:
                T = plain_text
                data.append(('plain', T))
            try:
                meta, input, output, i = extract_first_compute_cell(text)
                data.append(('compute', (meta,input,output)))
            except EOFError, msg:
                print msg
                break
            text = text[i:]

        ids = set([x[0]['id'] for typ, x in data if typ == 'compute' and  x[0].has_key('id')])
        used_ids = set([])

        cells = []
        for typ, T in data:
            if typ == 'plain':
                if len(T) > 0:
                    id = next_available_id(ids)
                    ids.add(id)
                    cells.append(self._new_text_cell(T, id=id))
                    used_ids.add(id)
            elif typ == 'compute':
                meta, input, output = T
                if not ignore_ids and meta.has_key('id'):
                    id = meta['id']
                    if id in used_ids:
                        # In this case don't reuse, since ids must be unique.
                        id = next_available_id(ids)
                        ids.add(id)
                    html = True
                else:
                    id = next_available_id(ids)
                    ids.add(id)
                    html = False
                used_ids.add(id)
                if hasattr(self, '__cells'):
                    C = self.get_cell_with_id(id = id)
                    if isinstance(C, TextCell):
                        C = self._new_cell(id)
                else:
                    C = self._new_cell(id)
                C.set_input_text(input)
                C.set_output_text(output, '')
                if html:
                    C.update_html_output(output)
                cells.append(C)

        if len(cells) == 0:   # there must be at least one cell.
            cells = [self._new_cell()]
        elif isinstance(cells[-1], TextCell):
            cells.append(self._new_cell())

        self.__cells = cells

        if not self.is_published():
            for c in self.cell_list():
                if c.is_interactive_cell():
                    if not c in self.__queue:
                        self.enqueue(c)

        # This *depends* on self.cell_list() being set!!
        self.set_cell_counter()



    ##########################################################
    # HTML rendering of the whole worksheet
    ##########################################################
    def html(self, include_title=True, do_print=False,
             confirm_before_leave=False, read_only=False):
        if self.is_published():
            try:
                return self.__html
            except AttributeError:
                s = self.html_worksheet_body(do_print=True)
                s += self.javascript_for_jsmath_rendering()
                self.__html = s
                return s

        s = ''

        s += self.html_worksheet_body(do_print=do_print)

        if do_print:
            s += self.javascript_for_jsmath_rendering()
        else:
            s += self.javascript_for_being_active_worksheet()

        if not do_print and confirm_before_leave:
            s += self.javascript_confirm_before_leave()

        return s

    def truncated_name(self, max=30):
        name = self.name()
        if len(name) > max:
            name = name[:max] + ' ...'
        return name

    def html_title(self, username='guest'):
        name = self.truncated_name()

        warn = self.warn_about_other_person_editing(username, WARN_THRESHOLD)

        s = ''
        s += '<div class="worksheet_title">'
        s += '<a id="worksheet_title" class="worksheet_title" onClick="rename_worksheet(); return false;" title="Click to rename this worksheet">%s</a>'%(name.replace('<','&lt;'))
        s += '<br>' + self.html_time_last_edited()
        if warn and username != 'guest' and not self.is_doc_worksheet():
            s += '&nbsp;&nbsp;<span class="pingdown">(Someone else is viewing this worksheet)</span>'
        s += '</div>'

        return s

    def is_doc_worksheet(self):
        try:
            return self.__is_doc_worksheet
        except AttributeError:
            return False

    def set_is_doc_worksheet(self, value):
        self.__is_doc_worksheet = value

    def html_save_discard_buttons(self):
        if self.is_doc_worksheet():
            return ''
        return """
        <button name="button_save" title="Save changes" onClick="save_worksheet();">Save</button><button title="Save changes and close window" onClick="save_worksheet_and_close();" name="button_save">Save & quit</button><button title="Discard changes to this worksheet" onClick="worksheet_discard();">Discard & quit</button>
        """

    def html_share_publish_buttons(self, select=None, backwards=False):
        if self.is_doc_worksheet():
            return ''
        def cls(x):
            if x == select:
                return "control-select"
            else:
                return "control"
        backwards = '../' if backwards else ''
        return """

        <a  title="Print this worksheet" class="usercontrol" onClick="print_worksheet()"><img border=0 src="/images/icon_print.gif" alt="Print">Print</a>
        <a class="%s" title="Interactively use this worksheet" onClick="edit_worksheet();">Worksheet</a>
        <a class="%s" title="Edit text version of this worksheet" href="%sedit">Edit</a>
        <a class="%s" title="View plain text version of this worksheet" href="%stext">Text</a>
        <a class="%s" href="%srevisions" title="View changes to this worksheet over time">Undo</a>
        <a class="%s" href="%sshare" title="Let others edit this worksheet">Share</a>
        <a class="%s" href="%spublish" title="Make this worksheet publicly viewable">Publish</a>
        """%(cls('use'),cls('edit'),backwards,cls('text'),backwards,cls('revisions'),backwards,cls('share'),backwards,cls('publish'),backwards)

    def html_data_options_list(self):
        D = self.attached_data_files()
        D.sort()
        x = '\n'.join(['<option value="datafile?name=%s">%s</option>'%(nm,nm) for nm in D])
        return x

    def html_file_menu(self):
##  <option title="Save this worksheet as an HTML web page" onClick="save_as('html');">Save as HTML (zipped) </option>
##  <option title="Save this worksheet to LaTeX format" onClick="save_as('latex');">Save as LaTeX (zipped) </option>
##  <option title="Save this worksheet as a PDF file" onClick="save_as('pdf');">Save as PDF</option>
##  <option title="Save this worksheet as a text file" onClick="save_as('text');">Save as Text</option>

        if self.is_doc_worksheet():
            system_select = ''
            pretty_print_check = ''
        else:
            system_select = self.notebook().html_system_select_form_element(self)
            pretty_print_check = self.notebook().html_pretty_print_check_form_element(self)

        data = self.html_data_options_list()

        return """
<select class="worksheet"  onchange="go_option(this);">
<option title="Select a file related function" value=""  selected>File...</option>
 <option title="Load a new worksheet stored in a file" value="upload_worksheet_button();">Upload worksheet from a file</option>
 <option title="Create a new worksheet" value="new_worksheet();">New worksheet</option>
 <option title="Save this worksheet to an sws file" value="download_worksheet('%s');">Download to a file</option>
 <option title="Print this worksheet" value="print_worksheet();">Print</option>
 <option title="Rename this worksheet" value="rename_worksheet();">Rename worksheet</option>
 <option title="Copy this worksheet" value="copy_worksheet();">Copy worksheet</option>
 <option title="Move this worksheet to the trash" value="delete_worksheet('%s');">Delete worksheet</option>
</select>

<select class="worksheet"  onchange="go_option(this);" >
 <option title="Select a worksheet function" value="" selected>Action...</option>
 <option title="Interrupt currently running calculations, if possible" value="interrupt();">Interrupt</option>
 <option title="Restart the worksheet process" value="restart_sage();">Restart worksheet</option>
 <option title="Quit the worksheet process" value="save_worksheet_and_close();">Save and quit worksheet</option>
 <option value="">---------------------------</option>
 <option title="Evaluate all input cells in the worksheet" value="evaluate_all();">Evaluate All</option>
 <option title="Hide all output" value="hide_all();">Hide All Output</option>
 <option title="Show all output" value="show_all();">Show All Output</option>
 <option title="Delete all output" value="delete_all_output();">Delete All Output</option>
 <option value="">---------------------------</option>
 <option title="Switch to single-cell mode" value="slide_mode();">One Cell Mode</option>
 <option title="Switch to multi-cell mode" value="cell_mode();">Multi Cell Mode</option>
 </select>

<select class="worksheet" onchange="handle_data_menu(this);" >
 <option title="Select an attached file" value="" selected>Data...</option>
 <option title="Upload or create a data file in a wide range of formats" value="__upload_data_file__">Upload or create file...</option>
 <option value="">--------------------</option>
%s
</select>

 %s
 %s
 """%(_notebook.clean_name(self.name()), self.filename(),
      data, system_select, pretty_print_check)
# <option title="Browse the data directory" value="data/">Browse data directory...</option>
# <option title="Browse the directory of output from cells" value="cells/">Browse cell output directories...</option>

# <option title="Configure this worksheet" value="worksheet_settings();">Worksheet settings</option>

    def html_menu(self):
        name = self.filename()

        menu = '&nbsp;'*3 + self.html_file_menu()

        filename = os.path.split(self.filename())[-1]
        download_name = _notebook.clean_name(self.name())

        #menu += '  </span>' #why is this here?  it isn't opened anywhere.

        return menu

    def html_worksheet_body(self, do_print, publish=False):
        n = len(self.cell_list())
        published = self.is_published() or publish

        s = '<div class="cell_input_active" id="cell_resizer"></div>'
        D = self.notebook().conf()
        ncols = D['word_wrap_cols']
        if not published:
            s += '<div class="worksheet_cell_list" id="worksheet_cell_list">\n'

        for i in range(n):
            cell = self.cell_list()[i]
            s += cell.html(ncols, do_print=do_print) + '\n'

        if not do_print and not published:
            s += '\n</div>\n'
            s += '\n<div class="insert_new_cell" id="insert_last_cell" onmousedown="insert_new_cell_after(cell_id_list[cell_id_list.length-1]);"> </div>\n'
            s += '<div class="worksheet_bottom_padding"></div>\n'
        return s

    def javascript_for_being_active_worksheet(self):
        s =  '<script type="text/javascript">cell_id_list=%s;</script>'%self.compute_cell_id_list()
        #s += 'for(i=0;i<cell_id_list.length;i++) prettify_cell(cell_id_list[i]);</script>\n'
        return s

    def javascript_for_jsmath_rendering(self):
        return '<script language=javascript>jsMath.ProcessBeforeShowing();</script>\n'

    def javascript_confirm_before_leave(self):
        return """<script type="text/javascript">
            window.onbeforeunload = confirmBrowseAway;
            function confirmBrowseAway()
            {
            return "Unsubmitted cells will be lost.";
            }
            </script>
            """



    ##########################################################
    # Last edited
    ##########################################################
    def last_edited(self):
        try:
            return self.__last_edited[0]
        except AttributeError:
            t = time.time()
            self.__last_edited = (t, self.owner())
            return t

    def date_edited(self):
        """
        Returns the date the worksheet was last edited if already recorded otherwise
        the current local time is recorded and returned.
        """
        try:
            return self.__date_edited[0]
        except AttributeError:
            t = time.localtime()
            self.__date_edited = (t, self.owner())
            return t

    def last_to_edit(self):
        try:
            return self.__last_edited[1]
        except AttributeError:
            return self.owner()

    def record_edit(self, user):
        self.__last_edited = (time.time(), user)
        self.__date_edited = (time.localtime(), user)
        self.autosave(user)

    def time_since_last_edited(self):
        return time.time() - self.last_edited()

    def warn_about_other_person_editing(self,username, threshold):
        r"""
        Check to see if another user besides username was the last to
        edit this worksheet during the last \var{threshold} seconds.  If
        so, return True and that user name.  If not, return False.

        INPUT:
           username -- user who would like to edit this file.
           threshold -- number of seconds, so if there was no activity on
                   this worksheet for this many seconds, then editing is
                   considered safe.
        """
        if self.time_since_last_edited() < threshold:
            user = self.last_to_edit()
            if user != username:
                return True, user
        False

    def ws_time_since_last_edited(self):
        t = self.time_since_last_edited()
        tm = convert_seconds_to_meaningful_time_span(t)
        return tm

    def html_time_since_last_edited(self):
        tm, who = time_since_last_edited
        return '<span class="lastedit">%s ago%s</span>'%(tm, self.last_to_edit)

    def html_time_last_edited(self):
        tm = convert_time_to_string(self.last_edited())
        who = self.last_to_edit()
        t = '<span class="lastedit">last edited on %s by %s</span>'%(tm, who)
        return t

    ##########################################################
    # Managing cells and groups of cells in this worksheet
    ##########################################################

    def cell_id_list(self):
        """
        Return a new list of the id's of cells in this worksheet.

        OUTPUT:
            a new list

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test Edit Save', 'admin')

        Now we set the worksheet to have two cells with the default id
        of 0 and another with id 10.
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n5\n}}}\n{{{id=10|\n2+8\n///\n10\n}}}')
            sage: W.cell_id_list()
            [0, 10]
        """
        return [C.id() for C in self.cell_list()]

    def compute_cell_id_list(self):
        return [C.id() for C in self.cell_list() if isinstance(C, Cell)]

    def cell_list(self):
        """
        Return a reference to the list of the all the cells in this worksheet.

        OUTPUT:
            list -- a list of cells

        NOTE: This function loads the cell list from disk (the file
        worksheet.txt) if it isn't available in memory.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test Edit Save', 'admin')
            sage: W.edit_save('Sage\n{{{\n2+3\n///\n5\n}}}\n{{{\n2+8\n///\n10\n}}}')
            sage: v = W.cell_list(); v
            [Cell 0; in=2+3, out=5, Cell 1; in=2+8, out=10]
            sage: v[0]
            Cell 0; in=2+3, out=5
        """
        try:
            return self.__cells
        except AttributeError:
            worksheet_txt = '%s/worksheet.txt'%self.__dir
            if not os.path.exists(worksheet_txt):
                #print "Creating new worksheet file %s"%worksheet_txt
                self.__cells = []
            else:
                #print "Loading worksheet %s"%worksheet_txt
                txt = open(worksheet_txt).read()
                self.edit_save(txt)
            return self.__cells

    def append_new_cell(self):
        """
        Create and append a new cell to the list of cells.

        OUTPUT:
            a new empty cell

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: W = nb.create_new_worksheet('Test Edit Save', 'admin')
            sage: W
            [Cell 0; in=, out=]
            sage: W.append_new_cell()
            Cell 1; in=, out=
            sage: W
            [Cell 0; in=, out=, Cell 1; in=, out=]
        """
        C = self._new_cell()
        self.cell_list().append(C)
        return C

    def new_cell_before(self, id, input=""):
        """
        Insert a new cell into the cell list before the cell
        with the given integer id.  If the id is not the
        id of any cell, inserts a new cell at the end of the
        cell list.

        INPUT:
            id -- integer
            input -- string

        OUTPUT:
            new cell with the given input text (empty by default).

        """
        cells = self.cell_list()
        for i in range(len(cells)):
            if cells[i].id() == id:
                C = self._new_cell(input=input)
                cells.insert(i, C)
                return C
        C = self._new_cell(input=input)
        cells.append(C)
        return C

    def new_cell_after(self, id, input=""):
        """
        Insert a new cell into the cell list after the cell
        with the given integer id.

        INPUT:
            id -- integer
            input -- string

        OUTPUT:
            new cell with the given input text (empty by default).

        """
        cells = self.cell_list()
        for i in range(len(cells)):
            if cells[i].id() == id:
                C = self._new_cell(input=input)
                cells.insert(i+1, C)
                return C
        C = self._new_cell(input=input)
        cells.append(C)
        return C

    def delete_cell_with_id(self, id):
        """
        Remove the cell with given id and return the cell before it.
        """
        cells = self.cell_list()
        for i in range(len(cells)):
            if cells[i].id() == id:

                # Delete this cell from the queued up calculation list:
                C = cells[i]
                if C in self.__queue and self.__queue[0] != C:
                    self.__queue.remove(C)

                # Delete this cell from the list of cells in this worksheet:
                del cells[i]
                if i > 0:
                    return cells[i-1].id()
                else:
                    break
        return cells[0].id()

    ##########################################################
    # Managing whether computing is happening: stop, start, clear, etc.
    ##########################################################
    def clear(self):
        self.__comp_is_running = False
        self.__queue = []
        self.__cells = [ ]
        for i in range(INITIAL_NUM_CELLS):
            self.append_new_cell()

    def computing(self):
        """
        Return whether or not a cell is currently being run in the
        worksheet Sage process.
        """
        try:
            return self.__comp_is_running
        except AttributeError:
            return False

    def set_not_computing(self):
        self.__comp_is_running = False
        self.__queue = []

    def quit(self):
        try:
            S = self.__sage
        except AttributeError:
            # no sage running anyways!
            return

        try:
            pid = S._expect.pid
            #print "PID = ", pid
            os.killpg(pid, 9)
            os.kill(pid, 9)
            S._expect = None
        except AttributeError, msg:
            print "WARNING: %s"%msg
        except Exception, msg:
            print msg
            print "WARNING: Error deleting Sage object!"

        try:
            os.kill(pid, 9)
        except:
            pass

        del self.__sage

        # We do this to avoid getting a stale Sage that uses old code.
        self.clear_queue()
        # We do this to avoid saving this worksheet's cells to disk repeatedly.
        self.save()
        del self.__cells

    def next_block_id(self):
        try:
            i = self.__next_block_id
        except AttributeError:
            i = 0
        i += 1
        self.__next_block_id = i
        return i

    def compute_process_has_been_started(self):
        """
        Return True precisely if the compute process has been started,
        irregardless of whether or not it is currently churning away
        on a computation.
        """
        try:
            S = self.__sage
            if S._expect is None:
                return False
        except AttributeError:
            return False
        return True

    def initialize_sage(self):
        self.delete_cell_input_files()
        object_directory = os.path.abspath(self.notebook().object_directory())
        S = self.sage()
        try:
            cmd = '__DIR__="%s/"; DIR=__DIR__; DATA="%s/"; '%(self.DIR(), os.path.abspath(self.data_directory()))
            cmd += '_support_.init(None, globals()); '
            S._send(cmd)   # non blocking
        except Exception, msg:
            print "ERROR initializing compute process:\n"
            print msg
            del self.__sage
            raise RuntimeError
        A = self.attached_files()
        for F in A.iterkeys():
            A[F] = 0  # expire all
        self._enqueue_auto_cells()
        return S

    def sage(self):
        """
        Return a started up copy of Sage initialized for computations.

        If this is a published worksheet, just return None, since published
        worksheets must not have any compute functionality.

        OUTPUT:
            a Sage interface
        """
        if self.is_published():
            return None
        try:
            S = self.__sage
            if S._expect is not None:
                return S
        except AttributeError:
            pass
        self.__sage = one_prestarted_sage(server = self.notebook().get_server(),
                                          ulimit = self.notebook().get_ulimit())
        self.__next_block_id = 0
        self.initialize_sage()

        # Check to see if the typeset/pretty print button is checked.
        # If so, send code to initialize the worksheet to have the
        # right pretty printing mode.
        if self.pretty_print():
            self.__sage._send('pretty_print_default(True);')

        return self.__sage

    def eval_asap_no_output(self, cmd, username=None):
        C = self._new_cell(hidden=True)
        C.set_asap(True)
        C.set_no_output(True)
        C.set_input_text(cmd)
        self.enqueue(C, username=username)


    def start_next_comp(self):
        if len(self.__queue) == 0:
            return

        if self.__comp_is_running:
            #self._record_that_we_are_computing()
            return

        C = self.__queue[0]

        if C.interrupted():
            # don't actually compute
            return

        D = C.directory()
        if not C.introspect():
            if C.is_interacting():
                I = C.interact
            else:
                I = C.input_text().strip()
            if I in ['restart', 'quit', 'exit']:
                self.restart_sage()
                S = self.system()
                if S is None: S = 'sage'
                C.set_output_text('Exited %s process'%S,'')
                return
            if not I.startswith('%timeit'):
                if I.startswith('%time'):
                    C.do_time()
                    I = after_first_word(I).lstrip()
                elif first_word(I) == 'time':
                    C.do_time()
                    I = after_first_word(I).lstrip()
        else:
            before_prompt, after_prompt = C.introspect()
            I = before_prompt

        S = self.sage()

        id = self.next_block_id()
        C.code_id = id

        # prevent directory disappear problems
        dir = self.directory()
        code_dir = '%s/code'%dir
        if not os.path.exists(code_dir):
            os.makedirs(code_dir)
        cell_dir = '%s/cells'%dir
        if not os.path.exists(cell_dir):
            os.makedirs(cell_dir)
        tmp = '%s/code/%s.py'%(dir, id)

        absD = os.path.abspath(D)
        input = 'os.chdir("%s")\n'%absD

        os.system('chmod -R a+rw "%s"'%absD)

        # This is useful mainly for interact -- it allows
        # a cell to know it's ID.
        input += 'sage.server.notebook.interact.SAGE_CELL_ID=%s\n'%(C.id())

        if C.time():
            input += '__SAGE_t__=cputime()\n__SAGE_w__=walltime()\n'

        # If the input ends in a question mark and is *not* a comment line,
        # then we introspect on it.
        Istrip = I.strip()
        if Istrip.endswith('?') and not Istrip.startswith('#'):
            C.set_introspect(I, '')
        I = I.replace('\\\n','')
        C._before_preparse = input + I
        input += self.preparse_input(I, C)

        try:
            compile(input, '', 'exec')
        except SyntaxError, msg:
            t = traceback.format_exc()
            s = 'File "<string>",'
            i = t.find(s)
            if i != -1:
                t = t[i+len(s):]
            i = t.find('\n')
            try:
                n = int(t[t[:i].rfind(' '):i])  # line number of the exception
                try:
                    t = 'Syntax Error:\n    %s'%C._before_preparse.split('\n')[n-1]
                except IndexError:
                    pass
                if False:
                    if i != -1:
                        t = t[i:]
                    v = [w for w in t.split('\n') if w]
                    t = '\n'.join(['Syntax Error:'] + v[0:-1])
                C.set_output_text(t, '')
                del self.__queue[0]
                return
            except ValueError:
                pass

        if C.time() and not C.introspect():
            input += 'print "CPU time: %.2f s,  Wall time: %.2f s"%(cputime(__SAGE_t__), walltime(__SAGE_w__))\n'

        input = self.synchronize(input)
        # Unfortunately, this has to go here at the beginning of the file until Python 2.6,
        # in order to support use of the with statement in the notebook.  Very annoying.
        input = 'from __future__ import with_statement\n' + input

        # This magic comment at the very start of the file allows utf8
        # characters in the file
        input = '# -*- coding: utf_8 -*-\n' + input

        open(tmp,'w').write(input)

        cmd = 'execfile("%s")\n'%os.path.abspath(tmp)
        # Signal an end (which would only be seen if there is an error.)
        cmd += 'print "\\x01r\\x01e%s"'%self.synchro()
        self.__comp_is_running = True
        try:
            S._send(cmd)
        except OSError, msg:
            self.restart_sage()
            C.set_output_text('The Sage compute process quit (possibly Sage crashed?).\nPlease retry your calculation.','')

    def check_comp(self, wait=0.2):
        """
        Check on currently computing cells in the queue.

        INPUT:
            wait -- float (default: 0.2); how long to wait for output.

        EXAMPLES:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\n{{{\n3^20\n}}}')
            sage: W.cell_list()[0].evaluate()
            sage: W.check_comp()     # random output -- depends on computer speed
            ('d', Cell 0; in=3^20, out=
            3486784401
            )
            sage: nb.delete()
        """
        if len(self.__queue) == 0:
            return 'e', None
        S = self.sage()
        C = self.__queue[0]
        if C.interrupted():
            self.__comp_is_running = False
            del self.__queue[0]
            return 'd', C

        try:
            done, out, new = S._so_far(wait=wait, alternate_prompt=SAGE_END+str(self.synchro()))
        except RuntimeError, msg:
            verbose("Computation was interrupted or failed. Restarting.\n%s"%msg)
            self.__comp_is_running = False
            self.start_next_comp()
            return 'w', C

        out = self.postprocess_output(out, C)
        if not done:
            # Still computing
            out = self._process_output(out)
            if not C.introspect():
                C.set_output_text(out, '')
            #self._record_that_we_are_computing()
            return 'w', C

        # Finished a computation.
        self.__comp_is_running = False
        del self.__queue[0]

        if C.is_no_output():
            # Clean up the temp directories associated to C, and do not set any output
            # text that C might have got.
            dir = self.directory()
            code_file = '%s/code/%s.py'%(dir, C.code_id)
            # NOTE -- this deletes the input file, which in the rare case when
            # the input defines a function and the user asks for the source of
            # that function, they wouldn't get it.
            os.unlink(code_file)
            cell_dir = '%s/cells/%s'%(dir, C.id())
            shutil.rmtree(cell_dir)
            return 'd', C

        out = self._process_output(out)
        if C.introspect():
            before_prompt, after_prompt = C.introspect()
            if len(before_prompt) == 0:
                return
            if before_prompt[-1] != '?':
                # completions
                c = self.best_completion(out, C._word_being_completed)
                C.set_changed_input_text(before_prompt + c + after_prompt)
                out = self.completions_html(C.id(), out)
                C.set_introspect_html(out, completing=True)
            else:
                C.set_introspect_html(out, completing=False)
        else:
            C.set_output_text(out, C.files_html(out), sage=self.sage())
            C.set_introspect_html('')

        return 'd', C

    def interrupt(self):
        """
        Interrupt all currently queued up calculations.

        OUTPUT:
            bool -- return True if no problems interrupting calculation
                    return False if the Sage interpreter had to be restarted.

        EXAMPLES:
        We create a worksheet and start a large factorization going:
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\n{{{\nfactor(2^997-1)\n}}}')
            sage: W.cell_list()[0].evaluate()

        It's running still
            sage: W.check_comp()
            ('w', Cell 0; in=factor(2^997-1), out=...)

        We interrupt it successfully.
            sage: W.interrupt()         # random -- could fail on heavily loaded machine
            True

        Now we check and nothing is computing.
            sage: W.check_comp()        # random -- could fail on heavily loaded machine
            ('e', None)

        Clean up.
            sage: nb.delete()
        """
        if len(self.__queue) == 0:
            # nothing to do
            return True

        success = False
        # stop the current computation in the running Sage
        try:
            S = self.__sage
        except AttributeError:
            pass
        else:
            success = S.interrupt(INTERRUPT_TRIES, quit_on_fail=False)

        if success:
            self.clear_queue()

        return success

    def clear_queue(self):
        # empty the queue
        for C in self.__queue:
            C.interrupt()
        self.__queue = []
        self.__comp_is_running = False

    def restart_sage(self):
        """
        Restart \sage kernel.
        """
        self.quit()

        self.__sage = initialized_sage(server = self.notebook().get_server(),
                                       ulimit = self.notebook().get_ulimit())
        self.initialize_sage()
        self.start_next_comp()


    def worksheet_command(self, cmd):
        return '/home/%s/%s'%(self.filename(), cmd)

    ##########################################################
    # Idle timeout
    ##########################################################
    def quit_if_idle(self, timeout):
        r"""
        Quit the worksheet process if it has been ``idle'' for more than \var{timeout} seconds,
        where idle is by definition that the worksheet has not reported back that it
        is actually computing.  I.e., an ignored worksheet process (since the user closed
        their browser) is also considered idle, even if code is running.
        """
        if self.time_idle() > timeout:
            print "Quitting ignored worksheet process for '%s'."%self.name()
            self.quit()

    def time_idle(self):
        return walltime() - self.last_compute_walltime()

    def last_compute_walltime(self):
        try:
            return self.__last_compute_walltime
        except AttributeError:
            t = walltime()
            self.__last_compute_walltime = t
            return t

    def _record_that_we_are_computing(self, username=None):
        self.__last_compute_walltime = walltime()
        if username:
            self.record_edit(username)

    def ping(self, username):
        if self.is_published():
            return
        self._record_that_we_are_computing(username)

    ##########################################################
    # Enqueuing cells
    ##########################################################
    def queue(self):
        return list(self.__queue)

    def queue_id_list(self):
        return [c.id() for c in self.__queue]


    def enqueue(self, C, username=None, next=False):
        r"""
        Queue up the cell C for evaluation in this worksheet.

        INPUT:
            C -- a Cell
            username -- the name of the user that is evaluating this
                        cell (mainly used for loging)

        NOTE: If \code{C.is_asap()} is True, then we put C as close to
        the beginning of the queue as possible, but after all asap cells.
        Otherwise, C goes at the end of the queue.
        """
        if self.is_published():
            return
        self._record_that_we_are_computing(username)
        if not isinstance(C, Cell):
            raise TypeError
        if C.worksheet() != self:
            raise ValueError, "C must be have self as worksheet."

        # Now enqueue the requested cell.
        if not (C in self.__queue):
            if C.is_asap():
                if self.computing():
                    i = 1
                else:
                    i = 0
                while i < len(self.__queue) and self.__queue[i].is_asap():
                    i += 1
                self.__queue.insert(i, C)
            else:
                self.__queue.append(C)
        self.start_next_comp()

    def _enqueue_auto_cells(self):
        for c in self.cell_list():
            if c.is_auto_cell():
                self.enqueue(c)

    def set_cell_counter(self):
        self.__next_id = 1 + max([C.id() for C in self.cell_list()])

    def _new_text_cell(self, plain_text, id=None):
        if id is None:
            id = self.__next_id
            self.__next_id += 1
        return TextCell(id, plain_text, self)

    def next_hidden_id(self):
        try:
            i = self.__next_hidden_id
            self.__next_hidden_id -= 1
        except AttributeError:
            i = -1
            self.__next_hidden_id = -2
        return i

    def _new_cell(self, id=None, hidden=False, input=''):
        if id is None:
            if hidden:
                id = self.next_hidden_id()
            else:
                id = self.__next_id
                self.__next_id += 1
        return Cell(id, input, '', self)

    def append(self, L):
        self.cell_list().append(L)

    ##########################################################
    # Accessing existing cells
    ##########################################################
    def __getitem__(self, n):
        try:
            return self.cell_list()[n]
        except IndexError:
            if n >= 0:  # this should never happen -- but for robustness we cover this case.
                for k in range(len(self.cell_list()),n+1):
                    self.cell_list().append(self._new_cell())
                return self.cell_list()[n]
            raise IndexError

    def get_cell_with_id(self, id):
        for c in self.cell_list():
            if c.id() == id:
                return c
        return self._new_cell(id)

    def synchronize(self, s):
        try:
            i = (self.__synchro + 1)%65536
        except AttributeError:
            i = 0
        self.__synchro = i
        return 'print "%s%s"\n'%(SAGE_BEGIN,i) + s + '\nprint "%s%s"\n'%(SAGE_END,i)

    def synchro(self):
        try:
            return self.__synchro
        except AttributeError:
            return 0

    def delete_cell_input_files(self):
        r"""
        Delete all the files \file{code_\%s.py} and \file{code_\%s.spyx} that are created
        when evaluating cells.  We do this when we first start the notebook
        to get rid of clutter.
        """
        D = self.directory() + '/code/'
        if os.path.exists(D):
            for X in os.listdir(D):
                os.unlink('%s/%s'%(D,X))
        else:
            os.makedirs(D)

    def check_cell(self, id):
        """
        Check the status on computation of the cell with given id.

        INPUT:
            id -- an integer

        OUTPUT:
            status -- a string, either 'd' (done) or 'w' (working)
            cell -- the cell with given id
        """
        cell = self.get_cell_with_id(id)

        if cell in self.__queue:
            status = 'w'
        else:
            status = 'd'
        return status, cell

    def is_last_id_and_previous_is_nonempty(self, id):
        if self.cell_list()[-1].id() != id:
            return False
        if len(self.cell_list()) == 1:
            return False
        if len(self.cell_list()[-2].output_text(ncols=0)) == 0:
            return False
        return True


    ##########################################################
    # (Tab) Completions
    ##########################################################
    def best_completion(self, s, word):
        completions = s.split()
        if len(completions) == 0:
            return ''
        n = len(word)
        i = n
        m = min([len(x) for x in completions])
        while i <= m:
            word = completions[0][:i]
            for w in completions[1:]:
                if w[:i] != word:
                    return w[n:i-1]
            i += 1
        return completions[0][n:m]

    def completions_html(self, id, s, cols=3):
        if 'no completions of' in s:
            return ''

        completions = s.split()

        n = len(completions)
        l = n/cols + n%cols

        if n == 1:
            return '' # don't show a window, just replace it

        rows = []
        for r in range(0,l):
            row = []
            for c in range(cols):
                try:
                    cell = completions[r + l*c]
                    row.append(cell)
                except:
                    pass
            rows.append(row)
        return format_completions_as_html(id, rows)

    ##########################################################
    # Processing of input and output to worksheet process.
    ##########################################################
    def preparse_input(self, input, C):
        C.set_is_html(False)
        introspect = C.introspect()
        if introspect:
            input = self.preparse_introspection_input(input, C, introspect)
        else:
            switched, input = self.check_for_system_switching(input, C)
            if not switched:
                input = self.preparse_nonswitched_input(input)
            input += '\n'
        return input

    def preparse_introspection_input(self, input, C, introspect):
        before_prompt, after_prompt = introspect
        i = 0
        while i < len(after_prompt):
            if after_prompt[i] == '?':
                if i < len(after_prompt)-1 and after_prompt[i+1] == '?':
                    i += 1
                before_prompt += after_prompt[:i+1]
                after_prompt = after_prompt[i+1:]
                C.set_introspect(before_prompt, after_prompt)
                break
            elif after_prompt[i] in ['"', "'", ' ', '\t', '\n']:
                break
            i += 1
        if before_prompt.endswith('??'):
            input = self._get_last_identifier(before_prompt[:-2])
            input = 'print _support_.source_code("%s", globals(), system="%s")'%(input, self.system())
        elif before_prompt.endswith('?'):
            input = self._get_last_identifier(before_prompt[:-1])
            input = 'print _support_.docstring("%s", globals(), system="%s")'%(input, self.system())
        else:
            input = self._get_last_identifier(before_prompt)
            C._word_being_completed = input
            input = 'print "\\n".join(_support_.completions("%s", globals(), system="%s"))'%(input, self.system())
        return input

    def preparse_nonswitched_input(self, input):
        input = ignore_prompts_and_output(input).rstrip()
        input = self.preparse(input)
        input = self.load_any_changed_attached_files(input)
        input = self.do_sage_extensions_preparsing(input)
        input = input.split('\n')

        # The following is all so the last line (or single lines)
        # will implicitly print as they should, unless they are
        # an assignment.   "display hook"  It's very complicated,
        # but it has to be...
        i = len(input)-1
        if i >= 0:
            while len(input[i]) > 0 and input[i][0] in ' \t':
                i -= 1
            t = '\n'.join(input[i:])
            if not t.startswith('def '):
                try:
                    compile(t+'\n', '', 'single')
                    t = t.replace("'", "\\u0027").replace('\n','\\u000a')
                    # IMPORTANT: If you change this line, also change
                    # the function format_exception in cell.py
                    input[i] = "exec compile(ur'%s' + '\\n', '', 'single')"%t
                    input = input[:i+1]
                except SyntaxError, msg:
                    pass
        input = '\n'.join(input)
        return input


    def _strip_synchro_from_start_of_output(self, s):
        z = SAGE_BEGIN+str(self.synchro())
        i = s.find(z)
        if i == -1:
            # did not find any synchronization info in the output stream
            j = s.find('Traceback')
            if j != -1:
                # Probably there was an error; better not hide it.
                return s[j:]
            else:
                # Maybe we just read too early -- supress displaying anything yet.
                return ''
        else:
            return s[i+len(z):]

    def _process_output(self, s):
        s = re.sub('\x08.','',s)
        s = self._strip_synchro_from_start_of_output(s)
        if SAGE_ERROR in s:
            i = s.rfind('>>>')
            if i >= 0:
                return s[:i-1]
        # Remove any control codes that might have not got stripped out.
        return s.replace(SAGE_BEGIN,'').replace(SAGE_END,'').replace(SC,'')

    def postprocess_output(self, out, C):
        i = out.find('\r\n')
        out = out[i+2:]
        #out = out.rstrip()
        if C.introspect():
            return out

        # this isn't needed anymore !
        # the python error message for list indices is not good enough.
        # out = out.replace('indices must be integers', 'indices must be of type Python int.\n(Hint: Use int(n) to make n into a Python int.)')

        out = out.replace("NameError: name 'os' is not defined", "NameError: name 'os' is not defined\nTHERE WAS AN ERROR LOADING THE SAGE LIBRARIES.  Try starting Sage from the command line to see what the error is.")

        try:
            tb = 'Traceback (most recent call last):'
            i = out.find(tb)
            if i != -1:
                t = '.py", line'
                j = out.find(t)
                z = out[j+5:].find(',')
                n = int(out[j+len(t):z + j+5])
                k = out[j:].find('\n')
                if k != -1:
                    k += j
                    l = out[k+1:].find('\n')
                    if l != -1:
                        l += k+1
                        I = C._before_preparse.split('\n')
                        out = out[:i + len(tb)+1] + '    ' + I[n-2] + out[l:]
        except (ValueError, IndexError), msg:
            pass
        return out

    def _get_last_identifier(self, s):
        return support.get_rightmost_identifier(s)

    def preparse(self, s):
        if sage.misc.interpreter.do_preparse:
            s = preparse_file(s, magic=False, do_time=True,
                              ignore_prompts=False)
        return s

    ##########################################################
    # Loading and attaching files
    ##########################################################
    def load_any_changed_attached_files(self, s):
        r"""
        Modify \var{s} by prepending any necessary load commands
        corresponding to attached files that have changed.
        """
        A = self.attached_files()
        init_sage = DOT_SAGE + 'init.sage'
        if not init_sage in A.keys() and os.path.exists(init_sage):
            A[init_sage] = 0

        # important that this is A.items() and not A.iteritems()
        # since we change A during the iteration.
        for F, tm in A.items():
            try:
                new_tm = os.path.getmtime(F)
            except OSError:
                del A[F]
            else:
                if new_tm > tm:
                    A[F] = new_tm
                    s = 'load %s\n'%F + s
        return s

    def attached_files(self):
        try:
            A = self.__attached
        except AttributeError:
            A = {}
            self.__attached = A

        return A

    def attach(self, filename):
        A = self.attached_files()
        try:
            A[filename] = os.path.getmtime(filename)
        except OSError:
            print "WARNING: File %s vanished"%filename
            pass

    def detach(self, filename):
        A = self.attached_files()
        try:
            A.pop(filename)
        except KeyError:
            pass

    def _normalized_filenames(self, L):
        i = L.find('#')
        if i != -1:
            L = L[:i]
        a = []
        OBJECTS = os.path.abspath(self.notebook().object_directory())
        for filename in L.split():
            filename = filename.strip('"').strip("'")
            if os.path.exists(OBJECTS + '/' + filename):
                filename = OBJECTS + '/' + filename
            elif os.path.exists(OBJECTS + '/' + filename + '.sobj'):
                filename = OBJECTS + '/' + filename + '.sobj'
            else:
                if len(filename) > 0 and filename[0] != '/':
                    filename = '%s/%s'%(self.DIR(), filename)
                if not filename.endswith('.py') and not filename.endswith('.sage') and \
                       not filename.endswith('.sobj') and not os.path.exists(filename):
                    if os.path.exists(filename + '.sage'):
                        filename = filename + '.sage'
                    elif os.path.exists(filename + '.py'):
                        filename = filename + '.py'
                    elif os.path.exists(filename + '.sobj'):
                        filename = filename + '.sobj'
            a.append(filename)
        return a

    def load_path(self):
        D = self.cells_directory()
        return [self.directory() + '/data/'] + [D + x for x in os.listdir(D)]

    def hunt_file(self, filename):
        if filename.lower().startswith('http://'):
            filename = remote_file.get_remote_file(filename)
        if not os.path.exists(filename):
            fn = os.path.split(filename)[-1]
            for D in self.load_path():
                t = D + '/' + fn
                if os.path.exists(t):
                    filename = t
                    break
                if os.path.exists(t + '.sobj'):
                    filename = t + '.sobj'
                    break
        return os.path.abspath(filename)

    def _load_file(self, filename, files_seen_so_far, this_file):
        if filename.endswith('.sobj'):
            name = os.path.splitext(filename)[0]
            name = os.path.split(name)[-1]
            return '%s = load("%s");'%(name, filename)

        if filename in files_seen_so_far:
            t = "print 'WARNING: Not loading %s -- would create recursive load'"%filename

        try:
            F = open(filename).read()
        except IOError:
            return "print 'Error loading %s -- file not found'"%filename
        else:
            filename_orig = filename
            filename = filename.rstrip('.txt')
            if filename.endswith('.py'):
                t = F
            elif filename.endswith('.spyx') or filename.endswith('.pyx'):
                cur = os.path.abspath(os.curdir)
                try:
                    mod, dir  = cython.cython(filename_orig, compile_message=True, use_cache=True)
                except (IOError, OSError, RuntimeError), msg:
                    return "print r'''Error compiling cython file:\n%s'''"%msg
                t  = "import sys\n"
                t += "sys.path.append('%s')\n"%dir
                t += "from %s import *\n"%mod
                return t
            elif filename.endswith('.sage'):
                t = self.preparse(F)
            else:
                t = "print 'Loading of file \"%s\" has type not implemented.'"%filename

        t = self.do_sage_extensions_preparsing(t,
                          files_seen_so_far + [this_file], filename)
        return t

    def _save_objects(self, s):
        s = s.replace(',',' ').replace('(',' ').replace(')',' ')
        v = s.split()
        return ';'.join(['save(%s,"%s")'%(x,x) for x in v])


    def do_sage_extensions_preparsing(self, s, files_seen_so_far=[], this_file=''):
        u = []
        for t in s.split('\n'):
            if t.startswith('load '):
                z = ''
                for filename in self._normalized_filenames(after_first_word(t)):
                    filename = self.hunt_file(filename)
                    z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t.startswith('attach '):
                z = ''
                for filename in self._normalized_filenames(after_first_word(t)):
                    filename = self.hunt_file(filename)
                    if not os.path.exists(filename):
                        z += "print 'Error attaching %s -- file not found'\n"%filename
                    else:
                        self.attach(filename)
                        z += self._load_file(filename, files_seen_so_far, this_file) + '\n'
                t = z

            elif t.startswith('detach '):
                filename = self.hunt_file(filename)
                for filename in self._normalized_filenames(after_first_word(t)):
                    self.detach(filename)
                t = ''

            elif t.startswith('save '):
                t = self._save_objects(after_first_word(t))

            u.append(t)

        return '\n'.join(u)

    def _eval_cmd(self, system, cmd, dir):
        cmd = cmd.replace("'", "\\u0027")
        return "print _support_.syseval(%s, ur'''%s''', '%s')"%(system, cmd, dir)

    ##########################################################
    # Parsing the %cython, %jsmath, %python, etc., extension.
    ##########################################################

    def cython_import(self, cmd, C):
        # Choice: Can use either C.relative_id() or self.next_block_id().
        # C.relative_id() has the advantage that block evals are cached, i.e.,
        # no need to recompile.  On the other hand tracebacks don't work if
        # you change a cell and create a new function in it.  Caching is
        # also annoying since the linked .c file disappears.
        id = self.next_block_id()
        # id = C.relative_id()
        spyx = os.path.abspath('%s/code/sage%s.spyx'%(self.directory(), id))
        if not (os.path.exists(spyx) and open(spyx).read() == cmd):
            open(spyx,'w').write(cmd)
        s  = '_support_.cython_import_all("%s", globals())'%spyx
        return s

    def check_for_system_switching(self, s, C):
        r"""
        Check for input cells that start with \code{\%foo},
        where \var{foo} is an object with an eval method.

        INPUT:
            s -- a string of the code from the cell to be executed
            C -- the cell object

        EXAMPLES:
        First, we set up a new notebook and worksheet.

            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')

        We first test running a native command in 'sage' mode and then a GAP cell
        within Sage mode.

            sage: W.edit_save('Sage\nsystem:sage\n{{{\n2+3\n}}}\n\n{{{\n%gap\nSymmetricGroup(5)\n}}}')
            sage: c0, c1 = W.cell_list()
            sage: W.check_for_system_switching(c0.input_text(), c0)
            (False, '2+3')
            sage: W.check_for_system_switching(c1.input_text(), c1)
            (True,
             "print _support_.syseval(gap, ur'''SymmetricGroup(5)''', '...')")

            sage: c0.evaluate()
            sage: W.check_comp()  #random output -- depends on the computer's speed
            ('d', Cell 0; in=2+3, out=
            5
            )
            sage: c1.evaluate()
            sage: W.check_comp()  #random output -- depends on the computer's speed
            ('d', Cell 1; in=%gap
            SymmetricGroup(5), out=
            Sym( [ 1 .. 5 ] )
            )

        Next, we run the same commands but from 'gap' mode.

            sage: W.edit_save('Sage\nsystem:gap\n{{{\n%sage\n2+3\n}}}\n\n{{{\nSymmetricGroup(5)\n}}}')
            sage: c0, c1 = W.cell_list()
            sage: W.check_for_system_switching(c0.input_text(), c0)
            (False, '2+3')
            sage: W.check_for_system_switching(c1.input_text(), c1)
            (True,
             "print _support_.syseval(gap, ur'''SymmetricGroup(5)''', '...')")
            sage: c0.evaluate()
            sage: W.check_comp()  #random output -- depends on the computer's speed
            ('d', Cell 0; in=%sage
            2+3, out=
            5
            )
            sage: c1.evaluate()
            sage: W.check_comp()  #random output -- depends on the computer's speed
            ('d', Cell 1; in=SymmetricGroup(5), out=
            Sym( [ 1 .. 5 ] )
            )

        """
        s = s.lstrip()
        if self.system() != 'sage':
            if len(s) == 0 or s[0] != "%":
                #Since no other system is specified, we return True
                #and evaluate the code using self.system()
                return True, self._eval_cmd(self.system(), s, os.path.abspath(C.directory()))
            elif s.startswith("%sage"):
                #Since the code we want to run is Sage code, we return
                #False.
                s = after_first_word(s).lstrip()
                return False, s
        else:
            #Since the system in Sage and there is not another
            #system specified, we return False.
            if len(s) == 0 or s[0] != '%':
                return False, s

        if s.startswith('%hide'):
            t = after_first_word(s).lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s.startswith('%save_server'):
            self.notebook().save()
            t = after_first_word(s).lstrip()
            if len(t) == 0 or t[0] != '%':
                return False, t
            s = t
        if s.startswith("%cython") or s.startswith("%pyrex") or s.startswith("%sagex"):  # a block of Cython code.
            return True, self.cython_import(after_first_word(s).lstrip(), C)

        # Determine system = the system doing the computation, e.g., %magma, and
        #           code_to_eval = the code to feed to the system via .eval.
        system = first_word(s)[1:]     # get rid of the percent sign.
        code_to_eval = after_first_word(s)
        cmd = self._eval_cmd(system, code_to_eval, os.path.abspath(C.directory()))
        if system == 'html':
            C.set_is_html(True)
        return True, cmd

    ##########################################################
    # List of attached files.
    ##########################################################
    def attached_html(self):
        s = ''
        div = '<div class="attached_filename" onClick="inspect_attached_file(\'%s\')">'
        A = self.attached_files()
        D = self.DIR()
        for F, tm in A.iteritems():
            # uncomment this to remove some absolute path info...
            # if F[:len(D)] == D: F = F[len(D)+1:]
            s += div%F + '%s</div>'%F
        return s

    ##########################################################
    # Showing and hiding all cells
    ##########################################################
    def show_all(self):
        for C in self.cell_list():
            try:
                C.set_cell_output_type('wrap')
            except AttributeError:   # for backwards compatibility
                pass

    def hide_all(self):
        for C in self.cell_list():
            try:
                C.set_cell_output_type('hidden')
            except AttributeError:
                pass

    def delete_all_output(self, username):
        """
        Delete all the output in all the worksheet cells.

        INPUT:
            username -- name of the user requesting the deletion.

        EXAMPLES:
        We create a new notebook, user, and a worksheet with one cell.
            sage: nb = sage.server.notebook.notebook.Notebook(tmp_dir())
            sage: nb.add_user('sage','sage','sage@sagemath.org',force=True)
            sage: W = nb.create_new_worksheet('Test', 'sage')
            sage: W.edit_save('Sage\nsystem:sage\n{{{\n2+3\n///\n5\n}}}')

        Notice that there is 1 cell with 5 in its output.
            sage: W.cell_list()
            [Cell 0; in=2+3, out=5]

        We now delete the output, observe that it is gone.
            sage: W.delete_all_output('sage')
            sage: W.cell_list()
            [Cell 0; in=2+3, out=]

        If an invalid user tries to delete all, a ValueError is raised.
            sage: W.delete_all_output('hacker')
            Traceback (most recent call last):
            ...
            ValueError: user 'hacker' not allowed to edit this worksheet

        Clean up.
            sage: nb.delete()
        """
        if not self.user_can_edit(username):
            raise ValueError, "user '%s' not allowed to edit this worksheet"%username
        self.save_snapshot(username)
        for C in self.cell_list():
            C.delete_output()


__internal_test1 = '''
def foo(x):
    "
    EXAMPLES:
        sage: 2+2
        4
    "
    return x
'''.lstrip()

__internal_test2 = '''
sage: 2 + 2
4
'''.lstrip()

def ignore_prompts_and_output(aString):
    r"""
    Given a string s that defines an input block of code,
    if the first line begins in \samp{sage:} (or \samp{>>>}),
    strip out all lines
    that don't begin in either \samp{sage:} (or \samp{>>>}) or \samp{...}, and
    remove all \samp{sage:} (or \samp{>>>}) and \samp{...} from the beginning
    of the remaining lines.

    TESTS:
        sage: test1 = sage.server.notebook.worksheet.__internal_test1
        sage: test1 == sage.server.notebook.worksheet.ignore_prompts_and_output(test1)
        True

        sage: test2 = sage.server.notebook.worksheet.__internal_test2
        sage: sage.server.notebook.worksheet.ignore_prompts_and_output(test2)
        '2 + 2\n'
    """
    s = aString.lstrip()
    is_example = s.startswith('sage:') or s.startswith('>>>')
    if not is_example:
        return aString # return original, not stripped copy
    new = ''
    lines = s.split('\n')
    for line in lines:
        line = line.lstrip()
        if line.startswith('sage:'):
            new += after_first_word(line).lstrip() + '\n'
        elif line.startswith('>>>'):
            new += after_first_word(line).lstrip() + '\n'
        elif line.startswith('...'):
            new += after_first_word(line) + '\n'
    return new

def extract_text_before_first_compute_cell(text):
    """
    OUTPUT:
        Everything in text up to the first \{\{\{.
    """
    i = text.find('{{{')
    if i == -1:
        return text
    return text[:i]

def extract_first_compute_cell(text):
    """
    INPUT:
        a block of wiki-like marked up text
    OUTPUT:
        meta -- meta information about the cell (as a dictionary)
        input -- string, the input text
        output -- string, the output text
        end -- integer, first position after \samp{\}\}\}} in text.
    """
    # Find the input block
    i = text.find('{{{')
    if i == -1:
        raise EOFError
    j = text[i:].find('\n')
    if j == -1:
        raise EOFError
    k = text[i:].find('|')
    if k != -1 and k < j:
        try:
            meta = dictify(text[i+3:i+k])
        except TypeError:
            meta = {}
        i += k + 1
    else:
        meta = {}
        i += 3

    j = text[i:].find('\n}}}')
    if j == -1:
        j = len(text)
    else:
        j += i
    k = text[i:].find('\n///')
    if k == -1 or k+i > j:
        input = text[i:j]
        output = ''
    else:
        input = text[i:i+k].strip()
        output = text[i+k+4:j].strip()

    return meta, input.strip(), output, j+4

def after_first_word(s):
    """
    Return everything after the first whitespace in the string s.
    Returns the empty string if there is nothing after the
    first whitespace.

    INPUT:
        s -- string
    OUTPUT:
        a string

    EXAMPLES:
        sage: from sage.server.notebook.worksheet import after_first_word
        sage: after_first_word("\%gap\n2+2\n")
        '2+2\n'
        sage: after_first_word("2+2")
        ''
    """
    i = whitespace.search(s)
    if i is None:
        return ''
    return s[i.start()+1:]

def first_word(s):
    """
    Returns everything before the first whitespace in the string s.
    If there is no whitespace, then the entire string s is returned.

    EXAMPLES:
        sage: from sage.server.notebook.worksheet import first_word
        sage: first_word("\%gap\n2+2\n")
        '\\%gap'
        sage: first_word("2+2")
        '2+2'
    """
    i = whitespace.search(s)
    if i is None:
        return s
    return s[:i.start()]



def format_completions_as_html(cell_id, completions):
    if len(completions) == 0:
        return ''
    lists = []

    # compute the width of each column
    column_width = []
    for i in range(len(completions[0])):
        column_width.append(max([len(x[i]) for x in completions if i < len(x)]))

    for i in range(len(completions)):
        row = completions[i]
        for j in range(len(row)):
            if len(lists) <= j:
                lists.append([])
            cell = """
<li id='completion%s_%s_%s' class='completion_menu_two'>
<a onClick='do_replacement(%s, "%s"); return false;'
   onMouseOver='this.focus(); select_replacement(%s,%s);'
>%s</a>
</li>"""%(cell_id, i, j, cell_id, row[j], i,j,
         row[j])
         #row[j] + '&nbsp;'*(column_width[j]-len(row[j])) )

            lists[j].append(cell)

    grid = "<ul class='completion_menu_one'>"
    for L in lists:
        s = "\n   ".join(L)
        grid += "\n <li class='completion_menu_one'>\n  <ul class='completion_menu_two'>\n%s\n  </ul>\n </li>"%s

    return grid + "\n</ul>"


def extract_name(text):
    # The first line is the title
    i = non_whitespace.search(text)
    if i is None:
        name = 'Untitled'
        n = 0
    else:
        i = i.start()
        j = text[i:].find('\n')
        if j != -1:
            name = text[i:i+j]
            n = j+1
        else:
            name = text[i:]
            n = len(text)-1
    return name.strip(), n

def extract_system(text):
    # If the first line is "system: ..." , then it is the system.  Otherwise the system is Sage.
    i = non_whitespace.search(text)
    if i is None:
        return 'sage', 0
    else:
        i = i.start()
        if not text[i:].startswith('system:'):
            return 'sage', 0
        j = text[i:].find('\n')
        if j != -1:
            system = text[i:i+j][7:].strip()
            n = j+1
        else:
            system = text[i:][7:].strip()
            n = len(text)-1
        return system, n


def dictify(s):
    """
    INPUT:
        s -- a string like 'in=5, out=7'
    OUTPUT:
        dict -- such as {'in':5, 'out':7}
    """
    w = []
    try:
        for v in s.split(','):
            a, b = v.strip().split('=')
            try:
                b = eval(b)
            except:
                pass
            w.append([a, b])
    except ValueError:
        return {}
    return dict(w)


def next_available_id(v):
    """
    Return smallest nonnegative integer not in v.
    """
    i = 0
    while i in v:
        i += 1
    return i


def convert_seconds_to_meaningful_time_span(t):
    if t < 60:
        s = int(t)
        if s == 1:
            return "1 second"
        return "%d seconds"%s
    if t < 3600:
        m = int(t/60)
        if m == 1:
            return "1 minute"
        return "%d minutes"%m
    if t < 3600*24:
        h = int(t/3600)
        if h == 1:
            return "1 hour"
        return "%d hours"%h
    d = int(t/(3600*24))
    if d == 1:
        return "1 day"
    return "%d days"%d


def convert_time_to_string(t):
    return time.strftime('%B %d, %Y %I:%M %p', time.localtime(float(t)))


def split_search_string_into_keywords(s):
    r"""
    The point of this function is to allow for searches like this:

    \begin{verbatim}
          "ws 7" foo bar  Modular  '"the" end'
    \end{verbatim}

    i.e., where search terms can be in quotes and the different quote
    types can be mixed.

    INPUT:
        s -- a string

    OUTPUT:
        list -- a list of strings
    """
    ans = []
    while len(s) > 0:
        word, i = _get_next(s, '"')
        if i != -1:
            ans.append(word)
            s = s[i:]
        word, j = _get_next(s, "'")
        if j != -1:
            ans.append(word)
            s = s[j:]
        if i == -1 and j == -1:
            break
    ans.extend(s.split())
    return ans


def _get_next(s, quote='"'):
    i = s.find(quote)
    if i != -1:
        j = s[i+1:].find(quote)
        if j != -1:
            return s[i+1:i+1+j].strip(), i+1+j
    return None, -1
