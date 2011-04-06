r"""
Sage Interface to the HG/Mercurial Revision Control System

These functions make setup and use of source control with Sage
easier, using the distributed Mercurial HG source control system.
To learn about Mercurial, see
http://www.selenic.com/mercurial/wiki/, in particular
http://mercurial.selenic.com/wiki/UnderstandingMercurial.
See also the `Sage Developer's Guide
<http://www.sagemath.org/doc/developer/>`_ and
http://wiki.sagemath.org/MercurialQueues for information about using
Mercurial with Sage.

Some useful commands:

-  Use ``hg_sage.diff()`` to view any changes made to the repository.

-  Use ``hg_sage.log()`` to see the change log for the repository.

-  Use ``hg_sage.serve()`` to start a web server for examining the
   repository.

-  Use ``hg_sage.commit()`` or ``hg_sage.record()`` to record any
   changes you've made to the repository.

-  Use ``hg_sage.export('tip')`` to produce a patch file for posting
   to the Sage trac server.

-  Use ``hg_sage.import_patch('file.patch')`` to import the Mercurial
   patch file ``file.patch``.

-  Use ``hg_sage.revert('file', rev=1234)`` reverts ``file`` to the
   contents it had in revision 1234.

-  Use ``hg_sage.rollback()`` to remove recorded patches without
   changing the working copy.

-  Use ``hg_sage.pull()`` to synchronize with the
   latest official stable Sage changesets.

If you want to use Mercurial queues, then type ``hg_sage.q[TAB]`` to
see the list of available methods.  Indeed, many Mercurial commands
are provided by methods here -- type ``hg_sage.[TAB]`` to get a full
list -- and you can execute any Mercurial command using
``hg_sage(COMMAND)``.  Finally, as listed, the above commands deal with
the Mercurial repository for the Sage library.  If you want to work
with other repositories distributed with Sage, this file provides the
following -- replace "hg_sage" with each of the following commands to
work with the given repository:

- ``hg_scripts`` -- the scripts repository (files in
  :file:`SAGE_ROOT/local/bin`)

- ``hg_sagenb`` -- the Sage notebook repository (files in
  :file:`SAGE_ROOT/devel/sagenb`)

- ``hg_root`` -- the Sage root repository (including files in
  :file:`SAGE_ROOT` and :file:`SAGE_ROOT/spkg`)

- ``hg_extcode`` -- the extcode repository (files in
  :file:`SAGE_ROOT/data/extcode`)
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2007 Jonathan Hanke <jonhanke@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, shutil

from   viewer import browser
from   misc   import tmp_filename, branch_current_hg, embedded
from   remote_file import get_remote_file as get_remote_file0
from   sage.server.misc import print_open_msg
from   subprocess import Popen
import re

sage_trac_re = re.compile('http[s]?://(sagetrac\.org|trac\.sagemath\.org)/sage_trac/attachment/ticket/[0-9]+/.*\.(patch|hg)')

def get_remote_file(f, **kwds):
    """
    Wrap the get_remote_file method to move the file if it ends in
    ?stuff, as happens with funny URLs from web servers.
    """
    g = get_remote_file0(f, **kwds)
    i = g.find('?')
    if i >= 0:
        h = g[:i]
        os.rename(g,h)
        return h
    return g

def pager():
    r"""
    Return a pager program, either 'cat' or 'less':
    'cat' if embedded in the notebook, 'less' otherwise.

    It is returned as a string suitable for a config option for the
    'hg' command.

    EXAMPLES::

        sage: sage.server.support.EMBEDDED_MODE=False
        sage: sage.misc.hg.pager()
        '--config pager.pager="LESS=\'R\' less"'
        sage: sage.server.support.EMBEDDED_MODE=True
        sage: sage.misc.hg.pager()
        '--config pager.pager=cat'
        sage: sage.server.support.EMBEDDED_MODE=False
    """
    if embedded():
        return '--config pager.pager=cat'
    else:
        return '--config pager.pager="LESS=\'R\' less"'

def color():
    """
    Color option for Mercurial.

    This is empty when called from the command-line, and it disables
    the "color" extension when called from the notebook.  According
    to the Mercurial docs, "color" is only used by the Mercurial
    commands diff, status, and qseries; however, it also seems to be
    used by a few other commands (like log, and qapplied, among
    others).  This function is used in :meth:`HG.diff`,
    :meth:`HG.log`, :meth:`HG.status`, :meth:`HG.qdiff`,
    :meth:`HG.qseries`, :meth:`HG.qapplied`, and
    :meth:`HG.qunapplied`.

    EXAMPLES::

        sage: sage.server.support.EMBEDDED_MODE=False
        sage: sage.misc.hg.color()
        ''
        sage: sage.server.support.EMBEDDED_MODE=True
        sage: sage.misc.hg.color()
        '--config color.mode=off'
        sage: sage.server.support.EMBEDDED_MODE=False
    """
    if embedded():
        return '--config color.mode=off'
    else:
        return ''

hg_docstring = r"""
This is an HG (Mercurial) repository.

To learn about Mercurial, see http://www.selenic.com/mercurial/wiki/.

This system is fully usable from both the command line and the Sage notebook.

Most commands are directly provided as member functions.  However,
you can use the full functionality of hg, i.e.,

``hg_%(obj_name)s("command line arguments")``

is *exactly* the same as typing::

        cd %(dir)s && hg command line arguments

"""


class HG:
    def __init__(self, dir, name, pull_url, push_url, target=None, cloneable=False, obj_name=''):
        """
        INPUT:

        - ``dir`` - directory that will contain the repository

        - ``name`` - a friendly name for the repository (only used for
          printing)

        - ``pull_url`` - a default URL to pull or record sends against
          (e.g., this could be a master repository on
          modular.math.washington.edu)

        - ``push_url`` - a default URL to push or record outgoing
          changes against (e.g., this could be a local repository on
          your favorite computer)

        - ``target`` - if the last part of dir is, e.g., sage-hg,
          create a symlink from sage-hg to target. If target=None,
          this symlink will not be created.

        TESTS::

            sage: 'scripts' in hg_scripts.__doc__
            True
        """
        self.__dir = os.path.abspath(dir)
        self.__name = name
        self.__pull_url = pull_url
        self.__push_url = push_url
        self.__initialized = False
        self.__target = target
        self.__cloneable = cloneable
        self.__obj_name = obj_name

        self.__doc__ = hg_docstring%{'obj_name':obj_name,
                                     'dir':self.__dir}

    def __repr__(self):
        """
        EXAMPLES::

            sage: hg_sage
            Hg repository 'Sage Library Source Code' in directory ...
        """
        return "Hg repository '%s' in directory %s"%(self.__name, self.__dir)


    def current_branch(self, print_flag=True):
        """
        Prints the current branch in the main Sage library.

        If ``print_flag`` is True, the default, then print the message
        "The current branch is NAME".  If False, return the string NAME.

        .. warning::

            This prints the current branch for the main Sage library
            repository, even if you call a command like
            "hg_scripts.current_branch()" which refers to a different
            repository.

        EXAMPLES::

            sage: hg_sage.current_branch()
            The current branch is: ...
        """
        branch_name = branch_current_hg()
        if print_flag:
            print "The current branch is: " + branch_name
        else:
            return branch_name

    def list_branches(self, print_flag=True):
        """
        Print all branches in the current Sage installation.

        If ``print_flag`` is True, the default, then print the message
        "Branches found:", followed by a list of the branches.  If
        False, return the list of the names of the branches.

        .. warning::

            This lists the branches for the main Sage library
            repository, even if you call a command like
            "hg_scripts.list_branches()" which refers to a different
            repository.

        EXAMPLES::

            sage: hg_sage.list_branches()
            Branches found:
            ...
            sage: 'main' in hg_sage.list_branches(print_flag=False)
            True
        """
        try:
            tmp_branch_list = [s[5:]  for s in os.listdir(SAGE_ROOT + "/devel")  if s.startswith("sage-")]
        except:
            raise RuntimeError, "Oops!  We had trouble...  Check that SAGE_ROOT gives the correct directory."

        if print_flag:
            print "Branches found:"
            for s in tmp_branch_list:
                print "    " + s
        else:
            return tmp_branch_list


    def status(self, debug=True):
        """
        Print the output of the command "hg status" for the repository.

        If ``debug`` is True, also print the full system command being
        executed.

        EXAMPLES::

            sage: hg_sage.status()
            Getting status of modified or unknown files:
            cd ... && hg status
            ...
            sage: hg_sage.status(debug=False)
            Getting status of modified or unknown files:
            ...
        """
        print("Getting status of modified or unknown files:")
        self('status %s' % (color(),), debug=debug)
        print "\n---\n"
        if self.__name == "Sage Library Source Code":
            b = branch_current_hg()
            if b == '': b='main'
            elif b[-1] == '/':
                b = b[:-1]
            print("Branch: %s"%b)

    def _changed_files(self):
        """
        EXAMPLES::

            sage: hg_sage._changed_files() # random
            False
        """
        out, err = self('status', interactive=False, debug=False)
        v = [x for x in out.split('\n') if (x.strip()[:1] != '?' and x.strip()[:1] != '!') and len(x) != 0]
        return len(v) > 0

    def _ensure_safe(self):
        """
        Ensure that the repository is in a safe state to have changes
        applied to it, i.e., that all changes to controlled files in the
        working directory are recorded.

        EXAMPLES:

            sage: hg_sage._ensure_safe()  # not tested
        """
        if self._changed_files():
            self.ci()
        if self._changed_files():
            raise RuntimeError, "Refusing to do operation since you still have unrecorded changes. You must check in all changes in your working repository first."

    def _warning(self):
        """
        Print a warning if the user has no .hgrc file.

        EXAMPLES::

            sage: hg_sage._warning() # random
        """
        from sage.plot.plot import DOCTEST_MODE
        if not os.path.exists(os.path.join(os.environ['HOME'], '.hgrc')) and not DOCTEST_MODE:
            print "\nWARNING:"
            print "Make sure to create a ~/.hgrc file:"
            print "-"*70
            print "[ui]"
            print "username = William Stein <wstein@gmail.com>"
            print "-"*70
            print "\n"

    def __call__(self, cmd=None, interactive=True, debug=True):
        """
        Run 'hg cmd' where cmd is an arbitrary string in the hg
        repository.

        INPUT:

        -  ``cmd`` - string, the hg command line (everything
           after 'hg')

        -  ``interactive`` - If True, runs using os.system, so
           user can interactively interact with hg, i.e., this is needed when
           you record changes because the editor pops up. If False, Popen is
           used to launch hg as a subprocess.

        - ``debug`` - if True, print the full system command being
          executed.

        OUTPUT:

        - If interactive is True, returns the exit code of the
          system call.

        - If interactive is False, returns the output and
          error text.

        - If cmd is not supplied, returns the output of the
          'status' command

        EXAMPLES::

            sage: hg_sage('hello') # not tested
            hg: unknown command 'hello'
            ...
            sage: hg_sage('status')  # not tested
            ...
        """
        self._warning()
        if cmd is None:
            cmd = 'status'
        s = 'cd "%s" && hg %s'%(self.__dir, cmd)
        if debug:
            print s
        if interactive:
            e = os.system(s)
            return e
        else:
            from subprocess import PIPE
            x = Popen(s, shell=True,
                      stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
            x.stdin.close()
            out = x.stdout.read()
            err = x.stderr.read()
            return out, err

    def serve(self, port=8200, address='localhost',
              open_viewer=True, options='', debug=True):
        """
        Start a web server for this repository.

        This server allows you to browse all files in the repository,
        see their changelogs, see who wrote any given line, etc.

        INPUT:

        -  ``port`` - port that the server will listen on

        -  ``address`` - (default: 'localhost') address to
           listen on

        -  ``open_viewer`` - boolean (default: True); whether
           to pop up the web page

        -  ``options`` - a string passed directly to hg's serve
           command.

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.serve()  # not tested
        """
        if open_viewer:
            cmd = 'sleep 1; %s http://%s:%s 1>&2 >/dev/null'%(browser(),
                                                              address, port)
            t = tmp_filename()
            open(t,'w').write(cmd)
            P = os.path.abspath(t)
            os.system('chmod +x %s; %s &'%(P, P))

        print_open_msg(address, port)
        self('serve --address %s --port %s  %s'%(address, port, options),
             debug=debug)
        print_open_msg(address, port)

    browse = serve

    def unbundle(self, bundle, update=True, options='', debug=True):
        """
        Apply patches from a hg patch to the repository.

        If the bundle is a .patch file, instead call the import_patch
        method. To see what is in a bundle before applying it, using
        self.incoming(bundle).

        INPUT:

        -  ``bundle`` - an hg bundle (created with the bundle
           command)

        -  ``update`` - if True (the default), update the
           working directory after unbundling.

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.unbundle('myhg.bundle')  # not tested
        """
        if bundle.startswith("http://") or bundle.startswith("https://"):
            if sage_trac_re.match(bundle):
                bundle = bundle.replace('sage_trac/attachment/', 'sage_trac/raw-attachment/')
            bundle = get_remote_file(bundle, verbose=True)
        if bundle[-6:] == '.patch':
            self.import_patch(bundle, options)
            return
        if bundle[-5:] == '.diff':
            return self.import_patch(bundle)
        self._ensure_safe()
        bundle = os.path.abspath(bundle)
        print "Unbundling bundle %s"%bundle
        if update:
            options = '-u'
        else:
            options = ''

        print "If you get an error 'abort: unknown parent'"
        print "this usually means either you need to do:"
        print "       hg_%s.pull()"%self.__obj_name
        print "or you're applying this patch to the wrong repository."
        self('unbundle %s "%s"'%(options, bundle), debug=debug)

    apply = unbundle

    def export(self, revs, filename=None, text=False, options='', debug=True):
        r"""
        Export patches with the changeset header and diffs for one or more
        revisions.

        If multiple revisions are given, one plain text unified diff file
        is generated for each one. These files should be applied using
        import_patch in order from smallest to largest revision number.
        The information shown in the changeset header is: author, changeset
        hash, parent and commit comment.

        .. note::

           If you are sending a patch to somebody using export and it
           depends on previous patches, make sure to include those
           revisions too! Alternatively, use the :meth:`.bundle`
           method, which includes enough information to patch against
           the default repository (but is an annoying and mysterious
           binary file).

        INPUT:

        -  ``revs`` - integer or list of integers (revision
           numbers); use the log() method to see these numbers.

        -  ``filename`` - (default: '%R.patch') The name of the file
           is given using a format string.  The formatting rules are
           as follows::

               %%   literal "%" character
               %H   changeset hash (40 bytes of hexadecimal)
               %N   number of patches being generated
               %R   changeset revision number
               %b   basename of the exporting repository
               %h   short-form changeset hash (12 bytes of hexadecimal)
               %n   zero-padded sequence number, starting at 1
               %r   zero-padded changeset revision number

        -  ``text`` - boolean (default False).  Setting this to be True
           has the same effect as passing the "-a" option below.

        -  ``options`` - string (default: '')

           - ``'-a'`` or ``'--text'`` - treat all files as text

           - ``'--switch-parent'`` -  diff against the second parent

           Without the ``-a`` option, export will avoid generating
           diffs of files it detects as binary. With ``-a``, export
           will generate a diff anyway, probably with undesirable
           results.

           With the ``--switch-parent`` option, the diff will be
           against the second parent. It can be useful to review a
           merge.

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.export('tip')  # not tested
            sage: hg_sage.export('tip', 'new.patch')  # not tested
            sage: hg_sage.export('tip', 'new.patch', options='--a')  # not tested
        """
        if filename is None:
            filename = '%R.patch'
        if not isinstance(revs, list):
            if revs == "tip":
                revs = [revs]
            else:
                revs = [int(revs)]
        if not isinstance(filename, str):
            raise TypeError, 'filename must be a string'
        if filename[-6:] != '.patch':
            filename += '.patch'
        full_path = os.path.abspath(filename)
        cwd = os.path.dirname(full_path)
        options += ' -o "%s"'%full_path
        if filename == '%R.patch':
            print "Output will be written to revision numbered file in the directory %s."%cwd
        else:
            print "Output will be written to '%s'"%full_path
        if text:
            options += ' -a'
        self('export %s %s'%(options, ' '.join([str(x) for x in revs])),
             debug=debug)

    def import_patch(self, filename, options='', debug=True):
        """
        Import an ordered set of patches from patch file, i.e., a plain
        text file created using the export command.

        If there are outstanding changes in the working directory,
        import_patch will abort unless given the -f flag.

        If imported patch was generated by the export command, user and
        description from patch override values from message headers and
        body. Values given as options with -m and -u override these.

        INPUT:

        -  ``filename`` - string

        -  ``options`` - string (default: '')::

            options: [-p NUM] [-b BASE] [-m MESSage] [-f] PATCH...

               -p --strip NUM      directory strip option for patch. This has
                                   the same meaning as the corresponding patch
                                   option (default: 1)
               -b --base PATH      base path
               -f --force          skip check for outstanding uncommitted changes
                  --no-commit      don't commit, just update the working directory
                  --exact          apply patch to the nodes from which it was
                                   generated
                  --import-branch  use any branch information in patch (implied
                                   by --exact)
               -m --message TEXT   use text as commit message
               -l --logfile FILE   read commit message from file
               -d --date DATE      record the specified date as commit date
               -u --user USER      record the specified user as committer
               -s --similarity SIMILARITY   guess renamed files by similarity (0<=s<=100)
                  --mq             operate on patch repository

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        ALIASES: patch

        EXAMPLES::

            sage: hg_sage.import_patch('trac_2001.patch')  # not tested
        """
        if filename.startswith("http://") or filename.startswith("https://"):
            filename = get_remote_file(filename, verbose=True)
        self._ensure_safe()
        self('import  %s "%s"'%(options, os.path.abspath(filename)),
             debug=debug)

    patch = import_patch

    def incoming(self, source, options='-p', debug=True):
        """
        Show new changesets found in the given source and display the
        corresponding diffs. This even works if the source is a bundle file
        (ends in .hg or .bundle). This is great because it lets you "see
        inside" the mysterious binary-only .hg files.

        Show new changesets found in the specified path/URL or the default
        pull location. These are the changesets that would be pulled if a
        pull was requested.

        For remote repository, using -bundle avoids downloading the
        changesets twice if the incoming is followed by a pull.

        See pull for valid source format details.

        ALIAS: inspect

        INPUT:

        -  ``filename`` - string, may be a URL

        -  ``options`` - (default: '-p')::

            string '[-p] [-n] [-M] [-r REV] ...'
              -M --no-merges       do not show merges
              -f --force           run even when remote repository is unrelated
              --style              display using template map file
              -n --newest-first    show newest record first
              --bundle             file to store the bundles into
              -p --patch           show patch
              -r --rev             a specific revision you would like to pull
              --template           display with template
              -e --ssh             specify ssh command to use
              --remotecmd          specify hg command to run on the remote side

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.incoming("http://hg.sagemath.org/sage-main") # not tested
            cd ... && hg incoming "http://hg.sagemath.org/sage-main"  --config pager.pager="LESS='R' less"
            comparing with http://hg.sagemath.org/sage-main
            searching for changes
            no changes found
        """
        if source.startswith("http://") or source.startswith("https://"):
            source = get_remote_file(source, verbose=True)
        if os.path.exists(source):
            source = os.path.abspath(source)
        if os.path.splitext(source)[1] in ['.hg', '.bundle']:
            source = 'bundle://%s'%source
        self('incoming %s "%s" %s'%(options, source, pager()), debug=debug)

    inspect = incoming


    def add(self, files, options='', debug=True):
        """
        Add the given list of files (or file) or directories to your HG
        repository. They must exist already.

        To see a list of files that haven't been added to the repository do
        self.status(). They will appear with an explanation point next
        them.

        Add needs to be called whenever you add a new file or directory to
        your project. Of course, it also needs to be called when you first
        create the project, to let hg know which files should be kept track
        of.

        INPUT:

        -  ``files`` - list or string; name of file or
           directory.

        -  ``options`` - string (e.g., '--dry-run')

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.add('module_list.pyc', options='--dry-run')
            Adding file module_list.pyc
            cd ... && hg add --dry-run "module_list.pyc"
        """
        if isinstance(files, str):
            if ' ' in files:
                files = files.split()
            else:
                files = [files]
        for file in files:
            print "Adding file %s"%file
            self('add %s "%s"'%(options, file), debug=debug)

    def remove(self, files, options='', debug=True):
        """
        Remove the given list of files (or file) or directories from your
        HG repository.

        INPUT:

        -  ``files`` - list or string; name of file or
           directory.

        -  ``options`` - string (e.g., '-f')

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.remove('sage/misc/remove_me.py') # not tested
            Removing file sage/misc/remove_me.py
            cd ... && hg rm "sage/misc/remove_me.py"
        """
        if isinstance(files, str):
            files = [files]
        for file in files:
            print "Removing file %s"%file
            self('rm %s "%s"'%(options, file), debug=debug)

    rm = remove

    def rename(self, src, dest, options='', debug=True):
        """
        Move (rename) the given file, from src to dest. This command takes
        effect in the next commit.

        INPUT:

        -  ``src, dest`` - strings that define a file, relative
           to self.dir()

        -  ``options``::

               -A --after    record a rename that has already occurred
               -f --force    forcibly copy over an existing managed file
               -n --dry-run  do not perform actions, just print output

        -  ``debug`` - boolean (default True); if True, print the full
           system command being executed.

        EXAMPLES::

            sage: hg_sage.rename('sage/misc/hg.py', 'sage/misc/hgnew.py', options='--dry-run')
            Moving sage/misc/hg.py --> sage/misc/hgnew.py
            cd ... && hg mv --dry-run "sage/misc/hg.py" "sage/misc/hgnew.py"
        """
        print "Moving %s --> %s"%(src,dest)
        self('mv %s "%s" "%s"'%(options, src,dest), debug=debug)

    move = rename
    mv = rename

    def log(self, branches=None, keyword=None, limit=None,
                  rev=None, merges=True, only_merges=False,
                  patch=None, template=False, include=None,
                  exclude=None, verbose=False, debug=True):
        """
        Display the change log for this repository. This is a list of
        changesets ordered by revision number.

        By default this command outputs: changeset id and hash, tags,
        non-trivial parents, user, date and time, and a summary for each
        commit.

        INPUT:

        -  ``branches`` - (string, default: None) show given
           branches

        -  ``keyword`` - (string, default: None) search for a
           keyword

        -  ``limit`` - (integer, default: None, or 20 in
           notebook mdoe) limit number of changes displayed

        -  ``rev`` - (integer) show the specified revision

        -  ``merges`` - (bool, default: True) whether or not
           to show merges

        -  ``only_merges`` - (bool, default: False) if true,
           show only merges

        -  ``patch`` - (string, default: None) show given
           patch

        -  ``template`` - (string, default: None) display with
           template

        -  ``include`` - (string, default: None) include names
           matching the given patterns

        -  ``exclude`` - (string, default: None) exclude names
           matching the given patterns

        -  ``verbose`` - (bool, default: False) If true, the
           list of changed files and full commit message is shown.

        -  ``debug`` - boolean (default True); if True, print the
           full system command being executed.

        EXAMPLES::

            sage: hg_sage.log() # not tested
            cd ... && hg log   --config pager.pager="LESS='R' less"
            ...
        """
        if embedded() and limit is None:
            limit = 20
        options = ''
        if branches:
            options += '-b %s '%branches
        if keyword:
            options += '-k "%s" '%keyword
        if limit:
            options += '-l %s '%limit
        if rev:
            options += '-r %s '%rev
        if not merges:
            options += '--no-merges '
        if only_merges:
            options += '-m '
        if patch:
            options += '-p "%s"'%patch
        if template:
            options += '--template'
        if include:
            options += '-I "%s"'%include
        if exclude:
            options += '-X "%s"'%exclude
        if verbose:
            options = '-v ' + options

        self('log %s %s %s'%(options, color(), pager()), debug=debug)

    changes = log
    history = log

    def diff(self, files='', rev=None, options='', debug=True):
        """
        Show differences between revisions for the specified files as a
        unified diff.

        By default this command tells you exactly what you have changed in
        your working repository since you last committed changes.

        INPUT:

        -  ``files`` - space separated list of files (relative
           to self.dir())

        -  ``rev`` - None or a list of integers.

        - ``options`` -- string (default '').  Some possibilities::

            -a --text                 treat all files as text
            -g --git                  use git extended diff format
            -p --show-function        show which function each change is in
            -I --include PATTERN [+]  include names matching the given patterns
            -X --exclude PATTERN [+]  exclude names matching the given patterns

        -  ``debug`` - boolean (default True); if True, print the
           full system command being executed.

        Differences between files are shown using the unified diff format.

        When two revision arguments are given, then changes are shown
        between those revisions. If only one revision is specified then
        that revision is compared to the working directory, and, when no
        revisions are specified, the working directory files are compared
        to its parent.

        EXAMPLES::

            sage: hg_sage.diff()
            cd ... && hg diff    --config pager.pager="LESS='R' less"

        To see the changes in this file since revision 10000:

            sage: hg_sage.diff('sage/misc/hg.py', rev=10000) # not tested
            cd ... && hg diff  -r 10000  sage/misc/hg.py  --config pager.pager="LESS='R' less"
            ...
        """
        if not rev is None:
            if not isinstance(rev, (list, tuple)):
                rev = [rev]
            extra_options = ' '.join(['-r %s'%r for r in rev]) + '  ' + files
        else:
            extra_options = files
        self('diff %s %s %s %s'%(options, extra_options, color(), pager()),
             debug=debug)

    what = diff

    def revert(self, files='', options='', rev=None, debug=True):
        """
        Revert files or dirs to their states as of some revision

        .. note::

            This command is most likely not what you are looking
            for. ``revert`` will partially overwrite content in the
            working directory without changing the working directory
            parents.  Use the method ``update(options='-r REV')`` to
            check out earlier revisions, or ``update(options='--clean
            .')`` to undo a merge which has added another parent.

        With no revision specified, revert the named files or
        directories to the contents they had in the parent of the
        working directory. This restores the contents of the affected
        files to an unmodified state and unschedules adds, removes,
        copies, and renames. If the working directory has two parents,
        you must explicitly specify a revision.

        Using the ``rev`` argument, revert the given files or
        directories to their contents as of a specific revision. This
        can be helpful to "roll back" some or all of an earlier change.
        Run ``hg_sage('help dates')`` for a list of formats valid for
        the ``-d/--date`` option.

        Revert modifies the working directory. It does not commit any
        changes, or change the parent of the working directory. If you
        revert to a revision other than the parent of the working
        directory, the reverted files will thus appear modified
        afterwards.

        If a file has been deleted, it is restored. If the executable
        mode of a file was changed, it is reset.

        If names are given, all files matching the names are
        reverted. If no arguments are given, no files are reverted.
        To revert all files in the repository, pass the argument
        ``options='--all'``.

        Modified files are saved with a .orig suffix before
        reverting. To disable these backups, use
        ``options='--no-backup'``.

        If ``debug`` is True, also print the full system command being
        executed.

        OPTIONS::

          -a --all                  revert all changes when no arguments given
          -d --date DATE            tipmost revision matching date
             --no-backup            do not save backup copies of files
          -I --include PATTERN [+]  include names matching the given patterns
          -X --exclude PATTERN [+]  exclude names matching the given patterns
          -n --dry-run              do not perform actions, just print output
             --mq                   operate on patch repository

        EXAMPLES::

            sage: hg_sage.revert('sage/misc/hg.py', rev=12000, options='--dry-run')
            cd ... && hg revert --dry-run -r 12000 sage/misc/hg.py
        """
        if not rev is None:
            options = options +' -r %s %s'%(rev, files)
        else:
            options = options + files
        self('revert %s'%options, debug=debug)

    def dir(self):
        """
        Return the directory where this repository is located.

        EXAMPLES::

            sage: os.path.realpath(hg_sage.dir()).startswith(os.path.realpath(os.environ['SAGE_ROOT']))
            True
        """
        return self.__dir

    def pull_url(self):
        """
        Return the default 'master url' for this repository.

        EXAMPLES::

            sage: hg_sage.pull_url()
            'http://hg.sagemath.org/sage-main/'
        """
        return self.__pull_url

    def push_url(self):
        """
        Return the default url for uploading this repository.

        EXAMPLES::

            sage: hg_sage.push_url()
            'http://hg.sagemath.org/sage-main/'
        """
        return self.__push_url


    def help(self, cmd='', debug=True):
        r"""
        Print a help message about ``cmd``, or if ``cmd`` is omitted,
        print a general Mercurial help message.

        If ``debug`` is True, also print the full system command being
        executed.

        If this hg object is called hg_sage, then you call a command using
        ``hg_sage('usual hg command line notation')``.  Type "hg_sage?" for
        more information.

        EXAMPLES::

            sage: hg_sage.help()  # not tested
            Mercurial Distributed SCM

            list of commands:
            ...
            sage: hg_sage.help('status')   # not tested
            hg status [OPTION]... [FILE]...

            aliases: st

            show changed files in the working directory
            ...
        """
        self('%s --help %s'%(cmd, pager()), debug=debug)

    def outgoing(self, url=None, opts='', debug=True):
        """
        Use this to find changsets that are in your branch, but not in the
        specified destination repository. If no destination is specified,
        the official repository is used. By default, push_url() is used.

        From the Mercurial documentation:

            Show changesets not found in the specified destination
            repository or the default push location.  These are the
            changesets that would be pushed if a push was requested.

            See push() for valid destination format details.

        INPUT:


        -  ``url`` - (Default: self.push_url())  the official
           repository

           - ``http://[user@]host[:port]/[path]``

           - ``https://[user@]host[:port]/[path]``

           - ``ssh://[user@]host[:port]/[path]``

           - local directory (starting with a /)

           - name of a branch (for hg_sage); no /'s

        - ``options`` - (Default: None)::

              -M --no-merges     do not show merges
              -f --force         run even when remote repository is unrelated
              -p --patch         show patch
              --style            display using template map file
              -r --rev           a specific revision you would like to push
              -n --newest-first  show newest record first
              --template         display with template
              -e --ssh           specify ssh command to use
              --remotecmd        specify hg command to run on the remote side

        - ``debug`` - if True, print the full system command being
          executed.

        EXAMPLES::

            sage: hg_sage.outgoing() # not tested
            cd ... && hg outgoing  http://hg.sagemath.org/sage-main/ --config pager.pager="LESS='R' less"
            comparing with http://hg.sagemath.org/sage-main/
            searching for changes
            ...
        """
        if url is None:
            url = self.__push_url

        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('outgoing %s %s %s' % (opts, url, pager()), debug=debug)

    def pull(self, url=None, options='-u', debug=True):
        """
        Pull all new patches from the repository at the given url, or use
        the default 'official' repository if no url is specified.

        INPUT:

        -  ``url`` - (Default: self.pull_url())  the official
           repository

           - ``http://[user@]host[:port]/[path]``

           - ``https://[user@]host[:port]/[path]``

           - ``ssh://[user@]host[:port]/[path]``

           - local directory (starting with a /)

           - name of a branch (for hg_sage); no /'s

        - ``options`` - (Default: '-u')::

              -u --update     update the working directory to tip after pull
              -e --ssh        specify ssh command to use
              -f --force      run even when remote repository is unrelated
              -r --rev        a specific revision you would like to pull
              --remotecmd     specify hg command to run on the remote side

        - ``debug`` - if True, print the full system command being
          executed.

        Some notes about using SSH with Mercurial:

        - SSH requires an accessible shell account on the destination
          machine and a copy of hg in the remote path or specified
          with as remotecmd.

        - path is relative to the remote user's home directory by
          default. Use an extra slash at the start of a path to
          specify an absolute path: ``ssh://example.com//tmp/repository``

        - Mercurial doesn't use its own compression via SSH; the right
          thing to do is to configure it in your /.ssh/ssh_config,
          e.g.::

              Host *.mylocalnetwork.example.com
                Compression off
              Host *
                Compression on

          Alternatively specify 'ssh -C' as your ssh command in your
          hgrc or with the -ssh command line option.

        EXAMPLES::

            sage: hg_sage.pull() # not tested
            cd ... && hg pull http://hg.sagemath.org/sage-main/
            ...
        """
        self._ensure_safe()

        if url is None:
            url = self.__pull_url
        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('pull %s %s'%(options, url), debug=debug)
        if self.__target == 'sage':
            print ""
            print "Now building the new Sage libraries"
            os.system('sage -b')
            print "You *MUST* restart Sage in order for the changes to take effect!"

        print "If it says use 'hg merge' above, then you should"
        print "type hg_%s.merge()."%self.__obj_name

    def push(self, url=None, options='', debug=True):
        """
        Push all new patches from the repository to the given destination.

        INPUT:

        -  ``url`` - (Default: self.push_url())  the official
           repository

           - ``http://[user@]host[:port]/[path]``

           - ``https://[user@]host[:port]/[path]``

           - ``ssh://[user@]host[:port]/[path]``

           - local directory (starting with a /)

           - name of a branch (for hg_sage); no /'s

        - ``options`` - (Default: '')::

              -e --ssh        specify ssh command to use
              -f --force      run even when remote repository is unrelated
              -r --rev        a specific revision you would like to pull
              --remotecmd     specify hg command to run on the remote side

        - ``debug`` - if True, print the full system command being
          executed.

        Some notes about using SSH with Mercurial:

        - SSH requires an accessible shell account on the destination
          machine and a copy of hg in the remote path or specified
          with as remotecmd.

        - path is relative to the remote user's home directory by
          default. Use an extra slash at the start of a path to
          specify an absolute path: ``ssh://example.com//tmp/repository``

        - Mercurial doesn't use its own compression via SSH; the right
          thing to do is to configure it in your /.ssh/ssh_config,
          e.g.::

              Host *.mylocalnetwork.example.com
                Compression off
              Host *
                Compression on

          Alternatively specify 'ssh -C' as your ssh command in your
          hgrc or with the -ssh command line option.

        EXAMPLES::

            sage: hg_sage.push() # not tested
            cd ... && hg push http://hg.sagemath.org/sage-main/
            ...
        """
        self._ensure_safe()

        if url is None:
            url = self.__push_url
        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('push %s %s'%(options, url), debug=debug)


    def merge(self, options='', debug=True):
        """
        Merge working directory with another revision

        Merge the contents of the current working directory and the
        requested revision. Files that changed between either parent are
        marked as changed for the next commit and a commit must be
        performed before any further updates are allowed.

        INPUT:

        -  ``options`` - default: ''::

             -f --force  force a merge with outstanding changes
             -r --rev    revision to merge

        - ``debug`` - if True, print the full system command being
          executed.

        EXAMPLES::

            sage: hg_sage.merge() # not tested
            cd ... && hg merge
        """
        self('merge %s'%options, debug=debug)

    def update(self, options='', debug=True):
        """
        update or merge working directory

        Update the working directory to the specified revision.

        If there are no outstanding changes in the working directory and
        there is a linear relationship between the current version and the
        requested version, the result is the requested version.

        To merge the working directory with another revision, use the merge
        command.

        By default, update will refuse to run if doing so would require
        merging or discarding local changes.

        aliases: up, checkout, co

        INPUT:

        -  ``options`` - string (default: '')::

            -C --clean  overwrite locally modified files
            -d --date   tipmost revision matching date
            -r --rev    revision

        - ``debug`` - if True, print the full system command being
          executed.

        EXAMPLES::

            sage: hg_sage.update() # not tested
            cd ... && hg update
        """
        self('update %s'%options, debug=debug)

    up = update
    checkout = update
    co = update

    def head(self, options='', debug=True):
        """
        show current repository heads

        Show all repository head changesets.

        Repository "heads" are changesets that don't have children
        changesets. They are where development generally takes place and
        are the usual targets for update and merge operations.

        INPUT:

        -  ``options`` - string (default: '')::

             -r --rev       show only heads which are descendants of rev
                --style     display using template map file
                --template  display with template

        - ``debug`` - if True, print the full system command being
          executed.

        EXAMPLES::

            sage: hg_sage.head() # random
            cd ... && hg head
            changeset:   15825:6ca08864b80c
            tag:         tip
            user:        Jeroen Demeyer <jdemeyer@cage.ugent.be>
            date:        Tue Jun 07 12:31:57 2011 +0000
            summary:     4.7.1.alpha2
        """
        self('head %s'%options, debug=debug)

    heads = head

    def switch(self, name=None):
        r"""
        Switch to a different branch. You must restart Sage after
        switching.

        Only available for ``hg_sage.``

        INPUT:

        -  ``name`` - name of a Sage branch (default: None)

        If the name is not given, this function returns a list of all
        branches.

        EXAMPLES::

            sage: hg_sage.switch() # random
            ['main']
            sage: hg_sage.switch('new') # not tested
            <BLANKLINE>
            ----------------------------------------------------------
            Building and installing modified Sage library files.
            <BLANKLINE>
            ...
        """
        if name is None:
            s = os.popen('ls -l %s/devel/ |grep sage-'%os.environ['SAGE_ROOT']).read()
            t = s.split('\n')
            v = []
            for X in t:
                i = X.rfind('sage-')
                n = X[i+5:]
                if n != '':
                    v.append(n)
            v = list(set(v))
            v.sort()
            return v
        os.system('sage -b "%s"'%name)

    def clone(self, name, rev=None):
        r"""
        Clone the current branch of the Sage library, and make it active.

        Only available for the ``hg_sage`` repository.

        Use ``hg_sage.switch('branch_name')`` to switch to a
        different branch. You must restart Sage after switching.

        INPUT:

        -  ``name`` - string

        -  ``rev`` - integer or None (default)


        If rev is None, clones the latest recorded version of the
        repository. This is very fast, e.g., about 30-60 seconds (including
        any build). If a specific revision is specified, cloning may take
        much longer (e.g., 5 minutes), since all Pyrex code has to be
        regenerated and compiled.

        EXAMPLES:

        Make a clone of the repository called testing. A copy of the
        current repository will be created in a directory sage-testing,
        then SAGE_ROOT/devel/sage will point to sage-testing, and when you
        next restart Sage that's the version you'll be using.

        ::

            sage: hg_sage.clone('testing')    # not tested
            ...

        Make a clone of the repository as it was at revision 1328.

        ::

            sage: hg_sage.clone('testing', 1328)    # not tested
            ...
        """
        if not self.__cloneable:
            raise RuntimeError, "only available for hg_sage"
        name = '_'.join(str(name).split())
        if rev is None:
            os.system('sage -clone %s'%name)
        else:
            os.system('sage -clone %s -r %s'%(name, int(rev)))

    def commit(self, files='', comment=None, options='', diff=True,
               debug=True):
        r"""
        Commit your changes to the repository.

        Quit out of the editor without saving to not record your changes.

        INPUT:

        - ``files`` - space separated string of file names (optional)
          If specified only those files are committed. The path must be
          absolute or relative to self.dir().

        - ``comment`` - optional changeset comment. If you don't give
           it you will be dumped into an editor. If you're using the
           Sage notebook, you *must* specify a comment.

        - ``options`` - string::

              -A --addremove  mark new/missing files as added/removed before committing
              -m --message    use <text> as commit message
              -l --logfile    read the commit message from <file>
              -d --date       record datecode as commit date
              -u --user       record user as committer
              -I --include    include names matching the given patterns
              -X --exclude    exclude names matching the given patterns

        - ``diff`` - (default: True) if True show diffs between your repository
          and your working repository before recording changes.

        - ``debug`` - if True, print the full system command being
          executed.

        .. note::

           If you create new files you should first add them with the
           add method.

        EXAMPLES::

            sage: hg_sage.commit('hg.py', comment='miscellaneous fixes') # not tested
            cd ... && hg commit -m "miscellaneous fixes" hg.py
        """
        if embedded() and comment is None:
            raise RuntimeError, "You're using the Sage notebook, so you *must* explicitly specify the comment in the commit command."
        if diff:
            self.diff(files)

        if isinstance(files, (list, tuple)):
            files = ' '.join([str(x) for x in files])

        if comment:
            self('commit %s -m "%s" %s '%(options, comment, files), debug=debug)
        else:
            self('commit %s %s'%(options, files), debug=debug)

    record = commit
    ci = commit

    def rollback(self, debug=True):
        """
        Remove recorded patches without changing the working copy.

        If ``debug`` is True, also print the full system command being
        executed.

        EXAMPLES::

            sage: hg_sage.rollback() # not tested
            cd ... && hg rollback
        """
        self('rollback', debug=debug)

    def bundle(self, filename, options='', url=None, base=None, to=None,
               debug=True):
        r"""
        Create an hg changeset bundle with the given filename against the
        repository at the given url (which is by default the 'official'
        Sage repository, unless push_url() is changed in a setup file).

        If you have internet access, it's best to just do
        ``hg_sage.bundle(filename)``. If you don't find a
        revision r that you and the person unbundling both have (by looking
        at ``hg_sage.log()``), then do
        ``hg_sage.bundle(filename, base=r)``.

        Use self.inspect('file.bundle') to inspect the resulting bundle.

        This is a file that you should probably send to William Stein
        (wstein@gmail.com), post to a web page, or send to sage-devel. It
        will be written to the current directory.

        INPUT:

        -  ``filename`` - output file in which to put bundle

        -  ``options`` - pass to hg

        -  ``url`` - url to bundle against (default:
           SAGE_SERVER, or push_url())

        -  ``base`` - a base changeset revision number to
           bundle against (doesn't require internet access)

        -  ``debug`` - if True, print the full system command being
           executed.

        EXAMPLES::

            sage: hg_sage.bundle('new-bundle') # not tested
            Writing to /.../new-bundle.hg
            cd ... && hg bundle tmphg http://hg.sagemath.org/sage-main/
            searching for changes
            133 changesets found
            Successfully created hg patch bundle /.../new-bundle.hg
        """
        if not base is None:
            url = ''
            options = '--base=%s %s'%(int(base), options)

        if url is None:
            url = self.__push_url

        # make sure that we don't accidentally create a file ending in '.hg.hg'
        if filename[-3:] == '.hg':
            filename = filename[:-3]
        # We write to a local tmp file, then move, since under
        # windows hg has a bug that makes it fail to write
        # to any filename that is at all complicated!
        filename = os.path.abspath(filename)
        if filename[-3:] != '.hg':
            filename += '.hg'
        print 'Writing to %s'%filename
        tmpfile = '%s/tmphg'%self.__dir
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        self('bundle %s tmphg %s'%(options, url), debug=debug)
        if os.path.exists(tmpfile):
            shutil.move(tmpfile, filename)
            print 'Successfully created hg patch bundle %s'%filename
            if not to is None:
                os.system('scp "%s" %s'%(filename, to))
        else:
            print 'Problem creating hg patch bundle %s'%filename

    send = bundle
    save = send

    # Mercurial queues

    def qseries(self, verbose=False, debug=True):
        """
        Mercurial queues: print the series file.

        If optional argument ``verbose`` is True, then also print the
        first line of each patch's header.

        If ``debug`` is True, also print the full system command being
        executed.

        EXAMPLES::

            sage: hg_sage.qseries()
            cd ... && hg qseries
        """
        options = "--summary" if verbose else ""
        self('qseries %s %s' % (options, color(),), debug=debug)

    def qapplied(self, verbose=False, debug=True):
        """
        Mercurial queues: print the patches in the queue which have
        been applied.

        If optional argument ``verbose`` is True, then also print the
        first line of each patch's header.

        If ``debug`` is True, also print the full system command being
        executed.

        EXAMPLES::

            sage: hg_sage.qapplied()
            cd ... && hg qapplied
        """
        options = "--summary" if verbose else ""
        self('qapplied %s %s' % (options, color(),), debug=debug)

    def qunapplied(self, verbose=False, debug=True):
        """
        Mercurial queues: print the patches in the queue which have
        not yet been applied.

        If optional argument ``verbose`` is True, then also print the
        first line of each patch's header.

        If ``debug`` is True, also print the full system command being
        executed.

        EXAMPLES::

            sage: hg_sage.qunapplied()
            cd ... && hg qunapplied
        """
        options = "--summary" if verbose else ""
        self('qunapplied %s %s' % (options, color(),), debug=debug)

    def qimport(self, filename, options='', debug=True):
        """
        Mercurial queues: insert the patch from ``filename`` in the queue.

        INPUT:

        - ``filename`` -- string
        - ``options`` -- string (default '')::

          -e --existing     import file in patch directory
          -n --name NAME    name of patch file
          -f --force        overwrite existing files
          -r --rev REV [+]  place existing revisions under mq control
          -g --git          use git extended diff format
          -P --push         qpush after importing

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qimport('old.patch') # not tested
            cd ... && hg qimport old.patch ...
        """
        if filename.startswith("http://") or filename.startswith("https://"):
            filename = get_remote_file(filename, verbose=True)
        self._ensure_safe()
        self('qimport %s %s' % (options, os.path.abspath(filename)),
             debug=debug)

    def qdelete(self, patches, keep_patch=False, debug=True):
        """
        Mercurial queues: delete the named patch from the queue.

        INPUT:

        - ``patches`` -- string, a patch or list of patches separated by
          white space

        - ``keep_patch`` -- (default False): if True, keep the patch
          files in the patch directory.

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qdelete('old.patch new.patch') # not tested
            cd ... && hg qdelete old.patch new.patch ...
        """
        options = "--keep" if keep_patch else ""
        self('qdelete %s %s' % (options, patches),
             debug=debug)

    qremove = qdelete

    def qpush(self, force=False, all=False, options='', debug=True):
        """
        Mercurial queues: push the next patch onto the stack.

        INPUT:

        - ``force`` -- boolean (default False): if True, apply even if the
          patch has rejects.

        - ``all`` -- boolean (default False): if True, apply all unapplied
          patches.

        - ``options`` -- string (default ''): extra options to pass to the command.

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qpush() # not tested
            cd ... && hg  qpush
        """
        extra_options = ''
        if force:
            extra_options += ' --force'
        if all:
            extra_options += ' --all'
        self('qpush %s %s' % (extra_options, options),
             debug=debug)

    def qpop(self, patch='', force=False, all=False, debug=True):
        """
        Mercurial queues: pop the top of the patch stack, or if given
        a patch, pop patches off of the stack until it is at the top.

        INPUT:

        - ``patch`` -- string (default ''): if nonempty, this should
          name an applied patch, and then the command will pop patches
          off the stack until ``patch`` is at the top.

        - ``force`` -- boolean (default False): if True, forget any
          changes to patched files.

        - ``all`` -- boolean (default False): if True, pop all
          patches.

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qpop() # not tested
            cd ... && hg qpop
        """
        extra_options = ''
        if force:
            extra_options += ' --force'
        if all:
            extra_options += ' --all'
        if patch:
            extra_options += ' %s' % patch
        self('qpop %s' % extra_options, debug=debug)

    def qrefresh(self, options='', debug=True):
        """
        Mercurial queues: update the current patch.

        INPUT:

        - ``options`` -- string (default '').  Some possibilities::

           -e --edit                 edit commit message
           -g --git                  use git extended diff format
           -I --include PATTERN [+]  include names matching the given patterns
           -X --exclude PATTERN [+]  exclude names matching the given patterns
           -m --message TEXT         use text as commit message

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qrefresh() # not tested
            cd ... && hg qrefresh
        """
        self('qrefresh %s' % options, debug=debug)

    def qdiff(self, options='', debug=True):
        """
        Mercurial queues: show a diff including the current patch and
        any more recent changes, thus showing what the current patch
        would become after :meth:`qrefresh`.

        Use :meth:`diff` if you only want to see the changes made
        since the last :meth:`qrefresh`.

        INPUT:

        - ``options`` -- string (default '').  Some possibilities::

           -a --text                 treat all files as text
           -g --git                  use git extended diff format
           -p --show-function        show which function each change is in
           -I --include PATTERN [+]  include names matching the given patterns
           -X --exclude PATTERN [+]  exclude names matching the given patterns

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qdiff()
            cd ... && hg qdiff ...
        """
        self('qdiff %s %s %s' % (options, color(), pager()), debug=debug)

    def qnew(self, patch, options='', debug=True):
        """
        Mercurial queues: create a new patch.

        INPUT:

        - ``patch`` -- string, the name of the patch.

        - ``options`` -- string (default '').  Some possibilities::

           -e --edit                 edit commit message
           -g --git                  use git extended diff format
           -I --include PATTERN [+]  include names matching the given patterns
           -X --exclude PATTERN [+]  exclude names matching the given patterns
           -m --message TEXT         use text as commit message

        - ``debug`` -- (default True): if True, print the full system
          command being executed.

        EXAMPLES::

            sage: hg_sage.qnew('my_new.patch') # not tested
            cd ... && hg qnew my_new.patch
        """
        self('qnew %s %s' % (patch, options), debug=debug)


##############################################################################
# Initialize the actual repositories.
##############################################################################

import misc

SAGE_ROOT = misc.SAGE_ROOT
DEFAULT_SERVER = "http://hg.sagemath.org"

SAGE_INCOMING_SERVER = os.getenv("SAGE_INCOMING_SERVER")
if SAGE_INCOMING_SERVER is None:
    try:
        SAGE_INCOMING_SERVER = os.environ['SAGE_HG_SERVER'].strip('/')
    except KeyError:
        #print "Falling back to a hard coded sage server in misc/hg.py"
        SAGE_INCOMING_SERVER = DEFAULT_SERVER

SAGE_OUTGOING_SERVER = os.getenv("SAGE_OUTGOING_SERVER")
if SAGE_OUTGOING_SERVER is None:
    SAGE_OUTGOING_SERVER = SAGE_INCOMING_SERVER

if (SAGE_INCOMING_SERVER == DEFAULT_SERVER):      ## Always uses the "main" branch on the default server.
    temp_in_branch_name = "main"
else:
    temp_in_branch_name = branch_current_hg()

if (SAGE_OUTGOING_SERVER == DEFAULT_SERVER):      ## Always uses the "main" branch on the default server.
    temp_out_branch_name = "main"
else:
    temp_out_branch_name = branch_current_hg()


if (SAGE_INCOMING_SERVER != DEFAULT_SERVER) or (SAGE_OUTGOING_SERVER != DEFAULT_SERVER):
    print "Non-default server settings detected:"
    print "    Incoming Server = " + SAGE_INCOMING_SERVER + ''.join(["  (default)"  \
                for i in range(1)  if (SAGE_INCOMING_SERVER == DEFAULT_SERVER)])
    print "    Outgoing Server = " + SAGE_OUTGOING_SERVER + ''.join(["  (default)"  \
                for i in range(1)  if (SAGE_OUTGOING_SERVER == DEFAULT_SERVER)])
    print


hg_sage    = HG('%s/devel/sage'%SAGE_ROOT,
                'Sage Library Source Code',
                    pull_url='%s/sage-%s/'%(SAGE_INCOMING_SERVER, temp_in_branch_name),
                    push_url='%s/sage-%s/'%(SAGE_OUTGOING_SERVER, temp_out_branch_name),
                cloneable=True,
                obj_name='sage')

hg_scripts = HG('%s/local/bin/'%SAGE_ROOT,
                'Sage Scripts',
                pull_url='%s/scripts-main/'%SAGE_INCOMING_SERVER,
                push_url='%s/scripts-main/'%SAGE_OUTGOING_SERVER,
                obj_name='scripts')

hg_extcode = HG('%s/data/extcode'%SAGE_ROOT,
                'Sage External System Code (e.g., PARI, MAGMA, etc.)',
                pull_url='%s/extcode-main/'%SAGE_INCOMING_SERVER,
                push_url='%s/extcode-main/'%SAGE_OUTGOING_SERVER,
                obj_name='extcode')


hg_examples = HG('%s/data/examples'%SAGE_ROOT,
                 'Sage Examples',
                 pull_url='%s/examples/'%SAGE_INCOMING_SERVER,
                 push_url='%s/examples/'%SAGE_OUTGOING_SERVER,
                 obj_name='examples')

hg_root = HG(SAGE_ROOT,
             'Sage Root',
             pull_url=SAGE_INCOMING_SERVER,
             push_url=SAGE_OUTGOING_SERVER,
             obj_name='root')

hg_sagenb = HG('%s/devel/sagenb' % SAGE_ROOT,
               'SageNB Source Code',
               pull_url='http://boxen.math.washington.edu:8123',
               push_url='http://boxen.math.washington.edu:8123',
               obj_name='sagenb')
