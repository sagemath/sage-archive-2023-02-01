"""
HG from SAGE.

These functions make setup and use of source control with SAGE easier, using
the distributed Mercurial HG source control system.  To learn about Mercurial,
see http://www.selenic.com/mercurial/wiki/.

This system should all be fully usable from the SAGE notebook (except
for merging, currently).
"""

########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import os, shutil

import sage.server.support
from   viewer import browser
from   misc   import tmp_filename, branch_current_hg

def pager():
    if sage.server.support.EMBEDDED_MODE:
        return 'cat'
    else:
        return 'more'

class HG:
    r"""
    This is an HG (Mercurial) repository.

    To learn about Mercurial, see http://www.selenic.com/mercurial/wiki/.

    This system should all be fully usable from the SAGE notebook.

    The few of the simplest and most useful commands are directly
    provided as member functions.  However, you can use the full
    functionality of hg by noting that typing, e.g.,
            \code{hg_sage("command line arguments")}
    is \emph{exactly} the same as typing
    \begin{verbatim}
            cd <SAGE_ROOT>/devel/sage/ && hg command line arguments
    \end{verbatim}
    """
    def __init__(self, dir, name, url, target=None, cloneable=False):
        """
        INPUT:
            dir -- directory that will contain the repository
            name -- a friendly name for the repository (only used for printing)
            url -- a default URL to pull or record sends against (e.g.,
                   this could be a master repository on modular.math.washington.edu)
            target -- if the last part of dir is, e.g., sage-hg,
                      create a symlink from sage-hg to target.
                      If target=None, this symlink will not be created.
        """
        self.__dir = os.path.abspath(dir)
        self.__name = name
        self.__url = url
        self.__initialized = False
        self.__target = target
        self.__cloneable = cloneable

    def __repr__(self):
        self.status()
        return "Hg repository '%s' in directory %s"%(self.__name, self.__dir)

    def status(self):
        print("Status of modified or unknown files:")
        self('status')
        print "\n---\n"
        if self.__name == "SAGE Library Source Code":
            b = branch_current_hg()
            if b == '': b='main'
            print("Branch: %s"%b)


    def _warning(self):
        if not os.path.exists(os.environ['HOME'] + '/.hgrc'):
            print "\nWARNING:"
            print "Make sure to create a ~/.hgrc file:"
            print "-"*70
            print "[ui]"
            print "username = William Stein <wstein@gmail.com>"
            print "-"*70
            print "\n"

    def __call__(self, cmd, check_initialized=True):
        """
        Run 'hg cmd' where cmd is an arbitrary string
        in the hg repository.
        """
        self._warning()
        s = 'cd "%s" && hg %s'%(self.__dir, cmd)
        print s
        return os.system(s)

    def serve(self, port=8200, open_viewer=False):
        """
        Start a web server for this repository.

        This server is very nice -- you can browse all files in the
        repository, see their changelogs, see who wrote any given
        line, etc.  Very nice.

        INPUT:
            port -- port that the server will listen on
            open_viewer -- boolean (default: False); whether to pop up the web page
        """
        print('Now serving repository on port %s'%port)
        print("Point your web browser at http://localhost:%s"%port)
        if open_viewer:
            cmd = 'sleep 1; %s http://%s:%s 1>&2 >/dev/null'%(browser(), 'localhost', port)
            t = tmp_filename()
            open(t,'w').write(cmd)
            os.system('source %s &'%(os.path.abspath(t)))
        self('serve --port %s'%port)

    browse = serve

    def unbundle(self, bundle, update=True, options=''):
        """
        Apply patches from a hg patch to the repository.

        If the bundle is a .patch file, instead call the import_patch method.

        INPUT:
             bundle -- an hg bundle (created with the bundle command)
             update -- if True (the default), update the working directory after unbundling.
        """
        if bundle[-6:] == '.patch':
            self.import_patch(bundle, options)
            return
        bundle = os.path.abspath(bundle)
        print "Unbundling bundle %s"%bundle
        if update:
            options = '-u'
        else:
            options = ''
        self('unbundle %s %s "%s"'%(options, options, bundle))

    apply = unbundle

    def export(self, revs, filename='%R.patch', text=False, options=''):
        r"""
        Export patches with the changeset header and diffs for one or
        more revisions.

        If multiple revisions are given, one plain text unified diff
        file is generated for each one.  These files should be applied
        using import_patch in order from smallest to largest revision
        number.  The information shown in the changeset header is:
        author, changeset hash, parent and commit comment.

        \note{If you are sending a patch to somebody using export and
        it depends on previous patches, make sure to include those
        revisions too!  Alternatively, use the bundle() method, which
        includes enough information to patch against the default
        repository (but is an annoying and mysterious binary file).}

        INPUT:
             revs -- integer or list of integers (revision numbers); use the log()
                     method to see these numbers.
             filename -- (default: '%s-%H.export') The name of the file is given using a format
                 string.  The formatting rules are as follows:
                    %%   literal "%" character
                    %H   changeset hash (40 bytes of hexadecimal)
                    %N   number of patches being generated
                    %R   changeset revision number
                    %b   basename of the exporting repository
                    %h   short-form changeset hash (12 bytes of hexadecimal)
                    %n   zero-padded sequence number, starting at 1
             options -- string (default: '')
                     -a --text           treat all files as text
                        --switch-parent  diff against the second parent
                    * Without the -a option, export will avoid
                      generating diffs of files it detects as
                      binary. With -a, export will generate a diff
                      anyway, probably with undesirable results.
                    * With the --switch-parent option, the diff will
                      be against the second parent. It can be useful
                      to review a merge.
        """
        if not isinstance(revs, list):
            revs = [int(revs)]
        if not isinstance(filename, str):
            raise TypeError, 'filename must be a string'
        if filename[-6:] != '.patch':
            filename += '.patch'
        options += ' -o "%s"'%(os.path.abspath(filename))
        if filename == '%R.patch':
            print "Output will be written to revision numbered file."%revs
        else:
            print "Output will be written to '%s'"%filename
        if text:
            options += ' -a'
        self('export %s %s'%(options, ' '.join([str(x) for x in revs])))

    def import_patch(self, filename, options=''):
        """
        Import an ordered set of patches from patch file, i.e., a plain
        text file created using the export command.

        If there are outstanding changes in the working directory, import
        will abort unless given the -f flag.

        If imported patch was generated by the export command, user
        and description from patch override values from message
        headers and body.  Values given as options with -m and -u
        override these.

        INPUT:
            filename  -- a string
            options -- a string
                options:  [-p NUM] [-b BASE] [-m MESSAGE] [-f] PATCH...
                 -p --strip    directory strip option for patch. This has the same
                               meaning as the corresponding patch option (default: 1)
                 -m --message  use <text> as commit message
                 -b --base     base path
                 -f --force    skip check for outstanding uncommitted changes
        """
        self('import "%s" %s'%(os.path.abspath(filename),options))

    def add(self, files, options=''):
        """
        Add the given list of files (or file) or directories
        to your HG repository.  They must exist already.

        To see a list of files that haven't been added to the
        repository do self.status().  They will appear with an
        explanation point next them.

        Add needs to be called whenever you add a new file or
        directory to your project.  Of course, it also needs to be
        called when you first create the project, to let hg know
        which files should be kept track of.

        INPUT:
            files -- list or string; name of file or directory.
            options -- string
        """
        if isinstance(files, str):
            files = [files]
        for file in files:
            file = os.path.abspath(file)
            print "Adding file %s"%file
            self('add %s "%s"'%(options, file))

    def remove(self, files, options=''):
        """
        Remove the given list of files (or file) or directories
        from your HG repository.

        INPUT:
            files -- list or string; name of file or directory.
            options -- string (e.g., '-f')
        """
        if isinstance(files, str):
            files = [files]
        for file in files:
            file = os.path.abspath(file)
            print "Removing file %s"%file
            self('rm %s "%s"'%(options, file))

    rm = remove

    def rename(self, src, dest, options=''):
        """
        Move (rename) the given file.

        INPUT:
            src, dest -- strings that define files, relative to self.dir()
            options --
                 -A --after    record a rename that has already occurred
                 -f --force    forcibly copy over an existing managed file
                 -I --include  include names matching the given patterns
                 -X --exclude  exclude names matching the given patterns
                 -n --dry-run  do not perform actions, just print output
        """
        if isinstance(files, str):
            files = [files]
        for file in files:
            file = os.path.abspath(file)
            print "Moving %s --> %s"%file
            self('mv %s "%s"'%(options, file))

    move = rename
    mv = rename

    def log(self, options=''):
        """
        Display the change log for this repository.  This is a list of
        all changesets ordered by revision number.
        """
        self('log %s | %s'%(options, pager()))

    changes = log


    def diff(self, files='', rev=None):
        """
        Show differences between revisions for the specified files as a unified diff.

        By default this command tells you exactly what you have
        changed in your working repository since you last commited
        changes.

        INPUT:
            files -- space separated list of files (relative to self.dir())
            rev -- None or a list of integers.

        Differences between files are shown using the unified diff format.

        When two revision arguments are given, then changes are shown
        between those revisions. If only one revision is specified then
        that revision is compared to the working directory, and, when no
        revisions are specified, the working directory files are compared
        to its parent.
        """
        if not rev is None:
            options = ' '.join(['-r %s'%r for r in rev]) + '  ' + files
        else:
            options = files
        self('diff %s | %s'%(options, pager()))

    what = diff

    def dir(self):
        """
        Return the directory where this repository is located.
        """
        return self.__dir

    def url(self):
        """
        Return the default 'master url' for this repository.
        """
        return self.__url

    def help(self, cmd=''):
        r"""
        Return help about the given command, or if cmd is omitted
        a list of commands.

        If this hg object is called hg_sage, then you
        call a command using
             \code{hg_sage('usual hg command line notation')}
        """
        self('%s --help | %s'%(cmd, pager()))

    def pull(self, url=None, options='-u'):
        """
        Pull all new patches from the repository at the given url,
        or use the default 'official' repository if no url is
        specified.

        INPUT:
            url:  default: self.url() -- the official repository
                   * http://[user@]host[:port]/[path]
                   * https://[user@]host[:port]/[path]
                   * ssh://[user@]host[:port]/[path]
                   * local directory (starting with a /)
                   * name of a branch (for hg_sage); no /'s
            options: (default: '-u')
                 -u --update     update the working directory to tip after pull
                 -e --ssh        specify ssh command to use
                 -f --force      run even when remote repository is unrelated
                 -r --rev        a specific revision you would like to pull
                 --remotecmd  specify hg command to run on the remote side

        Some notes about using SSH with Mercurial:
        - SSH requires an accessible shell account on the destination machine
          and a copy of hg in the remote path or specified with as remotecmd.
        - path is relative to the remote user's home directory by default.
          Use an extra slash at the start of a path to specify an absolute path:
            ssh://example.com//tmp/repository
        - Mercurial doesn't use its own compression via SSH; the right thing
          to do is to configure it in your ~/.ssh/ssh_config, e.g.:
            Host *.mylocalnetwork.example.com
              Compression off
            Host *
              Compression on
          Alternatively specify "ssh -C" as your ssh command in your hgrc or
          with the --ssh command line option.
        """
        if url is None:
            url = self.__url
        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('pull %s %s'%(options, url))
        if self.__target == 'sage':
            print ""
            print "Now building the new SAGE libraries"
            os.system('sage -b')
            print "You *MUST* restart SAGE in order for the changes to take effect!"

        print "If it says use 'hg merge' above, then you should"
        print "type hg_sage.merge(), where hg_sage is the name"
        print "of the repository you are using.  This might not"
        print "work with the notebook yet."

    def merge(self, options='-f'):
        """
        Merge working directory with another revision

        Merge the contents of the current working directory and the
        requested revision. Files that changed between either parent are
        marked as changed for the next commit and a commit must be
        performed before any further updates are allowed.

        INPUT:
            options -- default: '-f'
                 -b --branch  merge with head of a specific branch
                 -f --force   force a merge with outstanding changes
        """
        self('merge %s'%options)

    def update(self, options='-f'):
        """
        update or merge working directory

        Update the working directory to the specified revision.

        If there are no outstanding changes in the working directory and
        there is a linear relationship between the current version and the
        requested version, the result is the requested version.

        To merge the working directory with another revision, use the
        merge command.

        By default, update will refuse to run if doing so would require
        merging or discarding local changes.

        aliases: up, checkout, co

        INPUT:
            options -- string (default: '-f')
             -b --branch  checkout the head of a specific branch
             -C --clean   overwrite locally modified files
             -f --force   force a merge with outstanding changes
        """
        self('merge %s'%options)

    up = update
    checkout = update
    co = update

    def head(self, options=''):
        """
        show current repository heads

        Show all repository head changesets.

        Repository "heads" are changesets that don't have children
        changesets. They are where development generally takes place and
        are the usual targets for update and merge operations.

        INPUT:
            options -- string (default: '')
             -b --branches  show branches
                --style     display using template map file
             -r --rev       show only heads which are descendants of rev
                --template  display with template
        """
        self('head %s'%options)

    heads = head

    def switch(self, name=None):
        r"""
        Switch to a different branch.  You must restart SAGE after switching.

        Only available for \code{hg_sage.}

        INPUT:
            name -- name of a SAGE branch (default: None)

        If the name is not given, this function returns a list of all branches.
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
        Clone the current branch of the SAGE library, and make it active.

        Only available for \code{hg_sage.}

        Use \code{hg_sage.switch('branch_name')} to switch to a different branch.
        You must restart SAGE after switching.

        INPUT:
            name -- string
            rev -- integer or None (default)

        If rev is None, clones the latest recorded version of the repository.
        This is very fast, e.g., about 30-60 seconds (including any build).
        If a specific revision is specified, cloning may take much longer
        (e.g., 5 minutes), since all Pyrex code has to be regenerated and
        compiled.

        EXAMPLES:

        Make a clone of the repository called testing.  A copy of the
        current repository will be created in a directory sage-testing,
        then <SAGE_ROOT>/devel/sage will point to sage-testing, and
        when you next restart SAGE that's the version you'll be using.

            sage.: hg_sage.clone('testing')
            ...

        Make a clone of the repository as it was at revision 1328.
            sage.: hg_sage.clone('testing', 1328)
            ...
        """
        if not self.__cloneable:
            raise RuntimeError, "only available for hg_sage"
        name = '_'.join(str(name).split())
        if rev is None:
            os.system('sage -clone %s'%name)
        else:
            os.system('sage -clone %s -r %s'%(name, int(rev)))

    def commit(self, files='', comment=None, options=''):
        """
        Commit your changes to the repository.

        Quit out of the editor without saving to not record your
        changes.

        INPUT:
             files -- space separated string of file names (optional)
                      If specified only those files are commited.
                      The path must be absolute or relative to
                      self.dir().

             comment -- optional changeset comment.  If you don't give
                      it you will be dumped into an editor.  If you're
                      using the SAGE notebook, you *must* specify a comment.

             options -- string:
                 -A --addremove  mark new/missing files as added/removed before committing
                 -m --message    use <text> as commit message
                 -l --logfile    read the commit message from <file>
                 -d --date       record datecode as commit date
                 -u --user       record user as commiter
                 -I --include    include names matching the given patterns
                 -X --exclude    exclude names matching the given patterns

        \note{If you create new files you should first add them with the add method.}
        """
        if sage.server.support.EMBEDDED_MODE and comment is None:
            raise RuntimeError, "You're using the SAGE notebook, so you *must* explicitly specify the comment in the commit command."
        if comment:
            self('commit %s -m "%s" %s '%(options, comment, files))
        else:
            self('commit %s %s'%(options, files))

    record = commit
    ci = commit

    def rollback(self):
        """
        Remove recorded patches without changing the working copy.
        """
        self('rollback')

    def bundle(self, filename, options='', url=None):
        """
        Create an hg changeset bundle with the given filename against the
        repository at the given url (which is by default the
        'official' SAGE repository).

        This is a file that you should probably send to William Stein
        (wstein@gmail.com), post to a web page, or send to sage-devel.
        It will be written to the current directory.
        """
        if url is None:
            url = self.__url
        # We write to a local tmp file, then move, since unders
        # windows hg has a bug that makes it fail to write
        # to any filename that is at all complicated!
        filename = os.path.abspath(filename) + '.hg'
        print 'Writing to %s'%filename
        tmpfile = '%s/tmphg'%self.__dir
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        self('bundle %s tmphg %s'%(options, url))
        if os.path.exists(tmpfile):
            shutil.move(tmpfile, filename)
            print 'Successfully created hg patch bundle %s'%filename
        else:
            print 'Problem creating hg patch bundle %s'%filename

    send = bundle
    save = send


##############################################################################
# Initialize the actual repositories.
##############################################################################

import misc

SAGE_ROOT = misc.SAGE_ROOT
try:
    SAGE_SERVER = os.environ['SAGE_SERVER'] + '/hg/'
except KeyError:
    print "Falling back to a hard coded sage server in misc/hg.py"
    SAGE_SERVER = "http://sage.math.washington.edu/sage/hg/"

hg_sage    = HG('%s/devel/sage'%SAGE_ROOT,
                'SAGE Library Source Code',
                url='%s/sage-main'%SAGE_SERVER,
                cloneable=True)

hg_doc     = HG('%s/devel/doc'%SAGE_ROOT,
                'SAGE Documentation',
                url='%s/doc-main'%SAGE_SERVER)

hg_scripts = HG('%s/local/bin/'%SAGE_ROOT,
                'SAGE Scripts',
                url='%s/scripts-main'%SAGE_SERVER)

hg_extcode = HG('%s/data/extcode'%SAGE_ROOT,
                'SAGE External System Code (e.g., PARI, MAGMA, etc.)',
                url='%s/extcode-main'%SAGE_SERVER)
