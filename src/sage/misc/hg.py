r"""
SAGE Interface to the HG/Mercurial Revision Control System

These functions make setup and use of source control with SAGE easier, using
the distributed Mercurial HG source control system.  To learn about Mercurial,
see http://www.selenic.com/mercurial/wiki/.

This system should all be fully usable from the SAGE notebook (except
for merging, currently).
This system should all be mostly from the SAGE notebook.

\begin{itemize}
\item Use \code{hg_sage.record()} to record all of your changes.
\item Use \code{hg_sage.bundle('filename')} to bundle them up to send them.
\item Use \code{hg_sage.inspect('filename.hg')} to inspect a bundle.
\item Use \code{hg_sage.unbundle('filename.hg')} to import a bundle into your
      repository.
\item Use \code{hg_sage.pull()} to synchronize with the latest official
      stable SAGE changesets.
\end{itemize}
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
from   remote_file import get_remote_file
from   sage.server.misc import print_open_msg

def embedded():
    return sage.server.support.EMBEDDED_MODE

def pager():
    if embedded():
        return 'cat'
    else:
        return 'less'

class HG:
    r"""
    This is an HG (Mercurial) repository.

    To learn about Mercurial, see http://www.selenic.com/mercurial/wiki/.

    This system should all be fully usable from the SAGE notebook.

    Most commands are directly provided as member functions.  However,
    you can use the full functionality of hg, i.e.,
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
        return "Hg repository '%s' in directory %s"%(self.__name, self.__dir)

    def status(self):
        print("Getting status of modified or unknown files:")
        self('status')
        print "\n---\n"
        if self.__name == "SAGE Library Source Code":
            b = branch_current_hg()
            if b == '': b='main'
	    elif b[-1] == '/':
	        b = b[:-1]
            print("Branch: %s"%b)



    def _changed_files(self):
        out, err = self('status', interactive=False)
        v = [x for x in out.split('\n') if (x.strip()[:1] != '?' and x.strip()[:1] != '!') and len(x) != 0]
        return len(v) > 0

    def _ensure_safe(self):
        """
        Ensure that the repository is in a safe state to have changes
        applied to it, i.e., that all changes to controlled files in
        the working directory are recorded.
        """
        if self._changed_files():
            self.ci()
        if self._changed_files():
            raise RuntimeError, "Refusing to do operation since you still have unrecorded changes. You must check in all changes in your working repository first."

    def _warning(self):
        if not os.path.exists(os.environ['HOME'] + '/.hgrc'):
            print "\nWARNING:"
            print "Make sure to create a ~/.hgrc file:"
            print "-"*70
            print "[ui]"
            print "username = William Stein <wstein@gmail.com>"
            print "-"*70
            print "\n"

    def __call__(self, cmd=None, interactive=True):
        """
        Run 'hg cmd' where cmd is an arbitrary string
        in the hg repository.

        INPUT:
            cmd -- string, the hg command line (everything after 'hg')
            interactive -- If True, runs using os.system, so user can
                           interactively interact with hg, i.e., this
                           is needed when you record changes because
                           the editor pops up.
                           If False, popen3 is used to launch hg
                           as a subprocess.
        OUTPUT:
            * If interactive is True, returns the exit code of the system call.
            * If interactive is False, returns the output and error text.
            * If cmd is not supplied, returns the output of the 'status' command
        """
        self._warning()
        if cmd is None:
            cmd = 'status'
        s = 'cd "%s" && hg %s'%(self.__dir, cmd)
        print s
        if interactive:
            e = os.system(s)
            return e
        else:
            x = os.popen3(s)
            x[0].close()
            out = x[1].read()
            err = x[2].read()
            return out, err

    def serve(self, port=8200, address='localhost',
              open_viewer=True, options=''):
        """
        Start a web server for this repository.

        This server is very nice -- you can browse all files in the
        repository, see their changelogs, see who wrote any given
        line, etc.  Very nice.

        INPUT:
            port -- port that the server will listen on
            address --  (default: 'localhost') address to listen on
            open_viewer -- boolean (default: True); whether to pop up the web page
            options -- a string passed directly to hg's serve command.
        """
        if open_viewer:
            cmd = 'sleep 1; %s http://%s:%s 1>&2 >/dev/null'%(browser(),
                                                              address, port)
            t = tmp_filename()
            open(t,'w').write(cmd)
            P = os.path.abspath(t)
            os.system('chmod +x %s; %s &'%(P, P))

        print_open_msg(address, port)
        self('serve --address %s --port %s  %s'%(address, port, options))
        print_open_msg(address, port)

    browse = serve

    def unbundle(self, bundle, update=True, options=''):
        """
        Apply patches from a hg patch to the repository.

        If the bundle is a .patch file, instead call the import_patch method.
        To see what is in a bundle before applying it, using self.incoming(bundle).

        INPUT:
             bundle -- an hg bundle (created with the bundle command)
             update -- if True (the default), update the working directory after unbundling.
        """
        if bundle.startswith("http://") or bundle.startswith("https://"):
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
        print "this just means you need to do an x.pull(),"
        print "where x is the hg_ object you just called this method on."
        self('unbundle %s "%s"'%(options, bundle))

    apply = unbundle

    def export(self, revs, filename=None, text=False, options=''):
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
             filename -- (default: '%R.patch') The name of the file is given using a format
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
        if filename is None:
            filename = '%R.patch'
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

        ALIASES: patch
        """
        if filename.startswith("http://") or filename.startswith("https://"):
            filename = get_remote_file(filename, verbose=True)
        self._ensure_safe()
        self('import  %s "%s"'%(options, os.path.abspath(filename)))

    patch = import_patch

    def incoming(self, source, options=''):
        """
        Show new changesets found in the given source.  This even
        works if the source is a bundle file (ends in .hg or .bundle).

        Show new changesets found in the specified path/URL or the default
        pull location. These are the changesets that would be pulled if a pull
        was requested.

        For remote repository, using --bundle avoids downloading the changesets
        twice if the incoming is followed by a pull.

        See pull for valid source format details.

        ALIAS: inspect

        INPUT:
            filename -- string
            options -- string '[-p] [-n] [-M] [-r REV] ...'
                         -M --no-merges     do not show merges
                         -f --force         run even when remote repository is unrelated
                            --style         display using template map file
                         -n --newest-first  show newest record first
                            --bundle        file to store the bundles into
                         -p --patch         show patch
                         -r --rev           a specific revision you would like to pull
                            --template      display with template
                         -e --ssh           specify ssh command to use
                            --remotecmd     specify hg command to run on the remote side
        """
        if source.startswith("http://") or source.startswith("https://"):
            source = get_remote_file(source, verbose=True)
        if os.path.exists(source):
            source = os.path.abspath(source)
        if os.path.splitext(source)[1] in ['.hg', '.bundle']:
            source = 'bundle://%s'%source
        self('incoming %s "%s"'%(options, source))

    inspect = incoming


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
            print "Moving %s --> %s"%file
            self('mv %s "%s"'%(options, file))

    move = rename
    mv = rename

    def log(self, branches=None, keyword=None, limit=None,
                  rev=None, merges=False, only_merges=False,
                  patch=None, template=False, include=None,
                  exclude=None, verbose=False):
        """
        Display the change log for this repository.  This is a list of
        changesets ordered by revision number.

        By default this command outputs: changeset id and hash, tags,
        non-trivial parents, user, date and time, and a summary for each
        commit.

        INPUT:
            branches -- (string, default: None) show given branches
            keyword  -- (string, default: None) search for a keyword
            limit    -- (integer, default: None, or 20 in notebook mdoe)
                        limit number of changes displayed
            rev      -- (integer) show the specified revision
            merges   -- (bool, default: False) whether or not to show merges
            only_merges -- (bool, default: False) if true, show only merges
            patch    -- (string, default: None) show given patch
            template -- (string, default: None) display with template
            include  -- (string, default: None) include names matching the given patterns
            exclude  -- (string, default: None) exclude names matching the given patterns
            verbose  -- (bool, default: False) If true, the list of changed
                        files and full commit message is shown.
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
        if include:
            options += '-I "%s"'%include
        if exclude:
            options += '-X "%s"'%exclude
        if verbose:
            options = '-v ' + options

        self('log %s | %s'%(options, pager()))

    changes = log
    history = log

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
            if not isinstance(rev, (list, tuple)):
                rev = [rev]
            options = ' '.join(['-r %s'%r for r in rev]) + '  ' + files
        else:
            options = files
        self('diff %s | %s'%(options, pager()))

    what = diff

    def revert(self, files='', options='', rev=None):
        """
        Revert files or dirs to their states as of some revision

            With no revision specified, revert the named files or
            directories to the contents they had in the parent of the
            working directory.  This restores the contents of the
            affected files to an unmodified state.  If the working
            directory has two parents, you must explicitly specify the
            revision to revert to.

            Modified files are saved with a .orig suffix before
            reverting.  To disable these backups, use --no-backup.

            Using the -r option, revert the given files or directories
            to their contents as of a specific revision.  This can be
            helpful to 'roll back' some or all of a change that should
            not have been committed.

            Revert modifies the working directory.  It does not commit
            any changes, or change the parent of the working
            directory.  If you revert to a revision other than the
            parent of the working directory, the reverted files will
            thus appear modified afterwards.

            If a file has been deleted, it is recreated.  If the executable
            mode of a file was changed, it is reset.

            If names are given, all files matching the names are reverted.

            If no arguments are given, all files in the repository are
            reverted.

        OPTIONS:
            --no-backup  do not save backup copies of files
         -I --include    include names matching given patterns
         -X --exclude    exclude names matching given patterns
         -n --dry-run    do not perform actions, just print output
        """
        if not rev is None:
            options = ' -r %s %s'%(rev, files)
        else:
            options = files
        self('revert %s'%options)

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
        self._ensure_safe()

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

    def merge(self, options=''):
        """
        Merge working directory with another revision

        Merge the contents of the current working directory and the
        requested revision. Files that changed between either parent are
        marked as changed for the next commit and a commit must be
        performed before any further updates are allowed.

        INPUT:
            options -- default: ''
                'tip' -- tip
                 -b --branch  merge with head of a specific branch
                 -f --force   force a merge with outstanding changes
        """
        self('merge %s'%options)

    def update(self, options=''):
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
            options -- string (default: '')
             -b --branch  checkout the head of a specific branch
             -C --clean   overwrite locally modified files
             -f --force   force a merge with outstanding changes
        """
        self('update %s'%options)

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

        Only available for the \code{hg_sage} repository.

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

    def commit(self, files='', comment=None, options='', diff=True):
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

             diff -- (default: True) if True show diffs between your repository
                             and your working repository before recording changes.

        \note{If you create new files you should first add them with the add method.}
        """
        if sage.server.support.EMBEDDED_MODE and comment is None:
            raise RuntimeError, "You're using the SAGE notebook, so you *must* explicitly specify the comment in the commit command."
        if diff:
            self.diff(files)

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

    def bundle(self, filename, options='', url=None, base=None):
        r"""
        Create an hg changeset bundle with the given filename against the
        repository at the given url (which is by default the
        'official' SAGE repository).

        If you have internet access, it's best to just do
        \code{hg_sage.bundle(filename)}.  If you don't
        find a revision r that you and the person unbundling
        both have (by looking at \code{hg_sage.log()}), then
        do \code{hg_sage.bundle(filename, base=r)}.

        Use self.inspect('file.bundle') to inspect the resulting bundle.

        This is a file that you should probably send to William Stein
        (wstein@gmail.com), post to a web page, or send to sage-devel.
        It will be written to the current directory.

        INPUT:
            filename -- output file in which to put bundle
            options -- pass to hg
            url -- url to bundle against (default: SAGE_SERVER)
            base -- a base changeset revision number to bundle
                    against (doesn't require internet access)
        """
        if not base is None:
            url = ''
            options = '--base=%s %s'%(int(base), options)

        if url is None:
            url = self.__url

        # make sure that we don't accidentally create a file ending in '.hg.hg'
        if filename[-3:] == '.hg':
            filename = filename[:-3]
        # We write to a local tmp file, then move, since unders
        # windows hg has a bug that makes it fail to write
        # to any filename that is at all complicated!
        filename = os.path.abspath(filename)
        if filename[-3:] != '.hg':
            filename += '.hg'
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

hg_c_lib = HG('%s/devel/c_lib'%SAGE_ROOT,
                'SAGE C-library code',
                url='%s/extcode-main'%SAGE_SERVER)
