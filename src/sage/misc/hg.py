r"""
Sage Interface to the HG/Mercurial Revision Control System

These functions make setup and use of source control with Sage
easier, using the distributed Mercurial HG source control system.
To learn about Mercurial, see
http://www.selenic.com/mercurial/wiki/ , in particular
UnderstandingMercurial .

This system should all be fully usable from the Sage notebook
(except for merging, currently). This system should all be mostly
usable from the Sage notebook.


-  Use ``hg_sage.record()`` to record all of your
   changes.

-  Use ``hg_sage.bundle('filename')`` to bundle them
   up to send them.

-  Use ``hg_sage.inspect('filename.hg')`` to inspect a
   bundle.

-  Use ``hg_sage.unbundle('filename.hg')`` to import a
   bundle into your repository.

-  Use ``hg_sage.pull()`` to synchronize with the
   latest official stable Sage changesets.
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
import re

sage_trac_re = re.compile('http[s]?://(sagetrac\.org|trac\.sagemath\.org)/sage_trac/attachment/ticket/[0-9]+/.*\.(patch|hg)')

def get_remote_file(f, **kwds):
    """
    Wrap the get_remote_file method to move the file if it ends in
    ?stuff, as happens with funny urls from web servers.
    """
    g = get_remote_file0(f, **kwds)
    i = g.find('?')
    if i >= 0:
        h = g[:i]
        os.rename(g,h)
        return h
    return g

def pager():
    """
    Return a page program, which is either cat or less at present.

    Return cat if embedded in the notebook, and less otherwise.
    """
    if embedded():
        return 'cat'
    else:
        return 'less'


hg_docstring = r"""
This is an HG (Mercurial) repository.

To learn about Mercurial, see http://www.selenic.com/mercurial/wiki/.

This system should all be fully usable from the Sage notebook.

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
        return "Hg repository '%s' in directory %s"%(self.__name, self.__dir)


    def current_branch(self, print_flag=True):
        """
        Lists the current branch.
        """
        branch_name = branch_current_hg()
        if print_flag:
            print "The current branch is: " + branch_name
        else:
            return branch_name

    def list_branches(self, print_flag=True):
        """
        Print all branches in the current Sage installation.
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
        applied to it, i.e., that all changes to controlled files in the
        working directory are recorded.
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
        Run 'hg cmd' where cmd is an arbitrary string in the hg
        repository.

        INPUT:


        -  ``cmd`` - string, the hg command line (everything
           after 'hg')

        -  ``interactive`` - If True, runs using os.system, so
           user can interactively interact with hg, i.e., this is needed when
           you record changes because the editor pops up. If False, popen3 is
           used to launch hg as a subprocess.


        OUTPUT:

        - If interactive is True, returns the exit code of the
          system call.

        - If interactive is False, returns the output and
          error text.

        - If cmd is not supplied, returns the output of the
          'status' command
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

        This server is very nice - you can browse all files in the
        repository, see their changelogs, see who wrote any given line,
        etc. Very nice.

        INPUT:


        -  ``port`` - port that the server will listen on

        -  ``address`` - (default: 'localhost') address to
           listen on

        -  ``open_viewer`` - boolean (default: True); whether
           to pop up the web page

        -  ``options`` - a string passed directly to hg's serve
           command.
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

        If the bundle is a .patch file, instead call the import_patch
        method. To see what is in a bundle before applying it, using
        self.incoming(bundle).

        INPUT:


        -  ``bundle`` - an hg bundle (created with the bundle
           command)

        -  ``update`` - if True (the default), update the
           working directory after unbundling.
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
        self('unbundle %s "%s"'%(options, bundle))

    apply = unbundle

    def export(self, revs, filename=None, text=False, options=''):
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

        -  ``options`` - string (default: '')

           - ``'-a --text'`` - treat all files as text

           - ``'--switch-parent'`` -  diff against the second parent

           Without the ``-a`` option, export will avoid generating
           diffs of files it detects as binary. With ``-a``, export
           will generate a diff anyway, probably with undesirable
           results.

           With the ``--switch-parent`` option, the diff will be
           against the second parent. It can be useful to review a
           merge.
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

        If imported patch was generated by the export command, user and
        description from patch override values from message headers and
        body. Values given as options with -m and -u override these.

        INPUT:


        -  ``filename`` - string

        -  ``options`` - string (default: '')::

               options: [-p NUM] [-b BASE] [-m MESSage] [-f] PATCH...
                 -p --strip         directory strip option for patch. This has the same meaning as the corresponding patch option (default: 1)
                 -m --message       use text as commit message
                 -b --base          base path
                 -f --force         skip check for outstanding uncommitted changes

        ALIASES: patch
        """
        if filename.startswith("http://") or filename.startswith("https://"):
            filename = get_remote_file(filename, verbose=True)
        self._ensure_safe()
        self('import  %s "%s"'%(options, os.path.abspath(filename)))

    patch = import_patch

    def incoming(self, source, options='-p'):
        """
        Show new changesets found in the given source and display the
        corresponding diffs. This even works if the source is a bundle file
        (ends in .hg or .bundle). This is great because it lets you "see
        inside" the myserious binary-only .hg files.

        Show new changesets found in the specified path/URL or the default
        pull location. These are the changesets that would be pulled if a
        pull was requested.

        For remote repository, using -bundle avoids downloading the
        changesets twice if the incoming is followed by a pull.

        See pull for valid source format details.

        ALIAS: inspect

        INPUT:


        -  ``filename`` - string

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
        """
        if source.startswith("http://") or source.startswith("https://"):
            source = get_remote_file(source, verbose=True)
        if os.path.exists(source):
            source = os.path.abspath(source)
        if os.path.splitext(source)[1] in ['.hg', '.bundle']:
            source = 'bundle://%s'%source
        self('incoming %s "%s" | %s'%(options, source, pager()))

    inspect = incoming


    def add(self, files, options=''):
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

        -  ``options`` - string (e.g., '-dry-run')
        """
        if isinstance(files, str):
            if ' ' in files:
                files = files.split()
            else:
                files = [files]
        for file in files:
            print "Adding file %s"%file
            self('add %s "%s"'%(options, file))

    def remove(self, files, options=''):
        """
        Remove the given list of files (or file) or directories from your
        HG repository.

        INPUT:


        -  ``files`` - list or string; name of file or
           directory.

        -  ``options`` - string (e.g., '-f')
        """
        if isinstance(files, str):
            files = [files]
        for file in files:
            print "Removing file %s"%file
            self('rm %s "%s"'%(options, file))

    rm = remove

    def rename(self, src, dest, options=''):
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

        """
        print "Moving %s --> %s"%(src,dest)
        self('mv %s "%s" "%s"'%(options, src,dest))


    move = rename
    mv = rename

    def log(self, branches=None, keyword=None, limit=None,
                  rev=None, merges=True, only_merges=False,
                  patch=None, template=False, include=None,
                  exclude=None, verbose=False):
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

        self('log %s | %s'%(options, pager()))

    changes = log
    history = log

    def diff(self, files='', rev=None):
        """
        Show differences between revisions for the specified files as a
        unified diff.

        By default this command tells you exactly what you have changed in
        your working repository since you last commited changes.

        INPUT:


        -  ``files`` - space separated list of files (relative
           to self.dir())

        -  ``rev`` - None or a list of integers.


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

        With no revision specified, revert the named files or directories
        to the contents they had in the parent of the working directory.
        This restores the contents of the affected files to an unmodified
        state. If the working directory has two parents, you must
        explicitly specify the revision to revert to.

        Modified files are saved with a .orig suffix before reverting. To
        disable these backups, use -no-backup.

        Using the -r option, revert the given files or directories to their
        contents as of a specific revision. This can be helpful to 'roll
        back' some or all of a change that should not have been committed.

        Revert modifies the working directory. It does not commit any
        changes, or change the parent of the working directory. If you
        revert to a revision other than the parent of the working
        directory, the reverted files will thus appear modified
        afterwards.

        If a file has been deleted, it is recreated. If the executable mode
        of a file was changed, it is reset.

        If names are given, all files matching the names are reverted.

        If no arguments are given, all files in the repository are
        reverted.

        OPTIONS::

            --no-backup  do not save backup copies of files
         -I --include    include names matching given patterns
         -X --exclude    exclude names matching given patterns
         -n --dry-run    do not perform actions, just print output
        """
        if not rev is None:
            options = options +' -r %s %s'%(rev, files)
        else:
            options = options + files
        self('revert %s'%options)

    def dir(self):
        """
        Return the directory where this repository is located.
        """
        return self.__dir

    def pull_url(self):
        """
        Return the default 'master url' for this repository.
        """
        return self.__pull_url

    def push_url(self):
        """
        Return the default url for uploading this repository.
        """
        return self.__push_url


    def help(self, cmd=''):
        r"""
        Return help about the given command, or if cmd is omitted a list of
        commands.

        If this hg object is called hg_sage, then you call a command using
        ``hg_sage('usual hg command line notation')``
        """
        self('%s --help | %s'%(cmd, pager()))

    def outgoing(self, url=None, opts=''):
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

           - http://[user@]host[:port]/[path]

           - https://[user@]host[:port]/[path]

           - ssh://[user@]host[:port]/[path]

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

        """
        if url is None:
            url = self.__push_url

        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('outgoing %s %s | %s' % (opts, url, pager()))

    def pull(self, url=None, options='-u'):
        """
        Pull all new patches from the repository at the given url, or use
        the default 'official' repository if no url is specified.

        INPUT:

        -  ``url`` - (Default: self.push_url())  the official
           repository

           - http://[user@]host[:port]/[path]

           - https://[user@]host[:port]/[path]

           - ssh://[user@]host[:port]/[path]

           - local directory (starting with a /)

           - name of a branch (for hg_sage); no /'s

        - ``options`` - (Default: '-u')::

              -u --update     update the working directory to tip after pull
              -e --ssh        specify ssh command to use
              -f --force      run even when remote repository is unrelated
              -r --rev        a specific revision you would like to pull
              --remotecmd     specify hg command to run on the remote side


        Some notes about using SSH with Mercurial:

        - SSH requires an accessible shell account on the destination
          machine and a copy of hg in the remote path or specified
          with as remotecmd.

        - path is relative to the remote user's home directory by
          default. Use an extra slash at the start of a path to
          specify an absolute path: ssh://example.com//tmp/repository

        - Mercurial doesn't use its own compression via SSH; the right
          thing to do is to configure it in your /.ssh/ssh_config,
          e.g.::

              Host *.mylocalnetwork.example.com
                Compression off
              Host *
                Compression on

          Alternatively specify 'ssh -C' as your ssh command in your
          hgrc or with the -ssh command line option.
        """
        self._ensure_safe()

        if url is None:
            url = self.__pull_url
        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('pull %s %s'%(options, url))
        if self.__target == 'sage':
            print ""
            print "Now building the new SAGE libraries"
            os.system('sage -b')
            print "You *MUST* restart SAGE in order for the changes to take effect!"

        print "If it says use 'hg merge' above, then you should"
        print "type hg_%s.merge()."%self.__obj_name

    def push(self, url=None, options=''):
        """
        Push all new patches from the repository to the given destination.

        INPUT:

        -  ``url`` - (Default: self.push_url())  the official
           repository

           - http://[user@]host[:port]/[path]

           - https://[user@]host[:port]/[path]

           - ssh://[user@]host[:port]/[path]

           - local directory (starting with a /)

           - name of a branch (for hg_sage); no /'s

        - ``options`` - (Default: '-u')::

              -e --ssh        specify ssh command to use
              -f --force      run even when remote repository is unrelated
              -r --rev        a specific revision you would like to pull
              --remotecmd     specify hg command to run on the remote side


        Some notes about using SSH with Mercurial:

        - SSH requires an accessible shell account on the destination
          machine and a copy of hg in the remote path or specified
          with as remotecmd.

        - path is relative to the remote user's home directory by
          default. Use an extra slash at the start of a path to
          specify an absolute path: ssh://example.com//tmp/repository

        - Mercurial doesn't use its own compression via SSH; the right
          thing to do is to configure it in your /.ssh/ssh_config,
          e.g.::

              Host *.mylocalnetwork.example.com
                Compression off
              Host *
                Compression on

          Alternatively specify 'ssh -C' as your ssh command in your
          hgrc or with the -ssh command line option.

        """
        self._ensure_safe()

        if url is None:
            url = self.__push_url
        if not '/' in url:
            url = '%s/devel/sage-%s'%(SAGE_ROOT, url)

        self('push %s %s'%(options, url))


    def merge(self, options=''):
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

        """
        self('merge %s'%options)

    def update(self, options=''):
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


        -  ``options`` - string (default: '')::

             -r --rev       show only heads which are descendants of rev
                --style     display using template map file
                --template  display with template

        """
        self('head %s'%options)

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

    def commit(self, files='', comment=None, options='', diff=True):
        r"""
        Commit your changes to the repository.

        Quit out of the editor without saving to not record your changes.

        INPUT:

        - ``files`` - space separated string of file names (optional)
          If specified only those files are commited. The path must be
          absolute or relative to self.dir().

        - ``comment`` - optional changeset comment. If you don't give
           it you will be dumped into an editor. If you're using the
           Sage notebook, you *must* specify a comment.

        - ``options`` - string::

              -A --addremove  mark new/missing files as added/removed before committing
              -m --message    use <text> as commit message
              -l --logfile    read the commit message from <file>
              -d --date       record datecode as commit date
              -u --user       record user as commiter
              -I --include    include names matching the given patterns
              -X --exclude    exclude names matching the given patterns

        - ``diff`` - (default: True) if True show diffs between your repository
          and your working repository before recording changes.

        .. note::

           If you create new files you should first add them with the
           add method.
        """
        if embedded() and comment is None:
            raise RuntimeError, "You're using the SAGE notebook, so you *must* explicitly specify the comment in the commit command."
        if diff:
            self.diff(files)

        if isinstance(files, (list, tuple)):
            files = ' '.join([str(x) for x in files])

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

    def bundle(self, filename, options='', url=None, base=None, to=None):
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
        """
        if not base is None:
            url = ''
            options = '--base=%s %s'%(int(base), options)

        if url is None:
            url = self.__push_url

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
            if not to is None:
                os.system('scp "%s" %s'%(filename, to))
        else:
            print 'Problem creating hg patch bundle %s'%filename

    send = bundle
    save = send


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
