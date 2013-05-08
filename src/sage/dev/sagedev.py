"""
Sagedev
"""
import atexit
import collections
import ConfigParser as configparser
import email.utils
import os
import re
import shutil
import tempfile
import time

from datetime import datetime
from subprocess import call, check_call

from git_interface import GitInterface
from trac_interface import TracInterface
from user_interface import CmdLineInterface

from sage.env import DOT_SAGE
from sage.doctest import DOCTEST_MODE

# regular expressions to parse mercurial patches
HG_HEADER_REGEX = re.compile(r"^# HG changeset patch$")
HG_USER_REGEX = re.compile(r"^# User (.*)$")
HG_DATE_REGEX = re.compile(r"^# Date (\d+) (-?\d+)$")
HG_NODE_REGEX = re.compile(r"^# Node ID ([0-9a-f]+)$")
HG_PARENT_REGEX = re.compile(r"^# Parent +([0-9a-f]+)$")
HG_DIFF_REGEX = re.compile(r"^diff (?:-r [0-9a-f]+ ){1,2}(.*)$")
PM_DIFF_REGEX = re.compile(r"^(?:(?:\+\+\+)|(?:---)) [ab]/([^ ]*)(?: .*)?$")

# regular expressions to parse git patches -- at least those created by us
GIT_FROM_REGEX = re.compile(r"^From: (.*)$")
GIT_SUBJECT_REGEX = re.compile(r"^Subject: (.*)$")
GIT_DATE_REGEX = re.compile(r"^Date: (.*)$")
GIT_DIFF_REGEX = re.compile(r"^diff --git a/(.*) b/(.*)$") # this regex should work for our patches since we do not have spaces in file names

# regular expressions to determine whether a path was written for the new git
# repository of for the old hg repository
HG_PATH_REGEX = re.compile(r"^(?=sage/)|(?=doc/)|(?=module_list\.py)|(?=setup\.py)|(?=c_lib/)")
GIT_PATH_REGEX = re.compile(r"^(?=src/)")

class Config(collections.MutableMapping):
    """
    Wrapper around the ``devrc`` file storing the configuration for
    :class:`SageDev`.

    INPUT:

    - ``devrc`` -- a string (default: the absolute path of the ``devrc`` file in ``DOT_SAGE``)

    EXAMPLES::

        sage: dev._config
        Config('''
        [git]
        dot_git = ...
        doctest = ...
        [trac]
        username = doctest
        ''')

    """
    def __init__(self, devrc = os.path.join(DOT_SAGE, 'devrc')):
        """
        Initialization.
        """
        self._config = configparser.ConfigParser()
        self._devrc = devrc
        self._read_config()

    def __repr__(self):
        """
        Return a printable representation of this element.

        EXAMPLES::

            sage: repr(dev._config)
            "Config('''\n[git]\ndot_git = ...\ndoctest = ...\n[trac]\nusername = doctest\n''')"
        """
        return "Config('''\n"+"\n".join([ "[%s]\n"%s+"\n".join(["%s = %s"%(o,self[s][o]) for o in self[s] ]) for s in self ])+"\n''')"

    def _read_config(self):
        """
        Read the configuration from disk.

        EXAMPLES::

            sage: from sage.dev.sagedev import Config, doctest_config
            sage: c = doctest_config()
            sage: c._write_config()
            sage: c = Config(c._devrc)
            sage: c._read_config()
            sage: c
            Config('''
            [git]
            dot_git = ...
            doctest = ...
            [trac]
            username = doctest
            ''')
        """
        if os.path.exists(self._devrc):
            self._config.read(self._devrc)

    def _write_config(self):
        """
        Write the configuration to disk.

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: os.path.exists(c._devrc)
            False
            sage: c._write_config()
            sage: os.path.exists(c._devrc)
            True
        """
        with open(self._devrc, 'w') as F:
            self._config.write(F)
        # set the configuration file to read only by this user,
        # because it may contain the trac password
        os.chmod(self._devrc, 0600)

    def __getitem__(self, section):
        """
        Return the configurations in ``section``.

        EXAMPLES::

            sage: dev._config['git']
            IndexableForSection('''
            dot_git = ...
            doctest = ...
            ''')
            sage: dev._config['tig']
            Traceback (most recent call last):
            ...
            KeyError: 'tig'
        """
        if not section in self:
            raise KeyError(section)

        class IndexableForSection(collections.MutableMapping):
            def __init__(this, section):
                this._section = section
            def __repr__(this):
                return "IndexableForSection('''\n"+"\n".join(["%s = %s"%(o,this[o]) for o in this])+"\n''')"
            def __getitem__(this, option):
                try:
                    return self._config.get(this._section, option)
                except configparser.NoOptionError:
                    raise KeyError(option)
            def __iter__(this):
                return iter(self._config.options(this._section))
            def __setitem__(this, option, value):
                self._config.set(this._section, option, value)
            def getboolean(this, option):
                return self._config.getboolean(this._section, option)
            def _write_config(this):
                self._write_config()
            def __delitem__(this, option):
                self._config.remove_option(this._section, option)
            def __len__(this):
                return len(self._config.options(this._section))
            def __contains__(this, option):
                return option in self._config.options(this._section)

        return IndexableForSection(section)

    def __contains__(self, section):
        """
        returns true if section is in the configuration

        EXAMPLES::

            sage: 'git' in dev._config
            True
            sage: 'notgit' in dev._config
            False
        """
        return section in self._config.sections()

    def __iter__(self):
        """
        Return an iterator over the section names.

        EXAMPLES::

            sage: list(dev._config)
            ['git', 'trac']

        """
        return iter(self._config.sections())

    def __setitem__(self, section, dictionary):
        """
        Set ``section`` to ``dictionary``.

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: c['foo'] = {'foo':'foo'}
            sage: c['foo']['foo']
            'foo'

        """
        if self._config.has_section(section):
            self.remove_section(section)
        self._config.add_section(section)
        for option, value in dictionary.iteritems():
            self._config.set(section, option, value)

    def __len__(self):
        """
        get the number of sections in the configuration

        EXAMPLES::

            sage: len(dev._config)
            2
        """
        return len(self._config.sections())

    def __delitem__(self, section):
        """
        remove ``section`` from the configuration

        EXAMPLES::

            sage: from sage.dev.sagedev import doctest_config
            sage: c = doctest_config()
            sage: del c['git']
            sage: list(c)
            ['trac']
        """
        self._config.remove_section(section)

class SageDev(object):
    """
    The developer interface for sage.

    This class facilitates access to git and trac.

    EXAMPLES::

        sage: type(dev)
        <class 'sage.dev.sagedev.SageDev'>

    """
    def __init__(self, config = Config()):
        """
        Initialization.

        TESTS::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: type(SageDev(doctest_config()))
            <class 'sage.dev.sagedev.SageDev'>

        """
        self._config = config
        self._UI = CmdLineInterface()

        self.__git = None
        self.__trac = None

    ##
    ## Public interface
    ##

    def switch_ticket(self, ticket, branchname=None, offline=False):
        """
        Switch to a branch associated to ``ticket``.

        INPUT:

        - ``ticket`` -- an integer or a string. If a string, then
          switch to the branch ``ticket`` (``branchname`` must be
          ``None`` in this case).  If an integer, it must be the
          number of a ticket. If ``branchname`` is ``None``, then a
          new name for a ``branch`` is chosen. If ``branchname`` is
          the name of a branch that exists locally, then associate
          this branch to ``ticket``. Otherwise, switch to a new branch
          ``branchname``.

          If ``offline`` is ``False`` and a local branch was already
          associated to this ticket, then :meth:`remote_status` will
          be called on the ticket.  If no local branch was associated
          to the ticket, then the current version of that ticket will
          be downloaded from trac.

        - ``branchname`` -- a string or ``None`` (default: ``None``)

        - ``offline`` -- a boolean (default: ``False``)

        .. SEEALSO::

        - :meth:`download` -- Download the ticket from the remote
          server and merge in the changes.

        - :meth:`create_ticket` -- Make a new ticket.

        - :meth:`vanilla` -- Switch to a released version of Sage.
        """
        if isinstance(ticket, basestring):
            if branchname is not None:
                raise ValueError("Cannot specify two branch names")
        else:
            ticket = int(ticket)
        branch = self.git._ticket_to_branch(ticket)
        if branch is None:
            # ticket does not exist locally
            if isinstance(ticket, basestring):
                raise ValueError("Branch does not exist")
            elif not offline:
                if branchname is None:
                    branchname = "t/%s"%(ticket)
                # No branch associated to that ticket number
                tracbranch = self._trac_branch(ticket)
                if tracbranch is None:
                    # There is not yet a branch on trac for this ticket.
                    ref = "master"
                else:
                    ref = self._fetch(tracbranch)
                self.git.create_branch(branchname, ref)
                self.git._branch[ticket] = branchname
            else:
                raise ValueError("You cannot download a ticket while offline")
        else:
            if branchname is None:
                branchname = branch
            else:
                # User specified a branchname but there's already a branch for that ticket
                # ticket must be an int
                self.git._branch[ticket] = branchname
            self.git.switch_branch(branchname)
        if branch is not None and not offline:
            self.remote_status(branchname)

    def create_ticket(self, branchname=None, base="master", remote_branch=True):
        """
        Create a new ticket on trac.

        INPUT:

        - ``branchname`` -- a string or ``None`` (default: ``None``),
          the name of the local branch that will used for the new
          ticket; if ``None``, a name will be chosen automatically.

        - ``base`` -- a string or ``None`` (default: ``"master"``), a
          branch on which to base the ticket.  If ``None`` then the
          current ticket is used.  If the base is not ``"master"``
          then a dependency will be added.

        - ``remote_branch`` -- a string or boolean (default:
          ``True``), if a string, the name of the remote branch this
          branch should be tracking.  If ``True``, the remote branch
          is determined from the branchname.  If ``False`` no remote
          branch is stored.

        .. SEEALSO::

        - :meth:`switch_ticket` -- Switch to an already existing
          ticket.

        - :meth:`download` -- Download changes from an existing
          ticket.
        """
        ticketnum = self.trac.create_ticket_interactive()
        if ticketnum is None:
            # They didn't succeed.
            return
        if branchname is None:
            branchname = 't/%s'%(ticketnum)
        if base is None:
            base = self.git.current_branch()
            if base is None:
                raise ValueError("You cannot add a detached head as a dependency")
        self.git.create_branch(branchname, base, remote_branch)
        self._ticket[branchname] = ticketnum
        if base != "master":
            self.git._dependencies[branchname] = [base]
        self.git.switch_branch(branchname)

    def commit(self, message=None, interactive=False):
        """
        Create a commit from the pending changes on the current branch.

        INPUT:

        - ``message`` -- a string or ``None`` (default: ``None``), the message of
          the commit; if ``None``, prompt for a message.

        - ``interactive`` -- a boolean (default: ``False``), interactively select
          which part of the changes should be part of the commit. Prompts for the
          addition of untracked files even if ``interactive`` is ``False``.

        .. SEEALSO::

        - :meth:`upload` -- Upload changes to the remote server.  This
          is the next step once you've committed some changes.

        - :meth:`diff` -- Show changes that will be committed.
        """
        curticket = self.current_ticket()
        if curticket is None:
            curbranch = self.git.current_branch()
            if curbranch is None:
                raise ValueError("You may not commit in detached head mode.  Try switching to a branch.")
            prompt = "Are you sure you want to save your changes to branch %s"%(curbranch)
        else:
            prompt = "Are you sure you want to save your changes to ticket #%s"%(curticket)
        if self._UI.confirm(prompt):
            unknown_files = self.git.unknown_files()
            for file in unknown_files:
                if self._UI.confirm("Would you like to commit %s"%(file), default_yes=False):
                    self.git.add(file)
            kwds = {}
            if interactive:
                kwds['interactive'] = True
            else:
                kwds['all'] = True
            if message is None:
                message = self._UI.get_input("Please enter a commit message: ")
            kwds['m'] = message
            self.git.commit(**kwds)
        else:
            self._UI.show("If you want to commit these changes to another ticket use the switch_ticket() method")

    def upload(self, ticket=None, remote_branch=None, force=False, repository=None):
        """
        Upload the current branch to the sage repository.

        INPUT:

        - ``ticket`` -- an integer or ``None`` (default: ``None``), if an integer
          or if this branch is associated to a ticket, set the trac ticket to point
          to this branch.

        - ``remote_branch`` -- a string or ``None`` (default: ``None``), the remote
          branch to upload to; if ``None``, then a default is chosen (XXX: how?)

        - ``force`` -- a boolean (default: ``False``), whether to upload if this is
          not a fast-forward.

        .. SEEALSO::

        - :meth:`commit` -- Save changes to the local repository.

        - :meth:`download` -- Update a ticket with changes from the remote repository.
        """
        branch = self.git.current_branch()
        if branch is None: raise ValueError("You cannot upload a detached head.  See switch_ticket")
        oldticket = self.git._ticket[branch]
        if oldticket is None and ticket is None and remote_branch is None:
            self._UI.show("You don't have a ticket for this branch (%s)"%branch)
            return
        elif ticket is None:
            ticket = oldticket
        elif oldticket != ticket:
            if not self._UI.confirm("Are you sure you want to upload your changes to ticket %s instead of %s?"%(self._print(ticket), self._print(oldticket)), False):
                return
            self.git._ticket[branch] = ticket
        if ticket:
            ref = self._fetch(ticket)
            if not self.git.is_ancestor_of(ref, branch) and not force:
                if not self._UI.confirm("Changes not compatible with remote branch; consider downloading first.  Are you sure you want to continue?", False):
                    return
        remote_branch = remote_branch or self.git._local_to_remote_name(branch)
        if repository is None:
            repository = self.git._repo
        ref = self._fetch(remote_branch, repository=repository)
        if force or self.git.is_ancestor_of(ref, branch):
            self.git.push(repository, "%s:%s" % (branch, remote_branch))
        else:
            raise ValueError("The remote branch has changed; upload failed.  Consider downloading the changes.")
        if ticket:
            commit_id = self.git.branch_exists(branch)
            self.trac._set_branch(ticket, remote_branch, commit_id)
            git_deps = self._dependencies_as_tickets(branch)
            self.trac.set_dependencies(ticket, git_deps)

    def download(self, ticket=None, branchname=None, force=False, repository=None):
        """
        Download the changes made to a remote branch into a given
        ticket or the current branch.

        INPUT:

        - ``ticket`` -- an integer or ``None`` (default: ``None``).

          If an integer and there is a local branch corresponding to
          that ticket, switch to it.  Then merge the branch associated
          to the trac ticket ``ticket`` into that branch.

          If ``ticket`` is ``None`` and this branch is associated to a
          ticket and is not following a non-user remote branch, then
          also merge in the trac ticket branch. If this branch is
          following a non-user remote branch, then merge that branch
          instead.

        - ``branchname`` -- a string or ``None``, only used if there
          is no local branch already associated to ``ticket``.

        - ``force`` -- a boolean (default: ``False``), if ``False``,
          try to merge the remote branch into this branch; if
          ``False``, do not merge, but make this branch equal to the
          remote branch.

        .. SEEALSO::

        - :meth:`merge` -- Merge in local branches.

        - :meth:`upload` -- Upload changes to the remote server.

        - :meth:`switch_ticket` -- Switch to another ticket without
          updating.

        - :meth:`vanilla` -- Switch to a plain release (which is not a
          branch).

        - :meth:`import_patch` -- Import a patch into the current
          ticket.
        """
        if ticket is None:
            branch = self.git.current_branch()
            if branch is None:
                raise ValueError("Cannot download in detached head")
            ticket = self.git._ticket[branch]
        else:
            ticket = int(ticket)
            branch = self.git._ticket_to_branch(ticket)
            if branch is not None:
                self.git.switch_branch(branch)
        if branch is None:
            if branchname is None:
                branch = 't/%s'%(ticket)
            else:
                branch = branchname
            remote_branch = self._trac_branch(ticket)
        else:
            remote_branch = self._remote_pull_branch(branch)
        if remote_branch is None:
            raise ValueError("No remote branch associated to current branch")
        ref = self._fetch(remote_branch, repository=repository)
        if force:
            self.git.branch(branch, ref, f=True)
            overwrite_deps = True
        else:
            overwrite_deps = self.git.is_ancestor_of(branch, ref)
            self.merge(ref, create_dependency=False, download=False)
        if ticket is not None:
            trac_deps = self.trac.dependencies(ticket)
            git_deps = self._dependencies_as_tickets(branch)
            if overwrite_deps:
                self.git._dependencies[ticket] = trac_deps
            else:
                deps = trac_deps
                for d in git_deps:
                    if d not in deps:
                        deps.append(d)
                self.git._dependencies[ticket] = tuple(deps)

    def remote_status(self, ticket=None, quiet=False):
        """
        Show the remote status of ``ticket``.

        For tickets and remote branches, this shows the commit log of the branch on
        the trac ticket a summary of their difference to your related branches, and
        an overview of patchbot results (where applicable).

        INPUT:

        - ``ticket`` -- None, an integer, a string, or the special string "all"

        .. SEEALSO::

        - :meth:`local_tickets` -- Just shows local tickets without
          comparing them to the remote server.

        - :meth:`diff` -- Shows the actual differences on a given
          ticket.

        - :meth:`download` -- Merges in the changes on a given ticket
          from the remote server.

        - :meth:`upload` -- Pushes the changes on a given ticket to
          the remote server.
        """
        if ticket == "all":
            results = []
            for ticket, branch in self.local_tickets:
                results.add((ticket, branch, remote_status(ticket or branch, quiet=quiet)))
            if quiet:
                return results
        remote_branch = self._remote_pull_branch(ticket)
        if remote_branch is None:
            print ticket or "     ", branch, "not tracked remotely"
            return
        remote_ref = self._fetch(remote_branch)
        if isinstance(ticket, int):
            branch = self.git._branch[ticket]
        else:
            branch = ticket
        ahead, behind = self.git.read_output("rev-list", "--left-right", "%s..%s"%(branch, remote_ref), count=True).split()
        behind = int(behind)
        ahead = int(ahead)
        if quiet:
            return ahead, behind
        else:
            print ticket or "     ", branch, "ahead", ahead, "behind", behind


    def import_patch(self, patchname=None, url=None, local_file=None, diff_format=None, header_format=None, path_format=None):
        """
        Import a patch into the branch for the current ticket.

        If ``local_file`` is specified, apply the file it points to.

        Otherwise, apply the patch using :meth:`download_patch` and apply it.

        INPUT:

        - ``patchname`` -- a string or ``None`` (default: ``None``)

        - ``url`` -- a string or ``None`` (default: ``None``)

        - ``local_file`` -- a string or ``None`` (default: ``None``)

        .. SEEALSO::

        - :meth:`download_patch` -- This function downloads a patch to
          a local file.

        - :meth:`download` -- This function is used to merge in
          changes from a git branch rather than a patch.
        """
        ticketnum = self.current_ticket()
        if not self.git.reset_to_clean_state(): return
        if not self.git.reset_to_clean_working_directory(): return

        if not local_file:
            return self.import_patch(local_file = self.download_patch(ticketnum = ticketnum, patchname = patchname, url = url), diff_format=diff_format, header_format=header_format, path_format=path_format)
        else:
            if patchname or url:
                raise ValueError("If `local_file` is specified, `patchname`, and `url` must not be specified.")
            lines = open(local_file).read().splitlines()
            lines = self._rewrite_patch(lines, to_header_format="git", to_path_format="new", from_diff_format=diff_format, from_header_format=header_format, from_path_format=path_format)
            outfile = os.path.join(self.tmp_dir, "patch_new")
            open(outfile, 'w').writelines("\n".join(lines)+"\n")
            self._UI.show("Trying to apply reformatted patch `%s` ..."%outfile)
            shared_args = ["--ignore-whitespace",outfile]
            am_args = shared_args+["--resolvemsg=''"]
            am = self.git.am(*am_args)
            if am: # apply failed
                if not self._UI.confirm("The patch does not apply cleanly. Would you like to apply it anyway and create reject files for the parts that do not apply?", default_yes=False):
                    self._UI.show("Not applying patch.")
                    self.git.reset_to_clean_state(interactive=False)
                    return

                apply_args = shared_args + ["--reject"]
                apply = self.git.apply(*apply_args)
                if apply: # apply failed
                    if self._UI.get_input("The patch did not apply cleanly. Please integrate the `.rej` files that were created and resolve conflicts. When you did, type `resolved`. If you want to abort this process, type `abort`.",["resolved","abort"]) == "abort":
                        self.git.reset_to_clean_state(interactive=False)
                        self.git.reset_to_clean_working_directory(interactive=False)
                        return
                else:
                    self._UI.show("It seemed that the patch would not apply, but in fact it did.")

                self.git.add("--update")
                self.git.am("--resolved")

    def download_patch(self, ticketnum=None, patchname=None, url=None):
        """
        Download a patch to a temporary directory.

        If only ``ticketnum`` is specified and the ticket has only one
        attachment, download the patch attached to ``ticketnum``.

        If ``ticketnum`` and ``patchname`` are specified, download the
        patch ``patchname`` attached to ``ticketnum``.

        If ``url`` is specified, download ``url``.

        Raise an error on any other combination of parameters.

        INPUT:

        - ``ticketnum`` -- an int or an Integer or ``None`` (default:
          ``None``)

        - ``patchname`` -- a string or ``None`` (default: ``None``)

        - ``url`` -- a string or ``None`` (default: ``None``)

        OUTPUT:

        Returns the absolute file name of the returned file.

        .. SEEALSO::

        - :meth:`import_patch` -- also creates a commit on the current
          branch from the patch.
        """
        if url:
            if ticketnum or patchname:
                raise ValueError("If `url` is specifed, `ticketnum` and `patchname` must not be specified.")
            ret = os.path.join(self.tmp_dir,"patch")
            check_call(["wget","-r","-O",ret,url])
            return ret
        elif ticketnum:
            if patchname:
                return self.download_patch(url = self._config['trac']['server']+"raw-attachment/ticket/%s/%s"%(ticketnum,patchname))
            else:
                attachments = self.trac.attachment_names(ticketnum)
                if len(attachments) == 0:
                    raise ValueError("Ticket #%s has no attachments."%ticketnum)
                if len(attachments) == 1:
                    return self.download_patch(ticketnum = ticketnum, patchname = attachments[0])
                else:
                    raise ValueError("Ticket #%s has more than one attachment but parameter `patchname` is not present."%ticketnum)
        else:
            raise ValueError("If `url` is not specified, `ticketnum` must be specified")

    def diff(self, base="commit"):
        """
        Show how the current file system differs from ``base``.

        INPUT:

        - ``base`` -- show the differences against the latest
          ``'commit'`` (the default), against the branch ``'master'``
          (or any other branch name), or the merge of the
          ``'dependencies'`` of the current ticket (if the
          dependencies merge cleanly)

        .. SEEALSO::

        - :meth:`commit` -- record changes into the repository.

        - :meth:`local_tickets` -- list local tickets (you may want to
          commit your changes to a branch other than the current one).
        """
        if base == "commit":
            base = None
        if base == "dependencies":
            branch = self.git.current_branch()
            try:
                self.gather(self.trac.dependencies())
                self.git.diff("%s..%s"%(HEAD,branch))
            finally:
                self.git.checkout(branch)
        else:
            self.git.execute("diff", base)

    def prune_closed_tickets(self):
        """
        Remove branches for tickets that are already merged into master.

        .. SEEALSO::

        - :meth:`abandon_ticket` -- Abandon a single ticket or branch.
        """
        for branch in self.git.local_branches():
            if self.git.is_ancestor_of(branch, "master"):
                self._UI.show("Abandoning %s"%branch)
                self.git.abandon(branch)

    def abandon_ticket(self, ticket=None):
        """
        Abandon a ticket branch.

        INPUT:

        - ``ticket`` -- an integer or ``None`` (default: ``None``),
          remove the branch for ``ticket`` (or the current branch if
          ``None``). Also removes the users remote tracking branch.

        .. SEEALSO::

        - :meth:`prune_closed_tickets` -- abandon tickets that have
          been closed.

        - :meth:`local_tickets` -- list local tickets (by default only
          showing the non-abandoned ones).
        """
        if self._UI.confirm("Are you sure you want to delete your work on #%s?"%(ticketnum), default_yes=False):
            self.git.abandon(ticketnum)

    def gather(self, branchname, *tickets, **kwds):
        """
        Create a new branch ``branchname`` with ``tickets`` applied.

        INPUT:

        - ``branchname`` -- a string, the name of the new branch

        - ``tickets`` -- a list of integers or strings; for an
          integer, the branch on the trac ticket gets merged, for a
          string, that branch (or remote branch) gets merged.

        - ``create_dependency`` -- boolean (default ``True``, keyword
          only), whether to append the other ticket to the list of
          dependencies.  See :meth:`merge` for the consequences of
          having another branch as a dependency.

        - ``download`` -- boolean (default ``False``, keyword only),
          whether to download the most recent version of the other
          tickets before merging.

        .. SEEALSO::

        - :meth:`merge` -- merge into the current branch rather than
          creating a new one.

        - :meth:`show_dependencies` -- show the dependencies of a
          given branch.
        """
        create_dependencies = kwds.pop('create_dependencies',True)
        download = kwds.pop('download',False)
        if len(tickets) == 0:
            self._UI.show("Please include at least one input branch")
            return
        if self.git.branch_exists(branchname):
            if not self._UI.confirm("The %s branch already exists; do you want to merge into it?", default_yes=False):
                return
            self.git.execute_silent("checkout", branchname)
        else:
            self.switch_ticket(tickets[0], branchname=branchname, offline=not download)
            tickets = tickets[1:]
        for ticket in tickets:
            self.merge(ticket, create_dependency=create_dependencies,
                       download=download, message="Gathering %s into branch %s" % (ticket, branchname))

    def show_dependencies(self, ticket=None, all=False, _seen=None): # all = recursive
        """
        Show the dependencies of the given ticket.

        INPUT:

        - ``ticket`` -- string, int or None (default ``None``), the
          ticket for which dependencies are desired.  An int indicates
          a ticket number while a string indicates a branch name;
          ``None`` asks for the dependencies of the current ticket.

        - ``all`` -- boolean (default ``True``), whether to
          recursively list all tickets on which this ticket depends
          (in depth-first order), only including tickets that have a
          local branch.

        .. NOTE::

            Ticket dependencies are stored locally and only updated
            with respect to the remote server during :meth:`upload`
            and :meth:`download`.

        .. SEEALSO::

        - :meth:`TracInterface.dependencies` -- Query Trac to find
          dependencies.

        - :meth:`remote_status` -- will show the status of tickets
          with respect to the remote server.

        - :meth:`merge` -- Merge in changes from a dependency.

        - :meth:`diff` -- Show the changes in this branch over the
          dependencies.
        """
        if _seen is None:
            seen = []
        elif ticket in _seen:
            return
        else:
            seen = _seen
            seen.append(ticket)
        branchname = self.git._ticket_to_branch(ticket)
        if branchname is None:
            if _seen is not None: return
            raise ValueError("You must specify a valid ticket")
        dep = self._dependencies_as_tickets(branchname)
        if not all:
            self._UI.show("Ticket %s depends on %s"%(self._print(ticket), ", ".join([self._print(d) for d in dep])))
        else:
            for d in dep:
                self.show_dependencies(d, True, seen)
        if _seen is None:
            self._UI.show("Ticket %s depends on %s"%(self._print(ticket), ", ".join([self._print(d) for d in seen])))

    def merge(self, ticket="master", create_dependency=True, download=False, message=None):
        """
        Merge changes from another branch into the current branch.

        INPUT:

        - ``ticket`` -- string or int (default ``"master"``), a
          branch, ticket number or the current set of dependencies
          (indicated by the string ``"dependencies"``): the source of
          the changes to be merged.  If ``ticket = "dependencies"``
          then each updated dependency is merged in one by one,
          starting with the one listed first in the dependencies field
          on trac.  An int indicates a ticket number while a string
          indicates a branch name.

        - ``create_dependency`` -- boolean (default ``True``), whether
          to append the other ticket to the list of dependencies.

          Listing the other ticket as a dependency has the following
          consequences:

          - the other ticket must be positively reviewed and merged
            before this ticket may be merged into master.  The commits
            included from a dependency don't need to be reviewed in
            this ticket, whereas commits reviewed in this ticket from
            a non-dependency may make reviewing the other ticket
            easier.

          - you can more easily merge in future changes to
            dependencies.  So if you need a feature from another
            ticket it may be appropriate to create a dependency to
            that you may more easily benefit from others' work on that
            ticket.

          - if you depend on another ticket then you need to worry
            about the progress on that ticket.  If that ticket is
            still being actively developed then you may need to make
            many merges to keep up.

          Note that dependencies are stored locally and only updated
          with respect to the remote server during :meth:`upload` and
          :meth:`download`.

        - ``download`` -- boolean (default ``False``), whether to
          download the most recent version of the other ticket(s)
          before merging.

        .. SEEALSO::

        - :meth: `download` -- will download remote changes before
          merging.

        - :meth:`show_dependencies` -- see the current dependencies.

        - :meth:`GitInterface.merge` -- git's merge command has more
          options and can merge multiple branches at once.

        - :meth:`gather` -- creates a new branch to merge into rather
          than merging into the current branch.
        """
        curbranch = self.git.current_branch()
        if ticket == "dependencies":
            for d in self.git._dependencies[curbranch]:
                self.merge(d, False, download, message)
            return
        elif ticket is None:
            raise ValueError("You must specify a ticket to merge")
        ref = dep = None
        if download:
            remote_branch = self._remote_pull_branch(ticket)
            if remote_branch is not None:
                ref = self._fetch(remote_branch)
                dep = ticket
        if ref is None:
            dep = ref = self.git._ticket_to_branch(ticket)
        if ref is None:
            ref = ticket
            if isinstance(ticket, int):
                dep = ticket
        if create_dependency and dep and dep not in self.git._dependencies[curbranch]:
            self.git._dependencies[curbranch] += (dep,)
        if message is None:
            kwds = {}
        else:
            kwds = {'m':message}
        self.git.merge(ref, **kwds)

    def local_tickets(self, abandoned=False, quiet=False):
        """
        Print the tickets currently being worked on in your local
        repository.

        This function will show the branch names as well as the ticket
        numbers for all active tickets.  It will also show local
        branches that are not associated to ticket numbers.

        INPUT:

        - ``abandoned`` -- boolean (default ``False``), whether to show
          abandoned branches.

        - ``quite`` -- boolean (default ``False``), whether to show
          return the list of branches rather than printing them.

        .. SEEALSO::

        - :meth:`abandon_ticket` -- hide tickets from this method.

        - :meth:`remote_status` -- also show status compared to the
          trac server.

        - :meth:`current_ticket` -- get the current ticket.
        """
        raw_branches = self.git.read_output("branch").split()
        raw_branches.remove('*')
        branch_info = [(b, self.git._ticket[b]) for b in raw_branches
            if abandoned or not b.startswith("trash/")]
        if quiet:
            return branch_info
        else:
            for branch, ticket in branch_info:
                print ticket or '     ', '\t', branch

    def current_ticket(self, error=False):
        """
        Returns the current ticket as an int, or ``None`` if there is
        no current ticket.

        INPUT:

        - ``error`` -- boolean (default ``False``), whether to raise
          an error if there is no current ticket

        .. SEEALSO::

        - :meth:`local_tickets` -- show all local tickets.
        """
        curbranch = self.git.current_branch()
        if curbranch is not None and curbranch in self.git._ticket:
            return self.git._ticket[curbranch]
        if error: raise ValueError("You must specify a ticket")

    def vanilla(self, release="release"):
        """
        Returns to a basic release of Sage.

        INPUT:

        - ``release`` -- a string or decimal giving the release name.
          In fact, any tag, commit or branch will work.  If the tag
          does not exist locally an attempt to fetch it from the
          server will be made.

        Git equivalent::

            Checks out a given tag, commit or branch in detached head
            mode.

        .. SEEALSO::

        - :meth:`switch_ticket` -- switch to another branch, ready to
          develop on it.

        - :meth:`download` -- download a branch from the server and
          merge it.

        """
        if hasattr(release, 'literal'):
            release = release.literal
        release = str(release)
        if self._UI.confirm("Are you sure you want to revert to %s?"%(release)):
            self.git.switch_branch(release, detached = True)

    ##
    ## Auxilliary functions
    ##

    def _fetch(self, branch, repository=None):
        """
        Fetches ``branch`` from the remote repository, returning the name of
        the newly-updated local ref.

        INPUT:

        - ``branch`` -- name of a remote branch

        - ``repo`` -- name of a remote repository

        OUTPUT:

        The name of a newly created/updated local ref.

        """
        if repository is None:
            repository = self.git._repo
        local_ref = "refs/remotes/trac/%s" % branch
        self.git.fetch(repository, "%s:%s" % (branch, local_ref))
        return local_ref

    @property
    def tmp_dir(self):
        try:
            return self._tmp_dir
        except AttributeError:
            self._tmp_dir = tempfile.mkdtemp()
            atexit.register(shutil.rmtree, self._tmp_dir)
            return self._tmp_dir

    def _detect_patch_diff_format(self, lines):
        """
        Determine the format of the ``diff`` lines in ``lines``.

        INPUT:

        - ``lines`` -- a list of strings

        OUTPUT:

        Either ``git`` (for ``diff --git`` lines) or ``hg`` (for ``diff -r`` lines).

        .. NOTE::

            Most Sage developpers have configured mercurial to export
            patches in git format.

        EXAMPLES::

            sage: dev._detect_patch_diff_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"])
            'hg'
            sage: dev._detect_patch_diff_format(
            ....:     ["diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            'git'

            sage: import os.path
            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_diff_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'git'
            sage: dev._detect_patch_diff_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","diff.patch"
            ....:         )).read().splitlines())
            'hg'

        TESTS::

            sage: dev._detect_patch_diff_format(["# HG changeset patch"])
            Traceback (most recent call last):
            ...
            NotImplementedError: Failed to detect diff format.
            sage: dev._detect_patch_diff_format(
            ... ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py",
            ...  "diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            Traceback (most recent call last):
            ...
            ValueError: File appears to have mixed diff formats.

        """
        format = None
        regexs = { "hg" : HG_DIFF_REGEX, "git" : GIT_DIFF_REGEX }

        for line in lines:
            for name,regex in regexs.items():
                if regex.match(line):
                    if format is None:
                        format = name
                    if format != name:
                        raise ValueError("File appears to have mixed diff formats.")

        if format is None:
            raise NotImplementedError("Failed to detect diff format.")
        else:
            return format

    def _detect_patch_path_format(self, lines, diff_format = None):
        """
        Determine the format of the paths in the patch given in ``lines``.

        INPUT:

        - ``lines`` -- a list (or iterable) of strings

        - ``diff_format`` -- ``'hg'``,``'git'``, or ``None`` (default:
          ``None``), the format of the ``diff`` lines in the patch. If
          ``None``, the format will be determined by
          :meth:`_detect_patch_diff_format`.

        OUTPUT:

        A string, ``'new'`` (new repository layout) or ``'old'`` (old
        repository layout).

        EXAMPLES::

            sage: dev._detect_patch_path_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"])
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py"],
            ....:     diff_format="git")
            Traceback (most recent call last):
            ...
            NotImplementedError: Failed to detect path format.
            sage: dev._detect_patch_path_format(
            ....:     ["diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi"])
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi"])
            'new'

            sage: import os.path
            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_path_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'old'

        """
        lines = list(lines)
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        path_format = None

        if diff_format == "git":
            diff_regexs = (GIT_DIFF_REGEX, PM_DIFF_REGEX)
        elif diff_format == "hg":
            diff_regexs = (HG_DIFF_REGEX, PM_DIFF_REGEX)
        else:
            raise NotImplementedError(diff_format)

        regexs = { "old" : HG_PATH_REGEX, "new" : GIT_PATH_REGEX }

        for line in lines:
            for regex in diff_regexs:
                match = regex.match(line)
                if match:
                    for group in match.groups():
                        for name, regex in regexs.items():
                            if regex.match(group):
                                if path_format is None:
                                    path_format = name
                                if path_format != name:
                                    raise ValueError("File appears to have mixed path formats.")

        if path_format is None:
            raise NotImplementedError("Failed to detect path format.")
        else:
           return path_format

    def _rewrite_patch_diff_paths(self, lines, to_format, from_format=None, diff_format=None):
        """
        Rewrite the ``diff`` lines in ``lines`` to use ``to_format``.

        INPUT:

        - ``lines`` -- a list or iterable of strings

        - ``to_format`` -- ``'old'`` or ``'new'``

        - ``from_format`` -- ``'old'``, ``'new'``, or ``None`` (default:
          ``None``), the current formatting of the paths; detected
          automatically if ``None``

        - ``diff_format`` -- ``'git'``, ``'hg'``, or ``None`` (default:
          ``None``), the format of the ``diff`` lines; detected automatically
          if ``None``

        OUTPUT:

        A list of string, ``lines`` rewritten to conform to ``lines``.

        EXAMPLES:

        Paths in the old format::

            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="old")
            ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="old")
            ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="old", diff_format="git")
            ['--- a/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/sage/rings/padics/pow_computer_ext.pxd']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="new")
            ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="new")
            ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="new", diff_format="git")
            ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/src/sage/rings/padics/pow_computer_ext.pxd']

        Paths in the new format::

            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="old")
            ['diff -r 1492e39aff50 -r 5803166c5b11 sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="old")
            ['diff --git a/sage/rings/padics/FM_template.pxi b/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/src/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="old", diff_format="git")
            ['--- a/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/sage/rings/padics/pow_computer_ext.pxd']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py'],
            ....:     to_format="new")
            ['diff -r 1492e39aff50 -r 5803166c5b11 src/sage/schemes/elliptic_curves/ell_rational_field.py']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi'],
            ....:     to_format="new")
            ['diff --git a/src/sage/rings/padics/FM_template.pxi b/src/sage/rings/padics/FM_template.pxi']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
            ....:      '+++ b/src/sage/rings/padics/pow_computer_ext.pxd'],
            ....:     to_format="new", diff_format="git")
            ['--- a/src/sage/rings/padics/pow_computer_ext.pxd',
             '+++ b/src/sage/rings/padics/pow_computer_ext.pxd']

            sage: import os.path
            sage: from sage.env import SAGE_SRC
            sage: result = dev._rewrite_patch_diff_paths(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines(),
            ....:     to_format="new", diff_format="git")
            sage: len(result)
            2980
            sage: result[0]
            '#8703: Enumerated sets and data structure for ordered and binary trees'
            sage: result[12]
            'diff --git a/src/doc/en/reference/combinat/index.rst b/src/doc/en/reference/combinat/index.rst'

        """
        lines = list(lines)
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        if from_format is None:
            from_format = self._detect_patch_path_format(lines, diff_format=diff_format)

        if to_format == from_format:
            return lines

        def hg_path_to_git_path(path):
            if any([path.startswith(p) for p in "module_list.py","setup.py","c_lib/","sage/","doc/"]):
                return "src/%s"%path
            else:
                raise NotImplementedError("mapping hg path `%s`"%path)

        def git_path_to_hg_path(path):
            if any([path.startswith(p) for p in "src/module_list.py","src/setup.py","src/c_lib/","src/sage/","src/doc/"]):
                return path[4:]
            else:
                raise NotImplementedError("mapping git path `%s`"%path)

        def apply_replacements(lines, diff_regexs, replacement):
            ret = []
            for line in lines:
                for diff_regex in diff_regexs:
                    m = diff_regex.match(line)
                    if m:
                        line = line[:m.start(1)] + ("".join([ line[m.end(i-1):m.start(i)]+replacement(m.group(i)) for i in range(1,m.lastindex+1) ])) + line[m.end(m.lastindex):]
                ret.append(line)
            return ret

        diff_regex = None
        if diff_format == "hg":
            diff_regex = (HG_DIFF_REGEX, PM_DIFF_REGEX)
        elif diff_format == "git":
            diff_regex = (GIT_DIFF_REGEX, PM_DIFF_REGEX)
        else:
            raise NotImplementedError(diff_format)

        if from_format == "old":
            return self._rewrite_patch_diff_paths(apply_replacements(lines, diff_regex, hg_path_to_git_path), from_format="new", to_format=to_format, diff_format=diff_format)
        elif from_format == "new":
            if to_format == "old":
                return apply_replacements(lines, diff_regex, git_path_to_hg_path)
            else:
                raise NotImplementedError(to_format)
        else:
            raise NotImplementedError(from_format)

    def _detect_patch_header_format(self, lines):
        """
        Detect the format of the patch header in ``lines``.

        INPUT:

        - ``lines`` -- a list (or iterable) of strings

        OUTPUT:

        A string, ``'hg-export'`` (mercurial export header), ``'hg'``
        (mercurial header), ``'git'`` (git mailbox header), ``'diff'`` (no
        header)

        EXAMPLES::

            sage: dev._detect_patch_header_format(
            ... ['# HG changeset patch','# Parent 05fca316b08fe56c8eec85151d9a6dde6f435d46'])
            'hg'
            sage: dev._detect_patch_header_format(
            ... ['# HG changeset patch','# User foo@bar.com'])
            'hg-export'
            sage: dev._detect_patch_header_format(
            ... ['From: foo@bar'])
            'git'

            sage: import os.path
            sage: from sage.env import SAGE_SRC
            sage: dev._detect_patch_header_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'diff'
            sage: dev._detect_patch_header_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","diff.patch"
            ....:         )).read().splitlines())
            'diff'

        """
        lines = list(lines)
        if not lines:
            raise ValueError("patch is empty")

        if HG_HEADER_REGEX.match(lines[0]):
            if HG_USER_REGEX.match(lines[1]):
                return "hg-export"
            elif HG_PARENT_REGEX.match(lines[1]):
                return "hg"
        elif GIT_FROM_REGEX.match(lines[0]):
            return "git"
        return "diff"
        #raise NotImplementedError("Failed to determine patch header format.")

    def _detect_patch_modified_files(self, lines, diff_format = None):
        if diff_format is None:
            diff_format = self._detect_patch_diff_format(lines)

        if diff_format == "hg":
            regex = HG_DIFF_REGEX
        elif diff_format == "git":
            regex = GIT_DIFF_REGEX
        else:
            raise NotImplementedError(diff_format)

        ret = set()
        for line in lines:
            m = regex.match(line)
            if m:
                for group in m.groups():
                    split = group.split('/')
                    if split:
                        ret.add(split[-1])
        return list(ret)

    def _rewrite_patch_header(self, lines, to_format, from_format = None, diff_format = None):
        """
        Rewrite ``lines`` to match ``to_format``.

        INPUT:

        - ``lines`` -- a list of strings, the lines of the patch file

        - ``to_format`` -- one of ``'hg'``, ``'hg-export'``, ``'diff'``,
          ``'git'``, the format of the resulting patch file.

        - ``from_format`` -- one of ``None``, ``'hg'``, ``'hg-export'``, ``'diff'``, ``'git'``
          (default: ``None``), the format of the patch file.  The format is
          determined automatically if ``format`` is ``None``.

        OUTPUT:

        A list of lines, in the format specified by ``to_format``.

        Some sample patch files are in data/, in hg and git
        format. Since the translation is not perfect, the resulting
        file is also put there for comparison.

        EXAMPLES::

            sage: import os.path
            sage: from sage.env import SAGE_SRC
            sage: hg_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "data", "hg.patch")
            ....:     ).read().splitlines()
            sage: hg_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "data", "hg-output.patch")
            ....:     ).read().splitlines()
            sage: git_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "data", "git.patch")
            ....:     ).read().splitlines()
            sage: git_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "data", "git-output.patch")
            ....:     ).read().splitlines()

            sage: dev._rewrite_patch_header(git_lines, 'git') == git_lines
            True
            sage: dev._rewrite_patch_header(hg_lines, 'hg-export') == hg_lines
            True

            sage: dev._rewrite_patch_header(git_lines, 'hg-export') == hg_output_lines
            True
            sage: dev._rewrite_patch_header(hg_lines, 'git') == git_output_lines
            True

            sage: dev._rewrite_patch_header(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines(), 'git')[:5]
            ['From: "Unknown User" <unknown@sagemath.org>',
            'Subject: #8703: Enumerated sets and data structure for ordered and binary trees',
            'Date: ...',
            '',
            '- The Class Abstract[Labelled]Tree allows for inheritance from different']

        """
        lines = list(lines)
        if not lines:
            raise ValueError("empty patch file")

        if from_format is None:
            from_format = self._detect_patch_header_format(lines)

        if from_format == to_format:
            return lines

        def parse_header(lines, regexs, mandatory=False):
            header = {}
            i = 0
            for (key, regex) in regexs:
                if i > len(lines):
                    if mandatory:
                        raise ValueError("Malformed patch. Missing line for regular expression `%s`."%(regex.pattern))
                    else:
                        return
                match = regex.match(lines[i])
                if match is not None:
                    if len(match.groups()) > 0:
                        header[key] = match.groups()[0]
                    i += 1
                elif mandatory:
                    raise ValueError("Malformed patch. Line `%s` does not match regular expression `%s`."%(lines[i],regex.pattern))

            message = []
            for i in range(i,len(lines)):
                if lines[i].startswith("diff -"):
                    break
                else:
                    message.append(lines[i])

            header["message"] = message
            return header, lines[i:]

        if from_format == "git":
            header, diff = parse_header(lines, (("user", GIT_FROM_REGEX), ("subject", GIT_SUBJECT_REGEX), ("date", GIT_DATE_REGEX)),
                                        mandatory=True)

            if to_format == "hg-export":
                ret = []
                ret.append('# HG changeset patch')
                ret.append('# User %s'%(header["user"]))
                ret.append('# Date %s 00000'%int(time.mktime(email.utils.parsedate(header["date"])))) # this is not portable and the time zone is wrong
                ret.append('# Node ID 0000000000000000000000000000000000000000')
                ret.append('# Parent  0000000000000000000000000000000000000000')
                ret.append(header["subject"])
                ret.extend(header["message"])
                ret.extend(diff)
                return ret
            else:
                raise NotImplementedError(to_format)
        elif from_format in ["hg", "diff", "hg-export"]:
            header, diff = parse_header(lines,
                                        (("hg_header", HG_HEADER_REGEX),
                                         ("user", HG_USER_REGEX),
                                         ("date", HG_DATE_REGEX),
                                         ("node", HG_NODE_REGEX),
                                         ("parent", HG_PARENT_REGEX)))
            user    = header.get("user", '"Unknown User" <unknown@sagemath.org>')
            date    = email.utils.formatdate(int(header.get("date", time.time())))
            message = header.get("message", [])
            if message:
                subject = message[0]
                message = message[1:]
            else:
                subject = 'No Subject. Modified: %s'%(", ".join(self._detect_patch_modified_files(lines)))
            ret = []
            ret.append('From: %s'%user)
            ret.append('Subject: %s'%subject)
            ret.append('Date: %s'%date)
            ret.append('')
            if message and message != ['']: # avoid a double empty line
                ret.extend(message)
            ret.extend(diff)
            return self._rewrite_patch_header(ret, to_format=to_format, from_format="git", diff_format=diff_format)
        else:
            raise NotImplementedError(from_format)

    def _rewrite_patch(self, lines, to_path_format, to_header_format, from_diff_format=None, from_path_format=None, from_header_format=None):
        return self._rewrite_patch_diff_paths(self._rewrite_patch_header(lines, to_format=to_header_format, from_format=from_header_format, diff_format=from_diff_format), to_format=to_path_format, diff_format=from_diff_format, from_format=from_path_format)

    def _dependency_join(self, ticketnum=None):
        if ticketnum is None:
            ticketnum = self.current_ticket(error=True)
        for d in self.trac.dependencies(ticketnum):
            pass
        raise NotImplementedError

    @property
    def git(self):
        """
        A lazy property to provide a :class:`git_interface.GitInterface`.

        EXAMPLES::

            sage: g = dev.git; g
            GitInterface()
            sage: g is dev.git
            True

        """
        if self.__git is None:
            self.__git = GitInterface(self)
        return self.__git

    @property
    def trac(self):
        """
        A lazy property to provide a :class:`trac_interface.TracInterface`

        EXAMPLES::

            sage: t = dev.trac; t
             <sage.dev.trac_interface.TracInterface object at ...>
            sage: t is dev.trac
            True

        """
        if self.__trac is None:
            self.__trac = TracInterface(self)
        return self.__trac

    def __repr__(self):
        """
        Return a printable representation of this object.

        EXAMPLES::

            sage: dev # indirect doctest
            SageDev()

        """
        return "SageDev()"

    def _upload_ssh_key(self, keyfile, create_key_if_not_exists=True):
        """
        Upload ``keyfile`` to gitolite through the trac interface.

        INPUT:

        - ``keyfile`` -- the absolute path of the key file (default:
          ``~/.ssh/id_rsa``)

        - ``create_key_if_not_exists`` -- use ``ssh-keygen`` to create
          ``keyfile`` if ``keyfile`` or ``keyfile.pub`` does not exist.

        EXAMPLES::

            sage: import tempfile, os
            sage: tmp = tempfile.NamedTemporaryFile().name
            sage: dev._upload_ssh_key(tmp, create_key_if_not_exists = False)
            Traceback (most recent call last):
            ...
            IOError: [Errno 2] No such file or directory: ...
            sage: dev._upload_ssh_key(tmp, create_key_if_not_exists = True)
            Generating ssh key....
            Ssh key successfully generated
            sage: os.unlink(tmp)
        """
        cfg = self._config

        try:
            with open(keyfile, 'r') as F:
                pass
            with open(keyfile + '.pub', 'r') as F:
                pass
        except IOError:
            if create_key_if_not_exists:
                self._UI.show("Generating ssh key....")
                success = call(["ssh-keygen", "-q", "-f", keyfile, "-P", ""])
                if success == 0:
                    self._UI.show("Ssh key successfully generated")
                else:
                    raise RuntimeError("Ssh key generation failed.  Please create a key in `%s` and retry"%(keyfile))
            else:
                raise

        with open(keyfile + '.pub', 'r') as F:
            pubkey = F.read().strip()

        self.trac.sshkeys.setkeys(pubkey)

    def _trac_branch(self, ticket):
        D = self.trac._get_attributes(ticket)
        if 'branch' not in D: return None
        branch = D['branch']
        if branch: return branch

    def _remote_pull_branch(self, ticket):
        branchname = self.git._ticket_to_branch(ticket)
        remote_branch = self.git._remote[branchname]
        if remote_branch is None:
            userspace = True
        else:
            x = remote_branch.split('/')
            userspace = (x[0] == 'u' and x[1] == self.trac._username)
        if userspace and isinstance(ticket, int):
            remote_branch = self._trac_branch(ticket)
        return remote_branch

    def _print_ticket(self, ticket):
        if isinstance(ticket, int):
            return "#%s"%(ticket)
        ticket = self.git._ticket_to_branch(ticket)
        if ticket in self.git._ticket:
            return "#%s"%(self.git._ticket[ticket])
        return str(ticket)

    def _dependencies_as_tickets(self, branch):
        dep = self.git._dependencies[branch]
        dep = [d if isinstance(d, int) else self._ticket[d] for d in dep]
        dep = [d for d in dep if d]
        return dep

    def _save_uncommitted_changes(self):
        """
        Returns True if changes should be unstashed
        """
        curbranch = self.git.current_branch()
        if curbranch is None:
            options = ["new branch", "stash"]
        else:
            options = ["current branch", "new branch", "stash"]
        dest = self._UI.get_input("Where do you want to commit your changes?", options)
        if dest == "stash":
            self.git.stash()
        elif dest == "new branch":
            success = self.git.execute_silent("stash")
            if success != 0:
                raise RuntimeError("Failed to stash changes")
            return True
        elif dest == "current branch":
            self.commit()

    def _unstash_changes(self):
        success = self.git.execute_silent("stash", "apply")
        if success == 0:
            self.git.execute_silent("stash", "drop")
        else:
            self.git.execute_silent("reset", hard=True)
            self._UI.show("Changes did not apply cleanly to the new branch.  They are now in your stash")

def doctest_config():
    """
    creates a fake configuration used for doctesting

    EXAMPLE::

        sage: from sage.dev.sagedev import doctest_config
        sage: doctest_config()
        Config('''
        [git]
        dot_git = ...
        doctest = ...
        [trac]
        username = doctest
        ''')
    """
    ret = Config(devrc = tempfile.NamedTemporaryFile().name)
    tmp_dir = tempfile.mkdtemp()
    atexit.register(shutil.rmtree, tmp_dir)
    ret['git'] = {}
    ret['git']['dot_git'] = os.path.join(tmp_dir, '.git')
    ret['git']['doctest'] = tmp_dir
    ret['trac'] = {}
    ret['trac']['username'] = 'doctest'
    return ret

# default sagedev object
if DOCTEST_MODE:
    dev = SageDev(doctest_config())
else:
    dev = SageDev()
