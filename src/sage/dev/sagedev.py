r"""
Sagedev

TESTS::

    sage: from sage.dev.sagedev import SageDev
    sage: from sage.dev.git_interface import GitInterface
    sage: remote = GitInterface(...)
    sage: remote.init()
    sage: dev1 = SageDev(...)
    sage: dev2 = SageDev(...)

Collaboration::

    sage: # developer 1 creates a revision:
    sage: os.chdir(dev1.git._tmp_dir)
    sage: ticket = dev1.create_ticket()
    Created ticket #14366 (https://trac.sagemath.org/14366).
    Switched to branch 'ticket/14366'
    sage: with open('a_file', 'w') as f:
    ....:     f.write("revision 1")
    ....:
    sage: dev1.git.add('a_file')  # when not doctesting, you would do this interactively
    sage: dev1.commit(message="revision 1")
    Are you sure you want to save your changes to ticket #14366? [Yes/no]
    [ticket/14366 ...] revision 1
     1 file changed, 1 insertion(+)
     create mode 100644 a_file
    sage: dev1.diff()
    sage: dev1.upload()
    There does not seem to be a branch u/user1/ticket/14366 on the remote server yet. Do you want to create such a branch? [Yes/no]
    To ...
     * [new branch]      ticket/14366 -> u/user1/ticket/14366
    Ticket 14366 now refers to your branch u/user1/ticket/14366.
    sage: # developer 2 works on the ticket
    sage: os.chdir(dev2.git._tmp_dir)
    sage: dev2.switch_ticket(ticket)  # (or download??)
    Switched to branch 'ticket/14366'
    sage: dev2.download()
    sage: open('a_file').read()
    "revision 1"
    sage: with open('a_file, 'w') as f:
    ....:     f.write("revision 2a")
    ....:
    sage: dev2.commit()
    sage: dev2.upload()
    There does not seem to be a branch u/user2/ticket/14366 on the remote server yet. Do you want to create such a branch? [Yes/no]
    To ...
     * [new branch]      ticket/14366 -> u/user2/ticket/14366
    Ticket 14366 now refers to your branch u/user2/ticket/14366.

Merge conflicts (developer 1 is behind)::

    sage: $ developer 1 tries to work simultaneously
    sage: os.chdir(dev1.git._tmp_dir)
    sage: with open('a_file', 'w') as f:
    ....:     f.write("revision 2b")
    ....:
    sage: dev1.commit()
    sage: dev1.upload()
    ##### conflict!
    sage: dev1.download()
    #### resolve something
    sage: dev1.upload()

THINGS TO TEST:

    * not on a branch
    * in a merge / rebase /etc
    * pending changes
    * untracked files that are overwritten
    * no such remote branch
    * no such ticket
    * network errors, access errors
    * ticket<->branch association invalid
"""
import atexit
import email.utils
import os
import random
import re
import shutil
import tempfile
import time
import urllib2

from datetime import datetime
from subprocess import call, check_call

from git_interface import GitInterface
from trac_interface import TracInterface
from user_interface import UserInterface
from cmd_line_interface import CmdLineInterface
from config import Config
from saving_dict import SavingDict

from sage.env import DOT_SAGE, TRAC_SERVER_URI
from sage.doctest import DOCTEST_MODE

# regular expressions to parse mercurial patches
HG_HEADER_REGEX = re.compile(r"^# HG changeset patch$")
HG_USER_REGEX = re.compile(r"^# User (.*)$")
HG_DATE_REGEX = re.compile(r"^# Date (\d+) (-?\d+)$")
HG_NODE_REGEX = re.compile(r"^# Node ID ([0-9a-f]+)$")
HG_PARENT_REGEX = re.compile(r"^# Parent +([0-9a-f]+)$")
HG_DIFF_REGEX = re.compile(r"^diff (?:-r [0-9a-f]+ ){1,2}(.*)$")
PM_DIFF_REGEX = re.compile(r"^(?:(?:\+\+\+)|(?:---)) [ab]/([^ ]*)(?: .*)?$")
MV_DIFF_REGEX = re.compile(r"^rename (?:(?:to)|(?:from)) (.*)$")

# regular expressions to parse git patches -- at least those created by us
GIT_FROM_REGEX = re.compile(r"^From: (.*)$")
GIT_SUBJECT_REGEX = re.compile(r"^Subject: (.*)$")
GIT_DATE_REGEX = re.compile(r"^Date: (.*)$")
GIT_DIFF_REGEX = re.compile(r"^diff --git a/(.*) b/(.*)$") # this regex should work for our patches since we do not have spaces in file names

# regular expressions to determine whether a path was written for the new git
# repository of for the old hg repository
HG_PATH_REGEX = re.compile(r"^(?=sage/)|(?=doc/)|(?=module_list\.py)|(?=setup\.py)|(?=c_lib/)")
GIT_PATH_REGEX = re.compile(r"^(?=src/)")

# the name of the branch which holds the vanilla clone of sage - in the long
# run this should be "master", currently, "public/sage-git/master" contains some changes
# over "master" which have not been reviewed yet but which are needed to work
# using git
MASTER_BRANCH = "public/sage-git/master"

TracConnectionError = RuntimeError("could not connect with trac server")

class SageDev(object):
    r"""
    The developer interface for sage.

    This class facilitates access to git and trac.

    EXAMPLES::

        sage: type(dev)
        <class 'sage.dev.sagedev.SageDev'>
    """
    def __init__(self, config = Config(), UI=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: type(SageDev(doctest_config()))
            <class 'sage.dev.sagedev.SageDev'>
        """
        self._config = config
        if UI is not None:
            self._UI = UI
        else:
            self._UI = CmdLineInterface()

        self.__git = None
        self.__trac = None

        ticket_file = self._config.get('ticketfile',
                os.path.join(DOT_SAGE, 'branch_to_ticket'))
        branch_file = self._config.get('branchfile',
                os.path.join(DOT_SAGE, 'ticket_to_branch'))
        dependencies_file = self._config.get('dependenciesfile',
                os.path.join(DOT_SAGE, 'dependencies'))
        remote_branches_file = self._config.get('remotebranchesfile',
                os.path.join(DOT_SAGE, 'remote_branches'))

        self._ticket = SavingDict(ticket_file)
        self._branch = SavingDict(branch_file, paired=self._ticket)
        self._dependencies = SavingDict(dependencies_file, default=tuple)
        self._remote = SavingDict(remote_branches_file)

    @property
    def tmp_dir(self):
        r"""
        a lazy property to provide a temporary directory

        TESTS::

            sage: import os
            sage: os.path.isdir(dev.tmp_dir)
            True
        """
        try:
            return self._tmp_dir
        except AttributeError:
            self._tmp_dir = tempfile.mkdtemp()
            atexit.register(shutil.rmtree, self._tmp_dir)
            return self._tmp_dir

    @property
    def git(self):
        r"""
        a lazy property to provide a
        :class:`sage.dev.git_interface.GitInterface`

        TESTS::

            sage: g = dev.git; g
            GitInterface()
            sage: g is dev.git
            True
        """
        try:
            return self._git
        except AttributeError:
            self._git = GitInterface(self)
            return self._git

    @property
    def trac(self):
        r"""
        a lazy property to provide a
        :class:`sage.dev.trac_interface.TracInterface`

        TESTS::

            sage: t = dev.trac; t
            <sage.dev.trac_interface.TracInterface object at ...>
            sage: t is dev.trac
            True
        """
        try:
            return self._trac
        except AttributeError:
            self._trac = TracInterface(self)
            return self._trac

    def __repr__(self):
        r"""
        return a printable representation of this object

        TESTS::

            sage: dev # indirect doctest
            SageDev()
        """
        return "SageDev()"

    ##
    ## Public interface
    ##

    def switch_ticket(self, ticket, branchname=None):
        r"""
        switch to a branch associated to ``ticket``

        If there is already a local branch for ``ticket`` then
        :meth:`remote_status` will be called, otherwise the
        current branch attached to ``ticket`` will be downloaded
        from trac.

        INPUT:

        - ``ticket`` -- ticket to switch to

          If ``ticket`` is a string, then this method just switches
          to the branch ``ticket`` (in which case ``branchname``
          must me ``None``).

          If ``ticket`` is an integer, it must be the
          number of a ticket.

        - ``branchname`` -- branch to that stores changes for ``ticket``
          (default: ticket/``ticket``)

        .. SEEALSO::

        - :meth:`download` -- download the ticket from the remote
          server and merge in the changes

        - :meth:`create_ticket` -- make a new ticket

        - :meth:`vanilla` -- switch to a released version of Sage
        """
        try:
            ticket = int(ticket)
            if branchname is None:
                branchname = self._branch.get(ticket, 'ticket/%s'%ticket)
        except ValueError:
            if not isinstance(ticket, basestring):
                raise ValueError("bad ticket value")
            elif branchname is not None:
                raise ValueError("cannot specify two branch names")
            elif not self.git.branch_exists(ticket):
                raise ValueError("branch does not exist")
            branchname = ticket

        create = not self.git.branch_exists(branchname)
        if create:
            try:
                tracbranch = self._trac_branch(ticket)
                ref = self._fetch(tracbranch)
            except (urllib2.HTTPError, KeyError) as e:
                # there is not yet a branch on trac for this ticket
                tracbranch = None
                ref = MASTER_BRANCH
            self.git.create_branch(branchname, ref)

        self.git.switch_branch(branchname)

        if isinstance(ticket, int):
            self._branch.setdefault(ticket, branchname)
            try:
                if not create:
                    self.remote_status(branchname)
            except TracConnectionError:
                pass

    def edit_ticket(self, ticket=None):
        r"""
        Edit the description of ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or ``None`` (default: ``None``), the number
          of the ticket to edit. If ``None``, edit the :meth:`current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`add_comment`

        """
        if ticket is None:
            ticket = self.current_ticket()

        if ticket is None:
            raise ValueError("must specify a ticket")

        self.trac.edit_ticket(ticket)

    def add_comment(self, ticket=None):
        r"""
        Add a comment to ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or ``None`` (default: ``None``), the number
          of the ticket to edit. If ``None``, edit the :meth:`current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`edit_ticket`

        """
        if ticket is None:
            ticket = self.current_ticket()

        if ticket is None:
            raise ValueError("must specify a ticket")

        self.trac.add_comment(ticket)

    def browse_ticket(self, ticket=None):
        r"""
        start a webbrowser at the ticket page on Sage trac

        INPUT:
        - ``ticket`` -- an integer or ``None`` (default: ``None``), the number
          of the ticket to edit. If ``None``, edit the :meth:`current_ticket`.
        """
        from sage.misc.viewer import browser
        if ticket is None:
            ticket = self.current_ticket()

        if ticket is None:
            raise ValueError("must specify a ticket")

        browser_cmdline = browser() + ' ' + TRAC_SERVER_URI + '/ticket/' + str(ticket)
        os.system(browser_cmdline)

    def create_ticket(self,
            branchname=None, base=MASTER_BRANCH, remote_branch=None):
        r"""
        create a new ticket on trac and switch to a new local branch to
        work on said ticket

        INPUT:

        - ``branchname`` -- the name of the local branch that will be
          used for the new ticket (default: ticket/ticketnum)

        - ``base`` -- a branch on which to base the ticket
          (default: ``"master"``)

          If instead ``base`` is set to ``None`` then the current ticket
          is used. If the base is set to anything other than
          ``"master"`` then the corresponding dependency will be added.

        - ``remote_branch`` -- the name of the remote branch this
          ticket should be tracking (default: ``None``)

          alternatively you can set ``remote_branch`` to ``True`` and
          it will use the canonical ``remote_branch``

        .. SEEALSO::

        - :meth:`switch_ticket` -- Switch to an already existing ticket.

        - :meth:`download` -- Download changes from an existing ticket.
        """
        try:
            ticketnum = self.trac.create_ticket()
        except urllib2.HTTPError:
            raise TracConnectionError
        if ticketnum is None:
            # user aborted
            return
        if branchname is None:
            branchname = 'ticket/%s'%ticketnum
        if base is None:
            base = self.git.current_branch()
        self.git.create_branch(branchname, base, remote_branch)
        self._ticket[branchname] = ticketnum
        if base != MASTER_BRANCH:
            self._dependencies[branchname] = [base]
        self.git.switch_branch(branchname)
        return ticketnum

    def commit(self, message=None, interactive=False):
        r"""
        create a commit from the pending changes on the current branch

        This is most akin to mercurial's commit command, not git's.

        INPUT:

        - ``message`` -- the message of the commit (default: ``None``)

          if ``None``, prompt for a message

        - ``interactive`` -- if set, interactively select which part of the
          changes should be part of the commit

        .. SEEALSO::

        - :meth:`upload` -- Upload changes to the remote server.  This
          is the next step once you've committed some changes.

        - :meth:`diff` -- Show changes that will be committed.
        """
        curticket = self.current_ticket()

        prompt = "Are you sure you want to save your changes to "
        if curticket is None:
            curbranch = self.git.current_branch()
            prompt += "branch %s?"%curbranch
        else:
            prompt += "ticket %s?"%self._ticket_repr(curticket)
        if not self._UI.confirm(prompt):
            self._UI.show("If you want to commit these changes to another "+
                          "ticket use the switch_ticket() method")
            return

        for file in self.git.unknown_files():
            if self._UI.confirm("Would you like to commit %s"%file,
                    default_no=True):
                self.git.add(file)

        kwds = {}
        if interactive:
            kwds['patch'] = True
        else:
            kwds['all'] = True

        if message is not None:
            kwds['message'] = message

        self.git.commit(**kwds)

    def upload(self, ticket=None, remote_branch=None, force=False, repository=None):
        r"""
        Upload the current branch to the Sage repository

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

        - :meth:`download` -- Update a ticket with changes from the
          remote repository.
        """
        if repository is None:
            repository = self.git._repo
        branch = self.git.current_branch()
        try:
            oldticket = self._ticket[branch]
            if ticket is None:
                ticket = oldticket
            elif oldticket != ticket:
                if not self._UI.confirm("Are you sure you want to upload "+
                                        "your changes to ticket "+
                                        "%s "%self._ticket_repr(ticket)+
                                        "instead of "+
                                        "%s?"%self._ticket_repr(oldticket),
                                        default_no=True):
                    return
                self._ticket[branch] = ticket
        except KeyError:
            if ticket is None and remote_branch is None:
                raise ValueError("You don't have a ticket for this branch "+
                                 "(%s)"%branch)

        remote_branch = remote_branch or self.git._local_to_remote_name(branch)
        ref = None
        try:
            ref = self._fetch(remote_branch, repository=repository)
        except RuntimeError:
            if not self._UI.confirm("There does not seem to be a branch %s on the remote server yet. Do you want to create such a branch?"%remote_branch, default_no=False):
                return

        if ref and not (self.git.is_ancestor_of(ref, branch) or force):
            if self._UI.confirm("Changes not compatible with remote branch %s; consider downloading first. Are you sure you want to continue?"%remote_branch, default_no=True):
                force = True
            else: return

        self.git.push(repository, "%s:%s"%(branch, remote_branch), force=force)

        if ticket:
            ticket_branch = self.trac._get_attributes(ticket).get("branch", "").strip()
            if ticket_branch:
                ref = None
                try:
                    ref = self._fetch(ticket_branch, repository=repository)
                except RuntimeError: # no such branch
                    self._UI.show("The ticket %s refers to a non-existant branch %s - will overwrite the branch field on the ticket with %s"%(ticket,ticket_branch,remote_branch))

                if ref and not (self.git.is_ancestor_of(ref, branch) or force):
                    if not self._UI.show("Your changes would discard some of the commits on the current branch %s of the ticket %s. Download these changes first or use 'force' to overwrite them."%(ticket_branch,ticket)):
                        return

            git_deps = ", ".join(["#%s"%d for d in self._dependencies_as_tickets(branch)])
            self.trac.update(ticket, branch=remote_branch, dependencies=git_deps)
            self._UI.show("Ticket %s now refers to your branch %s."%(ticket,remote_branch))

    def download(self, ticket=None, branchname=None, force=False, repository=None):
        r"""
        download the changes made to a remote branch into a given
        ticket or the current branch

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
            try:
                ticket = self._ticket[branch]
            except KeyError:
                pass
        else:
            ticket = int(ticket)
            try:
                branch = self._branch[ticket]
            except KeyError:
                if branchname is None:
                    branch = 'ticket/%s'%ticket
                else:
                    branch = branchname
        remote_branch = self._remote_pull_branch(ticket or branch)
        ref = self._fetch(remote_branch, repository=repository)
        if force:
            self.git.branch(branch, ref, force=True)
            overwrite_deps = True
        else:
            overwrite_deps = self.git.is_ancestor_of(branch, ref)
            self.merge(ref, create_dependency=False, download=False)
        if ticket is not None:
            old_dependencies = ", ".join(map(str,self._dependencies[branch]))
            if old_dependencies == "": old_dependencies = "(no dependencies)"

            trac_deps = self.trac.dependencies(ticket)
            if overwrite_deps:
                self._dependencies[branch] = trac_deps
            else:
                deps = set(trac_deps)
                git_deps = self._dependencies_as_tickets(branch)
                deps.update(git_deps)
                self._dependencies[branch] = tuple(sorted(deps))

            new_dependencies = ", ".join(map(str,self._dependencies[branch]))
            if new_dependencies == "": new_dependencies = "(no dependencies)"

            if old_dependencies != new_dependencies:
                self._UI.show("WARNING: the dependencies of this ticket have changed from %s to %s"%(old_dependencies, new_dependencies))

    def remote_status(self, ticket=None, quiet=False):
        r"""
        show the remote status of ``ticket``

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
        def show(lines):
            lines = [list(str(l) for l in line) if not isinstance(line, basestring) else line
                              for line in lines]
            tabulated_lines = [line for line in lines if not isinstance(line, basestring)]
            if tabulated_lines:
                column_widths = [max(len(x) for x in col) for col in zip(*tabulated_lines)]
            to_display = []
            for line in lines:
                if isinstance(line, basestring):
                    to_display.append(line)
                else:
                    for i in xrange(len(line)):
                        line[i] += ' '*(column_widths[i]-len(line[i]))
                    line.insert(3, 'behind')
                    line.insert(2, 'ahead')
                    to_display.append(' '.join(line))
            self._UI.show('\n'.join(to_display))

        if ticket is None :
            ticket = self.current_ticket()

        if isinstance(ticket, int):
            branch = self._branch[ticket]
        else:
            branch = ticket

        if ticket == 'all':
            ret = (self.remote_status(ticket or branch, quiet=True)
                    for ticket, branch in self.local_tickets(quiet=True))
            if quiet:
                return tuple(ret)
            else:
                show(ret)
                return
        try:
            remote_branch = self._remote_pull_branch(ticket)
            remote_ref = self._fetch(remote_branch)
        except (KeyError, RuntimeError):
            ret = '%s not tracked remotely' % ticket
            if quiet:
                return ret
            else:
                show((ret,))
                return
        ahead, behind = self.git.read_output("rev-list",
                "%s...%s"%(branch, remote_ref),
                left_right=True, count=True).split()
        behind = int(behind)
        ahead = int(ahead)
        ret = (ticket or '     ', remote_branch, ahead, behind)
        if quiet:
            return (ticket or '     ', remote_branch, ahead, behind)
        else:
            show((ret,))

    def import_patch(self, patchname=None, url=None, local_file=None,
            diff_format=None, header_format=None, path_format=None):
        r"""
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
        if not self.git.reset_to_clean_state(): return
        if not self.git.reset_to_clean_working_directory(): return

        if not local_file:
            return self.import_patch(
                    local_file=self.download_patch(
                        patchname=patchname, url=url),
                    diff_format=diff_format,
                    header_format=header_format,
                    path_format=path_format)
        elif patchname or url:
            raise ValueError("if local_file is specified, patchname "+
                             "and url must not be specified")
        else:
            lines = open(local_file).read().splitlines()
            lines = self._rewrite_patch(lines, to_header_format="git",
                    to_path_format="new", from_diff_format=diff_format,
                    from_header_format=header_format,
                    from_path_format=path_format)

            outfile = tempfile.mkstemp(dir=self.tmp_dir)[1]
            open(outfile, 'w').writelines("\n".join(lines)+"\n")

            self._UI.show("Trying to apply reformatted patch `%s` ..."%outfile)
            try:
                self.git.am(outfile, ignore_whitespace=True, resolvemsg='')
            except GitError:
                if not self._UI.confirm("The patch does not apply cleanly. "+
                                        "Would you like to apply it anyway "+
                                        "and create reject files for the "+
                                        "parts that do not apply?",
                                        default_no=True):
                    self._UI.show("Not applying patch.")
                    self.git.reset_to_clean_state(interactive=False)
                    return

                try:
                    self.git.apply(output, ignore_whitespace=True, reject=True)
                except GitError:
                    if self._UI.select("The patch did not apply "+
                        "cleanly. Please integrate the `.rej` files that "+
                        "were created and resolve conflicts. After you do, "+
                        "type `resolved`. If you want to abort this process, "+
                        "type `abort`.", ("resolved","abort")) == "abort":
                        self.git.reset_to_clean_state(interactive=False)
                        self.git.reset_to_clean_working_directory(interactive=False)
                        return

                self._UI.show("It seemed that the patch would not apply, "+
                              "but in fact it did.")

                self.git.add(update=True)
                self.git.am(resolved=True)

    def download_patch(self, ticketnum=None, patchname=None, url=None):
        r"""
        download a patch to a temporary directory

        If only ``ticketnum`` is specified and the ticket has only one
        attachment, download the patch attached to ``ticketnum``.

        If ``ticketnum`` and ``patchname`` are specified, download the
        patch ``patchname`` attached to ``ticketnum``.

        If ``url`` is specified, download ``url``.

        If nothing is specified, and if the ''current'' ticket has only
        one attachment, download it.

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
            fd, ret = tempfile.mkstemp(dir=self.tmp_dir)
            os.close(fd)
            check_call(["wget","-r","--no-check-certificate", "-O",ret,url])
            return ret
        elif ticketnum:
            if patchname:
                return self.download_patch(url = TRAC_SERVER_URI+"/raw-attachment/ticket/%s/%s"%(ticketnum,patchname))
            else:
                attachments = self.trac.attachment_names(ticketnum)
                if len(attachments) == 0:
                    raise ValueError("Ticket %s has no attachments."%self._ticket_repr(ticketnum))
                if len(attachments) == 1:
                    return self.download_patch(ticketnum = ticketnum, patchname = attachments[0])
                else:
                    raise ValueError("Ticket %s has more than one attachment but parameter `patchname` is not present."%self._ticket_repr(ticketnum))
        elif not patchname:
            return self.download_patch(ticketnum=self.current_ticket())
        else:
            raise ValueError("If `url` is not specified, `ticketnum` must be specified")

    def diff(self, base=None):
        r"""
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
        r"""
        Remove branches for tickets that are already merged into master.

        .. SEEALSO::

        - :meth:`abandon_ticket` -- Abandon a single ticket or branch.
        """
        for branch in self.git.local_branches():
            if self.git.is_ancestor_of(branch, MASTER_BRANCH):
                self._UI.show("Abandoning %s"%branch)
                self.git.abandon(branch)

    def abandon_ticket(self, ticket=None):
        r"""
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
        if self._UI.confirm("Are you sure you want to delete your work on %s?"%self._ticket_repr(ticketnum), default_no=True):
            self.git.abandon(ticketnum)

    def gather(self, branchname, *tickets, **kwds):
        r"""
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
        create_dependencies = kwds.pop('create_dependencies', True)
        download = kwds.pop('download', False)
        if len(tickets) == 0:
            raise ValueError("must include at least one input branch")
        if self.git.branch_exists(branchname):
            if not self._UI.confirm("The branch %s already "%branchname+
                                    "exists; do you want to merge into it?",
                                    default_no=True):
                return
            self.git.execute_supersilent("checkout", branchname)
        else:
            self.switch_ticket(tickets[0],
                    branchname=branchname)
            tickets = tickets[1:]
        for ticket in tickets:
            self.merge(
                    ticket,
                    message="Gathering %s into "%self._ticket_repr(ticket) +
                            "branch %s"%branchname,
                    **kwds)

    def show_dependencies(self, ticket=None, all=False, _seen=None): # all = recursive
        r"""
        show the dependencies of the given ticket

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
        if ticket is None:
            ticket = self.current_ticket()
        try:
            branchname = self._branch[ticket]
        except KeyError:
            raise ValueError("you must specify a valid ticket")
        if _seen is None:
            seen = []
        elif ticket in _seen:
            return
        else:
            seen = _seen
            seen.append(ticket)
        dep = self._dependencies_as_tickets(branchname)
        if not all:
            self._UI.show("Ticket %s depends on %s"%(ticket, ", ".join([str(d) for d in dep])))
        else:
            for d in dep:
                self.show_dependencies(d, True, seen)
            if _seen is None:
                self._UI.show("Ticket %s depends on %s"%(ticket, ", ".join([str(d) for d in seen])))

    def merge(self, ticket=MASTER_BRANCH, create_dependency=True, download=False, message=None):
        r"""
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
            for dep in self._dependencies[curbranch]:
                self.merge(dep, False, download, message)
            return
        elif ticket is None:
            raise ValueError("you must specify a ticket to merge")
        ref = dep = None
        if download:
            remote_branch = self._remote_pull_branch(ticket)
            if remote_branch is not None:
                ref = self._fetch(remote_branch)
                dep = ticket
            else:
                raise ValueError("could not download branch for ticket %s - its `branch` field on trac is empty or invalid")
        if ref is None:
            try:
                dep = ref = self._branch[ticket]
            except KeyError:
                pass # ticket does not exists locally but we were not asked to download it
        if ref is None:
            ref = ticket

        if create_dependency:
            if dep is None:
                dep = int(ticket)
            if dep and dep not in self._dependencies[curbranch]:
                self._dependencies[curbranch] += (dep,)
                self._UI.show("recorded dependency on %s"%dep)
        if message is None:
            kwds = {}
        else:
            kwds = {'m':message}

        if ref is None:
            self._UI.show("Nothing has been merged because the branch for %s could not be found. "%ticket+ ("Probably the branch field for the ticket is empty or invalid." if download else "Probably the branch for the ticket does not exist locally, consider using '--download True'"))
            return

        self.git.merge(ref, **kwds)

    def local_tickets(self, abandoned=False, quiet=False):
        r"""
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
        branch_info = [(b, self._ticket[b]) for b in raw_branches
            if b in self._ticket and (abandoned or not b.startswith("trash/"))]
        if quiet:
            return branch_info
        else:
            self._UI.show('\n'.join([
                        "{0}\t{1}".format(ticket or "     ", branch)
                        for (branch, ticket) in branch_info]))

    def current_ticket(self, error=False):
        r"""
        Returns the current ticket as an int, or ``None`` if there is
        no current ticket.

        INPUT:

        - ``error`` -- boolean (default ``False``), whether to raise
          an error if there is no current ticket

        .. SEEALSO::

        - :meth:`local_tickets` -- show all local tickets.
        """
        curbranch = self.git.current_branch()
        if curbranch is not None and curbranch in self._ticket:
            return self._ticket[curbranch]
        if error: raise ValueError("You must specify a ticket")

    def vanilla(self, release="release"):
        r"""
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
        r"""
        fetches ``branch`` from the remote repository, returning the name of
        the newly-updated local ref

        INPUT:

        - ``branch`` -- name of a remote branch

        - ``repo`` -- name of a remote repository

        OUTPUT:

        The name of a newly created/updated local ref.
        """
        if repository is None:
            repository = self.git._repo
        local_ref = "refs/remotes/trac/%s"%branch
        self.git.execute_supersilent('fetch', repository, "+%s:%s"%(branch, local_ref))
        return local_ref

    def _detect_patch_diff_format(self, lines):
        r"""
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
        r"""
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
            sage: dev._detect_patch_path_format(
            ....:     ["rename to sage/rings/number_field/totallyreal.pyx"], diff_format='hg')
            'old'
            sage: dev._detect_patch_path_format(
            ....:     ["rename from src/sage/rings/number_field/totalyreal.pyx"], diff_format='git')
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
            diff_regexs = (GIT_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX)
        elif diff_format == "hg":
            diff_regexs = (HG_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX)
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
        r"""
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

            sage: dev._rewrite_patch_diff_paths(
            ....:     ['rename from sage/combinat/crystals/letters.py',
            ....:      'rename to sage/combinat/crystals/letters.pyx'],
            ....:     to_format="new", diff_format="hg")
            ['rename from src/sage/combinat/crystals/letters.py',
             'rename to src/sage/combinat/crystals/letters.pyx']
            sage: dev._rewrite_patch_diff_paths(
            ....:     ['rename from src/sage/combinat/crystals/letters.py',
            ....:      'rename to src/sage/combinat/crystals/letters.pyx'],
            ....:     to_format="old", diff_format="git")
            ['rename from sage/combinat/crystals/letters.py',
             'rename to sage/combinat/crystals/letters.pyx']

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
            diff_regex = (HG_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX)
        elif diff_format == "git":
            diff_regex = (GIT_DIFF_REGEX, PM_DIFF_REGEX, MV_DIFF_REGEX)
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
        r"""
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
        r"""
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

    def _upload_ssh_key(self, keyfile, create_key_if_not_exists=True):
        r"""
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
            sage: os.unlink(tmp+'.pub')
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
        branch = self.trac._get_attributes(ticket).get('branch')
        if branch:
            return branch
        raise KeyError("branch field not set for ticket %s on trac"%ticket)

    def _remote_pull_branch(self, ticket):
        if isinstance(ticket, basestring):
            try:
                ticket = self._branch[ticket]
            except KeyError:
                return self._remote[ticket]
        if isinstance(ticket, int):
            return self._trac_branch(ticket)
        raise ValueError("ticket(={value}) must be instance of basesting of int, but is instance of {type}"
                         .format(value = ticket, type = type(ticket)))

    def _ticket_repr(self, ticket):
        if isinstance(ticket, basestring):
            ticket = self._ticket_to_branch(ticket)
            try:
                ticket = self._ticket[ticket]
            except KeyError:
                return str(ticket)
        if isinstance(ticket, int):
            return "#%s"%ticket
        raise ValueError

    def _dependencies_as_tickets(self, branch):
        dep = self._dependencies[branch]
        dep = [d if isinstance(d, int) else self._ticket[d] for d in dep]
        dep = [d for d in dep if d]
        return dep

    def _reset_to_clean_state(self, interactive=True):
        states = self.git.get_state()
        if not states:
            return True
        if (interactive and not
                self._UI.confirm("Your repository is in an unclean state. It "+
                                 "seems you are in the middle of a merge of "+
                                 "some sort. To run this command you have to "+
                                 "reset your respository to a clean state. "+
                                 "Do you want me to reset your respository? "+
                                 "(This will discard any changes which are "+
                                 "not commited.)")):
            return False

    def _reset_to_clean_working_directory(self, interactive=True):
        if not self.git.has_uncommitted_changes():
            return True

        if (interactive and not
                self._UI.confirm("You have uncommited changes in your "+
                                 "working directory. To run this command you "+
                                 "have to discard your changes. Do you want "+
                                 "me to discard any changes which are not "+
                                 "commited?")):
            return False

        return self.git._reset_to_clean_working_directory()

# unused method
#    def save(self, interactive=True):
#        r"""
#        guided command for making a commit which includes all changes
#
#        EXAMPLES::
#
#            sage: from sage.dev.sagedev import SageDev, doctest_config
#            sage: git = SageDev(doctest_config()).git
#            sage: git.save(False)
#            [first_branch ...] doctesting message
#             2 files changed, 2 insertions(+)
#             create mode 100644 untracked_testfile1
#             create mode 100644 untracked_testfile2
#        """
#        if (interactive and
#                self._UI.confirm("Would you like to see a diff of the "+
#                                 "changes?", default_yes=False)):
#            self.execute("diff")
#        for F in self.unknown_files():
#            if (not interactive or
#                    self._UI.confirm("Would you like to start tracking "+
#                                     "%s?"%F)):
#                self.execute('add', F)
#        if interactive:
#            msg = self._UI.get_input("Please enter a commit message:")
#        else:
#            msg = 'doctesting message'
#        self.commit_all(m=msg)
#

    def _save_uncommitted_changes(self):
        r"""
        Returns True if changes should be unstashed
        """
        if not self._UI.confirm("You have uncommitted changes, would you "+
                                "like to save them?"):
            return
        try:
            curbranch = self.git.current_branch()
            options = ["current branch", "new branch", "stash"]
        except ValueError:
            options = ["new branch", "stash"]
        dest = self._UI.select("Where do you want to store your changes?", options)
        if dest == "stash":
            self.stash()
        elif dest == "new branch":
            self.execute_silent("stash")
            return True
        elif dest == "current branch":
            self.commit()

