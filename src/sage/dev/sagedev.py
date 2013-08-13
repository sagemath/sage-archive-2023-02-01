r"""
SageDev

This module provides :class:`SageDev`, the central object of the developer
scripts for sage.

AUTHORS:

- David Roe, Frej Drejhammar, Julian Rueth, Martin Raum, Nicolas M. Thiery, R.
  Andrew Ohana, Robert Bradshaw, Timo Kluck: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 David Roe <roed.math@gmail.com>
#                          Frej Drejhammar <frej.drejhammar@gmail.com>
#                          Julian Rueth <julian.rueth@fsfe.org>
#                          Martin Raum <martin@raum-brothers.eu>
#                          Nicolas M. Thiery <Nicolas.Thiery@u-psud.fr>
#                          R. Andrew Ohana <andrew.ohana@gmail.com>
#                          Robert Bradshaw <robertwb@gmail.com>
#                          Timo Kluck <tkluck@infty.nl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from user_interface_error import OperationCancelledError
from trac_error import TracConnectionError, TracInternalError, TracError
from git_error import GitError

import re
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

# regular expression to check validity of git options
GIT_BRANCH_REGEX = re.compile(r'^(?!.*/\.)(?!.*\.\.)(?!/)(?!.*//)(?!.*@\{)(?!.*\\)[^\040\177 ~^:?*[]+(?<!\.lock)(?<!/)(?<!\.)$') # http://stackoverflow.com/questions/12093748/how-do-i-check-for-valid-git-branch-names

# the name of the branch which holds the vanilla clone of sage - in the long
# run this should be "master", currently, "public/sage-git/master" contains some changes
# over "master" which have not been reviewed yet but which are needed to work
# using git
MASTER_BRANCH = "public/sage-git/master"

COMMIT_GUIDE=r"""
# Please type your commit message above.
# Lines starting with '#' are ignored.
#
# An empty file aborts the commit.
"""

class SageDev(object):
    r"""
    The developer interface for sage.

    This class facilitates access to git and trac.

    INPUT:

    - ``config`` -- a :class:`config.Config` or ``None`` (default: ``None``),
      the configuration of this object; the defaults uses the configuration
      stored in ``DOT_GIT/devrc``.

    - ``UI`` -- a :class:`user_interface.UserInterface` or ``None`` (default:
      ``None``), the default creates a
      :class:`cmd_line_interface.CmdLineInterface` from ``config['UI']``.

    - ``trac`` -- a :class:`trac_interface.TracInterface` or ``None`` (default:
      ``None``), the default creates a :class:`trac_interface.TracInterface`
      from ``config['trac']``.

    - ``git`` -- a :class:`git_interface.GitInterface` or ``None`` (default:
      ``None``), the default creates a :class:`git_interface.GitInterface` from
      ``config['git']``.

    EXAMPLES::

        sage: dev._sagedev
        SageDev()

    """
    def __init__(self, config=None, UI=None, trac=None, git=None):
        r"""
        Initialization.

        TESTS::

            sage: type(dev._sagedev)
            <class 'sage.dev.sagedev.SageDev'>

        """
        self.config = config
        if self.config is None:
            from config import Config
            self.config = Config()

        # create some empty config sections if they do not yet exist
        for section in ['UI','trac','git','sagedev']:
            if section not in self.config:
                self.config[section] = {}

        self._UI = UI
        if self._UI is None:
            from cmd_line_interface import CmdLineInterface
            self._UI = CmdLineInterface(self.config['UI'])

        self.trac = trac
        if self.trac is None:
            from trac_interface import TracInterface
            self.trac = TracInterface(self.config['trac'], self._UI)

        self.git = git
        if self.git is None:
            from git_interface import GitInterface
            self.git = GitInterface(self.config['git'], self._UI)

        # create some SavingDicts to store the relations between branches and tickets
        from sage.env import DOT_SAGE
        import os
        ticket_file = self.config['sagedev'].get('ticketfile', os.path.join(DOT_SAGE, 'branch_to_ticket'))
        branch_file = self.config['sagedev'].get('branchfile', os.path.join(DOT_SAGE, 'ticket_to_branch'))
        dependencies_file = self.config['sagedev'].get('dependenciesfile', os.path.join(DOT_SAGE, 'dependencies'))
        remote_branches_file = self.config['sagedev'].get('remotebranchesfile', os.path.join(DOT_SAGE, 'remote_branches'))

        # some people dislike double underscore fields; here you can very
        # seriously screw up your setup if you put something invalid into
        # these. Ideally these fields should only be touched by single
        # underscore methods such as _set_remote_branch which do some checking
        # on the parameters
        from saving_dict import SavingDict
        self.__branch_to_ticket = SavingDict(ticket_file)
        self.__ticket_to_branch = SavingDict(branch_file, paired=self.__branch_to_ticket)
        self.__ticket_dependencies = SavingDict(dependencies_file, default=tuple)
        self.__branch_to_remote_branch = SavingDict(remote_branches_file)

    @property
    def tmp_dir(self):
        r"""
        A lazy property to provide a temporary directory

        TESTS::

            sage: import os
            sage: os.path.isdir(dev._sagedev.tmp_dir)
            True

        """
        try:
            return self._tmp_dir
        except AttributeError:
            import tempfile
            self._tmp_dir = tempfile.mkdtemp()
            import atexit, shutil
            atexit.register(shutil.rmtree, self._tmp_dir)
            return self._tmp_dir

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: dev # indirect doctest
            SageDev()

        """
        return "SageDev()"

    def create_ticket(self, branch=None, base=MASTER_BRANCH, remote_branch=None):
        r"""
        Create a new ticket on trac and switch to a new local branch to work on
        said ticket.

        INPUT:

        - ``branch`` -- a string or ``None`` (default: ``None``), the name of
          the local branch that will be used for the new ticket; if ``None``,
          the branch will be called ``'ticket/ticket_number'``.

        - ``base`` -- a string or ``None``, a branch on which to base the
          ticket (default: the master branch ``'master'``), or a ticket; if
          ``base`` is set to ``None``, then the current ticket is used. If
          ``base`` is a ticket, then the corresponding dependency will be
          added.

        - ``remote_branch`` -- a string or ``None`` (default: ``None``), the
          branch to pull from and push to on trac's git server; if ``None``,
          then the default branch ``'u/username/ticket_number'`` will be used.

        OUTPUT:

        Returns the number of the newly created ticket as an int.

        .. SEEALSO::

            :meth:`switch_ticket`, :meth:`download`, :meth:`edit_ticket`

        TESTS:

        Set up a single user environment::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: dev._wrap("_dependencies_for_ticket")
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Create some tickets::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            1
            sage: UI.append("Summary: ticket2\ndescription")
            sage: dev.create_ticket()
            2
            sage: dev.git.commit(SUPER_SILENT, allow_empty=True, message="second commit")
            sage: dev.git.commit_for_branch('ticket/2') != dev.git.commit_for_branch('ticket/1')
            True

        Check that ``base`` works::

            sage: UI.append("Summary: ticket3\ndescription")
            sage: dev.create_ticket(base=2)
            3
            sage: dev.git.commit_for_branch('ticket/3') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(3)
            (2,)
            sage: UI.append("Summary: ticket4\ndescription")
            sage: dev.create_ticket(base='ticket/2')
            4
            sage: dev.git.commit_for_branch('ticket/4') == dev.git.commit_for_branch('ticket/2')
            True
            sage: dev._dependencies_for_ticket(4)
            ()
            sage: UI.append("Summary: ticket5\ndescription")

        In this example ``base`` does not exist::

            sage: dev.create_ticket(base=1000)
            ValueError: `1000` is not a valid ticket name or ticket does not exist on trac.

        In this example ``base`` does not exist locally::

            sage: dev.trac.create_ticket("summary5","description",{})
            5
            sage: dev.create_ticket(base=5)
            ValueError: Branch field is not set for ticket #5 on trac.

        This also fails if the internet connection is broken::

            sage: dev.trac._connected = False
            sage: dev.create_ticket(base=4)
            A network error ocurred, ticket creation aborted.
            Your command failed because no connection to trac could be established.
            sage: dev.trac._connected = True

        Creating a ticket when in detached HEAD state::

            sage: dev.git.checkout(SUPER_SILENT, 'HEAD', detach=True)
            sage: UI.append("Summary: ticket detached\ndescription")
            sage: dev.create_ticket()
            6
            sage: dev.git.current_branch()
            'ticket/6'

        Creating a ticket when in the middle of a merge::

            sage: dev.git.checkout(SUPER_SILENT, '-b','merge_branch')
            sage: with open('merge', 'w') as f: f.write("version 0")
            sage: dev.git.add(SUPER_SILENT, 'merge')
            sage: dev.git.commit(SUPER_SILENT, '-m','some change')
            sage: dev.git.checkout(SUPER_SILENT, 'ticket/6')
            sage: with open('merge', 'w') as f: f.write("version 1")
            sage: dev.git.add(SUPER_SILENT, 'merge')
            sage: dev.git.commit(SUPER_SILENT, '-m','conflicting change')
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.merge(SUPER_SILENT, 'merge_branch')
            ....: except GitError: pass
            sage: UI.append("n")
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard many changes which are not commited.) [yes/No] n
            Could not switch to branch ticket/7 because your working directory is not clean.
            sage: dev.git.reset_to_clean_state()

        Creating a ticket with uncommitted changes::

            sage: open('tracked', 'w').close()
            sage: dev.git.add(SUPER_SILENT, 'tracked')
            sage: UI.append("keep")
            sage: UI.append("Summary: ticket merge\ndescription")
            sage: dev.create_ticket()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] keep
            Could not switch to branch ticket/8 because your working directory is not clean.

        """
        dependencies = []

        if branch is not None:
            self._check_local_branch_name(branch, exists=False)

        if base is None:
            base = self._current_ticket()
        if base is None:
            raise SageDevValueError("currently on no ticket, `base` must not be None")
        if self._is_ticket_name(base):
            base = self._ticket_from_ticket_name(base)
            dependencies.append(base)
            base = self._local_branch_for_ticket(base, download_if_not_found=True)
        self._check_local_branch_name(base, exists=True)

        if remote_branch is not None:
            self._check_remote_branch_name(remote_branch, exists=any)

        # now that we have checked that the parameters are valid, let the user
        # interactively create a ticket
        try:
            ticket = self.trac.create_ticket_interactive()
        except OperationCancelledError:
            self._UI.info("Ticket creation aborted.")
            raise
        except TracConnectionError as e:
            self._UI.error("A network error ocurred, ticket creation aborted.")
            raise

        # note that the dependencies are not recorded on the newly created
        # ticket but only stored locally - a first push to trac will set the
        # dependencies

        if branch is None:
            branch = self._new_local_branch_for_ticket(ticket)
        if remote_branch is None:
            remote_branch = self._remote_branch_for_ticket(ticket)

        # create a new branch for the ticket
        self.git.branch(branch, base)
        self._UI.info("Branch {0} created from branch {1}.".format(branch, base))
        try:
            self._set_local_branch_for_ticket(ticket, branch)
            if dependencies:
                self._set_dependencies_for_ticket(ticket, dependencies)
                self._UI.info("Dependencies {0} recorded locally for ticket #{1}.".format(", ".join(['#'+str(dep) for dep in dependencies]), ticket))
            self._set_remote_branch_for_branch(branch, remote_branch)
            self._UI.info("Branch {0} will pull from/push to remote branch {1}. Use {2} to set a different remote branch.".format(branch, remote_branch, self._format_command("set_remote", {"branch":branch, "remote":"remote_branch"})))
        except:
            self._UI.info("An error ocurred. Deleting branch {0}.".format(branch))

            self.git.branch(branch, delete=True)
            self._set_dependencies_for_ticket(ticket, None)
            self._set_remote_branch_for_branch(branch, None)
            self._set_local_branch_for_ticket(ticket, None)

            raise

        # switch to the new branch
        self._UI.info("Now switching to your new branch {0}.".format(branch))
        try:
            self.switch_ticket(ticket)
        except:
            self._UI.info("Your ticket has been created on trac. However, an error ocurred while switching to your new branch {0}. Use {1} to manually switch to {0}.".format(branch,self._format_command("switch_ticket",str(ticket))))
            raise

        return ticket

    def switch_ticket(self, ticket, branch=None):
        r"""
        Switch to a branch associated to ``ticket``.

        If ``branch`` is an existing local branch, then ``ticket`` will be
        associated to it, and the working directory will be switched to
        ``branch``.

        Otherwise, if there is no local branch for ``ticket``, the branch
        specified on trac will be downloaded to ``branch``. If the trac ticket
        does not specify a branch yet, then a new one will be created from
        "master".

        INPUT:

        - ``ticket`` -- a string or an integer identifying a ticket

        - ``branch`` -- a string, the name of the local branch that stores
          changes for ``ticket`` (default: ticket/``ticket``)

        .. SEEALSO::

            :meth:`download`, :meth:`create_ticket`, :meth:`vanilla`

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config_alice = DoctestConfig('alice')
            sage: config_alice['trac']['password'] = 'secret'
            sage: alice = DoctestSageDevWrapper(config_alice, server)
            sage: alice._pull_master_branch()

            sage: config_bob = DoctestConfig('bob')
            sage: config_bob['trac']['password'] = 'secret'
            sage: bob = DoctestSageDevWrapper(config_bob, server)
            sage: bob._pull_master_branch()

        Alice tries to switch to ticket #1 which does not exist yet::

            sage: alice._chdir()
            sage: alice.switch_ticket(1)
            ValueError: `1` is not a valid ticket name or ticket does not exist on trac.

        Bob creates that ticket::

            sage: bob._chdir()
            sage: bob._UI.append("Summary: summary1\ndescription")
            sage: bob.create_ticket()
            1

        Now alice can switch to it, even though there is no branch on the
        ticket description::

            sage: alice._chdir()
            sage: alice.switch_ticket(1)

        If Bob commits something to the ticket, a ``switch_ticket`` by Alice
        does not take his changes into account::

            sage: bob._chdir()
            sage: bob.git.commit(SUPER_SILENT, allow_empty=True,message="empty commit")
            sage: bob._UI.append("y")
            sage: bob.upload()
            The branch u/bob/ticket/1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y

            sage: alice._chdir()
            sage: alice.switch_ticket(1)
            sage: alice.git.log('--pretty=%s')
            initial commit

        If Alice had not switched to that ticket before, she would of course
        see Bob's changes (this also checks that we can handle a corrupt ticket
        database and a detached HEAD)::

            sage: alice.git.checkout(SUPER_SILENT, 'HEAD', detach=True)
            sage: alice.git.branch(SUPER_SILENT, '-d','ticket/1')
            sage: alice.switch_ticket(1) # ticket #1 refers to the non-existant branch 'ticket/1'
            Ticket #1 refers to the non-existant local branch ticket/1. If you have not manually interacted with git, then this is a bug in sagedev. Removing the association from ticket #1 to branch ticket/1.
            sage: alice.git.current_branch()
            'ticket/1'
            sage: alice.git.log('--pretty=%s')
            empty commit
            initial commit

        Switching to a ticket with untracked files::

            sage: alice._UI.append("Summary: summary2\ndescription")
            sage: alice.create_ticket()
            2
            sage: alice.git.log('--pretty=%s')
            initial commit
            sage: open("untracked","w").close()
            sage: alice.switch_ticket(1)
            sage: alice.git.log('--pretty=%s')
            empty commit
            initial commit

        Switching to a ticket with untracked files which make a switch
        impossible::

            sage: alice.git.add(SUPER_SILENT, "untracked")
            sage: alice.git.commit(SUPER_SILENT, message="added untracked")
            sage: alice.switch_ticket(2)
            sage: open("untracked","w").close()
            sage: alice.switch_ticket(1)
            GitError: git exited with a non-zero exit code (1).
            This happened while executing `git checkout ticket/1`.
            git printed nothing to STDOUT.
            git printed the following to STDERR:
            error: The following untracked working tree files would be overwritten by checkout:
                untracked
            Please move or remove them before you can switch branches.
            Aborting

        Switching to a ticket with uncommited changes::

            sage: open("tracked","w").close()
            sage: alice.git.add(SUPER_SILENT, "tracked")
            sage: alice._UI.append('d')
            sage: alice.switch_ticket(2)
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] d

        """
        self._check_ticket_name(ticket, exists=True)

        ticket = self._ticket_from_ticket_name(ticket)

        if branch is None:
            if self._has_local_branch_for_ticket(ticket):
                branch = self._local_branch_for_ticket(ticket)
                self._UI.info("Switching to branch {0}.".format(branch))
                self.switch_branch(branch)
                return
            else:
                branch = self._new_local_branch_for_ticket(ticket)
                self._check_local_branch_name(branch, exists=False)

        self._check_local_branch_name(branch)

        if self._is_local_branch_name(branch, exists=True):
            # reset ticket to point to branch and checkout
            self._set_local_branch_for_ticket(ticket, branch)
            self._UI.info("Set local branch for ticket #{0} to {1}.".format(ticket, branch))
            self.switch_ticket(ticket, branch=None)
            return

        remote_branch = self.trac._branch_for_ticket(ticket)
        dependencies = self.trac.dependencies(ticket)
        if remote_branch is None: # branch field is not set on ticket
            self._UI.info("The branch field on ticket #{0} is not set. Creating a new branch {1} off the master branch {2}.".format(ticket, branch, MASTER_BRANCH))
            self.git.branch(branch, MASTER_BRANCH)
        else:
            try:
                self.download(remote_branch, branch)
            except:
                self._UI.error("Could not switch to ticket #{0} because the remote branch `{1}` for that ticket could not be downloaded.".format(ticket, remote_branch))
                raise

        try:
            self._set_local_branch_for_ticket(ticket, branch)
            self._set_dependencies_for_ticket(ticket, dependencies)
        except:
            self._UI.info("An error ocurred. Deleting branch {0}.".format(branch))
            self._set_local_branch_for_ticket(ticket, None)
            self.git.branch("-d",branch)
            raise

        self._UI.info("Switching to newly created branch {0}.".format(branch))
        self.switch_branch(branch)

    def switch_branch(self, branch):
        r"""
        Switch to the local branch ``branch``.

        INPUT:

        - ``branch`` - a string, the name of a local branch

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Create a few branches::

            sage: dev.git.branch("branch1")
            sage: dev.git.branch("branch2")

        Switch to a branch::

            sage: dev.switch_branch("branch1")
            sage: dev.git.current_branch()
            'branch1'

        The branch must exist::

            sage: dev.switch_branch("branch3")
            ValueError: Branch `branch3` does not exist locally.

        Switching branches with untracked files::

            sage: open("untracked","w").close()
            sage: dev.switch_branch("branch2")

        Switching branches with uncommitted changes::

            sage: open("tracked","w").close()
            sage: dev.git.add(SUPER_SILENT, "tracked")
            sage: dev.git.commit(SUPER_SILENT, message="added tracked")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("keep")
            sage: dev.switch_branch("branch1")
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] keep
            Could not switch to branch branch1 because your working directory is not clean.

        We can stash uncommitted changes::

            sage: UI.append("s")
            sage: dev.switch_branch("branch1")
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] s
            Your changes have been recorded on a new branch `stash/1`.

        And unstash the changes later::

            sage: dev.switch_branch('branch2')
            sage: dev.unstash()
            stash/1
            sage: dev.unstash('stash/1')

        Or we can just discard the changes::

            sage: UI.append("d")
            sage: dev.switch_branch("branch1")
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] d

        Switching branches when in the middle of a merge::

            sage: dev.git.checkout(SUPER_SILENT, '-b','merge_branch')
            sage: with open('merge', 'w') as f: f.write("version 0")
            sage: dev.git.add(SUPER_SILENT, 'merge')
            sage: dev.git.commit(SUPER_SILENT, '-m','some change')
            sage: dev.git.checkout(SUPER_SILENT, 'branch1')
            sage: with open('merge', 'w') as f: f.write("version 1")
            sage: dev.git.add(SUPER_SILENT, 'merge')
            sage: dev.git.commit(SUPER_SILENT, '-m','conflicting change')
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.merge(SUPER_SILENT, 'merge_branch')
            ....: except GitError: pass
            sage: UI.append('n')
            sage: dev.switch_branch('merge_branch')
            Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard many changes which are not commited.) [yes/No] n
            Could not switch to branch merge_branch because your working directory is not clean.
            sage: dev.git.reset_to_clean_state()

        Switching branches when in a detached HEAD::

            sage: dev.git.checkout(SUPER_SILENT, 'branch2', detach=True)
            sage: dev.switch_branch('branch1')

        With uncommitted changes::

            sage: dev.git.checkout(SUPER_SILENT, 'branch2', detach=True)
            sage: with open('tracked', 'w') as f: f.write("boo")
            sage: UI.append("discard")
            sage: dev.switch_branch('branch1')
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] discard

        Switching branches with untracked files that would be overwritten by
        the switch::

            sage: with open('tracked', 'w') as f: f.write("boo")
            sage: dev.switch_branch('branch2')
            GitError: git exited with a non-zero exit code (1).
            This happened while executing `git checkout branch2`.
            git printed nothing to STDOUT.
            git printed the following to STDERR:
            error: The following untracked working tree files would be overwritten by checkout:
                tracked
            Please move or remove them before you can switch branches.
            Aborting

        """
        self._check_local_branch_name(branch, exists=True)

        try:
            self.reset_to_clean_state()
            self.reset_to_clean_working_directory()
        except OperationCancelledError:
            self._UI.error("Could not switch to branch {0} because your working directory is not clean.".format(branch))
            raise

        try:
            from sage.dev.git_interface import SUPER_SILENT
            self.git.checkout(SUPER_SILENT, branch)
        except GitError as e:
            # the error message should be self explanatory
            raise

    def download(self, ticket_or_branch=None, branch=None):
        r"""
        Download ``ticket_or_branch`` to ``branch``.

        INPUT:

        - ``ticket_or_branch`` -- a string or an integer or ``None`` (default:
          ``None``), a ticket or a remote branch name; setting this to ``None``
          has the same effect as setting it to the :meth:`current_ticket`.

        - ``branch`` -- a string or ``None`` (default: ``None``), the branch to
          create or merge the changes into. If ``None``, then a new branch will
          be created unless there is already a branch for this ticket.

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config_alice = DoctestConfig('alice')
            sage: config_alice['trac']['password'] = 'secret'
            sage: alice = DoctestSageDevWrapper(config_alice, server)
            sage: alice._pull_master_branch()

            sage: config_bob = DoctestConfig('bob')
            sage: config_bob['trac']['password'] = 'secret'
            sage: bob = DoctestSageDevWrapper(config_bob, server)
            sage: bob._pull_master_branch()

        Alice creates ticket 1::

            sage: alice._chdir()
            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: ticket = alice.create_ticket()

        Bob attempts to download the ticket but fails because there is no
        branch for the ticket yet::

            sage: bob._chdir()
            sage: bob.download(ticket)
            ValueError: Branch field is not set for ticket #1 on trac.

        So, Bob starts to work on the ticket on a new branch::

            sage: bob.switch_ticket(ticket)

        Alice pushes a commit::

            sage: alice._chdir()
            sage: alice.git.commit(SUPER_SILENT, allow_empty=True, message="alice: empty commit")
            sage: alice._UI.append("y")
            sage: alice.upload()
            The branch u/alice/ticket/1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y

        Bob downloads the changes for ticket 1::

            sage: bob._chdir()
            sage: bob.download()
            sage: bob.git.log('--pretty=%s')
            alice: empty commit
            initial commit

        Bob commits a change::

            sage: open("bobs_file","w").close()
            sage: bob.git.add("bobs_file")
            sage: bob.git.commit(SUPER_SILENT, message="bob: added bobs_file")
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.upload()
            The branch u/bob/ticket/1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y
            I will now change the branch field of ticket #1 from its current value `u/alice/ticket/1` to `u/bob/ticket/1`. Is this what you want? [Yes/no] y

        Alice commits non-conflicting changes::

            sage: alice._chdir()
            sage: with open("alices_file","w") as f: f.write("1")
            sage: alice.git.add("alices_file")
            sage: alice.git.commit(SUPER_SILENT, message="alice: added alices_file")

        Alice can now download the changes by Bob without the need to merge
        manually::

            sage: alice.download()
            sage: alice.git.log('--pretty=%s')
            Merge branch 'u/bob/ticket/1' of /dev/shm/... into ticket/1
            alice: added alices_file
            bob: added bobs_file
            alice: empty commit
            initial commit

        Now, Bob commits some conflicting changes::

            sage: bob._chdir()
            sage: with open("alices_file","w") as f: f.write("2")
            sage: bob.git.add("alices_file")
            sage: bob.git.commit(SUPER_SILENT, message="bob: added alices_file")
            sage: bob._UI.append('y')
            sage: bob.upload()
            I will now upload the following new commits to the remote branch `u/bob/ticket/1`:
            ...: bob: added alices_file
            Is this what you want? [Yes/no] y

        Now, the download fails; one would have to use :meth:`merge`::

            sage: alice._chdir()
            sage: alice.download()
            GitError: git exited with a non-zero exit code (1).
            Pulling `u/bob/ticket/1` into `ticket/1` failed. Most probably this happened because this did not resolve as a fast-forward, i.e., there were conflicting changes. Maybe there are untracked files in your working directory which made the pull impossible.

        Undo the latest commit by alice, so we can download again::

            sage: alice.git.reset(SUPER_SILENT, 'HEAD~~', hard=True)
            sage: alice.download()
            sage: alice.git.log('--pretty=%s')
            bob: added alices_file
            bob: added bobs_file
            alice: empty commit
            initial commit

        Now, Alice creates an untracked file which makes a trivial merge
        impossible::

            sage: alice._chdir()
            sage: open("bobs_other_file","w").close()

            sage: bob._chdir()
            sage: open("bobs_other_file","w").close()
            sage: bob.git.add(SUPER_SILENT, "bobs_other_file")
            sage: bob.git.commit(SUPER_SILENT, message="bob: added bobs_other_file")
            sage: bob._UI.append('y')
            sage: bob.upload()
            I will now upload the following new commits to the remote branch `u/bob/ticket/1`:
            ...: bob: added bobs_other_file
            Is this what you want? [Yes/no] y

            sage: alice._chdir()
            sage: alice.download()
            GitError: git exited with a non-zero exit code (1).
            Pulling `u/bob/ticket/1` into `ticket/1` failed. Most probably this happened because this did not resolve as a fast-forward, i.e., there were conflicting changes. Maybe there are untracked files in your working directory which made the pull impossible.

        """
        if ticket_or_branch is None:
            ticket_or_branch = self._current_ticket()
            branch = self.git.current_branch()

        if ticket_or_branch is None:
            raise SageDevValueError("No `ticket_or_branch` specified to download.")

        if self._is_ticket_name(ticket_or_branch):
            ticket = self._ticket_from_ticket_name(ticket_or_branch)
            self._check_ticket_name(ticket, exists=True)

            remote_branch = self.trac._branch_for_ticket(ticket)
            if remote_branch is None:
                raise SageDevValueError("Branch field is not set for ticket #{0} on trac.".format(ticket))
            if branch is None:
                branch = self._new_local_branch_for_ticket(ticket)
            self._check_local_branch_name(branch)

        else:
            remote_branch = ticket_or_branch
            self._check_remote_branch_name(remote_branch)

            if branch is None:
                branch = remote_branch
            self._check_local_branch_name(branch)

        self._check_remote_branch_name(remote_branch, exists=True)

        self._UI.info("Fetching remote branch {0} into {1}.".format(remote_branch, branch))
        from git_error import DetachedHeadError
        try:
            current_branch = self.git.current_branch()
        except DetachedHeadError:
            current_branch = None

        if current_branch == branch:
            # we cannot fetch onto the current branch - we have to pull
            self.reset_to_clean_state()
            self.reset_to_clean_working_directory()

            try:
                from sage.dev.git_interface import SUPER_SILENT
                self.git.pull(SUPER_SILENT, self.git._repository, remote_branch)
            except GitError as e:
                # this might fail because the pull did not resolve as a
                # fast-forward or because there were untracked files around
                # that made a pull impossible
                # is there a way to find out?
                e.explain = "Pulling `{0}` into `{1}` failed. Most probably this happened because this did not resolve as a fast-forward, i.e., there were conflicting changes. Maybe there are untracked files in your working directory which made the pull impossible.".format(remote_branch, branch)
                e.advice =  "You can try to use {0} to resolve any conflicts manually.".format(self._format_command("merge",{"remote_branch":remote_branch}))
                raise
        else:
            try:
                from sage.dev.git_interface import SUPER_SILENT
                self.git.fetch(SUPER_SILENT, self.git._repository, "{0}:{1}".format(remote_branch, branch))
            except GitError as e:
                # there is not many scenarios in which this can fail - the most
                # likely being that branch already exists and this does not
                # resolve as a fast-forward; in any case, if the fetch fails,
                # then just nothing happened and we can abort the download
                # safely without a need to cleanup
                e.explain = "Fetching `{0}` into `{1}` failed.".format(remote_branch, branch)
                if self._is_local_branch_name(branch, exists=True):
                    e.explain += " Most probably this happened because the fetch did not resolve as a fast-forward, i.e., there were conflicting changes."
                    e.advice = "You can try to use {2} to switch to {1} and then use {3} to resolve these conflicts manually.".format(remote_branch, branch, self._format_command("switch-branch",branch), self._format_command("merge",{"remote_branch":remote_branch}))
                else:
                    # is there any advice one could give to the user?
                    pass
                raise

    def commit(self, message=None, interactive=False):
        r"""
        Create a commit from the pending changes on the current branch.

        This is most akin to mercurial's commit command, not git's.

        INPUT:

        - ``message`` -- the message of the commit (default: ``None``), if
          ``None``, prompt for a message.

        - ``interactive`` -- if set, interactively select which part of the
          changes should be part of the commit

        .. SEEALSO::

        - :meth:`upload` -- Upload changes to the remote server.  This
          is the next step once you've committed some changes.

        - :meth:`diff` -- Show changes that will be committed.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Commit an untracked file::

            sage: dev.git.checkout('-b','branch1')
            Switched to a new branch 'branch1'

            sage: open("tracked","w").close()
            sage: dev._UI.extend(["added tracked","y","y","y"])
            sage: dev.commit()
            The following files in your working directory are not tracked by git:
             tracked
            Do you want to add any of these files in this commit? [yes/No] y
            Do you want to add `tracked`? [yes/No] y
            Do you want to commit your changes to branch `branch1`? I will prompt you for a commit message if you do. [Yes/no] y
            [branch1 ...] added tracked
             0 files changed
             create mode ... tracked

        Commit a tracked file::

            sage: with open("tracked","w") as F: F.write("foo")
            sage: dev._UI.extend(["modified tracked","y"])
            sage: dev.commit()
            Do you want to commit your changes to branch `branch1`? I will prompt you for a commit message if you do. [Yes/no] y
            [branch1 ...] modified tracked
             1 file changed, 1 insertion(+)

        """
        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            self._UI.error("Cannot commit changes when not on any branch.")
            self._UI.info("Use {0} or {1} to switch to a branch.".format(self._format_command("switch_branch"), self._format_command("switch_ticket")))
            raise OperationCancelledError("cannot proceed in detached HEAD mode")

        # make sure the index is clean
        from git_interface import SUPER_SILENT
        self.git.reset(SUPER_SILENT)

        try:
            self._UI.info("Committing pending changes to branch `{0}`.".format(branch))

            try:
                untracked_files = self.git.untracked_files()
                if untracked_files:
                    if self._UI.confirm("The following files in your working directory are not tracked by git:\n{0}\nDo you want to add any of these files in this commit?".format("\n".join([" "+fname for fname in untracked_files])), default=False):
                        for file in untracked_files:
                            if self._UI.confirm("Do you want to add `{0}`?".format(file), default=False):
                                self.git.add(file)

                if interactive:
                    self.git.add(patch=True)
                else:
                    self.git.add(update=True)

                if not self._UI.confirm("Do you want to commit your changes to branch `{0}`? I will prompt you for a commit message if you do.".format(branch), default=True):
                    self._UI.info("If you want to commit to a different branch/ticket, run {0} or {1} first.".format(self._format_command("switch_branch"), self._format_command("switch_ticket")))
                    raise OperationCancelledError("user does not want to create a commit")

                from tempfile import NamedTemporaryFile
                commit_message = NamedTemporaryFile()
                if message: commit_message.write(message)
                commit_message.write(COMMIT_GUIDE)
                commit_message.close()

                self._UI.edit(commit_message.name)

                message = "\n".join([line for line in open(commit_message.name).read().splitlines() if not line.startswith("#")]).strip()

                if not message:
                    raise OperationCancelledError("empty commit message")

                self.git.commit(message=message)

            except OperationCancelledError:
                self._UI.info("Not creating a commit.")
                raise

        finally:
            # do not leave a non-clean index behind
            self.git.reset(SUPER_SILENT)

    def set_remote(self, branch_or_ticket, remote_branch):
        r"""
        Set the remote branch to push to for ``branch_or_ticket`` to
        ``remote_branch``.

        INPUT:

        - ``branch_or_ticket`` -- a string, the name of a local branch, or a
          string or an integer identifying a ticket or ``None``; if ``None``,
          the current branch is used.

        - ``remote_branch`` -- a string, the name of a remote branch (this
          branch may not exist yet)

        .. SEEALSO::

        - :meth:`upload` -- To upload changes after setting the remote branch

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: dev._wrap("_remote_branch_for_ticket")
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Create a new branch::

            sage: UI.append("Summary: ticket1\ndescription")
            sage: dev.create_ticket()
            1

        Modify the remote branch for this ticket's branch::

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev.set_remote('ticket/1', 'u/doctest/foo')
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/foo'
            sage: dev.set_remote('ticket/1', 'foo')
            The remote branch `foo` is not in your user scope. You might not have permission to push to that branch. Did you mean to set the remote branch to `u/doctest/foo`?
            sage: dev._remote_branch_for_ticket(1)
            'foo'
            sage: dev.set_remote('#1', 'u/doctest/foo')
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/foo'

        """
        if branch_or_ticket is None:
            from git_error import DetachedHeadError
            try:
                branch = self.git.current_branch()
            except DetachedHeadError:
                self._UI.error("`branch` must not be None because you are in detached HEAD state.")
                self._UI.info("Switch to a branch with `{0}` or specify branch explicitly.".format(self._format_command('switch_branch')))
                raise OperationCancelledError("detached head state")
        elif self._is_ticket_name(branch_or_ticket):
            ticket = self._ticket_from_ticket_name(branch_or_ticket)
            if not self._has_local_branch_for_ticket(ticket):
                self._UI.error("no local branch for ticket #{0} found. Cannot set remote branch for that ticket.".format(ticket))
                raise OperationCancelledError("no such ticket")
            branch = self._local_branch_for_ticket(ticket)
        else:
            branch = branch_or_ticket

        self._check_local_branch_name(branch, exists=True)
        self._check_remote_branch_name(remote_branch)

        if not remote_branch.startswith('u/{0}/'.format(self.trac._username)):
            self._UI.warning("The remote branch `{0}` is not in your user scope. You might not have permission to push to that branch. Did you mean to set the remote branch to `u/{1}/{0}`?".format(remote_branch, self.trac._username))

        self._set_remote_branch_for_branch(branch, remote_branch)

    def upload(self, ticket=None, remote_branch=None, force=False):
        r"""
        Upload the current branch to the Sage repository.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), if ``None`` and currently working on a ticket or
          if ``ticket`` specifies a ticket, then the branch on that ticket is
          set to ``remote_branch`` after the current branch has been uploaded there.

        - ``remote_branch`` -- a string or ``None`` (default: ``None``), the remote
          branch to upload to; if ``None``, then a default is chosen

        - ``force`` -- a boolean (default: ``False``), whether to upload if
          this is not a fast-forward.

        .. SEEALSO::

        - :meth:`commit` -- Save changes to the local repository.

        - :meth:`download` -- Update a ticket with changes from the remote
          repository.

        TESTS::

        Create a doctest setup with two users::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config_alice = DoctestConfig('alice')
            sage: config_alice['trac']['password'] = 'secret'
            sage: alice = DoctestSageDevWrapper(config_alice, server)
            sage: alice._pull_master_branch()

            sage: config_bob = DoctestConfig('bob')
            sage: config_bob['trac']['password'] = 'secret'
            sage: bob = DoctestSageDevWrapper(config_bob, server)
            sage: bob._pull_master_branch()

        Alice tries to upload to ticket 1 which does not exist yet::

            sage: alice._chdir()
            sage: alice.upload(ticket=1)
            ValueError: `1` is not a valid ticket name or ticket does not exist on trac.

        Alice creates ticket 1 and uploads some changes to it::

            sage: alice._UI.append("Summary: summary1\ndescription")
            sage: ticket = alice.create_ticket()
            sage: open("tracked", "w").close()
            sage: alice.git.add(SUPER_SILENT, "tracked")
            sage: alice.git.commit(SUPER_SILENT, message="alice: added tracked")
            sage: alice._UI.append("y")
            sage: alice.upload()
            The branch u/alice/ticket/1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y

        Now Bob can switch to that ticket and upload changes himself::

            sage: bob._chdir()
            sage: bob.switch_ticket(1)
            sage: with open("tracked", "w") as f: f.write("bob")
            sage: bob.git.add(SUPER_SILENT, "tracked")
            sage: bob.git.commit(SUPER_SILENT, message="bob: modified tracked")
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.upload()
            The branch u/bob/ticket/1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y
            I will now change the branch field of ticket #1 from its current value `u/alice/ticket/1` to `u/bob/ticket/1`. Is this what you want? [Yes/no] y

        Now Alice can download these changes::

            sage: alice._chdir()
            sage: alice.download()

        Alice and Bob make non-conflicting changes simultaneously::

            sage: with open("tracked", "w") as f: f.write("alice")
            sage: alice.git.add(SUPER_SILENT, "tracked")
            sage: alice.git.commit(SUPER_SILENT, message="alice: modified tracked")

            sage: bob._chdir()
            sage: open("tracked2", "w").close()
            sage: bob.git.add(SUPER_SILENT, "tracked2")
            sage: bob.git.commit(SUPER_SILENT, message="bob: added tracked2")

        After Alice uploaded her changes, Bob can not set the branch field anymore::

            sage: alice._chdir()
            sage: alice._UI.append("y")
            sage: alice._UI.append("y")
            sage: alice.upload()
            I will now upload the following new commits to the remote branch `u/alice/ticket/1`:
            ...: alice: modified tracked
            ...: bob: modified tracked
            Is this what you want? [Yes/no] y
            I will now change the branch field of ticket #1 from its current value `u/bob/ticket/1` to `u/alice/ticket/1`. Is this what you want? [Yes/no] y

            sage: bob._chdir()
            sage: bob._UI.append("y")
            sage: bob.upload()
            I will now upload the following new commits to the remote branch `u/bob/ticket/1`:
            ...: bob: added tracked2
            Is this what you want? [Yes/no] y
            Not setting the branch field for ticket #1 to `u/bob/ticket/1` because `u/bob/ticket/1` and the current value of the branch field `u/alice/ticket/1` have diverged.

        After merging the changes, this works again::

            sage: bob.download()
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.upload()
            I will now upload the following new commits to the remote branch `u/bob/ticket/1`:
            ...: Merge branch 'u/alice/ticket/1' of ... into ticket/1
            ...: alice: modified tracked
            Is this what you want? [Yes/no] y
            I will now change the branch field of ticket #1 from its current value `u/alice/ticket/1` to `u/bob/ticket/1`. Is this what you want? [Yes/no] y

        Check that ``ticket`` works::

            sage: bob.upload(2)
            ValueError: `2` is not a valid ticket name or ticket does not exist on trac.

        After creating the ticket, this works with a warning::

            sage: bob._UI.append("Summary: summary2\ndescription")
            sage: bob.create_ticket()
            2
            sage: bob.switch_ticket(1)
            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.upload(2)
            You are trying to push the branch `ticket/1` to `u/bob/ticket/2` for ticket #2. However, your local branch for ticket #2 seems to be `ticket/2`. Do you really want to proceed? [yes/No] y
            The branch u/bob/ticket/2 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y

        Check that ``remote_branch`` works::

            sage: bob._UI.append("y")
            sage: bob._UI.append("y")
            sage: bob.upload(remote_branch="u/bob/branch1")
            The branch u/bob/branch1 does not exist on the remote server yet. Do you want to create the branch? [Yes/no] y
            I will now change the branch field of ticket #1 from its current value `u/bob/ticket/1` to `u/bob/branch1`. Is this what you want? [Yes/no] y

        """
        from git_interface import SUPER_SILENT

        if ticket is None:
            ticket = self._current_ticket()
        if ticket is not None:
            ticket = self._ticket_from_ticket_name(ticket)
            self._check_ticket_name(ticket, exists=True)

        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            self._UI.error("Cannot upload while in detached HEAD state.")
            raise OperationCancelledError("cannot upload while in detached HEAD state")

        if remote_branch is None:
            if ticket:
                remote_branch = self._remote_branch_for_ticket(ticket)
                if remote_branch is None:
                    raise SageDevValueError("remote_branch must be specified since #{0} has no remote branch set.".format(ticket))
            else:
                remote_branch = self._remote_branch_for_branch(branch)
                if remote_branch is None:
                    raise SageDevValueError("remote_branch must be specified since the current branch has no remote branch set.")

        self._check_remote_branch_name(remote_branch)

        # whether the user already confirmed that he really wants to push and set the branch field
        user_confirmation = force

        if ticket is not None:
            if self._has_local_branch_for_ticket(ticket) and self._local_branch_for_ticket(ticket) == branch:
                pass
            elif self._has_local_branch_for_ticket(ticket) and self._local_branch_for_ticket(ticket) != branch:
                if user_confirmation or self._UI.confirm("You are trying to push the branch `{0}` to `{1}` for ticket #{2}. However, your local branch for ticket #{2} seems to be `{3}`. Do you really want to proceed?".format(branch, remote_branch, ticket, self._local_branch_for_ticket(ticket)), default=False):
                    self._UI.info("To permanently set the branch associated to ticket #{0} to `{1}`, use `{2}`.".format(ticket, branch, self._format_command("switch_ticket",ticket=ticket,branch=branch)))
                    user_confirmation = True
                else:
                    raise OperationCancelledError("user requsted")
            elif self._has_ticket_for_local_branch(branch) and self._ticket_for_local_branch(branch) != ticket:
                if user_confirmation or self._UI.confirm("You are trying to push the branch `{0}` to `{1}` for ticket #{2}. However, that branch is associated to ticket #{3}. Do you really want to proceed?".format(branch, remote_branch, ticket, self._ticket_for_local_branch(branch))):
                    self._UI.info("To permanently set the branch associated to ticket #{0} to `{1}`, use `{2}`. To create a new branch from `{1}` for #{0}, use `{3}` and `{4}`.".format(ticket, branch, self._format_command("switch_ticket",ticket=ticket,branch=branch), self._format_command("switch_ticket",ticket=ticket), self._format_command("merge", branch=branch)))
                    user_confirmation = True

        self._UI.info("Uploading your changes in `{0}` to `{1}`.".format(branch, remote_branch))
        try:
            remote_branch_exists = self._is_remote_branch_name(remote_branch, exists=True)
            if not remote_branch_exists:
                if not self._UI.confirm("The branch {0} does not exist on the remote server yet. Do you want to create the branch?".format(remote_branch), default=True):
                    raise OperationCancelledError("User did not want to create remote branch.")
            else:
                self.git.fetch(SUPER_SILENT, self.git._repository, remote_branch)

            # check whether force is necessary
            if remote_branch_exists and not self.git.is_child_of(branch, 'FETCH_HEAD'):
                if not force:
                    self._UI.error("Not uploading your changes because they would discard some of the commits on the remote branch `{0}`.".format(remote_branch))
                    self._UI.info("If this is really what you want, use `{0}` to upload your changes.".format(remote_branch, self._format_command("upload",ticket=ticket,remote_branch=remote_branch,force=True)))
                    raise OperationCancelledError("not a fast-forward")

            # check whether this is a nop
            if remote_branch_exists and not force and self.git.commit_for_branch(branch) == self.git.commit_for_branch('FETCH_HEAD'):
                self._UI.info("Not uploading your changes because the remote branch `{0}` is idential to your local branch `{1}`. Did you forget to commit your changes with `{2}`?".format(remote_branch, branch, self._format_command("commit")))
            else:
                try:
                    if not force:
                        if remote_branch_exists:
                            from git_interface import READ_OUTPUT
                            commits = self.git.log(READ_OUTPUT, "{0}..{1}".format('FETCH_HEAD', branch), '--pretty=%h: %s')
                            if not self._UI.confirm("I will now upload the following new commits to the remote branch `{0}`:\n{1}Is this what you want?".format(remote_branch, commits), default=True):
                                raise OperationCancelledError("user requested")

                    self.git.push(SUPER_SILENT, self.git._repository, "{0}:{1}".format(branch, remote_branch), force=force)
                except GitError as e:
                    # can we give any advice if this fails?
                    raise

            self._UI.info("Your changes in `{0}` have been uploaded to `{1}`.".format(branch, remote_branch))

        except OperationCancelledError:
            self._UI.info("Did not upload any changes.")
            raise


        if ticket:
            current_remote_branch = self.trac._branch_for_ticket(ticket)
            if current_remote_branch == remote_branch:
                self._UI.info("Not setting the branch field for ticket #{0} because it already points to your branch `{1}`.".format(ticket, remote_branch))
            else:
                self._UI.info("Setting the branch field of ticket #{0} to `{1}`.".format(ticket, remote_branch))

                if current_remote_branch is not None:
                    self.git.fetch(SUPER_SILENT, self.git._repository, current_remote_branch)
                    if force or self.git.is_ancestor_of('FETCH_HEAD', branch):
                        pass
                    else:
                        self._UI.error("Not setting the branch field for ticket #{0} to `{1}` because `{1}` and the current value of the branch field `{2}` have diverged.".format(ticket, remote_branch, current_remote_branch))
                        self._UI.info("If you really want to overwrite the branch field use `{0}`. Otherwise, you need to merge in the changes introduced by `{0}` by using `{1}`.".format(self._format_command("upload",ticket=ticket,remote_branch=remote_branch,force=True), self._format_command("download", ticket=ticket)))
                        raise OperationCancelledError("not a fast-forward")

                if current_remote_branch is not None and not force and not user_confirmation:
                    if not self._UI.confirm("I will now change the branch field of ticket #{0} from its current value `{1}` to `{2}`. Is this what you want?".format(ticket, current_remote_branch, remote_branch), default=True):
                        raise OperationCancelledError("user requested")

                attributes = self.trac._get_attributes(ticket)
                attributes['branch'] = remote_branch
                self.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)

        if ticket:
            old_dependencies = self.trac.dependencies(ticket)
            old_dependencies = ", ".join(["#"+dep for dep in old_dependencies])
            new_dependencies = self._dependencies_for_ticket(ticket)
            new_dependencies = ", ".join(["#"+dep for dep in new_dependencies])
            if old_dependencies != new_dependencies:
                self._UI.info("Uploading your dependencies for ticket #{0}: `{1}` => `{2}`".format(ticket, old_dependencies, new_dependencies))

                attributes = self.trac._get_attributes(ticket)
                attributes['dependencies'] = new_dependencies
                self.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)
            elif new_dependencies:
                self._UI.info("Not uploading your dependencies for ticket #{0} because the dependencies on trac are already up-to-date.".format(ticket))

    def reset_to_clean_state(self):
        r"""
        Reset the current working directory to a clean state.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Nothing happens if the directory is already clean::

            sage: dev.reset_to_clean_state()

        Bring the directory into a non-clean state::

            sage: dev.git.checkout(SUPER_SILENT, b="branch1")
            sage: with open("tracked", "w") as f: f.write("boo")
            sage: dev.git.add(SUPER_SILENT, "tracked")
            sage: dev.git.commit(SUPER_SILENT, message="added tracked")

            sage: dev.git.checkout(SUPER_SILENT, 'HEAD~')
            sage: dev.git.checkout(SUPER_SILENT, b="branch2")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: dev.git.add(SUPER_SILENT, "tracked")
            sage: dev.git.commit(SUPER_SILENT, message="added tracked")
            sage: from sage.dev.git_error import GitError
            sage: try:
            ....:     dev.git.merge(SUPER_SILENT, "branch1")
            ....: except GitError: pass
            sage: UI.append("n")
            sage: dev.reset_to_clean_state()
            Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard many changes which are not commited.) [yes/No] n
            sage: UI.append("y")
            sage: dev.reset_to_clean_state()
            Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard many changes which are not commited.) [yes/No] y
            sage: dev.reset_to_clean_state()

        A detached HEAD does not count as a non-clean state::

            sage: dev.git.checkout(SUPER_SILENT, 'HEAD', detach=True)
            sage: dev.reset_to_clean_state()

        """
        states = self.git.get_state()
        if not states:
            return
        if not self._UI.confirm("Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard many changes which are not commited.)", default=False):
            raise OperationCancelledError("User requested not to clean the current state.")

        self.git.reset_to_clean_state()

    def reset_to_clean_working_directory(self):
        r"""
        Drop any uncommitted changes in the working directory.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Check that nothing happens if there no changes::

            sage: dev.reset_to_clean_working_directory()

        Check that nothing happens if there are only untracked files::

            sage: open("untracked","w").close()
            sage: dev.reset_to_clean_working_directory()

        Uncommitted changes can simply be dropped::

            sage: open("tracked","w").close()
            sage: dev.git.add(SUPER_SILENT, "tracked")
            sage: dev.git.commit(SUPER_SILENT, message="added tracked")
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("discard")
            sage: dev.reset_to_clean_working_directory()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] discard
            sage: dev.reset_to_clean_working_directory()

        Uncommitted changes can be kept::

            sage: with open("tracked", "w") as f: f.write("foo")
            sage: UI.append("keep")
            sage: dev.reset_to_clean_working_directory()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] keep

        Or stashed::

            sage: UI.append("stash")
            sage: dev.reset_to_clean_working_directory()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] stash
            Your changes have been recorded on a new branch `stash/1`.
            sage: dev.reset_to_clean_working_directory()

        """
        try:
            self.reset_to_clean_state()
        except OperationCancelledError:
            self._UI.error("Can not clean the working directory unless in a clean state.")
            raise

        if not self.git.has_uncommitted_changes():
            return

        from git_interface import READ_OUTPUT
        files = "\n".join([line[2:] for line in self.git.status(READ_OUTPUT, porcelain=True).splitlines() if not line.startswith('?')])
        sel = self._UI.select("The following files in your working directory contain uncommitted changes:\n{0}\nDo you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later?".format(files), options=('discard','keep','stash'), default=1)
        if sel == 'discard':
            self.git.reset_to_clean_working_directory()
        elif sel == 'keep':
            raise OperationCancelledError("User requested not to clean the working directory.")
        elif sel == 'stash':
            from git_error import DetachedHeadError
            try:
                current_branch = self.git.current_branch()
            except DetachedHeadError:
                current_branch = None
                current_commit = self.git.current_commit()

            branch = self._new_local_branch_for_stash()
            try:
                try:
                    from git_interface import SUPER_SILENT
                    self.git.stash(SUPER_SILENT)
                    try:
                        self._UI.info("Creating a new branch `{0}` which contains your stashed changes.".format(branch))
                        self.git.stash(SUPER_SILENT,'branch',branch,'stash@{0}')
                        self._UI.info("Committing your changes to `{0}`.".format(branch))
                        self.git.commit(SUPER_SILENT,'-a',message="Changes stashed by reset_to_clean_working_directory()")
                    except:
                        self.git.stash(SUPER_SILENT, 'drop')
                        raise
                except:
                    if self._is_local_branch_name(branch, exists=True):
                        self.git.branch(SUPER_SILENT,"-D",branch)
                    raise
            finally:
                self.git.checkout(SUPER_SILENT, current_branch or current_commit)

            self._UI.show("Your changes have been recorded on a new branch `{0}`.".format(branch))
            self._UI.info("To recover your changes later use {1}.".format(branch, self._format_command("unstash",branch)))
        else:
            assert False

    def unstash(self, branch=None):
        r"""
        Unstash the changes recorded in ``branch``.

        INPUT:

        - ``branch`` -- the name of a local branch or ``None`` (default:
          ``None``), if ``None`` list all stashes.

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Create some stashes::

            sage: dev.unstash()
            (no stashes)
            sage: with open("tracked", "w") as f: f.write("foo")
            sage: dev.git.add("tracked")
            sage: UI.append("s")
            sage: dev.reset_to_clean_working_directory()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] s
            Your changes have been recorded on a new branch `stash/1`.
            sage: with open("tracked", "w") as f: f.write("boo")
            sage: dev.git.add("tracked")
            sage: UI.append("s")
            sage: dev.reset_to_clean_working_directory()
            The following files in your working directory contain uncommitted changes:
             tracked
            Do you want me to discard any changes which are not committed? Should the changes be kept? Or do you want to stash them for later? [discard/Keep/stash] s
            Your changes have been recorded on a new branch `stash/2`.
            sage: dev.unstash()
            stash/1
            stash/2

        Unstash a change::

            sage: dev.unstash("stash/1")

        Unstash something that is not a stash::

            sage: dev.unstash("HEAD")
            ValueError: `HEAD` is not a valid name for a stash.

        Unstash a conflicting change::

            sage: dev.unstash("stash/2")
            The changes recorded in `stash/2` do not apply cleanly to your working directory.

        """
        if branch is None:
            stashes = [stash for stash in self.git.local_branches() if self._is_stash_name(stash)]
            stashes.sort()
            stashes = "\n".join(stashes)
            stashes = stashes or "(no stashes)"
            self._UI.info("Use `{0}` to apply the changes recorded in the stash to your working directory where `name` is one of the following:\n{1}".format(self._format_command("unstash",branch="name"), stashes))
            self._UI.show(stashes)
            return

        self._check_stash_name(branch, exists=True)

        self.reset_to_clean_state()

        from git_interface import SUPER_SILENT
        try:
            self.git.cherry_pick(SUPER_SILENT, branch, no_commit=True)
        except GitError as e:
            self._UI.error("The changes recorded in `{0}` do not apply cleanly to your working directory.".format(branch))
            self._UI.info("You can try to resolve the conflicts manually with `{0}`.".format(self._format_command("merge", branch_or_ticket=branch)))
            raise OperationCancelledError("unstash failed")

        self.git.reset(SUPER_SILENT)

        self._UI.info("The changes recorded in `{0}` have been restored in your working directory. If you do not need the stash anymore, you can drop it with `{1}`.".format(branch, self._format_command("abandon",branch=branch)))

    def edit_ticket(self, ticket=None):
        r"""
        Edit the description of ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          edit the :meth:`_current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`add_comment`

        TESTS:

        Set up a single user for doctesting::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDevWrapper(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

        Create a ticket and edit it::

            sage: UI.append("Summary: summary1\ndescription")
            sage: dev.create_ticket()
            1
            sage: UI.append("Summary: summary1\ndescription...")
            sage: dev.edit_ticket()
            sage: dev.trac._get_attributes(1)
            {'description': 'description...', 'summary': 'summary1'}

        """
        if ticket is None:
            ticket = self._current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        self._check_ticket_name(ticket, exists=True)
        ticket = self._ticket_from_ticket_name(ticket)
        self.trac.edit_ticket_interactive(ticket)

    def add_comment(self, ticket=None):
        r"""
        Add a comment to ``ticket`` on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          edit the :meth:`_current_ticket`.

        .. SEEALSO::

            :meth:`create_ticket`, :meth:`edit_ticket`

        TESTS::

            TODO

        """
        if ticket is None:
            ticket = self.current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        ticket = self._ticket_from_ticket_name(ticket)
        self.trac.add_comment(ticket)

    def browse_ticket(self, ticket=None):
        r"""
        Start a webbrowser at the ticket page on trac.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``), the number of the ticket to edit. If ``None``,
          edit the :meth:`current_ticket`.

        .. SEEALSO::

            :meth:`edit_ticket`, :meth:`add_comment`

        TESTS::

            TODO

        """
        if ticket is None:
            ticket = self.current_ticket()

        if ticket is None:
            raise SageDevValueError("ticket must be specified if not currently on a ticket.")

        ticket = self._ticket_from_ticket_name(ticket)

        from sage.misc.viewer import browser
        from sage.env import TRAC_SERVER_URI
        browser_cmdline = browser() + ' ' + TRAC_SERVER_URI + '/ticket/' + str(ticket)
        import os
        os.system(browser_cmdline)

    def remote_status(self, ticket=None, quiet=False):
        r"""
        Show the remote status of ``ticket``.

        For tickets and remote branches, this shows the commit log of the branch on
        the trac ticket a summary of their difference to your related branches, and
        an overview of patchbot results (where applicable).

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          or ``'all'`` (default: ``None``), the number of the ticket to edit.
          If ``None``, edit the :meth:`current_ticket`.

        .. SEEALSO::

        - :meth:`local_tickets` -- Just shows local tickets without
          comparing them to the remote server.

        - :meth:`diff` -- Shows the actual differences on a given
          ticket.

        - :meth:`download` -- Merges in the changes on a given ticket
          from the remote server.

        - :meth:`upload` -- Pushes the changes on a given ticket to
          the remote server.

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
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

    def import_patch(self, patchname=None, url=None, local_file=None, diff_format=None, header_format=None, path_format=None):
        r"""
        Import a patch into the branch for the current ticket.

        If ``local_file`` is specified, apply the file it points to.

        Otherwise, download the patch using :meth:`download_patch` and apply
        it.

        INPUT:

        - ``patchname`` -- a string or ``None`` (default: ``None``)

        - ``url`` -- a string or ``None`` (default: ``None``)

        - ``local_file`` -- a string or ``None`` (default: ``None``)

        .. SEEALSO::

        - :meth:`download_patch` -- This function downloads a patch to
          a local file.

        - :meth:`download` -- This function is used to merge in
          changes from a git branch rather than a patch.

        TESTS::

            TODO

        """
        try:
            self.reset_to_clean_state()
            # TODO: clean untracked, otherwise the add below is not correct
            self.reset_to_clean_working_directory()
        except OperationCancelledError:
            self._UI.error("Cannot import patch. Your working directory is not in a clean state.")
            raise

        if not local_file:
            return self.import_patch(
                    local_file=self.download_patch(patchname=patchname, url=url),
                    diff_format=diff_format, header_format=header_format, path_format=path_format)
        elif patchname or url:
            raise SageDevValueError("if local_file is specified, patchname and url must not be specified")
        else:
            lines = open(local_file).read().splitlines()
            lines = self._rewrite_patch(lines, to_header_format="git",
                    to_path_format="new", from_diff_format=diff_format,
                    from_header_format=header_format,
                    from_path_format=path_format)

            import tempfile, os
            fd, outfile = tempfile.mkstemp(dir=self.tmp_dir)
            os.fdopen(fd, 'w').writelines("\n".join(lines)+"\n")

            self._UI.info("Trying to apply reformatted patch `%s`"%outfile)
            try:
                self.git.am(outfile, ignore_whitespace=True, resolvemsg='')
            except GitError:
                if not self._UI.confirm("The patch does not apply cleanly. Would you like to apply it anyway and create reject files for the parts that do not apply?", default=False):
                    self._UI.info("Not applying patch.")
                    self.git.reset_to_clean_state()
                    raise OperationCancelledError("User requested to cancel the apply.")

                try:
                    self.git.apply(output, ignore_whitespace=True, reject=True)
                except GitError:
                    if self._UI.select("The patch did not apply cleanly. Please integrate the `.rej` files that were created and resolve conflicts. After you do, type `resolved`. If you want to abort this process, type `abort`.", ("resolved","abort")) == "abort":
                        self.git.reset_to_clean_state()
                        self.git.reset_to_clean_working_directory()
                        raise OperationCancelledError("User requested to cancel the apply.")

                self._UI.show("It seemed that the patch would not apply, but in fact it did.")

                # TODO: add untracked (except .rej)
                self.git.add(update=True)
                self.git.am(resolved=True)

                # TODO: reset to clean to remove the .rej files

    def download_patch(self, ticket=None, patchname=None, url=None):
        r"""
        Download a patch to a temporary directory.

        If only ``ticket`` is specified and the ticket has only one
        attachment, download the patch attached to ``ticket``.

        If ``ticket`` and ``patchname`` are specified, download the
        patch ``patchname`` attached to ``ticket``.

        If ``url`` is specified, download ``url``.

        If nothing is specified, and if the ''current'' ticket has only
        one attachment, download it.

        Raise an error on any other combination of parameters.

        INPUT:

        - ``ticket`` -- an integer or string identifying a ticket or ``None``
          (default: ``None``)

        - ``patchname`` -- a string or ``None`` (default: ``None``)

        - ``url`` -- a string or ``None`` (default: ``None``)

        OUTPUT:

        Returns the absolute file name of the returned file.

        .. SEEALSO::

        - :meth:`import_patch` -- also creates a commit on the current branch
          from the patch.

        TESTS::

            TODO

        """
        if url:
            if ticket or patchname:
                raise ValueError("If `url` is specifed, `ticket` and `patchname` must not be specified.")
            import tempfile
            fd, ret = tempfile.mkstemp(dir=self.tmp_dir)
            import os
            os.close(fd)
            from subprocess import check_call
            check_call(["wget","-r","--no-check-certificate", "-O",ret,url])
            return ret
        elif ticket:
            ticket = self._ticket_from_ticket_name(ticket)

            if patchname:
                from sage.env import TRAC_SERVER_URI
                return self.download_patch(url = TRAC_SERVER_URI+"/raw-attachment/ticket/%s/%s"%(ticket,patchname))
            else:
                attachments = self.trac.attachment_names(ticket)
                if len(attachments) == 0:
                    raise SageDevValueError("Ticket #%s has no attachments."%ticket)
                if len(attachments) == 1:
                    return self.download_patch(ticket = ticket, patchname = attachments[0])
                else:
                    raise SageDevValueError("Ticket #%s has more than one attachment but parameter `patchname` is not present."%ticket)
        elif not patchname:
            return self.download_patch(ticket=self.current_ticket())
        else:
            raise SageDevValueError("If `url` is not specified, `ticket` must be specified")

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

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
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

            :meth:`abandon_ticket` -- Abandon a single ticket or branch.

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
        for branch in self.git.local_branches():
            if self.git.is_ancestor_of(branch, MASTER_BRANCH):
                self._UI.show("Abandoning %s"%branch)
                self.git.abandon(branch)

    def abandon(self, ticket_or_branch=None):
        r"""
        Abandon a ticket branch.

        INPUT:

        - ``ticket_or_branch`` -- an integer or string identifying a ticket or
          the name of a local branch or ``None`` (default: ``None``), remove
          the branch ``ticket_or_branch`` or the branch for the ticket
          ``ticket_or_branch`` (or the current branch if ``None``). Also
          removes the users remote tracking branch.

        .. SEEALSO::

        - :meth:`prune_closed_tickets` -- abandon tickets that have
          been closed.

        - :meth:`local_tickets` -- list local tickets (by default only
          showing the non-abandoned ones).

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
        if self._UI.confirm("Are you sure you want to delete your work on %s?"%self._ticket_repr(ticketnum), default=False):
            self.git.abandon(ticketnum)

    def gather(self, branchname, *tickets, **kwds):
        r"""
        Create a new branch ``branchname`` with ``tickets`` applied.

        INPUT:

        - ``branchname`` -- a string, the name of the new branch

        - ``tickets`` -- a list of integers and strings; for an integer or
          string identifying a ticket, the branch on the trac ticket gets
          merged, for the name of a local or remote branch, that branch gets
          merged.

        - ``create_dependency`` -- boolean (default: ``True``, keyword only),
          whether to append the other ticket to the list of dependencies.  See
          :meth:`merge` for the consequences of having another branch as a
          dependency.

        - ``download`` -- boolean (default: ``False``, keyword only), whether
          to download the most recent version of the other tickets before
          merging.

        .. SEEALSO::

        - :meth:`merge` -- merge into the current branch rather than
          creating a new one.

        - :meth:`show_dependencies` -- show the dependencies of a
          given branch.

        TEST::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
        create_dependencies = kwds.pop('create_dependencies', True)
        download = kwds.pop('download', False)
        if len(tickets) == 0:
            raise ValueError("must include at least one input branch")
        if self.git.commit_for_branch(branchname):
            if not self._UI.confirm("The branch %s already "%branchname+
                                    "exists; do you want to merge into it?",
                                    default=False):
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
        Show the dependencies of ``ticket``.

        INPUT:

        - ``ticket`` -- a string or integer identifying a ticket or ``None``
          (default: ``None``), the ticket for which dependencies are displayed.
          If ``None``, then the dependencies for the current ticket are
          displayed.

        - ``all`` -- boolean (default: ``True``), whether to recursively list
          all tickets on which this ticket depends (in depth-first order), only
          including tickets that have a local branch.

        .. NOTE::

            Ticket dependencies are stored locally and only updated with
            respect to the remote server during :meth:`upload` and
            :meth:`download`.

        .. SEEALSO::

        - :meth:`TracInterface.dependencies` -- Query Trac to find
          dependencies.

        - :meth:`remote_status` -- will show the status of tickets
          with respect to the remote server.

        - :meth:`merge` -- Merge in changes from a dependency.

        - :meth:`diff` -- Show the changes in this branch over the
          dependencies.

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
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
        TODO: fix this docstring

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
        raise NotImplementedError # the below does most probably not work anymore
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

        - ``abandoned`` -- boolean (default: ``False``), whether to show
          abandoned branches.

        - ``quite`` -- boolean (default: ``False``), whether to show return the
          list of branches rather than printing them.

        .. SEEALSO::

        - :meth:`abandon_ticket` -- hide tickets from this method.

        - :meth:`remote_status` -- also show status compared to the
          trac server.

        - :meth:`current_ticket` -- get the current ticket.

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
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

    def vanilla(self, release="release"):
        r"""
        Returns to an official release of Sage.

        INPUT:

        - ``release`` -- a string or decimal giving the release name.
          In fact, any tag, commit or branch will work.  If the tag
          does not exist locally an attempt to fetch it from the
          server will be made.

        Git equivalent::

            Checks out a given tag, commit or branch in detached head mode.

        .. SEEALSO::

        - :meth:`switch_ticket` -- switch to another branch, ready to
          develop on it.

        - :meth:`download` -- download a branch from the server and
          merge it.

        TESTS::

            TODO

        """
        raise NotImplementedError # the below most probably does not work anymore
        if hasattr(release, 'literal'):
            release = release.literal
        release = str(release)
        if self._UI.confirm("Are you sure you want to revert to %s?"%(release)):
            self.git.switch_branch(release, detached = True)

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

            sage: dev._wrap("_detect_patch_diff_format")
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
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'git'
            sage: dev._detect_patch_diff_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","diff.patch"
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
                        raise SageDevValueError("File appears to have mixed diff formats.")

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

            sage: dev._wrap("_detect_patch_path_format")
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
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
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
                                    raise SageDevValueError("File appears to have mixed path formats.")

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

            sage: dev._wrap("_rewrite_patch_diff_paths")
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
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
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

            sage: dev._wrap("_detect_patch_header_format")
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
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines())
            'diff'
            sage: dev._detect_patch_header_format(
            ....:     open(os.path.join(
            ....:             SAGE_SRC,"sage","dev","test","data","diff.patch"
            ....:         )).read().splitlines())
            'diff'
        """
        lines = list(lines)
        if not lines:
            raise SageDevValueError("patch is empty")

        if HG_HEADER_REGEX.match(lines[0]):
            if HG_USER_REGEX.match(lines[1]):
                return "hg-export"
            elif HG_PARENT_REGEX.match(lines[1]):
                return "hg"
        elif GIT_FROM_REGEX.match(lines[0]):
            return "git"

        return "diff"

    def _detect_patch_modified_files(self, lines, diff_format = None):
        #TODO: docstring
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
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "hg.patch")
            ....:     ).read().splitlines()
            sage: hg_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "hg-output.patch")
            ....:     ).read().splitlines()
            sage: git_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "git.patch")
            ....:     ).read().splitlines()
            sage: git_output_lines = open(
            ....:     os.path.join(SAGE_SRC, "sage", "dev", "test", "data", "git-output.patch")
            ....:     ).read().splitlines()

            sage: dev._wrap("_rewrite_patch_header")
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
            ....:             SAGE_SRC,"sage","dev","test","data","trac_8703-trees-fh.patch"
            ....:         )).read().splitlines(), 'git')[:5]
            ['From: "Unknown User" <unknown@sagemath.org>',
            'Subject: #8703: Enumerated sets and data structure for ordered and binary trees',
            'Date: ...',
            '',
            '- The Class Abstract[Labelled]Tree allows for inheritance from different']
        """
        import email.utils, time

        lines = list(lines)
        if not lines:
            raise SageDevValueError("empty patch file")

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
                        raise SageDevValueError("Malformed patch. Missing line for regular expression `%s`."%(regex.pattern))
                    else:
                        return
                match = regex.match(lines[i])
                if match is not None:
                    if len(match.groups()) > 0:
                        header[key] = match.groups()[0]
                    i += 1
                elif mandatory:
                    raise SageDevValueError("Malformed patch. Line `%s` does not match regular expression `%s`."%(lines[i],regex.pattern))

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
        #TODO: docstring
        return self._rewrite_patch_diff_paths(self._rewrite_patch_header(lines, to_format=to_header_format, from_format=from_header_format, diff_format=from_diff_format), to_format=to_path_format, diff_format=from_diff_format, from_format=from_path_format)

    # TODO: make sure this gets called before accessing the remote repository
    def _upload_ssh_key(self, keyfile, create_key_if_not_exists=True):
        r"""
        Upload ``keyfile.pub`` to gitolite through the trac interface.

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

        TESTS::

            TODO

        """
        raise NotImplementedError # the below does most probably not work anymore
        cfg = self._config

        try:
            with open(keyfile, 'r') as F:
                pass
            with open(keyfile + '.pub', 'r') as F:
                pass
        except IOError:
            if create_key_if_not_exists:
                self._UI.show("Generating ssh key....")
                from subprocess import call
                success = call(["ssh-keygen", "-q", "-f", keyfile, "-P", ""])
                if success == 0:
                    self._UI.show("Ssh key successfully generated")
                else:
                    raise RuntimeError("Ssh key generation failed.  Please create a key in `%s` and retry"%(keyfile))
            else:
                raise

        with open(keyfile + '.pub', 'r') as F:
            pubkey = F.read().strip()

        self.trac._authenticated_server_proxy.sshkeys.setkeys(pubkey)

    def _is_ticket_name(self, name, exists=False):
        r"""
        Return whether ``name`` is a valid ticket name, i.e., an integer.

        INPUT:

        - ``name`` -- a string or an int

        - ``exists`` -- a boolean (default: ``False``), if ``True``, return
          whether ``name`` is the name of an existing ticket

        EXAMPLES::

            sage: dev = dev._sagedev
            sage: dev._is_ticket_name(1000)
            True
            sage: dev._is_ticket_name("1000")
            True
            sage: dev._is_ticket_name("1 000")
            False
            sage: dev._is_ticket_name("#1000")
            True
            sage: dev._is_ticket_name("master")
            False
            sage: dev._is_ticket_name(1000, exists=True) # optional: internet
            True
            sage: dev._is_ticket_name(2^30, exists=True) # optional: internet
            False

        """
        if not isinstance(name, int):
            try:
                name = self._ticket_from_ticket_name(name)
            except SageDevValueError:
                return False

        if exists:
            try:
                self.trac._anonymous_server_proxy.ticket.get(name)
            except TracInternalError as e:
                if e.faultCode == 404: # ticket does not exist
                    return False
                raise
            except TracConnectionError as e:
                # if we cannot connect to trac, we assume that the ticket
                # exists; this makes more of the dev scripts usable in offline
                # scenarios
                pass

        return True

    def _check_ticket_name(self, name, exists=False):
        r"""
        Check that ``name`` is a valid ticket name.

        INPUT:

        - ``name`` -- a string or int

        - ``exists`` -- a boolean (default: ``False``), whether to check that
          the ticket exists on trac

        TESTS::

            sage: dev = dev._sagedev
            sage: dev._check_ticket_name(1000)
            sage: dev._check_ticket_name("1000")
            sage: dev._check_ticket_name("1 000")
            Traceback (most recent call last):
            ...
            SageDevValueError: `1 000` is not a valid ticket name.
            sage: dev._check_ticket_name("#1000")
            sage: dev._check_ticket_name("master")
            Traceback (most recent call last):
            ...
            SageDevValueError: `master` is not a valid ticket name.
            sage: dev._check_ticket_name(1000, exists=True) # optional: internet
            sage: dev._check_ticket_name(2^30, exists=True) # optional: internet
            Traceback (most recent call last):
            ...
            SageDevValueError: `1073741824` is not a valid ticket name or ticket does not exist on trac.

        """
        if not self._is_ticket_name(name, exists=exists):
            if exists:
                raise SageDevValueError("`{0}` is not a valid ticket name or ticket does not exist on trac.".format(name))
            else:
                raise SageDevValueError("`{0}` is not a valid ticket name.".format(name))

    def _ticket_from_ticket_name(self, name):
        r"""
        Return the ticket number for the ticket ``name``.

        EXAMPLES::

            sage: dev = dev._sagedev
            sage: dev._ticket_from_ticket_name("1000")
            1000
            sage: dev._ticket_from_ticket_name("#1000")
            1000
            sage: dev._ticket_from_ticket_name(1000)
            1000
            sage: dev._ticket_from_ticket_name(int(1000))
            1000
            sage: dev._ticket_from_ticket_name("1 000")
            Traceback (most recent call last):
            ...
            SageDevValueError: `1 000` is not a valid ticket name.

        """
        ticket = name
        if not isinstance(ticket, int):
            if isinstance(ticket, str) and ticket[0] == "#":
                ticket = ticket[1:]
            try:
                ticket = int(ticket)
            except ValueError:
                raise SageDevValueError("`{0}` is not a valid ticket name.".format(name))

        if ticket < 0:
            raise SageDevValueError("`{0}` is not a valid ticket name.".format(name))

        return ticket

    def _is_local_branch_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for a local branch.

        INPUT:

        - ``name`` -- a string

        - ``exists`` -- a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing local branch; if
          ``False``, check whether ``name`` is the name of a branch that does
          not exist yet.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._chdir()
            sage: dev._pull_master_branch()

            sage: dev._is_local_branch_name('')
            False
            sage: dev._is_local_branch_name('ticket/1')
            True
            sage: dev._is_local_branch_name('ticket/1', exists=True)
            False
            sage: dev._is_local_branch_name('ticket/1', exists=False)
            True
            sage: dev.git.branch('ticket/1')
            sage: dev._is_local_branch_name('ticket/1', exists=True)
            True
            sage: dev._is_local_branch_name('ticket/1', exists=False)
            False

        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")

        if not GIT_BRANCH_REGEX.match(name):
            return False

        if exists == True:
            return self.git.commit_for_branch(name) is not None
        elif exists == False:
            return self.git.commit_for_branch(name) is None
        elif exists is any:
            return True
        else:
            raise ValueError("exists")

    def _is_stash_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for a stash.

        INPUT:

        - ``name`` -- a string

        - ``exists`` - a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing stash; if
          ``False``, check whether ``name`` is the name of a stash that does
          not exist yet.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._chdir()
            sage: dev._pull_master_branch()

            sage: dev._is_stash_name("branch1")
            False
            sage: dev._is_stash_name("stash")
            False
            sage: dev._is_stash_name("stash/")
            False
            sage: dev._is_stash_name("stash/1")
            True
            sage: dev._is_stash_name("stash/1", exists=True)
            False

        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")

        if not name.startswith("stash/"):
            return False

        return self._is_local_branch_name(name, exists)

    def _check_stash_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a stash.

        INPUT:

        - ``name`` -- a string

        - ``exists`` - a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing stash; if
          ``False``, check whether ``name`` is the name of a stash that does
          not exist yet.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._chdir()
            sage: dev._pull_master_branch()

            sage: dev._check_stash_name("stash/1")
            sage: dev._check_stash_name("stash/1", exists=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: `stash/1` does not exist.
            sage: dev._check_stash_name("stash/1", exists=False)

        """
        if not self._is_stash_name(name):
            raise SageDevValueError("`{0}` is not a valid name for a stash.".format(name))
        if exists == True and not self._is_stash_name(name, exists):
            raise SageDevValueError("`{0}` does not exist.".format(name))
        elif exists == False and not self._is_stash_name(name, exists):
            raise SageDevValueError("`{0}` already exists, please choose a different name for the stash.")

    def _is_remote_branch_name(self, name, exists=any):
        r"""
        Return whether ``name`` is a valid name for a remote branch.

        INPUT:

        - ``name`` -- a string

        - ``exists`` -- a boolean or ``any`` (default: ``any``), if ``True``,
          check whether ``name`` is the name of an existing remote branch; if
          ``False``, check whether ``name`` is the name of a branch that does
          not exist yet.

        .. NOTE::

            Currently, this does not check whether name is in accordance with
            naming scheme configured on gitolite.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._chdir()
            sage: dev._pull_master_branch()

            sage: dev._is_remote_branch_name('')
            False
            sage: dev._is_remote_branch_name('ticket/1')
            True

            sage: dev._is_remote_branch_name('ticket/1', exists=True)
            False
            sage: dev._is_remote_branch_name('ticket/1', exists=False)
            True

        """
        if not isinstance(name, str):
            raise ValueError("name must be a string")

        if not GIT_BRANCH_REGEX.match(name):
            return False

        if exists is any:
            return True

        from git_error import GitError
        from git_interface import SUPER_SILENT
        try:
            self.git.ls_remote(SUPER_SILENT, self.git._repository, name, exit_code=True)
            remote_exists = True
        except GitError as e:
            if e.exit_code == 2:
                remote_exists = False
            else:
                raise

        if exists == True or exists == False:
            return remote_exists == exists
        else:
            raise ValueError("exists")

    def _check_local_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a local branch, raise a
        ``SageDevValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_local_branch_name`

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._check_local_branch_name('')
            Traceback (most recent call last):
            ...
            SageDevValueError: `` is not a valid name for a local branch.
            sage: dev._check_local_branch_name('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch `ticket/1` does not exist locally.
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            sage: dev.git.branch('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch `ticket/1` already exists, please choose a different name.

        """
        try:
            if not self._is_local_branch_name(name, exists=any):
                raise SageDevValueError("caught below")
        except SageDevValueError:
            raise SageDevValueError("`{0}` is not a valid name for a local branch.".format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_local_branch_name(name, exists=exists):
                raise SageDevValueError("Branch `{0}` does not exist locally.".format(name))
        elif exists == False:
            if not self._is_local_branch_name(name, exists=exists):
                raise SageDevValueError("Branch `{0}` already exists, please choose a different name.".format(name))
        else:
            assert False

    def _check_remote_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a remote branch, raise a
        ``SageDevValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_remote_branch_name`

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._check_remote_branch_name('')
            Traceback (most recent call last):
            ...
            SageDevValueError: `` is not a valid name for a remote branch.
            sage: dev._check_remote_branch_name('ticket/1')

            sage: dev._check_remote_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch `ticket/1` does not exist on the remote system.
            sage: dev._check_remote_branch_name('ticket/1', exists=False)

        """
        try:
            if not self._is_remote_branch_name(name, exists=any):
                raise SageDevValueError("caught below")
        except SageDevValueError:
            raise SageDevValueError("`{0}` is not a valid name for a remote branch.".format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_remote_branch_name(name, exists=exists):
                raise SageDevValueError("Branch `{0}` does not exist on the remote system.".format(name))
        elif exists == False:
            if not self._is_remote_branch_name(name, exists=exists):
                raise SageDevValueError("Branch `{0}` already exists, please choose a different name.".format(name))
        else:
            assert False

    def _remote_branch_for_ticket(self, ticket):
        r"""
        Return the name of the remote branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        .. NOTE:

            This does not take into account the ``branch`` field of the ticket
            on trac.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("#1")
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("1")
            'u/doctest/ticket/1'
            sage: dev._remote_branch_for_ticket("master")
            Traceback (most recent call last):
            ...
            SageDevValueError: `master` is not a valid ticket name.

            sage: UI.append("Summary: summary1\ndescription")
            sage: ticket = dev.create_ticket()

            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_ticket(1)
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'

        """
        ticket = self._ticket_from_ticket_name(ticket)

        default = "u/{0}/ticket/{1}".format(self.trac._username, ticket)

        try:
            branch = self._local_branch_for_ticket(ticket)
        except KeyError: # ticket has no branch yet
            return default

        ret = self._remote_branch_for_branch(branch)
        if ret is None:
            return default
        return ret

    def _has_local_branch_for_ticket(self, ticket):
        r"""
        Return whether there is a local branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._has_local_branch_for_ticket(1)
            False

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if ticket not in self.__ticket_to_branch:
            return False

        branch = self.__ticket_to_branch[ticket]
        if not self._is_local_branch_name(branch, exists=True):
            self._UI.warning("Ticket #{0} refers to the non-existant local branch {1}. If you have not manually interacted with git, then this is a bug in sagedev. Removing the association from ticket #{0} to branch {1}.".format(ticket, branch))
            del self.__ticket_to_branch[ticket]
            return False

        return True

    def _local_branch_for_ticket(self, ticket, download_if_not_found=False):
        r"""
        Return the name of the local branch for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or a string identifying a ticket

        - ``download_if_not_found`` -- a boolean (default: ``False``), whether
          to attempt to download a branch for ``ticket`` from trac if it does
          not exist locally

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config1 = DoctestConfig()
            sage: config1['trac']['password'] = 'secret'
            sage: dev1 = DoctestSageDev(config1, server)
            sage: dev1._pull_master_branch()

            sage: config2 = DoctestConfig('doctest2')
            sage: config2['trac']['password'] = 'secret'
            sage: dev2 = DoctestSageDev(config2, server)
            sage: dev2._pull_master_branch()

        If a local branch for the ticket exists, its name is returned::

            sage: dev1._chdir()
            sage: dev1._UI.append("Summary: ticket1\ndescription")
            sage: ticket = dev1.create_ticket()
            sage: dev1._local_branch_for_ticket(ticket)
            'ticket/1'

        If no local branch exists, the behaviour depends on ``download_if_not_found``::

            sage: dev2._chdir()
            sage: dev2._local_branch_for_ticket(ticket)
            Traceback (most recent call last):
            ...
            KeyError: 'No branch for ticket #1 in your repository.'
            sage: dev2._local_branch_for_ticket(ticket, download_if_not_found=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch field is not set for ticket #1 on trac.
            sage: attributes = dev1.trac._get_attributes(ticket)
            sage: attributes['branch'] = 'public/ticket/1'
            sage: dev1.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)
            'https://trac.sagemath.org/ticket/1#comment:1'
            sage: dev2._local_branch_for_ticket(ticket, download_if_not_found=True)
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch `public/ticket/1` does not exist on the remote system.

            sage: import os
            sage: os.chdir(server.git._config['src'])
            sage: server.git.branch('public/ticket/1')
            sage: dev2._chdir()
            sage: dev2._local_branch_for_ticket(ticket, download_if_not_found=True)
            'ticket/1'
            sage: dev2._local_branch_for_ticket(ticket)
            'ticket/1'

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if self._has_local_branch_for_ticket(ticket):
            return self.__ticket_to_branch[ticket]

        if not download_if_not_found:
            raise KeyError("No branch for ticket #{0} in your repository.".format(ticket))

        branch = self._new_local_branch_for_ticket(ticket)
        self.download(ticket, branch)
        self._set_local_branch_for_ticket(ticket, branch)
        return self._local_branch_for_ticket(ticket, download_if_not_found=False)

    def _new_local_branch_for_stash(self):
        r"""
        Return a new local branch name for a stash.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._new_local_branch_for_stash()
            'stash/1'
            sage: dev.git.branch('stash/1')
            sage: dev._new_local_branch_for_stash()
            'stash/2'

        """
        i = 0
        while True:
            i+=1
            branch = 'stash/{0}'.format(i)
            if self._is_stash_name(branch, exists=False):
                return branch

    def _new_local_branch_for_ticket(self, ticket):
        r"""
        Return a local branch name for ``ticket`` which does not exist yet.

        INPUT:

        - ``ticket`` -- a string or an int identifying a ticket

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._new_local_branch_for_ticket(1)
            'ticket/1'
            sage: dev.git.branch('ticket/1')
            sage: dev._new_local_branch_for_ticket(1)
            'ticket/1_'

        """
        ticket = self._ticket_from_ticket_name(ticket)

        branch = 'ticket/{0}'.format(ticket)

        while self._is_local_branch_name(branch, exists=True):
            branch = branch + "_"

        assert self._is_local_branch_name(branch, exists=False)

        return branch

    def _set_dependencies_for_ticket(self, ticket, dependencies):
        r"""
        Locally record ``dependencies`` for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or string identifying a ticket

        - ``dependencies`` -- an iterable of ticket numbers or ``None`` for no
          dependencies

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: UI.append("Summary: ticket1\ndescription")
            sage: ticket = dev.create_ticket()
            sage: dev._set_dependencies_for_ticket(ticket, [2, 3])
            sage: dev._dependencies_for_ticket(ticket)
            (2, 3)
            sage: dev._set_dependencies_for_ticket(ticket, None)
            sage: dev._dependencies_for_ticket(ticket)
            ()

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if dependencies is None:
            dependencies = []

        dependencies = [self._ticket_from_ticket_name(dep) for dep in dependencies]

        if not(dependencies):
            if ticket in self.__ticket_dependencies:
                del self.__ticket_dependencies[ticket]
            return

        if not self._has_local_branch_for_ticket(ticket):
            raise KeyError("no local branch for ticket #{0} found.".format(ticket))

        self.__ticket_dependencies[ticket] = tuple(sorted(dependencies))

    def _dependencies_for_ticket(self, ticket, download_if_not_found=False):
        r"""
        Return the locally recorded dependencies for ``ticket``.

        INPUT:

        - ``ticket`` -- an int or string identifying a ticket

        - ``download_if_not_found`` -- a boolean (default: ``False``), whether
          to take the information from trac if the ticket does not exist
          locally

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: UI.append("Summary: ticket1\ndescription")
            sage: ticket = dev.create_ticket()

            sage: dev._set_dependencies_for_ticket(ticket, [2, 3])
            sage: dev._dependencies_for_ticket(ticket)
            (2, 3)
            sage: dev._set_dependencies_for_ticket(ticket, None)
            sage: dev._dependencies_for_ticket(ticket)
            ()

            sage: dev._dependencies_for_ticket(2, download_if_not_found=True)
            Traceback (most recent call last):
            ...
            NotImplementedError

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if not self._has_local_branch_for_ticket(ticket):
            if download_if_not_found:
                raise NotImplementedError
            else:
                raise KeyError("no local branch for ticket #{0} found.".format(ticket))
        else:
            ret = self.__ticket_dependencies[ticket]

        return tuple(sorted([self._ticket_from_ticket_name(dep) for dep in ret]))

    def _set_remote_branch_for_branch(self, branch, remote_branch):
        r"""
        Set the remote branch of ``branch`` to ``remote_branch``.

        INPUT:

        - ``branch`` -- a string, a name of a local branch

        - ``remote_branch`` -- a string or ``None``, unset the remote branch if
          ``None``

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev.git.branch('ticket/1')

            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_ticket(1) # ticket/1 has not been set to be the branch for ticket #1
            'u/doctest/ticket/1'
            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            sage: dev._remote_branch_for_ticket(1)
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_ticket(1)
            'u/doctest/ticket/1'

        """
        self._check_local_branch_name(branch, exists=any)

        if remote_branch is None:
            if branch in self.__branch_to_remote_branch:
                del self.__branch_to_remote_branch[branch]
            return

        self._check_local_branch_name(branch, exists=True)
        self._check_remote_branch_name(remote_branch)

        self.__branch_to_remote_branch[branch] = remote_branch

    def _remote_branch_for_branch(self, branch):
        r"""
        Return the remote branch of ``branch`` or ``None`` if no remote branch is set.

        INPUT:

        - ``branch`` -- a string, the name of a local branch

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev.git.branch('ticket/1')

            sage: dev._remote_branch_for_branch('ticket/1') is None
            True
            sage: dev._set_remote_branch_for_branch("ticket/1", "public/1")
            sage: dev._remote_branch_for_branch('ticket/1')
            'public/1'
            sage: dev._set_remote_branch_for_branch("ticket/1", None)
            sage: dev._remote_branch_for_branch('ticket/1') is None
            True

        """
        self._check_local_branch_name(branch, exists=True)

        if branch in self.__branch_to_remote_branch:
            return self.__branch_to_remote_branch[branch]

        return None

    def _set_local_branch_for_ticket(self, ticket, branch):
        r"""
        Record that ``branch`` is the local branch associated to ``ticket``.

        INPUT:

        - ``ticket`` -- a string or int identifying a ticket

        - ``branch`` -- a string, the name of a local branch, or ``None`` to
          delete the association

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._local_branch_for_ticket(1)
            Traceback (most recent call last):
            ...
            KeyError: 'No branch for ticket #1 in your repository.'

            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            Traceback (most recent call last):
            ...
            SageDevValueError: Branch `ticket/1` does not exist locally.
            sage: dev.git.branch('ticket/1')
            sage: dev._set_local_branch_for_ticket(1, 'ticket/1')
            sage: dev._local_branch_for_ticket(1)
            'ticket/1'

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if branch is None:
            if ticket in self.__ticket_to_branch:
                del self.__ticket_to_branch[ticket]
            return

        self._check_local_branch_name(branch, exists=True)

        self.__ticket_to_branch[ticket] = branch

    def _format_command(self, command, *args, **kwargs):
        r"""
        Helper method for informational messages.

        OUTPUT:

        A command which the user can run from the command line/sage interactive
        shell to execute ``command`` with ``args`` and ``kwargs``.

        EXAMPLES::

            sage: dev._format_command('switch-ticket') # not tested (output depends on whether this test is run from within sage or not)
            'dev.switch_ticket()'
            sage: dev._format_command('switch-ticket',int(1)) # not tested
            'dev.switch_ticket(1)'

        """
        try:
            __IPYTHON__
        except NameError:
            args = [str(arg) for arg in args]
            kwargs = [ "--{0}={1}".format(str(key).replace("_","-"),kwargs[key]) for key in kwargs ]
            return "sage --dev {0} {1}".format(command.replace("_","-"), " ".join(args+kwargs))
        else:
            args = [str(arg) for arg in args]
            kwargs = [ "{0}={1}".format(str(key).replace("-","_"),kwargs[key]) for key in kwargs ]
            return "dev.{0}({1})".format(command.replace("-","_"), ", ".join(args+kwargs))

    def _current_ticket(self):
        r"""
        Return the ticket corresponding to the current branch or ``None`` if
        there is no ticket associated to that branch.

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: UI = dev._UI
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._current_ticket() is None
            True

            sage: UI.append("Summary: ticket1\ndescription")
            sage: ticket = dev.create_ticket()
            sage: dev._current_ticket()
            1

        """
        from git_error import DetachedHeadError
        try:
            branch = self.git.current_branch()
        except DetachedHeadError:
            return None

        if branch in self.__branch_to_ticket:
            return self.__branch_to_ticket[branch]

        return None

class SageDevValueError(ValueError):
    r"""
    A ``ValueError`` to indicate that the user supplied an invaid value.

    EXAMPLES::

        sage: dev.switch_ticket(-1)
        ValueError: `-1` is not a valid ticket name or ticket does not exist on trac.

    """
    def __init__(self, message):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.sagedev import SageDevValueError
            sage: type(SageDevValueError("message"))
            <class 'sage.dev.sagedev.SageDevValueError'>

        """
        ValueError.__init__(self, message)
