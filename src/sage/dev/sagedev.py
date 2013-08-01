r"""
SageDev

This module provides :class:`SageDev`, the central object of the developer
scripts for sage.

AUTHORS:

- TODO: add authors from github's history and trac's history

"""
#*****************************************************************************
#       Copyright (C) 2013 TODO
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
GIT_BRANCH_REGEX = re.compile(r'^(?!.*/\.)(?!.*\.\.)(?!/)(?!.*//)(?!.*@\{)(?!.*\\)[^\040\177 ~^:?*[]+/[^\040\177 ~^:?*[]+(?<!\.lock)(?<!/)(?<!\.)$') # http://stackoverflow.com/questions/12093748/how-do-i-check-for-valid-git-branch-names

# the name of the branch which holds the vanilla clone of sage - in the long
# run this should be "master", currently, "public/sage-git/master" contains some changes
# over "master" which have not been reviewed yet but which are needed to work
# using git
MASTER_BRANCH = "public/sage-git/master"

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

        sage: dev
        SageDev()

    """
    def __init__(self, config=None, UI=None, trac=None, git=None):
        r"""
        Initialization.

        TESTS::

            sage: type(dev)
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
            sage: os.path.isdir(dev.tmp_dir)
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
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
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
            Traceback (most recent call last):
            ...
            ValueError: `1000` is not a valid ticket name or ticket does not exist on trac.

        In this example ``base`` does not exist locally::

            sage: dev.trac.create_ticket("summary5","description",{})
            5
            sage: dev.create_ticket(base=5)
            Traceback (most recent call last):
            ...
            ValueError: Branch field is not set for ticket #5 on trac.

        This also fails if the internet connection is broken::

            sage: dev.trac._connected = False
            sage: dev.create_ticket(base=4)
            Traceback (most recent call last):
            ...
            TracConnectionError: Connection to trac server failed.

        """
        dependencies = []

        if branch is not None:
            self._check_local_branch_name(branch, exists=False)

        if base is None:
            base = self.current_ticket
        if base is None:
            raise ValueError("currently on no ticket, base must not be None")
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
            self._UI.info("An error ocurred while switching to your new branch {0}. Use {1} to manually switch to {0}.".format(branch,self._format_command("switch_ticket",str(ticket))))
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

        TESTS::

            TODO

        """
        ticket = self._ticket_from_ticket_name(ticket)

        if branch is None:
            if self._has_local_branch_for_ticket(ticket):
                branch = self._local_branch_for_ticket(ticket)
                self._UI.info("Switching to branch {0}.".format(branch))
                self.switch_branch(branch)
                return
            else:
                branch = self._new_local_branch_for_ticket(branch)
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
            except: #TODO
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
        #TODO
        from sage.dev.git_interface import SUPER_SILENT
        self.git.checkout(SUPER_SILENT, branch)

    def download(self, ticket_or_branch, branch=None):
        r"""
        Download ``ticket_or_branch`` to ``branch``.

        INPUT:

        - ``ticket_or_branch`` -- a string or an integer, a ticket or a remote
          branch name

        - ``branch`` -- a string or ``None`` (default: ``None``), the branch to
          create or merge the changes into. If ``None``, then a new branch will
          be created unless there is already a branch for this ticket.

        TESTS:

        Create a doctest setup with two users::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config1 = DoctestConfig()
            sage: config1['trac']['password'] = 'secret'
            sage: dev1 = DoctestSageDev(config1, server)
            sage: dev1._pull_master_branch()

            sage: config2 = DoctestConfig('doctest2')
            sage: config2['trac']['password'] = 'secret'
            sage: dev2 = DoctestSageDev(config2, server)
            sage: dev2._pull_master_branch()

        """
        if self._is_ticket_name(ticket_or_branch):
            ticket = self._ticket_from_ticket_name(ticket_or_branch)
            self._check_ticket_name(ticket, exists=True)

            remote_branch = self.trac._branch_for_ticket(ticket)
            if remote_branch is None:
                raise ValueError("Branch field is not set for ticket #{0} on trac.".format(ticket))
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
        try:
            from sage.dev.git_interface import SUPER_SILENT
            self.git.fetch(SUPER_SILENT, self.git._repository, "{0}:{1}".format(remote_branch, branch))
        except:
            # TODO
            raise

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
        #TODO
        from sage.dev.git_interface import SUPER_SILENT
        self.git.commit(SUPER_SILENT, all=True, message=message)

    def _is_ticket_name(self, name, exists=False):
        r"""
        Return whether ``name`` is a valid ticket name, i.e., an integer.

        INPUT:

        - ``name`` -- a string or an int

        - ``exists`` -- a boolean (default: ``False``), if ``True``, return
          whether ``name`` is the name of an existing ticket

        EXAMPLES::

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
            except ValueError:
                return False

        if name < 0:
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

            sage: dev._check_ticket_name(1000)
            sage: dev._check_ticket_name("1000")
            sage: dev._check_ticket_name("1 000")
            Traceback (most recent call last):
            ...
            ValueError: `1 000` is not a valid ticket name.
            sage: dev._check_ticket_name("#1000")
            sage: dev._check_ticket_name("master")
            Traceback (most recent call last):
            ...
            ValueError: `master` is not a valid ticket name.
            sage: dev._check_ticket_name(1000, exists=True) # optional: internet
            sage: dev._check_ticket_name(2^30, exists=True) # optional: internet
            Traceback (most recent call last):
            ...
            ValueError: `1073741824` is not a valid ticket name or ticket does not exist on trac.

        """
        if not self._is_ticket_name(name, exists=exists):
            if exists:
                raise ValueError("`{0}` is not a valid ticket name or ticket does not exist on trac.".format(name))
            else:
                raise ValueError("`{0}` is not a valid ticket name.".format(name))

    def _ticket_from_ticket_name(self, name):
        r"""
        Return the ticket number for the ticket ``name``.

        EXAMPLES::

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
            ValueError: `1 000` is not a valid ticket name.

        """
        if not isinstance(name, int):
            if isinstance(name, str) and name[0] == "#":
                name = name[1:]
            try:
                name = int(name)
            except ValueError:
                raise ValueError("`{0}` is not a valid ticket name.".format(name))

        return name

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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            raise ValueError

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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            raise ValueError

    def _check_local_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a local branch, raise a
        ``ValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_local_branch_name`

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._check_local_branch_name('')
            Traceback (most recent call last):
            ...
            ValueError: `` is not a valid branch name.
            sage: dev._check_local_branch_name('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            ValueError: Branch `ticket/1` does not exist locally.
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            sage: dev.git.branch('ticket/1')
            sage: dev._check_local_branch_name('ticket/1', exists=True)
            sage: dev._check_local_branch_name('ticket/1', exists=False)
            Traceback (most recent call last):
            ...
            ValueError: Branch `ticket/1` already exists, please choose a different name.

        """
        try:
            if not self._is_local_branch_name(name, exists=any):
                raise ValueError
        except ValueError:
            raise ValueError("`{0}` is not a valid branch name.".format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_local_branch_name(name, exists=exists):
                raise ValueError("Branch `{0}` does not exist locally.".format(name))
        elif exists == False:
            if not self._is_local_branch_name(name, exists=exists):
                raise ValueError("Branch `{0}` already exists, please choose a different name.".format(name))
        else:
            assert False

    def _check_remote_branch_name(self, name, exists=any):
        r"""
        Check whether ``name`` is a valid name for a remote branch, raise a
        ``ValueError`` if it is not.

        INPUT:

        same as for :meth:`_is_remote_branch_name`

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._check_remote_branch_name('')
            Traceback (most recent call last):
            ...
            ValueError: `` is not a valid branch name.
            sage: dev._check_remote_branch_name('ticket/1')

            sage: dev._check_remote_branch_name('ticket/1', exists=True)
            Traceback (most recent call last):
            ...
            ValueError: Branch `ticket/1` does not exist on the remote system.
            sage: dev._check_remote_branch_name('ticket/1', exists=False)

        """
        try:
            if not self._is_remote_branch_name(name, exists=any):
                raise ValueError
        except ValueError:
            raise ValueError("`{0}` is not a valid branch name.".format(name))

        if exists == any:
            return
        elif exists == True:
            if not self._is_remote_branch_name(name, exists=exists):
                raise ValueError("Branch `{0}` does not exist on the remote system.".format(name))
        elif exists == False:
            if not self._is_remote_branch_name(name, exists=exists):
                raise ValueError("Branch `{0}` already exists, please choose a different name.".format(name))
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
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
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
            ValueError: `master` is not a valid ticket name.

            sage: dev._UI.append("Summary: summary1\ndescription")
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
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._has_local_branch_for_ticket(1)
            False

        """
        return ticket in self.__ticket_to_branch

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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            ValueError: Branch field is not set for ticket #1 on trac.
            sage: attributes = dev1.trac._get_attributes(ticket)
            sage: attributes['branch'] = 'public/ticket/1'
            sage: dev1.trac._authenticated_server_proxy.ticket.update(ticket, "", attributes)
            'https://trac.sagemath.org/ticket/1#comment:1'
            sage: dev2._local_branch_for_ticket(ticket, download_if_not_found=True)
            Traceback (most recent call last):
            ...
            ValueError: Branch `public/ticket/1` does not exist on the remote system.

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

    def _new_local_branch_for_ticket(self, ticket):
        r"""
        Return a local branch name for ``ticket`` which does not exist yet.

        INPUT:

        - ``ticket`` -- a string or an int identifying a ticket

        TESTS::

            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._UI.append("Summary: ticket1\ndescription")
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
            sage: from sage.dev.git_interface import SUPER_SILENT
            sage: server = DoctestTracServer()
            sage: config = DoctestConfig()
            sage: config['trac']['password'] = 'secret'
            sage: dev = DoctestSageDev(config, server)
            sage: dev._pull_master_branch()
            sage: dev._chdir()

            sage: dev._UI.append("Summary: ticket1\ndescription")
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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            sage: from sage.dev.git_interface import SUPER_SILENT
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
            ValueError: Branch `ticket/1` does not exist locally.
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

            sage: dev._format_command('switch-ticket')
            'dev.switch_ticket()'
            sage: dev._format_command('switch-ticket',int(1))
            'dev.switch_ticket(1)'

        """
        try:
            __IPYTHON__
        except NameError:
            args = [str(arg) for arg in args]
            kwargs = [ "{0}={1}".format(str(key).replace("-","_"),kwargs[key]) for key in kwargs ]
            return "dev.{0}({1})".format(command.replace("-","_"), ", ".join(args+kwargs))
        else:
            args = [str(arg) for arg in args]
            kwargs = [ "--{0}={1}".format(str(key).replace("_","-"),kwargs[key]) for key in kwargs ]
            return "sage --dev {0} {1}".format(command.replace("_","-"), " ".join(args+kwargs))
