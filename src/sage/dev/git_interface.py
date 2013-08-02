r"""
Git Interface

This module provides a python interface to Sage's git repository.

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

import os

from sage.env import SAGE_DOT_GIT, SAGE_REPO_AUTHENTICATED, SAGE_SRC

from git_error import GitError, DetachedHeadError

SILENT = object()
SUPER_SILENT = object()
READ_OUTPUT = object()

class GitInterface(object):
    r"""
    A wrapper around the ``git`` command line tool.

    Most methods of this class correspond to actual git commands. Some add
    functionality which is not directly available in git. However, all of the
    methods should be non-interactive. If interaction is required the method
    should live in :class:`saged.dev.sagedev.SageDev`.

    EXAMPLES::

        sage: from sage.dev.config import Config
        sage: from sage.dev.user_interface import UserInterface
        sage: from sage.dev.git_interface import GitInterface
        sage: config = Config()
        sage: GitInterface(config, UserInterface(config))
        GitInterface()

    """
    def __init__(self, config, UI):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.config import Config
            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.git_interface import GitInterface
            sage: config = Config()
            sage: type(GitInterface(config, UserInterface(config)))
            <class 'sage.dev.git_interface.GitInterface'>

        """
        self._UI = UI
        self._config = config

        self._src = self._config.get('src', SAGE_SRC)
        self._dot_git = self._config.get('dot_git', SAGE_DOT_GIT)
        self._gitcmd = self._config.get('gitcmd', 'git')
        self._repository = self._config.get('repository', SAGE_REPO_AUTHENTICATED)

        if not os.path.exists(self._dot_git):
            raise ValueError("`%s` does not point to an existing directory."%self._dot_git)

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: from sage.dev.config import Config
            sage: from sage.dev.user_interface import UserInterface
            sage: from sage.dev.git_interface import GitInterface
            sage: config = Config()
            sage: repr(GitInterface(config, UserInterface(config)))
            'GitInterface()'

        """
        return "GitInterface()"

    def get_state(self):
        r"""
        Get the current state of merge/rebase/am/etc operations.

        OUTPUT:

        A tuple of strings which consists of any of the following:
        ``'rebase'``, ``'am'``, ``'rebase-i'``, ``'rebase-m'``, ``'merge'``,
        ``'bisect'``, ``'cherry-seq'``, ``'cherry'``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.add("file")
            sage: git.commit(SUPER_SILENT, "-m","initial commit")
            sage: git.checkout(SUPER_SILENT, "-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.commit(SUPER_SILENT, "-am","second commit")
            sage: git.checkout(SUPER_SILENT, "master")
            sage: git.checkout(SUPER_SILENT, "-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.commit(SUPER_SILENT, "-am","conflicting commit")

        A ``merge`` state::

            sage: git.checkout(SUPER_SILENT, "branch1")
            sage: git.merge(SUPER_SILENT, 'branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1)
            sage: git.get_state()
            ('merge',)
            sage: git.merge(SUPER_SILENT,abort=True)
            sage: git.get_state()
            ()

        A ``rebase`` state::

            sage: git.execute_supersilent('rebase', 'branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1)
            sage: git.get_state()
            ('rebase',)
            sage: git.rebase(SUPER_SILENT, abort=True)
            sage: git.get_state()
            ()

        A merge within an interactive rebase::

            sage: git.rebase(SUPER_SILENT, 'HEAD^', interactive=True, env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            sage: git.get_state()
            ('rebase-i',)
            sage: git.merge(SUPER_SILENT, 'branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1)
            sage: git.get_state()
            ('merge', 'rebase-i')
            sage: git.rebase(SUPER_SILENT, abort=True)
            sage: git.get_state()
            ()

        """
        # logic based on zsh's git backend for vcs_info
        opj = os.path.join
        p = lambda x: opj(self._dot_git, x)
        ret = []
        for d in map(p,("rebase-apply", "rebase", opj("..",".dotest"))):
            if os.path.isdir(d):
                if os.path.isfile(opj(d, 'rebasing')) and 'rebase' not in ret:
                    ret.append('rebase')
                if os.path.isfile(opj(d, 'applying')) and 'am' not in ret:
                    ret.append('am')
        for f in map(p, (opj('rebase-merge', 'interactive'),
                         opj('.dotest-merge', 'interactive'))):
            if os.path.isfile(f):
                ret.append('rebase-i')
                break
        else:
            for d in map(p, ('rebase-merge', '.dotest-merge')):
                if os.path.isdir(d):
                    ret.append('rebase-m')
                    break
        if os.path.isfile(p('MERGE_HEAD')):
            ret.append('merge')
        if os.path.isfile(p('BISECT_LOG')):
            ret.append('bisect')
        if os.path.isfile(p('CHERRY_PICK_HEAD')):
            if os.path.isdir(p('sequencer')):
                ret.append('cherry-seq')
            else:
                ret.append('cherry')
        # return in reverse order so reset operations are correctly ordered
        return tuple(reversed(ret))

    def reset_to_clean_state(self):
        r"""
        Get out of a merge/am/rebase/etc state.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.add("file")
            sage: git.commit(SUPER_SILENT, "-m","initial commit")
            sage: git.checkout(SUPER_SILENT, "-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.commit(SUPER_SILENT, "-am","second commit")
            sage: git.checkout(SUPER_SILENT, "master")
            sage: git.checkout(SUPER_SILENT, "-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.commit(SUPER_SILENT, "-am","conflicting commit")

        A merge within an interactive rebase::

            sage: git.checkout(SUPER_SILENT, "branch1")
            sage: git.rebase(SUPER_SILENT, 'HEAD^', interactive=True, env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            sage: git.get_state()
            ('rebase-i',)
            sage: git.merge(SUPER_SILENT, 'branch2')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (1)
            sage: git.get_state()
            ('merge', 'rebase-i')

        Get out of this state::

            sage: git.reset_to_clean_state()
            sage: git.get_state()
            ()

        """
        states = self.get_state()
        if not states:
            return

        state = states[0]
        if state.startswith('rebase'):
            self.execute_silent('rebase', abort=True)
        elif state == 'am':
            self.execute_silent('am', abort=True)
        elif state == 'merge':
            self.execute_silent('merge', abort=True)
        elif state == 'bisect':
            raise NotImplementedError(state)
        elif state.startswith('cherry'):
            self.execute_silent('cherry-pick', abort=True)
        else:
            raise RuntimeError("'%s' is not a valid state"%state)

        return self.reset_to_clean_state()

    def reset_to_clean_working_directory(self, remove_untracked_files=False, remove_untracked_directories=False, remove_ignored=False):
        r"""
        Reset any changes made to the working directory.

        INPUT:

        - ``remove_untracked_files`` -- a boolean (default: ``False``), whether
          to remove files which are not tracked by git

        - ``remove_untracked_directories`` -- a boolean (default: ``False``),
          whether to remove directories which are not tracked by git

        - ``remove_ignored`` -- a boolean (default: ``False``), whether to
          remove files directories which are ignored by git

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Set up some files/directories::

            sage: os.chdir(config['git']['src'])
            sage: open('tracked','w').close()
            sage: git.add(SUPER_SILENT, 'tracked')
            sage: with open('.gitignore','w') as f: f.write('ignored\nignored_dir')
            sage: git.add(SUPER_SILENT, '.gitignore')
            sage: git.commit(SUPER_SILENT, '-m', 'initial commit')

            sage: os.mkdir('untracked_dir')
            sage: open('untracked_dir/untracked','w').close()
            sage: open('untracked','w').close()
            sage: open('ignored','w').close()
            sage: os.mkdir('ignored_dir')
            sage: open('ignored_dir/untracked','w').close()
            sage: with open('tracked','w') as f: f.write('version 0')
            sage: git.status()
            # On branch master
            # Changes not staged for commit:
            #   (use "git add <file>..." to update what will be committed)
            #   (use "git checkout -- <file>..." to discard changes in working directory)
            #
            #   modified:   tracked
            #
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            #   untracked_dir/
            no changes added to commit (use "git add" and/or "git commit -a")

        Some invalid combinations of flags::

            sage: git.reset_to_clean_working_directory(remove_untracked_files = False, remove_untracked_directories = True)
            Traceback (most recent call last):
            ...
            ValueError: remove_untracked_directories only valid if remove_untracked_files is set
            sage: git.reset_to_clean_working_directory(remove_untracked_files = False, remove_ignored = True)
            Traceback (most recent call last):
            ...
            ValueError: remove_ignored only valid if remove_untracked_files is set

        Per default only the tracked modified files are reset to a clean
        state::

            sage: git.reset_to_clean_working_directory()
            sage: git.status()
            # On branch master
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked
            #   untracked_dir/
            nothing added to commit but untracked files present (use "git add" to track)

        Untracked items can be removed by setting the parameters::

            sage: git.reset_to_clean_working_directory(remove_untracked_files=True)
            Removing untracked
            Not removing untracked_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_untracked_directories=True)
            Removing untracked_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_ignored=True)
            Removing ignored
            Not removing ignored_dir/
            sage: git.reset_to_clean_working_directory(remove_untracked_files=True, remove_untracked_directories=True, remove_ignored=True)
            Removing ignored_dir/

        """
        if remove_untracked_directories and not remove_untracked_files:
            raise ValueError("remove_untracked_directories only valid if remove_untracked_files is set")
        if remove_ignored and not remove_untracked_files:
            raise ValueError("remove_ignored only valid if remove_untracked_files is set")

        self.reset(SILENT, hard=True)

        if remove_untracked_files:
            switches = ['-f']
            if remove_untracked_directories: switches.append("-d")
            if remove_ignored: switches.append("-x")
            self.clean(*switches)

    def _run_git(self, cmd, args, kwds, **ckwds):
        r"""
        Common implementation for :meth:`execute`, :meth:`execute_silent`,
        :meth:`execute_supersilent`, and :meth:`read_output`

        INPUT:

        - ``cmd`` - git command run

        - ``args`` - extra arguments for git

        - ``kwds`` - extra keywords for git

        - ``ckwds`` - Popen like keywords but with the following changes

          - ``stdout`` - if set to ``False`` will supress stdout

          - ``stderr`` - if set to ``False`` will supress stderr

        .. WARNING::

            This method does not raise an exception if the git call returns a
            non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git._run_git('status', (), {})
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)
            (0, None, None)
            sage: git._run_git('status', (), {}, stdout=False)
            (0, None, None)

        TESTS:

        Check that we refuse to touch the live source code in doctests::

            sage: dev.git.status()
            Traceback (most recent call last):
            ...
            AssertionError: possible attempt to work with the live repository/directory in a doctest - did you forget to dev._chdir()?

        """
        import sage.doctest
        import os
        assert not sage.doctest.DOCTEST_MODE or (self._dot_git != SAGE_DOT_GIT and self._repository != SAGE_REPO_AUTHENTICATED and os.path.abspath(os.getcwd()).startswith(self._src)), "possible attempt to work with the live repository/directory in a doctest - did you forget to dev._chdir()?"

        # not sure which commands could possibly create a commit object with
        # some crazy flags set - these commands should be safe
        if cmd not in [ "config", "diff", "grep", "log", "ls_remote", "remote", "reset", "show", "show_ref", "status", "symbolic_ref" ]:
            self._check_user_email()

        s = [self._gitcmd, "--git-dir=%s"%self._dot_git, cmd]

        env = ckwds.setdefault('env', dict(os.environ))
        env.update(kwds.pop('env', {}))

        for k, v in kwds.iteritems():
            if len(k) == 1:
                k = '-' + k
            else:
                k = '--' + k.replace('_', '-')
            if v is True:
                s.append(k)
            elif v is not False:
                s.extend((k, v))
        if args:
            s.extend(a for a in args if a is not None)

        s = [str(arg) for arg in s]

        from sage.dev.user_interface import INFO
        self._UI.show("[git] %s"%(" ".join(s)), log_level=INFO)

        if ckwds.get('dryrun', False):
            return s

        import subprocess
        devnull = open(os.devnull, 'w')
        if ckwds.get('stdout') is False:
            ckwds['stdout'] = devnull
        elif ckwds.get('stdout') is str:
            ckwds['stdout'] = subprocess.PIPE
        if ckwds.get('stderr') is False:
            ckwds['stderr'] = devnull
        elif ckwds.get('stderr') is str:
            ckwds['stderr'] = subprocess.PIPE
        process = subprocess.Popen(s, **ckwds)
        stdout, stderr = process.communicate()
        retcode = process.poll()
        return retcode, stdout, stderr

    def execute(self, cmd, *args, **kwds):
        r"""
        Run git.

        Raises an exception if git has non-zero exit code.

        INPUT:

        - ``cmd`` - git command run

        - ``args`` - extra arguments for git

        - ``kwds`` - extra keywords for git

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git.execute('status')
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)
            sage: git.execute_silent('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129)

        """
        exit_code = self._run_git(cmd, args, kwds)[0]
        if exit_code:
            raise GitError(exit_code)

    __call__ = execute

    def execute_silent(self, cmd, *args, **kwds):
        r"""
        Run git and supress its output to stdout.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git.execute_silent('status')
            sage: git.execute_silent('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129)

        """
        exit_code = self._run_git(cmd, args, kwds, stdout=False)[0]
        if exit_code:
            raise GitError(exit_code)

    def execute_supersilent(self, cmd, *args, **kwds):
        r"""
        Run git and supress its output to stdout and stderr.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git.execute_supersilent('status')
            sage: git.execute_supersilent('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129)

        """
        exit_code = self._run_git(cmd, args, kwds, stdout=False, stderr=False)[0]
        if exit_code:
            raise GitError(exit_code)

    def read_output(self, cmd, *args, **kwds):
        r"""
        Run git and return its output to stdout.

        Same input as :meth:`execute`.

        Raises an error if git returns a non-zero exit code.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])

            sage: git.read_output('status')
            '# On branch master\n#\n# Initial commit\n#\nnothing to commit (create/copy files and use "git add" to track)\n'
            sage: git.read_output('status',foo=True) # --foo is not a valid parameter
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (129)

        """
        exit_code, ret, _ = self._run_git(cmd, args, kwds, stdout=str, stderr=False)
        if exit_code:
            raise GitError(exit_code)
        return ret

    def is_child_of(self, a, b):
        r"""
        Return whether ``a`` is a child of ``b``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.add("file")
            sage: git.commit(SUPER_SILENT, "-m","initial commit")
            sage: git.checkout(SUPER_SILENT, "-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.commit(SUPER_SILENT, "-am","second commit")
            sage: git.checkout(SUPER_SILENT, "master")
            sage: git.checkout(SUPER_SILENT, "-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.commit(SUPER_SILENT, "-am","conflicting commit")

            sage: git.is_child_of('master', 'branch2')
            False
            sage: git.is_child_of('branch2', 'master')
            True
            sage: git.is_child_of('branch1', 'branch2')
            False
            sage: git.is_child_of('master', 'master')
            True

        """
        return self.is_ancestor_of(b, a)

    def is_ancestor_of(self, a, b):
        r"""
        Return whether ``a`` is an ancestor of ``b``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create two conflicting branches::

            sage: os.chdir(config['git']['src'])
            sage: with open("file","w") as f: f.write("version 0")
            sage: git.add("file")
            sage: git.commit(SUPER_SILENT, "-m","initial commit")
            sage: git.checkout(SUPER_SILENT, "-b","branch1")
            sage: with open("file","w") as f: f.write("version 1")
            sage: git.commit(SUPER_SILENT, "-am","second commit")
            sage: git.checkout(SUPER_SILENT, "master")
            sage: git.checkout(SUPER_SILENT, "-b","branch2")
            sage: with open("file","w") as f: f.write("version 2")
            sage: git.commit(SUPER_SILENT, "-am","conflicting commit")

            sage: git.is_ancestor_of('master', 'branch2')
            True
            sage: git.is_ancestor_of('branch2', 'master')
            False
            sage: git.is_ancestor_of('branch1', 'branch2')
            False
            sage: git.is_ancestor_of('master', 'master')
            True

        """
        return not self.rev_list(READ_OUTPUT, '{}..{}'.format(b, a)).splitlines()

    def has_uncommitted_changes(self):
        r"""
        Return whether there are uncommitted changes, i.e., whether there are
        modified files which are tracked by git.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        An untracked file does not count towards uncommited changes::

            sage: os.chdir(config['git']['src'])
            sage: open('untracked','w').close()
            sage: git.has_uncommitted_changes()
            False

        Once added to the index it does::

            sage: git.add('untracked')
            sage: git.has_uncommitted_changes()
            True
            sage: git.commit(SUPER_SILENT, '-m', 'tracking untracked')
            sage: git.has_uncommitted_changes()
            False
            sage: with open('untracked','w') as f: f.write('version 0')
            sage: git.has_uncommitted_changes()
            True

        """
        return bool([line for line in self.status(READ_OUTPUT, porcelain=True).splitlines() if not line.startswith('?')])

    def untracked_files(self):
        r"""
        Return a list of file names for files that are not tracked by git and
        not ignored.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        An untracked file::

            sage: os.chdir(config['git']['src'])
            sage: git.untracked_files()
            []
            sage: open('untracked','w').close()
            sage: git.untracked_files()
            ['untracked']

         Directories are not displayed here::

            sage: os.mkdir('untracked_dir')
            sage: git.untracked_files()
            ['untracked']
            sage: open('untracked_dir/untracked','w').close()
            sage: git.untracked_files()
            ['untracked', 'untracked_dir/untracked']

        """
        return self.read_output('ls-files', other=True, exclude_standard=True).splitlines()

    def local_branches(self):
        r"""
        Return a list of local branches sorted by last commit time.

        EXAMPLES::

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git.branch('branch1')
            sage: git.branch('branch2')

        Use this repository as a remote repository::

            sage: config2 = DoctestConfig()
            sage: git2 = GitInterface(config2["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config2['git']['src'])
            sage: git2.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git2.remote('add', 'git', config['git']['src'])
            sage: git2.fetch(SUPER_SILENT, 'git')
            sage: git2.checkout(SUPER_SILENT, "branch1")
            sage: git2.branch("-a")
            * branch1
              master
              remotes/git/branch1
              remotes/git/branch2
              remotes/git/master

            sage: git2.local_branches()
            ['branch1', 'master']
            sage: os.chdir(config['git']['src'])
            sage: git.local_branches()
            ['branch1', 'branch2', 'master']

        """
        result = self.for_each_ref(READ_OUTPUT, 'refs/heads/',
                    sort='-committerdate', format="%(refname)").splitlines()
        return [head[11:] for head in result]

    def current_branch(self):
        r"""
        Return the current branch

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git.commit(SILENT, '-m','second commit','--allow-empty')
            sage: git.branch('branch1')
            sage: git.branch('branch2')

            sage: git.current_branch()
            'master'
            sage: git.checkout(SUPER_SILENT, 'branch1')
            sage: git.current_branch()
            'branch1'

        If ``HEAD`` is detached::

            sage: git.checkout(SUPER_SILENT, 'master~')
            sage: git.current_branch()
            Traceback (most recent call last):
            ...
            DetachedHeadError: unexpectedly, git is in a detached HEAD state

        """
        try:
            return self.symbolic_ref(READ_OUTPUT, 'HEAD', short=True, quiet=True).strip()
        except GitError as e:
            if e.exit_code == 1:
               raise DetachedHeadError()
            raise

    def commit_for_branch(self, branch):
        r"""
        Return the commit id of the local ``branch``, or ``None`` if the branch
        does not exist

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git.branch('branch1')
            sage: git.branch('branch2')

        Check existence of branches::

            sage: git.commit_for_branch('branch1') # random output
            '087e1fdd0fe6f4c596f5db22bc54567b032f5d2b'
            sage: git.commit_for_branch('branch2') is not None
            True
            sage: git.commit_for_branch('branch3') is not None
            False

        """
        return self.commit_for_ref("refs/heads/%s"%branch)

    def commit_for_ref(self, ref):
        r"""
        Return the commit id of the ``ref``, or ``None`` if the ``ref`` does
        not exist.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git.branch('branch1')
            sage: git.branch('branch2')

        Check existence of branches::

            sage: git.commit_for_ref('refs/heads/branch1') # random output
            '087e1fdd0fe6f4c596f5db22bc54567b032f5d2b'
            sage: git.commit_for_ref('refs/heads/branch2') is not None
            True
            sage: git.commit_for_ref('refs/heads/branch3') is not None
            False

        """
        try:
            return self.show_ref(READ_OUTPUT, ref, hash=True, verify=True).strip()
        except GitError:
            return None

    def rename_branch(self, oldname, newname):
        r"""
        Rename ``oldname`` to ``newname``.

        EXAMPLES:

        Create a :class:`GitInterface` for doctesting::

            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))

        Create some branches::

            sage: os.chdir(config['git']['src'])
            sage: git.commit(SILENT, '-m','initial commit','--allow-empty')
            sage: git.branch('branch1')
            sage: git.branch('branch2')

        Rename some branches::

            sage: git.rename_branch('branch1', 'branch3')
            sage: git.rename_branch('branch2', 'branch3')
            Traceback (most recent call last):
            ...
            GitError: git returned with non-zero exit code (128)

        """
        self.branch(oldname, newname, move=True)

    def _check_user_email(self):
        r"""
        Make sure that a real name and an email are set for git. These will
        show up next to any commit that user creates.

        TESTS::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: del config['git']['user.name']
            sage: del config['git']['user.email']
            sage: UI = DoctestUserInterface(config["UI"])
            sage: git = GitInterface(config["git"], UI)
            sage: os.chdir(config['git']['src'])
            sage: UI.append("Doc Test")
            sage: UI.append("doc@test")
            sage: git._check_user_email()

        """
        try:
            self.config(SUPER_SILENT, "user.name")
        except GitError as e:
            if e.exit_code == 1:
                self._UI.normal("No real name has been set for git. This name shows up as the author for any commits you contribute to sage.")
                name = self._UI.question("Your real name:")
                self.git.config("user.name",name,local=True,add=True)
                self._UI.info("Your real name has been saved.")
            else:
                raise

        try:
            self.config(SUPER_SILENT, "user.email")
        except GitError as e:
            if e.exit_code == 1:
                self._UI.normal("No email address has been set for git. This email shows up as the author for any commits you contribute to sage.")
                email = self._UI.question("Your email address:")
                self.git.config("user.email",email,local=True,add=True)
                self._UI.info("Your email has been saved.")
            else:
                raise

for git_cmd_ in (
        "add",
        "am",
        "apply",
        "bisect",
        "branch",
        "config",
        "checkout",
        "clean",
        "clone",
        "commit",
        "diff",
        "fetch",
        "for_each_ref",
        "format_patch",
        "grep",
        "init",
        "log",
        "ls_remote",
        "merge",
        "mv",
        "pull",
        "push",
        "rebase",
        "remote",
        "reset",
        "rev_list",
        "rm",
        "show",
        "show_ref",
        "stash",
        "status",
        "symbolic_ref",
        "tag"
        ):
    def create_wrapper(git_cmd__):
        r"""
        Create a wrapper for `git_cmd__`.

        EXAMPLES::

            sage: import os
            sage: from sage.dev.git_interface import GitInterface, SILENT, SUPER_SILENT
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.user_interface import DoctestUserInterface
            sage: config = DoctestConfig()
            sage: git = GitInterface(config["git"], DoctestUserInterface(config["UI"]))
            sage: os.chdir(config['git']['src'])
            sage: git.status()
            # On branch master
            #
            # Initial commit
            #
            nothing to commit (create/copy files and use "git add" to track)

        """
        git_cmd = git_cmd__.replace("_","-")
        def meth(self, *args, **kwds):
            args = list(args)
            if len([arg for arg in args if arg in (SILENT, SUPER_SILENT, READ_OUTPUT)]) > 1:
                raise ValueError("at most one of SILENT, SUPER_SILENT, READ_OUTPUT allowed")
            if SILENT in args:
                args.remove(SILENT)
                return self.execute_silent(git_cmd, *args, **kwds)
            elif SUPER_SILENT in args:
                args.remove(SUPER_SILENT)
                return self.execute_supersilent(git_cmd, *args, **kwds)
            elif READ_OUTPUT in args:
                args.remove(READ_OUTPUT)
                return self.read_output(git_cmd, *args, **kwds)
            else:
                return self.execute(git_cmd, *args, **kwds)
        meth.__doc__ = r"""
        Call `git {0}`.

        If `args` contains ``SILENT``, then output to stdout is supressed.

        If `args` contains ``SUPER_SILENT``, then output to stdout and stderr
        is supressed.

        OUTPUT:

        Returns ``None`` unless `args` contains ``READ_OUTPUT``; in that case,
        the commands output to stdout is returned.

        See :meth:`execute` for more information.

        EXAMPLES:

            sage: dev.git.{1}() # not tested

        """.format(git_cmd, git_cmd__)
        return meth
    setattr(GitInterface, git_cmd_, create_wrapper(git_cmd_))
