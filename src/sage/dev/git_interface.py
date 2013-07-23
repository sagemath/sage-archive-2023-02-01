"""
Git Interface

This module provides a python interface to Sage's git repository.

TESTS::

    sage: from sage.dev.sagedev import SageDev, doctest_config
    sage: git = SageDev(doctest_config()).git
    sage: git.add("untracked_testfile1")
    0
    sage: git.rm('untracked_testfile1',force=True)
    rm 'untracked_testfile1'
    0
    sage: git.mv('testfile', 'new_testfile')
    0
    sage: git.checkout('HEAD', '--', 'testfile')
    0
    sage: git.reset('HEAD', '--', 'new_testfile')
    0
    sage: git.status()
    # On branch first_branch
    # Untracked files:
    #   (use "git add <file>..." to include in what will be committed)
    #
    #   new_testfile
    #   untracked_testfile2
    nothing added to commit but untracked files present (use "git add" to track)
    0
    sage: git.show('master')
    commit ...
    Author: doctest <doctest>
    Date:   Sat Mar 3 09:46:40 1973 +0000
    <BLANKLINE>
        add a testfile
    <BLANKLINE>
    diff --git a/testfile b/testfile
    new file mode 100644
    index 0000000...
    --- /dev/null
    +++ b/testfile
    @@ -0,0 +1 @@
    +this is a test file
    0
"""
import atexit
import os
import shutil
import subprocess
import tempfile

import sage.doctest
from sage.env import SAGE_DOT_GIT, SAGE_REPO_AUTHENTICATED, DOT_SAGE

def is_atomic_name(x):
    """
    returns true if x is a valid atomic branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import is_atomic_name
        sage: is_atomic_name(["branch"])
        True
        sage: is_atomic_name(["/branch"])
        False
        sage: is_atomic_name(["refs","heads","branch"])
        False
        sage: is_atomic_name([""])
        False
        sage: is_atomic_name(["u"])
        False
        sage: is_atomic_name(["1234"])
        True
    """
    if len(x) != 1:
        return False
    if '/' in x[0]:
        return False
    return x[0] not in ("t", "ticket", "u", "g", "abandoned","")

def is_ticket_name(x):
    """
    returns true if x is a valid ticket branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import is_ticket_name
        sage: is_ticket_name(["t", "1234"])
        True
        sage: is_ticket_name(["u", "doctest", "branch"])
        False
        sage: is_ticket_name(["padics", "feature"])
        False
        sage: is_ticket_name(["ticket", "myticket"])
        False
        sage: is_ticket_name(["ticket", "9876"])
        True
    """
    if len(x) != 2:
        return False
    if x[0] not in ('ticket', 't'):
        return False
    return x[1].isdigit()

def is_local_group_name(x):
    """
    returns true if x is a valid local group branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import is_local_group_name
        sage: is_local_group_name(["padic", "feature"])
        True
        sage: is_local_group_name(["g", "padic", "feature"])
        False
        sage: is_local_group_name(["padic", "feature", "1234"])
        False
        sage: is_local_group_name(["padic", "ticket", "1234"])
        True
    """
    if len(x) == 0:
        return False
    if not is_atomic_name(x[0:1]):
        return False
    if len(x) == 2:
        return is_atomic_name(x[1:])
    else:
        return is_ticket_name(x[1:])

def is_remote_group_name(x):
    """
    returns true if x is a valid remote group branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import is_remote_group_name
        sage: is_remote_group_name(["padic", "feature"])
        False
        sage: is_remote_group_name(["g", "padic", "feature"])
        True
        sage: is_remote_group_name(["g", "padic", "feature", "1234"])
        False
        sage: is_remote_group_name(["g", "padic", "ticket", "1234"])
        True
        sage: is_remote_group_name(["u", "doctest", "ticket", "1234"])
        False
    """
    if len(x) < 3:
        return False
    if x[0] != "g":
        return False
    return is_local_group_name(x[1:])

def is_remote_user_name(x):
    """
    returns true if x is a valid remote user branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import is_remote_user_name
        sage: is_remote_user_name(["u", "doctest", "ticket", "12345"])
        True
        sage: is_remote_user_name(["u", "doctest", "ticket", ""])
        False
        sage: is_remote_user_name(["u", "doctest"])
        False
        sage: is_remote_user_name(["g", "padic", "feature"])
        False
        sage: is_remote_user_name(["u", "doctest", "feature"])
        True
    """
    if len(x) < 3:
        return False
    if x[0] != "u":
        return False
    return all(x[1:])

def is_release_name(x):
    """
    returns true if x is a valid release name

    WARNING: this does not imply the existence of such a release

    EXAMPLES::

        sage: from sage.dev.git_interface import is_release_name
        sage: is_release_name(['5', '2', '7'])
        True
        sage: is_release_name(['6', '-2'])
        False
        sage: is_release_name(['6', 'beta0'])
        True
        sage: is_release_name(['7', 'rc'])
        False
        sage: is_release_name(['7', 'rc1'])
        True
    """
    for v in x[:-1]:
        try:
            if int(v) < 0:
                return False
        except ValueError:
            return False
    v = x[-1]
    if v.startswith('alpha'):
        v = v[5:]
    elif v.startswith('beta'):
        v = v[4:]
    elif v.startswith('rc'):
        v = v[2:]
    try:
        return int(v) >= 0
    except ValueError:
        return False

def normalize_ticket_name(x):
    """
    returns the normalized ticket branch name for x

    WARNING: it does not check to see if x is a valid ticket branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import normalize_ticket_name
        sage: normalize_ticket_name(["t", "12345"])
        'ticket/12345'
        sage: normalize_ticket_name(["cow", "goes", "moo"])
        'ticket/goes'
        sage: normalize_ticket_name(["branch"])
        Traceback (most recent call last):
        ...
        IndexError: list index out of range
    """
    return '/'.join(('ticket', x[1]))

class GitInterface(object):
    def __init__(self, sagedev):
        self._sagedev = sagedev
        self._UI = self._sagedev._UI

        self._config  = self._sagedev._config.get('git', {})
        self._dot_git = self._config.get('dot_git', SAGE_DOT_GIT)
        self._gitcmd  = self._config.get('gitcmd', 'git')
        self._repo    = self._config.get('repo', SAGE_REPO_AUTHENTICATED)
        self._doctest_mode = False

        if sage.doctest.DOCTEST_MODE:
            self._tmp_dir = tempfile.mkdtemp()
            atexit.register(shutil.rmtree, self._tmp_dir)
            self._dot_git = os.path.join(self._tmp_dir, '.git')
            self._doctest_mode = True
            self._prep_doctest_repo()

        if not os.path.exists(self._dot_git):
            raise ValueError("`%s` does not point to an existing directory."%self._dot_git)

    def _prep_doctest_repo(self):
        """
        creates a fake repository at self._tmp_dir for doctesting

        TESTS::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: import os
            sage: git = SageDev(doctest_config()).git  # indirect doctest
            sage: os.path.isfile(os.path.join(git._tmp_dir, 'testfile'))
            True
            sage: os.path.isfile(os.path.join(git._tmp_dir, 'untracked_testfile1'))
            True
            sage: os.path.isfile(os.path.join(git._tmp_dir, 'untracked_testfile2'))
            True
            sage: len(git.read_output('log','--oneline','first_branch').splitlines())
            2
            sage: len(git.read_output('log','--oneline','second_branch').splitlines())
            2
            sage: len(git.read_output('log','--oneline','master').splitlines())
            1
            sage: len(git.read_output('log','--oneline','test_tag').splitlines())
            1
        """
        os.chdir(self._tmp_dir)
        env = {s:'doctest' for s in
                ('GIT_COMMITTER_NAME', 'GIT_COMMITTER_EMAIL',
                 'GIT_AUTHOR_NAME',     'GIT_AUTHOR_EMAIL') }

        self.execute_silent('init')

        with open('testfile', 'w') as f:
            f.write('this is a test file\n')
        self.execute_silent('add', 'testfile')
        env['GIT_COMMITTER_DATE'] = env['GIT_AUTHOR_DATE'] = '100000000 +0000'
        self.execute_silent('commit', message="add a testfile", env=env)

        self.execute_silent('tag', 'test_tag')
        self.execute_silent('branch', 'first_branch')
        self.execute_silent('checkout', 'HEAD', quiet=True, detach=True)

        with open('testfile', 'w') as f:
            f.write('this test file has been edited\n')
        self.execute_silent('add', 'testfile')
        env['GIT_COMMITTER_DATE'] = env['GIT_AUTHOR_DATE'] = '200000000 +0000'
        self.execute_silent('commit', message="edit the testfile", env=env)

        self.execute_silent('branch', 'second_branch')
        self.execute_silent('checkout', 'first_branch', quiet=True)

        with open('testfile', 'w') as f:
            f.write('this test file has been edited differently\n')
        self.execute_silent('add', 'testfile')
        env['GIT_COMMITTER_DATE'] = env['GIT_AUTHOR_DATE'] = '300000000 +0000'
        self.execute_silent('commit',
                message="edit the testfile differently", env=env)

        with open('untracked_testfile1', 'w') as f:
            f.write('this test is untracked\n')
        with open('untracked_testfile2', 'w') as f:
            f.write('this test is also untracked\n')

    def __repr__(self):
        """
        TESTS::

            sage: repr(dev.git)
            'GitInterface()'
        """
        return "GitInterface()"

    #def released_sage_ver(self):
        # should return a string with the most recent released version
        # of Sage (in this branch's past?)
    #    raise NotImplementedError

    def get_state(self):
        """
        get the current state of merge/rebase/am/etc operations

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.execute_supersilent('merge', 'second_branch')
            1
            sage: git.get_state()
            ('merge',)
            sage: git.execute_silent('merge', abort=True)
            0
            sage: git.get_state()
            ()
            sage: git.execute_supersilent('rebase', 'second_branch')
            1
            sage: git.get_state()
            ('rebase',)
            sage: git.execute_silent('rebase', abort=True)
            0
            sage: git.get_state()
            ()
            sage: git.execute_supersilent('rebase', 'HEAD^', interactive=True,
            ....:     env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            0
            sage: git.get_state()
            ('rebase-i',)
            sage: git.execute_supersilent('merge', 'second_branch')
            1
            sage: git.get_state()
            ('merge', 'rebase-i')
            sage: git.execute_silent('rebase', abort=True)
            0
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

    def reset_to_clean_state(self, interactive=True):
        """
        gets out of a merge/am/rebase/etc state and returns True
        if successful

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.execute_supersilent('merge', 'second_branch')
            1
            sage: git.get_state()
            ('merge',)
            sage: git.reset_to_clean_state(False)
            True
            sage: git.get_state()
            ()
            sage: git.execute_supersilent('rebase', 'HEAD^', interactive=True,
            ....:     env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            0
            sage: git.execute_supersilent('merge', 'second_branch')
            1
            sage: git.get_state()
            ('merge', 'rebase-i')
            sage: git.reset_to_clean_state(False)
            True
            sage: git.get_state()
            ()
        """
        states = self.get_state()
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

        for state in states:
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

        if self.get_state():
            raise RuntimeError("failed to reset to clean state")
        return True

    def reset_to_clean_working_directory(self, interactive=True):
        r"""
        resets any changes made to the working directory and returns
        True if successful

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: import os
            sage: git = SageDev(doctest_config()).git
            sage: with open(os.path.join(git._tmp_dir, 'testfile'), 'w') as f:
            ....:     f.write('modified this file\n')
            sage: git.has_uncommitted_changes()
            True
            sage: git.reset_to_clean_working_directory(False)
            True
            sage: git.has_uncommitted_changes()
            False
        """
        if not self.has_uncommitted_changes():
            return True

        if (interactive and not
                self._UI.confirm("You have uncommited changes in your "+
                                 "working directory. To run this command you "+
                                 "have to discard your changes. Do you want "+
                                 "me to discard any changes which are not "+
                                 "commited?")):
            return False

        self.execute_silent('reset', hard=True)

        return True

    def _run_git(self, cmd, args, kwds, **ckwds):
        r"""
        common implementation for :meth:`execute`,
        :meth:`execute_silent`, :meth:`execute_supersilent`, and
        :meth:`read_output`

        INPUT:

        - ``cmd`` - git command run

        - ``args`` - extra arguments for git

        - ``kwds`` - extra keywords for git

        - ``ckwds`` - Popen like keywords but with the following changes

          - ``stdout`` - if set to `False` will supress stdout

          - ``stderr`` - if set to `False` will supress stderr

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: r = git._run_git('status', (), {})
            # On branch first_branch
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked_testfile1
            #   untracked_testfile2
            nothing added to commit but untracked files present (use "git add" to track)
            sage: r
            (0, None, None)
            sage: git._run_git('status', (), {}, stdout=False)
            (0, None, None)
            sage: git._run_git('rebase', ('HEAD^',),
            ....:     {'interactive':True,
            ....:      'env':{'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'}}, stdout=False)
            Stopped at ... edit the testfile differently
            You can amend the commit now, with
            <BLANKLINE>
                git commit --amend
            <BLANKLINE>
            Once you are satisfied with your changes, run
            <BLANKLINE>
                git rebase --continue
            <BLANKLINE>
            (0, None, None)
            sage: git._run_git('rebase', (), {'abort':True}, stdout=False, stderr=False)
            (0, None, None)
            sage: git._run_git('log', (), {'oneline':True}, stdout=str)
            (0, '... edit the testfile differently\n... add a testfile\n', None)
        """
        assert self._doctest_mode or not sage.doctest.DOCTEST_MODE, "running doctests which use git/trac is not supported from within a running session of sage"

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

        if not self._doctest_mode:
            self._UI.show("[git] %s"%(" ".join(s)))

        if ckwds.get('dryrun', False):
            return s

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
        """
        returns exit code of a git command

        INPUT:

        - ``cmd`` - git command to be run

        - ``args`` - extra arguments to supply to git command

        - ``kwds`` - extra arguments through keywords. E.g.

          ::

              graph=True                  => --graph
              message="this is a message" => --message 'this is a message'

        EXAMPLES::

            sage: r = dev.git.execute('status')
            # On branch first_branch
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked_testfile1
            #   untracked_testfile2
            nothing added to commit but untracked files present (use "git add" to track)
            sage: r
            0
            sage: _ = dev.git.execute('log', 'first_branch' , '--', 'testfile', graph=True)
            * commit ...
            | Author: doctest <doctest>
            | Date:   Thu Jul 5 05:20:00 1979 +0000
            |
            |     edit the testfile differently
            |
            * commit ...
              Author: doctest <doctest>
              Date:   Sat Mar 3 09:46:40 1973 +0000
            <BLANKLINE>
                  add a testfile
            sage: r = dev.git.execute('commit')
            # On branch first_branch
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked_testfile1
            #   untracked_testfile2
            nothing added to commit but untracked files present (use "git add" to track)
            sage: r
            1
        """
        return self._run_git(cmd, args, kwds)[0]

    __call__ = execute

    def execute_silent(self, cmd, *args, **kwds):
        """
        returns exit code of a git command while supressing stdout

        same input as :meth:`execute`

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.execute_silent('status')
            0
            sage: git.execute_silent('rebase', 'HEAD^', interactive=True,
            ....:     env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            Stopped at ... edit the testfile differently
            You can amend the commit now, with
            <BLANKLINE>
                git commit --amend
            <BLANKLINE>
            Once you are satisfied with your changes, run
            <BLANKLINE>
                git rebase --continue
            <BLANKLINE>
            0
        """
        return self._run_git(cmd, args, kwds, stdout=False)[0]

    def execute_supersilent(self, cmd, *args, **kwds):
        """
        returns exit code of a git command while supressing both stdout
        and stderr

        same input as :meth:`execute`

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.execute_supersilent('status')
            0
            sage: git.execute_supersilent('rebase', 'HEAD^', interactive=True,
            ....:     env={'GIT_SEQUENCE_EDITOR':'sed -i s+pick+edit+'})
            0
        """
        return self._run_git(cmd, args, kwds, stdout=False, stderr=False)[0]

    def read_output(self, cmd, *args, **kwds):
        r"""
        returns stdout of a git command

        same input as :meth:`execute`

        EXAMPLES::

            sage: dev.git.read_output('log', oneline=True)
            '... edit the testfile differently\n... add a testfile\n'
        """
        return self._run_git(cmd, args, kwds, stdout=str)[1]

    def is_child_of(self, a, b):
        """
        returns True if a is a child of b

        EXAMPLES::

            sage: dev.git.is_child_of('master', 'second_branch')
            False
            sage: dev.git.is_child_of('second_branch', 'master')
            True
            sage: dev.git.is_child_of('master', 'master')
            True
        """
        return self.is_ancestor_of(b, a)

    def is_ancestor_of(self, a, b):
        """
        returns True if a is an ancestor of b

        EXAMPLES::

            sage: dev.git.is_ancestor_of('master', 'second_branch')
            True
            sage: dev.git.is_ancestor_of('second_branch', 'master')
            False
            sage: dev.git.is_ancestor_of('master', 'master')
            True
        """
        revs = self.read_output('rev-list', '{}..{}'.format(b, a)).splitlines()
        return len(revs) == 0

    def has_uncommitted_changes(self):
        r"""
        returns True if there are uncommitted changes

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: import os
            sage: git = SageDev(doctest_config()).git
            sage: git.has_uncommitted_changes()
            False
            sage: with open(os.path.join(git._tmp_dir, 'testfile'), 'w') as f:
            ....:     f.write('modified this file\n')
            sage: git.has_uncommitted_changes()
            True
        """
        return self.execute('diff', quiet=True) != 0

    def commit_all(self, *args, **kwds):
        r"""
        commits all changes

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: import os
            sage: git = SageDev(doctest_config()).git
            sage: with open(os.path.join(git._tmp_dir, 'testfile'), 'w') as f:
            ....:     f.write('modified this file\n')
            sage: git.has_uncommitted_changes()
            True
            sage: git.commit_all(message="modified a file")
            [first_branch ...] modified a file
                 1 file changed, 1 insertion(+), 1 deletion(-)
            sage: git.has_uncommitted_changes()
            False
            sage: git.commit_all(message="made no changes")
            # On branch first_branch
            # Untracked files:
            #   (use "git add <file>..." to include in what will be committed)
            #
            #   untracked_testfile1
            #   untracked_testfile2
            nothing added to commit but untracked files present (use "git add" to track)
        """
        # if files are non-tracked and user doesn't want to add any of
        # them, there might be no changes being committed here....
        kwds['all'] = True
        self.execute("commit", *args, **kwds)

    def unknown_files(self):
        """
        returns the list of files that are not being tracked by git

        EXAMPLES::

            sage: dev.git.unknown_files()
            ['untracked_testfile1', 'untracked_testfile2']
        """
        return self.read_output('ls-files',
                        other=True, exclude_standard=True).splitlines()

    def save(self, interactive=True):
        """
        guided command for making a commit which includes all changes

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.save(False)
            [first_branch ...] doctesting message
             2 files changed, 2 insertions(+)
             create mode 100644 untracked_testfile1
             create mode 100644 untracked_testfile2
        """
        if (interactive and
                self._UI.confirm("Would you like to see a diff of the "+
                                 "changes?", default_yes=False)):
            self.execute("diff")
        for F in self.unknown_files():
            if (not interactive or
                    self._UI.confirm("Would you like to start tracking "+
                                     "%s?"%F)):
                self.execute('add', F)
        if interactive:
            msg = self._UI.get_input("Please enter a commit message:")
        else:
            msg = 'doctesting message'
        self.commit_all(m=msg)

    def local_branches(self):
        """
        return the list of the local branches sorted by last commit time

        EXAMPLES::

            sage: dev.git.local_branches()
            ['first_branch', 'second_branch', 'master']
        """
        result = self.read_output('for-each-ref', 'refs/heads/',
                    sort='-committerdate', format="%(refname)").splitlines()
        return [head[11:] for head in result]

    def current_branch(self):
        """
        return the current branch

        EXAMPLES::

            sage: dev.git.current_branch()
            'first_branch'
        """
        try:
            branch = self.read_output('symbolic-ref', 'HEAD').strip()
            if branch.startswith('refs/heads/'):
                return branch[11:]
            int(branch, 16)
        except ValueError:
            raise RuntimeError('HEAD is bizarre!')
        else:
            raise ValueError('HEAD is detached')

    def _branch_printname(self, branchname):
        """
        return branchname, where ticket branches are specially recognized

        EXAMPLES::

            sage: dev.git._branch_printname('first_branch')
            'first_branch'
            sage: dev.git._branch_printname('t/12345')
            '#12345'
            sage: dev.git._branch_printname('ticket/12345')
            '#12345'
        """
        if branchname.startswith('ticket/'):
            return '#' + branchname[7:]
        elif branchname.startswith('t/'):
            return '#' + branchname[2:]
        else:
            return branchname

    def _local_to_remote_name(self, branchname):
        """
        convert local branch name to 'canonical' remote branch name

        EXAMPLES::

            sage: dev.git._local_to_remote_name('padics/feature')
            'g/padics/feature'
            sage: dev.git._local_to_remote_name('t/12345')
            'u/doctest/ticket/12345'
            sage: dev.git._local_to_remote_name('ticket/12345')
            'u/doctest/ticket/12345'
            sage: dev.git._local_to_remote_name('u/doctest0/ticket/12345')
            'u/doctest0/ticket/12345'
            sage: dev.git._local_to_remote_name('some/local/project')
            'u/doctest/some/local/project'
        """
        if branchname in ("release", "beta", "master"):
            return branchname
        x = branchname.split('/')
        if is_local_group_name(x):
            if is_ticket_name(x[1:]):
                return '/'.join(('g', x[0], normalize_ticket_name(x[1:])))
            return '/'.join(['g']+x)
        elif is_ticket_name(x):
            return '/'.join(('u', self._sagedev.trac._username,
                normalize_ticket_name(x)))
        elif is_remote_user_name(x):
            return branchname
        else:
            return '/'.join(('u', self._sagedev.trac._username, branchname))

    def _remote_to_local_name(self, branchname):
        """
        convert remote branch name to 'canonical' local branch name

        EXAMPLES::

            sage: dev.git._remote_to_local_name('g/padics/feature')
            'padics/feature'
            sage: dev.git._remote_to_local_name('u/doctest/t/12345')
            'ticket/12345'
            sage: dev.git._remote_to_local_name('u/doctest/ticket/12345')
            'ticket/12345'
            sage: dev.git._remote_to_local_name('u/doctest0/ticket/12345')
            'u/doctest0/ticket/12345'
            sage: dev.git._remote_to_local_name('u/doctest/some/remote/project')
            'some/remote/project'
        """
        if branchname in ("release", "beta", "master"):
            return branchname
        x = branchname.split('/')
        if is_remote_group_name(x):
            if is_ticket_name(x[2:]):
                return '/'.join((x[1], normalize_ticket_name(x[2:])))
            return '/'.join(x[1:])
        elif is_remote_user_name(x):
            if x[1] != self._sagedev.trac._username:
                return branchname
            elif is_ticket_name(x[2:]):
                return normalize_ticket_name(x[2:])
            return '/'.join(x[2:])
        raise ValueError("not a valid remote branch name")

    def branch_exists(self, branch):
        """
        returns the commit id of the local branch, or ``None`` if
        branch does not exist

        EXAMPLES::

            sage: git = dev.git
            sage: git.branch_exists("first_branch")          # random
            '087e1fdd0fe6f4c596f5db22bc54567b032f5d2b'
            sage: int(git.branch_exists("first_branch"), 16) # random
            48484595766010161581593150175214386043155340587L
            sage: type(git.branch_exists("first_branch"))
            <type 'str'>
            sage: len(git.branch_exists("first_branch"))
            40
            sage: git.branch_exists("asdlkfjasdlf")
        """
        ref = "refs/heads/%s"%branch
        return self.ref_exists("refs/heads/%s"%branch)

    def ref_exists(self, ref):
        """
        returns the commit id of the ref, or ``None`` if
        branch does not exist

        EXAMPLES::

            sage: git = dev.git
            sage: git.ref_exists("refs/tags/test_tag")          # random
            'abdb32da3a1e50d4677e4760eda9433ac8b45414'
            sage: int(git.ref_exists("refs/tags/test_tag"), 16) # random
            981125714882459973971230819971657414365142078484L
            sage: type(git.ref_exists("refs/tags/test_tag"))
            <type 'str'>
            sage: len(git.ref_exists("refs/tags/test_tag"))
            40
            sage: git.ref_exists("refs/tags/asdlkfjasdlf")
        """
        # TODO: optimize and make this atomic :-)
        if self.execute("show-ref", ref, quiet=True, verify=True):
            return None
        else:
            return self.read_output("show-ref", ref,
                                hash=True, verify=True).strip()

    def create_branch(self, branchname, basebranch=None, remote_branch=True):
        """
        creates branch ``branchname`` based off of ``basebranch`` or the
        current branch if ``basebranch`` is ``None``

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.create_branch("third_branch")
            sage: git.branch_exists("third_branch") == git.branch_exists("first_branch")
            True
            sage: git.create_branch("fourth_branch", "second_branch")
            sage: git.branch_exists("fourth_branch") == git.branch_exists("second_branch")
            True
        """
        if branchname in ("t", "ticket", "all", "dependencies", "commit",
                "release", "beta", "master"):
            raise ValueError("bad branchname")
        if self.branch_exists(branchname):
            raise ValueError("branch already exists")

        if basebranch is None:
            ret = self.execute("branch", branchname)
        else:
            ret = self.execute("branch", branchname, basebranch)

        if remote_branch is True:
            remote_branch = self._local_to_remote_name(branchname)
        if remote_branch:
            self._sagedev._remote[branchname] = remote_branch

        if ret: # return non-zero exit codes
            return ret

    def rename_branch(self, oldname, newname):
        """
        renames branch ``oldname`` to ``newname``

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: bool(git.branch_exists("third_branch"))
            False
            sage: git.rename_branch("first_branch", "third_branch")
            sage: bool(git.branch_exists("first_branch"))
            False
            sage: bool(git.branch_exists("third_branch"))
            True
            sage: git.rename_branch("third_branch", "second_branch")
            fatal: A branch named 'second_branch' already exists.
        """
        self.execute("branch", oldname, newname, m=True)

    def fetch_project(self, group, branchname):
        raise NotImplementedError

    def switch_branch(self, branchname, detached = False):
        """
        switch to ``branchname`` in a detached state if ``detached`` is
        set to True

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.current_branch()
            'first_branch'
            sage: git.switch_branch('second_branch')
            Switched to branch 'second_branch'
            sage: git.current_branch()
            'second_branch'
            sage: git.branch_exists('third_branch')
            sage: git.switch_branch('third_branch')
            Switched to branch 'third_branch'
            sage: git.branch_exists('third_branch') # random
            '5249e7a56067e9f30244930192503d502558b6c3'
            sage: git.switch_branch('first_branch', detached=True)
            Note: checking out 'first_branch'.
            <BLANKLINE>
            You are in 'detached HEAD' state. You can look around, make experimental
            changes and commit them, and you can discard any commits you make in this
            state without impacting any branches by performing another checkout.
            <BLANKLINE>
            If you want to create a new branch to retain commits you create, you may
            do so (now or later) by using -b with the checkout command again. Example:
            <BLANKLINE>
              git checkout -b new_branch_name
            <BLANKLINE>
            HEAD is now at ... edit the testfile differently
        """
        move = None
        if self.has_uncommitted_changes():
            move = self._sagedev._save_uncommitted_changes()

        if not detached and self.branch_exists(branchname) is None:
            if self.create_branch(branchname) is not None:
                raise RuntimeError("could not create new branch")

        if self.execute("checkout", branchname, detach=detached) != 0:
            raise RuntimeError("failed to switch to new branch")

        if move:
            self._sagedev._unstash_changes()

    def vanilla(self, release=True):
        """
        switch to released version of sage
        """
        if release is False:
            release = "master"
        elif release is True:
            release = "release"
        else:
            release = str(release)
            if is_release_name(release.split('.')):
                self.execute('fetch', 'origin', tags=True)
                release = self.ref_exists('refs/tags/%s'%release)
                if release is None:
                    raise ValueError("was unable to find desired release")
        self.switch_branch(release)

    def abandon(self, branchname):
        """
        move branch to trash

        EXAMPLES::

            sage: from sage.dev.sagedev import SageDev, doctest_config
            sage: git = SageDev(doctest_config()).git
            sage: git.abandon('second_branch')
            sage: git.local_branches()
            ['first_branch', 'abandoned/second_branch', 'master']
        """
        trashname = "abandoned/" + branchname
        oldtrash = self.branch_exists(trashname)
        if oldtrash is not None:
            self._UI.show("Overwriting %s in trash"%oldtrash)
        self.rename_branch(branchname, trashname)
        # Need to delete remote branch (and have a hook move it to /g/abandoned/ and update the trac symlink)
        #remotename = self._remote[branchname]

def _git_cmd_wrapper(git_cmd):
    """
    creates a method for GitInterface that wraps a git command

    EXAMPLES::

        sage: from sage.dev.git_interface import _git_cmd_wrapper, GitInterface
        sage: cmd = _git_cmd_wrapper("ls-tree")
        sage: setattr(GitInterface, 'ls_tree', cmd)
        sage: r = dev.git.ls_tree('first_branch')
        100644 blob 13d5431ac2f249eb07313624ec8fa041ea0f34a4        testfile
        sage: r
        0
    """
    def meth(self, *args, **kwds):
        return self.execute(git_cmd.replace("_", "-"), *args, **kwds)
    meth.__doc__ = """
            direct call to \`git %s\`

            see :meth:`execute` for full documentation
            """%git_cmd
    return meth

for git_cmd in (
        "add",
        "am",
        "apply",
        "bisect",
        "branch",
        "checkout",
        "clean",
        "clone",
        "commit",
        "diff",
        "fetch",
        "format_patch",
        "grep",
        # "init",
        "log",
        "merge",
        "mv",
        "pull",
        "push",
        "rebase",
        "reset",
        "rm",
        "show",
        "stash",
        "status",
        # "tag"
        ):
    setattr(GitInterface, git_cmd, _git_cmd_wrapper(git_cmd))
