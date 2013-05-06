import collections
import cPickle
import functools
import os
import random
import types

from cStringIO import StringIO
from subprocess import call, check_output, CalledProcessError

from sage.env import SAGE_DOT_GIT, SAGE_REPO, DOT_SAGE

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

def normalize_ticket_name(x):
    """
    returns the normalized ticket branch name for x

    WARNING: it does not check to see if x is a valid ticket branch name

    EXAMPLES::

        sage: from sage.dev.git_interface import normalize_ticket_name
        sage: normalize_ticket_name(["t", "12345"])
        'ticket/12345'
        sage: normalize_ticket_name(["cow", "moo"])
        'ticket/moo'
        sage: normalize_ticket_name(["branch"])
        Traceback (most recent call last):
        ...
        IndexError: list index out of range
    """
    return '/'.join(('ticket', x[1]))

def load_dict_from_file(filename):
    if os.path.exists(filename):
        with open(filename) as F:
            s = F.read()
            unpickler = cPickle.Unpickler(StringIO(s))
            return unpickler.load()
    else:
        return {}

class SavingDict(collections.MutableMapping):
    def __init__(self, filename, values=None, default=None, paired=None):
        self._filename = filename
        if values is None:
            self._dict = load_dict_from_file(filename)
        else:
            self._dict = dict(values) # explicitly make copy
        if callable(default):
            self._default = default()
        else:
            self._default = default
        self._pairing = None
        if paired is not None:
            self.set_paired(paired)

    def unset_pairing(self):
        """
        unset any pairing that was constructed

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp1 = tempfile.NamedTemporaryFile().name
            sage: tmp2 = tempfile.NamedTemporaryFile().name
            sage: sd1= SavingDict(tmp1); sd2 = SavingDict(tmp2, paired=sd1)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: sd1.unset_pairing()
            sage: del sd1[0]; sd1[0]; sd2[1]
            0
            sage: os.unlink(tmp1); os.unlink(tmp2)
        """
        if self._pairing:
            with self._pairing as paired:
                paired.unset_pairing()
        self._pairing = None

    def set_paired(self, other):
        """
        set another dictionary to be updated with the reverse of this one
        and vice versa

        EXAMPLES::

            sage: from sage.dev.git_interface import SavingDict
            sage: import os, tempfile
            sage: tmp1 = tempfile.NamedTemporaryFile().name
            sage: tmp2 = tempfile.NamedTemporaryFile().name
            sage: sd1, sd2 = SavingDict(tmp1), SavingDict(tmp2)
            sage: sd1.set_paired(sd2)
            sage: sd1[0] = 1; sd1[0]; sd2[1]
            1
            0
            sage: del sd1[0]; sd1[0]; sd2[1]
            sage: sd2[2] = 3; sd1[3]; sd2[2]
            2
            3
            sage: sd2[2] = 4; sd1[3]
            sage: sd1[4]; sd2[2]
            2
            4
            sage: os.unlink(tmp1); os.unlink(tmp2)
        """
        if not isinstance(other, SavingDict):
            raise ValueError("other is not a SavingDict")

        self.unset_pairing()
        other.unset_pairing()

        class Pairing(object):
            def __init__(this, obj):
                this._obj = obj
            def __enter__(this):
                this._entered[0] += 1
                return other if this._obj is self else self
            def __exit__(this, type, value, tb):
                this._entered[0] -= 1
            def __nonzero__(this):
                return not this._entered[0]
        Pairing._entered = [0] # ints are immutable

        self._pairing, other._pairing = Pairing(self), Pairing(other)

    def _write(self):
        """
        writes self to disk
        """
        tmpfile = self._filename + '%016x'%(random.randrange(256**8))
        s = cPickle.dumps(self._dict, protocol=2)
        with open(tmpfile, 'wb') as F:
            F.write(s)
        # This move is atomic (the files are on the same filesystem)
        os.rename(tmpfile, self._filename)

    def __setitem__(self, key, value):
        if self._pairing:
            with self._pairing as paired:
                del paired[self[key]]
                paired[value] = key
        self._dict[key] = value
        self._write()

    def __delitem__(self, key):
        if self._pairing:
            with self._pairing as paired:
                del paired[self[key]]
        try:
            del self._dict[key]
        except KeyError:
            pass
        else:
            self._write()

    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            return self._default

    def __contains__(self, key):
        return key in self._dict

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._dict)

class authenticated(object):
    def __init__(self, func):
        self.func = func

    def __call__(self, git, *args, **kwargs):
        sshkeyfile = git._config.get('sshkeyfile',
                os.path.join(os.environ['HOME'], '.ssh', 'id_rsa'))

        if not "sshkey_set" in git._config or not git._config['sshkey_set']:
            git._sagedev._upload_ssh_key(sshkeyfile)
            git._config['sshkeyfile'] = sshkeyfile
            git._config['sshkey_set'] = True
            git._config._write_config()

        gitcmd = git._gitcmd

        ssh_git = "ssh -i '%s'" % sshkeyfile
        ssh_git = 'GIT_SSH="%s" ' % ssh_git
        try:
            saved_git_cmd = git._gitcmd
            git._gitcmd = ssh_git + gitcmd
            self.func(git, *args, **kwargs)
        finally:
            git._gitcmd = saved_git_cmd

class GitInterface(object):
    def __init__(self, sagedev):
        self._sagedev = sagedev
        self._UI = self._sagedev._UI

        self._config  = self._sagedev._config.get('git', {})
        self._doctest = self._config.get('doctest', None)
        self._dot_git = self._config.get('dot_git', SAGE_DOT_GIT)
        self._gitcmd  = self._config.get('gitcmd', 'git')
        self._repo    = self._config.get('repo', SAGE_REPO)

        if self._doctest is not None:
            # prepare the doctesting repository
            os.chdir(self._doctest)
            self.execute_silent('init')
            with open('testfile', 'w') as f:
                f.write('this is a test file\n')
            self.execute_silent('add', 'testfile')
            self.execute_silent('commit', '--message="add a testfile"')
            self.execute_silent('branch', 'first')
            with open('testfile', 'w') as f:
                f.write('this test file has been edited\n')
            self.execute_silent('add', 'testfile')
            self.execute_silent('commit', '--message="edit the testfile"')
            self.execute_silent('checkout', '-qb', 'second')

        if not os.path.exists(self._dot_git):
                raise ValueError("`%s` does not point to an existing directory."%self._dot_git)

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

    def __repr__(self):
        return "GitInterface()"

    def released_sage_ver(self):
        # should return a string with the most recent released version
        # of Sage (in this branch's past?)
        raise NotImplementedError

    def get_state(self):
        ret = []
        if os.path.exists(os.path.join(self._dot_git, "rebase-apply")):
            ret.append("am")
        return ret

    def reset_to_clean_state(self, interactive=True):
        states = self.get_state()
        if not states:
            return True
        for state in states:
            if state == "am":
                if interactive and not self._UI.confirm("Your repository is in an unclean state. It seems you are in the middle of a merge of some sort. To run this command you have to reset your respository to a clean state. Do you want me to reset your respository? (This will discard any changes which are not commited.)"):
                    return False

                self.execute_silent('am', '--abort')
                return self.reset_to_clean_state(interactive=interactive)
            else:
                raise NotImplementedError(state)

    def reset_to_clean_working_directory(self, interactive=True):
        if not self.has_uncommitted_changes():
            return True

        if interactive and not self._UI.confirm("You have uncommited changes in your working directory. To run this command you have to discard your changes. Do you want me to discard any changes which are not commited?"):
            return False

        self.execute_silent('reset', '--hard')

        return True

    def _clean_str(self, s):
        # for now, no error checking
        s = str(s)
        if " " in s:
            if "'" not in s:
                return "'" + s + "'"
            elif '"' not in s:
                return '"' + s + '"'
            else:
                raise RuntimeError("Quotes are too complicated")
        return s

    def _run_git(self, output_type, cmd, args, kwds):
        s = " ".join([self._gitcmd, "--git-dir=%s"%self._dot_git, cmd])
        dryrun = kwds.pop("dryrun", None)
        for k, v in kwds.iteritems():
            if len(k) == 1:
                k = ' -' + k
            else:
                k = ' --' + k
            if v is True:
                s += k
            else:
                s += k + " " + self._clean_str(v)
        if args:
            s += " " + " ".join([self._clean_str(a) for a in args if a is not None])
        if dryrun:
            return s
        else:
            if output_type == 'retval':
                return call(s, shell=True)
            elif output_type == 'retquiet':
                return call(s, shell=True, stdout=open(os.devnull, 'w'))
            elif output_type == 'stdout':
                return check_output(s, shell=True)

    def execute(self, cmd, *args, **kwds):
        return self._run_git('retval', cmd, args, kwds)

    def execute_silent(self, cmd, *args, **kwds):
        return self._run_git('retquiet', cmd, args, kwds)

    def read_output(self, cmd, *args, **kwds):
        return self._run_git('stdout', cmd, args, kwds)

    def is_child_of(self, a, b):
        """
        returns True if a is a child of b

        EXAMPLES::

            sage: sage_dev.git.is_child_of('first', 'second')
            False
            sage: sage_dev.git.is_child_of('second', 'first')
            True
            sage: sage_dev.git.is_child_of('master', 'master')
            True
        """
        return self.is_ancestor_of(b, a)

    def is_ancestor_of(self, a, b):
        """
        returns True if a is an ancestor of b

        EXAMPLES::

            sage: sage_dev.git.is_ancestor_of('first', 'second')
            True
            sage: sage_dev.git.is_ancestor_of('second', 'first')
            False
            sage: sage_dev.git.is_ancestor_of('master', 'master')
            True
        """
        revs = self.read_output('rev-list', '{}..{}'.format(b, a)).splitlines()
        return len(revs) == 0

    def has_uncommitted_changes(self):
        # Returns True if there are uncommitted changes
        return self.execute('diff', '--quiet') != 0

    def commit_all(self, *args, **kwds):
        # if files are non-tracked and user doesn't want to add any of
        # them, there might be no changes being committed here....
        kwds['a'] = True
        self.execute("commit", *args, **kwds)

    def unknown_files(self):
        # Should return a list of filenames of files that the user
        # might want to add
        status_output = self.read_output("status", porcelain=True)
        files = [line[3:] for line in status_output.splitlines()
                          if line[:2] == '??']
        return files

    def add_file(self, F):
        # Should add the file with filename F
        self.execute('add', "'{}'".format(F))

    def save(self):
        diff = self._UI.confirm("Would you like to see a diff of the changes?",
                               default_yes=False)
        if diff:
            self.execute("diff")
        added = self.files_added()
        for F in added:
            toadd = self._UI.confirm("Would you like to start tracking %s?"%F)
            if toadd:
                self.add_file(F)
        msg = self._UI.get_input("Please enter a commit message:")
        self.commit_all(m=msg)

    def local_branches(self):
        """
        Return the list of the local branches

        EXAMPLES::

            sage: sage_dev.git.local_branches()
            ['first', 'master', 'second']
        """
        result = self._run_git("stdout", "show-ref", ["--heads"], {}).split()
        result = [result[2*i+1][11:] for i in range(len(result)/2)]
        return result

    def current_branch(self):
        """
        Return the current branch

        EXAMPLES::

            sage: sage_dev.git.current_branch()
            'second'
        """
        try:
            branch = self.read_output('symbolic-ref', 'HEAD').strip()
            if branch.startswith('refs/heads/'):
                return branch[11:]
            int(branch, 16)
        except CalledProcessError:
            return None
        except ValueError:
            raise RuntimeError('HEAD is bizarre!')
        else:
            raise ValueError('HEAD is detached')

    def _ticket_to_branch(self, ticket):
        """
        Return the branch associated to an int or string.

        Returns ``None`` if no ticket is associated to this branch.
        """
        if ticket is None:
            ticket = self.current_branch()
        elif isinstance(ticket, int):
            ticket = self._branch[ticket]
        if ticket is not None and self.branch_exists(ticket):
            return ticket

    def _branch_to_ticketnum(self, branchname):
        """
        Return the ticket associated to this branch.

        Returns ``None`` if no local branch is associated to that
        ticket.
        """
        return self._ticket[branchname]

    def _branch_printname(self, branchname):
        """
        return branchname, where ticket branches are specially recognized

        EXAMPLES::

            sage: sage_dev.git._branch_printname('first')
            'first'
            sage: sage_dev.git._branch_printname('t/12345')
            '#12345'
            sage: sage_dev.git._branch_printname('ticket/12345')
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

            sage: sage_dev.git._local_to_remote_name('padics/feature')
            'g/padics/feature'
            sage: sage_dev.git._local_to_remote_name('t/12345')
            'u/doctest/ticket/12345'
            sage: sage_dev.git._local_to_remote_name('ticket/12345')
            'u/doctest/ticket/12345'
            sage: sage_dev.git._local_to_remote_name('u/doctest0/ticket/12345')
            'u/doctest0/ticket/12345'
            sage: sage_dev.git._local_to_remote_name('some/local/project')
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

            sage: sage_dev.git._remote_to_local_name('g/padics/feature')
            'padics/feature'
            sage: sage_dev.git._remote_to_local_name('u/doctest/t/12345')
            'ticket/12345'
            sage: sage_dev.git._remote_to_local_name('u/doctest/ticket/12345')
            'ticket/12345'
            sage: sage_dev.git._remote_to_local_name('u/doctest0/ticket/12345')
            'u/doctest0/ticket/12345'
            sage: sage_dev.git._remote_to_local_name('u/doctest/some/remote/project')
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
        Returns the commit id of the local branch, or None if branch does not exist.

        EXAMPLES::

            sage: git = sage_dev.git
            sage: git.branch_exists("master")    # random
            'c4512c860a162c962073a83fd08e984674dd4f44'
            sage: type(git.branch_exists("master"))
            <type 'str'>
            sage: len(git.branch_exists("master"))
            40
            sage: git.branch_exists("asdlkfjasdlf")
        """
        # TODO: optimize and make this atomic :-)
        ref = "refs/heads/%s"%branch
        if self.execute("show-ref", "--quiet", "--verify", ref):
            return None
        else:
            return self.read_output("show-ref", "--verify", ref).split()[0]

    def ref_exists(self, ref):
        raise NotImplementedError

    def create_branch(self, branchname, location=None, remote_branch=True):
        if branchname in ["t", "master", "all", "dependencies", "commit", "release"]:
            raise ValueError("Bad branchname")
        if self.branch_exists(branchname):
            raise ValueError("Branch already exists")
        move = None
        if self.has_uncommitted_changes():
            move = self._sagedev._save_uncommitted_changes()
        if location is None:
            self.branch(branchname)
        else:
            self.checkout(location, b = branchname)
        if remote_branch is True:
            remote_branch = self._local_to_remote_name(branchname)
        if remote_branch:
            self._remote[branchname] = remote_branch
        if move:
            self._sagedev._unstash_changes()

    def rename_branch(self, oldname, newname):
        self.execute("branch", oldname, newname, m=True)

    def fetch_project(self, group, branchname):
        raise NotImplementedError

    def switch_branch(self, branchname, detached = False):
        dest = move = None
        if self.has_uncommitted_changes():
            move = self._sagedev._save_uncommitted_changes()
        if detached:
            self.checkout(branchname, detach=True)
        else:
            success = self.execute_silent("checkout", branchname)
            if success != 0:
                success = self.execute_silent("branch", branchname)
                if success != 0:
                    raise RuntimeError("Could not create new branch")
                success = self.execute_silent("checkout", branchname)
                if success != 0:
                    raise RuntimeError("Failed to switch to new branch")
        if move:
            self._sagedev._unstash_changes()

    def vanilla(self, release=False):
        # switch to unstable branch in the past (release=False) or a
        # given named release
        if release is False:
            self.switch_branch("master")
        elif release is True:
            raise NotImplementedError
        else:
            release = self._validate_release_name(release)
            if not self.branch_exists(release):
                self.fetch_release(release)
            self.switch_branch(release)

    def abandon(self, branchname):
        """
        Move to trash/
        """
        trashname = "abandoned/" + branchname
        oldtrash = self.branch_exists(trashname)
        if oldtrash:
            self._UI.show("Overwriting %s in trash"(oldtrash))
        self.execute("branch", branchname, trashname, M=True)
        # Need to delete remote branch (and have a hook move it to /g/abandoned/ and update the trac symlink)
        #remotename = self._remote[branchname]

def git_cmd_wrapper(git_cmd):
    def f(self, *args, **kwds):
        return self.execute(git_cmd, *args, **kwds)
    return f

for git_cmd in ["add","am","apply","bisect","branch","checkout",
                "clean", "clone","commit","diff","fetch",
               "grep","init","log","merge",
               "mv","pull","push","rebase",
               "reset","rm","show","stash",
               "status","tag"]:
    setattr(GitInterface, git_cmd, git_cmd_wrapper(git_cmd))
