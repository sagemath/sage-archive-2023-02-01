r"""
Git errors

This module provides subclasses of ``RuntimeError`` to indicate error
conditions when calling git.

AUTHORS:

- Julian Rueth: initial version

"""
#*****************************************************************************
#       Copyright (C) 2013 Julian Rueth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
class GitError(RuntimeError):
    r"""
    Error raised when git exits with a non-zero exit code.

    EXAMPLES::

        sage: from sage.dev.git_error import GitError
        sage: raise GitError(128, "git foo", None, None)
        Traceback (most recent call last):
        ...
        GitError: git returned with non-zero exit code (128) for "git foo".

    """
    def __init__(self, exit_code, cmd, stdout, stderr, explain=None, advice=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.git_error import GitError
            sage: type(GitError(128, "git foo", None, None))
            <class 'sage.dev.git_error.GitError'>
        """
        self.exit_code = exit_code
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr
        self.explain = explain
        self.advice = advice

        msg = ['git returned with non-zero exit code ({0}) for "{1}".'.format(exit_code, cmd)]
        if stdout:
            msg.append("output to stdout:" + "\n".join( " " + l for l in stdout.splitlines() ))
        if stderr:
            msg.append("output to stderr:" + "\n".join( " " + l for l in stderr.splitlines() ))

        RuntimeError.__init__(self, "\n".join(msg))

class DetachedHeadError(RuntimeError):
    r"""
    Error raised when a git command can not be executed because the repository
    is in a detached HEAD state.

    EXAMPLES::

        sage: from sage.dev.git_error import DetachedHeadError
        sage: raise DetachedHeadError()
        Traceback (most recent call last):
        ...
        DetachedHeadError: unexpectedly, git is in a detached HEAD state

    """
    def __init__(self):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.git_error import DetachedHeadError
            sage: type(DetachedHeadError())
            <class 'sage.dev.git_error.DetachedHeadError'>
        """
        RuntimeError.__init__(self, "unexpectedly, git is in a detached HEAD state")

class InvalidStateError(RuntimeError):
    r"""
    Error raised when a git command can not be executed because the repository
    is not in a clean state.

    EXAMPLES::

        sage: from sage.dev.git_error import InvalidStateError
        sage: raise InvalidStateError()
        Traceback (most recent call last):
        ...
        InvalidStateError: unexpectedly, git is in an unclean state

    """
    def __init__(self):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.git_error import InvalidStateError
            sage: type(InvalidStateError())
            <class 'sage.dev.git_error.InvalidStateError'>
        """
        RuntimeError.__init__(self, "unexpectedly, git is in an unclean state")
