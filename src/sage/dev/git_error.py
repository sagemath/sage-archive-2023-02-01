r"""
Git errors

This module provides subclasses of ``RuntimeError`` to indicate error
conditions when calling git.

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
class GitError(RuntimeError):
    r"""
    Error raised when git exits with a non-zero exit code.

    EXAMPLES::

        sage: from sage.dev.git_error import GitError
        sage: raise GitError(128)
        Traceback (most recent call last):
        ...
        GitError: git returned with non-zero exit code (128)

    """
    def __init__(self, exit_code, cmd, stdout, stderr, advice=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.git_error import GitError
            sage: type(GitError(128))
            <class 'sage.dev.git_error.GitError'>

        """
        self.exit_code = exit_code
        self.cmd = cmd
        self.stdout = stdout
        self.stderr = stderr
        self.advice = advice
        RuntimeError.__init__(self, "git returned with non-zero exit code (%s)"%exit_code)

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
