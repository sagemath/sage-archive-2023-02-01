r"""
Wrapper for SageDev

This module provides a wrapper for :class:`sagedev.SageDev` to produce better error
handling.

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

obsolete_commands = {
    "add_comment": "comment",
    "switch_branch": "checkout",
    "download": "pull",
    "upload": "push",
    "switch_ticket": "checkout",
    "set_needs_work": "needs_work",
    "set_needs_review": "needs_review",
    "set_needs_info": "needs_info",
    "set_positive_review": "positive_review",
}

class SageDevWrapper(object):
    r"""
    Wrap a :class:`sagedev.SageDev` and its public methods.

    The need for this wrapper arises from the following problem:
    Some methods of :class:`sagedev.SageDev` call other public methods of
    :class:`sagedev.SageDev`, for example,
    :meth:`sagedev.SageDev.checkout` relies on
    :meth:`sagedev.SageDev.pull` in some cases. If an error occurs in
    :meth:`sagedev.SageDev.pull` such as an
    :class:`user_interface.OperationCancelledError`,
    :meth:`sagedev.SageDev.pull` must reraise that error so that
    :meth:`sagedev.SageDev.checkout` can do the necessary cleanup.
    However, if the user called :meth:`sagedev.SageDev.pull` directly, then
    the error should not be reraised: the user opted to cancel this operation
    and should not be bothered with an exception message. This wrapper takes
    care of getting this right.

    INPUT:

    - ``sagedev`` -- a :class:`sagedev.SageDev`

    EXAMPLES::

        sage: dev
        SageDev()

    TESTS::

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

    ``create_ticket`` silently fails for a wrapper::

        sage: UI.append("# abort")
        sage: dev.create_ticket()

    Without the wrapper an exception is raised::

        sage: UI.append("# abort")
        sage: dev._sagedev.create_ticket()
        Traceback (most recent call last):
        ...
        OperationCancelledError: ticket edit aborted

    """
    def __init__(self, sagedev):
        r"""
        Initialization.

        TESTS::

            sage: type(dev)
            <class 'sage.dev.sagedev_wrapper.SageDevWrapper'>

        """
        self._sagedev = sagedev

        self._wrap("abandon")
        self._wrap("comment")
        self._wrap("commit")
        self._wrap("create_ticket")
        self._wrap("diff")
        self._wrap("pull")
        self._wrap("download_patch")
        self._wrap("edit_ticket")
        self._wrap("gather")
        self._wrap("import_patch")
        self._wrap("local_tickets")
        self._wrap("merge")
        self._wrap("prune_closed_tickets")
        self._wrap("remote_status")
        self._wrap("reset_to_clean_working_directory")
        self._wrap("set_remote")
        self._wrap("show_dependencies")
        self._wrap("checkout")
        self._wrap("unstash")
        self._wrap("push")
        self._wrap("upload_ssh_key")
        self._wrap("vanilla")
        self._wrap("needs_work")
        self._wrap("needs_review")
        self._wrap("needs_info")
        self._wrap("positive_review")

        for old_command,new_command in obsolete_commands.items():
            self._obsolete(old_command, new_command)

        self.git = sagedev.git
        self.trac = sagedev.trac

    def _obsolete(self, old_method, new_method):
        r"""
        Create a method `old_method` which tells the user that `old_method`
        does not exist anymore and has been replaced by `new_method`.

        EXAMPLES::

            sage: dev.obsolete
            Traceback (most recent call last):
            ...
            AttributeError: 'SageDevWrapper' object has no attribute 'obsolete'
            sage: dev._obsolete("obsolete", "not_obsolete")
            sage: dev.obsolete
            <function wrapped at 0x...>

        """
        def wrap():
            from sage.misc.decorators import sage_wraps
            doc = "The command `{0}` does not exist anymore. Please use `{1}` instead.".format(self._sagedev._format_command(old_method), self._sagedev._format_command(new_method))
            def wrapped(*args, **kwargs):
                self._sagedev._UI.error(doc)
            wrapped.__doc__ = doc
            return wrapped

        setattr(self, old_method, wrap())

    def _wrap(self, method):
        r"""
        Create and register a wrapper for ``method``.

        EXAMPLES::

            sage: dev._local_branch_for_ticket
            Traceback (most recent call last):
            ...
            AttributeError: 'SageDevWrapper' object has no attribute '_local_branch_for_ticket'
            sage: dev._wrap("_local_branch_for_ticket")
            sage: dev._local_branch_for_ticket
            <function wrapped at 0x...>

        """
        from user_interface_error import OperationCancelledError
        from git_error import GitError, DetachedHeadError, InvalidStateError
        from trac_error import TracConnectionError
        from sagedev import SageDevValueError
        from user_interface import NORMAL, INFO, DEBUG
        def wrap(f):
            from sage.misc.decorators import sage_wraps
            @sage_wraps(f)
            def wrapped(*args, **kwargs):
                try:
                    return f(*args, **kwargs)
                except OperationCancelledError:
                    if self._sagedev._UI._config.get("log_level", NORMAL) >= DEBUG:
                        raise
                except GitError as e:
                    INFO_LEVEL = INFO if e.explain else NORMAL # show more info if the error was unexpected

                    self._sagedev._UI.error("GitError: git exited with a non-zero exit code ({0}).".format(e.exit_code))
                    self._sagedev._UI.show("This happened while executing `{0}`.".format(e.cmd), INFO_LEVEL)
                    self._sagedev._UI.info("I tried my best to put your working tree and repository back to its original state.")
                    if e.explain:
                        self._sagedev._UI.error(e.explain)
                    if e.advice:
                        self._sagedev._UI.info(e.advice)
                    if e.stdout is None:
                        pass
                    elif e.stdout.strip() == "":
                        self._sagedev._UI.show("git printed nothing to STDOUT.", INFO_LEVEL)
                    else:
                        self._sagedev._UI.show("git printed the following to STDOUT:\n{0}".format(e.stdout), INFO_LEVEL)
                    if e.stderr is None:
                        pass
                    elif e.stderr.strip() == "":
                        self._sagedev._UI.show("git printed nothing to STDERR.", INFO_LEVEL)
                    else:
                        self._sagedev._UI.show("git printed the following to STDERR:\n{0}".format(e.stderr), INFO_LEVEL)
                    if self._sagedev._UI._config.get("log_level", NORMAL) >= DEBUG:
                        raise
                except DetachedHeadError as e:
                    self._sagedev._UI.error("Unexpectedly your repository was found to be in a detached head state. This is probably a bug in sagedev.")
                    self._sagedev._UI.info("You can try to restore your repository to a clean state by running {0} and {1}.".format(self._sagedev._format_command("reset_to_clean_working_directory"), self._sagedev._format_command("switch",branch="master")))
                    raise
                except InvalidStateError as e:
                    self._sagedev._UI.error("Unexpectedly your repository was found to be in a non-clean state. This is probably a bug in sagedev.")
                    self._sagedev._UI.info("You can try to restore your repository to a clean state by running {0}.".format(self._sagedev._format_command("reset_to_clean_working_directory")))
                    raise
                except TracConnectionError as e:
                    self._sagedev._UI.error("Your command failed because no connection to trac could be established.")
                    if self._sagedev._UI._config.get("log_level", NORMAL) >= DEBUG:
                        raise
                except SageDevValueError as e:
                    self._sagedev._UI.error("ValueError: {0}".format(e.message))
                    if self._sagedev._UI._config.get("log_level", NORMAL) >= DEBUG:
                        raise
            return wrapped

        setattr(self, method, wrap(getattr(self._sagedev, method)))

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        TESTS::

            sage: repr(dev)
            'SageDev()'

        """
        return repr(self._sagedev)
