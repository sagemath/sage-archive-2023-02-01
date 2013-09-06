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

class SageDevWrapper(object):
    r"""
    Wrap a :class:`sagedev.SageDev` and its public methods.

    The need for this wrapper arises from the following problem:
    Some methods of :class:`sagedev.SageDev` call other public methods of
    :class:`sagedev.SageDev`, for example,
    :meth:`sagedev.SageDev.switch_ticket` relies on
    :meth:`sagedev.SageDev.download` in some cases. If an error occurs in
    :meth:`sagedev.SageDev.download` such as an
    :class:`user_interface.OperationCancelledError`,
    :meth:`sagedev.SageDev.download` must reraise that error so that
    :meth:`sagedev.SageDev.switch_ticket` can do the necessary cleanup.
    However, if the user called :meth:`sagedev.SageDev.download` directly, then
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
        self._wrap("add_comment")
        self._wrap("commit")
        self._wrap("create_ticket")
        self._wrap("diff")
        self._wrap("download")
        self._wrap("download_patch")
        self._wrap("edit_ticket")
        self._wrap("gather")
        self._wrap("import_patch")
        self._wrap("local_tickets")
        self._wrap("merge")
        self._wrap("prune_closed_tickets")
        self._wrap("remote_status")
        self._wrap("reset_to_clean_state")
        self._wrap("reset_to_clean_working_directory")
        self._wrap("set_remote")
        self._wrap("show_dependencies")
        self._wrap("switch_branch")
        self._wrap("switch_ticket")
        self._wrap("unstash")
        self._wrap("upload")
        self._wrap("upload_ssh_key")
        self._wrap("vanilla")
        self._wrap("set_needs_work")
        self._wrap("set_needs_review")
        self._wrap("set_needs_info")
        self._wrap("set_positive_review")

        self.git = sagedev.git
        self.trac = sagedev.trac

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
                    self._sagedev._UI.info("You can try to restore your repository to a clean state by running {0} and {1} and {2}.".format(self._sagedev._format_command("reset_to_clean_state"), self._sagedev._format_command("reset_to_clean_working_directory"), self._sagedev._format_command("vanilla")))
                    raise
                except InvalidStateError as e:
                    self._sagedev._UI.error("Unexpectedly your repository was found to be in a non-clean state. This is probably a bug in sagedev.")
                    self._sagedev._UI.info("You can try to restory your repository to a clean state by running {0} and {1}.".format(self._sagedev._format_command("reset_to_clean_state"), self._sagedev._format_command("reset_to_clean_working_directory")))
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
