r"""
SageDev objects for doctesting

This module provides special versions of :class:`sage.dev.sagedev.SageDev` and
:class:`sage.dev.sagedev_wrapper.SageDevWrapper` which are suitable for
doctesting.

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

import sage.dev.sagedev
import sage.dev.sagedev_wrapper

class DoctestSageDevWrapper(sage.dev.sagedev_wrapper.SageDevWrapper):
    r"""
    A :class:`sage.dev.sagedev_wrapper.SageDevWrapper` for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: DoctestSageDevWrapper(DoctestConfig(), DoctestTracServer())
        SageDev()

    """
    def __init__(self, config, trac_server):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.sagedev import DoctestSageDevWrapper
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: type(DoctestSageDevWrapper(DoctestConfig(), DoctestTracServer()))
            <class 'sage.dev.test.sagedev.DoctestSageDevWrapper'>

        """
        sagedev = DoctestSageDev(config, trac_server)
        sage.dev.sagedev_wrapper.SageDevWrapper.__init__(self, sagedev)

        self._UI = sagedev._UI

        self._wrap("_chdir")
        self._wrap("_pull_master_branch")

class DoctestSageDev(sage.dev.sagedev.SageDev):
    r"""
    A :class:`sage.dev.sagedev.SageDev` for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.sagedev import DoctestSageDev
        sage: from sage.dev.test.config import DoctestConfig
        sage: from sage.dev.test.trac_server import DoctestTracServer
        sage: DoctestSageDev(DoctestConfig(), DoctestTracServer())
        SageDev()

    """
    def __init__(self, config, trac_server):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: type(DoctestSageDev(DoctestConfig(), DoctestTracServer()))
            <class 'sage.dev.test.sagedev.DoctestSageDev'>

        """
        from user_interface import DoctestUserInterface
        UI = DoctestUserInterface(config['UI'])
        from trac_interface import DoctestTracInterface
        trac = DoctestTracInterface(config['trac'], UI, trac_server)
        from sage.dev.git_interface import GitInterface
        config['git']['repository_anonymous'] = config['git']['repository'] = trac_server.git._config['src']
        git = GitInterface(config['git'], UI)

        self._trac_server = trac_server

        sage.dev.sagedev.SageDev.__init__(self, config, UI, trac, git)

    def _pull_master_branch(self):
        r"""
        Pull the master branch of the repository of the
        :class:`trac_server.DoctestTracServer` into the local repository.

        EXAMPLES::

            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: dev = DoctestSageDev(DoctestConfig(), DoctestTracServer())
            sage: dev._pull_master_branch()

        """
        import os
        old_cwd = os.getcwd()
        self._chdir()
        try:
            from sage.dev.sagedev import MASTER_BRANCH
            if MASTER_BRANCH != "master":
                self.git.super_silent.chechkout("-b",MASTER_BRANCH)
                self.git.super_silent.checkout(MASTER_BRANCH)
            self.git.super_silent.pull(self._trac_server.git._config['src'], MASTER_BRANCH)
        finally:
            os.chdir(old_cwd)

    def _chdir(self):
        r"""
        Change the current working directory to the directory of the git
        repository of this object.

        EXAMPLES::

            sage: from sage.dev.test.sagedev import DoctestSageDev
            sage: from sage.dev.test.config import DoctestConfig
            sage: from sage.dev.test.trac_server import DoctestTracServer
            sage: dev = DoctestSageDev(DoctestConfig(), DoctestTracServer())
            sage: dev._chdir()

        """
        import os
        os.chdir(self.config['git']['src'])

def single_user_setup():
    r"""
    Create a typical single user setup for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.sagedev import single_user_setup
        sage: dev, config, UI, server = single_user_setup()

    """
    from trac_server import DoctestTracServer
    from config import DoctestConfig
    server = DoctestTracServer()
    config = DoctestConfig()
    config['trac']['password'] = 'secret'
    dev = DoctestSageDevWrapper(config, server)
    dev._pull_master_branch()
    dev._chdir()
    return dev, config, dev._UI, server

def two_user_setup():
    r"""
    Create a typical two user setup for doctesting.

    EXAMPLES::

        sage: from sage.dev.test.sagedev import two_user_setup
        sage: alice, alice_config, bob, bob_config, server = two_user_setup()

    """
    from trac_server import DoctestTracServer
    from config import DoctestConfig
    server = DoctestTracServer()
    config_alice = DoctestConfig('alice')
    config_alice['trac']['password'] = 'secret'
    alice = DoctestSageDevWrapper(config_alice, server)
    alice._pull_master_branch()

    config_bob = DoctestConfig('bob')
    config_bob['trac']['password'] = 'secret'
    bob = DoctestSageDevWrapper(config_bob, server)
    bob._pull_master_branch()

    return alice, config_alice, bob, config_bob, server
