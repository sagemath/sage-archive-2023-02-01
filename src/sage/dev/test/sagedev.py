import sage.dev.sagedev
import sage.dev.sagedev_wrapper

class DoctestSageDevWrapper(sage.dev.sagedev_wrapper.SageDevWrapper):
    def __init__(self, config, trac_server):
        sagedev = DoctestSageDev(config, trac_server)
        sage.dev.sagedev_wrapper.SageDevWrapper.__init__(self, sagedev)

        self._UI = sagedev._UI

        self._wrap("_chdir")
        self._wrap("_pull_master_branch")

class DoctestSageDev(sage.dev.sagedev.SageDev):
    def __init__(self, config, trac_server):
        from user_interface import DoctestUserInterface
        UI = DoctestUserInterface(config['UI'])
        from trac_interface import DoctestTracInterface
        trac = DoctestTracInterface(config['trac'], UI, trac_server)
        from sage.dev.git_interface import GitInterface
        config['git']['repository'] = trac_server.git._config['src']
        git = GitInterface(config['git'], UI)

        self._trac_server = trac_server

        sage.dev.sagedev.SageDev.__init__(self, config, UI, trac, git)

    def _pull_master_branch(self):
        from sage.dev.git_interface import SUPER_SILENT
        import os
        old_cwd = os.getcwd()
        self._chdir()
        try:
            from sage.dev.sagedev import MASTER_BRANCH
            self.git.fetch(SUPER_SILENT, self._trac_server.git._config['src'], "{0}:{0}".format(MASTER_BRANCH))
            self.git.checkout(SUPER_SILENT, MASTER_BRANCH)
        finally:
            os.chdir(old_cwd)

    def _chdir(self):
        import os
        os.chdir(self.config['git']['src'])

