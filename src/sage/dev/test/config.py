import atexit
import os

import sage.dev.config

class Config(sage.dev.config.Config):
    def __init__(self, trac_username = "doctest", repository=None):
        devrc = tempfile.mkstemp()[1]
        atexit.register(lambda: os.path.exists(devrc) or os.unlink(devrc))

        sage.dev.config.Config.__init__(devrc = devrc)

        self['trac'] = {'username': trac_username, 'password_timeout': '.5'}
        self['git'] = {}

        if remote:
            self['git']['repository'] = repository

        self._tmp_dir = tempfile.mkdtemp()
        atexit.register(shutil.rmtree, self._tmp_dir)
        self['git']['dot_git'] = self._tmp_dir

        from sage.dev.git_interface import GitInterface
        GitInterface(self).execute_silent("init")
