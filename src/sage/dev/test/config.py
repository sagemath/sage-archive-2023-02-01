import atexit
import os

import sage.dev.config

class Config(sage.dev.config.Config):
    def __init__(self, trac_username = "doctest", remote=None):
        devrc = tempfile.mkstemp()[1]
        atexit.register(lambda: os.path.exists(devrc) or os.unlink(devrc))

        sage.dev.config.Config.__init__(devrc = devrc)

        self['trac'] = {'username': trac_username, 'password_timeout': '.5'}
        if remote:
            self['git'] = {'repo': remote}
