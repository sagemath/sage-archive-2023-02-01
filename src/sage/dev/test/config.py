r"""
Configuration wrapper for doctesting

This module provides a wrapper for ``devrc`` which can be used for doctesting
without tampering with the user's ``devrc`` file.

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
from sage.dev.config import Config

class DoctestConfig(Config):
    r"""
    A :class:`sage.dev.config.Config` which lives in a temporary file and sets
    some sensible defaults for doctesting.

    This also initializes an empty git repository in a temporary directory.

    INPUT:

    - ``trac_username`` -- a string (default: ``'doctest'``), a (fake) username
    on trac

    - ``repository`` - a string or ``None`` (default: ``None``), a remote
    repository to push to and pull from

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: DoctestConfig()
        Config('''
        [trac]
        username = doctest
        [UI]
        log_level = 0
        [git]
        dot_git = ...
        ''')

    """
    def __init__(self, trac_username = "doctest", repository=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: type(DoctestConfig())
            <class 'sage.dev.test.config.DoctestConfig'>

        """
        import tempfile, atexit, shutil, os
        devrc = tempfile.mkstemp()[1]
        atexit.register(lambda: os.path.exists(devrc) or os.unlink(devrc))

        Config.__init__(self, devrc = devrc)

        self['trac'] = {'username': trac_username}
        self['UI'] = {'log_level': 0}
        self['git'] = {}

        if repository:
            self['git']['repository'] = repository

        self._tmp_dir = tempfile.mkdtemp()
        atexit.register(shutil.rmtree, self._tmp_dir)
        self['git']['dot_git'] = self._tmp_dir

        from sage.dev.git_interface import GitInterface, SILENT
        from sage.dev.test.user_interface import DoctestUserInterface
        GitInterface(self['git'], DoctestUserInterface(self["UI"])).init(SILENT, self._tmp_dir)
