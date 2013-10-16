r"""
Configuration wrapper for doctesting

This module provides a wrapper for ``devrc`` which can be used for doctesting
without tampering with the user's ``devrc`` file.

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

import os
import sage.dev.config


class DoctestConfig(sage.dev.config.Config):
    r"""
    A :class:`sage.dev.config.Config` which lives in a temporary file and sets
    some sensible defaults for doctesting.

    This also initializes an empty git repository in a temporary directory.

    INPUT:

    - ``trac_username`` -- a string (default: ``'doctest'``), a (fake)
      username on trac

    - ``repository`` - a string or ``None`` (default: ``None``), a
      remote repository to push to and pull from

    EXAMPLES::

        sage: from sage.dev.test.config import DoctestConfig
        sage: DoctestConfig()
        Config('''
        [trac]
        username = doctest
        ticket_cache = ...
        [UI]
        log_level = 1
        [git]
        ssh_key_set = True
        repository_anonymous = remote_repository_undefined
        repository = remote_repository_undefined
        src = ...
        dot_git = ...
        [sagedev]
        ''')
    """
    def __init__(self, trac_username="doctest", repository=None):
        r"""
        Initialization.

        TESTS::

            sage: from sage.dev.test.config import DoctestConfig
            sage: config = DoctestConfig()
            sage: type(config)
            <class 'sage.dev.test.config.DoctestConfig'>
        """
        from sage.dev.misc import tmp_dir, tmp_filename
        self._tmp_dir = tmp_dir()
        devrc = os.path.join(self._tmp_dir, 'devrc')
        sage.dev.config.Config.__init__(self, devrc=devrc)

        self['trac'] = {'username': trac_username}

        # Note: ConfigParser allows only string values
        from sage.dev.user_interface import INFO
        self['UI'] = {'log_level': str(INFO)}

        self['git'] = {'ssh_key_set': "True"}
        self['sagedev'] = {}

        self['git']['repository_anonymous'] = \
            self['git']['repository'] = \
            repository if repository else "remote_repository_undefined"

        self['trac']['ticket_cache'] = os.path.join(self._tmp_dir, "ticket_cache")
        repo = os.path.join(self._tmp_dir, 'repo')
        self['git']['src'] = repo
        self['git']['dot_git'] = os.path.join(repo, ".git")
        os.makedirs(self['git']['dot_git'])

        self['git']['user.name'] = trac_username
        self['git']['user.email'] = 'doc@test.test'

        from sage.dev.git_interface import GitInterface
        from sage.dev.test.user_interface import DoctestUserInterface
        old_cwd = os.getcwd()
        os.chdir(self['git']['src'])
        try:
            GitInterface(
                self['git'],
                DoctestUserInterface(self["UI"])
            ).silent.init(self['git']['src'])
        finally:
            os.chdir(old_cwd)
