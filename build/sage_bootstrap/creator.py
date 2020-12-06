# -*- coding: utf-8 -*-
"""
Package Creator
"""

#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_ROOT



class PackageCreator(object):

    def __init__(self, package_name):
        self.package_name = package_name
        self.path = os.path.join(SAGE_ROOT, 'build', 'pkgs', package_name)
        try:
            os.mkdir(self.path)
        except OSError:
            pass

    def set_version(self, version):
        """
        Write the version to ``package-version.txt``
        """
        with open(os.path.join(self.path, 'package-version.txt'), 'w+') as f:
            f.write(version)
            f.write('\n')
            
    def set_type(self, pkg_type):
        """
        Write the package type to ``type``
        """
        with open(os.path.join(self.path, 'type'), 'w+') as f:
            f.write(pkg_type)
            f.write('\n')
            
    def set_tarball(self, tarball, upstream_url):
        """
        Write the tarball name pattern to ``checksums.ini``
        """
        with open(os.path.join(self.path, 'checksums.ini'), 'w+') as f:
            f.write('tarball={0}'.format(tarball))
            f.write('\n')
            if upstream_url:
                f.write('upstream_url={0}'.format(upstream_url))
            f.write('\n')
            
    def set_description(self, description, license, upstream_contact):
        """
        Write the ``SPKG.rst`` file
        """
        with open(os.path.join(self.path, 'SPKG.rst'), 'w+') as f:
            def heading(title, char='-'):
                return '{0}\n{1}\n\n'.format(title, char * len(title))
            f.write(heading(self.package_name, '='))
            f.write(heading('Description'))
            if description:
                f.write('{0}\n\n'.format(description))
            f.write(heading('License'))
            if license:
                f.write('{0}\n\n'.format(license))
            f.write(heading('Upstream Contact'))
            if upstream_contact:
                f.write('{0}\n\n'.format(upstream_contact))

    def set_python_data_and_scripts(self):
        """
        Write the files ``dependencies``, ``spkg-install.in``, and ``install-requires.txt``
        """
        with open(os.path.join(self.path, 'dependencies'), 'w+') as f:
            f.write('$(PYTHON) | $(PYTHON_TOOLCHAIN)\n\n')
            f.write('----------\nAll lines of this file are ignored except the first.\n')
        with open(os.path.join(self.path, 'spkg-install.in'), 'w+') as f:
            f.write('cd src\nsdh_pip_install .\n')
        with open(os.path.join(self.path, 'install-requires.txt'), 'w+') as f:
            f.write('{0}\n'.format(self.package_name))
