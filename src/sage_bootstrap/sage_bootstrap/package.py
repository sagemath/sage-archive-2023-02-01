# -*- coding: utf-8 -*-
"""
Sage Packages
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import re
import os

import logging
log = logging.getLogger()

from sage_bootstrap.env import SAGE_ROOT



class Package(object):

    def __init__(self, package_name):
        """
        Sage Package

        A package is defined by a subdirectory of
        ``SAGE_ROOT/build/pkgs/``. The name of the package is the name
        of the subdirectory; The metadata of the package is contained
        in various files in the package directory. This class provides
        an abstraction to the metadata, you should never need to
        access the package directory directly.

        INPUT: 

        -- ``package_name`` -- string. Name of the package. The Sage
           convention is that all package names are lower case.
        """
        self.name = package_name
        self._init_checksum()
        self._init_version()

    def __eq__(self, other):
        return self.tarball == other.tarball
        
    @classmethod
    def all(cls):
        """
        Return all packages
        """
        base = os.path.join(SAGE_ROOT, 'build', 'pkgs')
        for subdir in os.listdir(base):
            path = os.path.join(base, subdir) 
            if not os.path.isdir(path):
                continue
            yield cls(subdir)

    @property
    def path(self):
        """
        Return the package directory
        """
        return os.path.join(SAGE_ROOT, 'build', 'pkgs', self.name)
            
    def _init_checksum(self):
        """
        Load the checksums from the appropriate ``checksums.ini`` file
        """
        checksums_ini = os.path.join(self.path, 'checksums.ini')
        assignment = re.compile('(?P<var>[a-zA-Z0-9]*)=(?P<value>.*)')
        result = dict()
        with open(checksums_ini, 'rt') as f:
            for line in f.readlines():
                match = assignment.match(line)
                if match is None:
                    continue
                var, value = match.groups()
                result[var] = value
        self.md5 = result['md5']
        self.sha1 = result['sha1']
        self.cksum = result['cksum']
        self.tarball_pattern = result['tarball']
        
    VERSION_PATCHLEVEL = re.compile('(?P<version>.*)\.p(?P<patchlevel>[0-9]+)')
    
    def _init_version(self):
        with open(os.path.join(self.path, 'package-version.txt')) as f:
            package_version = f.read().strip()
        match = self.VERSION_PATCHLEVEL.match(package_version)
        if match is None:
            self.version = package_version
            self.patchlevel = -1
        else:
            self.version = match.group('version')
            self.patchlevel = match.group('patchlevel')
        self.tarball = self.tarball_pattern.replace('VERSION', self.version)
        
        
