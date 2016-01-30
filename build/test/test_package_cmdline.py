# -*- coding: utf-8 -*-
"""
Test sage-download-file commandline utility
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
import unittest
import subprocess

from sage_bootstrap.env import SAGE_DISTFILES
from sage_bootstrap.download.mirror_list import MIRRORLIST_FILENAME
from sage_bootstrap.util import is_url
from sage_bootstrap.package import Package


EXECUTABLE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)),
    'bin',
    'sage-package',
)


class SagePackageTestCase(unittest.TestCase):

    def run_command(self, *args):
        proc = subprocess.Popen(
            args,
            stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        )
        stdout, stderr = proc.communicate()
        stdout = stdout.decode('utf-8')
        stderr = stderr.decode('utf-8')
        rc = proc.returncode
        return (rc, stdout, stderr)
    
    def test_config(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'config')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Configuration:\n'))
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_list(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'list')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue('configure' in stdout.splitlines())
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_name(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'name', pkg.tarball_filename)
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), 'configure')
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_tarball(self):
        pkg = Package('configure')
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'tarball', 'configure')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertEqual(stdout.rstrip(), pkg.tarball_filename)
        # Prints nothing to stderr
        self.assertEqual(stderr, '')

    def test_apropos(self):
        rc, stdout, stderr = self.run_command(EXECUTABLE, 'apropos', 'python')
        # returns successfully
        self.assertEqual(rc, 0)
        # Prints to stdout
        self.assertTrue(stdout.startswith('Did you mean:'))
        # Prints nothing to stderr
        self.assertEqual(stderr, '')
