# -*- coding: utf-8 -*-
"""
Commandline handling

Note that argparse is not part of Python 2.6, so we cannot rely on it here.
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

import os
import sys
import logging
log = logging.getLogger()



class CmdlineSubcommands(object):

    def __init__(self, argv=None):
        if argv is None:
            argv = sys.argv
        if len(argv) == 1:
            self.print_help()
            sys.exit(0)
        self.subcommand = argv[1]
        self.extra_args = argv[2:]

    def print_help(self):
        from textwrap import dedent
        print(dedent(self.__doc__).lstrip())
        
    def run(self):
        try:
            method = getattr(self, 'run_{0}'.format(self.subcommand))
        except AttributeError:
            log.error('unknown subcommand: {0}'.format(self.subcommand))
            sys.exit(1)
        try:
            method(*self.extra_args)
        except TypeError as err:
            log.error('invalid arguments to the {0} subcommand: {1}'
                      .format(self.subcommand, self.extra_args))
            sys.exit(1)

            
class SagePkgApplication(CmdlineSubcommands):
    """
    sage-pkg
    --------
    
    The package script is used to manage third-party tarballs.
    """
    
    def run_config(self):
        """
        Print the configuration
        """
        from sage_bootstrap.config import Configuration
        print(Configuration())

    def run_list(self):
        """
        Print a list of all available packages
        """
        pass

    def run_name(self, tarball_filename):
        """
        Find the package name given a tarball filename
        """
        pass

    def run_tarball(self, package_name):
        """
        Find the tarball filename given a package name
        """
        pass
    
