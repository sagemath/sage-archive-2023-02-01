# -*- coding: utf-8 -*-
"""
Set Up Logging

When using a script interactively, logging should go to stdout and be
human-readable. When using a script as part of a pipe (usually
involving tee), logging should go to stderr.
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

import sys
import os
import logging

logger = logging.getLogger()


default_formatter = logging.Formatter(
    '%(levelname)s [%(module)s.%(funcName)s:%(lineno)s]: %(message)s')

plain_formatter = logging.Formatter('%(message)s')


class ExcludeInfoFilter(logging.Filter):
    
    def filter(self, record):
        return record.levelno != logging.INFO
        
class OnlyInfoFilter(logging.Filter):
    
    def filter(self, record):
        return record.levelno == logging.INFO

        
def init_logger(config):
    level = getattr(logging, config.log.upper())
    logger.setLevel(level)
    ch_stderr = logging.StreamHandler(sys.stderr)
    ch_stderr.setLevel(logging.DEBUG)
    ch_stderr.setFormatter(default_formatter)
    if config.interactive:
        ch_stderr.addFilter(ExcludeInfoFilter())
        ch_stdout = logging.StreamHandler(sys.stdout)
        ch_stdout.setLevel(logging.DEBUG)
        ch_stdout.setFormatter(plain_formatter)
        ch_stdout.addFilter(OnlyInfoFilter())
        logger.addHandler(ch_stdout)
    logger.addHandler(ch_stderr)

