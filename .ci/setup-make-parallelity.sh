#!/bin/sh

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -ex

# Determine the number of threads that can run simultaneously on this system
# (we might not have nproc available.)
# Note that this value is incorrect for some CI providers (notably CircleCI:
# https://circleci.com/docs/2.0/configuration-reference/#resource_class) which
# provision fewer vCPUs than shown in /proc/cpuinfo. Also, setting this value
# too high can lead to RAM being insufficient, so it's best to set this
# variable manually in your CI configuration.
[[ -z "$NTHREADS" ]] && NTHREADS=`grep -E '^processor' /proc/cpuinfo | wc -l` || true
# Set -j and -l for make (though -l is probably stripped by Sage)
[[ -z "$MAKEOPTS" ]] && MAKEOPTS="-j $NTHREADS -l $((NTHREADS-1)).8" || true
# Not all parts of Sage seem to honor MAKEOPTS, so the current way of telling
# the system which concurrency to use, seems to be setting $MAKE.
[[ -z "$MAKE" ]] && MAKE="make $MAKEOPTS" || true
