#!/bin/sh

# This script gets called from CI to run minimal tests on the sagemath-dev image.
# This script expects a single argument, the full name of the docker image to
# test.

# ****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
# ****************************************************************************

set -exo pipefail

IMAGE="$1"

. .ci/setup-make-parallelity.sh

# Usage: timed_run limit args
# Runs $IMAGE with args and check that it terminates with a zero exit code in at most limit seconds.
function timed_run {
    START=`date +%s`
    docker run -e MAKE="$MAKE" "$IMAGE" "$2"
    END=`date +%s`
    TOTAL=$((END-START))
    echo "Checking that \"$2\" was fast…"
    [[ $TOTAL -lt $1 ]]
}

timed_run 60 true # runs make build
# TODO: Can't we get this faster than that?
timed_run $(( 1200/$NTHREADS )) make # runs make build and then make
