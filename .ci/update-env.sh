#!/bin/sh

# This script gets called from CI to establish the name of the current docker
# tag to build from the name of the branch/tag provided by CI.

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

# From the docker documentation: "A tag name must be valid ASCII and may
# contain lowercase and uppercase letters, digits, underscores, periods and
# dashes. A tag name may not start with a period or a dash and may contain a
# maximum of 128 characters."
export DOCKER_TAG=`echo $DOCKER_TAG | tr -d '[:space:]' | tr -c '[:alnum:]_.-' '-' | sed 's/^[-.]*//' | cut -c1-128`

[[ -z "$DOCKER_TAG" ]] && export DOCKER_TAG=none
[[ "$DOCKER_TAG" = "master" ]] && export DOCKER_TAG=latest

export DOCKER_IMAGE_CLI=${DOCKER_USER:-sagemath}/sagemath:$DOCKER_TAG
export DOCKER_IMAGE_DEV=${DOCKER_USER:-sagemath}/sagemath-dev:$DOCKER_TAG
