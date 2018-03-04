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

[[ -z "$DOCKER_TAG" ]] && DOCKER_TAG=none
[[ "$DOCKER_TAG" = "master" ]] && DOCKER_TAG=latest

DOCKER_IMAGE_CLI=${DOCKER_USER:-sagemath}/sagemath:$DOCKER_TAG
DOCKER_IMAGE_DEV=${DOCKER_USER:-sagemath}/sagemath-dev:$DOCKER_TAG
