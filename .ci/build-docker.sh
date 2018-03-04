#!/bin/sh

# This script gets called from CI to build several flavours of docker images
# which contain Sage.

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

# We speed up the build process by copying built artifacts from ARTIFACT_BASE
# during docker build. See /docker/Dockerfile for more details.
ARTIFACT_BASE=${ARTIFACT_BASE:-sagemath/sagemath-dev:develop}

# Seed our cache with $ARTIFACT_BASE if it exists
docker pull $ARTIFACT_BASE || true

function docker_build {
    time docker build -f docker/Dockerfile --build-arg "MAKE=${MAKE}" --build-arg ARTIFACT_BASE=$ARTIFACT_BASE $@
}

# We use a multi-stage build /docker/Dockerfile. For the caching to be
# effective, we populate the cache by building the build-time-dependencies and
# the make-all target. (Just building the last target is not enough as
# intermediate targets would be discarded from the cache and therefore the
# caching would fail for our actual builds below.)
docker_build --pull --target build-time-dependencies .
docker_build --pull --tag make-all --target make-all .

# Build the release image without build artifacts.
docker_build --target sagemath --tag "$DOCKER_IMAGE_CLI" .
# Build the developer image with the build artifacts intact.
# Note: It's important to build the dev image last because it might be tagged as ARTIFACT_BASE.
docker_build --target sagemath-dev --tag "$DOCKER_IMAGE_DEV" .
