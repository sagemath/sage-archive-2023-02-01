#! /usr/bin/env bash
## Write a Dockerfile to stdout that tests that the packages listed in the debian.txt files of standard spkg exist
## and satisfy the requirements tested by spkg-configure.m4
set -e
STRIP_COMMENTS="sed s/#.*//;"
SAGE_ROOT=.
SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SAGE_ROOT/build/pkgs/debian.txt))
# needed for bootstrap:
SYSTEM_PACKAGES+=" gettext autoconf automake libtool"
CONFIGURE_ARGS="--enable-option-checking "
for PKG_SCRIPTS in build/pkgs/*; do
    if [ -d $PKG_SCRIPTS ]; then
        PKG_BASE=$(basename $PKG_SCRIPTS)
        SYSTEM_PACKAGES_FILE=$PKG_SCRIPTS/debian.txt
        PKG_TYPE=$(cat $PKG_SCRIPTS/type)
        if [ -f $SYSTEM_PACKAGES_FILE -a "$PKG_TYPE" = "standard" ]; then
            SYSTEM_PACKAGES+=" "$(echo $(${STRIP_COMMENTS} $SYSTEM_PACKAGES_FILE))
            if [ -f $PKG_SCRIPTS/spkg-configure.m4 ]; then
                CONFIGURE_ARGS+="--with-system-$PKG_BASE=force "
            fi
        fi
    fi
done
cat <<EOF
ARG BASE_IMAGE=ubuntu:latest
FROM \${BASE_IMAGE}
RUN apt-get update &&  DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --yes $SYSTEM_PACKAGES
# Bootstrapping
RUN mkdir -p /sage
WORKDIR /sage
ADD Makefile VERSION.txt README.md bootstrap configure.ac sage ./
ADD m4 ./m4
ADD build ./build
ADD src/bin/sage-version.sh src/bin/sage-version.sh
RUN ./bootstrap
# Configuring
ADD src/ext src/ext
ADD src/bin src/bin
ADD src/Makefile.in src/Makefile.in
ARG EXTRA_CONFIGURE_ARGS
RUN ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (echo "********** configuring without forcing ***********"; cat config.log; ./configure --enable-build-as-root; cat config.log; exit 1)
RUN make -j4 toolchain V=0
# Compile something tricky
RUN make -j4 scipy
EOF
