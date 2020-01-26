#! /usr/bin/env bash
## Write a Dockerfile to stdout that tests that the packages listed in the debian.txt/fedora.txt files of standard spkg exist
## and satisfy the requirements tested by spkg-configure.m4
## This is called by $SAGE_ROOT/tox.ini
set -e
SYSTEM="${1:-debian}"
shopt -s extglob
TYPE_PATTERN="${2:-standard}"
WITH_SYSTEM_SPKG="${3:-yes}"
#
STRIP_COMMENTS="sed s/#.*//;"
SAGE_ROOT=.
SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SAGE_ROOT/build/pkgs/$SYSTEM{,-bootstrap}.txt))
CONFIGURE_ARGS="--enable-option-checking "
for PKG_SCRIPTS in build/pkgs/*; do
    if [ -d $PKG_SCRIPTS ]; then
        PKG_BASE=$(basename $PKG_SCRIPTS)
        SYSTEM_PACKAGES_FILE=$PKG_SCRIPTS/$SYSTEM.txt
        PKG_TYPE=$(cat $PKG_SCRIPTS/type)
        if [ -f $SYSTEM_PACKAGES_FILE ]; then
           PKG_SYSTEM_PACKAGES=$(echo $(${STRIP_COMMENTS} $SYSTEM_PACKAGES_FILE))
           if [ -n "PKG_SYSTEM_PACKAGES" ]; then
               case "$PKG_TYPE" in
                   $TYPE_PATTERN)
                       SYSTEM_PACKAGES+=" $PKG_SYSTEM_PACKAGES"
                       if [ -f $PKG_SCRIPTS/spkg-configure.m4 ]; then
                           CONFIGURE_ARGS+="--with-system-$PKG_BASE=${WITH_SYSTEM_SPKG} "
                       fi
                       ;;
               esac
           fi
        fi
    fi
done
case $SYSTEM in
    debian*|ubuntu*)
        cat <<EOF
ARG BASE_IMAGE=ubuntu:latest
FROM \${BASE_IMAGE}
RUN apt-get update &&  DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends --yes $SYSTEM_PACKAGES
EOF
        ;;
    fedora*|redhat*|centos*)
        cat <<EOF
ARG BASE_IMAGE=fedora:latest
FROM \${BASE_IMAGE}
RUN yum install -y $SYSTEM_PACKAGES
EOF
        ;;
    arch*)
        # https://hub.docker.com/_/archlinux/
        cat <<EOF
ARG BASE_IMAGE=archlinux:latest
FROM \${BASE_IMAGE}
RUN pacman -Syu --noconfirm $SYSTEM_PACKAGES
EOF
        ;;
    *)
        echo "Not implemented: package installation for SYSTEM=$SYSTEM"
        exit 1
        ;;
esac
cat <<EOF
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
EOF
if [ ${WITH_SYSTEM_SPKG} = "force" ]; then
    cat <<EOF
RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (echo "********** configuring without forcing ***********"; cat config.log; ./configure --enable-build-as-root; cat config.log; exit 1)
EOF
else
    cat <<EOF
RUN echo "****** Configuring: ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} *******"; ./configure --enable-build-as-root $CONFIGURE_ARGS \${EXTRA_CONFIGURE_ARGS} || (cat config.log; exit 1)
EOF
fi
cat <<EOF
# We first compile base-toolchain because otherwise lots of packages are missing their dependency on 'patch'
ARG NUMPROC=8
ENV MAKE="make -j\${NUMPROC}"
RUN make base-toolchain
# Avoid running the lengthy testsuite of the following.
RUN make cython
# Compile something tricky: Everything that uses BLAS.
ARG TARGETS="scipy cbc csdp fflas_ffpack gsl iml numpy r suitesparse cvxopt"
RUN SAGE_CHECK=yes MAKE="make -j4" make -k \${TARGETS}
EOF
