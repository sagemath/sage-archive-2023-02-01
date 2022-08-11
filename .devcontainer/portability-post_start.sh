#! /bin/sh
# The portability-.../devcontainer.json configurations run this script after the container is started.
#
# The script assumes that it is run from SAGE_ROOT.
#
# If "config.log" or "logs" are symlinks (for example, created by 'tox -e local-...',
# or after https://trac.sagemath.org/ticket/33262), remove them.
for f in config.log logs; do
    if [ -L $f ]; then
        rm -f $f
    fi
done
# Bootstrap, configure, and build the Sage distribution, reusing the Sage installation from the prebuilt image.
set -e
set -x
make configure
./configure --enable-build-as-root --prefix=/sage/local --with-sage-venv
make build V=0
