#! /bin/sh
# Run this script from SAGE_ROOT. Invoke with "--sudo" if sudo is needed.
#
# Install standard development tools - see SAGE_ROOT/build/pkgs/_develop
# This includes the prerequisites for VS Code remote containers,
export PATH=$(pwd)/build/bin:$PATH
SYSTEM=$(sage-guess-package-system)
eval $(sage-print-system-package-command $SYSTEM "$@" update)
eval $(sage-print-system-package-command $SYSTEM --yes "$@" --spkg install _prereq _develop $(head -n 1 build/pkgs/_develop/dependencies) $EXTRA_SAGE_PACKAGES)
