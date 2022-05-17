#! /bin/sh
# Run this script from SAGE_ROOT.
#
# Install standard development tools - see SAGE_ROOT/build/pkgs/_develop
# This includes the prerequisites for VS Code remote containers,
export PATH=$(pwd)/build/bin:$PATH
SYSTEM=$(sage-guess-package-system)
eval $(sage-print-system-package-command $SYSTEM --yes install $(sage-get-system-packages $SYSTEM _develop $(head -n 1 build/pkgs/_develop/dependencies)))
