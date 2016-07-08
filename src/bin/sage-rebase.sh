#! /bin/bash

# Author:
# * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
#
# Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories),
# but do not touch the ones already stored in the system database,
# and do not update it.
# Note that subsequent calls to 'rebaseall' will not update the Sage dlls.
#
# Invoke this script from a shell after going to the SAGE_LOCAL directory.

echo "Getting list of dlls. This may take a while..."
/bin/find . -name "*.dll" -o -name "*.so" -o -name "*.fas" > /tmp/sage-dlls.lst
echo "Now rebasing..."
/bin/rebase -O -T /tmp/sage-dlls.lst
