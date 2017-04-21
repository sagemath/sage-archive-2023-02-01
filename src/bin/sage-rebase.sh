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
SAGE_LOCAL=${1%/}

if [ -z "$SAGE_LOCAL" ]; then
    # Assume we are in $SAGE_LOCAL by default (the old behavior of this script)
    SAGE_LOCAL=.
fi

FINDFLAGS="-type f ( -name *.dll -o -name *.so -o -name *.fas ) -print"
FINDFLAGS="$FINDFLAGS -o -path "$SAGE_LOCAL"/var/tmp -prune"

echo "Getting list of dlls. This may take a while..."
/bin/find "$SAGE_LOCAL" $FINDFLAGS > /tmp/sage-dlls.lst
echo "Now rebasing..."
/bin/rebase -O -T /tmp/sage-dlls.lst
