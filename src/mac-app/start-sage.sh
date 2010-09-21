#!/bin/bash

# start-sage.sh
# Fluidium
#
# Created by Ivan Andrus on 16/1/10.
# Copyright 2010 Ivan Andrus. All rights reserved.

# Ensure we have enough arguments
if [ $# -lt 2 ]; then
    echo "usage: $0 SAGE_EXECUTABLE LOG"
    exit 1;
fi

# Ensure that we have the original sage and therefore the directory it's in should be SAGE_ROOT (right?)
SAGE_EXECUTABLE=`readlink -n "$1" 2> /dev/null` || \
SAGE_EXECUTABLE="$1"

# Strip the last section off -- this should be SAGE_ROOT.  We could
# just call `sage --root`, but I'm afraid this may not work (if there
# are spaces).  But if sageBinary is not set, the "sage" gets passed
# in, and we have no choice, but to call `sage --root`, and it should
# work in that case (assuming of course that sage is in PATH)
if [ "x$SAGE_EXECUTABLE" = "xsage" ]; then
    SAGE_ROOT=`sage --root`
else
    SAGE_ROOT="${SAGE_EXECUTABLE%/*}/"
fi
SAGE_LOG="$2"

# Work around spaces in the path (and perhaps cut down on location changes).
# We delete and recreate the symlink every time to ensure we are in the right place
rm -f /tmp/sage-mac-app
ln -s "$SAGE_ROOT" /tmp/sage-mac-app

# Move to a fake SAGE_ROOT -- without spaces
cd /tmp/sage-mac-app || exit 1

# Set (i.e. "source") SAGE_ROOT and all the other environment settings
echo Setting environment variables >> "$SAGE_LOG"
. ./local/bin/sage-env >> "$SAGE_LOG" 2>> "$SAGE_LOG"
export SAGE_ROOT

# Mac OS X app bundles are *intended* to be moved around, and/or given away
# So always run first the respective script handling this
# (This should also catch Intel vs. PPC or 32Bit vs. 64Bit conflicts - untested)
echo Checking install location >> "$SAGE_LOG"
./local/bin/sage-location >> "$SAGE_LOG" 2>> "$SAGE_LOG"

echo Checking existance of notebook directory >> "$SAGE_LOG"
if [ -d $DOT_SAGE/sage_notebook.sagenb ]; then
    echo Starting Notebook >> "$SAGE_LOG"
    ./sage --notebook >> "$SAGE_LOG" 2>> "$SAGE_LOG"
else
    echo Starting Notebook in Terminal >> "$SAGE_LOG"
    sage-native-execute osascript \
        -e 'tell app "Terminal"' \
        -e '    activate' \
        -e "    do script \"'$SAGE_ROOT'/sage --notebook\"" \
        -e 'end'
fi
exit 0
