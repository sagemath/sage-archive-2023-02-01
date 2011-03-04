#!/usr/bin/env bash

# Set SAGE_ROOT to the location of the sage install.
SAGE_ROOT="....."

CUR="`pwd`"   # save the current directory, so can change back after startup

if [ "x$CUR" = "x" ]; then
    echo "**************************************************************"
    echo " Error: the current working directory does not exist. Please  "
    echo " run sage again from an existing directory.                   "
    echo "**************************************************************"
    exit 1
fi

if [ "$SAGE_ROOT" = "....." ];  then
    SAGE_ROOT=`readlink -n "$0" 2> /dev/null` || \
    SAGE_ROOT=`realpath    "$0" 2> /dev/null` || \
    SAGE_ROOT="$0"

    SAGE_ROOT="${SAGE_ROOT%/*}/"

    if [ ! -f "$SAGE_ROOT/local/bin/sage-sage" ]; then
	echo "**************************************************************************"
	echo "You must compile Sage first using 'make' in the Sage root directory." >&2
	echo "(If you have already compiled Sage, you must set the SAGE_ROOT variable in "
	echo "the file '$0')".
	echo "**************************************************************************"
	exit 1
    fi
fi

# Make root absolute:
cd "$SAGE_ROOT"
SAGE_ROOT=`pwd`
export SAGE_ROOT
export CUR

"$SAGE_ROOT/local/bin/sage-sage" "$@"

# This should kill all children of this process too.
# Uncomment this if you have trouble with orphans.
# Note that you'll get an annoying "Killed" message
# whenever Sage exits.
# kill -9 -$$
