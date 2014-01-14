#!/bin/bash

# open-location.sh
# Sage
#
# Created by Ivan Andrus on 16/1/10.
# Copyright 2010 Ivan Andrus. All rights reserved.

if [ "x$1" = "x" ]; then
    echo "usage: $0 URL"
    echo "   open URL in SageNotebook"
    exit
fi

# Find path to this app
BROWSER=`readlink -n "$0" 2> /dev/null` || \
BROWSER=`realpath    "$0" 2> /dev/null` || \
BROWSER="$0"
BROWSER="${SAGE_BROWSER%.app*}.app"

# Tell it to open in this browser
open -a "$BROWSER" "$1"
