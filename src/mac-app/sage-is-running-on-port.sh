#!/bin/bash

# sage-is-running-on-port.sh
# Sage
#
# Created by Ivan Andrus on 16/9/10.
# Copyright 2010 __MyCompanyName__. All rights reserved.

WAIT=0
if [ "x$1" = "x--wait" ]; then
    WAIT=1
    shift
fi

# Set the pid file
PID_FILE=${1-~/.sage/sage_notebook.sagenb/sagenb.pid}

# Check the old pid file in case they are running an older version of Sage
if [ ! -e "$PID_FILE" ] && [ -e ~/.sage/sage_notebook.sagenb/twistd.pid ]; then
    PID_FILE=${1-~/.sage/sage_notebook.sagenb/twistd.pid}
fi

# Wait for a server to start (or at least for a pid file)
while [ "$WAIT" = 1 ] && [ ! -e "$PID_FILE" ]; do
    sleep 1
    # Check the old pid file in case they are running an older version of Sage
    if [ ! -e "$PID_FILE" ] && [ -e ~/.sage/sage_notebook.sagenb/twistd.pid ]; then
        PID_FILE=${1-~/.sage/sage_notebook.sagenb/twistd.pid}
    fi
done

# If the server is running, get the port for it (I hope this is right)
if [ -e "$PID_FILE" ] && kill -0 $(cat "$PID_FILE") 2>/dev/null; then
    grep -Eoe 'port=[0-9]+' ~/.sage/sage_notebook.sagenb/twistedconf.tac | cut -f 2 -d=
    exit 0
fi

# It's not running so return failure
echo 0
exit 1
