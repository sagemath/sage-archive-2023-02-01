#!/usr/bin/env bash

# Author: Gary Zablackis <gzabl@yahoo.com>

# Rebase all dlls in the current sage directory (and subdirectories)
# This is needed when installing SAGE binaries under Windows.

# It is currently not yet deployed, since I don't know when/how/what
# to do with it.  -- William Stein

echo "Getting list of dlls...This may take awhile..."

/bin/find -name *.dll > Sage-dlls.lst

rebaseall -T Sage-dlls.lst
