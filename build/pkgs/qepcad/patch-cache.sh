#!/bin/sh

# This file attempts to provide some sort of automation for converting
# patch files to the cached files that Sage uses at package install
# time.

# To use this, make sure that all your patches are in the patches
# directory and have a .patch extension.  Then:
# 1. Apply all your patches to the source
# 2. Run "patch-cache update-cache"
# 3. Replace your source with the original source

# In the spkg-install script, run "patch-cache update-files" to copy
# the modified files to the right place.

PATCHES=patches
PATCHES_CACHE=$PATCHES/cached
SRC=src

if [ $1 = update-cache ]; then
    rm $PATCHES_CACHE/*
    counter=1
    for p in $PATCHES/*.patch
    do
	for f in `cat $p | grep "+++" | cut -c5- | cut -f1`
	do
	    cachedfile=`basename $f`.$counter
	    cp $SRC/$f $PATCHES_CACHE/$cachedfile && \
		echo  $f > $PATCHES_CACHE/$cachedfile.filename
	    counter=`expr $counter + 1`
	done
    done

elif [ $1 = update-files ]; then
    for f in $PATCHES_CACHE/*.filename
    do
	cache_name=`basename $f .filename`
	real_name=`cat $f`
	echo Updating $real_name
	cp $PATCHES_CACHE/$cache_name $SRC/$real_name
    done
else
    echo "Error: wrong arguments; \nusage: $0 [update-cache | update-files]"
fi
