#! /bin/dash

# Authors:
# * Gary Zablackis <gzabl@yahoo.com>
# * Dmitrii Pasechnik <dimpase@gmail.com>
# * Jean-Pierre Flori <jean-pierre.flori@ssi.gouv.fr>
#
# Rebase all dlls in the SAGE_LOCAL directory (and its subdirectories)
# as well as the ones already stored in the system database,
# and update the database.
# This system-wide database is located in '/etc/rebase.db.i386' and
# includes the Cygwin system dlls.
#
# Invoke this script from a dash shell after going to the SAGE_LOCAL directory.
# Ensure that no other Cygwin processes are currently running.
# Note that you need write access to the system-wide rebase database
# (which usually means admin rights).

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
/bin/rebaseall -s dll -s exe -s so -s fas -T /tmp/sage-dlls.lst
