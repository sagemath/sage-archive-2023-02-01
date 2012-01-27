#!/usr/bin/env bash
#
# Determine whether the C compiler $CC $CFLAGS supports certain
# compiler flags.  "supported" simply means here that a very basic
# program compiles.
#
# Usage: testcflags.sh [CFLAG1] [CFLAG2]...
#
# This will test all given flags one by one (together with $CFLAGS) and
# return (on stdout) a space-separated list of supported flags.  Note
# that the order is important: if CFLAG1 works, it will be added while
# testing CFLAG2.
#
# The exit code is 0 if all given flags are supported, 1 if some flag
# is not supported.
#
# For example, running
#   $ testcflags.sh --foobar -Wall -march=native
# will return (on recent gcc versions)
#   -Wall -march=native
# with exit status 1 (because --foobar is not supported).
#
# A typical use would be:
#   $ CFLAGS="$CFLAGS `testcflags.sh -ffoo-bar`"
#
#
# AUTHORS:
#
# - Jeroen Demeyer (2012-01-27): initial version (#12367)
#
#*****************************************************************************
#       Copyright (C) 2012 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

if [ -z "$CC" ]; then
    echo >&2 "$0: set \$CC before running this script."
    exit 2
fi

cd /tmp
outfile=testcflags-$$
cfile=$outfile.c

cat >$cfile <<EOF
int main(int argc, char **argv)
{
    return 0;
}
EOF

status=0
NEWCFLAGS=""
for testflag in "$@"; do
    if $CC $CFLAGS $NEWCFLAGS "$testflag" $cfile -o $outfile 2>/dev/null; then
        # Success
        NEWCFLAGS="$NEWCFLAGS $testflag"
    else
        # Failure
        status=1
    fi
done

rm -f $TESTFILE $OUTFILE

echo $NEWCFLAGS
exit $status
