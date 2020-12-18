#!/usr/bin/env bash
# to be run from $SAGE_ROOT, with arguments sage-local-${{ env.PREVIOUS_STAGES }}.tar

if [ -z "$PREFIX" ]; then
    PREFIX=$(pwd)/local
fi

# Show all tar files
ls -l $*

# We specifically use the cygwin tar so that symlinks are saved/restored correctly on Windows.
for a in $*; do
    echo Extracting $a
    (cd / && tar xf -) < $a
    rm -f $a
done

# We set the installation records to the same mtime so that no rebuilds due to dependencies
# among these packages are triggered.
(cd "$PREFIX"/var/lib/sage/installed/ && touch .dummy && touch --reference=.dummy *)

# Show what has been built already.
ls -l "$PREFIX" "$PREFIX"/local/var/lib/sage/installed/
df -h

# Rebase!
src/bin/sage-rebase.sh "$PREFIX"
