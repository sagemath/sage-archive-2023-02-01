#!/bin/sh

export top_level=`pwd`

export saclib=$top_level/src/saclib2.2.0

echo "Patching SACLIB..."

cd $saclib
patch -p1 < $top_level/patches/saclib.patch

export qe=$top_level/src/qesource

echo "Patching QEPCADB..."
cd $qe
patch -p1 < $top_level/patches/cstring.patch
patch -p1 < $top_level/patches/make.patch
patch -p1 < $top_level/patches/qesource.patch

# CAD2D is not built until is fixed to compile under 64-bit
echo "Patching QEPCADB (cad2d)..."
cd $qe
patch -p1 <  $top_level/patches/ncurses-cad2d.patch
