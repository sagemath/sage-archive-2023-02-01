#!/usr/bin/env bash

# abort on error --- the spkg-install will catch it
set -e

VER_STRING="`sage_fortran --version`"
case "$VER_STRING" in
    G95*)
        echo "Using g95"
        cp "patches/setup.py.integrate"  "src/scipy/integrate/setup.py"
        cp "patches/setup.py.optimize" "src/scipy/optimize/setup.py"
        cp "patches/setup.py.special" "src/scipy/special/setup.py"
        cp "patches/setup.py.interpolate" "src/scipy/interpolate/setup.py"
        cp "patches/setup.py.odr"  "src/scipy/odr/setup.py"
        cp "patches/setup.py.stats"  "src/scipy/stats/setup.py"
        ;;
esac

# The following patch (optimize.py) is a temporary fix already included upstream.
cp "patches/optimize.py" "src/scipy/optimize/optimize.py"

# Fix an incorrect assert
cp "patches/mstats_basic.py" "src/scipy/stats/mstats_basic.py"
