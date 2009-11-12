#!/usr/bin/env perl
use File::Copy;

$ver_string = `sage_fortran --version`;
if ($ver_string =~ m/G95/)
{
    print "Using g95";
    copy("patches/setup.py.integrate", "src/scipy/integrate/setup.py");
    copy("patches/setup.py.optimize","src/scipy/optimize/setup.py");
    copy("patches/setup.py.special","src/scipy/special/setup.py");
    copy("patches/setup.py.interpolate","src/scipy/interpolate/setup.py");
    copy("patches/setup.py.odr", "src/scipy/odr/setup.py");
    copy("patches/setup.py.stats", "src/scipy/stats/setup.py");
}

# The following patch is a temporary fix already included upstream.
copy("patches/optimize.py","src/scipy/optimize/optimize.py");
copy("patches/mstats_basic.py","src/scipy/stats/mstats_basic.py");
