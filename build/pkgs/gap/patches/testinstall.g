#############################################################################
##
#W  testinstall.g               GAP library                      Frank Celler
##
##
#Y  Copyright (C)  1997,  Lehrstuhl D für Mathematik,  RWTH Aachen,  Germany
##
##  This file lists those files in the directory <F>tst</F> of the &GAP;
##  distribution that are recommended to be read after a &GAP; installation.
##
##  Each entry in the argument list of <C>RunStandardTests</C> is a pair that
##  consists of the filename (relative to the <F>tst</F> directory) and the
##  scaling factor that occurs in the <C>STOP_TEST</C> call at the end of the
##  test file.
##  <P/>
##  The documentation states the following:
##  <P/>
##  <#GAPDoc Label="[1]{testinstall.g}">
##  If you want to run a quick test of your &GAP; installation
##  (though this is not required), you can read in a test script
##  that exercises some &GAP;'s capabilities.
##  <P/>
##  <Log><![CDATA[
##  gap> Read( Filename( DirectoriesLibrary( "tst" ), "testinstall.g" ) );
##  ]]></Log>
##  <P/>
##  The test requires about 512MB of memory and runs about one
##  minute on an Intel Core 2 Duo / 2.53 GHz machine.
##  You will get a large number of lines with output about the progress
##  of the tests.
##  <#/GAPDoc>
##

Print( "You should start GAP4 using `gap -N -A -x 80 -r -m 100m -o 512m'.\n",
       "The more GAP4stones you get, the faster your system is.\n",
       "The runtime of the following tests (in general) increases.\n",
       "You should expect the test to take about one minute and show about\n",
       "100000 GAP4stones on an Intel Core 2 Duo / 2.53 GHz machine.\n",
       "The `next' time is an approximation of the running time ",
       "for the next file.\n\n" );

Reread( Filename( DirectoriesLibrary( "tst" ), "testutil.g" ) );

RunStandardTests( [
  [ "alghom.tst", 6000000 ],
  [ "algmat.tst", 180800000 ],
  [ "algsc.tst", 59000000 ],
  [ "combinat.tst", 5300000 ],
  [ "ctblfuns.tst", 3900000 ],
  [ "ctblmoli.tst", 98500000 ],
  [ "ctblmono.tst", 33400000 ],
  [ "cyclotom.tst", 900000 ],
  [ "ffe.tst", 3600000 ],
  [ "ffeconway.tst", 50200000 ],
  [ "gaussian.tst", 300000 ],
  [ "grpfp.tst", 146700000 ],
  [ "grpfree.tst", 700000 ],
  [ "grpmat.tst", 481000000 ],
  [ "grppcnrm.tst", 2333400000 ],
  [ "listgen.tst", 1000000 ],
  [ "mapping.tst", 37300000 ],
  [ "mgmring.tst", 1800000 ],
  [ "modfree.tst", 5800000 ],
  [ "onecohom.tst", 50600000 ],
  [ "oprt.tst", 2000000 ],
  [ "ratfun.tst", 800000 ],
  [ "relation.tst", 7700000 ],
  [ "rwspcgrp.tst", 59400000 ],
  [ "semicong.tst", 7800000 ],
  [ "semigrp.tst", 11200000 ],
  [ "semirel.tst", 10900000 ],
  [ "vspchom.tst", 10500000 ],
  [ "vspcmat.tst", 8400000 ],
  [ "vspcrow.tst", 47400000 ],
  [ "xgap.tst", 1206600000 ],
  [ "zlattice.tst", 100000 ],
  [ "zmodnz.tst", 2300000 ],
] );


#############################################################################
##
#E

