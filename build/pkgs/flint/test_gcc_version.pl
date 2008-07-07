#!/usr/bin/env perl
#
# This script tests if gcc is of version greater than
# 3.4.0. If so it exits with status 0 and no output.
# If not, it exits with status 1, and prints a warning
$ver_string=`gcc -dumpversion`;
$ver_string =~ m/(\d)\.(\d)/;



if($1 > 3)
   {
      print "Found gcc 4 or later\n";
      exit(0);

   }

if($1 == 3)
{
    if($2 > 3)
    {
        print "Found gcc 3.4.x\n";
	exit(0);
    }

}

print "WARNING: gcc version less than 3.4.0 \n";
exit(1);
