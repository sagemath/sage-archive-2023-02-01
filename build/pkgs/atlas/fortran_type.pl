#!/usr/bin/env perl

$output=`sage_fortran --version`;

if ($output =~ m/G95/)
{
print("g95\n");
}
else{
print("gfortran\n");
}
