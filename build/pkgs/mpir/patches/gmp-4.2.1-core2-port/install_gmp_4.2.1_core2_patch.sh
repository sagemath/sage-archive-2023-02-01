#!/bin/sh
#
#
# USAGE: install_patch.sh path/to/gmp-4.2.1
#
# Description: This script copies the updated files into the correct
# locations to patch GMP 4.2.1 to compile on Core 2 Mac OS X and Linux
# machines.  The script first tests to verify that it is running on
# a Core 2 processor before performing the patch.  If you want to
# force the patch manually, then just copy the following files:
#
#    amd64call.asm     --> gmp-4.2.1/tests
#    lahf_sahf_test.sh --> gmp-4.2.1/mpn/x86_64
#    add_n.asm         --> gmp-4.2.1/mpn/x86_64
#    sub_n.asm         --> gmp-4.2.1/mpn/x86_64
#    addmul_1.asm      --> gmp-4.2.1/mpn/x86_64
#    submul_1.asm      --> gmp-4.2.1/mpn/x86_64
#    mul_1.asm         --> gmp-4.2.1/mpn/x86_64
#
#
# This file and the rest are copyrighted by Jason Worth Martin
# and released under the Gnu LGPL.
#

if [ $# != 1 ] ; then
    echo
    echo "   USAGE: ./install_gmp_4.2.1_core2_patch.sh path/to/gmp-4.2.1";
    echo
    exit;
fi;

PATH_TO_GMP=$1;
cat > tmp_is_core2_cpu.c <<EOF
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>

/*
 * FUNCTION: uint32_t cpuid(uint32_t *output_array, uint32_t input)
 *
 * Places the input into eax and executes the CPUID instruction.
 * The value of eax is returned, the values of ebx,ecx,and edx
 * are placed on output_array[0], output_array[1], and output_array[2]
 * respectively.
 */
uint32_t cpuid(uint32_t *output_array, uint32_t input)
{
  register uint32_t eax asm ("eax");
  register uint32_t ebx asm ("ebx");
  register uint32_t ecx asm ("ecx");
  register uint32_t edx asm ("edx");

  ebx = 0;
  ecx = 0;
  edx = 0;

  eax = input;
  asm ("cpuid");
  output_array[0] = ebx;
  output_array[1] = ecx;
  output_array[2] = edx;
  return(eax);
}

main()
{
  uint32_t response[3];
  uint32_t input;
  uint32_t output;

  uint8_t steppingID, modelID, familyID, proc_type;
  uint8_t ext_modelID, ext_familyID;

  unsigned int Displayed_Family;
  unsigned int Displayed_Model;

  input = 1;
  output = cpuid(response,input);
  steppingID =   (uint8_t)( 0x0000000f & output);
  modelID =      (uint8_t)((0x000000f0 & output) >> 4);
  familyID =     (uint8_t)((0x00000f00 & output) >> 8);
  proc_type =    (uint8_t)((0x00003000 & output) >> 12);
  ext_modelID =  (uint8_t)((0x000f0000 & output) >> 16);
  ext_familyID = (uint8_t)((0x0ff00000 & output) >> 20);

  if (familyID != 0x0f)
    {
      Displayed_Family = familyID;
    }
  else
    {
      Displayed_Family = ext_familyID + familyID;
    }

  if ((familyID == 0x06) || (familyID == 0x0f) )
    {
      Displayed_Model = ((unsigned int)ext_modelID << 4) + modelID;
    }
  else
    {
      Displayed_Model = modelID;
    }
  if ( (Displayed_Model == 0xf) && (Displayed_Family == 0x6) )
    {
      printf("Yes");
    }
  else
    {
      printf("No");
    }
}
EOF
gcc -m64 tmp_is_core2_cpu.c -o tmp_is_core2_cpu > /dev/null 2>&1
if [ -x ./tmp_is_core2_cpu ]; then
    IS_CORE2_CPU=`./tmp_is_core2_cpu`;
    rm -f tmp_is_core2_cpu.c
    rm -f tmp_is_core2_cpu
else
    IS_CORE2_CPU=No;
    rm -f tmp_is_core2_cpu.c
fi
OS_NAME=`uname -s`
if [ $IS_CORE2_CPU == "Yes" ]; then
    echo "Detected Intel Core 2 CPU"
    if [ -r ${PATH_TO_GMP}/configure ]; then
	echo "Found GMP at ${PATH_TO_GMP}";
	if [ ${OS_NAME} == "Darwin" ] ; then
	    version=`grep "^ VERSION=" ${PATH_TO_GMP}/configure |\
             sed -E 's/.*VERSION=.([1-9.]*)./\1/'`;
	else
	    version=`grep "^ VERSION=" ${PATH_TO_GMP}/configure |\
             sed -r 's/.*VERSION=.([1-9.]*)./\1/'`;
	fi
	if [ $version == "4.2.1" ]; then
	    echo "Version 4.2.1 of GMP found."
	    echo "Copying patch files:"
	    echo "amd64call.asm     --> ${PATH_TO_GMP}/tests"
	    cp -f amd64call.asm ${PATH_TO_GMP}/tests
	    echo "lahf_sahf_test.sh --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f lahf_sahf_test.sh ${PATH_TO_GMP}/mpn/x86_64
	    echo "add_n.asm         --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f add_n.asm ${PATH_TO_GMP}/mpn/x86_64
	    echo "sub_n.asm         --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f sub_n.asm ${PATH_TO_GMP}/mpn/x86_64
	    echo "addmul_1.asm      --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f addmul_1.asm ${PATH_TO_GMP}/mpn/x86_64
	    echo "submul_1.asm      --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f submul_1.asm ${PATH_TO_GMP}/mpn/x86_64
	    echo "mul_1.asm         --> ${PATH_TO_GMP}/mpn/x86_64"
	    cp -f mul_1.asm ${PATH_TO_GMP}/mpn/x86_64
	    echo
	    if [ $OS_NAME == "Darwin" ] ; then
		echo
		echo
		echo "   *********************************************";
		echo "   * If compiling on Mac OS X, remember to set *";
		echo "   * CFLAGS=\"-m64 -fast\" and                   *";
		echo "   * build=x86_64-apple-darwin                 *";
		echo "   * when you configure.                       *"
		echo "   *********************************************";
		echo
	    fi
	else
	    echo "ERROR: Incorrect GMP Version.  Found GMP version: $version"
	    echo "       needed version 4.2.1"
	fi
    else
	echo "ERROR: Cannot read file ${PATH_TO_GMP}/configure"
	echo
	echo "USAGE: install_patch.sh path/to/gmp-4.2.1"
	echo "See the README file for details"
    fi;
else
    echo "Did not detect Intel Core 2 CPU.  Patch not installed"
fi;
