#!/usr/bin/env

import os
import sys

fortran = os.popen2(os.environ['SAGE_LOCAL']+'/bin/'+'which_fortran')[1].read()
SAGE_LOCAL_LIB = os.environ['SAGE_LOCAL']+'/lib'
SAGE_LOCAL_INCLUDE = os.environ['SAGE_LOCAL']+'/include'

if os.environ.has_key('SAGE_ATLAS_LIB'):
    ATLAS_LIB=os.environ['SAGE_ATLAS_LIB']
    if os.path.isdir(ATLAS_LIB):

        # Check for 4 static libraries.
        has_atlas = os.path.exists(ATLAS_LIB+'/lib//libatlas.a')
        has_lapack = os.path.exists(ATLAS_LIB+'/lib/liblapack.a')
        has_cblas = os.path.exists(ATLAS_LIB+'/lib/libcblas.a')
        has_f77blas = os.path.exists(ATLAS_LIB+'/lib/libf77blas.a')
        has_atlas_headers = os.path.exists(ATLAS_LIB+'/include/atlas')

        if has_atlas == True and has_lapack == True and has_cblas == True and has_f77blas == True and has_atlas_headers == True:

            s_gfortran = os.popen2('readelf -s ' +ATLAS_LIB+'/lib/libf77blas.so | grep gfortran')[1].read()
            s_g95 = os.popen2('readelf -s ' + ATLAS_LIB + '/lib/libf77blas.so | grep g95')[1].read()

            if s_gfortran !='' and not fortran.startswith('gfortran'):
                print "Symbols in lib77blas indicate it was build with gfortran \n"
                print "However SAGE is using a different fortran compiler \n"
                print "If you wish to use this blas library, make sure SAGE_FORTRAN points \n"
                print "to a fortran compiler compatible with this library. \n"
                sys.exit(2)

            if s_g95 !='' and not fortran.startswith('g95'):
                print "Symbols in lib77blas indicate it was build with g95 \n"
                print "However SAGE is using a different fortran compiler \n"
                print "If you wish to use this blas library, make sure SAGE_FORTRAN points \n"
                print "to a fortran compiler compatible with this library. \n"
                sys.exit(2)

            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libatlas.a '  + SAGE_LOCAL_LIB+'/libatlas.a')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libcblas.a '  + SAGE_LOCAL_LIB+ '/libcblas.a')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/liblapack.a '  + SAGE_LOCAL_LIB+'/liblapack.a')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libf77blas.a '  + SAGE_LOCAL_LIB+'/libf77blas.a')
            # We might as well include the shared libraries libatlas.so and libcblas.so, as these
            # build relieably on both Solaris and Linux, and might reduce the size of executable if
            # we link to them.
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libatlas.so '  + SAGE_LOCAL_LIB+'/libatlas.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libcblas.so '  + SAGE_LOCAL_LIB+ '/libcblas.so')
            # Buidling of liblapack.so and libf77blas.so is not reliable, but we will link to them
            # People will have to make sure if they do build them, that they work. But these changes
            # only impact SAGE_ATLAS_LIB, which currently is broken on Solaris and often broken on Linux.
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/liblapack.so '  + SAGE_LOCAL_LIB+'/liblapack.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libf77blas.so '  + SAGE_LOCAL_LIB+'/libf77blas.so')
            # Link to the ATLAS related header files.
            os.system(' ln -sf ' + ATLAS_LIB + '/include/atlas '  + SAGE_LOCAL_INCLUDE+'/')
            os.system(' ln -sf ' + ATLAS_LIB + '/include/cblas.h '  + SAGE_LOCAL_INCLUDE+'/cblas.h')
            os.system(' ln -sf ' + ATLAS_LIB + '/include/clapack.h '  + SAGE_LOCAL_INCLUDE+'/clapack.h')
            sys.exit(0)

        else:

            print("Unable to find one of liblapack.a, libcblas.a, libatlas.a or libf77blas.a")
            print("in supplied directory.")
            print("Set SAGE_ATLAS_LIB to the parent directory of liblapack.a, libcblas.a, libatlas.a and libf77blas.a,")
            print("to use existing ATLAS libraries, or unset SAGE_ATLAS_LIB to build ATLAS from source.")
            print("Then type make.")
            sys.exit(2)
    else:
            print("Unable to find one of liblapack.a, libcblas.a, libatlas.a or libf77blas.a")
            print("in supplied directory.")
            print("Set SAGE_ATLAS_LIB to the parent directory of liblapack.a, libcblas.a, libatlas.a and libf77blas.a,")
            print("to use existing ATLAS libraries, or unset SAGE_ATLAS_LIB to build ATLAS from source.")
            print("Then type make.")
            sys.exit(2)

else:
    sys.exit(1)
