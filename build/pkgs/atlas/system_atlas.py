#!/usr/bin/env

import os
import sys

fortran = os.popen2(os.environ['SAGE_LOCAL']+'/bin/'+'which_fortran')[1].read()
SAGE_LOCAL_LIB = os.environ['SAGE_LOCAL']+'/lib'
SAGE_LOCAL_INCLUDE = os.environ['SAGE_LOCAL']+'/include'

if os.environ.has_key('SAGE_ATLAS_LIB'):
    ATLAS_LIB=os.environ['SAGE_ATLAS_LIB']
    if os.path.isdir(ATLAS_LIB):

        has_atlas = os.path.exists(ATLAS_LIB+'/lib//libatlas.so')
        if sys.platform.startswith('sunos'):
            has_lapack = os.path.exists(ATLAS_LIB+'/lib/liblapack.a')
        else:
            has_lapack = os.path.exists(ATLAS_LIB+'/lib/liblapack.so')
        has_cblas = os.path.exists(ATLAS_LIB+'/lib/libcblas.so')
        has_f77blas = os.path.exists(ATLAS_LIB+'/lib/libf77blas.so')
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


            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libatlas.so '  + SAGE_LOCAL_LIB+'/libatlas.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libcblas.so '  + SAGE_LOCAL_LIB+ '/libcblas.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/liblapack.so '  + SAGE_LOCAL_LIB+'/liblapack.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/lib/libf77blas.so '  + SAGE_LOCAL_LIB+'/libf77blas.so')
            os.system(' ln -sf ' + ATLAS_LIB + '/include/atlas '  + SAGE_LOCAL_INCLUDE+'/')
            os.system(' ln -sf ' + ATLAS_LIB + '/include/cblas.h '  + SAGE_LOCAL_INCLUDE+'/cblas.h')
            os.system(' ln -sf ' + ATLAS_LIB + '/include/clapack.h '  + SAGE_LOCAL_INCLUDE+'/clapack.h')
            sys.exit(0)

        else:

            print("Unable to find one of liblapack.so, libcblas.so, libatlas.so, or libf77blas.so")
            print("in supplied directory.")
            print("Set SAGE_ATLAS_LIB to the lib directory of liblapack.so,libcblas.so,libatlas.so,libf77blas.so,")
            print("to use existing atlas libraries, or unset SAGE_ATLAS_LIB to build one from source.")
            print("Then type make.")
            sys.exit(2)
    else:
            print("Unable to find one of liblapack.so, libcblas.so, libatlas.so, or libf77blas.so")
            print("in supplied directory.")
            print("Set SAGE_ATLAS_LIB to the lib directory of liblapack.so,libcblas.so,libatlas.so,libf77blas.so,")
            print("to use existing atlas libraries, or unset SAGE_ATLAS_LIB to build one from source.")
            print("Then type make.")
            sys.exit(2)

else:
    sys.exit(1)
