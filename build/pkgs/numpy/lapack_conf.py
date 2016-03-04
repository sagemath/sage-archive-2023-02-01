#!/usr/bin/env python

import pkgconfig, os

conf_file=open('site.cfg', 'w')

conf_file.write('[DEFAULT]\n')
conf_file.write('library_dirs = '+ os.environ['SAGE_LOCAL']+ '/lib\n')
conf_file.write('include_dirs = '+ os.environ['SAGE_LOCAL']+ '/include\n')

if not (os.environ['UNAME'] == 'Darwin'):
    conf_file.write('[blas]\n')
    conf_file.write('include_dirs = '+ ':'.join(list(pkgconfig.parse('cblas')['include_dirs']))+'\n')
    conf_file.write('library_dirs = '+ ':'.join(list(pkgconfig.parse('blas cblas')['library_dirs']))+'\n')
    conf_file.write('blas_libs = '+ ','.join(list(pkgconfig.parse('blas cblas')['libraries']))+'\n')
    conf_file.write('[lapack]\n')
    conf_file.write('library_dirs = '+ ':'.join(list(pkgconfig.parse('lapack')['library_dirs']))+'\n')
    conf_file.write('lapack_libs = '+ ','.join(list(pkgconfig.parse('lapack')['libraries']))+'\n')

conf_file.close()
