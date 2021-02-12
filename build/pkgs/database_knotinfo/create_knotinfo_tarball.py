#!/usr/bin/python

r"""
Python script to create a tarball for the sage-package ``database_knotinfo``
in the given path. This utility should be used in case of a switch to a
new version of the data files (that is if the original files on the
KnotInfo LinkInfo web-page have changed). In that case an invocation of
``sage -package update database_knotinfo <new version>`` and
``sage -package fix-checksum database_knotinfo`` will be necessary.

..NOTE::

    This function demands the Python package ``pandas``, ``xlrd`` and
    ``xlsx2csv`` to be installed. If not you have to run::

        pip install pandas
        pip install xlrd
        pip install xlsx2csv

    before using this function.

INPUT:

- ``version`` -- string,  name of the new version to be created
  (by default date of the day of creation)
- ``path`` -- string of the path where the tarball should be stored
  (by default ``pwd``)

EXAMPLES::

    ~/sage $ build/pkgs/database_knotinfo/create_knotinfo_tarball.py 20210201 upstream
    src/
    src/knotinfo_data_complete.csv
    src/linkinfo_data_complete.csv
"""

import sys, os
from xlsx2csv import Xlsx2csv
from pandas import read_excel

##############################################################################
#       Copyright (C) 2021 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################


cmdline_args = sys.argv[1:]

version = None
path    = None

if len(cmdline_args) > 1:
    path    = cmdline_args[1]

if len(cmdline_args) > 0:
    version = cmdline_args[0]


if not version:
    from datetime import datetime
    version = str(datetime.today().date()).replace('-','')

if not path:
    path = os.environ['PWD']

path_temp = os.path.join(path, 'special_knotinfo_spkg_temp_dir')
path_src = os.path.join(path_temp, 'src')
os.makedirs(path_temp)
os.makedirs(path_src)

def convert(path_src, url, filen, reader):
    if reader == Xlsx2csv:
        excel = filen + '.xlsx'
    else:
        excel = filen + '.xls'
    csv   = filen + '.csv'
    inp   = os.path.join(url,      excel)
    out   = os.path.join(path_src, csv)
    if reader == Xlsx2csv:
        from six.moves.urllib.request import urlopen
        f = urlopen(inp)
        url_data = f.read()
        temp_file = os.path.join(path_temp, 'temp.xlsx')
        f = open(temp_file, 'wt')
        f.write(url_data)
        f.close()
        data = reader(temp_file, delimiter='|', skip_empty_lines=True)
        data.convert(out)
    else:
        data = reader(inp)
        data.to_csv(out, sep='|', index=False)

# first KnotInfo (using pandas and xlrd)
convert(path_src, 'https://knotinfo.math.indiana.edu/', 'knotinfo_data_complete', read_excel)

# now LinkInfo (using xlsx2csv)
convert(path_src, 'https://linkinfo.sitehost.iu.edu/', 'linkinfo_data_complete', Xlsx2csv)

tar_file = 'knotinfo-%s.tar.bz2' %version
path_tar = os.path.join(path_temp, tar_file)

os.system('cd %s; tar -cvjSf %s src' %(path_temp, tar_file))
os.system('mv %s %s; rm -rf %s' %(path_tar, path, path_temp))
