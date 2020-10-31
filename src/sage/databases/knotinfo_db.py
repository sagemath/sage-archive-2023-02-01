# -*- coding: utf-8 -*-
r"""
KontInfo Database

This module contains the class :class:`KnotInfoDataBase`  and auxilary classes for it
which serves as an interface to the lists of named knots and links provided at
https://knotinfo.math.indiana.edu/


AUTHORS:

- Sebastian Oehms August 2020: initial version
"""


##############################################################################
#       Copyright (C) 2020 Sebastian Oehms <seb.oehms@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
##############################################################################


import os
import csv
from enum import Enum

from sage.structure.sage_object import SageObject
from sage.misc.persist import save, load
from sage.misc.verbose import verbose
from sage.misc.cachefunc import cached_method
from sage.env import SAGE_SHARE, SAGE_ROOT


class KnotInfoColumnTypes(Enum):
    r"""
    Enum class to specify if a column from the table of knots and links provided by http://www.indiana.edu/~knotinfo
    is used for knots only, links only or both.

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoColumnTypes
        sage: [col_type for col_type in KnotInfoColumnTypes]
        [<KnotInfoColumnTypes.OnlyKnots: 'K'>,
        <KnotInfoColumnTypes.OnlyLinks: 'L'>,
        <KnotInfoColumnTypes.KnotsAndLinks: 'B'>]
    """

    OnlyKnots =     'K'       # column that is only used in the KnotInfo table
    OnlyLinks =     'L'       # column that is only used in the LinkInfo table
    KnotsAndLinks = 'B'       # column that is only used in both tables


class KnotInfoColumns(Enum):
    r"""
    Enum class to select a column from the table of knots and links provided by http://www.indiana.edu/~knotinfo

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase, KnotInfoColumns, KnotInfoColumnTypes
        sage: ki_db = KnotInfoDataBase()
        sage: KnotInfoColumns('Columns', ki_db.read_column_dict())
        <enum 'Columns'>
        sage: [col.column_name() for col in _ if col.column_type() == KnotInfoColumnTypes.OnlyLinks]  # optional - database_knotinfo
        ['Name - Unoriented',
         'Orientation',
         'Unoriented Rank',
         'PD Notation (vector)',
         'PD Notation (KnotTheory)',
         'Multivariable Alexander Polynomial',
         'HOMFLYPT Polynomial',
         'Unoriented',
         'Arc Notation',
         'Linking Matrix',
         'Rolfsen Name',
         'Components',
         'DT code',
         'Splitting Number',
         'Nullity',
         'Unlinking Number',
         'Weak Splitting Number']
    """
    def column_name(self):
        r"""
        Return the name of ``self`` displayed on the KnotInfo web-page.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase, KnotInfoColumns
            sage: ki_db = KnotInfoDataBase()
            sage: KIcols = KnotInfoColumns('Columns', ki_db.read_column_dict())
            sage: KIcols.dt_code.column_name()
            'DT code'
        """
        return self.value[0]

    def column_type(self):
        r"""
        Return the type of ``self``. That is an instance of ``Enum``
        :class:`KnotInfoColumnTypes`.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase, KnotInfoColumns, KnotInfoColumnTypes
            sage: ki_db = KnotInfoDataBase()
            sage: KIcols = KnotInfoColumns('Columns', ki_db.read_column_dict())
            sage: KIcols.homfly_polynomial.column_type()
            <KnotInfoColumnTypes.OnlyKnots: 'K'>
            sage: KIcols.homflypt_polynomial.column_type()
            <KnotInfoColumnTypes.OnlyLinks: 'L'>
            sage: KIcols.name.column_type()
            <KnotInfoColumnTypes.KnotsAndLinks: 'B'>
        """
        return self.value[1]


class KnotInfoFilename(Enum):
    r"""
    Enum for the different data files. The following choices are possible:

    - ``knots`` -- contains the the data from KnotInfo
    - ``links`` -- contains the the data for proper links from LinkInfo

    Examples::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase
        sage: ki_db = KnotInfoDataBase()
        sage: ki_db.filename
        <enum 'KnotInfoFilename'>
    """

    def url(self):
        """
        Return the URL to download the data from the web-page.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.url()
            'https://knotinfo.math.indiana.edu/'
        """
        if self == KnotInfoFilename.knots:
            return self.value[0]
        else:
            return self.value[0]

    def excel(self):
        """
        Return the Excel-file name to download the data from the web-page.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.excel()
            'knotinfo_data_complete.xls'
        """
        if self == KnotInfoFilename.knots:
            return '%s.xls' %(self.value[1])
        else:
            return '%s.xlsx' %(self.value[1])

    def csv(self):
        """
        Return the file name under which the data from the web-page
        are stored as csv file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.csv()
            'knotinfo_data_complete.csv'
        """
        return '%s.csv' %(self.value[1])

    def sobj_num_knots(self):
        """
        Return the file name under which the number of knots
        is stored as in python int in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_num_knots()
            'num_knots.sobj'
        """
        return 'num_knots.sobj'

    def sobj_row(self):
        """
        Return the file name under which the row-data of the csv-File
        is stored as python dictionary in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_row()
            'row_dict.sobj'
        """
        return 'row_dict.sobj'

    def sobj_column(self):
        """
        Return the file name under which the column-data of the csv-File
        is stored as python dictionary in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_column()
            'column_dict.sobj'
        """
        return 'column_dict.sobj'


    def sobj_data(self, column):
        """
        Return the file name under which the data of the given
        column is stored as python list in a sobj-file.

        Examples::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.filename.knots.sobj_data(ki_db.columns().braid_notation)
            'knotinfo_braid_notation'
        """
        if column.column_type() == KnotInfoColumnTypes.OnlyLinks:
            return 'linkinfo_%s' %(column.name)
        else:
            return 'knotinfo_%s' %(column.name)

    knots = ['https://knotinfo.math.indiana.edu/', 'knotinfo_data_complete']
    links = ['https://linkinfo.sitehost.iu.edu',   'linkinfo_data_complete']




#----------------------------------------------------------------------------------------------------------------------------
# Class to provide data for knots and links from the KnotInfo web-page
#----------------------------------------------------------------------------------------------------------------------------
class KnotInfoDataBase(SageObject):
    r"""
    Database interface to KnotInfo

    The original data are obtained from KnotInfo web-page (URL see the example below). In order
    to have these data installed during the build process as a sage-package they are converted
    as csv files into a tarball. This tarball has been created using the method :meth:`create_spkg_tarball`.

    EXAMPLES::

        sage: from sage.databases.knotinfo_db import KnotInfoDataBase
        sage: ki_db = KnotInfoDataBase()
        sage: ki_db.filename.knots
        <KnotInfoFilename.knots: ['https://knotinfo.math.indiana.edu/', 'knotinfo_data_complete']>
    """

    filename = KnotInfoFilename

    def __init__(self):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: from sage.env import SAGE_SHARE
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db._import_path
            '/home/sebastian/develop/sage/local/share/knotinfo'
        """
        self._package = 'knotinfo'
        version_file  = os.path.join(SAGE_ROOT, 'build/pkgs/%s/package-version.txt' %self._package)
        f = open(version_file)
        self._version = f.read().splitlines()[0]
        f.close()

        # some constants
        self._delimiter    = '|'
        self._names_column = 'name'
        self._components_column = 'components'
        self._knot_prefix  = 'K'
        self._import_path  = os.path.join(SAGE_SHARE, self._package)

        self._knot_list = None
        self._link_list = None
        self._available = None
        self._num_knots = None


    def _create_csv_file(self, filename, path_for_src):
        r"""
        Return the data fetched from the web-page as a csv file
        such that it can be parsed via pythons ``csv`` class.

        INPUT:

        - ``filename`` - instance of :class:`KnotInfoDataBase.filename`
        - ``path_for_src`` - string giving the pathname where to store
          the ``csv`` -files

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: import os
            sage: pwd = os.environ['PWD']
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db._create_csv_file(ki_db.filename.knots, pwd)
        """
        # TODO import directly from the internet page and convert to csv via pandoc
        return
        if not isinstance(filename, KnotInfoDataBase.filename):
            raise TypeError('File name must be an instance of enum %s' (KnotInfoDataBase.filename))

        import_file = '%s/%s' %(self._import_path, filename.csv())

        from six.moves.urllib.request import urlopen
        try:
            from urllib.error import HTTPError
        except ImportError:
            from urllib2 import HTTPError

        try:
            url = '%s/%s' %(filename.url(), filename.excel())
            url_data = urlopen(url).read().decode()
        except:
            pass

    def is_available(self):
        r"""
        Return wether the KnotInfo databases are installed or not.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.is_available()       # optional - database_knotinfo
            True
        """
        if not self._available:
            try:
                lib_path = self._import_path
                filename = self.filename.knots.sobj_num_knots()
                self._num_knots =  load('%s/%s' %(lib_path, filename))
                self._available = True
            except FileNotFoundError:
                self._available = False
                self._num_knots = len([v for v in row_demo_sample.values() if v[1]==1])
        return self._available


    def create_spkg_tarball(self, path_for_src=None):
        r"""
        Create a tarball for the sage-package ``knotinfo`` in the ``upstream`` directory. This
        utility should only be used by users who know what they do in case of a switch to a new
        version of the data files (that is if the original files on KnotInfo web-page have changed).
        In that case in invocation of ``sage -package update knotinfo <new version>`` and
        ``sage -package fix-checksum knotinfo`` will be necessary.

        INPUT:

        -- ``path_for_src`` - string of the path under which the source are stored in a
           subdirectory called ``src``

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.create_spkg_tarball()    # not tested (because of internet access)
        """
        if not path_for_src:
            path_for_src = os.environ['PWD']

        for filename in KnotInfoDataBase.filename:
            self._create_csv_file(filename, path_for_src)

        os.system('cd %s; tar -cvjSf %s/upstream/%s-%s.tar.bz2 src' %(path, SAGE_ROOT, self._package, self._version) )


    def version(self):
        r"""
        Return the current version.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.version()
            '20200713'
        """
        return self._version

    def knot_list(self):
        r"""
        Return the KnotInfo table rows as Python list.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.knot_list())  # not tested because its only used on installation
        """
        if self._knot_list:
            return self._knot_list

        print('Importing KnotInfo database from SPKG!')
        os.system('pwd')
        knot_csv = open('src/%s' %self.filename.knots.csv())
        knot_dict = csv.DictReader(knot_csv, delimiter=self._delimiter)
        self._knot_list = list(knot_dict)
        knot_csv.close()
        return self._knot_list


    def link_list(self):
        r"""
        Return the LinkInfo table rows as Python list.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.link_list())  # not tested because its only used on installation
        """
        if self._link_list:
            return self._link_list

        print('Importing LinkInfo database from SPKG!')
        link_csv = open('src/%s' %self.filename.links.csv())
        link_dict = csv.DictReader(link_csv, delimiter=self._delimiter)
        self._link_list = list(link_dict)
        link_csv.close()
        return self._link_list

    def create_col_dict_sobj(self):
        r"""
        Create ``sobj`` files containing the number of knots and a dictionary
        for the columns of the table.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.create_col_dict_sobj() # not tested because its only used on installation
        """
        knot_list = self.knot_list()
        knot_column_names = knot_list[0]
        len_knots = len(knot_list)

        link_list = self.link_list()
        link_column_names = link_list[0]

        from sage.misc.misc import sage_makedirs
        sage_makedirs(self._import_path)

        num_knots = len_knots - 1 
        save(num_knots, '%s/%s' %(self._import_path, self.filename.knots.sobj_num_knots()))

        column_dict = {}

        # ----------------------------------------------------------------
        # Columns that exist for knots and links
        # ----------------------------------------------------------------
        for col in knot_column_names:

            name = knot_column_names[col]
            if not name:
                # not of interest
                continue

            col_type = KnotInfoColumnTypes.OnlyKnots
            if col in link_column_names:
                col_type = KnotInfoColumnTypes.KnotsAndLinks
            column_dict[col] = [name, col_type]

        # ----------------------------------------------------------------
        # Columns that exist for links only
        # ----------------------------------------------------------------
        for col in link_column_names:

            name = link_column_names[col]
            if not name:
                # not of interest
                continue

            if col in knot_column_names:
                # already used
                continue

            col_type = KnotInfoColumnTypes.OnlyLinks
            column_dict[col] = [name, col_type]

        save(column_dict, '%s/%s' %(self._import_path, self.filename.knots.sobj_column()))



    def create_data_sobj(self):
        r"""
        Create ``sobj`` files containing the contents of the whole table.
        To each column there is created one file containing a list of
        strings giving the entries of the database table.

        The length of these lists depends on the type of the corresponding
        column. If a column is used in both tables (``KnotInfoColumnTypes.KnotsAndLinks``)
        the list of proper links is appended to the list of knots.
        In both other cases the lenght of the list corresponds to
        the number of listed knots and proper links respectively.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.create_data_sobj() # not tested because its only used on installation
        """
        knot_list = self.knot_list()
        link_list = self.link_list()
        len_knots = len(knot_list)
        len_links = len(link_list)

        row_dict = {}

        # ----------------------------------------------------------------
        # Columns that exist for knots and links
        # ----------------------------------------------------------------
        for col in self.columns():
            val_list = []

            if  col.column_type() != KnotInfoColumnTypes.OnlyLinks:
                for i in range(1 , len_knots):
                    if col.name == self._names_column:
                        row_dict[self._knot_prefix + knot_list[i][col.name]] = [i - 1 , 1]
                    else:
                        val_list.append(knot_list[i][col.name])

            if  col.column_type() != KnotInfoColumnTypes.OnlyKnots:
                for i in range(1 , len_links):
                    if col.name == self._names_column:
                        link_name = link_list[i][col.name]
                        link_name = link_name.replace('{', '_')
                        link_name = link_name.replace(',', '_')
                        link_name = link_name.replace('}', '')
 
                        num_comp = int(link_list[i][self._components_column])
                        row_dict[link_name] = [i + len_knots - 2 , num_comp]

                    else:
                        val_list.append(link_list[i][col.name])

            if val_list:
                save(val_list, '%s/%s' %(self._import_path, self.filename.knots.sobj_data(col)))

        save(row_dict,    '%s/%s' %(self._import_path, self.filename.knots.sobj_row()))


    @cached_method
    def columns(self):
        r"""
        Return the columns ot the databese table as instances of enum class
        :class:`KnotInfoColumns`.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: cols = ki_db.columns()
            sage: [col.column_name() for col in cols if col.name.startswith('pd')]   # optional - database_knotinfo
            ['PD Notation', 'PD Notation (vector)', 'PD Notation (KnotTheory)']
        """
        column_dict = self.read_column_dict()
        return KnotInfoColumns('Columns', column_dict)


    # -------------------------------------------------------------------------------------------------------------
    # read the dictionary for the column names from sobj-file
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read_column_dict(self):
        r"""
        Read the dictionary for the column names from the according sobj-file

        OUTPUT:

        A python dictionary containing the column names and types

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.read_column_dict())       # optional - database_knotinfo
            120
        """
        if not self.is_available():
            return column_demo_sample
        lib_path = self._import_path
        filename = self.filename.knots.sobj_column()
        return load('%s/%s' %(lib_path, filename))

    # -------------------------------------------------------------------------------------------------------------
    # read the dictionary for the row names that is the knot and link names from sobj-file
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read_row_dict(self):
        r"""
        Read the dictionary for the row names that is the knot and link names
        from the according sobj-file

        OUTPUT:

        A python dictionary containing the names of the knots and links
        together with their table index and the corresponding number of
        components

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: len(ki_db.read_row_dict())          # optional - database_knotinfo
            7166
        """
        if not self.is_available():
            return row_demo_sample
        lib_path = self._import_path
        filename = self.filename.knots.sobj_row()
        return load('%s/%s' %(lib_path, filename))

    # -------------------------------------------------------------------------------------------------------------
    # read the dictionary for the row names that is the knot and link names from sobj-file
    # -------------------------------------------------------------------------------------------------------------
    def read_num_knots(self):
        r"""
        Read the number of knots contained in the database (without
        proper links) from the according sobj-file.

        OUTPUT:

        Integer

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
            sage: ki_db.read_num_knots()              # optional - database_knotinfo
            2978
        """
        if not self._num_knots:
            self.is_available()
        return self._num_knots


    # -------------------------------------------------------------------------------------------------------------
    # read an sobj-file obtained from KnotInfo database
    # -------------------------------------------------------------------------------------------------------------
    @cached_method
    def read(self, column):
        r"""
        Access a column of KnotInfo / LinkInfo

        INPUT:

        ``column`` -- instance of enum :class:`KnotInfoColumns`
          to select the data to be read in

        OUTPUT:

        A python list containing the data corresponding to the column.

        EXAMPLES::

            sage: from sage.databases.knotinfo_db import KnotInfoDataBase
            sage: ki_db = KnotInfoDataBase()
        """
        if not isinstance(column, KnotInfoColumns):
            raise TypeError('column must be an instance of enum %s' (KnotInfoColumns))

        if not self.is_available():
            return data_demo_sample[column]

        lib_path = self._import_path
        if column.column_type() == KnotInfoColumnTypes.OnlyLinks:
            filename = self.filename.links.sobj_data(column)
        else:
            filename = self.filename.knots.sobj_data(column)

        verbose('loading data library %s ...' %(filename))
        res = load('%s/%s' %(lib_path, filename))
        verbose('... finished!')

        return res



column_demo_sample = {
    'name':                 ['Name',                 KnotInfoColumnTypes.KnotsAndLinks],
    'dt_notation':          ['DT Notation',          KnotInfoColumnTypes.OnlyKnots],
    'gauss_notation':       ['Gauss Notation',       KnotInfoColumnTypes.KnotsAndLinks],
    'pd_notation':          ['PD Notation',          KnotInfoColumnTypes.OnlyKnots],
    'pd_notation_vector':   ['PD Notation (vector)', KnotInfoColumnTypes.OnlyLinks],
    'crossing_number':      ['Crossing Number',      KnotInfoColumnTypes.KnotsAndLinks],
    'braid_index':          ['Braid Index',          KnotInfoColumnTypes.OnlyKnots],
    'braid_length':         ['Braid Length',         KnotInfoColumnTypes.OnlyKnots],
    'braid_notation':       ['Braid Notation',       KnotInfoColumnTypes.KnotsAndLinks],
    'alternating':          ['Alternating',          KnotInfoColumnTypes.KnotsAndLinks],
    'alexander_polynomial': ['Alexander',            KnotInfoColumnTypes.OnlyKnots],
    'jones_polynomial':     ['Jones',                KnotInfoColumnTypes.KnotsAndLinks],
    'conway_polynomial':    ['Conway',               KnotInfoColumnTypes.KnotsAndLinks],
    'homfly_polynomial':    ['HOMFLY',               KnotInfoColumnTypes.OnlyKnots],
    'kauffman_polynomial':  ['Kauffman',             KnotInfoColumnTypes.KnotsAndLinks],
    'symmetry_type':        ['Symmetry Type',        KnotInfoColumnTypes.OnlyKnots],
    'width':                ['Width',                KnotInfoColumnTypes.OnlyKnots],
    'homflypt_polynomial':  ['HOMFLYPT Polynomial',  KnotInfoColumnTypes.OnlyLinks],
    'arc_notation':         ['Arc Notation',         KnotInfoColumnTypes.OnlyLinks],
    'dt_code':              ['DT code',              KnotInfoColumnTypes.OnlyLinks]
}


row_demo_sample = {
    'K0_1':   [0, 1],
    'K3_1':   [1, 1],
    'K4_1':   [2, 1],
    'K5_1':   [3, 1],
    'K5_2':   [4, 1],
    'K6_1':   [5, 1],
    'K6_2':   [6, 1],
    'K6_3':   [7, 1],
    'K7_1':   [8, 1],
    'K7_2':   [9, 1],
    'L2a1_0': [10, 2],
    'L2a1_1': [11, 2],
    'L4a1_0': [12, 2],
    'L4a1_1': [13, 2],
    'L5a1_0': [14, 2],
    'L5a1_1': [15, 2],
    'L6a1_0': [16, 2],
    'L6a1_1': [17, 2],
    'L6a2_0': [18, 2],
    'L6a2_1': [19, 2]
}

db = KnotInfoDataBase()
dc = db.columns()

data_demo_sample = {
    dc.crossing_number: ['0', '3', '4', '5', '5', '6', '6', '6', '7', '7', '2', '2', '4', '4', '5', '5', '6', '6', '6', '6', '6'],
    dc.braid_notation: [
        '',
        '{1,1,1}',
        '{1,-2,1,-2}',
        '{1,1,1,1,1}',
        '{1,1,1,2,-1,2}',
        '{1,1,2,-1,-3,2,-3}',
        '{1,1,1,-2,1,-2}',
        '{1,1,-2,1,-2,-2}',
        '{1,1,1,1,1,1,1}',
        '{1,1,1,2,-1,2,3,-2,3}',
        '{2, {-1, -1}}',
        '{2, {1, 1}}',
        '{4, {1, -2, 3, -2, -1, -2, -3, -2}}',
        '{2, {1, 1, 1, 1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{3, {-1, 2, -1, 2, -1}}',
        '{4, {1, -2, 3, -2, 1, -2, -3, -2}}',
        '{4, {1, 2, 3, 2, 2, -1, 2, 2, -3, 2}}',
        '{4, {1, -2, -2, -2, 3, -2, -1, -2, -3, -2}}',
        '{4, {1, 2, -3, 2, -1, 2, 3, 2, 2, 2}}',
        '{2, {-1, -1, -1, -1, -1, -1}}'
        ],
    dc.braid_index: ['1', '2', '3', '2', '3', '4', '3', '3', '2', '4'],
    dc.braid_length: ['', '3', '4', '5', '6', '7', '6', '6', '7', '9'],
    dc.pd_notation: [
        '',
        '[[1,5,2,4],[3,1,4,6],[5,3,6,2]]',
        '[[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]',
        '[[2,8,3,7],[4,10,5,9],[6,2,7,1],[8,4,9,3],[10,6,1,5]]',
        '[[1,5,2,4],[3,9,4,8],[5,1,6,10],[7,3,8,2],[9,7,10,6]]',
        '[[1,7,2,6],[3,10,4,11],[5,3,6,2],[7,1,8,12],[9,4,10,5],[11,9,12,8]]',
        '[[1,8,2,9],[3,11,4,10],[5,1,6,12],[7,2,8,3],[9,7,10,6],[11,5,12,4]]',
        '[[4,2,5,1],[8,4,9,3],[12,9,1,10],[10,5,11,6],[6,11,7,12],[2,8,3,7]]',
        '[[1,9,2,8],[3,11,4,10],[5,13,6,12],[7,1,8,14],[9,3,10,2],[11,5,12,4],[13,7,14,6]]',
        '[[2,10,3,9],[4,14,5,13],[6,12,7,11],[8,2,9,1],[10,8,11,7],[12,6,13,5],[14,4,1,3]]',
        ],
    dc.pd_notation_vector: [
        '{{4, 1, 3, 2}, {2, 3, 1, 4}}',
        '{{4, 2, 3, 1}, {2, 4, 1, 3}}',
        '{{6, 1, 7, 2}, {8, 3, 5, 4}, {2, 5, 3, 6}, {4, 7, 1, 8}}',
        '{{6, 2, 7, 1}, {8, 4, 5, 3}, {2, 8, 3, 7}, {4, 6, 1, 5}}',
        '{{6, 1, 7, 2}, {10, 7, 5, 8}, {4, 5, 1, 6}, {2, 10, 3, 9}, {8, 4, 9, 3}}',
        '{{8, 2, 9, 1}, {10, 7, 5, 8}, {4, 10, 1, 9}, {2, 5, 3, 6}, {6, 3, 7, 4}}',
        '{{6, 1, 7, 2}, {10, 3, 11, 4}, {12, 8, 5, 7}, {8, 12, 9, 11}, {2, 5, 3, 6}, {4, 9, 1, 10}}',
        '{{10, 2, 11, 1}, {6, 4, 7, 3}, {12, 10, 5, 9}, {8, 6, 9, 5}, {2, 12, 3, 11}, {4, 8, 1, 7}}',
        '{{8, 1, 9, 2}, {12, 5, 7, 6}, {10, 3, 11, 4}, {4, 11, 5, 12}, {2, 7, 3, 8}, {6, 9, 1, 10}}',
        '{{10, 2, 11, 1}, {12, 6, 7, 5}, {8, 4, 9, 3}, {4, 8, 5, 7}, {2, 12, 3, 11}, {6, 10, 1, 9}}',
        '{{8, 1, 9, 2}, {2, 9, 3, 10}, {10, 3, 11, 4}, {12, 5, 7, 6}, {6, 7, 1, 8}, {4, 11, 5, 12}}'
        ],
    dc.dt_notation: [
        '',
        '[4, 6, 2]',
        '[4, 6, 8, 2]',
        '[6, 8, 10, 2, 4]',
        '[4, 8, 10, 2, 6]',
        '[4, 8, 12, 10, 2, 6]',
        '[4, 8, 10, 12, 2, 6]',
        '[4, 8, 10, 2, 12, 6]',
        '[8, 10, 12, 14, 2, 4, 6]',
        '[4, 10, 14, 12, 2, 8, 6]'
        ],
    dc.dt_code: [
        '[{4}, {2}]',
        '[{4}, {2}]',
        '[{6, 8}, {2, 4}]',
        '[{6, 8}, {4, 2}]',
        '[{6, 8}, {4, 10, 2}]',
        '[{8, 6}, {2, 10, 4}]',
        '[{6, 10}, {2, 12, 4, 8}]',
        '[{10, 6}, {8, 4, 12, 2}]',
        '[{8, 10, 12}, {2, 6, 4}]',
        '[{10, 8, 12}, {4, 6, 2}]',
        '[{8, 10, 12}, {6, 2, 4}]'
        ],
    dc.gauss_notation: [
        '',
        '{1, -2, 3, -1, 2, -3}',
        '{-1, 2, -3, 1, -4, 3, -2, 4}',
        '{-1, 2, -3, 4, -5, 1, -2, 3, -4, 5}',
        '{1, -2, 3, -1, 4, -5, 2, -3, 5, -4}',
        '{1, -2, 3, -4, 2, -1, 5, -6, 4, -3, 6, -5}',
        '{1, -2, 3, -4, 5, -6, 2, -1, 6, -3, 4, -5}',
        '{-1, 2, -3, 1, -4, 5, -2, 3, -6, 4, -5, 6}',
        '{1, -2, 3, -4, 5, -6, 7, -1, 2, -3, 4, -5, 6, -7}',
        '{-1, 2, -3, 4, -5, 6, -7, 1, -2, 7, -6, 5, -4, 3}',
        '{{1, -2}, {2, -1}}',
        '{{1, -2}, {2, -1}}',
        '{{1, -3, 2, -4}, {3, -1, 4, -2}}',
        '{{1, -3, 2, -4}, {4, -1, 3, -2}}',
        '{{1, -4, 5, -3}, {3, -1, 2, -5, 4, -2}}',
        '{{1, -4, 5, -3}, {4, -5, 2, -1, 3, -2}}',
        '{{1, -5, 2, -6}, {5, -1, 3, -4, 6, -2, 4, -3}}',
        '{{1, -5, 2, -6}, {4, -2, 6, -4, 3, -1, 5, -3}}',
        '{{1, -5, 3, -4, 2, -6}, {5, -1, 6, -3, 4, -2}}',
        '{{1, -5, 3, -4, 2, -6}, {4, -3, 6, -1, 5, -2}}',
        '{{1, -2, 3, -6, 4, -5}, {5, -1, 2, -3, 6, -4}}'
        ],
    dc.arc_notation: [
        '{{4, 2}, {3, 1}, {4, 2}, {1, 3}}',
        '{{2, 4}, {3, 1}, {2, 4}, {3, 1}}',
        '{{6, 4}, {3, 5}, {4, 2}, {1, 3}, {2, 6}, {5, 1}}',
        '{{3, 6}, {2, 5}, {6, 4}, {1, 3}, {5, 2}, {4, 1}}',
        '{{6, 2}, {1, 4}, {3, 5}, {4, 7}, {2, 6}, {7, 3}, {5, 1}}',
        '{{3, 5}, {6, 4}, {5, 2}, {7, 3}, {1, 6}, {2, 7}, {4, 1}}',
        '{{8, 4}, {3, 5}, {4, 2}, {6, 3}, {5, 7}, {1, 6}, {2, 8}, {7, 1}}',
        '{{2, 8}, {1, 7}, {8, 4}, {5, 3}, {4, 2}, {3, 6}, {7, 5}, {6, 1}}',
        '{{8, 3}, {2, 7}, {3, 1}, {4, 8}, {5, 2}, {6, 4}, {7, 5}, {1, 6}}',
        '{{3, 8}, {2, 7}, {8, 4}, {1, 3}, {5, 2}, {4, 6}, {7, 5}, {6, 1}}',
        '{{8, 2}, {1, 3}, {2, 4}, {3, 5}, {4, 6}, {5, 7}, {6, 8}, {7, 1}}'
        ],
    dc.alternating: ['Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'],
    dc.symmetry_type: [
        '',
        'reversible',
        'fully amphicheiral',
        'reversible',
        'reversible',
        'reversible',
        'reversible',
        'fully amphicheiral',
        'reversible',
        'reversible'
        ],
    dc.homfly_polynomial: [
        '',
        '(2*v^2-v^4)+ (v^2)*z^2',
        '(v^(-2)-1+ v^2)+ (-1)*z^2',
        '(3*v^4-2*v^6)+ (4*v^4-v^6)*z^2+ (v^4)*z^4',
        '(v^2+ v^4-v^6)+ (v^2+ v^4)*z^2',
        '(v^(-2)-v^2+ v^4)+ (-1-v^2)*z^2',
        '(2-2*v^2+ v^4)+ (1-3*v^2+ v^4)*z^2+ (-v^2)*z^4',
        '(-v^(-2)+ 3-v^2)+ (-v^(-2)+ 3-v^2)*z^2+ (1)*z^4',
        '(4*v^6-3*v^8)+ (10*v^6-4*v^8)*z^2+ (6*v^6-v^8)*z^4+ (v^6)*z^6',
        '(v^2+ v^6-v^8)+ (v^2+ v^4+ v^6)*z^2'
        ],
    dc.homflypt_polynomial: [
        '1/(v^3*z)-1/(v*z)-z/v',
        'v/z-v^3/z + v*z',
        '1/(v^5*z)-1/(v^3*z)-z/v^3-z/v',
        'v^3/z-v^5/z + 3*v^3*z-v^5*z + v^3*z^3',
        '1/(v*z)-v/z-z/v^3 + (2*z)/v-v*z + z^3/v',
        '1/(v*z)-v/z-z/v^3 + (2*z)/v-v*z + z^3/v',
        '1/(v^5*z)-1/(v^3*z)-(2*z)/v^3 + z/v-v*z + z^3/v',
        'v^3/z-v^5/z + 2*v^3*z + v^5*z-v^7*z + v^3*z^3 + v^5*z^3',
        '1/(v^7*z)-1/(v^5*z) + z/v^7-(2*z)/v^5-(2*z)/v^3-z^3/v^5-z^3/v^3',
        'v^5/z-v^7/z + 2*v^3*z + 2*v^5*z-v^7*z + v^3*z^3 + v^5*z^3',
        '1/(v^7*z)-1/(v^5*z) + (3*z)/v^7-(6*z)/v^5 + z^3/v^7-(5*z^3)/v^5-z^5/v^5'
        ]
}
