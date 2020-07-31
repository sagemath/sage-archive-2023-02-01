# -*- coding: utf-8 -*-
r"""
Cubic Hecke Database

This module contains the class :class:`CubicHeckeDataBase` which serves as an interface to
Ivan Marin's data files with respect to the cubic Hecke algebras. Furthermore, it contains
the class :class:`CubicHeckeFileCache` which enables :class:`CubicHeckeAlgebra` to keep
intermediate results of calculations in the file system.

AUTHORS:

- Sebastian Oehms May 2020: initial version
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
from enum import Enum

from sage.structure.sage_object import SageObject
from sage.misc.persist import db_save, db, save, load
from sage.misc.verbose import verbose
from sage.env import SAGE_SHARE, SAGE_ROOT
from sage.matrix.constructor import matrix, Matrix  # uppercase version used in Marin's file `MatricesRegH4.maple`
from sage.rings.integer_ring import ZZ
from sage.algebras.hecke_algebras.base_rings_of_definition.cubic_hecke_base_ring import CubicHeckeRingOfDefinition





#----------------------------------------------------------------------------------------------------------------------------
# functions to convert matrices and ring elements to and from flat python dictionaries in order to save matrices avoiding
# compatibility problems with older or newer sage versions and to save disc space
#----------------------------------------------------------------------------------------------------------------------------
# conversion of ring element to dictionary
#----------------------------------------------------------------------------------------------------------------------------
def convert_poly_to_dict_recursive(ring_elem):
    r"""
    Convert a ring element to a python dictionary recursively using the dict method of it.
    By recursion the dictionaries values are converted as well as long as they posses a
    ``dict`` method. If the values are sage Integers they are converted into python integers.

    INPUT:

    -- ``ring_elem`` - ring element to be converted into python dictionary

    OUTPUT:

    A python dictionary from which ``ring_elem`` can be reconstructed via element construction by recursion.
    The values of the dictionary may be dictionaries again if the parent of ring_elem has a base_ring
    different from itself.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import convert_poly_to_dict_recursive
        sage: L.<c>=LaurentPolynomialRing(ZZ, 'c')
        sage: P.<a,b> = L['a,b']; P
        Multivariate Polynomial Ring in a, b
          over Univariate Laurent Polynomial Ring in c over Integer Ring
        sage: elem = 5*b-3*a*~c; elem
        (-3*c^-1)*a + 5*b
        sage: convert_poly_to_dict_recursive(elem)
        {(0, 1): {0: 5}, (1, 0): {-1: -3}}
    """

    dict_res = {}
    if hasattr(ring_elem, 'dict'):
        ring_elem_dict = ring_elem.dict()
        for k in ring_elem_dict.keys():
            dict_res[k] = convert_poly_to_dict_recursive(ring_elem_dict[k])
    else:
        if ring_elem in ZZ:
            return int(ring_elem)
        return ring_elem
    return dict_res






#---------------------------------------------------------------------------------------------------------------------------
# conversion of matrix to dictionary
#---------------------------------------------------------------------------------------------------------------------------
def convert_mat_to_dict_recursive(mat):
    r"""
    Convert a matrix to a python dictionary using the dict method of it. Furthermore, the dictionaries
    values are converted as well by the convert_poly_to_dict_recursive function.

    INPUT:

    - ``mat`` -- matrix to be converted into python dictionary

    OUTPUT:

    A python dictionary from which mat can be reconstructed via element construction. The values of the
    dictionary may be dictionaries again if entries of the matrix have a dict method as well.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import convert_mat_to_dict_recursive
        sage: L.<c>=LaurentPolynomialRing(ZZ, 'c')
        sage: P.<a,b> = L['a,b']; P
        Multivariate Polynomial Ring in a, b
          over Univariate Laurent Polynomial Ring in c over Integer Ring
        sage: mat = matrix(P, [[2*a, -3], [c, 4*b*~c]]); mat
        [       2*a         -3]
        [         c (4*c^-1)*b]
        sage: convert_mat_to_dict_recursive(mat)
        {(0, 0): {(1, 0): {0: 2}},
        (0, 1): {(0, 0): {0: -3}},
        (1, 0): {(0, 0): {1: 1}},
        (1, 1): {(0, 1): {-1: 4}}}
    """

    mat_dict = {}
    mat_dict_temp = mat.dict()
    for k in mat_dict_temp.keys():
        mat_dict[k] = convert_poly_to_dict_recursive(mat_dict_temp[k]) 
    return mat_dict




class CubicHeckeDataFilename(Enum):
    r"""
    Enum for the different data files. The following choices are possible:

    - ``basis`` -- contains the basis for the cubic Hecke algebra up to 4 strands
    - ``regular_left`` -- contains the left regular representation matrices of the generators
    - ``regular_right`` -- contains the right regular representation matrices of the generators
    - ``irred_split`` -- contains representation matrices of the generators of the split irreducible
      representations


    Examples::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db.filename
        <enum 'CubicHeckeDataFilename'>
    """
    def download(self):
        """
        Return the file name to download the data from Ivan Marin's web-page.

        Examples::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.filename.basis.download()
            'baseH4.maple'
        """
        return self.value[0]
    def py(self):
        """
        Return the file name under which the data from Ivan Marin's web-page
        are stored as python file.

        Examples::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.filename.basis.py()
            'baseH4.maple.py'
        """
        return '%s.py' %(self.value[0])

    def sobj(self, nstrands=None):
        """
        Return the file name under which the data from Ivan Marin's web-page
        is converted into sobj-files.

        INPUT:

        - ``nstrands`` -- Integer number of strands of the underlying braid group
          if the data file depends on it. Otherwise use default None

        Examples::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.filename.basis.sobj()
            'monomial_basis.sobj'
            sage: cha_db.filename.basis.sobj()
            'monomial_basis.sobj'
            sage: cha_db.filename.irred_split.sobj(2)
            'irred_split_reprs_2.sobj'
            sage: cha_db.filename.regular_left.sobj(3)
            'regular_left_reprs_3.sobj'
        """
        if nstrands is None:
            return '%s.sobj' %(self.value[1])
        else:
            return '%s_%s.sobj' %(self.value[1], nstrands)

    basis         = ['baseH4.maple',             'monomial_basis']
    regular_left  = ['MatricesRegH4.maple',      'regular_left_reprs']
    regular_right = ['MatricesRegH4right.maple', 'regular_right_reprs']
    irred_split   = ['RepresentationsH25',       'irred_split_reprs']





#----------------------------------------------------------------------------------------------------------------------------
# Class to supply data for the basis and matrix representation for the cubic Hecke algebra
#----------------------------------------------------------------------------------------------------------------------------
class CubicHeckeDataBase(SageObject):
    r"""
    Database interface needed for :class:`CubicHeckeAlgebra`.

    The original data are obtained from Ivan Marin's web-page (URL see the example below). In order
    to have these data installed during the build process as a sage-package they are converted
    as python files into a tarball. This tarball has been created using the method :meth:`create_spkg_tarball`.

    EXAMPLES::

        sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
        sage: cha_db = CubicHeckeDataBase()
        sage: cha_db._url_marin
        'http://www.lamfa.u-picardie.fr/marin/softs/H4'
    """

    filename = CubicHeckeDataFilename

    def __init__(self):
        r"""
        Python constructor.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: from sage.env import SAGE_SHARE
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db._import_path_sobj == SAGE_SHARE + '/cubic_hecke_marin/sobj'
            True
        """
        self._url_marin              = 'http://www.lamfa.u-picardie.fr/marin/softs/H4'

        self._package = 'cubic_hecke_marin'
        version_file  = os.path.join(SAGE_ROOT, 'build/pkgs/%s/package-version.txt' %self._package)
        f = open(version_file)
        self._version = f.read().splitlines()[0]
        f.close()

        self._import_path      = os.path.join(SAGE_SHARE, self._package)
        self._import_path_py   = os.path.join(self._import_path, 'py')
        self._import_path_sobj = os.path.join(self._import_path, 'sobj')

        from sage.misc.misc import sage_makedirs
        sage_makedirs(self._import_path_py)
        sage_makedirs(self._import_path_sobj)

        self._data_library = {}

    def version(self):
        r"""
        Return the current version.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.version()
            '20200513'
        """
        return self._version

    def _create_python_file(self, filename):
        r"""
        Return the data fetched from Iwan Marin's homepage as a python file
        such that it can be loaded via `sage_eval`.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: load(cha_db._create_python_file(cha_db.filename.basis))    # not tested (because of internet access)
            Importing data for monomial_basis.sobj from http://www.lamfa.u-picardie.fr/marin/softs/H4/baseH4.maple
            sage: len(baseH4)                                                # not tested
            648
        """
        if not isinstance(filename, CubicHeckeDataBase.filename):
            raise TypeError('File name must be an instance of enum %s' (CubicHeckeDataBase.filename))

        import_file = '%s/%s' %(self._import_path_py, filename.py())

        # import directly from the internet page
        from six.moves.urllib.request import urlopen
        try:
            from urllib.error import HTTPError
        except ImportError:
            from urllib2 import HTTPError

        try:
            url = '%s/%s' %(self._url_marin, filename.download())
            url_data = urlopen(url).read().decode()
            print('Importing data for %s from %s' %(filename.sobj(), url))
            preparsed_data =url_data.replace(':=', '=').replace(';', '').replace('^', '**')
            f = open(import_file, 'wt')
            f.write(preparsed_data)
            f.close()
            return import_file
        except HTTPError:
            raise IOError('Data import file %s not found! Internet connection needed!' %(filename))

    def create_spkg_tarball(self):
        r"""
        Create a tarball for the sage-package ``cubic_heck_marin`` in the ``upstream`` directory. This
        utility should only be used by users who know what they do in case of a switch to a new
        version of the data files (that is if the original files on Iwan Marin's homepage have changed).
        In that case in invocation of ``sage -package fix-checksum cubic_hecke_marin`` will be necessary.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.create_spkg_tarball()    # not tested (because of internet access)
            Importing data for monomial_basis.sobj
               from http://www.lamfa.u-picardie.fr/marin/softs/H4/baseH4.maple
            Importing data for regular_left_reprs.sobj
               from http://www.lamfa.u-picardie.fr/marin/softs/H4/MatricesRegH4.maple
            Importing data for regular_right_reprs.sobj
               from http://www.lamfa.u-picardie.fr/marin/softs/H4/MatricesRegH4right.maple
            Importing data for irred_split_reprs.sobj
               from http://www.lamfa.u-picardie.fr/marin/softs/H4/RepresentationsH25
            py/
            py/MatricesRegH4.maple.py
            py/MatricesRegH4right.maple.py
            py/RepresentationsH25.py
            py/baseH4.maple.py
        """
        for filename in CubicHeckeDataBase.filename:
            self._create_python_file(filename)
        os.system('cd %s; tar -cvjSf %s/upstream/%s-%s.tar.bz2 py' %(self._import_path, SAGE_ROOT, self._package, self._version) )

    def import_data(self, filename, from_spkg=True):
        r"""
        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: load(cha_db.import_data(cha_db.filename.basis))
            sage: len(baseH4)
            648
        """
        if not isinstance(filename, CubicHeckeDataBase.filename):
            raise TypeError('File name must be an instance of enum %s' (CubicHeckeDataBase.filename))

        import_file = '%s/%s' %(self._import_path_py, filename.py())

        try:
            open(import_file)
            return import_file
        except IOError:
            if from_spkg:
                # import from the spkg tarball
                print('Importing cubic Hecke database from SPKG!')
                os.system('pwd')
                os.system('cp src/*.py %s' %(self._import_path_py))
                open(import_file)
                return import_file
            else:
                return self._create_python_file(filename)
                

    def create_static_db_marin_basis(self):
        r"""
        Create the basis of the cubic Hecke algebra according to the original
        data from Iwan Marin's home page.

        This method is called during the build procedure for the sage-package.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.create_static_db_marin_basis() 
        """
        global baseH4  # set by load
        load(self.import_data(self.filename.basis))

        basis_h1 = []
        basis_h2 = []
        basis_h3 = []
        basis_h4 = baseH4

        len_baseH4 = len(baseH4)

        ind_h1 = []
        ind_h2 = []
        ind_h3 = []
        ind_h4 = range(len_baseH4)

        for i in ind_h4:
            set_i = set(baseH4[i])
            if 3  not in set_i and -3  not in set_i:
                basis_h3.append( basis_h4[i] )
                ind_h3.append(i) 
                if 2  not in set_i and -2  not in set_i:
                    basis_h2.append( basis_h4[i] ) 
                    ind_h2.append(i) 
                    if 1  not in set_i and -1  not in set_i:
                        basis_h1.append( basis_h4[i] ) 
                        ind_h1.append(i) 

        # len_bas_h1 = len(basis_h1); len_bas_h2 = len(basis_h2); len_bas_h3 = len(basis_h3)

        basis = {1:[basis_h1, ind_h1], 2:[basis_h2, ind_h2], 3: [basis_h3, ind_h3], 4:[basis_h4, ind_h4]}
        save(basis, '%s/%s' %(self._import_path_sobj, self.filename.basis.sobj()) )


    def create_static_db_marin_regular(self, right=False):
        r"""
        Create the static data base for regular representations of the cubic
        Hecke algebra according to the original data from Iwan Marin's home page.

        This method is called during the build procedure for the sage-package.

        The invocations are not active in the doctest, since they cause a ``MemoryError``
        here. To refresh the data base you may call this method in a session. You will
        have to wait (maybe up to several minutes)!

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.create_static_db_marin_regular()             # not tested
            sage: cha_db.create_static_db_marin_regular(right=True)   # not tested
        """
        
        base_ring = CubicHeckeRingOfDefinition()
        global u, v, w, mm1, mm2, mm3, mm1I, mm2I, mm3I, reps      # set in load
        u, v, w = base_ring.gens_over_ground()

        if right == False:
            fname = self.filename.regular_left
        else:
            fname = self.filename.regular_right
        before = verbose('start loading %s' %fname)
        load(self.import_data(fname))
        before = verbose('end loading %s' %fname, t=before)

    
        def create_mat(ind_h, mat_h4):
            """
            Create restriction of regular representation of H4 to H1, H2 and H3
            """
 
            dim_mat = len(ind_h)
            mat = matrix(dim_mat, dim_mat, lambda i,j: mat_h4[ind_h[i], ind_h[j]])
            return convert_mat_to_dict_recursive(mat)            

        basis = self.read(self.filename.basis)
        ind_h1 = basis[1][1]
        ind_h2 = basis[2][1]
        ind_h3 = basis[3][1]
 
        representationH ={}
        representationH[0]  = [[create_mat(ind_h1, mm1)]]
        representationH[1]  = [[create_mat(ind_h2, mm1)]]
        representationH[2]  = [[create_mat(ind_h3, mm1), create_mat(ind_h3, mm2)]]
        representationH[3]  = [[convert_mat_to_dict_recursive(mm1), convert_mat_to_dict_recursive(mm2), convert_mat_to_dict_recursive(mm3)]]

        representationHI ={}
        representationHI[0] = [[create_mat(ind_h1, mm1I)]]
        representationHI[1] = [[create_mat(ind_h2, mm1I)]]
        representationHI[2] = [[create_mat(ind_h3, mm1I), create_mat(ind_h3, mm2I)]]
        representationHI[3] = [[convert_mat_to_dict_recursive(mm1I), convert_mat_to_dict_recursive(mm2I), convert_mat_to_dict_recursive(mm3I)]]
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import GenSign

        for i in range(4):
            if right == False:
                sobj_filename = '%s/%s' %(self._import_path_sobj, self.filename.regular_left.sobj(i+1))
            else:
                sobj_filename = '%s/%s' %(self._import_path_sobj, self.filename.regular_right.sobj(i+1))

            RegularMarinDict = {GenSign.pos:representationH[i], GenSign.neg:representationHI[i]}
            save(RegularMarinDict, sobj_filename )
        return


    def create_static_db_marin_split(self):
        r"""
        Create the static data base for split irreducible representations of the cubic
        Hecke algebra according to the original data from Iwan Marin's home page.

        This method is called during the build procedure for the sage-package.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: cha_db.create_static_db_marin_split()
        """
        # ------------------------------------------------------
        # Ivan Marin's data file uses a, b, c for the variables
        # corresponding to the eigenvalues of the cubic equation.
        # Therefore, we have to use them temporarily in that way
        # ------------------------------------------------------
        base_ring = CubicHeckeRingOfDefinition()
        extension_ring = base_ring.extension_ring()
        global a, b, c, j
        a, b, c = extension_ring.gens()
        j = extension_ring.cyclotomic_generator()

        load(self.import_data(self.filename.irred_split))

        a, b, c = base_ring.gens_over_ground() # now back to usual nameing
        cfs = [-c, b, -a, 1]
        cfse = [extension_ring(cf/c) for cf in cfs]

        def invert(matr):
           """
           Return inverse matrix for generators
           """

           matri = cfse[1]*matr.parent().one()
           matri += cfse[2]*matr
           matri += cfse[3]*matr**2
           d1, d2 = matr.dimensions()
           matrI = matrix(extension_ring, d1, d2, lambda i,j: extension_ring(matri[i,j]))
           return matrI

        # ------------------------------------------------------------------------------------------------
        # Restoring the split irreducibles from Iwan Marin's homepage
        # ------------------------------------------------------------------------------------------------

        anz_reps = len(reps)

        representation_h ={}
        representation_h[0]  = [[convert_mat_to_dict_recursive(Matrix(1,1,[extension_ring.one()]))]]
        representation_h[1]  = []
        representation_h[2]  = []
        representation_h[3]  = []

        representation_hI ={}
        representation_hI[0] = representation_h[0]
        representation_hI[1] = []
        representation_hI[2] = []
        representation_hI[3] = []
        for i in range( anz_reps ):
            repi = reps[i]
            if len(repi) != 3:
               raise RuntimeError( 'Error at position %d: three generators expected, got: %d' %( i, len(repi)))
            mt = []
            mtI = []
            for j in range(3):
                mat = matrix(repi[j])
                matI = invert(mat)
                mt.append(convert_mat_to_dict_recursive(mat))
                mtI.append(convert_mat_to_dict_recursive(matI))

            representation_h[3].append( mt )
            representation_hI[3].append( mtI )

            if i < 7:
                mt7 =  [ m for m in mt ]
                mt7I = [ m for m in mtI ]
                mt7.pop()
                mt7I.pop()
                representation_h[2].append( mt7 )
                representation_hI[2].append( mt7I )

            if i < 3:
                mt3 =  [ m for m in mt7 ]
                mt3I = [ m for m in mt7I ]
                mt3.pop()
                mt3I.pop()
                representation_h[1].append( mt3 )
                representation_hI[1].append( mt3I )

        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import GenSign
        for i in range(4):
            sobj_filename = '%s/%s' %(self._import_path_sobj, self.filename.irred_split.sobj(i+1))
            SplitIrredMarinDict = {GenSign.pos:representation_h[i], GenSign.neg:representation_hI[i]}
            save(SplitIrredMarinDict, sobj_filename)

        return

    # -------------------------------------------------------------------------------------------------------------
    # read from an sobj-file obtained from Ivan Marin's database
    # -------------------------------------------------------------------------------------------------------------
    def read(self, db_filename, nstrands=None):
        r"""
        Access various static data library.

        INPUT:

        ``db_filename`` -- instance of enum :class:`CubicHeckeDataBase.filename`
          to select the data to be read in

        OUTPUT:

        A dictionary containing the data corresponding to the db_filename.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: cha_db = CubicHeckeDataBase()
            sage: basis = cha_db.read(cha_db.filename.basis)
            sage: len(basis[3][0])
            24
        """
        if not isinstance(db_filename, CubicHeckeDataFilename):
            raise TypeError('db_filename must be an instance of enum %s' (CubicHeckeDataBase.filename))

        data_lib = self._data_library
        lib_path = self._import_path_sobj

        if (db_filename, nstrands) in data_lib.keys():
            return data_lib[(db_filename, nstrands)]

        verbose('loading data library %s ...' %(db_filename.sobj(nstrands=nstrands)))
        try:
            data_lib[(db_filename,nstrands)] = load('%s/%s' %(lib_path, db_filename.sobj(nstrands=nstrands)))
        except IOError:
            if db_filename == self.filename.basis:
                self.create_static_db_marin_basis()
            elif db_filename == self.filename.irred_split:
                self.create_static_db_marin_split()
            elif db_filename == self.filename.regular_right:
                self.create_static_db_marin_regular(right=True)
            else:
                self.create_static_db_marin_regular()
            data_lib[(db_filename, nstrands)] = load('%s/%s' %(lib_path, db_filename.sobj(nstrands=nstrands)))

        verbose('... finished!')

        return data_lib[(db_filename,nstrands)]


    # -------------------------------------------------------------------------------------------------------------
    # matrix_reprs_from_file_cache_
    # -------------------------------------------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, gen_ind, nstrands, ring_of_definition):
        r"""
        Return the matrix representations from the database.

        INPUT:

        - ``representation_type`` -- instance of enum :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation


        OUTPUT:

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeDataBase
            sage: CHA3 = algebras.CubicHecke(2)
            sage: GER = CHA3.extension_ring(generic=True)
            sage: cha_db = CHA3._database
            sage: rt = CHA3.repr_type
            sage: m1 =cha_db.read_matrix_representation(rt.SplitIrredMarin, 1, 3, GER)
            sage: len(m1)
            7
            sage: GBR = CHA3.base_ring(generic=True)
            sage: m1rl = cha_db.read_matrix_representation(rt.RegularLeft, 1, 3, GBR)
            sage: m1rl[0].dimensions()
            (24, 24)
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType, GenSign
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' (RepresentationType))

        num_rep = representation_type.number_of_representations(nstrands)
        rep_list = self.read(representation_type.data_filename(), nstrands=nstrands)
        if gen_ind > 0 :
            rep_list = [rep_list[GenSign.pos][i] for i in range(num_rep)]
            matrix_list = [matrix(ring_of_definition, rep[gen_ind-1 ], sparse=True) for rep in rep_list]
        else:
            # data of inverse of generators is stored under negative strand-index
            rep_list = [rep_list[GenSign.neg][i] for i in range(num_rep) ]
            matrix_list = [matrix(ring_of_definition, rep[-gen_ind-1 ], sparse=True) for rep in rep_list]
        for m in matrix_list: m.set_immutable()
        return matrix_list




class CubicHeckeFileCache(SageObject):
    """
    A class to cache calculations of the :class:`CubicHeckeAlgebra` in the local file system.
    """

    class section(Enum):
        r"""
        Enum for the different sections of file cache. The following choices are possible:

        - ``matrix_representations``  -- file cache for representation matrices of basis elements
        - ``braid_images``  -- file cache for images of braids
        - ``basis_extensions`` -- file cache for a dynamical growing basis used in the case of
          cubic Hecke algebras on more than 4 strands

        Examples::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(CHA2)
            sage: cha_fc.section
            <enum 'section'>
        """

        def filename(self, nstrands=None):
            r"""
            Return the file name under which the data of this file cache section
            is stored as an sobj-file.

            INPUT:

            - ``nstrands`` -- Integer number of strands of the underlying braid group
              if the data file depends on it. Otherwise use default None

            Examples::

                sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
                sage: CHA2 = algebras.CubicHecke(2)
                sage: cha_fc = CubicHeckeFileCache(CHA2)
                sage: cha_fc.section.matrix_representations.filename(2)
                'matrix_representations_2.sobj'
                sage: cha_fc.section.braid_images.filename(2)
                'braid_images_2.sobj'
            """
            if nstrands is None:
                return '%s.sobj' %(self.value)
            else:
                return '%s_%s.sobj' %(self.value, nstrands)

        matrix_representations  = 'matrix_representations'
        braid_images            = 'braid_images'
        basis_extensions        = 'basis_extensions'


    def __init__(self, num_strands):
        r"""
        Python constructor.

        INPUT:

        - ``cubic_hecke_algebra`` -- instance of :class:`CubicHeckeAlgebra`
          whose data should be cached in the file system.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc._file_cache_path == 'cubic_hecke'
            True
        """
        self._nstrands      = num_strands
        self._file_cache_path = 'cubic_hecke'
        self._data_library = {}

        from sage.misc.misc import sage_makedirs
        from sage.misc.persist import SAGE_DB
        sage_makedirs(os.path.join(SAGE_DB, self._file_cache_path))

    def reset_library(self, section=None):
        r"""
        Reset the file cache corresponding to the specified ``section``.

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
            sage: cha2_fc.reset_library(cha2_fc.section.matrix_representations)
            sage: data_mat = cha2_fc.read(cha2_fc.section.matrix_representations)
            sage: len(data_mat.keys())
            4
        """
        if section is None:
            for sec in self.section:
                self.reset_library(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' (CubicHeckeFileCache.section))
  
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        data_lib = self._data_library
        empty_dict = {}
        if  section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                new_dict={}
                empty_dict.update({rep_type:new_dict})
        elif  section == self.section.basis_extensions:
            empty_dict = []
        data_lib.update({section:empty_dict})


    def is_empty(self, section=None):
        r"""
        Return ``True`` if the cache of the given ``section`` is empty. 

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          to select the section of the file cache or ``None`` (default)
          meaning all sections
          

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library()
            sage: cha2_fc.is_empty()
            True
        """
        if section is None:
            return all(self.is_empty(section=sec) for sec in self.section)

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' (CubicHeckeFileCache.section))
  
        self.read(section)
        data_lib = self._data_library[section]
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if  section == self.section.matrix_representations:
            for rep_type in RepresentationType:
                if len(data_lib[rep_type]) > 0:
                    return False
            return True

        if  section == self.section.basis_extensions and self._nstrands > 4:
            # the new generators and their inverses are not counted
            # since they are added during initialization
            return len(data_lib) <= 2*(self._nstrands -4)
        return len(data_lib) == 0
       


    # -------------------------------------------------------------------------------------------------------------
    # save data file system
    # -------------------------------------------------------------------------------------------------------------
    def write(self, section=None):
        r"""
        Write data from memory to the file system.

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to.
          If omitted data of all sections is written to the file system

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.write(cha2_fc.section.braid_images)
        """
        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section is None:
            for sec in self.section:
                if sec in data_lib.keys():
                   self.write(section=sec)
            return

        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' (CubicHeckeFileCache.section))

        if section not in data_lib.keys():
            raise ValueError("No data for file %s in memory" %(section))

        db_save(data_lib[section], '%s/%s' %(lib_path, section.filename(self._nstrands)))


    # -------------------------------------------------------------------------------------------------------------
    # read from file system
    # -------------------------------------------------------------------------------------------------------------
    def read(self, section):
        r"""
        Read data into memory from the file system.
        

        INPUT:

        - ``section`` -- instance of enum :class:`CubicHeckeFileCache.section`
          specifying the section where the corresponding cached data belong to

        OUTPUT:

        Dictionary containing the data library corresponding to the section
        of file cache
        
        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: cha2_fc = CubicHeckeFileCache(2)
            sage: cha2_fc.reset_library(cha2_fc.section.braid_images)
            sage: cha2_fc.read(cha2_fc.section.braid_images)
            {}
        """
        if not isinstance(section, CubicHeckeFileCache.section):
            raise TypeError('section must be an instance of enum %s' (CubicHeckeFileCache.section))

        data_lib = self._data_library
        lib_path = self._file_cache_path

        if section in data_lib.keys():
            return data_lib[section]

        verbose('loading file cache %s ...' %(section))
        try:
            data_lib[section] = db('%s/%s' %(lib_path, section.filename(self._nstrands)))
            verbose('... finished!')
        except IOError:
            self.reset_library(section)
            verbose('... not found!')

        return data_lib[section]




    # -------------------------------------------------------------------------------------------------------------
    # matrix_reprs_from_file_cache_
    # -------------------------------------------------------------------------------------------------------------
    def read_matrix_representation(self, representation_type, monomial_tietze, ring_of_definition):
        r"""
        Return the matrix representations of the given monomial (in Tietze form)
        if it has been stored in the file cache before. Otherwise ``None`` is returned.

        INPUT:

        - ``representation_type`` -- instance of enum :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation

        - ``monomial_tietze`` -- tuple representing the braid in Tietze form

        - ``ring_of_definition`` -- instance of :class:`CubicHeckeRingOfDefinition` resp.
          :class:`CubicHeckeExtensionRing` (depending wether ``representation_type``
          is split or not)

        OUTPUT:

        Dictionary containing all matrix representations of ``self`` of the given representation_type
        which have been stored in the file cache.

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gt = g.Tietze()
            sage: rt = CHA2.repr_type
            sage: g.matrix(representation_type=rt.RegularLeft)
            [ 0 -v  1]
            [ 1  u  0]
            [ 0  w  0]
            sage: [_] == cha_fc.read_matrix_representation(rt.RegularLeft, gt, R)
            True
            sage: cha_fc.reset_library(cha_fc.section.matrix_representations)
            sage: cha_fc.write(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, gt, R) == None
            True
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' (RepresentationType))

        matrix_representations = self.read(self.section.matrix_representations)[representation_type]
        if monomial_tietze in matrix_representations.keys():
            matrix_list_dict = matrix_representations[monomial_tietze]
            matrix_list = [matrix(ring_of_definition, mat_dict, sparse=True) for mat_dict in matrix_list_dict]
            for m in matrix_list: m.set_immutable()
            return matrix_list
        return None




    # -------------------------------------------------------------------------------------------------------------
    # matrix_representation to file cache
    # -------------------------------------------------------------------------------------------------------------
    def write_matrix_representation(self, representation_type, monomial_tietze, matrix_list):
        r"""
        Write the matrix representation of a monomial to the file cache.

        INPUT:

        - ``representation_type`` -- instance of enum :class:`~sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep.RepresentationType`
          specifying the type of the representation

        - ``monomial_tietze`` -- tuple representing the braid in Tietze form

        - ``matrix_list`` -- list of matrices corresponding to the irreducible representations

        EXAMPLES::

            sage: CHA2 = algebras.CubicHecke(2)
            sage: R = CHA2.base_ring(generic=True)
            sage: cha_fc = CHA2._filecache
            sage: g, = CHA2.gens(); gi = ~g; git = gi.Tietze()
            sage: rt = CHA2.repr_type
            sage: m = gi.matrix(representation_type=rt.RegularRight)
            sage: cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            [
            [        0         1 (-w^-1)*u]
            [        0         0      w^-1]
            [        1         0  (w^-1)*v]
            ]
            sage: CHA2.reset_filecache(cha_fc.section.matrix_representations)
            sage: cha_fc.read_matrix_representation(rt.RegularLeft, git, R) == None
            True
            sage: cha_fc.write_matrix_representation(rt.RegularRight, git, [m])
            sage: [m] == cha_fc.read_matrix_representation(rt.RegularRight, git, R)
            True
        """
        from sage.algebras.hecke_algebras.matrix_representations.cubic_hecke_matrix_rep import RepresentationType
        if not isinstance(representation_type, RepresentationType):
            raise TypeError('representation_type must be an instance of enum %s' (RepresentationType))

        matrix_representations = self.read(self.section.matrix_representations)[representation_type]

        if monomial_tietze in matrix_representations.keys():
            # entry already registered
            return

        matrix_representation_dict = [convert_mat_to_dict_recursive(mat) for mat in list(matrix_list)]
        matrix_representations[monomial_tietze] = matrix_representation_dict

        self.write(self.section.matrix_representations)
        return

    # -------------------------------------------------------------------------------------------------------------
    # read braid images from file cache
    # -------------------------------------------------------------------------------------------------------------
    def read_braid_image(self, braid_tietze, ring_of_definition):
        r"""
        Return the list of pre calculated braid images from file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form

        - ``ring_of_definition`` -- instance of :class:`CubicHeckeRingOfDefinition`

        OUTPUT:

        A dictionary containing the pre calculated braid image of the given braid.

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b2 = b**2
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
            sage: b2_img = CHA2(b2); b2_img
            (-v) + u*c + w*c^-1
            sage: cha_fc.write_braid_image(b2.Tietze(), b2_img.to_vector())
            sage: cha_fc.read_braid_image(b2.Tietze(), ring_of_definition)
            (-v, u, w)
        """
        braid_images = self.read(self.section.braid_images)
        if braid_tietze in braid_images.keys():
            braid_image = braid_images[braid_tietze]
            result_list = [ring_of_definition(cf) for cf in list(braid_image)]
            from sage.modules.free_module_element import vector
            return vector(ring_of_definition, result_list)
        return None





    # -------------------------------------------------------------------------------------------------------------
    # braid image to_file cache
    # -------------------------------------------------------------------------------------------------------------
    def write_braid_image(self, braid_tietze, braid_image_vect):
        r"""
        Write the braid image of the given braid to the file cache.

        INPUT:

        - ``braid_tietze`` -- tuple representing the braid in Tietze form
        - ``braid_image_vect`` -- image of the given braid as a vector with respect
          to the basis of the cubic Hecke algebra

        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: ring_of_definition = CHA2.base_ring(generic=True)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: B2 = BraidGroup(2)
            sage: b, = B2.gens(); b3 = b**3
            sage: b3_img = CHA2(b3); b3_img
            (-u*v+w) + (u^2-v)*c + w*u*c^-1
            sage: cha_fc.write_braid_image(b3.Tietze(), b3_img.to_vector())
            sage: cha_fc.read_braid_image(b3.Tietze(), ring_of_definition)
            (-u*v + w, u^2 - v, w*u)
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.write(CubicHeckeFileCache.section.braid_images)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.braid_images)
            True
        """
        braid_images = self.read(self.section.braid_images)

        if braid_tietze in braid_images.keys():
            # entry already registered
            return

        braid_image_dict = [convert_poly_to_dict_recursive(cf) for cf in list(braid_image_vect)]
        braid_images[braid_tietze] = braid_image_dict

        self.write(self.section.braid_images)
        return


    # -------------------------------------------------------------------------------------------------------------
    # basis to file cache
    # -------------------------------------------------------------------------------------------------------------
    def update_basis_extensions(self, new_basis_extensions):
        r"""
        Update the file cache for basis extensions for cubic Hecke algebras on more than 4 strands
        according to the given ``new_basis_extensions``.

        INPUT:

        - ``new_basis_extensions`` -- list of additional (to the static basis) basis elements which should
          replace the former such list in the file.
 
        EXAMPLES::

            sage: from sage.databases.cubic_hecke_db import CubicHeckeFileCache
            sage: CHA2 = algebras.CubicHecke(2)
            sage: cha_fc = CubicHeckeFileCache(2)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
            sage: cha_fc.update_basis_extensions([(1,), (-1,)])
            sage: cha_fc.read(CubicHeckeFileCache.section.basis_extensions)
            [(1,), (-1,)]
            sage: cha_fc.reset_library(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.write(CubicHeckeFileCache.section.basis_extensions)
            sage: cha_fc.is_empty(CubicHeckeFileCache.section.basis_extensions)
            True
        """
        self._data_library.update({self.section.basis_extensions:new_basis_extensions})
        self.write(self.section.basis_extensions)
        return
