"""
Shared memory managers for F-symbol attributes.

This module provides an implementation for shared dictionary like
state attributes required by the orthogonal F-matrix solver.

Currently, the attributes only work when the base field of the FMatrix
factory is a cyclotomic field.
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

cimport cython
from cysignals.memory cimport sig_malloc
cimport numpy as np
from sage.combinat.root_system.poly_tup_engine cimport poly_to_tup, tup_fixes_sq
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.rings.polynomial.multi_polynomial_libsingular cimport MPolynomial_libsingular
from sage.rings.polynomial.polydict cimport ETuple

from multiprocessing import shared_memory
import numpy as np

cdef class KSHandler:
    def __init__(self, n_slots, field, use_mp=False, init_data={}, name=None):
        """
        Return a shared memory backed dict-like structure to manage the
        ``_ks`` attribute of an F-matrix factory object.

        This structure implements a representation of the known squares dictionary
        using a structured NumPy array backed by a contiguous shared memory
        object.

        The structure mimics a dictionary of ``(idx, known_sq)`` pairs. Each
        integer index corresponds to a variable and each ``known_sq`` is an
        element of the F-matrix factory's base cyclotomic field.

        Each cyclotomic coefficient is stored as a list of numerators and a
        list of denominators representing the rational coefficients. The
        structured array also maintains ``known`` attribute that indicates
        whether the structure contains an entry corresponding to the given index.

        The parent process should construct this object without a
        ``name`` attribute. Children processes use the ``name`` attribute,
        accessed via ``self.shm.name`` to attach to the shared memory block.

        INPUT:

        - ``n_slots`` -- The total number of F-symbols.
        - ``field`` -- F-matrix factory's base cyclotomic field.
        - ``use_mp`` -- a boolean indicating whether to construct a shared
          memory block to back ``self``.
        - ``name`` -- the name of a shared memory object
          (used by child processes for attaching).
        - ``init_data`` -- a dictionary or :class:`KSHandler` object containing
          known squares for initialization, e.g. from a solver checkpoint.

        .. NOTE::

            To properly dispose of shared memory resources,
            ``self.shm.unlink()`` must be called before exiting.

        .. WARNING::

            This structure does *not* cannot modify an entry that
            has already been set.

        EXAMPLES::

            sage: from sage.combinat.root_system.shm_managers import KSHandler
            sage: #Create shared data structure
            sage: f = FMatrix(FusionRing("A1",2), inject_variables=True)
            creating variables fx1..fx14
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: n = f._poly_ring.ngens()
            sage: ks = KSHandler(n,f._field,use_mp=True)
            sage: #In the same shell or in a different shell, attach to fvars
            sage: ks2 = KSHandler(n,f._field,name=ks.shm.name)
            sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
            sage: eqns = [fx1**2 - 4, fx3**2 + f._field.gen()**4 - 1/19*f._field.gen()**2]
            sage: ks.update([poly_to_tup(p) for p in eqns])
            sage: for idx, sq in ks.items():
            ....:     print("Index: {}, square: {}".format(idx, sq))
            ....:
            Index: 1, square: 4
            Index: 3, square: -zeta32^4 + 1/19*zeta32^2
            sage: ks.shm.unlink()
        """
        cdef int n, d
        self.field = field
        n = n_slots
        d = self.field.degree()
        ks_t = np.dtype([
            ('known', 'bool', (1,)),
            ('nums','i8',(d,)),
            ('denoms','u8',(d,))
        ])
        self.obj_cache = [None]*n
        if name is None:
            if use_mp:
                self.shm = shared_memory.SharedMemory(create=True,size=n*ks_t.itemsize)
                self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)
            else:
                self.ks_dat = np.ndarray((n,),dtype=ks_t)
            self.ks_dat['known'] = np.zeros((n,1),dtype='bool')
            self.ks_dat['nums'] = np.zeros((n,d),dtype='i8')
            self.ks_dat['denoms'] = np.ones((n,d),dtype='u8')
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)
        #Populate initializer data
        for idx, sq in init_data.items():
            self.setitem(idx,sq)

    @cython.nonecheck(False)
    @cython.wraparound(False)
    cdef NumberFieldElement_absolute get(self, int idx):
        """
        Retrieve the known square corresponding to the given index,
        if it exists.
        """
        if not self.ks_dat['known'][idx]:
            raise KeyError('Index {} does not correspond to a known square'.format(idx))
        if self.obj_cache[idx] is not None:
            return self.obj_cache[idx]
        cdef unsigned int i, d
        cdef Integer num, denom
        cdef Rational quo
        cdef list rat = list()
        cdef NumberFieldElement_absolute cyc_coeff
        d = self.field.degree()
        for i in range(d):
            num = Integer(self.ks_dat['nums'][idx,i])
            denom = Integer(self.ks_dat['denoms'][idx,i])
            quo = num / denom
            rat.append(quo)
        cyc_coeff = self.field(rat)
        self.obj_cache[idx] = cyc_coeff
        return cyc_coeff

    cpdef update(self, list eqns):
        r"""
        Update ```self``'s ``shared_memory``-backed dictionary of known
        squares. Keys are variable indices and corresponding values
        are the squares.

        EXAMPLES::

            sage: f = FMatrix(FusionRing("B5",1))
            sage: f._reset_solver_state()
            sage: for idx, sq in f._ks.items():
            ....:     k
            ....:
            sage: f.get_orthogonality_constraints()
            [fx0^2 - 1,
             fx1^2 - 1,
             fx2^2 - 1,
             fx3^2 - 1,
             fx4^2 - 1,
             fx5^2 - 1,
             fx6^2 - 1,
             fx7^2 - 1,
             fx8^2 - 1,
             fx9^2 - 1,
             fx10^2 + fx12^2 - 1,
             fx10*fx11 + fx12*fx13,
             fx10*fx11 + fx12*fx13,
             fx11^2 + fx13^2 - 1]
             sage: f.get_orthogonality_constraints(output=False)
             sage: f._ks.update(f.ideal_basis)
             sage: for idx, sq in f._ks.items():
             ....:     print(idx, "-->", sq)
             ....:
             0 --> 1
             1 --> 1
             2 --> 1
             3 --> 1
             4 --> 1
             5 --> 1
             6 --> 1
             7 --> 1
             8 --> 1
             9 --> 1

        .. WARNING::

            This method assumes every polynomial in ``eqns`` is *monic*.
        """
        cdef unsigned int i, idx
        cdef ETuple lm
        cdef list rhs
        cdef tuple eq_tup
        for i in range(len(eqns)):
            eq_tup = eqns[i]
            if tup_fixes_sq(eq_tup):
                rhs = [-v for v in eq_tup[-1][1]]
                #eq_tup is guaranteed univariate, so we extract variable idx from lm
                lm = eq_tup[0][0]
                idx = lm._data[0]
                self.setitem(idx, rhs)

    @cython.nonecheck(False)
    @cython.wraparound(False)
    @cython.infer_types(False)
    cdef setitem(self, int idx, rhs):
        """
        Create an entry corresponding to the given index.

        The ``rhs`` parameter may be a cyclotomic coefficient or its
        list/tuple representation.
        """
        cdef unsigned int i
        cdef long long num
        cdef unsigned long long denom
        self.ks_dat['known'][idx] = True
        if not isinstance(rhs, list):
            rhs = rhs._coefficients()
        for i in range(len(rhs)):
            num = <long long>(rhs[i].numerator())
            denom = <unsigned long long>(rhs[i].denominator())
            self.ks_dat['nums'][idx,i] = num
            self.ks_dat['denoms'][idx,i] = denom

    cdef bint contains(self, int idx):
        """
        Determine whether ``self`` contains entry corresponding to given
        ``idx``.
        """
        return self.ks_dat[idx]['known']

    def __eq__(self, KSHandler other):
        """
        Test for equality.

        TESTS::

            sage: f = FMatrix(FusionRing("C2",2))
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: from sage.combinat.root_system.shm_managers import KSHandler
            sage: n = f._poly_ring.ngens()
            sage: ks = KSHandler(n,f._field,use_mp=True,init_data=f._ks)
            sage: #In the same shell or in a different one, attach to shared memory handler
            sage: k2 = KSHandler(n,f._field,name=ks.shm.name)
            sage: ks == k2
            True
            sage: ks.shm.unlink()
        """
        ret = True
        for idx, sq in self.items():
            ret &= other.get(idx) == sq
        return ret

    def __reduce__(self):
        """
        Provide pickling / unpickling support for ``self.``

        TESTS::

            sage: f = FMatrix(FusionRing("A3",1))
            sage: f._reset_solver_state()
            sage: loads(dumps(f._ks)) == f._ks
            True
            sage: f.find_orthogonal_solution(verbose=False)    #long time
            sage: loads(dumps(f._ks)) == f._ks
            True
        """
        d = {i: sq for i, sq in self.items()}
        return make_KSHandler, (self.ks_dat.size,self.field,d)

    def items(self):
        """
        Iterate through existing entries using Python dict-style syntax.

        EXAMPLES::

            sage: f = FMatrix(FusionRing("A3",1))
            sage: f._reset_solver_state()
            sage: f.get_orthogonality_constraints(output=False)
            sage: f._ks.update(f.ideal_basis)
            sage: for idx, sq in f._ks.items():
            ....:     print("Index: {}, sq: {}".format(idx,sq))
            ....:
            Index: 0, sq: 1
            Index: 1, sq: 1
            Index: 2, sq: 1
            Index: 3, sq: 1
            Index: 4, sq: 1
            ...
            Index: 25, sq: 1
            Index: 26, sq: 1
        """
        cdef unsigned int i
        for i in range(self.ks_dat.size):
            if self.ks_dat['known'][i]:
                yield i, self.get(i)

def make_KSHandler(n_slots,field,init_data):
    """
    Provide pickling / unpickling support for :class:`KSHandler`.

    TESTS::

        sage: f = FMatrix(FusionRing("B4",1))
        sage: f._reset_solver_state()
        sage: loads(dumps(f._ks)) == f._ks                 #indirect doctest
        True
        sage: f.find_orthogonal_solution(verbose=False)    #long time
        sage: loads(dumps(f._ks)) == f._ks                 #indirect doctest
        True
    """
    return KSHandler(n_slots,field,init_data=init_data)

cdef class FvarsHandler:
    def __init__(self,n_slots,field,idx_to_sextuple,use_mp=False,name=None,init_data={},max_terms=20):
        """
        Return a shared memory backed dict-like structure to manage the
        ``_fvars`` attribute of an F-matrix factory object.

        This structure implements a representation of the F-symbols dictionary
        using a structured NumPy array backed by a contiguous shared memory
        object.

        The monomial data is stored in the ``exp_data`` structure. Monomial
        exponent data is stored contiguously and ``ticks`` are used to
        indicate different monomials.

        Coefficient data is stored in the ``coeff_nums`` and ``coeff_denom``
        arrays. The ``coeff_denom`` array stores the value
        ``d = coeff.denominator()`` for each cyclotomic coefficient. The
        ``coeff_nums`` array stores the values
        ``c.numerator() * d for c in coeff._coefficients()``, the abridged
        list representation of the cyclotomic coefficient ``coeff``.

        Each entry also has a boolean ``modified`` attribute, indicating
        whether it has been modified by the parent process. Entry retrieval
        is cached in each process, so each process must check whether
        entries have been modified before attempting retrieval.

        The parent process should construct this object without a
        ``name`` attribute. Children processes use the ``name`` attribute,
        accessed via ``self.shm.name`` to attach to the shared memory block.

        INPUT:

        - ``factory`` -- an F-matrix factory object.
        - ``name`` -- the name of a shared memory object
          (used by child processes for attaching).
        - ``max_terms`` -- maximum number of terms in each entry. Since
          we use contiguous C-style memory blocks, the size of the block
          must be known in advance.

        .. NOTE::

            To properly dispose of shared memory resources,
            ``self.shm.unlink()`` must be called before exiting.

        .. WARNING::

            The current data structure supports up to 2**16-1 entries,
            each with each monomial in each entry having at most 254
            nonzero terms. On average, each of the ``max_terms`` monomials
            can have at most 50 terms.

        EXAMPLES::

            sage: from sage.combinat.root_system.shm_managers import FvarsHandler
            sage: #Create shared data structure
            sage: f = FMatrix(FusionRing("A2",1), inject_variables=True)
            creating variables fx1..fx8
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7
            sage: fvars = FvarsHandler(8,f._field,f._idx_to_sextuple,use_mp=True)
            sage: #In the same shell or in a different shell, attach to fvars
            sage: fvars2 = FvarsHandler(8,f._field,f._idx_to_sextuple,name=fvars.shm.name)
            sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
            sage: fvars[f2, f1, f2, f2, f0, f0] = poly_to_tup(fx5**5)
            sage: f._tup_to_fpoly(fvars2[f2, f1, f2, f2, f0, f0])
            fx5^5
            sage: fvars.shm.unlink()
        """
        self.field = field
        cdef int d = self.field.degree()
        self.obj_cache = dict()
        self.fvars_t = np.dtype([
            ('modified','bool',(1,)),
            ('ticks', 'u1', (max_terms,)),
            ('exp_data', 'u2', (max_terms*50,)),
            ('coeff_nums','i4',(max_terms,d)),
            ('coeff_denom','u4',(max_terms,))
        ])
        self.sext_to_idx = {s: i for i, s in idx_to_sextuple.items()}
        self.ngens = n_slots
        if name is None:
            if use_mp:
                self.shm = shared_memory.SharedMemory(create=True,size=self.ngens*self.fvars_t.itemsize)
                self.fvars = np.ndarray((self.ngens,),dtype=self.fvars_t,buffer=self.shm.buf)
            else:
                self.fvars = np.ndarray((self.ngens,),dtype=self.fvars_t)
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.fvars = np.ndarray((self.ngens,),dtype=self.fvars_t,buffer=self.shm.buf)
        #Populate with initialziation data
        for sextuple, fvar in init_data.items():
            if isinstance(fvar, MPolynomial_libsingular):
                fvar = poly_to_tup(fvar)
            if isinstance(fvar, NumberFieldElement_absolute):
                fvar = ((ETuple({},self.ngens), fvar),)
            self[sextuple] = fvar

    cdef clear_modified(self):
        """
        Reset tagged entries modified by parent process.
        """
        self.fvars['modified'][:] = False

    @cython.nonecheck(False)
    @cython.wraparound(False)
    def __getitem__(self, sextuple):
        """
        Retrieve a record from the shared memory data structure by
        unflattening its representation and constructing relevant Python
        objects.

        This method returns a tuple of ``(ETuple, coeff)`` pairs,
        where ``coeff`` is an element of ``self.field``.

        EXAMPLES::

            sage: from sage.combinat.root_system.shm_managers import FvarsHandler
            sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
            sage: f = FMatrix(FusionRing("B7", 1), inject_variables=True)
            creating variables fx1..fx14
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
            sage: fvars = FvarsHandler(14,f._field,f._idx_to_sextuple,use_mp=True)
            sage: fvars[(f1, f2, f1, f2, f2, f2)] = poly_to_tup(1/8*fx0**15 - 23/79*fx2*fx13**3 - 799/2881*fx1*fx2**5*fx10)
            sage: fvars[f2, f2, f2, f2, f0, f0] = poly_to_tup(f._poly_ring.zero())
            sage: fvars[f2, f1, f2, f1, f2, f2] = poly_to_tup(-1/19*f._poly_ring.one())
            sage: s, t, r = (f1, f2, f1, f2, f2, f2), (f2, f2, f2, f2, f0, f0), (f2, f1, f2, f1, f2, f2)
            sage: f._tup_to_fpoly(fvars[s]) == 1/8*fx0**15 - 23/79*fx2*fx13**3 - 799/2881*fx1*fx2**5*fx10
            True
            sage: f._tup_to_fpoly(fvars[t]) == 0
            True
            sage: f._tup_to_fpoly(fvars[r]) == -1/19
            True
            sage: fvars.shm.unlink()

        .. NOTE::

            This method implements caching. Only the parent process is allowed
            to modify the shared fvars structure. Each process builds its own
            cache, so each process must update its cache before retrieving a
            modified entry, tagged via its ``modified`` property.
        """
        if not sextuple in self.sext_to_idx:
            raise KeyError('Invalid sextuple {}'.format(sextuple))
        cdef int idx = self.sext_to_idx[sextuple]
        # if idx in self.obj_cache:
        #     if self.fvars['modified'][idx]:
        #         del self.obj_cache[idx]
        #     else:
        #         return self.obj_cache[idx]
        cdef ETuple e = ETuple({}, self.ngens)
        cdef unsigned int cum, i, j, k, nnz
        cdef Integer d, num
        cdef list poly_tup, rats
        cdef Rational quo
        cdef tuple ret
        poly_tup = list()
        cum = 0
        for i in range(np.count_nonzero(self.fvars['ticks'][idx])):
            #Construct new ETuple for each monomial
            exp = e._new()
            #Handle constant coeff
            nnz = self.fvars['ticks'][idx,i] if self.fvars['ticks'][idx,i] < 255 else 0
            exp._nonzero = nnz
            exp._data = <int*>sig_malloc(sizeof(int)*nnz*2)
            for j in range(2*nnz):
                exp._data[j] = self.fvars['exp_data'][idx,cum]
                cum += 1

            #Construct cyclotomic field coefficient
            d = Integer(self.fvars['coeff_denom'][idx,i])
            rats = list()
            for k in range(self.field.degree()):
                num = Integer(self.fvars['coeff_nums'][idx,i,k])
                quo = num / d
                rats.append(quo)
            cyc_coeff = self.field(rats)

            poly_tup.append((exp, cyc_coeff))
        ret = tuple(poly_tup)
        # self.obj_cache[idx] = ret
        return ret

    @cython.nonecheck(False)
    @cython.wraparound(False)
    def __setitem__(self, sextuple, fvar):
        """
        Given a sextuple of labels and a tuple of ``(ETuple, cyc_coeff)`` pairs,
        create or overwrite an entry in the shared data structure
        corresponding to the given sextuple.

        EXAMPLES::

            sage: from sage.combinat.root_system.shm_managers import FvarsHandler
            sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
            sage: f = FMatrix(FusionRing("A3", 1), inject_variables=True)
            creating variables fx1..fx27
            Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13, fx14, fx15, fx16, fx17, fx18, fx19, fx20, fx21, fx22, fx23, fx24, fx25, fx26
            sage: fvars = FvarsHandler(27,f._field,f._idx_to_sextuple,use_mp=True)
            sage: fvars[(f3, f2, f1, f2, f1, f3)] = poly_to_tup(1/8*fx0**15 - 23/79*fx2*fx21**3 - 799/2881*fx1*fx2**5*fx10)
            sage: fvars[f3, f2, f3, f0, f1, f1] = poly_to_tup(f._poly_ring.zero())
            sage: fvars[f3, f3, f3, f1, f2, f2] = poly_to_tup(-1/19*f._poly_ring.one())
            sage: s, t, r = (f3, f2, f1, f2, f1, f3), (f3, f2, f3, f0, f1, f1), (f3, f3, f3, f1, f2, f2)
            sage: f._tup_to_fpoly(fvars[s]) == 1/8*fx0**15 - 23/79*fx2*fx21**3 - 799/2881*fx1*fx2**5*fx10
            True
            sage: f._tup_to_fpoly(fvars[t]) == 0
            True
            sage: f._tup_to_fpoly(fvars[r]) == -1/19
            True
            sage: fvars.shm.unlink()
        """
        cdef unsigned int cum, i, j, k, idx
        cdef unsigned long denom
        cdef long c
        cdef ETuple exp
        cdef NumberFieldElement_absolute coeff
        idx = self.sext_to_idx[sextuple]
        #Clear entry before inserting
        self.fvars[idx] = np.zeros((1,), dtype=self.fvars_t)
        cum = 0
        i = 0
        for exp, coeff in fvar:
            #Handle constant coefficient
            if exp._nonzero > 0:
                self.fvars['ticks'][idx,i] = exp._nonzero
            else:
                self.fvars['ticks'][idx,i] = -1
            for j in range(2*exp._nonzero):
                self.fvars['exp_data'][idx,cum] = exp._data[j]
                cum += 1
            denom = coeff.denominator()
            self.fvars['coeff_denom'][idx,i] = denom
            for k, r in enumerate(coeff._coefficients()):
                c = r * denom
                self.fvars['coeff_nums'][idx,i,k] = c
            i += 1
        self.fvars['modified'][idx] = True

    def __reduce__(self):
        """
        Provide pickling / unpickling support for ``self.``

        TESTS::

            sage: f = FMatrix(FusionRing("F4",1))
            sage: from sage.combinat.root_system.shm_managers import FvarsHandler
            sage: n = f._poly_ring.ngens()
            sage: fvars = FvarsHandler(n,f._field,f._idx_to_sextuple,init_data=f._fvars,use_mp=True)
            sage: for s, fvar in loads(dumps(fvars)).items():
            ....:     assert f._fvars[s] == f._tup_to_fpoly(fvar)
            ....:
            sage: fvars.shm.unlink()
        """
        n = self.fvars.size
        idx_map = {i: s for s, i in self.sext_to_idx.items()}
        d = {s: fvar for s, fvar in self.items()}
        return make_FvarsHandler, (n,self.field,idx_map,d)

    def items(self):
        """
        Iterates through key-value pairs in the data structure as if it
        were a Python dict.

        As in a Python dict, the key-value pairs are yielded in no particular
        order.

        EXAMPLES::

            sage: f = FMatrix(FusionRing("G2", 1), inject_variables=True)
            creating variables fx1..fx5
            Defining fx0, fx1, fx2, fx3, fx4
            sage: f.start_worker_pool()
            sage: for sextuple, fvar in f._shared_fvars.items():
            ....:     if sextuple == (f1, f1, f1, f1, f1, f1):
            ....:         f._tup_to_fpoly(fvar)
            ....:
            fx4
            sage: f.shutdown_worker_pool()
        """
        for sextuple in self.sext_to_idx:
            yield sextuple, self[sextuple]

def make_FvarsHandler(n,field,idx_map,init_data):
    """
    Provide pickling / unpickling support for :class:`FvarsHandler`.

    TESTS::

        sage: f = FMatrix(FusionRing("G2",1))
        sage: from sage.combinat.root_system.shm_managers import FvarsHandler
        sage: n = f._poly_ring.ngens()
        sage: fvars = FvarsHandler(n,f._field,f._idx_to_sextuple,init_data=f._fvars,use_mp=True)
        sage: for s, fvar in loads(dumps(fvars)).items():        # indirect doctest
        ....:     assert f._fvars[s] == f._tup_to_fpoly(fvar)
        ....:
        sage: fvars.shm.unlink()
    """
    return FvarsHandler(n,field,idx_map,init_data=init_data)
