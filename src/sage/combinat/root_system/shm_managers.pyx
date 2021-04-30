"""
Shared memory managers for F-symbol attributes
"""
# ****************************************************************************
#  Copyright (C) 2021 Guillermo Aboumrad <gh_willieab>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from cysignals.memory cimport sig_malloc
cimport numpy as np
from sage.rings.polynomial.polydict cimport ETuple
from sage.rings.number_field.number_field_base cimport NumberField

from multiprocessing import shared_memory
from sage.misc.cachefunc import cached_method
import numpy as np

class KSHandler():
    # cdef field
    # cdef np.ndarray ks_dat
    # cdef public shm

    def __init__(self, factory, name=None):
        self.field = factory._field
        n = factory._poly_ring.ngens()
        d = self.field.degree()
        ks_t = np.dtype([
            ('known', 'bool', (1,)),
            ('nums','i8',(d,)),
            ('denoms','u8',(d,))
        ])
        if name is None:
            self.shm = shared_memory.SharedMemory(create=True,size=n*ks_t.itemsize)
            self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)
            self.ks_dat['known'] = np.zeros((n,1),dtype='bool')
            self.ks_dat['nums'] = np.zeros((n,d),dtype='i4')
            self.ks_dat['denoms'] = np.ones((n,d),dtype='u4')
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)

    @cached_method
    def __getitem__(self, idx):
        if not self.ks_dat['known'][idx]:
            raise KeyError('Index {} does not correspond to a known square'.format(idx))
        cdef unsigned int i, d
        cdef list rat = list()
        d = self.field.degree()
        for i in range(d):
            rat.append(self.ks_dat['nums'][idx][i] / self.ks_dat['denoms'][idx][i])
        return self.field(rat)

    def __setitem__(self, idx, rhs):
        self.ks_dat['known'][idx] = True
        cdef unsigned int i
        for i in range(len(rhs)):
            self.ks_dat['nums'][idx][i] = rhs[i].numerator()
            self.ks_dat['denoms'][idx][i] = rhs[i].denominator()

    def __contains__(self, idx):
        return self.ks_dat['known'][idx]

    def __len__(self):
        """
        Compute the number of known squares.

        """
        return self.ks_dat['known'].sum()

    def __eq__(self, other):
        """
        Test for equality.
        """
        return np.all(self.ks_dat == other.ks_dat)

    def items(self):
        cdef unsigned int i
        for i in range(self.ks_dat.size):
            if self.ks_dat['known'][i]:
                yield i, self[i]

cdef class FvarsHandler():
    cdef dict sext_to_idx, obj_cache
    cdef unsigned int ngens
    cdef fvars_t
    cdef NumberField field
    cdef public np.ndarray fvars
    cdef public shm

    def __init__(self, factory, name=None, max_terms=20):
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
            sage: fvars = FvarsHandler(f)
            sage: #In the same shell or in a different shell, attach to fvars
            sage: fvars2 = FvarsHandler(f, name=fvars.shm.name)
            sage: from sage.combinat.root_system.poly_tup_engine import poly_to_tup
            sage: fvars[f2, f1, f2, f2, f0, f0] = poly_to_tup(fx5**5)
            sage: f._tup_to_fpoly(fvars2[f2, f1, f2, f2, f0, f0])
            fx5^5
            sage: fvars.shm.unlink()
        """
        n = len(factory._fvars)
        d = factory._field.degree()
        self.obj_cache = dict()
        self.fvars_t = np.dtype([
            ('modified','bool',(1,)),
            ('ticks', 'u1', (max_terms,)),
            ('exp_data', 'u2', (max_terms*50,)),
            ('coeff_nums','i4',(max_terms,d)),
            ('coeff_denom','u4',(max_terms,))
        ])
        self.sext_to_idx = {s: i for i, s in factory._idx_to_sextuple.items()}
        self.ngens = factory._poly_ring.ngens()
        self.field = factory._field
        if name is None:
            self.shm = shared_memory.SharedMemory(create=True,size=n*self.fvars_t.itemsize)
            self.fvars = np.ndarray((n,),dtype=self.fvars_t,buffer=self.shm.buf)
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.fvars = np.ndarray((n,),dtype=self.fvars_t,buffer=self.shm.buf)

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
            sage: fvars = FvarsHandler(f)
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
        cdef int cum, denom, i, j, k, idx
        cdef ETuple exp
        idx = self.sext_to_idx[sextuple]
        #Clear entry before inserting
        self.fvars[idx] = np.zeros((1,), dtype=self.fvars_t)
        cum = 0
        i = 0
        for exp, coeff in fvar:
            #Handle constant coefficient
            if exp._nonzero > 0:
                self.fvars['ticks'][idx][i] = exp._nonzero
            else:
                self.fvars['ticks'][idx][i] = -1
            for j in range(2*exp._nonzero):
                self.fvars['exp_data'][idx][cum] = exp._data[j]
                cum += 1
            denom = coeff.denominator()
            self.fvars['coeff_denom'][idx][i] = denom
            for k, c in enumerate(coeff._coefficients()):
                self.fvars['coeff_nums'][idx][i,k] = c * denom
            i += 1
        self.fvars['modified'][idx] = True

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
            sage: fvars = FvarsHandler(f)
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
        if idx in self.obj_cache:
            if self.fvars['modified'][idx]:
                del self.obj_cache[idx]
            else:
                return self.obj_cache[idx]
        cdef ETuple e = ETuple({}, self.ngens)
        cdef unsigned int cum, d, i, j, nnz
        cdef list poly_tup = list()
        cdef tuple ret
        cum = 0
        for i in range(np.count_nonzero(self.fvars['ticks'][idx])):
            #Construct new ETuple for each monomial
            exp = e._new()
            #Handle constant coeff
            nnz = self.fvars['ticks'][idx][i] if self.fvars['ticks'][idx][i] < 255 else 0
            exp._nonzero = nnz
            exp._data = <int*>sig_malloc(sizeof(int)*nnz*2)
            for j in range(2*nnz):
                exp._data[j] = self.fvars['exp_data'][idx][cum]
                cum += 1

            #Construct cyclotomic field coefficient
            d = self.fvars['coeff_denom'][idx][i]
            cyc_coeff = self.field([num / d for num in self.fvars['coeff_nums'][idx][i]])

            poly_tup.append((exp, cyc_coeff))
        ret = tuple(poly_tup)
        self.obj_cache[idx] = ret
        return ret

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
            sage: p = f.get_worker_pool()
            sage: for sextuple, fvar in f._shared_fvars.items():
            ....:     if sextuple == (f1, f1, f1, f1, f1, f1):
            ....:         f._tup_to_fpoly(fvar)
            ....:
            fx4
            sage: f.shutdown_worker_pool(p)
        """
        for sextuple in self.sext_to_idx:
            yield sextuple, self[sextuple]
