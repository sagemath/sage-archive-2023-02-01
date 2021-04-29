from cysignals.memory cimport sig_malloc
from sage.rings.polynomial.polydict cimport ETuple

from multiprocessing import shared_memory
import numpy as np
from sage.misc.cachefunc import cached_method

class KSHandler():
    def __init__(self, factory, name=None):
        self.field = factory._field
        self.index = 0
        n = len(factory._fvars)
        d = self.field.degree()
        ks_t = np.dtype([
            ('known', 'bool', (1,)),
            ('nums','i4',(d,)),
            ('denoms','u4',(d,))
        ])
        if name is None:
            self.shm = shared_memory.SharedMemory(create=True,size=n*ks_t.itemsize)
            self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)
            self.reset()
        else:
            self.shm = shared_memory.SharedMemory(name=name)
            self.ks_dat = np.ndarray((n,),dtype=ks_t,buffer=self.shm.buf)

    @cached_method
    def __getitem__(self, idx):
        if not self.ks_dat['known'][idx]:
            raise KeyError('Index {} does not correspond to a known square'.format(idx))
        cdef list rat = list()
        for i in range(self.field.degree()):
            rat.append(self.ks_dat['nums'][idx][i] / self.ks_dat['denoms'][idx][i])
        return self.field(rat)

    def __setitem__(self, idx, rhs):
        self.ks_dat['known'][idx] = True
        cdef int i
        for i in range(len(rhs)):
            self.ks_dat['nums'][idx][i] = -rhs[i].numerator()
            self.ks_dat['denoms'][idx][i] = rhs[i].denominator()

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index == self.ks_dat.size:
            raise StopIteration

        #Skip indices that are not known
        while not self.ks_dat['known'][self.index]:
            self.index += 1
            if self.index == self.ks_dat.size:
                raise StopIteration

        self.index += 1
        return self[self.index-1]

    def __contains__(self, idx):
        return self.ks_dat['known'][idx]

    def __len__(self):
        return sum(self.ks_dat['known'])

    def __eq__(self, other):
        return np.all(self.ks_dat == other.ks_dat)

    def items(self):
        for v in self:
            yield self.index-1, v

    def reset(self):
        n = self.ks_dat.size
        d = self.field.degree()
        self.__getitem__.clear_cache()
        self.ks_dat['known'] = np.zeros((n,1),dtype='bool')
        self.ks_dat['nums'] = np.zeros((n,d),dtype='i4')
        self.ks_dat['denoms'] = np.ones((n,d),dtype='u4')

class FvarsHandler():
    # cdef dict sext_to_idx
    # cdef int ngens
    # cdef field, fvars, fvars_t, shm

    def __init__(self, factory, name=None, max_terms=20):
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

    #Given a tuple of labels and a tuple of (ETuple, cyc_coeff) pairs, insert into
    #shared dictionary
    #WARNING: current data structure supports up to 2**16-1 entries,
    #each with each monomial in each entry having at most 254 nonzero terms.
    def __setitem__(self, sextuple, fvar):
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

    #Retrieve a record by unflattening and constructing relevant Python objects
    #Only the parent is allowed to modify the shared _fvars object, so to implement caching
    #we must tag modified objects and delete cache if present before retrieval
    def __getitem__(self, sextuple):
        if not sextuple in self.sext_to_idx:
            raise KeyError('Invalid sextuple {}'.format(sextuple))
        cdef int idx = self.sext_to_idx[sextuple]
        if idx in self.obj_cache:
            if self.fvars['modified'][idx]:
                del self.obj_cache[idx]
            else:
                return self.obj_cache[idx]
        cdef ETuple e = ETuple({}, self.ngens)
        cdef int cum, d, i, j
        cdef list poly_tup = list()
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
        cdef tuple ret = tuple(poly_tup)
        self.obj_cache[idx] = ret
        return ret

    def __len__(self):
        return self.fvars.size

    def __repr__(self):
        return str({sextuple: self[sextuple] for sextuple in self.sext_to_idx})

    def items(self):
        for sextuple in self.sext_to_idx:
            yield sextuple, self[sextuple]
