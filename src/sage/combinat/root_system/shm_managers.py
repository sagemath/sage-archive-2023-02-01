from multiprocessing import shared_memory
import numpy as np
from sage.misc.cachefunc import cached_method

class KSHandler():
    def __init__(self, factory, names=None):
        self.field = factory._field
        self.index = 0
        n = len(factory._fvars)
        d = self.field.degree()
        coeff_dtype = np.dtype([('nums','i4',(d,)), ('denoms','u4',(d,))])
        if names is None:
            self.ks_shm = shared_memory.SharedMemory(create=True,size=n)
            self.ks = np.ndarray((n,),dtype='bool',buffer=self.ks_shm.buf)
            self.ks[:] = np.zeros((n,),dtype='bool')
            self.coeff_shm = shared_memory.SharedMemory(create=True,size=2*n*d*4)
            self.coeffs = np.ndarray((n,),dtype=coeff_dtype,buffer=self.coeff_shm.buf)
            self.coeffs['nums'][:] = np.zeros((n,d),dtype='i4')
            self.coeffs['denoms'][:] = np.ones((n,d),dtype='u4')
        else:
            self.ks_shm = shared_memory.SharedMemory(name=names[0])
            self.ks = np.ndarray((n,),dtype='bool',buffer=self.ks_shm.buf)
            self.coeff_shm = shared_memory.SharedMemory(name=names[1])
            self.coeffs = np.ndarray((n,),dtype=coeff_dtype,buffer=self.coeff_shm.buf)

    @cached_method
    def __getitem__(self, idx):
        if not self.ks[idx]:
            raise KeyError('Index {} does not correspond to a known square'.format(idx))
        rat = list()
        for i in range(self.field.degree()):
            rat.append(self.coeffs['nums'][idx][i] / self.coeffs['denoms'][idx][i])
        return self.field(rat)

    def __setitem__(self, idx, rhs):
        self.ks[idx] = True
        for i, c in enumerate(rhs):
            self.coeffs['nums'][idx][i] = -c.numerator()
            self.coeffs['denoms'][idx][i] = c.denominator()

    def __iter__(self):
        self.index = 0
        return self

    def __next__(self):
        if self.index == len(self.ks):
            raise StopIteration

        #Skip indices that are not known
        while not self.ks[self.index]:
            self.index += 1
            if self.index == len(self.ks):
                raise StopIteration

        self.index += 1
        return self[self.index-1]

    def __contains__(self, idx):
        return self.ks[idx]

    def __len__(self):
        return sum(self.ks)

    def __eq__(self, other):
        return np.all(self.ks == other.ks) and np.all(self.coeffs == other.coeffs)

    def items(self):
        for v in self:
            yield self.index-1, v

    def reset(self):
        n = len(self.ks)
        d = self.field.degree()
        self.ks[:] = np.zeros((n,),dtype='bool')
        self.__getitem__.clear_cache()
        self.coeffs['nums'][:] = np.zeros((n,d),dtype='i4')
        self.coeffs['denoms'][:] = np.ones((n,d),dtype='u4')

    def unlink(self):
        self.ks_shm.unlink()
        self.coeff_shm.unlink()

    def get_names(self):
        return (self.ks_shm.name,self.coeff_shm.name)
