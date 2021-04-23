from multiprocessing import shared_memory

class KSHandler():
    def __init__(self, factory, names=None):
        n = len(factory._fvars)
        self.field = factory._field
        self.index = 0
        if names is None:
            self.ks = shared_memory.ShareableList([False]*n)
            self.nums = [shared_memory.ShareableList([int(0)]*self.field.degree()) for i in range(n)]
            self.denoms = shared_memory.ShareableList([int(0)]*n)
        else:
            self.ks = shared_memory.ShareableList(name=names[0])
            self.nums = [shared_memory.ShareableList(name=names[1][i]) for i in range(n)]
            self.denoms = shared_memory.ShareableList(name=names[2])

    def __getitem__(self,idx):
        if not self.ks[idx]:
            raise KeyError('Index {} does not correspond to a known square'.format(idx))
        denom = self.denoms[idx]
        return self.field([num / denom for num in self.nums[idx]])

    def __setitem__(self,idx,rhs):
        self.ks[idx] = True
        self.denoms[idx] = int(rhs.denominator())
        for i, c in enumerate(rhs._coefficients()):
            self.nums[idx][i] = int(c * self.denoms[idx])

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

        denom = self.denoms[self.index]
        ret = self.field([num / denom for num in self.nums[self.index]])
        self.index += 1
        return ret

    def __contains__(self,idx):
        return self.ks[idx]

    def __len__(self):
        return sum(self.ks)

    def items(self):
        for v in self:
            yield self.index-1, v

    def reset(self):
        n = len(self.ks)
        for i in range(n):
            self.ks[i] = False
            self.denoms[i] = 0
            for j in range(len(self.nums[i])):
                self.nums[i][j] = 0

    def unlink(self):
        self.ks.shm.unlink()
        for num_list in self.nums:
            num_list.shm.unlink()
        self.denoms.shm.unlink()

    def get_names(self):
        ks_name = self.ks.shm.name
        nums_names = [num_list.shm.name for num_list in self.nums]
        denom_name = self.denoms.shm.name
        return (ks_name,nums_names,denom_name)
