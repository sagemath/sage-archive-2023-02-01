from sage.misc.misc import prod
from sage.combinat.family import Family
from root_lattice_realization import RootLatticeRealization

class WeightLatticeRealization (RootLatticeRealization):

    def check(self):
        RootLatticeRealization.check(self)
        Lambda     = self.fundamental_weights()
        alphacheck = self.simple_coroots()

        for i in self.index_set():
            assert(Lambda[i].is_dominant())
            for j in self.index_set():
                assert(Lambda[j].scalar(alphacheck[i]) == (1 if i==j else 0))

        assert(self.rho().is_dominant())
        assert(self.highest_root().is_dominant())

    # Should this be a method or an attribute?
    # same question for the roots, ...
    def fundamental_weights(self):
        """
        Returns the family $(\Lambda_i)_{i\in I}$ of the fundamental weights
        """
        if not hasattr(self,"_fundamental_weights"):
            self._fundamental_weights = Family(self.index_set(),
                                               self.fundamental_weight,
                                               name = "Lambda")
        return self._fundamental_weights

    def rho(self):
        """
        EXAMPLES:
            sage: RootSystem(['A',3]).ambient_lattice().rho()
            (3, 2, 1, 0)
        """
        return sum(self.fundamental_weights())

    # Should it be a method of highest_weight?
    def weyl_dimension(self, highest_weight):
        """
        EXAMPLES:
            sage: RootSystem(['A',3]).ambient_lattice().weyl_dimension([2,1,0,0])
            20
        """
        highest_weight = self(highest_weight)
        assert(highest_weight.is_dominant())
        rho = self.rho()
        n = prod([(rho+highest_weight).dot_product(x) for x in self.positive_roots()])
        d = prod([ rho.dot_product(x) for x in self.positive_roots()])
        return n/d

