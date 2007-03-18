from finite_subgroup import FiniteSubgroup


class CuspidalSubgroup(FiniteSubgroup):
    def _repr_(self):
        return "Cuspidal subgroup of %s"%self.abelian_variety()

    def _generators(self):
        """
        Return a list of vectors that define elements of the rational
        homology that generate this finite subgroup.

        EXAMPLES:
            sage: J = J0(37)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, 0, 0, 1/3)]
            sage: J = J0(43)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, -1/7, 0, 1/7, 0, 2/7)]
            sage: J = J0(22)
            sage: C = J.cuspidal_subgroup()
            sage: C._generators()
            [(0, 0, 0, -1/5),
             (4/5, -1/5, 1/5, -2/5),
             (-1/5, 4/5, 1/5, -2/5),
             (-1/5, -1/5, 6/5, -2/5),
             (-1/5, -1/5, 1/5, 3/5),
             (-1/5, -1/5, 1/5, -1/5),
             (0, 0, 0, 1/5)]
            sage: J = J1(13)
            sage: C = J.cuspidal_subgroup()
            sage: len(C._generators())
            15
            sage: C._generators()[:3]
            [(0, 1/19, 1/19, -1/19), (-2/19, 0, 0, 1/19), (2/19, 0, 0, -1/19)]
        """
        A = self.abelian_variety()
        # Now take the image of the integral structure on the ambient
        # modular symbols space.
        M = A.modular_symbols()
        I = M.integral_period_mapping().matrix()
        Amb = M.ambient_module()
        J = Amb.integral_structure().basis_matrix()
        R = (J * I).rows()
        return [x for x in R if x.denominator() != 1]


class RationalCuspidalSubgroup(FiniteSubgroup):
    def _repr_(self):
        return "Rational cuspidal subgroup of %s"%self.abelian_variety()
