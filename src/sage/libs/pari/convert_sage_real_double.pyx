from cypari2.convert cimport new_gen_from_double

cpdef Gen new_gen_from_real_double_element(RealDoubleElement self):

    return new_gen_from_double(self._value)
