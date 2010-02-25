cdef extern from "linbox/field/modular.h":
    # double

    cdef cppclass ModDoubleFieldElement "LinBox::Modular<double>::Element":
        pass

    cdef cppclass ModDoubleField "LinBox::Modular<double>":
        ModDoubleField(int modulus)
        ModDoubleFieldElement init(ModDoubleFieldElement res, int v)
        ModDoubleFieldElement init(ModDoubleFieldElement res, double v)
        ModDoubleFieldElement inv(  ModDoubleFieldElement x, ModDoubleFieldElement y)
        ModDoubleFieldElement neg(  ModDoubleFieldElement x, ModDoubleFieldElement y)
        ModDoubleFieldElement mul(  ModDoubleFieldElement r, ModDoubleFieldElement x, ModDoubleFieldElement y)
        ModDoubleFieldElement mulin(ModDoubleFieldElement x, ModDoubleFieldElement y)
        ModDoubleFieldElement addin(ModDoubleFieldElement x, ModDoubleFieldElement y)
        ModDoubleFieldElement invin(ModDoubleFieldElement y)
        ModDoubleFieldElement negin(ModDoubleFieldElement y)
        int characteristic(int c)
        bint isZero(ModDoubleFieldElement x)


    # float
    cdef cppclass ModFloatFieldElement "LinBox::Modular<float>::Element":
        pass

    cdef cppclass ModFloatField "LinBox::Modular<float>":
        ModFloatField(int modulus)
        ModFloatFieldElement init(ModFloatFieldElement res, int v)
        ModFloatFieldElement init(ModFloatFieldElement res, double v)
        ModFloatFieldElement inv(  ModFloatFieldElement x, ModFloatFieldElement y)
        ModFloatFieldElement neg(  ModFloatFieldElement x, ModFloatFieldElement y)
        ModFloatFieldElement mul(  ModFloatFieldElement r, ModFloatFieldElement x, ModFloatFieldElement y)
        ModFloatFieldElement mulin(ModFloatFieldElement x, ModFloatFieldElement y)
        ModFloatFieldElement addin(ModFloatFieldElement x, ModFloatFieldElement y)
        ModFloatFieldElement invin(ModFloatFieldElement y)
        ModFloatFieldElement negin(ModFloatFieldElement y)
        int characteristic(int c)
        bint isZero(ModFloatFieldElement x)


