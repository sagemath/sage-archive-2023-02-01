cdef extern from "linbox/matrix/blas-matrix.h" namespace "LinBox":
    cdef cppclass BlasMatrix[T]:
        BlasMatrix(size_t nrows, size_t ncols)
        void setEntry(size_t i, size_t j, T)
        T &getEntry(size_t i, size_t j)

from sage.libs.linbox.modular cimport ModDoubleField, ModDoubleFieldElement, ModFloatField, ModFloatFieldElement

cdef extern from "linbox/algorithms/echelon-form.h":
    cdef cppclass EchelonFormDomainDouble "LinBox::EchelonFormDomain<LinBox::Modular<double> >":
        EchelonFormDomainDouble(ModDoubleField)
        int rowReducedEchelon(BlasMatrix[ModDoubleFieldElement], BlasMatrix[ModDoubleFieldElement])

    cdef cppclass EchelonFormDomainFloat "LinBox::EchelonFormDomain<LinBox::Modular<float> >":
        EchelonFormDomainFloat(ModFloatField)
        int rowReducedEchelon(BlasMatrix[ModFloatFieldElement], BlasMatrix[ModFloatFieldElement])

